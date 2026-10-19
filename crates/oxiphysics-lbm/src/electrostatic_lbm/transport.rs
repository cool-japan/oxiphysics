// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! D2Q5 lattice-Boltzmann drift-diffusion solver for charge carriers.

use super::{
    MAX_MOVING_WEIGHT, STANDARD_MOVING_WEIGHT, constants, drift_equilibrium, stream_periodic,
    weighted_lattice_diffusivity,
};
use crate::error::{Error, Result};

/// LBM solver for the 2-D drift-diffusion equation of one carrier species.
///
/// Solves ∂n/∂t + ∇·(n **v** − D ∇n) = S with the drift velocity
/// **v** = s·μ·**E** (s = [`ChargeTransportLbm::charge_sign`]: −1 for
/// electrons, +1 for holes / positive ions) and the Einstein relation
/// D = μ·V_T.
///
/// The BGK rate ω and the weight `w` of each moving D2Q5 population
/// ([`ChargeTransportLbm::lattice_weight`], 1/6 by default) fix the lattice
/// diffusivity `D_lu = 2w·(1/ω − ½)`; the Einstein relation then fixes the
/// physical time step
/// `Δt = D_lu·Δx²/(μ·V_T)` (see [`ChargeTransportLbm::time_step`]) and the
/// lattice drift velocity `u = v·Δt/Δx = s·E·Δx·D_lu/V_T`.
///
/// The drift enters through an exponentially fitted equilibrium (moving
/// populations built from `e^{±a/2}`, see the module documentation):
/// for weak fields it is the usual `w_q·n·(1 + 3 c_q·u)`, and for strong
/// fields its zero-current steady state still reproduces the Boltzmann ratio
/// `n ∝ exp(−s·φ/V_T)` exactly, so minority carriers stay accurate across
/// potential drops of several V_T per cell.  Drops beyond ≈ 2.6 V_T per cell
/// (for w = 1/6) are limited to keep every population non-negative.  At
/// ω = 1 the post-collision state is that non-negative equilibrium, so the
/// scheme preserves positivity unconditionally; ω > 1 (over-relaxation) can
/// produce transient negative populations next to steep density jumps.
///
/// Streaming is periodic; cells pinned with
/// [`ChargeTransportLbm::set_fixed_density`] act as Dirichlet (contact) nodes.
pub struct ChargeTransportLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions f\[y * nx * 5 + x * 5 + q\].
    pub f: Vec<f64>,
    /// Carrier density n at each cell.
    pub n: Vec<f64>,
    /// Carrier mobility μ (m²/V·s), magnitude.
    pub mobility: f64,
    /// Thermal voltage V_T = kT/q (V).
    pub vt: f64,
    /// BGK relaxation rate.
    pub omega: f64,
    /// Weight `w` of each of the four moving D2Q5 populations (rest weight
    /// `1 − 4w`), in (0, 1/4]; 1/6 is the standard lattice.
    pub lattice_weight: f64,
    /// Ex field (V/m).
    pub ex: Vec<f64>,
    /// Ey field (V/m).
    pub ey: Vec<f64>,
    /// Source/recombination term: carrier density added per lattice time
    /// step (a physical rate R in m⁻³·s⁻¹ corresponds to
    /// `R·`[`ChargeTransportLbm::time_step`]).
    pub source: Vec<f64>,
    /// Cell spacing (m).
    pub dx: f64,
    /// Sign of the carrier charge: −1 for electrons, +1 for holes/ions.
    pub charge_sign: f64,
    /// Dirichlet constraints: `Some(n)` pins the carrier density of a cell.
    pub fixed: Vec<Option<f64>>,
    /// Streaming buffer reused between steps.
    scratch: Vec<f64>,
}

impl ChargeTransportLbm {
    /// Create a new charge-transport solver for positive carriers
    /// (`charge_sign = +1`; use [`ChargeTransportLbm::with_charge_sign`] for
    /// electrons).
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `mobility` - carrier mobility μ (m²/V·s)
    /// * `temperature` - lattice temperature (K)
    /// * `omega` - BGK relaxation rate
    /// * `dx` - cell spacing (m)
    pub fn new(nx: usize, ny: usize, mobility: f64, temperature: f64, omega: f64, dx: f64) -> Self {
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let n = nx * ny;
        Self {
            nx,
            ny,
            f: vec![0.0; n * 5],
            n: vec![0.0; n],
            mobility,
            vt,
            omega,
            lattice_weight: STANDARD_MOVING_WEIGHT,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            source: vec![0.0; n],
            dx,
            charge_sign: 1.0,
            fixed: vec![None; n],
            scratch: vec![0.0; n * 5],
        }
    }

    /// Set the sign of the carrier charge (−1 electrons, +1 holes/ions).
    pub fn with_charge_sign(mut self, sign: f64) -> Self {
        self.charge_sign = sign;
        self
    }

    /// Check that the physical parameters define a valid lattice mapping.
    ///
    /// # Errors
    /// [`Error::InvalidGridDimensions`] for an empty grid and
    /// [`Error::ParameterOutOfBounds`] for a non-physical mobility, thermal
    /// voltage, cell spacing, BGK rate or lattice weight.
    pub fn validate(&self) -> Result<()> {
        if self.nx == 0 || self.ny == 0 {
            return Err(Error::InvalidGridDimensions {
                nx: self.nx,
                ny: self.ny,
            });
        }
        let checks = [
            (
                "mobility",
                self.mobility,
                0.0,
                f64::MAX,
                self.mobility >= 0.0,
            ),
            ("thermal_voltage", self.vt, 0.0, f64::MAX, self.vt > 0.0),
            ("dx", self.dx, 0.0, f64::MAX, self.dx > 0.0),
            (
                "omega",
                self.omega,
                0.0,
                2.0,
                self.omega > 0.0 && self.omega < 2.0,
            ),
            (
                "lattice_weight",
                self.lattice_weight,
                0.0,
                MAX_MOVING_WEIGHT,
                self.lattice_weight > 0.0 && self.lattice_weight <= MAX_MOVING_WEIGHT,
            ),
        ];
        for (name, value, min, max, ok) in checks {
            if !(ok && value.is_finite()) {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min,
                    max,
                });
            }
        }
        Ok(())
    }

    /// Diffusion coefficient via Einstein relation D = μ·V_T.
    pub fn diffusion_coeff(&self) -> f64 {
        self.mobility * self.vt
    }

    /// Lattice diffusivity `D_lu = 2w·(1/ω − ½)`.
    pub fn lattice_diffusivity(&self) -> f64 {
        weighted_lattice_diffusivity(self.omega, self.lattice_weight)
    }

    /// Physical duration of one lattice step, `Δt = D_lu·Δx²/(μ·V_T)` (s).
    ///
    /// Returns `f64::INFINITY` for immobile carriers (μ·V_T ≤ 0).
    pub fn time_step(&self) -> f64 {
        let d = self.diffusion_coeff();
        if d > 0.0 {
            self.lattice_diffusivity() * self.dx * self.dx / d
        } else {
            f64::INFINITY
        }
    }

    /// Choose the lattice weight so that one lattice step lasts `dt` seconds
    /// at the current ω.
    ///
    /// Used to advance several species with different mobilities on a common
    /// physical clock without changing their relaxation rate (in particular
    /// keeping the positivity-preserving ω = 1).
    ///
    /// # Errors
    /// [`Error::ParameterOutOfBounds`] if `dt`, the mobility, the thermal
    /// voltage or the spacing is not positive and finite, if ω ∉ (0, 2), or if
    /// the required weight exceeds 1/4 (time step too long for this ω).
    pub fn set_time_step(&mut self, dt: f64) -> Result<()> {
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(Error::ParameterOutOfBounds {
                name: "dt".to_string(),
                value: dt,
                min: 0.0,
                max: f64::MAX,
            });
        }
        let d = self.diffusion_coeff();
        if !(d > 0.0 && d.is_finite()) {
            return Err(Error::ParameterOutOfBounds {
                name: "diffusion_coeff".to_string(),
                value: d,
                min: 0.0,
                max: f64::MAX,
            });
        }
        if !(self.dx > 0.0 && self.dx.is_finite()) {
            return Err(Error::ParameterOutOfBounds {
                name: "dx".to_string(),
                value: self.dx,
                min: 0.0,
                max: f64::MAX,
            });
        }
        if !(self.omega > 0.0 && self.omega < 2.0) {
            return Err(Error::ParameterOutOfBounds {
                name: "omega".to_string(),
                value: self.omega,
                min: 0.0,
                max: 2.0,
            });
        }
        let d_lu = d * dt / (self.dx * self.dx);
        let weight = d_lu / (2.0 * (1.0 / self.omega - 0.5));
        if !(weight > 0.0 && weight <= MAX_MOVING_WEIGHT) {
            return Err(Error::ParameterOutOfBounds {
                name: "lattice_weight".to_string(),
                value: weight,
                min: 0.0,
                max: MAX_MOVING_WEIGHT,
            });
        }
        self.lattice_weight = weight;
        Ok(())
    }

    /// Lattice drift velocity `u = s·E·Δx·D_lu/V_T` for the field (ex, ey)
    /// (weak-field limit of the fitted equilibrium).
    pub fn lattice_drift(&self, ex: f64, ey: f64) -> (f64, f64) {
        let k = self.charge_sign * self.dx * self.lattice_diffusivity() / self.vt;
        (k * ex, k * ey)
    }

    /// Set carrier density at cell (x, y).
    pub fn set_density(&mut self, x: usize, y: usize, val: f64) {
        let i = y * self.nx + x;
        self.n[i] = val;
    }

    /// Pin the carrier density of cell (x, y) to `val` (Dirichlet/contact).
    pub fn set_fixed_density(&mut self, x: usize, y: usize, val: f64) {
        let i = y * self.nx + x;
        self.n[i] = val;
        self.fixed[i] = Some(val);
        let feq = drift_equilibrium(val, 0.0, 0.0, self.omega, self.lattice_weight);
        self.f[i * 5..i * 5 + 5].copy_from_slice(&feq);
    }

    /// Remove the Dirichlet condition of cell (x, y).
    pub fn release_density(&mut self, x: usize, y: usize) {
        let i = y * self.nx + x;
        self.fixed[i] = None;
    }

    /// Set electric field components at cell (x, y).
    pub fn set_field(&mut self, x: usize, y: usize, ex: f64, ey: f64) {
        let i = y * self.nx + x;
        self.ex[i] = ex;
        self.ey[i] = ey;
    }

    /// Initialise distributions from density.
    pub fn init(&mut self) {
        for i in 0..self.nx * self.ny {
            if let Some(val) = self.fixed[i] {
                self.n[i] = val;
            }
            let feq = drift_equilibrium(self.n[i], 0.0, 0.0, self.omega, self.lattice_weight);
            self.f[i * 5..i * 5 + 5].copy_from_slice(&feq);
        }
    }

    /// Perform one BGK + streaming step.
    ///
    /// After the step [`ChargeTransportLbm::n`] holds the density of the
    /// streamed state.  Immobile carriers (μ = 0, hence D = 0 by the Einstein
    /// relation) neither drift nor diffuse; only the source is applied.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let cells = nx * ny;
        if self.diffusion_coeff() <= 0.0 {
            for i in 0..cells {
                self.n[i] = match self.fixed[i] {
                    Some(val) => val,
                    None => self.n[i] + self.source[i],
                };
                let feq = drift_equilibrium(self.n[i], 0.0, 0.0, self.omega, self.lattice_weight);
                self.f[i * 5..i * 5 + 5].copy_from_slice(&feq);
            }
            return;
        }
        let omega = self.omega;
        let weight = self.lattice_weight;
        // drift strength per cell in thermal units: a = s·E·Δx/V_T
        let k = self.charge_sign * self.dx / self.vt;

        // collision
        for i in 0..cells {
            let ax = k * self.ex[i];
            let ay = k * self.ey[i];
            let base = i * 5;
            match self.fixed[i] {
                Some(val) => {
                    self.n[i] = val;
                    let feq = drift_equilibrium(val, ax, ay, omega, weight);
                    self.f[base..base + 5].copy_from_slice(&feq);
                }
                None => {
                    let n_loc: f64 = self.f[base..base + 5].iter().sum();
                    self.n[i] = n_loc;
                    let feq = drift_equilibrium(n_loc, ax, ay, omega, weight);
                    let src = drift_equilibrium(self.source[i], 0.0, 0.0, omega, weight);
                    for q in 0..5_usize {
                        self.f[base + q] += omega * (feq[q] - self.f[base + q]) + src[q];
                    }
                }
            }
        }

        // streaming
        if self.scratch.len() != self.f.len() {
            self.scratch = vec![0.0; self.f.len()];
        }
        stream_periodic(&self.f, nx, ny, &mut self.scratch);
        std::mem::swap(&mut self.f, &mut self.scratch);

        // density of the streamed state
        for i in 0..cells {
            self.n[i] = match self.fixed[i] {
                Some(val) => val,
                None => self.f[i * 5..i * 5 + 5].iter().sum(),
            };
        }
    }

    /// Run `n_steps` drift-diffusion steps.
    pub fn run(&mut self, n_steps: usize) {
        self.init();
        for _ in 0..n_steps {
            self.step();
        }
    }
}
