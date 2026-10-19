// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Coupled Poisson + drift-diffusion LBM model of a 1-D p-n junction.

use super::{ChargeTransportLbm, PoissonLbm, constants, neg_gradient};
use crate::error::{Error, Result};

/// Default number of Poisson LBM relaxation sweeps per coupled step.
const DEFAULT_POISSON_SUBSTEPS: usize = 8;

/// Coupled Poisson + drift-diffusion LBM for an abrupt p-n junction.
///
/// The device occupies an `nx × 1` grid: acceptors `N_a` fill the cells
/// `x < nx/2` (p side) and donors `N_d` the cells `x ≥ nx/2` (n side).  Both
/// end cells are ohmic contacts: carrier densities take their local
/// charge-neutral equilibrium values and the potential is pinned to
/// `φ = V_T·ln(n/n_i)` (intrinsic level as reference), so that the contacts
/// differ by the built-in voltage.
///
/// Each [`SemiconductorLbm::step`]
/// 1. assembles ρ = q·(p − n + N_d⁺ − N_a⁻) from the current carriers,
/// 2. relaxes the Poisson LBM ([`SemiconductorLbm::poisson_substeps`] sweeps,
///    warm-started from the previous potential),
/// 3. feeds **E** = −∇φ to both carrier solvers, and
/// 4. advances electrons (charge sign −1) and holes (+1) by one common
///    physical time step: both relax at ω = 1 (positivity preserving); the
///    faster species uses the standard D2Q5 weights and the slower one the
///    reduced moving weight that matches the same Δt through the Einstein
///    relation.
///
/// The coupled iteration converges to the self-consistent equilibrium, in
/// which the interior potential drop equals the built-in voltage
/// V_bi = V_T · ln(N_a · N_d / n_i²), a depletion region forms around the
/// junction and n·p = n_i² holds everywhere.  The grid spacing must resolve
/// the Debye length (Δx ≲ L_D/3) so that the potential drop per cell stays of
/// order V_T.
pub struct SemiconductorLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Poisson solver.
    pub poisson: PoissonLbm,
    /// Electron transport.
    pub electrons: ChargeTransportLbm,
    /// Hole transport.
    pub holes: ChargeTransportLbm,
    /// Intrinsic carrier concentration (m⁻³).
    pub ni: f64,
    /// Thermal voltage (V).
    pub vt: f64,
    /// Acceptor doping density N_a (m⁻³) — left half of device.
    pub na: f64,
    /// Donor doping density N_d (m⁻³) — right half of device.
    pub nd: f64,
    /// Cell spacing (m).
    pub dx: f64,
    /// Poisson LBM relaxation sweeps per coupled step.
    pub poisson_substeps: usize,
}

impl SemiconductorLbm {
    /// Create a 1-D p-n junction modelled on an `nx × 1` grid.
    ///
    /// # Arguments
    /// * `nx` - number of grid cells (≥ 4)
    /// * `na` - acceptor doping (m⁻³)
    /// * `nd` - donor doping (m⁻³)
    /// * `ni` - intrinsic carrier density (m⁻³)
    /// * `temperature` - device temperature (K)
    /// * `mobility_e` - electron mobility (m²/V·s)
    /// * `mobility_h` - hole mobility (m²/V·s)
    /// * `dx` - cell spacing (m)
    ///
    /// # Errors
    /// [`Error::InvalidGridDimensions`] if `nx < 4`, and
    /// [`Error::ParameterOutOfBounds`] if a doping, the intrinsic density,
    /// the temperature, a mobility or `dx` is not positive and finite.
    pub fn new(
        nx: usize,
        na: f64,
        nd: f64,
        ni: f64,
        temperature: f64,
        mobility_e: f64,
        mobility_h: f64,
        dx: f64,
    ) -> Result<Self> {
        if nx < 4 {
            return Err(Error::InvalidGridDimensions { nx, ny: 1 });
        }
        let positive = [
            ("na", na),
            ("nd", nd),
            ("ni", ni),
            ("temperature", temperature),
            ("mobility_e", mobility_e),
            ("mobility_h", mobility_h),
            ("dx", dx),
        ];
        for (name, value) in positive {
            if !(value > 0.0 && value.is_finite()) {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: 0.0,
                    max: f64::MAX,
                });
            }
        }
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let poisson = PoissonLbm::new(nx, 1, constants::SILICON_EPSILON_R, 1.0).with_spacing(dx);
        let mut electrons =
            ChargeTransportLbm::new(nx, 1, mobility_e, temperature, 1.0, dx).with_charge_sign(-1.0);
        let mut holes =
            ChargeTransportLbm::new(nx, 1, mobility_h, temperature, 1.0, dx).with_charge_sign(1.0);
        // Common clock: the faster carrier has D_lu = 1/6 (ω = 1, w = 1/6);
        // the slower one keeps ω = 1 with a reduced moving weight.
        let dt = dx * dx / (6.0 * mobility_e.max(mobility_h) * vt);
        electrons.set_time_step(dt)?;
        holes.set_time_step(dt)?;
        Ok(Self {
            nx,
            ny: 1,
            poisson,
            electrons,
            holes,
            ni,
            vt,
            na,
            nd,
            dx,
            poisson_substeps: DEFAULT_POISSON_SUBSTEPS,
        })
    }

    /// Built-in potential V_bi = V_T · ln(N_a · N_d / n_i²).
    pub fn builtin_voltage(&self) -> f64 {
        self.vt * (self.na * self.nd / (self.ni * self.ni)).ln()
    }

    /// Physical duration of one coupled step (s).
    pub fn time_step(&self) -> f64 {
        self.electrons.time_step()
    }

    /// Net ionised doping N_d⁺ − N_a⁻ (m⁻³) of cell `x`.
    pub fn net_doping(&self, x: usize) -> f64 {
        if x < self.nx / 2 { -self.na } else { self.nd }
    }

    /// Charge-neutral equilibrium (n, p) for net doping `c`, computed without
    /// cancellation for the minority carrier.
    fn neutral_carriers(&self, c: f64) -> (f64, f64) {
        let ni2 = self.ni * self.ni;
        let root = (0.25 * c * c + ni2).sqrt();
        if c >= 0.0 {
            let n = 0.5 * c + root;
            (n, ni2 / n)
        } else {
            let p = -0.5 * c + root;
            (ni2 / p, p)
        }
    }

    /// Initialise the junction.
    ///
    /// Every cell starts at its local charge-neutral equilibrium (majority
    /// carrier ≈ doping, minority = n_i²/majority).  The two end cells become
    /// ohmic contacts whose potentials φ = V_T·ln(n/n_i) differ by the
    /// built-in voltage; since the initial charge vanishes, the initial
    /// potential is the Laplace solution between the contacts (a linear
    /// ramp).  The depletion region and its space charge then develop
    /// self-consistently during [`SemiconductorLbm::step`].
    pub fn init_junction(&mut self) {
        let last = self.nx - 1;
        let mut contact_phi = [0.0_f64; 2];
        for x in 0..self.nx {
            let (n, p) = self.neutral_carriers(self.net_doping(x));
            self.electrons.n[x] = n;
            self.holes.n[x] = p;
            if x == 0 || x == last {
                let phi = self.vt * (n / self.ni).ln();
                contact_phi[usize::from(x == last)] = phi;
                self.electrons.set_fixed_density(x, 0, n);
                self.holes.set_fixed_density(x, 0, p);
                self.poisson.set_phi(x, 0, phi);
            }
        }
        let [phi_p, phi_n] = contact_phi;
        for x in 0..self.nx {
            self.poisson.phi[x] = phi_p + (phi_n - phi_p) * x as f64 / last as f64;
        }
        self.electrons.init();
        self.holes.init();
        self.update_charge();
        self.poisson.init_from_phi();
        self.update_fields();
    }

    /// Assemble ρ = q·(p − n + N_d⁺ − N_a⁻) into the Poisson solver.
    fn update_charge(&mut self) {
        for x in 0..self.nx {
            let c = self.net_doping(x);
            self.poisson.rho_charge[x] =
                constants::ELEM_CHARGE * (self.holes.n[x] - self.electrons.n[x] + c);
        }
    }

    /// Feed **E** = −∇φ to both carrier solvers.
    fn update_fields(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        for y in 0..ny {
            for x in 0..nx {
                let (ex, ey) = neg_gradient(&self.poisson.phi, nx, ny, self.dx, x, y, false);
                let i = y * nx + x;
                self.electrons.ex[i] = ex;
                self.electrons.ey[i] = ey;
                self.holes.ex[i] = ex;
                self.holes.ey[i] = ey;
            }
        }
    }

    /// Run one coupled Poisson + transport step.
    pub fn step(&mut self) {
        self.update_charge();
        for _ in 0..self.poisson_substeps {
            self.poisson.step();
        }
        self.update_fields();
        self.electrons.step();
        self.holes.step();
    }

    /// Step until the coupled state is stationary or `max_steps` is reached.
    ///
    /// The residual of a step is the largest of `|Δφ|/V_T`, `|Δn|/n` and
    /// `|Δp|/p` over all cells.  Returns the number of steps taken and the
    /// last residual; the caller inspects the residual to decide convergence.
    pub fn solve_equilibrium(&mut self, max_steps: usize, tol: f64) -> (usize, f64) {
        let mut residual = f64::MAX;
        let mut steps = 0;
        let mut phi_old = self.poisson.phi.clone();
        let mut n_old = self.electrons.n.clone();
        let mut p_old = self.holes.n.clone();
        while steps < max_steps && residual > tol {
            phi_old.copy_from_slice(&self.poisson.phi);
            n_old.copy_from_slice(&self.electrons.n);
            p_old.copy_from_slice(&self.holes.n);
            self.step();
            residual = 0.0;
            for i in 0..self.nx * self.ny {
                let d_phi = (self.poisson.phi[i] - phi_old[i]).abs() / self.vt;
                let d_n = relative_change(self.electrons.n[i], n_old[i]);
                let d_p = relative_change(self.holes.n[i], p_old[i]);
                residual = residual.max(d_phi).max(d_n).max(d_p);
            }
            steps += 1;
        }
        (steps, residual)
    }
}

/// Relative change `|new − old|/|old|` (absolute change when `old` is zero).
fn relative_change(new: f64, old: f64) -> f64 {
    let diff = (new - old).abs();
    if old != 0.0 { diff / old.abs() } else { diff }
}
