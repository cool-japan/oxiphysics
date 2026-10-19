// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Two-fluid (electron + ion) electrostatic drift-diffusion plasma LBM.

use super::{ChargeTransportLbm, PoissonLbm, constants, neg_gradient};
use crate::error::{Error, Result};

/// Default number of Poisson LBM relaxation sweeps per coupled step.
const DEFAULT_POISSON_SUBSTEPS: usize = 16;

/// Two-fluid plasma LBM (electron fluid + ion fluid) in the drift-diffusion
/// limit.
///
/// Electrons (charge sign −1) and singly charged ions (+1) each evolve with a
/// D2Q5 drift-diffusion LBM ([`ChargeTransportLbm`]); both relax at ω = 1 and
/// advance on one common physical clock (the more mobile species uses the
/// standard D2Q5 weights, the other a reduced moving weight matching the same
/// Δt through the Einstein relation D = μ·V_T).
/// They are coupled through the self-consistent field **E** = −∇φ + **E**_ext
/// where φ solves ∇²φ = −e(n_i − n_e)/ε with the Poisson LBM.
///
/// Without Dirichlet cells in [`PlasmaLbm::poisson`] the domain is periodic;
/// the periodic Poisson problem is then solved for the charge relative to a
/// uniform neutralising background (the mean charge density is removed),
/// which is the standard closure for a periodic plasma.
pub struct PlasmaLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Electron fluid LBM (charge transport, charge sign −1).
    pub electrons: ChargeTransportLbm,
    /// Ion fluid LBM (charge transport, charge sign +1).
    pub ions: ChargeTransportLbm,
    /// Poisson solver for self-consistent field.
    pub poisson: PoissonLbm,
    /// Cell spacing (m).
    pub dx: f64,
    /// Uniform externally applied field (Ex, Ey) in V/m.
    pub external_field: [f64; 2],
    /// Poisson LBM relaxation sweeps per coupled step.
    pub poisson_substeps: usize,
}

impl PlasmaLbm {
    /// Create a new plasma LBM simulation.
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `mobility_e` - electron mobility (m²/V·s)
    /// * `mobility_i` - ion mobility (m²/V·s)
    /// * `temperature` - plasma temperature (K)
    /// * `dx` - cell spacing (m)
    /// * `epsilon_r` - relative permittivity of background medium
    ///
    /// # Errors
    /// [`Error::InvalidGridDimensions`] for an empty grid and
    /// [`Error::ParameterOutOfBounds`] if a mobility is negative, both are
    /// zero, or the temperature, spacing or permittivity is not positive.
    pub fn new(
        nx: usize,
        ny: usize,
        mobility_e: f64,
        mobility_i: f64,
        temperature: f64,
        dx: f64,
        epsilon_r: f64,
    ) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidGridDimensions { nx, ny });
        }
        let checks = [
            ("mobility_e", mobility_e, mobility_e >= 0.0),
            ("mobility_i", mobility_i, mobility_i >= 0.0),
            ("temperature", temperature, temperature > 0.0),
            ("dx", dx, dx > 0.0),
            ("epsilon_r", epsilon_r, epsilon_r > 0.0),
        ];
        for (name, value, ok) in checks {
            if !(ok && value.is_finite()) {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: 0.0,
                    max: f64::MAX,
                });
            }
        }
        let mu_max = mobility_e.max(mobility_i);
        if mu_max <= 0.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "mobility_e.max(mobility_i)".to_string(),
                value: mu_max,
                min: f64::MIN_POSITIVE,
                max: f64::MAX,
            });
        }
        let mut electrons = ChargeTransportLbm::new(nx, ny, mobility_e, temperature, 1.0, dx)
            .with_charge_sign(-1.0);
        let mut ions =
            ChargeTransportLbm::new(nx, ny, mobility_i, temperature, 1.0, dx).with_charge_sign(1.0);
        let vt = electrons.vt;
        // Common clock: the more mobile species has D_lu = 1/6 (ω = 1,
        // w = 1/6); the other keeps ω = 1 with a reduced moving weight.
        let dt = dx * dx / (6.0 * mu_max * vt);
        if mobility_e > 0.0 {
            electrons.set_time_step(dt)?;
        }
        if mobility_i > 0.0 {
            ions.set_time_step(dt)?;
        }
        let poisson = PoissonLbm::new(nx, ny, epsilon_r, 1.0).with_spacing(dx);
        Ok(Self {
            nx,
            ny,
            electrons,
            ions,
            poisson,
            dx,
            external_field: [0.0, 0.0],
            poisson_substeps: DEFAULT_POISSON_SUBSTEPS,
        })
    }

    /// Physical duration of one coupled step (s).
    pub fn time_step(&self) -> f64 {
        self.electrons.time_step().min(self.ions.time_step())
    }

    /// Assemble ρ = e·(n_i − n_e) into the Poisson solver (relative to the
    /// neutralising background when the domain is periodic).
    fn update_charge(&mut self) {
        let cells = self.nx * self.ny;
        for i in 0..cells {
            self.poisson.rho_charge[i] =
                constants::ELEM_CHARGE * (self.ions.n[i] - self.electrons.n[i]);
        }
        if !self.poisson.has_dirichlet() {
            let mean = self.poisson.rho_charge.iter().sum::<f64>() / cells as f64;
            for rho in &mut self.poisson.rho_charge {
                *rho -= mean;
            }
        }
    }

    /// Feed **E** = −∇φ + **E**_ext to both species.
    fn update_fields(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let periodic = !self.poisson.has_dirichlet();
        let [ex_ext, ey_ext] = self.external_field;
        for y in 0..ny {
            for x in 0..nx {
                let (ex, ey) = neg_gradient(&self.poisson.phi, nx, ny, self.dx, x, y, periodic);
                let i = y * nx + x;
                self.electrons.ex[i] = ex + ex_ext;
                self.electrons.ey[i] = ey + ey_ext;
                self.ions.ex[i] = ex + ex_ext;
                self.ions.ey[i] = ey + ey_ext;
            }
        }
    }

    /// One coupled step: charge assembly, Poisson relaxation, field update,
    /// then transport of both species.
    pub fn step(&mut self) {
        self.update_charge();
        for _ in 0..self.poisson_substeps {
            self.poisson.step();
        }
        self.update_fields();
        self.electrons.step();
        self.ions.step();
    }

    /// Run `n_steps` plasma steps.
    pub fn run(&mut self, n_steps: usize) {
        self.electrons.init();
        self.ions.init();
        for _ in 0..n_steps {
            self.step();
        }
    }
}
