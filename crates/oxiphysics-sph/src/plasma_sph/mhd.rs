// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoothed-particle magnetohydrodynamics: the [`MhdSph`] solver.
//!
//! # Discretisation
//!
//! Every sum runs over the neighbours `j` of particle `i`.  The pair gradient
//! is the symmetrised cubic-spline gradient
//!
//! `∇ᵢWᵢⱼ = ½ [W'(rᵢⱼ, hᵢ) + W'(rᵢⱼ, hⱼ)] (xᵢ − xⱼ)/rᵢⱼ`,
//!
//! with `W' ≤ 0`, so that it points from particle `i` towards particle `j`
//! (`∇ⱼWⱼᵢ = −∇ᵢWᵢⱼ`), and `Vⱼ = mⱼ/ρⱼ` is the neighbour volume.  With
//! `vᵢⱼ = vᵢ − vⱼ` the estimators are
//!
//! | quantity | SPH estimate |
//! |---|---|
//! | `∇×B` | `Σⱼ Vⱼ ∇ᵢWᵢⱼ × (Bⱼ − Bᵢ)` |
//! | `∇·B` | `(1/ρᵢ) Σⱼ mⱼ (Bⱼ − Bᵢ)·∇ᵢWᵢⱼ` |
//! | `DB/Dt = (B·∇)v − B(∇·v)` | `Σⱼ Vⱼ [Bᵢ (vᵢⱼ·∇ᵢWᵢⱼ) − vᵢⱼ (Bᵢ·∇ᵢWᵢⱼ)]` |
//! | `∇²B` | `2 Σⱼ Vⱼ (Bᵢ − Bⱼ) (rᵢⱼ·∇ᵢWᵢⱼ)/(rᵢⱼ² + 0.01 hᵢ²)` (Brookshaw) |
//! | `−∇P/ρ` | `−Σⱼ mⱼ (Pᵢ/ρᵢ² + Pⱼ/ρⱼ²) ∇ᵢWᵢⱼ` |
//! | `du/dt = −(P/ρ)∇·v` | `(Pᵢ/ρᵢ²) Σⱼ mⱼ vᵢⱼ·∇ᵢWᵢⱼ` |
//!
//! The momentum equation is `dv/dt = (∇×B)×B/(μ₀ρ) − ∇P/ρ` and the magnetic
//! field obeys the Lagrangian induction equation
//! `dB/dt = (B·∇)v − B(∇·v) + (η/μ₀)∇²B`: the particles move with the flow,
//! so the advection term of `∂B/∂t = ∇×(v×B)` is carried by the particle
//! motion and must not be counted a second time.  The gas is ideal,
//! `P = (γ − 1) ρ u`, and its specific internal energy `u` follows the
//! adiabatic compression/expansion work, so the symmetric pressure force and
//! the energy equation conserve `Σ mᵢ (½ vᵢ² + uᵢ)`.
//!
//! `∇·B = 0` is maintained with the hyperbolic/parabolic (mixed GLM) cleaning of
//! Dedner et al. (2002) in the SPH form of Tricco & Price (2012):
//!
//! `dB/dt = −∇ψ`, `dψ/dt = −c_h² ∇·B − ψ σ c_h/h`,
//!
//! with the gradient `∇ψᵢ = ρᵢ Σⱼ mⱼ (ψᵢ/ρᵢ² + ψⱼ/ρⱼ²) ∇ᵢWᵢⱼ`, which is the
//! discrete adjoint of the divergence above.
//!
//! # References
//! - Price (2012): Smoothed particle hydrodynamics and magnetohydrodynamics
//! - Tricco & Price (2012): Constrained hyperbolic divergence cleaning for SPMHD
//! - Dedner et al. (2002): Hyperbolic divergence cleaning for the MHD equations
//! - Brookshaw (1985): A method of calculating radiative heat diffusion in
//!   particle simulations

use super::{
    K_B, M_PROTON, MU_0, PlasmaParticle, add3, cross3, cubic_kernel, cubic_kernel_grad, dot3, len3,
    scale3, sub3,
};

/// Regularisation of the Brookshaw Laplacian denominator in units of `h²`.
const LAPLACIAN_EPS2: f64 = 0.01;

/// Default Dedner damping ratio σ (Tricco & Price 2012).
const DEFAULT_DEDNER_SIGMA: f64 = 0.8;

/// Default Dedner cleaning speed c_h \[m s⁻¹\].
const DEFAULT_DEDNER_CH: f64 = 1e4;

/// One interacting pair `(i, j)` of the particles of an [`MhdSph`].
struct Pair<'a> {
    /// The neighbour particle `j`.
    pj: &'a PlasmaParticle,
    /// Separation `x_i − x_j` \[m\].
    r_ij: [f64; 3],
    /// Symmetrised kernel gradient `∇_i W_ij` \[m⁻⁴\].
    grad_w: [f64; 3],
}

/// Symmetrised cubic-spline gradient `∇ᵢWᵢⱼ` for the separation `r_ij = xᵢ − xⱼ`.
///
/// It is the mean of the gradients with the smoothing lengths of both
/// particles, so it is exactly antisymmetric under `i ↔ j` and the pair forces
/// conserve momentum also for variable smoothing lengths.  It vanishes outside
/// the support radius `2 max(hᵢ, hⱼ)`.
fn pair_gradient(r_ij: [f64; 3], h_i: f64, h_j: f64) -> [f64; 3] {
    let grad_i = cubic_kernel_grad(r_ij, h_i);
    if h_i == h_j {
        return grad_i;
    }
    scale3(add3(grad_i, cubic_kernel_grad(r_ij, h_j)), 0.5)
}

/// Particle volume `m/ρ` \[m³\] (zero for a vanishing density).
fn volume(p: &PlasmaParticle) -> f64 {
    if p.density > 0.0 {
        p.mass / p.density
    } else {
        0.0
    }
}

/// Magnetohydrodynamics SPH solver.
///
/// Implements (see the [module documentation](self) for the estimators)
/// - the Lagrangian ideal and resistive induction equation
///   `dB/dt = (B·∇)v − B(∇·v) + (η/μ₀)∇²B`,
/// - the Lorentz force per unit mass `(1/(μ₀ρ)) (∇×B)×B`,
/// - the gas pressure force `−∇P/ρ` with `P = (γ − 1) ρ u` and the internal
///   energy equation `du/dt = −(P/ρ)∇·v`,
/// - Dedner hyperbolic/parabolic divergence cleaning for `∇·B = 0`,
/// - the explicit time integration of all of the above.
///
/// All smoothing-length dependent quantities use the per-particle
/// `smoothing_length`; the cleaning step is explicit, so it requires
/// `c_h · dt` to be well below the smoothing length (about `0.2 h`).
#[derive(Debug, Clone)]
pub struct MhdSph {
    /// List of plasma particles.
    pub particles: Vec<PlasmaParticle>,
    /// Resistivity η \[Ω m\]; the magnetic diffusivity is η/μ₀ \[m² s⁻¹\].
    pub resistivity: f64,
    /// Dedner cleaning speed c_h \[m s⁻¹\].
    pub dedner_ch: f64,
    /// Dedner damping ratio σ: the cleaning scalar ψ of a particle with
    /// smoothing length `h` decays on the time scale `h/(σ c_h)`.
    pub dedner_sigma: f64,
    /// Adiabatic index γ.
    pub gamma: f64,
}

impl MhdSph {
    /// Create a new MHD-SPH solver.
    pub fn new(resistivity: f64, gamma: f64) -> Self {
        Self {
            particles: Vec::new(),
            resistivity,
            dedner_ch: DEFAULT_DEDNER_CH,
            dedner_sigma: DEFAULT_DEDNER_SIGMA,
            gamma,
        }
    }

    /// Add a particle to the simulation.
    ///
    /// A particle that has no internal energy yet gets the specific internal
    /// energy `u = k_B T / ((γ − 1) m_p)` of its temperature, which gives the
    /// pressure `P = (γ − 1) ρ u = ρ k_B T / m_p` that
    /// [`PlasmaParticle::new`] reports.
    pub fn add_particle(&mut self, mut p: PlasmaParticle) {
        if p.internal_energy <= 0.0 && p.temperature > 0.0 && self.gamma > 1.0 {
            p.internal_energy = K_B * p.temperature / ((self.gamma - 1.0) * M_PROTON);
        }
        self.particles.push(p);
    }

    /// The neighbours `j ≠ i` of particle `i` that interact with it.
    ///
    /// Pairs outside the kernel support have a vanishing gradient, contribute
    /// nothing to any estimator and are skipped.  An out-of-range `i` has no
    /// neighbours.
    fn pairs(&self, i: usize) -> impl Iterator<Item = Pair<'_>> {
        let pi = self.particles.get(i);
        self.particles
            .iter()
            .enumerate()
            .filter_map(move |(j, pj)| {
                let pi = pi?;
                if j == i {
                    return None;
                }
                let r_ij = sub3(pi.position, pj.position);
                let grad_w = pair_gradient(r_ij, pi.smoothing_length, pj.smoothing_length);
                if grad_w == [0.0; 3] {
                    return None;
                }
                Some(Pair { pj, r_ij, grad_w })
            })
    }

    /// Recompute the thermodynamic state of every particle from the particle
    /// distribution.
    ///
    /// - density: the summation `ρᵢ = Σⱼ mⱼ W(rᵢⱼ, hᵢ)` over all particles
    ///   including `i` itself,
    /// - pressure: the ideal-gas law `Pᵢ = (γ − 1) ρᵢ uᵢ`,
    /// - temperature: `Tᵢ = (γ − 1) uᵢ m_p / k_B`.
    ///
    /// [`euler_step`](Self::euler_step) calls this first; call it after
    /// building or editing the particle set to evaluate the operators of this
    /// solver on a consistent state.
    pub fn update_thermodynamics(&mut self) {
        let gamma_m1 = (self.gamma - 1.0).max(0.0);
        let densities: Vec<f64> = self
            .particles
            .iter()
            .map(|pi| {
                self.particles
                    .iter()
                    .map(|pj| {
                        let r = len3(sub3(pi.position, pj.position));
                        pj.mass * cubic_kernel(r, pi.smoothing_length)
                    })
                    .sum()
            })
            .collect();
        for (p, rho) in self.particles.iter_mut().zip(densities) {
            p.density = rho;
            p.pressure = gamma_m1 * rho * p.internal_energy;
            if gamma_m1 > 0.0 {
                p.temperature = gamma_m1 * p.internal_energy * M_PROTON / K_B;
            }
        }
    }

    /// SPH curl of the magnetic field at particle `i`, `∇×B` \[T m⁻¹\]:
    /// `Σⱼ (mⱼ/ρⱼ) ∇ᵢWᵢⱼ × (Bⱼ − Bᵢ)`.
    ///
    /// For `B = (0, x, 0)` this is `(0, 0, 1)`.
    pub fn curl_b(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        let mut curl = [0.0_f64; 3];
        for pair in self.pairs(i) {
            let db = sub3(pair.pj.magnetic_field, pi.magnetic_field);
            curl = add3(curl, scale3(cross3(pair.grad_w, db), volume(pair.pj)));
        }
        curl
    }

    /// SPH divergence of the magnetic field at particle `i`, `∇·B` \[T m⁻¹\]:
    /// `(1/ρᵢ) Σⱼ mⱼ (Bⱼ − Bᵢ)·∇ᵢWᵢⱼ`.
    pub fn divergence_b(&self, i: usize) -> f64 {
        let Some(pi) = self.particles.get(i) else {
            return 0.0;
        };
        if pi.density <= 0.0 {
            return 0.0;
        }
        let sum: f64 = self
            .pairs(i)
            .map(|pair| {
                pair.pj.mass * dot3(sub3(pair.pj.magnetic_field, pi.magnetic_field), pair.grad_w)
            })
            .sum();
        sum / pi.density
    }

    /// Ideal induction equation at particle `i`: the Lagrangian rate of change
    /// `dB/dt = (B·∇)v − B(∇·v)` \[T s⁻¹\] of the field carried by the
    /// particle.
    ///
    /// SPH estimate `Σⱼ (mⱼ/ρⱼ) [Bᵢ (vᵢⱼ·∇ᵢWᵢⱼ) − vᵢⱼ (Bᵢ·∇ᵢWᵢⱼ)]` with
    /// `vᵢⱼ = vᵢ − vⱼ`.  It equals the Eulerian `∇×(v×B)` plus the advection
    /// `(v·∇)B`, which the particle motion supplies: a uniform field carried
    /// by a uniform flow is constant, a shear `v = (y, 0, 0)` of a field
    /// `B = (0, 1, 0)` winds it up at `dBₓ/dt = +1`, and a field perpendicular
    /// to a compression `v = α x̂` changes at `−α B` (flux freezing).
    pub fn induction_rhs(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        let mut db_dt = [0.0_f64; 3];
        for pair in self.pairs(i) {
            let v_ij = sub3(pi.velocity, pair.pj.velocity);
            let stretch = scale3(v_ij, dot3(pi.magnetic_field, pair.grad_w));
            let compress = scale3(pi.magnetic_field, dot3(v_ij, pair.grad_w));
            db_dt = add3(db_dt, scale3(sub3(compress, stretch), volume(pair.pj)));
        }
        db_dt
    }

    /// Lorentz acceleration of particle `i`, `a = (1/(μ₀ ρᵢ)) (∇×B)×B`
    /// \[m s⁻²\] (the force per unit mass).
    ///
    /// Uses the SPH curl [`curl_b`](Self::curl_b); for a field `B = (0, 1 + x, 0)`
    /// the magnetic pressure gradient pushes the plasma towards `−x`.
    pub fn lorentz_force(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        if pi.density <= 0.0 {
            return [0.0; 3];
        }
        let f = cross3(self.curl_b(i), pi.magnetic_field);
        scale3(f, 1.0 / (MU_0 * pi.density))
    }

    /// Gas pressure acceleration of particle `i`, `−∇P/ρ` \[m s⁻²\]:
    /// `−Σⱼ mⱼ (Pᵢ/ρᵢ² + Pⱼ/ρⱼ²) ∇ᵢWᵢⱼ` with the stored pressures and densities
    /// (see [`update_thermodynamics`](Self::update_thermodynamics)).
    pub fn pressure_acceleration(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        if pi.density <= 0.0 {
            return [0.0; 3];
        }
        let p_i = pi.pressure / (pi.density * pi.density);
        let mut accel = [0.0_f64; 3];
        for pair in self.pairs(i) {
            if pair.pj.density <= 0.0 {
                continue;
            }
            let p_j = pair.pj.pressure / (pair.pj.density * pair.pj.density);
            accel = sub3(accel, scale3(pair.grad_w, pair.pj.mass * (p_i + p_j)));
        }
        accel
    }

    /// Rate of change of the specific internal energy of particle `i`,
    /// `du/dt = −(P/ρ) ∇·v` \[J kg⁻¹ s⁻¹\]:
    /// `(Pᵢ/ρᵢ²) Σⱼ mⱼ vᵢⱼ·∇ᵢWᵢⱼ`.  Positive in compression.
    pub fn internal_energy_rate(&self, i: usize) -> f64 {
        let Some(pi) = self.particles.get(i) else {
            return 0.0;
        };
        if pi.density <= 0.0 {
            return 0.0;
        }
        let sum: f64 = self
            .pairs(i)
            .map(|pair| pair.pj.mass * dot3(sub3(pi.velocity, pair.pj.velocity), pair.grad_w))
            .sum();
        pi.pressure / (pi.density * pi.density) * sum
    }

    /// Gradient of the cleaning scalar at particle `i`,
    /// `∇ψᵢ = ρᵢ Σⱼ mⱼ (ψᵢ/ρᵢ² + ψⱼ/ρⱼ²) ∇ᵢWᵢⱼ` \[T m⁻¹\].
    fn psi_gradient(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        if pi.density <= 0.0 {
            return [0.0; 3];
        }
        let psi_i = pi.psi / (pi.density * pi.density);
        let mut grad = [0.0_f64; 3];
        for pair in self.pairs(i) {
            if pair.pj.density <= 0.0 {
                continue;
            }
            let psi_j = pair.pj.psi / (pair.pj.density * pair.pj.density);
            grad = add3(grad, scale3(pair.grad_w, pair.pj.mass * (psi_i + psi_j)));
        }
        scale3(grad, pi.density)
    }

    /// Apply one step of Dedner divergence cleaning (Tricco & Price 2012).
    ///
    /// The field is cleaned with the scalar `ψ` of the previous step,
    /// `B ← B − dt ∇ψ`, and `ψ` is then sourced by the divergence of the
    /// cleaned field and damped,
    /// `dψ/dt = −c_h² ∇·B − ψ/τ` with `τ = h/(σ c_h)`.  The damping is integrated
    /// exactly (`ψ ← ψ e^{−dt/τ} − c_h² ∇·B τ (1 − e^{−dt/τ})`), so it stays
    /// stable when `dt ≫ τ`; the hyperbolic coupling is explicit and needs
    /// `c_h · dt ≲ 0.2 h`.
    ///
    /// Uses the stored particle densities.  A non-positive or non-finite `dt`
    /// does nothing.
    pub fn dedner_cleaning_step(&mut self, dt: f64) {
        if !(dt.is_finite() && dt > 0.0) {
            return;
        }
        let n = self.particles.len();
        let grad_psi: Vec<[f64; 3]> = (0..n).map(|i| self.psi_gradient(i)).collect();
        for (p, grad) in self.particles.iter_mut().zip(grad_psi) {
            p.magnetic_field = sub3(p.magnetic_field, scale3(grad, dt));
        }
        let div_b: Vec<f64> = (0..n).map(|i| self.divergence_b(i)).collect();
        let c_h = self.dedner_ch;
        let sigma = self.dedner_sigma;
        for (p, div) in self.particles.iter_mut().zip(div_b) {
            let rate = if p.smoothing_length > 0.0 {
                sigma * c_h / p.smoothing_length
            } else {
                0.0
            };
            let x = rate * dt;
            let decay = (-x).exp();
            // (1 − e^{−x}) / x, the exact weight of a constant source over dt.
            let source_weight = if x > 1e-12 { -(-x).exp_m1() / x } else { 1.0 };
            p.psi = p.psi * decay - c_h * c_h * div * dt * source_weight;
        }
    }

    /// Compute the resistive diffusion term `(η/μ₀) ∇²B` \[T s⁻¹\] at particle `i`.
    ///
    /// Uses the Brookshaw Laplacian
    /// `2 Σⱼ (mⱼ/ρⱼ) (Bᵢ − Bⱼ) (rᵢⱼ·∇ᵢWᵢⱼ)/(rᵢⱼ² + 0.01 hᵢ²)`, which is exact for
    /// linear fields (zero) and gives `∇²Bₓ ≈ 2` for `Bₓ = x²`, and the magnetic
    /// diffusivity `η/μ₀` of the resistivity `η` \[Ω m\].
    pub fn resistive_diffusion(&self, i: usize) -> [f64; 3] {
        let Some(pi) = self.particles.get(i) else {
            return [0.0; 3];
        };
        let eps2 = LAPLACIAN_EPS2 * pi.smoothing_length * pi.smoothing_length;
        let mut lap_b = [0.0_f64; 3];
        for pair in self.pairs(i) {
            let r2 = dot3(pair.r_ij, pair.r_ij);
            let coeff = 2.0 * volume(pair.pj) * dot3(pair.r_ij, pair.grad_w) / (r2 + eps2);
            lap_b = add3(
                lap_b,
                scale3(sub3(pi.magnetic_field, pair.pj.magnetic_field), coeff),
            );
        }
        scale3(lap_b, self.resistivity / MU_0)
    }

    /// Advance the simulation by one explicit Euler step of duration `dt`.
    ///
    /// 1. [`update_thermodynamics`](Self::update_thermodynamics): density
    ///    summation and the ideal-gas pressure.
    /// 2. All rates are evaluated from that one state: the acceleration
    ///    `(∇×B)×B/(μ₀ρ) − ∇P/ρ`, the field rate `dB/dt` of the induction
    ///    equation plus resistive diffusion, and `du/dt`.
    /// 3. Velocity and magnetic field are advanced.
    /// 4. The internal energy is advanced with `du/dt` evaluated for the
    ///    mean of the old and the new velocity.  Together with the symmetric
    ///    pressure force this makes the work done by the pressure exactly equal
    ///    to the heat gained by the gas, so a pure gas conserves
    ///    `Σ mᵢ (½ vᵢ² + uᵢ)` to rounding error.
    /// 5. [`dedner_cleaning_step`](Self::dedner_cleaning_step) cleans the
    ///    advanced field on the not yet moved particle distribution.
    /// 6. The positions are advanced with the new velocities.
    ///
    /// The `accel` of every particle holds the acceleration of this step; the
    /// `density` describes the distribution at the beginning of the step, and
    /// `pressure` and `temperature` the gas state after the energy update.
    pub fn euler_step(&mut self, dt: f64) {
        self.update_thermodynamics();
        let rates: Vec<([f64; 3], [f64; 3], f64)> = (0..self.particles.len())
            .map(|i| {
                (
                    add3(self.lorentz_force(i), self.pressure_acceleration(i)),
                    add3(self.induction_rhs(i), self.resistive_diffusion(i)),
                    self.internal_energy_rate(i),
                )
            })
            .collect();

        let mut heating_old = Vec::with_capacity(rates.len());
        for (p, (accel, db_dt, du_dt)) in self.particles.iter_mut().zip(rates) {
            p.accel = accel;
            p.velocity = add3(p.velocity, scale3(accel, dt));
            p.magnetic_field = add3(p.magnetic_field, scale3(db_dt, dt));
            heating_old.push(du_dt);
        }

        let heating_new: Vec<f64> = (0..self.particles.len())
            .map(|i| self.internal_energy_rate(i))
            .collect();
        let gamma_m1 = (self.gamma - 1.0).max(0.0);
        for (p, (du_old, du_new)) in self
            .particles
            .iter_mut()
            .zip(heating_old.into_iter().zip(heating_new))
        {
            p.internal_energy = (p.internal_energy + 0.5 * (du_old + du_new) * dt).max(0.0);
            p.pressure = gamma_m1 * p.density * p.internal_energy;
            if gamma_m1 > 0.0 {
                p.temperature = gamma_m1 * p.internal_energy * M_PROTON / K_B;
            }
        }

        self.dedner_cleaning_step(dt);

        for p in self.particles.iter_mut() {
            p.position = add3(p.position, scale3(p.velocity, dt));
        }
    }

    /// Total kinetic energy of all particles \[J\].
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Total thermal (internal) energy Σ mᵢ uᵢ \[J\].
    pub fn total_thermal_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.mass * p.internal_energy)
            .sum()
    }

    /// Total magnetic energy Σ B²/(2μ₀) V_i \[J\].
    pub fn total_magnetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let b2 = dot3(p.magnetic_field, p.magnetic_field);
                b2 / (2.0 * MU_0) * volume(p)
            })
            .sum()
    }
}
