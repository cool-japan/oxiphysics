// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Two-fluid (Landau) model of superfluid helium: a lattice-Boltzmann normal
//! fluid, an inviscid superfluid and the mutual friction that couples them.
//!
//! See [`SuperfluidLbm`].

use super::superfluid_fraction;
use crate::error::Error;

/// Number of populations of the D1Q3 lattice of the normal fluid.
const Q: usize = 3;

/// Courant number `(|v| + c_s) h / Δx` of the explicit superfluid sub-steps.
const EULER_CFL: f64 = 0.4;

/// Upper bound on the number of explicit superfluid sub-steps per lattice step.
const MAX_EULER_SUBSTEPS: usize = 4096;

/// Densities at or below this value are treated as vacuum.
const RHO_FLOOR: f64 = f64::MIN_POSITIVE;

/// Largest ratio `|v|/c` for which the rest population of the D1Q3 equilibrium,
/// `(2/3) ρ (1 − 3 v²/(2 c²))`, stays positive: `sqrt(2/3)`.
const MAX_NORMAL_SPEED_RATIO: f64 = 0.816_496_580_927_726;

/// D1Q3 equilibrium populations for density `rho` and velocity `u`.
///
/// The lattice velocities are `e = [0, +c, -c]` with weights `[2/3, 1/6, 1/6]`
/// and `c_s² = c²/3`, and
///
/// `f_q = w_q ρ [1 + e_q u/c_s² + (e_q u)²/(2 c_s⁴) − u²/(2 c_s²)]`,
///
/// which reproduces `Σ f_q = ρ`, `Σ e_q f_q = ρ u` and `Σ e_q² f_q = ρ (c_s² + u²)`.
fn d1q3_equilibrium(rho: f64, u: f64, c: f64) -> [f64; Q] {
    let m = u / c;
    let quad = 3.0 * m * m;
    [
        2.0 / 3.0 * rho * (1.0 - 0.5 * quad),
        rho / 6.0 * (1.0 + 3.0 * m + quad),
        rho / 6.0 * (1.0 - 3.0 * m + quad),
    ]
}

/// Minmod slope limiter: the argument of smaller magnitude if both have the
/// same sign, zero otherwise.
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// Velocity `mom / rho` of an isothermal state (zero in vacuum).
fn velocity_of(rho: f64, mom: f64) -> f64 {
    if rho > RHO_FLOOR { mom / rho } else { 0.0 }
}

/// Work arrays of the finite-volume evaluation of the isothermal Euler
/// equations, allocated once per lattice step.
struct EulerScratch {
    velocity: Vec<f64>,
    slope_rho: Vec<f64>,
    slope_velocity: Vec<f64>,
    flux: Vec<[f64; 2]>,
}

impl EulerScratch {
    fn new(n: usize) -> Self {
        Self {
            velocity: vec![0.0; n],
            slope_rho: vec![0.0; n],
            slope_velocity: vec![0.0; n],
            flux: vec![[0.0; 2]; n],
        }
    }

    /// Write `−∂ₓF` of the isothermal Euler equations
    ///
    /// `∂ρ/∂t + ∂ₓ(ρ v) = 0`, `∂(ρ v)/∂t + ∂ₓ(ρ v² + c_s² ρ) = 0`
    ///
    /// for the periodic cell averages `(rho, mom)` into `(d_rho, d_mom)`.
    ///
    /// The interface states are the MUSCL reconstruction of the primitive
    /// variables `(ρ, v)` with minmod-limited slopes and the interface flux is
    /// the local Lax-Friedrichs (Rusanov) flux with the maximum signal speed
    /// `|v| + c_s`.  The scheme is conservative, so the total mass and
    /// momentum change only through the update itself.
    fn rhs(
        &mut self,
        rho: &[f64],
        mom: &[f64],
        cs: f64,
        dx: f64,
        d_rho: &mut [f64],
        d_mom: &mut [f64],
    ) {
        let Self {
            velocity,
            slope_rho,
            slope_velocity,
            flux,
        } = self;
        let n = rho.len();
        let cs2 = cs * cs;
        for (v, (&r, &m)) in velocity.iter_mut().zip(rho.iter().zip(mom)) {
            *v = velocity_of(r, m);
        }
        for i in 0..n {
            let west = (i + n - 1) % n;
            let east = (i + 1) % n;
            slope_rho[i] = minmod(rho[i] - rho[west], rho[east] - rho[i]);
            slope_velocity[i] = minmod(velocity[i] - velocity[west], velocity[east] - velocity[i]);
        }
        // Flux through the interface between cell `i` and cell `i + 1`.
        for i in 0..n {
            let east = (i + 1) % n;
            let rho_l = rho[i] + 0.5 * slope_rho[i];
            let v_l = velocity[i] + 0.5 * slope_velocity[i];
            let rho_r = rho[east] - 0.5 * slope_rho[east];
            let v_r = velocity[east] - 0.5 * slope_velocity[east];
            let speed = v_l.abs().max(v_r.abs()) + cs;
            let flux_l = [rho_l * v_l, rho_l * v_l * v_l + cs2 * rho_l];
            let flux_r = [rho_r * v_r, rho_r * v_r * v_r + cs2 * rho_r];
            flux[i] = [
                0.5 * (flux_l[0] + flux_r[0]) - 0.5 * speed * (rho_r - rho_l),
                0.5 * (flux_l[1] + flux_r[1]) - 0.5 * speed * (rho_r * v_r - rho_l * v_l),
            ];
        }
        for i in 0..n {
            let west = (i + n - 1) % n;
            d_rho[i] = -(flux[i][0] - flux[west][0]) / dx;
            d_mom[i] = -(flux[i][1] - flux[west][1]) / dx;
        }
    }
}

/// Two-fluid model of He-II on a periodic one-dimensional lattice: a normal
/// fluid and a superfluid component coupled by mutual friction.
///
/// Landau's decomposition splits the total density into `ρ = ρ_s + ρ_n`; the
/// superfluid fraction `ρ_s/ρ` approaches 1 at `T = 0` and vanishes at the
/// critical temperature.  Both components are isothermal fluids with the
/// lattice sound speed `c_s = (Δx/Δt)/√3` (partial pressures `p_i = c_s² ρ_i`),
/// and each one has its own continuity and momentum equation:
///
/// * The **normal fluid** is a D1Q3 BGK lattice Boltzmann fluid.  Its
///   populations are collided with relaxation time `tau_n` (in units of Δt,
///   `tau_n > 1/2`) and streamed one site per step, which gives the kinematic
///   viscosity `ν_n = c_s² (τ_n − ½) Δt` ([`normal_viscosity`](Self::normal_viscosity)).
///   The fields `rho_n` and `vel_n` are the density and velocity moments of the
///   populations; only the non-equilibrium part of the populations (which has
///   neither mass nor momentum) is stored besides them, so the fields stay the
///   single source of truth for the state.
/// * The **superfluid** is inviscid and irrotational.  Its mass and momentum
///   obey the isothermal Euler equations
///   `∂ρ_s/∂t + ∂ₓ(ρ_s v_s) = 0`, `∂(ρ_s v_s)/∂t + ∂ₓ(ρ_s v_s² + c_s² ρ_s) = 0`,
///   which are advanced on the lattice spacing Δx by a conservative
///   second-order finite-volume scheme (MUSCL reconstruction with the minmod
///   limiter, Rusanov flux, SSP-RK2 time stepping, sub-cycled to the Courant
///   limit).
/// * **Mutual friction** exchanges momentum between the components,
///   `ρ_n ∂v_n/∂t = F`, `ρ_s ∂v_s/∂t = −F` with `F = α ρ_s (v_s − v_n)`.
///   The slip velocity `v_s − v_n` therefore relaxes exponentially,
///   `∂(v_s − v_n)/∂t = −κ (v_s − v_n)` with `κ = α ρ / ρ_n`
///   ([`slip_decay_rate`](Self::slip_decay_rate)), while the momentum density
///   `ρ_s v_s + ρ_n v_n` is conserved.
///
/// One [`step`](Self::step) is a Strang splitting: half a step of mutual
/// friction (integrated exactly, so it is unconditionally stable and conserves
/// momentum to rounding error), the transport step of both components (normal
/// fluid collision and streaming, superfluid Euler update), and the second
/// half step of mutual friction.  The temperature only enters through the
/// initial density split of [`new`](Self::new); the model has no thermal
/// (second-sound) dynamics.
///
/// All derived quantities use lattice-consistent units: `Δx` and `Δt` set the
/// lattice speed `c = Δx/Δt`, and the velocities, densities and `α` are in the
/// same units as `Δx/Δt` and the initial densities.
#[derive(Debug, Clone)]
pub struct SuperfluidLbm {
    /// Number of lattice sites.
    pub nx: usize,
    /// Superfluid density ρ_s(x).
    pub rho_s: Vec<f64>,
    /// Normal fluid density ρ_n(x).
    pub rho_n: Vec<f64>,
    /// Superfluid velocity field v_s(x).
    pub vel_s: Vec<f64>,
    /// Normal fluid velocity field v_n(x).
    pub vel_n: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
    /// Critical temperature T_c (K).
    pub t_critical: f64,
    /// Mutual friction coefficient α (the slip velocity decays at the rate
    /// `α ρ/ρ_n`).
    pub alpha: f64,
    /// BGK relaxation time of the normal component, in units of Δt (> 1/2).
    pub tau_n: f64,
    /// Lattice spacing Δx.
    pub dx: f64,
    /// Time step Δt.
    pub dt: f64,
    /// Non-equilibrium part `f_q − f_q^eq(ρ_n, v_n)` of the normal-fluid
    /// populations (no mass, no momentum).
    pub(super) f_neq: Vec<[f64; Q]>,
}

impl SuperfluidLbm {
    /// Create a two-fluid LBM lattice at rest.
    ///
    /// The unit total density is split into `ρ_s = f_s(T)` and
    /// `ρ_n = 1 − f_s(T)` with the superfluid fraction of
    /// [`superfluid_fraction`].
    ///
    /// # Arguments
    /// * `nx`         – number of lattice sites
    /// * `temperature`– current temperature (K)
    /// * `t_critical` – critical temperature (K)
    /// * `alpha`      – mutual friction coefficient
    /// * `tau_n`      – normal fluid relaxation time (units of Δt, > 1/2)
    /// * `dx`         – lattice spacing
    /// * `dt`         – time step
    pub fn new(
        nx: usize,
        temperature: f64,
        t_critical: f64,
        alpha: f64,
        tau_n: f64,
        dx: f64,
        dt: f64,
    ) -> Self {
        let fs = superfluid_fraction(temperature, t_critical);
        let rho_s = vec![fs; nx];
        let rho_n = vec![1.0 - fs; nx];
        Self {
            nx,
            rho_s,
            rho_n,
            vel_s: vec![0.0; nx],
            vel_n: vec![0.0; nx],
            temperature,
            t_critical,
            alpha,
            tau_n,
            dx,
            dt,
            f_neq: vec![[0.0; Q]; nx],
        }
    }

    /// Superfluid fraction at current temperature.
    pub fn superfluid_fraction(&self) -> f64 {
        superfluid_fraction(self.temperature, self.t_critical)
    }

    /// Total density at site i.
    pub fn total_density(&self, i: usize) -> f64 {
        self.rho_s[i] + self.rho_n[i]
    }

    /// Lattice speed `c = Δx/Δt`.
    pub fn lattice_speed(&self) -> f64 {
        self.dx / self.dt
    }

    /// Isothermal sound speed `c_s = c/√3` of both components.
    pub fn sound_speed(&self) -> f64 {
        self.lattice_speed() / 3.0_f64.sqrt()
    }

    /// Kinematic viscosity `ν_n = c_s² (τ_n − ½) Δt` of the normal fluid.
    ///
    /// In one dimension the viscous stress is longitudinal, `2 ρ_n ν_n ∂ₓv_n`,
    /// so sound waves of wavenumber `k` are damped at the rate `ν_n k²`.
    pub fn normal_viscosity(&self) -> f64 {
        let cs = self.sound_speed();
        cs * cs * (self.tau_n - 0.5) * self.dt
    }

    /// Relaxation rate `κ = α ρ / ρ_n` of the slip velocity `v_s − v_n` at
    /// site `i` caused by mutual friction.
    ///
    /// The rate is zero where one of the components is absent, since there is
    /// nothing to exchange momentum with.
    pub fn slip_decay_rate(&self, i: usize) -> f64 {
        let (rho_s, rho_n) = (self.rho_s[i], self.rho_n[i]);
        if rho_s <= RHO_FLOOR || rho_n <= RHO_FLOOR {
            return 0.0;
        }
        self.alpha * (rho_s + rho_n) / rho_n
    }

    /// Total mass `Σ (ρ_s + ρ_n) Δx` of both components.
    pub fn total_mass(&self) -> f64 {
        let sum: f64 = self.rho_s.iter().chain(&self.rho_n).sum();
        sum * self.dx
    }

    /// Total momentum `Σ (ρ_s v_s + ρ_n v_n) Δx` of both components.
    pub fn total_momentum(&self) -> f64 {
        let superfluid: f64 = self.rho_s.iter().zip(&self.vel_s).map(|(r, v)| r * v).sum();
        let normal: f64 = self.rho_n.iter().zip(&self.vel_n).map(|(r, v)| r * v).sum();
        (superfluid + normal) * self.dx
    }

    /// Check that the lattice is well formed and the parameters are physical.
    fn validate(&self) -> Result<(), Error> {
        if self.nx == 0 {
            return Err(Error::InvalidGridDimensions { nx: 0, ny: 1 });
        }
        for (field, len) in [
            ("rho_s", self.rho_s.len()),
            ("rho_n", self.rho_n.len()),
            ("vel_s", self.vel_s.len()),
            ("vel_n", self.vel_n.len()),
        ] {
            if len != self.nx {
                return Err(Error::ComponentCountMismatch {
                    expected: self.nx,
                    found: len,
                    field: field.to_string(),
                });
            }
        }
        if !self.tau_n.is_finite() || self.tau_n <= 0.5 {
            return Err(Error::UnstableRelaxationTime { tau: self.tau_n });
        }
        for (name, value) in [("dx", self.dx), ("dt", self.dt)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: f64::MIN_POSITIVE,
                    max: f64::MAX,
                });
            }
        }
        if !self.alpha.is_finite() || self.alpha < 0.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "alpha".to_string(),
                value: self.alpha,
                min: 0.0,
                max: f64::MAX,
            });
        }
        self.check_fields()?;
        self.check_normal_speed()
    }

    /// Check that all fields are finite and the densities non-negative.
    fn check_fields(&self) -> Result<(), Error> {
        for (field, values) in [
            ("rho_s", &self.rho_s),
            ("rho_n", &self.rho_n),
            ("vel_s", &self.vel_s),
            ("vel_n", &self.vel_n),
        ] {
            if let Some((index, &value)) = values.iter().enumerate().find(|(_, v)| !v.is_finite()) {
                return Err(Error::NonFiniteField {
                    field: field.to_string(),
                    index,
                    value,
                });
            }
        }
        for (name, values) in [("rho_s", &self.rho_s), ("rho_n", &self.rho_n)] {
            if let Some(&value) = values.iter().find(|v| **v < 0.0) {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: 0.0,
                    max: f64::MAX,
                });
            }
        }
        Ok(())
    }

    /// Check that the normal fluid stays within the speed range of the D1Q3
    /// lattice, `|v_n| < sqrt(2/3) c`, where its equilibrium populations are
    /// positive.
    fn check_normal_speed(&self) -> Result<(), Error> {
        let c = self.lattice_speed();
        let limit = MAX_NORMAL_SPEED_RATIO * c;
        for (&rho, &v) in self.rho_n.iter().zip(&self.vel_n) {
            if rho > RHO_FLOOR && v.abs() >= limit {
                let cs = self.sound_speed();
                return Err(Error::MachNumberExceeded {
                    mach: v.abs() / cs,
                    limit: limit / cs,
                });
            }
        }
        Ok(())
    }

    /// Relax the slip velocity of every site by mutual friction over a time `h`.
    ///
    /// The local system `ρ_n v_n' = F`, `ρ_s v_s' = −F` with
    /// `F = α ρ_s (v_s − v_n)` conserves the momentum density
    /// `ρ_s v_s + ρ_n v_n`, i.e. the centre-of-mass velocity `v_cm`, and the
    /// slip velocity decays as `exp(−κ h)`; both velocities are then
    /// `v_s = v_cm + (ρ_n/ρ) slip`, `v_n = v_cm − (ρ_s/ρ) slip`.
    fn apply_mutual_friction(&mut self, h: f64) {
        if self.alpha == 0.0 {
            return;
        }
        for i in 0..self.nx {
            let (rho_s, rho_n) = (self.rho_s[i], self.rho_n[i]);
            if rho_s <= RHO_FLOOR || rho_n <= RHO_FLOOR {
                continue;
            }
            let rho = rho_s + rho_n;
            let slip = self.vel_s[i] - self.vel_n[i];
            let v_cm = (rho_s * self.vel_s[i] + rho_n * self.vel_n[i]) / rho;
            let slip_after = slip * (-self.alpha * rho / rho_n * h).exp();
            self.vel_s[i] = v_cm + rho_n / rho * slip_after;
            self.vel_n[i] = v_cm - rho_s / rho * slip_after;
        }
    }

    /// Collision and streaming of the normal-fluid populations.
    ///
    /// The populations are `f_q = f_q^eq(ρ_n, v_n) + f_neq_q`.  The BGK
    /// collision relaxes them towards the equilibrium with rate `1/τ_n`, the
    /// populations moving with `±c` are shifted by one site (periodic), and
    /// the density, velocity and non-equilibrium part are recomputed from the
    /// streamed populations.
    fn step_normal_fluid(&mut self) {
        let n = self.nx;
        let c = self.lattice_speed();
        let retained = 1.0 - 1.0 / self.tau_n;
        let post: Vec<[f64; Q]> = (0..n)
            .map(|i| {
                let feq = d1q3_equilibrium(self.rho_n[i], self.vel_n[i], c);
                std::array::from_fn(|q| feq[q] + retained * self.f_neq[i][q])
            })
            .collect();
        for (i, post_i) in post.iter().enumerate() {
            let west = (i + n - 1) % n;
            let east = (i + 1) % n;
            // The +c population arrives from the west neighbour, the -c
            // population from the east neighbour.
            let f = [post_i[0], post[west][1], post[east][2]];
            let rho = f[0] + f[1] + f[2];
            let velocity = if rho > RHO_FLOOR {
                c * (f[1] - f[2]) / rho
            } else {
                self.vel_n[i]
            };
            let feq = d1q3_equilibrium(rho, velocity, c);
            self.rho_n[i] = rho;
            self.vel_n[i] = velocity;
            self.f_neq[i] = std::array::from_fn(|q| f[q] - feq[q]);
        }
    }

    /// Advance the superfluid by one time step `Δt` of the isothermal Euler
    /// equations (see the type documentation).
    ///
    /// # Errors
    ///
    /// [`Error::MachNumberExceeded`] if the superfluid speed is so large that
    /// the explicit scheme would need more than `MAX_EULER_SUBSTEPS`
    /// sub-steps.
    fn step_superfluid(&mut self) -> Result<(), Error> {
        let n = self.nx;
        let cs = self.sound_speed();
        let mut rho = self.rho_s.clone();
        let mut mom: Vec<f64> = rho.iter().zip(&self.vel_s).map(|(r, v)| r * v).collect();

        let v_max = rho
            .iter()
            .zip(&self.vel_s)
            .filter(|(r, _)| **r > RHO_FLOOR)
            .fold(0.0_f64, |acc, (_, v)| acc.max(v.abs()));
        let substeps = ((v_max + cs) * self.dt / (EULER_CFL * self.dx)).ceil();
        if substeps > MAX_EULER_SUBSTEPS as f64 {
            return Err(Error::MachNumberExceeded {
                mach: v_max / cs,
                limit: (MAX_EULER_SUBSTEPS as f64 * EULER_CFL * self.dx / self.dt - cs) / cs,
            });
        }
        let substeps = (substeps as usize).max(1);
        let h = self.dt / substeps as f64;

        let mut scratch = EulerScratch::new(n);
        let mut d_rho = vec![0.0; n];
        let mut d_mom = vec![0.0; n];
        let mut stage_rho = vec![0.0; n];
        let mut stage_mom = vec![0.0; n];
        for _ in 0..substeps {
            // SSP-RK2 (Heun): U1 = U + h L(U), U+ = (U + U1 + h L(U1)) / 2.
            scratch.rhs(&rho, &mom, cs, self.dx, &mut d_rho, &mut d_mom);
            for i in 0..n {
                stage_rho[i] = rho[i] + h * d_rho[i];
                stage_mom[i] = mom[i] + h * d_mom[i];
            }
            scratch.rhs(&stage_rho, &stage_mom, cs, self.dx, &mut d_rho, &mut d_mom);
            for i in 0..n {
                rho[i] = 0.5 * (rho[i] + stage_rho[i] + h * d_rho[i]);
                mom[i] = 0.5 * (mom[i] + stage_mom[i] + h * d_mom[i]);
            }
        }

        for i in 0..n {
            self.rho_s[i] = rho[i];
            if rho[i] > RHO_FLOOR {
                self.vel_s[i] = mom[i] / rho[i];
            }
        }
        Ok(())
    }

    /// Advance the two-fluid state by one time step `Δt`.
    ///
    /// The step is a Strang splitting (see the type documentation): half a
    /// step of mutual friction, the transport of both components, and half a
    /// step of mutual friction.  It conserves the mass of each component and
    /// the total momentum to rounding error.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidGridDimensions`] for an empty lattice,
    /// * [`Error::ComponentCountMismatch`] if a field does not have `nx`
    ///   entries,
    /// * [`Error::UnstableRelaxationTime`] if `tau_n` is not larger than `1/2`,
    /// * [`Error::ParameterOutOfBounds`] if Δx or Δt is not positive and
    ///   finite, `alpha` is negative or not finite, or a density is negative,
    /// * [`Error::NonFiniteField`] if a field holds NaN or an infinity (also
    ///   checked after the update),
    /// * [`Error::MachNumberExceeded`] if the normal fluid reaches the speed
    ///   `sqrt(2/3) Δx/Δt` beyond which the D1Q3 equilibrium has negative
    ///   populations, or the superfluid speed is too large for the explicit
    ///   update.
    ///
    /// The state is not modified when an error is returned.
    pub fn try_step(&mut self) -> Result<(), Error> {
        self.validate()?;
        if self.f_neq.len() != self.nx {
            self.f_neq = vec![[0.0; Q]; self.nx];
        }
        let backup = self.clone();
        let result = self.advance();
        if result.is_err() {
            *self = backup;
        }
        result
    }

    /// The splitting sequence of [`try_step`](Self::try_step) on a validated
    /// lattice.
    fn advance(&mut self) -> Result<(), Error> {
        let half = 0.5 * self.dt;
        self.apply_mutual_friction(half);
        self.check_normal_speed()?;
        self.step_normal_fluid();
        self.step_superfluid()?;
        self.apply_mutual_friction(half);
        self.check_fields()?;
        self.check_normal_speed()
    }

    /// Advance the two-fluid state by one time step.
    ///
    /// Equivalent to [`try_step`](Self::try_step); an invalid configuration
    /// leaves the state untouched.
    pub fn step(&mut self) {
        // The configuration error, if any, is reported by `try_step`.
        let _ = self.try_step();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Helium-II-like lattice with unit spacing and time step, so `c = 1`.
    fn unit_lattice(nx: usize, temperature: f64, alpha: f64, tau_n: f64) -> SuperfluidLbm {
        SuperfluidLbm::new(nx, temperature, 2.17, alpha, tau_n, 1.0, 1.0)
    }

    /// Slip velocity `v_s − v_n` at site `i`.
    fn slip(sf: &SuperfluidLbm, i: usize) -> f64 {
        sf.vel_s[i] - sf.vel_n[i]
    }

    /// Bit patterns of a field, so that a `NaN` entry compares equal to itself.
    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    #[test]
    fn mutual_friction_relaxes_the_slip_velocity_exponentially() {
        // Uniform counter-flow: the slip must decay as exp(-alpha rho t / rho_n),
        // while the momentum density rho_s v_s + rho_n v_n stays constant.
        let mut sf = SuperfluidLbm::new(8, 1.0, 2.17, 0.05, 0.8, 1.0, 0.05);
        sf.vel_s.fill(0.3);
        sf.vel_n.fill(-0.1);
        let kappa = sf.slip_decay_rate(0);
        let rho = sf.total_density(0);
        assert!((kappa - 0.05 * rho / sf.rho_n[0]).abs() < 1e-14);
        let slip0 = slip(&sf, 0);
        let momentum0 = sf.rho_s[0] * sf.vel_s[0] + sf.rho_n[0] * sf.vel_n[0];
        let mut previous = slip0.abs();
        for step in 1..=60 {
            sf.step();
            let t = step as f64 * sf.dt;
            let expected = slip0 * (-kappa * t).exp();
            for i in 0..sf.nx {
                let s = slip(&sf, i);
                assert!(
                    (s - expected).abs() <= 1e-9 * slip0.abs(),
                    "step {step}: slip {s} vs {expected}"
                );
                let momentum = sf.rho_s[i] * sf.vel_s[i] + sf.rho_n[i] * sf.vel_n[i];
                assert!((momentum - momentum0).abs() < 1e-12);
            }
            // The slip never grows: the original sign error multiplied it by
            // 1 + alpha dt (1 + rho_s/rho_n) in every step.
            assert!(slip(&sf, 0).abs() < previous);
            previous = slip(&sf, 0).abs();
        }
        let decay = slip(&sf, 0) / slip0;
        assert!(decay > 0.0 && decay < 0.2, "decay factor {decay}");
    }

    #[test]
    fn mutual_friction_is_stable_for_arbitrarily_large_steps() {
        // kappa dt is of the order of 10^3 here; an explicit update would be
        // wildly unstable, the exact exponential just locks both fluids to the
        // centre-of-mass velocity.
        let mut sf = SuperfluidLbm::new(6, 1.0, 2.17, 10.0, 0.9, 100.0, 10.0);
        sf.vel_s.fill(0.4);
        sf.vel_n.fill(-0.2);
        let v_cm = (sf.rho_s[0] * 0.4 + sf.rho_n[0] * -0.2) / sf.total_density(0);
        sf.step();
        for i in 0..sf.nx {
            assert!((sf.vel_s[i] - v_cm).abs() < 1e-12, "{}", sf.vel_s[i]);
            assert!((sf.vel_n[i] - v_cm).abs() < 1e-12, "{}", sf.vel_n[i]);
        }
    }

    #[test]
    fn friction_conserves_the_centre_of_mass_velocity_for_unequal_densities() {
        // Superfluid-dominated mixture: the heavy component barely changes.
        let mut sf = unit_lattice(4, 0.4, 2.0, 0.7);
        sf.vel_s.fill(0.02);
        sf.vel_n.fill(0.0);
        let (rho_s, rho_n) = (sf.rho_s[0], sf.rho_n[0]);
        assert!(rho_s > 50.0 * rho_n);
        let v_cm = rho_s * 0.02 / (rho_s + rho_n);
        sf.step();
        // alpha rho / rho_n is large: the normal fluid is dragged to the
        // superfluid, which keeps almost all of its velocity.
        for i in 0..sf.nx {
            assert!((sf.vel_n[i] - v_cm).abs() < 1e-9);
            assert!((sf.vel_s[i] - v_cm).abs() < 1e-9);
            assert!(sf.vel_s[i] > 0.0199 && sf.vel_s[i] < 0.02);
        }
    }

    #[test]
    fn no_friction_without_a_partner_component() {
        // T = 0: no normal fluid, T > Tc: no superfluid.  Nothing to exchange.
        for temperature in [0.0, 3.0] {
            let mut sf = unit_lattice(8, temperature, 5.0, 0.8);
            sf.vel_s.fill(0.1);
            sf.vel_n.fill(-0.05);
            assert_eq!(sf.slip_decay_rate(3), 0.0);
            sf.step();
            if temperature == 0.0 {
                assert!(sf.vel_s.iter().all(|v| (v - 0.1).abs() < 1e-12));
            } else {
                assert!(sf.vel_n.iter().all(|v| (v + 0.05).abs() < 1e-12));
            }
        }
    }

    #[test]
    fn bgk_collision_relaxes_the_non_equilibrium_part_with_the_relaxation_time() {
        // A non-equilibrium population with neither mass nor momentum
        // (second moment only) relaxes by exactly 1 - 1/tau per step on a
        // uniform lattice, where streaming changes nothing.
        for &tau in &[0.6, 0.8, 1.0, 1.5, 2.0] {
            let mut sf = unit_lattice(8, 3.0, 0.0, tau);
            sf.rho_n.fill(0.8);
            sf.vel_n.fill(0.02);
            let perturbation = [-2.0e-3, 1.0e-3, 1.0e-3];
            sf.f_neq.fill(perturbation);
            sf.step();
            for (q, &part) in perturbation.iter().enumerate() {
                let expected = (1.0 - 1.0 / tau) * part;
                for site in &sf.f_neq {
                    assert!((site[q] - expected).abs() < 1e-15, "tau {tau} q {q}");
                }
            }
            assert!(sf.rho_n.iter().all(|r| (r - 0.8).abs() < 1e-15));
            assert!(sf.vel_n.iter().all(|v| (v - 0.02).abs() < 1e-15));
        }
    }

    #[test]
    fn normal_fluid_populations_carry_the_moments_of_the_fields() {
        let c = 2.5;
        for &(rho, u) in &[(1.0, 0.0), (0.7, 0.3), (2.0, -0.9), (1.3, 4.0)] {
            let f = d1q3_equilibrium(rho, u, c);
            let mass: f64 = f.iter().sum();
            let momentum = c * (f[1] - f[2]);
            let flux = c * c * (f[1] + f[2]);
            assert!((mass - rho).abs() < 1e-14);
            assert!((momentum - rho * u).abs() < 1e-14);
            assert!((flux - rho * (c * c / 3.0 + u * u)).abs() < 1e-13);
            // The moments are polynomial identities that hold at any speed;
            // the populations are positive up to |u| = sqrt(2/3) c.
            if u.abs() < MAX_NORMAL_SPEED_RATIO * c {
                assert!(f.iter().all(|x| *x > 0.0), "rho {rho} u {u}");
            } else {
                assert!(f[0] < 0.0);
            }
        }
    }

    #[test]
    fn normal_fluid_speeds_beyond_the_lattice_limit_are_rejected() {
        // c = 1: sqrt(2/3) = 0.8165 is the largest speed with a positive
        // equilibrium; the check applies to the input and to the result of
        // the friction half step.
        let mut sf = unit_lattice(8, 1.0, 0.0, 0.8);
        sf.vel_n[2] = 0.9;
        let err = sf.try_step().expect_err("supersonic normal fluid");
        assert!(matches!(err, Error::MachNumberExceeded { .. }), "{err:?}");
        assert_eq!(sf.vel_n[2], 0.9);

        // A fast superfluid drags a light normal fluid beyond the limit: the
        // input state is valid (the normal fluid is at rest), but the friction
        // half step locks both fluids to the centre-of-mass velocity
        // 0.9 * 2.0 + 0.1 * 0.0 = 1.8 > sqrt(2/3).
        let mut sf = SuperfluidLbm::new(8, 1.0, 2.17, 50.0, 0.8, 1.0, 1.0);
        sf.rho_s.fill(0.9);
        sf.rho_n.fill(0.1);
        sf.vel_s.fill(2.0);
        sf.vel_n.fill(0.0);
        sf.validate().expect("the input state itself is valid");
        let before = sf.clone();
        let err = sf
            .try_step()
            .expect_err("friction pushes v_n past the limit");
        assert!(matches!(err, Error::MachNumberExceeded { .. }), "{err:?}");
        // The rejected step leaves the whole state untouched.
        assert_eq!(bits(&sf.rho_s), bits(&before.rho_s));
        assert_eq!(bits(&sf.rho_n), bits(&before.rho_n));
        assert_eq!(bits(&sf.vel_s), bits(&before.vel_s));
        assert_eq!(bits(&sf.vel_n), bits(&before.vel_n));
        assert_eq!(sf.f_neq, before.f_neq);
    }

    /// Fourier amplitude `2/N sum_i rho_i cos(k x_i)` of a density profile.
    fn cosine_amplitude(rho: &[f64], k: f64) -> f64 {
        let n = rho.len() as f64;
        2.0 / n
            * rho
                .iter()
                .enumerate()
                .map(|(i, r)| r * (k * i as f64).cos())
                .sum::<f64>()
    }

    #[test]
    fn normal_fluid_sound_waves_follow_the_viscous_dispersion_relation() {
        // Linearised isothermal Navier-Stokes in 1-D:
        //   d_t rho' + d_x u' = 0,  d_t u' + c_s^2 d_x rho' = 2 nu d_xx u'
        // gives omega = -i nu k^2 +- k sqrt(c_s^2 - nu^2 k^2), so a standing
        // density wave evolves as
        //   a exp(-G t) [cos(w t) + (G/w) sin(w t)],  G = nu k^2.
        let nx = 128;
        let tau = 0.8;
        let mut sf = unit_lattice(nx, 3.0, 0.0, tau);
        let nu = sf.normal_viscosity();
        assert!((nu - (tau - 0.5) / 3.0).abs() < 1e-15);
        let cs = sf.sound_speed();
        assert!((cs - 1.0 / 3.0_f64.sqrt()).abs() < 1e-15);
        let k = 2.0 * PI / nx as f64;
        let amplitude = 1e-4;
        for (i, rho) in sf.rho_n.iter_mut().enumerate() {
            *rho = 1.0 + amplitude * (k * i as f64).cos();
        }
        let gamma = nu * k * k;
        let omega = k * (cs * cs - nu * nu * k * k).sqrt();
        let steps = 400;
        let mut worst = 0.0_f64;
        for step in 1..=steps {
            sf.step();
            let t = step as f64;
            let expected = amplitude
                * (-gamma * t).exp()
                * ((omega * t).cos() + gamma / omega * (omega * t).sin());
            let measured = cosine_amplitude(&sf.rho_n, k);
            worst = worst.max((measured - expected).abs() / amplitude);
        }
        assert!(worst < 0.01, "worst relative deviation {worst}");
        // Mass is conserved to rounding.
        let mass: f64 = sf.rho_n.iter().sum();
        let mass0 = nx as f64;
        assert!((mass - mass0).abs() < 1e-11);
    }

    #[test]
    fn normal_fluid_damping_grows_with_the_relaxation_time() {
        // Larger tau_n means larger viscosity: a standing wave of the same
        // wavelength keeps less of its amplitude envelope.
        let envelope = |tau: f64| -> f64 {
            let nx = 64;
            let mut sf = unit_lattice(nx, 3.0, 0.0, tau);
            let k = 2.0 * PI / nx as f64;
            for (i, rho) in sf.rho_n.iter_mut().enumerate() {
                *rho = 1.0 + 1e-4 * (k * i as f64).cos();
            }
            let mut peak = 0.0_f64;
            for step in 0..600 {
                sf.step();
                if step >= 400 {
                    peak = peak.max(cosine_amplitude(&sf.rho_n, k).abs());
                }
            }
            peak
        };
        let thin = envelope(0.55);
        let thick = envelope(1.2);
        assert!(thin > thick, "peak {thin} (tau 0.55) vs {thick} (tau 1.2)");
        assert!(thin / thick > 1.5);
    }

    #[test]
    fn superfluid_density_pulse_splits_into_two_sound_pulses() {
        let nx = 400;
        let mut sf = unit_lattice(nx, 0.0, 0.0, 0.8);
        let centre = 200.0;
        let width = 6.0;
        let bump = 0.01;
        for (i, rho) in sf.rho_s.iter_mut().enumerate() {
            let x = i as f64 - centre;
            *rho = 1.0 + bump * (-0.5 * x * x / (width * width)).exp();
        }
        let mass0: f64 = sf.rho_s.iter().sum();
        let steps = 120;
        for _ in 0..steps {
            sf.step();
        }
        let travelled = sf.sound_speed() * steps as f64;
        let peak_of = |range: std::ops::Range<usize>| -> (usize, f64) {
            range
                .map(|i| (i, sf.rho_s[i] - 1.0))
                .fold(
                    (0, f64::MIN),
                    |best, cur| if cur.1 > best.1 { cur } else { best },
                )
        };
        let (right_at, right_peak) = peak_of(200..nx);
        let (left_at, left_peak) = peak_of(0..200);
        assert!(
            (right_at as f64 - (centre + travelled)).abs() <= 1.5,
            "right pulse at {right_at}, expected {}",
            centre + travelled
        );
        assert!(
            (left_at as f64 - (centre - travelled)).abs() <= 1.5,
            "left pulse at {left_at}, expected {}",
            centre - travelled
        );
        // Each pulse carries half of the initial amplitude (up to numerical
        // diffusion of the limited scheme).
        assert!(right_peak > 0.8 * 0.5 * bump && right_peak < 1.05 * 0.5 * bump);
        assert!((right_peak - left_peak).abs() < 1e-9);
        // Conservative: mass preserved, no net momentum for a symmetric start.
        let mass: f64 = sf.rho_s.iter().sum();
        assert!((mass - mass0).abs() < 1e-10);
        let momentum: f64 = sf.rho_s.iter().zip(&sf.vel_s).map(|(r, v)| r * v).sum();
        assert!(momentum.abs() < 1e-10);
        assert!(sf.rho_s.iter().all(|r| *r > 0.99));
    }

    #[test]
    fn superfluid_update_uses_the_lattice_spacing_and_time_step() {
        // The dynamics only depend on c = dx/dt and on alpha dt: the same
        // physical problem on a lattice of spacing 2 and time step 2 (and half
        // the friction coefficient) reproduces the unit lattice site by site.
        let build = |dx: f64, dt: f64, alpha: f64| -> SuperfluidLbm {
            let mut sf = SuperfluidLbm::new(48, 1.6, 2.17, alpha, 0.75, dx, dt);
            for i in 0..48 {
                let phase = 2.0 * PI * i as f64 / 48.0;
                sf.rho_s[i] *= 1.0 + 0.04 * phase.sin();
                sf.rho_n[i] *= 1.0 + 0.04 * (2.0 * phase).cos();
                sf.vel_s[i] = 0.03 * phase.cos();
                sf.vel_n[i] = -0.02 * (3.0 * phase).sin();
            }
            sf
        };
        let mut unit = build(1.0, 1.0, 0.2);
        let mut scaled = build(2.0, 2.0, 0.1);
        for _ in 0..80 {
            unit.step();
            scaled.step();
        }
        for i in 0..48 {
            assert!((unit.rho_s[i] - scaled.rho_s[i]).abs() < 1e-12);
            assert!((unit.rho_n[i] - scaled.rho_n[i]).abs() < 1e-12);
            assert!((unit.vel_s[i] - scaled.vel_s[i]).abs() < 1e-12);
            assert!((unit.vel_n[i] - scaled.vel_n[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn nonuniform_flow_conserves_mass_and_momentum() {
        let nx = 64;
        let mut sf = unit_lattice(nx, 1.5, 0.3, 0.7);
        for i in 0..nx {
            let phase = 2.0 * PI * i as f64 / nx as f64;
            sf.rho_s[i] *= 1.0 + 0.05 * phase.sin();
            sf.rho_n[i] *= 1.0 + 0.08 * phase.cos();
            sf.vel_s[i] = 0.06 * (2.0 * phase).sin() + 0.02;
            sf.vel_n[i] = -0.04 * phase.cos();
        }
        let mass_s0: f64 = sf.rho_s.iter().sum();
        let mass_n0: f64 = sf.rho_n.iter().sum();
        let momentum0 = sf.total_momentum();
        assert!(momentum0.abs() > 0.1, "start with real momentum");
        let mut slip0 = 0.0_f64;
        for i in 0..nx {
            slip0 = slip0.max(slip(&sf, i).abs());
        }
        for _ in 0..300 {
            sf.step();
        }
        assert!((sf.rho_s.iter().sum::<f64>() - mass_s0).abs() < 1e-10);
        assert!((sf.rho_n.iter().sum::<f64>() - mass_n0).abs() < 1e-10);
        assert!((sf.total_momentum() - momentum0).abs() < 1e-10);
        assert!(sf.rho_s.iter().all(|r| *r > 0.0));
        assert!(sf.rho_n.iter().all(|r| *r > 0.0));
        let slip_end = (0..nx).fold(0.0_f64, |m, i| m.max(slip(&sf, i).abs()));
        assert!(
            slip_end < 0.2 * slip0,
            "friction must damp the slip: {slip_end} vs {slip0}"
        );
        assert!(sf.vel_s.iter().chain(&sf.vel_n).all(|v| v.is_finite()));
    }

    #[test]
    fn fields_modified_between_steps_stay_consistent_with_the_populations() {
        // The fields are the state: kicking the normal fluid between two steps
        // changes the momentum by exactly the kick and nothing else.
        let nx = 32;
        let mut sf = unit_lattice(nx, 1.2, 0.1, 0.9);
        for i in 0..nx {
            sf.vel_n[i] = 0.02 * (2.0 * PI * i as f64 / nx as f64).sin();
        }
        for _ in 0..25 {
            sf.step();
        }
        let before = sf.total_momentum();
        for v in sf.vel_n.iter_mut() {
            *v += 0.03;
        }
        let kick = sf.rho_n.iter().sum::<f64>() * 0.03;
        assert!((sf.total_momentum() - (before + kick)).abs() < 1e-12);
        let expected = sf.total_momentum();
        for _ in 0..25 {
            sf.step();
        }
        assert!((sf.total_momentum() - expected).abs() < 1e-10);
    }

    #[test]
    fn invalid_configurations_are_reported_and_leave_the_state_untouched() {
        let base = unit_lattice(8, 1.0, 0.2, 0.8);
        let assert_untouched = |mut sf: SuperfluidLbm, check: &dyn Fn(&Error) -> bool| {
            let before = sf.clone();
            let err = sf.try_step().expect_err("configuration must be rejected");
            assert!(check(&err), "{err:?}");
            assert_eq!(bits(&sf.rho_s), bits(&before.rho_s));
            assert_eq!(bits(&sf.rho_n), bits(&before.rho_n));
            assert_eq!(bits(&sf.vel_s), bits(&before.vel_s));
            assert_eq!(bits(&sf.vel_n), bits(&before.vel_n));
            assert_eq!(sf.f_neq, before.f_neq);
            sf.step();
            assert_eq!(bits(&sf.vel_s), bits(&before.vel_s));
            assert_eq!(bits(&sf.vel_n), bits(&before.vel_n));
        };

        for tau in [0.5, 0.2, 0.0, f64::NAN, f64::INFINITY] {
            let mut sf = base.clone();
            sf.tau_n = tau;
            assert_untouched(sf, &|e| matches!(e, Error::UnstableRelaxationTime { .. }));
        }
        for (dx, dt) in [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (1.0, f64::NAN)] {
            let mut sf = base.clone();
            sf.dx = dx;
            sf.dt = dt;
            assert_untouched(sf, &|e| matches!(e, Error::ParameterOutOfBounds { .. }));
        }
        let mut sf = base.clone();
        sf.alpha = -0.1;
        assert_untouched(sf, &|e| matches!(e, Error::ParameterOutOfBounds { .. }));
        let mut sf = base.clone();
        sf.rho_s[2] = -1e-3;
        assert_untouched(sf, &|e| matches!(e, Error::ParameterOutOfBounds { .. }));
        let mut sf = base.clone();
        sf.vel_n[5] = f64::NAN;
        assert_untouched(sf, &|e| matches!(e, Error::NonFiniteField { index: 5, .. }));
        let mut sf = base.clone();
        sf.vel_s.pop();
        assert_untouched(sf, &|e| matches!(e, Error::ComponentCountMismatch { .. }));
        let mut sf = base.clone();
        sf.rho_n.push(0.1);
        assert_untouched(sf, &|e| matches!(e, Error::ComponentCountMismatch { .. }));
        let mut empty = unit_lattice(0, 1.0, 0.2, 0.8);
        assert!(matches!(
            empty.try_step().expect_err("no sites"),
            Error::InvalidGridDimensions { .. }
        ));
        empty.step();
    }

    #[test]
    fn superfluid_speeds_beyond_the_explicit_limit_are_rejected() {
        let mut sf = unit_lattice(16, 0.0, 0.0, 0.8);
        sf.vel_s[4] = 1.0e6;
        let before = sf.clone();
        let err = sf.try_step().expect_err("Mach number far too large");
        assert!(matches!(err, Error::MachNumberExceeded { .. }), "{err:?}");
        assert_eq!(sf.vel_s, before.vel_s);
        assert_eq!(sf.rho_s, before.rho_s);
    }

    #[test]
    fn existing_constructor_contract_is_unchanged() {
        let sf = SuperfluidLbm::new(16, 1.0, 2.17, 0.1, 0.6, 1.0, 0.01);
        let fs = superfluid_fraction(1.0, 2.17);
        assert_eq!(sf.nx, 16);
        assert!((sf.superfluid_fraction() - fs).abs() < 1e-15);
        for i in 0..16 {
            assert!((sf.total_density(i) - 1.0).abs() < 1e-15);
            assert!((sf.rho_s[i] - fs).abs() < 1e-15);
        }
        assert!((sf.lattice_speed() - 100.0).abs() < 1e-12);
        assert!((sf.total_mass() - 16.0).abs() < 1e-12);
        assert_eq!(sf.total_momentum(), 0.0);
    }
}
