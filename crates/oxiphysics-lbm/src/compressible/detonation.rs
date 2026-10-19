// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! One-dimensional reactive-Euler detonation solver.
//!
//! Solves the 1-D reactive Euler equations for a single-step Arrhenius
//! mixture with one adiabatic exponent γ (ideal gas, `R = 1`, `T = p/ρ`):
//!
//! ```text
//! ∂t ρ    + ∂x(ρu)          = 0
//! ∂t (ρu) + ∂x(ρu² + p)     = 0
//! ∂t E    + ∂x((E + p)u)    = ρ q w
//! ∂t (ρλ) + ∂x(ρλu)         = ρ w,      w = A exp(−Ea/T)(1 − λ)
//! ```
//!
//! with `E = p/(γ−1) + ρu²/2`.  The transport part uses a first-order
//! finite-volume scheme with the Rusanov (local Lax-Friedrichs) interface
//! flux and transmissive (zero-gradient) ends; the reaction is integrated
//! exactly over each sub-step at frozen temperature (Godunov splitting), so
//! the released heat `ρ q Δλ` is added to `E` without stiffness limits.
//! Each call to [`DetonationLbm::step`] advances by `dt`, sub-cycling to
//! respect a CFL number of [`DetonationLbm::CFL`].

/// Detonation wave simulation (1-D reactive Euler, Chapman-Jouguet/ZND).
///
/// The unburnt mixture is at rest with `ρ = rho0`, `p = p0`, `λ = 0`.  The
/// detonation is initiated by a driver region at the left end filled with
/// the burnt Chapman-Jouguet state (see [`DetonationLbm::ignite`]); the front
/// then develops a ZND structure and propagates at the CJ speed `dcj`.
#[derive(Debug, Clone)]
pub struct DetonationLbm {
    /// Number of cells.
    pub n: usize,
    /// Density profile.
    pub rho: Vec<f64>,
    /// Velocity profile.
    pub vel: Vec<f64>,
    /// Pressure profile.
    pub pressure: Vec<f64>,
    /// Reaction progress variable (0=unburnt, 1=burnt).
    pub progress: Vec<f64>,
    /// Heat release per unit mass Q.
    pub heat_release: f64,
    /// Activation energy Ea (in units of `T = p/ρ`).
    pub activation_energy: f64,
    /// Pre-exponential factor A.
    pub reaction_rate_a: f64,
    /// Gamma.
    pub gamma: f64,
    /// Chapman-Jouguet detonation velocity for the unburnt state.
    pub dcj: f64,
    /// Grid spacing.
    pub dx: f64,
    /// Time step advanced by one call to [`DetonationLbm::step`].
    pub dt: f64,
    /// Unburnt density.
    pub rho0: f64,
    /// Unburnt pressure.
    pub p0: f64,
    /// Elapsed simulation time.
    pub time: f64,
    /// Time-integrated net mass flux through the two open ends
    /// (`∫ (F_left − F_right) dt`); `Σρ dx − initial mass` equals this value.
    pub boundary_mass_inflow: f64,
}

/// Conserved variables `[ρ, ρu, E, ρλ]`.
type Conserved = [f64; 4];

impl DetonationLbm {
    /// CFL number used for internal sub-cycling.
    pub const CFL: f64 = 0.45;

    /// Create a new detonation simulation.
    ///
    /// The unburnt state is `ρ0 = 1`, `p0 = 1`, `u = 0`; the leftmost
    /// `max(2, n/20)` cells are ignited with the CJ burnt state.
    pub fn new(
        n: usize,
        dx: f64,
        dt: f64,
        gamma: f64,
        heat_release: f64,
        activation_energy: f64,
        reaction_rate_a: f64,
    ) -> Self {
        let rho0 = 1.0;
        let p0 = 1.0;
        let dcj = Self::cj_velocity_exact(gamma, heat_release, p0, rho0);
        let mut sim = Self {
            n,
            rho: vec![rho0; n],
            vel: vec![0.0; n],
            pressure: vec![p0; n],
            progress: vec![0.0; n],
            heat_release,
            activation_energy,
            reaction_rate_a,
            gamma,
            dcj,
            dx,
            dt,
            rho0,
            p0,
            time: 0.0,
            boundary_mass_inflow: 0.0,
        };
        sim.ignite(2usize.max(n / 20));
        sim
    }

    /// Chapman-Jouguet detonation velocity in the strong-detonation limit
    /// `D ≈ sqrt(2(γ²−1)Q)` (valid when `Q ≫ c0²`).
    pub fn cj_velocity(gamma: f64, q: f64) -> f64 {
        (2.0 * (gamma * gamma - 1.0) * q).sqrt()
    }

    /// Exact CJ Mach number `M = sqrt(1 + H) + sqrt(H)`,
    /// `H = (γ²−1)Q/(2c0²)`, `c0² = γ p0/ρ0` (single-γ mixture).
    pub fn cj_mach(gamma: f64, q: f64, p0: f64, rho0: f64) -> f64 {
        let c0_sq = gamma * p0 / rho0;
        let h = (gamma * gamma - 1.0) * q / (2.0 * c0_sq);
        (1.0 + h).sqrt() + h.sqrt()
    }

    /// Exact CJ detonation velocity `D = M_cj c0` for the unburnt state
    /// `(p0, ρ0)`.
    pub fn cj_velocity_exact(gamma: f64, q: f64, p0: f64, rho0: f64) -> f64 {
        Self::cj_mach(gamma, q, p0, rho0) * (gamma * p0 / rho0).sqrt()
    }

    /// CJ pressure ratio `p_cj/p0 = (1 + γ M_cj²)/(γ + 1)`.
    ///
    /// Follows from the Rayleigh line together with the CJ condition that
    /// the burnt gas leaves the front sonically.  (The normal-shock ratio
    /// `(2γM² − (γ−1))/(γ+1)` is the von Neumann spike, not the CJ state.)
    pub fn cj_pressure_ratio(gamma: f64, m_cj: f64) -> f64 {
        (1.0 + gamma * m_cj * m_cj) / (gamma + 1.0)
    }

    /// Burnt CJ state `(ρ, u, p)` behind a detonation running into the
    /// unburnt state at rest.
    pub fn cj_state(&self) -> (f64, f64, f64) {
        let m = Self::cj_mach(self.gamma, self.heat_release, self.p0, self.rho0);
        let g = self.gamma;
        let p_cj = self.p0 * Self::cj_pressure_ratio(g, m);
        let rho_cj = self.rho0 * (g + 1.0) * m * m / (1.0 + g * m * m);
        let d = m * (g * self.p0 / self.rho0).sqrt();
        let u_cj = d * (1.0 - self.rho0 / rho_cj);
        (rho_cj, u_cj, p_cj)
    }

    /// Fill the leftmost `cells` cells with the burnt CJ state (λ = 1).
    pub fn ignite(&mut self, cells: usize) {
        let (rho_cj, u_cj, p_cj) = self.cj_state();
        let count = cells.min(self.n);
        for i in 0..count {
            self.rho[i] = rho_cj;
            self.vel[i] = u_cj;
            self.pressure[i] = p_cj;
            self.progress[i] = 1.0;
        }
    }

    /// Arrhenius rate constant `A exp(−Ea/T)`.
    fn rate_constant(&self, temp: f64) -> f64 {
        self.reaction_rate_a * (-self.activation_energy / temp.max(0.01)).exp()
    }

    /// Arrhenius reaction rate `w = A exp(−Ea/T)(1 − λ)`.
    pub fn reaction_rate(&self, temp: f64, progress: f64) -> f64 {
        if progress >= 1.0 {
            return 0.0;
        }
        self.rate_constant(temp) * (1.0 - progress)
    }

    /// Total mass `Σ ρ dx`.
    pub fn total_mass(&self) -> f64 {
        self.rho.iter().sum::<f64>() * self.dx
    }

    /// Position of the leading shock: the rightmost location where the
    /// pressure crosses `p0 + (p_cj − p0)/2`, linearly interpolated between
    /// cell centres `(i + 1/2) dx`.  `None` if no cell exceeds the threshold.
    pub fn shock_position(&self) -> Option<f64> {
        let (_, _, p_cj) = self.cj_state();
        let threshold = self.p0 + 0.5 * (p_cj - self.p0);
        let i = self.pressure.iter().rposition(|&p| p > threshold)?;
        let x_i = (i as f64 + 0.5) * self.dx;
        match self.pressure.get(i + 1) {
            Some(&p_next) => {
                let p_i = self.pressure[i];
                let frac = (p_i - threshold) / (p_i - p_next);
                Some(x_i + frac * self.dx)
            }
            None => Some(x_i),
        }
    }

    fn conserved(&self, i: usize) -> Conserved {
        let rho = self.rho[i];
        let u = self.vel[i];
        let p = self.pressure[i];
        [
            rho,
            rho * u,
            p / (self.gamma - 1.0) + 0.5 * rho * u * u,
            rho * self.progress[i],
        ]
    }

    /// Primitive `(ρ, u, p, λ)` from conserved variables.
    fn primitive(&self, c: &Conserved) -> (f64, f64, f64, f64) {
        let rho = c[0];
        let u = c[1] / rho;
        let p = (self.gamma - 1.0) * (c[2] - 0.5 * c[1] * u);
        (rho, u, p, c[3] / rho)
    }

    /// Physical flux and maximal signal speed of a state.
    fn flux_and_speed(&self, c: &Conserved) -> (Conserved, f64) {
        let (rho, u, p, lam) = self.primitive(c);
        let sound = (self.gamma * p.max(0.0) / rho).sqrt();
        (
            [rho * u, rho * u * u + p, (c[2] + p) * u, rho * lam * u],
            u.abs() + sound,
        )
    }

    /// Rusanov interface flux between left and right states.
    fn rusanov(&self, left: &Conserved, right: &Conserved) -> Conserved {
        let (fl, sl) = self.flux_and_speed(left);
        let (fr, sr) = self.flux_and_speed(right);
        let s = sl.max(sr);
        let mut flux = [0.0f64; 4];
        for k in 0..4 {
            flux[k] = 0.5 * (fl[k] + fr[k]) - 0.5 * s * (right[k] - left[k]);
        }
        flux
    }

    fn max_signal_speed(&self) -> f64 {
        (0..self.n)
            .map(|i| {
                let sound = (self.gamma * self.pressure[i].max(0.0) / self.rho[i]).sqrt();
                self.vel[i].abs() + sound
            })
            .fold(0.0f64, f64::max)
    }

    /// Advance the transport and reaction by one sub-step `h`.
    fn advance(&mut self, h: f64) {
        let n = self.n;
        if n == 0 {
            return;
        }
        let states: Vec<Conserved> = (0..n).map(|i| self.conserved(i)).collect();
        // Interface fluxes F_{i-1/2}, i = 0..=n, with transmissive ghosts.
        let fluxes: Vec<Conserved> = (0..=n)
            .map(|i| {
                let left = &states[i.saturating_sub(1)];
                let right = &states[i.min(n - 1)];
                self.rusanov(left, right)
            })
            .collect();
        self.boundary_mass_inflow += h * (fluxes[0][0] - fluxes[n][0]);
        let ratio = h / self.dx;
        for i in 0..n {
            let mut c = states[i];
            for k in 0..4 {
                c[k] -= ratio * (fluxes[i + 1][k] - fluxes[i][k]);
            }
            let (rho, u, p, lam) = self.primitive(&c);
            // Exact integration of dλ/dt = k(T)(1 − λ) at frozen T.
            let lam_old = lam.clamp(0.0, 1.0);
            let k_rate = if lam_old < 1.0 {
                self.rate_constant(p / rho)
            } else {
                0.0
            };
            let lam_new = 1.0 - (1.0 - lam_old) * (-k_rate * h).exp();
            let released = rho * self.heat_release * (lam_new - lam_old);
            let energy = c[2] + released;
            self.rho[i] = rho;
            self.vel[i] = u;
            self.pressure[i] = (self.gamma - 1.0) * (energy - 0.5 * rho * u * u);
            self.progress[i] = lam_new;
        }
    }

    /// Advance the detonation by `dt`, sub-cycling at CFL ≤ [`Self::CFL`].
    pub fn step(&mut self) {
        let mut remaining = self.dt;
        while remaining > 0.0 {
            let speed = self.max_signal_speed();
            let dt_cfl = Self::CFL * self.dx / speed;
            let h = if dt_cfl.is_finite() && dt_cfl > 0.0 {
                remaining.min(dt_cfl)
            } else {
                remaining
            };
            self.advance(h);
            remaining -= h;
        }
        self.time += self.dt;
    }

    /// Induction zone length estimate.
    pub fn induction_length(&self) -> f64 {
        let dcj = self.dcj;
        let ea = self.activation_energy;
        let a = self.reaction_rate_a;
        dcj / (a * (-ea).exp())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cj_pressure_ratio_rayleigh_line() {
        let r = DetonationLbm::cj_pressure_ratio(1.4, 5.0);
        assert!((r - 15.0).abs() < 1e-12, "ratio = {r}");
        let r1 = DetonationLbm::cj_pressure_ratio(1.4, 1.0);
        assert!((r1 - 1.0).abs() < 1e-12);
        // Not a constant: depends on Mach number.
        assert!(DetonationLbm::cj_pressure_ratio(1.4, 3.0) > 1.0);
    }

    #[test]
    fn cj_state_satisfies_rankine_hugoniot_and_sonic_condition() {
        let det = DetonationLbm::new(10, 0.1, 0.01, 1.4, 10.0, 10.0, 10.0);
        let g = det.gamma;
        let (rho1, u1, p1) = det.cj_state();
        let d = det.dcj;
        let (rho0, p0) = (det.rho0, det.p0);
        // mass, momentum, energy across the front (front frame)
        let w1 = d - u1;
        assert!((rho0 * d - rho1 * w1).abs() < 1e-12);
        assert!((p0 + rho0 * d * d - (p1 + rho1 * w1 * w1)).abs() < 1e-10);
        let h = |p: f64, r: f64| g / (g - 1.0) * p / r;
        let lhs = h(p0, rho0) + 0.5 * d * d + det.heat_release;
        let rhs = h(p1, rho1) + 0.5 * w1 * w1;
        assert!((lhs - rhs).abs() < 1e-10);
        // CJ: burnt gas leaves sonically
        assert!((w1 - (g * p1 / rho1).sqrt()).abs() < 1e-12);
        // exact value for q=10, p0=rho0=1, gamma=1.4
        assert!((d - 4.680_870_149_618_411).abs() < 1e-9);
    }

    #[test]
    fn new_handles_tiny_grids() {
        let det0 = DetonationLbm::new(0, 0.1, 0.01, 1.4, 10.0, 10.0, 10.0);
        assert!(det0.rho.is_empty());
        let mut det1 = DetonationLbm::new(1, 0.1, 0.01, 1.4, 10.0, 10.0, 10.0);
        assert_eq!(det1.progress, vec![1.0]);
        det1.step();
        assert!(det1.rho[0].is_finite());
    }

    /// Detonation parameters with negligible auto-ignition of the unburnt
    /// gas (Ea/T0 = 20: rate ~5e-7 at T0 = 1) and a post-shock reaction
    /// zone of about six cells at the von Neumann state (T ≈ 4).
    fn resolved_detonation() -> DetonationLbm {
        DetonationLbm::new(600, 0.1, 0.05, 1.4, 10.0, 20.0, 250.0)
    }

    #[test]
    fn density_and_velocity_evolve() {
        let mut det = resolved_detonation();
        let rho_init = det.rho.clone();
        for _ in 0..20 {
            det.step();
        }
        let changed_rho = det
            .rho
            .iter()
            .zip(rho_init.iter())
            .filter(|(a, b)| (*a - *b).abs() > 1e-3)
            .count();
        assert!(changed_rho > 5, "density never changed: {changed_rho}");
        let moving = det.vel.iter().skip(40).filter(|v| v.abs() > 1e-2).count();
        assert!(moving > 3, "no flow induced ahead of the driver");
        assert!(det.rho.iter().all(|r| *r > 0.0 && r.is_finite()));
        assert!(det.pressure.iter().all(|p| *p > 0.0 && p.is_finite()));
    }

    #[test]
    fn total_mass_conserved_including_boundary_flux() {
        let mut det = resolved_detonation();
        let m0 = det.total_mass();
        for _ in 0..60 {
            det.step();
        }
        let m1 = det.total_mass();
        let err = (m1 - m0 - det.boundary_mass_inflow).abs();
        assert!(err < 1e-10 * m0, "mass drift {err}");
    }

    #[test]
    fn unburnt_gas_does_not_auto_ignite() {
        let mut det = resolved_detonation();
        for _ in 0..40 {
            det.step();
        }
        let front = det.shock_position().unwrap_or(0.0);
        let ahead = ((front / det.dx) as usize + 10).min(det.n);
        let max_ahead = det.progress[ahead..].iter().fold(0.0f64, |m, &l| m.max(l));
        assert!(max_ahead < 1e-3, "unburnt gas reacted: {max_ahead}");
    }

    #[test]
    fn front_speed_matches_cj_velocity() {
        let mut det = resolved_detonation();
        // Let the ZND structure develop.
        while det.time < 3.0 {
            det.step();
        }
        let (t1, x1) = (det.time, det.shock_position().unwrap_or(0.0));
        while det.time < 9.0 {
            det.step();
        }
        let (t2, x2) = (det.time, det.shock_position().unwrap_or(0.0));
        let speed = (x2 - x1) / (t2 - t1);
        let rel = (speed - det.dcj).abs() / det.dcj;
        assert!(
            rel < 0.15,
            "front speed {speed} vs D_cj {} (rel {rel})",
            det.dcj
        );
        // The front is a reacting detonation, not an inert shock: the gas
        // behind it is burnt.
        let behind = ((x2 / det.dx) as usize).saturating_sub(30);
        assert!(det.progress[behind] > 0.95);
    }
}
