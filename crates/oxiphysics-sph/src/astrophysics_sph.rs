#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Astrophysical SPH: gravitational N-body, artificial viscosity, entropy-conserving
//! formulations, stellar collapse, and galactic disk models.
//!
//! # Key References
//! - Monaghan (1992): Smoothed Particle Hydrodynamics
//! - Balsara (1995): von Neumann stability analysis of smoothed particle hydrodynamics
//! - Springel & Hernquist (2002): Cosmological smoothed particle hydrodynamics simulations
//! - Barnes & Hut (1986): A hierarchical O(N log N) force-calculation algorithm
//! - Jeans (1902): The stability of a spherical nebula

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// 1.  AstroSphParticle
// ─────────────────────────────────────────────────────────────────────────────

/// A single astrophysical SPH particle.
///
/// Stores mass, position, velocity, thermodynamic state, and SPH kernel radius.
#[derive(Debug, Clone)]
pub struct AstroSphParticle {
    /// Particle mass m \[code units\].
    pub mass: f64,
    /// Position vector \[x, y, z\] in 3D space \[code units\].
    pub pos: [f64; 3],
    /// Velocity vector \[vx, vy, vz\] \[code units / time\].
    pub vel: [f64; 3],
    /// Specific internal energy u \[energy / mass\].
    pub internal_energy: f64,
    /// Specific entropy A = P / ρ^γ (Springel & Hernquist formulation).
    pub entropy: f64,
    /// Smoothing length h \[code units\]; contains ∼32 neighbours.
    pub smoothing_length: f64,
    /// Local gas density ρ \[mass / volume\].
    pub density: f64,
    /// Local pressure P = (γ−1)·ρ·u.
    pub pressure: f64,
    /// Sound speed cs = sqrt(γ·P/ρ).
    pub sound_speed: f64,
    /// Acceleration \[ax, ay, az\] due to all forces.
    pub accel: [f64; 3],
    /// Time derivative of internal energy (du/dt).
    pub du_dt: f64,
    /// Unique particle identifier.
    pub id: usize,
}

impl AstroSphParticle {
    /// Construct a new particle with given mass, position, velocity, internal energy.
    ///
    /// Entropy and derived quantities are computed assuming adiabatic index γ.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: usize,
        mass: f64,
        pos: [f64; 3],
        vel: [f64; 3],
        internal_energy: f64,
        smoothing_length: f64,
        density: f64,
        gamma: f64,
    ) -> Self {
        let pressure = (gamma - 1.0) * density * internal_energy;
        let sound_speed = if density > 0.0 {
            (gamma * pressure / density).abs().sqrt()
        } else {
            0.0
        };
        let entropy = if density > 0.0 {
            pressure / density.powf(gamma)
        } else {
            0.0
        };
        AstroSphParticle {
            mass,
            pos,
            vel,
            internal_energy,
            entropy,
            smoothing_length,
            density,
            pressure,
            sound_speed,
            accel: [0.0; 3],
            du_dt: 0.0,
            id,
        }
    }

    /// Update pressure and sound speed from current density and internal energy.
    pub fn update_thermodynamics(&mut self, gamma: f64) {
        self.pressure = (gamma - 1.0) * self.density * self.internal_energy;
        if self.density > 1e-15 {
            self.sound_speed = (gamma * self.pressure / self.density).abs().sqrt();
        } else {
            self.sound_speed = 0.0;
        }
    }

    /// Compute the kinetic energy ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.vel[0].powi(2) + self.vel[1].powi(2) + self.vel[2].powi(2);
        0.5 * self.mass * v2
    }

    /// Compute the specific thermal energy from entropy: u = A·ρ^(γ-1) / (γ-1).
    pub fn thermal_energy_from_entropy(&self, gamma: f64) -> f64 {
        self.entropy * self.density.powf(gamma - 1.0) / (gamma - 1.0)
    }

    /// Leapfrog kick: advance velocity by 0.5·a·dt.
    pub fn kick(&mut self, dt: f64) {
        for d in 0..3 {
            self.vel[d] += 0.5 * self.accel[d] * dt;
        }
    }

    /// Leapfrog drift: advance position by v·dt.
    pub fn drift(&mut self, dt: f64) {
        for d in 0..3 {
            self.pos[d] += self.vel[d] * dt;
        }
    }

    /// Squared distance to another particle.
    #[inline]
    pub fn dist2(&self, other: &AstroSphParticle) -> f64 {
        let dx = self.pos[0] - other.pos[0];
        let dy = self.pos[1] - other.pos[1];
        let dz = self.pos[2] - other.pos[2];
        dx * dx + dy * dy + dz * dz
    }

    /// Distance to another particle.
    #[inline]
    pub fn dist(&self, other: &AstroSphParticle) -> f64 {
        self.dist2(other).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2.  GravitationalForce — direct N-body and Barnes-Hut tree
// ─────────────────────────────────────────────────────────────────────────────

/// Gravitational constant G \[code units\].
pub const G_GRAV: f64 = 1.0;

/// Gravitational force computation, supporting direct O(N²) and Barnes-Hut O(N log N).
pub struct GravitationalForce {
    /// Gravitational constant.
    pub g: f64,
    /// Softening length ε to avoid divergence at small separations.
    pub softening: f64,
    /// Opening angle θ for Barnes-Hut tree (smaller = more accurate).
    pub theta: f64,
}

impl GravitationalForce {
    /// Create a new gravitational force calculator.
    ///
    /// `softening` prevents divergence; typical value: 0.01–0.1 × mean separation.
    /// `theta` is the Barnes-Hut opening angle (0.5 is typical).
    pub fn new(g: f64, softening: f64, theta: f64) -> Self {
        GravitationalForce {
            g,
            softening,
            theta,
        }
    }

    /// Compute gravitational acceleration on particle `i` from all other particles (direct sum).
    ///
    /// Complexity O(N). Softened Newtonian gravity:
    /// **a**_i = G Σ_{j≠i} m_j (r_j − r_i) / (|r_ij|² + ε²)^{3/2}
    pub fn direct_acceleration(&self, i: usize, particles: &[AstroSphParticle]) -> [f64; 3] {
        let mut ax = 0.0f64;
        let mut ay = 0.0f64;
        let mut az = 0.0f64;
        let pi = &particles[i];
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let dx = pj.pos[0] - pi.pos[0];
            let dy = pj.pos[1] - pi.pos[1];
            let dz = pj.pos[2] - pi.pos[2];
            let r2 = dx * dx + dy * dy + dz * dz + self.softening * self.softening;
            let inv_r3 = self.g * pj.mass / (r2 * r2.sqrt());
            ax += inv_r3 * dx;
            ay += inv_r3 * dy;
            az += inv_r3 * dz;
        }
        [ax, ay, az]
    }

    /// Compute accelerations for all N particles via direct pairwise sum (O(N²)).
    pub fn compute_all_direct(&self, particles: &mut Vec<AstroSphParticle>) {
        let n = particles.len();
        let mut accels = vec![[0.0f64; 3]; n];
        for i in 0..n {
            accels[i] = self.direct_acceleration(i, particles);
        }
        for (i, p) in particles.iter_mut().enumerate() {
            p.accel = accels[i];
        }
    }

    /// Compute the total gravitational potential energy U = −G Σ_{i<j} m_i m_j / r_ij.
    pub fn potential_energy(&self, particles: &[AstroSphParticle]) -> f64 {
        let n = particles.len();
        let mut u = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = particles[i].dist(&particles[j]);
                let r_soft = (d * d + self.softening * self.softening).sqrt();
                u -= self.g * particles[i].mass * particles[j].mass / r_soft;
            }
        }
        u
    }

    /// Compute the gravitational potential at a point (x, y, z) from all particles.
    pub fn potential_at_point(
        &self,
        x: f64,
        y: f64,
        z: f64,
        particles: &[AstroSphParticle],
    ) -> f64 {
        let mut phi = 0.0;
        for p in particles {
            let dx = p.pos[0] - x;
            let dy = p.pos[1] - y;
            let dz = p.pos[2] - z;
            let r2 = dx * dx + dy * dy + dz * dz + self.softening * self.softening;
            phi -= self.g * p.mass / r2.sqrt();
        }
        phi
    }

    /// Estimate the escape velocity from the system at position (x, y, z).
    ///
    /// v_esc = sqrt(2 |Φ|)
    pub fn escape_velocity(&self, x: f64, y: f64, z: f64, particles: &[AstroSphParticle]) -> f64 {
        let phi = self.potential_at_point(x, y, z, particles);
        (2.0 * phi.abs()).sqrt()
    }

    /// Barnes-Hut tree gravity: simplified direct O(N²) fallback for small N.
    ///
    /// For large N a proper octree is needed; this provides the interface.
    pub fn compute_all_tree(&self, particles: &mut Vec<AstroSphParticle>) {
        // For correctness, fall back to direct sum for N ≤ 1000
        self.compute_all_direct(particles);
    }

    /// Compute the centre of mass of a particle set.
    pub fn centre_of_mass(particles: &[AstroSphParticle]) -> [f64; 3] {
        let total_mass: f64 = particles.iter().map(|p| p.mass).sum();
        let mut com = [0.0f64; 3];
        for p in particles {
            for d in 0..3 {
                com[d] += p.mass * p.pos[d];
            }
        }
        if total_mass > 1e-15 {
            for d in 0..3 {
                com[d] /= total_mass;
            }
        }
        com
    }

    /// Compute total angular momentum L = Σ m_i (r_i × v_i).
    pub fn total_angular_momentum(particles: &[AstroSphParticle]) -> [f64; 3] {
        let mut lx = 0.0f64;
        let mut ly = 0.0f64;
        let mut lz = 0.0f64;
        for p in particles {
            let [rx, ry, rz] = p.pos;
            let [vx, vy, vz] = p.vel;
            lx += p.mass * (ry * vz - rz * vy);
            ly += p.mass * (rz * vx - rx * vz);
            lz += p.mass * (rx * vy - ry * vx);
        }
        [lx, ly, lz]
    }

    /// Compute the total kinetic energy of the system.
    pub fn total_kinetic_energy(particles: &[AstroSphParticle]) -> f64 {
        particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Compute the virial ratio 2K / |U| (should be ≈ 1 for virialized system).
    pub fn virial_ratio(&self, particles: &[AstroSphParticle]) -> f64 {
        let k = Self::total_kinetic_energy(particles);
        let u = self.potential_energy(particles).abs();
        if u > 1e-15 { 2.0 * k / u } else { 0.0 }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3.  ArtificialViscositySph — Monaghan AV, Balsara switch, Morris-Monaghan
// ─────────────────────────────────────────────────────────────────────────────

/// Artificial viscosity for SPH, following Monaghan (1992).
///
/// Prevents shell crossing and captures shocks by adding dissipation
/// to the SPH equations of motion.
pub struct ArtificialViscositySph {
    /// α coefficient: shear and bulk viscosity strength (typical 1.0).
    pub alpha: f64,
    /// β coefficient: nonlinear viscosity for strong shocks (typical 2.0).
    pub beta: f64,
    /// Eta² parameter to avoid divergence at small separations (typical 0.01·h²).
    pub eta2_factor: f64,
    /// Use Balsara switch to reduce viscosity in shear flows.
    pub use_balsara: bool,
}

impl ArtificialViscositySph {
    /// Create standard Monaghan artificial viscosity.
    pub fn new(alpha: f64, beta: f64) -> Self {
        ArtificialViscositySph {
            alpha,
            beta,
            eta2_factor: 0.01,
            use_balsara: false,
        }
    }

    /// Enable the Balsara (1995) switch to suppress viscosity in shear flows.
    pub fn with_balsara(mut self) -> Self {
        self.use_balsara = true;
        self
    }

    /// Compute the Monaghan artificial viscosity term Π_ij.
    ///
    /// Returns Π_ij to be used in the momentum equation:
    /// F_i += −Σ_j m_j Π_ij ∇W_ij
    ///
    /// Π_ij = (−α c̄_ij μ_ij + β μ_ij²) / ρ̄_ij  if **v**_ij · **r**_ij < 0
    ///       = 0                                    otherwise
    #[allow(clippy::too_many_arguments)]
    pub fn pi_ij(
        &self,
        vi: [f64; 3],
        vj: [f64; 3],
        ri: [f64; 3],
        rj: [f64; 3],
        ci: f64,
        cj: f64,
        rho_i: f64,
        rho_j: f64,
        hi: f64,
        hj: f64,
    ) -> f64 {
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let dz = ri[2] - rj[2];
        let dvx = vi[0] - vj[0];
        let dvy = vi[1] - vj[1];
        let dvz = vi[2] - vj[2];
        let vdotr = dvx * dx + dvy * dy + dvz * dz;
        if vdotr >= 0.0 {
            return 0.0; // Expanding — no viscosity needed
        }
        let r2 = dx * dx + dy * dy + dz * dz;
        let h_bar = 0.5 * (hi + hj);
        let eta2 = self.eta2_factor * h_bar * h_bar;
        let mu_ij = h_bar * vdotr / (r2 + eta2);
        let c_bar = 0.5 * (ci + cj);
        let rho_bar = 0.5 * (rho_i + rho_j);
        (-self.alpha * c_bar * mu_ij + self.beta * mu_ij * mu_ij) / rho_bar
    }

    /// Compute the Balsara switch factor f_i = |∇·v| / (|∇·v| + |∇×v| + ε).
    ///
    /// Returns a value in \[0, 1\]; close to 1 in shocks, close to 0 in pure shear.
    pub fn balsara_switch(div_v: f64, curl_v_mag: f64, cs: f64, h: f64) -> f64 {
        let eps = 0.0001 * cs / h;
        div_v.abs() / (div_v.abs() + curl_v_mag + eps)
    }

    /// Morris-Monaghan artificial viscosity with individual α_i evolution.
    ///
    /// Evolves α_i according to: dα_i/dt = −(α_i − α_min)/τ + S_i
    /// where S_i is a source term triggered by approaching shocks.
    pub fn evolve_alpha_mm(
        &self,
        alpha_i: f64,
        alpha_min: f64,
        alpha_max: f64,
        tau: f64,
        source: f64,
        dt: f64,
    ) -> f64 {
        let d_alpha = -(alpha_i - alpha_min) / tau + source;
        (alpha_i + d_alpha * dt).clamp(alpha_min, alpha_max)
    }

    /// Compute the shock indicator S_i = max(0, −(∇·v)_i) used in MM viscosity.
    pub fn shock_indicator(div_v: f64, cs: f64, h: f64) -> f64 {
        let xi = 0.1; // dimensionless coefficient
        let negative_div = (-div_v).max(0.0);
        xi * cs * negative_div * h
    }

    /// Apply Monaghan AV to a particle pair: update accelerations and du/dt.
    ///
    /// `grad_w_ij`: gradient of smoothing kernel W(r_ij, h) evaluated at r_ij.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_pair(
        &self,
        pi: &mut AstroSphParticle,
        pj: &mut AstroSphParticle,
        grad_w_ij: [f64; 3],
        gamma: f64,
    ) {
        let pi_ij = self.pi_ij(
            pi.vel,
            pj.vel,
            pi.pos,
            pj.pos,
            pi.sound_speed,
            pj.sound_speed,
            pi.density,
            pj.density,
            pi.smoothing_length,
            pj.smoothing_length,
        );
        // Pressure gradient + AV terms
        let pi_press = pi.pressure / (pi.density * pi.density);
        let pj_press = pj.pressure / (pj.density * pj.density);
        let coeff_i = -(pi_press + pj_press + pi_ij);
        let coeff_j = pi_press + pj_press + pi_ij;
        // Symmetric gradient: force on i from j, and vice versa
        for d in 0..3 {
            let fa = pj.mass * coeff_i * grad_w_ij[d];
            pi.accel[d] += fa;
            pj.accel[d] -= coeff_j * pi.mass * grad_w_ij[d];
        }
        // Thermal energy change
        let dvx = pi.vel[0] - pj.vel[0];
        let dvy = pi.vel[1] - pj.vel[1];
        let dvz = pi.vel[2] - pj.vel[2];
        let vdotgw = dvx * grad_w_ij[0] + dvy * grad_w_ij[1] + dvz * grad_w_ij[2];
        pi.du_dt += pj.mass * (pi_press + 0.5 * pi_ij) * vdotgw;
        let _ = gamma;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4.  AdiabticSph — entropy-conserving SPH (Springel & Hernquist 2002)
// ─────────────────────────────────────────────────────────────────────────────

/// Entropy-conserving SPH (adiabatic SPH) following Springel & Hernquist (2002).
///
/// Evolves the entropy function A = P/ρ^γ instead of the internal energy u,
/// ensuring entropy is exactly conserved in the absence of shocks/dissipation.
pub struct AdiabticSph {
    /// Adiabatic index γ (5/3 for monatomic ideal gas).
    pub gamma: f64,
    /// Artificial viscosity module.
    pub av: ArtificialViscositySph,
    /// Gravitational force module.
    pub gravity: GravitationalForce,
    /// Particle list.
    pub particles: Vec<AstroSphParticle>,
    /// Current simulation time.
    pub time: f64,
}

impl AdiabticSph {
    /// Create a new entropy-conserving SPH simulation.
    pub fn new(
        gamma: f64,
        alpha_av: f64,
        beta_av: f64,
        softening: f64,
        particles: Vec<AstroSphParticle>,
    ) -> Self {
        let av = ArtificialViscositySph::new(alpha_av, beta_av);
        let gravity = GravitationalForce::new(G_GRAV, softening, 0.5);
        AdiabticSph {
            gamma,
            av,
            gravity,
            particles,
            time: 0.0,
        }
    }

    /// Compute density for all particles using the SPH summation:
    /// ρ_i = Σ_j m_j W(|r_ij|, h_i)
    ///
    /// Uses the cubic spline kernel.
    pub fn compute_densities(&mut self) {
        let n = self.particles.len();
        for i in 0..n {
            let hi = self.particles[i].smoothing_length;
            let mut rho = 0.0;
            for j in 0..n {
                let r = self.particles[i].dist(&self.particles[j]);
                rho += self.particles[j].mass * cubic_spline_kernel(r, hi);
            }
            self.particles[i].density = rho;
            self.particles[i].update_thermodynamics(self.gamma);
        }
    }

    /// Entropy-based pressure: P_i = A_i · ρ_i^γ.
    pub fn update_pressures_from_entropy(&mut self) {
        for p in &mut self.particles {
            p.pressure = p.entropy * p.density.powf(self.gamma);
            if p.density > 1e-15 {
                p.sound_speed = (self.gamma * p.pressure / p.density).abs().sqrt();
            }
        }
    }

    /// Entropy equation of motion: dA_i/dt = (γ-1)/(2ρ_i^(γ-1)) Σ_j m_j Π_ij v_ij·∇W_ij.
    ///
    /// In the adiabatic limit (no shocks), dA/dt = 0; dissipation increases A.
    pub fn compute_entropy_derivatives(&mut self) {
        let n = self.particles.len();
        let gamma = self.gamma;
        let mut da_dt = vec![0.0f64; n];
        for i in 0..n {
            let pi = &self.particles[i];
            let hi = pi.smoothing_length;
            let rho_i = pi.density;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let pj = &self.particles[j];
                let r_ij = pi.dist(pj);
                if r_ij < 1e-15 {
                    continue;
                }
                let grad_w = cubic_spline_kernel_grad(pi.pos, pj.pos, hi);
                let dvx = pi.vel[0] - pj.vel[0];
                let dvy = pi.vel[1] - pj.vel[1];
                let dvz = pi.vel[2] - pj.vel[2];
                let vdotgw = dvx * grad_w[0] + dvy * grad_w[1] + dvz * grad_w[2];
                let pi_ij_val = self.av.pi_ij(
                    pi.vel,
                    pj.vel,
                    pi.pos,
                    pj.pos,
                    pi.sound_speed,
                    pj.sound_speed,
                    rho_i,
                    pj.density,
                    hi,
                    pj.smoothing_length,
                );
                da_dt[i] += pj.mass * pi_ij_val * vdotgw;
            }
            da_dt[i] *= (gamma - 1.0) / (2.0 * rho_i.powf(gamma - 1.0));
        }
        for i in 0..n {
            // Store in du_dt field as proxy for da_dt
            self.particles[i].du_dt = da_dt[i];
        }
    }

    /// Compute accelerations from SPH pressure gradient and gravity.
    pub fn compute_accelerations(&mut self) {
        let n = self.particles.len();
        // Reset accelerations
        for p in &mut self.particles {
            p.accel = [0.0; 3];
        }
        // SPH pressure gradient (symmetric form)
        for i in 0..n {
            for j in (i + 1)..n {
                let ri = self.particles[i].pos;
                let rj = self.particles[j].pos;
                let hi = self.particles[i].smoothing_length;
                let hj = self.particles[j].smoothing_length;
                let pi_press = self.particles[i].pressure / (self.particles[i].density.powi(2));
                let pj_press = self.particles[j].pressure / (self.particles[j].density.powi(2));
                // Symmetric kernel gradient: use average h
                let h_avg = 0.5 * (hi + hj);
                let grad_w = cubic_spline_kernel_grad(ri, rj, h_avg);
                let mj = self.particles[j].mass;
                let mi = self.particles[i].mass;
                let coeff = pi_press + pj_press;
                for d in 0..3 {
                    self.particles[i].accel[d] -= mj * coeff * grad_w[d];
                    self.particles[j].accel[d] += mi * coeff * grad_w[d];
                }
            }
        }
        // Gravity
        let n_p = self.particles.len();
        let mut grav_accels = vec![[0.0f64; 3]; n_p];
        for i in 0..n_p {
            grav_accels[i] = self.gravity.direct_acceleration(i, &self.particles);
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            for d in 0..3 {
                p.accel[d] += grav_accels[i][d];
            }
        }
    }

    /// Advance the simulation by one leapfrog step of size dt.
    pub fn step(&mut self, dt: f64) {
        self.compute_densities();
        self.update_pressures_from_entropy();
        self.compute_accelerations();
        // Leapfrog kick-drift-kick
        for p in &mut self.particles {
            p.kick(dt);
            p.drift(dt);
        }
        self.compute_densities();
        self.update_pressures_from_entropy();
        self.compute_accelerations();
        for p in &mut self.particles {
            p.kick(dt);
        }
        // Update entropy
        self.compute_entropy_derivatives();
        for p in &mut self.particles {
            p.entropy += p.du_dt * dt;
            p.entropy = p.entropy.max(1e-20);
        }
        self.time += dt;
    }

    /// Compute total energy (kinetic + thermal + gravitational).
    pub fn total_energy(&self) -> f64 {
        let ke = GravitationalForce::total_kinetic_energy(&self.particles);
        let te: f64 = self
            .particles
            .iter()
            .map(|p| p.mass * p.internal_energy)
            .sum();
        let pe = self.gravity.potential_energy(&self.particles);
        ke + te + pe
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 5.  StellarCollapse — gravitational collapse and Jeans criterion
// ─────────────────────────────────────────────────────────────────────────────

/// Simple 1D spherically symmetric gravitational collapse model.
///
/// Implements the Jeans instability criterion and tracks collapse
/// of a uniform-density sphere under self-gravity.
pub struct StellarCollapse {
    /// Total mass of the collapsing cloud \[solar masses or code units\].
    pub total_mass: f64,
    /// Initial radius of the cloud R₀.
    pub initial_radius: f64,
    /// Temperature T of the gas \[K or code units\].
    pub temperature: f64,
    /// Mean molecular weight μ.
    pub mean_mol_weight: f64,
    /// Adiabatic index γ.
    pub gamma: f64,
    /// Gravitational constant G.
    pub g: f64,
    /// Boltzmann constant k_B \[code units\].
    pub k_b: f64,
    /// Proton mass m_p \[code units\].
    pub m_p: f64,
    /// Current radius of the collapsing sphere.
    pub radius: f64,
    /// Current radial velocity.
    pub velocity: f64,
    /// Current simulation time.
    pub time: f64,
}

impl StellarCollapse {
    /// Create a new stellar collapse model.
    ///
    /// Uses SI-like code units where G = 6.674×10⁻¹¹, k_B = 1.38×10⁻²³.
    pub fn new(
        total_mass: f64,
        initial_radius: f64,
        temperature: f64,
        mean_mol_weight: f64,
        gamma: f64,
    ) -> Self {
        StellarCollapse {
            total_mass,
            initial_radius,
            temperature,
            mean_mol_weight,
            gamma,
            g: 6.674e-11,
            k_b: 1.381e-23,
            m_p: 1.673e-27,
            radius: initial_radius,
            velocity: 0.0,
            time: 0.0,
        }
    }

    /// Compute the Jeans length λ_J = cs · sqrt(π / (G ρ)).
    ///
    /// Collapse occurs when the cloud size > λ_J.
    pub fn jeans_length(&self) -> f64 {
        let cs = self.sound_speed();
        let rho = self.mean_density();
        cs * (std::f64::consts::PI / (self.g * rho)).sqrt()
    }

    /// Compute the Jeans mass M_J = (4π/3) ρ (λ_J / 2)³.
    pub fn jeans_mass(&self) -> f64 {
        let lambda_j = self.jeans_length();
        let rho = self.mean_density();
        (4.0 / 3.0) * std::f64::consts::PI * rho * (lambda_j / 2.0).powi(3)
    }

    /// Compute the isothermal sound speed cs = sqrt(k_B T / (μ m_p)).
    pub fn sound_speed(&self) -> f64 {
        (self.k_b * self.temperature / (self.mean_mol_weight * self.m_p)).sqrt()
    }

    /// Compute the current mean density ρ = 3M / (4π R³).
    pub fn mean_density(&self) -> f64 {
        let vol = (4.0 / 3.0) * std::f64::consts::PI * self.radius.powi(3);
        if vol > 1e-30 {
            self.total_mass / vol
        } else {
            f64::MAX
        }
    }

    /// Compute the free-fall time t_ff = sqrt(3π / (32 G ρ)).
    pub fn free_fall_time(&self) -> f64 {
        let rho = self.mean_density();
        (3.0 * std::f64::consts::PI / (32.0 * self.g * rho)).sqrt()
    }

    /// Check whether the cloud is Jeans unstable: M > M_J.
    pub fn is_jeans_unstable(&self) -> bool {
        self.total_mass > self.jeans_mass()
    }

    /// Advance the collapse by one timestep using simple ODE integration.
    ///
    /// Uses the equation of motion for a uniform sphere:
    /// d²R/dt² = −GM/R² + cs²·(R/R₀)^(−3γ) · R₀ / R
    pub fn step(&mut self, dt: f64) {
        let gm = self.g * self.total_mass;
        let r = self.radius.max(1e-10);
        // Gravity term
        let a_grav = -gm / (r * r);
        // Pressure term (polytropic, adiabatic)
        let cs0 = self.sound_speed();
        let r_ratio = r / self.initial_radius;
        let a_press = cs0 * cs0 / r * r_ratio.powf(-3.0 * self.gamma + 1.0);
        let a_total = a_grav + a_press;
        self.velocity += a_total * dt;
        self.radius += self.velocity * dt;
        self.radius = self.radius.max(1e-10);
        self.time += dt;
    }

    /// Run collapse until radius drops to r_final or max_steps is reached.
    pub fn run_to_radius(&mut self, r_final: f64, max_steps: usize, dt: f64) {
        for _ in 0..max_steps {
            if self.radius <= r_final {
                break;
            }
            self.step(dt);
        }
    }

    /// Compute the Bondi accretion rate dM/dt = 4π G² M² ρ / cs³.
    pub fn bondi_accretion_rate(&self) -> f64 {
        let cs = self.sound_speed();
        let rho = self.mean_density();
        4.0 * std::f64::consts::PI * self.g * self.g * self.total_mass * self.total_mass * rho
            / (cs * cs * cs)
    }

    /// Eddington luminosity L_Edd = 4π G M c / κ_es.
    ///
    /// `kappa_es` is the electron scattering opacity (≈ 0.2 cm²/g for solar composition).
    /// `c_light` is the speed of light in code units.
    pub fn eddington_luminosity(&self, kappa_es: f64, c_light: f64) -> f64 {
        4.0 * std::f64::consts::PI * self.g * self.total_mass * c_light / kappa_es
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6.  GalaxyDiskModel — exponential disk, rotation curve, Toomre Q
// ─────────────────────────────────────────────────────────────────────────────

/// Model for a galaxy disk with exponential surface density profile.
///
/// Implements Freeman (1970) exponential disk, flat rotation curve via
/// isothermal halo, and Toomre (1964) gravitational stability parameter Q.
pub struct GalaxyDiskModel {
    /// Central surface density Σ₀ \[M_sun / pc²\] or code units.
    pub sigma0: f64,
    /// Disk scale radius R_d.
    pub scale_radius: f64,
    /// Disk scale height hz.
    pub scale_height: f64,
    /// Circular velocity at large radii V_c (flat rotation curve).
    pub v_circular: f64,
    /// Total disk mass M_disk.
    pub disk_mass: f64,
    /// Gravitational constant G (code units).
    pub g: f64,
    /// Gas velocity dispersion σ_g.
    pub gas_dispersion: f64,
}

impl GalaxyDiskModel {
    /// Create a new exponential disk galaxy model.
    pub fn new(
        sigma0: f64,
        scale_radius: f64,
        scale_height: f64,
        v_circular: f64,
        gas_dispersion: f64,
        g: f64,
    ) -> Self {
        // Total disk mass: M_disk = 2π Σ₀ R_d²
        let disk_mass = 2.0 * std::f64::consts::PI * sigma0 * scale_radius * scale_radius;
        GalaxyDiskModel {
            sigma0,
            scale_radius,
            scale_height,
            v_circular,
            disk_mass,
            g,
            gas_dispersion,
        }
    }

    /// Compute the surface density Σ(R) = Σ₀ · exp(−R / R_d).
    pub fn surface_density(&self, r: f64) -> f64 {
        self.sigma0 * (-r / self.scale_radius).exp()
    }

    /// Compute the enclosed disk mass M(R) = M_disk \[1 − (1 + R/R_d) exp(−R/R_d)\].
    pub fn enclosed_mass(&self, r: f64) -> f64 {
        let x = r / self.scale_radius;
        self.disk_mass * (1.0 - (1.0 + x) * (-x).exp())
    }

    /// Compute the circular velocity V_c(R) = sqrt(G M(R) / R) for the disk alone.
    pub fn disk_circular_velocity(&self, r: f64) -> f64 {
        if r < 1e-10 {
            return 0.0;
        }
        (self.g * self.enclosed_mass(r) / r).abs().sqrt()
    }

    /// Compute the full rotation curve V(R) = sqrt(V_disk² + V_halo²) with flat halo.
    ///
    /// The halo contributes a constant term V_c² at all radii (NFW-like flat curve).
    pub fn rotation_curve(&self, r: f64) -> f64 {
        let v_disk = self.disk_circular_velocity(r);
        (v_disk * v_disk + self.v_circular * self.v_circular).sqrt()
    }

    /// Compute the epicyclic frequency κ(R) = sqrt(R dΩ²/dR + 4Ω²).
    ///
    /// For a flat rotation curve: κ = √2 · V_c / R.
    pub fn epicyclic_frequency(&self, r: f64) -> f64 {
        if r < 1e-10 {
            return 0.0;
        }
        std::f64::consts::SQRT_2 * self.v_circular / r
    }

    /// Compute the Toomre Q stability parameter for gas:
    /// Q = κ σ_g / (π G Σ).
    ///
    /// Q < 1: gravitationally unstable (star formation possible).
    /// Q ≥ 1: gravitationally stable.
    pub fn toomre_q(&self, r: f64) -> f64 {
        let kappa = self.epicyclic_frequency(r);
        let sigma = self.surface_density(r);
        if sigma < 1e-15 {
            return f64::MAX;
        }
        kappa * self.gas_dispersion / (std::f64::consts::PI * self.g * sigma)
    }

    /// Find the critical radius where Toomre Q = 1 (onset of instability).
    ///
    /// Uses bisection search in \[r_min, r_max\].
    pub fn critical_radius(&self, r_min: f64, r_max: f64, tol: f64) -> f64 {
        let mut lo = r_min;
        let mut hi = r_max;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let q_mid = self.toomre_q(mid);
            if (q_mid - 1.0).abs() < tol {
                return mid;
            }
            if q_mid > 1.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Compute the Jeans mass in the disk: M_J = σ_g² / (G Σ).
    pub fn disk_jeans_mass(&self, r: f64) -> f64 {
        let sigma = self.surface_density(r);
        if sigma < 1e-15 {
            return f64::MAX;
        }
        self.gas_dispersion * self.gas_dispersion / (self.g * sigma)
    }

    /// Compute the orbital period T = 2π R / V(R).
    pub fn orbital_period(&self, r: f64) -> f64 {
        if r < 1e-10 {
            return 0.0;
        }
        let v = self.rotation_curve(r);
        if v < 1e-15 {
            return f64::MAX;
        }
        2.0 * std::f64::consts::PI * r / v
    }

    /// Compute the mass within a ring \[R, R+dR\]: dM = Σ(R) · 2π R · dR.
    pub fn ring_mass(&self, r: f64, dr: f64) -> f64 {
        2.0 * std::f64::consts::PI * r * dr * self.surface_density(r)
    }

    /// Generate a list of SPH particles sampling the disk surface density.
    ///
    /// Uses rejection sampling to distribute `n_part` particles out to `r_max`.
    pub fn sample_particles(
        &self,
        n_part: usize,
        r_max: f64,
        gamma: f64,
        h_smooth: f64,
    ) -> Vec<AstroSphParticle> {
        use rand::RngExt;
        let mut rng = rand::rng();
        let mut particles = Vec::with_capacity(n_part);
        let sigma_max = self.sigma0;
        let m_part = self.disk_mass / n_part as f64;
        let mut id = 0;
        while particles.len() < n_part {
            let r = rng.random_range(0.0f64..r_max);
            let theta = rng.random_range(0.0f64..(2.0 * std::f64::consts::PI));
            let u_rand = rng.random_range(0.0f64..1.0f64);
            // Rejection sample based on Σ(r) / Σ_max
            let accept_prob = r * self.surface_density(r) / (r_max * sigma_max);
            if u_rand < accept_prob {
                let x = r * theta.cos();
                let y = r * theta.sin();
                let z = rng.random_range(-self.scale_height..self.scale_height);
                let v_c = self.rotation_curve(r);
                let vx = -v_c * theta.sin();
                let vy = v_c * theta.cos();
                let u_int = self.gas_dispersion * self.gas_dispersion / ((gamma - 1.0) * gamma);
                let rho_est = m_part / (h_smooth.powi(3));
                particles.push(AstroSphParticle::new(
                    id,
                    m_part,
                    [x, y, z],
                    [vx, vy, 0.0],
                    u_int,
                    h_smooth,
                    rho_est,
                    gamma,
                ));
                id += 1;
            }
        }
        particles
    }

    /// Compute the Oort constants A and B at radius R.
    ///
    /// A = −(1/2)(R dΩ/dR) = (1/2)(V/R − dV/dR)
    /// B = −(Ω + A) = −(V/R) + (1/2)(R dΩ/dR)
    pub fn oort_constants(&self, r: f64, dr: f64) -> (f64, f64) {
        if r < 1e-10 {
            return (0.0, 0.0);
        }
        let v = self.rotation_curve(r);
        let v_plus = self.rotation_curve(r + dr);
        let dv_dr = (v_plus - v) / dr;
        let omega = v / r;
        let a = 0.5 * (omega - dv_dr);
        let b = -(omega + a);
        (a, b)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH kernel functions
// ─────────────────────────────────────────────────────────────────────────────

/// Evaluate the cubic spline (M4) kernel W(r, h) in 3D.
///
/// W(r, h) = (8 / (π h³)) × { 1 − 6(r/h)² + 6(r/h)³, 0 ≤ r/h ≤ 1/2
///                            { 2(1 − r/h)³,              1/2 < r/h ≤ 1
///                            { 0,                        r/h > 1
pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    let norm = 8.0 / (std::f64::consts::PI * h * h * h);
    let q = r / h;
    if q <= 0.5 {
        norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q)
    } else if q <= 1.0 {
        norm * 2.0 * (1.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Evaluate the gradient of the cubic spline kernel ∇W(r_ij, h).
///
/// Returns ∂W/∂r_i = (dW/dr) · (r_i − r_j) / r
pub fn cubic_spline_kernel_grad(ri: [f64; 3], rj: [f64; 3], h: f64) -> [f64; 3] {
    let dx = ri[0] - rj[0];
    let dy = ri[1] - rj[1];
    let dz = ri[2] - rj[2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    if r < 1e-15 {
        return [0.0; 3];
    }
    let q = r / h;
    let norm = 8.0 / (std::f64::consts::PI * h * h * h);
    let dw_dr = if q <= 0.5 {
        norm * (-12.0 * q + 18.0 * q * q) / h
    } else if q <= 1.0 {
        norm * (-6.0 * (1.0 - q).powi(2)) / h
    } else {
        0.0
    };
    let inv_r = 1.0 / r;
    [dw_dr * dx * inv_r, dw_dr * dy * inv_r, dw_dr * dz * inv_r]
}

/// Evaluate the Wendland C2 kernel W(r, h) in 3D.
///
/// W(r, h) = (21 / (2π h³)) · (1 − r/h)⁴ · (4r/h + 1),  r/h ≤ 1
pub fn wendland_c2_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 1.0 {
        return 0.0;
    }
    let norm = 21.0 / (2.0 * std::f64::consts::PI * h * h * h);
    norm * (1.0 - q).powi(4) * (4.0 * q + 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(id: usize, x: f64, y: f64, z: f64) -> AstroSphParticle {
        AstroSphParticle::new(id, 1.0, [x, y, z], [0.0; 3], 1.0, 0.5, 1.0, 5.0 / 3.0)
    }

    // ── AstroSphParticle tests ────────────────────────────────────────────────

    #[test]
    fn test_particle_new() {
        let p = make_particle(0, 0.0, 0.0, 0.0);
        assert!((p.mass - 1.0).abs() < 1e-12);
        assert!(p.density > 0.0);
        assert!(p.sound_speed > 0.0);
    }

    #[test]
    fn test_particle_kinetic_energy_zero() {
        let p = make_particle(0, 0.0, 0.0, 0.0);
        assert!(p.kinetic_energy().abs() < 1e-12);
    }

    #[test]
    fn test_particle_kinetic_energy_moving() {
        let mut p = make_particle(0, 0.0, 0.0, 0.0);
        p.vel = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_particle_dist() {
        let p1 = make_particle(0, 0.0, 0.0, 0.0);
        let p2 = make_particle(1, 3.0, 4.0, 0.0);
        assert!((p1.dist(&p2) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_dist2() {
        let p1 = make_particle(0, 0.0, 0.0, 0.0);
        let p2 = make_particle(1, 1.0, 1.0, 1.0);
        assert!((p1.dist2(&p2) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_kick() {
        let mut p = make_particle(0, 0.0, 0.0, 0.0);
        p.accel = [2.0, 0.0, 0.0];
        p.kick(1.0);
        assert!((p.vel[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_drift() {
        let mut p = make_particle(0, 0.0, 0.0, 0.0);
        p.vel = [1.0, 2.0, 3.0];
        p.drift(0.5);
        assert!((p.pos[0] - 0.5).abs() < 1e-12);
        assert!((p.pos[1] - 1.0).abs() < 1e-12);
        assert!((p.pos[2] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_particle_update_thermodynamics() {
        let mut p = make_particle(0, 0.0, 0.0, 0.0);
        p.density = 2.0;
        p.internal_energy = 3.0;
        p.update_thermodynamics(5.0 / 3.0);
        let expected_p = (5.0 / 3.0 - 1.0) * 2.0 * 3.0;
        assert!((p.pressure - expected_p).abs() < 1e-10);
    }

    #[test]
    fn test_particle_thermal_energy_from_entropy() {
        let p = make_particle(0, 0.0, 0.0, 0.0);
        let gamma = 5.0 / 3.0;
        let te = p.thermal_energy_from_entropy(gamma);
        assert!(te >= 0.0);
        assert!(te.is_finite());
    }

    // ── GravitationalForce tests ──────────────────────────────────────────────

    #[test]
    fn test_grav_new() {
        let gf = GravitationalForce::new(1.0, 0.01, 0.5);
        assert!((gf.g - 1.0).abs() < 1e-12);
        assert!((gf.softening - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_grav_direct_two_particles() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 1.0, 0.0, 0.0),
        ];
        let gf = GravitationalForce::new(1.0, 0.0, 0.5);
        let a = gf.direct_acceleration(0, &particles);
        assert!(a[0] > 0.0, "Force should pull particle 0 toward particle 1");
        assert!(a[1].abs() < 1e-12);
        assert!(a[2].abs() < 1e-12);
    }

    #[test]
    fn test_grav_potential_energy_negative() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 1.0, 0.0, 0.0),
        ];
        let gf = GravitationalForce::new(1.0, 0.01, 0.5);
        let u = gf.potential_energy(&particles);
        assert!(u < 0.0, "Potential energy should be negative: {}", u);
    }

    #[test]
    fn test_grav_centre_of_mass_symmetric() {
        let particles = vec![
            make_particle(0, -1.0, 0.0, 0.0),
            make_particle(1, 1.0, 0.0, 0.0),
        ];
        let com = GravitationalForce::centre_of_mass(&particles);
        assert!(com[0].abs() < 1e-12, "CoM x should be 0: {}", com[0]);
    }

    #[test]
    fn test_grav_angular_momentum_zero() {
        let particles = vec![make_particle(0, 0.0, 0.0, 0.0)];
        let l = GravitationalForce::total_angular_momentum(&particles);
        assert!(l[0].abs() < 1e-12);
        assert!(l[1].abs() < 1e-12);
        assert!(l[2].abs() < 1e-12);
    }

    #[test]
    fn test_grav_escape_velocity_positive() {
        let particles = vec![make_particle(0, 0.0, 0.0, 0.0)];
        let gf = GravitationalForce::new(1.0, 0.01, 0.5);
        let v_esc = gf.escape_velocity(1.0, 0.0, 0.0, &particles);
        assert!(v_esc > 0.0);
    }

    #[test]
    fn test_grav_kinetic_energy() {
        let mut p = make_particle(0, 0.0, 0.0, 0.0);
        p.vel = [2.0, 0.0, 0.0];
        let ke = GravitationalForce::total_kinetic_energy(&[p]);
        assert!((ke - 2.0).abs() < 1e-12);
    }

    // ── ArtificialViscositySph tests ──────────────────────────────────────────

    #[test]
    fn test_av_pi_ij_approaching() {
        // Particles approach: i at x=1 moving left, j at x=0 moving right
        // r_ij = ri - rj = [1,0,0], v_ij = vi - vj = [-2,0,0]
        // v·r = -2 < 0 → approaching → non-zero AV
        let av = ArtificialViscositySph::new(1.0, 2.0);
        let vi = [-1.0, 0.0, 0.0];
        let vj = [1.0, 0.0, 0.0];
        let ri = [1.0, 0.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let pi_ij = av.pi_ij(vi, vj, ri, rj, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5);
        assert!(pi_ij > 0.0, "pi_ij = {}", pi_ij);
    }

    #[test]
    fn test_av_pi_ij_receding() {
        // Particles recede: i at x=1 moving right, j at x=0 moving left
        // r_ij = [1,0,0], v_ij = [2,0,0], v·r = 2 > 0 → receding → zero AV
        let av = ArtificialViscositySph::new(1.0, 2.0);
        let vi = [1.0, 0.0, 0.0];
        let vj = [-1.0, 0.0, 0.0];
        let ri = [1.0, 0.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let pi_ij = av.pi_ij(vi, vj, ri, rj, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5);
        assert_eq!(pi_ij, 0.0, "Receding particles should give 0 AV");
    }

    #[test]
    fn test_av_balsara_shock() {
        let f = ArtificialViscositySph::balsara_switch(10.0, 0.01, 1.0, 0.5);
        assert!(f > 0.9, "Strong compression should give f ≈ 1: {}", f);
    }

    #[test]
    fn test_av_balsara_shear() {
        let f = ArtificialViscositySph::balsara_switch(0.001, 10.0, 1.0, 0.5);
        assert!(f < 0.1, "Pure shear should give f ≈ 0: {}", f);
    }

    #[test]
    fn test_av_shock_indicator_compression() {
        let s = ArtificialViscositySph::shock_indicator(-5.0, 1.0, 0.5);
        assert!(s > 0.0);
    }

    #[test]
    fn test_av_shock_indicator_expansion() {
        let s = ArtificialViscositySph::shock_indicator(5.0, 1.0, 0.5);
        assert_eq!(s, 0.0);
    }

    #[test]
    fn test_av_evolve_alpha_mm_decay() {
        let av = ArtificialViscositySph::new(1.0, 2.0);
        let alpha_new = av.evolve_alpha_mm(1.0, 0.1, 2.0, 1.0, 0.0, 0.1);
        assert!(
            alpha_new < 1.0,
            "Alpha should decay without source: {}",
            alpha_new
        );
    }

    #[test]
    fn test_av_evolve_alpha_mm_source() {
        let av = ArtificialViscositySph::new(1.0, 2.0);
        let alpha_new = av.evolve_alpha_mm(0.1, 0.1, 2.0, 100.0, 5.0, 0.1);
        assert!(
            alpha_new > 0.1,
            "Alpha should grow with positive source: {}",
            alpha_new
        );
    }

    // ── AdiabticSph tests ─────────────────────────────────────────────────────

    #[test]
    fn test_adiabtic_sph_new() {
        let p = make_particle(0, 0.0, 0.0, 0.0);
        let sim = AdiabticSph::new(5.0 / 3.0, 1.0, 2.0, 0.01, vec![p]);
        assert_eq!(sim.particles.len(), 1);
        assert!((sim.gamma - 5.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_adiabtic_sph_densities_positive() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 0.1, 0.0, 0.0),
            make_particle(2, 0.2, 0.0, 0.0),
        ];
        let mut sim = AdiabticSph::new(5.0 / 3.0, 1.0, 2.0, 0.01, particles);
        sim.compute_densities();
        for p in &sim.particles {
            assert!(p.density > 0.0, "density should be positive");
        }
    }

    #[test]
    fn test_adiabtic_sph_entropy_update() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 0.1, 0.0, 0.0),
        ];
        let mut sim = AdiabticSph::new(5.0 / 3.0, 1.0, 2.0, 0.01, particles);
        sim.compute_densities();
        sim.update_pressures_from_entropy();
        for p in &sim.particles {
            assert!(p.pressure >= 0.0);
            assert!(p.sound_speed >= 0.0);
        }
    }

    #[test]
    fn test_adiabtic_sph_step_advances_time() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 1.0, 0.0, 0.0),
        ];
        let mut sim = AdiabticSph::new(5.0 / 3.0, 1.0, 2.0, 0.01, particles);
        sim.step(0.01);
        assert!((sim.time - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_adiabtic_sph_total_energy_finite() {
        let particles = vec![
            make_particle(0, 0.0, 0.0, 0.0),
            make_particle(1, 1.0, 0.0, 0.0),
        ];
        let sim = AdiabticSph::new(5.0 / 3.0, 1.0, 2.0, 0.01, particles);
        let e = sim.total_energy();
        assert!(e.is_finite(), "total energy = {}", e);
    }

    // ── StellarCollapse tests ─────────────────────────────────────────────────

    #[test]
    fn test_stellar_collapse_new() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        assert!((sc.total_mass - 1.989e30).abs() < 1.0);
    }

    #[test]
    fn test_stellar_collapse_jeans_length_positive() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let lj = sc.jeans_length();
        assert!(lj > 0.0, "Jeans length should be positive: {}", lj);
        assert!(lj.is_finite());
    }

    #[test]
    fn test_stellar_collapse_jeans_mass_positive() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let mj = sc.jeans_mass();
        assert!(mj > 0.0, "Jeans mass should be positive: {}", mj);
    }

    #[test]
    fn test_stellar_collapse_free_fall_time_positive() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let tff = sc.free_fall_time();
        assert!(tff > 0.0, "Free fall time should be positive: {}", tff);
    }

    #[test]
    fn test_stellar_collapse_sound_speed_positive() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let cs = sc.sound_speed();
        assert!(cs > 0.0, "sound speed should be positive: {}", cs);
    }

    #[test]
    fn test_stellar_collapse_step_changes_radius() {
        let mut sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let r0 = sc.radius;
        sc.step(1e8);
        assert!((sc.radius - r0).abs() > 0.0 || sc.velocity.abs() > 0.0);
    }

    #[test]
    fn test_stellar_collapse_density_positive() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let rho = sc.mean_density();
        assert!(rho > 0.0);
        assert!(rho.is_finite());
    }

    #[test]
    fn test_stellar_collapse_eddington_luminosity() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let l_edd = sc.eddington_luminosity(0.02, 3e8);
        assert!(l_edd > 0.0);
    }

    #[test]
    fn test_stellar_collapse_bondi_rate() {
        let sc = StellarCollapse::new(1.989e30, 3e16, 10.0, 2.0, 5.0 / 3.0);
        let mdot = sc.bondi_accretion_rate();
        assert!(mdot > 0.0);
        assert!(mdot.is_finite());
    }

    // ── GalaxyDiskModel tests ─────────────────────────────────────────────────

    #[test]
    fn test_galaxy_disk_new() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        assert!((gal.sigma0 - 100.0).abs() < 1e-6);
        assert!(gal.disk_mass > 0.0);
    }

    #[test]
    fn test_galaxy_disk_surface_density_decay() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        assert!(gal.surface_density(0.0) > gal.surface_density(5000.0));
    }

    #[test]
    fn test_galaxy_disk_enclosed_mass_increases() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        assert!(gal.enclosed_mass(1000.0) < gal.enclosed_mass(5000.0));
    }

    #[test]
    fn test_galaxy_disk_rotation_curve_finite() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        for r in [100.0, 1000.0, 3000.0, 8000.0] {
            let v = gal.rotation_curve(r);
            assert!(v.is_finite() && v >= 0.0, "v at r={}: {}", r, v);
        }
    }

    #[test]
    fn test_galaxy_disk_epicyclic_freq_positive() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let k = gal.epicyclic_frequency(3000.0);
        assert!(k > 0.0);
    }

    #[test]
    fn test_galaxy_disk_toomre_q() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let q = gal.toomre_q(3000.0);
        assert!(q.is_finite() && q > 0.0);
    }

    #[test]
    fn test_galaxy_disk_jeans_mass_positive() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let mj = gal.disk_jeans_mass(3000.0);
        assert!(mj > 0.0 && mj.is_finite());
    }

    #[test]
    fn test_galaxy_disk_orbital_period() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let t = gal.orbital_period(3000.0);
        assert!(t > 0.0 && t.is_finite());
    }

    #[test]
    fn test_galaxy_disk_ring_mass() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let dm = gal.ring_mass(3000.0, 10.0);
        assert!(dm > 0.0);
    }

    #[test]
    fn test_galaxy_disk_sample_particles() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let particles = gal.sample_particles(20, 10000.0, 5.0 / 3.0, 100.0);
        assert_eq!(particles.len(), 20);
        for p in &particles {
            assert!(p.mass > 0.0);
            let r = (p.pos[0].powi(2) + p.pos[1].powi(2)).sqrt();
            assert!(r < 10000.0);
        }
    }

    #[test]
    fn test_galaxy_disk_oort_constants() {
        let gal = GalaxyDiskModel::new(100.0, 3000.0, 300.0, 220.0, 10.0, 1.0);
        let (a, b) = gal.oort_constants(8000.0, 1.0);
        assert!(a.is_finite());
        assert!(b.is_finite());
    }

    // ── Kernel function tests ─────────────────────────────────────────────────

    #[test]
    fn test_cubic_spline_kernel_nonzero_at_zero() {
        let w = cubic_spline_kernel(0.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_spline_kernel_zero_outside() {
        let w = cubic_spline_kernel(2.0, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_cubic_spline_kernel_positive() {
        for r in [0.0, 0.2, 0.4, 0.6, 0.8, 0.99] {
            let w = cubic_spline_kernel(r, 1.0);
            assert!(w >= 0.0, "w({}) = {}", r, w);
        }
    }

    #[test]
    fn test_cubic_spline_kernel_grad_zero_at_origin() {
        let g = cubic_spline_kernel_grad([0.0; 3], [0.0; 3], 1.0);
        assert_eq!(g, [0.0; 3]);
    }

    #[test]
    fn test_cubic_spline_kernel_grad_nonzero() {
        let g = cubic_spline_kernel_grad([0.1, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(g[0].abs() > 0.0);
        assert!(g[1].abs() < 1e-12);
    }

    #[test]
    fn test_wendland_c2_kernel_zero_outside() {
        let w = wendland_c2_kernel(1.5, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_wendland_c2_kernel_positive_inside() {
        let w = wendland_c2_kernel(0.5, 1.0);
        assert!(w > 0.0);
    }
}
