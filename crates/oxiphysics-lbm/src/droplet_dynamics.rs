// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Droplet dynamics for LBM multiphase flow.
//!
//! Implements droplet tracking, coalescence detection, breakup criteria,
//! and deformation analysis using the Weber number for droplets in LBM
//! phase-field or multiphase simulations.

use std::f64::consts::PI;

/// A single droplet in the simulation domain.
///
/// Stores the centroid, equivalent radius, velocity, and interface tension.
#[derive(Debug, Clone)]
pub struct Droplet {
    /// Droplet centroid \[x, y, z\] in lattice units.
    pub center: [f64; 3],
    /// Equivalent sphere radius in lattice units.
    pub radius: f64,
    /// Droplet velocity \[vx, vy, vz\] in lattice units per time step.
    pub velocity: [f64; 3],
    /// Surface tension coefficient σ (force per unit length).
    pub surface_tension: f64,
}

impl Droplet {
    /// Create a new droplet with the given parameters.
    pub fn new(center: [f64; 3], radius: f64, velocity: [f64; 3], surface_tension: f64) -> Self {
        Self {
            center,
            radius,
            velocity,
            surface_tension,
        }
    }

    /// Droplet volume (sphere approximation) in lattice units³.
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }

    /// Droplet mass given a bulk density ρ.
    pub fn mass(&self, density: f64) -> f64 {
        density * self.volume()
    }

    /// Kinetic energy of the droplet given density ρ.
    pub fn kinetic_energy(&self, density: f64) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        0.5 * self.mass(density) * v2
    }

    /// Speed (magnitude of velocity vector).
    pub fn speed(&self) -> f64 {
        (self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2)).sqrt()
    }

    /// Weber number We = ρ u² R / σ.
    ///
    /// Quantifies the ratio of inertial forces to surface tension forces.
    pub fn weber_number(&self, density: f64) -> f64 {
        let u2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        if self.surface_tension == 0.0 {
            return f64::INFINITY;
        }
        density * u2 * self.radius / self.surface_tension
    }

    /// Capillary number Ca = μ u / σ.
    pub fn capillary_number(&self, viscosity: f64) -> f64 {
        if self.surface_tension == 0.0 {
            return f64::INFINITY;
        }
        viscosity * self.speed() / self.surface_tension
    }

    /// Centre-to-centre distance to another droplet.
    pub fn distance_to(&self, other: &Droplet) -> f64 {
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        let dz = self.center[2] - other.center[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Returns `true` if the surfaces of the two droplets are closer than `gap`.
    pub fn is_near(&self, other: &Droplet, gap: f64) -> bool {
        self.distance_to(other) < self.radius + other.radius + gap
    }
}

/// Model that determines how two touching droplets merge.
#[derive(Debug, Clone, PartialEq)]
pub enum CoalescenceModel {
    /// Droplets merge instantly upon contact.
    Immediate,
    /// Coalescence is delayed by thin-film drainage; time scales with Ca.
    FilmDrainage {
        /// Drainage time constant τ (time steps).
        drainage_time: f64,
    },
    /// Arbitrary user-defined delay (time steps).
    Delayed {
        /// Fixed delay before merging (time steps).
        delay_steps: u64,
    },
}

impl CoalescenceModel {
    /// Returns the coalescence delay in time steps given the capillary number Ca.
    pub fn delay(&self, capillary_number: f64) -> f64 {
        match self {
            CoalescenceModel::Immediate => 0.0,
            CoalescenceModel::FilmDrainage { drainage_time } => drainage_time * capillary_number,
            CoalescenceModel::Delayed { delay_steps } => *delay_steps as f64,
        }
    }
}

/// Criteria for droplet breakup.
///
/// A droplet breaks when *either* the Weber number or the capillary number
/// exceeds its respective critical threshold.
#[derive(Debug, Clone)]
pub struct BreakupCriteria {
    /// Critical Weber number above which breakup occurs.
    pub critical_weber: f64,
    /// Critical capillary number above which breakup occurs.
    pub critical_capillary: f64,
}

impl BreakupCriteria {
    /// Create breakup criteria with given critical We and Ca.
    pub fn new(critical_weber: f64, critical_capillary: f64) -> Self {
        Self {
            critical_weber,
            critical_capillary,
        }
    }

    /// Returns `true` if the droplet should break given density and viscosity.
    pub fn should_break(&self, droplet: &Droplet, density: f64, viscosity: f64) -> bool {
        droplet.weber_number(density) > self.critical_weber
            || droplet.capillary_number(viscosity) > self.critical_capillary
    }
}

/// Tracks droplets found in a phase-field or density field.
///
/// Droplets are identified by thresholding the order-parameter field and
/// labelling connected regions (simplified flood-fill approach).
#[derive(Debug, Clone)]
pub struct DropletTracker {
    /// Currently tracked droplets.
    pub droplets: Vec<Droplet>,
    /// Phase-field threshold separating droplet interior from background.
    pub phi_threshold: f64,
    /// Default surface tension assigned to newly detected droplets.
    pub default_surface_tension: f64,
    /// Minimum radius (lattice units) — smaller blobs are discarded.
    pub min_radius: f64,
}

impl DropletTracker {
    /// Create a new tracker.
    pub fn new(phi_threshold: f64, default_surface_tension: f64, min_radius: f64) -> Self {
        Self {
            droplets: Vec::new(),
            phi_threshold,
            default_surface_tension,
            min_radius,
        }
    }

    /// Update the droplet list from a 2D phase field (row-major, size nx×ny).
    ///
    /// Uses a simplified labelling that finds the bounding centroid and
    /// equivalent radius for each connected component above `phi_threshold`.
    pub fn update_from_phase_field_2d(
        &mut self,
        phi: &[f64],
        nx: usize,
        ny: usize,
        velocity_field: Option<&[[f64; 3]]>,
    ) {
        let mut visited = vec![false; nx * ny];
        self.droplets.clear();

        for start in 0..(nx * ny) {
            if visited[start] || phi[start] < self.phi_threshold {
                continue;
            }
            // BFS flood fill
            let mut queue = vec![start];
            let mut cells: Vec<usize> = Vec::new();
            visited[start] = true;

            while let Some(idx) = queue.pop() {
                cells.push(idx);
                let x = idx % nx;
                let y = idx / nx;
                let neighbors = [
                    if x > 0 { Some(idx - 1) } else { None },
                    if x + 1 < nx { Some(idx + 1) } else { None },
                    if y > 0 { Some(idx - nx) } else { None },
                    if y + 1 < ny { Some(idx + nx) } else { None },
                ];
                for nb in neighbors.into_iter().flatten() {
                    if !visited[nb] && phi[nb] >= self.phi_threshold {
                        visited[nb] = true;
                        queue.push(nb);
                    }
                }
            }

            // Compute centroid and equivalent radius
            let n = cells.len() as f64;
            let cx = cells.iter().map(|&i| (i % nx) as f64).sum::<f64>() / n;
            let cy = cells.iter().map(|&i| (i / nx) as f64).sum::<f64>() / n;
            let area = n;
            let radius = (area / PI).sqrt();

            if radius < self.min_radius {
                continue;
            }

            // Average velocity from field if provided
            let vel = if let Some(vf) = velocity_field {
                let sum: [f64; 3] = cells.iter().fold([0.0; 3], |mut acc, &i| {
                    if i < vf.len() {
                        acc[0] += vf[i][0];
                        acc[1] += vf[i][1];
                        acc[2] += vf[i][2];
                    }
                    acc
                });
                [sum[0] / n, sum[1] / n, sum[2] / n]
            } else {
                [0.0; 3]
            };

            self.droplets.push(Droplet::new(
                [cx, cy, 0.0],
                radius,
                vel,
                self.default_surface_tension,
            ));
        }
    }

    /// Detect pairs of droplets that are close enough to coalesce.
    ///
    /// Returns indices of droplet pairs whose surfaces are within `gap` of
    /// each other.
    pub fn detect_coalescence_candidates(&self, gap: f64) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..self.droplets.len() {
            for j in (i + 1)..self.droplets.len() {
                if self.droplets[i].is_near(&self.droplets[j], gap) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }

    /// Merge two droplets (volume-conserving) and return the merged droplet.
    ///
    /// The merged centre is the volume-weighted average of the two centres.
    pub fn merge(&self, a: usize, b: usize) -> Droplet {
        let da = &self.droplets[a];
        let db = &self.droplets[b];
        let va = da.volume();
        let vb = db.volume();
        let v_total = va + vb;
        let r_new = (3.0 * v_total / (4.0 * PI)).cbrt();
        let cx = (va * da.center[0] + vb * db.center[0]) / v_total;
        let cy = (va * da.center[1] + vb * db.center[1]) / v_total;
        let cz = (va * da.center[2] + vb * db.center[2]) / v_total;
        // Momentum-conserving velocity (equal density assumed)
        let vx = (va * da.velocity[0] + vb * db.velocity[0]) / v_total;
        let vy = (va * da.velocity[1] + vb * db.velocity[1]) / v_total;
        let vz = (va * da.velocity[2] + vb * db.velocity[2]) / v_total;
        let sigma = (da.surface_tension + db.surface_tension) / 2.0;
        Droplet::new([cx, cy, cz], r_new, [vx, vy, vz], sigma)
    }

    /// Return the count of currently tracked droplets.
    pub fn count(&self) -> usize {
        self.droplets.len()
    }
}

/// Simulates droplet deformation under external flow using the Weber number.
///
/// Uses the Taylor deformation parameter D = (L - B) / (L + B) where L and B
/// are the major and minor semi-axes of the deformed droplet.
#[derive(Debug, Clone)]
pub struct DropletDynamics {
    /// The droplet being simulated.
    pub droplet: Droplet,
    /// Continuous-phase (outer fluid) density ρ_c.
    pub outer_density: f64,
    /// Continuous-phase dynamic viscosity μ_c.
    pub outer_viscosity: f64,
    /// Droplet-phase dynamic viscosity μ_d.
    pub inner_viscosity: f64,
    /// Current Taylor deformation parameter D ∈ \[-1, 1\].
    pub deformation: f64,
    /// Breakup criteria applied during stepping.
    pub breakup_criteria: BreakupCriteria,
}

impl DropletDynamics {
    /// Create a new DropletDynamics instance.
    pub fn new(
        droplet: Droplet,
        outer_density: f64,
        outer_viscosity: f64,
        inner_viscosity: f64,
        breakup_criteria: BreakupCriteria,
    ) -> Self {
        Self {
            droplet,
            outer_density,
            outer_viscosity,
            inner_viscosity,
            deformation: 0.0,
            breakup_criteria,
        }
    }

    /// Viscosity ratio λ = μ_d / μ_c.
    pub fn viscosity_ratio(&self) -> f64 {
        self.inner_viscosity / self.outer_viscosity
    }

    /// Equilibrium Taylor deformation (small-deformation theory):
    /// D_eq = We * (19λ + 16) / (16λ + 16).
    pub fn equilibrium_deformation(&self) -> f64 {
        let we = self.droplet.weber_number(self.outer_density);
        let lambda = self.viscosity_ratio();
        we * (19.0 * lambda + 16.0) / (16.0 * (lambda + 1.0))
    }

    /// Advance the deformation towards its equilibrium value over one time step.
    ///
    /// Uses a relaxation model: dD/dt = (D_eq - D) / τ_deform.
    pub fn step(&mut self, dt: f64, relaxation_time: f64) {
        let d_eq = self.equilibrium_deformation();
        self.deformation += dt / relaxation_time * (d_eq - self.deformation);
    }

    /// Returns `true` if the current We or Ca exceeds the breakup threshold.
    pub fn will_break(&self) -> bool {
        self.breakup_criteria
            .should_break(&self.droplet, self.outer_density, self.outer_viscosity)
    }

    /// Apply an external shear velocity increment to the droplet.
    pub fn apply_shear(&mut self, shear_rate: f64, dt: f64) {
        // Simple shear: u_x += shear_rate * R * dt
        self.droplet.velocity[0] += shear_rate * self.droplet.radius * dt;
    }
}

/// Compute the natural oscillation frequency of a free droplet (Rayleigh mode).
///
/// f = (1/2π) * sqrt(n(n-1)(n+2)σ / (ρ_in + ρ_out/(n+1)) R³)
///
/// For the dominant mode n=2 this simplifies to:
/// f = (1/2π) * sqrt(8σ / (ρ R³))
///
/// # Arguments
/// * `radius`           – droplet radius \[m or lattice units\]
/// * `surface_tension`  – surface tension σ \[N/m or LU equivalent\]
/// * `density`          – droplet bulk density ρ \[kg/m³ or LU equivalent\]
pub fn droplet_oscillation_freq(radius: f64, surface_tension: f64, density: f64) -> f64 {
    if radius <= 0.0 || density <= 0.0 || surface_tension < 0.0 {
        return 0.0;
    }
    let omega2 = 8.0 * surface_tension / (density * radius.powi(3));
    omega2.sqrt() / (2.0 * PI)
}

/// Compute the Laplace pressure jump across a spherical droplet interface.
///
/// ΔP = 2σ / R
pub fn laplace_pressure(radius: f64, surface_tension: f64) -> f64 {
    if radius <= 0.0 {
        return 0.0;
    }
    2.0 * surface_tension / radius
}

/// Compute the Bond number Bo = ρ g R² / σ (ratio of gravity to surface tension).
pub fn bond_number(density: f64, gravity: f64, radius: f64, surface_tension: f64) -> f64 {
    if surface_tension == 0.0 {
        return f64::INFINITY;
    }
    density * gravity * radius.powi(2) / surface_tension
}

/// Compute the Ohnesorge number Oh = μ / sqrt(ρ σ R).
///
/// Relates viscous to inertial and surface-tension forces.
pub fn ohnesorge_number(viscosity: f64, density: f64, surface_tension: f64, radius: f64) -> f64 {
    let denom = (density * surface_tension * radius).sqrt();
    if denom == 0.0 {
        return f64::INFINITY;
    }
    viscosity / denom
}

// ---------------------------------------------------------------------------
// DropletLBM — diffuse-interface droplet using a phase-field order parameter
// ---------------------------------------------------------------------------

/// Diffuse-interface droplet simulation using a phase-field (order-parameter) LBM.
///
/// The phase field φ ∈ \[0,1\] distinguishes liquid (φ≈1) from gas (φ≈0).
/// The distribution function `f_dist` holds the D2Q9 populations for the
/// flow solver (9 velocities per node) and `phi` stores the scalar field.
#[derive(Debug, Clone)]
pub struct DropletLBM {
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Phase-field (order parameter) φ, length nx*ny.
    pub phi: Vec<f64>,
    /// Distribution function populations, length nx*ny*9.
    pub f_dist: Vec<f64>,
    /// Surface tension coefficient σ.
    pub sigma: f64,
    /// Liquid density ρ_l.
    pub rho_l: f64,
    /// Gas density ρ_g.
    pub rho_g: f64,
    /// Relaxation time τ for the BGK collision.
    pub tau: f64,
}

impl DropletLBM {
    /// Create a new `DropletLBM` with equilibrium initialisation.
    ///
    /// The phase field is initialised to `rho_g` everywhere; individual
    /// droplets can be added with [`Self::init_circular_droplet`].
    pub fn new(nx: usize, ny: usize, sigma: f64, rho_l: f64, rho_g: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            phi: vec![rho_g; n],
            f_dist: vec![0.0; n * 9],
            sigma,
            rho_l,
            rho_g,
            tau: 1.0,
        }
    }

    /// Row-major index for node (i, j).
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Index into the distribution function array for velocity direction `q`.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }

    /// Initialise the phase field for a circular droplet centred at (cx, cy)
    /// with radius `r`, using a hyperbolic-tangent interface profile.
    pub fn init_circular_droplet(&mut self, cx: f64, cy: f64, r: f64) {
        let w = 2.0_f64; // interface width in lattice units
        for j in 0..self.ny {
            for i in 0..self.nx {
                let dx = i as f64 - cx;
                let dy = j as f64 - cy;
                let dist = (dx * dx + dy * dy).sqrt();
                let phi_val = 0.5 * (1.0 + ((r - dist) / w).tanh());
                let idx = self.index(i, j);
                let phi_mapped = self.rho_g + (self.rho_l - self.rho_g) * phi_val;
                self.phi[idx] = phi_mapped;
            }
        }
    }
}

/// Young-Laplace pressure jump across a 2D circular droplet interface.
///
/// In 2D the curvature is 1/R, so ΔP = σ/R.
///
/// # Arguments
/// * `sigma`  – Surface tension (N/m or lattice units).
/// * `radius` – Droplet radius (m or lattice units).
pub fn surface_tension_pressure_jump(sigma: f64, radius: f64) -> f64 {
    if radius <= 0.0 {
        return 0.0;
    }
    sigma / radius
}

/// Young's equation: contact angle from interfacial energies.
///
/// cos θ = (γ_sg − γ_sl) / γ_lg
///
/// Returns the contact angle in radians, clamped to \[0, π\].
///
/// # Arguments
/// * `gamma_sg` – Solid–gas interfacial energy.
/// * `gamma_sl` – Solid–liquid interfacial energy.
/// * `gamma_lg` – Liquid–gas interfacial energy.
pub fn contact_angle_from_energy(gamma_sg: f64, gamma_sl: f64, gamma_lg: f64) -> f64 {
    if gamma_lg.abs() < f64::EPSILON {
        return 0.0;
    }
    let cos_theta = (gamma_sg - gamma_sl) / gamma_lg;
    cos_theta.clamp(-1.0, 1.0).acos()
}

/// Mass-weighted velocity of the droplet phase.
///
/// `u = Σ φ_i u_i / Σ φ_i`
///
/// Returns `[0.0, 0.0]` when the total mass is negligible.
///
/// # Arguments
/// * `phi` – Phase field (nx*ny).
/// * `ux`  – x-velocity field (nx*ny).
/// * `uy`  – y-velocity field (nx*ny).
/// * `nx`  – Grid width.
/// * `ny`  – Grid height.
pub fn droplet_velocity(phi: &[f64], ux: &[f64], uy: &[f64], nx: usize, ny: usize) -> [f64; 2] {
    let n = nx * ny;
    let mut sum_phi = 0.0_f64;
    let mut sum_px = 0.0_f64;
    let mut sum_py = 0.0_f64;
    for i in 0..n {
        sum_phi += phi[i];
        sum_px += phi[i] * ux[i];
        sum_py += phi[i] * uy[i];
    }
    if sum_phi < f64::EPSILON {
        return [0.0, 0.0];
    }
    [sum_px / sum_phi, sum_py / sum_phi]
}

/// Estimate the droplet radius from the total liquid volume.
///
/// In 2D: volume ≈ π R², so R = sqrt(A_liquid / π).
///
/// Cells with φ > (ρ_l + ρ_g)/2 are counted as liquid.
///
/// # Arguments
/// * `phi`   – Phase field.
/// * `nx`    – Grid width.
/// * `ny`    – Grid height.
/// * `rho_l` – Liquid density (threshold midpoint).
/// * `rho_g` – Gas density.
pub fn droplet_radius_estimate(phi: &[f64], nx: usize, ny: usize, rho_l: f64, rho_g: f64) -> f64 {
    let threshold = (rho_l + rho_g) * 0.5;
    let area = phi.iter().filter(|&&p| p > threshold).count() as f64;
    let _ = (nx, ny); // grid dims provided for interface consistency
    (area / PI).sqrt()
}

/// Compute the centre of mass of the droplet phase.
///
/// Returns the (x, y) centroid weighted by φ.
///
/// # Arguments
/// * `phi` – Phase field (nx*ny).
/// * `nx`  – Grid width.
/// * `ny`  – Grid height.
pub fn droplet_center_of_mass(phi: &[f64], nx: usize, ny: usize) -> [f64; 2] {
    let mut sum_phi = 0.0_f64;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    for j in 0..ny {
        for i in 0..nx {
            let p = phi[j * nx + i];
            sum_phi += p;
            cx += p * i as f64;
            cy += p * j as f64;
        }
    }
    if sum_phi < f64::EPSILON {
        return [0.0, 0.0];
    }
    [cx / sum_phi, cy / sum_phi]
}

/// Weber number: ratio of inertial to surface tension forces.
///
/// We = ρ v² R / σ
///
/// # Arguments
/// * `rho`   – Fluid density.
/// * `v`     – Characteristic velocity.
/// * `r`     – Droplet radius.
/// * `sigma` – Surface tension.
pub fn weber_number(rho: f64, v: f64, r: f64, sigma: f64) -> f64 {
    if sigma == 0.0 {
        return f64::INFINITY;
    }
    rho * v * v * r / sigma
}

/// Capillary number: ratio of viscous to surface tension forces.
///
/// Ca = μ v / σ
///
/// # Arguments
/// * `mu`    – Dynamic viscosity.
/// * `v`     – Characteristic velocity.
/// * `sigma` – Surface tension.
pub fn capillary_number(mu: f64, v: f64, sigma: f64) -> f64 {
    if sigma == 0.0 {
        return f64::INFINITY;
    }
    mu * v / sigma
}

/// Ohnesorge number: ratio of viscous forces to inertial and surface tension.
///
/// Oh = μ / sqrt(ρ σ R)
///
/// # Arguments
/// * `mu`    – Dynamic viscosity.
/// * `rho`   – Fluid density.
/// * `r`     – Droplet radius.
/// * `sigma` – Surface tension.
pub fn ohnesorge_number_lbm(mu: f64, rho: f64, r: f64, sigma: f64) -> f64 {
    let denom = (rho * sigma * r).sqrt();
    if denom < f64::EPSILON {
        return f64::INFINITY;
    }
    mu / denom
}

/// Bond number: ratio of gravitational to surface tension forces.
///
/// Bo = ρ g R² / σ
///
/// # Arguments
/// * `rho`   – Fluid density.
/// * `g`     – Gravitational acceleration.
/// * `r`     – Droplet radius.
/// * `sigma` – Surface tension.
pub fn bond_number_lbm(rho: f64, g: f64, r: f64, sigma: f64) -> f64 {
    if sigma == 0.0 {
        return f64::INFINITY;
    }
    rho * g * r * r / sigma
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_droplet(r: f64, sigma: f64) -> Droplet {
        Droplet::new([0.0, 0.0, 0.0], r, [1.0, 0.0, 0.0], sigma)
    }

    // --- Droplet basics ---

    #[test]
    fn test_droplet_volume() {
        let d = make_droplet(1.0, 0.1);
        let expected = 4.0 / 3.0 * PI;
        assert!((d.volume() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_mass() {
        let d = make_droplet(1.0, 0.1);
        let rho = 2.0;
        assert!((d.mass(rho) - rho * d.volume()).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_speed() {
        let d = Droplet::new([0.0; 3], 1.0, [3.0, 4.0, 0.0], 0.1);
        assert!((d.speed() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_kinetic_energy() {
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 0.1);
        let rho = 1.0;
        let ke = d.kinetic_energy(rho);
        let expected = 0.5 * d.mass(rho);
        assert!((ke - expected).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_weber_zero_velocity() {
        let d = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        assert_eq!(d.weber_number(1.0), 0.0);
    }

    #[test]
    fn test_droplet_weber_nonzero() {
        let d = Droplet::new([0.0; 3], 2.0, [1.0, 0.0, 0.0], 1.0);
        // We = 1.0 * 1.0 * 2.0 / 1.0 = 2.0
        assert!((d.weber_number(1.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_weber_zero_sigma() {
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 0.0);
        assert!(d.weber_number(1.0).is_infinite());
    }

    #[test]
    fn test_droplet_capillary() {
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 1.0);
        // Ca = mu * u / sigma = 0.5 * 1 / 1 = 0.5
        assert!((d.capillary_number(0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_distance() {
        let a = Droplet::new([0.0, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([3.0, 4.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_is_near_true() {
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([2.0, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        // distance = 2.0, r+r = 2.0, gap = 0.1 → is_near
        assert!(a.is_near(&b, 0.1));
    }

    #[test]
    fn test_droplet_is_near_false() {
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([10.0, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!(!a.is_near(&b, 0.1));
    }

    // --- CoalescenceModel ---

    #[test]
    fn test_coalescence_immediate() {
        let m = CoalescenceModel::Immediate;
        assert_eq!(m.delay(0.5), 0.0);
    }

    #[test]
    fn test_coalescence_film_drainage() {
        let m = CoalescenceModel::FilmDrainage {
            drainage_time: 10.0,
        };
        assert!((m.delay(0.3) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_coalescence_delayed() {
        let m = CoalescenceModel::Delayed { delay_steps: 5 };
        assert!((m.delay(0.0) - 5.0).abs() < 1e-12);
    }

    // --- BreakupCriteria ---

    #[test]
    fn test_breakup_no_break_low_we() {
        let crit = BreakupCriteria::new(10.0, 10.0);
        let d = Droplet::new([0.0; 3], 1.0, [0.5, 0.0, 0.0], 1.0);
        assert!(!crit.should_break(&d, 1.0, 1.0));
    }

    #[test]
    fn test_breakup_breaks_high_we() {
        let crit = BreakupCriteria::new(1.0, 100.0);
        // We = 1.0 * 100.0 * 1.0 / 0.1 = 1000 > 1.0
        let d = Droplet::new([0.0; 3], 1.0, [10.0, 0.0, 0.0], 0.1);
        assert!(crit.should_break(&d, 1.0, 0.1));
    }

    #[test]
    fn test_breakup_breaks_high_ca() {
        let crit = BreakupCriteria::new(1000.0, 0.5);
        // Ca = mu * u / sigma = 1.0 * 1.0 / 1.0 = 1.0 > 0.5
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 1.0);
        assert!(crit.should_break(&d, 1.0, 1.0));
    }

    // --- DropletTracker ---

    #[test]
    fn test_tracker_empty_field() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        let phi = vec![0.0; 100];
        tracker.update_from_phase_field_2d(&phi, 10, 10, None);
        assert_eq!(tracker.count(), 0);
    }

    #[test]
    fn test_tracker_single_droplet() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        // 5×5 patch of 1.0 in a 20×20 field
        let mut phi = vec![0.0f64; 20 * 20];
        for y in 7..12usize {
            for x in 7..12usize {
                phi[y * 20 + x] = 1.0;
            }
        }
        tracker.update_from_phase_field_2d(&phi, 20, 20, None);
        assert_eq!(tracker.count(), 1);
        let d = &tracker.droplets[0];
        assert!((d.center[0] - 9.0).abs() < 1.0);
        assert!((d.center[1] - 9.0).abs() < 1.0);
    }

    #[test]
    fn test_tracker_two_droplets() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        let mut phi = vec![0.0f64; 40 * 20];
        // Droplet 1: centre ~(5,5)
        for y in 3..8usize {
            for x in 3..8usize {
                phi[y * 40 + x] = 1.0;
            }
        }
        // Droplet 2: centre ~(30,10)
        for y in 7..13usize {
            for x in 27..33usize {
                phi[y * 40 + x] = 1.0;
            }
        }
        tracker.update_from_phase_field_2d(&phi, 40, 20, None);
        assert_eq!(tracker.count(), 2);
    }

    #[test]
    fn test_tracker_small_blob_discarded() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 5.0);
        let mut phi = vec![0.0f64; 20 * 20];
        // 2×2 blob — radius ~ 0.8, below min_radius 5.0
        phi[5 * 20 + 5] = 1.0;
        phi[5 * 20 + 6] = 1.0;
        phi[6 * 20 + 5] = 1.0;
        phi[6 * 20 + 6] = 1.0;
        tracker.update_from_phase_field_2d(&phi, 20, 20, None);
        assert_eq!(tracker.count(), 0);
    }

    #[test]
    fn test_tracker_detect_coalescence() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        // Manually place two close droplets
        tracker
            .droplets
            .push(Droplet::new([0.0, 0.0, 0.0], 2.0, [0.0; 3], 0.07));
        tracker
            .droplets
            .push(Droplet::new([3.5, 0.0, 0.0], 2.0, [0.0; 3], 0.07));
        let pairs = tracker.detect_coalescence_candidates(0.0);
        // distance = 3.5 < 2+2=4 → candidate
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 1));
    }

    #[test]
    fn test_tracker_no_coalescence_far() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        tracker
            .droplets
            .push(Droplet::new([0.0, 0.0, 0.0], 1.0, [0.0; 3], 0.07));
        tracker
            .droplets
            .push(Droplet::new([20.0, 0.0, 0.0], 1.0, [0.0; 3], 0.07));
        let pairs = tracker.detect_coalescence_candidates(0.0);
        assert_eq!(pairs.len(), 0);
    }

    #[test]
    fn test_tracker_merge_volume_conservation() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        tracker
            .droplets
            .push(Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 0.07));
        tracker
            .droplets
            .push(Droplet::new([3.0, 0.0, 0.0], 2.0, [0.0; 3], 0.07));
        let merged = tracker.merge(0, 1);
        let v0 = tracker.droplets[0].volume();
        let v1 = tracker.droplets[1].volume();
        assert!((merged.volume() - (v0 + v1)).abs() < 1e-8);
    }

    #[test]
    fn test_tracker_merge_momentum_conservation() {
        let mut tracker = DropletTracker::new(0.5, 0.07, 1.0);
        tracker
            .droplets
            .push(Droplet::new([0.0; 3], 1.0, [2.0, 0.0, 0.0], 0.07));
        tracker
            .droplets
            .push(Droplet::new([3.0, 0.0, 0.0], 1.0, [0.0; 3], 0.07));
        let merged = tracker.merge(0, 1);
        let v0 = tracker.droplets[0].volume();
        let v1 = tracker.droplets[1].volume();
        let expected_vx = (v0 * 2.0 + v1 * 0.0) / (v0 + v1);
        assert!((merged.velocity[0] - expected_vx).abs() < 1e-10);
    }

    // --- DropletDynamics ---

    #[test]
    fn test_dynamics_viscosity_ratio() {
        let d = make_droplet(1.0, 0.07);
        let crit = BreakupCriteria::new(10.0, 10.0);
        let dd = DropletDynamics::new(d, 1.0, 0.1, 0.5, crit);
        assert!((dd.viscosity_ratio() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_dynamics_equilibrium_deformation_zero_velocity() {
        let d = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.07);
        let crit = BreakupCriteria::new(10.0, 10.0);
        let dd = DropletDynamics::new(d, 1.0, 0.1, 0.1, crit);
        assert!((dd.equilibrium_deformation() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dynamics_step_approaches_equilibrium() {
        let d = Droplet::new([0.0; 3], 1.0, [0.1, 0.0, 0.0], 0.07);
        let crit = BreakupCriteria::new(10.0, 10.0);
        let mut dd = DropletDynamics::new(d, 1.0, 0.1, 0.1, crit);
        let d_eq = dd.equilibrium_deformation();
        for _ in 0..1000 {
            dd.step(0.01, 1.0);
        }
        assert!((dd.deformation - d_eq).abs() < 0.01);
    }

    #[test]
    fn test_dynamics_will_break_false() {
        let d = Droplet::new([0.0; 3], 1.0, [0.1, 0.0, 0.0], 1.0);
        let crit = BreakupCriteria::new(100.0, 100.0);
        let dd = DropletDynamics::new(d, 1.0, 0.01, 0.01, crit);
        assert!(!dd.will_break());
    }

    #[test]
    fn test_dynamics_will_break_true() {
        let d = Droplet::new([0.0; 3], 1.0, [100.0, 0.0, 0.0], 0.01);
        let crit = BreakupCriteria::new(1.0, 100.0);
        let dd = DropletDynamics::new(d, 1.0, 0.1, 0.1, crit);
        assert!(dd.will_break());
    }

    #[test]
    fn test_dynamics_shear_increases_velocity() {
        let d = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let crit = BreakupCriteria::new(10.0, 10.0);
        let mut dd = DropletDynamics::new(d, 1.0, 0.1, 0.1, crit);
        dd.apply_shear(1.0, 1.0);
        assert!(dd.droplet.velocity[0] > 0.0);
    }

    // --- Standalone functions ---

    #[test]
    fn test_oscillation_freq_positive() {
        let f = droplet_oscillation_freq(1e-3, 0.072, 1000.0);
        assert!(f > 0.0);
    }

    #[test]
    fn test_oscillation_freq_zero_radius() {
        assert_eq!(droplet_oscillation_freq(0.0, 0.072, 1000.0), 0.0);
    }

    #[test]
    fn test_oscillation_freq_zero_density() {
        assert_eq!(droplet_oscillation_freq(1e-3, 0.072, 0.0), 0.0);
    }

    #[test]
    fn test_oscillation_freq_negative_radius() {
        assert_eq!(droplet_oscillation_freq(-1.0, 0.072, 1000.0), 0.0);
    }

    #[test]
    fn test_oscillation_freq_scales_with_radius() {
        let f1 = droplet_oscillation_freq(1.0, 1.0, 1.0);
        let f2 = droplet_oscillation_freq(2.0, 1.0, 1.0);
        // ω ∝ R^{-3/2}, so f2 < f1
        assert!(f2 < f1);
    }

    #[test]
    fn test_oscillation_freq_scales_with_sigma() {
        let f1 = droplet_oscillation_freq(1.0, 1.0, 1.0);
        let f2 = droplet_oscillation_freq(1.0, 4.0, 1.0);
        assert!((f2 / f1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_laplace_pressure_sphere() {
        let dp = laplace_pressure(1.0, 0.5);
        assert!((dp - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_laplace_pressure_zero_radius() {
        assert_eq!(laplace_pressure(0.0, 0.5), 0.0);
    }

    #[test]
    fn test_bond_number_zero_sigma() {
        assert!(bond_number(1.0, 9.81, 1.0, 0.0).is_infinite());
    }

    #[test]
    fn test_bond_number_value() {
        let bo = bond_number(1000.0, 9.81, 1e-3, 0.072);
        // Bo = 1000 * 9.81 * 1e-6 / 0.072 ≈ 0.136
        assert!(bo > 0.0 && bo < 1.0);
    }

    #[test]
    fn test_ohnesorge_number() {
        let oh = ohnesorge_number(0.001, 1000.0, 0.072, 1e-3);
        assert!(oh > 0.0);
    }

    #[test]
    fn test_ohnesorge_zero_denom() {
        let oh = ohnesorge_number(0.001, 0.0, 0.072, 1e-3);
        assert!(oh.is_infinite());
    }

    // --- DropletLBM ---

    #[test]
    fn test_droplet_lbm_new_sizes() {
        let d = DropletLBM::new(10, 8, 0.01, 1.0, 0.1);
        assert_eq!(d.phi.len(), 80);
        assert_eq!(d.f_dist.len(), 720);
    }

    #[test]
    fn test_droplet_lbm_index() {
        let d = DropletLBM::new(10, 8, 0.01, 1.0, 0.1);
        assert_eq!(d.index(3, 2), 23);
    }

    #[test]
    fn test_droplet_lbm_fi() {
        let d = DropletLBM::new(10, 8, 0.01, 1.0, 0.1);
        assert_eq!(d.fi(1, 0, 3), 12); // (0*10+1)*9 + 3 = 12
    }

    #[test]
    fn test_droplet_lbm_init_circular_center_high() {
        let mut d = DropletLBM::new(20, 20, 0.01, 1.0, 0.01);
        d.init_circular_droplet(10.0, 10.0, 5.0);
        let idx = d.index(10, 10);
        // centre should be close to rho_l
        assert!(d.phi[idx] > 0.5);
    }

    #[test]
    fn test_droplet_lbm_init_circular_outside_low() {
        let mut d = DropletLBM::new(20, 20, 0.01, 1.0, 0.01);
        d.init_circular_droplet(10.0, 10.0, 3.0);
        let idx = d.index(0, 0);
        // corner far from droplet should be near rho_g
        assert!(d.phi[idx] < 0.5);
    }

    #[test]
    fn test_droplet_lbm_phi_range() {
        let mut d = DropletLBM::new(16, 16, 0.01, 1.0, 0.0);
        d.init_circular_droplet(8.0, 8.0, 4.0);
        for &p in &d.phi {
            assert!((0.0..=(1.0 + 1e-12)).contains(&p), "p={p}");
        }
    }

    // --- surface_tension_pressure_jump ---

    #[test]
    fn test_surface_tension_pressure_jump_basic() {
        let dp = surface_tension_pressure_jump(0.5, 1.0);
        assert!((dp - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_surface_tension_pressure_jump_zero_radius() {
        assert_eq!(surface_tension_pressure_jump(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_surface_tension_pressure_jump_increases_with_sigma() {
        let dp1 = surface_tension_pressure_jump(1.0, 1.0);
        let dp2 = surface_tension_pressure_jump(2.0, 1.0);
        assert!(dp2 > dp1);
    }

    #[test]
    fn test_surface_tension_pressure_jump_decreases_with_radius() {
        let dp1 = surface_tension_pressure_jump(1.0, 1.0);
        let dp2 = surface_tension_pressure_jump(1.0, 2.0);
        assert!(dp2 < dp1);
    }

    // --- contact_angle_from_energy ---

    #[test]
    fn test_contact_angle_hydrophilic() {
        // gamma_sg > gamma_sl → cos θ > 0 → θ < pi/2
        let theta = contact_angle_from_energy(1.0, 0.5, 1.0);
        assert!(theta < PI / 2.0);
    }

    #[test]
    fn test_contact_angle_hydrophobic() {
        // gamma_sg < gamma_sl → cos θ < 0 → θ > pi/2
        let theta = contact_angle_from_energy(0.5, 1.0, 1.0);
        assert!(theta > PI / 2.0);
    }

    #[test]
    fn test_contact_angle_zero_gamma_lg() {
        let theta = contact_angle_from_energy(1.0, 0.5, 0.0);
        assert_eq!(theta, 0.0);
    }

    #[test]
    fn test_contact_angle_range() {
        let theta = contact_angle_from_energy(0.7, 0.3, 1.0);
        assert!((0.0..=PI).contains(&theta));
    }

    // --- weber_number (standalone) ---

    #[test]
    fn test_weber_number_standalone() {
        let we = weber_number(1.0, 1.0, 1.0, 1.0);
        assert!((we - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_weber_number_zero_sigma() {
        assert!(weber_number(1.0, 1.0, 1.0, 0.0).is_infinite());
    }

    #[test]
    fn test_weber_number_scales_with_velocity() {
        let we1 = weber_number(1.0, 1.0, 1.0, 1.0);
        let we2 = weber_number(1.0, 2.0, 1.0, 1.0);
        assert!((we2 - 4.0 * we1).abs() < 1e-12);
    }

    // --- capillary_number (standalone) ---

    #[test]
    fn test_capillary_number_standalone() {
        let ca = capillary_number(0.5, 2.0, 1.0);
        assert!((ca - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_capillary_number_zero_sigma() {
        assert!(capillary_number(1.0, 1.0, 0.0).is_infinite());
    }

    // --- ohnesorge_number_lbm ---

    #[test]
    fn test_ohnesorge_number_lbm_positive() {
        let oh = ohnesorge_number_lbm(0.001, 1000.0, 1e-3, 0.072);
        assert!(oh > 0.0);
    }

    #[test]
    fn test_ohnesorge_number_lbm_zero_denom() {
        let oh = ohnesorge_number_lbm(0.001, 0.0, 1e-3, 0.072);
        assert!(oh.is_infinite());
    }

    // --- bond_number_lbm ---

    #[test]
    fn test_bond_number_lbm_positive() {
        let bo = bond_number_lbm(1000.0, 9.81, 1e-3, 0.072);
        assert!(bo > 0.0);
    }

    #[test]
    fn test_bond_number_lbm_zero_sigma() {
        assert!(bond_number_lbm(1.0, 9.81, 1.0, 0.0).is_infinite());
    }

    // --- droplet_velocity ---

    #[test]
    fn test_droplet_velocity_uniform() {
        let phi = vec![1.0; 4];
        let ux = vec![2.0; 4];
        let uy = vec![3.0; 4];
        let v = droplet_velocity(&phi, &ux, &uy, 2, 2);
        assert!((v[0] - 2.0).abs() < 1e-12);
        assert!((v[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_droplet_velocity_zero_phi() {
        let phi = vec![0.0; 4];
        let ux = vec![1.0; 4];
        let uy = vec![1.0; 4];
        let v = droplet_velocity(&phi, &ux, &uy, 2, 2);
        assert_eq!(v, [0.0, 0.0]);
    }

    // --- droplet_radius_estimate ---

    #[test]
    fn test_droplet_radius_estimate_full() {
        // All cells are liquid → area = 100 → R = sqrt(100/pi)
        let phi = vec![1.0; 100];
        let r = droplet_radius_estimate(&phi, 10, 10, 1.0, 0.0);
        let expected = (100.0_f64 / PI).sqrt();
        assert!((r - expected).abs() < 1e-10);
    }

    #[test]
    fn test_droplet_radius_estimate_empty() {
        let phi = vec![0.0; 100];
        let r = droplet_radius_estimate(&phi, 10, 10, 1.0, 0.0);
        assert_eq!(r, 0.0);
    }

    // --- droplet_center_of_mass ---

    #[test]
    fn test_droplet_center_of_mass_uniform() {
        // Uniform phi → COM = geometric center
        let phi = vec![1.0; 9];
        let com = droplet_center_of_mass(&phi, 3, 3);
        assert!((com[0] - 1.0).abs() < 1e-10);
        assert!((com[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_droplet_center_of_mass_zero_phi() {
        let phi = vec![0.0; 4];
        let com = droplet_center_of_mass(&phi, 2, 2);
        assert_eq!(com, [0.0, 0.0]);
    }
}
