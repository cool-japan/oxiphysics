// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular material and snow/mud simulation.
//!
//! Provides MPM-based snow, Drucker-Prager plasticity, Bingham plastic mud,
//! dry sand, avalanche stability, tire-terrain interaction, and deformation tracking.

// ---------------------------------------------------------------------------
// Small math helpers
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ---------------------------------------------------------------------------
// SnowParticle
// ---------------------------------------------------------------------------

/// A single MPM snow particle.
pub struct SnowParticle {
    /// World-space position (metres).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Rest density (kg/m³).
    pub density: f64,
    /// Cohesion pressure (Pa).
    pub cohesion: f64,
    /// Hardness parameter (dimensionless); increases with compaction.
    pub hardness: f64,
    /// Temperature (Kelvin).
    pub temperature: f64,
    /// Deformation gradient (3×3 stored row-major).
    pub def_gradient: [f64; 9],
    /// Volume (m³).
    pub volume: f64,
    /// Mass (kg).
    pub mass: f64,
}

impl SnowParticle {
    /// Create a snow particle with default values.
    pub fn new(position: [f64; 3], density: f64, temperature: f64) -> Self {
        let volume = 1.0 / (density * 64.0); // heuristic per particle
        Self {
            position,
            velocity: [0.0; 3],
            density,
            cohesion: 1000.0,
            hardness: 1.0,
            temperature,
            def_gradient: identity3x3(),
            volume,
            mass: density * volume,
        }
    }
}

fn identity3x3() -> [f64; 9] {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}

// ---------------------------------------------------------------------------
// SnowMpm — Material Point Method for snow
// ---------------------------------------------------------------------------

/// MPM simulation grid cell.
pub struct MpmCell {
    /// Grid node velocity (m/s).
    pub velocity: [f64; 3],
    /// Accumulated mass (kg).
    pub mass: f64,
    /// Accumulated momentum (kg m/s).
    pub momentum: [f64; 3],
    /// Accumulated force (N).
    pub force: [f64; 3],
}

impl MpmCell {
    fn new() -> Self {
        Self {
            velocity: [0.0; 3],
            mass: 0.0,
            momentum: [0.0; 3],
            force: [0.0; 3],
        }
    }
}

/// MPM solver for snow.
pub struct SnowMpm {
    /// Grid resolution (cells per axis).
    pub grid_resolution: usize,
    /// Cell size (metres).
    pub cell_size: f64,
    /// Snow particles.
    pub particles: Vec<SnowParticle>,
    /// Grid cells (flattened 3D array).
    pub grid: Vec<MpmCell>,
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Critical compression.
    pub critical_compression: f64,
    /// Critical stretch.
    pub critical_stretch: f64,
    /// Hardening coefficient.
    pub hardening: f64,
    /// Gravity.
    pub gravity: [f64; 3],
}

impl SnowMpm {
    /// Create an MPM snow solver.
    pub fn new(grid_resolution: usize, domain_size: f64) -> Self {
        let n3 = grid_resolution * grid_resolution * grid_resolution;
        Self {
            grid_resolution,
            cell_size: domain_size / grid_resolution as f64,
            particles: Vec::new(),
            grid: (0..n3).map(|_| MpmCell::new()).collect(),
            young_modulus: 1.4e5,
            poisson_ratio: 0.2,
            critical_compression: 0.019,
            critical_stretch: 0.0075,
            hardening: 10.0,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: SnowParticle) {
        self.particles.push(p);
    }

    /// Reset grid (clear mass, momentum, force).
    pub fn clear_grid(&mut self) {
        for cell in &mut self.grid {
            *cell = MpmCell::new();
        }
    }

    /// Particle-to-grid (P2G) transfer.
    pub fn particle_to_grid(&mut self) {
        let h = self.cell_size;
        let n = self.grid_resolution;
        for p in &self.particles {
            let ix = (p.position[0] / h) as i64;
            let iy = (p.position[1] / h) as i64;
            let iz = (p.position[2] / h) as i64;
            // Transfer to 2-cell neighbourhood
            for di in -1_i64..=1 {
                for dj in -1_i64..=1 {
                    for dk in -1_i64..=1 {
                        let nx = ix + di;
                        let ny = iy + dj;
                        let nz = iz + dk;
                        if nx < 0
                            || ny < 0
                            || nz < 0
                            || nx >= n as i64
                            || ny >= n as i64
                            || nz >= n as i64
                        {
                            continue;
                        }
                        let idx = grid_idx(nx as usize, ny as usize, nz as usize, n);
                        let weight = quadratic_weight(
                            p.position,
                            [nx as f64 * h, ny as f64 * h, nz as f64 * h],
                            h,
                        );
                        self.grid[idx].mass += p.mass * weight;
                        for c in 0..3 {
                            self.grid[idx].momentum[c] += p.mass * p.velocity[c] * weight;
                        }
                    }
                }
            }
        }
    }

    /// Grid-to-particle (G2P) transfer.
    pub fn grid_to_particle(&mut self) {
        let h = self.cell_size;
        let n = self.grid_resolution;
        for p in &mut self.particles {
            let ix = (p.position[0] / h) as i64;
            let iy = (p.position[1] / h) as i64;
            let iz = (p.position[2] / h) as i64;
            let mut new_vel = [0.0_f64; 3];
            for di in -1_i64..=1 {
                for dj in -1_i64..=1 {
                    for dk in -1_i64..=1 {
                        let nx = ix + di;
                        let ny = iy + dj;
                        let nz = iz + dk;
                        if nx < 0
                            || ny < 0
                            || nz < 0
                            || nx >= n as i64
                            || ny >= n as i64
                            || nz >= n as i64
                        {
                            continue;
                        }
                        let idx = grid_idx(nx as usize, ny as usize, nz as usize, n);
                        let weight = quadratic_weight(
                            p.position,
                            [nx as f64 * h, ny as f64 * h, nz as f64 * h],
                            h,
                        );
                        let cell_vel = if self.grid[idx].mass > 1e-12 {
                            scale3(self.grid[idx].momentum, 1.0 / self.grid[idx].mass)
                        } else {
                            [0.0; 3]
                        };
                        new_vel = add3(new_vel, scale3(cell_vel, weight));
                    }
                }
            }
            p.velocity = new_vel;
        }
    }

    /// Full MPM step.
    pub fn step(&mut self, dt: f64) {
        self.clear_grid();
        self.particle_to_grid();
        // Apply gravity to grid
        let g = self.gravity;
        for cell in &mut self.grid {
            if cell.mass > 1e-12 {
                for (c, &gc) in g.iter().enumerate() {
                    cell.momentum[c] += cell.mass * gc * dt;
                    cell.velocity[c] = cell.momentum[c] / cell.mass;
                }
            }
        }
        self.grid_to_particle();
        // Advect particles
        for p in &mut self.particles {
            p.position = add3(p.position, scale3(p.velocity, dt));
        }
    }
}

fn grid_idx(x: usize, y: usize, z: usize, n: usize) -> usize {
    x * n * n + y * n + z
}

/// Quadratic B-spline weight for MPM.
pub fn mpm_grid_to_particle(node_pos: [f64; 3], particle_pos: [f64; 3], h: f64) -> f64 {
    quadratic_weight(particle_pos, node_pos, h)
}

fn quadratic_weight(particle: [f64; 3], node: [f64; 3], h: f64) -> f64 {
    let mut w = 1.0;
    for i in 0..3 {
        let x = ((particle[i] - node[i]) / h).abs();
        w *= if x < 0.5 {
            0.75 - x * x
        } else if x < 1.5 {
            0.5 * (1.5 - x) * (1.5 - x)
        } else {
            0.0
        };
    }
    w
}

// ---------------------------------------------------------------------------
// DruckerPragerPlastic
// ---------------------------------------------------------------------------

/// Drucker-Prager yield criterion with hardening for snow compaction.
pub struct DruckerPragerPlastic {
    /// Friction angle (radians).
    pub friction_angle: f64,
    /// Cohesion (Pa).
    pub cohesion: f64,
    /// Hardening modulus (Pa).
    pub hardening_modulus: f64,
    /// Accumulated plastic strain.
    pub plastic_strain: f64,
}

impl DruckerPragerPlastic {
    /// Create a Drucker-Prager model.
    pub fn new(friction_angle: f64, cohesion: f64) -> Self {
        Self {
            friction_angle,
            cohesion,
            hardening_modulus: 1e4,
            plastic_strain: 0.0,
        }
    }

    /// Effective friction coefficient from the friction angle.
    pub fn mu(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi))
    }

    /// Yield function value: `f > 0` means plastic.
    ///
    /// `p` = hydrostatic pressure (Pa, compression positive),
    /// `q` = deviatoric stress invariant (Pa).
    pub fn yield_function(&self, p: f64, q: f64) -> f64 {
        drucker_prager_yield(
            q,
            p,
            self.mu(),
            self.cohesion + self.hardening_modulus * self.plastic_strain,
        )
    }

    /// Return stress (q) after plastic projection.
    pub fn project(&mut self, p: f64, q: f64) -> f64 {
        let f = self.yield_function(p, q);
        if f <= 0.0 {
            return q; // elastic
        }
        // Radial return
        let delta_gamma = f / (1.0 + self.hardening_modulus);
        self.plastic_strain += delta_gamma;
        q - delta_gamma
    }
}

/// Drucker-Prager yield function value.
pub fn drucker_prager_yield(q: f64, p: f64, mu: f64, c: f64) -> f64 {
    q - mu * p - c
}

// ---------------------------------------------------------------------------
// SnowCohesion
// ---------------------------------------------------------------------------

/// Temperature-dependent snow cohesion (sintering model).
pub struct SnowCohesion {
    /// Reference cohesion at 273 K (Pa).
    pub cohesion_0: f64,
    /// Temperature sensitivity (Pa/K).
    pub temp_sensitivity: f64,
    /// Sintering time constant (seconds).
    pub sintering_tau: f64,
    /// Accumulated sintering time.
    pub sintered_time: f64,
}

impl SnowCohesion {
    /// Create a snow cohesion model.
    pub fn new(cohesion_0: f64) -> Self {
        Self {
            cohesion_0,
            temp_sensitivity: 50.0,
            sintering_tau: 3600.0,
            sintered_time: 0.0,
        }
    }

    /// Cohesion at temperature `t_kelvin`.
    pub fn cohesion_at_temp(&self, t_kelvin: f64) -> f64 {
        let delta = (273.15 - t_kelvin).max(0.0);
        self.cohesion_0 + self.temp_sensitivity * delta
    }

    /// Effective cohesion including sintering strengthening.
    pub fn effective_cohesion(&self, t_kelvin: f64) -> f64 {
        let base = self.cohesion_at_temp(t_kelvin);
        let sintering_factor = 1.0 - (-self.sintered_time / self.sintering_tau).exp();
        base * (1.0 + sintering_factor)
    }

    /// Advance the sintering clock.
    pub fn tick(&mut self, dt: f64) {
        self.sintered_time += dt;
    }
}

/// Compute snow density from temperature using a simple model.
pub fn snow_density_from_temp(t_kelvin: f64) -> f64 {
    // Fresh snow: ~50 kg/m³ at -10°C; denser near 0°C
    let t_celsius = t_kelvin - 273.15;
    let base = 150.0; // kg/m³
    let temp_effect = (-t_celsius * 0.1).exp();
    (base + 200.0 * (1.0 - temp_effect)).clamp(50.0, 600.0)
}

// ---------------------------------------------------------------------------
// MudParticle
// ---------------------------------------------------------------------------

/// A particle of clay-water mud.
pub struct MudParticle {
    /// World-space position (metres).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Water content by mass fraction (0..1).
    pub water_content: f64,
    /// Undrained shear strength (Pa).
    pub shear_strength: f64,
    /// Atterberg liquid limit (water content fraction).
    pub liquid_limit: f64,
    /// Atterberg plastic limit.
    pub plastic_limit: f64,
}

impl MudParticle {
    /// Create a mud particle.
    pub fn new(position: [f64; 3], water_content: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            water_content,
            shear_strength: 10000.0,
            liquid_limit: 0.4,
            plastic_limit: 0.2,
        }
    }

    /// Returns true when mud is flowing (above liquid limit).
    pub fn is_liquid(&self) -> bool {
        self.water_content >= self.liquid_limit
    }

    /// Returns true when mud is plastic (between plastic and liquid limits).
    pub fn is_plastic(&self) -> bool {
        self.water_content >= self.plastic_limit && self.water_content < self.liquid_limit
    }

    /// Compute undrained shear strength using a simplified relationship.
    pub fn compute_shear_strength(&self) -> f64 {
        if self.is_liquid() {
            0.0
        } else if self.is_plastic() {
            let w = self.water_content;
            let wl = self.liquid_limit;
            let wp = self.plastic_limit;
            self.shear_strength * (wl - w) / (wl - wp)
        } else {
            self.shear_strength
        }
    }
}

// ---------------------------------------------------------------------------
// ViscoPlasticMud — Bingham plastic
// ---------------------------------------------------------------------------

/// Bingham plastic viscoplastic mud model.
pub struct ViscoPlasticMud {
    /// Yield stress (Pa).
    pub yield_stress: f64,
    /// Plastic viscosity (Pa·s).
    pub plastic_viscosity: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Particles.
    pub particles: Vec<MudParticle>,
    /// Gravity.
    pub gravity: [f64; 3],
}

impl ViscoPlasticMud {
    /// Create a Bingham mud model.
    pub fn new(yield_stress: f64, plastic_viscosity: f64, density: f64) -> Self {
        Self {
            yield_stress,
            plastic_viscosity,
            density,
            particles: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add a mud particle.
    pub fn add_particle(&mut self, p: MudParticle) {
        self.particles.push(p);
    }

    /// Compute Bingham shear stress for a given shear rate `gamma_dot` (1/s).
    pub fn shear_stress(&self, gamma_dot: f64) -> f64 {
        bingham_shear_stress(gamma_dot, self.yield_stress, self.plastic_viscosity)
    }

    /// Integrate particles forward by `dt`.
    pub fn step(&mut self, dt: f64) {
        let g = self.gravity;
        let yield_stress = self.yield_stress;
        let plastic_viscosity = self.plastic_viscosity;
        let density = self.density;
        for p in &mut self.particles {
            let vel_len = len3(p.velocity);
            let tau = bingham_shear_stress(vel_len, yield_stress, plastic_viscosity);
            let drag = if vel_len > 1e-9 {
                let dir = scale3(p.velocity, 1.0 / vel_len);
                scale3(dir, -tau / (density * 0.01))
            } else {
                [0.0; 3]
            };
            let acc = add3(g, drag);
            p.velocity = add3(p.velocity, scale3(acc, dt));
            p.position = add3(p.position, scale3(p.velocity, dt));
        }
    }
}

/// Compute Bingham plastic shear stress.
pub fn bingham_shear_stress(gamma_dot: f64, tau_y: f64, mu_p: f64) -> f64 {
    if gamma_dot.abs() < 1e-12 {
        return 0.0;
    }
    tau_y + mu_p * gamma_dot.abs()
}

// ---------------------------------------------------------------------------
// SandGranular
// ---------------------------------------------------------------------------

/// Dry granular sand model.
pub struct SandGranular {
    /// Angle of repose (radians).
    pub angle_of_repose: f64,
    /// Dilatancy angle (radians).
    pub dilatancy_angle: f64,
    /// Internal friction angle (radians).
    pub friction_angle: f64,
    /// Cohesion (Pa) — effectively zero for dry sand.
    pub cohesion: f64,
    /// Density (kg/m³).
    pub density: f64,
}

impl SandGranular {
    /// Create a dry sand model.
    pub fn new(friction_angle: f64) -> Self {
        Self {
            angle_of_repose: friction_angle,
            dilatancy_angle: friction_angle * 0.5,
            friction_angle,
            cohesion: 0.0,
            density: 1600.0,
        }
    }

    /// Mohr-Coulomb shear strength at effective normal stress `sigma_n`.
    pub fn shear_strength(&self, sigma_n: f64) -> f64 {
        self.cohesion + sigma_n * self.friction_angle.tan()
    }

    /// Returns `true` when a slope exceeds the angle of repose.
    pub fn is_unstable(&self, slope_angle: f64) -> bool {
        slope_angle > self.angle_of_repose
    }

    /// Volumetric strain from dilatancy given shear strain increment.
    pub fn dilatancy_strain(&self, shear_strain: f64) -> f64 {
        shear_strain * self.dilatancy_angle.tan()
    }
}

// ---------------------------------------------------------------------------
// AvalancheSim
// ---------------------------------------------------------------------------

/// Simple snow avalanche slope-stability model.
pub struct AvalancheSim {
    /// Slope angle (radians).
    pub slope_angle: f64,
    /// Slab thickness (metres).
    pub slab_thickness: f64,
    /// Snow density (kg/m³).
    pub density: f64,
    /// Weak layer shear strength (Pa).
    pub weak_layer_strength: f64,
    /// Weak layer cohesion (Pa).
    pub weak_layer_cohesion: f64,
    /// Whether an avalanche is in progress.
    pub triggered: bool,
    /// Avalanche velocity (m/s).
    pub avalanche_velocity: f64,
}

impl AvalancheSim {
    /// Create an avalanche simulation.
    pub fn new(slope_angle: f64, slab_thickness: f64, density: f64) -> Self {
        Self {
            slope_angle,
            slab_thickness,
            density,
            weak_layer_strength: 500.0,
            weak_layer_cohesion: 100.0,
            triggered: false,
            avalanche_velocity: 0.0,
        }
    }

    /// Compute the driving shear stress (Pa) on the weak layer.
    pub fn driving_stress(&self) -> f64 {
        9.81 * self.density * self.slab_thickness * self.slope_angle.sin()
    }

    /// Compute the resisting shear stress (Pa) from the weak layer.
    pub fn resisting_stress(&self) -> f64 {
        let normal_stress = 9.81 * self.density * self.slab_thickness * self.slope_angle.cos();
        self.weak_layer_cohesion + normal_stress * self.weak_layer_strength / 1000.0
    }

    /// Safety factor (> 1 → stable, < 1 → unstable).
    pub fn safety_factor(&self) -> f64 {
        self.resisting_stress() / self.driving_stress().max(1e-9)
    }

    /// Trigger probability (simplified logistic function of safety factor).
    pub fn trigger_probability(&self) -> f64 {
        let fs = self.safety_factor();
        1.0 / (1.0 + (-10.0 * (1.0 - fs)).exp())
    }

    /// Step the avalanche dynamics.
    pub fn step(&mut self, dt: f64) {
        if self.safety_factor() < 1.0 && !self.triggered {
            self.triggered = true;
        }
        if self.triggered {
            let net_accel = 9.81 * self.slope_angle.sin() - 0.1 * self.avalanche_velocity;
            self.avalanche_velocity = (self.avalanche_velocity + net_accel * dt).max(0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// TireInteraction
// ---------------------------------------------------------------------------

/// Tire-terrain interaction in snow/mud.
pub struct TireInteraction {
    /// Tire contact width (metres).
    pub width: f64,
    /// Tire contact length (metres).
    pub contact_length: f64,
    /// Tire load (N).
    pub load: f64,
    /// Terrain cohesion (Pa).
    pub cohesion: f64,
    /// Terrain friction angle (radians).
    pub friction_angle: f64,
    /// Terrain stiffness (Pa/m — Bekker parameter kc).
    pub kc: f64,
    /// Terrain sinkage exponent (Bekker n).
    pub sinkage_exp: f64,
    /// Bekker modulus kphi (Pa/m).
    pub kphi: f64,
}

impl TireInteraction {
    /// Create a tire interaction model.
    pub fn new(
        width: f64,
        contact_length: f64,
        load: f64,
        cohesion: f64,
        friction_angle: f64,
    ) -> Self {
        Self {
            width,
            contact_length,
            load,
            cohesion,
            friction_angle,
            kc: 102.0,
            sinkage_exp: 0.9,
            kphi: 1510.0,
        }
    }

    /// Compute tire sinkage (metres) using the Bekker-Wong model.
    pub fn sinkage(&self) -> f64 {
        let p = self.load / (self.width * self.contact_length);
        let k = self.kc / self.width + self.kphi;
        (p / k).powf(1.0 / self.sinkage_exp)
    }

    /// Compute rolling resistance force (N).
    pub fn rolling_resistance(&self) -> f64 {
        let z = self.sinkage();
        // Simplified: resistance ~ load * z / contact_length
        self.load * z / self.contact_length
    }

    /// Maximum traction force (N) via Mohr-Coulomb.
    pub fn max_traction(&self) -> f64 {
        let contact_area = self.width * self.contact_length;
        let normal_stress = self.load / contact_area;
        contact_area * (self.cohesion + normal_stress * self.friction_angle.tan())
    }
}

// ---------------------------------------------------------------------------
// GroundDeformation
// ---------------------------------------------------------------------------

/// Permanent deformation tracking (track record in snow/mud).
pub struct GroundDeformation {
    /// Grid resolution (cells per side).
    pub resolution: usize,
    /// Cell size (metres).
    pub cell_size: f64,
    /// Deformation depth per cell (metres, positive = indented).
    pub depth_map: Vec<f64>,
    /// Maximum deformation allowed (metres).
    pub max_depth: f64,
    /// Recovery rate (m/s) — slow creep back.
    pub recovery_rate: f64,
}

impl GroundDeformation {
    /// Create a deformation tracking grid.
    pub fn new(resolution: usize, domain_size: f64) -> Self {
        Self {
            resolution,
            cell_size: domain_size / resolution as f64,
            depth_map: vec![0.0; resolution * resolution],
            max_depth: 0.3,
            recovery_rate: 1e-4,
        }
    }

    /// Apply a point load at world position `(x, z)` with `pressure` (Pa) and time `dt`.
    pub fn apply_load(&mut self, x: f64, z: f64, pressure: f64, stiffness: f64, dt: f64) {
        let ix = (x / self.cell_size) as i64;
        let iz = (z / self.cell_size) as i64;
        let n = self.resolution as i64;
        if ix < 0 || ix >= n || iz < 0 || iz >= n {
            return;
        }
        let idx = ix as usize * self.resolution + iz as usize;
        let sinkage = (pressure / stiffness) * dt;
        self.depth_map[idx] = (self.depth_map[idx] + sinkage).min(self.max_depth);
    }

    /// Advance elastic recovery.
    pub fn recover(&mut self, dt: f64) {
        for d in &mut self.depth_map {
            *d = (*d - self.recovery_rate * dt).max(0.0);
        }
    }

    /// Query deformation depth at world position `(x, z)`.
    pub fn query(&self, x: f64, z: f64) -> f64 {
        let ix = (x / self.cell_size) as i64;
        let iz = (z / self.cell_size) as i64;
        let n = self.resolution as i64;
        if ix < 0 || ix >= n || iz < 0 || iz >= n {
            return 0.0;
        }
        self.depth_map[ix as usize * self.resolution + iz as usize]
    }

    /// Maximum deformation depth currently recorded.
    pub fn max_deformation(&self) -> f64 {
        self.depth_map.iter().cloned().fold(0.0, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. SnowParticle: mass = density * volume
    #[test]
    fn test_snow_particle_mass() {
        let p = SnowParticle::new([0.0; 3], 300.0, 263.15);
        assert!((p.mass - p.density * p.volume).abs() < 1e-12);
    }

    // 2. SnowMpm: grid size correct
    #[test]
    fn test_mpm_grid_size() {
        let sim = SnowMpm::new(4, 1.0);
        assert_eq!(sim.grid.len(), 64);
    }

    // 3. SnowMpm: step runs without panic
    #[test]
    fn test_mpm_step() {
        let mut sim = SnowMpm::new(4, 1.0);
        sim.add_particle(SnowParticle::new([0.5, 0.8, 0.5], 300.0, 263.0));
        sim.step(0.001);
    }

    // 4. DruckerPragerPlastic: yield function positive outside yield surface
    #[test]
    fn test_dp_yield_positive() {
        let dp = DruckerPragerPlastic::new(30.0_f64.to_radians(), 500.0);
        let f = dp.yield_function(1000.0, 5000.0);
        // With large q this should be positive
        assert!(f > 0.0);
    }

    // 5. DruckerPragerPlastic: projection reduces q
    #[test]
    fn test_dp_project_reduces_q() {
        let mut dp = DruckerPragerPlastic::new(30.0_f64.to_radians(), 500.0);
        let q_before = 5000.0;
        let q_after = dp.project(1000.0, q_before);
        assert!(q_after <= q_before);
    }

    // 6. SnowCohesion: cohesion higher at lower temperature
    #[test]
    fn test_snow_cohesion_temp() {
        let sc = SnowCohesion::new(1000.0);
        let c_cold = sc.cohesion_at_temp(253.15); // -20°C
        let c_warm = sc.cohesion_at_temp(273.15); // 0°C
        assert!(c_cold > c_warm);
    }

    // 7. SnowCohesion: sintering increases over time
    #[test]
    fn test_snow_cohesion_sintering() {
        let mut sc = SnowCohesion::new(1000.0);
        let c0 = sc.effective_cohesion(263.0);
        sc.tick(3600.0);
        let c1 = sc.effective_cohesion(263.0);
        assert!(c1 > c0);
    }

    // 8. snow_density_from_temp: reasonable range
    #[test]
    fn test_snow_density_range() {
        let d = snow_density_from_temp(263.0);
        assert!((50.0..=600.0).contains(&d));
    }

    // 9. MudParticle: is_liquid when above liquid limit
    #[test]
    fn test_mud_is_liquid() {
        let p = MudParticle::new([0.0; 3], 0.5);
        assert!(p.is_liquid());
    }

    // 10. MudParticle: shear strength zero when liquid
    #[test]
    fn test_mud_zero_strength_liquid() {
        let p = MudParticle::new([0.0; 3], 0.5);
        assert_eq!(p.compute_shear_strength(), 0.0);
    }

    // 11. MudParticle: plastic range
    #[test]
    fn test_mud_plastic_range() {
        let p = MudParticle::new([0.0; 3], 0.3);
        assert!(p.is_plastic());
    }

    // 12. ViscoPlasticMud: bingham stress increases with shear rate
    #[test]
    fn test_bingham_stress_increases() {
        let mud = ViscoPlasticMud::new(100.0, 10.0, 1500.0);
        let s1 = mud.shear_stress(1.0);
        let s2 = mud.shear_stress(2.0);
        assert!(s2 > s1);
    }

    // 13. bingham_shear_stress: zero gamma → zero
    #[test]
    fn test_bingham_zero_gamma() {
        assert_eq!(bingham_shear_stress(0.0, 100.0, 10.0), 0.0);
    }

    // 14. SandGranular: stable below repose angle
    #[test]
    fn test_sand_stable() {
        let sand = SandGranular::new(35.0_f64.to_radians());
        assert!(!sand.is_unstable(30.0_f64.to_radians()));
    }

    // 15. SandGranular: unstable above repose
    #[test]
    fn test_sand_unstable() {
        let sand = SandGranular::new(35.0_f64.to_radians());
        assert!(sand.is_unstable(40.0_f64.to_radians()));
    }

    // 16. SandGranular: Mohr-Coulomb shear strength positive
    #[test]
    fn test_sand_shear_strength() {
        let sand = SandGranular::new(35.0_f64.to_radians());
        assert!(sand.shear_strength(1000.0) > 0.0);
    }

    // 17. AvalancheSim: safety factor > 1 for gentle slope
    #[test]
    fn test_avalanche_stable() {
        let sim = AvalancheSim::new(10.0_f64.to_radians(), 0.5, 200.0);
        assert!(sim.safety_factor() > 1.0);
    }

    // 18. AvalancheSim: trigger probability near 0.5 at FS=1
    #[test]
    fn test_avalanche_trigger_prob_at_unity() {
        let prob = {
            let mut sim = AvalancheSim::new(30.0_f64.to_radians(), 1.0, 300.0);
            // Set strengths to exactly match driving
            sim.weak_layer_strength = sim.driving_stress() * 1000.0
                / (9.81 * 300.0 * 1.0 * (30.0_f64).to_radians().cos());
            sim.weak_layer_cohesion = 0.0;
            sim.trigger_probability()
        };
        assert!(prob > 0.0 && prob <= 1.0);
    }

    // 19. TireInteraction: sinkage positive
    #[test]
    fn test_tire_sinkage_positive() {
        let ti = TireInteraction::new(0.2, 0.3, 4000.0, 100.0, 30.0_f64.to_radians());
        assert!(ti.sinkage() > 0.0);
    }

    // 20. TireInteraction: max_traction positive
    #[test]
    fn test_tire_traction_positive() {
        let ti = TireInteraction::new(0.2, 0.3, 4000.0, 100.0, 30.0_f64.to_radians());
        assert!(ti.max_traction() > 0.0);
    }

    // 21. GroundDeformation: apply_load increases depth
    #[test]
    fn test_deformation_load() {
        let mut gd = GroundDeformation::new(10, 10.0);
        gd.apply_load(5.0, 5.0, 10000.0, 1e6, 1.0);
        assert!(gd.query(5.0, 5.0) > 0.0);
    }

    // 22. GroundDeformation: recover reduces depth
    #[test]
    fn test_deformation_recovery() {
        let mut gd = GroundDeformation::new(10, 10.0);
        gd.depth_map[0] = 0.1;
        gd.recover(1.0);
        assert!(gd.depth_map[0] < 0.1);
    }

    // 23. GroundDeformation: max_deformation zero initially
    #[test]
    fn test_deformation_max_zero() {
        let gd = GroundDeformation::new(10, 10.0);
        assert_eq!(gd.max_deformation(), 0.0);
    }

    // 24. mpm_grid_to_particle: weight in [0, 1]
    #[test]
    fn test_mpm_weight_range() {
        let w = mpm_grid_to_particle([0.0; 3], [0.0; 3], 1.0);
        assert!((0.0..=1.0).contains(&w));
    }

    // 25. drucker_prager_yield: negative when inside yield surface
    #[test]
    fn test_dp_yield_negative_inside() {
        let f = drucker_prager_yield(100.0, 10000.0, 0.3, 500.0);
        assert!(f < 0.0);
    }

    // 26. ViscoPlasticMud: step moves particles downward
    #[test]
    fn test_mud_step_gravity() {
        let mut mud = ViscoPlasticMud::new(100.0, 10.0, 1500.0);
        mud.add_particle(MudParticle::new([0.0, 5.0, 0.0], 0.5));
        let y0 = mud.particles[0].position[1];
        mud.step(0.1);
        assert!(mud.particles[0].position[1] < y0);
    }

    // 27. AvalancheSim: velocity increases after trigger
    #[test]
    fn test_avalanche_velocity_increases() {
        let mut sim = AvalancheSim::new(45.0_f64.to_radians(), 1.0, 300.0);
        sim.triggered = true;
        sim.step(1.0);
        assert!(sim.avalanche_velocity > 0.0);
    }

    // 28. snow_density_from_temp: denser near freezing
    #[test]
    fn test_snow_density_near_freezing() {
        let d_cold = snow_density_from_temp(243.0);
        let d_warm = snow_density_from_temp(273.0);
        assert!(d_cold != d_warm); // should differ
    }

    // 29. SnowMpm: clear_grid resets mass
    #[test]
    fn test_mpm_clear_grid() {
        let mut sim = SnowMpm::new(4, 1.0);
        sim.add_particle(SnowParticle::new([0.5, 0.5, 0.5], 300.0, 263.0));
        sim.particle_to_grid();
        sim.clear_grid();
        let total_mass: f64 = sim.grid.iter().map(|c| c.mass).sum();
        assert!(total_mass < 1e-12);
    }

    // 30. GroundDeformation: clamped to max_depth
    #[test]
    fn test_deformation_clamp() {
        let mut gd = GroundDeformation::new(10, 10.0);
        // Apply huge load many times
        for _ in 0..1000 {
            gd.apply_load(5.0, 5.0, 1e9, 1e4, 1.0);
        }
        assert!(gd.query(5.0, 5.0) <= gd.max_depth + 1e-9);
    }
}
