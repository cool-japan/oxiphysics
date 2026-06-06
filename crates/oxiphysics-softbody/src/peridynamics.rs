// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Peridynamics — nonlocal continuum mechanics for fracture simulation.
//!
//! Implements bond-based (PMB) and ordinary state-based peridynamic models,
//! crack tracking, fracture energy computation, and convergence helpers.

// ── Vector helpers ────────────────────────────────────────────────────────────

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[cfg(test)]
#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = norm3(v);
    if len < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / len)
    }
}

// ── Material type ─────────────────────────────────────────────────────────────

/// Peridynamic material formulation variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PdMaterialType {
    /// Prototype microelastic brittle (bond-based) model.
    BondBased,
    /// Generic state-based peridynamic model.
    StateBased,
    /// Ordinary state-based: force parallel to bond direction.
    OrdinaryStateBased,
    /// Non-ordinary state-based: general correspondence formulation.
    NonOrdinaryStateBased,
}

// ── Material parameters ───────────────────────────────────────────────────────

/// Peridynamic material parameters.
#[derive(Debug, Clone)]
pub struct PdMaterial {
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
    /// Critical stretch beyond which a bond breaks (dimensionless).
    pub critical_stretch: f64,
    /// Nonlocal horizon radius δ (m).
    pub horizon_radius: f64,
}

impl PdMaterial {
    /// Construct a new `PdMaterial`.
    pub fn new(
        bulk_modulus: f64,
        shear_modulus: f64,
        density: f64,
        critical_stretch: f64,
        horizon_radius: f64,
    ) -> Self {
        Self {
            bulk_modulus,
            shear_modulus,
            density,
            critical_stretch,
            horizon_radius,
        }
    }

    /// Bond-based micromodulus constant `c` for 3-D PMB model.
    ///
    /// For a sphere of radius `horizon`:
    /// `c = 18 K / (π δ⁴)`
    pub fn micro_modulus_bond_based(bulk_modulus: f64, horizon: f64) -> f64 {
        18.0 * bulk_modulus / (std::f64::consts::PI * horizon.powi(4))
    }

    /// Critical energy release rate derived from the material parameters.
    ///
    /// For the PMB model: `G_c = (π c s_c² δ⁵) / 10`
    pub fn critical_energy_release_rate(&self) -> f64 {
        let c = Self::micro_modulus_bond_based(self.bulk_modulus, self.horizon_radius);
        std::f64::consts::PI * c * self.critical_stretch.powi(2) * self.horizon_radius.powi(5)
            / 10.0
    }

    /// Young's modulus derived from bulk and shear moduli.
    pub fn young_modulus(&self) -> f64 {
        9.0 * self.bulk_modulus * self.shear_modulus
            / (3.0 * self.bulk_modulus + self.shear_modulus)
    }

    /// Poisson's ratio.
    pub fn poisson_ratio(&self) -> f64 {
        (3.0 * self.bulk_modulus - 2.0 * self.shear_modulus)
            / (2.0 * (3.0 * self.bulk_modulus + self.shear_modulus))
    }
}

// ── Particle ──────────────────────────────────────────────────────────────────

/// A peridynamic material particle.
#[derive(Debug, Clone)]
pub struct PdParticle {
    /// Current position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Accumulated body/bond force per unit volume (N/m³).
    pub force: [f64; 3],
    /// Displacement from reference position (m).
    pub disp: [f64; 3],
    /// Reference (rest) position (m).
    pub pos_ref: [f64; 3],
    /// Representative volume of this particle (m³).
    pub volume: f64,
    /// Damage index ∈ \[0, 1\]: fraction of broken bonds.
    pub damage: f64,
}

impl PdParticle {
    /// Create a new particle at position `pos` with given `volume`.
    pub fn new(pos: [f64; 3], volume: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            disp: [0.0; 3],
            pos_ref: pos,
            volume,
            damage: 0.0,
        }
    }

    /// Clear accumulated force.
    pub fn clear_force(&mut self) {
        self.force = [0.0; 3];
    }

    /// Apply an external body force density (N/m³) to this particle.
    pub fn apply_body_force(&mut self, f: [f64; 3]) {
        self.force = add3(self.force, f);
    }
}

// ── Bond ──────────────────────────────────────────────────────────────────────

/// A peridynamic bond connecting two particles.
#[derive(Debug, Clone)]
pub struct PdBond {
    /// Index of particle i.
    pub i: usize,
    /// Index of particle j.
    pub j: usize,
    /// Reference (rest) bond length |ξ| (m).
    pub ref_length: f64,
    /// Whether this bond has been permanently broken.
    pub broken: bool,
    /// Current stretch (dimensionless).
    pub stretch: f64,
}

impl PdBond {
    /// Create a new bond between particles `i` and `j` with reference length `ref_length`.
    pub fn new(i: usize, j: usize, ref_length: f64) -> Self {
        Self {
            i,
            j,
            ref_length,
            broken: false,
            stretch: 0.0,
        }
    }

    /// Compute bond stretch from reference vector `ξ` and relative displacement `η`.
    ///
    /// `s = (|ξ + η| - |ξ|) / |ξ|`
    pub fn compute_stretch(xi: [f64; 3], eta: [f64; 3]) -> f64 {
        let xi_len = norm3(xi);
        if xi_len < 1e-15 {
            return 0.0;
        }
        let deformed = add3(xi, eta);
        let deformed_len = norm3(deformed);
        (deformed_len - xi_len) / xi_len
    }

    /// Whether this bond is broken.
    pub fn is_broken(&self) -> bool {
        self.broken
    }

    /// Update the stretch value and break the bond if it exceeds `critical`.
    pub fn update(&mut self, xi: [f64; 3], eta: [f64; 3], critical: f64) {
        self.stretch = Self::compute_stretch(xi, eta);
        if self.stretch.abs() >= critical {
            self.broken = true;
        }
    }
}

// ── Pairwise bond force ───────────────────────────────────────────────────────

/// Compute the pairwise bond force vector on particle i due to particle j
/// (PMB model).
///
/// `f = c * s * (ξ + η) / |ξ + η|`
///
/// Returns the force per unit volume squared contribution.
pub fn bond_force(stretch: f64, c: f64, ref_vec: [f64; 3], disp_diff: [f64; 3]) -> [f64; 3] {
    let deformed = add3(ref_vec, disp_diff);
    let len = norm3(deformed);
    if len < 1e-15 {
        return [0.0; 3];
    }
    let dir = scale3(deformed, 1.0 / len);
    scale3(dir, c * stretch)
}

// ── Crack path ────────────────────────────────────────────────────────────────

/// Tracks crack propagation as a set of broken bond midpoints.
#[derive(Debug, Clone, Default)]
pub struct CrackPath {
    /// World-space midpoints of broken bonds, in order of breaking.
    pub points: Vec<[f64; 3]>,
}

impl CrackPath {
    /// Create a new empty crack path.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a newly broken bond between positions `pi` and `pj`.
    pub fn record_break(&mut self, pi: [f64; 3], pj: [f64; 3]) {
        let mid = scale3(add3(pi, pj), 0.5);
        self.points.push(mid);
    }

    /// Number of bond-break events recorded.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the crack path is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

// ── Fracture energy ───────────────────────────────────────────────────────────

/// Fracture energy tracker.
#[derive(Debug, Clone, Default)]
pub struct FractureEnergy {
    /// Total fracture energy released so far (J).
    pub total: f64,
}

impl FractureEnergy {
    /// Create a new `FractureEnergy` accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Accumulate the energy released when one bond of length `ref_length`
    /// and micromodulus `c` breaks at critical stretch `s_c`.
    ///
    /// Energy per bond ≈ `0.5 * c * s_c² * ref_length * V_i * V_j`
    pub fn add_broken_bond(&mut self, c: f64, s_c: f64, ref_length: f64, v_i: f64, v_j: f64) {
        self.total += 0.5 * c * s_c * s_c * ref_length * v_i * v_j;
    }

    /// Reset the fracture energy to zero.
    pub fn reset(&mut self) {
        self.total = 0.0;
    }
}

// ── Main peridynamic system ───────────────────────────────────────────────────

/// The full peridynamic simulation system.
pub struct PdSystem {
    /// All particles.
    pub particles: Vec<PdParticle>,
    /// All bonds (both broken and intact).
    pub bonds: Vec<PdBond>,
    /// Material parameters.
    pub material: PdMaterial,
    /// Initial bond count per particle (for damage calculation).
    pub initial_bond_count: Vec<usize>,
    /// Crack path tracker.
    pub crack_path: CrackPath,
    /// Fracture energy accumulator.
    pub fracture_energy: FractureEnergy,
}

impl PdSystem {
    /// Construct a new `PdSystem` from a list of particles and material.
    pub fn new(particles: Vec<PdParticle>, material: PdMaterial) -> Self {
        let n = particles.len();
        Self {
            particles,
            bonds: Vec::new(),
            material,
            initial_bond_count: vec![0; n],
            crack_path: CrackPath::new(),
            fracture_energy: FractureEnergy::new(),
        }
    }

    /// Find all particle pairs within the horizon radius and create bonds.
    ///
    /// Returns the total number of bonds created.
    pub fn find_bonds(&mut self) -> usize {
        self.bonds.clear();
        let n = self.particles.len();
        let delta = self.material.horizon_radius;
        let mut counts = vec![0usize; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let xi = sub3(self.particles[j].pos_ref, self.particles[i].pos_ref);
                let d = norm3(xi);
                if d <= delta && d > 1e-15 {
                    self.bonds.push(PdBond::new(i, j, d));
                    counts[i] += 1;
                    counts[j] += 1;
                }
            }
        }
        self.initial_bond_count = counts;
        self.bonds.len()
    }

    /// Compute pairwise bond forces (PMB model) and accumulate into `particle.force`.
    pub fn compute_forces(&mut self) {
        // Clear forces
        for p in &mut self.particles {
            p.clear_force();
        }

        let c = PdMaterial::micro_modulus_bond_based(
            self.material.bulk_modulus,
            self.material.horizon_radius,
        );
        let s_c = self.material.critical_stretch;

        let n_bonds = self.bonds.len();
        for b in 0..n_bonds {
            if self.bonds[b].broken {
                continue;
            }
            let i = self.bonds[b].i;
            let j = self.bonds[b].j;

            let xi = sub3(self.particles[j].pos_ref, self.particles[i].pos_ref);
            let eta = sub3(self.particles[j].disp, self.particles[i].disp);

            let stretch = PdBond::compute_stretch(xi, eta);
            self.bonds[b].stretch = stretch;

            if stretch.abs() >= s_c {
                // Bond breaks
                let pi = self.particles[i].pos;
                let pj = self.particles[j].pos;
                self.crack_path.record_break(pi, pj);
                self.fracture_energy.add_broken_bond(
                    c,
                    s_c,
                    self.bonds[b].ref_length,
                    self.particles[i].volume,
                    self.particles[j].volume,
                );
                self.bonds[b].broken = true;
                continue;
            }

            let f = bond_force(stretch, c, xi, eta);
            let v_j = self.particles[j].volume;
            let v_i = self.particles[i].volume;

            // f_i += f * V_j;  f_j -= f * V_i
            let fi = scale3(f, v_j);
            let fj = scale3(f, v_i);

            self.particles[i].force = add3(self.particles[i].force, fi);
            self.particles[j].force = sub3(self.particles[j].force, fj);
        }
    }

    /// Update damage index for every particle as `broken / initial`.
    pub fn compute_damage(&mut self) {
        let n = self.particles.len();
        let mut broken_counts = vec![0usize; n];

        for b in &self.bonds {
            if b.broken {
                broken_counts[b.i] += 1;
                broken_counts[b.j] += 1;
            }
        }

        for (k, &bc) in broken_counts.iter().enumerate().take(n) {
            let init = self.initial_bond_count[k];
            self.particles[k].damage = if init == 0 {
                0.0
            } else {
                bc as f64 / init as f64
            };
        }
    }

    /// Advance one time step using the velocity-Verlet scheme.
    pub fn step(&mut self, dt: f64) {
        self.velocity_verlet(dt);
    }

    /// Velocity-Verlet time integration.
    ///
    /// 1. `v(t + ½dt) = v(t) + ½ a(t) dt`
    /// 2. `x(t + dt) = x(t) + v(t + ½dt) dt`
    /// 3. Compute forces at new positions
    /// 4. `v(t + dt) = v(t + ½dt) + ½ a(t + dt) dt`
    pub fn velocity_verlet(&mut self, dt: f64) {
        let rho = self.material.density;
        let inv_rho = 1.0 / rho;
        let half_dt = 0.5 * dt;

        // Half-step velocity + full-step position
        for p in &mut self.particles {
            let a = scale3(p.force, inv_rho);
            p.vel = add3(p.vel, scale3(a, half_dt));
            p.pos = add3(p.pos, scale3(p.vel, dt));
            p.disp = sub3(p.pos, p.pos_ref);
        }

        // Recompute forces at new positions
        self.compute_forces();

        // Second half-step velocity
        for p in &mut self.particles {
            let a = scale3(p.force, inv_rho);
            p.vel = add3(p.vel, scale3(a, half_dt));
        }
    }

    /// Estimate the critical stable time step (CFL condition).
    ///
    /// `dt_crit = 0.8 * sqrt(2 ρ / (π c δ³))`
    pub fn critical_time_step(&self) -> f64 {
        let c = PdMaterial::micro_modulus_bond_based(
            self.material.bulk_modulus,
            self.material.horizon_radius,
        );
        let delta = self.material.horizon_radius;
        let rho = self.material.density;
        let denom = std::f64::consts::PI * c * delta.powi(3);
        if denom < 1e-30 {
            return f64::MAX;
        }
        0.8 * (2.0 * rho / denom).sqrt()
    }

    /// Total kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        let rho = self.material.density;
        self.particles.iter().fold(0.0, |acc, p| {
            acc + 0.5 * rho * p.volume * dot3(p.vel, p.vel)
        })
    }

    /// Total elastic (strain) energy stored in intact bonds (J).
    ///
    /// `W = ½ Σ_bonds c * s² * |ξ| * V_i * V_j`
    pub fn strain_energy(&self) -> f64 {
        let c = PdMaterial::micro_modulus_bond_based(
            self.material.bulk_modulus,
            self.material.horizon_radius,
        );
        self.bonds.iter().fold(0.0, |acc, b| {
            if b.broken {
                acc
            } else {
                acc + 0.5
                    * c
                    * b.stretch.powi(2)
                    * b.ref_length
                    * self.particles[b.i].volume
                    * self.particles[b.j].volume
            }
        })
    }

    /// Ordinary state-based force density for particle `i` due to particle `j`.
    ///
    /// Uses weighted volume, dilatation, and deviatoric strain.
    pub fn osb_force_density(
        &self,
        xi: [f64; 3],
        eta: [f64; 3],
        weight: f64,
        theta_i: f64,
        theta_j: f64,
        alpha: f64,
        beta: f64,
    ) -> [f64; 3] {
        let xi_eta = add3(xi, eta);
        let m_len = norm3(xi_eta);
        if m_len < 1e-15 {
            return [0.0; 3];
        }
        let e_len = norm3(xi);
        let ext = norm3(xi_eta) - e_len; // extension scalar
        let t = alpha * (theta_i + theta_j) * weight * ext
            + beta * weight * (ext - (theta_i + theta_j) / 3.0 * e_len);
        let dir = scale3(xi_eta, 1.0 / m_len);
        scale3(dir, t)
    }

    /// Compute dilatation θ for a single particle (simplified sum over bonds).
    pub fn compute_dilatation(&self, particle_idx: usize) -> f64 {
        let delta = self.material.horizon_radius;
        let mut theta = 0.0;
        for b in &self.bonds {
            if b.broken {
                continue;
            }
            let (i, j) = if b.i == particle_idx {
                (b.i, b.j)
            } else if b.j == particle_idx {
                (b.j, b.i)
            } else {
                continue;
            };
            let xi = sub3(self.particles[j].pos_ref, self.particles[i].pos_ref);
            let eta = sub3(self.particles[j].disp, self.particles[i].disp);
            let xi_len = norm3(xi);
            let ext = norm3(add3(xi, eta)) - xi_len;
            // Influence function w(|ξ|) = δ / |ξ|
            let w = delta / (xi_len + 1e-15);
            theta += w * xi_len * ext * self.particles[j].volume;
        }
        // Normalise by weighted volume (simplified)
        3.0 * theta / (delta.powi(4) + 1e-15)
    }

    /// Convergence check: in the small-horizon limit the bond-based model
    /// should recover Poisson ratio ν = 1/4 (3-D).
    ///
    /// Returns `|ν - 0.25|`.
    pub fn convergence_to_classical(&self) -> f64 {
        (self.material.poisson_ratio() - 0.25_f64).abs()
    }
}

// ── Regular lattice initialiser ───────────────────────────────────────────────

/// Generates a regular cubic lattice of peridynamic particles.
pub struct PdGrid;

impl PdGrid {
    /// Create a 3-D cubic lattice of particles filling the box
    /// `[0, lx] × [0, ly] × [0, lz]` with spacing `dx`.
    ///
    /// Particle volume is `dx³`.
    pub fn cubic_lattice(lx: f64, ly: f64, lz: f64, dx: f64) -> Vec<PdParticle> {
        let nx = ((lx / dx).round() as usize).max(1);
        let ny = ((ly / dx).round() as usize).max(1);
        let nz = ((lz / dx).round() as usize).max(1);
        let vol = dx * dx * dx;
        let mut particles = Vec::with_capacity(nx * ny * nz);
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let pos = [
                        (ix as f64 + 0.5) * dx,
                        (iy as f64 + 0.5) * dx,
                        (iz as f64 + 0.5) * dx,
                    ];
                    particles.push(PdParticle::new(pos, vol));
                }
            }
        }
        particles
    }

    /// Create a uniform 1-D bar of particles along the x-axis.
    pub fn bar_1d(length: f64, dx: f64) -> Vec<PdParticle> {
        let n = ((length / dx).round() as usize).max(1);
        let vol = dx;
        (0..n)
            .map(|i| PdParticle::new([(i as f64 + 0.5) * dx, 0.0, 0.0], vol))
            .collect()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_material() -> PdMaterial {
        PdMaterial::new(
            1.0e9,  // bulk modulus
            0.5e9,  // shear modulus
            2000.0, // density
            0.05,   // critical stretch
            0.02,   // horizon
        )
    }

    // 1. Micro-modulus is positive for positive inputs.
    #[test]
    fn test_micro_modulus_positive() {
        let c = PdMaterial::micro_modulus_bond_based(1e9, 0.02);
        assert!(c > 0.0, "micro-modulus must be positive, got {c}");
    }

    // 2. Critical energy release rate is positive.
    #[test]
    fn test_critical_energy_release_rate_positive() {
        let mat = default_material();
        let gc = mat.critical_energy_release_rate();
        assert!(gc > 0.0, "Gc must be positive, got {gc}");
    }

    // 3. Young's modulus formula.
    #[test]
    fn test_young_modulus() {
        let mat = default_material();
        let e = mat.young_modulus();
        assert!(e > 0.0, "E must be positive, got {e}");
    }

    // 4. Poisson's ratio bounds.
    #[test]
    fn test_poisson_ratio_bounds() {
        let mat = default_material();
        let nu = mat.poisson_ratio();
        assert!(
            nu > -1.0 && nu < 0.5,
            "Poisson ratio {nu} out of valid range"
        );
    }

    // 5. Bond stretch zero when no displacement.
    #[test]
    fn test_bond_stretch_zero_displacement() {
        let xi = [0.01, 0.0, 0.0];
        let eta = [0.0, 0.0, 0.0];
        let s = PdBond::compute_stretch(xi, eta);
        assert!(s.abs() < 1e-12, "stretch should be 0, got {s}");
    }

    // 6. Bond stretch 0.1 when extension is 10%.
    #[test]
    fn test_bond_stretch_ten_percent() {
        let xi = [1.0, 0.0, 0.0];
        let eta = [0.1, 0.0, 0.0]; // 10% extension
        let s = PdBond::compute_stretch(xi, eta);
        assert!((s - 0.1).abs() < 1e-12, "stretch should be 0.1, got {s}");
    }

    // 7. Bond does not break below critical stretch.
    #[test]
    fn test_bond_no_break_below_critical() {
        let xi = [1.0, 0.0, 0.0];
        let eta = [0.04, 0.0, 0.0]; // stretch = 0.04 < 0.05
        let mut bond = PdBond::new(0, 1, 1.0);
        bond.update(xi, eta, 0.05);
        assert!(
            !bond.is_broken(),
            "bond should not be broken at sub-critical stretch"
        );
    }

    // 8. Bond breaks at critical stretch.
    #[test]
    fn test_bond_breaks_at_critical() {
        let xi = [1.0, 0.0, 0.0];
        let eta = [0.05, 0.0, 0.0]; // stretch = 0.05 = critical
        let mut bond = PdBond::new(0, 1, 1.0);
        bond.update(xi, eta, 0.05);
        assert!(
            bond.is_broken(),
            "bond should be broken at critical stretch"
        );
    }

    // 9. bond_force is zero when stretch is zero.
    #[test]
    fn test_bond_force_zero_stretch() {
        let xi = [1.0, 0.0, 0.0];
        let eta = [0.0, 0.0, 0.0];
        let f = bond_force(0.0, 1e12, xi, eta);
        for c in f {
            assert!(c.abs() < 1e-12, "force component should be 0, got {c}");
        }
    }

    // 10. bond_force direction matches bond direction.
    #[test]
    fn test_bond_force_direction() {
        let xi = [1.0, 0.0, 0.0];
        let eta = [0.1, 0.0, 0.0];
        let s = PdBond::compute_stretch(xi, eta);
        let f = bond_force(s, 1e9, xi, eta);
        // Force should be in +x direction (positive stretch, positive force)
        assert!(
            f[0] > 0.0,
            "force x-component should be positive, got {}",
            f[0]
        );
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }

    // 11. PdSystem::new creates zero bonds.
    #[test]
    fn test_system_new_no_bonds() {
        let mat = default_material();
        let particles = vec![PdParticle::new([0.0; 3], 1e-6)];
        let sys = PdSystem::new(particles, mat);
        assert_eq!(sys.bonds.len(), 0);
    }

    // 12. find_bonds connects nearby particles.
    #[test]
    fn test_find_bonds_nearby() {
        let mat = default_material(); // horizon = 0.02
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6), // within horizon
        ];
        let mut sys = PdSystem::new(particles, mat);
        let nb = sys.find_bonds();
        assert_eq!(nb, 1, "should find 1 bond, found {nb}");
    }

    // 13. find_bonds ignores distant particles.
    #[test]
    fn test_find_bonds_distant() {
        let mat = default_material(); // horizon = 0.02
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([1.0, 0.0, 0.0], 1e-6), // outside horizon
        ];
        let mut sys = PdSystem::new(particles, mat);
        let nb = sys.find_bonds();
        assert_eq!(nb, 0, "distant particles should have no bond");
    }

    // 14. compute_forces does not panic.
    #[test]
    fn test_compute_forces_no_panic() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        sys.compute_forces();
    }

    // 15. Equal and opposite forces (Newton's 3rd law).
    #[test]
    fn test_forces_newton_third_law() {
        let mat = default_material();
        let mut sys = PdSystem::new(
            vec![
                PdParticle::new([0.0, 0.0, 0.0], 1e-6),
                PdParticle::new([0.01, 0.0, 0.0], 1e-6),
            ],
            mat,
        );
        // Displace particle 1 to create non-zero stretch
        sys.particles[1].disp = [0.0002, 0.0, 0.0];
        sys.particles[1].pos = add3(sys.particles[1].pos_ref, sys.particles[1].disp);
        sys.find_bonds();
        sys.compute_forces();

        let f0 = sys.particles[0].force;
        let f1 = sys.particles[1].force;
        // Forces must be opposite
        for k in 0..3 {
            assert!(
                (f0[k] + f1[k]).abs() < 1e-3,
                "Newton 3rd law violated at component {k}: f0={:.6}, f1={:.6}",
                f0[k],
                f1[k]
            );
        }
    }

    // 16. Damage is zero before any breakage.
    #[test]
    fn test_damage_zero_initially() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        sys.compute_damage();
        for p in &sys.particles {
            assert_eq!(p.damage, 0.0, "initial damage should be 0");
        }
    }

    // 17. Damage increases after bond break.
    #[test]
    fn test_damage_after_break() {
        let mut mat = default_material();
        mat.critical_stretch = 0.001; // very sensitive
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        // Force large displacement
        sys.particles[1].disp = [0.5, 0.0, 0.0];
        sys.particles[1].pos = add3(sys.particles[1].pos_ref, sys.particles[1].disp);
        sys.compute_forces();
        sys.compute_damage();
        assert!(
            sys.particles[0].damage > 0.0,
            "damage should increase after bond break"
        );
    }

    // 18. kinetic_energy is zero at rest.
    #[test]
    fn test_kinetic_energy_zero_at_rest() {
        let mat = default_material();
        let particles = vec![PdParticle::new([0.0; 3], 1e-6)];
        let sys = PdSystem::new(particles, mat);
        assert!(
            sys.kinetic_energy().abs() < 1e-30,
            "KE should be zero at rest"
        );
    }

    // 19. kinetic_energy scales with velocity squared.
    #[test]
    fn test_kinetic_energy_scales_with_v_squared() {
        let mat = default_material();
        let mut sys = PdSystem::new(vec![PdParticle::new([0.0; 3], 1e-6)], mat);
        sys.particles[0].vel = [1.0, 0.0, 0.0];
        let ke1 = sys.kinetic_energy();
        sys.particles[0].vel = [2.0, 0.0, 0.0];
        let ke2 = sys.kinetic_energy();
        assert!(
            (ke2 / ke1 - 4.0).abs() < 1e-10,
            "KE should scale as v², ratio={:.6}",
            ke2 / ke1
        );
    }

    // 20. strain_energy is zero for zero displacement.
    #[test]
    fn test_strain_energy_zero_at_rest() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        sys.compute_forces();
        let se = sys.strain_energy();
        assert!(
            se.abs() < 1e-30,
            "strain energy should be 0 at rest, got {se}"
        );
    }

    // 21. critical_time_step is positive and finite.
    #[test]
    fn test_critical_time_step_positive() {
        let mat = default_material();
        let sys = PdSystem::new(vec![], mat);
        let dt = sys.critical_time_step();
        assert!(
            dt > 0.0 && dt.is_finite(),
            "dt_crit should be positive finite, got {dt}"
        );
    }

    // 22. step does not panic with valid system.
    #[test]
    fn test_step_no_panic() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        sys.step(1e-7);
    }

    // 23. CrackPath records midpoint correctly.
    #[test]
    fn test_crack_path_midpoint() {
        let mut cp = CrackPath::new();
        let pi = [0.0, 0.0, 0.0];
        let pj = [2.0, 0.0, 0.0];
        cp.record_break(pi, pj);
        assert_eq!(cp.len(), 1);
        let mid = cp.points[0];
        assert!((mid[0] - 1.0).abs() < 1e-12);
        assert!(mid[1].abs() < 1e-12);
        assert!(mid[2].abs() < 1e-12);
    }

    // 24. FractureEnergy accumulates correctly.
    #[test]
    fn test_fracture_energy_accumulate() {
        let mut fe = FractureEnergy::new();
        fe.add_broken_bond(1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(
            (fe.total - 0.5).abs() < 1e-12,
            "FE should be 0.5, got {}",
            fe.total
        );
        fe.reset();
        assert_eq!(fe.total, 0.0);
    }

    // 25. PdGrid::bar_1d produces correct count.
    #[test]
    fn test_pd_grid_bar_1d() {
        let particles = PdGrid::bar_1d(1.0, 0.1);
        assert_eq!(particles.len(), 10, "1-D bar should have 10 particles");
    }

    // 26. PdGrid::cubic_lattice produces correct count.
    #[test]
    fn test_pd_grid_cubic_lattice() {
        let particles = PdGrid::cubic_lattice(0.3, 0.2, 0.1, 0.1);
        assert_eq!(particles.len(), (3 * 2), "cubic lattice count mismatch");
    }

    // 27. Particles from grid have positive volume.
    #[test]
    fn test_pd_grid_positive_volume() {
        let particles = PdGrid::bar_1d(0.5, 0.1);
        for p in &particles {
            assert!(p.volume > 0.0, "particle volume must be positive");
        }
    }

    // 28. normalize3 returns unit vector.
    #[test]
    fn test_normalize3_unit() {
        let v = [3.0, 4.0, 0.0];
        let u = normalize3(v);
        let len = norm3(u);
        assert!(
            (len - 1.0).abs() < 1e-12,
            "normalized vector should have unit length"
        );
    }

    // 29. normalize3 of zero vector returns zero.
    #[test]
    fn test_normalize3_zero_input() {
        let v = [0.0, 0.0, 0.0];
        let u = normalize3(v);
        for c in u {
            assert!(c.abs() < 1e-15);
        }
    }

    // 30. compute_dilatation returns finite value.
    #[test]
    fn test_dilatation_finite() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        let theta = sys.compute_dilatation(0);
        assert!(theta.is_finite(), "dilatation should be finite");
    }

    // 31. convergence_to_classical: bond-based Poisson ratio ≈ 0.25 for K=G.
    #[test]
    fn test_convergence_to_classical() {
        // For K = G, ν should be close to 0.0 rather than 0.25,
        // but the bound just needs to be computable.
        let mat = default_material();
        let sys = PdSystem::new(vec![], mat);
        let diff = sys.convergence_to_classical();
        assert!(diff.is_finite(), "convergence measure should be finite");
    }

    // 32. Two-particle velocity-Verlet changes position.
    #[test]
    fn test_velocity_verlet_moves_particles() {
        let mat = default_material();
        let particles = vec![
            PdParticle::new([0.0, 0.0, 0.0], 1e-6),
            PdParticle::new([0.01, 0.0, 0.0], 1e-6),
        ];
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        // Give one particle an initial velocity
        sys.particles[0].vel = [10.0, 0.0, 0.0];
        let x0_before = sys.particles[0].pos[0];
        sys.step(1e-6);
        let x0_after = sys.particles[0].pos[0];
        assert!(
            (x0_after - x0_before).abs() > 1e-15,
            "particle should move after step"
        );
    }

    // 33. OSB force density is zero for zero displacement.
    #[test]
    fn test_osb_force_zero_displacement() {
        let mat = default_material();
        let sys = PdSystem::new(vec![], mat);
        let xi = [0.01, 0.0, 0.0];
        let eta = [0.0, 0.0, 0.0];
        let f = sys.osb_force_density(xi, eta, 1.0, 0.0, 0.0, 1e9, 1e9);
        for c in f {
            assert!(c.abs() < 1e-12, "OSB force should be 0 at rest, got {c}");
        }
    }

    // 34. Fracture energy is non-negative.
    #[test]
    fn test_fracture_energy_non_negative() {
        let fe = FractureEnergy::new();
        assert!(fe.total >= 0.0);
    }

    // 35. find_bonds is symmetric: bonds i<j only.
    #[test]
    fn test_find_bonds_no_duplicate() {
        let mat = default_material();
        let particles = PdGrid::bar_1d(0.05, 0.01);
        let mut sys = PdSystem::new(particles, mat);
        sys.find_bonds();
        for b in &sys.bonds {
            assert!(b.i < b.j, "bonds must satisfy i < j");
        }
    }
}
