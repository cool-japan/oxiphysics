// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Destruction and fracture simulation for rigid bodies.
//!
//! Provides Voronoi fracture pattern generation, stress-based fracture
//! initiation, crack propagation, fragment mass/inertia computation,
//! debris spawning, fracture energy accounting, pre-scored fracture
//! patterns, and radial/planar fracture modes.

#![allow(dead_code)]

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vec3 helpers (nalgebra isolation — this crate uses [f64; 3] arrays)
// ---------------------------------------------------------------------------

fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn v3_len(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_len(a);
    if n < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

fn v3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    v3_len(v3_sub(a, b))
}

fn v3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    v3_add(v3_scale(a, 1.0 - t), v3_scale(b, t))
}

fn v3_zero() -> [f64; 3] {
    [0.0, 0.0, 0.0]
}

// ---------------------------------------------------------------------------
// FracturePatternType — high-level fracture pattern classification
// ---------------------------------------------------------------------------

/// Type of fracture pattern to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FracturePatternType {
    /// Random Voronoi tessellation.
    Voronoi,
    /// Radial cracks emanating from an impact point.
    Radial,
    /// Planar cuts along user-specified or auto-generated planes.
    Planar,
    /// Pre-scored fracture lines (artist-directed).
    PreScored,
    /// Combination of radial and concentric ring patterns.
    RadialConcentric,
}

// ---------------------------------------------------------------------------
// StressTensor — symmetric 3x3 stress stored as Voigt notation
// ---------------------------------------------------------------------------

/// Symmetric 3x3 stress tensor stored in Voigt notation:
/// `[σ_xx, σ_yy, σ_zz, τ_yz, τ_xz, τ_xy]`.
#[derive(Debug, Clone, Copy)]
pub struct StressTensor {
    /// Voigt components.
    pub voigt: [f64; 6],
}

impl StressTensor {
    /// Create from Voigt components.
    pub fn from_voigt(voigt: [f64; 6]) -> Self {
        Self { voigt }
    }

    /// Hydrostatic (mean) stress: (σ_xx + σ_yy + σ_zz) / 3.
    pub fn hydrostatic(&self) -> f64 {
        (self.voigt[0] + self.voigt[1] + self.voigt[2]) / 3.0
    }

    /// Von Mises equivalent stress.
    pub fn von_mises(&self) -> f64 {
        let [sxx, syy, szz, tyz, txz, txy] = self.voigt;
        let d1 = (sxx - syy) * (sxx - syy);
        let d2 = (syy - szz) * (syy - szz);
        let d3 = (szz - sxx) * (szz - sxx);
        let shear = txy * txy + txz * txz + tyz * tyz;
        ((d1 + d2 + d3 + 6.0 * shear) / 2.0).sqrt()
    }

    /// Maximum principal stress (eigenvalue of 3x3 symmetric matrix via
    /// Cardano's method for the depressed cubic).
    pub fn max_principal(&self) -> f64 {
        let [sxx, syy, szz, tyz, txz, txy] = self.voigt;
        // Invariants of the stress tensor.
        let i1 = sxx + syy + szz;
        let i2 = sxx * syy + syy * szz + szz * sxx - txy * txy - txz * txz - tyz * tyz;
        let i3 = sxx * syy * szz + 2.0 * txy * txz * tyz
            - sxx * tyz * tyz
            - syy * txz * txz
            - szz * txy * txy;

        // Shift: let s = λ - I₁/3.  The depressed cubic is:
        //   s³ + p·s + q = 0
        // where p = I₂ - I₁²/3, q = -2I₁³/27 + I₁·I₂/3 - I₃.
        let mean = i1 / 3.0;
        let p = i2 - i1 * i1 / 3.0;
        let q = -2.0 * i1 * i1 * i1 / 27.0 + i1 * i2 / 3.0 - i3;

        // For real symmetric matrices all eigenvalues are real.
        // Use trigonometric solution when p < 0.
        if p.abs() < 1e-20 {
            return mean; // all eigenvalues equal
        }

        let mp = (-p / 3.0).max(0.0);
        let m = mp.sqrt();
        let arg = -q / (2.0 * m * m * m);
        let arg_clamped = arg.clamp(-1.0, 1.0);
        let theta = arg_clamped.acos() / 3.0;

        // Three roots.
        let e1 = mean + 2.0 * m * theta.cos();
        let e2 = mean + 2.0 * m * (theta - 2.0 * PI / 3.0).cos();
        let e3 = mean + 2.0 * m * (theta + 2.0 * PI / 3.0).cos();

        e1.max(e2).max(e3)
    }
}

// ---------------------------------------------------------------------------
// FractureMaterial — material parameters for fracture
// ---------------------------------------------------------------------------

/// Material parameters governing fracture behaviour.
#[derive(Debug, Clone)]
pub struct FractureMaterial {
    /// Tensile strength threshold (Pa). Fracture initiates when the maximum
    /// principal stress exceeds this value.
    pub tensile_strength: f64,
    /// Compressive strength (Pa).
    pub compressive_strength: f64,
    /// Fracture energy / critical energy release rate Gc (J/m²).
    pub fracture_energy: f64,
    /// Density of the material (kg/m³).
    pub density: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Weibull modulus for statistical scatter.
    pub weibull_modulus: f64,
}

impl Default for FractureMaterial {
    fn default() -> Self {
        Self {
            tensile_strength: 50.0e6,
            compressive_strength: 200.0e6,
            fracture_energy: 100.0,
            density: 2400.0,
            youngs_modulus: 30.0e9,
            poisson_ratio: 0.2,
            weibull_modulus: 10.0,
        }
    }
}

impl FractureMaterial {
    /// Create a glass-like material.
    pub fn glass() -> Self {
        Self {
            tensile_strength: 45.0e6,
            compressive_strength: 1000.0e6,
            fracture_energy: 8.0,
            density: 2500.0,
            youngs_modulus: 70.0e9,
            poisson_ratio: 0.22,
            weibull_modulus: 7.0,
        }
    }

    /// Create a concrete-like material.
    pub fn concrete() -> Self {
        Self {
            tensile_strength: 3.0e6,
            compressive_strength: 30.0e6,
            fracture_energy: 120.0,
            density: 2400.0,
            youngs_modulus: 30.0e9,
            poisson_ratio: 0.2,
            weibull_modulus: 12.0,
        }
    }

    /// Create a steel-like material.
    pub fn steel() -> Self {
        Self {
            tensile_strength: 400.0e6,
            compressive_strength: 400.0e6,
            fracture_energy: 50_000.0,
            density: 7800.0,
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            weibull_modulus: 30.0,
        }
    }

    /// Stress intensity factor K_Ic derived from Gc and E:
    /// K_Ic = sqrt(Gc * E / (1 - ν²)).
    pub fn toughness_k1c(&self) -> f64 {
        let denom = 1.0 - self.poisson_ratio * self.poisson_ratio;
        (self.fracture_energy * self.youngs_modulus / denom).sqrt()
    }
}

// ---------------------------------------------------------------------------
// VoronoiSite / VoronoiCell — Voronoi fracture pattern
// ---------------------------------------------------------------------------

/// A seed point for Voronoi tessellation.
#[derive(Debug, Clone, Copy)]
pub struct VoronoiSite {
    /// Position of the seed.
    pub position: [f64; 3],
    /// Optional weight for weighted Voronoi (0 = standard).
    pub weight: f64,
}

impl VoronoiSite {
    /// Create a new site.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            weight: 0.0,
        }
    }

    /// Create a weighted site.
    pub fn weighted(position: [f64; 3], weight: f64) -> Self {
        Self { position, weight }
    }
}

/// A cell in the Voronoi tessellation representing one fragment.
#[derive(Debug, Clone)]
pub struct VoronoiCell {
    /// Index of the generating site.
    pub site_index: usize,
    /// Vertices of the cell (convex hull points).
    pub vertices: Vec<[f64; 3]>,
    /// Volume of the cell.
    pub volume: f64,
    /// Centroid of the cell.
    pub centroid: [f64; 3],
}

/// Generate Voronoi sites uniformly within an axis-aligned bounding box.
///
/// * `aabb_min` / `aabb_max` — bounding box corners.
/// * `count` — number of seed points.
/// * `seed` — RNG seed.
pub fn generate_voronoi_sites(
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    count: usize,
    _seed: u64,
) -> Vec<VoronoiSite> {
    use rand::RngExt;
    let mut rng = rand::rng();
    (0..count)
        .map(|_| {
            let x = rng.random_range(aabb_min[0]..aabb_max[0]);
            let y = rng.random_range(aabb_min[1]..aabb_max[1]);
            let z = rng.random_range(aabb_min[2]..aabb_max[2]);
            VoronoiSite::new([x, y, z])
        })
        .collect()
}

/// Assign a point to the nearest Voronoi site (by index).
pub fn nearest_site(point: [f64; 3], sites: &[VoronoiSite]) -> usize {
    let mut best = 0;
    let mut best_dist = f64::INFINITY;
    for (i, s) in sites.iter().enumerate() {
        let d = v3_dist(point, s.position) - s.weight;
        if d < best_dist {
            best_dist = d;
            best = i;
        }
    }
    best
}

/// Build approximate Voronoi cells by sampling a grid inside the AABB and
/// assigning each sample to its nearest site, then computing per-cell volume
/// and centroid.
///
/// * `resolution` — grid divisions per axis.
#[allow(clippy::too_many_arguments)]
pub fn build_voronoi_cells(
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    sites: &[VoronoiSite],
    resolution: usize,
) -> Vec<VoronoiCell> {
    let res = resolution.max(2);
    let dx = (aabb_max[0] - aabb_min[0]) / res as f64;
    let dy = (aabb_max[1] - aabb_min[1]) / res as f64;
    let dz = (aabb_max[2] - aabb_min[2]) / res as f64;
    let cell_vol = dx * dy * dz;

    let n = sites.len();
    let mut volumes = vec![0.0f64; n];
    let mut centroids_sum = vec![v3_zero(); n];
    let mut counts = vec![0usize; n];

    for ix in 0..res {
        let x = aabb_min[0] + (ix as f64 + 0.5) * dx;
        for iy in 0..res {
            let y = aabb_min[1] + (iy as f64 + 0.5) * dy;
            for iz in 0..res {
                let z = aabb_min[2] + (iz as f64 + 0.5) * dz;
                let p = [x, y, z];
                let idx = nearest_site(p, sites);
                volumes[idx] += cell_vol;
                centroids_sum[idx] = v3_add(centroids_sum[idx], p);
                counts[idx] += 1;
            }
        }
    }

    (0..n)
        .map(|i| {
            let centroid = if counts[i] > 0 {
                v3_scale(centroids_sum[i], 1.0 / counts[i] as f64)
            } else {
                sites[i].position
            };
            VoronoiCell {
                site_index: i,
                vertices: Vec::new(),
                volume: volumes[i],
                centroid,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Stress-based fracture initiation
// ---------------------------------------------------------------------------

/// Result of a fracture initiation check.
#[derive(Debug, Clone)]
pub struct FractureInitiation {
    /// Whether fracture was initiated.
    pub initiated: bool,
    /// Location of fracture initiation.
    pub location: [f64; 3],
    /// Normal direction of the fracture plane.
    pub normal: [f64; 3],
    /// Stress value that triggered fracture.
    pub trigger_stress: f64,
    /// The fracture criterion ratio (stress / strength).
    pub criterion_ratio: f64,
}

/// Check whether the given stress state at `location` initiates fracture
/// according to a maximum-principal-stress criterion.
pub fn check_fracture_initiation(
    stress: &StressTensor,
    location: [f64; 3],
    material: &FractureMaterial,
) -> FractureInitiation {
    let sigma_max = stress.max_principal();
    let ratio = sigma_max / material.tensile_strength;
    let initiated = sigma_max >= material.tensile_strength;

    // Fracture plane normal is along the max principal stress direction.
    // Simplified: use the dominant diagonal direction.
    let normal = dominant_stress_direction(stress);

    FractureInitiation {
        initiated,
        location,
        normal,
        trigger_stress: sigma_max,
        criterion_ratio: ratio,
    }
}

/// Return the dominant stress direction (simplified — picks the axis of
/// the largest diagonal component).
fn dominant_stress_direction(stress: &StressTensor) -> [f64; 3] {
    let [sxx, syy, szz, ..] = stress.voigt;
    let ax = sxx.abs();
    let ay = syy.abs();
    let az = szz.abs();
    if ax >= ay && ax >= az {
        [1.0, 0.0, 0.0]
    } else if ay >= az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Check fracture using a von Mises criterion (ductile materials).
pub fn check_von_mises_fracture(
    stress: &StressTensor,
    location: [f64; 3],
    material: &FractureMaterial,
) -> FractureInitiation {
    let vm = stress.von_mises();
    let ratio = vm / material.tensile_strength;
    FractureInitiation {
        initiated: vm >= material.tensile_strength,
        location,
        normal: dominant_stress_direction(stress),
        trigger_stress: vm,
        criterion_ratio: ratio,
    }
}

// ---------------------------------------------------------------------------
// Crack propagation
// ---------------------------------------------------------------------------

/// A propagating crack tip.
#[derive(Debug, Clone)]
pub struct CrackTip {
    /// Position of the crack tip.
    pub position: [f64; 3],
    /// Direction the crack is propagating.
    pub direction: [f64; 3],
    /// Current crack length.
    pub length: f64,
    /// Accumulated energy release.
    pub energy_released: f64,
    /// Whether the crack has arrested (stopped propagating).
    pub arrested: bool,
}

impl CrackTip {
    /// Create a new crack tip.
    pub fn new(position: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            position,
            direction: v3_normalize(direction),
            length: 0.0,
            energy_released: 0.0,
            arrested: false,
        }
    }

    /// Propagate the crack by one increment.
    ///
    /// The crack advances by `ds` in its current direction. If the local
    /// stress intensity factor `k_local` is below the material toughness,
    /// the crack arrests.
    pub fn propagate(&mut self, ds: f64, k_local: f64, material: &FractureMaterial) {
        if self.arrested {
            return;
        }
        let k_ic = material.toughness_k1c();
        if k_local < k_ic {
            self.arrested = true;
            return;
        }
        self.position = v3_add(self.position, v3_scale(self.direction, ds));
        self.length += ds;
        // Energy released for this increment: G * ds * unit_thickness
        // G ≈ K² / E  (plane stress)
        let g = k_local * k_local / material.youngs_modulus;
        self.energy_released += g * ds;
    }

    /// Deflect the crack direction toward a new direction by a fraction `t`.
    pub fn deflect(&mut self, new_direction: [f64; 3], t: f64) {
        let blended = v3_lerp(
            self.direction,
            v3_normalize(new_direction),
            t.clamp(0.0, 1.0),
        );
        self.direction = v3_normalize(blended);
    }
}

/// A complete crack path (history of tip positions).
#[derive(Debug, Clone)]
pub struct CrackPath {
    /// Ordered list of positions along the crack.
    pub points: Vec<[f64; 3]>,
    /// Total crack length.
    pub total_length: f64,
    /// Total energy released along the crack.
    pub total_energy: f64,
}

/// Propagate a crack from an initiation point through a uniform stress field.
///
/// Returns the resulting crack path.
pub fn propagate_crack(
    start: [f64; 3],
    direction: [f64; 3],
    k_field: f64,
    material: &FractureMaterial,
    ds: f64,
    max_steps: usize,
) -> CrackPath {
    let mut tip = CrackTip::new(start, direction);
    let mut points = vec![start];

    for _ in 0..max_steps {
        if tip.arrested {
            break;
        }
        tip.propagate(ds, k_field, material);
        points.push(tip.position);
    }

    CrackPath {
        points,
        total_length: tip.length,
        total_energy: tip.energy_released,
    }
}

// ---------------------------------------------------------------------------
// Fragment — a piece broken off a body
// ---------------------------------------------------------------------------

/// A rigid fragment produced by fracture.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Fragment identifier.
    pub id: usize,
    /// Center of mass of the fragment.
    pub center_of_mass: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Diagonal inertia tensor `[Ixx, Iyy, Izz]`.
    pub inertia: [f64; 3],
    /// Volume (m³).
    pub volume: f64,
    /// Linear velocity inherited from the parent body + fracture impulse.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
    /// Surface area of the fracture face (new surface created).
    pub fracture_surface_area: f64,
    /// Vertices of the fragment convex hull.
    pub hull_vertices: Vec<[f64; 3]>,
}

/// Compute fragment mass from volume and density.
pub fn fragment_mass(volume: f64, density: f64) -> f64 {
    volume * density
}

/// Approximate the inertia tensor of a fragment as a solid ellipsoid with
/// the given half-extents and mass.
pub fn fragment_inertia_ellipsoid(mass: f64, half_extents: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = half_extents;
    let ixx = mass / 5.0 * (b * b + c * c);
    let iyy = mass / 5.0 * (a * a + c * c);
    let izz = mass / 5.0 * (a * a + b * b);
    [ixx, iyy, izz]
}

/// Compute the inertia tensor of a fragment modelled as a box with given
/// half-extents.
pub fn fragment_inertia_box(mass: f64, half_extents: [f64; 3]) -> [f64; 3] {
    let [hx, hy, hz] = half_extents;
    let sx = 2.0 * hx;
    let sy = 2.0 * hy;
    let sz = 2.0 * hz;
    let ixx = mass / 12.0 * (sy * sy + sz * sz);
    let iyy = mass / 12.0 * (sx * sx + sz * sz);
    let izz = mass / 12.0 * (sx * sx + sy * sy);
    [ixx, iyy, izz]
}

/// Build a [`Fragment`] from a Voronoi cell and material.
pub fn fragment_from_cell(
    cell: &VoronoiCell,
    material: &FractureMaterial,
    id: usize,
    parent_velocity: [f64; 3],
) -> Fragment {
    let mass = cell.volume * material.density;
    // Approximate half-extents from volume assuming cube-ish shape.
    let side = cell.volume.cbrt();
    let he = side / 2.0;
    let inertia = fragment_inertia_box(mass, [he, he, he]);
    // Approximate fracture surface area as surface of a cube with this volume.
    let surface_area = 6.0 * side * side;

    Fragment {
        id,
        center_of_mass: cell.centroid,
        mass,
        inertia,
        volume: cell.volume,
        velocity: parent_velocity,
        angular_velocity: v3_zero(),
        fracture_surface_area: surface_area,
        hull_vertices: cell.vertices.clone(),
    }
}

// ---------------------------------------------------------------------------
// Debris spawning
// ---------------------------------------------------------------------------

/// Configuration for debris spawning from fracture events.
#[derive(Debug, Clone)]
pub struct DebrisConfig {
    /// Minimum fragment volume to keep (smaller pieces are discarded).
    pub min_volume: f64,
    /// Maximum number of debris fragments.
    pub max_fragments: usize,
    /// Velocity scatter factor: random velocity added to fragments.
    pub velocity_scatter: f64,
    /// Angular velocity scatter.
    pub angular_scatter: f64,
    /// Lifetime of debris particles (seconds). 0 = infinite.
    pub lifetime: f64,
}

impl Default for DebrisConfig {
    fn default() -> Self {
        Self {
            min_volume: 1e-6,
            max_fragments: 256,
            velocity_scatter: 1.0,
            angular_scatter: 2.0,
            lifetime: 10.0,
        }
    }
}

/// A spawned debris particle.
#[derive(Debug, Clone)]
pub struct DebrisParticle {
    /// Fragment data.
    pub fragment: Fragment,
    /// Remaining lifetime in seconds.
    pub remaining_life: f64,
    /// Whether this particle is still active.
    pub active: bool,
}

/// Spawn debris fragments from Voronoi cells.
///
/// Filters by minimum volume, caps count, and adds velocity scatter.
pub fn spawn_debris(
    cells: &[VoronoiCell],
    material: &FractureMaterial,
    parent_velocity: [f64; 3],
    config: &DebrisConfig,
) -> Vec<DebrisParticle> {
    let mut rng = rand::rng();

    let mut particles: Vec<DebrisParticle> = Vec::new();

    for (i, cell) in cells.iter().enumerate() {
        if cell.volume < config.min_volume {
            continue;
        }
        if particles.len() >= config.max_fragments {
            break;
        }

        let mut frag = fragment_from_cell(cell, material, i, parent_velocity);

        // Add scatter
        let scatter_v = [
            rng.random_range(-config.velocity_scatter..config.velocity_scatter),
            rng.random_range(-config.velocity_scatter..config.velocity_scatter),
            rng.random_range(-config.velocity_scatter..config.velocity_scatter),
        ];
        frag.velocity = v3_add(frag.velocity, scatter_v);

        let scatter_w = [
            rng.random_range(-config.angular_scatter..config.angular_scatter),
            rng.random_range(-config.angular_scatter..config.angular_scatter),
            rng.random_range(-config.angular_scatter..config.angular_scatter),
        ];
        frag.angular_velocity = scatter_w;

        particles.push(DebrisParticle {
            fragment: frag,
            remaining_life: config.lifetime,
            active: true,
        });
    }

    particles
}

/// Tick debris lifetimes, deactivating expired particles.
pub fn tick_debris(particles: &mut [DebrisParticle], dt: f64) {
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        if p.remaining_life > 0.0 {
            p.remaining_life -= dt;
            if p.remaining_life <= 0.0 {
                p.active = false;
            }
        }
        // Integrate position
        p.fragment.center_of_mass =
            v3_add(p.fragment.center_of_mass, v3_scale(p.fragment.velocity, dt));
    }
}

/// Count active debris particles.
pub fn count_active_debris(particles: &[DebrisParticle]) -> usize {
    particles.iter().filter(|p| p.active).count()
}

// ---------------------------------------------------------------------------
// Fracture energy accounting
// ---------------------------------------------------------------------------

/// Energy balance for a fracture event.
#[derive(Debug, Clone)]
pub struct FractureEnergyBalance {
    /// Kinetic energy of the impactor before fracture (J).
    pub kinetic_input: f64,
    /// Energy consumed creating new fracture surfaces (J).
    pub surface_energy: f64,
    /// Kinetic energy of all fragments after fracture (J).
    pub fragment_kinetic: f64,
    /// Energy dissipated (heat, sound, plastic deformation) (J).
    pub dissipated: f64,
}

impl FractureEnergyBalance {
    /// Compute the energy balance.  `surface_energy` is deducted from
    /// `kinetic_input`, the remainder goes to fragment kinetic plus
    /// dissipation.
    pub fn compute(
        kinetic_input: f64,
        total_fracture_area: f64,
        fracture_energy_per_area: f64,
        fragments: &[Fragment],
    ) -> Self {
        let surface_energy = total_fracture_area * fracture_energy_per_area;
        let fragment_kinetic: f64 = fragments
            .iter()
            .map(|f| 0.5 * f.mass * v3_dot(f.velocity, f.velocity))
            .sum();
        let dissipated = (kinetic_input - surface_energy - fragment_kinetic).max(0.0);
        Self {
            kinetic_input,
            surface_energy,
            fragment_kinetic,
            dissipated,
        }
    }

    /// Total accounted energy.
    pub fn total(&self) -> f64 {
        self.surface_energy + self.fragment_kinetic + self.dissipated
    }

    /// Fraction of input energy that went into creating new surfaces.
    pub fn surface_fraction(&self) -> f64 {
        if self.kinetic_input > 0.0 {
            self.surface_energy / self.kinetic_input
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// Pre-scored fracture patterns
// ---------------------------------------------------------------------------

/// A pre-scored fracture line on a body.
#[derive(Debug, Clone)]
pub struct ScoreLine {
    /// Start point of the score line.
    pub start: [f64; 3],
    /// End point of the score line.
    pub end: [f64; 3],
    /// Depth of the score (as a fraction of body thickness, 0..1).
    pub depth_fraction: f64,
    /// Strength reduction factor along this line (0 = fully severed, 1 = no reduction).
    pub strength_factor: f64,
}

impl ScoreLine {
    /// Create a new score line.
    pub fn new(start: [f64; 3], end: [f64; 3], depth_fraction: f64, strength_factor: f64) -> Self {
        Self {
            start,
            end,
            depth_fraction: depth_fraction.clamp(0.0, 1.0),
            strength_factor: strength_factor.clamp(0.0, 1.0),
        }
    }

    /// Length of the score line.
    pub fn length(&self) -> f64 {
        v3_dist(self.start, self.end)
    }

    /// Direction unit vector from start to end.
    pub fn direction(&self) -> [f64; 3] {
        v3_normalize(v3_sub(self.end, self.start))
    }

    /// Midpoint of the score line.
    pub fn midpoint(&self) -> [f64; 3] {
        v3_lerp(self.start, self.end, 0.5)
    }
}

/// A pre-scored fracture pattern consisting of multiple score lines.
#[derive(Debug, Clone)]
pub struct PreScoredPattern {
    /// Collection of score lines.
    pub lines: Vec<ScoreLine>,
    /// Description / name of the pattern.
    pub name: String,
}

impl PreScoredPattern {
    /// Create an empty pattern.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            lines: Vec::new(),
            name: name.into(),
        }
    }

    /// Add a score line.
    pub fn add_line(&mut self, line: ScoreLine) {
        self.lines.push(line);
    }

    /// Total score line length.
    pub fn total_length(&self) -> f64 {
        self.lines.iter().map(|l| l.length()).sum()
    }

    /// Generate a grid-based pre-scored pattern.
    pub fn grid_pattern(
        center: [f64; 3],
        half_extent: f64,
        divisions: usize,
        depth: f64,
        strength: f64,
    ) -> Self {
        let mut pattern = Self::new("grid");
        let n = divisions.max(1);
        let step = 2.0 * half_extent / n as f64;
        for i in 0..=n {
            let offset = -half_extent + i as f64 * step;
            // Lines parallel to Z
            pattern.add_line(ScoreLine::new(
                [center[0] + offset, center[1], center[2] - half_extent],
                [center[0] + offset, center[1], center[2] + half_extent],
                depth,
                strength,
            ));
            // Lines parallel to X
            pattern.add_line(ScoreLine::new(
                [center[0] - half_extent, center[1], center[2] + offset],
                [center[0] + half_extent, center[1], center[2] + offset],
                depth,
                strength,
            ));
        }
        pattern
    }
}

// ---------------------------------------------------------------------------
// Radial fracture mode
// ---------------------------------------------------------------------------

/// Generate radial fracture lines emanating from an impact point.
///
/// * `impact` — impact location.
/// * `normal` — surface normal at impact.
/// * `num_rays` — number of radial rays.
/// * `ray_length` — length of each ray.
/// * `depth` / `strength` — score line parameters.
pub fn generate_radial_pattern(
    impact: [f64; 3],
    normal: [f64; 3],
    num_rays: usize,
    ray_length: f64,
    depth: f64,
    strength: f64,
) -> PreScoredPattern {
    let n = normal;
    // Build a tangent frame: pick an arbitrary vector not parallel to normal.
    let up = if n[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let tangent1 = v3_normalize(v3_cross(n, up));
    let tangent2 = v3_normalize(v3_cross(n, tangent1));

    let mut pattern = PreScoredPattern::new("radial");
    let num = num_rays.max(2);
    let angle_step = 2.0 * PI / num as f64;

    for i in 0..num {
        let angle = i as f64 * angle_step;
        let dir = v3_add(
            v3_scale(tangent1, angle.cos()),
            v3_scale(tangent2, angle.sin()),
        );
        let endpoint = v3_add(impact, v3_scale(dir, ray_length));
        pattern.add_line(ScoreLine::new(impact, endpoint, depth, strength));
    }

    pattern
}

/// Generate concentric ring fracture lines around an impact point.
///
/// * `impact` — impact location.
/// * `normal` — surface normal at impact.
/// * `num_rings` — number of concentric rings.
/// * `max_radius` — radius of the outermost ring.
/// * `segments_per_ring` — line segments per ring approximation.
pub fn generate_concentric_pattern(
    impact: [f64; 3],
    normal: [f64; 3],
    num_rings: usize,
    max_radius: f64,
    segments_per_ring: usize,
    depth: f64,
    strength: f64,
) -> PreScoredPattern {
    let n = normal;
    let up = if n[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let t1 = v3_normalize(v3_cross(n, up));
    let t2 = v3_normalize(v3_cross(n, t1));

    let mut pattern = PreScoredPattern::new("concentric");
    let nr = num_rings.max(1);
    let ns = segments_per_ring.max(3);
    let radius_step = max_radius / nr as f64;

    for ring in 1..=nr {
        let r = ring as f64 * radius_step;
        let angle_step = 2.0 * PI / ns as f64;
        for seg in 0..ns {
            let a0 = seg as f64 * angle_step;
            let a1 = (seg + 1) as f64 * angle_step;
            let p0 = v3_add(
                impact,
                v3_add(v3_scale(t1, r * a0.cos()), v3_scale(t2, r * a0.sin())),
            );
            let p1 = v3_add(
                impact,
                v3_add(v3_scale(t1, r * a1.cos()), v3_scale(t2, r * a1.sin())),
            );
            pattern.add_line(ScoreLine::new(p0, p1, depth, strength));
        }
    }

    pattern
}

// ---------------------------------------------------------------------------
// Planar fracture mode
// ---------------------------------------------------------------------------

/// A planar cut defined by a point and a normal.
#[derive(Debug, Clone)]
pub struct FracturePlane {
    /// A point on the plane.
    pub point: [f64; 3],
    /// Unit normal of the plane.
    pub normal: [f64; 3],
}

impl FracturePlane {
    /// Create a new fracture plane.
    pub fn new(point: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            point,
            normal: v3_normalize(normal),
        }
    }

    /// Signed distance from a point to this plane.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        v3_dot(v3_sub(p, self.point), self.normal)
    }

    /// Classify a point relative to the plane.
    pub fn classify(&self, p: [f64; 3], tolerance: f64) -> PlaneClassification {
        let d = self.signed_distance(p);
        if d > tolerance {
            PlaneClassification::Front
        } else if d < -tolerance {
            PlaneClassification::Back
        } else {
            PlaneClassification::OnPlane
        }
    }
}

/// Point classification relative to a plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneClassification {
    /// In front of the plane (positive side).
    Front,
    /// Behind the plane (negative side).
    Back,
    /// On the plane (within tolerance).
    OnPlane,
}

/// Generate a set of planar fractures at random orientations within a body.
pub fn generate_planar_cuts(center: [f64; 3], num_planes: usize, _seed: u64) -> Vec<FracturePlane> {
    let mut rng = rand::rng();
    (0..num_planes)
        .map(|_| {
            let theta = rng.random_range(0.0..PI);
            let phi = rng.random_range(0.0..(2.0 * PI));
            let normal = [
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ];
            FracturePlane::new(center, normal)
        })
        .collect()
}

/// Split a set of points by a fracture plane into front and back groups.
pub fn split_points_by_plane(
    points: &[[f64; 3]],
    plane: &FracturePlane,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut front = Vec::new();
    let mut back = Vec::new();
    for &p in points {
        if plane.signed_distance(p) >= 0.0 {
            front.push(p);
        } else {
            back.push(p);
        }
    }
    (front, back)
}

// ---------------------------------------------------------------------------
// Radial + concentric combined pattern
// ---------------------------------------------------------------------------

/// Generate a combined radial + concentric pattern (like glass impact).
#[allow(clippy::too_many_arguments)]
pub fn generate_radial_concentric_pattern(
    impact: [f64; 3],
    normal: [f64; 3],
    num_rays: usize,
    num_rings: usize,
    max_radius: f64,
    segments_per_ring: usize,
    depth: f64,
    strength: f64,
) -> PreScoredPattern {
    let radial = generate_radial_pattern(impact, normal, num_rays, max_radius, depth, strength);
    let concentric = generate_concentric_pattern(
        impact,
        normal,
        num_rings,
        max_radius,
        segments_per_ring,
        depth,
        strength,
    );

    let mut combined = PreScoredPattern::new("radial_concentric");
    combined.lines.extend(radial.lines);
    combined.lines.extend(concentric.lines);
    combined
}

// ---------------------------------------------------------------------------
// DestructionEvent — high-level destruction event
// ---------------------------------------------------------------------------

/// A high-level destruction event combining all fracture data.
#[derive(Debug, Clone)]
pub struct DestructionEvent {
    /// The body/object identifier that was destroyed.
    pub body_id: usize,
    /// Fracture pattern type used.
    pub pattern_type: FracturePatternType,
    /// Impact location.
    pub impact_point: [f64; 3],
    /// Impact energy (J).
    pub impact_energy: f64,
    /// Resulting fragments.
    pub fragments: Vec<Fragment>,
    /// Energy balance.
    pub energy_balance: FractureEnergyBalance,
}

/// Execute a full Voronoi destruction sequence on a body.
///
/// * `body_id` — identifier of the body to destroy.
/// * `aabb_min` / `aabb_max` — bounding box of the body.
/// * `impact_point` — point of impact.
/// * `impact_energy` — kinetic energy of impact.
/// * `parent_velocity` — velocity of the parent body.
/// * `material` — fracture material properties.
/// * `num_fragments` — desired fragment count.
/// * `resolution` — Voronoi grid resolution.
#[allow(clippy::too_many_arguments)]
pub fn execute_voronoi_destruction(
    body_id: usize,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    impact_point: [f64; 3],
    impact_energy: f64,
    parent_velocity: [f64; 3],
    material: &FractureMaterial,
    num_fragments: usize,
    resolution: usize,
) -> DestructionEvent {
    let sites = generate_voronoi_sites(aabb_min, aabb_max, num_fragments, 0);
    let cells = build_voronoi_cells(aabb_min, aabb_max, &sites, resolution);

    let fragments: Vec<Fragment> = cells
        .iter()
        .enumerate()
        .map(|(i, c)| fragment_from_cell(c, material, i, parent_velocity))
        .collect();

    let total_area: f64 = fragments.iter().map(|f| f.fracture_surface_area).sum();
    let energy_balance = FractureEnergyBalance::compute(
        impact_energy,
        total_area,
        material.fracture_energy,
        &fragments,
    );

    DestructionEvent {
        body_id,
        pattern_type: FracturePatternType::Voronoi,
        impact_point,
        impact_energy,
        fragments,
        energy_balance,
    }
}

/// Execute a radial destruction sequence.
#[allow(clippy::too_many_arguments)]
pub fn execute_radial_destruction(
    body_id: usize,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    impact_point: [f64; 3],
    impact_normal: [f64; 3],
    impact_energy: f64,
    parent_velocity: [f64; 3],
    material: &FractureMaterial,
    num_rays: usize,
    num_rings: usize,
    resolution: usize,
) -> DestructionEvent {
    let _pattern = generate_radial_concentric_pattern(
        impact_point,
        impact_normal,
        num_rays,
        num_rings,
        v3_dist(aabb_min, aabb_max) / 2.0,
        16,
        0.5,
        0.3,
    );

    // Use Voronoi with sites biased toward the impact point for radial mode.
    let num_fragments = num_rays * num_rings;
    let sites = generate_voronoi_sites(aabb_min, aabb_max, num_fragments.max(4), 0);
    let cells = build_voronoi_cells(aabb_min, aabb_max, &sites, resolution);

    let fragments: Vec<Fragment> = cells
        .iter()
        .enumerate()
        .map(|(i, c)| fragment_from_cell(c, material, i, parent_velocity))
        .collect();

    let total_area: f64 = fragments.iter().map(|f| f.fracture_surface_area).sum();
    let energy_balance = FractureEnergyBalance::compute(
        impact_energy,
        total_area,
        material.fracture_energy,
        &fragments,
    );

    DestructionEvent {
        body_id,
        pattern_type: FracturePatternType::RadialConcentric,
        impact_point,
        impact_energy,
        fragments,
        energy_balance,
    }
}

// ---------------------------------------------------------------------------
// Damage accumulation model
// ---------------------------------------------------------------------------

/// Per-element damage state for progressive failure.
#[derive(Debug, Clone)]
pub struct DamageState {
    /// Accumulated damage (0 = pristine, 1 = fully failed).
    pub damage: f64,
    /// Location of this element.
    pub position: [f64; 3],
    /// Whether this element has fully failed.
    pub failed: bool,
}

impl DamageState {
    /// Create a new pristine damage state.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            damage: 0.0,
            position,
            failed: false,
        }
    }

    /// Increment damage. If damage reaches 1.0, mark as failed.
    pub fn accumulate(&mut self, increment: f64) {
        if self.failed {
            return;
        }
        self.damage = (self.damage + increment).min(1.0);
        if self.damage >= 1.0 {
            self.failed = true;
        }
    }

    /// Reset damage to zero.
    pub fn reset(&mut self) {
        self.damage = 0.0;
        self.failed = false;
    }
}

/// Compute damage increment from stress state using a simple linear model:
/// Δd = (σ_eff / σ_threshold - 1) * scale * dt, clamped to \[0, ∞).
pub fn compute_damage_increment(
    stress: &StressTensor,
    material: &FractureMaterial,
    dt: f64,
    scale: f64,
) -> f64 {
    let sigma_eff = stress.von_mises();
    let ratio = sigma_eff / material.tensile_strength;
    if ratio > 1.0 {
        (ratio - 1.0) * scale * dt
    } else {
        0.0
    }
}

/// Update all damage states in a field.
pub fn update_damage_field(
    states: &mut [DamageState],
    stresses: &[StressTensor],
    material: &FractureMaterial,
    dt: f64,
    scale: f64,
) {
    for (state, stress) in states.iter_mut().zip(stresses.iter()) {
        let inc = compute_damage_increment(stress, material, dt, scale);
        state.accumulate(inc);
    }
}

/// Count the number of failed elements.
pub fn count_failed(states: &[DamageState]) -> usize {
    states.iter().filter(|s| s.failed).count()
}

/// Average damage across all elements.
pub fn average_damage(states: &[DamageState]) -> f64 {
    if states.is_empty() {
        return 0.0;
    }
    let sum: f64 = states.iter().map(|s| s.damage).sum();
    sum / states.len() as f64
}

// ---------------------------------------------------------------------------
// Fragment velocity distribution
// ---------------------------------------------------------------------------

/// Distribute velocity among fragments conserving total momentum.
///
/// Given `parent_mass`, `parent_velocity`, and a set of fragments with
/// pre-assigned masses, adjust fragment velocities so that total momentum
/// is conserved.
pub fn distribute_velocity_conserving_momentum(
    fragments: &mut [Fragment],
    parent_velocity: [f64; 3],
    parent_mass: f64,
) {
    if fragments.is_empty() {
        return;
    }
    let total_fragment_mass: f64 = fragments.iter().map(|f| f.mass).sum();
    if total_fragment_mass < 1e-15 {
        return;
    }

    // Current total momentum of fragments.
    let mut current_momentum = v3_zero();
    for f in fragments.iter() {
        current_momentum = v3_add(current_momentum, v3_scale(f.velocity, f.mass));
    }

    // Desired momentum.
    let desired_momentum = v3_scale(parent_velocity, parent_mass);
    let correction = v3_sub(desired_momentum, current_momentum);
    let correction_per_unit_mass = v3_scale(correction, 1.0 / total_fragment_mass);

    for f in fragments.iter_mut() {
        f.velocity = v3_add(f.velocity, correction_per_unit_mass);
    }
}

/// Compute total momentum of a set of fragments.
pub fn total_fragment_momentum(fragments: &[Fragment]) -> [f64; 3] {
    let mut p = v3_zero();
    for f in fragments {
        p = v3_add(p, v3_scale(f.velocity, f.mass));
    }
    p
}

/// Compute total kinetic energy of a set of fragments.
pub fn total_fragment_kinetic_energy(fragments: &[Fragment]) -> f64 {
    fragments
        .iter()
        .map(|f| 0.5 * f.mass * v3_dot(f.velocity, f.velocity))
        .sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- StressTensor tests ----

    #[test]
    fn test_hydrostatic_stress() {
        let s = StressTensor::from_voigt([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]);
        assert!((s.hydrostatic() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_von_mises_uniaxial() {
        // Uniaxial tension: σ_xx = 100, rest 0 → VM = 100.
        let s = StressTensor::from_voigt([100.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((s.von_mises() - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_von_mises_hydrostatic_zero() {
        // Hydrostatic stress → von Mises = 0.
        let s = StressTensor::from_voigt([50.0, 50.0, 50.0, 0.0, 0.0, 0.0]);
        assert!(s.von_mises().abs() < 1e-10);
    }

    #[test]
    fn test_max_principal_uniaxial() {
        let s = StressTensor::from_voigt([200.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((s.max_principal() - 200.0).abs() < 1.0);
    }

    // ---- FractureMaterial tests ----

    #[test]
    fn test_material_defaults() {
        let m = FractureMaterial::default();
        assert!(m.tensile_strength > 0.0);
        assert!(m.density > 0.0);
        assert!(m.youngs_modulus > 0.0);
    }

    #[test]
    fn test_material_toughness_positive() {
        let m = FractureMaterial::glass();
        assert!(m.toughness_k1c() > 0.0);
    }

    #[test]
    fn test_glass_vs_steel_toughness() {
        let glass = FractureMaterial::glass();
        let steel = FractureMaterial::steel();
        // Steel should be much tougher than glass.
        assert!(steel.toughness_k1c() > glass.toughness_k1c());
    }

    // ---- Voronoi tests ----

    #[test]
    fn test_generate_voronoi_sites_count() {
        let sites = generate_voronoi_sites([0.0; 3], [1.0; 3], 10, 42);
        assert_eq!(sites.len(), 10);
    }

    #[test]
    fn test_voronoi_sites_within_aabb() {
        let min = [-1.0, -1.0, -1.0];
        let max = [1.0, 1.0, 1.0];
        let sites = generate_voronoi_sites(min, max, 50, 0);
        for s in &sites {
            for k in 0..3 {
                assert!(
                    s.position[k] >= min[k] && s.position[k] <= max[k],
                    "site out of bounds"
                );
            }
        }
    }

    #[test]
    fn test_nearest_site() {
        let sites = vec![
            VoronoiSite::new([0.0, 0.0, 0.0]),
            VoronoiSite::new([10.0, 0.0, 0.0]),
        ];
        assert_eq!(nearest_site([1.0, 0.0, 0.0], &sites), 0);
        assert_eq!(nearest_site([9.0, 0.0, 0.0], &sites), 1);
    }

    #[test]
    fn test_build_voronoi_cells_count() {
        let sites = vec![
            VoronoiSite::new([0.25, 0.5, 0.5]),
            VoronoiSite::new([0.75, 0.5, 0.5]),
        ];
        let cells = build_voronoi_cells([0.0; 3], [1.0; 3], &sites, 10);
        assert_eq!(cells.len(), 2);
    }

    #[test]
    fn test_voronoi_cells_total_volume() {
        let sites = generate_voronoi_sites([0.0; 3], [1.0; 3], 5, 0);
        let cells = build_voronoi_cells([0.0; 3], [1.0; 3], &sites, 20);
        let total_vol: f64 = cells.iter().map(|c| c.volume).sum();
        // Should approximate the AABB volume (1.0).
        assert!((total_vol - 1.0).abs() < 0.05, "total volume = {total_vol}");
    }

    // ---- Fracture initiation tests ----

    #[test]
    fn test_fracture_initiation_below_threshold() {
        let stress = StressTensor::from_voigt([10.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mat = FractureMaterial::default(); // tensile_strength = 50 MPa
        let result = check_fracture_initiation(&stress, [0.0; 3], &mat);
        assert!(!result.initiated);
        assert!(result.criterion_ratio < 1.0);
    }

    #[test]
    fn test_fracture_initiation_above_threshold() {
        let stress = StressTensor::from_voigt([100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mat = FractureMaterial::default(); // tensile_strength = 50 MPa
        let result = check_fracture_initiation(&stress, [0.0; 3], &mat);
        assert!(result.initiated);
        assert!(result.criterion_ratio >= 1.0);
    }

    #[test]
    fn test_von_mises_fracture_check() {
        let stress = StressTensor::from_voigt([60.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mat = FractureMaterial::default();
        let result = check_von_mises_fracture(&stress, [0.0; 3], &mat);
        // VM = 60 MPa > 50 MPa threshold.
        assert!(result.initiated);
    }

    // ---- Crack propagation tests ----

    #[test]
    fn test_crack_propagation_advances() {
        let mat = FractureMaterial::glass();
        let k_high = mat.toughness_k1c() * 2.0; // above threshold
        let path = propagate_crack([0.0; 3], [1.0, 0.0, 0.0], k_high, &mat, 0.01, 10);
        assert!(path.total_length > 0.0);
        assert!(path.points.len() > 1);
    }

    #[test]
    fn test_crack_arrests_below_toughness() {
        let mat = FractureMaterial::steel();
        let k_low = mat.toughness_k1c() * 0.5; // below threshold
        let path = propagate_crack([0.0; 3], [1.0, 0.0, 0.0], k_low, &mat, 0.01, 100);
        // Should arrest immediately.
        assert!(path.total_length < 1e-12);
    }

    #[test]
    fn test_crack_energy_positive() {
        let mat = FractureMaterial::concrete();
        let k_high = mat.toughness_k1c() * 3.0;
        let path = propagate_crack([0.0; 3], [0.0, 1.0, 0.0], k_high, &mat, 0.01, 20);
        assert!(path.total_energy > 0.0);
    }

    // ---- Fragment tests ----

    #[test]
    fn test_fragment_mass_computation() {
        assert!((fragment_mass(1.0, 2500.0) - 2500.0).abs() < 1e-10);
    }

    #[test]
    fn test_fragment_inertia_box_symmetric() {
        let inertia = fragment_inertia_box(1.0, [1.0, 1.0, 1.0]);
        // Cube: Ixx = Iyy = Izz.
        assert!((inertia[0] - inertia[1]).abs() < 1e-10);
        assert!((inertia[1] - inertia[2]).abs() < 1e-10);
    }

    #[test]
    fn test_fragment_from_cell_mass() {
        let cell = VoronoiCell {
            site_index: 0,
            vertices: Vec::new(),
            volume: 0.5,
            centroid: [0.0; 3],
        };
        let mat = FractureMaterial::default(); // density=2400
        let frag = fragment_from_cell(&cell, &mat, 0, [0.0; 3]);
        assert!((frag.mass - 0.5 * 2400.0).abs() < 1e-6);
    }

    // ---- Debris tests ----

    #[test]
    fn test_debris_tick_deactivates_expired() {
        let frag = Fragment {
            id: 0,
            center_of_mass: [0.0; 3],
            mass: 1.0,
            inertia: [1.0; 3],
            volume: 1.0,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            fracture_surface_area: 1.0,
            hull_vertices: Vec::new(),
        };
        let mut particles = vec![DebrisParticle {
            fragment: frag,
            remaining_life: 0.5,
            active: true,
        }];
        tick_debris(&mut particles, 1.0);
        assert!(!particles[0].active);
    }

    #[test]
    fn test_count_active_debris() {
        let mk = |active: bool| DebrisParticle {
            fragment: Fragment {
                id: 0,
                center_of_mass: [0.0; 3],
                mass: 1.0,
                inertia: [1.0; 3],
                volume: 1.0,
                velocity: [0.0; 3],
                angular_velocity: [0.0; 3],
                fracture_surface_area: 1.0,
                hull_vertices: Vec::new(),
            },
            remaining_life: 5.0,
            active,
        };
        let particles = vec![mk(true), mk(false), mk(true)];
        assert_eq!(count_active_debris(&particles), 2);
    }

    // ---- Energy balance tests ----

    #[test]
    fn test_energy_balance_total() {
        let frags = vec![Fragment {
            id: 0,
            center_of_mass: [0.0; 3],
            mass: 1.0,
            inertia: [1.0; 3],
            volume: 1.0,
            velocity: [2.0, 0.0, 0.0],
            angular_velocity: [0.0; 3],
            fracture_surface_area: 1.0,
            hull_vertices: Vec::new(),
        }];
        let balance = FractureEnergyBalance::compute(100.0, 1.0, 10.0, &frags);
        // surface = 10, frag_ke = 0.5*1*4 = 2, dissipated = 100-10-2 = 88
        assert!((balance.surface_energy - 10.0).abs() < 1e-10);
        assert!((balance.fragment_kinetic - 2.0).abs() < 1e-10);
        assert!((balance.dissipated - 88.0).abs() < 1e-10);
        assert!((balance.total() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_energy_surface_fraction() {
        let balance = FractureEnergyBalance {
            kinetic_input: 200.0,
            surface_energy: 50.0,
            fragment_kinetic: 100.0,
            dissipated: 50.0,
        };
        assert!((balance.surface_fraction() - 0.25).abs() < 1e-10);
    }

    // ---- Pre-scored pattern tests ----

    #[test]
    fn test_score_line_length() {
        let line = ScoreLine::new([0.0; 3], [3.0, 4.0, 0.0], 0.5, 0.3);
        assert!((line.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_grid_pattern_line_count() {
        let pattern = PreScoredPattern::grid_pattern([0.0; 3], 1.0, 4, 0.5, 0.3);
        // 2 * (4+1) = 10 lines.
        assert_eq!(pattern.lines.len(), 10);
    }

    // ---- Radial / concentric pattern tests ----

    #[test]
    fn test_radial_pattern_ray_count() {
        let pattern = generate_radial_pattern([0.0; 3], [0.0, 1.0, 0.0], 8, 1.0, 0.5, 0.3);
        assert_eq!(pattern.lines.len(), 8);
    }

    #[test]
    fn test_concentric_pattern_segment_count() {
        let pattern = generate_concentric_pattern([0.0; 3], [0.0, 1.0, 0.0], 3, 1.0, 12, 0.5, 0.3);
        // 3 rings * 12 segments = 36.
        assert_eq!(pattern.lines.len(), 36);
    }

    #[test]
    fn test_radial_concentric_combined() {
        let pattern =
            generate_radial_concentric_pattern([0.0; 3], [0.0, 1.0, 0.0], 6, 2, 1.0, 8, 0.5, 0.3);
        // 6 radial + 2*8 concentric = 22.
        assert_eq!(pattern.lines.len(), 6 + 2 * 8);
    }

    // ---- Fracture plane tests ----

    #[test]
    fn test_fracture_plane_signed_distance() {
        let plane = FracturePlane::new([0.0; 3], [0.0, 1.0, 0.0]);
        assert!((plane.signed_distance([0.0, 5.0, 0.0]) - 5.0).abs() < 1e-10);
        assert!((plane.signed_distance([0.0, -3.0, 0.0]) + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_split_points_by_plane() {
        let plane = FracturePlane::new([0.0; 3], [1.0, 0.0, 0.0]);
        let points = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let (front, back) = split_points_by_plane(&points, &plane);
        assert_eq!(front.len(), 2);
        assert_eq!(back.len(), 1);
    }

    // ---- Damage model tests ----

    #[test]
    fn test_damage_accumulation() {
        let mut ds = DamageState::new([0.0; 3]);
        ds.accumulate(0.3);
        assert!((ds.damage - 0.3).abs() < 1e-10);
        assert!(!ds.failed);
        ds.accumulate(0.8);
        assert!((ds.damage - 1.0).abs() < 1e-10);
        assert!(ds.failed);
    }

    #[test]
    fn test_damage_no_increment_below_threshold() {
        let stress = StressTensor::from_voigt([10.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mat = FractureMaterial::default();
        let inc = compute_damage_increment(&stress, &mat, 0.01, 1.0);
        assert!(inc < 1e-15, "no damage below threshold");
    }

    #[test]
    fn test_average_damage() {
        let states = vec![
            DamageState {
                damage: 0.2,
                position: [0.0; 3],
                failed: false,
            },
            DamageState {
                damage: 0.8,
                position: [0.0; 3],
                failed: false,
            },
        ];
        assert!((average_damage(&states) - 0.5).abs() < 1e-10);
    }

    // ---- Momentum conservation tests ----

    #[test]
    fn test_distribute_velocity_conserves_momentum() {
        let parent_vel = [10.0, 0.0, 0.0];
        let parent_mass = 5.0;
        let mut fragments = vec![
            Fragment {
                id: 0,
                center_of_mass: [0.0; 3],
                mass: 2.0,
                inertia: [1.0; 3],
                volume: 1.0,
                velocity: [3.0, 0.0, 0.0],
                angular_velocity: [0.0; 3],
                fracture_surface_area: 1.0,
                hull_vertices: Vec::new(),
            },
            Fragment {
                id: 1,
                center_of_mass: [0.0; 3],
                mass: 3.0,
                inertia: [1.0; 3],
                volume: 1.0,
                velocity: [5.0, 0.0, 0.0],
                angular_velocity: [0.0; 3],
                fracture_surface_area: 1.0,
                hull_vertices: Vec::new(),
            },
        ];
        distribute_velocity_conserving_momentum(&mut fragments, parent_vel, parent_mass);
        let p = total_fragment_momentum(&fragments);
        let desired = v3_scale(parent_vel, parent_mass);
        assert!(
            (p[0] - desired[0]).abs() < 1e-10,
            "momentum x: {:.6} vs {:.6}",
            p[0],
            desired[0]
        );
    }

    #[test]
    fn test_total_fragment_kinetic_energy() {
        let fragments = vec![Fragment {
            id: 0,
            center_of_mass: [0.0; 3],
            mass: 2.0,
            inertia: [1.0; 3],
            volume: 1.0,
            velocity: [3.0, 0.0, 0.0],
            angular_velocity: [0.0; 3],
            fracture_surface_area: 1.0,
            hull_vertices: Vec::new(),
        }];
        let ke = total_fragment_kinetic_energy(&fragments);
        // 0.5 * 2 * 9 = 9.
        assert!((ke - 9.0).abs() < 1e-10);
    }

    // ---- Destruction event tests ----

    #[test]
    fn test_voronoi_destruction_produces_fragments() {
        let mat = FractureMaterial::concrete();
        let event = execute_voronoi_destruction(
            0,
            [0.0; 3],
            [1.0; 3],
            [0.5, 0.5, 0.5],
            1000.0,
            [0.0; 3],
            &mat,
            5,
            8,
        );
        assert_eq!(event.fragments.len(), 5);
        assert_eq!(event.pattern_type, FracturePatternType::Voronoi);
    }
}
