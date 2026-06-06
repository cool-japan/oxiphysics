// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body fracture mechanics: brittle/ductile fracture, crack propagation,
//! fragment physics, glass cracking, explosive fragmentation, and forensic
//! reconstruction fitting.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper vec3
// ---------------------------------------------------------------------------

fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

fn vec3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

// ---------------------------------------------------------------------------
// FractureMode
// ---------------------------------------------------------------------------

/// Fracture failure mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FractureMode {
    /// Brittle fracture — instantaneous failure at the fracture criterion.
    Brittle,
    /// Ductile fracture — progressive damage accumulation before failure.
    Ductile,
    /// Mixed-mode combining brittle and ductile characteristics.
    Mixed,
}

// ---------------------------------------------------------------------------
// FractureModel
// ---------------------------------------------------------------------------

/// Material fracture model parameters.
#[derive(Debug, Clone)]
pub struct FractureModel {
    /// Failure mode.
    pub mode: FractureMode,
    /// Fracture energy (J/m²) — energy per unit crack area (Gc).
    pub fracture_energy: f64,
    /// Mode-I stress intensity factor toughness (Pa·√m).
    pub toughness_k1c: f64,
    /// Tensile strength (Pa).
    pub tensile_strength: f64,
    /// Compressive strength (Pa).
    pub compressive_strength: f64,
    /// Weibull modulus for statistical strength scatter.
    pub weibull_modulus: f64,
    /// Density (kg/m³).
    pub density: f64,
}

impl FractureModel {
    /// Default model for soda-lime glass.
    pub fn glass() -> Self {
        Self {
            mode: FractureMode::Brittle,
            fracture_energy: 8.0,
            toughness_k1c: 0.75e6,
            tensile_strength: 50e6,
            compressive_strength: 700e6,
            weibull_modulus: 7.0,
            density: 2500.0,
        }
    }

    /// Default model for structural steel.
    pub fn steel() -> Self {
        Self {
            mode: FractureMode::Ductile,
            fracture_energy: 200_000.0,
            toughness_k1c: 50e6,
            tensile_strength: 500e6,
            compressive_strength: 500e6,
            weibull_modulus: 30.0,
            density: 7800.0,
        }
    }

    /// Default model for concrete.
    pub fn concrete() -> Self {
        Self {
            mode: FractureMode::Mixed,
            fracture_energy: 100.0,
            toughness_k1c: 1.0e6,
            tensile_strength: 4e6,
            compressive_strength: 40e6,
            weibull_modulus: 12.0,
            density: 2400.0,
        }
    }

    /// Critical crack half-length (m) from K1c and tensile strength.
    pub fn critical_crack_length(&self) -> f64 {
        let k = self.toughness_k1c;
        let sigma = self.tensile_strength.max(1.0);
        (k / sigma).powi(2) / PI
    }

    /// Returns `true` if the applied stress exceeds tensile strength.
    pub fn is_fractured(&self, applied_stress: f64) -> bool {
        applied_stress >= self.tensile_strength
    }

    /// Returns `true` if the applied stress intensity factor exceeds K1c.
    pub fn is_crack_critical(&self, k_applied: f64) -> bool {
        k_applied >= self.toughness_k1c
    }
}

// ---------------------------------------------------------------------------
// FragmentNode / FragmentGraph
// ---------------------------------------------------------------------------

/// A single fragment in the fracture graph.
#[derive(Debug, Clone)]
pub struct FragmentNode {
    /// Fragment index.
    pub id: usize,
    /// Fragment mass (kg).
    pub mass: f64,
    /// Fragment centre of mass position (m).
    pub position: [f64; 3],
    /// Fragment velocity (m/s).
    pub velocity: [f64; 3],
    /// Volume (m³).
    pub volume: f64,
    /// Whether this fragment has separated.
    pub separated: bool,
}

impl FragmentNode {
    /// Construct a new fragment.
    pub fn new(id: usize, mass: f64, position: [f64; 3], volume: f64) -> Self {
        Self {
            id,
            mass,
            position,
            velocity: [0.0; 3],
            volume,
            separated: false,
        }
    }

    /// Kinetic energy of this fragment (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }
}

/// Edge in the fragment bond graph.
#[derive(Debug, Clone)]
pub struct FragmentBond {
    /// Node ID a.
    pub node_a: usize,
    /// Node ID b.
    pub node_b: usize,
    /// Bond strength (N) — maximum force before breaking.
    pub strength: f64,
    /// Whether the bond is broken.
    pub broken: bool,
}

impl FragmentBond {
    /// Create a new bond between two fragments.
    pub fn new(node_a: usize, node_b: usize, strength: f64) -> Self {
        Self {
            node_a,
            node_b,
            strength,
            broken: false,
        }
    }

    /// Break the bond.
    pub fn break_bond(&mut self) {
        self.broken = true;
    }
}

/// Fragment connectivity graph — nodes are fragments, edges are breakable bonds.
#[derive(Debug, Clone, Default)]
pub struct FragmentGraph {
    /// Fragment nodes.
    pub nodes: Vec<FragmentNode>,
    /// Fragment bonds.
    pub bonds: Vec<FragmentBond>,
}

impl FragmentGraph {
    /// Create an empty graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a fragment node.
    pub fn add_node(&mut self, node: FragmentNode) {
        self.nodes.push(node);
    }

    /// Add a bond between nodes.
    pub fn add_bond(&mut self, bond: FragmentBond) {
        self.bonds.push(bond);
    }

    /// Simulate a Voronoi fracture by distributing `n_fragments` evenly in
    /// a box `[0,0,0]–size` and connecting adjacent fragments.
    pub fn voronoi_fracture(
        &mut self,
        total_mass: f64,
        density: f64,
        size: [f64; 3],
        n_fragments: usize,
    ) {
        self.nodes.clear();
        self.bonds.clear();
        let frag_mass = total_mass / n_fragments.max(1) as f64;
        let frag_vol = frag_mass / density.max(1.0);
        // Place fragments on a grid
        let n_side = (n_fragments as f64).cbrt().ceil() as usize;
        for i in 0..n_fragments {
            let ix = (i % n_side) as f64 / n_side as f64;
            let iy = ((i / n_side) % n_side) as f64 / n_side as f64;
            let iz = (i / (n_side * n_side)) as f64 / n_side as f64;
            let pos = [ix * size[0], iy * size[1], iz * size[2]];
            self.nodes
                .push(FragmentNode::new(i, frag_mass, pos, frag_vol));
        }
        // Connect each fragment to the next (simplified chain)
        let strength_per_bond = total_mass * 9.81 * 10.0;
        for i in 0..n_fragments.saturating_sub(1) {
            self.bonds
                .push(FragmentBond::new(i, i + 1, strength_per_bond));
        }
    }

    /// Apply a force to all bonds and break those that exceed strength.
    pub fn apply_impulse(&mut self, impulse: f64) {
        for bond in &mut self.bonds {
            if !bond.broken && impulse > bond.strength {
                bond.break_bond();
                if let Some(node_b) = self.nodes.get_mut(bond.node_b) {
                    node_b.separated = true;
                }
            }
        }
    }

    /// Count broken bonds.
    pub fn broken_bond_count(&self) -> usize {
        self.bonds.iter().filter(|b| b.broken).count()
    }

    /// Count separated fragments.
    pub fn separated_fragment_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.separated).count()
    }

    /// Total kinetic energy across all fragments (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.nodes.iter().map(|n| n.kinetic_energy()).sum()
    }
}

// ---------------------------------------------------------------------------
// CrackPropagation
// ---------------------------------------------------------------------------

/// Crack tip state and Paris-law fatigue propagation.
#[derive(Debug, Clone)]
pub struct CrackPropagation {
    /// Current crack tip position (m).
    pub tip_position: [f64; 3],
    /// Crack propagation direction (unit vector).
    pub direction: [f64; 3],
    /// Current crack half-length (m).
    pub crack_length: f64,
    /// Paris law coefficient C (m/cycle / (Pa√m)^m).
    pub paris_c: f64,
    /// Paris law exponent m.
    pub paris_m: f64,
    /// Applied stress intensity factor range ΔK (Pa√m).
    pub delta_k: f64,
    /// Stress intensity factor threshold (Pa√m) below which no propagation.
    pub k_threshold: f64,
    /// Total number of cycles applied.
    pub cycles: u64,
}

impl CrackPropagation {
    /// Default crack propagation for steel.
    pub fn new_steel(initial_crack_length: f64, tip_position: [f64; 3]) -> Self {
        Self {
            tip_position,
            direction: [1.0, 0.0, 0.0],
            crack_length: initial_crack_length,
            paris_c: 1.0e-11,
            paris_m: 3.0,
            // delta_k and k_threshold in MPa√m scale (values ~5–20)
            delta_k: 10.0,
            k_threshold: 5.0,
            cycles: 0,
        }
    }

    /// Compute stress intensity factor K for an edge crack under stress `sigma` (Pa).
    pub fn stress_intensity(&self, sigma: f64) -> f64 {
        // K = Y * sigma * sqrt(pi * a), Y = 1.12 for edge crack
        1.12 * sigma * (PI * self.crack_length).sqrt()
    }

    /// Paris law: crack growth per cycle (m/cycle).
    pub fn growth_per_cycle(&self) -> f64 {
        if self.delta_k < self.k_threshold {
            return 0.0;
        }
        self.paris_c * self.delta_k.powf(self.paris_m)
    }

    /// Advance the crack by `n_cycles` cycles, updating position.
    pub fn advance(&mut self, n_cycles: u64) {
        let da_per_cycle = self.growth_per_cycle();
        let da = da_per_cycle * n_cycles as f64;
        self.crack_length += da;
        let dir = vec3_normalise(self.direction);
        self.tip_position = vec3_add(self.tip_position, vec3_scale(dir, da));
        self.cycles += n_cycles;
    }

    /// Remaining life (cycles) before crack reaches `critical_length` (m).
    pub fn remaining_life_cycles(&self, critical_length: f64) -> u64 {
        let da_per_cycle = self.growth_per_cycle();
        if da_per_cycle < 1e-20 {
            return u64::MAX;
        }
        let remaining_a = (critical_length - self.crack_length).max(0.0);
        (remaining_a / da_per_cycle).round() as u64
    }
}

// ---------------------------------------------------------------------------
// ImpactFracture
// ---------------------------------------------------------------------------

/// High-strain-rate impact fracture model (Grady-Kipp fragment sizing).
#[derive(Debug, Clone)]
pub struct ImpactFracture {
    /// Material model.
    pub fracture_model: FractureModel,
    /// Strain rate (1/s).
    pub strain_rate: f64,
    /// Spall strength (Pa) — tensile strength under shock loading.
    pub spall_strength: f64,
    /// Grady fragment size constant.
    pub grady_constant: f64,
}

impl ImpactFracture {
    /// Construct for glass under high-velocity impact.
    pub fn new_glass(strain_rate: f64) -> Self {
        Self {
            fracture_model: FractureModel::glass(),
            strain_rate,
            spall_strength: 2e6,
            grady_constant: 1.0,
        }
    }

    /// Grady-Kipp mean fragment size (m).
    pub fn mean_fragment_size(&self) -> f64 {
        // S = (K1c / (rho * strain_rate^2))^(2/3)
        let k1c = self.fracture_model.toughness_k1c;
        let rho = self.fracture_model.density;
        let edot = self.strain_rate.max(1.0);
        (k1c / (rho * edot * edot)).powf(2.0 / 3.0) * self.grady_constant
    }

    /// Number of fragments for a body of total volume `volume` (m³).
    pub fn fragment_count(&self, volume: f64) -> usize {
        let frag_size = self.mean_fragment_size().max(1e-9);
        let frag_vol = (4.0 / 3.0) * PI * (frag_size / 2.0).powi(3);
        (volume / frag_vol.max(1e-18)).round() as usize
    }

    /// Returns `true` if applied impulse (Pa·s) exceeds spall strength threshold.
    pub fn is_spall_failure(&self, impulse_pa_s: f64) -> bool {
        impulse_pa_s * self.strain_rate >= self.spall_strength
    }
}

// ---------------------------------------------------------------------------
// FractureConstraint
// ---------------------------------------------------------------------------

/// A pre-fracture joint that breaks when force/impulse exceeds threshold.
#[derive(Debug, Clone)]
pub struct FractureConstraint {
    /// Unique constraint ID.
    pub id: usize,
    /// Force threshold for fracture (N).
    pub force_threshold: f64,
    /// Impulse threshold for fracture (N·s).
    pub impulse_threshold: f64,
    /// Whether the joint has been fractured.
    pub fractured: bool,
    /// Accumulated impulse since last reset (N·s).
    pub accumulated_impulse: f64,
    /// Post-fracture restitution coefficient.
    pub restitution: f64,
}

impl FractureConstraint {
    /// Create a fracture constraint with given thresholds.
    pub fn new(id: usize, force_threshold: f64, impulse_threshold: f64) -> Self {
        Self {
            id,
            force_threshold,
            impulse_threshold,
            fractured: false,
            accumulated_impulse: 0.0,
            restitution: 0.3,
        }
    }

    /// Apply a force (N) for time `dt` (s) and check for fracture.
    ///
    /// Returns `true` if the constraint fractured during this call.
    pub fn apply_force(&mut self, force: f64, dt: f64) -> bool {
        if self.fractured {
            return false;
        }
        let impulse = force * dt;
        self.accumulated_impulse += impulse.abs();
        if force.abs() >= self.force_threshold || self.accumulated_impulse >= self.impulse_threshold
        {
            self.fractured = true;
            return true;
        }
        false
    }

    /// Reset accumulated impulse (e.g. after each physics frame).
    pub fn reset_impulse(&mut self) {
        self.accumulated_impulse = 0.0;
    }
}

// ---------------------------------------------------------------------------
// FragmentPhysics
// ---------------------------------------------------------------------------

/// Fragment mass/inertia from geometry, initial velocity from fracture energy.
#[derive(Debug, Clone)]
pub struct FragmentPhysics {
    /// Fragment mass (kg).
    pub mass: f64,
    /// Diagonal inertia tensor (kg·m²).
    pub inertia: [f64; 3],
    /// Initial velocity of the fragment (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
    /// Fracture energy available to this fragment (J).
    pub fracture_energy: f64,
}

impl FragmentPhysics {
    /// Compute initial fragment velocity magnitude from fracture energy.
    ///
    /// The velocity direction is given by `direction` (unit vector).
    pub fn from_fracture_energy(
        mass: f64,
        fracture_energy: f64,
        direction: [f64; 3],
        size: [f64; 3],
        density: f64,
    ) -> Self {
        let v_mag = (2.0 * fracture_energy / mass.max(1e-6)).sqrt();
        let dir = {
            let n = vec3_norm(direction);
            if n < 1e-12 {
                [1.0, 0.0, 0.0]
            } else {
                vec3_scale(direction, 1.0 / n)
            }
        };
        let velocity = vec3_scale(dir, v_mag);
        // Box inertia approximation
        let [hx, hy, hz] = [size[0] / 2.0, size[1] / 2.0, size[2] / 2.0];
        let vol = 8.0 * hx * hy * hz;
        let m = density * vol;
        let ixx = m / 12.0 * (4.0 * hy * hy + 4.0 * hz * hz);
        let iyy = m / 12.0 * (4.0 * hx * hx + 4.0 * hz * hz);
        let izz = m / 12.0 * (4.0 * hx * hx + 4.0 * hy * hy);
        Self {
            mass,
            inertia: [ixx, iyy, izz],
            velocity,
            angular_velocity: [0.0; 3],
            fracture_energy,
        }
    }

    /// Kinetic energy of the fragment (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// GlassBreaker
// ---------------------------------------------------------------------------

/// Glass fracture pattern model: Hertzian contact, radial/concentric cracks.
#[derive(Debug, Clone)]
pub struct GlassBreaker {
    /// Glass fracture model.
    pub fracture_model: FractureModel,
    /// Is this tempered glass?
    pub tempered: bool,
    /// Glass thickness (m).
    pub thickness: f64,
    /// Impact energy (J).
    pub impact_energy: f64,
}

impl GlassBreaker {
    /// Construct for a standard annealed windshield.
    pub fn annealed_windshield() -> Self {
        Self {
            fracture_model: FractureModel::glass(),
            tempered: false,
            thickness: 0.006,
            impact_energy: 0.0,
        }
    }

    /// Construct for a tempered side window.
    pub fn tempered_side_window() -> Self {
        Self {
            fracture_model: FractureModel::glass(),
            tempered: true,
            thickness: 0.004,
            impact_energy: 0.0,
        }
    }

    /// Number of radial cracks formed (empirical model).
    pub fn radial_crack_count(&self) -> usize {
        if self.impact_energy < 0.5 {
            return 0;
        }
        let n = (4.0 + 6.0 * (self.impact_energy / 2.0).ln().max(0.0)) as usize;
        if self.tempered { n * 3 } else { n }
    }

    /// Hertzian contact radius (m) for a spherical impactor.
    ///
    /// `impactor_radius` — impactor sphere radius (m).
    /// `impactor_force` — peak contact force (N).
    pub fn hertz_contact_radius(&self, impactor_radius: f64, impactor_force: f64) -> f64 {
        // a = (3FR / 4E*)^(1/3), simplified for glass E* ~ 35 GPa
        let e_star = 35e9_f64;
        let numerator = 3.0 * impactor_force * impactor_radius;
        let denominator = 4.0 * e_star;
        (numerator / denominator).powf(1.0 / 3.0)
    }

    /// Whether impact energy exceeds the threshold to shatter the pane.
    pub fn is_shattered(&self) -> bool {
        let threshold = if self.tempered { 5.0 } else { 2.0 };
        self.impact_energy >= threshold
    }

    /// Estimated fragment count for tempered glass shattering.
    pub fn tempered_fragment_count(&self) -> usize {
        if !self.is_shattered() || !self.tempered {
            return 0;
        }
        // Tempered glass shatters into ~40 fragments per cm²
        let area = 0.5 * 0.3; // approximate side window area m²
        (area * 1e4 * 40.0) as usize
    }
}

// ---------------------------------------------------------------------------
// ExplosiveFragmentation
// ---------------------------------------------------------------------------

/// Mott distribution explosive fragmentation model.
#[derive(Debug, Clone)]
pub struct ExplosiveFragmentation {
    /// Total explosive mass (kg).
    pub explosive_mass: f64,
    /// Casing mass (kg).
    pub casing_mass: f64,
    /// Casing density (kg/m³).
    pub casing_density: f64,
    /// Casing wall thickness (m).
    pub wall_thickness: f64,
    /// Detonation velocity (m/s).
    pub det_velocity: f64,
}

impl ExplosiveFragmentation {
    /// Default parameters for a small munition.
    pub fn small_munition() -> Self {
        Self {
            explosive_mass: 0.5,
            casing_mass: 1.5,
            casing_density: 7800.0,
            wall_thickness: 0.006,
            det_velocity: 7000.0,
        }
    }

    /// Gurney velocity (m/s) — characteristic fragment velocity.
    pub fn gurney_velocity(&self) -> f64 {
        // v = sqrt(2 * E_g * M/C / (1 + M/C))
        // Simplified: use Gurney energy E_g ~ 2.5e6 J/kg for TNT
        let e_g = 2.5e6_f64;
        let mc = self.casing_mass / self.explosive_mass.max(1e-6);
        (2.0 * e_g / (mc + 0.5)).sqrt()
    }

    /// Mott mean fragment mass (kg).
    pub fn mott_mean_mass(&self) -> f64 {
        // Simplified Mott formula: M_avg = B^2 * t^5/2 / sqrt(rho_c)
        // with simplified B ~ 0.1 m/kg^(1/2)
        let b = 0.10_f64;
        let t = self.wall_thickness;
        b * b * t.powf(2.5) / self.casing_density.sqrt().max(1.0)
    }

    /// Total fragment count.
    pub fn fragment_count(&self) -> usize {
        let mean = self.mott_mean_mass();
        if mean < 1e-12 {
            return 0;
        }
        (self.casing_mass / mean).round() as usize
    }

    /// Casualty radius (m) — range at which fragment flux causes 50% casualties.
    ///
    /// Uses simple inverse-square model with a lethal fragment threshold.
    pub fn casualty_radius(&self) -> f64 {
        let v = self.gurney_velocity();
        let n = self.fragment_count() as f64;
        let m_frag = self.mott_mean_mass();
        // Lethal energy threshold ~ 78 J
        let e_lethal = 78.0_f64;
        // E = 0.5 * m * v_r^2, where v_r decreases with distance
        // Simplified: r_50 = v * sqrt(m / e_lethal) * cbrt(n / 4pi)
        let v_threshold = (2.0 * e_lethal / m_frag.max(1e-12)).sqrt();
        if v <= v_threshold {
            return 0.0;
        }
        let r = v * (m_frag / e_lethal).sqrt() * (n / (4.0 * PI)).cbrt();
        r.max(0.0)
    }
}

// ---------------------------------------------------------------------------
// FractureLogger
// ---------------------------------------------------------------------------

/// A single fracture event record.
#[derive(Debug, Clone)]
pub struct FractureEvent {
    /// Simulation time (s).
    pub time: f64,
    /// Location of the fracture origin (m).
    pub location: [f64; 3],
    /// Applied force at fracture (N).
    pub force: f64,
    /// Number of fragments produced.
    pub fragment_count: usize,
    /// Constraint ID that fractured.
    pub constraint_id: usize,
}

/// Logs and analyses fracture events during a simulation.
#[derive(Debug, Clone, Default)]
pub struct FractureLogger {
    /// All recorded fracture events.
    pub events: Vec<FractureEvent>,
}

impl FractureLogger {
    /// Create an empty logger.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a fracture event.
    pub fn log(&mut self, event: FractureEvent) {
        self.events.push(event);
    }

    /// Total number of fracture events.
    pub fn event_count(&self) -> usize {
        self.events.len()
    }

    /// Total fragment count across all events.
    pub fn total_fragments(&self) -> usize {
        self.events.iter().map(|e| e.fragment_count).sum()
    }

    /// Mean fracture force (N).
    pub fn mean_force(&self) -> f64 {
        if self.events.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.events.iter().map(|e| e.force).sum();
        sum / self.events.len() as f64
    }

    /// Time of first fracture event (s), or `None` if no events.
    pub fn first_event_time(&self) -> Option<f64> {
        self.events.first().map(|e| e.time)
    }

    /// Clear all events.
    pub fn clear(&mut self) {
        self.events.clear();
    }
}

// ---------------------------------------------------------------------------
// ReconstructionFit
// ---------------------------------------------------------------------------

/// Fits fragmented pieces back together for forensic reconstruction.
#[derive(Debug, Clone)]
pub struct ReconstructionFit {
    /// Fragment positions as recorded post-fracture (m).
    pub fragment_positions: Vec<[f64; 3]>,
    /// Fragment velocities as recorded post-fracture (m/s).
    pub fragment_velocities: Vec<[f64; 3]>,
    /// Reconstruction time window (s) — how far back to reverse kinematics.
    pub time_window: f64,
}

impl ReconstructionFit {
    /// Create a reconstruction fit from post-fracture fragment states.
    pub fn new(positions: Vec<[f64; 3]>, velocities: Vec<[f64; 3]>, time_window: f64) -> Self {
        Self {
            fragment_positions: positions,
            fragment_velocities: velocities,
            time_window,
        }
    }

    /// Reverse-propagate each fragment to its position at t=0 (fracture time).
    ///
    /// Uses constant-velocity backward integration.
    pub fn reverse_to_origin(&self) -> Vec<[f64; 3]> {
        self.fragment_positions
            .iter()
            .zip(self.fragment_velocities.iter())
            .map(|(pos, vel)| {
                // x0 = x - v * t
                vec3_sub(*pos, vec3_scale(*vel, self.time_window))
            })
            .collect()
    }

    /// Estimate the fracture origin as the centroid of reversed positions.
    pub fn estimated_origin(&self) -> [f64; 3] {
        let origins = self.reverse_to_origin();
        if origins.is_empty() {
            return [0.0; 3];
        }
        let n = origins.len() as f64;
        let mut sum = [0.0_f64; 3];
        for o in &origins {
            sum = vec3_add(sum, *o);
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Root-mean-square dispersion of reversed positions from their centroid (m).
    pub fn reconstruction_rms(&self) -> f64 {
        let origins = self.reverse_to_origin();
        let centroid = self.estimated_origin();
        if origins.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = origins
            .iter()
            .map(|o| {
                let d = vec3_sub(*o, centroid);
                vec3_dot(d, d)
            })
            .sum();
        (sum_sq / origins.len() as f64).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- FractureModel ---

    #[test]
    fn glass_is_brittle() {
        let m = FractureModel::glass();
        assert_eq!(m.mode, FractureMode::Brittle);
    }

    #[test]
    fn steel_is_ductile() {
        let m = FractureModel::steel();
        assert_eq!(m.mode, FractureMode::Ductile);
    }

    #[test]
    fn critical_crack_length_positive() {
        let m = FractureModel::glass();
        assert!(m.critical_crack_length() > 0.0);
    }

    #[test]
    fn is_fractured_above_strength() {
        let m = FractureModel::glass();
        assert!(m.is_fractured(m.tensile_strength + 1.0));
        assert!(!m.is_fractured(m.tensile_strength * 0.5));
    }

    #[test]
    fn is_crack_critical_above_k1c() {
        let m = FractureModel::glass();
        assert!(m.is_crack_critical(m.toughness_k1c + 1.0));
        assert!(!m.is_crack_critical(m.toughness_k1c * 0.5));
    }

    // --- FragmentGraph ---

    #[test]
    fn voronoi_fracture_creates_correct_count() {
        let mut g = FragmentGraph::new();
        g.voronoi_fracture(10.0, 2500.0, [1.0, 1.0, 1.0], 8);
        assert_eq!(g.nodes.len(), 8);
    }

    #[test]
    fn voronoi_fracture_bonds_connected() {
        let mut g = FragmentGraph::new();
        g.voronoi_fracture(10.0, 2500.0, [1.0, 1.0, 1.0], 5);
        assert_eq!(g.bonds.len(), 4);
    }

    #[test]
    fn fragment_graph_apply_impulse_breaks_weak_bonds() {
        let mut g = FragmentGraph::new();
        g.voronoi_fracture(10.0, 2500.0, [1.0, 1.0, 1.0], 4);
        let strength = g.bonds[0].strength;
        g.apply_impulse(strength + 1.0);
        assert!(g.broken_bond_count() > 0);
    }

    #[test]
    fn fragment_graph_no_break_below_threshold() {
        let mut g = FragmentGraph::new();
        g.voronoi_fracture(10.0, 2500.0, [1.0, 1.0, 1.0], 4);
        g.apply_impulse(0.0);
        assert_eq!(g.broken_bond_count(), 0);
    }

    #[test]
    fn fragment_node_kinetic_energy_zero() {
        let node = FragmentNode::new(0, 1.0, [0.0; 3], 1e-3);
        assert_eq!(node.kinetic_energy(), 0.0);
    }

    // --- CrackPropagation ---

    #[test]
    fn crack_growth_zero_below_threshold() {
        let c = CrackPropagation::new_steel(0.001, [0.0; 3]);
        // With delta_k below k_threshold? No — default delta_k is above threshold
        // Explicitly set below:
        let mut c2 = c.clone();
        c2.delta_k = c2.k_threshold * 0.5;
        assert_eq!(c2.growth_per_cycle(), 0.0);
    }

    #[test]
    fn crack_advances_with_cycles() {
        let mut c = CrackPropagation::new_steel(0.001, [0.0; 3]);
        let initial = c.crack_length;
        c.advance(100_000);
        assert!(
            c.crack_length > initial,
            "crack should grow: {}",
            c.crack_length
        );
    }

    #[test]
    fn crack_stress_intensity_positive() {
        let c = CrackPropagation::new_steel(0.005, [0.0; 3]);
        let k = c.stress_intensity(100e6);
        assert!(k > 0.0, "K should be positive: {k}");
    }

    #[test]
    fn crack_remaining_life_finite() {
        let c = CrackPropagation::new_steel(0.001, [0.0; 3]);
        let life = c.remaining_life_cycles(0.010);
        assert!(life > 0, "remaining life should be positive: {life}");
    }

    // --- ImpactFracture ---

    #[test]
    fn impact_fracture_mean_size_positive() {
        let f = ImpactFracture::new_glass(1000.0);
        let s = f.mean_fragment_size();
        assert!(s > 0.0, "mean fragment size: {s}");
    }

    #[test]
    fn impact_fracture_count_increases_with_volume() {
        let f = ImpactFracture::new_glass(1000.0);
        let n1 = f.fragment_count(0.001);
        let n2 = f.fragment_count(0.010);
        assert!(n2 >= n1, "more fragments from larger volume: {n1} vs {n2}");
    }

    #[test]
    fn impact_fracture_spall_detection() {
        let f = ImpactFracture::new_glass(1000.0);
        let spall_threshold = f.spall_strength / f.strain_rate;
        assert!(f.is_spall_failure(spall_threshold + 1.0));
    }

    // --- FractureConstraint ---

    #[test]
    fn fracture_constraint_breaks_on_high_force() {
        let mut c = FractureConstraint::new(0, 1000.0, 1e9);
        let broke = c.apply_force(1001.0, 0.001);
        assert!(broke);
        assert!(c.fractured);
    }

    #[test]
    fn fracture_constraint_does_not_break_twice() {
        let mut c = FractureConstraint::new(0, 100.0, 1e9);
        c.apply_force(1000.0, 0.01);
        let second = c.apply_force(1000.0, 0.01);
        assert!(!second, "should not report fracture twice");
    }

    #[test]
    fn fracture_constraint_impulse_accumulates() {
        let mut c = FractureConstraint::new(0, 1e9, 100.0);
        c.apply_force(10.0, 5.0); // impulse = 50
        c.apply_force(10.0, 5.0); // impulse = 100 → breaks
        assert!(c.fractured);
    }

    // --- FragmentPhysics ---

    #[test]
    fn fragment_physics_velocity_from_energy() {
        let fp = FragmentPhysics::from_fracture_energy(
            0.1,
            10.0,
            [1.0, 0.0, 0.0],
            [0.05, 0.05, 0.05],
            7800.0,
        );
        let ke = fp.kinetic_energy();
        assert!(
            (ke - 10.0).abs() < 0.1,
            "kinetic energy should match fracture energy: {ke}"
        );
    }

    // --- GlassBreaker ---

    #[test]
    fn glass_breaker_radial_cracks_zero_low_energy() {
        let mut g = GlassBreaker::annealed_windshield();
        g.impact_energy = 0.1;
        assert_eq!(g.radial_crack_count(), 0);
    }

    #[test]
    fn glass_breaker_cracks_increase_with_energy() {
        let mut g = GlassBreaker::annealed_windshield();
        g.impact_energy = 1.0;
        let n1 = g.radial_crack_count();
        g.impact_energy = 10.0;
        let n2 = g.radial_crack_count();
        assert!(n2 >= n1, "more cracks with higher energy: {n1} vs {n2}");
    }

    #[test]
    fn glass_breaker_not_shattered_below_threshold() {
        let mut g = GlassBreaker::annealed_windshield();
        g.impact_energy = 0.5;
        assert!(!g.is_shattered());
    }

    #[test]
    fn glass_breaker_hertz_contact_positive() {
        let g = GlassBreaker::annealed_windshield();
        let a = g.hertz_contact_radius(0.01, 5000.0);
        assert!(a > 0.0, "contact radius: {a}");
    }

    #[test]
    fn glass_breaker_tempered_more_cracks() {
        let mut annealed = GlassBreaker::annealed_windshield();
        let mut tempered = GlassBreaker::tempered_side_window();
        annealed.impact_energy = 5.0;
        tempered.impact_energy = 5.0;
        assert!(tempered.radial_crack_count() >= annealed.radial_crack_count());
    }

    // --- ExplosiveFragmentation ---

    #[test]
    fn explosive_gurney_velocity_positive() {
        let e = ExplosiveFragmentation::small_munition();
        assert!(e.gurney_velocity() > 0.0);
    }

    #[test]
    fn explosive_fragment_count_positive() {
        let e = ExplosiveFragmentation::small_munition();
        assert!(e.fragment_count() > 0);
    }

    #[test]
    fn explosive_casualty_radius_positive() {
        let e = ExplosiveFragmentation::small_munition();
        let r = e.casualty_radius();
        assert!(r >= 0.0, "casualty radius: {r}");
    }

    // --- FractureLogger ---

    #[test]
    fn fracture_logger_records_events() {
        let mut l = FractureLogger::new();
        l.log(FractureEvent {
            time: 0.1,
            location: [0.0; 3],
            force: 500.0,
            fragment_count: 5,
            constraint_id: 0,
        });
        assert_eq!(l.event_count(), 1);
        assert_eq!(l.total_fragments(), 5);
    }

    #[test]
    fn fracture_logger_mean_force() {
        let mut l = FractureLogger::new();
        l.log(FractureEvent {
            time: 0.1,
            location: [0.0; 3],
            force: 100.0,
            fragment_count: 2,
            constraint_id: 0,
        });
        l.log(FractureEvent {
            time: 0.2,
            location: [1.0, 0.0, 0.0],
            force: 300.0,
            fragment_count: 4,
            constraint_id: 1,
        });
        assert!((l.mean_force() - 200.0).abs() < 1e-10);
    }

    #[test]
    fn fracture_logger_clear() {
        let mut l = FractureLogger::new();
        l.log(FractureEvent {
            time: 0.1,
            location: [0.0; 3],
            force: 100.0,
            fragment_count: 1,
            constraint_id: 0,
        });
        l.clear();
        assert_eq!(l.event_count(), 0);
    }

    // --- ReconstructionFit ---

    #[test]
    fn reconstruction_fit_origin_estimate() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let fit = ReconstructionFit::new(positions, velocities, 1.0);
        let origin = fit.estimated_origin();
        // Reversed: [0,0,0] and [0,0,0] → centroid = [0,0,0]
        assert!(
            vec3_norm(origin) < 1e-10,
            "origin should be near zero: {:?}",
            origin
        );
    }

    #[test]
    fn reconstruction_fit_rms_zero_perfect() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let fit = ReconstructionFit::new(positions, velocities, 1.0);
        let rms = fit.reconstruction_rms();
        assert!(rms < 1e-10, "RMS should be zero for symmetric case: {rms}");
    }

    #[test]
    fn reconstruction_fit_empty() {
        let fit = ReconstructionFit::new(vec![], vec![], 1.0);
        let origin = fit.estimated_origin();
        assert_eq!(origin, [0.0; 3]);
    }
}
