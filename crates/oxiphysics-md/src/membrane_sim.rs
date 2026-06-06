// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lipid bilayer membrane simulation.
//!
//! Provides coarse-grained lipid models, bilayer construction, MD stepping,
//! and analysis (order parameter, area per lipid, bending modulus, tension).

use std::f64::consts::PI;

/// Type of lipid molecule in the bilayer.
#[derive(Debug, Clone, PartialEq)]
pub enum LipidType {
    /// Dipalmitoylphosphatidylcholine — saturated, forms ordered gel phase.
    Dppc,
    /// Dioleoylphosphatidylcholine — unsaturated, forms fluid phase.
    Dopc,
    /// Cholesterol — modulates membrane fluidity and order.
    Cholesterol,
    /// Custom lipid with a user-specified name.
    Custom(String),
}

impl LipidType {
    /// Returns the typical area per lipid in nm² at 300 K.
    pub fn reference_area(&self) -> f64 {
        match self {
            LipidType::Dppc => 0.48,
            LipidType::Dopc => 0.72,
            LipidType::Cholesterol => 0.38,
            LipidType::Custom(_) => 0.60,
        }
    }

    /// Returns the typical hydrophobic tail length in nm.
    pub fn hydrophobic_length(&self) -> f64 {
        match self {
            LipidType::Dppc => 1.67,
            LipidType::Dopc => 1.49,
            LipidType::Cholesterol => 1.70,
            LipidType::Custom(_) => 1.50,
        }
    }
}

/// A single lipid molecule in the coarse-grained model.
///
/// Position is the headgroup bead; orientation points from head toward tail.
#[derive(Debug, Clone)]
pub struct Lipid {
    /// Position of the headgroup bead \[x, y, z\] in nm.
    pub position: [f64; 3],
    /// Unit orientation vector (head → tail direction).
    pub orientation: [f64; 3],
    /// Tail order parameter S₂ = 1.5 cos²θ - 0.5 (range -0.5 to 1.0).
    pub tail_order: f64,
    /// Area per lipid in nm² (assigned during analysis).
    pub area_per_lipid: f64,
    /// Type of this lipid.
    pub lipid_type: LipidType,
}

impl Lipid {
    /// Create a new lipid with default tail_order and area_per_lipid.
    pub fn new(position: [f64; 3], orientation: [f64; 3], lipid_type: LipidType) -> Self {
        let apl = lipid_type.reference_area();
        Self {
            position,
            orientation: normalise(orientation),
            tail_order: 0.0,
            area_per_lipid: apl,
            lipid_type,
        }
    }

    /// Tilt angle (radians) between orientation and the membrane normal \[0,0,1\].
    pub fn tilt_angle(&self) -> f64 {
        let cos_theta = self.orientation[2].clamp(-1.0, 1.0);
        cos_theta.abs().acos()
    }

    /// Tail order parameter from the current orientation relative to z-axis.
    pub fn compute_order(&self) -> f64 {
        let cos_theta = self.orientation[2];
        1.5 * cos_theta * cos_theta - 0.5
    }
}

/// Normalise a 3-vector (returns unchanged if zero).
fn normalise(v: [f64; 3]) -> [f64; 3] {
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag < 1e-15 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / mag, v[1] / mag, v[2] / mag]
    }
}

/// Mechanical properties of a membrane patch.
#[derive(Debug, Clone)]
pub struct MembraneProperties {
    /// Bending modulus κ in k_B T.
    pub bending_modulus: f64,
    /// Lateral tension γ in mN/m.
    pub tension: f64,
    /// Area compressibility modulus K_A in mN/m.
    pub area_compressibility: f64,
}

impl MembraneProperties {
    /// Default properties for a typical phospholipid bilayer.
    pub fn default_phospholipid() -> Self {
        Self {
            bending_modulus: 20.0,       // k_BT
            tension: 0.0,                // zero-tension equilibrium
            area_compressibility: 230.0, // mN/m
        }
    }

    /// Create custom membrane properties.
    pub fn new(bending_modulus: f64, tension: f64, area_compressibility: f64) -> Self {
        Self {
            bending_modulus,
            tension,
            area_compressibility,
        }
    }
}

/// A transmembrane protein embedded in the bilayer.
#[derive(Debug, Clone)]
pub struct ProteinInsertion {
    /// Centre of mass of the protein \[x, y, z\] in nm.
    pub position: [f64; 3],
    /// Tilt angle of the protein helix axis relative to membrane normal (radians).
    pub tilt_angle: f64,
    /// Hydrophobic thickness of the transmembrane domain in nm.
    pub hydrophobic_thickness: f64,
    /// Radius of the transmembrane cylinder in nm.
    pub radius: f64,
    /// Mismatch between protein hydrophobic thickness and bilayer thickness.
    pub hydrophobic_mismatch: f64,
}

impl ProteinInsertion {
    /// Create a protein insertion with a given tilt and hydrophobic thickness.
    pub fn new(
        position: [f64; 3],
        tilt_angle: f64,
        hydrophobic_thickness: f64,
        bilayer_thickness: f64,
        radius: f64,
    ) -> Self {
        Self {
            position,
            tilt_angle,
            hydrophobic_thickness,
            radius,
            hydrophobic_mismatch: hydrophobic_thickness - bilayer_thickness,
        }
    }

    /// Deformation energy cost of inserting this protein (harmonic model).
    ///
    /// E = 0.5 * K_A * mismatch² / A_protein
    pub fn deformation_energy(&self, area_compressibility: f64) -> f64 {
        let area = PI * self.radius * self.radius;
        if area == 0.0 {
            return 0.0;
        }
        0.5 * area_compressibility * self.hydrophobic_mismatch.powi(2) / area
    }
}

/// Membrane fluidity characterisation.
#[derive(Debug, Clone)]
pub struct MembraneFluidity {
    /// Lateral diffusion coefficient D in nm²/ns.
    pub diffusion_coefficient: f64,
    /// Flip-flop rate (leaflet exchange rate) in ns⁻¹.
    pub flip_flop_rate: f64,
    /// Temperature at which these values were measured (K).
    pub temperature: f64,
}

impl MembraneFluidity {
    /// Create fluidity parameters.
    pub fn new(diffusion_coefficient: f64, flip_flop_rate: f64, temperature: f64) -> Self {
        Self {
            diffusion_coefficient,
            flip_flop_rate,
            temperature,
        }
    }

    /// Mean squared displacement at time t: MSD = 4 D t.
    pub fn msd(&self, t: f64) -> f64 {
        4.0 * self.diffusion_coefficient * t
    }

    /// Half-time for flip-flop (ns): t_{1/2} = ln(2) / k_ff.
    pub fn flip_flop_half_time(&self) -> f64 {
        if self.flip_flop_rate == 0.0 {
            f64::INFINITY
        } else {
            2_f64.ln() / self.flip_flop_rate
        }
    }
}

/// A lipid bilayer membrane patch.
///
/// Contains two leaflets (upper: z > 0, lower: z < 0). Periodic boundary
/// conditions are applied in x and y; z is the membrane normal direction.
#[derive(Debug, Clone)]
pub struct Membrane {
    /// All lipid molecules (upper leaflet first, then lower).
    pub lipids: Vec<Lipid>,
    /// Number of lipids per leaflet.
    pub lipids_per_leaflet: usize,
    /// Bilayer thickness in nm (distance between leaflet headgroup planes).
    pub bilayer_thickness: f64,
    /// Total projected area of the patch in nm².
    pub area: f64,
    /// Box size \[Lx, Ly\] in nm.
    pub box_xy: [f64; 2],
    /// Mechanical properties of this membrane.
    pub properties: MembraneProperties,
    /// Embedded proteins.
    pub proteins: Vec<ProteinInsertion>,
    /// Current simulation time step.
    pub step_count: u64,
}

impl Membrane {
    /// Build a flat bilayer patch with `n_per_leaflet` lipids per leaflet.
    ///
    /// Lipids are placed on a square grid. Upper leaflet at z = +thickness/2,
    /// lower leaflet at z = -thickness/2 with inverted orientation.
    pub fn build_bilayer(
        n_per_leaflet: usize,
        lipid_type: LipidType,
        box_xy: [f64; 2],
        bilayer_thickness: f64,
        properties: MembraneProperties,
    ) -> Self {
        let mut lipids = Vec::with_capacity(2 * n_per_leaflet);
        let cols = (n_per_leaflet as f64).sqrt().ceil() as usize;
        let dx = box_xy[0] / cols as f64;
        let dy = box_xy[1] / n_per_leaflet.div_ceil(cols) as f64;
        let z_upper = bilayer_thickness / 2.0;
        let z_lower = -bilayer_thickness / 2.0;

        for i in 0..n_per_leaflet {
            let col = i % cols;
            let row = i / cols;
            let x = (col as f64 + 0.5) * dx;
            let y = (row as f64 + 0.5) * dy;
            lipids.push(Lipid::new(
                [x, y, z_upper],
                [0.0, 0.0, -1.0], // head up, tail down
                lipid_type.clone(),
            ));
        }
        for i in 0..n_per_leaflet {
            let col = i % cols;
            let row = i / cols;
            let x = (col as f64 + 0.5) * dx;
            let y = (row as f64 + 0.5) * dy;
            lipids.push(Lipid::new(
                [x, y, z_lower],
                [0.0, 0.0, 1.0], // head down, tail up
                lipid_type.clone(),
            ));
        }

        let area = box_xy[0] * box_xy[1];
        Self {
            lipids,
            lipids_per_leaflet: n_per_leaflet,
            bilayer_thickness,
            area,
            box_xy,
            properties,
            proteins: Vec::new(),
            step_count: 0,
        }
    }

    /// Perform one MD step: apply thermal noise to orientations and positions.
    ///
    /// Each lipid undergoes Brownian translational diffusion and rotational
    /// fluctuation scaled by temperature `kT` (in units of energy).
    pub fn step_md(&mut self, kt: f64, dt: f64, diffusion: f64) {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        // Deterministic pseudo-random noise using step count + index hash
        for (i, lipid) in self.lipids.iter_mut().enumerate() {
            let mut h = DefaultHasher::new();
            (self.step_count, i, 0u8).hash(&mut h);
            let bits = h.finish();
            // Map hash to [-1, 1]
            let noise_x = ((bits & 0xFFFF) as f64 / 32767.5) - 1.0;
            let noise_y = (((bits >> 16) & 0xFFFF) as f64 / 32767.5) - 1.0;
            let noise_z = (((bits >> 32) & 0xFFFF) as f64 / 32767.5) - 1.0;

            let amp = (2.0 * diffusion * dt / kt).sqrt() * kt.sqrt();
            lipid.position[0] += amp * noise_x;
            lipid.position[1] += amp * noise_y;
            // Apply PBC in xy
            lipid.position[0] = lipid.position[0].rem_euclid(self.box_xy[0]);
            lipid.position[1] = lipid.position[1].rem_euclid(self.box_xy[1]);

            // Small rotational fluctuation
            let rot_amp = (kt / 100.0).sqrt();
            let new_orient = [
                lipid.orientation[0] + rot_amp * noise_x,
                lipid.orientation[1] + rot_amp * noise_y,
                lipid.orientation[2] + rot_amp * noise_z * 0.1,
            ];
            lipid.orientation = normalise(new_orient);
            lipid.tail_order = lipid.compute_order();
        }
        self.step_count += 1;
    }

    /// Insert a protein into the membrane.
    pub fn insert_protein(&mut self, protein: ProteinInsertion) {
        self.proteins.push(protein);
    }

    /// Recompute the bilayer thickness from the mean z-positions of the two leaflets.
    pub fn update_bilayer_thickness(&mut self) {
        if self.lipids_per_leaflet == 0 {
            return;
        }
        let upper_z: f64 = self.lipids[..self.lipids_per_leaflet]
            .iter()
            .map(|l| l.position[2])
            .sum::<f64>()
            / self.lipids_per_leaflet as f64;
        let lower_z: f64 = self.lipids[self.lipids_per_leaflet..]
            .iter()
            .map(|l| l.position[2])
            .sum::<f64>()
            / self.lipids_per_leaflet as f64;
        self.bilayer_thickness = (upper_z - lower_z).abs();
    }

    /// Count of all lipids (both leaflets).
    pub fn total_lipids(&self) -> usize {
        self.lipids.len()
    }
}

/// Compute the deuterium order parameter S_CD for a lipid assembly.
///
/// S_CD = ⟨1.5 cos²θ - 0.5⟩, where θ is the angle between the C–D bond
/// (approximated by the lipid orientation) and the bilayer normal (z-axis).
///
/// Returns a value in \[-0.5, 1.0\]; a fluid bilayer gives ~0.2–0.3.
pub fn compute_order_parameter(lipids: &[Lipid]) -> f64 {
    if lipids.is_empty() {
        return 0.0;
    }
    let sum: f64 = lipids.iter().map(|l| l.compute_order()).sum();
    sum / lipids.len() as f64
}

/// Compute the area per lipid for a membrane by dividing the projected box area
/// by the number of lipids in one leaflet.
pub fn compute_area_per_lipid(membrane: &Membrane) -> f64 {
    if membrane.lipids_per_leaflet == 0 {
        return 0.0;
    }
    membrane.area / membrane.lipids_per_leaflet as f64
}

/// Compute the deuterium order parameter profile along the membrane normal.
///
/// Bins lipids into `n_bins` z-slices and returns (z_centres, S_CD_per_bin).
pub fn order_parameter_profile(
    lipids: &[Lipid],
    z_min: f64,
    z_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dz = (z_max - z_min) / n_bins as f64;
    let mut counts = vec![0usize; n_bins];
    let mut sums = vec![0.0f64; n_bins];
    for l in lipids {
        let z = l.position[2];
        if z < z_min || z >= z_max {
            continue;
        }
        let bin = ((z - z_min) / dz) as usize;
        let bin = bin.min(n_bins - 1);
        counts[bin] += 1;
        sums[bin] += l.compute_order();
    }
    let z_centres: Vec<f64> = (0..n_bins).map(|i| z_min + (i as f64 + 0.5) * dz).collect();
    let s_cd: Vec<f64> = sums
        .iter()
        .zip(counts.iter())
        .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
        .collect();
    (z_centres, s_cd)
}

/// Estimate membrane bending rigidity from area fluctuation spectrum (simplified).
///
/// κ ≈ k_B T * q^4 * A / ⟨|h(q)|²⟩  (Helfrich model, single-mode approximation)
///
/// Here we use a simplified formula: κ = k_BT * A_proj / (8π²σ_h²)
/// where σ_h is the standard deviation of lipid z-positions.
pub fn estimate_bending_modulus(lipids: &[Lipid], kbt: f64, projected_area: f64) -> f64 {
    if lipids.len() < 2 {
        return 0.0;
    }
    let n = lipids.len() as f64;
    let mean_z = lipids.iter().map(|l| l.position[2]).sum::<f64>() / n;
    let var_z = lipids
        .iter()
        .map(|l| (l.position[2] - mean_z).powi(2))
        .sum::<f64>()
        / n;
    if var_z < 1e-30 {
        return f64::INFINITY;
    }
    kbt * projected_area / (8.0 * PI * PI * var_z)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_lipid(z: f64) -> Lipid {
        Lipid::new([0.0, 0.0, z], [0.0, 0.0, 1.0], LipidType::Dppc)
    }

    fn build_test_membrane(n: usize) -> Membrane {
        Membrane::build_bilayer(
            n,
            LipidType::Dppc,
            [10.0, 10.0],
            4.0,
            MembraneProperties::default_phospholipid(),
        )
    }

    // --- LipidType ---

    #[test]
    fn test_lipid_type_reference_area() {
        assert!((LipidType::Dppc.reference_area() - 0.48).abs() < 1e-10);
        assert!((LipidType::Dopc.reference_area() - 0.72).abs() < 1e-10);
        assert!((LipidType::Cholesterol.reference_area() - 0.38).abs() < 1e-10);
    }

    #[test]
    fn test_lipid_type_hydrophobic_length() {
        assert!(LipidType::Dppc.hydrophobic_length() > 0.0);
        assert!(LipidType::Dopc.hydrophobic_length() > 0.0);
    }

    #[test]
    fn test_lipid_type_custom() {
        let lt = LipidType::Custom("POPE".into());
        assert_eq!(lt.reference_area(), 0.60);
    }

    // --- Lipid struct ---

    #[test]
    fn test_lipid_creation_normalises_orientation() {
        let l = Lipid::new([0.0; 3], [0.0, 0.0, 5.0], LipidType::Dppc);
        let mag =
            (l.orientation[0].powi(2) + l.orientation[1].powi(2) + l.orientation[2].powi(2)).sqrt();
        assert!((mag - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lipid_tilt_angle_zero() {
        let l = Lipid::new([0.0; 3], [0.0, 0.0, 1.0], LipidType::Dppc);
        assert!(l.tilt_angle() < 1e-10);
    }

    #[test]
    fn test_lipid_tilt_angle_90() {
        let l = Lipid::new([0.0; 3], [1.0, 0.0, 0.0], LipidType::Dppc);
        assert!((l.tilt_angle() - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lipid_compute_order_parallel() {
        let l = Lipid::new([0.0; 3], [0.0, 0.0, 1.0], LipidType::Dppc);
        // cos θ = 1 → S₂ = 1.5 - 0.5 = 1.0
        assert!((l.compute_order() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lipid_compute_order_perpendicular() {
        let l = Lipid::new([0.0; 3], [1.0, 0.0, 0.0], LipidType::Dppc);
        // cos θ = 0 → S₂ = -0.5
        assert!((l.compute_order() - (-0.5)).abs() < 1e-12);
    }

    #[test]
    fn test_lipid_area_per_lipid_set() {
        let l = Lipid::new([0.0; 3], [0.0, 0.0, 1.0], LipidType::Dopc);
        assert!((l.area_per_lipid - 0.72).abs() < 1e-10);
    }

    // --- Membrane construction ---

    #[test]
    fn test_membrane_total_lipids() {
        let m = build_test_membrane(16);
        assert_eq!(m.total_lipids(), 32);
    }

    #[test]
    fn test_membrane_lipids_per_leaflet() {
        let m = build_test_membrane(9);
        assert_eq!(m.lipids_per_leaflet, 9);
    }

    #[test]
    fn test_membrane_upper_leaflet_positive_z() {
        let m = build_test_membrane(4);
        for l in &m.lipids[..4] {
            assert!(l.position[2] > 0.0, "upper leaflet z = {}", l.position[2]);
        }
    }

    #[test]
    fn test_membrane_lower_leaflet_negative_z() {
        let m = build_test_membrane(4);
        for l in &m.lipids[4..] {
            assert!(l.position[2] < 0.0, "lower leaflet z = {}", l.position[2]);
        }
    }

    #[test]
    fn test_membrane_bilayer_thickness() {
        let m = build_test_membrane(4);
        assert!((m.bilayer_thickness - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_membrane_area() {
        let m = build_test_membrane(4);
        assert!((m.area - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_membrane_step_md_changes_positions() {
        let mut m = build_test_membrane(4);
        let z0 = m.lipids[0].position[0];
        m.step_md(1.0, 0.01, 0.1);
        // After one step positions should differ (noise is deterministic but nonzero)
        let changed = m.lipids.iter().any(|l| (l.position[0] - z0).abs() > 1e-15);
        assert!(changed);
    }

    #[test]
    fn test_membrane_step_md_pbc_x() {
        let mut m = build_test_membrane(4);
        // Force lipid 0 near the boundary
        m.lipids[0].position[0] = 9.9;
        m.step_md(1000.0, 1.0, 10.0);
        for l in &m.lipids {
            assert!(
                l.position[0] >= 0.0 && l.position[0] < m.box_xy[0],
                "PBC violation in x: {}",
                l.position[0]
            );
            assert!(
                l.position[1] >= 0.0 && l.position[1] < m.box_xy[1],
                "PBC violation in y: {}",
                l.position[1]
            );
        }
    }

    #[test]
    fn test_membrane_insert_protein() {
        let mut m = build_test_membrane(4);
        let p = ProteinInsertion::new([5.0, 5.0, 0.0], 0.0, 3.8, m.bilayer_thickness, 0.5);
        m.insert_protein(p);
        assert_eq!(m.proteins.len(), 1);
    }

    #[test]
    fn test_membrane_update_thickness() {
        let mut m = build_test_membrane(4);
        m.update_bilayer_thickness();
        // Upper at +2, lower at -2 → thickness = 4
        assert!(m.bilayer_thickness > 0.0);
    }

    // --- MembraneProperties ---

    #[test]
    fn test_membrane_properties_default() {
        let p = MembraneProperties::default_phospholipid();
        assert!(p.bending_modulus > 0.0);
        assert_eq!(p.tension, 0.0);
        assert!(p.area_compressibility > 0.0);
    }

    // --- ProteinInsertion ---

    #[test]
    fn test_protein_mismatch() {
        let p = ProteinInsertion::new([0.0; 3], 0.0, 4.0, 3.5, 0.5);
        assert!((p.hydrophobic_mismatch - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_protein_deformation_energy_positive() {
        let p = ProteinInsertion::new([0.0; 3], 0.0, 4.0, 3.5, 0.5);
        let e = p.deformation_energy(230.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_protein_deformation_energy_zero_mismatch() {
        let p = ProteinInsertion::new([0.0; 3], 0.0, 3.5, 3.5, 0.5);
        assert_eq!(p.deformation_energy(230.0), 0.0);
    }

    // --- MembraneFluidity ---

    #[test]
    fn test_fluidity_msd() {
        let f = MembraneFluidity::new(1.0, 0.001, 300.0);
        assert!((f.msd(1.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_fluidity_half_time() {
        let f = MembraneFluidity::new(1.0, 1.0, 300.0);
        assert!((f.flip_flop_half_time() - 2_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn test_fluidity_half_time_zero_rate() {
        let f = MembraneFluidity::new(1.0, 0.0, 300.0);
        assert!(f.flip_flop_half_time().is_infinite());
    }

    // --- compute_order_parameter ---

    #[test]
    fn test_order_parameter_empty() {
        assert_eq!(compute_order_parameter(&[]), 0.0);
    }

    #[test]
    fn test_order_parameter_perfect_alignment() {
        let lipids: Vec<Lipid> = (0..10).map(|i| flat_lipid(i as f64)).collect();
        let s = compute_order_parameter(&lipids);
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_order_parameter_perpendicular() {
        let lipids: Vec<Lipid> = (0..10)
            .map(|i| Lipid::new([0.0, 0.0, i as f64], [1.0, 0.0, 0.0], LipidType::Dppc))
            .collect();
        let s = compute_order_parameter(&lipids);
        assert!((s - (-0.5)).abs() < 1e-12);
    }

    #[test]
    fn test_order_parameter_range() {
        let m = build_test_membrane(9);
        let s = compute_order_parameter(&m.lipids);
        assert!((-0.5..=1.0).contains(&s));
    }

    // --- compute_area_per_lipid ---

    #[test]
    fn test_area_per_lipid_zero_leaflet() {
        let m = Membrane {
            lipids: vec![],
            lipids_per_leaflet: 0,
            bilayer_thickness: 4.0,
            area: 100.0,
            box_xy: [10.0, 10.0],
            properties: MembraneProperties::default_phospholipid(),
            proteins: vec![],
            step_count: 0,
        };
        assert_eq!(compute_area_per_lipid(&m), 0.0);
    }

    #[test]
    fn test_area_per_lipid_value() {
        let m = build_test_membrane(25);
        let apl = compute_area_per_lipid(&m);
        // area = 10*10 = 100 nm², n = 25, apl = 4.0
        assert!((apl - 4.0).abs() < 1e-10);
    }

    // --- order_parameter_profile ---

    #[test]
    fn test_order_profile_bins() {
        let lipids: Vec<Lipid> = (0..20).map(|i| flat_lipid(i as f64 * 0.1)).collect();
        let (zc, s) = order_parameter_profile(&lipids, 0.0, 2.0, 4);
        assert_eq!(zc.len(), 4);
        assert_eq!(s.len(), 4);
    }

    #[test]
    fn test_order_profile_empty_bin_zero() {
        let lipids = vec![flat_lipid(0.5)];
        let (_, s) = order_parameter_profile(&lipids, 0.0, 2.0, 4);
        // Only bin 1 is populated; others should be 0.0
        assert_eq!(s[2], 0.0);
    }

    // --- estimate_bending_modulus ---

    #[test]
    fn test_bending_modulus_flat_infinite() {
        // All lipids at same z → var_z = 0 → κ = ∞
        let lipids: Vec<Lipid> = (0..10).map(|_| flat_lipid(0.0)).collect();
        let kappa = estimate_bending_modulus(&lipids, 1.0, 100.0);
        assert!(kappa.is_infinite());
    }

    #[test]
    fn test_bending_modulus_fluctuating() {
        let lipids: Vec<Lipid> = (0..100)
            .map(|i| flat_lipid((i as f64) * 0.01 - 0.5))
            .collect();
        let kappa = estimate_bending_modulus(&lipids, 1.0, 100.0);
        assert!(kappa > 0.0 && kappa.is_finite());
    }

    #[test]
    fn test_bending_modulus_too_few_lipids() {
        let lipids = vec![flat_lipid(0.0)];
        let kappa = estimate_bending_modulus(&lipids, 1.0, 100.0);
        assert_eq!(kappa, 0.0);
    }
}
