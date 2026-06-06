// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Surgical simulation: tissue cutting, needle insertion, organ deformation.
//!
//! Provides models for:
//! - Surgical tools (scalpel, needle, grasper)
//! - Tissue cutting with force model
//! - Needle insertion into layered tissue
//! - Grasper jaw mechanics
//! - Soft tissue deformation under surgical loads
//! - Vessel damage and bleeding model
//! - Laparoscopic trocar port kinematics (remote center of motion)
//! - Robotic surgery (da Vinci) wrist kinematics

/// 3D vector type for surgical simulation \[x, y, z\].
pub type Vec3 = [f64; 3];

/// Surgical tool type classification.
#[derive(Debug, Clone, PartialEq)]
pub enum ToolType {
    /// Scalpel for incision
    Scalpel,
    /// Needle for injection or suturing
    Needle,
    /// Grasper for tissue manipulation
    Grasper,
    /// Electrocautery device
    Cautery,
    /// Stapler for tissue closure
    Stapler,
}

/// Surgical tool with tip position and velocity.
pub struct SurgicalTool {
    /// Tool type
    pub tool_type: ToolType,
    /// Tip position in world space \[m\]
    pub tip_position: Vec3,
    /// Tip velocity \[m/s\]
    pub tip_velocity: Vec3,
    /// Tool shaft direction (unit vector)
    pub shaft_direction: Vec3,
    /// Tool diameter \[m\]
    pub diameter: f64,
    /// Applied force at tip \[N\]
    pub applied_force: Vec3,
}

impl SurgicalTool {
    /// Create a new surgical tool.
    ///
    /// # Arguments
    /// * `tool_type` - Type of surgical tool
    /// * `tip_position` - Initial tip position \[m\]
    /// * `shaft_direction` - Unit vector along shaft
    /// * `diameter` - Tool diameter \[m\]
    pub fn new(
        tool_type: ToolType,
        tip_position: Vec3,
        shaft_direction: Vec3,
        diameter: f64,
    ) -> Self {
        Self {
            tool_type,
            tip_position,
            tip_velocity: [0.0; 3],
            shaft_direction,
            diameter,
            applied_force: [0.0; 3],
        }
    }

    /// Update tip position and velocity.
    pub fn update_position(&mut self, new_position: Vec3, dt: f64) {
        let old = self.tip_position;
        self.tip_velocity = [
            (new_position[0] - old[0]) / dt,
            (new_position[1] - old[1]) / dt,
            (new_position[2] - old[2]) / dt,
        ];
        self.tip_position = new_position;
    }

    /// Tip speed \[m/s\].
    pub fn tip_speed(&self) -> f64 {
        vec3_norm(self.tip_velocity)
    }

    /// Set applied force at tool tip.
    pub fn set_force(&mut self, force: Vec3) {
        self.applied_force = force;
    }
}

/// Tissue cutting model for soft tissue meshes.
///
/// Models edge-tool intersection, remeshing, and cutting force.
pub struct TissueCutting {
    /// Tensile strength of tissue σ_c \[Pa\]
    pub tensile_strength: f64,
    /// Cut depth achieved so far \[m\]
    pub cut_depth: f64,
    /// Total cut area so far \[m²\]
    pub cut_area: f64,
    /// Number of edges split
    pub edges_split: usize,
}

impl TissueCutting {
    /// Create a new tissue cutting model.
    ///
    /// # Arguments
    /// * `tensile_strength` - Tissue tensile strength \[Pa\]
    pub fn new(tensile_strength: f64) -> Self {
        Self {
            tensile_strength,
            cut_depth: 0.0,
            cut_area: 0.0,
            edges_split: 0,
        }
    }

    /// Cutting force model: F = σ_c * A.
    ///
    /// # Arguments
    /// * `cut_area` - Area being cut \[m²\]
    pub fn cutting_force(&self, cut_area: f64) -> f64 {
        self.tensile_strength * cut_area
    }

    /// Test if a tool ray intersects a triangle edge (ray vs. edge test).
    ///
    /// Returns the intersection parameter t ∈ \[0,1\] along the edge, or None.
    ///
    /// # Arguments
    /// * `tool_pos` - Tool position
    /// * `tool_dir` - Tool direction (unit vector)
    /// * `edge_a` - First edge vertex
    /// * `edge_b` - Second edge vertex
    /// * `tolerance` - Distance tolerance for intersection \[m\]
    pub fn ray_edge_intersection(
        tool_pos: Vec3,
        tool_dir: Vec3,
        edge_a: Vec3,
        edge_b: Vec3,
        tolerance: f64,
    ) -> Option<f64> {
        // Parametric: P(t) = edge_a + t*(edge_b - edge_a)
        // Find t where |P(t) - tool_pos - s*tool_dir| is minimized
        let ab = vec3_sub(edge_b, edge_a);
        let ao = vec3_sub(tool_pos, edge_a);
        let ab_len2 = vec3_dot(ab, ab);
        if ab_len2 < 1e-30 {
            return None;
        }
        // Project ao onto ab to find closest point on edge
        let t = vec3_dot(ao, ab) / ab_len2;
        if !(0.0..=1.0).contains(&t) {
            return None;
        }
        let pt = vec3_add(edge_a, vec3_scale(ab, t));
        let diff = vec3_sub(pt, tool_pos);
        // Distance from edge point to tool ray
        let along = vec3_dot(diff, tool_dir);
        let closest = vec3_sub(diff, vec3_scale(tool_dir, along));
        let dist = vec3_norm(closest);
        if dist < tolerance { Some(t) } else { None }
    }

    /// Perform an edge split remesh operation.
    ///
    /// Records the split in the cutting model statistics.
    ///
    /// # Arguments
    /// * `edge_a` - First vertex of edge \[m\]
    /// * `edge_b` - Second vertex of edge \[m\]
    /// * `split_t` - Parametric split position \[0,1\]
    pub fn split_edge(&mut self, edge_a: Vec3, edge_b: Vec3, split_t: f64) -> Vec3 {
        self.edges_split += 1;
        let ab = vec3_sub(edge_b, edge_a);
        vec3_add(edge_a, vec3_scale(ab, split_t))
    }

    /// Simulate cutting a region of given area and update internal state.
    pub fn perform_cut(&mut self, area: f64, depth: f64) {
        self.cut_area += area;
        self.cut_depth += depth;
    }
}

/// Needle insertion model for layered tissue.
///
/// Models puncture force, stiffness along path, and tissue displacement.
pub struct NeedleInsertion {
    /// Needle tip diameter \[m\]
    pub tip_diameter: f64,
    /// Skin puncture strength σ_p \[Pa\]
    pub puncture_strength: f64,
    /// Tissue layers (stiffness \[Pa/m\], thickness \[m\])
    pub layers: Vec<(f64, f64)>,
    /// Current insertion depth \[m\]
    pub depth: f64,
    /// Whether needle has punctured initial skin
    pub has_punctured: bool,
}

impl NeedleInsertion {
    /// Create needle insertion model.
    ///
    /// # Arguments
    /// * `tip_diameter` - Needle tip diameter \[m\]
    /// * `puncture_strength` - Skin puncture strength \[Pa\]
    pub fn new(tip_diameter: f64, puncture_strength: f64) -> Self {
        Self {
            tip_diameter,
            puncture_strength,
            layers: Vec::new(),
            depth: 0.0,
            has_punctured: false,
        }
    }

    /// Add a tissue layer.
    ///
    /// # Arguments
    /// * `stiffness` - Layer stiffness \[Pa/m\]
    /// * `thickness` - Layer thickness \[m\]
    pub fn add_layer(&mut self, stiffness: f64, thickness: f64) {
        self.layers.push((stiffness, thickness));
    }

    /// Skin puncture force: F_p = σ_p * π * (d_tip/2)².
    pub fn puncture_force(&self) -> f64 {
        let r = self.tip_diameter / 2.0;
        self.puncture_strength * std::f64::consts::PI * r * r
    }

    /// Insertion stiffness at current depth.
    ///
    /// Returns the spring constant of the layer currently being traversed.
    pub fn stiffness_at_depth(&self, depth: f64) -> f64 {
        let mut accumulated = 0.0;
        for &(k, t) in &self.layers {
            accumulated += t;
            if depth <= accumulated {
                return k;
            }
        }
        // Beyond all layers: return last layer stiffness or zero
        self.layers.last().map(|&(k, _)| k).unwrap_or(0.0)
    }

    /// Insertion force at current depth (force = stiffness × depth within layer).
    pub fn insertion_force(&self, depth: f64) -> f64 {
        if depth <= 0.0 {
            return 0.0;
        }
        if !self.has_punctured && depth > 0.0 {
            return self.puncture_force();
        }
        let k = self.stiffness_at_depth(depth);
        k * depth
    }

    /// Advance needle insertion by given increment.
    pub fn advance(&mut self, delta: f64) {
        if !self.has_punctured {
            self.has_punctured = true;
        }
        self.depth += delta;
    }

    /// Tissue push-through displacement (fraction of needle diameter).
    pub fn tissue_displacement(&self, depth: f64) -> f64 {
        0.1 * self.tip_diameter * (depth / 0.01).min(1.0)
    }
}

/// Grasper jaw model for tissue holding.
///
/// Models grip force vs. jaw angle/displacement.
pub struct GrasperModel {
    /// Jaw spring stiffness k \[N/m\]
    pub jaw_stiffness: f64,
    /// Maximum jaw opening \[m\]
    pub max_opening: f64,
    /// Current jaw displacement from closed position \[m\]
    pub jaw_displacement: f64,
    /// Friction coefficient between jaw and tissue
    pub friction_coeff: f64,
}

impl GrasperModel {
    /// Create a grasper model.
    ///
    /// # Arguments
    /// * `jaw_stiffness` - Jaw spring stiffness \[N/m\]
    /// * `max_opening` - Maximum jaw opening \[m\]
    /// * `friction_coeff` - Jaw-tissue friction coefficient
    pub fn new(jaw_stiffness: f64, max_opening: f64, friction_coeff: f64) -> Self {
        Self {
            jaw_stiffness,
            max_opening,
            jaw_displacement: 0.0,
            friction_coeff,
        }
    }

    /// Grip force: F = k * Δ_jaw.
    pub fn grip_force(&self) -> f64 {
        self.jaw_stiffness * self.jaw_displacement
    }

    /// Maximum grip force (at full closure from max_opening).
    pub fn max_grip_force(&self) -> f64 {
        self.jaw_stiffness * self.max_opening
    }

    /// Check if tissue can be held: μ * N ≥ F_required.
    ///
    /// # Arguments
    /// * `required_force` - Force required to hold tissue \[N\]
    pub fn can_hold(&self, required_force: f64) -> bool {
        let normal_force = self.grip_force();
        self.friction_coeff * normal_force >= required_force
    }

    /// Set jaw displacement (clamped to \[0, max_opening\]).
    pub fn set_displacement(&mut self, displacement: f64) {
        self.jaw_displacement = displacement.clamp(0.0, self.max_opening);
    }

    /// Jaw angle from displacement (assuming half-jaw arm length L).
    ///
    /// θ = arcsin(Δ / (2*L))
    pub fn jaw_angle(&self, jaw_arm_length: f64) -> f64 {
        (self.jaw_displacement / (2.0 * jaw_arm_length)).asin()
    }
}

/// Soft tissue deformation under surgical loads.
///
/// Linear viscoelastic: σ = E*ε + η*ε̇
pub struct TissueDeformation {
    /// Young's modulus E \[Pa\]
    pub elastic_modulus: f64,
    /// Viscosity η \[Pa·s\]
    pub viscosity: f64,
    /// Current strain ε
    pub strain: f64,
    /// Current strain rate \[1/s\]
    pub strain_rate: f64,
}

impl TissueDeformation {
    /// Create tissue deformation model.
    ///
    /// # Arguments
    /// * `elastic_modulus` - Young's modulus E \[Pa\]
    /// * `viscosity` - Viscosity η \[Pa·s\]
    pub fn new(elastic_modulus: f64, viscosity: f64) -> Self {
        Self {
            elastic_modulus,
            viscosity,
            strain: 0.0,
            strain_rate: 0.0,
        }
    }

    /// Stress: σ = E*ε + η*ε̇.
    pub fn stress(&self) -> f64 {
        self.elastic_modulus * self.strain + self.viscosity * self.strain_rate
    }

    /// Deformation under applied stress (quasi-static: ε = σ/E).
    pub fn deformation_from_stress(&self, applied_stress: f64) -> f64 {
        applied_stress / self.elastic_modulus
    }

    /// Apply load and update strain.
    ///
    /// # Arguments
    /// * `applied_stress` - Applied stress \[Pa\]
    /// * `dt` - Time step \[s\]
    pub fn apply_load(&mut self, applied_stress: f64, dt: f64) {
        let new_strain = applied_stress / self.elastic_modulus;
        self.strain_rate = (new_strain - self.strain) / dt;
        self.strain = new_strain;
    }

    /// Tissue stiffness under compression is higher (tissue nonlinearity factor).
    pub fn effective_stiffness(&self, compression_factor: f64) -> f64 {
        self.elastic_modulus * (1.0 + compression_factor * self.strain.abs())
    }
}

/// Vessel damage and bleeding model.
///
/// Pressure-threshold model for vessel rupture and bleeding rate.
pub struct BleedingModel {
    /// Vessel rupture pressure threshold \[Pa\]
    pub rupture_threshold: f64,
    /// Vessel cross-sectional area \[m²\]
    pub vessel_area: f64,
    /// Blood viscosity \[Pa·s\]
    pub blood_viscosity: f64,
    /// Vessel resistance \[Pa·s/m³\]
    pub vessel_resistance: f64,
    /// Whether vessel is ruptured
    pub is_ruptured: bool,
}

impl BleedingModel {
    /// Create bleeding model.
    ///
    /// # Arguments
    /// * `rupture_threshold` - Pressure to rupture vessel \[Pa\]
    /// * `vessel_diameter` - Vessel diameter \[m\]
    /// * `blood_viscosity` - Blood dynamic viscosity \[Pa·s\]
    pub fn new(rupture_threshold: f64, vessel_diameter: f64, blood_viscosity: f64) -> Self {
        let r = vessel_diameter / 2.0;
        let vessel_area = std::f64::consts::PI * r * r;
        // Simplified Poiseuille resistance for unit length
        let vessel_resistance = 8.0 * blood_viscosity / (std::f64::consts::PI * r.powi(4));
        Self {
            rupture_threshold,
            vessel_area,
            blood_viscosity,
            vessel_resistance,
            is_ruptured: false,
        }
    }

    /// Check if pressure causes rupture.
    pub fn check_rupture(&mut self, pressure: f64) -> bool {
        if pressure >= self.rupture_threshold {
            self.is_ruptured = true;
        }
        self.is_ruptured
    }

    /// Blood flow rate from damaged vessel: Q = ΔP / R (Hagen-Poiseuille).
    ///
    /// Returns 0 if vessel is not ruptured.
    pub fn flow_rate(&self, pressure_drop: f64) -> f64 {
        if !self.is_ruptured {
            return 0.0;
        }
        (pressure_drop / self.vessel_resistance).max(0.0)
    }

    /// Blood loss volume over time interval.
    ///
    /// # Arguments
    /// * `pressure_drop` - Driving pressure \[Pa\]
    /// * `time` - Duration \[s\]
    pub fn blood_loss(&self, pressure_drop: f64, time: f64) -> f64 {
        self.flow_rate(pressure_drop) * time
    }

    /// Repair vessel (stop bleeding).
    pub fn repair(&mut self) {
        self.is_ruptured = false;
    }
}

/// Laparoscopic trocar port kinematics model.
///
/// Implements the remote center of motion (RCM) constraint at the trocar point.
/// The tool rotates about the fulcrum (trocar insertion point).
pub struct LaparoscopicTool {
    /// Trocar insertion point (fulcrum) position \[m\]
    pub trocar_point: Vec3,
    /// Tool tip position \[m\]
    pub tip_position: Vec3,
    /// Tool length from trocar to tip \[m\]
    pub tool_length: f64,
    /// Current azimuth angle \[rad\]
    pub azimuth: f64,
    /// Current elevation angle \[rad\]
    pub elevation: f64,
}

impl LaparoscopicTool {
    /// Create laparoscopic tool with trocar constraint.
    ///
    /// # Arguments
    /// * `trocar_point` - Fulcrum point (trocar insertion) \[m\]
    /// * `tool_length` - Length of instrument from trocar to tip \[m\]
    pub fn new(trocar_point: Vec3, tool_length: f64) -> Self {
        let tip = [
            trocar_point[0],
            trocar_point[1] - tool_length,
            trocar_point[2],
        ];
        Self {
            trocar_point,
            tip_position: tip,
            tool_length,
            azimuth: 0.0,
            elevation: -std::f64::consts::FRAC_PI_2,
        }
    }

    /// Update tip position based on azimuth and elevation angles.
    ///
    /// The tool rotates about the trocar point (remote center of motion).
    pub fn set_angles(&mut self, azimuth: f64, elevation: f64) {
        self.azimuth = azimuth;
        self.elevation = elevation;
        // Compute tip position from spherical coordinates about trocar
        let dx = self.tool_length * elevation.cos() * azimuth.sin();
        let dy = self.tool_length * elevation.sin();
        let dz = self.tool_length * elevation.cos() * azimuth.cos();
        self.tip_position = [
            self.trocar_point[0] + dx,
            self.trocar_point[1] + dy,
            self.trocar_point[2] + dz,
        ];
    }

    /// Verify RCM constraint: distance from trocar to tip = tool_length.
    pub fn verify_rcm_constraint(&self) -> f64 {
        let d = vec3_sub(self.tip_position, self.trocar_point);
        vec3_norm(d)
    }

    /// Is the RCM constraint satisfied within tolerance?
    pub fn rcm_satisfied(&self, tolerance: f64) -> bool {
        (self.verify_rcm_constraint() - self.tool_length).abs() < tolerance
    }

    /// Shaft direction (unit vector from trocar to tip).
    pub fn shaft_direction(&self) -> Vec3 {
        let d = vec3_sub(self.tip_position, self.trocar_point);
        let len = vec3_norm(d);
        if len < 1e-15 {
            [0.0, 1.0, 0.0]
        } else {
            [d[0] / len, d[1] / len, d[2] / len]
        }
    }
}

/// Robotic surgery (da Vinci inspired) wrist kinematic model.
///
/// Models a 2-DOF wrist at the end of the instrument shaft.
pub struct RoboticSurgery {
    /// Instrument shaft tip position \[m\]
    pub shaft_tip: Vec3,
    /// Shaft direction unit vector
    pub shaft_direction: Vec3,
    /// Wrist pitch angle \[rad\]
    pub wrist_pitch: f64,
    /// Wrist yaw angle \[rad\]
    pub wrist_yaw: f64,
    /// Length of wrist link \[m\]
    pub wrist_length: f64,
    /// Maximum wrist angle \[rad\]
    pub max_wrist_angle: f64,
}

impl RoboticSurgery {
    /// Create robotic surgery instrument model.
    ///
    /// # Arguments
    /// * `shaft_tip` - Position of shaft tip (proximal wrist) \[m\]
    /// * `shaft_direction` - Unit vector along shaft
    /// * `wrist_length` - Length of wrist segment \[m\]
    /// * `max_wrist_angle` - Maximum wrist deflection angle \[rad\]
    pub fn new(
        shaft_tip: Vec3,
        shaft_direction: Vec3,
        wrist_length: f64,
        max_wrist_angle: f64,
    ) -> Self {
        Self {
            shaft_tip,
            shaft_direction,
            wrist_pitch: 0.0,
            wrist_yaw: 0.0,
            wrist_length,
            max_wrist_angle,
        }
    }

    /// Set wrist angles (clamped to ±max_wrist_angle).
    pub fn set_wrist(&mut self, pitch: f64, yaw: f64) {
        self.wrist_pitch = pitch.clamp(-self.max_wrist_angle, self.max_wrist_angle);
        self.wrist_yaw = yaw.clamp(-self.max_wrist_angle, self.max_wrist_angle);
    }

    /// Compute instrument tip position given wrist angles.
    ///
    /// Simple 2D model: pitch rotates in XY plane, yaw rotates in XZ plane.
    pub fn instrument_tip(&self) -> Vec3 {
        // Build local rotation: combine pitch (around Z) and yaw (around Y)
        let cp = self.wrist_pitch.cos();
        let sp = self.wrist_pitch.sin();
        let cy = self.wrist_yaw.cos();
        let sy = self.wrist_yaw.yaw_factor();

        // Simplified: tip offset in shaft frame
        let local_dir = [sy * cp, sp, cy * cp];

        // Rotate local direction by shaft direction
        let tip_offset = vec3_scale(local_dir, self.wrist_length);
        // For simplicity, use shaft as Z axis
        [
            self.shaft_tip[0]
                + self.shaft_direction[0] * self.wrist_length * cp * cy
                + tip_offset[0] * 0.1,
            self.shaft_tip[1]
                + self.shaft_direction[1] * self.wrist_length * cp * cy
                + tip_offset[1] * 0.1,
            self.shaft_tip[2]
                + self.shaft_direction[2] * self.wrist_length * cp * cy
                + tip_offset[2] * 0.1,
        ]
    }

    /// Distance from shaft tip to instrument tip.
    pub fn reach(&self) -> f64 {
        let tip = self.instrument_tip();
        vec3_norm(vec3_sub(tip, self.shaft_tip))
    }

    /// Check if wrist angles are within workspace.
    pub fn in_workspace(&self) -> bool {
        self.wrist_pitch.abs() <= self.max_wrist_angle
            && self.wrist_yaw.abs() <= self.max_wrist_angle
    }
}

/// Extension trait for yaw factor computation.
trait YawFactor {
    /// Returns sin of self (used for yaw computation).
    fn yaw_factor(self) -> f64;
}

impl YawFactor for f64 {
    fn yaw_factor(self) -> f64 {
        self.sin()
    }
}

// -------------------------------------------------------------------------
// Vector utility functions
// -------------------------------------------------------------------------

/// Add two Vec3.
fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract b from a.
fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale Vec3 by scalar.
fn vec3_scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm.
fn vec3_norm(a: Vec3) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Cutting force proportional to cut area
    #[test]
    fn test_cutting_force_proportional_to_area() {
        let tc = TissueCutting::new(1e6);
        let f1 = tc.cutting_force(1e-6);
        let f2 = tc.cutting_force(2e-6);
        assert!(
            (f2 / f1 - 2.0).abs() < EPS,
            "Cutting force should be proportional to area: ratio={:.6}",
            f2 / f1
        );
    }

    // 2. Cutting force = σ_c * A
    #[test]
    fn test_cutting_force_formula() {
        let sigma_c = 5e5;
        let area = 3e-6;
        let tc = TissueCutting::new(sigma_c);
        let f = tc.cutting_force(area);
        assert!(
            (f - sigma_c * area).abs() < EPS,
            "Cutting force should be σ_c*A"
        );
    }

    // 3. Ray-edge intersection detects hit
    #[test]
    fn test_ray_edge_intersection_hit() {
        let tool_pos = [0.0, 0.0, 0.0];
        let tool_dir = [1.0, 0.0, 0.0];
        let edge_a = [1.0, -0.001, 0.0];
        let edge_b = [1.0, 0.001, 0.0];
        let result = TissueCutting::ray_edge_intersection(tool_pos, tool_dir, edge_a, edge_b, 0.01);
        assert!(result.is_some(), "Ray should intersect edge");
    }

    // 4. Ray-edge intersection misses non-intersecting edge
    #[test]
    fn test_ray_edge_intersection_miss() {
        let tool_pos = [0.0, 0.0, 0.0];
        let tool_dir = [1.0, 0.0, 0.0];
        let edge_a = [1.0, 1.0, 0.0];
        let edge_b = [1.0, 2.0, 0.0];
        let result = TissueCutting::ray_edge_intersection(tool_pos, tool_dir, edge_a, edge_b, 0.01);
        assert!(result.is_none(), "Ray should miss distant edge");
    }

    // 5. Edge split generates midpoint
    #[test]
    fn test_edge_split_midpoint() {
        let mut tc = TissueCutting::new(1e6);
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let mid = tc.split_edge(a, b, 0.5);
        assert!(
            (mid[0] - 1.0).abs() < EPS && mid[1].abs() < EPS && mid[2].abs() < EPS,
            "Split midpoint should be [1,0,0], got {:?}",
            mid
        );
    }

    // 6. Edges split count increments
    #[test]
    fn test_edge_split_count() {
        let mut tc = TissueCutting::new(1e6);
        tc.split_edge([0.0; 3], [1.0, 0.0, 0.0], 0.5);
        tc.split_edge([0.0; 3], [1.0, 0.0, 0.0], 0.3);
        assert_eq!(tc.edges_split, 2, "Should have split 2 edges");
    }

    // 7. Needle: puncture force = σ_p * π * r²
    #[test]
    fn test_needle_puncture_force() {
        let sigma_p = 1e5;
        let d = 1e-3;
        let ni = NeedleInsertion::new(d, sigma_p);
        let fp = ni.puncture_force();
        let expected = sigma_p * std::f64::consts::PI * (d / 2.0).powi(2);
        assert!(
            (fp - expected).abs() < 1e-12,
            "Puncture force formula mismatch: {:.6e} vs {:.6e}",
            fp,
            expected
        );
    }

    // 8. Needle: insertion force increases with depth (after puncture)
    #[test]
    fn test_needle_force_increases_with_depth() {
        let mut ni = NeedleInsertion::new(1e-3, 1e5);
        ni.add_layer(1e4, 0.05);
        ni.has_punctured = true;
        let f1 = ni.insertion_force(0.01);
        let f2 = ni.insertion_force(0.02);
        assert!(f2 > f1, "Insertion force should increase with depth");
    }

    // 9. Needle: stiffness at zero depth is first layer stiffness
    #[test]
    fn test_needle_stiffness_at_zero_depth() {
        let mut ni = NeedleInsertion::new(1e-3, 1e5);
        ni.add_layer(2e4, 0.05);
        let k = ni.stiffness_at_depth(0.01);
        assert!(
            (k - 2e4).abs() < EPS,
            "Stiffness at small depth should be first layer k"
        );
    }

    // 10. Needle: advance updates depth
    #[test]
    fn test_needle_advance_updates_depth() {
        let mut ni = NeedleInsertion::new(1e-3, 1e5);
        ni.advance(0.01);
        ni.advance(0.005);
        assert!(
            (ni.depth - 0.015).abs() < EPS,
            "Depth should accumulate correctly"
        );
    }

    // 11. Grasper: grip force proportional to jaw displacement
    #[test]
    fn test_grasper_force_proportional_to_displacement() {
        let grasper = GrasperModel::new(100.0, 0.02, 0.5);
        let mut g1 = GrasperModel::new(100.0, 0.02, 0.5);
        let mut g2 = GrasperModel::new(100.0, 0.02, 0.5);
        g1.set_displacement(0.005);
        g2.set_displacement(0.010);
        let f1 = g1.grip_force();
        let f2 = g2.grip_force();
        assert!(
            (f2 / f1 - 2.0).abs() < EPS,
            "Grip force should be proportional to jaw displacement: ratio={:.6}",
            f2 / f1
        );
        let _ = grasper; // suppress unused warning
    }

    // 12. Grasper: grip force = k * Δ_jaw
    #[test]
    fn test_grasper_grip_force_formula() {
        let mut g = GrasperModel::new(200.0, 0.03, 0.4);
        g.set_displacement(0.01);
        let f = g.grip_force();
        assert!((f - 200.0 * 0.01).abs() < EPS, "Grip force = k*Δ");
    }

    // 13. Grasper: can_hold returns true when friction sufficient
    #[test]
    fn test_grasper_can_hold_true() {
        let mut g = GrasperModel::new(500.0, 0.05, 0.6);
        g.set_displacement(0.02); // F = 500*0.02 = 10 N, μ*N = 0.6*10 = 6 N
        assert!(g.can_hold(5.0), "Should be able to hold 5 N");
    }

    // 14. Grasper: can_hold returns false when friction insufficient
    #[test]
    fn test_grasper_can_hold_false() {
        let mut g = GrasperModel::new(100.0, 0.05, 0.2);
        g.set_displacement(0.001); // F = 0.1 N, μ*N = 0.02 N
        assert!(
            !g.can_hold(5.0),
            "Should not hold 5 N with insufficient grip"
        );
    }

    // 15. Grasper: displacement clamped to max_opening
    #[test]
    fn test_grasper_displacement_clamped() {
        let mut g = GrasperModel::new(100.0, 0.02, 0.5);
        g.set_displacement(0.1); // exceeds max_opening
        assert!(
            (g.jaw_displacement - 0.02).abs() < EPS,
            "Displacement should be clamped to max_opening"
        );
    }

    // 16. Trocar: RCM constraint satisfied after set_angles
    #[test]
    fn test_trocar_rcm_constraint() {
        let mut tool = LaparoscopicTool::new([0.0, 0.0, 0.0], 0.3);
        tool.set_angles(0.5, -0.5);
        assert!(
            tool.rcm_satisfied(1e-10),
            "RCM constraint should be satisfied: dist={:.6}",
            tool.verify_rcm_constraint()
        );
    }

    // 17. Trocar: RCM constraint satisfied at different angles
    #[test]
    fn test_trocar_rcm_multiple_angles() {
        let mut tool = LaparoscopicTool::new([1.0, 2.0, 3.0], 0.25);
        for &az in &[0.0, 0.3, -0.3, 1.0] {
            for &el in &[-0.5, -0.3, 0.0, 0.3] {
                tool.set_angles(az, el);
                assert!(
                    tool.rcm_satisfied(1e-10),
                    "RCM violated at az={az:.3}, el={el:.3}: dist={:.6}",
                    tool.verify_rcm_constraint()
                );
            }
        }
    }

    // 18. Trocar: shaft direction is unit vector
    #[test]
    fn test_trocar_shaft_direction_unit() {
        let mut tool = LaparoscopicTool::new([0.0, 0.0, 0.0], 0.3);
        tool.set_angles(0.4, -0.6);
        let dir = tool.shaft_direction();
        let len = vec3_norm(dir);
        assert!(
            (len - 1.0).abs() < 1e-10,
            "Shaft direction should be unit vector"
        );
    }

    // 19. Tissue deformation proportional to load
    #[test]
    fn test_tissue_deformation_proportional_to_load() {
        let tissue = TissueDeformation::new(1e4, 100.0);
        let d1 = tissue.deformation_from_stress(100.0);
        let d2 = tissue.deformation_from_stress(200.0);
        assert!(
            (d2 / d1 - 2.0).abs() < EPS,
            "Deformation should be proportional to stress"
        );
    }

    // 20. Tissue deformation: σ = E*ε at steady state
    #[test]
    fn test_tissue_stress_at_steady_state() {
        let mut tissue = TissueDeformation::new(1e4, 100.0);
        tissue.apply_load(500.0, 0.01);
        tissue.apply_load(500.0, 0.01); // apply again to reach steady state
        let s = tissue.stress();
        assert!(
            s > 0.0,
            "Tissue stress should be positive under load, got {:.6}",
            s
        );
    }

    // 21. Bleeding: no flow below rupture threshold
    #[test]
    fn test_bleeding_no_flow_below_threshold() {
        let model = BleedingModel::new(1e4, 1e-3, 3e-3);
        let flow = model.flow_rate(5000.0); // below 10000 threshold
        assert!(flow.abs() < EPS, "No bleeding below rupture threshold");
    }

    // 22. Bleeding: rupture occurs above threshold
    #[test]
    fn test_bleeding_rupture_above_threshold() {
        let mut model = BleedingModel::new(1e4, 1e-3, 3e-3);
        model.check_rupture(1.5e4);
        assert!(model.is_ruptured, "Vessel should rupture above threshold");
    }

    // 23. Bleeding: flow > 0 after rupture
    #[test]
    fn test_bleeding_flow_after_rupture() {
        let mut model = BleedingModel::new(1e4, 1e-3, 3e-3);
        model.check_rupture(2e4);
        let flow = model.flow_rate(1000.0);
        assert!(flow > 0.0, "Blood should flow after rupture");
    }

    // 24. Bleeding: repair stops bleeding
    #[test]
    fn test_bleeding_repair() {
        let mut model = BleedingModel::new(1e4, 1e-3, 3e-3);
        model.check_rupture(2e4);
        model.repair();
        let flow = model.flow_rate(1000.0);
        assert!(flow.abs() < EPS, "No flow after repair, got {:.6e}", flow);
    }

    // 25. Bleeding: blood loss = flow_rate × time
    #[test]
    fn test_bleeding_blood_loss() {
        let mut model = BleedingModel::new(1e4, 2e-3, 3e-3);
        model.check_rupture(2e4);
        let t = 5.0;
        let dp = 1000.0;
        let volume = model.blood_loss(dp, t);
        let expected = model.flow_rate(dp) * t;
        assert!(
            (volume - expected).abs() < EPS,
            "Blood loss should be flow_rate * time"
        );
    }

    // 26. Tool: tip speed is zero initially
    #[test]
    fn test_tool_tip_speed_zero_initially() {
        let tool = SurgicalTool::new(ToolType::Scalpel, [0.0; 3], [0.0, 0.0, 1.0], 1e-3);
        assert!(
            tool.tip_speed().abs() < EPS,
            "Initial tip speed should be zero"
        );
    }

    // 27. Tool: update position updates velocity
    #[test]
    fn test_tool_update_position_velocity() {
        let mut tool = SurgicalTool::new(ToolType::Needle, [0.0; 3], [0.0, 0.0, 1.0], 5e-4);
        tool.update_position([0.01, 0.0, 0.0], 0.1);
        assert!(
            (tool.tip_velocity[0] - 0.1).abs() < 1e-10,
            "Velocity should be Δx/dt = 0.1 m/s"
        );
    }

    // 28. Robotic surgery: wrist angles clamped
    #[test]
    fn test_robotic_wrist_angles_clamped() {
        let mut robot =
            RoboticSurgery::new([0.0; 3], [0.0, 0.0, 1.0], 0.05, std::f64::consts::FRAC_PI_4);
        robot.set_wrist(2.0, -2.0); // exceeds π/4
        assert!(
            robot.wrist_pitch.abs() <= std::f64::consts::FRAC_PI_4 + EPS,
            "Wrist pitch should be clamped"
        );
        assert!(
            robot.wrist_yaw.abs() <= std::f64::consts::FRAC_PI_4 + EPS,
            "Wrist yaw should be clamped"
        );
    }

    // 29. Robotic surgery: reach is nonzero at neutral pose
    #[test]
    fn test_robotic_reach_nonzero() {
        let robot =
            RoboticSurgery::new([0.0; 3], [0.0, 0.0, 1.0], 0.05, std::f64::consts::FRAC_PI_4);
        let r = robot.reach();
        assert!(r > 0.0, "Reach should be nonzero");
    }

    // 30. Robotic surgery: in_workspace at zero angles
    #[test]
    fn test_robotic_in_workspace_zero() {
        let robot =
            RoboticSurgery::new([0.0; 3], [0.0, 0.0, 1.0], 0.05, std::f64::consts::FRAC_PI_4);
        assert!(robot.in_workspace(), "Zero angles should be in workspace");
    }

    // 31. Tissue effective stiffness increases with compression
    #[test]
    fn test_tissue_effective_stiffness_increases() {
        let mut tissue = TissueDeformation::new(1e4, 100.0);
        tissue.apply_load(1000.0, 0.01);
        let k0 = tissue.effective_stiffness(0.0);
        let k1 = tissue.effective_stiffness(1.0);
        assert!(
            k1 >= k0,
            "Effective stiffness should increase with compression"
        );
    }

    // 32. Needle tissue displacement increases with depth
    #[test]
    fn test_needle_tissue_displacement_increases() {
        let ni = NeedleInsertion::new(1e-3, 1e5);
        let d1 = ni.tissue_displacement(0.005);
        let d2 = ni.tissue_displacement(0.01);
        assert!(d2 >= d1, "Tissue displacement should increase with depth");
    }

    // 33. Cutting perform_cut accumulates area and depth
    #[test]
    fn test_cutting_accumulates() {
        let mut tc = TissueCutting::new(1e6);
        tc.perform_cut(1e-6, 1e-3);
        tc.perform_cut(2e-6, 2e-3);
        assert!(
            (tc.cut_area - 3e-6).abs() < EPS,
            "Cut area should accumulate"
        );
        assert!(
            (tc.cut_depth - 3e-3).abs() < EPS,
            "Cut depth should accumulate"
        );
    }

    // 34. Grasper jaw angle increases with displacement
    #[test]
    fn test_grasper_jaw_angle_increases() {
        let mut g = GrasperModel::new(100.0, 0.05, 0.5);
        g.set_displacement(0.005);
        let angle1 = g.jaw_angle(0.05);
        g.set_displacement(0.010);
        let angle2 = g.jaw_angle(0.05);
        assert!(
            angle2 > angle1,
            "Jaw angle should increase with displacement"
        );
    }

    // 35. Vessel bleeding only above pressure threshold
    #[test]
    fn test_bleeding_threshold_boundary() {
        let threshold = 8000.0;
        let mut model = BleedingModel::new(threshold, 1e-3, 3e-3);
        model.check_rupture(threshold - 1.0);
        assert!(!model.is_ruptured, "Should not rupture below threshold");
        model.check_rupture(threshold);
        assert!(model.is_ruptured, "Should rupture at threshold");
    }
}
