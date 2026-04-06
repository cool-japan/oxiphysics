// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin shell wrinkling mechanics.
//!
//! Implements wrinkling state detection, tension-field theory, wrinkle geometry
//! (wavelength / amplitude), sheet draping under gravity, and critical buckling
//! stress for thin membranes.
//!
//! All quantities are in SI units (Pa, m, N/m …) unless otherwise stated.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// ThinShell
// ---------------------------------------------------------------------------

/// Elastic thin shell (Kirchhoff-Love plate) descriptor.
#[derive(Clone, Debug)]
pub struct ThinShell {
    /// Shell thickness (m).
    pub thickness: f64,
    /// Young's modulus (Pa).
    pub elastic_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
}

impl ThinShell {
    /// Construct a new thin shell.
    pub fn new(thickness: f64, elastic_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            thickness,
            elastic_modulus,
            poisson_ratio,
        }
    }

    /// Bending stiffness D = E t³ / (12 (1 − ν²))  \[N·m\].
    pub fn bending_stiffness(&self) -> f64 {
        let t = self.thickness;
        self.elastic_modulus * t * t * t / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }

    /// Membrane (extensional) stiffness N = E t  \[N/m\].
    pub fn membrane_stiffness(&self) -> f64 {
        self.elastic_modulus * self.thickness
    }
}

// ---------------------------------------------------------------------------
// WrinklingState
// ---------------------------------------------------------------------------

/// Wrinkling state of a thin membrane under in-plane loads.
#[derive(Clone, Debug, PartialEq)]
pub enum WrinklingState {
    /// Both principal stresses are tensile: no wrinkling.
    Taut,
    /// Both principal stresses are compressive: sheet is slack/buckled.
    Slack,
    /// One principal stress tensile, one compressive: wrinkles form.
    Wrinkled,
}

// ---------------------------------------------------------------------------
// WrinklingModel
// ---------------------------------------------------------------------------

/// Wrinkling model for a thin shell element under prescribed principal strains.
#[derive(Clone, Debug)]
pub struct WrinklingModel {
    /// The underlying thin shell.
    pub shell: ThinShell,
    /// Principal strains \[ε₁, ε₂\] where ε₁ ≥ ε₂ by convention.
    pub principal_strains: [f64; 2],
}

impl WrinklingModel {
    /// Construct a new wrinkling model.
    pub fn new(shell: ThinShell, principal_strains: [f64; 2]) -> Self {
        Self {
            shell,
            principal_strains,
        }
    }

    /// Principal stresses from plane-stress constitutive law.
    ///
    /// σ₁ = E/(1−ν²) (ε₁ + ν ε₂),  σ₂ = E/(1−ν²) (ε₂ + ν ε₁)
    fn principal_stresses(&self) -> [f64; 2] {
        let e = self.shell.elastic_modulus;
        let nu = self.shell.poisson_ratio;
        let fac = e / (1.0 - nu * nu);
        let e1 = self.principal_strains[0];
        let e2 = self.principal_strains[1];
        [fac * (e1 + nu * e2), fac * (e2 + nu * e1)]
    }

    /// Determine the wrinkling state based on principal stresses.
    pub fn determine_state(&self) -> WrinklingState {
        let [s1, s2] = self.principal_stresses();
        if s1 >= 0.0 && s2 >= 0.0 {
            WrinklingState::Taut
        } else if s1 <= 0.0 && s2 <= 0.0 {
            WrinklingState::Slack
        } else {
            WrinklingState::Wrinkled
        }
    }

    /// Wrinkle half-wavelength (m) based on bending vs membrane stiffness ratio.
    ///
    /// λ/2 = π (D / N)^(1/2) where D is bending stiffness and N membrane stiffness.
    pub fn wrinkle_wavelength(&self) -> f64 {
        let d = self.shell.bending_stiffness();
        let n = self.shell.membrane_stiffness();
        PI * (d / n).sqrt()
    }

    /// Wrinkle amplitude (m) for a given tension N_tension (N/m) along the
    /// wrinkle direction.
    ///
    /// A = t / sqrt(1 + 4 D N_tension / (π² D_max²))
    /// Simplified practical formula: A ≈ t * sqrt(N_crit / N_tension).
    pub fn wrinkle_amplitude(&self, tension: f64) -> f64 {
        if tension <= 0.0 {
            return 0.0;
        }
        let d = self.shell.bending_stiffness();
        let n_crit = PI * PI * d / (self.wrinkle_wavelength().powi(2));
        self.shell.thickness * (n_crit / tension).sqrt().min(1.0)
    }
}

// ---------------------------------------------------------------------------
// TensionField
// ---------------------------------------------------------------------------

/// 2-D in-plane stress tensor stored column-major as \[σ₁₁, σ₂₁, σ₁₂, σ₂₂\].
///
/// Layout: `[[0,2\],[1,3]]` → `tensor[i + 2*j]` for row i, col j.
#[derive(Clone, Debug)]
pub struct TensionField {
    /// Stress tensor components \[σ₁₁, σ₂₁, σ₁₂, σ₂₂\] (Pa).
    pub stress_tensor: [f64; 4],
}

impl TensionField {
    /// Construct a tension field from a 2×2 stress tensor.
    pub fn new(stress_tensor: [f64; 4]) -> Self {
        Self { stress_tensor }
    }

    /// Isotropic tension field at pressure `p` (Pa): σ = p I₂.
    pub fn isotropic(p: f64) -> Self {
        Self::new([p, 0.0, 0.0, p])
    }

    /// Compute the two principal stresses via eigenvalues of the 2×2 tensor.
    ///
    /// For a symmetric tensor `[[a,b\],[b,d]]`:
    /// `σ_{1,2} = (a+d)/2 ± sqrt(((a-d)/2)² + b²)`
    pub fn compute_principal_stresses(&self) -> [f64; 2] {
        let a = self.stress_tensor[0];
        let b = self.stress_tensor[1]; // σ₂₁ (= σ₁₂ for symmetric tensors)
        let d = self.stress_tensor[3];
        let mean = (a + d) * 0.5;
        let half_diff = (a - d) * 0.5;
        let r = (half_diff * half_diff + b * b).sqrt();
        [mean + r, mean - r]
    }

    /// Modified tension-field stress tensor: compressive principal stresses
    /// are zeroed out (wrinkle constraint removes compression).
    ///
    /// Returns the modified 2×2 tensor in the same column-major layout.
    pub fn modified_tension_field(&self) -> [f64; 4] {
        let [s1, s2] = self.compute_principal_stresses();
        let s1_mod = s1.max(0.0);
        let s2_mod = s2.max(0.0);

        // Reconstruct in principal frame and rotate back — for a diagonal case.
        // Full reconstruction via spectral decomposition of 2×2.
        let a = self.stress_tensor[0];
        let b = self.stress_tensor[1];
        let d = self.stress_tensor[3];

        let half_diff = (a - d) * 0.5;
        let r = (half_diff * half_diff + b * b).sqrt();

        if r < 1e-15 {
            // Already diagonal.
            return [s1_mod, 0.0, 0.0, s2_mod];
        }

        // Angle of principal frame
        let cos2 = half_diff / r;
        let _sin2 = b / r;
        let cos_t = ((1.0 + cos2) * 0.5).sqrt();
        let sin_t = if b >= 0.0 {
            ((1.0 - cos2) * 0.5).sqrt()
        } else {
            -((1.0 - cos2) * 0.5).sqrt()
        };

        // Rotation: σ_orig = R · diag(s1_mod, s2_mod) · Rᵀ
        let s11 = s1_mod * cos_t * cos_t + s2_mod * sin_t * sin_t;
        let s12 = (s1_mod - s2_mod) * cos_t * sin_t;
        let s22 = s1_mod * sin_t * sin_t + s2_mod * cos_t * cos_t;

        [s11, s12, s12, s22]
    }
}

// ---------------------------------------------------------------------------
// SheetDraping
// ---------------------------------------------------------------------------

/// Simple mass-spring sheet draping under gravity.
///
/// Particles are connected by springs to their grid neighbors; fixed points
/// are pinned and do not move.
#[derive(Clone, Debug)]
pub struct SheetDraping {
    /// Vertex positions (m).
    pub vertices: Vec<[f64; 3]>,
    /// Indices of fixed (pinned) vertices.
    pub fix_points: Vec<usize>,
    /// Gravity acceleration vector (m/s²), e.g. `[0, -9.81, 0]`.
    pub gravity: [f64; 3],
    /// Grid width (number of vertices along X).
    pub nx: usize,
    /// Grid height (number of vertices along Y).
    pub ny: usize,
    /// Spring stiffness (N/m).
    pub stiffness: f64,
}

impl SheetDraping {
    /// Construct a flat rectangular sheet of `nx × ny` vertices with spacing
    /// `dx` (m) along X and `dy` (m) along Y, centred at the origin in XZ.
    pub fn new_flat(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        gravity: [f64; 3],
        stiffness: f64,
    ) -> Self {
        let mut vertices = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                vertices.push([i as f64 * dx, 0.0, j as f64 * dy]);
            }
        }
        Self {
            vertices,
            fix_points: Vec::new(),
            gravity,
            nx,
            ny,
            stiffness,
        }
    }

    /// Pin vertex at grid position (ix, iy).
    pub fn pin(&mut self, ix: usize, iy: usize) {
        self.fix_points.push(iy * self.nx + ix);
    }

    /// Compute drape by iterating a simple relaxation (spring + gravity).
    ///
    /// Returns the final vertex positions after `iterations` steps.
    pub fn compute_drape(&self, iterations: usize) -> Vec<[f64; 3]> {
        let n = self.vertices.len();
        let mut pos = self.vertices.clone();
        // Use a smaller dt to keep the explicit integration stable.
        // Stability criterion: dt < sqrt(m/k).
        let mass = 0.01_f64; // kg per vertex
        let dt = (0.1 * (mass / self.stiffness.max(1e-30)).sqrt()).min(0.002);
        let damping = 0.95_f64;

        let mut vel = vec![[0.0_f64; 3]; n];

        // Build neighbor list (4-connected grid).
        let neighbors: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|idx| {
                let ix = idx % self.nx;
                let iy = idx / self.nx;
                let mut nbrs = Vec::new();
                let rest = {
                    // rest length = distance from vertices array
                    let p = self.vertices[idx];
                    move |other: usize| -> f64 {
                        let q = self.vertices[other];
                        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2))
                            .sqrt()
                    }
                };
                if ix > 0 {
                    let j = idx - 1;
                    nbrs.push((j, rest(j)));
                }
                if ix + 1 < self.nx {
                    let j = idx + 1;
                    nbrs.push((j, rest(j)));
                }
                if iy > 0 {
                    let j = idx - self.nx;
                    nbrs.push((j, rest(j)));
                }
                if iy + 1 < self.ny {
                    let j = idx + self.nx;
                    nbrs.push((j, rest(j)));
                }
                nbrs
            })
            .collect();

        let fixed: std::collections::HashSet<usize> = self.fix_points.iter().copied().collect();

        for _ in 0..iterations {
            // Compute forces.
            let mut force = vec![[0.0_f64; 3]; n];
            for i in 0..n {
                // Gravity.
                force[i][0] += mass * self.gravity[0];
                force[i][1] += mass * self.gravity[1];
                force[i][2] += mass * self.gravity[2];
                // Spring forces.
                for &(j, rest_len) in &neighbors[i] {
                    let dx = pos[j][0] - pos[i][0];
                    let dy_v = pos[j][1] - pos[i][1];
                    let dz = pos[j][2] - pos[i][2];
                    let dist = (dx * dx + dy_v * dy_v + dz * dz).sqrt().max(1e-15);
                    let f_mag = self.stiffness * (dist - rest_len);
                    force[i][0] += f_mag * dx / dist;
                    force[i][1] += f_mag * dy_v / dist;
                    force[i][2] += f_mag * dz / dist;
                }
            }
            // Integrate.
            for i in 0..n {
                if fixed.contains(&i) {
                    continue;
                }
                vel[i][0] = (vel[i][0] + dt * force[i][0] / mass) * damping;
                vel[i][1] = (vel[i][1] + dt * force[i][1] / mass) * damping;
                vel[i][2] = (vel[i][2] + dt * force[i][2] / mass) * damping;
                pos[i][0] += dt * vel[i][0];
                pos[i][1] += dt * vel[i][1];
                pos[i][2] += dt * vel[i][2];
            }
        }
        pos
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Critical buckling stress (Pa) for a thin plate under in-plane compression.
///
/// Classical formula:  σ_cr = k π² E / (12 (1−ν²)) * (t/L)²
/// where k ≈ 4 for simply-supported edges (4-edge pinned).
pub fn critical_buckling_stress(e: f64, nu: f64, t: f64, l: f64) -> f64 {
    let k = 4.0_f64; // buckling coefficient for SS edges
    k * PI * PI * e / (12.0 * (1.0 - nu * nu)) * (t / l).powi(2)
}

/// Wrinkling number (dimensionless) characterising the severity of wrinkling.
///
/// w = (σ_max − σ_min) / (σ_max + σ_min + ε)
/// Returns 0 when both stresses are equal (no shear / biaxial) and approaches
/// 1 when the compressive stress dominates.
pub fn wrinkling_number(sigma_min: f64, sigma_max: f64) -> f64 {
    let denom = (sigma_max + sigma_min).abs() + 1e-30;
    (sigma_max - sigma_min).abs() / denom
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── ThinShell bending stiffness ──────────────────────────────────────────

    #[test]
    fn test_bending_stiffness_formula() {
        // D = E t³ / (12 (1−ν²))
        let shell = ThinShell::new(0.001, 200e9, 0.3);
        let d = shell.bending_stiffness();
        let expected = 200e9 * 0.001_f64.powi(3) / (12.0 * (1.0 - 0.3_f64.powi(2)));
        assert!(
            (d - expected).abs() < EPS * expected.abs() + EPS,
            "D={d}, expected={expected}"
        );
    }

    #[test]
    fn test_bending_stiffness_scales_cubic_with_thickness() {
        let shell1 = ThinShell::new(1e-3, 70e9, 0.33);
        let shell2 = ThinShell::new(2e-3, 70e9, 0.33);
        let ratio = shell2.bending_stiffness() / shell1.bending_stiffness();
        assert!(
            (ratio - 8.0).abs() < 1e-6,
            "D ∝ t³: ratio should be 8, got {ratio}"
        );
    }

    #[test]
    fn test_membrane_stiffness_formula() {
        let shell = ThinShell::new(2e-3, 100e9, 0.25);
        let n = shell.membrane_stiffness();
        let expected = 100e9 * 2e-3;
        assert!((n - expected).abs() < EPS);
    }

    #[test]
    fn test_membrane_stiffness_scales_linear_with_thickness() {
        let s1 = ThinShell::new(1e-3, 70e9, 0.3);
        let s2 = ThinShell::new(3e-3, 70e9, 0.3);
        let ratio = s2.membrane_stiffness() / s1.membrane_stiffness();
        assert!((ratio - 3.0).abs() < EPS);
    }

    #[test]
    fn test_thin_shell_positive_stiffness() {
        let s = ThinShell::new(1e-3, 200e9, 0.3);
        assert!(s.bending_stiffness() > 0.0);
        assert!(s.membrane_stiffness() > 0.0);
    }

    // ── WrinklingModel state ─────────────────────────────────────────────────

    #[test]
    fn test_wrinkling_state_taut_biaxial_tension() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        // Both strains positive → both stresses positive → Taut.
        let model = WrinklingModel::new(shell, [1e-3, 0.5e-3]);
        assert_eq!(model.determine_state(), WrinklingState::Taut);
    }

    #[test]
    fn test_wrinkling_state_slack_biaxial_compression() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        // Both strains negative → both stresses negative → Slack.
        let model = WrinklingModel::new(shell, [-1e-3, -0.5e-3]);
        assert_eq!(model.determine_state(), WrinklingState::Slack);
    }

    #[test]
    fn test_wrinkling_state_wrinkled_mixed() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        // ε₁ large positive, ε₂ large negative → Wrinkled.
        let model = WrinklingModel::new(shell, [5e-3, -5e-3]);
        assert_eq!(model.determine_state(), WrinklingState::Wrinkled);
    }

    #[test]
    fn test_wrinkle_wavelength_positive() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        let model = WrinklingModel::new(shell, [0.0, 0.0]);
        let lam = model.wrinkle_wavelength();
        assert!(lam > 0.0, "Wavelength must be positive, got {lam}");
    }

    #[test]
    fn test_wrinkle_wavelength_scales_with_thickness() {
        // λ ∝ t (since D ∝ t³ and N ∝ t, so sqrt(D/N) ∝ t).
        let s1 = ThinShell::new(1e-3, 200e9, 0.3);
        let s2 = ThinShell::new(2e-3, 200e9, 0.3);
        let lam1 = WrinklingModel::new(s1, [0.0, 0.0]).wrinkle_wavelength();
        let lam2 = WrinklingModel::new(s2, [0.0, 0.0]).wrinkle_wavelength();
        let ratio = lam2 / lam1;
        assert!((ratio - 2.0).abs() < 1e-6, "λ ∝ t: ratio={ratio}");
    }

    #[test]
    fn test_wrinkle_amplitude_zero_tension() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        let model = WrinklingModel::new(shell, [0.0, 0.0]);
        assert!((model.wrinkle_amplitude(0.0) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_wrinkle_amplitude_positive_tension() {
        let shell = ThinShell::new(1e-3, 200e9, 0.3);
        let model = WrinklingModel::new(shell, [1e-3, -1e-3]);
        let a = model.wrinkle_amplitude(1e4);
        assert!(a >= 0.0, "Amplitude must be non-negative");
    }

    #[test]
    fn test_wrinkling_state_clone() {
        let s = WrinklingState::Wrinkled;
        assert_eq!(s.clone(), WrinklingState::Wrinkled);
    }

    // ── TensionField ─────────────────────────────────────────────────────────

    #[test]
    fn test_principal_stresses_isotropic() {
        // Isotropic: σ = p I → both principal stresses = p.
        let tf = TensionField::isotropic(1000.0);
        let [s1, s2] = tf.compute_principal_stresses();
        assert!((s1 - 1000.0).abs() < 1e-6, "s1={s1}");
        assert!((s2 - 1000.0).abs() < 1e-6, "s2={s2}");
    }

    #[test]
    fn test_principal_stresses_uniaxial() {
        // [[σ, 0], [0, 0]] → principals are σ and 0.
        let tf = TensionField::new([500.0, 0.0, 0.0, 0.0]);
        let [s1, s2] = tf.compute_principal_stresses();
        let smax = s1.max(s2);
        let smin = s1.min(s2);
        assert!((smax - 500.0).abs() < 1e-6, "s_max={smax}");
        assert!(smin.abs() < 1e-6, "s_min={smin}");
    }

    #[test]
    fn test_principal_stresses_shear() {
        // Pure shear [[0, τ], [τ, 0]] → principals ±τ.
        let tau = 200.0_f64;
        let tf = TensionField::new([0.0, tau, tau, 0.0]);
        let [s1, s2] = tf.compute_principal_stresses();
        assert!((s1.abs() - tau).abs() < 1e-6, "s1={s1}");
        assert!((s2.abs() - tau).abs() < 1e-6, "s2={s2}");
        // They should be opposite sign.
        assert!((s1 + s2).abs() < 1e-6, "s1+s2={}", s1 + s2);
    }

    #[test]
    fn test_modified_tension_field_removes_compression() {
        // Tensor with one negative principal → modified should zero that component.
        let tf = TensionField::new([1000.0, 0.0, 0.0, -500.0]);
        let modified = tf.modified_tension_field();
        let tf_mod = TensionField::new(modified);
        let [s1, s2] = tf_mod.compute_principal_stresses();
        assert!(s1 >= -1e-6, "Modified s1 must be non-negative: {s1}");
        assert!(s2 >= -1e-6, "Modified s2 must be non-negative: {s2}");
    }

    #[test]
    fn test_modified_tension_field_isotropic_tension_unchanged() {
        // Fully tensile tensor should be unchanged.
        let tf = TensionField::isotropic(500.0);
        let orig = tf.compute_principal_stresses();
        let modified = tf.modified_tension_field();
        let tf_mod = TensionField::new(modified);
        let modp = tf_mod.compute_principal_stresses();
        assert!(
            (orig[0] - modp[0]).abs() < 1e-4,
            "Tensile stress must be unchanged"
        );
        assert!((orig[1] - modp[1]).abs() < 1e-4);
    }

    #[test]
    fn test_tension_field_new_roundtrip() {
        let t = [100.0, 50.0, 50.0, 200.0];
        let tf = TensionField::new(t);
        assert_eq!(tf.stress_tensor, t);
    }

    // ── SheetDraping ─────────────────────────────────────────────────────────

    #[test]
    fn test_sheet_draping_no_gravity_no_movement() {
        // Zero gravity → no movement from flat initial positions.
        let sheet = SheetDraping::new_flat(3, 3, 0.1, 0.1, [0.0, 0.0, 0.0], 100.0);
        let initial = sheet.vertices.clone();
        let result = sheet.compute_drape(50);
        for (i, (a, b)) in initial.iter().zip(result.iter()).enumerate() {
            let d2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
            assert!(d2.sqrt() < 1e-6, "Vertex {i} moved without gravity");
        }
    }

    #[test]
    fn test_sheet_draping_gravity_center_sags() {
        // 3x3 sheet pinned at corners, gravity downward → centre must sag.
        let mut sheet = SheetDraping::new_flat(3, 3, 0.1, 0.1, [0.0, -9.81, 0.0], 200.0);
        sheet.pin(0, 0);
        sheet.pin(2, 0);
        sheet.pin(0, 2);
        sheet.pin(2, 2);
        let result = sheet.compute_drape(200);
        let centre_y = result[4][1]; // centre vertex
        assert!(
            centre_y < -1e-4,
            "Centre should sag under gravity, y={centre_y}"
        );
    }

    #[test]
    fn test_sheet_draping_fixed_points_dont_move() {
        let mut sheet = SheetDraping::new_flat(3, 3, 0.1, 0.1, [0.0, -9.81, 0.0], 200.0);
        sheet.pin(0, 0);
        sheet.pin(2, 2);
        let initial = sheet.vertices.clone();
        let result = sheet.compute_drape(100);
        // Fixed vertices 0 and 8 must not move.
        for &idx in &sheet.fix_points {
            let d2 = (initial[idx][0] - result[idx][0]).powi(2)
                + (initial[idx][1] - result[idx][1]).powi(2)
                + (initial[idx][2] - result[idx][2]).powi(2);
            assert!(d2.sqrt() < 1e-10, "Fixed vertex {idx} should not move");
        }
    }

    #[test]
    fn test_sheet_draping_vertex_count_preserved() {
        let sheet = SheetDraping::new_flat(4, 5, 0.1, 0.1, [0.0, -9.81, 0.0], 100.0);
        let result = sheet.compute_drape(10);
        assert_eq!(result.len(), 20);
    }

    // ── Critical buckling stress ─────────────────────────────────────────────

    #[test]
    fn test_critical_buckling_stress_positive() {
        let sigma = critical_buckling_stress(200e9, 0.3, 1e-3, 0.1);
        assert!(sigma > 0.0, "Buckling stress must be positive, got {sigma}");
    }

    #[test]
    fn test_critical_buckling_stress_scales_inverse_square_length() {
        // σ_cr ∝ (t/L)² → doubling L quarters σ_cr.
        let s1 = critical_buckling_stress(200e9, 0.3, 1e-3, 0.1);
        let s2 = critical_buckling_stress(200e9, 0.3, 1e-3, 0.2);
        let ratio = s1 / s2;
        assert!((ratio - 4.0).abs() < 1e-6, "σ_cr ∝ (1/L)²: ratio={ratio}");
    }

    #[test]
    fn test_critical_buckling_stress_scales_cubic_thickness() {
        // σ_cr ∝ t² → doubling t quadruples σ_cr.
        let s1 = critical_buckling_stress(200e9, 0.3, 1e-3, 0.1);
        let s2 = critical_buckling_stress(200e9, 0.3, 2e-3, 0.1);
        let ratio = s2 / s1;
        assert!((ratio - 4.0).abs() < 1e-6, "σ_cr ∝ t²: ratio={ratio}");
    }

    // ── Wrinkling number ─────────────────────────────────────────────────────

    #[test]
    fn test_wrinkling_number_equal_stresses_zero() {
        // σ_max = σ_min → no wrinkling → number should be 0.
        let wn = wrinkling_number(100.0, 100.0);
        assert!(
            wn.abs() < EPS,
            "Equal stresses → wrinkling_number=0, got {wn}"
        );
    }

    #[test]
    fn test_wrinkling_number_positive_range() {
        let wn = wrinkling_number(-50.0, 200.0);
        assert!((0.0..=2.0).contains(&wn), "wrinkling_number={wn}");
    }

    #[test]
    fn test_wrinkling_number_large_difference() {
        // Large σ_max, small σ_min → wrinkling_number close to 1.
        let wn = wrinkling_number(0.0, 1e6);
        assert!(
            wn > 0.9,
            "Large shear → wrinkling_number should be ~1, got {wn}"
        );
    }
}
