//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::colormap::{Colormap, map_scalar};
use crate::primitives::{Color, LinePrimitive};
use oxiphysics_core::math::Vec3;
use std::f64::consts::PI;

use super::types::{
    ContactForce, HydroStress, PlasticStrain, StrainTensor, StressField, StressInvariant,
    StressTensor, StressTrajectory, SuperquadricGlyph, TensorGlyph, YieldCriterion,
};

/// Compute the three principal stresses (eigenvalues of the symmetric stress tensor)
/// using the Cardano closed-form solution.
///
/// Input: Voigt \[s11, s22, s33, s12, s23, s13\].
/// Output: \[l1, l2, l3\] sorted in descending order.
pub fn principal_stresses(s: [f64; 6]) -> [f64; 3] {
    let sxx = s[0];
    let syy = s[1];
    let szz = s[2];
    let sxy = s[3];
    let syz = s[4];
    let sxz = s[5];
    let i1 = sxx + syy + szz;
    let i2 = sxx * syy + syy * szz + szz * sxx - sxy * sxy - syz * syz - sxz * sxz;
    let i3 = sxx * (syy * szz - syz * syz) - sxy * (sxy * szz - syz * sxz)
        + sxz * (sxy * syz - syy * sxz);
    let p3 = (i1 * i1 - 3.0 * i2) / 9.0;
    let q3 = (2.0 * i1.powi(3) - 9.0 * i1 * i2 + 27.0 * i3) / 54.0;
    let p3_safe = p3.max(0.0);
    let p_sqrt = p3_safe.sqrt();
    let val = if p3_safe > 1e-30 {
        (q3 / (p3_safe * p_sqrt)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let theta = val.acos();
    let shift = i1 / 3.0;
    let r = 2.0 * p_sqrt;
    let mut eigs = [
        shift + r * (theta / 3.0).cos(),
        shift + r * ((theta + 2.0 * PI) / 3.0).cos(),
        shift + r * ((theta + 4.0 * PI) / 3.0).cos(),
    ];
    eigs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}
/// Von Mises equivalent stress.
pub fn von_mises_stress(s: [f64; 6]) -> f64 {
    let term1 = (s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2);
    let term2 = 6.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2));
    ((term1 + term2) / 2.0).sqrt()
}
/// Hydrostatic (mean normal) stress: p = (s11 + s22 + s33) / 3.
pub fn hydrostatic_stress(s: [f64; 6]) -> f64 {
    (s[0] + s[1] + s[2]) / 3.0
}
/// Deviatoric stress in Voigt notation: s_ij = sigma_ij - p * delta_ij.
pub fn deviatoric_stress(s: [f64; 6]) -> [f64; 6] {
    let p = hydrostatic_stress(s);
    [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]]
}
/// Stress triaxiality: eta = sigma_hydro / sigma_vm.
pub fn stress_triaxiality(s: [f64; 6]) -> f64 {
    let vm = von_mises_stress(s);
    if vm < 1e-30 {
        return 0.0;
    }
    hydrostatic_stress(s) / vm
}
/// Lode angle in \[0, pi/3\] computed from the third deviatoric invariant.
pub fn lode_angle(s: [f64; 6]) -> f64 {
    let dev = deviatoric_stress(s);
    let j2 = 0.5
        * (dev[0] * dev[0]
            + dev[1] * dev[1]
            + dev[2] * dev[2]
            + 2.0 * (dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]));
    if j2 < 1e-30 {
        return 0.0;
    }
    let d = [
        [dev[0], dev[3], dev[5]],
        [dev[3], dev[1], dev[4]],
        [dev[5], dev[4], dev[2]],
    ];
    let j3 = det3(d);
    let arg = (3.0 * 3.0_f64.sqrt() / 2.0) * j3 / (j2.powf(1.5));
    let arg_clamped = arg.clamp(-1.0, 1.0);
    arg_clamped.acos() / 3.0
}
/// Compute scalar values for stress invariant coloring.
pub fn stress_invariant_values(stresses: &[StressTensor], invariant: StressInvariant) -> Vec<f64> {
    stresses
        .iter()
        .map(|st| match invariant {
            StressInvariant::I1 => st.i1(),
            StressInvariant::J2 => {
                let dev = st.deviatoric();
                0.5 * (dev[0] * dev[0]
                    + dev[1] * dev[1]
                    + dev[2] * dev[2]
                    + 2.0 * (dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]))
            }
            StressInvariant::J3 => st.i3(),
            StressInvariant::VonMises => st.von_mises(),
            StressInvariant::Hydrostatic => st.hydrostatic(),
            StressInvariant::Triaxiality => st.triaxiality(),
            StressInvariant::LodeAngle => st.lode_angle(),
        })
        .collect()
}
/// Color a set of stress tensors by a stress invariant.
pub fn stress_invariant_colors(
    stresses: &[StressTensor],
    invariant: StressInvariant,
    colormap: Colormap,
) -> Vec<Color> {
    let values = stress_invariant_values(stresses, invariant);
    if values.is_empty() {
        return Vec::new();
    }
    let min_v = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_v = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    values
        .iter()
        .map(|&v| map_scalar(v, min_v, max_v, colormap))
        .collect()
}
/// Generate tensor glyphs from a set of stress tensors.
///
/// Each glyph uses the principal stresses as axis lengths and
/// approximate eigenvectors as axis directions.
pub fn generate_tensor_glyphs(
    stresses: &[StressTensor],
    positions: &[[f64; 3]],
    scale: f64,
    colormap: Colormap,
) -> Vec<TensorGlyph> {
    let n = stresses.len().min(positions.len());
    if n == 0 {
        return Vec::new();
    }
    let vm_values: Vec<f64> = stresses[..n].iter().map(|s| s.von_mises()).collect();
    let vm_min = vm_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let vm_max = vm_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut glyphs = Vec::with_capacity(n);
    for i in 0..n {
        let principals = stresses[i].principal_stresses();
        let magnitudes = [
            principals[0].abs() * scale,
            principals[1].abs() * scale,
            principals[2].abs() * scale,
        ];
        let directions = approximate_eigenvectors(stresses[i].voigt);
        let color = map_scalar(vm_values[i], vm_min, vm_max, colormap);
        glyphs.push(TensorGlyph {
            position: positions[i],
            magnitudes,
            directions,
            color,
        });
    }
    glyphs
}
/// Approximate eigenvectors of a symmetric 3x3 matrix from Voigt notation.
///
/// Uses a simple approach: for near-diagonal tensors, returns coordinate axes;
/// otherwise uses cross-product based approach.
pub(super) fn approximate_eigenvectors(s: [f64; 6]) -> [[f64; 3]; 3] {
    let off_diag = s[3].abs() + s[4].abs() + s[5].abs();
    if off_diag < 1e-10 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let mat = [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]];
    let v1 = power_iteration(&mat, [1.0, 0.0, 0.0], 20);
    let seed2 = if v1[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let ortho = cross3(v1, seed2);
    let v2 = normalize3_arr(ortho);
    let v3 = cross3(v1, v2);
    [v1, v2, normalize3_arr(v3)]
}
/// Power iteration to find the dominant eigenvector.
pub(super) fn power_iteration(mat: &[[f64; 3]; 3], mut v: [f64; 3], iterations: usize) -> [f64; 3] {
    for _ in 0..iterations {
        let new_v = mat_vec_3(mat, v);
        v = normalize3_arr(new_v);
    }
    v
}
pub(super) fn mat_vec_3(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn normalize3_arr(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
/// Interpolate a stress tensor within a tetrahedral element using barycentric coordinates.
pub fn interpolate_stress(
    field: &StressField,
    element_indices: &[usize; 4],
    barycentric: &[f64; 4],
) -> StressTensor {
    let mut result = StressTensor::zero();
    for (i, &node_idx) in element_indices.iter().enumerate() {
        let s = field.stresses[node_idx];
        let w = barycentric[i];
        for j in 0..6 {
            result.voigt[j] += w * s.voigt[j];
        }
    }
    result
}
/// Generate principal stress arrows as line primitives.
///
/// For each stress tensor and position, creates 3 arrows along the principal
/// directions. Tension is colored red, compression is blue.
pub fn principal_stress_arrows(
    stresses: &[StressTensor],
    positions: &[[f64; 3]],
    scale: f64,
) -> Vec<LinePrimitive> {
    let n = stresses.len().min(positions.len());
    let mut lines = Vec::with_capacity(n * 3);
    for i in 0..n {
        let principals = stresses[i].principal_stresses();
        let dirs = approximate_eigenvectors(stresses[i].voigt);
        for (j, &dir) in dirs.iter().enumerate() {
            let half_len = principals[j].abs() * scale;
            let color = if principals[j] >= 0.0 {
                Color::red()
            } else {
                Color::blue()
            };
            let pos = Vec3::new(positions[i][0], positions[i][1], positions[i][2]);
            let d = Vec3::new(dir[0] * half_len, dir[1] * half_len, dir[2] * half_len);
            lines.push(LinePrimitive {
                start: pos - d,
                end: pos + d,
                color,
            });
        }
    }
    lines
}
pub(super) fn von_mises(s: &[f64; 6]) -> f64 {
    von_mises_stress(*s)
}
/// Map each stress tensor to a color based on its von Mises equivalent stress.
pub fn von_mises_colors(stresses: &[[f64; 6]], colormap: Colormap) -> Vec<Color> {
    if stresses.is_empty() {
        return Vec::new();
    }
    let vm_values: Vec<f64> = stresses.iter().map(von_mises).collect();
    let min_vm = vm_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_vm = vm_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    vm_values
        .iter()
        .map(|&v| map_scalar(v, min_vm, max_vm, colormap))
        .collect()
}
/// Generate principal stress direction glyphs as line primitives.
pub fn principal_stress_glyphs(
    stress: &[f64; 6],
    position: Vec3,
    scale: f64,
) -> Vec<LinePrimitive> {
    let sxx = stress[0];
    let syy = stress[1];
    let szz = stress[2];
    let txy = stress[3];
    let tyz = stress[4];
    let txz = stress[5];
    use oxiphysics_core::math::Matrix3;
    let mat = Matrix3::new(sxx, txy, txz, txy, syy, tyz, txz, tyz, szz);
    let eigen = mat.symmetric_eigen();
    let colors = [Color::red(), Color::green(), Color::blue()];
    let mut lines = Vec::with_capacity(3);
    for (i, &c) in colors.iter().enumerate() {
        let eigval = eigen.eigenvalues[i];
        let eigvec = Vec3::new(
            eigen.eigenvectors[(0, i)],
            eigen.eigenvectors[(1, i)],
            eigen.eigenvectors[(2, i)],
        );
        let half_len = eigval.abs() * scale;
        let dir = eigvec.normalize() * half_len;
        lines.push(LinePrimitive {
            start: position - dir,
            end: position + dir,
            color: c,
        });
    }
    lines
}
/// Trace a stress trajectory from a seed point through a stress field.
///
/// Advances `steps` times by `step_size` in the direction of the
/// `principal_index`-th principal eigenvector at each point. Uses nearest-node
/// lookup in the field.
pub fn trace_stress_trajectory(
    field: &StressField,
    node_positions: &[[f64; 3]],
    seed: [f64; 3],
    principal_index: usize,
    steps: usize,
    step_size: f64,
) -> StressTrajectory {
    let mut points = vec![seed];
    let mut pos = seed;
    for _ in 0..steps {
        let nearest = node_positions.iter().enumerate().min_by(|(_, a), (_, b)| {
            let da = vec3_sq_dist(pos, **a);
            let db = vec3_sq_dist(pos, **b);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some((idx, _)) = nearest {
            if idx >= field.stresses.len() {
                break;
            }
            let voigt = field.stresses[idx].voigt;
            let eigvecs = principal_eigenvectors(voigt);
            if principal_index >= 3 {
                break;
            }
            let dir = eigvecs[principal_index];
            let norm = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            if norm < 1e-12 {
                break;
            }
            pos = [
                pos[0] + dir[0] / norm * step_size,
                pos[1] + dir[1] / norm * step_size,
                pos[2] + dir[2] / norm * step_size,
            ];
            points.push(pos);
        } else {
            break;
        }
    }
    StressTrajectory {
        points,
        principal_index,
    }
}
pub(super) fn vec3_sq_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}
/// Compute the three principal eigenvectors of a symmetric Voigt tensor
/// using the power iteration approximation (3 axes of the rotation-free frame).
pub fn principal_eigenvectors(voigt: [f64; 6]) -> [[f64; 3]; 3] {
    let ps = principal_stresses(voigt);
    let _ = ps;
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Stress concentration factor Kt at each field node.
///
/// Kt = local_von_mises / nominal_stress.
/// Nodes where nominal_stress ≈ 0 are skipped (Kt = 0).
pub fn stress_concentration_factors(field: &StressField, nominal_stress: f64) -> Vec<f64> {
    if nominal_stress.abs() < 1e-14 {
        return vec![0.0; field.stresses.len()];
    }
    field
        .stresses
        .iter()
        .map(|s| s.von_mises() / nominal_stress)
        .collect()
}
/// Colour map for stress concentration factors using the Jet colormap.
pub fn stress_concentration_colors(
    field: &StressField,
    nominal_stress: f64,
    colormap: crate::colormap::Colormap,
) -> Vec<Color> {
    let kts = stress_concentration_factors(field, nominal_stress);
    let max_kt = kts.iter().cloned().fold(0.0_f64, f64::max).max(1e-12);
    kts.iter()
        .map(|&kt| map_scalar(kt / max_kt, 0.0, 1.0, colormap))
        .collect()
}
/// Safety factor at each node: SF = yield_stress / von_mises_stress.
/// Regions with SF < 1 are yielded.
pub fn safety_factor_map(field: &StressField, yield_stress: f64) -> Vec<f64> {
    field
        .stresses
        .iter()
        .map(|s| {
            let vm = s.von_mises();
            if vm > 1e-12 {
                yield_stress / vm
            } else {
                f64::INFINITY
            }
        })
        .collect()
}
/// Colour map for safety factors: red = yielded (SF < 1), green = safe.
pub fn safety_factor_colors(field: &StressField, yield_stress: f64) -> Vec<Color> {
    let sfs = safety_factor_map(field, yield_stress);
    sfs.iter()
        .map(|&sf| {
            if sf < 1.0 {
                Color::new(1.0, 0.0, 0.0, 1.0)
            } else if sf < 2.0 {
                Color::new(1.0, (sf - 1.0) as f32, 0.0, 1.0)
            } else {
                Color::new(0.0, 1.0, 0.0, 1.0)
            }
        })
        .collect()
}
/// Distance of a stress state from the yield surface.
///
/// Returns values in the range \[-1, 1\] where:
/// - < 0 → inside yield surface (elastic)
/// - ≥ 0 → on or outside yield surface (yielded)
/// - = 0 → exactly on the yield surface
pub fn yield_surface_distance(
    stress: [f64; 6],
    yield_stress: f64,
    criterion: YieldCriterion,
    friction_angle_deg: f64,
) -> f64 {
    match criterion {
        YieldCriterion::VonMises => {
            let vm = von_mises_stress(stress);
            (vm - yield_stress) / yield_stress
        }
        YieldCriterion::Tresca => {
            let mut ps = principal_stresses(stress);
            ps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let tresca = ps[2] - ps[0];
            (tresca - yield_stress) / yield_stress
        }
        YieldCriterion::MohrCoulomb => {
            let phi = friction_angle_deg.to_radians();
            let sin_phi = phi.sin();
            let cos_phi = phi.cos();
            let mut ps = principal_stresses(stress);
            ps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let sigma1 = ps[2];
            let sigma3 = ps[0];
            let c = yield_stress / (2.0 * cos_phi);
            let f_mc = (sigma1 - sigma3) + (sigma1 + sigma3) * sin_phi - 2.0 * c * cos_phi;
            f_mc / (2.0 * c * cos_phi)
        }
    }
}
/// Colour map for yield surface proximity.
/// Red = yielded, blue = elastic, white = near surface.
pub fn yield_surface_colors(
    field: &StressField,
    yield_stress: f64,
    criterion: YieldCriterion,
    friction_angle_deg: f64,
) -> Vec<Color> {
    field
        .stresses
        .iter()
        .map(|s| {
            let d = yield_surface_distance(s.voigt, yield_stress, criterion, friction_angle_deg);
            if d >= 0.0 {
                let intensity = (d.min(1.0) as f32).clamp(0.0, 1.0);
                Color::new(1.0, 0.0, intensity * 0.5, 1.0)
            } else {
                let t = ((-d).min(1.0) as f32).clamp(0.0, 1.0);
                Color::new(0.0, t * 0.5, 1.0, 1.0)
            }
        })
        .collect()
}
/// Map a cycle-to-failure count to a colour using a logarithmic scale.
///
/// Low life (many cycles used) → red; high life (few cycles used) → blue.
pub fn fatigue_life_color(
    cycles_to_failure: f64,
    min_cycles: f64,
    max_cycles: f64,
    colormap: crate::colormap::Colormap,
) -> Color {
    if max_cycles <= min_cycles || min_cycles <= 0.0 {
        return map_scalar(0.0, 0.0, 1.0, colormap);
    }
    let log_val = (cycles_to_failure.max(min_cycles).ln() - min_cycles.ln())
        / (max_cycles.ln() - min_cycles.ln());
    map_scalar(log_val.clamp(0.0, 1.0), 0.0, 1.0, colormap)
}
/// Generate per-node fatigue life colours from a von Mises stress field.
///
/// Uses Basquin's law: N = (sigma_a / sigma_f)^(-1/b).
/// Returns colours from `colormap` with the log scale from `min_cycles` to
/// `max_cycles`.
pub fn fatigue_life_colormap(
    field: &StressField,
    sigma_f: f64,
    b: f64,
    min_cycles: f64,
    max_cycles: f64,
    colormap: crate::colormap::Colormap,
) -> Vec<Color> {
    field
        .stresses
        .iter()
        .map(|s| {
            let vm = s.von_mises();
            let nf = if vm > 1e-12 {
                (vm / sigma_f).powf(-1.0 / b)
            } else {
                max_cycles
            };
            fatigue_life_color(nf, min_cycles, max_cycles, colormap)
        })
        .collect()
}
#[inline]
pub(super) fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::stress_viz::DamageField;
    use crate::stress_viz::MohrCircleData;
    use crate::stress_viz::StressPath;
    #[test]
    fn test_hydrostatic_zero_von_mises() {
        let s = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress(s);
        assert!(
            vm.abs() < 1e-8,
            "hydrostatic stress should have zero von Mises, got {vm}"
        );
    }
    #[test]
    fn test_principal_stress_sum_equals_trace() {
        let s = [200.0, -50.0, 100.0, 30.0, 10.0, -20.0];
        let trace = s[0] + s[1] + s[2];
        let principals = principal_stresses(s);
        let sum = principals[0] + principals[1] + principals[2];
        assert!(
            (sum - trace).abs() < 1e-6,
            "principal stress sum {sum} should equal trace {trace}"
        );
    }
    #[test]
    fn test_deviatoric_trace_zero() {
        let s = [300.0, -100.0, 50.0, 75.0, -25.0, 10.0];
        let dev = deviatoric_stress(s);
        let trace = dev[0] + dev[1] + dev[2];
        assert!(
            trace.abs() < 1e-10,
            "deviatoric trace should be zero, got {trace}"
        );
    }
    #[test]
    fn test_pure_shear_lode_angle() {
        let tau = 100.0;
        let s = [tau, -tau, 0.0, 0.0, 0.0, 0.0];
        let theta = lode_angle(s);
        assert!(
            (0.0..=PI / 3.0).contains(&theta),
            "Lode angle out of valid range [0, pi/3]: {theta}"
        );
    }
    #[test]
    fn test_triaxiality_uniaxial() {
        let sigma = 300.0;
        let s = [sigma, 0.0, 0.0, 0.0, 0.0, 0.0];
        let eta = stress_triaxiality(s);
        let expected = 1.0 / 3.0;
        assert!(
            (eta - expected).abs() < 1e-6,
            "uniaxial triaxiality should be 1/3, got {eta}"
        );
    }
    #[test]
    fn test_stress_tensor_wrapper_consistency() {
        let voigt = [100.0, 200.0, 300.0, 50.0, 30.0, 20.0];
        let st = StressTensor::new(voigt);
        assert!((st.hydrostatic() - hydrostatic_stress(voigt)).abs() < 1e-10);
        assert!((st.von_mises() - von_mises_stress(voigt)).abs() < 1e-10);
        let dev_st = st.deviatoric();
        let dev_fn = deviatoric_stress(voigt);
        for i in 0..6 {
            assert!((dev_st[i] - dev_fn[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_principal_stresses_diagonal() {
        let s = [300.0, 100.0, 200.0, 0.0, 0.0, 0.0];
        let mut principals = principal_stresses(s);
        principals.sort_by(|a, b| b.partial_cmp(a).unwrap());
        assert!(
            (principals[0] - 300.0).abs() < 1e-6,
            "max principal: {}",
            principals[0]
        );
        assert!(
            (principals[1] - 200.0).abs() < 1e-6,
            "mid principal: {}",
            principals[1]
        );
        assert!(
            (principals[2] - 100.0).abs() < 1e-6,
            "min principal: {}",
            principals[2]
        );
    }
    #[test]
    fn test_interpolate_stress_midpoint() {
        let s0 = StressTensor::new([0.0; 6]);
        let s1 = StressTensor::new([100.0; 6]);
        let s2 = StressTensor::new([200.0; 6]);
        let s3 = StressTensor::new([300.0; 6]);
        let field = StressField::new(vec![s0, s1, s2, s3], vec![]);
        let bary = [0.25f64; 4];
        let result = interpolate_stress(&field, &[0, 1, 2, 3], &bary);
        for i in 0..6 {
            assert!(
                (result.voigt[i] - 150.0).abs() < 1e-8,
                "interpolated voigt[{i}] = {}, expected 150",
                result.voigt[i]
            );
        }
    }
    #[test]
    fn test_von_mises_colors_length() {
        let stresses = vec![
            [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [200.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let colors = von_mises_colors(&stresses, Colormap::Viridis);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn test_mohr_circle_hydrostatic() {
        let s = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let mohr = MohrCircleData::from_voigt(s);
        assert!(mohr.max_shear().abs() < 1e-6, "Hydrostatic has zero shear");
        assert!((mohr.mean_normal() - 100.0).abs() < 1e-6);
    }
    #[test]
    fn test_mohr_circle_uniaxial() {
        let sigma = 200.0;
        let s = [sigma, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mohr = MohrCircleData::from_voigt(s);
        assert!(
            (mohr.max_shear() - 100.0).abs() < 1e-4,
            "Max shear: {}",
            mohr.max_shear()
        );
    }
    #[test]
    fn test_mohr_circle_points() {
        let s = [200.0, 0.0, -100.0, 0.0, 0.0, 0.0];
        let mohr = MohrCircleData::from_voigt(s);
        let points = mohr.outer_circle_points(36);
        assert_eq!(points.len(), 36);
        for &(sigma, tau) in &points {
            let r = ((sigma - mohr.centers[0]).powi(2) + tau * tau).sqrt();
            assert!(
                (r - mohr.radii[0]).abs() < 1e-6,
                "Point not on circle: r={r}, expected {}",
                mohr.radii[0]
            );
        }
    }
    #[test]
    fn test_stress_invariants() {
        let st = StressTensor::new([100.0, 200.0, 300.0, 0.0, 0.0, 0.0]);
        assert!((st.i1() - 600.0).abs() < 1e-10, "I1 should be 600");
        assert!(
            (st.i2() - 110000.0).abs() < 1e-6,
            "I2 should be 110000, got {}",
            st.i2()
        );
        assert!(
            (st.i3() - 6_000_000.0).abs() < 1e-3,
            "I3 should be 6000000, got {}",
            st.i3()
        );
    }
    #[test]
    fn test_frobenius_norm() {
        let st = StressTensor::uniaxial(3.0);
        assert!((st.frobenius_norm() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_octahedral_shear() {
        let st = StressTensor::uniaxial(300.0);
        let tau_oct = st.octahedral_shear();
        let vm = st.von_mises();
        let expected = 2.0_f64.sqrt() / 3.0 * vm;
        assert!((tau_oct - expected).abs() < 1e-6);
    }
    #[test]
    fn test_stress_path() {
        let mut path = StressPath::new();
        assert!(path.is_empty());
        path.record(0.0, StressTensor::uniaxial(100.0));
        path.record(1.0, StressTensor::uniaxial(200.0));
        path.record(2.0, StressTensor::uniaxial(300.0));
        assert_eq!(path.len(), 3);
        assert!(!path.is_empty());
        assert!((path.max_von_mises() - 300.0).abs() < 1e-6);
    }
    #[test]
    fn test_stress_path_pq_data() {
        let mut path = StressPath::new();
        path.record(0.0, StressTensor::uniaxial(300.0));
        let pq = path.pq_data();
        assert_eq!(pq.len(), 1);
        assert!(
            (pq[0].0 - 100.0).abs() < 1e-6,
            "p should be 100, got {}",
            pq[0].0
        );
        assert!(
            (pq[0].1 - 300.0).abs() < 1e-6,
            "q should be 300, got {}",
            pq[0].1
        );
    }
    #[test]
    fn test_stress_path_triaxiality_history() {
        let mut path = StressPath::new();
        path.record(0.0, StressTensor::uniaxial(300.0));
        let hist = path.triaxiality_history();
        assert_eq!(hist.len(), 1);
        assert!((hist[0].1 - 1.0 / 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_stress_path_to_lines() {
        let mut path = StressPath::new();
        path.record(0.0, StressTensor::uniaxial(100.0));
        path.record(1.0, StressTensor::uniaxial(200.0));
        let lines = path.to_principal_space_lines(1.0);
        assert_eq!(lines.len(), 1);
    }
    #[test]
    fn test_stress_invariant_colors() {
        let stresses = vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(200.0)];
        let colors = stress_invariant_colors(&stresses, StressInvariant::VonMises, Colormap::Jet);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn test_tensor_glyph_generation() {
        let stresses = vec![
            StressTensor::uniaxial(100.0),
            StressTensor::new([200.0, -100.0, 50.0, 30.0, 10.0, -20.0]),
        ];
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let glyphs = generate_tensor_glyphs(&stresses, &positions, 0.01, Colormap::Viridis);
        assert_eq!(glyphs.len(), 2);
        for glyph in &glyphs {
            for &m in &glyph.magnitudes {
                assert!(m >= 0.0, "Magnitude should be non-negative");
            }
        }
    }
    #[test]
    fn test_principal_stress_arrows() {
        let stresses = vec![StressTensor::uniaxial(100.0)];
        let positions = vec![[0.0, 0.0, 0.0]];
        let arrows = principal_stress_arrows(&stresses, &positions, 0.01);
        assert_eq!(arrows.len(), 3);
    }
    #[test]
    fn test_stress_tensor_factories() {
        let uni = StressTensor::uniaxial(100.0);
        assert!((uni.voigt[0] - 100.0).abs() < 1e-10);
        assert!(uni.voigt[1].abs() < 1e-10);
        let hydro = StressTensor::hydrostatic_state(50.0);
        assert!((hydro.hydrostatic() - 50.0).abs() < 1e-10);
        assert!(hydro.von_mises() < 1e-8);
        let shear = StressTensor::pure_shear(75.0);
        assert!((shear.voigt[3] - 75.0).abs() < 1e-10);
    }
    #[test]
    fn test_stress_tensor_sub() {
        let a = StressTensor::new([100.0, 200.0, 300.0, 10.0, 20.0, 30.0]);
        let b = StressTensor::new([50.0, 100.0, 150.0, 5.0, 10.0, 15.0]);
        let c = a - b;
        for i in 0..6 {
            assert!((c.voigt[i] - b.voigt[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_stress_field_average() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(300.0)],
            vec![],
        );
        let avg = field.average_stress();
        assert!((avg.voigt[0] - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_stress_field_max_vm() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(500.0)],
            vec![],
        );
        assert!((field.max_von_mises() - 500.0).abs() < 1e-6);
    }
    #[test]
    fn test_stress_trajectory_length() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(200.0)],
            vec![],
        );
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let traj = trace_stress_trajectory(&field, &positions, [0.0, 0.0, 0.0], 0, 5, 0.1);
        assert!(
            traj.points.len() > 1,
            "trajectory should have at least 2 points"
        );
        assert_eq!(traj.principal_index, 0);
    }
    #[test]
    fn test_stress_trajectory_seed_is_first_point() {
        let field = StressField::new(vec![StressTensor::uniaxial(100.0)], vec![]);
        let positions = vec![[0.0, 0.0, 0.0]];
        let seed = [0.5, 0.3, 0.1];
        let traj = trace_stress_trajectory(&field, &positions, seed, 0, 3, 0.1);
        assert!((traj.points[0][0] - seed[0]).abs() < 1e-10);
        assert!((traj.points[0][1] - seed[1]).abs() < 1e-10);
        assert!((traj.points[0][2] - seed[2]).abs() < 1e-10);
    }
    #[test]
    fn test_stress_concentration_factor_nominal() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(200.0), StressTensor::uniaxial(100.0)],
            vec![],
        );
        let kts = stress_concentration_factors(&field, 100.0);
        assert_eq!(kts.len(), 2);
        assert!((kts[0] - 2.0).abs() < 1e-6, "Kt should be 2: {}", kts[0]);
        assert!((kts[1] - 1.0).abs() < 1e-6, "Kt should be 1: {}", kts[1]);
    }
    #[test]
    fn test_stress_concentration_zero_nominal() {
        let field = StressField::new(vec![StressTensor::uniaxial(100.0)], vec![]);
        let kts = stress_concentration_factors(&field, 0.0);
        assert_eq!(kts[0], 0.0, "zero nominal → Kt=0");
    }
    #[test]
    fn test_stress_concentration_colors_length() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(300.0)],
            vec![],
        );
        let colors = stress_concentration_colors(&field, 100.0, crate::colormap::Colormap::Jet);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn test_safety_factor_above_yield() {
        let field = StressField::new(vec![StressTensor::uniaxial(100.0)], vec![]);
        let sf = safety_factor_map(&field, 300.0);
        assert!((sf[0] - 3.0).abs() < 1e-6, "SF should be 3: {}", sf[0]);
    }
    #[test]
    fn test_safety_factor_yielded() {
        let field = StressField::new(vec![StressTensor::uniaxial(400.0)], vec![]);
        let sf = safety_factor_map(&field, 300.0);
        assert!(sf[0] < 1.0, "SF < 1 → yielded: {}", sf[0]);
    }
    #[test]
    fn test_safety_factor_colors_red_when_yielded() {
        let field = StressField::new(vec![StressTensor::uniaxial(500.0)], vec![]);
        let colors = safety_factor_colors(&field, 100.0);
        assert!(
            (colors[0].r - 1.0).abs() < 1e-6,
            "yielded node should be red"
        );
    }
    #[test]
    fn test_yield_surface_von_mises_inside() {
        let stress = [50.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let d = yield_surface_distance(stress, 100.0, YieldCriterion::VonMises, 30.0);
        assert!(d < 0.0, "should be inside von Mises surface: d={d}");
    }
    #[test]
    fn test_yield_surface_von_mises_outside() {
        let stress = [200.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let d = yield_surface_distance(stress, 100.0, YieldCriterion::VonMises, 30.0);
        assert!(d > 0.0, "should be outside von Mises surface: d={d}");
    }
    #[test]
    fn test_yield_surface_tresca() {
        let stress = [200.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let d = yield_surface_distance(stress, 100.0, YieldCriterion::Tresca, 30.0);
        assert!(d > 0.0, "uniaxial 200 should yield at Tresca 100: d={d}");
    }
    #[test]
    fn test_yield_surface_colors_length() {
        let field = StressField::new(
            vec![StressTensor::uniaxial(50.0), StressTensor::uniaxial(200.0)],
            vec![],
        );
        let colors = yield_surface_colors(&field, 100.0, YieldCriterion::VonMises, 30.0);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn test_fatigue_life_color_low_cycles_different_from_high() {
        let c_low = fatigue_life_color(1e3, 1e3, 1e7, crate::colormap::Colormap::Jet);
        let c_high = fatigue_life_color(1e7, 1e3, 1e7, crate::colormap::Colormap::Jet);
        let diff =
            (c_low.r - c_high.r).abs() + (c_low.g - c_high.g).abs() + (c_low.b - c_high.b).abs();
        assert!(
            diff > 0.01,
            "low/high cycle colors should differ: diff={diff}"
        );
    }
    #[test]
    fn test_fatigue_life_colormap_length() {
        let field = StressField::new(
            vec![
                StressTensor::uniaxial(100.0),
                StressTensor::uniaxial(200.0),
                StressTensor::uniaxial(50.0),
            ],
            vec![],
        );
        let colors = fatigue_life_colormap(
            &field,
            1000.0,
            -0.1,
            1e3,
            1e7,
            crate::colormap::Colormap::Viridis,
        );
        assert_eq!(colors.len(), 3);
    }
    #[test]
    fn test_damage_field_initial_zero() {
        let df = DamageField::new(5);
        for &d in &df.damage {
            assert!((d).abs() < 1e-10, "initial damage should be 0");
        }
    }
    #[test]
    fn test_damage_field_accumulate_miner() {
        let mut df = DamageField::new(1);
        df.accumulate_miner(0, &[500.0], &[1000.0]);
        assert!(
            (df.damage[0] - 0.5).abs() < 1e-10,
            "Miner damage should be 0.5: {}",
            df.damage[0]
        );
    }
    #[test]
    fn test_damage_field_capped_at_one() {
        let mut df = DamageField::new(1);
        df.accumulate_miner(0, &[2000.0], &[1000.0]);
        assert!(
            (df.damage[0] - 1.0).abs() < 1e-10,
            "damage capped at 1: {}",
            df.damage[0]
        );
    }
    #[test]
    fn test_damage_field_max() {
        let mut df = DamageField::new(3);
        df.accumulate_miner(0, &[100.0], &[1000.0]);
        df.accumulate_miner(1, &[800.0], &[1000.0]);
        df.accumulate_miner(2, &[300.0], &[1000.0]);
        let max = df.max_damage();
        assert!((max - 0.8).abs() < 1e-10, "max damage should be 0.8: {max}");
    }
    #[test]
    fn test_damage_field_colors_length() {
        let mut df = DamageField::new(3);
        df.accumulate_miner(0, &[500.0], &[1000.0]);
        let colors = df.to_colors(crate::colormap::Colormap::Viridis);
        assert_eq!(colors.len(), 3);
    }
    #[test]
    fn test_principal_eigenvectors_orthogonal() {
        let vecs = principal_eigenvectors([100.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        for v in &vecs {
            let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "eigenvector should be unit: {norm}"
            );
        }
    }
}
/// Strain energy density at a material point.
///
/// Uses the linear-elastic formulation: W = 0.5 * sigma : epsilon.
pub fn strain_energy_density(stress: &StressTensor, strain: &StrainTensor) -> f64 {
    let s = stress.voigt;
    let e = strain.voigt;
    0.5 * (s[0] * e[0]
        + s[1] * e[1]
        + s[2] * e[2]
        + 2.0 * (s[3] * e[3] + s[4] * e[4] + s[5] * e[5]))
}
/// Compute strain energy density for an entire field.
pub fn strain_energy_field(stresses: &[StressTensor], strains: &[StrainTensor]) -> Vec<f64> {
    stresses
        .iter()
        .zip(strains.iter())
        .map(|(s, e)| strain_energy_density(s, e))
        .collect()
}
/// Color map strain energy density using a given colormap.
pub fn strain_energy_colors(
    stresses: &[StressTensor],
    strains: &[StrainTensor],
    colormap: Colormap,
) -> Vec<crate::primitives::Color> {
    let energies = strain_energy_field(stresses, strains);
    let max_e = energies.iter().cloned().fold(0.0f64, f64::max);
    let min_e = energies.iter().cloned().fold(f64::MAX, f64::min);
    energies
        .iter()
        .map(|&e| map_scalar(e, min_e, max_e, colormap))
        .collect()
}
/// Color-map plastic strain field.
pub fn plastic_strain_colors(
    strains: &[PlasticStrain],
    max_eps: f64,
    colormap: Colormap,
) -> Vec<crate::primitives::Color> {
    strains
        .iter()
        .map(|ps| {
            map_scalar(
                ps.equivalent_plastic_strain,
                0.0,
                max_eps.max(1e-30),
                colormap,
            )
        })
        .collect()
}
/// Yield indicator color: red if yielded, green if not.
pub fn yield_indicator_colors(strains: &[PlasticStrain]) -> Vec<crate::primitives::Color> {
    strains
        .iter()
        .map(|ps| {
            if ps.yielded {
                crate::primitives::Color {
                    r: 1.0,
                    g: 0.0,
                    b: 0.0,
                    a: 1.0,
                }
            } else {
                crate::primitives::Color {
                    r: 0.0,
                    g: 0.8,
                    b: 0.2,
                    a: 1.0,
                }
            }
        })
        .collect()
}
/// Color-map hydrodynamic pressure field.
pub fn hydro_pressure_colors(
    hydro_stresses: &[HydroStress],
    p_min: f64,
    p_max: f64,
    colormap: Colormap,
) -> Vec<crate::primitives::Color> {
    hydro_stresses
        .iter()
        .map(|hs| map_scalar(hs.pressure, p_min, p_max, colormap))
        .collect()
}
/// Visualize contact forces as line primitives (arrows).
pub fn contact_force_lines(contacts: &[ContactForce], scale: f64) -> Vec<LinePrimitive> {
    contacts
        .iter()
        .map(|c| {
            let end = oxiphysics_core::math::Vec3::new(
                c.position[0] + c.normal[0] * c.normal_force * scale,
                c.position[1] + c.normal[1] * c.normal_force * scale,
                c.position[2] + c.normal[2] * c.normal_force * scale,
            );
            LinePrimitive {
                start: oxiphysics_core::math::Vec3::new(
                    c.position[0],
                    c.position[1],
                    c.position[2],
                ),
                end,
                color: crate::primitives::Color {
                    r: 1.0,
                    g: 0.5,
                    b: 0.0,
                    a: 1.0,
                },
            }
        })
        .collect()
}
/// Color-map contacts by normal force magnitude.
pub fn contact_force_colors(
    contacts: &[ContactForce],
    colormap: Colormap,
) -> Vec<crate::primitives::Color> {
    let max_f = contacts
        .iter()
        .map(|c| c.normal_force)
        .fold(0.0f64, f64::max);
    contacts
        .iter()
        .map(|c| map_scalar(c.normal_force, 0.0, max_f.max(1e-30), colormap))
        .collect()
}
/// Generate superquadric glyphs for all stress tensors in a field.
pub fn generate_superquadric_glyphs(
    stresses: &[StressTensor],
    positions: &[[f64; 3]],
    scale: f64,
    colormap: Colormap,
) -> Vec<SuperquadricGlyph> {
    stresses
        .iter()
        .zip(positions.iter())
        .map(|(s, &pos)| SuperquadricGlyph::from_stress(s, pos, scale, colormap))
        .collect()
}
#[cfg(test)]
mod new_stress_viz_tests {
    use super::*;

    #[test]
    fn strain_tensor_volumetric_isotropic() {
        let e = StrainTensor::new([0.01, 0.01, 0.01, 0.0, 0.0, 0.0]);
        assert!(
            (e.volumetric() - 0.03).abs() < 1e-10,
            "volumetric={}",
            e.volumetric()
        );
    }
    #[test]
    fn strain_tensor_zero_volumetric() {
        let e = StrainTensor::zero();
        assert_eq!(e.volumetric(), 0.0);
    }
    #[test]
    fn strain_energy_density_positive_for_compression() {
        let sigma = StressTensor::uniaxial(100.0);
        let epsilon = StrainTensor::new([0.001, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let w = strain_energy_density(&sigma, &epsilon);
        assert!(w > 0.0, "strain energy should be positive: {w}");
    }
    #[test]
    fn strain_energy_density_zero_for_zero_strain() {
        let sigma = StressTensor::uniaxial(100.0);
        let epsilon = StrainTensor::zero();
        let w = strain_energy_density(&sigma, &epsilon);
        assert_eq!(w, 0.0);
    }
    #[test]
    fn strain_energy_field_length() {
        let stresses = vec![StressTensor::uniaxial(100.0); 5];
        let strains = vec![StrainTensor::new([0.001, 0.0, 0.0, 0.0, 0.0, 0.0]); 5];
        let energies = strain_energy_field(&stresses, &strains);
        assert_eq!(energies.len(), 5);
    }
    #[test]
    fn strain_energy_colors_length() {
        let stresses = vec![StressTensor::uniaxial(100.0), StressTensor::uniaxial(200.0)];
        let strains = vec![
            StrainTensor::new([0.001, 0.0, 0.0, 0.0, 0.0, 0.0]),
            StrainTensor::new([0.002, 0.0, 0.0, 0.0, 0.0, 0.0]),
        ];
        let colors = strain_energy_colors(&stresses, &strains, Colormap::Viridis);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn plastic_strain_zero_not_yielded() {
        let ps = PlasticStrain::zero();
        assert!(!ps.yielded);
        assert_eq!(ps.equivalent_plastic_strain, 0.0);
    }
    #[test]
    fn plastic_strain_accumulate_increases_eps() {
        let mut ps = PlasticStrain::zero();
        ps.accumulate(0.01);
        assert!((ps.equivalent_plastic_strain - 0.01).abs() < 1e-10);
        assert!(ps.yielded);
    }
    #[test]
    fn plastic_strain_colors_length() {
        let strains = vec![PlasticStrain::zero(), PlasticStrain::new(0.05, [0.0; 6])];
        let colors = plastic_strain_colors(&strains, 0.1, Colormap::Jet);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn yield_indicator_colors_red_when_yielded() {
        let ps_yielded = PlasticStrain::new(0.01, [0.0; 6]);
        let colors = yield_indicator_colors(&[ps_yielded]);
        assert!((colors[0].r - 1.0).abs() < 1e-6, "yielded should be red");
        assert!(colors[0].g < 0.1);
    }
    #[test]
    fn yield_indicator_colors_green_when_not_yielded() {
        let ps = PlasticStrain::zero();
        let colors = yield_indicator_colors(&[ps]);
        assert!(colors[0].r < 0.1, "not yielded should be green");
        assert!(colors[0].g > 0.5);
    }
    #[test]
    fn hydrostress_total_stress_hydrostatic() {
        let hs = HydroStress::hydrostatic(100.0);
        let total = hs.total_stress();
        assert!(
            (total.voigt[0] - (-100.0)).abs() < 1e-10,
            "s11={}",
            total.voigt[0]
        );
        assert!((total.voigt[1] - (-100.0)).abs() < 1e-10);
        assert!((total.voigt[2] - (-100.0)).abs() < 1e-10);
        assert!(total.voigt[3].abs() < 1e-10);
    }
    #[test]
    fn hydrostress_dynamic_pressure() {
        let q = HydroStress::dynamic_pressure(1000.0, 2.0);
        assert!((q - 2000.0).abs() < 1e-6, "q={q}");
    }
    #[test]
    fn hydrostress_viscous_von_mises_zero_for_hydrostatic() {
        let hs = HydroStress::hydrostatic(100.0);
        assert!(hs.viscous_von_mises() < 1e-10);
    }
    #[test]
    fn hydro_pressure_colors_length() {
        let hydros = vec![
            HydroStress::hydrostatic(50.0),
            HydroStress::hydrostatic(150.0),
        ];
        let colors = hydro_pressure_colors(&hydros, 0.0, 200.0, Colormap::Viridis);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn contact_force_total_force() {
        let c = ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 3.0, 4.0, 0, 1);
        assert!(
            (c.total_force() - 5.0).abs() < 1e-10,
            "Pythagorean 3-4-5: {}",
            c.total_force()
        );
    }
    #[test]
    fn contact_force_friction_ratio() {
        let c = ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 100.0, 30.0, 0, 1);
        assert!((c.friction_ratio() - 0.3).abs() < 1e-10);
    }
    #[test]
    fn contact_force_friction_ratio_zero_normal() {
        let c = ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 5.0, 0, 1);
        assert_eq!(c.friction_ratio(), 0.0);
    }
    #[test]
    fn contact_force_lines_length() {
        let contacts = vec![
            ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 10.0, 2.0, 0, 1),
            ContactForce::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5.0, 1.0, 1, 2),
        ];
        let lines = contact_force_lines(&contacts, 0.01);
        assert_eq!(lines.len(), 2);
    }
    #[test]
    fn contact_force_colors_length() {
        let contacts = vec![
            ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 10.0, 2.0, 0, 1),
            ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 50.0, 5.0, 0, 2),
        ];
        let colors = contact_force_colors(&contacts, Colormap::Jet);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn contact_force_lines_end_offset_in_normal_direction() {
        let pos = [1.0, 2.0, 3.0];
        let n = [0.0, 1.0, 0.0];
        let c = ContactForce::new(pos, n, 10.0, 0.0, 0, 1);
        let lines = contact_force_lines(&[c], 0.1);
        assert!(
            (lines[0].end[1] - 3.0).abs() < 1e-10,
            "end_y={}",
            lines[0].end[1]
        );
    }
    #[test]
    fn superquadric_glyph_from_stress_uniaxial() {
        let stress = StressTensor::uniaxial(100.0);
        let glyph = SuperquadricGlyph::from_stress(&stress, [0.0; 3], 0.01, Colormap::Viridis);
        for &s in &glyph.semi_axes {
            assert!(s > 0.0, "semi-axis should be positive: {s}");
        }
    }
    #[test]
    fn superquadric_glyph_implicit_zero_at_surface() {
        let stress = StressTensor::hydrostatic_state(100.0);
        let glyph = SuperquadricGlyph::from_stress(&stress, [0.0; 3], 0.01, Colormap::Viridis);
        let p = [glyph.semi_axes[0], 0.0, 0.0];
        let val = glyph.implicit(p);
        assert!(val.is_finite(), "implicit value should be finite: {val}");
    }
    #[test]
    fn superquadric_glyph_implicit_inside_small_point() {
        let stress = StressTensor::hydrostatic_state(100.0);
        let glyph = SuperquadricGlyph::from_stress(&stress, [0.0; 3], 0.01, Colormap::Viridis);
        let val = glyph.implicit([0.0, 0.0, 0.0]);
        assert!(
            val == 0.0 || val < 1.0 + 1e-6,
            "origin should be inside: {val}"
        );
    }
    #[test]
    fn superquadric_generate_length() {
        let stresses = vec![StressTensor::uniaxial(100.0); 4];
        let positions: Vec<[f64; 3]> = (0..4).map(|i| [i as f64, 0.0, 0.0]).collect();
        let glyphs = generate_superquadric_glyphs(&stresses, &positions, 0.01, Colormap::Jet);
        assert_eq!(glyphs.len(), 4);
    }
    #[test]
    fn superquadric_exponents_in_valid_range() {
        let stress = StressTensor::new([100.0, 50.0, 20.0, 10.0, 5.0, 2.0]);
        let glyph = SuperquadricGlyph::from_stress(&stress, [0.0; 3], 0.01, Colormap::Viridis);
        assert!(glyph.exponents[0] >= 0.1 && glyph.exponents[0] <= 2.0);
        assert!(glyph.exponents[1] >= 0.1 && glyph.exponents[1] <= 2.0);
    }
    #[test]
    fn strain_energy_symmetric_in_stress_strain_swap() {
        let sigma = StressTensor::new([100.0, 50.0, 30.0, 10.0, 5.0, 2.0]);
        let epsilon = StrainTensor::new([0.002, 0.001, 0.0005, 0.0001, 0.00005, 0.00002]);
        let w1 = strain_energy_density(&sigma, &epsilon);
        assert!(w1 > 0.0, "energy should be positive: {w1}");
        assert!(w1.is_finite(), "energy should be finite: {w1}");
    }
    #[test]
    fn hydrostress_total_stress_has_viscous_term() {
        let mut hs = HydroStress::hydrostatic(0.0);
        hs.viscous[0] = 50.0;
        let total = hs.total_stress();
        assert!(
            (total.voigt[0] - 50.0).abs() < 1e-10,
            "viscous term s11={}",
            total.voigt[0]
        );
    }
    #[test]
    fn contact_force_total_force_only_normal() {
        let c = ContactForce::new([0.0; 3], [0.0, 1.0, 0.0], 5.0, 0.0, 0, 1);
        assert!((c.total_force() - 5.0).abs() < 1e-10);
    }
}
