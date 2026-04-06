//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Eigen3;

/// A symmetric 3×3 tensor stored in Voigt notation: `[T_xx, T_yy, T_zz, T_xy, T_yz, T_xz]`.
pub type Tensor3 = [f64; 6];
/// A 3-D vector.
pub type Vec3f = [f64; 3];
/// A 3×3 matrix stored in row-major order.
pub type Mat3 = [f64; 9];
/// Convert a Voigt tensor to full 3×3 row-major matrix.
pub fn voigt_to_mat3(t: &Tensor3) -> Mat3 {
    [t[0], t[3], t[5], t[3], t[1], t[4], t[5], t[4], t[2]]
}
/// Convert a 3×3 full matrix to Voigt form (assuming symmetry).
pub fn mat3_to_voigt(m: &Mat3) -> Tensor3 {
    [m[0], m[4], m[8], m[1], m[5], m[2]]
}
/// Matrix–vector product for a 3×3 row-major matrix.
pub fn mat3_vec(m: &Mat3, v: &Vec3f) -> Vec3f {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}
/// Matrix–matrix product for two 3×3 row-major matrices.
pub fn mat3_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}
/// Transpose of a 3×3 row-major matrix.
pub fn mat3_transpose(m: &Mat3) -> Mat3 {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}
/// Identity 3×3 matrix.
pub fn mat3_identity() -> Mat3 {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}
/// Dot product of two 3-D vectors.
pub fn dot3(a: &Vec3f, b: &Vec3f) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// L2 norm of a 3-D vector.
pub fn norm3(v: &Vec3f) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3-D vector; returns zero vector if norm < ε.
pub fn normalize3(v: &Vec3f) -> Vec3f {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Cross product of two 3-D vectors.
pub fn cross3(a: &Vec3f, b: &Vec3f) -> Vec3f {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Scale a 3-D vector.
pub fn scale3(v: &Vec3f, s: f64) -> Vec3f {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Add two 3-D vectors.
pub fn add3(a: &Vec3f, b: &Vec3f) -> Vec3f {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-D vectors.
pub fn sub3(a: &Vec3f, b: &Vec3f) -> Vec3f {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Compute eigendecomposition of a symmetric 3×3 tensor using Jacobi iteration.
///
/// Returns eigenvalues in descending order with corresponding eigenvectors.
pub fn eig3(t: &Tensor3) -> Eigen3 {
    let mut a = voigt_to_mat3(t);
    let mut v = mat3_identity();
    pub(super) const MAX_ITER: usize = 64;
    pub(super) const EPS: f64 = 1e-12;
    for _ in 0..MAX_ITER {
        let mut max_val = 0.0_f64;
        let mut p = 0usize;
        let mut q = 1usize;
        for i in 0..3 {
            for j in (i + 1)..3 {
                let v_ij = a[i * 3 + j].abs();
                if v_ij > max_val {
                    max_val = v_ij;
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < EPS {
            break;
        }
        let a_pp = a[p * 3 + p];
        let a_qq = a[q * 3 + q];
        let a_pq = a[p * 3 + q];
        let theta = 0.5 * (a_qq - a_pp) / a_pq;
        let t_val = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            1.0 / (theta - (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t_val * t_val).sqrt();
        let s = t_val * c;
        let mut a_new = a;
        a_new[p * 3 + p] = a_pp - t_val * a_pq;
        a_new[q * 3 + q] = a_qq + t_val * a_pq;
        a_new[p * 3 + q] = 0.0;
        a_new[q * 3 + p] = 0.0;
        for r in 0..3 {
            if r != p && r != q {
                let a_rp = a[r * 3 + p];
                let a_rq = a[r * 3 + q];
                a_new[r * 3 + p] = c * a_rp - s * a_rq;
                a_new[p * 3 + r] = a_new[r * 3 + p];
                a_new[r * 3 + q] = s * a_rp + c * a_rq;
                a_new[q * 3 + r] = a_new[r * 3 + q];
            }
        }
        a = a_new;
        let mut v_new = v;
        for r in 0..3 {
            let v_rp = v[r * 3 + p];
            let v_rq = v[r * 3 + q];
            v_new[r * 3 + p] = c * v_rp - s * v_rq;
            v_new[r * 3 + q] = s * v_rp + c * v_rq;
        }
        v = v_new;
    }
    let mut vals = [a[0], a[4], a[8]];
    let mut vecs: [Vec3f; 3] = [[v[0], v[3], v[6]], [v[1], v[4], v[7]], [v[2], v[5], v[8]]];
    for i in 0..3 {
        for j in (i + 1)..3 {
            if vals[j] > vals[i] {
                vals.swap(i, j);
                vecs.swap(i, j);
            }
        }
    }
    let vecs_norm = [
        normalize3(&vecs[0]),
        normalize3(&vecs[1]),
        normalize3(&vecs[2]),
    ];
    Eigen3 {
        values: vals,
        vectors: vecs_norm,
    }
}
/// Linear interpolation of two Voigt tensors (not SPD-aware).
pub(super) fn voigt_lerp(a: &Tensor3, b: &Tensor3, t: f64) -> Tensor3 {
    std::array::from_fn(|i| a[i] * (1.0 - t) + b[i] * t)
}
/// Fractional anisotropy from three eigenvalues.
///
/// `FA = sqrt(3/2) * ||λ - λ̄I|| / ||λ||`  where `λ̄ = (λ₁+λ₂+λ₃)/3`.
pub fn fractional_anisotropy(eigenvalues: &[f64; 3]) -> f64 {
    let [l1, l2, l3] = *eigenvalues;
    let mean = (l1 + l2 + l3) / 3.0;
    let numerator = ((l1 - mean).powi(2) + (l2 - mean).powi(2) + (l3 - mean).powi(2)).sqrt();
    let denominator = (l1.powi(2) + l2.powi(2) + l3.powi(2)).sqrt();
    if denominator < 1e-15 {
        return 0.0;
    }
    (3.0_f64 / 2.0).sqrt() * numerator / denominator
}
/// Mean diffusivity (trace/3) from three eigenvalues.
pub fn mean_diffusivity(eigenvalues: &[f64; 3]) -> f64 {
    (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]) / 3.0
}
/// Axial diffusivity (largest eigenvalue).
pub fn axial_diffusivity(eigenvalues: &[f64; 3]) -> f64 {
    eigenvalues[0]
}
/// Radial diffusivity (average of two smaller eigenvalues).
pub fn radial_diffusivity(eigenvalues: &[f64; 3]) -> f64 {
    (eigenvalues[1] + eigenvalues[2]) / 2.0
}
/// Relative anisotropy.
pub fn relative_anisotropy(eigenvalues: &[f64; 3]) -> f64 {
    let md = mean_diffusivity(eigenvalues);
    if md.abs() < 1e-15 {
        return 0.0;
    }
    let [l1, l2, l3] = *eigenvalues;
    let num = ((l1 - md).powi(2) + (l2 - md).powi(2) + (l3 - md).powi(2)).sqrt();
    num / (3.0_f64.sqrt() * md)
}
/// Volume ratio (product of eigenvalues / cube of mean diffusivity).
pub fn volume_ratio(eigenvalues: &[f64; 3]) -> f64 {
    let md = mean_diffusivity(eigenvalues);
    if md.abs() < 1e-15 {
        return 0.0;
    }
    let prod = eigenvalues[0] * eigenvalues[1] * eigenvalues[2];
    prod / (md.powi(3))
}
pub(super) fn spectral_color(t: f64) -> (u8, u8, u8) {
    let r = (((t - 0.5) * 4.0).clamp(0.0, 1.0) * 255.0) as u8;
    let g = ((1.0 - (t * 2.0 - 1.0).abs().clamp(0.0, 1.0)) * 255.0) as u8;
    let b = ((1.0 - t * 2.0).clamp(0.0, 1.0) * 255.0) as u8;
    (r, g, b)
}
/// Von Mises equivalent stress from a Voigt stress tensor.
pub fn von_mises_stress(t: &Tensor3) -> f64 {
    let [sxx, syy, szz, sxy, syz, sxz] = *t;
    let term1 = (sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2);
    let term2 = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);
    (0.5 * (term1 + term2)).sqrt()
}
/// Compute the tensor invariants I₁, I₂, I₃ (principal invariants).
pub fn tensor_invariants(t: &Tensor3) -> [f64; 3] {
    let eig = eig3(t);
    let [l1, l2, l3] = eig.values;
    let i1 = l1 + l2 + l3;
    let i2 = l1 * l2 + l2 * l3 + l3 * l1;
    let i3 = l1 * l2 * l3;
    [i1, i2, i3]
}
/// Compute the Lode angle (normalized: 0 = triaxial tension, π/3 = triaxial compression).
pub fn lode_angle(t: &Tensor3) -> f64 {
    let eig = eig3(t);
    let [s1, s2, s3] = eig.values;
    let j2 = ((s1 - s2).powi(2) + (s2 - s3).powi(2) + (s3 - s1).powi(2)) / 6.0;
    let j3 =
        (s1 - (s1 + s2 + s3) / 3.0) * (s2 - (s1 + s2 + s3) / 3.0) * (s3 - (s1 + s2 + s3) / 3.0);
    if j2.abs() < 1e-15 {
        return 0.0;
    }
    let sin_3theta = -1.5 * (3.0_f64).sqrt() * j3 / j2.powf(1.5);
    sin_3theta.clamp(-1.0, 1.0).asin() / 3.0
}
/// Compute the Frobenius norm of a Voigt tensor.
pub fn frobenius_norm(t: &Tensor3) -> f64 {
    let m = voigt_to_mat3(t);
    m.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Check if a tensor is positive definite (all eigenvalues > 0).
pub fn is_positive_definite(t: &Tensor3) -> bool {
    let eig = eig3(t);
    eig.values[2] > 0.0
}
/// Check if a tensor is positive semi-definite (all eigenvalues ≥ 0).
pub fn is_positive_semidefinite(t: &Tensor3) -> bool {
    let eig = eig3(t);
    eig.values[2] >= -1e-12
}
/// Construct a random SPD tensor for testing.
#[allow(dead_code)]
pub(super) fn random_spd(seed: u64) -> Tensor3 {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    let l1 = next() * 2.0 + 0.5;
    let l2 = next() * 1.5 + 0.3;
    let l3 = next() * 1.0 + 0.1;
    [l1, l2, l3, 0.0, 0.0, 0.0]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor_visualization::DiffusionTensorImaging;
    use crate::tensor_visualization::DtiVoxel;
    use crate::tensor_visualization::HyperstreamlineTracer;
    use crate::tensor_visualization::StressTensorField;
    use crate::tensor_visualization::TensorGlyphBuilder;
    use crate::tensor_visualization::TensorInterpolation;
    use crate::tensor_visualization::TensorTopology;
    use crate::tensor_visualization::TopologyPointType;
    use std::f64::consts::PI;
    pub(super) const EPS: f64 = 1e-6;
    pub(super) const EPS_LOOSE: f64 = 1e-4;
    fn isotropic(lambda: f64) -> Tensor3 {
        [lambda, lambda, lambda, 0.0, 0.0, 0.0]
    }
    fn diag(l1: f64, l2: f64, l3: f64) -> Tensor3 {
        [l1, l2, l3, 0.0, 0.0, 0.0]
    }
    #[test]
    fn test_voigt_mat3_roundtrip() {
        let t = [1.0, 2.0, 3.0, 0.4, 0.5, 0.6_f64];
        let m = voigt_to_mat3(&t);
        let t2 = mat3_to_voigt(&m);
        for i in 0..6 {
            assert!((t[i] - t2[i]).abs() < EPS, "roundtrip failed at index {i}");
        }
    }
    #[test]
    fn test_mat3_identity_multiply() {
        let id = mat3_identity();
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0_f64];
        let result = mat3_mul(&id, &a);
        for i in 0..9 {
            assert!((result[i] - a[i]).abs() < EPS);
        }
    }
    #[test]
    fn test_normalize3_unit_vector() {
        let v = [3.0, 4.0, 0.0_f64];
        let n = normalize3(&v);
        let len = norm3(&n);
        assert!(
            (len - 1.0).abs() < EPS,
            "normalized vector should have unit length"
        );
    }
    #[test]
    fn test_normalize3_zero_vector() {
        let v = [0.0_f64; 3];
        let n = normalize3(&v);
        assert_eq!(n, [0.0; 3]);
    }
    #[test]
    fn test_dot3() {
        let a = [1.0, 0.0, 0.0_f64];
        let b = [0.0, 1.0, 0.0_f64];
        assert!(
            (dot3(&a, &b)).abs() < EPS,
            "orthogonal vectors should have zero dot product"
        );
        assert!((dot3(&a, &a) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_cross3_orthogonality() {
        let a = [1.0, 0.0, 0.0_f64];
        let b = [0.0, 1.0, 0.0_f64];
        let c = cross3(&a, &b);
        assert!((c[0]).abs() < EPS);
        assert!((c[1]).abs() < EPS);
        assert!((c[2] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_eig3_isotropic_eigenvalues() {
        let t = isotropic(3.0);
        let eig = eig3(&t);
        for val in eig.values {
            assert!(
                (val - 3.0).abs() < EPS_LOOSE,
                "isotropic tensor should have equal eigenvalues {val}"
            );
        }
    }
    #[test]
    fn test_eig3_diagonal_tensor_sorted() {
        let t = diag(5.0, 2.0, 8.0);
        let eig = eig3(&t);
        assert!(eig.values[0] >= eig.values[1]);
        assert!(eig.values[1] >= eig.values[2]);
        assert!((eig.values[0] - 8.0).abs() < EPS_LOOSE);
        assert!((eig.values[1] - 5.0).abs() < EPS_LOOSE);
        assert!((eig.values[2] - 2.0).abs() < EPS_LOOSE);
    }
    #[test]
    fn test_eig3_eigenvectors_orthonormal() {
        let t = [4.0, 3.0, 2.0, 1.0, 0.5, 0.25_f64];
        let eig = eig3(&t);
        for i in 0..3 {
            let n = norm3(&eig.vectors[i]);
            assert!((n - 1.0).abs() < EPS_LOOSE, "eigenvector {i} has norm {n}");
        }
        let d01 = dot3(&eig.vectors[0], &eig.vectors[1]).abs();
        let d02 = dot3(&eig.vectors[0], &eig.vectors[2]).abs();
        let d12 = dot3(&eig.vectors[1], &eig.vectors[2]).abs();
        assert!(d01 < EPS_LOOSE, "eigenvectors 0,1 not orthogonal: {d01}");
        assert!(d02 < EPS_LOOSE, "eigenvectors 0,2 not orthogonal: {d02}");
        assert!(d12 < EPS_LOOSE, "eigenvectors 1,2 not orthogonal: {d12}");
    }
    #[test]
    fn test_eig3_reconstruction() {
        let t = [3.0, 2.0, 1.0, 0.5, 0.3, 0.2_f64];
        let eig = eig3(&t);
        let m_orig = voigt_to_mat3(&t);
        let mut m_rec = [0.0_f64; 9];
        for k in 0..3 {
            let v = &eig.vectors[k];
            for i in 0..3 {
                for j in 0..3 {
                    m_rec[i * 3 + j] += eig.values[k] * v[i] * v[j];
                }
            }
        }
        for i in 0..9 {
            assert!(
                (m_orig[i] - m_rec[i]).abs() < 1e-3,
                "reconstruction failed at [{i}]: {} vs {}",
                m_orig[i],
                m_rec[i]
            );
        }
    }
    #[test]
    fn test_glyph_builder_semi_axes_positive() {
        let t = diag(4.0, 2.0, 1.0);
        let builder = TensorGlyphBuilder::new();
        let glyph = builder.build(&t, [0.0; 3]);
        for ax in glyph.semi_axes {
            assert!(ax > 0.0, "semi-axis must be positive");
        }
    }
    #[test]
    fn test_glyph_builder_scale() {
        let t = diag(1.0, 1.0, 1.0);
        let builder = TensorGlyphBuilder {
            scale: 2.0,
            ..TensorGlyphBuilder::new()
        };
        let glyph = builder.build(&t, [0.0; 3]);
        for ax in glyph.semi_axes {
            assert!(ax >= 1.9, "scale=2 should produce semi-axes ~2, got {ax}");
        }
    }
    #[test]
    fn test_glyph_builder_clamp() {
        let t = diag(100.0, 100.0, 100.0);
        let builder = TensorGlyphBuilder {
            max_axis: 5.0,
            ..TensorGlyphBuilder::new()
        };
        let glyph = builder.build(&t, [0.0; 3]);
        for ax in glyph.semi_axes {
            assert!(ax <= 5.0 + EPS, "semi-axis clamped to max_axis");
        }
    }
    #[test]
    fn test_glyph_builder_field() {
        let tensors: Vec<(Tensor3, Vec3f)> = (0..5)
            .map(|i| (diag(1.0 + i as f64, 1.0, 0.5), [i as f64, 0.0, 0.0]))
            .collect();
        let builder = TensorGlyphBuilder::new();
        let glyphs = builder.build_field(&tensors);
        assert_eq!(glyphs.len(), 5);
    }
    #[test]
    fn test_glyph_builder_grid() {
        let nx = 2;
        let ny = 2;
        let nz = 2;
        let field: Vec<Tensor3> = (0..nx * ny * nz).map(|_| diag(2.0, 1.0, 0.5)).collect();
        let builder = TensorGlyphBuilder::new();
        let glyphs = builder.build_grid(&field, nx, ny, nz, [0.0; 3], [1.0; 3]);
        assert_eq!(glyphs.len(), nx * ny * nz);
    }
    #[test]
    fn test_fa_isotropic_is_zero() {
        let fa = fractional_anisotropy(&[3.0, 3.0, 3.0]);
        assert!(fa.abs() < EPS, "isotropic FA should be 0, got {fa}");
    }
    #[test]
    fn test_fa_linear_is_near_one() {
        let fa = fractional_anisotropy(&[10.0, 0.1, 0.1]);
        assert!(fa > 0.9, "linear diffusion should have FA ~1, got {fa}");
    }
    #[test]
    fn test_fa_range() {
        let vals = [3.0, 1.5, 0.5_f64];
        let fa = fractional_anisotropy(&vals);
        assert!((0.0..=1.0).contains(&fa), "FA out of [0,1]: {fa}");
    }
    #[test]
    fn test_mean_diffusivity() {
        let md = mean_diffusivity(&[3.0, 2.0, 1.0]);
        assert!((md - 2.0).abs() < EPS);
    }
    #[test]
    fn test_axial_radial_diffusivity() {
        let ad = axial_diffusivity(&[5.0, 2.0, 1.0]);
        assert!((ad - 5.0).abs() < EPS);
        let rd = radial_diffusivity(&[5.0, 2.0, 1.0]);
        assert!((rd - 1.5).abs() < EPS);
    }
    #[test]
    fn test_relative_anisotropy() {
        let ra = relative_anisotropy(&[3.0, 3.0, 3.0]);
        assert!(ra.abs() < EPS, "isotropic RA should be 0");
    }
    #[test]
    fn test_log_exp_roundtrip() {
        let t = diag(4.0, 2.0, 1.0);
        let log_t = TensorInterpolation::log(&t);
        let t_rec = TensorInterpolation::exp(&log_t);
        for i in 0..6 {
            assert!(
                (t[i] - t_rec[i]).abs() < 1e-4,
                "log-exp roundtrip failed at [{i}]"
            );
        }
    }
    #[test]
    fn test_lerp_endpoints() {
        let t0 = diag(1.0, 1.0, 1.0);
        let t1 = diag(4.0, 4.0, 4.0);
        let t_start = TensorInterpolation::lerp(&t0, &t1, 0.0);
        let t_end = TensorInterpolation::lerp(&t0, &t1, 1.0);
        assert!(
            (t_start[0] - 1.0).abs() < 1e-3,
            "lerp at 0 should return t0"
        );
        assert!((t_end[0] - 4.0).abs() < 1e-3, "lerp at 1 should return t1");
    }
    #[test]
    fn test_lerp_midpoint_spd() {
        let t0 = diag(1.0, 1.0, 1.0);
        let t1 = diag(4.0, 4.0, 4.0);
        let tmid = TensorInterpolation::lerp(&t0, &t1, 0.5);
        assert!(
            is_positive_definite(&tmid),
            "interpolated tensor must be SPD"
        );
    }
    #[test]
    fn test_weighted_mean_equal_weights() {
        let t0 = diag(2.0, 2.0, 2.0);
        let t1 = diag(8.0, 8.0, 8.0);
        let mean = TensorInterpolation::weighted_mean(&[t0, t1], &[1.0, 1.0]);
        assert!(
            (mean[0] - 4.0).abs() < 1e-3,
            "equal-weight mean should be geometric mean ~4, got {}",
            mean[0]
        );
    }
    #[test]
    fn test_geodesic_distance_self_zero() {
        let t = diag(3.0, 2.0, 1.0);
        let d = TensorInterpolation::distance(&t, &t);
        assert!(d.abs() < EPS, "distance to self should be 0, got {d}");
    }
    #[test]
    fn test_geodesic_distance_positive() {
        let t0 = diag(1.0, 1.0, 1.0);
        let t1 = diag(4.0, 2.0, 1.0);
        let d = TensorInterpolation::distance(&t0, &t1);
        assert!(
            d > 0.0,
            "distance between different tensors should be positive"
        );
    }
    #[test]
    fn test_dti_process_length() {
        let tensors: Vec<Tensor3> = (0..10)
            .map(|i| diag(1.0 + i as f64 * 0.1, 1.0, 0.5))
            .collect();
        let dti = DiffusionTensorImaging::new();
        let voxels = dti.process(&tensors);
        assert_eq!(voxels.len(), 10);
    }
    #[test]
    fn test_dti_fa_color_map_length() {
        let tensors = vec![diag(5.0, 1.0, 0.5); 8];
        let dti = DiffusionTensorImaging::new();
        let voxels = dti.process(&tensors);
        let colors = dti.fa_color_map(&voxels);
        assert_eq!(colors.len(), 8);
    }
    #[test]
    fn test_dti_grayscale_isotropic() {
        let tensors = vec![diag(1.0, 1.0, 1.0)];
        let dti = DiffusionTensorImaging::new();
        let voxels = dti.process(&tensors);
        let gray = dti.fa_grayscale_map(&voxels);
        assert_eq!(gray[0], 0, "isotropic voxel should be black in FA map");
    }
    #[test]
    fn test_dti_tracking_mask() {
        let tensors = vec![diag(5.0, 0.5, 0.4), diag(1.0, 1.0, 1.0)];
        let dti = DiffusionTensorImaging::new();
        let voxels = dti.process(&tensors);
        let mask = dti.tracking_mask(&voxels);
        assert!(mask[0], "high-FA voxel should be tracked");
        assert!(!mask[1], "isotropic voxel should not be tracked");
    }
    #[test]
    fn test_dti_statistics() {
        let tensors: Vec<Tensor3> = (0..20)
            .map(|i| diag(1.0 + i as f64 * 0.5, 1.0, 0.5))
            .collect();
        let dti = DiffusionTensorImaging::new();
        let voxels = dti.process(&tensors);
        let stats = dti.statistics(&voxels);
        assert_eq!(stats.n_voxels, 20);
        assert!(stats.mean_fa >= 0.0 && stats.mean_fa <= 1.0);
        assert!(stats.max_fa >= stats.min_fa);
    }
    #[test]
    fn test_dti_voxel_from_tensor() {
        let t = diag(3.0, 1.5, 0.5);
        let v = DtiVoxel::from_tensor(t);
        assert!(v.fa >= 0.0 && v.fa <= 1.0);
        assert!(v.md > 0.0);
        assert!(v.ad >= v.rd);
    }
    #[test]
    fn test_von_mises_uniaxial() {
        let t = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let vm = von_mises_stress(&t);
        assert!(
            (vm - 100.0).abs() < EPS_LOOSE,
            "uniaxial von Mises should be 100, got {vm}"
        );
    }
    #[test]
    fn test_von_mises_hydrostatic_zero() {
        let t = isotropic(50.0);
        let vm = von_mises_stress(&t);
        assert!(
            vm.abs() < EPS_LOOSE,
            "hydrostatic von Mises should be 0, got {vm}"
        );
    }
    #[test]
    fn test_stress_decompose_hydrostatic() {
        let t = isotropic(10.0);
        let field = StressTensorField::new();
        let ps = field.decompose(&t);
        assert!((ps.hydrostatic - 10.0).abs() < EPS_LOOSE);
    }
    #[test]
    fn test_stress_decompose_deviatoric_traceless() {
        let t = [100.0, 50.0, 20.0, 10.0, 5.0, 3.0_f64];
        let field = StressTensorField::new();
        let ps = field.decompose(&t);
        let trace = ps.deviatoric[0] + ps.deviatoric[1] + ps.deviatoric[2];
        assert!(
            trace.abs() < EPS_LOOSE,
            "deviatoric tensor must be traceless, got {trace}"
        );
    }
    #[test]
    fn test_stress_trajectories_count() {
        let tensors = vec![diag(100.0, 0.0, -50.0); 5];
        let positions: Vec<Vec3f> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let field = StressTensorField::new();
        let lines = field.principal_trajectories(&tensors, &positions);
        assert_eq!(lines.len(), 15);
    }
    #[test]
    fn test_stress_von_mises_field() {
        let tensors = vec![diag(100.0, 0.0, 0.0); 4];
        let field = StressTensorField::new();
        let vm = field.von_mises_field(&tensors);
        assert_eq!(vm.len(), 4);
        for v in &vm {
            assert!((*v - 100.0).abs() < EPS_LOOSE);
        }
    }
    #[test]
    fn test_stress_hydrostatic_field() {
        let tensors = vec![isotropic(30.0); 3];
        let field = StressTensorField::new();
        let hyd = field.hydrostatic_field(&tensors);
        for h in &hyd {
            assert!((*h - 30.0).abs() < EPS, "hydrostatic should be 30, got {h}");
        }
    }
    #[test]
    fn test_topology_classify_regular() {
        let t = diag(3.0, 2.0, 1.0);
        let topo = TensorTopology::new();
        let pt = topo.classify_point(&t, [0.0; 3]);
        assert_eq!(pt.kind, TopologyPointType::Regular);
    }
    #[test]
    fn test_topology_classify_isotropic() {
        let t = isotropic(2.0);
        let topo = TensorTopology::new();
        let pt = topo.classify_point(&t, [0.0; 3]);
        assert_eq!(pt.kind, TopologyPointType::Isotropic);
    }
    #[test]
    fn test_topology_classify_planar_degenerate() {
        let t = diag(5.0, 1.0, 1.0);
        let topo = TensorTopology::new();
        let pt = topo.classify_point(&t, [0.0; 3]);
        assert_eq!(pt.kind, TopologyPointType::PlanarDegenerate);
    }
    #[test]
    fn test_topology_classify_linear_degenerate() {
        let t = diag(3.0, 3.0, 1.0);
        let topo = TensorTopology::new();
        let pt = topo.classify_point(&t, [0.0; 3]);
        assert_eq!(pt.kind, TopologyPointType::LinearDegenerate);
    }
    #[test]
    fn test_topology_find_degenerate_points() {
        let t_iso = isotropic(2.0);
        let t_reg = diag(3.0, 2.0, 1.0);
        let topo = TensorTopology::new();
        let pts = vec![
            topo.classify_point(&t_iso, [0.0; 3]),
            topo.classify_point(&t_reg, [1.0; 3]),
        ];
        let deg = topo.find_degenerate_points(&pts);
        assert_eq!(deg.len(), 1);
    }
    #[test]
    fn test_topology_degeneracy_measure_isotropic() {
        let t = isotropic(3.0);
        let topo = TensorTopology::new();
        let m = topo.degeneracy_measure(&t);
        assert!(
            m < 0.1,
            "isotropic tensor should have low degeneracy measure, got {m}"
        );
    }
    #[test]
    fn test_topology_classify_field() {
        let tensors = vec![diag(3.0, 2.0, 1.0), isotropic(2.0), diag(5.0, 1.0, 1.0)];
        let positions: Vec<Vec3f> = (0..3).map(|i| [i as f64, 0.0, 0.0]).collect();
        let topo = TensorTopology::new();
        let pts = topo.classify_field(&tensors, &positions);
        assert_eq!(pts.len(), 3);
        assert_eq!(pts[0].kind, TopologyPointType::Regular);
        assert_eq!(pts[1].kind, TopologyPointType::Isotropic);
        assert_eq!(pts[2].kind, TopologyPointType::PlanarDegenerate);
    }
    #[test]
    fn test_hyperstreamline_forward_moves() {
        let field_fn = |_p: Vec3f| -> Tensor3 { diag(5.0, 1.0, 0.5) };
        let tracer = HyperstreamlineTracer {
            max_steps: 10,
            bidirectional: false,
            min_fa: 0.05,
            ..HyperstreamlineTracer::new()
        };
        let sl = tracer.trace([0.0; 3], &field_fn);
        assert!(
            sl.points.len() > 1,
            "streamline should have more than 1 point"
        );
    }
    #[test]
    fn test_hyperstreamline_fa_above_threshold() {
        let field_fn = |_p: Vec3f| -> Tensor3 { diag(5.0, 1.0, 0.5) };
        let tracer = HyperstreamlineTracer {
            max_steps: 10,
            bidirectional: false,
            min_fa: 0.05,
            ..HyperstreamlineTracer::new()
        };
        let sl = tracer.trace([0.0; 3], &field_fn);
        for fa in &sl.fa {
            assert!(*fa >= 0.0 && *fa <= 1.0, "FA out of bounds: {fa}");
        }
    }
    #[test]
    fn test_hyperstreamline_stops_on_isotropic() {
        let high_fa_field = |_p: Vec3f| -> Tensor3 { diag(5.0, 1.0, 0.5) };
        let tracer = HyperstreamlineTracer {
            max_steps: 20,
            bidirectional: false,
            min_fa: 0.5,
            ..HyperstreamlineTracer::new()
        };
        let sl = tracer.trace([0.0; 3], &high_fa_field);
        assert!(!sl.points.is_empty());
    }
    #[test]
    fn test_hyperstreamline_bidirectional_longer() {
        let field_fn = |_p: Vec3f| -> Tensor3 { diag(5.0, 1.0, 0.5) };
        let tracer_bi = HyperstreamlineTracer {
            max_steps: 5,
            bidirectional: true,
            min_fa: 0.05,
            ..HyperstreamlineTracer::new()
        };
        let tracer_uni = HyperstreamlineTracer {
            max_steps: 5,
            bidirectional: false,
            min_fa: 0.05,
            ..HyperstreamlineTracer::new()
        };
        let sl_bi = tracer_bi.trace([0.0; 3], &field_fn);
        let sl_uni = tracer_uni.trace([0.0; 3], &field_fn);
        assert!(sl_bi.points.len() >= sl_uni.points.len());
    }
    #[test]
    fn test_hyperstreamline_many_seeds() {
        let field_fn = |_p: Vec3f| -> Tensor3 { diag(3.0, 1.0, 0.5) };
        let tracer = HyperstreamlineTracer {
            max_steps: 5,
            bidirectional: false,
            min_fa: 0.05,
            ..HyperstreamlineTracer::new()
        };
        let seeds: Vec<Vec3f> = (0..4).map(|i| [i as f64, 0.0, 0.0]).collect();
        let lines = tracer.trace_many(&seeds, &field_fn);
        assert_eq!(lines.len(), 4);
    }
    #[test]
    fn test_tensor_invariants_diagonal() {
        let t = diag(3.0, 2.0, 1.0);
        let [i1, i2, i3] = tensor_invariants(&t);
        assert!(
            (i1 - 6.0).abs() < EPS_LOOSE,
            "I1 should be trace=6, got {i1}"
        );
        assert!(
            (i2 - 11.0).abs() < EPS_LOOSE,
            "I2 should be 3*2+2*1+3*1=11, got {i2}"
        );
        assert!((i3 - 6.0).abs() < EPS_LOOSE, "I3 should be det=6, got {i3}");
    }
    #[test]
    fn test_frobenius_norm_identity() {
        let t = isotropic(1.0);
        let norm = frobenius_norm(&t);
        assert!(
            (norm - 3.0_f64.sqrt()).abs() < EPS_LOOSE,
            "Frobenius norm of I should be sqrt(3), got {norm}"
        );
    }
    #[test]
    fn test_is_positive_definite() {
        assert!(is_positive_definite(&diag(1.0, 1.0, 1.0)));
        assert!(!is_positive_definite(&diag(1.0, 1.0, -0.1)));
    }
    #[test]
    fn test_is_positive_semidefinite() {
        assert!(is_positive_semidefinite(&diag(1.0, 0.0, 1.0)));
        assert!(!is_positive_semidefinite(&diag(1.0, 1.0, -1.0)));
    }
    #[test]
    fn test_lode_angle_range() {
        let t = [100.0, 50.0, 20.0, 10.0, 5.0, 3.0_f64];
        let theta = lode_angle(&t);
        assert!(
            (-PI / 6.0 - 0.1..=PI / 6.0 + 0.1).contains(&theta),
            "Lode angle out of range: {theta}"
        );
    }
    #[test]
    fn test_volume_ratio_isotropic() {
        let vr = volume_ratio(&[2.0, 2.0, 2.0]);
        assert!((vr - 1.0).abs() < EPS, "isotropic VR should be 1, got {vr}");
    }
}
