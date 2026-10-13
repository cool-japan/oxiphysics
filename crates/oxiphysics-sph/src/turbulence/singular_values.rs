//! Singular values of 3×3 matrices and the Nicoud σ operator.
//!
//! Shared by [`crate::turbulence::SigmaModel`] and
//! [`crate::turbulence_sph::SigmaModel`] so that both LES closures use one
//! (tested) singular-value solver.
//!
//! Reference: Nicoud, Baya Toda, Cabrit, Bose & Lee (2011), "Using singular
//! values to build a subgrid-scale model for large eddy simulations",
//! Physics of Fluids 23, 085106.

use super::functions::Mat3;

/// Convergence threshold of the Jacobi sweeps: a column pair is considered
/// orthogonal when `|a_p · a_q| <= JACOBI_TOLERANCE |a_p| |a_q|`.
const JACOBI_TOLERANCE: f64 = 1.0e-15;

/// Upper bound on the number of Jacobi sweeps (3×3 matrices converge in about five).
const MAX_JACOBI_SWEEPS: usize = 30;

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Singular values `σ1 ≥ σ2 ≥ σ3 ≥ 0` of a real 3×3 matrix `g`.
///
/// One-sided Jacobi (Hestenes) method: plane rotations of the columns of `g`
/// until they are mutually orthogonal; the singular values are then the column
/// norms.  Unlike the eigenvalues of the Gram matrix `gᵀg` (obtained from the
/// characteristic cubic) this keeps full accuracy where the σ model needs it
/// most — at repeated or vanishing singular values, which lose half of the
/// digits (`σ3 ≈ 1e-8 σ1` instead of `0`) when `gᵀg` is formed explicitly.  A
/// rank-deficient `g` (two-dimensional flow, pure shear, solid rotation)
/// therefore yields `σ3` at rounding level (exactly `0` for axis-aligned
/// states), and repeated singular values agree to rounding level.
pub(crate) fn singular_values_3x3(g: &Mat3) -> [f64; 3] {
    let mut columns = [
        [g[0][0], g[1][0], g[2][0]],
        [g[0][1], g[1][1], g[2][1]],
        [g[0][2], g[1][2], g[2][2]],
    ];
    for _ in 0..MAX_JACOBI_SWEEPS {
        let mut rotated = false;
        for (p, q) in [(0_usize, 1_usize), (0, 2), (1, 2)] {
            let (col_p, col_q) = (columns[p], columns[q]);
            let alpha = dot(&col_p, &col_p);
            let beta = dot(&col_q, &col_q);
            let gamma = dot(&col_p, &col_q);
            if gamma.abs() <= JACOBI_TOLERANCE * (alpha * beta).sqrt() {
                continue;
            }
            rotated = true;
            // Rotation that annihilates the inner product of the two columns.
            let zeta = (beta - alpha) / (2.0 * gamma);
            let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = c * t;
            columns[p] = std::array::from_fn(|k| c * col_p[k] - s * col_q[k]);
            columns[q] = std::array::from_fn(|k| s * col_p[k] + c * col_q[k]);
        }
        if !rotated {
            break;
        }
    }
    let mut sigma = columns.map(|column| dot(&column, &column).sqrt());
    sigma.sort_by(|a, b| b.total_cmp(a));
    sigma
}

/// Nicoud's σ operator `D_σ = σ3 (σ1 − σ2)(σ2 − σ3) / σ1²` for singular
/// values sorted as `σ1 ≥ σ2 ≥ σ3 ≥ 0`.
///
/// It vanishes for a vanishing gradient, for two-dimensional or
/// two-component flow (`σ3 = 0`), for axisymmetric/isotropic expansion
/// (`σ1 = σ2` or `σ2 = σ3`), for pure shear and for solid-body rotation.
pub(crate) fn sigma_operator(sigma: [f64; 3]) -> f64 {
    let [s1, s2, s3] = sigma;
    if s1 <= 1.0e-30 || s1.is_nan() {
        return 0.0;
    }
    (s3 * (s1 - s2) * (s2 - s3) / (s1 * s1)).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::turbulence::SigmaModel;
    use crate::turbulence_sph::SigmaModel as SphSigmaModel;

    fn diag(a: f64, b: f64, c: f64) -> Mat3 {
        [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
    }

    fn det3(a: &Mat3) -> f64 {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    }

    fn matmul(a: &Mat3, b: &Mat3) -> Mat3 {
        let mut c = [[0.0; 3]; 3];
        for (c_row, a_row) in c.iter_mut().zip(a.iter()) {
            for (j, c_ij) in c_row.iter_mut().enumerate() {
                *c_ij = a_row
                    .iter()
                    .zip(b.iter())
                    .map(|(a_ik, b_row)| a_ik * b_row[j])
                    .sum();
            }
        }
        c
    }

    fn transpose(a: &Mat3) -> Mat3 {
        [
            [a[0][0], a[1][0], a[2][0]],
            [a[0][1], a[1][1], a[2][1]],
            [a[0][2], a[1][2], a[2][2]],
        ]
    }

    /// Proper rotation about the axis `(x, y, z)` (normalised internally) by `angle`.
    fn rotation(axis: [f64; 3], angle: f64) -> Mat3 {
        let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let (x, y, z) = (axis[0] / n, axis[1] / n, axis[2] / n);
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ]
    }

    /// `Q1 diag(a, b, c) Q2ᵀ` with two fixed, generic rotations.
    fn rotated_diagonal(a: f64, b: f64, c: f64) -> Mat3 {
        let q1 = rotation([1.0, 2.0, 3.0], 0.9);
        let q2 = rotation([-0.5, 0.1, 0.7], 2.4);
        matmul(&matmul(&q1, &diag(a, b, c)), &transpose(&q2))
    }

    /// Deterministic pseudo-random matrix with entries in `(-1, 1)`.
    fn pseudo_random_matrix(seed: u64) -> Mat3 {
        let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        };
        let mut m = [[0.0; 3]; 3];
        for row in m.iter_mut() {
            for v in row.iter_mut() {
                *v = next();
            }
        }
        m
    }

    fn frobenius_sq(a: &Mat3) -> f64 {
        a.iter().flatten().map(|v| v * v).sum()
    }

    fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol * (1.0 + expected.abs()),
            "{what}: got {actual}, expected {expected}"
        );
    }

    #[test]
    fn singular_values_of_rank_two_diagonal_gradient() {
        // G = gᵀg = diag(4, 1, 0): σ = (2, 1, 0).  The previous cubic had the wrong
        // sign of q and returned λ = (3.33, 2.33, -0.67) for this matrix.
        assert_eq!(singular_values_3x3(&diag(2.0, 1.0, 0.0)), [2.0, 1.0, 0.0]);
    }

    #[test]
    fn singular_values_of_full_rank_diagonal_gradient() {
        assert_eq!(singular_values_3x3(&diag(1.0, -2.0, 3.0)), [3.0, 2.0, 1.0]);
        // The order of the diagonal entries does not matter.
        assert_eq!(singular_values_3x3(&diag(-3.0, 1.0, 2.0)), [3.0, 2.0, 1.0]);
    }

    #[test]
    fn singular_values_satisfy_trace_and_determinant_identities() {
        for seed in 1..=60u64 {
            let g = pseudo_random_matrix(seed);
            let s = singular_values_3x3(&g);
            assert!(
                s[0] >= s[1] && s[1] >= s[2] && s[2] >= 0.0,
                "ordering {s:?}"
            );
            let sum_sq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
            assert_close(sum_sq, frobenius_sq(&g), 1e-13, "Σσ² = ‖g‖²_F");
            assert_close(s[0] * s[1] * s[2], det3(&g).abs(), 1e-13, "Πσ = |det g|");
            // The largest singular value bounds ‖g v‖ for every unit vector v.
            for axis in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
                let gv: f64 = g
                    .iter()
                    .map(|row| {
                        let component: f64 = row.iter().zip(axis.iter()).map(|(a, b)| a * b).sum();
                        component * component
                    })
                    .sum();
                assert!(gv.sqrt() <= s[0] * (1.0 + 1e-13), "seed {seed}");
            }
        }
    }

    #[test]
    fn singular_values_are_invariant_under_orthogonal_transforms() {
        let g = pseudo_random_matrix(77);
        let reference = singular_values_3x3(&g);
        let q1 = rotation([1.0, 2.0, 3.0], 0.9);
        let q2 = rotation([-0.5, 0.1, 0.7], 2.4);
        let rotated = matmul(&matmul(&q1, &g), &transpose(&q2));
        let s = singular_values_3x3(&rotated);
        for (value, expected) in s.iter().zip(reference.iter()) {
            assert_close(*value, *expected, 1e-13, "σ_k");
        }
    }

    #[test]
    fn repeated_and_vanishing_singular_values_keep_full_accuracy() {
        // Rotated matrices with a repeated singular value: no loss of digits
        // (the characteristic cubic of gᵀg only resolves them to about 1e-8).
        let s = singular_values_3x3(&rotated_diagonal(4.0, 1.0, 1.0));
        assert_close(s[0], 4.0, 1e-13, "σ1");
        assert_close(s[1], 1.0, 1e-13, "σ2");
        assert_close(s[2], 1.0, 1e-13, "σ3");
        let s = singular_values_3x3(&rotated_diagonal(3.0, 3.0, 0.5));
        assert_close(s[0], 3.0, 1e-13, "σ1");
        assert_close(s[1], 3.0, 1e-13, "σ2");
        assert_close(s[2], 0.5, 1e-13, "σ3");
        // Rank deficiency: the vanishing singular values stay at rounding level (not 1e-8).
        let s = singular_values_3x3(&rotated_diagonal(3.0, 2.0, 0.0));
        assert_close(s[0], 3.0, 1e-13, "σ1");
        assert_close(s[1], 2.0, 1e-13, "σ2");
        assert!(
            s[2] < 3e-14,
            "σ3 = {} must vanish for a rank-2 gradient",
            s[2]
        );
        let s = singular_values_3x3(&rotated_diagonal(5.0, 0.0, 0.0));
        assert_close(s[0], 5.0, 1e-13, "σ1");
        assert!(s[1] < 5e-14 && s[2] < 5e-14, "{s:?}");
    }

    #[test]
    fn singular_values_of_rotated_planar_flow_have_zero_third_value() {
        // A two-dimensional strain rotated into a generic orientation stays rank 2.
        let q = rotation([0.4, -0.7, 0.2], 1.3);
        let planar = [[1.5, 0.4, 0.0], [-0.2, -0.9, 0.0], [0.0, 0.0, 0.0]];
        let g = matmul(&matmul(&q, &planar), &transpose(&q));
        let s = singular_values_3x3(&g);
        assert!(
            s[2] < 1e-13 * s[0],
            "σ3 = {} must vanish for a planar flow",
            s[2]
        );
        assert!(s[1] > 0.1, "σ2 = {} should be O(1)", s[1]);
    }

    #[test]
    fn singular_values_of_zero_and_non_finite_input_terminate() {
        assert_eq!(singular_values_3x3(&[[0.0; 3]; 3]), [0.0, 0.0, 0.0]);
        let s = singular_values_3x3(&[[f64::NAN; 3]; 3]);
        assert!(s.iter().all(|v| v.is_nan()));
        let s = singular_values_3x3(&[[1.0, f64::INFINITY, 0.0], [0.0; 3], [0.0; 3]]);
        assert_eq!(sigma_operator(s), 0.0);
    }

    #[test]
    fn sigma_operator_known_value() {
        // σ = (3, 2, 1): D_σ = 1 · (3 − 2)(2 − 1) / 9 = 1/9.
        assert_close(sigma_operator([3.0, 2.0, 1.0]), 1.0 / 9.0, 1e-15, "D_σ");
    }

    #[test]
    fn sigma_operator_vanishes_in_degenerate_states() {
        assert_eq!(sigma_operator([0.0, 0.0, 0.0]), 0.0);
        assert_eq!(sigma_operator([2.0, 1.0, 0.0]), 0.0, "two-component flow");
        assert_eq!(sigma_operator([1.0, 1.0, 1.0]), 0.0, "isotropic expansion");
        assert_eq!(
            sigma_operator([2.0, 2.0, 1.0]),
            0.0,
            "axisymmetric (σ1 = σ2)"
        );
        assert_eq!(
            sigma_operator([2.0, 1.0, 1.0]),
            0.0,
            "axisymmetric (σ2 = σ3)"
        );
        assert_eq!(
            sigma_operator([f64::NAN, 1.0, 1.0]),
            0.0,
            "non-finite input"
        );
    }

    #[test]
    fn sigma_model_matches_closed_form_for_diagonal_gradient() {
        // g = diag(1, 2, 3), C_σ Δ = 1: ν_t = 1/9.  The sign error in the cubic gave 0.069.
        let model = SigmaModel::new(1.0, 1.0);
        assert_close(
            model.eddy_viscosity(diag(1.0, 2.0, 3.0)),
            1.0 / 9.0,
            1e-14,
            "ν_t",
        );
        // The filter scale enters squared.
        let scaled = SigmaModel::new(1.35, 0.1);
        let l = 1.35 * 0.1;
        assert_close(
            scaled.eddy_viscosity(diag(1.0, 2.0, 3.0)),
            l * l / 9.0,
            1e-14,
            "(C_σ Δ)² D_σ",
        );
    }

    #[test]
    fn sigma_model_vanishes_for_shear_rotation_and_planar_flows() {
        let model = SigmaModel::new(1.35, 0.1);
        let shear = [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let rotation_rate = [[0.0, -2.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let planar_strain = diag(3.0, -3.0, 0.0);
        let isotropic = diag(1.0, 1.0, 1.0);
        for (name, g) in [
            ("pure shear", shear),
            ("solid rotation", rotation_rate),
            ("planar strain", planar_strain),
            ("isotropic expansion", isotropic),
            ("zero gradient", [[0.0; 3]; 3]),
        ] {
            assert_eq!(model.eddy_viscosity(g), 0.0, "{name}");
        }
        // The same states in a rotated frame (g -> Q g Qᵀ) stay at rounding level.
        let q = rotation([0.3, 0.8, -0.5], 1.1);
        for (name, g) in [
            ("pure shear", shear),
            ("solid rotation", rotation_rate),
            ("planar strain", planar_strain),
            ("isotropic expansion", isotropic),
        ] {
            let rotated = matmul(&matmul(&q, &g), &transpose(&q));
            let nu_t = model.eddy_viscosity(rotated);
            assert!(nu_t.abs() < 1e-14, "{name} (rotated): ν_t = {nu_t}");
        }
    }

    #[test]
    fn sigma_model_is_frame_invariant_and_homogeneous() {
        let model = SigmaModel::new(1.35, 0.1);
        let g = [[0.9, 0.3, -0.2], [0.1, -0.4, 0.6], [0.5, -0.7, 0.2]];
        let reference = model.eddy_viscosity(g);
        assert!(
            reference > 0.0,
            "a generic 3-D gradient must produce ν_t > 0"
        );
        let q = rotation([0.6, 0.2, -0.9], 0.8);
        let rotated = matmul(&matmul(&q, &g), &transpose(&q));
        assert_close(model.eddy_viscosity(rotated), reference, 1e-12, "rotated g");
        let doubled = g.map(|row| row.map(|v| 2.0 * v));
        assert_close(
            model.eddy_viscosity(doubled),
            2.0 * reference,
            1e-12,
            "D_σ(2g) = 2 D_σ(g)",
        );
    }

    #[test]
    fn both_sigma_models_agree() {
        let delta = 0.05;
        let c_sigma = 1.35;
        let tensor_model = SigmaModel::new(c_sigma, delta);
        let sph_model = SphSigmaModel::new(c_sigma);
        for seed in 1..=25u64 {
            let g = pseudo_random_matrix(seed);
            let a = tensor_model.eddy_viscosity(g);
            let b = sph_model.eddy_viscosity(&g, delta);
            assert_close(
                a,
                b,
                1e-14,
                "turbulence::SigmaModel vs turbulence_sph::SigmaModel",
            );
        }
        let g = diag(1.0, 2.0, 3.0);
        assert_close(
            SphSigmaModel::new(1.0).eddy_viscosity(&g, 1.0),
            1.0 / 9.0,
            1e-14,
            "turbulence_sph::SigmaModel for diag(1,2,3)",
        );
        // Shear and rotation: zero for the root-exported model as well (it used to
        // return a Smagorinsky-like 0.1 |S| (C h)² there).
        let shear = [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let rotation_rate = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        assert_eq!(sph_model.eddy_viscosity(&shear, delta), 0.0);
        assert_eq!(sph_model.eddy_viscosity(&rotation_rate, delta), 0.0);
    }
}
