//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// 3×3 matrix type (row-major).
pub type Mat3 = [[f64; 3]; 3];
/// Identity 3×3 matrix.
pub const IDENTITY3: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// Add two 3×3 matrices.
pub fn mat3_add(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}
/// Subtract two 3×3 matrices.
pub fn mat3_sub(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    c
}
/// Scale a 3×3 matrix by a scalar.
pub fn mat3_scale(a: &Mat3, s: f64) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}
/// Matrix product `a · b` of two 3×3 matrices.
pub fn mat3_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Matrix-vector product: a * v (3D).
pub fn mat3_vec(a: &Mat3, v: [f64; 3]) -> [f64; 3] {
    let mut r = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i] += a[i][j] * v[j];
        }
    }
    r
}
/// Trace of a 3×3 matrix.
pub fn mat3_trace(a: &Mat3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}
/// Dot product of two 3D vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean norm of a 3D vector.
#[inline]
pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3D vector (returns zero vector if near-zero).
#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Outer product n⊗n for a 3D unit vector.
pub fn outer3(a: [f64; 3], b: [f64; 3]) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i] * b[j];
        }
    }
    c
}
/// Compute eigenvalues of a symmetric 3×3 matrix (Jacobi iteration, 2D path for efficiency).
///
/// Returns sorted eigenvalues \[λ_min, λ_mid, λ_max\].
pub fn symmetric_eigenvalues_3x3(a: &Mat3) -> [f64; 3] {
    let p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    if p1 < 1e-30 {
        let mut eigs = [a[0][0], a[1][1], a[2][2]];
        eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        return eigs;
    }
    let q = mat3_trace(a) / 3.0;
    let b00 = a[0][0] - q;
    let b11 = a[1][1] - q;
    let b22 = a[2][2] - q;
    let p2 = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    if p < 1e-30 {
        return [q, q, q];
    }
    let b: Mat3 = [
        [(a[0][0] - q) / p, a[0][1] / p, a[0][2] / p],
        [a[1][0] / p, (a[1][1] - q) / p, a[1][2] / p],
        [a[2][0] / p, a[2][1] / p, (a[2][2] - q) / p],
    ];
    let det_b = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    let r = (det_b * 0.5).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
    let eig2 = 3.0 * q - eig1 - eig3;
    let mut eigs = [eig3, eig2, eig1];
    eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}
/// Cross product of two 3D vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Unit eigenvector of the symmetric matrix `a` belonging to its largest
/// eigenvalue (the direction of the maximum principal value).
///
/// With `λ_max` from [`symmetric_eigenvalues_3x3`], the matrix
/// `M = A − λ_max I` is negative semi-definite and its null space is exactly
/// the `λ_max` eigenspace. For a simple `λ_max` the null space is one
/// dimensional and orthogonal to every row of `M`, so the cross product of two
/// independent rows is the eigenvector; the cross product with the largest
/// norm is taken because it is the best conditioned. When `λ_max` is repeated
/// the rows of `M` are parallel, all cross products vanish, and any unit vector
/// orthogonal to a row lies in the eigenspace (one is returned). For `M ≈ 0`
/// (`A` is a multiple of the identity) every direction is an eigenvector and
/// `[1, 0, 0]` is returned.
///
/// The sign of an eigenvector is arbitrary; the result is oriented so that its
/// component of largest magnitude is positive, which makes it deterministic.
pub fn max_eigenvector_3x3(a: &Mat3) -> [f64; 3] {
    let lambda_max = symmetric_eigenvalues_3x3(a)[2];
    if !lambda_max.is_finite() {
        return [1.0, 0.0, 0.0];
    }
    let m = mat3_sub(a, &mat3_scale(&IDENTITY3, lambda_max));
    let rows = [m[0], m[1], m[2]];
    let row_norms = [norm3(rows[0]), norm3(rows[1]), norm3(rows[2])];
    let mut largest = 0;
    for k in 1..3 {
        if row_norms[k] > row_norms[largest] {
            largest = k;
        }
    }
    let scale = row_norms[largest];
    if scale.is_nan() || scale < 1e-30 {
        return [1.0, 0.0, 0.0];
    }
    let crosses = [
        cross3(rows[0], rows[1]),
        cross3(rows[0], rows[2]),
        cross3(rows[1], rows[2]),
    ];
    let mut best = crosses[0];
    let mut best_norm = norm3(best);
    for c in &crosses[1..] {
        let c_norm = norm3(*c);
        if c_norm > best_norm {
            best = *c;
            best_norm = c_norm;
        }
    }
    let direction = if best_norm > 1e-12 * scale * scale {
        normalize3(best)
    } else {
        // Repeated λ_max: pick a unit vector orthogonal to the dominant row by
        // crossing it with the coordinate axis it is least aligned with.
        let row = rows[largest];
        let mut axis = 0;
        for k in 1..3 {
            if row[k].abs() < row[axis].abs() {
                axis = k;
            }
        }
        let mut e = [0.0_f64; 3];
        e[axis] = 1.0;
        normalize3(cross3(row, e))
    };
    let mut pivot = 0;
    for k in 1..3 {
        if direction[k].abs() > direction[pivot].abs() {
            pivot = k;
        }
    }
    if direction[pivot] < 0.0 {
        [-direction[0], -direction[1], -direction[2]]
    } else {
        direction
    }
}
/// Wendland C2 kernel value in 3D.
///
/// `W(r, h) = (21/(16πh³)) * (1 - q)⁴ * (1 + 4q)` with `q = r/(2h)` for
/// `r < 2h`, and 0 beyond. The support radius is `2h` and the constant is the
/// canonical `21/(2πH³)` for support `H = 2h`, so `∫ W dV = 1`.
pub fn wendland_c2_3d(r: f64, h: f64) -> f64 {
    if h < 1e-30 {
        return 0.0;
    }
    let q = r / (2.0 * h);
    if q >= 1.0 {
        return 0.0;
    }
    let factor = 21.0 / (16.0 * PI * h * h * h);
    let t = 1.0 - q;
    factor * t * t * t * t * (1.0 + 4.0 * q)
}
/// Wendland C2 kernel gradient magnitude dW/dr (same normalisation as
/// [`wendland_c2_3d`]).
pub fn wendland_c2_3d_grad(r: f64, h: f64) -> f64 {
    if h < 1e-30 || r < 1e-30 {
        return 0.0;
    }
    let q = r / (2.0 * h);
    if q >= 1.0 {
        return 0.0;
    }
    let factor = 21.0 / (16.0 * PI * h * h * h);
    let t = 1.0 - q;
    let dw_dq = t * t * t * (-4.0 * (1.0 + 4.0 * q) + 4.0 * t);
    factor * dw_dq / (2.0 * h)
}
pub(super) trait CrackSurfaceHelper {
    fn smooth_h_or_default(&self, default: f64) -> f64;
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fracture_sph::types::*;
    fn make_particle(pos: [f64; 3], vel: [f64; 3]) -> SphFractureParticle {
        SphFractureParticle::new(pos, vel, 1e-3, 0.01, 2500.0)
    }
    fn make_damage_model() -> ContinuumDamageModel {
        ContinuumDamageModel::new_isotropic(1e-4, 50.0, 70e9, 0.3, 0.01)
    }
    fn make_fracture_energy() -> FractureEnergy {
        FractureEnergy::new(50.0, 70e9, 1e6, 1e-5, true)
    }
    fn make_simulation(n: usize) -> SphFractureSimulation {
        let particles: Vec<SphFractureParticle> = (0..n)
            .map(|i| make_particle([i as f64 * 0.01, 0.0, 0.0], [0.0; 3]))
            .collect();
        SphFractureSimulation::new(
            particles,
            make_damage_model(),
            make_fracture_energy(),
            1e-7,
            0.01,
        )
    }
    #[test]
    fn test_particle_new_intact() {
        let p = make_particle([0.0; 3], [0.0; 3]);
        assert!((p.damage - 0.0).abs() < 1e-14);
        assert!(p.is_intact(0.5));
    }
    #[test]
    fn test_particle_effective_stress_zero_damage() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress[0][0] = 1e6;
        p.stress[1][1] = 2e6;
        let eff = p.effective_stress();
        assert!((eff[0][0] - 1e6).abs() < 1.0);
        assert!((eff[1][1] - 2e6).abs() < 1.0);
    }
    #[test]
    fn test_particle_effective_stress_full_damage() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress[0][0] = 1e6;
        p.damage = 1.0;
        let eff = p.effective_stress();
        for row in &eff {
            for &val in row {
                assert!(val.abs() < 1e-10, "Fully damaged stress should be ~0");
            }
        }
    }
    #[test]
    fn test_particle_hydrostatic_pressure() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress[0][0] = 3e6;
        p.stress[1][1] = 3e6;
        p.stress[2][2] = 3e6;
        let ph = p.hydrostatic_pressure();
        assert!((ph + 3e6).abs() < 1.0, "hydrostatic pressure: {}", ph);
    }
    #[test]
    fn test_particle_von_mises_zero_for_isotropic() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress[0][0] = 1e6;
        p.stress[1][1] = 1e6;
        p.stress[2][2] = 1e6;
        p.update_deviatoric();
        let vm = p.von_mises_stress();
        assert!(vm < 1e-3, "VM stress should be ~0 for isotropic: {}", vm);
    }
    #[test]
    fn test_particle_distance_to() {
        let p1 = make_particle([0.0, 0.0, 0.0], [0.0; 3]);
        let p2 = make_particle([3.0, 4.0, 0.0], [0.0; 3]);
        assert!((p1.distance_to(&p2) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_particle_max_principal_stress_positive() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress[0][0] = 5e6;
        p.stress[1][1] = 2e6;
        p.stress[2][2] = 1e6;
        let sigma_max = p.max_principal_stress();
        assert!(
            (sigma_max - 5e6).abs() < 1e3,
            "max principal: {}",
            sigma_max
        );
    }
    #[test]
    fn test_particle_is_not_intact_when_damaged() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.damage = 0.8;
        assert!(!p.is_intact(0.5));
    }
    #[test]
    fn test_eigenvalues_diagonal() {
        let a: Mat3 = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        let eigs = symmetric_eigenvalues_3x3(&a);
        let expected = [1.0, 2.0, 3.0];
        for (e, &ex) in eigs.iter().zip(expected.iter()) {
            assert!((e - ex).abs() < 1e-6, "eigenvalue {} != {}", e, ex);
        }
    }
    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) - 14.0).abs() < 1e-14);
    }
    #[test]
    fn test_norm3() {
        assert!((norm3([3.0, 4.0, 0.0]) - 5.0).abs() < 1e-14);
    }
    #[test]
    fn test_normalize3_unit_vector() {
        let v = normalize3([3.0, 0.0, 4.0]);
        assert!((norm3(v) - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_normalize3_zero_vector() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert!((norm3(v)).abs() < 1e-14);
    }
    #[test]
    fn test_outer3_identity_property() {
        let e1 = [1.0, 0.0, 0.0];
        let o = outer3(e1, e1);
        assert!((o[0][0] - 1.0).abs() < 1e-14);
        assert!((o[0][1]).abs() < 1e-14);
    }
    #[test]
    fn test_mat3_trace() {
        let a: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        assert!((mat3_trace(&a) - 15.0).abs() < 1e-14);
    }
    /// `mat3_mul` is the row-by-column product: a worked example, the identity
    /// as neutral element, `(A B) v = A (B v)` and the non-commutativity that
    /// distinguishes `Ω σ − σ Ω` from zero.
    #[test]
    fn test_mat3_mul_is_the_row_by_column_product() {
        let a: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b: Mat3 = [[0.0, 1.0, 0.0], [-1.0, 0.0, 2.0], [3.0, 0.0, 1.0]];
        let ab = mat3_mul(&a, &b);
        assert_eq!(ab, [[7.0, 1.0, 7.0], [13.0, 4.0, 16.0], [19.0, 7.0, 25.0]]);
        assert_eq!(mat3_mul(&a, &IDENTITY3), a);
        assert_eq!(mat3_mul(&IDENTITY3, &a), a);
        let v = [1.0, -2.0, 0.5];
        let left = mat3_vec(&ab, v);
        let right = mat3_vec(&a, mat3_vec(&b, v));
        for (l, r) in left.iter().zip(right.iter()) {
            assert!((l - r).abs() < 1e-12, "(AB)v = {left:?}, A(Bv) = {right:?}");
        }
        assert_ne!(mat3_mul(&b, &a), ab);
    }
    #[test]
    fn test_damage_model_modulus() {
        let dm = make_damage_model();
        let e_eff = dm.effective_modulus(0.5);
        assert!((e_eff - 35e9).abs() < 1e3, "E_eff = {}", e_eff);
    }
    #[test]
    fn test_damage_model_zero_damage_full_modulus() {
        let dm = make_damage_model();
        assert!((dm.effective_modulus(0.0) - dm.youngs_modulus).abs() < 1.0);
    }
    #[test]
    fn test_damage_model_full_damage_zero_modulus() {
        let dm = make_damage_model();
        assert!(dm.effective_modulus(1.0).abs() < 1e-5);
    }
    #[test]
    fn test_damage_rate_zero_for_zero_strain_rate() {
        let dm = make_damage_model();
        let rate = dm.damage_rate_isotropic(0.1, 0.0, 1e-4);
        assert!((rate).abs() < 1e-14);
    }
    #[test]
    fn test_damage_rate_zero_for_full_damage() {
        let dm = make_damage_model();
        let rate = dm.damage_rate_isotropic(1.0, 1000.0, 1e-4);
        assert!((rate).abs() < 1e-14);
    }
    #[test]
    fn test_regularized_damage_single_neighbor() {
        let dm = make_damage_model();
        let d = dm.regularized_damage(&[0.5], &[0.0]);
        assert!((d - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_damage_model_bulk_modulus() {
        let dm = make_damage_model();
        let k = dm.bulk_modulus();
        let expected = 70e9 / (3.0 * (1.0 - 0.6));
        assert!(
            (k - expected).abs() < 1e3,
            "K = {}, expected {}",
            k,
            expected
        );
    }
    #[test]
    fn test_damage_model_shear_modulus() {
        let dm = make_damage_model();
        let g = dm.shear_modulus();
        let expected = 70e9 / (2.0 * 1.3);
        assert!((g - expected).abs() < 1e3);
    }
    #[test]
    fn test_damage_model_critical_energy_density() {
        let dm = make_damage_model();
        let wc = dm.critical_energy_density();
        assert!(wc > 0.0);
    }
    #[test]
    fn test_crack_new() {
        let c = CrackSurface::new([0.0; 3], 1e6);
        assert!((c.k_ic - 1e6).abs() < 1.0);
        assert!(!c.arrested);
        assert_eq!(c.crack_path.len(), 1);
    }
    #[test]
    fn test_crack_advance() {
        let mut c = CrackSurface::new([0.0, 0.0, 0.0], 1e6);
        c.propagation_dir = [1.0, 0.0, 0.0];
        c.advance(0.1);
        assert!((c.tip_position[0] - 0.1).abs() < 1e-12);
        assert_eq!(c.crack_path.len(), 2);
    }
    #[test]
    fn test_crack_arrest() {
        let mut c = CrackSurface::new([0.0; 3], 1e6);
        c.arrest();
        assert!(c.arrested);
        c.advance(1.0);
        assert!((c.crack_length).abs() < 1e-14);
    }
    #[test]
    fn test_crack_should_propagate_below_kic() {
        let mut c = CrackSurface::new([0.0; 3], 1e6);
        c.k_i = 0.5e6;
        c.k_ii = 0.0;
        assert!(!c.should_propagate());
    }
    #[test]
    fn test_crack_should_propagate_above_kic() {
        let mut c = CrackSurface::new([0.0; 3], 1e6);
        c.k_i = 2e6;
        c.k_ii = 0.0;
        assert!(c.should_propagate());
    }
    #[test]
    fn test_crack_mts_angle_zero_for_mode_i() {
        let mut c = CrackSurface::new([0.0; 3], 1e6);
        c.k_i = 1e6;
        c.k_ii = 0.0;
        let theta = c.mts_propagation_angle();
        assert!(theta.abs() < 1e-10, "Mode-I should give zero MTS angle");
    }
    #[test]
    fn test_crack_compute_k_factors() {
        let mut c = CrackSurface::new([0.0; 3], 1e6);
        c.compute_k_factors(1e6, 0.0, 0.01);
        assert!(c.k_i > 0.0);
        assert!((c.k_ii).abs() < 1e-10);
    }
    #[test]
    fn test_crack_nearest_particle() {
        let particles = vec![
            make_particle([1.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.1, 0.0, 0.0], [0.0; 3]),
            make_particle([5.0, 0.0, 0.0], [0.0; 3]),
        ];
        let c = CrackSurface::new([0.0, 0.0, 0.0], 1e6);
        let idx = c.nearest_particle_index(&particles);
        assert_eq!(idx, 1);
    }
    #[test]
    fn test_fracture_energy_kic_plane_stress() {
        let fe = make_fracture_energy();
        let kic = fe.k_ic(0.3);
        let expected = (50.0 * 70e9_f64).sqrt();
        assert!(
            (kic - expected).abs() < 1.0,
            "K_Ic = {}, expected {}",
            kic,
            expected
        );
    }
    #[test]
    fn test_fracture_energy_griffith_length() {
        let fe = make_fracture_energy();
        let a = fe.griffith_crack_length(1e6);
        assert!(a > 0.0);
        assert!(a.is_finite());
    }
    #[test]
    fn test_cohesive_traction_linear_at_zero() {
        let fe = make_fracture_energy();
        let t = fe.cohesive_traction_linear(0.0);
        assert!((t - fe.cohesive_traction).abs() < 1.0);
    }
    #[test]
    fn test_cohesive_traction_linear_at_delta_c() {
        let fe = make_fracture_energy();
        let t = fe.cohesive_traction_linear(fe.critical_separation);
        assert!(t.abs() < 1e-10);
    }
    #[test]
    fn test_cohesive_traction_exponential_positive() {
        let fe = make_fracture_energy();
        let t = fe.cohesive_traction_exponential(fe.critical_separation * 0.5);
        assert!(t > 0.0);
    }
    #[test]
    fn test_cohesive_energy_positive() {
        let fe = make_fracture_energy();
        assert!(fe.cohesive_energy() > 0.0);
    }
    #[test]
    fn test_griffith_criterion_met() {
        let fe = make_fracture_energy();
        assert!(!fe.griffith_criterion_met(fe.g_c * 0.5));
        assert!(fe.griffith_criterion_met(fe.g_c));
    }
    #[test]
    fn test_dugdale_process_zone_positive() {
        let fe = make_fracture_energy();
        let rp = fe.dugdale_process_zone(0.3);
        assert!(rp > 0.0);
    }
    #[test]
    fn test_fragment_tracking_single_cluster() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.005, 0.0, 0.0], [0.0; 3]),
            make_particle([0.01, 0.0, 0.0], [0.0; 3]),
        ];
        let mut ft = FragmentTracking::new(0.02, 0.99);
        let n = ft.label_fragments(&particles);
        assert_eq!(n, 1, "All intact particles → 1 fragment");
    }
    #[test]
    fn test_fragment_tracking_two_clusters() {
        let mut particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.005, 0.0, 0.0], [0.0; 3]),
            make_particle([1.0, 0.0, 0.0], [0.0; 3]),
            make_particle([1.005, 0.0, 0.0], [0.0; 3]),
        ];
        particles[2].damage = 0.0;
        let mut ft = FragmentTracking::new(0.02, 0.99);
        let n = ft.label_fragments(&particles);
        assert_eq!(n, 2, "Two separated clusters → 2 fragments");
    }
    #[test]
    fn test_fragment_tracking_isolated_damaged() {
        let mut particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.005, 0.0, 0.0], [0.0; 3]),
        ];
        particles[0].damage = 1.0;
        let mut ft = FragmentTracking::new(0.02, 0.99);
        let n = ft.label_fragments(&particles);
        assert_eq!(n, 1);
    }
    #[test]
    fn test_fragment_sizes_non_empty() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.005, 0.0, 0.0], [0.0; 3]),
        ];
        let mut ft = FragmentTracking::new(0.02, 0.99);
        ft.label_fragments(&particles);
        let sizes = ft.fragment_sizes();
        assert!(!sizes.is_empty());
    }
    #[test]
    fn test_fragment_masses_positive() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([0.005, 0.0, 0.0], [0.0; 3]),
        ];
        let mut ft = FragmentTracking::new(0.02, 0.99);
        ft.label_fragments(&particles);
        let masses = ft.fragment_masses(&particles);
        for &m in &masses {
            assert!(m > 0.0);
        }
    }
    #[test]
    fn test_fragment_centroids_computed() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], [0.0; 3]),
            make_particle([1.0, 0.0, 0.0], [0.0; 3]),
        ];
        let mut ft = FragmentTracking::new(0.5, 0.99);
        ft.label_fragments(&particles);
        let centroids = ft.fragment_centroids(&particles);
        assert!(!centroids.is_empty());
    }
    #[test]
    fn test_grady_kipp_mean_size_positive() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let s = gk.mean_fragment_size(1e4);
        assert!(s > 0.0 && s.is_finite());
    }
    #[test]
    fn test_grady_kipp_size_decreases_with_strain_rate() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let s1 = gk.mean_fragment_size(1e3);
        let s2 = gk.mean_fragment_size(1e6);
        assert!(s2 < s1, "Higher strain rate → smaller fragments");
    }
    #[test]
    fn test_grady_energy_size_positive() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let s = gk.grady_energy_fragment_size(1e4, 50.0);
        assert!(s > 0.0 && s.is_finite());
    }
    #[test]
    fn test_weibull_flaw_density_increases_with_strain() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let n1 = gk.weibull_flaw_density(1e-4);
        let n2 = gk.weibull_flaw_density(1e-3);
        assert!(n2 > n1);
    }
    #[test]
    fn test_weibull_flaw_density_zero_at_zero_strain() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        assert!((gk.weibull_flaw_density(0.0)).abs() < 1e-30);
    }
    #[test]
    fn test_fragment_cdf_zero_at_zero() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let p = gk.fragment_cdf(0.0, 1e4);
        assert!((p).abs() < 1e-14);
    }
    #[test]
    fn test_fragment_cdf_approaches_one() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let p = gk.fragment_cdf(1e3, 1e4);
        assert!(p > 0.99);
    }
    #[test]
    fn test_fragmentation_time_positive() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let t = gk.fragmentation_time(1e4);
        assert!(t > 0.0 && t.is_finite());
    }
    #[test]
    fn test_flaw_activation_strain_roundtrip() {
        let gk = ImpactFragmentation::new(1e-6, 3.0, 2500.0, 5000.0, 8.0, 1e15);
        let n = 1e12_f64;
        let eps = gk.flaw_activation_strain(n);
        let n2 = gk.weibull_flaw_density(eps);
        assert!((n2 - n).abs() / n < 1e-6, "n2 = {}, n = {}", n2, n);
    }
    #[test]
    fn test_simulation_new() {
        let sim = make_simulation(5);
        assert_eq!(sim.particles.len(), 5);
        assert!((sim.time).abs() < 1e-30);
        assert_eq!(sim.step_count, 0);
    }
    #[test]
    fn test_simulation_step_increments_time() {
        let mut sim = make_simulation(5);
        let dt = sim.dt;
        sim.step();
        assert!((sim.time - dt).abs() < 1e-20);
    }
    #[test]
    fn test_simulation_run_n_steps() {
        let mut sim = make_simulation(5);
        sim.run(10);
        assert_eq!(sim.step_count, 10);
    }
    #[test]
    fn test_simulation_damage_non_decreasing() {
        let mut sim = make_simulation(10);
        for p in sim.particles.iter_mut() {
            p.stress[0][0] = 1e9;
        }
        let d0: Vec<f64> = sim.particles.iter().map(|p| p.damage).collect();
        sim.step();
        for (i, p) in sim.particles.iter().enumerate() {
            assert!(
                p.damage >= d0[i] - 1e-12,
                "damage decreased at particle {}",
                i
            );
        }
    }
    #[test]
    fn test_simulation_max_damage_non_negative() {
        let sim = make_simulation(5);
        assert!(sim.max_damage() >= 0.0);
    }
    #[test]
    fn test_simulation_mean_damage_initially_zero() {
        let sim = make_simulation(5);
        assert!((sim.mean_damage()).abs() < 1e-14);
    }
    #[test]
    fn test_simulation_kinetic_energy_non_negative() {
        let mut sim = make_simulation(5);
        for p in sim.particles.iter_mut() {
            p.velocity = [1.0, 0.0, 0.0];
        }
        assert!(sim.kinetic_energy() > 0.0);
    }
    #[test]
    fn test_simulation_add_crack() {
        let mut sim = make_simulation(5);
        let crack = CrackSurface::new([0.0; 3], 1e6);
        sim.add_crack(crack);
        assert_eq!(sim.cracks.len(), 1);
    }
    #[test]
    fn test_simulation_identify_fragments_at_least_one() {
        let mut sim = make_simulation(5);
        let n = sim.identify_fragments();
        assert!(n > 0 || sim.particles.is_empty());
    }
    #[test]
    fn test_simulation_fully_damaged_count_zero_initially() {
        let sim = make_simulation(5);
        assert_eq!(sim.fully_damaged_count(), 0);
    }
    #[test]
    fn test_simulation_positions_change_with_nonzero_velocity() {
        let mut sim = make_simulation(3);
        for p in sim.particles.iter_mut() {
            p.velocity = [1.0, 0.0, 0.0];
        }
        let pos_before: Vec<[f64; 3]> = sim.particles.iter().map(|p| p.position).collect();
        sim.step();
        for (i, p) in sim.particles.iter().enumerate() {
            assert!((p.position[0] - pos_before[i][0]).abs() > 1e-20);
        }
    }
    #[test]
    fn test_wendland_kernel_zero_outside_support() {
        let w = wendland_c2_3d(3.0, 1.0);
        assert!((w).abs() < 1e-30);
    }
    #[test]
    fn test_wendland_kernel_positive_inside() {
        let w = wendland_c2_3d(0.5, 1.0);
        assert!(w > 0.0);
    }
    #[test]
    fn test_wendland_kernel_grad_zero_outside() {
        let g = wendland_c2_3d_grad(3.0, 1.0);
        assert!((g).abs() < 1e-30);
    }
    #[test]
    fn test_eigenvalues_sum_equals_trace() {
        let a: Mat3 = [[4.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 2.0]];
        let eigs = symmetric_eigenvalues_3x3(&a);
        let sum: f64 = eigs.iter().sum();
        let tr = mat3_trace(&a);
        assert!(
            (sum - tr).abs() < 1e-6,
            "sum of eigenvalues {} != trace {}",
            sum,
            tr
        );
    }
    /// Midpoint rule for `∫_0^{2h} 4π r² W(r) dr`.
    fn kernel_volume_integral(h: f64) -> f64 {
        let steps = 20_000;
        let dr = 2.0 * h / steps as f64;
        (0..steps)
            .map(|k| {
                let r = (k as f64 + 0.5) * dr;
                4.0 * PI * r * r * wendland_c2_3d(r, h) * dr
            })
            .sum()
    }
    #[test]
    fn test_wendland_c2_3d_is_normalised() {
        for h in [0.01, 1.0] {
            let integral = kernel_volume_integral(h);
            assert!(
                (integral - 1.0).abs() < 1e-4,
                "∫W dV = {integral} for h = {h}, expected 1"
            );
        }
    }
    #[test]
    fn test_wendland_c2_3d_central_value() {
        let h = 0.5;
        let expected = 21.0 / (16.0 * PI * h * h * h);
        assert!((wendland_c2_3d(0.0, h) - expected).abs() < 1e-12 * expected);
    }
    #[test]
    fn test_wendland_c2_3d_grad_is_derivative_of_kernel() {
        for h in [0.01, 1.0] {
            for frac in [0.2, 0.7, 1.3, 1.9] {
                let r = frac * h;
                let d = 1e-6 * h;
                let numeric = (wendland_c2_3d(r + d, h) - wendland_c2_3d(r - d, h)) / (2.0 * d);
                let analytic = wendland_c2_3d_grad(r, h);
                assert!(
                    (numeric - analytic).abs() < 1e-6 * analytic.abs(),
                    "dW/dr at r = {frac} h, h = {h}: numeric {numeric}, analytic {analytic}"
                );
            }
        }
    }
    /// Stress rate of a lattice in uniaxial strain rate: with a correctly
    /// normalised kernel the SPH velocity gradient reproduces `ε̇` and the
    /// stress increment is `(λ + 2μ) ε̇ dt` (an 8x over-normalised kernel
    /// gives eight times that).
    #[test]
    fn test_update_stress_uniaxial_strain_rate_matches_hooke() {
        let dx = 0.01;
        let n_side = 7;
        let strain_rate = 10.0;
        let dt = 1e-6;
        let density = 2500.0;
        let mass = density * dx * dx * dx;
        let mut particles = Vec::new();
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    let pos = [ix as f64 * dx, iy as f64 * dx, iz as f64 * dx];
                    let vel = [strain_rate * pos[0], 0.0, 0.0];
                    particles.push(SphFractureParticle::new(pos, vel, mass, 0.015, density));
                }
            }
        }
        let mut sim = SphFractureSimulation::new(
            particles,
            make_damage_model(),
            make_fracture_energy(),
            dt,
            0.015,
        );
        sim.update_stress();
        let centre = (n_side * n_side * n_side - 1) / 2;
        let e = 70e9;
        let nu = 0.3;
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let sigma_xx = sim.particles[centre].stress[0][0];
        let sigma_yy = sim.particles[centre].stress[1][1];
        let expected_xx = (lambda + 2.0 * mu) * strain_rate * dt;
        let expected_yy = lambda * strain_rate * dt;
        assert!(
            (sigma_xx - expected_xx).abs() < 0.03 * expected_xx,
            "σ_xx = {sigma_xx}, expected {expected_xx}"
        );
        assert!(
            (sigma_yy - expected_yy).abs() < 0.03 * expected_yy,
            "σ_yy = {sigma_yy}, expected {expected_yy}"
        );
    }
    #[test]
    fn test_max_eigenvector_diagonal_picks_largest_eigenvalue() {
        let a: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let v = max_eigenvector_3x3(&a);
        assert!(v[0].abs() < 1e-12 && v[1].abs() < 1e-12, "v = {v:?}");
        assert!((v[2].abs() - 1.0).abs() < 1e-12, "v = {v:?}");
        let b: Mat3 = [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]];
        let w = max_eigenvector_3x3(&b);
        assert!(w[0].abs() < 1e-12 && w[2].abs() < 1e-12, "w = {w:?}");
        assert!((w[1].abs() - 1.0).abs() < 1e-12, "w = {w:?}");
    }
    #[test]
    fn test_max_eigenvector_rotated_matrix() {
        // R diag(1, 2, 5) Rᵀ with R a 30° rotation about x: the dominant
        // eigenvector is R e_z = (0, -sin 30°, cos 30°).
        let theta = 30.0_f64.to_radians();
        let (s, c) = theta.sin_cos();
        let r: Mat3 = [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]];
        let d: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]];
        let rt: Mat3 = [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]];
        let mul = |x: &Mat3, y: &Mat3| -> Mat3 {
            let mut out = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        out[i][j] += x[i][k] * y[k][j];
                    }
                }
            }
            out
        };
        let a = mul(&mul(&r, &d), &rt);
        let v = max_eigenvector_3x3(&a);
        let expected = [0.0, -s, c];
        let sign = dot3(v, expected).signum();
        for k in 0..3 {
            assert!(
                (sign * v[k] - expected[k]).abs() < 1e-10,
                "v = {v:?}, expected ±{expected:?}"
            );
        }
    }
    #[test]
    fn test_max_eigenvector_repeated_maximum_lies_in_eigenspace() {
        // λ_max = 3 is double (x-y plane); any unit vector there is valid.
        let a: Mat3 = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]];
        let v = max_eigenvector_3x3(&a);
        assert!((norm3(v) - 1.0).abs() < 1e-12, "v = {v:?}");
        assert!(v[2].abs() < 1e-12, "v must lie in the x-y plane: {v:?}");
        // Isotropic stress: every direction is an eigenvector.
        let iso: Mat3 = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let w = max_eigenvector_3x3(&iso);
        assert!((norm3(w) - 1.0).abs() < 1e-12, "w = {w:?}");
    }
    #[test]
    fn test_max_eigenvector_is_deterministic_in_sign() {
        let a: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let v = max_eigenvector_3x3(&a);
        assert!(v[2] > 0.0, "largest component must be positive: {v:?}");
        let b: Mat3 = [[1.0, 0.5, 0.0], [0.5, -4.0, 0.0], [0.0, 0.0, 0.5]];
        let w = max_eigenvector_3x3(&b);
        let pivot = (0..3).fold(0, |m, k| if w[k].abs() > w[m].abs() { k } else { m });
        assert!(w[pivot] > 0.0, "largest component must be positive: {w:?}");
    }
    #[test]
    fn test_max_eigenvector_satisfies_eigen_equation_for_random_matrices() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        for _ in 0..200 {
            let (a00, a11, a22) = (next() * 1e6, next() * 1e6, next() * 1e6);
            let (a01, a02, a12) = (next() * 1e6, next() * 1e6, next() * 1e6);
            let a: Mat3 = [[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]];
            let lambda_max = symmetric_eigenvalues_3x3(&a)[2];
            let v = max_eigenvector_3x3(&a);
            assert!((norm3(v) - 1.0).abs() < 1e-12, "not a unit vector: {v:?}");
            let av = mat3_vec(&a, v);
            let residual = norm3([
                av[0] - lambda_max * v[0],
                av[1] - lambda_max * v[1],
                av[2] - lambda_max * v[2],
            ]);
            assert!(
                residual < 1e-8 * lambda_max.abs().max(1.0),
                "A v != λ_max v for {a:?}: residual {residual}, λ_max {lambda_max}"
            );
        }
    }
    /// The damage update must orient the crack plane along the maximum
    /// principal stress, not the minimum one.
    #[test]
    fn test_update_damage_aligns_crack_normal_with_max_principal_stress() {
        let mut p = make_particle([0.0; 3], [0.0; 3]);
        p.stress = [[1.0e6, 0.0, 0.0], [0.0, 5.0e8, 0.0], [0.0, 0.0, 2.0e6]];
        p.damage = 0.2;
        let mut sim = SphFractureSimulation::new(
            vec![p],
            make_damage_model(),
            make_fracture_energy(),
            1e-7,
            0.01,
        );
        sim.update_damage();
        let p = &sim.particles[0];
        assert!(p.damage > 0.01, "damage {} must be active", p.damage);
        let n = p.crack_normal;
        assert!(
            n[0].abs() < 1e-12 && n[2].abs() < 1e-12 && (n[1].abs() - 1.0).abs() < 1e-12,
            "crack normal {n:?} must be ±e_y (σ_yy is the max principal stress)"
        );
    }
}
