//! Germano–Lilly dynamic procedure for the Smagorinsky coefficient on particle data.
//!
//! The Smagorinsky closure `τ_ij − δ_ij τ_kk/3 = −2 C Δ² |S̄| S̄_ij` has the
//! squared constant `C = C_s²`.  The dynamic procedure (Germano et al., 1991;
//! Lilly, 1992) determines `C` from the resolved flow itself with a second,
//! coarser *test filter* (`^`, width `Δ̂ = α Δ`):
//!
//! ```text
//! L_ij = (ū_i ū_j)^ − û_i û_j                     (resolved, "Leonard" stress)
//! M_ij = 2 Δ² [ (|S̄| S̄_ij)^ − α² |Ŝ| Ŝ_ij ]        (Smagorinsky model difference)
//! L_ij − δ_ij L_kk/3 = C M_ij                      (Germano identity)
//! C    = ⟨L^d_ij M_ij⟩ / ⟨M_ij M_ij⟩               (Lilly's least squares)
//! ```
//!
//! On particles the test filter is the box filter over the particle itself and
//! its neighbour list (the same stencil as
//! [`box_filter_scalar`](super::functions::box_filter_scalar)), and `⟨·⟩` is the
//! average of the local contractions over that stencil, which keeps the
//! coefficient well conditioned.  The ratio `α = Δ̂/Δ` should match the extent
//! of the stencil: a top-hat filter whose per-axis second moment is `⟨δx²⟩` has
//! the width `Δ̂ = sqrt(12 ⟨δx²⟩)`.
//!
//! Unlike a strain-only formulation, `L_ij` depends on the *velocity* field
//! `ū`: it vanishes for a uniform velocity and it changes when only the
//! velocity changes at fixed strain.
//!
//! References: Germano, Piomelli, Moin & Cabot (1991), Phys. Fluids A 3, 1760;
//! Lilly (1992), Phys. Fluids A 4, 633.

use super::functions::Mat3;

/// Symmetric 3×3 tensor in Voigt order `(xx, yy, zz, xy, yz, xz)`.
type Sym6 = [f64; 6];

/// Relative threshold below which `⟨M:M⟩` is considered to have vanished.
///
/// It is relative to `(Δ² ⟨|S̄|²⟩)²`, the natural size of `M:M`, so that slow
/// flows are not mistaken for degenerate ones, while a vanishing strain or a
/// test filter equal to the grid filter (`α = 1`) are.
const DEGENERATE_RELATIVE_TOLERANCE: f64 = 1.0e-24;

fn sym6(m: &Mat3) -> Sym6 {
    [
        m[0][0],
        m[1][1],
        m[2][2],
        0.5 * (m[0][1] + m[1][0]),
        0.5 * (m[1][2] + m[2][1]),
        0.5 * (m[0][2] + m[2][0]),
    ]
}

/// `a_i a_j` as a symmetric tensor.
fn outer(a: &[f64; 3]) -> Sym6 {
    [
        a[0] * a[0],
        a[1] * a[1],
        a[2] * a[2],
        a[0] * a[1],
        a[1] * a[2],
        a[0] * a[2],
    ]
}

/// Double contraction `A : B`.
fn double_dot(a: &Sym6, b: &Sym6) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
}

/// Strain magnitude `|S| = sqrt(2 S_ij S_ij)`.
fn magnitude(s: &Sym6) -> f64 {
    (2.0 * double_dot(s, s)).sqrt()
}

/// Deviatoric part `A − tr(A) I/3`.
fn deviator(a: &Sym6) -> Sym6 {
    let mean = (a[0] + a[1] + a[2]) / 3.0;
    [a[0] - mean, a[1] - mean, a[2] - mean, a[3], a[4], a[5]]
}

/// `acc += weight · x`.
fn accumulate(acc: &mut Sym6, x: &Sym6, weight: f64) {
    for (a, v) in acc.iter_mut().zip(x.iter()) {
        *a += weight * v;
    }
}

/// Particle indices of the test-filter stencil of `k`: `k` itself followed by
/// its neighbours.  Neighbour indices that are out of range or equal to `k`
/// are skipped.
fn stencil(neighbors: &[Vec<usize>], k: usize, n: usize) -> impl Iterator<Item = usize> + '_ {
    std::iter::once(k).chain(
        neighbors
            .get(k)
            .into_iter()
            .flatten()
            .copied()
            .filter(move |&j| j < n && j != k),
    )
}

/// Local contractions of the Germano identity at one particle.
#[derive(Clone, Copy, Default)]
struct LocalTerms {
    /// `L^d : M^d`
    lm: f64,
    /// `M^d : M^d`
    mm: f64,
    /// `(Δ² ⟨|S̄|²⟩)²`, the natural size of `M:M`
    scale: f64,
}

/// Evaluate `L`, `M` and their contractions at particle `k` from the test-filter
/// stencil of `k`.
fn local_terms(
    u: &[[f64; 3]],
    s: &[Mat3],
    neighbors: &[Vec<usize>],
    k: usize,
    delta: f64,
    filter_ratio: f64,
) -> LocalTerms {
    let n = u.len().min(s.len());
    if k >= n {
        return LocalTerms::default();
    }
    let mut count = 0.0_f64;
    let mut u_hat = [0.0_f64; 3];
    let mut uu_hat: Sym6 = [0.0; 6];
    let mut s_hat: Sym6 = [0.0; 6];
    let mut strain_flux_hat: Sym6 = [0.0; 6];
    let mut magnitude_sq_hat = 0.0_f64;
    for j in stencil(neighbors, k, n) {
        let sj = sym6(&s[j]);
        let mag = magnitude(&sj);
        count += 1.0;
        for (acc, v) in u_hat.iter_mut().zip(u[j].iter()) {
            *acc += v;
        }
        accumulate(&mut uu_hat, &outer(&u[j]), 1.0);
        accumulate(&mut s_hat, &sj, 1.0);
        accumulate(&mut strain_flux_hat, &sj, mag);
        magnitude_sq_hat += mag * mag;
    }
    let inv = 1.0 / count;
    for v in uu_hat
        .iter_mut()
        .chain(s_hat.iter_mut())
        .chain(strain_flux_hat.iter_mut())
    {
        *v *= inv;
    }
    for v in &mut u_hat {
        *v *= inv;
    }
    magnitude_sq_hat *= inv;

    // L_ij = (u_i u_j)^ − û_i û_j
    let mut leonard = uu_hat;
    accumulate(&mut leonard, &outer(&u_hat), -1.0);

    // M_ij = 2 Δ² [ (|S| S_ij)^ − α² |Ŝ| Ŝ_ij ]
    let delta_sq = delta * delta;
    let test_weight = filter_ratio * filter_ratio * magnitude(&s_hat);
    let mut model_diff = strain_flux_hat;
    accumulate(&mut model_diff, &s_hat, -test_weight);
    for v in &mut model_diff {
        *v *= 2.0 * delta_sq;
    }

    let leonard_dev = deviator(&leonard);
    let model_dev = deviator(&model_diff);
    let strain_scale = delta_sq * magnitude_sq_hat;
    LocalTerms {
        lm: double_dot(&leonard_dev, &model_dev),
        mm: double_dot(&model_dev, &model_dev),
        scale: strain_scale * strain_scale,
    }
}

/// `⟨L:M⟩ / ⟨M:M⟩` for the averaged contractions of `terms`, or `0` when the
/// denominator has vanished (zero strain, `α = 1`, non-finite input).
fn lilly_ratio(terms: impl Iterator<Item = LocalTerms>) -> f64 {
    let (mut lm, mut mm, mut scale) = (0.0_f64, 0.0_f64, 0.0_f64);
    for t in terms {
        lm += t.lm;
        mm += t.mm;
        scale += t.scale;
    }
    if !mm.is_finite() || mm <= DEGENERATE_RELATIVE_TOLERANCE * scale {
        return 0.0;
    }
    let ratio = lm / mm;
    if ratio.is_finite() { ratio } else { 0.0 }
}

/// Unclipped Germano–Lilly ratio `⟨L:M⟩ / ⟨M:M⟩` at particle `i`.
///
/// * `u` – resolved velocity `ū` of every particle,
/// * `s` – resolved strain-rate tensor `S̄_ij` of every particle,
/// * `neighbors` – neighbour list of every particle (without the particle itself); the
///   test-filter stencil of particle `k` is `k` plus `neighbors[k]`,
/// * `delta` – grid filter width `Δ` (particle spacing),
/// * `test_filter_ratio` – `α = Δ̂/Δ` (typically 2).
///
/// The result is the least-squares value of `C_s²` over the stencil of `i`; it
/// is *not* restricted to be positive.  Clip it with the bounds of the model
/// (see `DynamicSmagorinsky::dynamic_cs_sq`).  Returns `0` when `M` vanishes
/// (no strain variation to calibrate against) or when `i` is out of range.
/// Evaluating one particle costs `O(k²)` for `k` neighbours; use
/// [`germano_lilly_ratio_field`] for all particles.
pub fn germano_lilly_ratio(
    u: &[[f64; 3]],
    s: &[Mat3],
    neighbors: &[Vec<usize>],
    i: usize,
    delta: f64,
    test_filter_ratio: f64,
) -> f64 {
    let n = u.len().min(s.len());
    if i >= n {
        return 0.0;
    }
    lilly_ratio(
        stencil(neighbors, i, n).map(|k| local_terms(u, s, neighbors, k, delta, test_filter_ratio)),
    )
}

/// [`germano_lilly_ratio`] for every particle (`O(N k)`).
///
/// The returned vector has `min(u.len(), s.len())` entries.
pub fn germano_lilly_ratio_field(
    u: &[[f64; 3]],
    s: &[Mat3],
    neighbors: &[Vec<usize>],
    delta: f64,
    test_filter_ratio: f64,
) -> Vec<f64> {
    let n = u.len().min(s.len());
    let terms: Vec<LocalTerms> = (0..n)
        .map(|k| local_terms(u, s, neighbors, k, delta, test_filter_ratio))
        .collect();
    (0..n)
        .map(|i| lilly_ratio(stencil(neighbors, i, n).map(|k| terms[k])))
        .collect()
}

/// Restrict `value` to `[lower, upper]` without panicking for inverted bounds
/// (the upper bound wins).
pub(super) fn clip_to_bounds(value: f64, lower: f64, upper: f64) -> f64 {
    value.max(lower).min(upper)
}

#[cfg(test)]
mod tests {
    use super::super::functions::{box_filter_scalar, box_filter_vector};
    use super::super::types::{DynamicSmagorinsky, DynamicSmagorinskyLes};
    use super::*;

    /// Index of lattice node `(ix, iy, iz)` in an `n³` lattice.
    fn node(n: usize, ix: usize, iy: usize, iz: usize) -> usize {
        (ix * n + iy) * n + iz
    }

    /// Positions of an `n³` lattice with spacing `dx`, centred on the origin.
    fn lattice_positions(n: usize, dx: f64) -> Vec<[f64; 3]> {
        let half = (n as f64 - 1.0) / 2.0;
        let mut points = Vec::with_capacity(n * n * n);
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    points.push([
                        (ix as f64 - half) * dx,
                        (iy as f64 - half) * dx,
                        (iz as f64 - half) * dx,
                    ]);
                }
            }
        }
        points
    }

    /// 26-point (Chebyshev distance 1) neighbour lists of an `n³` lattice.
    fn chebyshev_neighbors(n: usize) -> Vec<Vec<usize>> {
        let mut lists = Vec::with_capacity(n * n * n);
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let mut list = Vec::new();
                    for ox in [-1_isize, 0, 1] {
                        for oy in [-1_isize, 0, 1] {
                            for oz in [-1_isize, 0, 1] {
                                if (ox, oy, oz) == (0, 0, 0) {
                                    continue;
                                }
                                if let (Some(jx), Some(jy), Some(jz)) = (
                                    ix.checked_add_signed(ox),
                                    iy.checked_add_signed(oy),
                                    iz.checked_add_signed(oz),
                                ) && jx < n
                                    && jy < n
                                    && jz < n
                                {
                                    list.push(node(n, jx, jy, jz));
                                }
                            }
                        }
                    }
                    lists.push(list);
                }
            }
        }
        lists
    }

    fn apply(m: &Mat3, x: &[f64; 3]) -> [f64; 3] {
        [
            m[0][0] * x[0] + m[0][1] * x[1] + m[0][2] * x[2],
            m[1][0] * x[0] + m[1][1] * x[1] + m[1][2] * x[2],
            m[2][0] * x[0] + m[2][1] * x[1] + m[2][2] * x[2],
        ]
    }

    fn frobenius_sq(m: &Mat3) -> f64 {
        m.iter().flatten().map(|v| v * v).sum()
    }

    /// Deterministic pseudo-random numbers in `[-1, 1)`.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        }
    }

    fn random_velocity(count: usize, rng: &mut Lcg) -> Vec<[f64; 3]> {
        (0..count)
            .map(|_| [rng.next_f64(), rng.next_f64(), rng.next_f64()])
            .collect()
    }

    fn random_strain(count: usize, rng: &mut Lcg) -> Vec<Mat3> {
        (0..count)
            .map(|_| {
                let (xx, yy) = (rng.next_f64(), rng.next_f64());
                let zz = -(xx + yy);
                let (xy, yz, xz) = (rng.next_f64(), rng.next_f64(), rng.next_f64());
                [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
            })
            .collect()
    }

    /// Reference evaluation of `⟨L:M⟩/⟨M:M⟩` built from the public box filters.
    fn reference_ratio(
        u: &[[f64; 3]],
        s: &[Mat3],
        neighbors: &[Vec<usize>],
        i: usize,
        delta: f64,
        alpha: f64,
    ) -> f64 {
        let n = u.len();
        let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
        let mag: Vec<f64> = s.iter().map(|m| (2.0 * frobenius_sq(m)).sqrt()).collect();
        let u_hat = box_filter_vector(u, neighbors);
        let filtered = |values: Vec<f64>| box_filter_scalar(&values, neighbors);
        let uu: Vec<Vec<f64>> = pairs
            .iter()
            .map(|&(a, b)| filtered(u.iter().map(|v| v[a] * v[b]).collect()))
            .collect();
        let strain_flux: Vec<Vec<f64>> = pairs
            .iter()
            .map(|&(a, b)| filtered((0..n).map(|j| mag[j] * s[j][a][b]).collect()))
            .collect();
        let s_hat: Vec<Vec<f64>> = pairs
            .iter()
            .map(|&(a, b)| filtered(s.iter().map(|m| m[a][b]).collect()))
            .collect();
        let mut lm = vec![0.0; n];
        let mut mm = vec![0.0; n];
        for k in 0..n {
            let mut leonard = [[0.0; 3]; 3];
            let mut model = [[0.0; 3]; 3];
            let mut s_full = [[0.0; 3]; 3];
            for (c, &(a, b)) in pairs.iter().enumerate() {
                s_full[a][b] = s_hat[c][k];
                s_full[b][a] = s_hat[c][k];
            }
            let s_hat_mag = (2.0 * frobenius_sq(&s_full)).sqrt();
            for (c, &(a, b)) in pairs.iter().enumerate() {
                let l_ab = uu[c][k] - u_hat[k][a] * u_hat[k][b];
                let m_ab = 2.0
                    * delta
                    * delta
                    * (strain_flux[c][k] - alpha * alpha * s_hat_mag * s_hat[c][k]);
                leonard[a][b] = l_ab;
                leonard[b][a] = l_ab;
                model[a][b] = m_ab;
                model[b][a] = m_ab;
            }
            let l_mean = (leonard[0][0] + leonard[1][1] + leonard[2][2]) / 3.0;
            let m_mean = (model[0][0] + model[1][1] + model[2][2]) / 3.0;
            for a in 0..3 {
                leonard[a][a] -= l_mean;
                model[a][a] -= m_mean;
            }
            for (l_row, m_row) in leonard.iter().zip(model.iter()) {
                for (l, m) in l_row.iter().zip(m_row.iter()) {
                    lm[k] += l * m;
                    mm[k] += m * m;
                }
            }
        }
        let lm_avg = box_filter_scalar(&lm, neighbors);
        let mm_avg = box_filter_scalar(&mm, neighbors);
        lm_avg[i] / mm_avg[i]
    }

    /// Closed-form `C_s²` at the centre of a 5³ lattice for `u = B x`, `S̄ = A`
    /// (`A` symmetric and traceless, `B = A + W` with `W` antisymmetric):
    /// `L = (2/3) dx² B Bᵀ`, `M = 2 Δ² (1 − α²) |A| A`.
    fn closed_form_ratio(a: &Mat3, b: &Mat3, dx: f64, delta: f64, alpha: f64) -> f64 {
        let mut bbt_dot_a = 0.0;
        for (i, a_row) in a.iter().enumerate() {
            for (j, a_ij) in a_row.iter().enumerate() {
                let bbt_ij: f64 = (0..3).map(|k| b[i][k] * b[j][k]).sum();
                bbt_dot_a += bbt_ij * a_ij;
            }
        }
        let a_sq = frobenius_sq(a);
        let m_scale = 2.0 * delta * delta * (1.0 - alpha * alpha) * (2.0 * a_sq).sqrt();
        (2.0 / 3.0) * dx * dx * bbt_dot_a / (m_scale * a_sq)
    }

    /// Fields of the linear flow `u = B x`, `S̄ = A` on the 5³ lattice and the centre index.
    fn linear_flow(a: &Mat3, b: &Mat3, dx: f64) -> (Vec<[f64; 3]>, Vec<Mat3>, usize) {
        let n = 5;
        let positions = lattice_positions(n, dx);
        let u: Vec<[f64; 3]> = positions.iter().map(|x| apply(b, x)).collect();
        let s = vec![*a; u.len()];
        (u, s, node(n, 2, 2, 2))
    }

    fn assert_close(actual: f64, expected: f64, rel: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= rel * expected.abs().max(1e-300),
            "{what}: got {actual}, expected {expected}"
        );
    }

    #[test]
    fn constant_velocity_field_gives_zero_coefficient() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let u = vec![[0.8, -0.3, 1.7]; count];
        let mut rng = Lcg(3);
        // Constant velocity => no resolved stress, whatever the strain field is.
        let s = random_strain(count, &mut rng);
        let model = DynamicSmagorinsky::new(0.1, 2.0);
        for i in 0..count {
            let c = model.dynamic_cs_sq(&u, &s, &neighbors, i);
            assert!(c.abs() < 1e-12, "node {i}: C_s² = {c}");
        }
        // Zero strain => M vanishes: nothing to calibrate against.
        let zero = vec![[[0.0; 3]; 3]; count];
        let u_random = random_velocity(count, &mut rng);
        assert_eq!(model.dynamic_cs_sq(&u_random, &zero, &neighbors, 62), 0.0);
        assert!(
            model
                .dynamic_cs_sq_field(&u_random, &zero, &neighbors)
                .iter()
                .all(|&c| c == 0.0)
        );
    }

    #[test]
    fn different_random_fields_give_different_coefficients() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let model = DynamicSmagorinsky {
            c_max: 1.0e3,
            c_min: -1.0e3,
            ..DynamicSmagorinsky::new(0.1, 2.0)
        };
        let mut rng = Lcg(2024);
        let u_a = random_velocity(count, &mut rng);
        let s_a = random_strain(count, &mut rng);
        let u_b = random_velocity(count, &mut rng);
        let s_b = random_strain(count, &mut rng);
        let c_a = model.dynamic_cs_sq(&u_a, &s_a, &neighbors, 62);
        let c_b = model.dynamic_cs_sq(&u_b, &s_b, &neighbors, 62);
        assert!(c_a.is_finite() && c_b.is_finite());
        assert!(
            (c_a - c_b).abs() > 1e-6 * c_a.abs().max(c_b.abs()),
            "{c_a} vs {c_b}"
        );
        // At a fixed strain field the coefficient still depends on the velocity.
        let c_a_s = model.dynamic_cs_sq(&u_a, &s_b, &neighbors, 62);
        let c_b_s = model.dynamic_cs_sq(&u_b, &s_b, &neighbors, 62);
        assert!(
            (c_a_s - c_b_s).abs() > 1e-6 * c_a_s.abs().max(c_b_s.abs()),
            "{c_a_s} vs {c_b_s}"
        );
    }

    #[test]
    fn coefficient_field_is_not_a_constant() {
        // The strain-only formulation returned the same value (the upper clipping bound)
        // for every particle of every flow.
        let n = 6;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(404);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);

        // Unclipped (very wide bounds): essentially every particle gets its own value.
        let wide = DynamicSmagorinsky {
            c_min: -1.0e3,
            c_max: 1.0e3,
            ..DynamicSmagorinsky::new(0.1, 2.0)
        };
        let mut values = wide.dynamic_cs_sq_field(&u, &s, &neighbors);
        values.sort_by(|a, b| a.total_cmp(b));
        values.dedup();
        assert!(
            values.len() > count / 2,
            "only {} distinct coefficients for {count} particles",
            values.len()
        );

        // With the default bounds the field is clipped, but not saturated everywhere.
        let model = DynamicSmagorinsky::new(0.1, 2.0);
        let field = model.dynamic_cs_sq_field(&u, &s, &neighbors);
        assert!(field.iter().all(|&c| c >= model.c_min && c <= model.c_max));
        assert!(
            field.iter().any(|&c| c < model.c_max),
            "not every particle may sit on the upper clipping bound"
        );
        assert!(
            field.iter().any(|&c| c > model.c_min && c < model.c_max),
            "some particles must keep their unclipped value"
        );
    }

    #[test]
    fn linear_flow_matches_the_closed_form_and_lies_inside_the_clipping_bounds() {
        let dx = 0.1;
        let delta = 0.1;
        let alpha = 2.0;
        let model = DynamicSmagorinsky::new(delta, alpha);

        // Axisymmetric extension diag(1, 1, -2): C_s² > 0.  The strain-only formula
        // returned the upper clipping bound for every flow.
        let a: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -2.0]];
        let (u, s, centre) = linear_flow(&a, &a, dx);
        let neighbors = chebyshev_neighbors(5);
        let expected = closed_form_ratio(&a, &a, dx, delta, alpha);
        assert!(expected > 0.0 && expected < model.c_max, "{expected}");
        let ratio = germano_lilly_ratio(&u, &s, &neighbors, centre, delta, alpha);
        assert_close(ratio, expected, 1e-12, "unclipped ratio");
        let cs_sq = model.dynamic_cs_sq(&u, &s, &neighbors, centre);
        assert_close(cs_sq, expected, 1e-12, "C_s²");
        assert!(cs_sq > model.c_min && cs_sq < model.c_max);

        // Reversed flow (u -> -u, S -> -S): the sign of tr(A³) flips, so the least-squares
        // value is negative and is clipped to the lower bound.
        let a_rev = a.map(|row| row.map(|v| -v));
        let (u_rev, s_rev, _) = linear_flow(&a_rev, &a_rev, dx);
        let expected_rev = closed_form_ratio(&a_rev, &a_rev, dx, delta, alpha);
        assert_close(expected_rev, -expected, 1e-12, "closed form is odd in A");
        let ratio_rev = germano_lilly_ratio(&u_rev, &s_rev, &neighbors, centre, delta, alpha);
        assert_close(ratio_rev, expected_rev, 1e-12, "unclipped reversed ratio");
        assert_eq!(
            model.dynamic_cs_sq(&u_rev, &s_rev, &neighbors, centre),
            model.c_min
        );

        // Generic symmetric traceless strain plus a rotation: the rotation changes L
        // (the velocity) but not S̄, and the closed form follows it.
        let a_gen: Mat3 = [[1.0, 0.4, -0.3], [0.4, -0.6, 0.5], [-0.3, 0.5, -0.4]];
        let w: Mat3 = [[0.0, -0.7, 0.2], [0.7, 0.0, -0.4], [-0.2, 0.4, 0.0]];
        let b_gen: Mat3 = [
            [
                a_gen[0][0] + w[0][0],
                a_gen[0][1] + w[0][1],
                a_gen[0][2] + w[0][2],
            ],
            [
                a_gen[1][0] + w[1][0],
                a_gen[1][1] + w[1][1],
                a_gen[1][2] + w[1][2],
            ],
            [
                a_gen[2][0] + w[2][0],
                a_gen[2][1] + w[2][1],
                a_gen[2][2] + w[2][2],
            ],
        ];
        let wide = DynamicSmagorinsky {
            c_min: -1.0e3,
            c_max: 1.0e3,
            ..DynamicSmagorinsky::new(delta, alpha)
        };
        let (u_pure, s_pure, _) = linear_flow(&a_gen, &a_gen, dx);
        let (u_rot, s_rot, _) = linear_flow(&a_gen, &b_gen, dx);
        assert_eq!(s_pure, s_rot, "the strain field is the same");
        let c_pure = wide.dynamic_cs_sq(&u_pure, &s_pure, &neighbors, centre);
        let c_rot = wide.dynamic_cs_sq(&u_rot, &s_rot, &neighbors, centre);
        assert_close(
            c_pure,
            closed_form_ratio(&a_gen, &a_gen, dx, delta, alpha),
            1e-11,
            "pure strain",
        );
        assert_close(
            c_rot,
            closed_form_ratio(&a_gen, &b_gen, dx, delta, alpha),
            1e-11,
            "strain + rotation",
        );
        assert!(
            (c_pure - c_rot).abs() > 1e-3 * c_pure.abs(),
            "the velocity (not only the strain) must matter: {c_pure} vs {c_rot}"
        );
    }

    #[test]
    fn matches_a_reference_built_from_the_public_box_filters() {
        let n = 6;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(99);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let delta = 0.07;
        let alpha = 1.8;
        let field = germano_lilly_ratio_field(&u, &s, &neighbors, delta, alpha);
        assert_eq!(field.len(), count);
        for i in [
            node(n, 2, 2, 2),
            node(n, 3, 2, 4),
            node(n, 0, 0, 0),
            node(n, 5, 3, 1),
            77,
        ] {
            let reference = reference_ratio(&u, &s, &neighbors, i, delta, alpha);
            assert_close(field[i], reference, 1e-11, "field version");
            let single = germano_lilly_ratio(&u, &s, &neighbors, i, delta, alpha);
            assert_close(single, reference, 1e-11, "per-particle version");
        }
    }

    #[test]
    fn coefficient_is_invariant_under_scaling_and_galilean_shifts() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(5);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let reference = germano_lilly_ratio(&u, &s, &neighbors, 62, 0.1, 2.0);
        // u -> λ u, S̄ -> λ S̄ : L and M both scale with λ².
        let lambda = 3.7;
        let u_scaled: Vec<[f64; 3]> = u.iter().map(|v| v.map(|c| lambda * c)).collect();
        let s_scaled: Vec<Mat3> = s.iter().map(|m| m.map(|r| r.map(|c| lambda * c))).collect();
        let scaled = germano_lilly_ratio(&u_scaled, &s_scaled, &neighbors, 62, 0.1, 2.0);
        assert_close(scaled, reference, 1e-11, "scale invariance");
        // u -> u + U : the Leonard stress is a covariance and does not see a uniform shift.
        let shift = [4.0, -2.5, 0.75];
        let u_shifted: Vec<[f64; 3]> = u
            .iter()
            .map(|v| [v[0] + shift[0], v[1] + shift[1], v[2] + shift[2]])
            .collect();
        let shifted = germano_lilly_ratio(&u_shifted, &s, &neighbors, 62, 0.1, 2.0);
        assert_close(shifted, reference, 1e-9, "Galilean invariance");
    }

    #[test]
    fn isolated_particles_and_malformed_inputs_are_handled() {
        let u = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let s = [[[0.5, 0.2, 0.0], [0.2, -0.5, 0.0], [0.0, 0.0, 0.0]]; 2];
        // No neighbours: the stencil is the particle itself, so L = 0 and C = 0.
        let isolated: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
        assert_eq!(germano_lilly_ratio(&u, &s, &isolated, 0, 0.1, 2.0), 0.0);
        // Missing lists, self references and out-of-range indices are ignored.
        let malformed = [vec![0, 1, 9, 9]];
        let value = germano_lilly_ratio(&u, &s, &malformed, 0, 0.1, 2.0);
        assert!(value.is_finite());
        assert_eq!(germano_lilly_ratio(&u, &s, &malformed, 5, 0.1, 2.0), 0.0);
        let field = germano_lilly_ratio_field(&u, &s, &malformed, 0.1, 2.0);
        assert_eq!(field.len(), 2);
        assert!(field.iter().all(|c| c.is_finite()));
        // Mismatched array lengths use the common prefix.
        assert_eq!(
            germano_lilly_ratio_field(&u, &s[..1], &malformed, 0.1, 2.0).len(),
            1
        );
        assert!(germano_lilly_ratio_field(&[], &[], &[], 0.1, 2.0).is_empty());
    }

    #[test]
    fn test_filter_equal_to_the_grid_filter_is_degenerate() {
        // α = 1 makes M vanish identically.
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(8);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let ratio = germano_lilly_ratio(&u, &s, &neighbors, 62, 0.1, 1.0);
        assert!(ratio.is_finite());
        // M = 2Δ²((|S|S)^ − |Ŝ|Ŝ) does not vanish for a varying strain, but must for a uniform one.
        let uniform = vec![s[0]; count];
        assert_eq!(
            germano_lilly_ratio(&u, &uniform, &neighbors, 62, 0.1, 1.0),
            0.0
        );
    }

    #[test]
    fn clipping_never_panics_and_respects_the_bounds() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(17);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let unclipped = germano_lilly_ratio(&u, &s, &neighbors, 62, 0.1, 2.0);
        let capped = DynamicSmagorinsky {
            c_min: 0.0,
            c_max: 0.5 * unclipped.abs().max(1e-6),
            ..DynamicSmagorinsky::new(0.1, 2.0)
        };
        let c = capped.dynamic_cs_sq(&u, &s, &neighbors, 62);
        assert!(c >= capped.c_min && c <= capped.c_max, "{c}");
        // Inverted bounds must not panic (f64::clamp would).
        let inverted = DynamicSmagorinsky {
            c_min: 1.0,
            c_max: 0.1,
            ..DynamicSmagorinsky::new(0.1, 2.0)
        };
        assert_eq!(inverted.dynamic_cs_sq(&u, &s, &neighbors, 62), 0.1);
        let field = capped.dynamic_cs_sq_field(&u, &s, &neighbors);
        assert!(
            field
                .iter()
                .all(|&c| c >= capped.c_min && c <= capped.c_max)
        );
    }

    #[test]
    fn eddy_viscosity_and_sgs_stress_follow_the_dynamic_coefficient() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(31);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let delta = 0.05;
        let model = DynamicSmagorinsky {
            c_max: 10.0,
            ..DynamicSmagorinsky::new(delta, 2.0)
        };
        let cs_field = model.dynamic_cs_sq_field(&u, &s, &neighbors);
        let nu_field = model.eddy_viscosity_field(&u, &s, &neighbors);
        assert_eq!(nu_field.len(), count);
        let mut positive = 0;
        for i in 0..count {
            let mag = (2.0 * frobenius_sq(&s[i])).sqrt();
            assert_close(
                nu_field[i],
                cs_field[i] * delta * delta * mag,
                1e-14,
                "ν_t = C_s² Δ² |S̄|",
            );
            assert!(nu_field[i] >= 0.0, "ν_t must not be negative");
            if nu_field[i] > 0.0 {
                positive += 1;
            }
        }
        assert!(positive > 0, "random fields should activate the model");

        let density = 1000.0;
        let i = 62;
        let tau = model.sgs_stress_with_cs_sq(cs_field[i], &s[i], density);
        let mut dissipation = 0.0;
        for (a, row) in tau.iter().enumerate() {
            for (b, &t_ab) in row.iter().enumerate() {
                assert_close(t_ab, tau[b][a], 1e-14, "τ symmetric");
                assert_close(
                    t_ab,
                    -2.0 * density * nu_field[i] * s[i][a][b],
                    1e-14,
                    "τ = -2ρ ν_t S̄",
                );
                dissipation += t_ab * s[i][a][b];
            }
        }
        assert!(
            dissipation <= 0.0,
            "the SGS stress must drain resolved energy: {dissipation}"
        );
    }

    #[test]
    fn both_dynamic_models_share_the_same_procedure() {
        let n = 5;
        let neighbors = chebyshev_neighbors(n);
        let count = n * n * n;
        let mut rng = Lcg(12);
        let u = random_velocity(count, &mut rng);
        let s = random_strain(count, &mut rng);
        let wide = DynamicSmagorinsky {
            c_min: -1.0e3,
            c_max: 1.0e3,
            ..DynamicSmagorinsky::new(0.1, 2.0)
        };
        let les = DynamicSmagorinskyLes {
            cs_sq_min: -1.0e3,
            cs_sq_max: 1.0e3,
            ..DynamicSmagorinskyLes::new(0.1, 2.0)
        };
        for i in [0, 31, 62, 100] {
            assert_eq!(
                wide.dynamic_cs_sq(&u, &s, &neighbors, i),
                les.dynamic_cs_sq(&u, &s, &neighbors, i)
            );
        }
        assert_eq!(
            wide.dynamic_cs_sq_field(&u, &s, &neighbors),
            les.dynamic_cs_sq_field(&u, &s, &neighbors)
        );
        // The default clipping bounds differ (0.25 vs 0.04) and are honoured by each model.
        let default_les = DynamicSmagorinskyLes::default();
        for c in default_les.dynamic_cs_sq_field(&u, &s, &neighbors) {
            assert!((0.0..=0.04).contains(&c), "{c}");
        }
        let nu = les.eddy_viscosity_field(&u, &s, &neighbors);
        let cs = les.dynamic_cs_sq_field(&u, &s, &neighbors);
        for i in 0..count {
            let mag = (2.0 * frobenius_sq(&s[i])).sqrt();
            assert_close(nu[i], cs[i] * 0.1 * 0.1 * mag, 1e-14, "Les ν_t");
        }
        let tau = les.sgs_stress_with_cs_sq(cs[62], &s[62], 1000.0);
        assert_close(
            tau[0][1],
            -2.0 * 1000.0 * nu[62] * s[62][0][1],
            1e-14,
            "Les τ_xy",
        );
    }
}
