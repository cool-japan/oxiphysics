//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    GeodesicState, LeviCivitaConnection, MetricTensorND, OneForm, RicciTensor, RiemannTensor,
    RiemannianMetric, SecondFundamentalForm, So3, TwoForm,
};

/// A 3-vector.
pub type Vec3 = [f64; 3];
/// A 3×3 matrix stored row-major.
pub type Mat3 = [[f64; 3]; 3];
/// A 4×4 matrix stored row-major.
pub type Mat4 = [[f64; 4]; 4];
/// A unit quaternion \[w, x, y, z\].
pub type Quat = [f64; 4];
/// Add two 3-vectors.
pub fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
pub fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
pub fn scale3(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product of two 3-vectors.
pub fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
pub fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Norm of a 3-vector.
pub fn norm3(a: Vec3) -> f64 {
    dot3(a, a).sqrt()
}
/// Normalize a 3-vector (returns zero vector if near zero).
pub fn normalize3(a: Vec3) -> Vec3 {
    let n = norm3(a);
    if n < 1e-14 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(a, 1.0 / n)
    }
}
/// 3×3 identity matrix.
pub fn mat3_identity() -> Mat3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// 3×3 zero matrix.
pub fn mat3_zero() -> Mat3 {
    [[0.0; 3]; 3]
}
/// Matrix–vector multiply: M·v.
pub fn mat3_mul_vec3(m: Mat3, v: Vec3) -> Vec3 {
    [dot3(m[0], v), dot3(m[1], v), dot3(m[2], v)]
}
/// Matrix–matrix multiply: A·B.
pub fn mat3_mul(a: Mat3, b: Mat3) -> Mat3 {
    let bt = mat3_transpose(b);
    [
        [dot3(a[0], bt[0]), dot3(a[0], bt[1]), dot3(a[0], bt[2])],
        [dot3(a[1], bt[0]), dot3(a[1], bt[1]), dot3(a[1], bt[2])],
        [dot3(a[2], bt[0]), dot3(a[2], bt[1]), dot3(a[2], bt[2])],
    ]
}
/// Matrix transpose.
pub fn mat3_transpose(m: Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Matrix trace.
pub fn mat3_trace(m: Mat3) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}
/// Frobenius norm of a 3×3 matrix.
pub fn mat3_frobenius_norm(m: Mat3) -> f64 {
    let mut s = 0.0;
    for row in &m {
        for &v in row {
            s += v * v;
        }
    }
    s.sqrt()
}
/// Add two 3×3 matrices.
pub fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut r = mat3_zero();
    for (ri, (ai, bi)) in r.iter_mut().zip(a.iter().zip(b.iter())) {
        for (rij, (aij, bij)) in ri.iter_mut().zip(ai.iter().zip(bi.iter())) {
            *rij = *aij + *bij;
        }
    }
    r
}
/// Scale a 3×3 matrix by scalar.
pub fn mat3_scale(m: Mat3, s: f64) -> Mat3 {
    let mut r = mat3_zero();
    for (ri, mi) in r.iter_mut().zip(m.iter()) {
        for (rij, mij) in ri.iter_mut().zip(mi.iter()) {
            *rij = *mij * s;
        }
    }
    r
}
/// 3×3 determinant.
pub fn mat3_det(m: Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Skew-symmetric (hat) operator: converts a 3-vector ω to \[ω\]×.
pub fn skew3(v: Vec3) -> Mat3 {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}
/// Vee (inverse hat) operator: extracts 3-vector from skew-symmetric matrix.
pub fn vee3(m: Mat3) -> Vec3 {
    [m[2][1], m[0][2], m[1][0]]
}
/// Wedge product of two 1-forms: α ∧ β.
pub fn wedge(alpha: OneForm, beta: OneForm) -> TwoForm {
    let a = alpha.components;
    let b = beta.components;
    TwoForm {
        dxdy: a[0] * b[1] - a[1] * b[0],
        dydz: a[1] * b[2] - a[2] * b[1],
        dzdx: a[2] * b[0] - a[0] * b[2],
    }
}
/// Parallel transport of a vector along a curve on SO(3).
///
/// `initial_vec` – vector to be transported.
/// `angular_velocity` – angular velocity of the frame along the curve (rad/s).
/// `dt` – time step.
///
/// Returns the transported vector after one step.
pub fn parallel_transport_step(initial_vec: Vec3, angular_velocity: Vec3, dt: f64) -> Vec3 {
    let omega = scale3(angular_velocity, dt);
    let rot = So3::exp(omega);
    rot.rotate(initial_vec)
}
/// Parallel transport along a sequence of angular velocities.
pub fn parallel_transport(initial_vec: Vec3, angular_velocities: &[Vec3], dt: f64) -> Vec3 {
    let mut v = initial_vec;
    for &omega in angular_velocities {
        v = parallel_transport_step(v, omega, dt);
    }
    v
}
/// Riemannian exponential map on SO(3) at base point R₀.
///
/// `r0` – base point (rotation).
/// `tangent` – tangent vector in T_{R₀}SO(3) (Lie algebra element).
///
/// Returns the geodesic endpoint.
pub fn so3_exp_map(r0: &So3, tangent: Vec3) -> So3 {
    let step = So3::exp(tangent);
    r0.compose(&step)
}
/// Riemannian logarithm map on SO(3): log_{R₀}(R₁).
///
/// Returns the tangent vector pointing from R₀ to R₁ in T_{R₀}SO(3).
pub fn so3_log_map(r0: &So3, r1: &So3) -> Vec3 {
    let rel = r0.inverse().compose(r1);
    rel.log()
}
/// Geodesic on SO(3) between two rotations, sampled at parameter t ∈ \[0,1\].
pub fn so3_geodesic(r0: &So3, r1: &So3, t: f64) -> So3 {
    r0.slerp(r1, t)
}
/// Numerically approximate Christoffel symbols Γ^k_ij for a surface
/// defined by a mapping r(u, v) via finite differences.
///
/// `r_fn` – parametric surface function r(u, v) → R³.
/// `u`, `v` – evaluation point.
/// `h` – finite difference step size.
///
/// Returns Christoffel symbols of the second kind as a 2×2×2 array:
/// `gamma[k][i][j]` = Γ^k_ij with index 0 ↔ u and index 1 ↔ v.  They are
/// obtained from the first-kind symbols Γ_{ij,l} = r_ij · r_l by raising the
/// last index with the inverse first fundamental form, Γ^k_ij = g^{kl} Γ_{ij,l}.
pub fn christoffel_symbols_numerical<F>(r_fn: &F, u: f64, v: f64, h: f64) -> [[[f64; 2]; 2]; 2]
where
    F: Fn(f64, f64) -> Vec3,
{
    let ru = scale3(sub3(r_fn(u + h, v), r_fn(u - h, v)), 1.0 / (2.0 * h));
    let rv = scale3(sub3(r_fn(u, v + h), r_fn(u, v - h)), 1.0 / (2.0 * h));
    let r0 = r_fn(u, v);
    let ruu = scale3(
        sub3(add3(r_fn(u + h, v), r_fn(u - h, v)), scale3(r0, 2.0)),
        1.0 / (h * h),
    );
    let rvv = scale3(
        sub3(add3(r_fn(u, v + h), r_fn(u, v - h)), scale3(r0, 2.0)),
        1.0 / (h * h),
    );
    let ruv = scale3(
        sub3(
            add3(r_fn(u + h, v + h), r_fn(u - h, v - h)),
            add3(r_fn(u + h, v - h), r_fn(u - h, v + h)),
        ),
        1.0 / (4.0 * h * h),
    );
    let g = RiemannianMetric::from_tangents(ru, rv);
    let (g11, g12, g22) = g.inverse();
    // Christoffel symbols of the first kind of a parametrised surface:
    // Γ_{ij,k} = r_{ij} · r_k.  Spelled out in terms of the fundamental form
    // they are Γ_{11,1} = E_u/2, Γ_{11,2} = F_u - E_v/2, Γ_{12,1} = E_v/2,
    // Γ_{12,2} = G_u/2, Γ_{22,1} = F_v - G_u/2 and Γ_{22,2} = G_v/2.
    let c111 = dot3(ruu, ru);
    let c112 = dot3(ruu, rv);
    let c121 = dot3(ruv, ru);
    let c122 = dot3(ruv, rv);
    let c221 = dot3(rvv, ru);
    let c222 = dot3(rvv, rv);
    let gamma_1_11 = g11 * c111 + g12 * c112;
    let gamma_1_12 = g11 * c121 + g12 * c122;
    let gamma_1_22 = g11 * c221 + g12 * c222;
    let gamma_2_11 = g12 * c111 + g22 * c112;
    let gamma_2_12 = g12 * c121 + g22 * c122;
    let gamma_2_22 = g12 * c221 + g22 * c222;
    [
        [[gamma_1_11, gamma_1_12], [gamma_1_12, gamma_1_22]],
        [[gamma_2_11, gamma_2_12], [gamma_2_12, gamma_2_22]],
    ]
}
/// Integrate a geodesic on a surface by one step (Euler).
///
/// `state` – current geodesic state.
/// `gamma` – Christoffel symbols at current parameter point.
/// `ds` – arc-length step.
pub fn geodesic_step(state: GeodesicState, gamma: [[[f64; 2]; 2]; 2], ds: f64) -> GeodesicState {
    let u = state.u;
    let v = state.v;
    let du = state.du;
    let dv = state.dv;
    let ddu =
        -(gamma[0][0][0] * du * du + 2.0 * gamma[0][0][1] * du * dv + gamma[0][1][1] * dv * dv);
    let ddv =
        -(gamma[1][0][0] * du * du + 2.0 * gamma[1][0][1] * du * dv + gamma[1][1][1] * dv * dv);
    GeodesicState {
        u: u + du * ds,
        v: v + dv * ds,
        du: du + ddu * ds,
        dv: dv + ddv * ds,
    }
}
/// Compute a geodesic curve on a surface given initial conditions.
///
/// Returns list of parameter pairs (u, v).
pub fn geodesic_curve<F>(
    r_fn: &F,
    u0: f64,
    v0: f64,
    du0: f64,
    dv0: f64,
    num_steps: usize,
    ds: f64,
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64) -> Vec3,
{
    let mut state = GeodesicState::new(u0, v0, du0, dv0);
    let mut curve = Vec::with_capacity(num_steps);
    curve.push((state.u, state.v));
    let h = ds * 1e-3;
    for _ in 0..num_steps {
        let gamma = christoffel_symbols_numerical(r_fn, state.u, state.v, h);
        state = geodesic_step(state, gamma, ds);
        curve.push((state.u, state.v));
    }
    curve
}
/// Numerically integrate Gaussian curvature over a patch of a surface.
///
/// Uses a simple uniform grid for integration.
/// Returns the integral ∬K dA.
pub fn gauss_bonnet_integral<F>(
    r_fn: &F,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    nu: usize,
    nv: usize,
) -> f64
where
    F: Fn(f64, f64) -> Vec3,
{
    let du = (u_max - u_min) / nu as f64;
    let dv = (v_max - v_min) / nv as f64;
    let h = du.min(dv) * 1e-3;
    let mut integral = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = u_min + (i as f64 + 0.5) * du;
            let v = v_min + (j as f64 + 0.5) * dv;
            let ru = scale3(sub3(r_fn(u + h, v), r_fn(u - h, v)), 0.5 / h);
            let rv = scale3(sub3(r_fn(u, v + h), r_fn(u, v - h)), 0.5 / h);
            let normal_vec = cross3(ru, rv);
            let area = norm3(normal_vec);
            if area < 1e-15 {
                continue;
            }
            let unit_normal = scale3(normal_vec, 1.0 / area);
            let r0 = r_fn(u, v);
            let ruu = scale3(
                sub3(add3(r_fn(u + h, v), r_fn(u - h, v)), scale3(r0, 2.0)),
                1.0 / (h * h),
            );
            let rvv = scale3(
                sub3(add3(r_fn(u, v + h), r_fn(u, v - h)), scale3(r0, 2.0)),
                1.0 / (h * h),
            );
            let ruv = scale3(
                sub3(
                    add3(r_fn(u + h, v + h), r_fn(u - h, v - h)),
                    add3(r_fn(u + h, v - h), r_fn(u - h, v + h)),
                ),
                0.25 / (h * h),
            );
            let g_metric = RiemannianMetric::from_tangents(ru, rv);
            let sff = SecondFundamentalForm::from_derivatives(ruu, ruv, rvv, unit_normal);
            let k = sff.gaussian_curvature(&g_metric);
            integral += k * area * du * dv;
        }
    }
    integral
}
/// Compute the matrix exponential `exp(A)` of a general 3×3 matrix.
///
/// The exponential is evaluated by the scaling-and-squaring method with Padé
/// approximants (Higham; Al-Mohy & Higham) as implemented by
/// `nalgebra::Matrix3::exp`.  The scaling exponent grows only like
/// `log2(‖A‖)`, so matrices with large norm are handled to near full double
/// precision: for example `exp(diag(50, 0, 0))` reproduces `diag(e^50, 1, 1)`
/// and the exponential of a skew-symmetric matrix of any norm stays
/// orthogonal to rounding level.
///
/// The result overflows to `±inf`/`NaN` exactly as `f64::exp` does once the
/// true entries leave the representable range.  If any entry of `a` is not
/// finite, every entry of the result is `NaN`.
pub fn mat3_exp(a: Mat3) -> Mat3 {
    if a.iter().flatten().any(|value| !value.is_finite()) {
        return [[f64::NAN; 3]; 3];
    }
    let matrix = nalgebra::Matrix3::new(
        a[0][0], a[0][1], a[0][2], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2],
    );
    let exp = matrix.exp();
    [
        [exp[(0, 0)], exp[(0, 1)], exp[(0, 2)]],
        [exp[(1, 0)], exp[(1, 1)], exp[(1, 2)]],
        [exp[(2, 0)], exp[(2, 1)], exp[(2, 2)]],
    ]
}
/// Compute matrix logarithm of a rotation matrix (SO(3)) via Rodrigues formula.
pub fn so3_log_matrix(r: Mat3) -> Mat3 {
    let trace = mat3_trace(r);
    let cos_theta = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    if theta.abs() < 1e-12 {
        return mat3_zero();
    }
    let rt = mat3_transpose(r);

    mat3_scale(
        mat3_add(r, mat3_scale(rt, -1.0)),
        theta / (2.0 * theta.sin()),
    )
}
/// Lie bracket \[A, B\] = AB - BA for 3×3 matrices (elements of gl(3)).
pub fn lie_bracket_mat3(a: Mat3, b: Mat3) -> Mat3 {
    let ab = mat3_mul(a, b);
    let ba = mat3_mul(b, a);
    mat3_add(ab, mat3_scale(ba, -1.0))
}
/// Lie bracket for so(3) elements (skew-symmetric matrices): \[ω̂₁, ω̂₂\] = ω̂₁×ω̂₂ hat.
///
/// Equivalently: \[ω₁\]× × \[ω₂\]× = \[ω₁ × ω₂\]×.
pub fn so3_lie_bracket(omega1: Vec3, omega2: Vec3) -> Vec3 {
    cross3(omega1, omega2)
}
/// Adjoint action of SO(3) on so(3): Ad_R(ω) = R ω.
pub fn so3_adjoint(r: &So3, omega: Vec3) -> Vec3 {
    r.rotate(omega)
}
/// Coadjoint action of SO(3) on so(3)*: Ad*_R(μ) = R⁻ᵀ μ = R μ.
pub fn so3_coadjoint(r: &So3, mu: Vec3) -> Vec3 {
    r.rotate(mu)
}
/// Small adjoint (ad): ad_ω₁(ω₂) = ω₁ × ω₂.
pub fn so3_small_adjoint(omega1: Vec3, omega2: Vec3) -> Vec3 {
    cross3(omega1, omega2)
}
/// Maximum dimension supported for general metric tensor operations.
pub const MAX_DIM: usize = 4;
/// Ricci scalar R = g^{mu nu} R_{mu nu}.
pub fn ricci_scalar(metric: &MetricTensorND, ricci: &RicciTensor) -> f64 {
    let ginv = metric.inverse();
    let dim = metric.dim;
    let mut scalar = 0.0;
    for (mu, row) in ginv.iter().enumerate().take(dim) {
        for (nu, &g_val) in row.iter().enumerate().take(dim) {
            scalar += g_val * ricci.components[mu][nu];
        }
    }
    scalar
}
/// Compute the Ricci scalar directly from a metric function at a point.
pub fn ricci_scalar_from_metric_fn<F>(
    metric_fn: &F,
    point: &[f64; MAX_DIM],
    dim: usize,
    h: f64,
) -> f64
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let ricci = RicciTensor::from_metric_fn(metric_fn, point, dim, h);
    let g_here = metric_fn(point);
    let metric = MetricTensorND { dim, g: g_here };
    ricci_scalar(&metric, &ricci)
}
/// Einstein tensor G_{mu nu} = R_{mu nu} - 0.5 g_{mu nu} R.
pub fn einstein_tensor(metric: &MetricTensorND, ricci: &RicciTensor) -> [[f64; MAX_DIM]; MAX_DIM] {
    let r = ricci_scalar(metric, ricci);
    let dim = metric.dim;
    let mut g_tensor = [[0.0_f64; MAX_DIM]; MAX_DIM];
    for (mu, row) in g_tensor.iter_mut().enumerate().take(dim) {
        for (nu, cell) in row.iter_mut().enumerate().take(dim) {
            *cell = ricci.components[mu][nu] - 0.5 * metric.g[mu][nu] * r;
        }
    }
    g_tensor
}
/// Kretschner scalar K = R_{abcd} R^{abcd} (fully contracted Riemann tensor).
pub fn kretschner_scalar<F>(metric_fn: &F, point: &[f64; MAX_DIM], dim: usize, h: f64) -> f64
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let riemann = RiemannTensor::from_metric_fn(metric_fn, point, dim, h);
    let g_here = metric_fn(point);
    let metric = MetricTensorND { dim, g: g_here };
    let ginv = metric.inverse();
    let mut r_lower = [[[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM];
    for (a, ra) in r_lower.iter_mut().enumerate().take(dim) {
        for (b, rab) in ra.iter_mut().enumerate().take(dim) {
            for (c, rabc) in rab.iter_mut().enumerate().take(dim) {
                for (d, cell) in rabc.iter_mut().enumerate().take(dim) {
                    let mut val = 0.0;
                    for (e, &ge) in g_here[a].iter().enumerate().take(dim) {
                        val += ge * riemann.components[e][b][c][d];
                    }
                    *cell = val;
                }
            }
        }
    }
    let mut kretschner = 0.0;
    for a in 0..dim {
        for b in 0..dim {
            for c in 0..dim {
                for d in 0..dim {
                    let mut r_upper = 0.0;
                    for e in 0..dim {
                        for f in 0..dim {
                            for gg in 0..dim {
                                for hh in 0..dim {
                                    r_upper += ginv[a][e]
                                        * ginv[b][f]
                                        * ginv[c][gg]
                                        * ginv[d][hh]
                                        * r_lower[e][f][gg][hh];
                                }
                            }
                        }
                    }
                    kretschner += r_lower[a][b][c][d] * r_upper;
                }
            }
        }
    }
    kretschner
}
/// Covariant derivative of a vector field.
///
/// Given V^mu and its partial derivatives dV^mu/dx^nu at a point,
/// returns (nabla_nu V)^mu = dV^mu/dx^nu + Gamma^mu_{nu lambda} V^lambda.
pub fn covariant_derivative_vector(
    conn: &LeviCivitaConnection,
    v: &[f64; MAX_DIM],
    dv: &[[f64; MAX_DIM]; MAX_DIM],
) -> [[f64; MAX_DIM]; MAX_DIM] {
    let dim = conn.dim;
    let mut result = [[0.0_f64; MAX_DIM]; MAX_DIM];
    for mu in 0..dim {
        for nu in 0..dim {
            let mut val = dv[mu][nu];
            for (lambda, &vl) in v.iter().enumerate().take(dim) {
                val += conn.christoffel[mu][nu][lambda] * vl;
            }
            result[mu][nu] = val;
        }
    }
    result
}
/// Covariant derivative of a covector (1-form) field.
///
/// Given omega_mu and its partial derivatives d omega_mu/dx^nu,
/// returns (nabla_nu omega)_mu = d omega_mu/dx^nu - Gamma^lambda_{nu mu} omega_lambda.
pub fn covariant_derivative_covector(
    conn: &LeviCivitaConnection,
    omega: &[f64; MAX_DIM],
    domega: &[[f64; MAX_DIM]; MAX_DIM],
) -> [[f64; MAX_DIM]; MAX_DIM] {
    let dim = conn.dim;
    let mut result = [[0.0_f64; MAX_DIM]; MAX_DIM];
    for mu in 0..dim {
        for nu in 0..dim {
            let mut val = domega[mu][nu];
            for (lambda, &wl) in omega.iter().enumerate().take(dim) {
                val -= conn.christoffel[lambda][nu][mu] * wl;
            }
            result[mu][nu] = val;
        }
    }
    result
}
/// Configuration for geodesic integration.
#[derive(Debug, Clone, Copy)]
pub struct GeodesicConfig {
    /// Manifold dimension.
    pub dim: usize,
    /// Time step for RK4 integration.
    pub dt: f64,
    /// Number of integration steps.
    pub steps: usize,
    /// Finite-difference step for Christoffel symbol computation.
    pub h_christoffel: f64,
}

/// Solve the geodesic equation in N dimensions.
///
/// x''^\mu + Gamma^\mu_{\alpha\beta} x'^\alpha x'^\beta = 0
///
/// Uses RK4 integration.
///
/// Returns a vector of (position, velocity) pairs along the geodesic.
pub fn geodesic_equation_solve<F>(
    metric_fn: &F,
    x0: &[f64; MAX_DIM],
    v0: &[f64; MAX_DIM],
    cfg: &GeodesicConfig,
) -> Vec<([f64; MAX_DIM], [f64; MAX_DIM])>
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let dim = cfg.dim;
    let dt = cfg.dt;
    let steps = cfg.steps;
    let h_christoffel = cfg.h_christoffel;
    let mut trajectory = Vec::with_capacity(steps + 1);
    let mut x = *x0;
    let mut v = *v0;
    trajectory.push((x, v));
    for _ in 0..steps {
        let accel = |pos: &[f64; MAX_DIM], vel: &[f64; MAX_DIM]| -> [f64; MAX_DIM] {
            let conn = LeviCivitaConnection::from_metric_fn(metric_fn, pos, dim, h_christoffel);
            let mut a = [0.0; MAX_DIM];
            for (mu, cell) in a.iter_mut().enumerate().take(dim) {
                let mut val = 0.0;
                for alpha in 0..dim {
                    for beta in 0..dim {
                        val += conn.christoffel[mu][alpha][beta] * vel[alpha] * vel[beta];
                    }
                }
                *cell = -val;
            }
            a
        };
        let a1 = accel(&x, &v);
        let mut x1 = [0.0; MAX_DIM];
        let mut v1 = [0.0; MAX_DIM];
        for i in 0..dim {
            x1[i] = x[i] + 0.5 * dt * v[i];
            v1[i] = v[i] + 0.5 * dt * a1[i];
        }
        let a2 = accel(&x1, &v1);
        let mut x2 = [0.0; MAX_DIM];
        let mut v2 = [0.0; MAX_DIM];
        for i in 0..dim {
            x2[i] = x[i] + 0.5 * dt * v1[i];
            v2[i] = v[i] + 0.5 * dt * a2[i];
        }
        let a3 = accel(&x2, &v2);
        let mut x3 = [0.0; MAX_DIM];
        let mut v3 = [0.0; MAX_DIM];
        for i in 0..dim {
            x3[i] = x[i] + dt * v2[i];
            v3[i] = v[i] + dt * a3[i];
        }
        let a4 = accel(&x3, &v3);
        for i in 0..dim {
            x[i] += dt / 6.0 * (v[i] + 2.0 * v1[i] + 2.0 * v2[i] + v3[i]);
            v[i] += dt / 6.0 * (a1[i] + 2.0 * a2[i] + 2.0 * a3[i] + a4[i]);
        }
        trajectory.push((x, v));
    }
    trajectory
}
/// Parallel transport a vector along a curve in N dimensions.
///
/// `curve` is a sequence of coordinate points.
/// `v0` is the initial vector to transport.
/// `metric_fn` returns the metric at any point.
///
/// Returns the transported vector at each point of the curve.
pub fn parallel_transport_nd<F>(
    metric_fn: &F,
    dim: usize,
    curve: &[[f64; MAX_DIM]],
    v0: &[f64; MAX_DIM],
    h_christoffel: f64,
) -> Vec<[f64; MAX_DIM]>
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let mut result = Vec::with_capacity(curve.len());
    let mut v = *v0;
    result.push(v);
    for i in 1..curve.len() {
        let conn =
            LeviCivitaConnection::from_metric_fn(metric_fn, &curve[i - 1], dim, h_christoffel);
        let mut dx = [0.0; MAX_DIM];
        for k in 0..dim {
            dx[k] = curve[i][k] - curve[i - 1][k];
        }
        let mut dv = [0.0; MAX_DIM];
        for (mu, dv_mu) in dv.iter_mut().enumerate().take(dim) {
            for (alpha, &v_alpha) in v.iter().enumerate().take(dim) {
                for (beta, &dx_beta) in dx.iter().enumerate().take(dim) {
                    *dv_mu -= conn.christoffel[mu][alpha][beta] * v_alpha * dx_beta;
                }
            }
        }
        for mu in 0..dim {
            v[mu] += dv[mu];
        }
        result.push(v);
    }
    result
}
/// Check if a vector field is a Killing vector field.
///
/// A Killing vector field satisfies: nabla_mu xi_nu + nabla_nu xi_mu = 0.
///
/// `xi_fn` returns the Killing vector at a point.
/// Returns the maximum violation of the Killing equation.
pub fn killing_equation_violation<F, G>(
    metric_fn: &F,
    xi_fn: &G,
    point: &[f64; MAX_DIM],
    dim: usize,
    h: f64,
) -> f64
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
    G: Fn(&[f64; MAX_DIM]) -> [f64; MAX_DIM],
{
    let conn = LeviCivitaConnection::from_metric_fn(metric_fn, point, dim, h);
    let g_here = metric_fn(point);
    let xi = xi_fn(point);
    let mut xi_lower = [0.0; MAX_DIM];
    for mu in 0..dim {
        for nu in 0..dim {
            xi_lower[mu] += g_here[mu][nu] * xi[nu];
        }
    }
    let mut dxi_lower = [[0.0_f64; MAX_DIM]; MAX_DIM];
    for nu in 0..dim {
        let mut p_plus = *point;
        let mut p_minus = *point;
        p_plus[nu] += h;
        p_minus[nu] -= h;
        let xi_plus = xi_fn(&p_plus);
        let xi_minus = xi_fn(&p_minus);
        let g_plus = metric_fn(&p_plus);
        let g_minus = metric_fn(&p_minus);
        for mu in 0..dim {
            let mut xi_low_plus = 0.0;
            let mut xi_low_minus = 0.0;
            for rho in 0..dim {
                xi_low_plus += g_plus[mu][rho] * xi_plus[rho];
                xi_low_minus += g_minus[mu][rho] * xi_minus[rho];
            }
            dxi_lower[mu][nu] = (xi_low_plus - xi_low_minus) / (2.0 * h);
        }
    }
    let nabla_xi = covariant_derivative_covector(&conn, &xi_lower, &dxi_lower);
    let mut max_viol = 0.0_f64;
    for (mu, row_mu) in nabla_xi.iter().enumerate().take(dim) {
        for (nu, _) in row_mu.iter().enumerate().take(dim) {
            let viol = (nabla_xi[nu][mu] + nabla_xi[mu][nu]).abs();
            max_viol = max_viol.max(viol);
        }
    }
    max_viol
}
/// Weyl tensor C^rho_{sigma mu nu} for dim >= 3.
///
/// C = R - (2/(n-2))(g * Ric - g * g * R/(n-1)) where * denotes Kulkarni-Nomizu product
/// Simplified direct formula:
/// C^rho_{sigma mu nu} = R^rho_{sigma mu nu}
///   - 1/(n-2) (delta^rho_mu R_{sigma nu} - delta^rho_nu R_{sigma mu}
///     + g_{sigma nu} R^rho_mu - g_{sigma mu} R^rho_nu)
///   + R/((n-1)(n-2)) (delta^rho_mu g_{sigma nu} - delta^rho_nu g_{sigma mu})
pub fn weyl_tensor<F>(
    metric_fn: &F,
    point: &[f64; MAX_DIM],
    dim: usize,
    h: f64,
) -> [[[[f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM]
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    assert!(dim >= 3, "Weyl tensor only defined for dim >= 3");
    let riemann = RiemannTensor::from_metric_fn(metric_fn, point, dim, h);
    let ricci = RicciTensor::from_riemann(&riemann);
    let g_here = metric_fn(point);
    let metric = MetricTensorND { dim, g: g_here };
    let ginv = metric.inverse();
    let r = ricci_scalar(&metric, &ricci);
    let mut ricci_mixed = [[0.0_f64; MAX_DIM]; MAX_DIM];
    for (rho, rm_row) in ricci_mixed.iter_mut().enumerate().take(dim) {
        for (mu, cell) in rm_row.iter_mut().enumerate().take(dim) {
            for (sigma, &ginv_val) in ginv[rho].iter().enumerate().take(dim) {
                *cell += ginv_val * ricci.components[sigma][mu];
            }
        }
    }
    let n = dim as f64;
    let mut weyl = [[[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM];
    let delta = |i: usize, j: usize| -> f64 { if i == j { 1.0 } else { 0.0 } };
    for rho in 0..dim {
        for sigma in 0..dim {
            for mu in 0..dim {
                for nu in 0..dim {
                    let ricci_term = delta(rho, mu) * ricci.components[sigma][nu]
                        - delta(rho, nu) * ricci.components[sigma][mu]
                        + g_here[sigma][nu] * ricci_mixed[rho][mu]
                        - g_here[sigma][mu] * ricci_mixed[rho][nu];
                    let scalar_term = r
                        * (delta(rho, mu) * g_here[sigma][nu] - delta(rho, nu) * g_here[sigma][mu]);
                    weyl[rho][sigma][mu][nu] = riemann.components[rho][sigma][mu][nu]
                        - ricci_term / (n - 2.0)
                        + scalar_term / ((n - 1.0) * (n - 2.0));
                }
            }
        }
    }
    weyl
}
/// Configuration for geodesic deviation (Jacobi equation) integration.
#[derive(Debug, Clone, Copy)]
pub struct DeviationConfig {
    /// Manifold dimension.
    pub dim: usize,
    /// Time step (must match the geodesic time step).
    pub dt: f64,
    /// Finite-difference step for Riemann tensor computation.
    pub h: f64,
}

/// Geodesic deviation equation (Jacobi equation).
///
/// Given a geodesic (x(t), v(t)) and an initial deviation vector xi(0) and its derivative,
/// computes the evolution of the deviation vector.
///
/// D^2 xi^mu / dt^2 = -R^mu_{alpha beta gamma} v^alpha xi^beta v^gamma
pub fn geodesic_deviation<F>(
    metric_fn: &F,
    geodesic: &[([f64; MAX_DIM], [f64; MAX_DIM])],
    xi0: &[f64; MAX_DIM],
    dxi0: &[f64; MAX_DIM],
    cfg: &DeviationConfig,
) -> Vec<[f64; MAX_DIM]>
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let dim = cfg.dim;
    let dt = cfg.dt;
    let h = cfg.h;
    let mut result = Vec::with_capacity(geodesic.len());
    let mut xi = *xi0;
    let mut dxi = *dxi0;
    result.push(xi);
    for i in 1..geodesic.len() {
        let (pos, vel) = &geodesic[i - 1];
        let riemann = RiemannTensor::from_metric_fn(metric_fn, pos, dim, h);
        let mut ddxi = [0.0; MAX_DIM];
        for (mu, cell) in ddxi.iter_mut().enumerate().take(dim) {
            for alpha in 0..dim {
                for (beta, &xi_b) in xi.iter().enumerate().take(dim) {
                    for gamma in 0..dim {
                        *cell -= riemann.components[mu][alpha][beta][gamma]
                            * vel[alpha]
                            * xi_b
                            * vel[gamma];
                    }
                }
            }
        }
        for mu in 0..dim {
            xi[mu] += dxi[mu] * dt;
            dxi[mu] += ddxi[mu] * dt;
        }
        result.push(xi);
    }
    result
}
/// Compute the Gauss-Bonnet integrand in N=2 dimensions.
///
/// For a 2D surface with metric g_{ij}, the Gauss-Bonnet theorem states:
/// integral(K * sqrt(det g) du dv) = 2*pi*chi
///
/// where K is the Gaussian curvature and chi is the Euler characteristic.
pub fn gauss_bonnet_2d_integrand<F>(metric_fn: &F, point: &[f64; MAX_DIM], h: f64) -> f64
where
    F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
{
    let g_here = metric_fn(point);
    let metric = MetricTensorND { dim: 2, g: g_here };
    let r = ricci_scalar_from_metric_fn(metric_fn, point, 2, h);
    let k = r / 2.0;
    k * metric.volume_element()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::differential_geometry::Se3;
    use crate::differential_geometry::Se3Algebra;
    use std::f64::consts::PI;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_so3_identity_exp_zero() {
        let r = So3::exp([0.0, 0.0, 0.0]);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((r.mat[i][j] - expected).abs() < EPS);
            }
        }
    }
    #[test]
    fn test_so3_exp_log_roundtrip() {
        let omega = [0.3, -0.5, 0.2];
        let r = So3::exp(omega);
        let omega2 = r.log();
        for i in 0..3 {
            assert!(
                (omega[i] - omega2[i]).abs() < 1e-10,
                "Component {} mismatch",
                i
            );
        }
    }
    #[test]
    fn test_so3_orthogonality() {
        let r = So3::exp([0.5, -0.3, 0.7]);
        let rt = mat3_transpose(r.mat);
        let rrt = mat3_mul(r.mat, rt);
        let i = mat3_identity();
        for row in 0..3 {
            for col in 0..3 {
                assert!((rrt[row][col] - i[row][col]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_so3_determinant_one() {
        let r = So3::exp([1.0, 0.5, -0.3]);
        let d = mat3_det(r.mat);
        assert!((d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_so3_inverse_compose_identity() {
        let r = So3::exp([0.4, -0.2, 0.8]);
        let r_inv = r.inverse();
        let prod = r.compose(&r_inv);
        let omega = prod.log();
        assert!(norm3(omega) < 1e-10);
    }
    #[test]
    fn test_so3_slerp_endpoints() {
        let r0 = So3::identity();
        let r1 = So3::exp([0.3, 0.1, -0.2]);
        let s0 = r0.slerp(&r1, 0.0);
        let omega_s0 = s0.log();
        assert!(norm3(omega_s0) < 1e-10);
    }
    #[test]
    fn test_so3_slerp_midpoint() {
        let r0 = So3::identity();
        let omega = [0.6, 0.0, 0.0];
        let r1 = So3::exp(omega);
        let mid = r0.slerp(&r1, 0.5);
        let log_mid = mid.log();
        assert!((log_mid[0] - 0.3).abs() < 1e-10);
    }
    #[test]
    fn test_so3_quaternion_roundtrip() {
        let r = So3::exp([0.5, -0.3, 0.1]);
        let q = r.to_quaternion();
        let r2 = So3::from_quaternion(q);
        for i in 0..3 {
            for j in 0..3 {
                assert!((r.mat[i][j] - r2.mat[i][j]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_so3_geodesic_distance_zero_to_self() {
        let r = So3::exp([0.3, -0.1, 0.5]);
        assert!(r.dist(&r) < 1e-10);
    }
    #[test]
    fn test_so3_adjoint_equals_rotation() {
        let r = So3::exp([0.2, 0.3, -0.1]);
        let v = [1.0, 0.0, 0.0];
        let adj = so3_adjoint(&r, v);
        let rot = r.rotate(v);
        for i in 0..3 {
            assert!((adj[i] - rot[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_se3_identity() {
        let id = Se3::identity();
        let p = [1.0, 2.0, 3.0];
        let t = id.transform_point(p);
        for i in 0..3 {
            assert!((t[i] - p[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_se3_exp_log_roundtrip() {
        let twist = Se3Algebra::new([0.1, -0.2, 0.3], [1.0, -0.5, 0.2]);
        let se3 = Se3::exp(twist);
        let twist2 = se3.log();
        for i in 0..3 {
            assert!((twist.omega[i] - twist2.omega[i]).abs() < 1e-8);
        }
    }
    #[test]
    fn test_se3_inverse() {
        let se3 = Se3::from_rt(So3::exp([0.3, 0.1, -0.2]), [1.0, 2.0, 3.0]);
        let inv = se3.inverse();
        let composed = se3.compose(&inv);
        let t = composed.translation;
        for &ti in t.iter() {
            assert!(ti.abs() < 1e-10);
        }
    }
    #[test]
    fn test_se3_compose_translation() {
        let t1 = Se3::from_rt(So3::identity(), [1.0, 0.0, 0.0]);
        let t2 = Se3::from_rt(So3::identity(), [2.0, 0.0, 0.0]);
        let t12 = t1.compose(&t2);
        assert!((t12.translation[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_se3_mat4_last_row() {
        let se3 = Se3::from_rt(So3::exp([0.1, 0.2, -0.3]), [5.0, -2.0, 1.0]);
        let m = se3.to_mat4();
        assert!((m[3][0]).abs() < 1e-10);
        assert!((m[3][1]).abs() < 1e-10);
        assert!((m[3][2]).abs() < 1e-10);
        assert!((m[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_metric_sphere_area_element() {
        let r = 2.0_f64;
        let theta = std::f64::consts::PI / 2.0;
        let g = RiemannianMetric::new(r * r, 0.0, r * r * theta.sin() * theta.sin());
        let area = g.area_element();
        assert!((area - r * r).abs() < 1e-10);
    }
    #[test]
    fn test_metric_flat_plane() {
        let g = RiemannianMetric::new(1.0, 0.0, 1.0);
        assert!((g.det() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_metric_length_orthogonal() {
        let g = RiemannianMetric::new(4.0, 0.0, 9.0);
        assert!((g.length(1.0, 0.0) - 2.0).abs() < 1e-10);
        assert!((g.length(0.0, 1.0) - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_gauss_curvature_sphere() {
        let r = 3.0_f64;
        let ru = [r, 0.0, 0.0];
        let rv = [0.0, r, 0.0];
        let ruu = [-r, 0.0, 0.0];
        let ruv = [0.0, 0.0, 0.0];
        let rvv = [0.0, -r, 0.0];
        let n = [0.0, 0.0, 1.0];
        let g = RiemannianMetric::from_tangents(ru, rv);
        let sff = SecondFundamentalForm::from_derivatives(ruu, ruv, rvv, n);
        let k = sff.gaussian_curvature(&g);
        let expected = (sff.l * sff.n - sff.m * sff.m) / g.det();
        assert!((k - expected).abs() < 1e-10);
    }
    #[test]
    fn test_mean_curvature_flat_surface() {
        let g = RiemannianMetric::new(1.0, 0.0, 1.0);
        let sff = SecondFundamentalForm {
            l: 0.0,
            m: 0.0,
            n: 0.0,
        };
        assert_eq!(sff.mean_curvature(&g), 0.0);
    }
    #[test]
    fn test_one_form_evaluation() {
        let alpha = OneForm::new([1.0, 2.0, 3.0]);
        let v = [4.0, 5.0, 6.0];
        assert!((alpha.evaluate(v) - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_two_form_evaluation_antisymmetry() {
        let alpha = OneForm::new([1.0, 0.0, 0.0]);
        let beta = OneForm::new([0.0, 1.0, 0.0]);
        let form = wedge(alpha, beta);
        let u = [1.0, 0.0, 0.0];
        let v = [0.0, 1.0, 0.0];
        let val1 = form.evaluate(u, v);
        let val2 = form.evaluate(v, u);
        assert!((val1 + val2).abs() < 1e-10);
    }
    #[test]
    fn test_wedge_product_basis() {
        let ex = OneForm::new([1.0, 0.0, 0.0]);
        let ey = OneForm::new([0.0, 1.0, 0.0]);
        let form = wedge(ex, ey);
        assert!((form.dxdy - 1.0).abs() < 1e-10);
        assert!((form.dydz).abs() < 1e-10);
        assert!((form.dzdx).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_transport_zero_rotation_identity() {
        let v = [1.0, 0.0, 0.0];
        let transported = parallel_transport(v, &[[0.0, 0.0, 0.0]; 10], 0.01);
        for i in 0..3 {
            assert!((transported[i] - v[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_parallel_transport_preserves_length() {
        let v = [1.0, 2.0, 3.0];
        let init_len = norm3(v);
        let omegas: Vec<Vec3> = (0..20).map(|i| [0.01 * i as f64 * 0.1, 0.0, 0.0]).collect();
        let transported = parallel_transport(v, &omegas, 0.01);
        let final_len = norm3(transported);
        assert!((final_len - init_len).abs() < 1e-8);
    }
    #[test]
    fn test_mat3_exp_identity_at_zero() {
        let a = mat3_zero();
        let e = mat3_exp(a);
        let i = mat3_identity();
        for row in 0..3 {
            for col in 0..3 {
                assert!((e[row][col] - i[row][col]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_mat3_exp_so3_agrees() {
        let omega = [0.3, -0.2, 0.1];
        let r_rodrigues = So3::exp(omega);
        let skew = skew3(omega);
        let r_exp = mat3_exp(skew);
        for (i, (rod_row, exp_row)) in r_rodrigues.mat.iter().zip(r_exp.iter()).enumerate() {
            for (j, (&rod_val, &exp_val)) in rod_row.iter().zip(exp_row.iter()).enumerate() {
                assert!((rod_val - exp_val).abs() < 1e-8, "mismatch at [{i}][{j}]");
            }
        }
    }
    #[test]
    fn test_lie_bracket_antisymmetric() {
        let a = skew3([1.0, 0.0, 0.0]);
        let b = skew3([0.0, 1.0, 0.0]);
        let ab = lie_bracket_mat3(a, b);
        let ba = lie_bracket_mat3(b, a);
        for (ab_row, ba_row) in ab.iter().zip(ba.iter()) {
            for (&ab_val, &ba_val) in ab_row.iter().zip(ba_row.iter()) {
                assert!((ab_val + ba_val).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_so3_lie_bracket_cross_product() {
        let o1 = [1.0, 0.0, 0.0];
        let o2 = [0.0, 1.0, 0.0];
        let bracket = so3_lie_bracket(o1, o2);
        let cross = cross3(o1, o2);
        for i in 0..3 {
            assert!((bracket[i] - cross[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_geodesic_curve_flat_plane() {
        let r_fn = |u: f64, v: f64| [u, v, 0.0];
        let curve = geodesic_curve(&r_fn, 0.0, 0.0, 1.0, 0.0, 10, 0.1);
        assert_eq!(curve.len(), 11);
        for (u, v) in &curve {
            assert!(
                v.abs() < 1e-8,
                "v should be ~0 on flat plane geodesic, got {}",
                v
            );
            let _ = u;
        }
    }
    #[test]
    fn test_gauss_bonnet_sphere() {
        let r = 1.0_f64;
        let r_fn = |u: f64, v: f64| {
            let x = r * u.sin() * v.cos();
            let y = r * u.sin() * v.sin();
            let z = r * u.cos();
            [x, y, z]
        };
        let integral = gauss_bonnet_integral(&r_fn, 0.05, PI - 0.05, 0.0, 2.0 * PI, 20, 40);
        assert!(
            (integral - 4.0 * PI).abs() < 0.5,
            "Gauss-Bonnet integral = {}, expected ~4π",
            integral
        );
    }
    #[test]
    fn test_skew_vee_roundtrip() {
        let v = [1.0, 2.0, 3.0];
        let m = skew3(v);
        let v2 = vee3(m);
        for i in 0..3 {
            assert!((v[i] - v2[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_se3_interp_fraction_zero() {
        let se3_a = Se3::from_rt(So3::exp([0.1, -0.2, 0.3]), [1.0, 2.0, 3.0]);
        let se3_b = Se3::from_rt(So3::exp([0.4, 0.1, -0.5]), [4.0, -1.0, 2.0]);
        let interp = se3_a.interp(&se3_b, 0.0);
        let t = interp.translation;
        let t_a = se3_a.translation;
        for i in 0..3 {
            assert!((t[i] - t_a[i]).abs() < 1e-8);
        }
    }
    /// Euclidean 2D metric function.
    fn flat_2d_metric(_p: &[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM] {
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 1.0;
        g[1][1] = 1.0;
        g
    }
    /// Polar coordinates metric: ds^2 = dr^2 + r^2 d\theta^2.
    fn polar_metric(p: &[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM] {
        let r = p[0];
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 1.0;
        g[1][1] = r * r;
        g
    }
    /// Sphere metric (2D): ds^2 = R^2 (dtheta^2 + sin^2(theta) dphi^2).
    fn sphere_2d_metric(p: &[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM] {
        let theta = p[0];
        let r = 1.0;
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = r * r;
        g[1][1] = r * r * theta.sin() * theta.sin();
        g
    }
    /// Euclidean 3D metric function.
    fn flat_3d_metric(_p: &[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM] {
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 1.0;
        g[1][1] = 1.0;
        g[2][2] = 1.0;
        g
    }
    #[test]
    fn test_metric_tensor_nd_euclidean() {
        let m = MetricTensorND::euclidean(3);
        assert!((m.determinant() - 1.0).abs() < 1e-12);
        let inv = m.inverse();
        for (i, inv_row) in inv.iter().enumerate().take(3) {
            for (j, &val) in inv_row.iter().enumerate().take(3) {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_metric_tensor_nd_minkowski() {
        let m = MetricTensorND::minkowski();
        assert!((m.determinant() - (-1.0)).abs() < 1e-12);
        let inv = m.inverse();
        assert!((inv[0][0] - (-1.0)).abs() < 1e-12);
        assert!((inv[1][1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_metric_tensor_nd_inverse_roundtrip() {
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 2.0;
        g[0][1] = 0.5;
        g[1][0] = 0.5;
        g[1][1] = 3.0;
        let m = MetricTensorND::new(2, &g);
        let inv = m.inverse();
        for (i, g_row) in g.iter().enumerate().take(2) {
            for j in 0..2 {
                let val: f64 = g_row
                    .iter()
                    .zip(inv.iter())
                    .map(|(&g_ik, inv_k)| g_ik * inv_k[j])
                    .take(2)
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-10,
                    "g * g^-1 [{i}][{j}] = {val}, expected {expected}"
                );
            }
        }
    }
    #[test]
    fn test_metric_raise_lower_roundtrip() {
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 2.0;
        g[0][1] = 0.5;
        g[1][0] = 0.5;
        g[1][1] = 3.0;
        let m = MetricTensorND::new(2, &g);
        let mut v = [0.0; MAX_DIM];
        v[0] = 1.0;
        v[1] = 2.0;
        let lowered = m.lower_index(&v);
        let raised = m.raise_index(&lowered);
        for i in 0..2 {
            assert!(
                (raised[i] - v[i]).abs() < 1e-10,
                "Raise/lower roundtrip failed at {i}"
            );
        }
    }
    #[test]
    fn test_metric_inner_product_euclidean() {
        let m = MetricTensorND::euclidean(3);
        let mut u = [0.0; MAX_DIM];
        let mut v = [0.0; MAX_DIM];
        u[0] = 1.0;
        u[1] = 2.0;
        u[2] = 3.0;
        v[0] = 4.0;
        v[1] = 5.0;
        v[2] = 6.0;
        let ip = m.inner_product(&u, &v);
        assert!((ip - 32.0).abs() < 1e-12);
    }
    #[test]
    fn test_metric_volume_element() {
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        g[0][0] = 4.0;
        g[1][1] = 9.0;
        let m = MetricTensorND::new(2, &g);
        assert!((m.volume_element() - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_levi_civita_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let conn = LeviCivitaConnection::from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        for sigma in 0..2 {
            for mu in 0..2 {
                for nu in 0..2 {
                    assert!(
                        conn.christoffel[sigma][mu][nu].abs() < 1e-6,
                        "Christoffel[{sigma}][{mu}][{nu}] = {} in flat space",
                        conn.christoffel[sigma][mu][nu]
                    );
                }
            }
        }
    }
    #[test]
    fn test_levi_civita_polar_coords() {
        let r = 2.0;
        let point = [r, 1.0, 0.0, 0.0];
        let conn = LeviCivitaConnection::from_metric_fn(&polar_metric, &point, 2, 1e-5);
        assert!(
            (conn.christoffel[0][1][1] - (-r)).abs() < 0.01,
            "Gamma^r_{{theta theta}} = {}, expected {}",
            conn.christoffel[0][1][1],
            -r
        );
        assert!(
            (conn.christoffel[1][0][1] - 1.0 / r).abs() < 0.01,
            "Gamma^theta_{{r theta}} = {}, expected {}",
            conn.christoffel[1][0][1],
            1.0 / r
        );
    }
    #[test]
    fn test_riemann_tensor_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let riemann = RiemannTensor::from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        for rho in 0..2 {
            for sigma in 0..2 {
                for mu in 0..2 {
                    for nu in 0..2 {
                        assert!(
                            riemann.components[rho][sigma][mu][nu].abs() < 1e-3,
                            "R[{rho}][{sigma}][{mu}][{nu}] = {} in flat space",
                            riemann.components[rho][sigma][mu][nu]
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_riemann_tensor_sphere() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let riemann = RiemannTensor::from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        let mut max_component = 0.0_f64;
        for rho in 0..2 {
            for sigma in 0..2 {
                for mu in 0..2 {
                    for nu in 0..2 {
                        max_component =
                            max_component.max(riemann.components[rho][sigma][mu][nu].abs());
                    }
                }
            }
        }
        assert!(
            max_component > 0.1,
            "Sphere should have non-zero Riemann tensor"
        );
    }
    #[test]
    fn test_riemann_bianchi_identity() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let riemann = RiemannTensor::from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        let viol = riemann.bianchi_identity_violation();
        assert!(
            viol < 0.5,
            "Bianchi identity violation = {viol}, should be small"
        );
    }
    #[test]
    fn test_ricci_tensor_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let ricci = RicciTensor::from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        for mu in 0..2 {
            for nu in 0..2 {
                assert!(
                    ricci.components[mu][nu].abs() < 1e-3,
                    "Ricci[{mu}][{nu}] = {} in flat space",
                    ricci.components[mu][nu]
                );
            }
        }
    }
    #[test]
    fn test_ricci_tensor_symmetry() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let ricci = RicciTensor::from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        let viol = ricci.symmetry_violation();
        assert!(
            viol < 0.1,
            "Ricci tensor should be symmetric, violation = {viol}"
        );
    }
    #[test]
    fn test_ricci_scalar_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let r = ricci_scalar_from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        assert!(r.abs() < 0.01, "Ricci scalar in flat 2D = {r}, expected 0");
    }
    #[test]
    fn test_ricci_scalar_sphere() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let r = ricci_scalar_from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        assert!(
            (r - 2.0).abs() < 0.5,
            "Ricci scalar on unit sphere = {r}, expected ~2"
        );
    }
    #[test]
    fn test_einstein_tensor_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let ricci = RicciTensor::from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        let g_here = flat_2d_metric(&point);
        let metric = MetricTensorND { dim: 2, g: g_here };
        let g_tensor = einstein_tensor(&metric, &ricci);
        for (mu, row) in g_tensor.iter().enumerate() {
            for (nu, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < 0.01,
                    "Einstein[{mu}][{nu}] = {} in flat space",
                    val
                );
            }
        }
    }
    #[test]
    fn test_covariant_derivative_flat_space() {
        let point = [0.5, 0.5, 0.0, 0.0];
        let conn = LeviCivitaConnection::from_metric_fn(&flat_2d_metric, &point, 2, 1e-4);
        let mut v = [0.0; MAX_DIM];
        v[0] = 1.0;
        v[1] = 2.0;
        let mut dv = [[0.0; MAX_DIM]; MAX_DIM];
        dv[0][0] = 0.5;
        dv[1][1] = -0.3;
        let nabla = covariant_derivative_vector(&conn, &v, &dv);
        assert!((nabla[0][0] - 0.5).abs() < 0.01);
        assert!((nabla[1][1] - (-0.3)).abs() < 0.01);
    }
    #[test]
    fn test_geodesic_equation_flat_space() {
        let mut x0 = [0.0; MAX_DIM];
        let mut v0 = [0.0; MAX_DIM];
        x0[0] = 0.0;
        x0[1] = 0.0;
        v0[0] = 1.0;
        v0[1] = 0.5;
        let gcfg = GeodesicConfig {
            dim: 2,
            dt: 0.1,
            steps: 10,
            h_christoffel: 1e-4,
        };
        let traj = geodesic_equation_solve(&flat_2d_metric, &x0, &v0, &gcfg);
        assert_eq!(traj.len(), 11);
        let (final_x, _) = &traj[10];
        assert!(
            (final_x[0] - 1.0).abs() < 0.05,
            "x = {}, expected ~1.0",
            final_x[0]
        );
        assert!(
            (final_x[1] - 0.5).abs() < 0.05,
            "y = {}, expected ~0.5",
            final_x[1]
        );
    }
    #[test]
    fn test_parallel_transport_nd_flat() {
        let curve: Vec<[f64; MAX_DIM]> = (0..10)
            .map(|i| {
                let mut p = [0.0; MAX_DIM];
                p[0] = i as f64 * 0.1;
                p[1] = i as f64 * 0.05;
                p
            })
            .collect();
        let mut v0 = [0.0; MAX_DIM];
        v0[0] = 1.0;
        v0[1] = 0.0;
        let transported = parallel_transport_nd(&flat_2d_metric, 2, &curve, &v0, 1e-4);
        let last = transported.last().unwrap();
        assert!((last[0] - 1.0).abs() < 0.01);
        assert!(last[1].abs() < 0.01);
    }
    #[test]
    fn test_killing_vector_rotation() {
        let xi_fn = |p: &[f64; MAX_DIM]| -> [f64; MAX_DIM] {
            let mut xi = [0.0; MAX_DIM];
            xi[0] = -p[1];
            xi[1] = p[0];
            xi
        };
        let point = [1.0, 1.0, 0.0, 0.0];
        let viol = killing_equation_violation(&flat_2d_metric, &xi_fn, &point, 2, 1e-4);
        assert!(
            viol < 0.01,
            "Rotation Killing vector violation = {viol}, should be ~0"
        );
    }
    #[test]
    fn test_killing_vector_translation() {
        let xi_fn = |_p: &[f64; MAX_DIM]| -> [f64; MAX_DIM] {
            let mut xi = [0.0; MAX_DIM];
            xi[0] = 1.0;
            xi
        };
        let point = [2.0, 3.0, 0.0, 0.0];
        let viol = killing_equation_violation(&flat_2d_metric, &xi_fn, &point, 2, 1e-4);
        assert!(
            viol < 0.01,
            "Translation Killing vector violation = {viol}, should be ~0"
        );
    }
    #[test]
    fn test_non_killing_vector() {
        let xi_fn = |p: &[f64; MAX_DIM]| -> [f64; MAX_DIM] {
            let mut xi = [0.0; MAX_DIM];
            xi[0] = p[0];
            xi[1] = p[1];
            xi
        };
        let point = [1.0, 1.0, 0.0, 0.0];
        let viol = killing_equation_violation(&flat_2d_metric, &xi_fn, &point, 2, 1e-4);
        assert!(
            viol > 0.5,
            "Dilation should NOT be a Killing vector, violation = {viol}"
        );
    }
    #[test]
    fn test_geodesic_deviation_flat_space() {
        let x0 = [0.0; MAX_DIM];
        let mut v0 = [0.0; MAX_DIM];
        v0[0] = 1.0;
        let gcfg = GeodesicConfig {
            dim: 2,
            dt: 0.1,
            steps: 10,
            h_christoffel: 1e-4,
        };
        let traj = geodesic_equation_solve(&flat_2d_metric, &x0, &v0, &gcfg);
        let mut xi0 = [0.0; MAX_DIM];
        xi0[1] = 1.0;
        let dxi0 = [0.0; MAX_DIM];
        let dcfg = DeviationConfig {
            dim: 2,
            dt: 0.1,
            h: 1e-4,
        };
        let deviation = geodesic_deviation(&flat_2d_metric, &traj, &xi0, &dxi0, &dcfg);
        let last = deviation.last().unwrap();
        assert!(
            (last[1] - 1.0).abs() < 0.05,
            "Deviation should remain ~1.0, got {}",
            last[1]
        );
    }
    #[test]
    fn test_gauss_bonnet_2d_flat() {
        let point = [1.0, 1.0, 0.0, 0.0];
        let integrand = gauss_bonnet_2d_integrand(&flat_2d_metric, &point, 1e-4);
        assert!(
            integrand.abs() < 0.01,
            "Flat space GB integrand = {integrand}"
        );
    }
    #[test]
    fn test_levi_civita_connection_symmetry() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let conn = LeviCivitaConnection::from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        for sigma in 0..2 {
            for mu in 0..2 {
                for nu in 0..2 {
                    let diff =
                        (conn.christoffel[sigma][mu][nu] - conn.christoffel[sigma][nu][mu]).abs();
                    assert!(
                        diff < 0.01,
                        "Christoffel not symmetric: [{sigma}][{mu}][{nu}] diff = {diff}"
                    );
                }
            }
        }
    }
    #[test]
    fn test_weyl_tensor_3d_flat() {
        let point = [1.0, 1.0, 1.0, 0.0];
        let weyl = weyl_tensor(&flat_3d_metric, &point, 3, 1e-4);
        for (rho, w_rho) in weyl.iter().enumerate() {
            for (sigma, w_sig) in w_rho.iter().enumerate() {
                for (mu, w_mu) in w_sig.iter().enumerate() {
                    for (nu, &val) in w_mu.iter().enumerate() {
                        assert!(
                            val.abs() < 0.1,
                            "Weyl[{rho}][{sigma}][{mu}][{nu}] = {} in flat 3D",
                            val
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_riemann_antisymmetry() {
        let theta = 1.0;
        let point = [theta, 1.0, 0.0, 0.0];
        let riemann = RiemannTensor::from_metric_fn(&sphere_2d_metric, &point, 2, 1e-4);
        for rho in 0..2 {
            for sigma in 0..2 {
                for mu in 0..2 {
                    for nu in 0..2 {
                        let sum = riemann.components[rho][sigma][mu][nu]
                            + riemann.components[rho][sigma][nu][mu];
                        assert!(
                            sum.abs() < 0.1,
                            "Riemann not antisymmetric in last pair: [{rho}][{sigma}][{mu}][{nu}], sum = {sum}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_metric_determinant_4d() {
        let m = MetricTensorND::minkowski();
        assert!((m.determinant() - (-1.0)).abs() < 1e-12);
    }
    // ── mat3_exp (scaling exponent, accuracy, non-finite input) ──────────────
    /// Independent reference: `exp(A)` by a long Taylor series (modest ‖A‖ only).
    fn mat3_exp_taylor_reference(a: Mat3) -> Mat3 {
        let mut result = mat3_identity();
        let mut term = mat3_identity();
        for k in 1..=80u32 {
            term = mat3_scale(mat3_mul(term, a), 1.0 / f64::from(k));
            result = mat3_add(result, term);
        }
        result
    }
    fn mat3_max_abs_diff(a: Mat3, b: Mat3) -> f64 {
        a.iter()
            .flatten()
            .zip(b.iter().flatten())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max)
    }
    #[test]
    fn test_mat3_exp_large_diagonal_does_not_overflow() {
        // The old scaling exponent grew linearly with the norm: for ‖A‖ >= 44.4
        // `2_u64.pow(s)` overflowed (panic in debug, NaN in release).
        let e = mat3_exp([[50.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let e50 = 50.0_f64.exp();
        assert!(
            ((e[0][0] - e50) / e50).abs() < 1e-10,
            "exp(50) entry = {}, expected {e50}",
            e[0][0]
        );
        assert!((e[1][1] - 1.0).abs() < 1e-10, "e[1][1] = {}", e[1][1]);
        assert!((e[2][2] - 1.0).abs() < 1e-10, "e[2][2] = {}", e[2][2]);
        for (i, j) in [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)] {
            assert!(
                e[i][j].abs() < 1e-8,
                "off-diagonal e[{i}][{j}] = {}",
                e[i][j]
            );
        }
    }
    #[test]
    fn test_mat3_exp_large_skew_is_orthogonal_and_matches_rodrigues() {
        // |ω| = 100: a very long rotation, norm of the generator ≈ 141.
        let omega = [100.0 / 3.0_f64.sqrt(); 3];
        let r = mat3_exp(skew3(omega));
        let rrt = mat3_mul(r, mat3_transpose(r));
        assert!(
            mat3_max_abs_diff(rrt, mat3_identity()) < 1e-9,
            "R Rᵀ deviates from I by {}",
            mat3_max_abs_diff(rrt, mat3_identity())
        );
        assert!((mat3_det(r) - 1.0).abs() < 1e-9, "det = {}", mat3_det(r));
        let rodrigues = So3::exp(omega).mat;
        assert!(
            mat3_max_abs_diff(r, rodrigues) < 1e-9,
            "exp(skew) differs from Rodrigues by {}",
            mat3_max_abs_diff(r, rodrigues)
        );
    }
    #[test]
    fn test_mat3_exp_general_matrix_matches_taylor_reference() {
        let a = [[0.3, -0.7, 0.2], [1.1, 0.4, -0.5], [0.0, 0.9, -0.8]];
        let e = mat3_exp(a);
        let reference = mat3_exp_taylor_reference(a);
        assert!(
            mat3_max_abs_diff(e, reference) < 1e-12,
            "diff to Taylor reference = {}",
            mat3_max_abs_diff(e, reference)
        );
    }
    #[test]
    fn test_mat3_exp_inverse_property() {
        let a = [[1.5, -2.0, 0.7], [0.3, 2.2, -1.4], [-0.9, 0.5, 1.1]];
        let e = mat3_exp(a);
        let e_inv = mat3_exp(mat3_scale(a, -1.0));
        let prod = mat3_mul(e, e_inv);
        assert!(
            mat3_max_abs_diff(prod, mat3_identity()) < 1e-9,
            "exp(A) exp(-A) deviates from I by {}",
            mat3_max_abs_diff(prod, mat3_identity())
        );
    }
    #[test]
    fn test_mat3_exp_tiny_matrix_keeps_first_order_term() {
        // A² = 0, so exp(A) = I + A exactly; the old `norm < 1e-12` shortcut
        // returned the bare identity and dropped the linear term.
        let a = [[0.0, 1e-14, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let e = mat3_exp(a);
        assert!(
            ((e[0][1] - 1e-14) / 1e-14).abs() < 1e-12,
            "linear term lost: e[0][1] = {}",
            e[0][1]
        );
        assert!((e[0][0] - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_mat3_exp_non_finite_input_gives_nan() {
        let infinite = mat3_exp([[f64::INFINITY, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        assert!(infinite.iter().flatten().all(|v| v.is_nan()));
        let not_a_number = mat3_exp([[0.0, f64::NAN, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        assert!(not_a_number.iter().flatten().all(|v| v.is_nan()));
    }
    // ── christoffel_symbols_numerical (first-kind index placement) ───────────
    #[test]
    fn test_christoffel_numerical_polar_plane() {
        // Polar coordinates r(u, v) = (v cos u, v sin u, 0): u is the angle and
        // v the radius, so Γ^u_uv = 1/v and Γ^v_uu = -v; all else vanishes.
        let r_fn = |u: f64, v: f64| [v * u.cos(), v * u.sin(), 0.0];
        let (u, v) = (0.3, 2.0);
        let gamma = christoffel_symbols_numerical(&r_fn, u, v, 1e-4);
        assert!(
            (gamma[0][0][1] - 1.0 / v).abs() < 1e-5,
            "Γ^u_uv = {}, expected {}",
            gamma[0][0][1],
            1.0 / v
        );
        assert!(
            (gamma[0][1][0] - 1.0 / v).abs() < 1e-5,
            "Γ^u_vu = {}, expected {}",
            gamma[0][1][0],
            1.0 / v
        );
        assert!(
            (gamma[1][0][0] + v).abs() < 1e-5,
            "Γ^v_uu = {}, expected {}",
            gamma[1][0][0],
            -v
        );
        for (k, i, j) in [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)] {
            assert!(
                gamma[k][i][j].abs() < 1e-5,
                "Γ^{k}_{i}{j} = {} should vanish for the polar plane",
                gamma[k][i][j]
            );
        }
    }
    #[test]
    fn test_christoffel_numerical_unit_sphere() {
        // r(u, v) = (sin u cos v, sin u sin v, cos u): u is the polar angle.
        // Γ^u_vv = -sin u cos u and Γ^v_uv = cot u.
        let r_fn = |u: f64, v: f64| [u.sin() * v.cos(), u.sin() * v.sin(), u.cos()];
        let (u, v) = (0.9, 0.3);
        let gamma = christoffel_symbols_numerical(&r_fn, u, v, 1e-4);
        let expected_uvv = -u.sin() * u.cos();
        let expected_vuv = u.cos() / u.sin();
        assert!(
            (gamma[0][1][1] - expected_uvv).abs() < 1e-5,
            "Γ^u_vv = {}, expected {expected_uvv}",
            gamma[0][1][1]
        );
        assert!(
            (gamma[1][0][1] - expected_vuv).abs() < 1e-5,
            "Γ^v_uv = {}, expected {expected_vuv}",
            gamma[1][0][1]
        );
        assert!(gamma[0][0][0].abs() < 1e-5);
        assert!(gamma[0][0][1].abs() < 1e-5);
        assert!(gamma[1][0][0].abs() < 1e-5);
        assert!(gamma[1][1][1].abs() < 1e-5);
    }
    #[test]
    fn test_christoffel_numerical_skewed_parametrisation() {
        // z = u v has F = u v != 0, r_uu = r_vv = 0 and r_uv = (0, 0, 1), so the
        // only non-zero first-kind symbols are Γ_{12,1} = v and Γ_{12,2} = u.
        // With D = 1 + u² + v²:  Γ^u_uv = v / D, Γ^v_uv = u / D, the rest 0.
        let r_fn = |u: f64, v: f64| [u, v, u * v];
        let (u, v) = (0.7, -0.4);
        let d = 1.0 + u * u + v * v;
        let gamma = christoffel_symbols_numerical(&r_fn, u, v, 1e-3);
        assert!(
            (gamma[0][0][1] - v / d).abs() < 1e-6,
            "Γ^u_uv = {}, expected {}",
            gamma[0][0][1],
            v / d
        );
        assert!(
            (gamma[0][1][0] - v / d).abs() < 1e-6,
            "Γ^u_vu = {}, expected {}",
            gamma[0][1][0],
            v / d
        );
        assert!(
            (gamma[1][0][1] - u / d).abs() < 1e-6,
            "Γ^v_uv = {}, expected {}",
            gamma[1][0][1],
            u / d
        );
        for (k, i, j) in [(0, 0, 0), (0, 1, 1), (1, 0, 0), (1, 1, 1)] {
            assert!(
                gamma[k][i][j].abs() < 1e-6,
                "Γ^{k}_{i}{j} = {} should vanish",
                gamma[k][i][j]
            );
        }
    }
    #[test]
    fn test_christoffel_numerical_matches_metric_connection() {
        // Independent cross-check: the connection built from the metric
        // components g = [[E, F], [F, G]] of z = u v must agree with the
        // surface-based symbols.
        let r_fn = |u: f64, v: f64| [u, v, u * v];
        let metric_fn = |p: &[f64; MAX_DIM]| {
            let mut g = [[0.0; MAX_DIM]; MAX_DIM];
            g[0][0] = 1.0 + p[1] * p[1];
            g[0][1] = p[0] * p[1];
            g[1][0] = p[0] * p[1];
            g[1][1] = 1.0 + p[0] * p[0];
            g
        };
        let (u, v) = (-0.35, 0.8);
        let surface = christoffel_symbols_numerical(&r_fn, u, v, 1e-3);
        let conn = LeviCivitaConnection::from_metric_fn(&metric_fn, &[u, v, 0.0, 0.0], 2, 1e-5);
        for (k, (surface_k, conn_k)) in surface.iter().zip(conn.christoffel.iter()).enumerate() {
            for (i, (surface_i, conn_i)) in surface_k.iter().zip(conn_k.iter()).enumerate() {
                for (j, (from_surface, from_metric)) in
                    surface_i.iter().zip(conn_i.iter()).enumerate()
                {
                    assert!(
                        (from_surface - from_metric).abs() < 1e-5,
                        "Γ^{k}_{i}{j}: surface {from_surface} vs metric {from_metric}"
                    );
                }
            }
        }
    }
    #[test]
    fn test_geodesic_curve_polar_plane_is_a_straight_line() {
        // In polar coordinates r(u, v) = (v cos u, v sin u, 0) geodesics are the
        // straight lines of the plane.  Start at (1, 0) heading tangentially
        // (u' = 1, v' = 0): the exact geodesic is x = 1, y = s, i.e.
        // u = atan(s), v = sqrt(1 + s²).  The swapped first-kind symbols bent
        // this path into a curve that fell towards the origin (end point
        // ≈ (0.46, 0.44) instead of (1, 1) after s = 1).
        let r_fn = |u: f64, v: f64| [v * u.cos(), v * u.sin(), 0.0];
        let ds = 0.01;
        let steps = 100;
        let curve = geodesic_curve(&r_fn, 0.0, 1.0, 1.0, 0.0, steps, ds);
        assert_eq!(curve.len(), steps + 1);
        for (k, &(u, v)) in curve.iter().enumerate() {
            let x = v * u.cos();
            assert!(
                (x - 1.0).abs() < 2e-2,
                "step {k}: x = {x} should stay on the line x = 1"
            );
        }
        let (u_end, v_end) = curve[steps];
        assert!(
            (v_end * u_end.sin() - 1.0).abs() < 2e-2,
            "arc length s = 1 along the line gives y = 1, got {}",
            v_end * u_end.sin()
        );
        assert!(
            (v_end - 2.0_f64.sqrt()).abs() < 2e-2,
            "radius should grow to sqrt(2), got {v_end}"
        );
    }
}
