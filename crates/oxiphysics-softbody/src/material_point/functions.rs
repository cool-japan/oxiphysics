//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions_2::g2p_transfer_with_dt;
use super::types::{MaterialType, MpmGrid, MpmParticle, SandParams, SnowParams};

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Length of a 3-vector.
#[inline]
pub fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3-vector (returns zero vector if length < 1e-15).
#[inline]
pub fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}
/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Identity 3×3 matrix.
#[inline]
pub fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Zero 3×3 matrix.
#[inline]
pub fn mat3_zero() -> [[f64; 3]; 3] {
    [[0.0; 3]; 3]
}
/// Multiply two 3×3 matrices.
#[inline]
pub fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}
/// Multiply a 3×3 matrix by a 3-vector.
#[inline]
pub fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// Transpose a 3×3 matrix.
#[inline]
pub fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Trace of a 3×3 matrix.
#[inline]
pub fn mat3_trace(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}
/// Determinant of a 3×3 matrix.
#[inline]
pub fn mat3_det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Add two 3×3 matrices.
#[inline]
pub fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][j] + b[i][j];
        }
    }
    r
}
/// Scale a 3×3 matrix.
#[inline]
pub fn mat3_scale(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = m[i][j] * s;
        }
    }
    r
}
/// Outer product of two 3-vectors: a ⊗ b.
#[inline]
pub fn outer3(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i] * b[j];
        }
    }
    r
}
/// Inverse of a 3×3 matrix. Returns identity if singular.
#[inline]
pub fn mat3_inv(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = mat3_det(m);
    if det.abs() < 1e-15 {
        return mat3_identity();
    }
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}
/// Polar decomposition: F = R * S.
/// Returns (R, S) where R is rotation, S is symmetric positive definite.
/// Uses iterative method (Higham's algorithm).
pub fn polar_decomp(f: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut r = f;
    for _ in 0..20 {
        let r_t = mat3_transpose(r);
        let r_t_inv = mat3_inv(r_t);
        let r_next_raw = mat3_add(r, r_t_inv);
        r = mat3_scale(r_next_raw, 0.5);
        let diff = mat3_add(r, mat3_scale(mat3_transpose(r), -1.0));
        let err: f64 = diff.iter().flat_map(|row| row.iter()).map(|x| x * x).sum();
        if err < 1e-20 {
            break;
        }
    }
    let s = mat3_mul(mat3_transpose(r), f);
    (r, s)
}
/// SVD of a 3×3 matrix (approximate via polar decomposition).
/// Returns (U, sigma_vec, V^T).
pub fn mat3_svd(m: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3]) {
    let (u, s) = polar_decomp(m);
    let sigma = [s[0][0], s[1][1], s[2][2]];
    let vt = mat3_transpose(mat3_mul(mat3_transpose(u), m));
    (u, sigma, vt)
}
/// Quadratic B-spline kernel weight N(x).
/// Defined for |x| in \[0, 1.5\].
#[inline]
pub fn bspline_quadratic(x: f64) -> f64 {
    let x = x.abs();
    if x < 0.5 {
        0.75 - x * x
    } else if x < 1.5 {
        let t = 1.5 - x;
        0.5 * t * t
    } else {
        0.0
    }
}
/// Derivative of quadratic B-spline kernel.
#[inline]
pub fn bspline_quadratic_grad(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 0.5 {
        -2.0 * x
    } else if ax < 1.5 {
        let sign = if x >= 0.0 { 1.0 } else { -1.0 };
        sign * (ax - 1.5)
    } else {
        0.0
    }
}
/// Cubic B-spline kernel weight.
#[inline]
pub fn bspline_cubic(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 1.0 {
        0.5 * ax * ax * ax - ax * ax + 2.0 / 3.0
    } else if ax < 2.0 {
        let t = 2.0 - ax;
        t * t * t / 6.0
    } else {
        0.0
    }
}
/// Derivative of cubic B-spline kernel.
#[inline]
pub fn bspline_cubic_grad(x: f64) -> f64 {
    let ax = x.abs();
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    if ax < 1.0 {
        sign * (1.5 * ax * ax - 2.0 * ax)
    } else if ax < 2.0 {
        let t = 2.0 - ax;
        -sign * t * t / 2.0
    } else {
        0.0
    }
}
/// Neo-Hookean Piola-Kirchhoff stress P.
///
/// P = μ (F - F^{-T}) + λ ln(J) F^{-T}
pub fn neo_hookean_pk1(f: [[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let f_inv_t = mat3_transpose(mat3_inv(f));
    let log_j = if j > 1e-10 { j.ln() } else { -10.0 };
    let term1 = mat3_scale(mat3_add(f, mat3_scale(f_inv_t, -1.0)), mu);
    let term2 = mat3_scale(f_inv_t, lambda * log_j);
    mat3_add(term1, term2)
}
/// Neo-Hookean Cauchy stress σ.
///
/// σ = (1/J) * P * F^T
pub fn neo_hookean_cauchy(f: [[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 3]; 3] {
    let j = mat3_det(f).abs();
    if j < 1e-15 {
        return mat3_zero();
    }
    let pk1 = neo_hookean_pk1(f, mu, lambda);
    let ft = mat3_transpose(f);
    let pk1_ft = mat3_mul(pk1, ft);
    mat3_scale(pk1_ft, 1.0 / j)
}
/// Fixed corotated constitutive model (Stomakhin et al. 2012).
///
/// P = 2μ (F - R) + λ (J - 1) J F^{-T}
pub fn fixed_corotated_pk1(f: [[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 3]; 3] {
    let (r, _s) = polar_decomp(f);
    let j = mat3_det(f);
    let f_inv_t = mat3_transpose(mat3_inv(f));
    let fe_minus_r = mat3_add(f, mat3_scale(r, -1.0));
    let term1 = mat3_scale(fe_minus_r, 2.0 * mu);
    let term2 = mat3_scale(f_inv_t, lambda * (j - 1.0) * j);
    mat3_add(term1, term2)
}
/// Apply von Mises yield criterion to the elastic deformation gradient.
///
/// Projects Fe onto the yield surface if the von Mises stress exceeds yield_stress.
/// Returns the corrected Fe.
pub fn von_mises_return_mapping(
    fe: [[f64; 3]; 3],
    mu: f64,
    lambda: f64,
    yield_stress: f64,
) -> [[f64; 3]; 3] {
    let cauchy = neo_hookean_cauchy(fe, mu, lambda);
    let trace_c = mat3_trace(cauchy);
    let mut deviatoric = cauchy;
    let mean_stress = trace_c / 3.0;
    for (i, row) in deviatoric.iter_mut().enumerate() {
        row[i] -= mean_stress;
    }
    let norm_sq: f64 = deviatoric
        .iter()
        .flat_map(|row| row.iter())
        .map(|x| x * x)
        .sum();
    let norm_dev = norm_sq.sqrt();
    let vm_stress = (1.5_f64).sqrt() * norm_dev;
    if vm_stress <= yield_stress {
        return fe;
    }
    let scale = yield_stress / vm_stress;
    let (r, _s) = polar_decomp(fe);
    let u_matrix = mat3_mul(mat3_transpose(r), fe);
    let trace_u = mat3_trace(u_matrix);
    let mut u_dev = u_matrix;
    let mean_u = trace_u / 3.0;
    for (i, row) in u_dev.iter_mut().enumerate() {
        row[i] -= mean_u;
    }
    let u_corrected = {
        let mut uc = mat3_zero();
        for (i, uc_row) in uc.iter_mut().enumerate() {
            for (j, uc_ij) in uc_row.iter_mut().enumerate() {
                *uc_ij = mean_u * (if i == j { 1.0 } else { 0.0 }) + u_dev[i][j] * scale;
            }
        }
        uc
    };
    mat3_mul(r, u_corrected)
}
/// Snow hardening model (Stomakhin et al. 2013).
///
/// Returns (mu_e, lambda_e) scaled by hardening exponent.
pub fn snow_hardening_params(
    fe: [[f64; 3]; 3],
    fp: [[f64; 3]; 3],
    mu0: f64,
    lambda0: f64,
    hardening: f64,
    theta_c: f64,
    theta_s: f64,
) -> (f64, f64) {
    let je = mat3_det(fe).abs();
    let jp = mat3_det(fp).abs();
    let je_clamped = je.max(1e-6).min(1.0 - theta_c + theta_s);
    let exp_h = (hardening * (1.0 - jp)).exp();
    let mu_e = mu0 * exp_h;
    let lambda_e = lambda0 * exp_h;
    let _ = je_clamped;
    (mu_e, lambda_e)
}
/// Snow plastic update (Stomakhin et al. 2013).
///
/// Clamps singular values of Fe to \[1 - theta_c, 1 + theta_s\].
pub fn snow_plastic_update(
    fe: [[f64; 3]; 3],
    _fp: [[f64; 3]; 3],
    f_total: [[f64; 3]; 3],
    theta_c: f64,
    theta_s: f64,
) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let (u, sigma, vt) = mat3_svd(fe);
    let sigma_clamped = [
        sigma[0].max(1.0 - theta_c).min(1.0 + theta_s),
        sigma[1].max(1.0 - theta_c).min(1.0 + theta_s),
        sigma[2].max(1.0 - theta_c).min(1.0 + theta_s),
    ];
    let diag_s = [
        [sigma_clamped[0], 0.0, 0.0],
        [0.0, sigma_clamped[1], 0.0],
        [0.0, 0.0, sigma_clamped[2]],
    ];
    let fe_new = mat3_mul(mat3_mul(u, diag_s), vt);
    let fp_new = mat3_mul(mat3_inv(fe_new), f_total);
    (fe_new, fp_new)
}
/// Drucker-Prager yield criterion for sand.
///
/// Applies the DP yield surface: f = ||dev(ε)|| - A * tr(ε)
/// where A is related to the friction angle.
pub fn drucker_prager_return_mapping(
    fe: [[f64; 3]; 3],
    mu: f64,
    lambda: f64,
    friction_angle_deg: f64,
) -> [[f64; 3]; 3] {
    let friction = (friction_angle_deg * std::f64::consts::PI / 180.0).tan();
    let alpha = friction / (3.0_f64).sqrt();
    let (u, sigma, vt) = mat3_svd(fe);
    let epsilon = [sigma[0].ln(), sigma[1].ln(), sigma[2].ln()];
    let tr_eps = epsilon[0] + epsilon[1] + epsilon[2];
    let eps_dev = [
        epsilon[0] - tr_eps / 3.0,
        epsilon[1] - tr_eps / 3.0,
        epsilon[2] - tr_eps / 3.0,
    ];
    let norm_dev = len3(eps_dev);
    if tr_eps > 0.0 {
        return mat3_identity();
    }
    let dp_f = norm_dev + alpha * tr_eps;
    if dp_f <= 0.0 {
        return fe;
    }
    let delta_gamma = dp_f;
    let eps_corrected = [
        epsilon[0] - delta_gamma * (eps_dev[0] / (norm_dev + 1e-15) + alpha / 3.0),
        epsilon[1] - delta_gamma * (eps_dev[1] / (norm_dev + 1e-15) + alpha / 3.0),
        epsilon[2] - delta_gamma * (eps_dev[2] / (norm_dev + 1e-15) + alpha / 3.0),
    ];
    let _ = (mu, lambda);
    let sigma_new = [
        eps_corrected[0].exp(),
        eps_corrected[1].exp(),
        eps_corrected[2].exp(),
    ];
    let diag_s = [
        [sigma_new[0], 0.0, 0.0],
        [0.0, sigma_new[1], 0.0],
        [0.0, 0.0, sigma_new[2]],
    ];
    mat3_mul(mat3_mul(u, diag_s), vt)
}
/// Compute B-spline quadratic weights and gradients for a particle.
///
/// Returns (weights\[3\]\[3\], grad_weights\[3\]\[3\]) for each dimension × 3 stencil.
pub fn compute_weights(base: [i64; 3], frac: [f64; 3], h: f64) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut w = [[0.0f64; 3]; 3];
    let mut dw = [[0.0f64; 3]; 3];
    for d in 0..3 {
        for k in 0..3 {
            let x = frac[d] - (k as f64 - 1.0);
            w[d][k] = bspline_quadratic(x);
            dw[d][k] = bspline_quadratic_grad(x) / h;
        }
        let _ = base[d];
    }
    (w, dw)
}
/// APIC P2G transfer: scatter particle momentum and mass to grid.
pub fn p2g_transfer(grid: &mut MpmGrid, particles: &[MpmParticle], dt: f64, gravity: [f64; 3]) {
    let h = grid.h;
    let inv_h = 1.0 / h;
    let inv_h2 = inv_h * inv_h;
    for p in particles {
        if !p.active {
            continue;
        }
        let (base, frac) = grid.base_and_frac(p.position);
        let (w, _dw) = compute_weights(base, frac, h);
        let stress = match p.material {
            MaterialType::NeoHookean => neo_hookean_pk1(p.deformation_gradient, p.mu0, p.lambda0),
            MaterialType::VonMises => neo_hookean_pk1(p.fe, p.mu0, p.lambda0),
            MaterialType::Snow => {
                let (mu_e, lambda_e) =
                    snow_hardening_params(p.fe, p.fp, p.mu0, p.lambda0, p.hardening, 0.025, 0.0075);
                fixed_corotated_pk1(p.fe, mu_e, lambda_e)
            }
            MaterialType::Sand => neo_hookean_pk1(p.fe, p.mu0, p.lambda0),
            MaterialType::Fluid => {
                let j = mat3_det(p.deformation_gradient);
                let pressure = p.lambda0 * (j - 1.0);
                let mut s = mat3_zero();
                for (i, row) in s.iter_mut().enumerate() {
                    row[i] = -pressure;
                }
                s
            }
            MaterialType::Rigid => mat3_zero(),
        };
        let vol_stress = mat3_scale(stress, p.volume0);
        for (di, w0_di) in w[0].iter().enumerate() {
            for (dj, w1_dj) in w[1].iter().enumerate() {
                for (dk, w2_dk) in w[2].iter().enumerate() {
                    let ix = base[0] + di as i64 - 1;
                    let iy = base[1] + dj as i64 - 1;
                    let iz = base[2] + dk as i64 - 1;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let w_ijk = w0_di * w1_dj * w2_dk;
                    let node_idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    let xip = sub3(
                        grid.node_pos(ix as usize, iy as usize, iz as usize),
                        p.position,
                    );
                    grid.nodes[node_idx].mass += w_ijk * p.mass;
                    let c_xip = mat3_vec(p.affine_c, xip);
                    let mom_contrib = scale3(add3(p.velocity, c_xip), p.mass * w_ijk);
                    grid.nodes[node_idx].momentum =
                        add3(grid.nodes[node_idx].momentum, mom_contrib);
                    let grad_w = scale3(xip, w_ijk * inv_h2 * 4.0);
                    let f_force = mat3_vec(vol_stress, grad_w);
                    grid.nodes[node_idx].force =
                        add3(grid.nodes[node_idx].force, scale3(f_force, -1.0));
                    let g_force = scale3(gravity, p.mass * w_ijk);
                    grid.nodes[node_idx].force = add3(grid.nodes[node_idx].force, g_force);
                }
            }
        }
    }
    for node in &mut grid.nodes {
        node.normalize_velocity();
        node.apply_force(dt);
    }
}
/// APIC G2P transfer: gather velocity and affine C matrix from grid to particles.
pub fn g2p_transfer_dt(grid: &MpmGrid, particles: &mut [MpmParticle], _dt: f64) {
    let h = grid.h;
    let inv_h = 1.0 / h;
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        let (base, frac) = grid.base_and_frac(p.position);
        let (w, _dw) = compute_weights(base, frac, h);
        let mut new_vel = [0.0f64; 3];
        let mut new_c = mat3_zero();
        for (di, w0_di) in w[0].iter().enumerate() {
            for (dj, w1_dj) in w[1].iter().enumerate() {
                for (dk, w2_dk) in w[2].iter().enumerate() {
                    let ix = base[0] + di as i64 - 1;
                    let iy = base[1] + dj as i64 - 1;
                    let iz = base[2] + dk as i64 - 1;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let w_ijk = w0_di * w1_dj * w2_dk;
                    let node_idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    let vi = grid.nodes[node_idx].velocity;
                    new_vel = add3(new_vel, scale3(vi, w_ijk));
                    let xip = sub3(
                        grid.node_pos(ix as usize, iy as usize, iz as usize),
                        p.position,
                    );
                    let outer = outer3(vi, xip);
                    let scaled_outer = mat3_scale(outer, w_ijk * 4.0 * inv_h * inv_h);
                    new_c = mat3_add(new_c, scaled_outer);
                }
            }
        }
        p.velocity = new_vel;
        p.affine_c = new_c;
    }
}
/// Update particle deformation gradient: F ← (I + dt * C) * F.
pub fn update_deformation_gradient(
    particles: &mut [MpmParticle],
    dt: f64,
    theta_c: f64,
    theta_s: f64,
) {
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        let i_plus_dt_c = {
            let dt_c = mat3_scale(p.affine_c, dt);
            mat3_add(mat3_identity(), dt_c)
        };
        let f_new = mat3_mul(i_plus_dt_c, p.deformation_gradient);
        p.deformation_gradient = f_new;
        match p.material {
            MaterialType::Snow => {
                let fe_trial = mat3_mul(i_plus_dt_c, p.fe);
                let (fe_new, fp_new) = snow_plastic_update(fe_trial, p.fp, f_new, theta_c, theta_s);
                p.fe = fe_new;
                p.fp = fp_new;
            }
            MaterialType::VonMises => {
                let fe_trial = mat3_mul(i_plus_dt_c, p.fe);
                p.fe = von_mises_return_mapping(fe_trial, p.mu0, p.lambda0, p.yield_param);
                p.fp = mat3_mul(mat3_inv(p.fe), f_new);
            }
            MaterialType::Sand => {
                let fe_trial = mat3_mul(i_plus_dt_c, p.fe);
                p.fe = drucker_prager_return_mapping(fe_trial, p.mu0, p.lambda0, p.yield_param);
                p.fp = mat3_mul(mat3_inv(p.fe), f_new);
            }
            _ => {
                p.fe = f_new;
            }
        }
    }
}
/// Advect particles: x ← x + dt * v.
pub fn advect_particles(particles: &mut [MpmParticle], dt: f64) {
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        p.position = add3(p.position, scale3(p.velocity, dt));
    }
}
/// Perform one full MPM time step.
pub fn mpm_step(
    grid: &mut MpmGrid,
    particles: &mut [MpmParticle],
    dt: f64,
    gravity: [f64; 3],
    theta_c: f64,
    theta_s: f64,
    boundary_thickness: usize,
) {
    grid.reset();
    p2g_transfer(grid, particles, dt, gravity);
    grid.apply_boundary_conditions(boundary_thickness);
    grid.apply_slip_floor_bc(1);
    g2p_transfer(grid, particles);
    update_deformation_gradient(particles, dt, theta_c, theta_s);
    advect_particles(particles, dt);
}
/// Compute the signed distance field (level-set) from a particle cloud.
///
/// Uses a simple particle splatting approach.
pub fn compute_level_set(
    particles: &[MpmParticle],
    grid: &MpmGrid,
    kernel_radius: f64,
) -> Vec<f64> {
    let n = grid.num_nodes();
    let mut phi = vec![f64::MAX; n];
    let [nx, ny, nz] = grid.dims;
    for p in particles {
        if !p.active {
            continue;
        }
        let (base, _frac) = grid.base_and_frac(p.position);
        let r_cells = (kernel_radius / grid.h) as i64 + 2;
        for di in -r_cells..=r_cells {
            for dj in -r_cells..=r_cells {
                for dk in -r_cells..=r_cells {
                    let ix = base[0] + di;
                    let iy = base[1] + dj;
                    let iz = base[2] + dk;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let xg = grid.node_pos(ix as usize, iy as usize, iz as usize);
                    let d = len3(sub3(xg, p.position));
                    let idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    if d - kernel_radius < phi[idx] {
                        phi[idx] = d - kernel_radius;
                    }
                }
            }
        }
    }
    let _ = (nx, ny, nz);
    phi
}
/// Apply level-set coupling: push particles out of solid boundaries.
pub fn apply_level_set_coupling(
    particles: &mut [MpmParticle],
    phi: &[f64],
    grid: &MpmGrid,
    normal_restitution: f64,
) {
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        let (base, frac) = grid.base_and_frac(p.position);
        let (w, _) = compute_weights(base, frac, grid.h);
        let mut phi_p = 0.0;
        let grad_phi = [0.0f64; 3];
        for (di, w0_di) in w[0].iter().enumerate() {
            for (dj, w1_dj) in w[1].iter().enumerate() {
                for (dk, w2_dk) in w[2].iter().enumerate() {
                    let ix = base[0] + di as i64 - 1;
                    let iy = base[1] + dj as i64 - 1;
                    let iz = base[2] + dk as i64 - 1;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let w_ijk = w0_di * w1_dj * w2_dk;
                    let idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    phi_p += w_ijk * phi[idx];
                }
            }
        }
        if phi_p < 0.0 {
            p.level_set_phi = phi_p;
            let n_hat = norm3(grad_phi);
            let v_n = dot3(p.velocity, n_hat);
            if v_n < 0.0 {
                let v_reflect = scale3(n_hat, -normal_restitution * v_n);
                p.velocity = add3(p.velocity, v_reflect);
            }
            p.position = add3(p.position, scale3(n_hat, -phi_p));
        }
    }
}
/// Simple pressure projection for fluid MPM.
///
/// Solves the incompressibility constraint using Jacobi iteration.
pub fn pressure_projection(grid: &mut MpmGrid, dt: f64, density0: f64, num_iterations: usize) {
    let [nx, ny, nz] = grid.dims;
    let h = grid.h;
    let mut pressure = vec![0.0f64; nx * ny * nz];
    let inv_h = 1.0 / h;
    let inv_dt = 1.0 / dt;
    let inv_rho_h = inv_h / density0;
    for _iter in 0..num_iterations {
        let pressure_old = pressure.clone();
        for ix in 1..nx - 1 {
            for iy in 1..ny - 1 {
                for iz in 1..nz - 1 {
                    let idx = grid.idx(ix, iy, iz);
                    if grid.nodes[idx].mass < 1e-15 {
                        continue;
                    }
                    let vxp = grid.nodes[grid.idx(ix + 1, iy, iz)].velocity[0];
                    let vxm = grid.nodes[grid.idx(ix - 1, iy, iz)].velocity[0];
                    let vyp = grid.nodes[grid.idx(ix, iy + 1, iz)].velocity[1];
                    let vym = grid.nodes[grid.idx(ix, iy - 1, iz)].velocity[1];
                    let vzp = grid.nodes[grid.idx(ix, iy, iz + 1)].velocity[2];
                    let vzm = grid.nodes[grid.idx(ix, iy, iz - 1)].velocity[2];
                    let div_v = ((vxp - vxm) + (vyp - vym) + (vzp - vzm)) * 0.5 * inv_h;
                    let rhs = -density0 * div_v * inv_dt;
                    let pp = pressure_old[grid.idx(ix + 1, iy, iz)]
                        + pressure_old[grid.idx(ix - 1, iy, iz)]
                        + pressure_old[grid.idx(ix, iy + 1, iz)]
                        + pressure_old[grid.idx(ix, iy - 1, iz)]
                        + pressure_old[grid.idx(ix, iy, iz + 1)]
                        + pressure_old[grid.idx(ix, iy, iz - 1)];
                    pressure[idx] = (pp - h * h * rhs) / 6.0;
                }
            }
        }
        let _ = inv_rho_h;
    }
    for ix in 1..nx - 1 {
        for iy in 1..ny - 1 {
            for iz in 1..nz - 1 {
                let idx = grid.idx(ix, iy, iz);
                if grid.nodes[idx].mass < 1e-15 {
                    continue;
                }
                let dp_dx = (pressure[grid.idx(ix + 1, iy, iz)]
                    - pressure[grid.idx(ix - 1, iy, iz)])
                    * 0.5
                    * inv_h;
                let dp_dy = (pressure[grid.idx(ix, iy + 1, iz)]
                    - pressure[grid.idx(ix, iy - 1, iz)])
                    * 0.5
                    * inv_h;
                let dp_dz = (pressure[grid.idx(ix, iy, iz + 1)]
                    - pressure[grid.idx(ix, iy, iz - 1)])
                    * 0.5
                    * inv_h;
                let dp = [dp_dx, dp_dy, dp_dz];
                let inv_rho = if grid.nodes[idx].mass > 1e-15 {
                    1.0 / density0
                } else {
                    0.0
                };
                grid.nodes[idx].velocity = sub3(grid.nodes[idx].velocity, scale3(dp, dt * inv_rho));
            }
        }
    }
}
/// Initialize a uniform block of particles (for simulation setup).
pub fn init_particle_block(
    min_corner: [f64; 3],
    max_corner: [f64; 3],
    particle_spacing: f64,
    mass_per_particle: f64,
    material: MaterialType,
    mu0: f64,
    lambda0: f64,
) -> Vec<MpmParticle> {
    let mut particles = Vec::new();
    let vol0 = particle_spacing * particle_spacing * particle_spacing;
    let mut x = min_corner[0];
    while x <= max_corner[0] {
        let mut y = min_corner[1];
        while y <= max_corner[1] {
            let mut z = min_corner[2];
            while z <= max_corner[2] {
                particles.push(MpmParticle::new(
                    [x, y, z],
                    [0.0; 3],
                    mass_per_particle,
                    vol0,
                    material,
                    mu0,
                    lambda0,
                ));
                z += particle_spacing;
            }
            y += particle_spacing;
        }
        x += particle_spacing;
    }
    particles
}
/// Initialize a sphere of particles.
pub fn init_particle_sphere(
    center: [f64; 3],
    radius: f64,
    particle_spacing: f64,
    mass_per_particle: f64,
    material: MaterialType,
    mu0: f64,
    lambda0: f64,
) -> Vec<MpmParticle> {
    let block = init_particle_block(
        [center[0] - radius, center[1] - radius, center[2] - radius],
        [center[0] + radius, center[1] + radius, center[2] + radius],
        particle_spacing,
        mass_per_particle,
        material,
        mu0,
        lambda0,
    );
    block
        .into_iter()
        .filter(|p| len3(sub3(p.position, center)) <= radius)
        .collect()
}
/// Neo-Hookean elastic energy density ψ(F).
///
/// ψ = μ/2 (tr(F^T F) - 3) - μ ln(J) + λ/2 (ln J)²
pub fn neo_hookean_energy(f: [[f64; 3]; 3], mu: f64, lambda: f64) -> f64 {
    let ft_f = mat3_mul(mat3_transpose(f), f);
    let trace = mat3_trace(ft_f);
    let j = mat3_det(f);
    if j <= 0.0 {
        return f64::MAX;
    }
    let log_j = j.ln();
    0.5 * mu * (trace - 3.0) - mu * log_j + 0.5 * lambda * log_j * log_j
}
/// Fixed corotated elastic energy density.
pub fn fixed_corotated_energy(f: [[f64; 3]; 3], mu: f64, lambda: f64) -> f64 {
    let (r, _) = polar_decomp(f);
    let f_minus_r = mat3_add(f, mat3_scale(r, -1.0));
    let norm_sq: f64 = f_minus_r
        .iter()
        .flat_map(|row| row.iter())
        .map(|x| x * x)
        .sum();
    let j = mat3_det(f);
    mu * norm_sq + 0.5 * lambda * (j - 1.0) * (j - 1.0)
}
/// FLIP blending: mix PIC and FLIP velocity updates.
///
/// v_new = (1 - flip_ratio) * v_pic + flip_ratio * (v_old + delta_v_flip)
pub fn g2p_flip(
    grid: &MpmGrid,
    grid_old_velocity: &[[f64; 3]],
    particles: &mut [MpmParticle],
    flip_ratio: f64,
) {
    let h = grid.h;
    let inv_h = 1.0 / h;
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        let (base, frac) = grid.base_and_frac(p.position);
        let (w, _) = compute_weights(base, frac, h);
        let mut v_pic = [0.0f64; 3];
        let mut delta_v = [0.0f64; 3];
        let mut new_c = mat3_zero();
        for (di, w0_di) in w[0].iter().enumerate() {
            for (dj, w1_dj) in w[1].iter().enumerate() {
                for (dk, w2_dk) in w[2].iter().enumerate() {
                    let ix = base[0] + di as i64 - 1;
                    let iy = base[1] + dj as i64 - 1;
                    let iz = base[2] + dk as i64 - 1;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let w_ijk = w0_di * w1_dj * w2_dk;
                    let node_idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    let vi_new = grid.nodes[node_idx].velocity;
                    let vi_old = grid_old_velocity[node_idx];
                    v_pic = add3(v_pic, scale3(vi_new, w_ijk));
                    delta_v = add3(delta_v, scale3(sub3(vi_new, vi_old), w_ijk));
                    let xip = sub3(
                        grid.node_pos(ix as usize, iy as usize, iz as usize),
                        p.position,
                    );
                    let outer = outer3(vi_new, xip);
                    new_c = mat3_add(new_c, mat3_scale(outer, w_ijk * 4.0 * inv_h * inv_h));
                }
            }
        }
        let v_flip = add3(p.velocity, delta_v);
        p.velocity = add3(scale3(v_pic, 1.0 - flip_ratio), scale3(v_flip, flip_ratio));
        p.affine_c = new_c;
    }
}
/// Create a snow particle with default snow parameters.
pub fn make_snow_particle(
    position: [f64; 3],
    velocity: [f64; 3],
    volume0: f64,
    params: &SnowParams,
) -> MpmParticle {
    let mass = params.rho0 * volume0;
    let mut p = MpmParticle::new(
        position,
        velocity,
        mass,
        volume0,
        MaterialType::Snow,
        params.mu0,
        params.lambda0,
    );
    p.hardening = params.hardening;
    p
}
/// Create a sand particle with Drucker-Prager plasticity.
pub fn make_sand_particle(
    position: [f64; 3],
    velocity: [f64; 3],
    volume0: f64,
    params: &SandParams,
) -> MpmParticle {
    let mass = params.rho0 * volume0;
    let mut p = MpmParticle::new(
        position,
        velocity,
        mass,
        volume0,
        MaterialType::Sand,
        params.mu0,
        params.lambda0,
    );
    p.yield_param = params.friction_angle;
    p
}
/// Compute the CFL-safe time step for MPM.
///
/// dt_cfl = CFL * h / v_max
pub fn compute_cfl_dt(particles: &[MpmParticle], h: f64, cfl: f64) -> f64 {
    let mut v_max = 1e-6;
    for p in particles {
        if !p.active {
            continue;
        }
        let v = len3(p.velocity);
        if v > v_max {
            v_max = v;
        }
    }
    cfl * h / v_max
}
/// G2P without dt (pure gather; dt used only by advection step).
pub fn g2p_transfer(grid: &MpmGrid, particles: &mut [MpmParticle]) {
    g2p_transfer_with_dt(grid, particles, 1.0);
}
