#![allow(clippy::needless_range_loop)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    BcEnforcementMethod, ConvergencePoint, ConvergenceRate, ConvergenceStudySummary,
    CouplingInterface, CouplingMethod, CouplingNode, EfgMaterial, EfgResult, ErrorIndicator,
    EssentialBc, KernelType, MeshfreeNode, MlsShape, RbfInterpolant, RbfType, RkpmShape,
    StressState, SupportDomain, TriCell,
};

/// 2-D point alias.
pub(super) type Point2 = [f64; 2];
/// 3-D point alias.
pub(super) type Point3 = [f64; 3];
/// Euclidean distance in 2-D.
#[allow(dead_code)]
pub(super) fn dist2(a: Point2, b: Point2) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}
/// Euclidean distance in 3-D.
#[allow(dead_code)]
pub(super) fn dist3(a: Point3, b: Point3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}
/// Dot product of two slices.
#[allow(dead_code)]
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Solve 2x2 linear system A x = b. Returns `None` if singular.
#[allow(dead_code)]
pub(super) fn solve2(a: [[f64; 2]; 2], b: [f64; 2]) -> Option<[f64; 2]> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-14 {
        return None;
    }
    Some([
        (b[0] * a[1][1] - b[1] * a[0][1]) / det,
        (a[0][0] * b[1] - a[1][0] * b[0]) / det,
    ])
}
/// Solve 3x3 linear system A x = b. Returns `None` if singular.
#[allow(dead_code)]
pub(super) fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv = 1.0 / det;
    let x0 = inv
        * (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            + b[1] * (a[0][2] * a[2][1] - a[0][1] * a[2][2])
            + b[2] * (a[0][1] * a[1][2] - a[0][2] * a[1][1]));
    let x1 = inv
        * (b[0] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
            + b[1] * (a[0][0] * a[2][2] - a[0][2] * a[2][0])
            + b[2] * (a[0][2] * a[1][0] - a[0][0] * a[1][2]));
    let x2 = inv
        * (b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
            + b[1] * (a[0][1] * a[2][0] - a[0][0] * a[2][1])
            + b[2] * (a[0][0] * a[1][1] - a[0][1] * a[1][0]));
    Some([x0, x1, x2])
}
/// Invert a symmetric 3x3 matrix. Returns `None` if singular.
#[allow(dead_code)]
pub(super) fn invert3(a: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv = 1.0 / det;
    Some([
        [
            inv * (a[1][1] * a[2][2] - a[1][2] * a[2][1]),
            inv * (a[0][2] * a[2][1] - a[0][1] * a[2][2]),
            inv * (a[0][1] * a[1][2] - a[0][2] * a[1][1]),
        ],
        [
            inv * (a[1][2] * a[2][0] - a[1][0] * a[2][2]),
            inv * (a[0][0] * a[2][2] - a[0][2] * a[2][0]),
            inv * (a[0][2] * a[1][0] - a[0][0] * a[1][2]),
        ],
        [
            inv * (a[1][0] * a[2][1] - a[1][1] * a[2][0]),
            inv * (a[0][1] * a[2][0] - a[0][0] * a[2][1]),
            inv * (a[0][0] * a[1][1] - a[0][1] * a[1][0]),
        ],
    ])
}
/// Multiply 3x3 matrix by 3-vector.
#[allow(dead_code)]
pub(super) fn mat_vec3(a: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2],
        a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2],
        a[2][0] * v[0] + a[2][1] * v[1] + a[2][2] * v[2],
    ]
}
/// Multiply 3x3 matrices C = A * B.
#[allow(dead_code)]
pub(super) fn mat_mul3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Evaluate a scalar kernel function at normalised distance `s = r / h`.
///
/// Returns `(w, dw_ds)` — value and derivative with respect to `s`.
#[allow(dead_code)]
pub fn kernel_eval(s: f64, ktype: KernelType) -> (f64, f64) {
    if !(0.0..=1.0).contains(&s) {
        return (0.0, 0.0);
    }
    match ktype {
        KernelType::CubicSpline => {
            if s <= 0.5 {
                let s2 = s * s;
                let w = 2.0 / 3.0 - 4.0 * s2 + 4.0 * s2 * s;
                let dw = -8.0 * s + 12.0 * s2;
                (w, dw)
            } else {
                let t = 1.0 - s;
                let w = (4.0 / 3.0) * t * t * t;
                let dw = -4.0 * t * t;
                (w, dw)
            }
        }
        KernelType::QuarticSpline => {
            let t = 1.0 - s * s;
            let w = t * t;
            let dw = -4.0 * s * t;
            (w, dw)
        }
        KernelType::Gaussian => {
            let beta = 4.0;
            let e1 = (-beta * s * s).exp();
            let e2 = (-beta).exp();
            let w = (e1 - e2) / (1.0 - e2);
            let dw = -2.0 * beta * s * e1 / (1.0 - e2);
            (w, dw)
        }
        KernelType::WendlandC2 => {
            let t = 1.0 - s;
            let w = t.powi(4) * (1.0 + 4.0 * s);
            let dw = -20.0 * s * t.powi(3);
            (w, dw)
        }
        KernelType::WendlandC4 => {
            let t = 1.0 - s;
            let w = t.powi(6) * (1.0 + 6.0 * s + 35.0 / 3.0 * s * s);
            let dw = t.powi(5)
                * (-6.0 * (1.0 + 6.0 * s + 35.0 / 3.0 * s * s) + t * (6.0 + 70.0 / 3.0 * s));
            (w, dw)
        }
    }
}
/// Compute MLS shape functions at `pt` using a linear basis `[1, x, y]`.
///
/// # Arguments
/// * `pt` — evaluation point
/// * `domain` — support domain (supplies nodes and neighbour search)
/// * `ktype` — kernel function type
///
/// Returns `None` if moment matrix is singular.
#[allow(dead_code)]
pub fn mls_shape_2d(pt: Point2, domain: &SupportDomain, ktype: KernelType) -> Option<MlsShape> {
    let nbrs = domain.neighbours(pt);
    let n = nbrs.len();
    if n < 3 {
        return None;
    }
    let mut mom = [[0.0_f64; 3]; 3];
    let mut dmom_dx = [[0.0_f64; 3]; 3];
    let mut dmom_dy = [[0.0_f64; 3]; 3];
    let mut weights = Vec::with_capacity(n);
    let mut dw_dx = Vec::with_capacity(n);
    let mut dw_dy = Vec::with_capacity(n);
    let mut basis_vals = Vec::with_capacity(n);
    for &idx in &nbrs {
        let node = &domain.nodes[idx];
        let h = node.support_radius * domain.dilation;
        let dx = pt[0] - node.pos[0];
        let dy = pt[1] - node.pos[1];
        let r = (dx * dx + dy * dy).sqrt().max(1e-15);
        let s = r / h;
        let (w, dw_ds) = kernel_eval(s.min(1.0), ktype);
        let ds_dr = 1.0 / h;
        let dr_dx = dx / r;
        let dr_dy = dy / r;
        let dwx = dw_ds * ds_dr * dr_dx;
        let dwy = dw_ds * ds_dr * dr_dy;
        weights.push(w);
        dw_dx.push(dwx);
        dw_dy.push(dwy);
        let p = [1.0, dx, dy];
        basis_vals.push(p);
        for a in 0..3 {
            for b in 0..3 {
                mom[a][b] += w * p[a] * p[b];
                dmom_dx[a][b] += dwx * p[a] * p[b];
                dmom_dy[a][b] += dwy * p[a] * p[b];
            }
        }
        dmom_dx[0][1] += w;
        dmom_dx[1][0] += w;
        dmom_dx[1][1] += 2.0 * w * 1.0;
    }
    dmom_dx = [[0.0; 3]; 3];
    dmom_dy = [[0.0; 3]; 3];
    for k in 0..n {
        let p = basis_vals[k];
        let w = weights[k];
        let wx = dw_dx[k];
        let wy = dw_dy[k];
        let dp_dx = [0.0, 1.0, 0.0];
        let dp_dy = [0.0, 0.0, 1.0];
        for a in 0..3 {
            for b in 0..3 {
                dmom_dx[a][b] += wx * p[a] * p[b] + w * dp_dx[a] * p[b] + w * p[a] * dp_dx[b];
                dmom_dy[a][b] += wy * p[a] * p[b] + w * dp_dy[a] * p[b] + w * p[a] * dp_dy[b];
            }
        }
    }
    let inv_m = invert3(mom)?;
    let p0 = [1.0, 0.0, 0.0];
    let gamma = mat_vec3(inv_m, p0);
    let mut phi = Vec::with_capacity(n);
    let mut dphi_dx_vec = Vec::with_capacity(n);
    let mut dphi_dy_vec = Vec::with_capacity(n);
    let dinv_dx = {
        let tmp = mat_mul3(inv_m, dmom_dx);
        let res = mat_mul3(tmp, inv_m);
        let mut neg = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg[i][j] = -res[i][j];
            }
        }
        neg
    };
    let dinv_dy = {
        let tmp = mat_mul3(inv_m, dmom_dy);
        let res = mat_mul3(tmp, inv_m);
        let mut neg = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg[i][j] = -res[i][j];
            }
        }
        neg
    };
    let dgamma_dx = mat_vec3(dinv_dx, p0);
    let dgamma_dy = mat_vec3(dinv_dy, p0);
    for k in 0..n {
        let p = basis_vals[k];
        let w = weights[k];
        let wx = dw_dx[k];
        let wy = dw_dy[k];
        let dp_dx = [0.0, 1.0, 0.0];
        let dp_dy = [0.0, 0.0, 1.0];
        let gp: f64 = gamma.iter().zip(p.iter()).map(|(a, b)| a * b).sum();
        phi.push(gp * w);
        let dgp_dx: f64 = dgamma_dx
            .iter()
            .zip(p.iter())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            + gamma
                .iter()
                .zip(dp_dx.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
        dphi_dx_vec.push(dgp_dx * w + gp * wx);
        let dgp_dy: f64 = dgamma_dy
            .iter()
            .zip(p.iter())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            + gamma
                .iter()
                .zip(dp_dy.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
        dphi_dy_vec.push(dgp_dy * w + gp * wy);
    }
    Some(MlsShape {
        phi,
        dphi_dx: dphi_dx_vec,
        dphi_dy: dphi_dy_vec,
        neighbours: nbrs,
    })
}
/// Evaluate an RBF at distance `r` with shape parameter `c`.
#[allow(dead_code)]
pub fn rbf_eval(r: f64, c: f64, rtype: RbfType) -> f64 {
    match rtype {
        RbfType::Multiquadric => (r * r + c * c).sqrt(),
        RbfType::InverseMultiquadric => 1.0 / (r * r + c * c).sqrt(),
        RbfType::Gaussian => (-(c * c) * r * r).exp(),
        RbfType::ThinPlateSpline => {
            if r.abs() < 1e-15 {
                0.0
            } else {
                r * r * r.abs().ln()
            }
        }
        RbfType::Polyharmonic3 => r * r * r,
    }
}
/// Solve a dense linear system (Gaussian elimination with partial pivoting).
///
/// Overwrites `a` (n x n stored row-major in flat vec) and `b`.
#[allow(dead_code)]
pub(super) fn dense_solve(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
    for col in 0..n {
        let mut max_val = a[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let v = a[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-14 {
            return false;
        }
        if max_row != col {
            for j in 0..n {
                a.swap(col * n + j, max_row * n + j);
            }
            b.swap(col, max_row);
        }
        let pivot = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            for j in col..n {
                let v = a[col * n + j];
                a[row * n + j] -= factor * v;
            }
            let bc = b[col];
            b[row] -= factor * bc;
        }
    }
    for col in (0..n).rev() {
        let mut s = b[col];
        for j in (col + 1)..n {
            s -= a[col * n + j] * b[j];
        }
        b[col] = s / a[col * n + col];
    }
    true
}
/// Build an RBF interpolant from scattered 2-D data without polynomial
/// augmentation (pure RBF).
///
/// # Arguments
/// * `centres` — data point positions
/// * `values` — data point scalar values
/// * `c` — shape parameter
/// * `rtype` — RBF type
#[allow(dead_code)]
pub fn rbf_interpolant(
    centres: &[Point2],
    values: &[f64],
    c: f64,
    rtype: RbfType,
) -> Option<RbfInterpolant> {
    let n = centres.len();
    if n == 0 || n != values.len() {
        return None;
    }
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let r = dist2(centres[i], centres[j]);
            a[i * n + j] = rbf_eval(r, c, rtype);
        }
    }
    let mut b = values.to_vec();
    if !dense_solve(&mut a, &mut b, n) {
        return None;
    }
    Some(RbfInterpolant {
        centres: centres.to_vec(),
        coeffs: b,
        c,
        rtype,
    })
}
/// Plane-stress constitutive matrix D (3x3).
#[allow(dead_code)]
pub(super) fn plane_stress_d(young: f64, poisson: f64) -> [[f64; 3]; 3] {
    let c = young / (1.0 - poisson * poisson);
    [
        [c, c * poisson, 0.0],
        [c * poisson, c, 0.0],
        [0.0, 0.0, c * (1.0 - poisson) / 2.0],
    ]
}
/// Assemble global stiffness for 2-D plane-stress EFG using triangular
/// background cells with 1-point Gauss quadrature.
///
/// Returns the stiffness matrix as a dense row-major flat vector and the
/// force vector.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn efg_assemble_2d(
    domain: &SupportDomain,
    cells: &[TriCell],
    mat: &EfgMaterial,
    ktype: KernelType,
    body_force: Point2,
) -> Option<(Vec<f64>, Vec<f64>)> {
    let ndof = domain.node_count() * 2;
    let mut stiffness = vec![0.0; ndof * ndof];
    let mut force = vec![0.0; ndof];
    let d_mat = plane_stress_d(mat.young, mat.poisson);
    for cell in cells {
        let area = cell.area();
        let cen = cell.centroid();
        let shape = mls_shape_2d(cen, domain, ktype)?;
        let nn = shape.neighbours.len();
        for a in 0..nn {
            let ia = shape.neighbours[a];
            let b_a = [
                [shape.dphi_dx[a], 0.0],
                [0.0, shape.dphi_dy[a]],
                [shape.dphi_dy[a], shape.dphi_dx[a]],
            ];
            force[ia * 2] += shape.phi[a] * body_force[0] * area * mat.thickness;
            force[ia * 2 + 1] += shape.phi[a] * body_force[1] * area * mat.thickness;
            for b in 0..nn {
                let ib = shape.neighbours[b];
                let b_b = [
                    [shape.dphi_dx[b], 0.0],
                    [0.0, shape.dphi_dy[b]],
                    [shape.dphi_dy[b], shape.dphi_dx[b]],
                ];
                for r in 0..2 {
                    for cc in 0..2 {
                        let mut val = 0.0;
                        for s in 0..3 {
                            for t in 0..3 {
                                val += b_a[s][r] * d_mat[s][t] * b_b[t][cc];
                            }
                        }
                        stiffness[(ia * 2 + r) * ndof + ib * 2 + cc] += val * area * mat.thickness;
                    }
                }
            }
        }
    }
    Some((stiffness, force))
}
/// Compute RKPM shape functions at `pt` using a linear correction function
/// and kernel `ktype`.
///
/// The RKPM correction ensures zeroth- and first-order consistency:
///   Psi_I(x) = C(x; x-x_I) * Phi_a(x-x_I)
/// where Phi_a is the kernel and C is a correction polynomial.
#[allow(dead_code)]
pub fn rkpm_shape_2d(pt: Point2, domain: &SupportDomain, ktype: KernelType) -> Option<RkpmShape> {
    let nbrs = domain.neighbours(pt);
    let n = nbrs.len();
    if n < 3 {
        return None;
    }
    let mut mom = [[0.0_f64; 3]; 3];
    let mut dmom_dx = [[0.0_f64; 3]; 3];
    let mut dmom_dy = [[0.0_f64; 3]; 3];
    pub(super) struct NodeData {
        w: f64,
        dwx: f64,
        dwy: f64,
        p: [f64; 3],
    }
    let mut ndata = Vec::with_capacity(n);
    for &idx in &nbrs {
        let node = &domain.nodes[idx];
        let h = node.support_radius * domain.dilation;
        let dx = pt[0] - node.pos[0];
        let dy = pt[1] - node.pos[1];
        let r = (dx * dx + dy * dy).sqrt().max(1e-15);
        let s = (r / h).min(1.0);
        let (w, dw_ds) = kernel_eval(s, ktype);
        let ds_dr = 1.0 / h;
        let dr_dx = dx / r;
        let dr_dy = dy / r;
        let dwx = dw_ds * ds_dr * dr_dx;
        let dwy = dw_ds * ds_dr * dr_dy;
        let p = [1.0, dx, dy];
        ndata.push(NodeData { w, dwx, dwy, p });
        let dp_dx = [0.0, 1.0, 0.0];
        let dp_dy = [0.0, 0.0, 1.0];
        for a in 0..3 {
            for b in 0..3 {
                mom[a][b] += w * p[a] * p[b];
                dmom_dx[a][b] += dwx * p[a] * p[b] + w * dp_dx[a] * p[b] + w * p[a] * dp_dx[b];
                dmom_dy[a][b] += dwy * p[a] * p[b] + w * dp_dy[a] * p[b] + w * p[a] * dp_dy[b];
            }
        }
    }
    let inv_m = invert3(mom)?;
    let p0 = [1.0, 0.0, 0.0];
    let dinv_dx = {
        let tmp = mat_mul3(inv_m, dmom_dx);
        let res = mat_mul3(tmp, inv_m);
        let mut neg = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg[i][j] = -res[i][j];
            }
        }
        neg
    };
    let dinv_dy = {
        let tmp = mat_mul3(inv_m, dmom_dy);
        let res = mat_mul3(tmp, inv_m);
        let mut neg = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg[i][j] = -res[i][j];
            }
        }
        neg
    };
    let gamma = mat_vec3(inv_m, p0);
    let dgamma_dx = mat_vec3(dinv_dx, p0);
    let dgamma_dy = mat_vec3(dinv_dy, p0);
    let mut psi = Vec::with_capacity(n);
    let mut dpsi_dx_v = Vec::with_capacity(n);
    let mut dpsi_dy_v = Vec::with_capacity(n);
    for k in 0..n {
        let nd = &ndata[k];
        let gp: f64 = gamma.iter().zip(nd.p.iter()).map(|(a, b)| a * b).sum();
        psi.push(gp * nd.w);
        let dp_dx = [0.0, 1.0, 0.0];
        let dp_dy = [0.0, 0.0, 1.0];
        let dgp_dx: f64 = dgamma_dx
            .iter()
            .zip(nd.p.iter())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            + gamma
                .iter()
                .zip(dp_dx.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
        dpsi_dx_v.push(dgp_dx * nd.w + gp * nd.dwx);
        let dgp_dy: f64 = dgamma_dy
            .iter()
            .zip(nd.p.iter())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            + gamma
                .iter()
                .zip(dp_dy.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
        dpsi_dy_v.push(dgp_dy * nd.w + gp * nd.dwy);
    }
    Some(RkpmShape {
        psi,
        dpsi_dx: dpsi_dx_v,
        dpsi_dy: dpsi_dy_v,
        neighbours: nbrs,
    })
}
/// Generate a regular grid of triangular integration cells over a rectangular
/// domain `[x0, x1] x [y0, y1]` with `nx` by `ny` divisions.
///
/// Each rectangular sub-cell is split into two triangles.
#[allow(dead_code)]
pub fn generate_tri_cells(
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    nx: usize,
    ny: usize,
) -> Vec<TriCell> {
    let hx = (x1 - x0) / nx as f64;
    let hy = (y1 - y0) / ny as f64;
    let mut cells = Vec::with_capacity(nx * ny * 2);
    for i in 0..nx {
        for j in 0..ny {
            let lx = x0 + i as f64 * hx;
            let ly = y0 + j as f64 * hy;
            let rx = lx + hx;
            let uy = ly + hy;
            cells.push(TriCell {
                vertices: [[lx, ly], [rx, ly], [lx, uy]],
            });
            cells.push(TriCell {
                vertices: [[rx, ly], [rx, uy], [lx, uy]],
            });
        }
    }
    cells
}
/// Generate a Voronoi-style cell list from a set of nodes (simple 2-D).
///
/// This is a simplified version that creates a triangulation by connecting
/// each node to its nearest neighbours. Returns triangle cells.
#[allow(dead_code)]
pub fn generate_voronoi_cells(nodes: &[Point2], _max_neighbours: usize) -> Vec<TriCell> {
    let n = nodes.len();
    if n < 3 {
        return Vec::new();
    }
    let mut cells = Vec::new();
    for i in 0..n {
        let mut dists: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, dist2(nodes[i], nodes[j])))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        if dists.len() >= 2 {
            let j = dists[0].0;
            let k = dists[1].0;
            if i < j && i < k {
                cells.push(TriCell {
                    vertices: [nodes[i], nodes[j], nodes[k]],
                });
            }
        }
    }
    cells
}
/// Apply essential BCs via the penalty method to a dense stiffness matrix
/// and force vector.
///
/// Adds `alpha * I` to the diagonal and `alpha * g` to the force at
/// constrained DOFs.
#[allow(dead_code)]
pub fn apply_penalty_bc(
    stiffness: &mut [f64],
    force: &mut [f64],
    ndof: usize,
    bcs: &[EssentialBc],
    penalty: f64,
) {
    for bc in bcs {
        let idx = bc.node * 2 + bc.dof;
        if idx < ndof {
            stiffness[idx * ndof + idx] += penalty;
            force[idx] += penalty * bc.value;
        }
    }
}
/// Apply essential BCs via the Lagrange multiplier method.
///
/// Augments the system from (K, f) of size n to size (n + m) where m is
/// the number of constraints. Returns the augmented system.
#[allow(dead_code)]
pub fn apply_lagrange_bc(
    stiffness: &[f64],
    force: &[f64],
    ndof: usize,
    bcs: &[EssentialBc],
) -> (Vec<f64>, Vec<f64>) {
    let m = bcs.len();
    let n_aug = ndof + m;
    let mut k_aug = vec![0.0; n_aug * n_aug];
    let mut f_aug = vec![0.0; n_aug];
    for i in 0..ndof {
        for j in 0..ndof {
            k_aug[i * n_aug + j] = stiffness[i * ndof + j];
        }
        f_aug[i] = force[i];
    }
    for (c, bc) in bcs.iter().enumerate() {
        let idx = bc.node * 2 + bc.dof;
        if idx < ndof {
            k_aug[(ndof + c) * n_aug + idx] = 1.0;
            k_aug[idx * n_aug + ndof + c] = 1.0;
            f_aug[ndof + c] = bc.value;
        }
    }
    (k_aug, f_aug)
}
/// Apply essential BCs via Nitsche's method.
///
/// Modifies the stiffness matrix and force vector in place.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn apply_nitsche_bc(
    stiffness: &mut [f64],
    force: &mut [f64],
    ndof: usize,
    bcs: &[EssentialBc],
    nitsche_param: f64,
) {
    for bc in bcs {
        let idx = bc.node * 2 + bc.dof;
        if idx < ndof {
            stiffness[idx * ndof + idx] += nitsche_param;
            force[idx] += nitsche_param * bc.value;
        }
    }
}
/// Compute stress at a single evaluation point from nodal displacements
/// using MLS shape function derivatives and a plane-stress constitutive
/// matrix.
#[allow(dead_code)]
pub fn recover_stress_at_point(
    pt: Point2,
    displacements: &[f64],
    domain: &SupportDomain,
    ktype: KernelType,
    mat: &EfgMaterial,
) -> Option<StressState> {
    let shape = mls_shape_2d(pt, domain, ktype)?;
    let d_mat = plane_stress_d(mat.young, mat.poisson);
    let mut strain = [0.0_f64; 3];
    for (k, &idx) in shape.neighbours.iter().enumerate() {
        let ux = displacements[idx * 2];
        let uy = displacements[idx * 2 + 1];
        strain[0] += shape.dphi_dx[k] * ux;
        strain[1] += shape.dphi_dy[k] * uy;
        strain[2] += shape.dphi_dy[k] * ux + shape.dphi_dx[k] * uy;
    }
    let sx = d_mat[0][0] * strain[0] + d_mat[0][1] * strain[1];
    let sy = d_mat[1][0] * strain[0] + d_mat[1][1] * strain[1];
    let txy = d_mat[2][2] * strain[2];
    Some(StressState {
        sigma_xx: sx,
        sigma_yy: sy,
        tau_xy: txy,
    })
}
/// Perform nodal stress recovery by averaging stresses from nearby
/// integration points.
#[allow(dead_code)]
pub fn nodal_stress_recovery(
    domain: &SupportDomain,
    cells: &[TriCell],
    displacements: &[f64],
    ktype: KernelType,
    mat: &EfgMaterial,
) -> Vec<StressState> {
    let n = domain.node_count();
    let mut stresses = vec![StressState::default(); n];
    let mut counts = vec![0usize; n];
    for cell in cells {
        let cen = cell.centroid();
        if let Some(s) = recover_stress_at_point(cen, displacements, domain, ktype, mat) {
            let nbrs = domain.neighbours(cen);
            for idx in nbrs {
                stresses[idx].sigma_xx += s.sigma_xx;
                stresses[idx].sigma_yy += s.sigma_yy;
                stresses[idx].tau_xy += s.tau_xy;
                counts[idx] += 1;
            }
        }
    }
    for i in 0..n {
        if counts[i] > 0 {
            let c = counts[i] as f64;
            stresses[i].sigma_xx /= c;
            stresses[i].sigma_yy /= c;
            stresses[i].tau_xy /= c;
        }
    }
    stresses
}
/// Superconvergent Patch Recovery (SPR) for meshfree stress recovery.
///
/// Fits a polynomial to integration point stresses in a patch around each
/// node, then evaluates at the node position. Uses a linear fit
/// `sigma = a0 + a1*x + a2*y` within each patch.
#[allow(dead_code)]
pub fn spr_stress_recovery(
    domain: &SupportDomain,
    cells: &[TriCell],
    displacements: &[f64],
    ktype: KernelType,
    mat: &EfgMaterial,
) -> Vec<StressState> {
    let n = domain.node_count();
    let mut stresses = vec![StressState::default(); n];
    let mut ip_stresses = Vec::new();
    let mut ip_positions = Vec::new();
    for cell in cells {
        let cen = cell.centroid();
        if let Some(s) = recover_stress_at_point(cen, displacements, domain, ktype, mat) {
            ip_stresses.push(s);
            ip_positions.push(cen);
        }
    }
    for node_idx in 0..n {
        let node_pos = domain.nodes[node_idx].pos;
        let h = domain.nodes[node_idx].support_radius * domain.dilation;
        let mut patch_pts = Vec::new();
        let mut patch_stresses = Vec::new();
        for (k, pos) in ip_positions.iter().enumerate() {
            if dist2(node_pos, *pos) < h {
                patch_pts.push(*pos);
                patch_stresses.push(ip_stresses[k]);
            }
        }
        if patch_pts.len() < 3 {
            if !patch_stresses.is_empty() {
                let c = patch_stresses.len() as f64;
                for ps in &patch_stresses {
                    stresses[node_idx].sigma_xx += ps.sigma_xx;
                    stresses[node_idx].sigma_yy += ps.sigma_yy;
                    stresses[node_idx].tau_xy += ps.tau_xy;
                }
                stresses[node_idx].sigma_xx /= c;
                stresses[node_idx].sigma_yy /= c;
                stresses[node_idx].tau_xy /= c;
            }
            continue;
        }
        let m = patch_pts.len();
        let mut ptp = [[0.0_f64; 3]; 3];
        let mut pts_xx = [0.0_f64; 3];
        let mut pts_yy = [0.0_f64; 3];
        let mut pts_xy = [0.0_f64; 3];
        for k in 0..m {
            let px = patch_pts[k][0] - node_pos[0];
            let py = patch_pts[k][1] - node_pos[1];
            let p = [1.0, px, py];
            for a in 0..3 {
                for b in 0..3 {
                    ptp[a][b] += p[a] * p[b];
                }
                pts_xx[a] += p[a] * patch_stresses[k].sigma_xx;
                pts_yy[a] += p[a] * patch_stresses[k].sigma_yy;
                pts_xy[a] += p[a] * patch_stresses[k].tau_xy;
            }
        }
        if let Some(a) = solve3(ptp, pts_xx) {
            stresses[node_idx].sigma_xx = a[0];
        }
        if let Some(a) = solve3(ptp, pts_yy) {
            stresses[node_idx].sigma_yy = a[0];
        }
        if let Some(a) = solve3(ptp, pts_xy) {
            stresses[node_idx].tau_xy = a[0];
        }
    }
    stresses
}
/// Compute error indicators based on stress gradient for adaptive refinement.
///
/// Nodes with high stress gradients are candidates for refinement (adding
/// new particles nearby).
#[allow(dead_code)]
pub fn compute_error_indicators(
    domain: &SupportDomain,
    stresses: &[StressState],
) -> Vec<ErrorIndicator> {
    let n = domain.node_count();
    let mut indicators = Vec::with_capacity(n);
    for i in 0..n {
        let nbrs = domain.neighbours(domain.nodes[i].pos);
        let vm_i = stresses[i].von_mises();
        let mut max_grad = 0.0_f64;
        for &j in &nbrs {
            if j == i {
                continue;
            }
            let d = dist2(domain.nodes[i].pos, domain.nodes[j].pos);
            if d > 1e-15 {
                let grad = (stresses[j].von_mises() - vm_i).abs() / d;
                max_grad = max_grad.max(grad);
            }
        }
        indicators.push(ErrorIndicator {
            node: i,
            error: max_grad,
        });
    }
    indicators
}
/// Refine the meshfree domain by inserting new particles at mid-points
/// between nodes that exceed the error threshold.
///
/// Returns a new set of nodes (original plus inserted).
#[allow(dead_code)]
pub fn adaptive_refine(
    domain: &SupportDomain,
    indicators: &[ErrorIndicator],
    threshold: f64,
) -> Vec<MeshfreeNode> {
    let mut new_nodes: Vec<MeshfreeNode> = domain.nodes.clone();
    for ind in indicators {
        if ind.error > threshold {
            let i = ind.node;
            let nbrs = domain.neighbours(domain.nodes[i].pos);
            for &j in &nbrs {
                if j <= i {
                    continue;
                }
                let mid = [
                    (domain.nodes[i].pos[0] + domain.nodes[j].pos[0]) / 2.0,
                    (domain.nodes[i].pos[1] + domain.nodes[j].pos[1]) / 2.0,
                ];
                let avg_h =
                    (domain.nodes[i].support_radius + domain.nodes[j].support_radius) / 2.0 * 0.75;
                new_nodes.push(MeshfreeNode {
                    pos: mid,
                    support_radius: avg_h,
                    value: (domain.nodes[i].value + domain.nodes[j].value) / 2.0,
                });
            }
        }
    }
    new_nodes
}
/// Compute the density of nodes in a region and flag under-resolved areas.
#[allow(dead_code)]
pub fn node_density_map(
    domain: &SupportDomain,
    grid_x: usize,
    grid_y: usize,
    bounds: [f64; 4],
) -> Vec<f64> {
    let nx = grid_x.max(1);
    let ny = grid_y.max(1);
    let hx = (bounds[2] - bounds[0]) / nx as f64;
    let hy = (bounds[3] - bounds[1]) / ny as f64;
    let mut density = vec![0.0; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            let pt = [
                bounds[0] + (i as f64 + 0.5) * hx,
                bounds[1] + (j as f64 + 0.5) * hy,
            ];
            density[i * ny + j] = domain.neighbours(pt).len() as f64;
        }
    }
    density
}
/// Set up bridging domain coupling between a meshfree domain and FEM mesh.
///
/// Identifies nodes in the overlap region and assigns blending weights.
#[allow(dead_code)]
pub fn setup_bridging_domain(
    meshfree_domain: &SupportDomain,
    fem_nodes: &[Point2],
    overlap_x_min: f64,
    overlap_x_max: f64,
) -> CouplingInterface {
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1e6);
    let width = (overlap_x_max - overlap_x_min).max(1e-15);
    for (i, node) in meshfree_domain.nodes.iter().enumerate() {
        if node.pos[0] >= overlap_x_min && node.pos[0] <= overlap_x_max {
            let w = (node.pos[0] - overlap_x_min) / width;
            interface.add_node(CouplingNode {
                pos: node.pos,
                meshfree_idx: Some(i),
                fem_idx: None,
                blend_weight: w,
            });
        }
    }
    for (j, &pos) in fem_nodes.iter().enumerate() {
        if pos[0] >= overlap_x_min && pos[0] <= overlap_x_max {
            let w = (pos[0] - overlap_x_min) / width;
            interface.add_node(CouplingNode {
                pos,
                meshfree_idx: None,
                fem_idx: Some(j),
                blend_weight: 1.0 - w,
            });
        }
    }
    interface
}
/// Solve a dense linear system K u = f with the assembled EFG matrices.
#[allow(dead_code)]
pub fn solve_dense(stiffness: &[f64], force: &[f64], ndof: usize) -> Option<Vec<f64>> {
    let mut a = stiffness.to_vec();
    let mut b = force.to_vec();
    if dense_solve(&mut a, &mut b, ndof) {
        Some(b)
    } else {
        None
    }
}
/// Run a complete 2-D plane-stress EFG analysis.
///
/// # Arguments
/// * `domain` — support domain with nodes
/// * `cells` — background integration cells
/// * `mat` — material properties
/// * `ktype` — kernel function type
/// * `body_force` — body force per unit volume
/// * `bcs` — essential boundary conditions
/// * `bc_method` — boundary condition enforcement method
/// * `penalty` — penalty parameter (used for Penalty and Nitsche methods)
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn efg_analysis_2d(
    domain: &SupportDomain,
    cells: &[TriCell],
    mat: &EfgMaterial,
    ktype: KernelType,
    body_force: Point2,
    bcs: &[EssentialBc],
    bc_method: BcEnforcementMethod,
    penalty: f64,
) -> Option<EfgResult> {
    let (mut stiffness, mut force) = efg_assemble_2d(domain, cells, mat, ktype, body_force)?;
    let ndof = domain.node_count() * 2;
    match bc_method {
        BcEnforcementMethod::Penalty => {
            apply_penalty_bc(&mut stiffness, &mut force, ndof, bcs, penalty);
            let u = solve_dense(&stiffness, &force, ndof)?;
            Some(EfgResult {
                displacements: u,
                num_nodes: domain.node_count(),
            })
        }
        BcEnforcementMethod::LagrangeMultiplier => {
            let (k_aug, f_aug) = apply_lagrange_bc(&stiffness, &force, ndof, bcs);
            let n_aug = ndof + bcs.len();
            let u_aug = solve_dense(&k_aug, &f_aug, n_aug)?;
            Some(EfgResult {
                displacements: u_aug[..ndof].to_vec(),
                num_nodes: domain.node_count(),
            })
        }
        BcEnforcementMethod::Nitsche => {
            apply_nitsche_bc(&mut stiffness, &mut force, ndof, bcs, penalty);
            let u = solve_dense(&stiffness, &force, ndof)?;
            Some(EfgResult {
                displacements: u,
                num_nodes: domain.node_count(),
            })
        }
    }
}
/// Compute the condition number estimate of a dense matrix (ratio of
/// max to min diagonal element magnitudes — rough estimate).
#[allow(dead_code)]
pub fn condition_estimate(matrix: &[f64], n: usize) -> f64 {
    let mut max_diag = 0.0_f64;
    let mut min_diag = f64::MAX;
    for i in 0..n {
        let v = matrix[i * n + i].abs();
        if v > max_diag {
            max_diag = v;
        }
        if v < min_diag && v > 1e-30 {
            min_diag = v;
        }
    }
    if min_diag < 1e-30 {
        f64::INFINITY
    } else {
        max_diag / min_diag
    }
}
/// Compute the L2 norm of a vector.
#[allow(dead_code)]
pub fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Compute the infinity norm of a vector.
#[allow(dead_code)]
pub fn linf_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0_f64, f64::max)
}
/// Generate a regular grid of meshfree nodes over `[x0,x1] x [y0,y1]`.
#[allow(dead_code)]
pub fn generate_regular_nodes(
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    nx: usize,
    ny: usize,
    support_factor: f64,
) -> Vec<MeshfreeNode> {
    let hx = (x1 - x0) / (nx.max(1) - 1).max(1) as f64;
    let hy = (y1 - y0) / (ny.max(1) - 1).max(1) as f64;
    let h = hx.max(hy) * support_factor;
    let mut nodes = Vec::with_capacity(nx * ny);
    for i in 0..nx {
        for j in 0..ny {
            nodes.push(MeshfreeNode {
                pos: [x0 + i as f64 * hx, y0 + j as f64 * hy],
                support_radius: h,
                value: 0.0,
            });
        }
    }
    nodes
}
/// Generate a set of scattered nodes using a random perturbation of a
/// regular grid.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn generate_scattered_nodes(
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    nx: usize,
    ny: usize,
    support_factor: f64,
    perturbation: f64,
) -> Vec<MeshfreeNode> {
    use rand::RngExt;
    let hx = (x1 - x0) / (nx.max(1) - 1).max(1) as f64;
    let hy = (y1 - y0) / (ny.max(1) - 1).max(1) as f64;
    let h = hx.max(hy) * support_factor;
    let mut rng = rand::rng();
    let mut nodes = Vec::with_capacity(nx * ny);
    for i in 0..nx {
        for j in 0..ny {
            let px: f64 = rng.random_range(-perturbation..perturbation) * hx;
            let py: f64 = rng.random_range(-perturbation..perturbation) * hy;
            nodes.push(MeshfreeNode {
                pos: [x0 + i as f64 * hx + px, y0 + j as f64 * hy + py],
                support_radius: h,
                value: 0.0,
            });
        }
    }
    nodes
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshfree_fem::*;
    fn make_test_domain() -> SupportDomain {
        let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 3, 3, 2.5);
        SupportDomain::new(nodes, 1.0)
    }
    fn make_fine_domain() -> SupportDomain {
        let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 5, 5, 2.5);
        SupportDomain::new(nodes, 1.0)
    }
    #[test]
    fn test_dist2_basic() {
        let d = dist2([0.0, 0.0], [3.0, 4.0]);
        assert!((d - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_dist3_basic() {
        let d = dist3([0.0, 0.0, 0.0], [1.0, 2.0, 2.0]);
        assert!((d - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_dot_product() {
        let v = dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]);
        assert!((v - 32.0).abs() < 1e-12);
    }
    #[test]
    fn test_solve2() {
        let sol = solve2([[2.0, 3.0], [1.0, 1.0]], [8.0, 3.0]).unwrap();
        assert!((sol[0] - 1.0).abs() < 1e-12);
        assert!((sol[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_solve2_singular() {
        assert!(solve2([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0]).is_none());
    }
    #[test]
    fn test_solve3() {
        let a = [[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 2.0]];
        let b = [5.0, 10.0, 5.0];
        let x = solve3(a, b).unwrap();
        for i in 0..3 {
            let row_sum: f64 = (0..3).map(|j| a[i][j] * x[j]).sum();
            assert!((row_sum - b[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_invert3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = invert3(id).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_mat_vec3() {
        let a = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let v = [1.0, 2.0, 3.0];
        let r = mat_vec3(a, v);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 4.0).abs() < 1e-12);
        assert!((r[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_mat_mul3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let c = mat_mul3(id, a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_kernel_cubic_spline_at_zero() {
        let (w, _dw) = kernel_eval(0.0, KernelType::CubicSpline);
        assert!(w > 0.5, "Cubic spline at 0 should be > 0.5, got {w}");
    }
    #[test]
    fn test_kernel_cubic_spline_at_one() {
        let (w, _dw) = kernel_eval(1.0, KernelType::CubicSpline);
        assert!(w.abs() < 1e-12, "Cubic spline at 1 should be ~0, got {w}");
    }
    #[test]
    fn test_kernel_quartic_spline() {
        let (w, _dw) = kernel_eval(0.0, KernelType::QuarticSpline);
        assert!((w - 1.0).abs() < 1e-12);
        let (w1, _) = kernel_eval(1.0, KernelType::QuarticSpline);
        assert!(w1.abs() < 1e-12);
    }
    #[test]
    fn test_kernel_gaussian() {
        let (w, _) = kernel_eval(0.0, KernelType::Gaussian);
        assert!(
            (w - 1.0).abs() < 1e-6,
            "Gaussian at 0 should be ~1, got {w}"
        );
    }
    #[test]
    fn test_kernel_wendland_c2() {
        let (w, _) = kernel_eval(0.0, KernelType::WendlandC2);
        assert!((w - 1.0).abs() < 1e-12);
        let (w1, _) = kernel_eval(1.0, KernelType::WendlandC2);
        assert!(w1.abs() < 1e-12);
    }
    #[test]
    fn test_kernel_wendland_c4() {
        let (w, _) = kernel_eval(0.0, KernelType::WendlandC4);
        assert!((w - 1.0).abs() < 1e-12);
        let (w1, _) = kernel_eval(1.0, KernelType::WendlandC4);
        assert!(w1.abs() < 1e-12);
    }
    #[test]
    fn test_kernel_outside_support() {
        for kt in [
            KernelType::CubicSpline,
            KernelType::QuarticSpline,
            KernelType::Gaussian,
            KernelType::WendlandC2,
            KernelType::WendlandC4,
        ] {
            let (w, dw) = kernel_eval(1.5, kt);
            assert!((w).abs() < 1e-12);
            assert!((dw).abs() < 1e-12);
        }
    }
    #[test]
    fn test_support_domain_neighbours() {
        let domain = make_test_domain();
        let nbrs = domain.neighbours([0.5, 0.5]);
        assert!(!nbrs.is_empty(), "Centre of domain should have neighbours");
    }
    #[test]
    fn test_support_domain_sorted() {
        let domain = make_test_domain();
        let pairs = domain.neighbours_sorted([0.5, 0.5]);
        for w in pairs.windows(2) {
            assert!(w[0].1 <= w[1].1, "Should be sorted by distance");
        }
    }
    #[test]
    fn test_generate_regular_nodes() {
        let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 4, 4, 2.0);
        assert_eq!(nodes.len(), 16);
    }
    #[test]
    fn test_generate_scattered_nodes() {
        let nodes = generate_scattered_nodes(0.0, 0.0, 1.0, 1.0, 4, 4, 2.0, 0.1);
        assert_eq!(nodes.len(), 16);
    }
    #[test]
    fn test_tri_cell_area() {
        let cell = TriCell {
            vertices: [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        };
        assert!((cell.area() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_tri_cell_centroid() {
        let cell = TriCell {
            vertices: [[0.0, 0.0], [3.0, 0.0], [0.0, 3.0]],
        };
        let c = cell.centroid();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_generate_tri_cells() {
        let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 2, 2);
        assert_eq!(cells.len(), 8);
    }
    #[test]
    fn test_rbf_eval_gaussian() {
        let v = rbf_eval(0.0, 1.0, RbfType::Gaussian);
        assert!((v - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rbf_eval_multiquadric() {
        let v = rbf_eval(0.0, 2.0, RbfType::Multiquadric);
        assert!((v - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_rbf_eval_tps_zero() {
        let v = rbf_eval(0.0, 1.0, RbfType::ThinPlateSpline);
        assert!(v.abs() < 1e-12);
    }
    #[test]
    fn test_rbf_interpolant_constant() {
        let centres = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let values = vec![5.0, 5.0, 5.0, 5.0];
        let interp = rbf_interpolant(&centres, &values, 1.0, RbfType::Gaussian).unwrap();
        let v = interp.eval([0.5, 0.5]);
        assert!(
            (v - 5.0).abs() < 2.0,
            "RBF constant interpolation at centre should be ~5, got {v}"
        );
    }
    #[test]
    fn test_mls_shape_partition_of_unity() {
        let domain = make_fine_domain();
        let pt = [0.5, 0.5];
        let shape = mls_shape_2d(pt, &domain, KernelType::QuarticSpline).unwrap();
        let sum: f64 = shape.phi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "MLS shape functions should sum to 1, got {sum}"
        );
    }
    #[test]
    fn test_mls_shape_derivative_sum_zero() {
        let domain = make_fine_domain();
        let pt = [0.5, 0.5];
        let shape = mls_shape_2d(pt, &domain, KernelType::QuarticSpline).unwrap();
        let sum_dx: f64 = shape.dphi_dx.iter().sum();
        let sum_dy: f64 = shape.dphi_dy.iter().sum();
        assert!(
            sum_dx.abs() < 1e-4,
            "Sum of dphi/dx should be ~0, got {sum_dx}"
        );
        assert!(
            sum_dy.abs() < 1e-4,
            "Sum of dphi/dy should be ~0, got {sum_dy}"
        );
    }
    #[test]
    fn test_rkpm_shape_partition_of_unity() {
        let domain = make_fine_domain();
        let pt = [0.5, 0.5];
        let shape = rkpm_shape_2d(pt, &domain, KernelType::QuarticSpline).unwrap();
        let sum: f64 = shape.psi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "RKPM shape functions should sum to 1, got {sum}"
        );
    }
    #[test]
    fn test_plane_stress_d_symmetric() {
        let d = plane_stress_d(200e9, 0.3);
        assert!((d[0][1] - d[1][0]).abs() < 1e-3);
    }
    #[test]
    fn test_apply_penalty_bc() {
        let ndof = 4;
        let mut k = vec![0.0; ndof * ndof];
        let mut f = vec![0.0; ndof];
        let bcs = vec![EssentialBc {
            node: 0,
            dof: 0,
            value: 1.0,
        }];
        apply_penalty_bc(&mut k, &mut f, ndof, &bcs, 1e10);
        assert!(k[0] > 1e9);
        assert!(f[0] > 1e9);
    }
    #[test]
    fn test_apply_lagrange_bc() {
        let ndof = 4;
        let k = vec![1.0; ndof * ndof];
        let f = vec![0.0; ndof];
        let bcs = vec![EssentialBc {
            node: 0,
            dof: 0,
            value: 2.0,
        }];
        let (k_aug, f_aug) = apply_lagrange_bc(&k, &f, ndof, &bcs);
        assert_eq!(k_aug.len(), 25);
        assert_eq!(f_aug.len(), 5);
        assert!((f_aug[4] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_apply_nitsche_bc() {
        let ndof = 4;
        let mut k = vec![0.0; ndof * ndof];
        let mut f = vec![0.0; ndof];
        let bcs = vec![EssentialBc {
            node: 1,
            dof: 1,
            value: 3.0,
        }];
        apply_nitsche_bc(&mut k, &mut f, ndof, &bcs, 1e8);
        let idx = 2 + 1;
        assert!(k[idx * ndof + idx] > 1e7);
        assert!((f[idx] - 3e8).abs() < 1e-3);
    }
    #[test]
    fn test_stress_von_mises() {
        let s = StressState {
            sigma_xx: 100.0,
            sigma_yy: 0.0,
            tau_xy: 0.0,
        };
        assert!((s.von_mises() - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_stress_principal() {
        let s = StressState {
            sigma_xx: 100.0,
            sigma_yy: 0.0,
            tau_xy: 0.0,
        };
        let (p1, p2) = s.principal();
        assert!((p1 - 100.0).abs() < 1e-10);
        assert!(p2.abs() < 1e-10);
    }
    #[test]
    fn test_condition_estimate() {
        let mat = vec![10.0, 0.0, 0.0, 1.0];
        let cond = condition_estimate(&mat, 2);
        assert!((cond - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_l2_norm() {
        assert!((l2_norm(&[3.0, 4.0]) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_linf_norm() {
        assert!((linf_norm(&[-5.0, 3.0, 4.0]) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_dense_solve_identity() {
        let mut a = vec![1.0, 0.0, 0.0, 1.0];
        let mut b = vec![3.0, 7.0];
        assert!(dense_solve(&mut a, &mut b, 2));
        assert!((b[0] - 3.0).abs() < 1e-12);
        assert!((b[1] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_coupling_interface_new() {
        let ci = CouplingInterface::new(CouplingMethod::Direct, 1e6);
        assert_eq!(ci.num_coupling_nodes(), 0);
        assert_eq!(ci.method, CouplingMethod::Direct);
    }
    #[test]
    fn test_coupling_add_node() {
        let mut ci = CouplingInterface::new(CouplingMethod::Direct, 1e6);
        ci.add_node(CouplingNode {
            pos: [0.5, 0.5],
            meshfree_idx: Some(0),
            fem_idx: Some(0),
            blend_weight: 0.5,
        });
        assert_eq!(ci.num_coupling_nodes(), 1);
    }
    #[test]
    fn test_coupling_penalty_matrix() {
        let mut ci = CouplingInterface::new(CouplingMethod::Direct, 1e6);
        ci.add_node(CouplingNode {
            pos: [0.5, 0.5],
            meshfree_idx: Some(0),
            fem_idx: Some(0),
            blend_weight: 0.5,
        });
        let k = ci.compute_penalty_coupling(4, 4);
        assert_eq!(k.len(), 64);
    }
    #[test]
    fn test_bridging_domain_setup() {
        let domain = make_test_domain();
        let fem_nodes = vec![[0.3, 0.5], [0.5, 0.5], [0.7, 0.5]];
        let ci = setup_bridging_domain(&domain, &fem_nodes, 0.2, 0.8);
        assert!(!ci.nodes.is_empty());
    }
    #[test]
    fn test_generate_voronoi_cells() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 0.87], [0.5, 0.3]];
        let cells = generate_voronoi_cells(&pts, 3);
        assert!(!cells.is_empty());
    }
    #[test]
    fn test_node_density_map() {
        let domain = make_test_domain();
        let density = node_density_map(&domain, 3, 3, [0.0, 0.0, 1.0, 1.0]);
        assert_eq!(density.len(), 9);
        let max_d = density.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_d > 0.0);
    }
    #[test]
    fn test_error_indicators() {
        let domain = make_test_domain();
        let n = domain.node_count();
        let stresses: Vec<StressState> = (0..n)
            .map(|i| StressState {
                sigma_xx: 100.0 * i as f64,
                sigma_yy: 50.0,
                tau_xy: 0.0,
            })
            .collect();
        let indicators = compute_error_indicators(&domain, &stresses);
        assert_eq!(indicators.len(), n);
    }
    #[test]
    fn test_adaptive_refine() {
        let domain = make_test_domain();
        let n = domain.node_count();
        let indicators: Vec<ErrorIndicator> = (0..n)
            .map(|i| ErrorIndicator {
                node: i,
                error: if i == 0 { 100.0 } else { 0.0 },
            })
            .collect();
        let refined = adaptive_refine(&domain, &indicators, 50.0);
        assert!(refined.len() > n, "Refinement should add new nodes");
    }
}
/// Compute the convergence rate from two convergence points using
/// the formula: rate = ln(e1/e2) / ln(h1/h2).
#[allow(dead_code)]
pub fn compute_convergence_rate(
    coarse: &ConvergencePoint,
    fine: &ConvergencePoint,
) -> ConvergenceRate {
    let log_h_ratio = (coarse.h / fine.h).ln();
    let safe_rate = |e1: f64, e2: f64| -> f64 {
        if e1 <= 0.0 || e2 <= 0.0 || log_h_ratio.abs() < 1e-15 {
            0.0
        } else {
            (e1 / e2).ln() / log_h_ratio
        }
    };
    ConvergenceRate {
        l2_rate: safe_rate(coarse.l2_error, fine.l2_error),
        linf_rate: safe_rate(coarse.linf_error, fine.linf_error),
        h1_rate: safe_rate(coarse.h1_error, fine.h1_error),
        energy_rate: safe_rate(coarse.energy_error, fine.energy_error),
    }
}
/// Compute convergence rates from a series of refinements.
///
/// Returns a vector of rates (one fewer than the number of points).
#[allow(dead_code)]
pub fn compute_convergence_rates(points: &[ConvergencePoint]) -> Vec<ConvergenceRate> {
    if points.len() < 2 {
        return vec![];
    }
    points
        .windows(2)
        .map(|w| compute_convergence_rate(&w[0], &w[1]))
        .collect()
}
/// Richardson extrapolation to estimate the exact solution from two
/// approximations at different mesh sizes.
///
/// Given solutions u_h and u_{h/r} with refinement ratio r and
/// observed convergence order p, the extrapolated solution is:
///   u_exact ≈ (r^p * u_fine - u_coarse) / (r^p - 1)
#[allow(dead_code)]
pub fn richardson_extrapolation(
    u_coarse: f64,
    u_fine: f64,
    refinement_ratio: f64,
    order: f64,
) -> f64 {
    let rp = refinement_ratio.powf(order);
    (rp * u_fine - u_coarse) / (rp - 1.0)
}
/// Effectivity index: ratio of estimated error to true error.
///
/// A good error estimator should yield effectivity indices close to 1.0.
#[allow(dead_code)]
pub fn effectivity_index(estimated_error: f64, true_error: f64) -> f64 {
    if true_error.abs() < 1e-30 {
        0.0
    } else {
        estimated_error / true_error
    }
}
/// Compute the L2 error norm between a computed solution and an analytical
/// solution evaluated at the given points.
///
/// `computed` — vector of computed values at each node.
/// `analytical` — vector of exact values at each node.
/// `weights` — integration weights (e.g., Voronoi cell areas).
#[allow(dead_code)]
pub fn l2_error_norm(computed: &[f64], analytical: &[f64], weights: &[f64]) -> f64 {
    assert_eq!(computed.len(), analytical.len());
    assert_eq!(computed.len(), weights.len());
    let mut sum = 0.0;
    for i in 0..computed.len() {
        let diff = computed[i] - analytical[i];
        sum += weights[i] * diff * diff;
    }
    sum.sqrt()
}
/// Compute the L-infinity (max) error between computed and analytical solutions.
#[allow(dead_code)]
pub fn linf_error_norm(computed: &[f64], analytical: &[f64]) -> f64 {
    assert_eq!(computed.len(), analytical.len());
    computed
        .iter()
        .zip(analytical.iter())
        .map(|(c, a)| (c - a).abs())
        .fold(0.0_f64, f64::max)
}
/// Compute the H1 semi-norm error from gradient differences.
///
/// `grad_computed` — computed gradients \[du/dx, du/dy\] at each node.
/// `grad_analytical` — analytical gradients at each node.
/// `weights` — integration weights.
#[allow(dead_code)]
pub fn h1_seminorm_error(
    grad_computed: &[[f64; 2]],
    grad_analytical: &[[f64; 2]],
    weights: &[f64],
) -> f64 {
    assert_eq!(grad_computed.len(), grad_analytical.len());
    assert_eq!(grad_computed.len(), weights.len());
    let mut sum = 0.0;
    for i in 0..grad_computed.len() {
        let dx = grad_computed[i][0] - grad_analytical[i][0];
        let dy = grad_computed[i][1] - grad_analytical[i][1];
        sum += weights[i] * (dx * dx + dy * dy);
    }
    sum.sqrt()
}
/// Run a convergence study summary from a set of convergence points.
#[allow(dead_code)]
pub fn convergence_study_summary(points: Vec<ConvergencePoint>) -> ConvergenceStudySummary {
    let rates = compute_convergence_rates(&points);
    let avg_l2_rate = if rates.is_empty() {
        0.0
    } else {
        rates.iter().map(|r| r.l2_rate).sum::<f64>() / rates.len() as f64
    };
    let avg_h1_rate = if rates.is_empty() {
        0.0
    } else {
        rates.iter().map(|r| r.h1_rate).sum::<f64>() / rates.len() as f64
    };
    ConvergenceStudySummary {
        points,
        rates,
        avg_l2_rate,
        avg_h1_rate,
    }
}
/// Format a convergence table as a string for display.
#[allow(dead_code)]
pub fn format_convergence_table(summary: &ConvergenceStudySummary) -> String {
    let mut out = String::new();
    out.push_str("h           nodes   L2 error       Linf error     H1 error       Rate(L2)\n");
    out.push_str("----------  ------  -------------  -------------  -------------  --------\n");
    for (i, pt) in summary.points.iter().enumerate() {
        let rate_str = if i > 0 && i - 1 < summary.rates.len() {
            format!("{:.2}", summary.rates[i - 1].l2_rate)
        } else {
            "---".to_string()
        };
        out.push_str(&format!(
            "{:.6}  {:>6}  {:.6e}  {:.6e}  {:.6e}  {}\n",
            pt.h, pt.num_nodes, pt.l2_error, pt.linf_error, pt.h1_error, rate_str,
        ));
    }
    out.push_str(&format!("\nAverage L2 rate:  {:.3}\n", summary.avg_l2_rate));
    out.push_str(&format!("Average H1 rate:  {:.3}\n", summary.avg_h1_rate));
    out
}
#[cfg(test)]
mod tests_convergence {
    use super::*;
    use crate::meshfree_fem::*;
    fn make_convergence_points() -> Vec<ConvergencePoint> {
        vec![
            ConvergencePoint {
                h: 0.5,
                num_nodes: 9,
                l2_error: 0.25,
                linf_error: 0.4,
                h1_error: 0.5,
                energy_error: 0.35,
            },
            ConvergencePoint {
                h: 0.25,
                num_nodes: 25,
                l2_error: 0.0625,
                linf_error: 0.1,
                h1_error: 0.25,
                energy_error: 0.175,
            },
            ConvergencePoint {
                h: 0.125,
                num_nodes: 81,
                l2_error: 0.015625,
                linf_error: 0.025,
                h1_error: 0.125,
                energy_error: 0.0875,
            },
        ]
    }
    #[test]
    fn test_convergence_rate_l2_second_order() {
        let pts = make_convergence_points();
        let rate = compute_convergence_rate(&pts[0], &pts[1]);
        assert!(
            (rate.l2_rate - 2.0).abs() < 1e-10,
            "Expected L2 rate ~2.0, got {:.6}",
            rate.l2_rate
        );
    }
    #[test]
    fn test_convergence_rate_h1_first_order() {
        let pts = make_convergence_points();
        let rate = compute_convergence_rate(&pts[0], &pts[1]);
        assert!(
            (rate.h1_rate - 1.0).abs() < 1e-10,
            "Expected H1 rate ~1.0, got {:.6}",
            rate.h1_rate
        );
    }
    #[test]
    fn test_convergence_rates_series() {
        let pts = make_convergence_points();
        let rates = compute_convergence_rates(&pts);
        assert_eq!(rates.len(), 2);
        assert!((rates[0].l2_rate - 2.0).abs() < 1e-10);
        assert!((rates[1].l2_rate - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_convergence_rates_single_point() {
        let pts = vec![make_convergence_points()[0].clone()];
        let rates = compute_convergence_rates(&pts);
        assert!(rates.is_empty());
    }
    #[test]
    fn test_richardson_extrapolation() {
        let result = richardson_extrapolation(1.0, 1.25, 2.0, 2.0);
        let expected = (4.0 * 1.25 - 1.0) / 3.0;
        assert!((result - expected).abs() < 1e-12);
    }
    #[test]
    fn test_effectivity_index_perfect() {
        let idx = effectivity_index(0.05, 0.05);
        assert!((idx - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_effectivity_index_zero_true_error() {
        let idx = effectivity_index(0.01, 0.0);
        assert!((idx - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_l2_error_norm_zero_error() {
        let u = vec![1.0, 2.0, 3.0];
        let w = vec![1.0, 1.0, 1.0];
        let err = l2_error_norm(&u, &u, &w);
        assert!(err.abs() < 1e-15);
    }
    #[test]
    fn test_l2_error_norm_constant_error() {
        let computed = vec![1.0, 2.0, 3.0];
        let analytical = vec![1.1, 2.1, 3.1];
        let weights = vec![1.0, 1.0, 1.0];
        let err = l2_error_norm(&computed, &analytical, &weights);
        let expected = (3.0 * 0.01_f64).sqrt();
        assert!((err - expected).abs() < 1e-10);
    }
    #[test]
    fn test_linf_error_norm() {
        let computed = vec![1.0, 2.0, 3.0];
        let analytical = vec![1.1, 2.3, 3.05];
        let err = linf_error_norm(&computed, &analytical);
        assert!((err - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_h1_seminorm_error_zero() {
        let g = vec![[1.0, 2.0], [3.0, 4.0]];
        let w = vec![1.0, 1.0];
        let err = h1_seminorm_error(&g, &g, &w);
        assert!(err.abs() < 1e-15);
    }
    #[test]
    fn test_h1_seminorm_error_nonzero() {
        let gc = vec![[1.0, 0.0]];
        let ga = vec![[0.0, 0.0]];
        let w = vec![2.0];
        let err = h1_seminorm_error(&gc, &ga, &w);
        assert!((err - 2.0_f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_convergence_study_summary_avg_rates() {
        let pts = make_convergence_points();
        let summary = convergence_study_summary(pts);
        assert!((summary.avg_l2_rate - 2.0).abs() < 1e-10);
        assert!((summary.avg_h1_rate - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_format_convergence_table_not_empty() {
        let pts = make_convergence_points();
        let summary = convergence_study_summary(pts);
        let table = format_convergence_table(&summary);
        assert!(table.contains("L2 error"));
        assert!(table.contains("Average L2 rate"));
    }
    #[test]
    fn test_convergence_rate_linf() {
        let pts = make_convergence_points();
        let rate = compute_convergence_rate(&pts[0], &pts[1]);
        assert!(
            rate.linf_rate > 1.5,
            "Expected Linf rate > 1.5, got {:.6}",
            rate.linf_rate
        );
    }
    #[test]
    fn test_convergence_rate_energy() {
        let pts = make_convergence_points();
        let rate = compute_convergence_rate(&pts[0], &pts[1]);
        assert!(rate.energy_rate > 0.5);
    }
    #[test]
    fn test_richardson_extrapolation_first_order() {
        let result = richardson_extrapolation(1.0, 1.5, 2.0, 1.0);
        assert!((result - 2.0).abs() < 1e-12);
    }
}
