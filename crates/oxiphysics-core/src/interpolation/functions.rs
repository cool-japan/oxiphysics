//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BSplineBasis, RbfKernel};

/// Linear interpolation between two scalars.
///
/// Returns `a` at `t = 0` and `b` at `t = 1`.
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
/// Component-wise linear interpolation between two 3-vectors.
pub fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        lerp(a[0], b[0], t),
        lerp(a[1], b[1], t),
        lerp(a[2], b[2], t),
    ]
}
/// Bilinear interpolation on a unit square.
///
/// `f00` = value at (0,0), `f10` at (1,0), `f01` at (0,1), `f11` at (1,1).
/// `u` and `v` are parameters in \[0, 1\].
pub fn bilinear(f00: f64, f10: f64, f01: f64, f11: f64, u: f64, v: f64) -> f64 {
    lerp(lerp(f00, f10, u), lerp(f01, f11, u), v)
}
/// Trilinear interpolation on a unit cube.
///
/// `values` is ordered as `[f000, f100, f010, f110, f001, f101, f011, f111]`
/// where the index bits are (w-bit, v-bit, u-bit).
/// `u`, `v`, `w` are parameters in \[0, 1\].
pub fn trilinear(values: &[f64; 8], u: f64, v: f64, w: f64) -> f64 {
    let c00 = lerp(values[0], values[1], u);
    let c10 = lerp(values[2], values[3], u);
    let c01 = lerp(values[4], values[5], u);
    let c11 = lerp(values[6], values[7], u);
    let c0 = lerp(c00, c10, v);
    let c1 = lerp(c01, c11, v);
    lerp(c0, c1, w)
}
/// Dot product of two quaternions.
pub fn quat_dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}
/// Euclidean norm of a quaternion.
pub fn quat_norm(q: [f64; 4]) -> f64 {
    quat_dot(q, q).sqrt()
}
/// Normalize a quaternion to unit length.
pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let n = quat_norm(q);
    if n < f64::EPSILON {
        return [1.0, 0.0, 0.0, 0.0];
    }
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}
/// Normalized linear interpolation (nlerp) of two quaternions.
///
/// Fast approximation of slerp; constant angular velocity is not guaranteed.
pub fn quat_nlerp(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    let b = if quat_dot(a, b) < 0.0 {
        [-b[0], -b[1], -b[2], -b[3]]
    } else {
        b
    };
    let raw = [
        lerp(a[0], b[0], t),
        lerp(a[1], b[1], t),
        lerp(a[2], b[2], t),
        lerp(a[3], b[3], t),
    ];
    quat_normalize(raw)
}
/// Spherical linear interpolation (slerp) between two unit quaternions.
///
/// Falls back to `quat_nlerp` when the quaternions are nearly parallel
/// (`|dot| > 0.9995`).
pub fn quat_slerp(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    let mut dot = quat_dot(a, b);
    let b = if dot < 0.0 {
        dot = -dot;
        [-b[0], -b[1], -b[2], -b[3]]
    } else {
        b
    };
    if dot > 0.9995 {
        return quat_nlerp(a, b, t);
    }
    let theta_0 = dot.acos();
    let theta = theta_0 * t;
    let sin_theta = theta.sin();
    let sin_theta_0 = theta_0.sin();
    let s0 = (theta_0 - theta).sin() / sin_theta_0;
    let s1 = sin_theta / sin_theta_0;
    [
        s0 * a[0] + s1 * b[0],
        s0 * a[1] + s1 * b[1],
        s0 * a[2] + s1 * b[2],
        s0 * a[3] + s1 * b[3],
    ]
}
/// SQUAD — smooth quaternion spline interpolation.
///
/// Interpolates between `q1` and `q2` using auxiliary control points derived
/// from `q0` and `q3`. `t` ∈ \[0, 1\].
pub fn quat_squad(q0: [f64; 4], q1: [f64; 4], q2: [f64; 4], q3: [f64; 4], t: f64) -> [f64; 4] {
    let s1 = squad_inner(q0, q1, q2);
    let s2 = squad_inner(q1, q2, q3);
    let u = 2.0 * t * (1.0 - t);
    quat_slerp(quat_slerp(q1, q2, t), quat_slerp(s1, s2, t), u)
}
/// Compute the inner quadrangle point for SQUAD.
pub(super) fn squad_inner(q_prev: [f64; 4], q: [f64; 4], q_next: [f64; 4]) -> [f64; 4] {
    let inv_q_prev = quat_slerp(q, q_prev, 1.0);
    let inv_q_next = quat_slerp(q, q_next, 1.0);
    let mid = quat_nlerp(inv_q_prev, inv_q_next, 0.5);
    quat_slerp(q, mid, -0.25)
}
/// Cubic Hermite interpolation (scalar).
///
/// Blends position and tangent values: H00·p0 + H10·m0 + H01·p1 + H11·m1.
pub fn hermite(p0: f64, m0: f64, p1: f64, m1: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
}
/// Derivative of cubic Hermite at local parameter `t`.
pub fn hermite_deriv(p0: f64, m0: f64, p1: f64, m1: f64, t: f64) -> f64 {
    let t2 = t * t;
    let dh00 = 6.0 * t2 - 6.0 * t;
    let dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    let dh01 = -6.0 * t2 + 6.0 * t;
    let dh11 = 3.0 * t2 - 2.0 * t;
    dh00 * p0 + dh10 * m0 + dh01 * p1 + dh11 * m1
}
/// Cubic Hermite interpolation (3-vector).
pub fn hermite3(p0: [f64; 3], m0: [f64; 3], p1: [f64; 3], m1: [f64; 3], t: f64) -> [f64; 3] {
    [
        hermite(p0[0], m0[0], p1[0], m1[0], t),
        hermite(p0[1], m0[1], p1[1], m1[1], t),
        hermite(p0[2], m0[2], p1[2], m1[2], t),
    ]
}
/// Catmull-Rom segment between `p1` and `p2` using finite-difference tangents
/// derived from `p0` and `p3`.
///
/// At `t = 0` the result is `p1`; at `t = 1` it is `p2`.
pub fn catmull_rom(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], t: f64) -> [f64; 3] {
    let m1 = [
        (p2[0] - p0[0]) * 0.5,
        (p2[1] - p0[1]) * 0.5,
        (p2[2] - p0[2]) * 0.5,
    ];
    let m2 = [
        (p3[0] - p1[0]) * 0.5,
        (p3[1] - p1[1]) * 0.5,
        (p3[2] - p1[2]) * 0.5,
    ];
    hermite3(p1, m1, p2, m2, t)
}
/// Helper: distance between two 3-D points.
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Bilinear interpolation on a regular 2-D grid of values.
///
/// `grid[row][col]` holds the value at `(x_min + col*dx, y_min + row*dy)`.
/// Returns the linearly interpolated value at `(x, y)`.
/// Clamps `(x, y)` to the grid boundary.
pub fn bilinear_grid(
    grid: &[Vec<f64>],
    x_min: f64,
    y_min: f64,
    dx: f64,
    dy: f64,
    x: f64,
    y: f64,
) -> f64 {
    let rows = grid.len();
    if rows == 0 {
        return 0.0;
    }
    let cols = grid[0].len();
    if cols == 0 {
        return 0.0;
    }
    let fx = ((x - x_min) / dx).clamp(0.0, (cols - 1) as f64);
    let fy = ((y - y_min) / dy).clamp(0.0, (rows - 1) as f64);
    let c0 = fx.floor() as usize;
    let c1 = (c0 + 1).min(cols - 1);
    let r0 = fy.floor() as usize;
    let r1 = (r0 + 1).min(rows - 1);
    let u = fx - c0 as f64;
    let v = fy - r0 as f64;
    bilinear(grid[r0][c0], grid[r0][c1], grid[r1][c0], grid[r1][c1], u, v)
}
/// Parameters for a regular 3-D rectilinear grid.
#[derive(Debug, Clone, Copy)]
pub struct Grid3Params {
    /// Minimum x coordinate.
    pub x_min: f64,
    /// Minimum y coordinate.
    pub y_min: f64,
    /// Minimum z coordinate.
    pub z_min: f64,
    /// Grid spacing in x.
    pub dx: f64,
    /// Grid spacing in y.
    pub dy: f64,
    /// Grid spacing in z.
    pub dz: f64,
}

/// Trilinear interpolation on a regular 3-D grid.
///
/// `grid[z][y][x]` holds the value at grid index `(x, y, z)`.
/// Grid coordinates start at `(params.x_min, params.y_min, params.z_min)`
/// with steps `(params.dx, params.dy, params.dz)`.
/// Clamps to the grid boundary.
pub fn trilinear_grid(grid: &[Vec<Vec<f64>>], params: &Grid3Params, x: f64, y: f64, z: f64) -> f64 {
    let Grid3Params {
        x_min,
        y_min,
        z_min,
        dx,
        dy,
        dz,
    } = *params;
    let nz = grid.len();
    if nz == 0 {
        return 0.0;
    }
    let ny = grid[0].len();
    if ny == 0 {
        return 0.0;
    }
    let nx = grid[0][0].len();
    if nx == 0 {
        return 0.0;
    }
    let fi = ((x - x_min) / dx).clamp(0.0, (nx - 1) as f64);
    let fj = ((y - y_min) / dy).clamp(0.0, (ny - 1) as f64);
    let fk = ((z - z_min) / dz).clamp(0.0, (nz - 1) as f64);
    let x0 = fi.floor() as usize;
    let x1 = (x0 + 1).min(nx - 1);
    let y0 = fj.floor() as usize;
    let y1 = (y0 + 1).min(ny - 1);
    let z0 = fk.floor() as usize;
    let z1 = (z0 + 1).min(nz - 1);
    let u = fi - x0 as f64;
    let v = fj - y0 as f64;
    let w = fk - z0 as f64;
    let vals = [
        grid[z0][y0][x0],
        grid[z0][y0][x1],
        grid[z0][y1][x0],
        grid[z0][y1][x1],
        grid[z1][y0][x0],
        grid[z1][y0][x1],
        grid[z1][y1][x0],
        grid[z1][y1][x1],
    ];
    trilinear(&vals, u, v, w)
}
/// Cubic polynomial weight for the Catmull-Rom / Keys kernel.
pub(super) fn bicubic_weight(t: f64) -> f64 {
    let at = t.abs();
    if at < 1.0 {
        1.5 * at * at * at - 2.5 * at * at + 1.0
    } else if at < 2.0 {
        -0.5 * at * at * at + 2.5 * at * at - 4.0 * at + 2.0
    } else {
        0.0
    }
}
/// Bicubic interpolation using a 4×4 stencil on a regular 2-D grid.
///
/// `grid[row][col]` is indexed as `(y, x)`.  Grid spacing is assumed uniform
/// with spacing 1 (normalized coordinates); `x` and `y` are in `[0, cols-1]`
/// and `[0, rows-1]` respectively.  Out-of-bounds accesses are clamped.
pub fn bicubic(grid: &[Vec<f64>], x: f64, y: f64) -> f64 {
    let rows = grid.len() as isize;
    if rows == 0 {
        return 0.0;
    }
    let cols = grid[0].len() as isize;
    if cols == 0 {
        return 0.0;
    }
    let get = |r: isize, c: isize| -> f64 {
        let r = r.clamp(0, rows - 1) as usize;
        let c = c.clamp(0, cols - 1) as usize;
        grid[r][c]
    };
    let xi = x.floor() as isize;
    let yi = y.floor() as isize;
    let dx = x - xi as f64;
    let dy = y - yi as f64;
    let mut result = 0.0;
    for m in -1_isize..=2 {
        let wy = bicubic_weight(dy - m as f64);
        for n in -1_isize..=2 {
            let wx = bicubic_weight(dx - n as f64);
            result += wx * wy * get(yi + m, xi + n);
        }
    }
    result
}
/// Compute the barycentric coordinates `(u, v, w)` of point `p` with respect
/// to triangle `(a, b, c)` in 2-D.
///
/// Returns `(u, v, w)` such that `p = u·a + v·b + w·c` and `u + v + w = 1`.
/// If the triangle is degenerate (zero area) all coordinates are `1/3`.
pub fn barycentric_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> (f64, f64, f64) {
    let denom = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
    if denom.abs() < f64::EPSILON {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }
    let u = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / denom;
    let v = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / denom;
    let w = 1.0 - u - v;
    (u, v, w)
}
/// Interpolate a scalar at point `p` inside triangle `(a, b, c)` with values
/// `(fa, fb, fc)` using barycentric coordinates.
pub fn barycentric_interp_2d(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    fa: f64,
    fb: f64,
    fc: f64,
) -> f64 {
    let (u, v, w) = barycentric_2d(p, a, b, c);
    u * fa + v * fb + w * fc
}
/// Natural neighbor interpolation using a simplified area-based Sibson weight.
///
/// For each query point `p`, computes inverse-distance-squared weights to the
/// `k` nearest supplied data points (a practical stand-in for true Sibson
/// weights when no Voronoi library is available).
///
/// `points` — 2-D scatter positions.
/// `values` — scalar value at each point.
/// `p`      — query location.
/// `k`      — number of nearest neighbours to use (clamped to `n`).
pub fn natural_neighbor_interp(points: &[[f64; 2]], values: &[f64], p: [f64; 2], k: usize) -> f64 {
    let n = points.len();
    if n == 0 || values.len() != n {
        return 0.0;
    }
    let k = k.min(n).max(1);
    let mut dists: Vec<(usize, f64)> = points
        .iter()
        .enumerate()
        .map(|(i, q)| {
            let dx = p[0] - q[0];
            let dy = p[1] - q[1];
            (i, dx * dx + dy * dy)
        })
        .collect();
    dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    if dists[0].1 < f64::EPSILON * f64::EPSILON {
        return values[dists[0].0];
    }
    let mut weight_sum = 0.0;
    let mut weighted_val = 0.0;
    for &(idx, d2) in dists.iter().take(k) {
        let w = 1.0 / d2;
        weight_sum += w;
        weighted_val += w * values[idx];
    }
    if weight_sum < f64::EPSILON {
        return 0.0;
    }
    weighted_val / weight_sum
}
/// Evaluate a NURBS curve using the Cox-de Boor algorithm and return the
/// point together with its first derivative via automatic differentiation.
///
/// `t` is clamped to the valid knot span `[knots[degree\], knots[n - degree - 1 + 1]]`.
/// Returns `(point, tangent)` as two 3-D arrays.
pub fn nurbs_evaluate_with_derivative(
    control_points: &[[f64; 3]],
    weights: &[f64],
    degree: usize,
    knots: &[f64],
    t: f64,
) -> ([f64; 3], [f64; 3]) {
    let n = control_points.len();
    if n == 0 || weights.len() != n || knots.len() < n + degree + 1 {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut point = [0.0_f64; 3];
    let mut deriv = [0.0_f64; 3];
    let mut denom = 0.0_f64;
    let mut denom_deriv = 0.0_f64;
    for i in 0..n {
        let b = BSplineBasis::basis(i, degree, t, knots);
        let db = BSplineBasis::basis_derivative(i, degree, t, knots);
        let w = weights[i];
        let wb = w * b;
        let wdb = w * db;
        for k in 0..3 {
            point[k] += wb * control_points[i][k];
            deriv[k] += wdb * control_points[i][k];
        }
        denom += wb;
        denom_deriv += wdb;
    }
    if denom.abs() < f64::EPSILON {
        return ([0.0; 3], [0.0; 3]);
    }
    let inv_d = 1.0 / denom;
    let pt = [point[0] * inv_d, point[1] * inv_d, point[2] * inv_d];
    let tang = [
        (deriv[0] - denom_deriv * pt[0]) * inv_d,
        (deriv[1] - denom_deriv * pt[1]) * inv_d,
        (deriv[2] - denom_deriv * pt[2]) * inv_d,
    ];
    (pt, tang)
}
/// Thin-plate spline RBF: `r² · ln(r)` (convention: 0 at r=0).
pub fn rbf_thin_plate_spline(r: f64) -> f64 {
    if r < f64::EPSILON {
        return 0.0;
    }
    r * r * r.ln()
}
/// Evaluate a 2-D RBF interpolant using the thin-plate spline kernel.
///
/// `centers` — 2-D scatter positions.
/// `weights` — fitted weights (same length as `centers`).
/// `p`        — query point.
pub fn rbf_tps_evaluate(centers: &[[f64; 2]], weights: &[f64], p: [f64; 2]) -> f64 {
    centers
        .iter()
        .zip(weights.iter())
        .map(|(c, &w)| {
            let dx = p[0] - c[0];
            let dy = p[1] - c[1];
            let r = (dx * dx + dy * dy).sqrt();
            w * rbf_thin_plate_spline(r)
        })
        .sum()
}
/// Fit a 2-D RBF interpolant with the thin-plate spline kernel.
///
/// Returns fitted weights or an error if the system is singular.
pub fn rbf_tps_fit(centers: &[[f64; 2]], values: &[f64]) -> Result<Vec<f64>, String> {
    let n = centers.len();
    if n == 0 || values.len() != n {
        return Err("centers and values must be non-empty and same length".into());
    }
    let mut mat: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let dx = centers[i][0] - centers[j][0];
                    let dy = centers[i][1] - centers[j][1];
                    let r = (dx * dx + dy * dy).sqrt();
                    rbf_thin_plate_spline(r)
                })
                .collect()
        })
        .collect();
    let mut rhs = values.to_vec();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&a, &b| {
                mat[a][col]
                    .abs()
                    .partial_cmp(&mat[b][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .ok_or("RBF TPS matrix is singular")?;
        if mat[pivot][col].abs() < 1e-14 {
            return Err("RBF TPS matrix is singular or near-singular".into());
        }
        mat.swap(col, pivot);
        rhs.swap(col, pivot);
        let diag = mat[col][col];
        for cell in mat[col][col..].iter_mut() {
            *cell /= diag;
        }
        rhs[col] /= diag;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = mat[row][col];
            let col_vals: Vec<f64> = mat[col][col..].to_vec();
            for (cell, &cv) in mat[row][col..].iter_mut().zip(col_vals.iter()) {
                *cell -= cv * f;
            }
            let rv = rhs[col] * f;
            rhs[row] -= rv;
        }
    }
    Ok(rhs)
}
/// Monotone cubic (PCHIP / Fritsch-Carlson) interpolation.
///
/// Given data points `(xs[i], ys[i])` sorted by `xs`, computes the
/// Fritsch-Carlson derivative estimates at each knot, then evaluates the
/// piecewise Hermite cubic at query point `x`.
///
/// The Fritsch-Carlson method modifies the local finite-difference slopes so
/// that the interpolant is monotone within each interval where the data is
/// monotone.  This avoids the spurious oscillations that can occur with
/// unconstrained cubic splines.
///
/// # Panics
/// Panics if `xs` is empty or `xs.len() != ys.len()`.
pub fn monotone_cubic(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    assert!(
        n > 0 && n == ys.len(),
        "xs and ys must be non-empty and equal length"
    );
    if n == 1 {
        return ys[0];
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let delta: Vec<f64> = (0..n - 1)
        .map(|i| {
            let h = xs[i + 1] - xs[i];
            if h.abs() < f64::EPSILON {
                0.0
            } else {
                (ys[i + 1] - ys[i]) / h
            }
        })
        .collect();
    let mut d: Vec<f64> = vec![0.0; n];
    d[0] = delta[0];
    d[n - 1] = delta[n - 2];
    for i in 1..n - 1 {
        d[i] = (delta[i - 1] + delta[i]) / 2.0;
    }
    for i in 0..n - 1 {
        if delta[i].abs() < f64::EPSILON {
            d[i] = 0.0;
            d[i + 1] = 0.0;
            continue;
        }
        let alpha = d[i] / delta[i];
        let beta = d[i + 1] / delta[i];
        let tau = alpha * alpha + beta * beta;
        if tau > 9.0 {
            let scale = 3.0 / tau.sqrt();
            d[i] = scale * alpha * delta[i];
            d[i + 1] = scale * beta * delta[i];
        }
    }
    let seg = (0..n - 1).rev().find(|&i| x >= xs[i]).unwrap_or(0);
    let h = xs[seg + 1] - xs[seg];
    if h.abs() < f64::EPSILON {
        return ys[seg];
    }
    let t = (x - xs[seg]) / h;
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * ys[seg] + h10 * h * d[seg] + h01 * ys[seg + 1] + h11 * h * d[seg + 1]
}
/// Evaluate an RBF interpolant at query point `x`.
///
/// # Arguments
/// * `centers` - Slice of centre vectors (each of the same dimensionality).
/// * `weights` - Fitted RBF weights (one per centre).
/// * `x`       - Query point vector (same dimensionality as centres).
/// * `kernel`  - The RBF kernel to use.
///
/// # Returns
/// The interpolated value `sum_i weights[i] * phi(||x - centers[i]||)`.
pub fn rbf_interpolate(centers: &[Vec<f64>], weights: &[f64], x: &[f64], kernel: RbfKernel) -> f64 {
    centers
        .iter()
        .zip(weights.iter())
        .map(|(c, &w)| {
            let r = c
                .iter()
                .zip(x.iter())
                .map(|(ci, xi)| (ci - xi).powi(2))
                .sum::<f64>()
                .sqrt();
            w * kernel.eval(r)
        })
        .sum()
}
/// Fit RBF weights by solving the interpolation system Φ w = f.
///
/// Uses Gaussian elimination with partial pivoting.
/// Returns an error if the kernel matrix is singular.
pub fn rbf_fit(
    centers: &[Vec<f64>],
    values: &[f64],
    kernel: RbfKernel,
) -> Result<Vec<f64>, String> {
    let n = centers.len();
    if n == 0 || values.len() != n {
        return Err("centers and values must be non-empty and the same length".into());
    }
    let mut mat: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let r = centers[i]
                        .iter()
                        .zip(centers[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    kernel.eval(r)
                })
                .collect()
        })
        .collect();
    let mut rhs = values.to_vec();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&a, &b| {
                mat[a][col]
                    .abs()
                    .partial_cmp(&mat[b][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .ok_or("RBF matrix is empty")?;
        if mat[pivot][col].abs() < 1e-14 {
            return Err("RBF kernel matrix is singular or near-singular".into());
        }
        mat.swap(col, pivot);
        rhs.swap(col, pivot);
        let diag = mat[col][col];
        for cell in mat[col][col..].iter_mut() {
            *cell /= diag;
        }
        rhs[col] /= diag;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = mat[row][col];
            let col_vals: Vec<f64> = mat[col][col..].to_vec();
            for (cell, &cv) in mat[row][col..].iter_mut().zip(col_vals.iter()) {
                *cell -= cv * f;
            }
            let rv = rhs[col] * f;
            rhs[row] -= rv;
        }
    }
    Ok(rhs)
}
/// Floater-Hormann family of barycentric rational interpolants.
///
/// Computes the interpolant of blending parameter `d` at query point `x`.
///
/// The Floater-Hormann weights are:
/// ```text
/// w_i = (-1)^i * sum_{j=max(0,i-d)}^{min(i,n-1-d)} C(d, i-j)
/// ```
/// giving an interpolant of rational type with no poles on the real line
/// and convergence rate O(h^(d+1)).
///
/// Setting `d = 0` gives the Berrut interpolant; larger `d` gives higher
/// accuracy on smooth functions.
///
/// # Arguments
/// * `xs` - Distinct, sorted abscissae.
/// * `ys` - Function values at `xs`.
/// * `d`  - Blending parameter in `0..xs.len()`.
/// * `x`  - Query point.
pub fn barycentric_rational(xs: &[f64], ys: &[f64], d: usize, x: f64) -> f64 {
    let n = xs.len();
    assert!(
        n > 0 && n == ys.len(),
        "xs and ys must be non-empty and equal length"
    );
    assert!(d < n, "blending parameter d must be < n");
    for (i, &xi) in xs.iter().enumerate() {
        if (x - xi).abs() < f64::EPSILON * 100.0 {
            return ys[i];
        }
    }
    let mut weights = vec![0.0_f64; n];
    for (i, w) in weights.iter_mut().enumerate() {
        let j_lo = i.saturating_sub(d);
        let j_hi = (i + 1).min(n - d);
        let sign = if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 };
        let mut s = 0.0_f64;
        for j in j_lo..j_hi {
            let k = i - j;
            let mut c = 1.0_f64;
            for l in 0..k {
                c *= (d - l) as f64 / (l + 1) as f64;
            }
            s += c;
        }
        *w = sign * s;
    }
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for i in 0..n {
        let t = weights[i] / (x - xs[i]);
        num += t * ys[i];
        den += t;
    }
    num / den
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AkimaSpline;
    use crate::BSplineCurve;
    use crate::CatmullRomSpline;
    use crate::HermiteSpline;
    use crate::MonotoneCubicSpline;
    use crate::NaturalCubicSpline;
    use crate::NurbsCurve;
    use crate::RBFInterpolation;
    pub(super) const EPS: f64 = 1e-9;
    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }
    fn approx3(a: [f64; 3], b: [f64; 3]) -> bool {
        approx(a[0], b[0]) && approx(a[1], b[1]) && approx(a[2], b[2])
    }
    fn approx_eq(a: f64, b: f64) -> bool {
        approx(a, b)
    }
    fn approx_eq4(a: [f64; 4], b: [f64; 4]) -> bool {
        approx(a[0], b[0]) && approx(a[1], b[1]) && approx(a[2], b[2]) && approx(a[3], b[3])
    }
    #[test]
    fn test_lerp_at_zero_is_a() {
        assert!(approx(lerp(3.0, 7.0, 0.0), 3.0));
    }
    #[test]
    fn test_lerp_at_one_is_b() {
        assert!(approx(lerp(3.0, 7.0, 1.0), 7.0));
    }
    #[test]
    fn test_lerp_midpoint() {
        assert!(approx(lerp(0.0, 1.0, 0.5), 0.5));
    }
    #[test]
    fn test_bilinear_corners() {
        assert!(approx(bilinear(1.0, 2.0, 3.0, 4.0, 0.0, 0.0), 1.0));
        assert!(approx(bilinear(1.0, 2.0, 3.0, 4.0, 1.0, 0.0), 2.0));
        assert!(approx(bilinear(1.0, 2.0, 3.0, 4.0, 0.0, 1.0), 3.0));
        assert!(approx(bilinear(1.0, 2.0, 3.0, 4.0, 1.0, 1.0), 4.0));
    }
    #[test]
    fn test_trilinear_corners() {
        let vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        assert!(approx(trilinear(&vals, 0.0, 0.0, 0.0), 1.0));
        assert!(approx(trilinear(&vals, 1.0, 1.0, 1.0), 8.0));
    }
    #[test]
    fn test_slerp_t0() {
        let a = quat_normalize([1.0, 0.0, 0.0, 0.0]);
        let b = quat_normalize([0.0, 1.0, 0.0, 0.0]);
        assert!(approx_eq4(quat_slerp(a, b, 0.0), a));
    }
    #[test]
    fn test_slerp_t1() {
        let a = quat_normalize([1.0, 0.0, 0.0, 0.0]);
        let b = quat_normalize([0.0, 1.0, 0.0, 0.0]);
        assert!(approx_eq4(quat_slerp(a, b, 1.0), b));
    }
    #[test]
    fn test_hermite_spline_endpoints() {
        let pts = vec![1.0, 3.0, 2.0];
        let tans = vec![0.5, -0.5, 1.0];
        let sp = HermiteSpline::new(pts, tans).unwrap();
        assert!(approx(sp.evaluate(0.0), 1.0), "start = first point");
        assert!(approx(sp.evaluate(2.0), 2.0), "end = last point");
    }
    #[test]
    fn test_hermite_spline_passes_through_interior() {
        let pts = vec![0.0, 1.0, 0.0];
        let tans = vec![1.0, 0.0, -1.0];
        let sp = HermiteSpline::new(pts, tans).unwrap();
        assert!(approx(sp.evaluate(1.0), 1.0));
    }
    #[test]
    fn test_hermite_spline_error_mismatch() {
        let result = HermiteSpline::new(vec![1.0, 2.0], vec![0.0]);
        assert!(result.is_err());
    }
    #[test]
    fn test_catmull_rom_spline_passes_through_points() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 3.0, 0.0],
        ];
        let sp = CatmullRomSpline::new(pts.clone(), 0.5);
        for (i, &p) in pts.iter().enumerate() {
            let ev = sp.evaluate(i as f64);
            assert!(approx3(ev, p), "t={i}: expected {p:?}, got {ev:?}");
        }
    }
    #[test]
    fn test_catmull_rom_arc_length_positive() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let sp = CatmullRomSpline::new(pts, 0.5);
        let l = sp.arc_length(1000);
        assert!(l > 0.0, "arc length should be positive");
    }
    #[test]
    fn test_natural_cubic_spline_interpolates_at_knots() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 1.0, 0.0, 1.0, 0.0];
        let sp = NaturalCubicSpline::fit(xs.clone(), ys.clone()).unwrap();
        for (xi, yi) in xs.iter().zip(ys.iter()) {
            let v = sp.evaluate(*xi);
            assert!((v - yi).abs() < 1e-10, "at x={xi}: expected {yi}, got {v}");
        }
    }
    #[test]
    fn test_natural_cubic_spline_derivative_continuous() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 4.0, 9.0];
        let sp = NaturalCubicSpline::fit(xs, ys).unwrap();
        for &x in &[1.0_f64, 2.0_f64] {
            let dl = sp.derivative(x - 1e-7);
            let dr = sp.derivative(x + 1e-7);
            assert!(
                (dl - dr).abs() < 1e-5,
                "derivative jump at x={x}: left={dl}, right={dr}"
            );
        }
    }
    #[test]
    fn test_natural_cubic_spline_integral() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 2.0];
        let sp = NaturalCubicSpline::fit(xs, ys).unwrap();
        let ig = sp.integral(0.0, 2.0);
        assert!((ig - 2.0).abs() < 1e-10, "integral = {ig}");
    }
    #[test]
    fn test_bspline_degree1_piecewise_linear() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::new(pts.clone(), 1);
        let p0 = curve.evaluate(0.0);
        assert!(approx3(p0, pts[0]), "start: {p0:?}");
        let pn = curve.evaluate(1.0);
        assert!(approx3(pn, *pts.last().unwrap()), "end: {pn:?}");
    }
    #[test]
    fn test_bspline_endpoints_clamped() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = BSplineCurve::new(pts.clone(), 3);
        let p0 = curve.evaluate(0.0);
        let pn = curve.evaluate(1.0);
        assert!(approx3(p0, pts[0]), "clamped start: {p0:?}");
        assert!(approx3(pn, *pts.last().unwrap()), "clamped end: {pn:?}");
    }
    #[test]
    fn test_nurbs_unit_weights_equals_bspline() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let degree = 3;
        let weights = vec![1.0; pts.len()];
        let nurbs = NurbsCurve::new(pts.clone(), weights, degree);
        let bspline = BSplineCurve::new(pts, degree);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let pn = nurbs.evaluate(t);
            let pb = bspline.evaluate(t);
            assert!(approx3(pn, pb), "t={t}: nurbs={pn:?}, bspline={pb:?}");
        }
    }
    #[test]
    fn test_rbf_exact_at_data_points() {
        let points: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let values = vec![1.0, 2.0, 3.0, 4.0];
        let rbf = RBFInterpolation::fit(&points, &values, 1.0).unwrap();
        for (p, &v) in points.iter().zip(values.iter()) {
            let ev = rbf.evaluate(*p);
            assert!((ev - v).abs() < 1e-6, "at {p:?}: expected {v}, got {ev}");
        }
    }
    #[test]
    fn test_bilinear_grid_at_corners() {
        let grid = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        assert!(approx(
            bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0),
            1.0
        ));
        assert!(approx(
            bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0),
            2.0
        ));
        assert!(approx(
            bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0),
            3.0
        ));
        assert!(approx(
            bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
            4.0
        ));
    }
    #[test]
    fn test_bilinear_grid_midpoint() {
        let grid = vec![vec![0.0, 0.0], vec![0.0, 4.0]];
        let v = bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, 0.5, 0.5);
        assert!((v - 1.0).abs() < 1e-9, "midpoint should be 1.0, got {v}");
    }
    #[test]
    fn test_bilinear_grid_clamping() {
        let grid = vec![vec![5.0, 6.0], vec![7.0, 8.0]];
        let v = bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0);
        assert!(approx(v, 5.0), "clamped to (0,0) should be 5.0, got {v}");
    }
    #[test]
    fn test_trilinear_grid_corners() {
        let grid = vec![
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            vec![vec![5.0, 6.0], vec![7.0, 8.0]],
        ];
        let gp = Grid3Params {
            x_min: 0.0,
            y_min: 0.0,
            z_min: 0.0,
            dx: 1.0,
            dy: 1.0,
            dz: 1.0,
        };
        assert!(approx(trilinear_grid(&grid, &gp, 0.0, 0.0, 0.0), 1.0));
        assert!(approx(trilinear_grid(&grid, &gp, 1.0, 1.0, 1.0), 8.0));
    }
    #[test]
    fn test_bicubic_at_integer_coords() {
        let grid = vec![
            vec![1.0, 2.0, 3.0, 4.0],
            vec![5.0, 6.0, 7.0, 8.0],
            vec![9.0, 10.0, 11.0, 12.0],
            vec![13.0, 14.0, 15.0, 16.0],
        ];
        let v = bicubic(&grid, 1.0, 1.0);
        assert!(
            (v - 6.0).abs() < 1e-6,
            "bicubic at (1,1) should be ~6, got {v}"
        );
    }
    #[test]
    fn test_bicubic_interior_bounded() {
        let grid = vec![
            vec![0.0, 1.0, 0.0],
            vec![1.0, 4.0, 1.0],
            vec![0.0, 1.0, 0.0],
        ];
        let v = bicubic(&grid, 1.0, 1.0);
        assert!(v > 0.0, "bicubic center should be positive, got {v}");
    }
    #[test]
    fn test_barycentric_2d_centroid() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        let p = [1.0 / 3.0, 1.0 / 3.0];
        let (u, v, w) = barycentric_2d(p, a, b, c);
        assert!((u - 1.0 / 3.0).abs() < 1e-9);
        assert!((v - 1.0 / 3.0).abs() < 1e-9);
        assert!((w - 1.0 / 3.0).abs() < 1e-9);
        assert!((u + v + w - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_barycentric_2d_vertex_a() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        let (u, v, w) = barycentric_2d(a, a, b, c);
        assert!((u - 1.0).abs() < 1e-9);
        assert!(v.abs() < 1e-9);
        assert!(w.abs() < 1e-9);
    }
    #[test]
    fn test_barycentric_interp_2d_linear() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        let p = [0.5, 0.0];
        let v = barycentric_interp_2d(p, a, b, c, 0.0, 1.0, 2.0);
        assert!((v - 0.5).abs() < 1e-9, "expected 0.5, got {v}");
    }
    #[test]
    fn test_natural_neighbor_at_data_point() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let values = vec![3.0, 5.0, 7.0];
        let v = natural_neighbor_interp(&points, &values, [0.0, 0.0], 3);
        assert!(
            (v - 3.0).abs() < 1e-9,
            "exact data point should return its value, got {v}"
        );
    }
    #[test]
    fn test_natural_neighbor_midpoint() {
        let points = vec![[0.0, 0.0], [2.0, 0.0]];
        let values = vec![0.0, 4.0];
        let v = natural_neighbor_interp(&points, &values, [1.0, 0.0], 2);
        assert!((v - 2.0).abs() < 1e-9, "midpoint should be 2.0, got {v}");
    }
    #[test]
    fn test_monotone_cubic_spline_interpolates_knots() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 4.0, 9.0];
        let sp = MonotoneCubicSpline::fit(xs.clone(), ys.clone()).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            let v = sp.evaluate(x);
            assert!((v - y).abs() < 1e-9, "at x={x}: expected {y}, got {v}");
        }
    }
    #[test]
    fn test_monotone_cubic_spline_monotone() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let sp = MonotoneCubicSpline::fit(xs, ys).unwrap();
        let mut prev = sp.evaluate(0.0);
        for i in 1..=40 {
            let x = i as f64 * 0.1;
            let cur = sp.evaluate(x);
            assert!(
                cur >= prev - 1e-9,
                "monotone cubic should be non-decreasing at x={x}"
            );
            prev = cur;
        }
    }
    #[test]
    fn test_monotone_cubic_spline_clamped() {
        let xs = vec![0.0, 1.0];
        let ys = vec![2.0, 5.0];
        let sp = MonotoneCubicSpline::fit(xs, ys).unwrap();
        assert!(
            approx(sp.evaluate(-1.0), 2.0),
            "should clamp to first value"
        );
        assert!(approx(sp.evaluate(2.0), 5.0), "should clamp to last value");
    }
    #[test]
    fn test_akima_spline_interpolates_knots() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = vec![0.0, 1.0, 8.0, 27.0, 64.0, 125.0];
        let sp = AkimaSpline::fit(xs.clone(), ys.clone()).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            let v = sp.evaluate(x);
            assert!((v - y).abs() < 1e-6, "at x={x}: expected {y}, got {v}");
        }
    }
    #[test]
    fn test_akima_spline_clamped() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![1.0, 3.0, 2.0];
        let sp = AkimaSpline::fit(xs, ys).unwrap();
        assert!(approx(sp.evaluate(-5.0), 1.0));
        assert!(approx(sp.evaluate(10.0), 2.0));
    }
    #[test]
    fn test_akima_spline_error_on_unsorted() {
        assert!(AkimaSpline::fit(vec![1.0, 0.0], vec![0.0, 1.0]).is_err());
    }
    #[test]
    fn test_nurbs_with_derivative_endpoint() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let weights = vec![1.0; 4];
        let degree = 3;
        let knots = BSplineCurve::uniform_knots(4, degree);
        let (pt, _tang) = nurbs_evaluate_with_derivative(&pts, &weights, degree, &knots, 0.0);
        assert!(
            approx3(pt, pts[0]),
            "t=0 should be first control point, got {pt:?}"
        );
    }
    #[test]
    fn test_nurbs_with_derivative_at_one() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let weights = vec![1.0; 4];
        let degree = 3;
        let knots = BSplineCurve::uniform_knots(4, degree);
        let (pt, _tang) = nurbs_evaluate_with_derivative(&pts, &weights, degree, &knots, 1.0);
        assert!(
            approx3(pt, *pts.last().unwrap()),
            "t=1 should be last control point, got {pt:?}"
        );
    }
    #[test]
    fn test_rbf_tps_fit_and_evaluate_exact() {
        let centers = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let values = vec![0.0, 1.0, 2.0, 3.0];
        let weights = rbf_tps_fit(&centers, &values).unwrap();
        for (c, &v) in centers.iter().zip(values.iter()) {
            let ev = rbf_tps_evaluate(&centers, &weights, *c);
            assert!((ev - v).abs() < 1e-5, "at {c:?}: expected {v}, got {ev}");
        }
    }
    #[test]
    fn test_rbf_tps_kernel_zero_at_origin() {
        assert!(rbf_thin_plate_spline(0.0).abs() < 1e-15);
    }
    #[test]
    fn test_lerp_endpoints_legacy() {
        assert!(approx_eq(lerp(3.0, 7.0, 0.0), 3.0));
        assert!(approx_eq(lerp(3.0, 7.0, 1.0), 7.0));
    }
    #[test]
    fn test_hermite_endpoints_legacy() {
        assert!(approx_eq(hermite(1.0, 0.5, 3.0, 0.5, 0.0), 1.0));
        assert!(approx_eq(hermite(1.0, 0.5, 3.0, 0.5, 1.0), 3.0));
    }
    #[test]
    fn test_catmull_rom_endpoints_legacy() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [2.0, 1.0, 0.0];
        let p3 = [3.0, 0.0, 0.0];
        assert!(approx3(catmull_rom(p0, p1, p2, p3, 0.0), p1));
        assert!(approx3(catmull_rom(p0, p1, p2, p3, 1.0), p2));
    }
    #[test]
    fn test_lerp_extrapolation() {
        assert!(approx(lerp(0.0, 10.0, 2.0), 20.0));
        assert!(approx(lerp(0.0, 10.0, -1.0), -10.0));
    }
    #[test]
    fn test_lerp3_component_wise() {
        let a = [1.0, 2.0, 3.0];
        let b = [3.0, 4.0, 5.0];
        let mid = lerp3(a, b, 0.5);
        assert!(approx3(mid, [2.0, 3.0, 4.0]));
    }
    #[test]
    fn test_bilinear_centre() {
        let v = bilinear(1.0, 3.0, 5.0, 7.0, 0.5, 0.5);
        assert!(approx(v, (1.0 + 3.0 + 5.0 + 7.0) / 4.0));
    }
    #[test]
    fn test_trilinear_unit_cube_centre() {
        let vals = [0.0_f64; 8];
        assert!(approx(trilinear(&vals, 0.5, 0.5, 0.5), 0.0));
    }
    #[test]
    fn test_natural_cubic_spline_extrapolation_left() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 4.0];
        let sp = NaturalCubicSpline::fit(xs, ys).unwrap();
        let v = sp.evaluate(-1.0);
        assert!(v.is_finite(), "extrapolated value should be finite");
    }
    #[test]
    fn test_natural_cubic_spline_extrapolation_right() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 4.0];
        let sp = NaturalCubicSpline::fit(xs, ys).unwrap();
        let v = sp.evaluate(10.0);
        assert!(v.is_finite(), "extrapolated value should be finite");
    }
    #[test]
    fn test_natural_cubic_spline_monotone_input() {
        let xs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x * x).collect();
        let sp = NaturalCubicSpline::fit(xs.clone(), ys).unwrap();
        for &x in &[0.5_f64, 1.5, 2.5, 3.5] {
            let v = sp.evaluate(x);
            let expected = x * x;
            assert!(
                (v - expected).abs() < 0.2,
                "x={x}: spline={v}, parabola={expected}"
            );
        }
    }
    #[test]
    fn test_natural_cubic_spline_too_few_points_error() {
        assert!(NaturalCubicSpline::fit(vec![1.0], vec![1.0]).is_err());
    }
    #[test]
    fn test_nurbs_midpoint_inside_hull() {
        let pts = vec![[0.0, 0.0, 0.0], [2.0, 2.0, 0.0], [4.0, 0.0, 0.0]];
        let weights = vec![1.0; 3];
        let nurbs = NurbsCurve::new(pts, weights, 2);
        let mid = nurbs.evaluate(0.5);
        assert!(mid[0] > 0.0 && mid[0] < 4.0, "x out of range: {}", mid[0]);
        assert!(mid[1] >= 0.0, "y should be non-negative: {}", mid[1]);
    }
    #[test]
    fn test_nurbs_evaluate_multiple_interior_points() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0; 3];
        let nurbs = NurbsCurve::new(pts, weights, 2);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let p = nurbs.evaluate(t);
            assert!(
                p[0].is_finite() && p[1].is_finite(),
                "t={t}: non-finite point {p:?}"
            );
        }
    }
    #[test]
    fn test_bspline_evaluates_midpoint() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::new(pts, 2);
        let mid = curve.evaluate(0.5);
        assert!(mid[0] > 0.0 && mid[0] < 2.0, "mid[0]={}", mid[0]);
    }
    #[test]
    fn test_bspline_evaluates_many_t() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::new(pts, 2);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let p = curve.evaluate(t);
            assert!(
                p[0].is_finite() && p[1].is_finite(),
                "t={t}: non-finite {p:?}"
            );
        }
    }
    #[test]
    fn test_rbf_interpolation_constant() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let vals = vec![5.0, 5.0, 5.0];
        let rbf = RBFInterpolation::fit(&pts, &vals, 1.0).unwrap();
        for (p, &v) in pts.iter().zip(vals.iter()) {
            let ev = rbf.evaluate(*p);
            assert!(
                (ev - v).abs() < 1e-4,
                "RBF at data point: expected {v}, got {ev}"
            );
        }
    }
    #[test]
    fn test_pchip_constant_data() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![7.0, 7.0, 7.0, 7.0];
        let sp = MonotoneCubicSpline::fit(xs, ys).unwrap();
        for &x in &[0.3_f64, 1.0, 1.7, 2.5] {
            assert!(approx(sp.evaluate(x), 7.0));
        }
    }
    #[test]
    fn test_pchip_decreasing_data_is_monotone() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![4.0, 3.0, 2.0, 1.0];
        let sp = MonotoneCubicSpline::fit(xs, ys).unwrap();
        let mut prev = sp.evaluate(0.0);
        for i in 1..=30 {
            let x = i as f64 * 0.1;
            let cur = sp.evaluate(x);
            assert!(
                cur <= prev + 1e-9,
                "PCHIP should be non-increasing at x={x}"
            );
            prev = cur;
        }
    }
    #[test]
    fn test_akima_spline_linear_data() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * x + 1.0).collect();
        let sp = AkimaSpline::fit(xs, ys.clone()).unwrap();
        for &x in &[0.5_f64, 1.5, 2.5, 3.5, 4.5] {
            let v = sp.evaluate(x);
            let expected = 2.0 * x + 1.0;
            assert!(
                (v - expected).abs() < 1e-6,
                "Akima linear at x={x}: {v} != {expected}"
            );
        }
    }
    #[test]
    fn test_akima_interpolates_all_knots() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = vec![1.0, -1.0, 2.0, 0.5, 3.0, 1.5];
        let sp = AkimaSpline::fit(xs.clone(), ys.clone()).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            let v = sp.evaluate(x);
            assert!(
                (v - y).abs() < 1e-8,
                "Akima at x={x}: expected {y}, got {v}"
            );
        }
    }
    #[test]
    fn test_hermite_midpoint() {
        let v = hermite(0.0, 0.0, 2.0, 0.0, 0.5);
        assert!((v - 1.0).abs() < 0.1, "Hermite midpoint approx: {v}");
    }
    #[test]
    fn test_hermite3_endpoints() {
        let p0 = [0.0, 0.0, 0.0];
        let m0 = [1.0, 0.0, 0.0];
        let p1 = [1.0, 1.0, 0.0];
        let m1 = [0.0, 1.0, 0.0];
        let r0 = hermite3(p0, m0, p1, m1, 0.0);
        let r1 = hermite3(p0, m0, p1, m1, 1.0);
        assert!(approx3(r0, p0), "hermite3 at t=0: {r0:?}");
        assert!(approx3(r1, p1), "hermite3 at t=1: {r1:?}");
    }
    #[test]
    fn test_quat_normalize_unit() {
        let q = [3.0, 1.0, 2.0, 1.0];
        let n = quat_normalize(q);
        let norm = quat_norm(n);
        assert!((norm - 1.0).abs() < 1e-12, "normalized quat norm={norm}");
    }
    #[test]
    fn test_quat_nlerp_t_half() {
        let a = [1.0, 0.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0, 0.0];
        let mid = quat_nlerp(a, b, 0.5);
        let norm = quat_norm(mid);
        assert!(
            (norm - 1.0).abs() < 1e-10,
            "nlerp result should be unit: {norm}"
        );
    }
    #[test]
    fn test_slerp_unit_result() {
        let a = quat_normalize([1.0, 0.5, 0.0, 0.0]);
        let b = quat_normalize([0.0, 0.0, 1.0, 0.5]);
        for t_int in 0..=10 {
            let t = t_int as f64 / 10.0;
            let r = quat_slerp(a, b, t);
            let norm = quat_norm(r);
            assert!((norm - 1.0).abs() < 1e-10, "slerp at t={t}: norm={norm}");
        }
    }
    #[test]
    fn test_bilinear_grid_linear_x() {
        let grid = vec![vec![0.0, 1.0], vec![0.0, 1.0]];
        for i in 0..=10 {
            let x = i as f64 / 10.0;
            let v = bilinear_grid(&grid, 0.0, 0.0, 1.0, 1.0, x, 0.5);
            assert!(
                (v - x).abs() < 1e-12,
                "bilinear_grid linear x at x={x}: {v}"
            );
        }
    }
    #[test]
    fn test_trilinear_grid_z_linear() {
        let grid = vec![
            vec![vec![0.0, 0.0], vec![0.0, 0.0]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        ];
        let gp = Grid3Params {
            x_min: 0.0,
            y_min: 0.0,
            z_min: 0.0,
            dx: 1.0,
            dy: 1.0,
            dz: 1.0,
        };
        for i in 0..=10 {
            let z = i as f64 / 10.0;
            let v = trilinear_grid(&grid, &gp, 0.5, 0.5, z);
            assert!((v - z).abs() < 1e-12, "trilinear z-linear at z={z}: {v}");
        }
    }
    #[test]
    fn test_monotone_cubic_interpolates_knots() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 4.0, 9.0];
        for (xi, yi) in xs.iter().zip(ys.iter()) {
            let v = monotone_cubic(&xs, &ys, *xi);
            assert!(
                (v - yi).abs() < 1e-10,
                "pchip at x={xi}: got {v}, expected {yi}"
            );
        }
    }
    #[test]
    fn test_monotone_cubic_linear_data_exact() {
        let xs: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * x + 1.0).collect();
        for i in 0..=20 {
            let x = i as f64 * 0.25;
            let expected = 2.0 * x + 1.0;
            let v = monotone_cubic(&xs, &ys, x);
            assert!(
                (v - expected).abs() < 1e-9,
                "monotone_cubic linear at x={x}: got {v}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_monotone_cubic_monotone_on_monotone_data() {
        let xs = vec![0.0, 1.0, 2.0, 4.0, 8.0];
        let ys = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let mut prev = monotone_cubic(&xs, &ys, 0.0);
        for i in 1..=40 {
            let x = i as f64 * 0.2;
            if x > 8.0 {
                break;
            }
            let v = monotone_cubic(&xs, &ys, x);
            assert!(
                v >= prev - 1e-10,
                "monotone_cubic should be non-decreasing at x={x}: {prev} -> {v}"
            );
            prev = v;
        }
    }
    #[test]
    fn test_monotone_cubic_extrapolates_to_boundary() {
        let xs = vec![1.0, 2.0, 3.0];
        let ys = vec![10.0, 20.0, 30.0];
        assert!((monotone_cubic(&xs, &ys, 0.0) - 10.0).abs() < 1e-12);
        assert!((monotone_cubic(&xs, &ys, 5.0) - 30.0).abs() < 1e-12);
    }
    #[test]
    fn test_rbf_gaussian_interpolates_centres() {
        let centers = vec![vec![0.0_f64], vec![1.0], vec![2.0]];
        let values = vec![1.0, 2.0, 3.0];
        let kernel = RbfKernel::Gaussian(1.0);
        let weights = rbf_fit(&centers, &values, kernel).expect("RBF fit should succeed");
        for (c, &v) in centers.iter().zip(values.iter()) {
            let pred = rbf_interpolate(&centers, &weights, c, kernel);
            assert!(
                (pred - v).abs() < 1e-8,
                "RBF at centre {}: pred={pred}, expected={v}",
                c[0]
            );
        }
    }
    #[test]
    fn test_rbf_multiquadric_fit_1d() {
        let centers: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let values: Vec<f64> = centers.iter().map(|c| c[0] * c[0]).collect();
        let kernel = RbfKernel::Multiquadric(0.5);
        let weights = rbf_fit(&centers, &values, kernel).expect("fit should succeed");
        for (c, &v) in centers.iter().zip(values.iter()) {
            let pred = rbf_interpolate(&centers, &weights, c, kernel);
            assert!(
                (pred - v).abs() < 1e-6,
                "MQ RBF at {}: pred={pred}, expected={v}",
                c[0]
            );
        }
    }
    #[test]
    fn test_rbf_fit_empty_returns_error() {
        let centers: Vec<Vec<f64>> = vec![];
        let values: Vec<f64> = vec![];
        let result = rbf_fit(&centers, &values, RbfKernel::Gaussian(1.0));
        assert!(result.is_err(), "empty input should return error");
    }
    #[test]
    fn test_barycentric_rational_interpolates_nodes() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys: Vec<f64> = xs.iter().map(|&x| x * x - x + 1.0).collect();
        for (xi, yi) in xs.iter().zip(ys.iter()) {
            let v = barycentric_rational(&xs, &ys, 2, *xi);
            assert!(
                (v - yi).abs() < 1e-8,
                "barycentric_rational at x={xi}: got {v}, expected {yi}"
            );
        }
    }
    #[test]
    fn test_barycentric_rational_d0_is_berrut() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![1.0, 4.0, 2.0, 5.0];
        for (xi, yi) in xs.iter().zip(ys.iter()) {
            let v = barycentric_rational(&xs, &ys, 0, *xi);
            assert!(
                (v - yi).abs() < 1e-8,
                "Berrut (d=0) at x={xi}: got {v}, expected {yi}"
            );
        }
    }
    #[test]
    fn test_barycentric_rational_smooth_function() {
        let n = 8;
        let xs: Vec<f64> = (0..n)
            .map(|i| i as f64 * std::f64::consts::PI / (n - 1) as f64)
            .collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x.sin()).collect();
        for i in 0..n - 1 {
            let x = (xs[i] + xs[i + 1]) / 2.0;
            let exact = x.sin();
            let v = barycentric_rational(&xs, &ys, 3, x);
            assert!(
                (v - exact).abs() < 0.01,
                "barycentric_rational sin at x={x:.4}: got {v:.6}, exact={exact:.6}"
            );
        }
    }
}
