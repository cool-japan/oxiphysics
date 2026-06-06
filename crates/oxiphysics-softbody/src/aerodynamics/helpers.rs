//! Internal helper vector math functions for aerodynamics.

#[inline]
pub(crate) fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
pub(crate) fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub(crate) fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub(crate) fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
pub(crate) fn v3_normalize(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = v3_norm(a);
    if n < 1e-12 {
        None
    } else {
        Some(v3_scale(a, 1.0 / n))
    }
}

/// Biot-Savart induced velocity due to a finite vortex filament
/// from point `a` to point `b` with circulation `gamma` at field point `p`.
pub(crate) fn biot_savart_filament(a: [f64; 3], b: [f64; 3], p: [f64; 3], gamma: f64) -> [f64; 3] {
    let r1 = v3_sub(p, a);
    let r2 = v3_sub(p, b);
    let r0 = v3_sub(b, a);
    let cross = v3_cross(r1, r2);
    let cross_sq = v3_dot(cross, cross);
    let core_radius = 1e-6;
    if cross_sq < core_radius * core_radius {
        return [0.0; 3];
    }
    let r1_mag = v3_norm(r1);
    let r2_mag = v3_norm(r2);
    if r1_mag < 1e-12 || r2_mag < 1e-12 {
        return [0.0; 3];
    }
    let dot_factor = v3_dot(
        r0,
        v3_sub(v3_scale(r1, 1.0 / r1_mag), v3_scale(r2, 1.0 / r2_mag)),
    );
    let coeff = gamma / (4.0 * std::f64::consts::PI * cross_sq) * dot_factor;
    v3_scale(cross, coeff)
}

/// Alias for `biot_savart_filament` with clearer argument names.
pub(crate) fn biot_savart_velocity(
    r1: [f64; 3],
    r2: [f64; 3],
    gamma: f64,
    eval: [f64; 3],
) -> [f64; 3] {
    biot_savart_filament(r1, r2, eval, gamma)
}

/// Gaussian elimination with partial pivoting.
pub(crate) fn gauss_solve_opt(n: usize, a: &mut [f64], b: &[f64]) -> Option<Vec<f64>> {
    let mut x = b.to_vec();
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = a[col * n + col].abs();
        for row in (col + 1)..n {
            let val = a[row * n + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None;
        }
        if max_row != col {
            for k in 0..n {
                a.swap(col * n + k, max_row * n + k);
            }
            x.swap(col, max_row);
        }
        let pivot = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            for k in col..n {
                let val = a[col * n + k];
                a[row * n + k] -= factor * val;
            }
            x[row] -= factor * x[col];
        }
    }
    for row in (0..n).rev() {
        let mut sum = x[row];
        for k in (row + 1)..n {
            sum -= a[row * n + k] * x[k];
        }
        let diag = a[row * n + row];
        if diag.abs() < 1e-30 {
            return None;
        }
        x[row] = sum / diag;
    }
    Some(x)
}

/// Smooth sigmoid blending function centred at `centre` with half-width `w`.
///
/// Returns 0 when `x << centre` and 1 when `x >> centre`.
#[inline]
pub(crate) fn stall_blend(x: f64, centre: f64, w: f64) -> f64 {
    let t = (x - centre) / w.max(1e-10);
    1.0 / (1.0 + (-t).exp())
}

/// Linear interpolation lookup in a table of (xs, ys) pairs.
pub(crate) fn interp_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if xs.is_empty() || ys.is_empty() {
        return 0.0;
    }
    if x <= xs[0] {
        return ys[0];
    }
    let last = xs.len() - 1;
    if x >= xs[last] {
        return ys[last.min(ys.len() - 1)];
    }
    for i in 0..last {
        if x >= xs[i] && x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] + t * (ys[i + 1] - ys[i]);
        }
    }
    ys[last.min(ys.len() - 1)]
}
