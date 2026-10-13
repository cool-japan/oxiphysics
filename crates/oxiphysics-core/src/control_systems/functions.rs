//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use nalgebra::{DMatrix, DVector};

use super::response::{C64, SisoRealization, polynomial_from_roots, real_polynomial, sample_count};
use super::types::{Complex, PidController, StateSpaceModel, TransferFunction};

/// Evaluate a polynomial with coefficients `[a_n, a_{n-1}, ..., a_0]` at complex z.
pub fn poly_eval_complex(coeffs: &[f64], z: Complex) -> Complex {
    let mut result = Complex::new(0.0, 0.0);
    for &c in coeffs {
        result = result.mul(&z).add(&Complex::new(c, 0.0));
    }
    result
}
/// Evaluate polynomial at real x using Horner's method.
pub fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
    let mut r = 0.0;
    for &c in coeffs {
        r = r * x + c;
    }
    r
}
/// Multiply two polynomials (convolution of coefficient arrays).
pub fn poly_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let n = a.len() + b.len() - 1;
    let mut result = vec![0.0; n];
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            result[i + j] += ai * bj;
        }
    }
    result
}
/// Add two polynomials (may differ in length).
pub fn poly_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len().max(b.len());
    let mut result = vec![0.0; n];
    let offset_a = n - a.len();
    let offset_b = n - b.len();
    for (i, &v) in a.iter().enumerate() {
        result[i + offset_a] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        result[i + offset_b] += v;
    }
    result
}
/// Derivative of polynomial coefficients `[a_n, ..., a_0]`.
pub fn poly_deriv(coeffs: &[f64]) -> Vec<f64> {
    if coeffs.len() <= 1 {
        return vec![0.0];
    }
    let deg = coeffs.len() - 1;
    let mut d = Vec::with_capacity(deg);
    for (i, &c) in coeffs.iter().enumerate().take(deg) {
        d.push(c * (deg - i) as f64);
    }
    d
}
/// Find polynomial roots using Durand-Kerner method.
///
/// Returns approximate complex roots for polynomials of any degree.
pub fn poly_roots(coeffs: &[f64]) -> Vec<Complex> {
    let mut start = 0;
    while start < coeffs.len() && coeffs[start].abs() < 1e-300 {
        start += 1;
    }
    let coeffs = &coeffs[start..];
    if coeffs.len() <= 1 {
        return vec![];
    }
    let n = coeffs.len() - 1;
    if n == 0 {
        return vec![];
    }
    let lead = coeffs[0];
    if lead.abs() < 1e-300 {
        return vec![];
    }
    let norm: Vec<f64> = coeffs.iter().map(|&c| c / lead).collect();
    if n == 1 {
        return vec![Complex::new(-norm[1], 0.0)];
    }
    if n == 2 {
        let a = 1.0;
        let b = norm[1];
        let c = norm[2];
        let disc = b * b - 4.0 * a * c;
        if disc >= 0.0 {
            let sq = disc.sqrt();
            return vec![
                Complex::new((-b + sq) / (2.0 * a), 0.0),
                Complex::new((-b - sq) / (2.0 * a), 0.0),
            ];
        } else {
            let sq = (-disc).sqrt();
            return vec![
                Complex::new(-b / (2.0 * a), sq / (2.0 * a)),
                Complex::new(-b / (2.0 * a), -sq / (2.0 * a)),
            ];
        }
    }
    let mut roots: Vec<Complex> = (0..n)
        .map(|k| {
            let angle = 2.0 * PI * k as f64 / n as f64 + 0.3;
            Complex::from_polar(0.5 + 0.1 * k as f64, angle)
        })
        .collect();
    for _iter in 0..200 {
        let mut max_delta = 0.0_f64;
        for i in 0..n {
            let mut denom = Complex::new(1.0, 0.0);
            for j in 0..n {
                if j != i {
                    denom = denom.mul(&roots[i].sub(&roots[j]));
                }
            }
            let val = poly_eval_complex(&norm, roots[i]);
            if denom.norm_sq() < 1e-300 {
                continue;
            }
            let delta = val.div(&denom);
            roots[i] = roots[i].sub(&delta);
            max_delta = max_delta.max(delta.abs());
        }
        if max_delta < 1e-12 {
            break;
        }
    }
    roots
}
/// Routh-Hurwitz stability test.
///
/// Returns `true` if all roots of the polynomial have negative real parts.
/// The polynomial is given as `[a_n, a_{n-1}, ..., a_0]` (descending order).
pub fn routh_hurwitz_stable(coeffs: &[f64]) -> bool {
    if coeffs.is_empty() {
        return true;
    }
    let n = coeffs.len();
    if n == 1 {
        return coeffs[0] > 0.0;
    }
    let sign = coeffs[0].signum();
    for &c in coeffs {
        if c * sign <= 0.0 {
            return false;
        }
    }
    let rows = n;
    let cols = n.div_ceil(2);
    let mut table = vec![vec![0.0; cols]; rows];
    for j in 0..cols {
        if 2 * j < n {
            table[0][j] = coeffs[2 * j];
        }
        if 2 * j + 1 < n {
            table[1][j] = coeffs[2 * j + 1];
        }
    }
    for i in 2..rows {
        let pivot = table[i - 1][0];
        if pivot.abs() < 1e-15 {
            return false;
        }
        for j in 0..cols - 1 {
            let a = table[i - 2][j + 1];
            let b = if j + 1 < cols {
                table[i - 1][j + 1]
            } else {
                0.0
            };
            table[i][j] = (pivot * a - table[i - 2][0] * b) / pivot;
        }
    }
    let first_sign = table[0][0].signum();
    for row in &table {
        if row[0] * first_sign < 0.0 {
            return false;
        }
        if row[0].abs() < 1e-15 {
            return false;
        }
    }
    true
}
/// Advance state by one Euler step: x_{k+1} = x_k + dt*(A*x_k + B*u_k).
pub fn ss_step(model: &StateSpaceModel, state: &[f64], input: &[f64], dt: f64) -> Vec<f64> {
    let dx: Vec<f64> = model
        .a
        .iter()
        .zip(model.b.iter())
        .map(|(a_row, b_row)| {
            let ax: f64 = a_row
                .iter()
                .zip(state.iter())
                .map(|(aij, sj)| aij * sj)
                .sum();
            let bu: f64 = b_row
                .iter()
                .zip(input.iter())
                .map(|(bij, uj)| bij * uj)
                .sum();
            ax + bu
        })
        .collect();
    state
        .iter()
        .zip(dx.iter())
        .map(|(si, dxi)| si + dt * dxi)
        .collect()
}
/// Compute output: y = C*x + D*u.
pub fn ss_output(model: &StateSpaceModel, state: &[f64], input: &[f64]) -> Vec<f64> {
    model
        .c
        .iter()
        .zip(model.d.iter())
        .map(|(c_row, d_row)| {
            let cx: f64 = c_row
                .iter()
                .zip(state.iter())
                .map(|(cij, sj)| cij * sj)
                .sum();
            let du: f64 = d_row
                .iter()
                .zip(input.iter())
                .map(|(dij, uj)| dij * uj)
                .sum();
            cx + du
        })
        .collect()
}
/// Simulate state-space model over time steps, returning state history.
pub fn ss_simulate(
    model: &StateSpaceModel,
    x0: &[f64],
    inputs: &[Vec<f64>],
    dt: f64,
) -> Vec<Vec<f64>> {
    let mut history = Vec::with_capacity(inputs.len() + 1);
    let mut x = x0.to_vec();
    history.push(x.clone());
    for u in inputs {
        x = ss_step(model, &x, u, dt);
        history.push(x.clone());
    }
    history
}
/// Update PID controller and return the control output.
///
/// `setpoint` is the desired value, `measured` is the process variable,
/// `dt` is the time step.
pub fn pid_update(ctrl: &mut PidController, setpoint: f64, measured: f64, dt: f64) -> f64 {
    let error = setpoint - measured;
    let p_term = ctrl.kp * (ctrl.setpoint_weight_b * setpoint - measured);
    let anti_windup_correction = ctrl.anti_windup_gain
        * (ctrl.prev_unclipped - ctrl.prev_unclipped.clamp(ctrl.out_min, ctrl.out_max));
    ctrl.integral += error * dt - anti_windup_correction * dt;
    let i_term = ctrl.ki * ctrl.integral;
    let deriv_input = ctrl.setpoint_weight_c * setpoint - measured;
    let raw_deriv = if dt > 0.0 {
        (deriv_input - ctrl.prev_error) / dt
    } else {
        0.0
    };
    let alpha = ctrl.deriv_filter_coeff;
    let filtered_deriv = alpha * ctrl.prev_deriv_filtered + (1.0 - alpha) * raw_deriv;
    let d_term = ctrl.kd * filtered_deriv;
    ctrl.prev_error = deriv_input;
    ctrl.prev_deriv_filtered = filtered_deriv;
    let output_unclipped = p_term + i_term + d_term;
    ctrl.prev_unclipped = output_unclipped;
    output_unclipped.clamp(ctrl.out_min, ctrl.out_max)
}
/// Reset PID controller state.
pub fn pid_reset(ctrl: &mut PidController) {
    ctrl.integral = 0.0;
    ctrl.prev_error = 0.0;
    ctrl.prev_deriv_filtered = 0.0;
    ctrl.prev_unclipped = 0.0;
}
/// Generate logarithmically spaced frequencies.
pub fn logspace_omega(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![start];
    }
    let log_start = start.ln();
    let log_end = end.ln();
    (0..n)
        .map(|i| {
            let frac = i as f64 / (n - 1) as f64;
            (log_start + frac * (log_end - log_start)).exp()
        })
        .collect()
}
/// Bode plot data: returns (frequencies, magnitude_dB, phase_deg).
pub fn bode_data(
    tf: &TransferFunction,
    omega_start: f64,
    omega_end: f64,
    n_points: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let omegas = logspace_omega(omega_start, omega_end, n_points);
    let mut mag_db = Vec::with_capacity(n_points);
    let mut phase_deg = Vec::with_capacity(n_points);
    for &w in &omegas {
        let h = tf.eval_jw(w);
        let m = h.abs();
        mag_db.push(20.0 * m.max(1e-300).log10());
        phase_deg.push(h.arg() * 180.0 / PI);
    }
    (omegas, mag_db, phase_deg)
}
/// Compute gain margin and phase margin from Bode data.
///
/// Returns `(gain_margin_dB, phase_margin_deg)`.
pub fn stability_margins(
    tf: &TransferFunction,
    omega_start: f64,
    omega_end: f64,
    n_points: usize,
) -> (f64, f64) {
    let omegas = logspace_omega(omega_start, omega_end, n_points);
    let mut gain_margin_db = f64::INFINITY;
    let mut phase_margin_deg = f64::INFINITY;
    for &w in &omegas {
        let h = tf.eval_jw(w);
        let mag = h.abs();
        let phase = h.arg() * 180.0 / PI;
        if (mag - 1.0).abs() < 0.05 {
            let pm = 180.0 + phase;
            if pm.abs() < phase_margin_deg.abs() {
                phase_margin_deg = pm;
            }
        }
        if (phase + 180.0).abs() < 5.0 {
            let gm = -20.0 * mag.max(1e-300).log10();
            if gm.abs() < gain_margin_db.abs() {
                gain_margin_db = gm;
            }
        }
    }
    (gain_margin_db, phase_margin_deg)
}
/// Nyquist plot data: returns (real_parts, imag_parts).
pub fn nyquist_data(
    tf: &TransferFunction,
    omega_start: f64,
    omega_end: f64,
    n_points: usize,
) -> (Vec<f64>, Vec<f64>) {
    let omegas = logspace_omega(omega_start, omega_end, n_points);
    let mut re_parts = Vec::with_capacity(n_points);
    let mut im_parts = Vec::with_capacity(n_points);
    for &w in &omegas {
        let h = tf.eval_jw(w);
        re_parts.push(h.re);
        im_parts.push(h.im);
    }
    (re_parts, im_parts)
}
/// Check Nyquist stability: counts encirclements of -1+0j.
///
/// Returns the number of clockwise encirclements.
pub fn nyquist_encirclements(
    tf: &TransferFunction,
    omega_start: f64,
    omega_end: f64,
    n_points: usize,
) -> i32 {
    let omegas = logspace_omega(omega_start, omega_end, n_points);
    let mut total_angle = 0.0;
    let mut prev_angle = 0.0_f64;
    for (k, &w) in omegas.iter().enumerate() {
        let h = tf.eval_jw(w);
        let shifted = Complex::new(h.re + 1.0, h.im);
        let angle = shifted.arg();
        if k > 0 {
            let mut d = angle - prev_angle;
            if d > PI {
                d -= 2.0 * PI;
            }
            if d < -PI {
                d += 2.0 * PI;
            }
            total_angle += d;
        }
        prev_angle = angle;
    }
    (total_angle / (2.0 * PI)).round() as i32
}
/// Compute root locus: poles of 1 + K*H(s) for varying gain K.
///
/// Returns a vector of (gain, poles) pairs.
pub fn root_locus(tf: &TransferFunction, gains: &[f64]) -> Vec<(f64, Vec<Complex>)> {
    let mut results = Vec::with_capacity(gains.len());
    for &k in gains {
        let scaled_num: Vec<f64> = tf.num.iter().map(|&c| c * k).collect();
        let cl_poly = poly_add(&tf.den, &scaled_num);
        let poles = poly_roots(&cl_poly);
        results.push((k, poles));
    }
    results
}
/// Find the break-away/break-in points on the real axis for root locus.
///
/// These are points where dK/ds = 0 on the real axis.
pub fn root_locus_breakpoints(tf: &TransferFunction) -> Vec<f64> {
    let num_d = poly_deriv(&tf.num);
    let den_d = poly_deriv(&tf.den);
    let p1 = poly_mul(&num_d, &tf.den);
    let p2 = poly_mul(&tf.num, &den_d);
    let neg_p2: Vec<f64> = p2.iter().map(|&c| -c).collect();
    let result_poly = poly_add(&p1, &neg_p2);
    let roots = poly_roots(&result_poly);
    roots
        .iter()
        .filter(|r| r.im.abs() < 1e-6)
        .map(|r| r.re)
        .collect()
}
/// Discretize a state-space model using Zero-Order Hold (ZOH).
///
/// Approximation: Ad = I + A*T + (A*T)^2/2, Bd = (I*T + A*T^2/2)*B
/// where T is the sampling period.
pub fn discretize_zoh(model: &StateSpaceModel, dt: f64) -> StateSpaceModel {
    let n = model.state_dim();
    // a2[i][j] = (A^2)[i][j]
    let a2: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    model.a[i]
                        .iter()
                        .zip(model.a.iter())
                        .map(|(&aik, ak_row)| aik * ak_row[j])
                        .sum()
                })
                .collect()
        })
        .collect();
    let ad: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let base = model.a[i][j] * dt + a2[i][j] * dt * dt / 2.0;
                    if i == j { base + 1.0 } else { base }
                })
                .collect()
        })
        .collect();
    let bd: Vec<Vec<f64>> = model
        .a
        .iter()
        .zip(model.b.iter())
        .map(|(a_row, b_row)| {
            b_row
                .iter()
                .enumerate()
                .map(|(j, &bij)| {
                    let a_b_term: f64 = a_row
                        .iter()
                        .zip(model.b.iter())
                        .map(|(&aik, bk_row)| aik * bk_row[j])
                        .sum();
                    dt * bij + a_b_term * dt * dt / 2.0
                })
                .collect()
        })
        .collect();
    StateSpaceModel::new(ad, bd, model.c.clone(), model.d.clone())
}
/// Discretize a transfer function using Tustin (bilinear) transform.
///
/// Substitutes s = (2/T)*(z-1)/(z+1) into H(s).
/// Returns discrete-time transfer function H(z).
pub fn discretize_tustin(tf: &TransferFunction, dt: f64) -> TransferFunction {
    let n_num = tf.num.len();
    let n_den = tf.den.len();
    let order = (n_num.max(n_den)).saturating_sub(1);
    let c = 2.0 / dt;
    let z_minus_1 = vec![1.0, -1.0];
    let z_plus_1 = vec![1.0, 1.0];
    let mut pow_zm1 = vec![vec![1.0]; order + 1];
    let mut pow_zp1 = vec![vec![1.0]; order + 1];
    for k in 1..=order {
        pow_zm1[k] = poly_mul(&pow_zm1[k - 1], &z_minus_1);
        pow_zp1[k] = poly_mul(&pow_zp1[k - 1], &z_plus_1);
    }
    let mut num_z = vec![0.0; order + 1];
    for (i, &coeff) in tf.num.iter().rev().enumerate() {
        let ck = c.powi(i as i32);
        let term = poly_mul(&pow_zm1[i], &pow_zp1[order - i]);
        let expanded = term.iter().map(|&t| t * coeff * ck).collect::<Vec<_>>();
        if expanded.len() > num_z.len() {
            num_z.resize(expanded.len(), 0.0);
        }
        let offset = num_z.len() - expanded.len();
        for (j, &v) in expanded.iter().enumerate() {
            num_z[j + offset] += v;
        }
    }
    let mut den_z = vec![0.0; order + 1];
    for (i, &coeff) in tf.den.iter().rev().enumerate() {
        let ck = c.powi(i as i32);
        let term = poly_mul(&pow_zm1[i], &pow_zp1[order - i]);
        let expanded = term.iter().map(|&t| t * coeff * ck).collect::<Vec<_>>();
        if expanded.len() > den_z.len() {
            den_z.resize(expanded.len(), 0.0);
        }
        let offset = den_z.len() - expanded.len();
        for (j, &v) in expanded.iter().enumerate() {
            den_z[j + offset] += v;
        }
    }
    let lead = den_z[0];
    if lead.abs() > 1e-300 {
        for v in &mut num_z {
            *v /= lead;
        }
        for v in &mut den_z {
            *v /= lead;
        }
    }
    TransferFunction::new(num_z, den_z)
}
/// Multiply n-by-k matrix A by k-by-m matrix B, returning n-by-m.
pub fn mat_mul(a: &[f64], b: &[f64], n: usize, k: usize, m: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * m];
    for i in 0..n {
        for j in 0..m {
            for l in 0..k {
                c[i * m + j] += a[i * k + l] * b[l * m + j];
            }
        }
    }
    c
}
/// Transpose n-by-m matrix to m-by-n.
pub fn mat_transpose(a: &[f64], n: usize, m: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * m];
    for i in 0..n {
        for j in 0..m {
            t[j * n + i] = a[i * m + j];
        }
    }
    t
}
/// Identity matrix of size n.
pub fn mat_eye(n: usize) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}
/// Add two n-by-m matrices.
pub fn mat_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}
/// Subtract two n-by-m matrices: a - b.
pub fn mat_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x - y).collect()
}
/// Scale matrix by scalar.
pub fn mat_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|&x| x * s).collect()
}
/// Invert a small matrix (up to 3x3) using cofactors.
pub fn mat_inv_small(a: &[f64], n: usize) -> Vec<f64> {
    if n == 1 {
        if a[0].abs() < 1e-300 {
            return vec![0.0];
        }
        return vec![1.0 / a[0]];
    }
    if n == 2 {
        let det = a[0] * a[3] - a[1] * a[2];
        if det.abs() < 1e-300 {
            return vec![0.0; 4];
        }
        return vec![a[3] / det, -a[1] / det, -a[2] / det, a[0] / det];
    }
    if n == 3 {
        let det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
            + a[2] * (a[3] * a[7] - a[4] * a[6]);
        if det.abs() < 1e-300 {
            return vec![0.0; 9];
        }
        let inv_det = 1.0 / det;
        return vec![
            (a[4] * a[8] - a[5] * a[7]) * inv_det,
            (a[2] * a[7] - a[1] * a[8]) * inv_det,
            (a[1] * a[5] - a[2] * a[4]) * inv_det,
            (a[5] * a[6] - a[3] * a[8]) * inv_det,
            (a[0] * a[8] - a[2] * a[6]) * inv_det,
            (a[2] * a[3] - a[0] * a[5]) * inv_det,
            (a[3] * a[7] - a[4] * a[6]) * inv_det,
            (a[1] * a[6] - a[0] * a[7]) * inv_det,
            (a[0] * a[4] - a[1] * a[3]) * inv_det,
        ];
    }
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = a[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col * 2 * n + col].abs();
        for row in (col + 1)..n {
            let v = aug[row * 2 * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-300 {
            return vec![0.0; n * n];
        }
        if max_row != col {
            for j in 0..(2 * n) {
                aug.swap(col * 2 * n + j, max_row * 2 * n + j);
            }
        }
        let pivot = aug[col * 2 * n + col];
        for j in 0..(2 * n) {
            aug[col * 2 * n + j] /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row * 2 * n + col];
                for j in 0..(2 * n) {
                    let v = aug[col * 2 * n + j];
                    aug[row * 2 * n + j] -= factor * v;
                }
            }
        }
    }
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * 2 * n + n + j];
        }
    }
    inv
}
/// Compute the controllability matrix \[B, AB, A^2B, ...\].
///
/// `a` is n-by-n, `b` is n-by-m (flat row-major). Returns n-by-(n*m).
pub fn controllability_matrix(a: &[f64], b: &[f64], n: usize, m: usize) -> Vec<f64> {
    let mut result = Vec::with_capacity(n * n * m);
    let mut ab = b.to_vec();
    for _k in 0..n {
        result.extend_from_slice(&ab);
        ab = mat_mul(a, &ab, n, n, m);
    }
    let cols = n * m;
    let mut matrix = vec![0.0; n * cols];
    for k in 0..n {
        for i in 0..n {
            for j in 0..m {
                matrix[i * cols + k * m + j] = result[k * n * m + i * m + j];
            }
        }
    }
    matrix
}
/// Compute the observability matrix \[C; CA; CA^2; ...\].
///
/// `a` is n-by-n, `c` is p-by-n (flat row-major). Returns (n*p)-by-n.
pub fn observability_matrix(a: &[f64], c: &[f64], n: usize, p: usize) -> Vec<f64> {
    let rows = n * p;
    let mut result = vec![0.0; rows * n];
    let mut ca = c.to_vec();
    for k in 0..n {
        for i in 0..p {
            for j in 0..n {
                result[(k * p + i) * n + j] = ca[i * n + j];
            }
        }
        ca = mat_mul(&ca, a, p, n, n);
    }
    result
}
/// Rank of an `n`-by-`m` row-major matrix by Gaussian elimination.
///
/// The columns are processed left to right while a *separate* pivot-row
/// counter tracks the rows already used: a column without a usable pivot is
/// skipped without consuming a row, and columns beyond `min(n, m)` are still
/// examined.  Rank-deficient and wide matrices (such as controllability
/// matrices `[B, AB, ...]` whose first columns vanish) are therefore ranked
/// correctly.  A pivot is accepted when its magnitude is at least
/// `1e-10 * max(1, max |a_ij|)`.
///
/// Entries missing from a slice shorter than `n * m` are treated as zero and
/// surplus entries are ignored.
pub fn matrix_rank(mat: &[f64], n: usize, m: usize) -> usize {
    let Some(size) = n.checked_mul(m) else {
        return 0;
    };
    if size == 0 {
        return 0;
    }
    let mut work: Vec<f64> = mat
        .iter()
        .copied()
        .chain(std::iter::repeat(0.0))
        .take(size)
        .collect();
    let scale = work.iter().fold(1.0_f64, |acc, v| acc.max(v.abs()));
    let tol = 1e-10 * scale;
    let mut pivot_row = 0;
    for col in 0..m {
        if pivot_row == n {
            break;
        }
        let mut max_val = 0.0_f64;
        let mut max_row = pivot_row;
        for row in pivot_row..n {
            let v = work[row * m + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < tol {
            continue;
        }
        if max_row != pivot_row {
            for j in 0..m {
                work.swap(pivot_row * m + j, max_row * m + j);
            }
        }
        let pivot = work[pivot_row * m + col];
        for row in (pivot_row + 1)..n {
            let factor = work[row * m + col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for j in col..m {
                let v = work[pivot_row * m + j];
                work[row * m + j] -= factor * v;
            }
        }
        pivot_row += 1;
    }
    pivot_row
}
/// Check if the system is controllable (rank of controllability matrix == n).
pub fn is_controllable(a: &[f64], b: &[f64], n: usize, m: usize) -> bool {
    let cm = controllability_matrix(a, b, n, m);
    matrix_rank(&cm, n, n * m) == n
}
/// Check if the system is observable (rank of observability matrix == n).
pub fn is_observable(a: &[f64], c: &[f64], n: usize, p: usize) -> bool {
    let om = observability_matrix(a, c, n, p);
    matrix_rank(&om, n * p, n) == n
}
/// Compute the LQR cost: J = x'Qx + u'Ru.
pub fn lqr_cost(x: &[f64], u: &[f64], q: &[Vec<f64>], r: &[Vec<f64>]) -> f64 {
    let n = x.len();
    let m = u.len();
    let mut cost = 0.0;
    for i in 0..n {
        for j in 0..n {
            cost += x[i] * q[i][j] * x[j];
        }
    }
    for i in 0..m {
        for j in 0..m {
            cost += u[i] * r[i][j] * u[j];
        }
    }
    cost
}
/// Design a lead compensator for a desired phase margin boost at crossover frequency.
///
/// Returns `(zero_freq, pole_freq)` for the compensator C(s) = (s + z) / (s + p).
pub fn lead_compensator(crossover_freq: f64, phase_margin_deg: f64) -> (f64, f64) {
    let phi_max = phase_margin_deg * PI / 180.0;
    let sin_phi = phi_max.sin();
    let alpha = ((1.0 - sin_phi) / (1.0 + sin_phi)).max(1e-10);
    let wm = crossover_freq;
    let zero = wm * alpha.sqrt();
    let pole = wm / alpha.sqrt();
    (zero, pole)
}
/// Design a lag compensator for steady-state error reduction.
///
/// Returns `(zero_freq, pole_freq)`.
pub fn lag_compensator(crossover_freq: f64, dc_gain_boost: f64) -> (f64, f64) {
    let beta = dc_gain_boost.max(1.0);
    let wm = crossover_freq;
    let zero = wm / 10.0;
    let pole = zero / beta;
    (zero, pole)
}
/// Observer gain via Ackermann's formula for a single-output system.
///
/// Returns the gain `L` (length `n`) for which the eigenvalues of `A - L C`
/// are exactly `desired_poles`, using the dual Ackermann formula
///
/// ```text
/// L = phi(A) * O^-1 * e_n,
/// O = [C; C A; ...; C A^(n-1)],   phi(s) = prod_i (s - p_i),
/// ```
///
/// where `e_n = [0, ..., 0, 1]^T`.  `a` is `n x n`, `c` is the single output
/// row (`1 x n`) and `desired_poles` holds one `[re, im]` entry per pole, so a
/// complex pole must be listed together with its conjugate.
///
/// Returns `None` when
/// * `c` does not have exactly one row of length `n`, `a` is not square, or
///   `desired_poles.len() != n`,
/// * any input is not finite,
/// * the pole set is not closed under complex conjugation (`phi` would have
///   complex coefficients, so no real gain exists), or
/// * the pair `(A, C)` is not observable (`O` is singular).
///
/// For `n == 0` the empty gain `Some(vec![])` is returned when `c` and
/// `desired_poles` are empty as well.
pub fn observer_gain_ackermann(
    a: &[Vec<f64>],
    c: &[Vec<f64>],
    desired_poles: &[[f64; 2]],
) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return (c.is_empty() && desired_poles.is_empty()).then(Vec::new);
    }
    if a.iter().any(|row| row.len() != n)
        || c.len() != 1
        || c[0].len() != n
        || desired_poles.len() != n
    {
        return None;
    }
    if a.iter()
        .flatten()
        .chain(c[0].iter())
        .chain(desired_poles.iter().flatten())
        .any(|v| !v.is_finite())
    {
        return None;
    }

    // Desired characteristic polynomial phi(s) = prod (s - p_i), real coefficients.
    let roots: Vec<C64> = desired_poles.iter().map(|p| C64::new(p[0], p[1])).collect();
    let phi = polynomial_from_roots(C64::new(1.0, 0.0), &roots);
    let magnitude = phi.iter().map(|c| c.re.abs()).fold(1.0_f64, f64::max);
    if phi.iter().any(|c| c.im.abs() > 1e-9 * magnitude) {
        return None;
    }

    // phi(A) by Horner's scheme: ((I A + phi_1 I) A + phi_2 I) ... + phi_n I.
    let a_mat = DMatrix::from_fn(n, n, |i, j| a[i][j]);
    let identity = DMatrix::<f64>::identity(n, n);
    let mut phi_a = identity.clone();
    for coefficient in phi.iter().skip(1) {
        phi_a = &phi_a * &a_mat + &identity * coefficient.re;
    }

    // v = O^-1 e_n through an LU solve (never an explicit inverse).
    let a_flat: Vec<f64> = a.iter().flatten().copied().collect();
    let observability = observability_matrix(&a_flat, &c[0], n, 1);
    if matrix_rank(&observability, n, n) < n {
        return None;
    }
    let mut e_n = DVector::<f64>::zeros(n);
    e_n[n - 1] = 1.0;
    let v = DMatrix::from_row_slice(n, n, &observability)
        .lu()
        .solve(&e_n)?;

    let gain: Vec<f64> = (phi_a * v).iter().copied().collect();
    gain.iter().all(|g| g.is_finite()).then_some(gain)
}
/// Check if all poles have negative real part (stable).
pub fn pole_placement_check(eigs: &[[f64; 2]]) -> bool {
    eigs.iter().all(|e| e[0] < 0.0)
}
/// Approximate the H-infinity norm (peak of |H(jw)|) over given frequencies.
pub fn h_infinity_norm_bound(tf: &TransferFunction, omegas: &[f64]) -> f64 {
    let mut max_mag = 0.0_f64;
    for &w in omegas {
        let h = tf.eval_jw(w);
        max_mag = max_mag.max(h.abs());
    }
    max_mag
}
/// Ziegler-Nichols PID tuning from ultimate gain and period.
///
/// Returns `(Kp, Ki, Kd)`.
pub fn ziegler_nichols_pid(ku: f64, tu: f64) -> (f64, f64, f64) {
    let kp = 0.6 * ku;
    let ki = 2.0 * kp / tu;
    let kd = kp * tu / 8.0;
    (kp, ki, kd)
}
/// Cohen-Coon PID tuning from plant parameters.
///
/// `k` is the static gain, `tau` is the time constant, `theta` is the dead time.
/// Returns `(Kp, Ki, Kd)`.
pub fn cohen_coon_pid(k: f64, tau: f64, theta: f64) -> (f64, f64, f64) {
    let r = theta / tau;
    let kp = (1.0 / k) * (tau / theta) * (1.33 + r / 4.0);
    let ti = theta * (32.0 + 6.0 * r) / (13.0 + 8.0 * r);
    let td = theta * 4.0 / (11.0 + 2.0 * r);
    let ki = kp / ti;
    let kd = kp * td;
    (kp, ki, kd)
}
/// Check Lyapunov stability for a 2x2 system by examining eigenvalues of A.
///
/// Returns `true` if all eigenvalues have negative real part.
pub fn lyapunov_stable_2x2(a: &[f64; 4]) -> bool {
    let tr = a[0] + a[3];
    let det = a[0] * a[3] - a[1] * a[2];
    let disc = tr * tr - 4.0 * det;
    if disc >= 0.0 {
        let sq = disc.sqrt();
        let l1 = (tr + sq) / 2.0;
        let l2 = (tr - sq) / 2.0;
        l1 < 0.0 && l2 < 0.0
    } else {
        tr < 0.0
    }
}
/// Solve continuous Lyapunov equation A*P + P*A' + Q = 0 for 1x1 system.
///
/// Returns P\[0,0\].
pub fn solve_lyapunov_1x1(a: f64, q: f64) -> f64 {
    if a.abs() < 1e-300 {
        return 0.0;
    }
    -q / (2.0 * a)
}
/// Step response of a transfer function (via state-space simulation).
///
/// The transfer function is realised in controllable canonical form and
/// discretised with the exact zero-order hold, so the returned samples equal
/// the true unit-step response `y(k * dt)` (up to rounding) for any pole
/// pattern, including repeated, complex, unstable and integrating poles.
///
/// Returns `(time_vec, output_vec)` with `ceil(t_final / dt)` samples at
/// `t_k = k * dt`.  Both vectors are empty when no sampled response exists:
/// non-positive or non-finite `dt` / `t_final`, a transfer function that is
/// improper (numerator degree above the denominator degree), has a zero
/// denominator or non-finite coefficients.
pub fn step_response(tf: &TransferFunction, t_final: f64, dt: f64) -> (Vec<f64>, Vec<f64>) {
    sampled_response(tf, t_final, dt, SisoRealization::step_samples)
}
/// Impulse response of a transfer function (via state-space simulation).
///
/// Uses the same exact realisation as [`step_response`]: the unit impulse
/// moves the state to `B` and the system then evolves freely, so the samples
/// equal the true impulse response `g(k * dt)` up to rounding.  A direct
/// feedthrough term `D * delta(t)` of a biproper transfer function is a Dirac
/// impulse that cannot be represented by samples and is omitted.
///
/// Returns `(time_vec, output_vec)`; both are empty in the same situations as
/// for [`step_response`].
pub fn impulse_response(tf: &TransferFunction, t_final: f64, dt: f64) -> (Vec<f64>, Vec<f64>) {
    sampled_response(tf, t_final, dt, SisoRealization::impulse_samples)
}
/// Shared driver of [`step_response`] and [`impulse_response`].
fn sampled_response(
    tf: &TransferFunction,
    t_final: f64,
    dt: f64,
    sampler: fn(&SisoRealization, f64, usize) -> Option<Vec<f64>>,
) -> (Vec<f64>, Vec<f64>) {
    let Some(count) = sample_count(t_final, dt) else {
        return (Vec::new(), Vec::new());
    };
    let Some(realization) =
        SisoRealization::new(&real_polynomial(&tf.num), &real_polynomial(&tf.den))
    else {
        return (Vec::new(), Vec::new());
    };
    match sampler(&realization, dt, count) {
        Some(outputs) => ((0..count).map(|k| k as f64 * dt).collect(), outputs),
        None => (Vec::new(), Vec::new()),
    }
}
/// Compute the H2 norm of a stable SISO system via numerical integration.
///
/// H2 = sqrt( (1/2pi) * integral |H(jw)|^2 dw )
pub fn h2_norm(tf: &TransferFunction, omega_max: f64, n_points: usize) -> f64 {
    let dw = omega_max / n_points as f64;
    let mut integral = 0.0;
    for k in 0..n_points {
        let w = (k as f64 + 0.5) * dw;
        let h = tf.eval_jw(w);
        integral += h.norm_sq() * dw;
    }
    (integral / PI).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::control::KalmanFilter;
    use crate::control::LqrController;
    use crate::control_systems::MpcController;
    #[test]
    fn test_complex_abs() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.abs() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_complex_arg() {
        let z = Complex::new(1.0, 1.0);
        assert!((z.arg() - PI / 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_complex_mul() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a.mul(&b);
        assert!((c.re - (-5.0)).abs() < 1e-10);
        assert!((c.im - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_complex_div() {
        let a = Complex::new(1.0, 0.0);
        let b = Complex::new(0.0, 1.0);
        let c = a.div(&b);
        assert!((c.re).abs() < 1e-10);
        assert!((c.im - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_complex_from_polar() {
        let z = Complex::from_polar(2.0, PI / 2.0);
        assert!((z.re).abs() < 1e-10);
        assert!((z.im - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_poly_eval_constant() {
        assert!((poly_eval(&[5.0], 10.0) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_poly_eval_linear() {
        assert!((poly_eval(&[2.0, 1.0], 3.0) - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_poly_mul_degree() {
        let a = vec![1.0, 1.0];
        let b = vec![1.0, -1.0];
        let c = poly_mul(&a, &b);
        assert_eq!(c.len(), 3);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1]).abs() < 1e-10);
        assert!((c[2] + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_poly_add_different_lengths() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0];
        let c = poly_add(&a, &b);
        assert_eq!(c.len(), 3);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 6.0).abs() < 1e-10);
        assert!((c[2] - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_poly_deriv() {
        let d = poly_deriv(&[3.0, 2.0, 1.0]);
        assert_eq!(d.len(), 2);
        assert!((d[0] - 6.0).abs() < 1e-10);
        assert!((d[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_dc_gain_first_order() {
        let tf = TransferFunction::first_order(3.0, 2.0);
        assert!((tf.dc_gain() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_dc_gain_second_order() {
        let tf = TransferFunction::second_order(5.0, 0.7);
        assert!((tf.dc_gain() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_is_stable_first_order() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        assert!(tf.is_stable());
    }
    #[test]
    fn test_tf_is_unstable() {
        let tf = TransferFunction::new(vec![1.0, 1.0], vec![1.0, -1.0]);
        assert!(!tf.is_stable());
    }
    #[test]
    fn test_tf_series() {
        let h1 = TransferFunction::first_order(2.0, 1.0);
        let h2 = TransferFunction::first_order(3.0, 1.0);
        let hs = TransferFunction::series(&h1, &h2);
        assert!((hs.dc_gain() - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_poles_first_order() {
        let tf = TransferFunction::first_order(1.0, 2.0);
        let poles = tf.poles();
        assert_eq!(poles.len(), 1);
        assert!((poles[0].re - (-0.5)).abs() < 1e-6);
    }
    #[test]
    fn test_tf_poles_second_order_underdamped() {
        let tf = TransferFunction::second_order(10.0, 0.3);
        let poles = tf.poles();
        assert_eq!(poles.len(), 2);
        for p in &poles {
            assert!(p.re < 0.0, "poles should have negative real part");
        }
    }
    #[test]
    fn test_routh_hurwitz_stable_first_order() {
        assert!(routh_hurwitz_stable(&[1.0, 2.0]));
    }
    #[test]
    fn test_routh_hurwitz_unstable_negative_coeff() {
        assert!(!routh_hurwitz_stable(&[1.0, -1.0]));
    }
    #[test]
    fn test_routh_hurwitz_stable_second_order() {
        assert!(routh_hurwitz_stable(&[1.0, 3.0, 2.0]));
    }
    #[test]
    fn test_routh_hurwitz_stable_third_order() {
        assert!(routh_hurwitz_stable(&[1.0, 6.0, 11.0, 6.0]));
    }
    #[test]
    fn test_routh_hurwitz_unstable_third_order() {
        assert!(!routh_hurwitz_stable(&[1.0, 1.0, 1.0, 10.0]));
    }
    #[test]
    fn test_ss_step_integrator() {
        let model = StateSpaceModel::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        );
        let new_state = ss_step(&model, &[0.0], &[2.0], 0.5);
        assert!((new_state[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ss_output_passthrough() {
        let model = StateSpaceModel::new(
            vec![vec![0.0]],
            vec![vec![0.0]],
            vec![vec![0.0]],
            vec![vec![2.0]],
        );
        let y = ss_output(&model, &[0.0], &[3.0]);
        assert!((y[0] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_ss_step_2x2() {
        let model = StateSpaceModel::new(
            vec![vec![0.0, 1.0], vec![-1.0, 0.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        );
        let new_state = ss_step(&model, &[1.0, 0.0], &[0.0], 0.01);
        assert_eq!(new_state.len(), 2);
    }
    #[test]
    fn test_ss_simulate_length() {
        let model = StateSpaceModel::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        );
        let inputs = vec![vec![1.0]; 10];
        let history = ss_simulate(&model, &[0.0], &inputs, 0.1);
        assert_eq!(history.len(), 11);
    }
    #[test]
    fn test_ss_dimensions() {
        let model = StateSpaceModel::new(
            vec![vec![0.0, 1.0], vec![-1.0, 0.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            vec![vec![0.0], vec![0.0]],
        );
        assert_eq!(model.state_dim(), 2);
        assert_eq!(model.input_dim(), 1);
        assert_eq!(model.output_dim(), 2);
    }
    #[test]
    fn test_pid_proportional_only() {
        let mut ctrl = PidController::new(2.0, 0.0, 0.0);
        let out = pid_update(&mut ctrl, 5.0, 3.0, 0.1);
        assert!((out - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_pid_clamped() {
        let mut ctrl = PidController::new(100.0, 0.0, 0.0).with_limits(-1.0, 1.0);
        let out = pid_update(&mut ctrl, 10.0, 0.0, 0.1);
        assert!((out - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pid_reset_clears() {
        let mut ctrl = PidController::new(1.0, 1.0, 0.0);
        pid_update(&mut ctrl, 1.0, 0.0, 0.5);
        pid_reset(&mut ctrl);
        assert!(ctrl.integral.abs() < 1e-12);
        assert!(ctrl.prev_error.abs() < 1e-12);
    }
    #[test]
    fn test_pid_setpoint_weighting() {
        let mut ctrl = PidController::new(2.0, 0.0, 0.0).with_setpoint_weights(0.5, 0.0);
        let out = pid_update(&mut ctrl, 10.0, 0.0, 0.1);
        assert!((out - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_pid_derivative_filter() {
        let mut ctrl = PidController::new(0.0, 0.0, 1.0).with_deriv_filter(0.5);
        let _out1 = pid_update(&mut ctrl, 0.0, 0.0, 0.1);
        let _out2 = pid_update(&mut ctrl, 1.0, 0.0, 0.1);
        assert!(ctrl.prev_deriv_filtered.is_finite());
    }
    #[test]
    fn test_pid_anti_windup() {
        let mut ctrl = PidController::new(1.0, 10.0, 0.0)
            .with_limits(-1.0, 1.0)
            .with_anti_windup(1.0);
        for _ in 0..100 {
            pid_update(&mut ctrl, 100.0, 0.0, 0.1);
        }
        assert!(
            ctrl.integral.abs() < 1e6,
            "anti-windup should limit integral growth"
        );
    }
    #[test]
    fn test_logspace_length() {
        let omegas = logspace_omega(0.1, 1000.0, 50);
        assert_eq!(omegas.len(), 50);
    }
    #[test]
    fn test_logspace_endpoints() {
        let omegas = logspace_omega(1.0, 100.0, 10);
        assert!((omegas[0] - 1.0).abs() < 1e-10);
        assert!((omegas[9] - 100.0).abs() < 1e-6);
    }
    #[test]
    fn test_bode_data_length() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (f, m, p) = bode_data(&tf, 0.01, 100.0, 50);
        assert_eq!(f.len(), 50);
        assert_eq!(m.len(), 50);
        assert_eq!(p.len(), 50);
    }
    #[test]
    fn test_nyquist_data_length() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (re, im) = nyquist_data(&tf, 0.01, 100.0, 50);
        assert_eq!(re.len(), 50);
        assert_eq!(im.len(), 50);
    }
    #[test]
    fn test_root_locus_returns_correct_count() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let gains: Vec<f64> = (0..10).map(|i| i as f64 * 0.5).collect();
        let rl = root_locus(&tf, &gains);
        assert_eq!(rl.len(), 10);
    }
    #[test]
    fn test_root_locus_zero_gain_gives_open_loop_poles() {
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 2.0]);
        let rl = root_locus(&tf, &[0.0]);
        assert_eq!(rl.len(), 1);
        let poles = &rl[0].1;
        assert_eq!(poles.len(), 1);
        assert!((poles[0].re - (-2.0)).abs() < 1e-6);
    }
    #[test]
    fn test_root_locus_breakpoints() {
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 4.0, 3.0]);
        let bp = root_locus_breakpoints(&tf);
        assert!(!bp.is_empty(), "should find breakpoint");
        let has_near_minus2 = bp.iter().any(|&x| (x + 2.0).abs() < 0.5);
        assert!(has_near_minus2, "breakpoint near -2 expected, got {:?}", bp);
    }
    #[test]
    fn test_discretize_zoh_dimensions() {
        let model = StateSpaceModel::new(
            vec![vec![0.0, 1.0], vec![-2.0, -3.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        );
        let disc = discretize_zoh(&model, 0.01);
        assert_eq!(disc.a.len(), 2);
        assert_eq!(disc.a[0].len(), 2);
        assert_eq!(disc.b.len(), 2);
        assert_eq!(disc.b[0].len(), 1);
    }
    #[test]
    fn test_discretize_zoh_identity_for_zero_a() {
        let model = StateSpaceModel::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        );
        let disc = discretize_zoh(&model, 0.1);
        assert!((disc.a[0][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_discretize_tustin_preserves_dc() {
        let tf = TransferFunction::first_order(2.0, 1.0);
        let disc = discretize_tustin(&tf, 0.01);
        let h_z1 = poly_eval(&disc.num, 1.0) / poly_eval(&disc.den, 1.0);
        assert!(
            (h_z1 - 2.0).abs() < 0.01,
            "Tustin should preserve DC gain, got {:.6}",
            h_z1
        );
    }
    #[test]
    fn test_controllable_system() {
        let a = vec![0.0, 1.0, -2.0, -3.0];
        let b = vec![0.0, 1.0];
        assert!(is_controllable(&a, &b, 2, 1));
    }
    #[test]
    fn test_uncontrollable_system() {
        let a = vec![1.0, 0.0, 0.0, 2.0];
        let b = vec![1.0, 0.0];
        assert!(!is_controllable(&a, &b, 2, 1));
    }
    #[test]
    fn test_observable_system() {
        let a = vec![0.0, 1.0, -2.0, -3.0];
        let c = vec![1.0, 0.0];
        assert!(is_observable(&a, &c, 2, 1));
    }
    #[test]
    fn test_matrix_rank() {
        let m = vec![1.0, 0.0, 0.0, 1.0];
        assert_eq!(matrix_rank(&m, 2, 2), 2);
    }
    #[test]
    fn test_matrix_rank_deficient() {
        let m = vec![1.0, 2.0, 2.0, 4.0];
        assert_eq!(matrix_rank(&m, 2, 2), 1);
    }
    #[test]
    fn test_lqr_zero_state_zero_control() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let b = vec![vec![1.0], vec![0.0]];
        let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let r = vec![vec![1.0]];
        if let Some(lqr) = LqrController::solve(&a, &b, q, r) {
            let u = lqr.compute(&[0.0, 0.0]);
            assert_eq!(u[0], 0.0);
        }
    }
    #[test]
    fn test_lqr_cost_zero() {
        let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let r = vec![vec![1.0]];
        let cost = lqr_cost(&[0.0, 0.0], &[0.0], &q, &r);
        assert!(cost.abs() < 1e-12);
    }
    #[test]
    fn test_lqr_cost_non_negative() {
        let q = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let r = vec![vec![1.0]];
        let cost = lqr_cost(&[0.5, -0.3], &[0.1], &q, &r);
        assert!(cost >= 0.0);
    }
    #[test]
    fn test_kalman_initial_state() {
        let x0 = vec![0.0];
        let p0 = vec![vec![1.0]];
        let kf = KalmanFilter::new(x0, p0);
        assert_eq!(kf.x_hat[0], 0.0);
    }
    #[test]
    fn test_kalman_prediction_preserves_state() {
        let x0 = vec![5.0];
        let p0 = vec![vec![1.0]];
        let mut kf = KalmanFilter::new(x0, p0);
        // A = I, B = [0], u = [0], Q = [0.01]
        let a = vec![vec![1.0]];
        let b = vec![vec![0.0]];
        let q = vec![vec![0.01]];
        kf.predict(&a, &b, &[0.0], &q);
        assert!((kf.x_hat[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_kalman_update_moves_toward_measurement() {
        let x0 = vec![0.0];
        let p0 = vec![vec![1.0]];
        let mut kf = KalmanFilter::new(x0, p0);
        let a = vec![vec![1.0]];
        let b = vec![vec![0.0]];
        let q = vec![vec![0.01]];
        kf.predict(&a, &b, &[0.0], &q);
        let h = vec![vec![1.0]];
        let r = vec![vec![0.1]];
        let _ = kf.update(&h, &r, &[10.0]);
        assert!(kf.x_hat[0] > 0.0);
    }
    #[test]
    fn test_kalman_predict_then_update() {
        let x0 = vec![0.0];
        let p0 = vec![vec![1.0]];
        let mut kf = KalmanFilter::new(x0, p0);
        let a = vec![vec![1.0]];
        let b = vec![vec![0.0]];
        let q = vec![vec![0.01]];
        let h = vec![vec![1.0]];
        let r = vec![vec![0.1]];
        kf.predict(&a, &b, &[0.0], &q);
        let _ = kf.update(&h, &r, &[5.0]);
        assert_eq!(kf.x_hat.len(), 1);
        assert!(kf.x_hat[0] > 0.0);
    }
    #[test]
    fn test_lead_compensator_zero_lt_pole() {
        let (z, p) = lead_compensator(10.0, 45.0);
        assert!(z < p);
    }
    #[test]
    fn test_lead_compensator_positive() {
        let (z, p) = lead_compensator(5.0, 30.0);
        assert!(z > 0.0 && p > 0.0);
    }
    #[test]
    fn test_lag_compensator_pole_lt_zero() {
        let (z, p) = lag_compensator(10.0, 5.0);
        assert!(p < z);
    }
    #[test]
    fn test_observer_gain_1d() {
        let a = vec![vec![2.0]];
        let c = vec![vec![1.0]];
        let desired = [[-1.0, 0.0]];
        let l = observer_gain_ackermann(&a, &c, &desired).expect("observable 1-D system");
        assert_eq!(l.len(), 1);
        assert!((l[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_observer_gain_empty() {
        let l = observer_gain_ackermann(&[], &[], &[]);
        assert_eq!(l, Some(Vec::new()));
    }
    /// Coefficients `[1, c_{n-1}, ..., c_0]` of `det(sI - M)` (Faddeev-LeVerrier).
    fn characteristic_polynomial(m: &DMatrix<f64>) -> Vec<f64> {
        let n = m.nrows();
        let identity = DMatrix::<f64>::identity(n, n);
        let mut coefficients = vec![1.0];
        let mut mk = DMatrix::<f64>::zeros(n, n);
        for k in 1..=n {
            mk = m * &mk + &identity * coefficients[k - 1];
            coefficients.push(-(m * &mk).trace() / k as f64);
        }
        coefficients
    }
    fn closed_loop_matrix(a: &[Vec<f64>], c: &[Vec<f64>], l: &[f64]) -> DMatrix<f64> {
        let n = a.len();
        DMatrix::from_fn(n, n, |i, j| a[i][j] - l[i] * c[0][j])
    }
    #[test]
    fn test_observer_gain_double_integrator_places_poles() {
        // A = [[0, 1], [0, 0]], C = [1, 0], poles {-1, -2}: the old diagonal
        // heuristic gave L = [1, 2] (eig(A - LC) = -0.5 +- 1.32j); the correct
        // gain is L = [3, 2].
        let a = vec![vec![0.0, 1.0], vec![0.0, 0.0]];
        let c = vec![vec![1.0, 0.0]];
        let l = observer_gain_ackermann(&a, &c, &[[-1.0, 0.0], [-2.0, 0.0]])
            .expect("observable system");
        assert!((l[0] - 3.0).abs() < 1e-12, "L = {l:?}");
        assert!((l[1] - 2.0).abs() < 1e-12, "L = {l:?}");
        let cp = characteristic_polynomial(&closed_loop_matrix(&a, &c, &l));
        for (got, want) in cp.iter().zip([1.0, 3.0, 2.0]) {
            assert!((got - want).abs() < 1e-12, "char poly {cp:?}");
        }
    }
    #[test]
    fn test_observer_gain_complex_pole_pair() {
        // Poles -1 +- 2j: phi(s) = s^2 + 2 s + 5, phi(A) = [[5, 2], [0, 5]], O = I
        // => L = phi(A) e_2 = [2, 5].
        let a = vec![vec![0.0, 1.0], vec![0.0, 0.0]];
        let c = vec![vec![1.0, 0.0]];
        let l = observer_gain_ackermann(&a, &c, &[[-1.0, 2.0], [-1.0, -2.0]])
            .expect("conjugate pair is realisable");
        assert!(
            (l[0] - 2.0).abs() < 1e-12 && (l[1] - 5.0).abs() < 1e-12,
            "L = {l:?}"
        );
    }
    #[test]
    fn test_observer_gain_three_states_matches_requested_spectrum() {
        let a = vec![
            vec![1.0, 2.0, 0.0],
            vec![0.0, -1.0, 3.0],
            vec![1.0, 0.0, 2.0],
        ];
        let c = vec![vec![1.0, 0.0, 1.0]];
        let poles = [[-2.0, 0.0], [-3.0, 1.0], [-3.0, -1.0]];
        let l = observer_gain_ackermann(&a, &c, &poles).expect("observable system");
        let closed_loop = closed_loop_matrix(&a, &c, &l);
        // (s + 2)(s^2 + 6 s + 10) = s^3 + 8 s^2 + 22 s + 20.
        let cp = characteristic_polynomial(&closed_loop);
        for (got, want) in cp.iter().zip([1.0, 8.0, 22.0, 20.0]) {
            assert!((got - want).abs() < 1e-9, "char poly {cp:?}");
        }
        let eigenvalues = closed_loop.complex_eigenvalues();
        for pole in &poles {
            let nearest = eigenvalues
                .iter()
                .map(|e| ((e.re - pole[0]).powi(2) + (e.im - pole[1]).powi(2)).sqrt())
                .fold(f64::MAX, f64::min);
            assert!(nearest < 1e-9, "pole {pole:?} missing, nearest {nearest}");
        }
    }
    #[test]
    fn test_observer_gain_rejects_unrealisable_requests() {
        let a = vec![vec![0.0, 1.0], vec![0.0, 0.0]];
        let c = vec![vec![1.0, 0.0]];
        let real_poles = [[-1.0, 0.0], [-2.0, 0.0]];
        // Wrong number of poles.
        assert!(observer_gain_ackermann(&a, &c, &real_poles[..1]).is_none());
        // Two output rows: Ackermann's formula is single-output only.
        let c2 = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(observer_gain_ackermann(&a, &c2, &real_poles).is_none());
        // Output row of the wrong length / ragged A.
        assert!(observer_gain_ackermann(&a, &[vec![1.0]], &real_poles).is_none());
        assert!(observer_gain_ackermann(&[vec![0.0, 1.0], vec![0.0]], &c, &real_poles).is_none());
        // Unobservable pair: x2 is not seen through y = x1 when A = diag(1, 2).
        let diag = vec![vec![1.0, 0.0], vec![0.0, 2.0]];
        assert!(observer_gain_ackermann(&diag, &c, &real_poles).is_none());
        // A lone complex pole has no real realisation.
        assert!(observer_gain_ackermann(&[vec![2.0]], &[vec![1.0]], &[[-1.0, 1.0]]).is_none());
        // Non-finite data.
        assert!(observer_gain_ackermann(&a, &c, &[[f64::NAN, 0.0], [-2.0, 0.0]]).is_none());
    }
    // ── matrix_rank / controllability (pivot row decoupled from column) ──────
    #[test]
    fn test_matrix_rank_wide_matrix_with_late_pivot_columns() {
        // Controllability matrix [B, AB] of A = [[0, 1], [0, 0]], B = [[0, 0], [0, 1]].
        // The old column-index pivot row reported rank 1 here.
        let m = vec![0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0];
        assert_eq!(matrix_rank(&m, 2, 4), 2);
    }
    #[test]
    fn test_matrix_rank_leading_zero_columns() {
        // [[0, 1, 0], [0, 0, 1]] has rank 2 (the old code returned 0).
        let m = vec![0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(matrix_rank(&m, 2, 3), 2);
    }
    #[test]
    fn test_matrix_rank_assorted_shapes() {
        assert_eq!(matrix_rank(&[0.0; 6], 2, 3), 0);
        assert_eq!(matrix_rank(&[], 0, 3), 0);
        // Tall 4x2 with interleaved zero rows.
        assert_eq!(
            matrix_rank(&[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], 4, 2),
            2
        );
        // Tall rank-1: rows are multiples of (1, 2).
        assert_eq!(matrix_rank(&[1.0, 2.0, 2.0, 4.0, 3.0, 6.0], 3, 2), 1);
        // Classic 3x3 of rank 2.
        let m = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        assert_eq!(matrix_rank(&m, 3, 3), 2);
        // Row swap needed: first column starts with zero.
        assert_eq!(matrix_rank(&[0.0, 1.0, 1.0, 0.0], 2, 2), 2);
    }
    #[test]
    fn test_matrix_rank_is_scale_aware() {
        // Rank-1 matrix of large entries: elimination residue ~1e-6 must not
        // be mistaken for a second pivot.
        let big = 1e8;
        let m = vec![big, 2.0 * big, 3.0 * big, 6.000_000_000_001 * big];
        assert_eq!(matrix_rank(&m, 2, 2), 1);
        // A genuinely small but well separated second pivot still counts.
        assert_eq!(matrix_rank(&[1.0, 0.0, 0.0, 1e-6], 2, 2), 2);
    }
    #[test]
    fn test_matrix_rank_short_slice_is_zero_padded() {
        assert_eq!(matrix_rank(&[1.0, 0.0, 0.0], 2, 2), 1);
    }
    #[test]
    fn test_is_controllable_with_two_inputs_and_late_pivots() {
        // The auditor's example: the controllability matrix is [[0,0,0,1],[0,1,0,0]].
        let a = vec![0.0, 1.0, 0.0, 0.0];
        let b = vec![0.0, 0.0, 0.0, 1.0];
        assert!(is_controllable(&a, &b, 2, 2));
        // A single input acting on x2 drives the chain x2 -> x1: controllable.
        let b_single = vec![0.0, 1.0];
        assert!(is_controllable(&a, &b_single, 2, 1));
        // A single input acting on x1 never reaches x2: not controllable.
        let b_dead = vec![1.0, 0.0];
        assert!(!is_controllable(&a, &b_dead, 2, 1));
    }
    #[test]
    fn test_is_observable_tall_observability_matrix() {
        // Two outputs: O = [C; CA] is 4x2.
        let a = vec![0.0, 1.0, 0.0, 0.0];
        let c = vec![0.0, 1.0, 0.0, 1.0];
        // Both rows measure x2 only; x1 never influences the output.
        assert!(!is_observable(&a, &c, 2, 2));
        let c_good = vec![1.0, 0.0, 0.0, 0.0];
        assert!(is_observable(&a, &c_good, 2, 2));
    }
    #[test]
    fn test_pole_placement_stable() {
        let eigs = [[-1.0, 0.0], [-2.0, 0.5]];
        assert!(pole_placement_check(&eigs));
    }
    #[test]
    fn test_pole_placement_unstable() {
        let eigs = [[-1.0, 0.0], [0.1, 0.0]];
        assert!(!pole_placement_check(&eigs));
    }
    #[test]
    fn test_h_infinity_norm_positive() {
        let tf = TransferFunction::first_order(2.0, 1.0);
        let omegas = logspace_omega(0.01, 100.0, 50);
        let norm = h_infinity_norm_bound(&tf, &omegas);
        assert!(norm > 0.0);
    }
    #[test]
    fn test_h2_norm_positive() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let norm = h2_norm(&tf, 100.0, 1000);
        assert!(norm > 0.0);
    }
    #[test]
    fn test_mat_inv_2x2() {
        let a = vec![1.0, 0.0, 0.0, 2.0];
        let inv = mat_inv_small(&a, 2);
        assert!((inv[0] - 1.0).abs() < 1e-10);
        assert!((inv[3] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_mpc_returns_correct_length() {
        let a = vec![1.0, 0.1, 0.0, 1.0];
        let b = vec![0.005, 0.1];
        let q = vec![1.0, 0.0, 0.0, 1.0];
        let r = vec![0.1];
        let mpc = MpcController {
            a,
            b,
            n: 2,
            m: 1,
            horizon: 3,
            q,
            r,
        };
        let u = mpc.compute(&[1.0, 0.0], &[0.0, 0.0]);
        assert_eq!(u.len(), 1);
    }
    #[test]
    fn test_ziegler_nichols() {
        let (kp, ki, kd) = ziegler_nichols_pid(10.0, 2.0);
        assert!((kp - 6.0).abs() < 1e-10);
        assert!((ki - 6.0).abs() < 1e-10);
        assert!((kd - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_cohen_coon_positive_gains() {
        let (kp, ki, kd) = cohen_coon_pid(1.0, 5.0, 1.0);
        assert!(kp > 0.0);
        assert!(ki > 0.0);
        assert!(kd > 0.0);
    }
    #[test]
    fn test_lyapunov_stable_2x2() {
        assert!(lyapunov_stable_2x2(&[-1.0, 0.0, 0.0, -2.0]));
    }
    #[test]
    fn test_lyapunov_unstable_2x2() {
        assert!(!lyapunov_stable_2x2(&[1.0, 0.0, 0.0, -1.0]));
    }
    #[test]
    fn test_solve_lyapunov_1x1() {
        let p = solve_lyapunov_1x1(-2.0, 1.0);
        assert!((p - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_step_response_length() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (t, y) = step_response(&tf, 5.0, 0.01);
        assert_eq!(t.len(), y.len());
        assert!(t.len() > 400);
    }
    #[test]
    fn test_impulse_response_length() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (t, y) = impulse_response(&tf, 5.0, 0.01);
        assert_eq!(t.len(), y.len());
        assert!(t.len() > 400);
    }
    #[test]
    fn test_stability_margins_finite() {
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 3.0, 3.0, 1.0]);
        let (gm, pm) = stability_margins(&tf, 0.01, 100.0, 500);
        assert!(gm.is_finite());
        assert!(pm.is_finite());
    }
    #[test]
    fn test_nyquist_encirclements_stable_system() {
        let tf = TransferFunction::first_order(0.5, 1.0);
        let enc = nyquist_encirclements(&tf, 0.01, 100.0, 500);
        assert_eq!(enc, 0);
    }
    #[test]
    fn test_feedback_dc_gain() {
        let plant = TransferFunction::first_order(10.0, 1.0);
        let ctrl = TransferFunction::new(vec![1.0], vec![1.0]);
        let cl = TransferFunction::feedback(&plant, &ctrl);
        let dc = cl.dc_gain();
        assert!(
            (dc - 10.0 / 11.0).abs() < 0.1,
            "expected ~0.909, got {:.6}",
            dc
        );
    }
}
