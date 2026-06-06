//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Simplified linear MPC controller.
pub struct MpcController {
    /// System matrix A (n-by-n flat).
    pub a: Vec<f64>,
    /// Input matrix B (n-by-m flat).
    pub b: Vec<f64>,
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Prediction horizon.
    pub horizon: usize,
    /// State weight Q (n-by-n flat).
    pub q: Vec<f64>,
    /// Input weight R (m-by-m flat).
    pub r: Vec<f64>,
}
impl MpcController {
    /// Compute the MPC control action for current state toward target.
    ///
    /// Uses a simplified one-step LQR-like approach.
    pub fn compute(&self, state: &[f64], target: &[f64]) -> Vec<f64> {
        let _error: Vec<f64> = state
            .iter()
            .zip(target.iter())
            .map(|(&s, &t)| s - t)
            .collect();
        let lqr = LqrController::solve(
            &self.a,
            &self.b,
            &self.q,
            &self.r,
            LqrDims {
                n: self.n,
                m: self.m,
                max_iter: self.horizon * 10,
            },
        );
        let error_state: Vec<f64> = state
            .iter()
            .zip(target.iter())
            .map(|(&s, &t)| s - t)
            .collect();
        lqr.control(&error_state)
    }
}
/// PID controller with anti-windup, derivative filter, and setpoint weighting.
///
/// Transfer function form:
///   u(t) = Kp * (b*r - y) + Ki * integral(e) + Kd * d/dt(c*r - y)
///
/// where `b` is the proportional setpoint weight and `c` is the derivative
/// setpoint weight.
#[derive(Debug, Clone)]
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Integral accumulator.
    pub integral: f64,
    /// Previous error (for derivative).
    pub prev_error: f64,
    /// Previous filtered derivative value.
    pub prev_deriv_filtered: f64,
    /// Output lower limit.
    pub out_min: f64,
    /// Output upper limit.
    pub out_max: f64,
    /// Derivative filter coefficient (0 = no filter, closer to 1 = more filtering).
    pub deriv_filter_coeff: f64,
    /// Setpoint weighting for proportional term (0 to 1).
    pub setpoint_weight_b: f64,
    /// Setpoint weighting for derivative term (0 to 1).
    pub setpoint_weight_c: f64,
    /// Anti-windup: back-calculation gain.
    pub anti_windup_gain: f64,
    /// Previous unclipped output (for anti-windup).
    pub prev_unclipped: f64,
}
impl PidController {
    /// Create a PID controller with given gains.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            prev_deriv_filtered: 0.0,
            out_min: f64::NEG_INFINITY,
            out_max: f64::INFINITY,
            deriv_filter_coeff: 0.0,
            setpoint_weight_b: 1.0,
            setpoint_weight_c: 1.0,
            anti_windup_gain: 0.0,
            prev_unclipped: 0.0,
        }
    }
    /// Set output limits and return self.
    pub fn with_limits(mut self, min: f64, max: f64) -> Self {
        self.out_min = min;
        self.out_max = max;
        self
    }
    /// Set derivative filter coefficient (0..1) and return self.
    pub fn with_deriv_filter(mut self, alpha: f64) -> Self {
        self.deriv_filter_coeff = alpha.clamp(0.0, 0.999);
        self
    }
    /// Set setpoint weights and return self.
    pub fn with_setpoint_weights(mut self, b: f64, c: f64) -> Self {
        self.setpoint_weight_b = b;
        self.setpoint_weight_c = c;
        self
    }
    /// Set anti-windup back-calculation gain and return self.
    pub fn with_anti_windup(mut self, gain: f64) -> Self {
        self.anti_windup_gain = gain;
        self
    }
}
/// Dimensions and iteration limit for LQR solving.
#[derive(Debug, Clone, Copy)]
pub struct LqrDims {
    /// State dimension `n`.
    pub n: usize,
    /// Input dimension `m`.
    pub m: usize,
    /// Maximum Riccati iterations.
    pub max_iter: usize,
}

/// LQR controller: solve the discrete algebraic Riccati equation iteratively.
///
/// Minimizes J = sum(x'Qx + u'Ru).
pub struct LqrController {
    /// Optimal gain matrix K (m-by-n, flat).
    pub k: Vec<f64>,
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
}
impl LqrController {
    /// Solve LQR by iterating the Riccati equation.
    ///
    /// `a` is n-by-n, `b` is n-by-m, `q` is n-by-n, `r` is m-by-m.
    /// Dimensions and iteration limit are given via `dims`.
    pub fn solve(a: &[f64], b: &[f64], q: &[f64], r: &[f64], dims: LqrDims) -> Self {
        let n = dims.n;
        let m = dims.m;
        let max_iter = dims.max_iter;
        let bt = mat_transpose(b, n, m);
        let mut p = q.to_vec();
        for _iter in 0..max_iter {
            let pb = mat_mul(&p, b, n, n, m);
            let bpb = mat_mul(&bt, &pb, m, n, m);
            let s = mat_add(r, &bpb);
            let s_inv = mat_inv_small(&s, m);
            let pa = mat_mul(&p, a, n, n, n);
            let bpa = mat_mul(&bt, &pa, m, n, n);
            let k = mat_mul(&s_inv, &bpa, m, m, n);
            let at = mat_transpose(a, n, n);
            let ap = mat_mul(&at, &p, n, n, n);
            let apa = mat_mul(&ap, a, n, n, n);
            let apb = mat_mul(&ap, b, n, n, m);
            let apbk = mat_mul(&apb, &k, n, m, n);
            let p_new_1 = mat_add(q, &apa);
            let p_new = mat_sub(&p_new_1, &apbk);
            let diff: f64 = p_new
                .iter()
                .zip(p.iter())
                .map(|(&a, &b)| (a - b).abs())
                .sum();
            p = p_new;
            if diff < 1e-12 {
                break;
            }
        }
        let pb = mat_mul(&p, b, n, n, m);
        let bpb = mat_mul(&bt, &pb, m, n, m);
        let s = mat_add(r, &bpb);
        let s_inv = mat_inv_small(&s, m);
        let pa = mat_mul(&p, a, n, n, n);
        let bpa = mat_mul(&bt, &pa, m, n, n);
        let k = mat_mul(&s_inv, &bpa, m, m, n);
        Self { k, n, m }
    }
    /// Compute optimal control: u = -K * x.
    pub fn control(&self, x: &[f64]) -> Vec<f64> {
        let mut u = vec![0.0; self.m];
        for (i, ui) in u.iter_mut().enumerate() {
            for (j, &xj) in x.iter().enumerate() {
                *ui -= self.k[i * self.n + j] * xj;
            }
        }
        u
    }
}
/// SISO transfer function H(s) = num(s) / den(s).
///
/// Coefficients are in descending power order: `[a_n, a_{n-1}, ..., a_0]`.
#[derive(Debug, Clone)]
pub struct TransferFunction {
    /// Numerator polynomial coefficients.
    pub num: Vec<f64>,
    /// Denominator polynomial coefficients.
    pub den: Vec<f64>,
}
impl TransferFunction {
    /// Create a transfer function from numerator and denominator polynomials.
    pub fn new(num: Vec<f64>, den: Vec<f64>) -> Self {
        Self { num, den }
    }
    /// First-order system: K / (tau*s + 1).
    pub fn first_order(gain: f64, time_constant: f64) -> Self {
        Self {
            num: vec![gain],
            den: vec![time_constant, 1.0],
        }
    }
    /// Second-order system: wn^2 / (s^2 + 2*zeta*wn*s + wn^2).
    pub fn second_order(omega_n: f64, zeta: f64) -> Self {
        Self {
            num: vec![omega_n * omega_n],
            den: vec![1.0, 2.0 * zeta * omega_n, omega_n * omega_n],
        }
    }
    /// DC gain: H(0) = num(0) / den(0).
    pub fn dc_gain(&self) -> f64 {
        let n0 = self.num.last().copied().unwrap_or(0.0);
        let d0 = self.den.last().copied().unwrap_or(1.0);
        if d0.abs() < 1e-300 {
            return f64::INFINITY;
        }
        n0 / d0
    }
    /// Evaluate H(jw) at a given frequency.
    pub fn eval_jw(&self, omega: f64) -> Complex {
        let s = Complex::jw(omega);
        let num = poly_eval_complex(&self.num, s);
        let den = poly_eval_complex(&self.den, s);
        num.div(&den)
    }
    /// Return the poles (roots of denominator) via companion matrix eigenvalues (2x2 only).
    /// For higher order, returns approximate roots via Durand-Kerner iteration.
    pub fn poles(&self) -> Vec<Complex> {
        poly_roots(&self.den)
    }
    /// Return the zeros (roots of numerator).
    pub fn zeros(&self) -> Vec<Complex> {
        poly_roots(&self.num)
    }
    /// Check BIBO stability: all poles must have negative real part.
    pub fn is_stable(&self) -> bool {
        let poles = self.poles();
        poles.iter().all(|p| p.re < 0.0)
    }
    /// Series connection: H1(s) * H2(s).
    pub fn series(h1: &TransferFunction, h2: &TransferFunction) -> TransferFunction {
        TransferFunction {
            num: poly_mul(&h1.num, &h2.num),
            den: poly_mul(&h1.den, &h2.den),
        }
    }
    /// Feedback connection: H(s) / (1 + H(s)*G(s)) for negative feedback.
    pub fn feedback(plant: &TransferFunction, controller: &TransferFunction) -> TransferFunction {
        let open_num = poly_mul(&plant.num, &controller.num);
        let open_den = poly_mul(&plant.den, &controller.den);
        let cl_den = poly_add(&open_den, &open_num);
        TransferFunction {
            num: open_num,
            den: cl_den,
        }
    }
}
/// State-space model: dx/dt = Ax + Bu, y = Cx + Du.
///
/// Matrices stored as `Vec<Vec`f64`>` (row-major).
#[derive(Debug, Clone)]
pub struct StateSpaceModel {
    /// System matrix (n x n).
    pub a: Vec<Vec<f64>>,
    /// Input matrix (n x m).
    pub b: Vec<Vec<f64>>,
    /// Output matrix (p x n).
    pub c: Vec<Vec<f64>>,
    /// Feedthrough matrix (p x m).
    pub d: Vec<Vec<f64>>,
}
impl StateSpaceModel {
    /// Create a state-space model.
    pub fn new(a: Vec<Vec<f64>>, b: Vec<Vec<f64>>, c: Vec<Vec<f64>>, d: Vec<Vec<f64>>) -> Self {
        Self { a, b, c, d }
    }
    /// Return the state dimension (n).
    pub fn state_dim(&self) -> usize {
        self.a.len()
    }
    /// Return the input dimension (m).
    pub fn input_dim(&self) -> usize {
        if self.b.is_empty() {
            0
        } else {
            self.b[0].len()
        }
    }
    /// Return the output dimension (p).
    pub fn output_dim(&self) -> usize {
        self.c.len()
    }
}
/// Linear Kalman filter for state estimation.
///
/// State model: x_{k+1} = F*x_k + w_k, z_k = H*x_k + v_k
/// where w_k ~ N(0,Q), v_k ~ N(0,R).
pub struct KalmanFilter {
    /// State transition matrix (n-by-n, flat).
    pub f: Vec<f64>,
    /// Observation matrix (m-by-n, flat).
    pub h: Vec<f64>,
    /// Process noise covariance (n-by-n, flat).
    pub q: Vec<f64>,
    /// Measurement noise covariance (m-by-m, flat).
    pub r: Vec<f64>,
    /// State estimate (n).
    pub x: Vec<f64>,
    /// Error covariance (n-by-n, flat).
    pub p: Vec<f64>,
    /// State dimension.
    pub n: usize,
    /// Measurement dimension.
    pub m: usize,
}
impl KalmanFilter {
    /// Create a new Kalman filter.
    pub fn new(f: Vec<f64>, h: Vec<f64>, q: Vec<f64>, r: Vec<f64>, n: usize, m: usize) -> Self {
        let x = vec![0.0; n];
        let p = mat_eye(n);
        Self {
            f,
            h,
            q,
            r,
            x,
            p,
            n,
            m,
        }
    }
    /// Predict step: x = F*x, P = F*P*F' + Q.
    pub fn predict(&mut self) {
        let new_x = mat_mul(&self.f, &self.x, self.n, self.n, 1);
        let fp = mat_mul(&self.f, &self.p, self.n, self.n, self.n);
        let ft = mat_transpose(&self.f, self.n, self.n);
        let fpf = mat_mul(&fp, &ft, self.n, self.n, self.n);
        self.p = mat_add(&fpf, &self.q);
        self.x = new_x;
    }
    /// Update step: incorporate measurement z.
    pub fn update(&mut self, z: &[f64]) {
        let hx = mat_mul(&self.h, &self.x, self.m, self.n, 1);
        let y: Vec<f64> = z.iter().zip(hx.iter()).map(|(&a, &b)| a - b).collect();
        let hp = mat_mul(&self.h, &self.p, self.m, self.n, self.n);
        let ht = mat_transpose(&self.h, self.m, self.n);
        let hph = mat_mul(&hp, &ht, self.m, self.n, self.m);
        let s = mat_add(&hph, &self.r);
        let ph = mat_mul(&self.p, &ht, self.n, self.n, self.m);
        let s_inv = mat_inv_small(&s, self.m);
        let k = mat_mul(&ph, &s_inv, self.n, self.m, self.m);
        let ky = mat_mul(&k, &y, self.n, self.m, 1);
        self.x = mat_add(&self.x, &ky);
        let kh = mat_mul(&k, &self.h, self.n, self.m, self.n);
        let ikh = mat_sub(&mat_eye(self.n), &kh);
        self.p = mat_mul(&ikh, &self.p, self.n, self.n, self.n);
    }
    /// Run a full predict-update cycle and return the current state estimate.
    pub fn step(&mut self, z: &[f64]) -> Vec<f64> {
        self.predict();
        self.update(z);
        self.x.clone()
    }
}
/// A complex number (real, imaginary).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}
impl Complex {
    /// Create a new complex number.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    /// Magnitude |z|.
    pub fn abs(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    /// Angle (argument) in radians.
    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }
    /// Multiply two complex numbers.
    pub fn mul(&self, other: &Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
    /// Divide two complex numbers.
    pub fn div(&self, other: &Complex) -> Complex {
        let denom = other.re * other.re + other.im * other.im;
        if denom < 1e-300 {
            return Complex::new(0.0, 0.0);
        }
        Complex {
            re: (self.re * other.re + self.im * other.im) / denom,
            im: (self.im * other.re - self.re * other.im) / denom,
        }
    }
    /// Add two complex numbers.
    pub fn add(&self, other: &Complex) -> Complex {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
    /// Subtract two complex numbers.
    pub fn sub(&self, other: &Complex) -> Complex {
        Complex {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
    /// Complex conjugate.
    pub fn conj(&self) -> Complex {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }
    /// Magnitude squared.
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    /// Convert from polar form (magnitude, angle_rad).
    pub fn from_polar(mag: f64, angle: f64) -> Self {
        Self {
            re: mag * angle.cos(),
            im: mag * angle.sin(),
        }
    }
    /// Pure imaginary: jw.
    pub fn jw(omega: f64) -> Self {
        Self { re: 0.0, im: omega }
    }
    /// Scale by a real factor.
    pub fn scale(&self, k: f64) -> Complex {
        Complex {
            re: self.re * k,
            im: self.im * k,
        }
    }
}
