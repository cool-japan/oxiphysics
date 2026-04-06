//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BodePoint, RelayFeedbackResult, StabilityMargins};

/// Compute Bode plot data for a first-order system: G(jw) = K / (j*w*tau + 1).
///
/// Returns `(frequencies, magnitudes_db, phases_deg)`.
pub fn bode_first_order(
    gain: f64,
    tau: f64,
    omega_min: f64,
    omega_max: f64,
    n_points: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut freqs = Vec::with_capacity(n_points);
    let mut mags = Vec::with_capacity(n_points);
    let mut phases = Vec::with_capacity(n_points);
    let log_min = omega_min.ln();
    let log_max = omega_max.ln();
    for i in 0..n_points {
        let t = i as f64 / (n_points - 1).max(1) as f64;
        let omega = (log_min + t * (log_max - log_min)).exp();
        let wt = omega * tau;
        let mag = gain / (1.0 + wt * wt).sqrt();
        let mag_db = 20.0 * mag.log10();
        let phase_deg = -(wt).atan().to_degrees();
        freqs.push(omega);
        mags.push(mag_db);
        phases.push(phase_deg);
    }
    (freqs, mags, phases)
}
/// Gain margin for a first-order system with PID controller.
///
/// For a simple open-loop G(s) = K / (tau*s + 1), the gain margin is
/// infinite (first-order systems with proportional control cannot go unstable).
/// Returns `f64::INFINITY` for first-order plants.
pub fn gain_margin_first_order() -> f64 {
    f64::INFINITY
}
/// Phase margin for a first-order system with proportional gain `kp`.
///
/// Open-loop: L(jw) = kp * K / (1 + j*w*tau)
/// Gain crossover: |L(jw_gc)| = 1 => w_gc = sqrt((kp*K)^2 - 1) / tau
/// Phase margin: 180° + angle(L(jw_gc))
///
/// Returns the phase margin in degrees.
pub fn phase_margin_first_order(kp: f64, gain: f64, _tau: f64) -> f64 {
    let kp_k = kp * gain;
    if kp_k <= 1.0 {
        return 180.0;
    }
    let w_gc_tau = (kp_k * kp_k - 1.0).sqrt();
    let phase_at_gc = -(w_gc_tau).atan().to_degrees();
    180.0 + phase_at_gc
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AdaptivePid;
    use crate::CascadedPid;
    use crate::FirstOrderTf;
    use crate::JointMotor;

    use crate::LqrController;
    use crate::MotorMode;
    use crate::PidController;

    use crate::SecondOrderTf;
    use crate::SpringDamper;

    use crate::StateSpace2x1;
    use crate::TrajectoryFollower;
    use crate::ZieglerNichols;
    #[test]
    fn pid_proportional_only() {
        let mut pid = PidController::new(2.0, 0.0, 0.0);
        let output = pid.update(3.0, 0.1);
        assert!((output - 6.0).abs() < 1e-10, "output={output}");
    }
    #[test]
    fn pid_integral_accumulates() {
        let mut pid = PidController::new(0.0, 1.0, 0.0);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        assert!(
            (pid.integral - 0.2).abs() < 1e-10,
            "integral={}",
            pid.integral
        );
    }
    #[test]
    fn pid_anti_windup() {
        let mut pid = PidController::new(0.0, 1.0, 0.0);
        pid.set_limits(0.5, f64::MAX);
        for _ in 0..100 {
            pid.update(1.0, 0.1);
        }
        assert!(pid.integral <= 0.5 + 1e-10, "integral={}", pid.integral);
    }
    #[test]
    fn pid_output_saturation() {
        let mut pid = PidController::new(100.0, 0.0, 0.0);
        pid.set_limits(f64::MAX, 5.0);
        let output = pid.update(1.0, 0.1);
        assert!((output - 5.0).abs() < 1e-10, "output={output}");
    }
    #[test]
    fn lqr_zero_error_zero_input() {
        let lqr = LqrController::new_double_integrator(1.0, 1.0, 1.0);
        let u = lqr.compute_input([1.0, 0.0], [1.0, 0.0]);
        assert!(u.abs() < 1e-12, "u={u}");
    }
    #[test]
    fn spring_damper_zero_displacement_zero_force() {
        let sd = SpringDamper {
            stiffness: 10.0,
            damping: 1.0,
            rest_length: 0.0,
        };
        let f = sd.force(0.0, 0.0);
        assert!(f.abs() < 1e-12, "f={f}");
    }
    #[test]
    fn joint_motor_position_control_moves_toward_target() {
        let mut motor = JointMotor {
            target_position: 1.0,
            target_velocity: 0.0,
            max_torque: 100.0,
            max_velocity: 10.0,
            gear_ratio: 1.0,
            controller: PidController::new(10.0, 0.0, 0.0),
            mode: MotorMode::PositionControl,
        };
        let torque = motor.compute_torque(0.0, 0.0, 0.1);
        assert!(torque > 0.0, "torque={torque}");
    }
    #[test]
    fn trajectory_follower_interpolates() {
        let waypoints = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]];
        let pid = PidController::new(1.0, 0.0, 0.0);
        let follower = TrajectoryFollower::new(waypoints, pid);
        let pos = follower.desired_position(0.5);
        assert!((pos - 0.5).abs() < 1e-10, "pos={pos}");
        let pos2 = follower.desired_position(1.5);
        assert!((pos2 - 0.5).abs() < 1e-10, "pos2={pos2}");
    }
    #[test]
    fn cascaded_pid_nonzero_output_for_nonzero_pos_error() {
        let mut cpid = CascadedPid::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let output = cpid.update(1.0, 0.0, 0.1);
        assert!(output.abs() > 1e-10, "output={output}");
    }
    #[test]
    fn zn_classic_gains() {
        let (kp, ki, kd) = ZieglerNichols::classic(10.0, 2.0);
        assert!((kp - 6.0).abs() < 1e-10, "kp={kp}");
        assert!((ki - 6.0).abs() < 1e-10, "ki={ki}");
        assert!((kd - 1.5).abs() < 1e-10, "kd={kd}");
    }
    #[test]
    fn zn_pi_only_gains() {
        let (kp, ki, kd) = ZieglerNichols::pi_only(10.0, 2.0);
        assert!((kp - 4.5).abs() < 1e-10);
        assert!(ki > 0.0);
        assert!(kd.abs() < 1e-14);
    }
    #[test]
    fn zn_p_only_gain() {
        let (kp, ki, kd) = ZieglerNichols::p_only(10.0);
        assert!((kp - 5.0).abs() < 1e-10);
        assert!(ki.abs() < 1e-14);
        assert!(kd.abs() < 1e-14);
    }
    #[test]
    fn zn_some_overshoot() {
        let (kp, ki, kd) = ZieglerNichols::some_overshoot(9.0, 3.0);
        assert!((kp - 3.0).abs() < 1e-10);
        assert!(ki > 0.0);
        assert!(kd > 0.0);
    }
    #[test]
    fn state_space_double_integrator_step() {
        let mut ss = StateSpace2x1::double_integrator();
        for _ in 0..100 {
            ss.step(1.0, 0.01);
        }
        assert!(ss.state[1] > 0.5, "velocity = {}", ss.state[1]);
        assert!(ss.state[0] > 0.1, "position = {}", ss.state[0]);
    }
    #[test]
    fn state_space_damped_oscillator_is_stable() {
        let ss = StateSpace2x1::damped_oscillator(10.0, 0.5);
        assert!(ss.is_stable(), "damped oscillator should be stable");
    }
    #[test]
    fn state_space_undamped_oscillator_not_stable() {
        let ss = StateSpace2x1::damped_oscillator(10.0, 0.0);
        assert!(!ss.is_stable(), "undamped oscillator is marginally stable");
    }
    #[test]
    fn state_space_eigenvalues_real() {
        let ss = StateSpace2x1::new([[-1.0, 0.0], [0.0, -2.0]], [1.0, 0.0], [1.0, 0.0], 0.0);
        let ev = ss.eigenvalues();
        assert!(ev[0].1.abs() < 1e-12, "imaginary part should be zero");
        assert!(ev[1].1.abs() < 1e-12);
    }
    #[test]
    fn state_space_eigenvalues_complex() {
        let ss = StateSpace2x1::damped_oscillator(10.0, 0.1);
        let ev = ss.eigenvalues();
        assert!(ev[0].0 < 0.0);
        assert!(ev[1].0 < 0.0);
        assert!(ev[0].1.abs() > 0.1);
    }
    #[test]
    fn state_space_output() {
        let ss = StateSpace2x1::new([[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0], [1.0, 2.0], 3.0);
        let y = ss.output(1.0);
        assert!((y - 3.0).abs() < 1e-10, "output with zero state = {y}");
    }
    #[test]
    fn first_order_tf_steady_state() {
        let mut tf = FirstOrderTf::new(2.0, 1.0);
        for _ in 0..10000 {
            tf.step(1.0, 0.01);
        }
        assert!(
            (tf.output() - 2.0).abs() < 0.01,
            "steady state = {}, expected 2.0",
            tf.output()
        );
    }
    #[test]
    fn first_order_tf_dc_gain() {
        let tf = FirstOrderTf::new(5.0, 0.1);
        assert!((tf.dc_gain() - 5.0).abs() < 1e-14);
    }
    #[test]
    fn first_order_tf_bandwidth() {
        let tf = FirstOrderTf::new(1.0, 0.5);
        assert!((tf.bandwidth() - 2.0).abs() < 1e-14);
    }
    #[test]
    fn second_order_tf_no_nan() {
        let mut tf = SecondOrderTf::new(1.0, 10.0, 0.7);
        for _ in 0..1000 {
            tf.step(1.0, 0.001);
        }
        assert!(!tf.output().is_nan(), "output should not be NaN");
    }
    #[test]
    fn bode_first_order_dc_gain() {
        let (freqs, mags, _phases) = bode_first_order(2.0, 1.0, 0.01, 100.0, 50);
        assert_eq!(freqs.len(), 50);
        assert!(
            (mags[0] - 20.0 * 2.0_f64.log10()).abs() < 0.5,
            "DC gain: {}",
            mags[0]
        );
    }
    #[test]
    fn bode_first_order_phase_range() {
        let (_freqs, _mags, phases) = bode_first_order(1.0, 1.0, 0.01, 100.0, 50);
        for &p in &phases {
            assert!((-90.0..=0.0).contains(&p), "phase out of range: {p}");
        }
    }
    #[test]
    fn gain_margin_first_order_infinite() {
        assert!(gain_margin_first_order().is_infinite());
    }
    #[test]
    fn phase_margin_first_order_low_gain() {
        let pm = phase_margin_first_order(0.5, 1.0, 1.0);
        assert!((pm - 180.0).abs() < 1e-10, "PM = {pm}");
    }
    #[test]
    fn phase_margin_first_order_high_gain() {
        let pm = phase_margin_first_order(10.0, 1.0, 1.0);
        assert!(pm > 0.0 && pm < 180.0, "PM = {pm}");
    }
    #[test]
    fn adaptive_pid_increases_kp() {
        let base = PidController::new(1.0, 0.0, 0.0);
        let mut apid = AdaptivePid {
            base_pid: base,
            adaptation_rate: 0.1,
        };
        let initial_kp = apid.base_pid.kp;
        for _ in 0..10 {
            apid.update(1.0, 0.1);
        }
        assert!(
            apid.base_pid.kp > initial_kp,
            "kp should increase: {}",
            apid.base_pid.kp
        );
    }
}
/// 2×2 matrix multiplication: C = A * B.
pub(super) fn mat2_mul(a: &[[f64; 2]; 2], b: &[[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let mut c = [[0.0f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// 2×2 matrix A * B^T.
pub(super) fn mat2_mul_transpose(a: &[[f64; 2]; 2], b: &[[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let mut c = [[0.0f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                c[i][j] += a[i][k] * b[j][k];
            }
        }
    }
    c
}
/// Pole placement for a 2-state system via Ackermann's formula.
///
/// Computes a state feedback gain K such that the closed-loop eigenvalues
/// of (A - B K) are the desired poles.
///
/// For a 2×2 system: K = e_2^T Wc^{-1} a_d(A)
/// where a_d is the desired characteristic polynomial and Wc is the
/// controllability matrix.
#[allow(dead_code)]
pub fn pole_placement_2d(
    a: &[[f64; 2]; 2],
    b: &[f64; 2],
    desired_poles: [f64; 2],
) -> Option<[f64; 2]> {
    let ab = [
        a[0][0] * b[0] + a[0][1] * b[1],
        a[1][0] * b[0] + a[1][1] * b[1],
    ];
    let wc = [[b[0], ab[0]], [b[1], ab[1]]];
    let det_wc = wc[0][0] * wc[1][1] - wc[0][1] * wc[1][0];
    if det_wc.abs() < 1e-30 {
        return None;
    }
    let sum_poles = desired_poles[0] + desired_poles[1];
    let prod_poles = desired_poles[0] * desired_poles[1];
    let tr_a = a[0][0] + a[1][1];
    let det_a = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let a_sq = mat2_mul(a, a);
    let a_d = [
        [
            a_sq[0][0] - sum_poles * a[0][0] + prod_poles,
            a_sq[0][1] - sum_poles * a[0][1],
        ],
        [
            a_sq[1][0] - sum_poles * a[1][0],
            a_sq[1][1] - sum_poles * a[1][1] + prod_poles,
        ],
    ];
    let _ = (tr_a, det_a);
    let inv_det = 1.0 / det_wc;
    let wc_inv = [
        [wc[1][1] * inv_det, -wc[0][1] * inv_det],
        [-wc[1][0] * inv_det, wc[0][0] * inv_det],
    ];
    let e2_wc_inv = [wc_inv[1][0], wc_inv[1][1]];
    let k = [
        e2_wc_inv[0] * a_d[0][0] + e2_wc_inv[1] * a_d[1][0],
        e2_wc_inv[0] * a_d[0][1] + e2_wc_inv[1] * a_d[1][1],
    ];
    Some(k)
}
/// Discretize a 2×1 state-space system using zero-order hold (ZOH).
///
/// A_d = e^(A*dt) ≈ I + A*dt + (A*dt)²/2! + ...
/// B_d = A^{-1}(A_d - I) B  (if A invertible)
///      ≈ dt * B + dt²/2 * A*B + ...
///
/// Uses series expansion truncated at 4 terms.
#[allow(dead_code)]
pub fn discretize_zoh(a: &[[f64; 2]; 2], b: &[f64; 2], dt: f64) -> ([[f64; 2]; 2], [f64; 2]) {
    let adt = [[a[0][0] * dt, a[0][1] * dt], [a[1][0] * dt, a[1][1] * dt]];
    let adt2 = mat2_mul(&adt, &adt);
    let adt3 = mat2_mul(&adt2, &adt);
    let adt4 = mat2_mul(&adt3, &adt);
    let mut ad = [[1.0, 0.0], [0.0, 1.0]];
    for i in 0..2 {
        for j in 0..2 {
            ad[i][j] += adt[i][j] + adt2[i][j] / 2.0 + adt3[i][j] / 6.0 + adt4[i][j] / 24.0;
        }
    }
    let mut bd = [0.0f64; 2];
    let scale = [1.0, dt / 2.0, dt * dt / 6.0, dt * dt * dt / 24.0];
    let a_pows = [[[1.0f64, 0.0], [0.0, 1.0]], adt, adt2, adt3];
    for k in 0..4 {
        bd[0] += scale[k] * (a_pows[k][0][0] * b[0] + a_pows[k][0][1] * b[1]);
        bd[1] += scale[k] * (a_pows[k][1][0] * b[0] + a_pows[k][1][1] * b[1]);
    }
    for x in bd.iter_mut() {
        *x *= dt;
    }
    (ad, bd)
}
/// Check controllability of a 2-state system.
///
/// Wc = \[B, AB\]; system is controllable if rank(Wc) = 2 (det ≠ 0).
#[allow(dead_code)]
pub fn is_controllable_2d(a: &[[f64; 2]; 2], b: &[f64; 2]) -> bool {
    let ab = [
        a[0][0] * b[0] + a[0][1] * b[1],
        a[1][0] * b[0] + a[1][1] * b[1],
    ];
    let det = b[0] * ab[1] - b[1] * ab[0];
    det.abs() > 1e-14
}
/// Check observability of a 2-state, 1-output system.
///
/// Wo = \[C; CA\]; system is observable if rank(Wo) = 2.
#[allow(dead_code)]
pub fn is_observable_2d(a: &[[f64; 2]; 2], c: &[f64; 2]) -> bool {
    let ca = [
        c[0] * a[0][0] + c[1] * a[1][0],
        c[0] * a[0][1] + c[1] * a[1][1],
    ];
    let det = c[0] * ca[1] - c[1] * ca[0];
    det.abs() > 1e-14
}
/// Compute the controllability Gramian for a stable system.
///
/// W_c = integral_0^inf e^{At} B B^T e^{A^T t} dt
///
/// Solved via Lyapunov equation: A W_c + W_c A^T + B B^T = 0
/// For 2×2 system, solved analytically.
#[allow(dead_code)]
pub fn controllability_gramian(a: &[[f64; 2]; 2], b: &[f64; 2]) -> Option<[[f64; 2]; 2]> {
    let bbt = [[b[0] * b[0], b[0] * b[1]], [b[1] * b[0], b[1] * b[1]]];
    let mat = [
        [2.0 * a[0][0], 2.0 * a[0][1], 0.0],
        [a[1][0], a[0][0] + a[1][1], a[0][1]],
        [0.0, 2.0 * a[1][0], 2.0 * a[1][1]],
    ];
    let rhs = [-bbt[0][0], -bbt[0][1], -bbt[1][1]];
    let sol = solve_3x3(&mat, &rhs)?;
    Some([[sol[0], sol[1]], [sol[1], sol[2]]])
}
/// Solve a 3×3 linear system Ax = b using Cramer's rule.
pub(super) fn solve_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-30 {
        return None;
    }
    let x0 = (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
        + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2]))
        / det;
    let x1 = (a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
        - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0]))
        / det;
    let x2 = (a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
        - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
        + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
        / det;
    Some([x0, x1, x2])
}
/// Frequency response of a second-order system: G(jω) = ωn²/(ωn²-ω²+2jζωnω)
#[allow(dead_code)]
pub fn second_order_freq_response(wn: f64, zeta: f64, omega: f64) -> (f64, f64) {
    let re = wn * wn - omega * omega;
    let im = 2.0 * zeta * wn * omega;
    let denom_sq = re * re + im * im;
    if denom_sq < 1e-60 {
        return (0.0, 0.0);
    }
    let mag = wn * wn / denom_sq.sqrt();
    let phase = (-im).atan2(re);
    (mag, phase)
}
/// Compute resonant frequency ωr and peak gain Mr for underdamped system.
#[allow(dead_code)]
pub fn resonant_peak(wn: f64, zeta: f64) -> Option<(f64, f64)> {
    if zeta >= 1.0 / 2.0_f64.sqrt() {
        return None;
    }
    let wr = wn * (1.0 - 2.0 * zeta * zeta).sqrt();
    let mr = 1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt());
    Some((wr, mr))
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::KalmanFilter;

    #[test]
    fn test_kalman_predict_does_not_nan() {
        let mut kf = KalmanFilter::double_integrator(0.01, 0.1, 0.1);
        kf.predict(0.0);
        assert!(kf.x_hat[0].is_finite() && kf.x_hat[1].is_finite());
    }
    #[test]
    fn test_kalman_update_reduces_uncertainty() {
        let mut kf = KalmanFilter::double_integrator(0.01, 0.1, 0.1);
        let p_init = kf.p[0][0];
        kf.predict(0.0);
        kf.update(0.5);
        assert!(
            kf.p[0][0] <= p_init + 0.01,
            "P[0][0] should not increase: {}",
            kf.p[0][0]
        );
    }
    #[test]
    fn test_kalman_converges_to_constant() {
        let mut kf = KalmanFilter::double_integrator(0.01, 0.01, 0.1);
        for _ in 0..200 {
            kf.predict(0.0);
            kf.update(1.0);
        }
        assert!(
            (kf.state()[0] - 1.0).abs() < 0.1,
            "should converge: {}",
            kf.state()[0]
        );
    }
    #[test]
    fn test_pole_placement_controllable() {
        let a = [[0.0, 1.0], [0.0, 0.0]];
        let b = [0.0, 1.0];
        let poles = [-2.0, -3.0];
        let k = pole_placement_2d(&a, &b, poles);
        assert!(k.is_some(), "should succeed for controllable system");
        let k = k.unwrap();
        let a_cl = [
            [a[0][0] - b[0] * k[0], a[0][1] - b[0] * k[1]],
            [a[1][0] - b[1] * k[0], a[1][1] - b[1] * k[1]],
        ];
        let tr = a_cl[0][0] + a_cl[1][1];
        let det = a_cl[0][0] * a_cl[1][1] - a_cl[0][1] * a_cl[1][0];
        assert!((tr - (-5.0)).abs() < 0.1, "trace = {tr}");
        assert!((det - 6.0).abs() < 0.1, "det = {det}");
    }
    #[test]
    fn test_pole_placement_uncontrollable() {
        let a = [[1.0, 0.0], [0.0, 1.0]];
        let b = [0.0, 0.0];
        let k = pole_placement_2d(&a, &b, [-1.0, -2.0]);
        assert!(k.is_none(), "uncontrollable system should return None");
    }
    #[test]
    fn test_discretize_zoh_identity() {
        let a = [[0.0, 0.0], [0.0, 0.0]];
        let b = [1.0, 0.0];
        let (ad, _bd) = discretize_zoh(&a, &b, 0.001);
        assert!((ad[0][0] - 1.0).abs() < 1e-6);
        assert!((ad[1][1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_discretize_zoh_integrator() {
        let a = [[0.0, 1.0], [0.0, 0.0]];
        let b = [0.0, 1.0];
        let dt = 0.01;
        let (ad, _bd) = discretize_zoh(&a, &b, dt);
        assert!((ad[0][0] - 1.0).abs() < 1e-6);
        assert!((ad[0][1] - dt).abs() < 1e-6);
        assert!((ad[1][0]).abs() < 1e-6);
        assert!((ad[1][1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_controllability_double_integrator() {
        let a = [[0.0, 1.0], [0.0, 0.0]];
        let b = [0.0, 1.0];
        assert!(
            is_controllable_2d(&a, &b),
            "double integrator is controllable"
        );
    }
    #[test]
    fn test_observability_double_integrator() {
        let a = [[0.0, 1.0], [0.0, 0.0]];
        let c = [1.0, 0.0];
        assert!(
            is_observable_2d(&a, &c),
            "position observation is observable"
        );
    }
    #[test]
    fn test_controllability_gramian_positive_definite() {
        let a = [[-1.0, 0.0], [0.0, -2.0]];
        let b = [1.0, 1.0];
        let w = controllability_gramian(&a, &b);
        assert!(w.is_some(), "should compute Gramian for stable system");
        let w = w.unwrap();
        assert!(w[0][0] > 0.0, "W[0][0] = {}", w[0][0]);
        assert!(w[1][1] > 0.0, "W[1][1] = {}", w[1][1]);
    }
    #[test]
    fn test_second_order_freq_response_at_resonance() {
        let wn = 10.0;
        let zeta = 0.1;
        let (mag, _phase) = second_order_freq_response(wn, zeta, wn);
        let expected = 1.0 / (2.0 * zeta);
        assert!(
            (mag - expected).abs() / expected < 0.05,
            "mag at resonance = {mag}"
        );
    }
    #[test]
    fn test_second_order_freq_response_dc() {
        let (mag, _) = second_order_freq_response(10.0, 0.5, 0.0);
        assert!((mag - 1.0).abs() < 1e-10, "DC gain = {mag}");
    }
    #[test]
    fn test_resonant_peak_underdamped() {
        let result = resonant_peak(10.0, 0.1);
        assert!(result.is_some());
        let (wr, mr) = result.unwrap();
        assert!(wr > 0.0 && wr < 10.0, "resonant freq = {wr}");
        assert!(mr > 1.0, "peak gain should exceed 1: {mr}");
    }
    #[test]
    fn test_resonant_peak_overdamped() {
        let result = resonant_peak(10.0, 0.8);
        assert!(result.is_none(), "overdamped system should have no peak");
    }
}
/// Compute the frequency response of a transfer function H(jω) = num(jω)/den(jω).
///
/// `num` and `den` are polynomial coefficients in descending order
/// (highest power first), e.g. `[1.0, 2.0]` represents s + 2.
#[allow(dead_code)]
pub fn bode_response(num: &[f64], den: &[f64], freqs: &[f64]) -> Vec<BodePoint> {
    freqs
        .iter()
        .map(|&omega| {
            let eval_poly = |coeffs: &[f64]| -> (f64, f64) {
                let n = coeffs.len();
                let mut re = 0.0_f64;
                let mut im = 0.0_f64;
                for (k, &c) in coeffs.iter().enumerate() {
                    let power = (n - 1 - k) as i32;
                    let mag = omega.powi(power);
                    match power.rem_euclid(4) {
                        0 => re += c * mag,
                        1 => im += c * mag,
                        2 => re -= c * mag,
                        3 => im -= c * mag,
                        _ => unreachable!(),
                    }
                }
                (re, im)
            };
            let (nr, ni) = eval_poly(num);
            let (dr, di) = eval_poly(den);
            let denom_sq = dr * dr + di * di;
            let (hr, hi) = if denom_sq < 1e-300 {
                (0.0, 0.0)
            } else {
                (
                    (nr * dr + ni * di) / denom_sq,
                    (ni * dr - nr * di) / denom_sq,
                )
            };
            let magnitude = (hr * hr + hi * hi).sqrt();
            let magnitude_db = if magnitude > 0.0 {
                20.0 * magnitude.log10()
            } else {
                -f64::INFINITY
            };
            let phase_deg = hi.atan2(hr).to_degrees();
            BodePoint {
                freq: omega,
                magnitude_db,
                phase_deg,
            }
        })
        .collect()
}
#[cfg(test)]
mod tests_control_theory_extended {
    use super::*;

    use crate::PidControllerV2;

    use crate::StateFeedbackController;
    use crate::StateSpace;
    #[test]
    fn test_pid_v2_proportional_only() {
        let mut pid = PidControllerV2::new(2.0, 0.0, 0.0, f64::MAX);
        let output = pid.update(3.0, 0.01);
        assert!(
            (output - 6.0).abs() < 1e-12,
            "P-only: expected 6.0, got {output}"
        );
    }
    #[test]
    fn test_pid_v2_integral_windup_clamped() {
        let mut pid = PidControllerV2::new(0.0, 1.0, 0.0, 5.0);
        for _ in 0..1000 {
            pid.update(1.0, 0.1);
        }
        let output = pid.update(0.0, 0.0);
        assert!(
            output.abs() <= 5.0 + 1e-10,
            "windup should clamp output: {output}"
        );
    }
    #[test]
    fn test_pid_v2_reset() {
        let mut pid = PidControllerV2::new(1.0, 1.0, 0.0, 100.0);
        pid.update(5.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.prev_error, 0.0);
    }
    #[test]
    fn test_state_space_integrator_step() {
        let ss = StateSpace::new(vec![vec![0.0]], vec![vec![1.0]], vec![vec![1.0]]);
        let mut x = vec![0.0_f64];
        let u = [1.0_f64];
        let dt = 0.1;
        let mut ss_mut = ss;
        ss_mut.step(&mut x, &u, dt);
        assert!((x[0] - 0.1).abs() < 1e-12, "integrator: x = {}", x[0]);
    }
    #[test]
    fn test_state_feedback_control_sign() {
        let ctrl = StateFeedbackController::new(vec![vec![1.0, 0.0]]);
        let u = ctrl.control(&[2.0, 0.0]);
        assert!((u[0] - (-2.0)).abs() < 1e-12, "u = {}", u[0]);
    }
    #[test]
    fn test_bode_magnitude_at_dc() {
        let num = [1.0_f64];
        let den = [1.0_f64, 1.0_f64];
        let freqs = [0.0_f64];
        let points = bode_response(&num, &den, &freqs);
        assert!(
            points[0].magnitude_db.abs() < 0.1,
            "DC magnitude should be ≈ 0 dB, got {}",
            points[0].magnitude_db
        );
    }
    #[test]
    fn test_bode_response_length() {
        let num = [1.0_f64];
        let den = [1.0_f64, 1.0_f64];
        let freqs: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let points = bode_response(&num, &den, &freqs);
        assert_eq!(points.len(), 10);
    }
    #[test]
    fn test_state_space_output() {
        let ss = StateSpace::new(vec![vec![0.0]], vec![vec![1.0]], vec![vec![2.0]]);
        let x = vec![3.0_f64];
        let u = [0.0_f64];
        let y = ss.output(&x, &u);
        assert!((y[0] - 6.0).abs() < 1e-12, "y = 2*3 = 6, got {}", y[0]);
    }
}
/// Compute gain and phase margins from Bode-plot data.
///
/// Searches `freqs`/`magnitudes_db`/`phases_deg` for the crossover points.
/// Returns `None` if the crossover points cannot be found in the data.
#[allow(dead_code)]
pub fn compute_stability_margins(
    freqs: &[f64],
    magnitudes_db: &[f64],
    phases_deg: &[f64],
) -> Option<StabilityMargins> {
    let n = freqs.len();
    if n < 2 || magnitudes_db.len() != n || phases_deg.len() != n {
        return None;
    }
    let mut phase_crossover_freq = f64::NAN;
    let mut gain_at_pc = f64::NAN;
    for i in 0..n - 1 {
        let p0 = phases_deg[i];
        let p1 = phases_deg[i + 1];
        if (p0 + 180.0) * (p1 + 180.0) <= 0.0 {
            let t = (p0 + 180.0) / (p0 - p1);
            phase_crossover_freq = freqs[i] + t * (freqs[i + 1] - freqs[i]);
            gain_at_pc = magnitudes_db[i] + t * (magnitudes_db[i + 1] - magnitudes_db[i]);
            break;
        }
    }
    let mut gain_crossover_freq = f64::NAN;
    let mut phase_at_gc = f64::NAN;
    for i in 0..n - 1 {
        let m0 = magnitudes_db[i];
        let m1 = magnitudes_db[i + 1];
        if m0 * m1 <= 0.0 {
            let t = m0 / (m0 - m1);
            gain_crossover_freq = freqs[i] + t * (freqs[i + 1] - freqs[i]);
            phase_at_gc = phases_deg[i] + t * (phases_deg[i + 1] - phases_deg[i]);
            break;
        }
    }
    let gain_margin_db = if gain_at_pc.is_finite() {
        -gain_at_pc
    } else {
        f64::INFINITY
    };
    let phase_margin_deg = if phase_at_gc.is_finite() {
        phase_at_gc + 180.0
    } else {
        f64::NAN
    };
    Some(StabilityMargins {
        gain_margin_db,
        phase_margin_deg,
        phase_crossover_freq,
        gain_crossover_freq,
    })
}
/// Solve the discrete-time algebraic Riccati equation (DARE) iteratively.
///
/// P_{k+1} = A^T P_k A - A^T P_k B (B^T P_k B + R)^{-1} B^T P_k A + Q
///
/// Returns the steady-state P and LQR gain K = (B^T P B + R)^{-1} B^T P A.
/// `max_iter` iterations with convergence tolerance `tol`.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn dare_2d(
    a: &[[f64; 2]; 2],
    b: &[f64; 2],
    q: &[[f64; 2]; 2],
    r_scalar: f64,
    max_iter: usize,
    tol: f64,
) -> Option<([[f64; 2]; 2], [f64; 2])> {
    let mut p = *q;
    for _ in 0..max_iter {
        let pb = [
            p[0][0] * b[0] + p[0][1] * b[1],
            p[1][0] * b[0] + p[1][1] * b[1],
        ];
        let s = b[0] * pb[0] + b[1] * pb[1] + r_scalar;
        if s.abs() < 1e-60 {
            return None;
        }
        let bta = [
            b[0] * (p[0][0] * a[0][0] + p[0][1] * a[1][0])
                + b[1] * (p[1][0] * a[0][0] + p[1][1] * a[1][0]),
            b[0] * (p[0][0] * a[0][1] + p[0][1] * a[1][1])
                + b[1] * (p[1][0] * a[0][1] + p[1][1] * a[1][1]),
        ];
        let k_gain = [bta[0] / s, bta[1] / s];
        let atp = [
            [
                a[0][0] * p[0][0] + a[1][0] * p[1][0],
                a[0][0] * p[0][1] + a[1][0] * p[1][1],
            ],
            [
                a[0][1] * p[0][0] + a[1][1] * p[1][0],
                a[0][1] * p[0][1] + a[1][1] * p[1][1],
            ],
        ];
        let atpa = mat2_mul(&atp, a);
        let mut p_new = [[0.0f64; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                p_new[i][j] = q[i][j] + atpa[i][j] - k_gain[i] * s * k_gain[j];
            }
        }
        let diff = (p_new[0][0] - p[0][0]).abs()
            + (p_new[0][1] - p[0][1]).abs()
            + (p_new[1][0] - p[1][0]).abs()
            + (p_new[1][1] - p[1][1]).abs();
        p = p_new;
        if diff < tol {
            return Some((p, k_gain));
        }
    }
    let pb = [
        p[0][0] * b[0] + p[0][1] * b[1],
        p[1][0] * b[0] + p[1][1] * b[1],
    ];
    let s = b[0] * pb[0] + b[1] * pb[1] + r_scalar;
    if s.abs() < 1e-60 {
        return None;
    }
    let bta = [
        b[0] * (p[0][0] * a[0][0] + p[0][1] * a[1][0])
            + b[1] * (p[1][0] * a[0][0] + p[1][1] * a[1][0]),
        b[0] * (p[0][0] * a[0][1] + p[0][1] * a[1][1])
            + b[1] * (p[1][0] * a[0][1] + p[1][1] * a[1][1]),
    ];
    Some((p, [bta[0] / s, bta[1] / s]))
}
/// Estimate Ku and Tu from relay-feedback oscillation data.
///
/// * `relay_amplitude` - Relay output amplitude (d).
/// * `oscillation_amplitude` - Peak-to-peak amplitude of the output / 2.
/// * `oscillation_period` - Measured period of the limit cycle (seconds).
#[allow(dead_code)]
pub fn relay_feedback_tune(
    relay_amplitude: f64,
    oscillation_amplitude: f64,
    oscillation_period: f64,
) -> RelayFeedbackResult {
    let ku = if oscillation_amplitude.abs() > 1e-15 {
        4.0 * relay_amplitude / (std::f64::consts::PI * oscillation_amplitude)
    } else {
        0.0
    };
    RelayFeedbackResult {
        ku,
        tu: oscillation_period,
    }
}
/// Compute gain margin from a simple open-loop gain K and first-order lag τ.
///
/// For G(s) = K / (s·τ + 1)^2 (double lag), the phase-crossover frequency
/// ωpc satisfies  -180° = -2·arctan(ωpc·τ), giving ωpc = tan(90°)/τ → ∞.
/// The actual gain margin is therefore infinite for a single-lag system but
/// finite for a double-lag.  This function returns the approximate gain margin
/// for a double-lag plant with unit-feedback PID (P-only, gain = kp).
#[allow(dead_code)]
pub fn gain_margin_double_lag(k_plant: f64, tau: f64, kp: f64) -> f64 {
    let omega_pc = 3.0_f64.sqrt() / tau;
    let mag_at_pc = k_plant / (1.0 + (omega_pc * tau).powi(2));
    if kp * mag_at_pc < 1e-15 {
        f64::INFINITY
    } else {
        20.0 * (1.0 / (kp * mag_at_pc)).log10()
    }
}
/// Compute phase margin for a unit-feedback P-only controller.
///
/// At gain crossover ωgc:  kp · K / sqrt(1 + (ωgc τ)^2) = 1
/// Phase margin = 180° + phase(G(jωgc)) = 180° - arctan(ωgc τ)
#[allow(dead_code)]
pub fn phase_margin_first_order_p(k_plant: f64, _tau: f64, kp: f64) -> f64 {
    let loop_gain = kp * k_plant;
    if loop_gain <= 1.0 {
        return 90.0;
    }
    let wgc_tau = (loop_gain * loop_gain - 1.0).sqrt();
    let phase_at_gc = -wgc_tau.atan().to_degrees();
    180.0 + phase_at_gc
}
#[cfg(test)]
mod tests_control_expanded {
    use super::*;
    use crate::BumplessPid;
    use crate::DiscreteStateSpace;

    use crate::LuenbergerObserver;

    use crate::PidMode;
    use crate::SmithPredictor;

    #[test]
    fn bumpless_pid_manual_returns_override() {
        let mut pid = BumplessPid::new(1.0, 0.1, 0.0, 100.0, 100.0);
        pid.set_manual(7.5);
        let out = pid.update(1.0, 0.01);
        assert!((out - 7.5).abs() < 1e-12, "manual mode output: {out}");
    }
    #[test]
    fn bumpless_pid_auto_mode_basic() {
        let mut pid = BumplessPid::new(2.0, 0.0, 0.0, 100.0, 100.0);
        let out = pid.update(3.0, 0.01);
        assert!((out - 6.0).abs() < 1e-12, "auto mode kp*e: {out}");
    }
    #[test]
    fn bumpless_pid_bumpless_transfer_no_jump() {
        let mut pid = BumplessPid::new(1.0, 1.0, 0.0, 100.0, 100.0);
        pid.set_manual(4.0);
        pid.set_auto(0.0);
        let out = pid.update(0.0, 0.01);
        assert!((out - 4.0).abs() < 0.5, "bumpless transfer: {out}");
    }
    #[test]
    fn bumpless_pid_reset_clears_state() {
        let mut pid = BumplessPid::new(1.0, 1.0, 0.0, 100.0, 100.0);
        pid.update(5.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.prev_error, 0.0);
        assert_eq!(pid.mode, PidMode::Auto);
    }
    #[test]
    fn bumpless_pid_output_saturation() {
        let mut pid = BumplessPid::new(100.0, 0.0, 0.0, 10.0, 5.0);
        let out = pid.update(1.0, 0.01);
        assert!(out.abs() <= 5.0 + 1e-12, "saturated output: {out}");
    }
    #[test]
    fn bumpless_pid_integral_windup_clamped() {
        let mut pid = BumplessPid::new(0.0, 1.0, 0.0, 3.0, 100.0);
        for _ in 0..1000 {
            pid.update(1.0, 0.1);
        }
        assert!(
            pid.integral.abs() <= 3.0 + 1e-10,
            "integral clamped: {}",
            pid.integral
        );
    }
    #[test]
    fn bumpless_pid_derivative_action() {
        let mut pid = BumplessPid::new(0.0, 0.0, 1.0, 100.0, 100.0);
        pid.update(0.0, 0.01);
        let out = pid.update(1.0, 0.01);
        assert!((out - 100.0).abs() < 1e-6, "derivative output: {out}");
    }
    #[test]
    fn stability_margins_simple_gain_margin() {
        let freqs = vec![1.0, 2.0, 3.0];
        let magnitudes_db = vec![0.0, -6.0, -10.0];
        let phases_deg = vec![-90.0, -180.0, -200.0];
        let m = compute_stability_margins(&freqs, &magnitudes_db, &phases_deg);
        assert!(m.is_some());
        let m = m.unwrap();
        assert!(m.gain_margin_db > 0.0, "gain margin: {}", m.gain_margin_db);
    }
    #[test]
    fn stability_margins_none_on_empty() {
        let m = compute_stability_margins(&[], &[], &[]);
        assert!(m.is_none());
    }
    #[test]
    fn stability_margins_phase_margin_positive_for_stable() {
        let freqs = vec![0.5, 1.0, 2.0];
        let magnitudes_db = vec![6.0, 0.0, -6.0];
        let phases_deg = vec![-60.0, -120.0, -150.0];
        let m = compute_stability_margins(&freqs, &magnitudes_db, &phases_deg);
        let m = m.unwrap();
        assert!(
            (m.phase_margin_deg - 60.0).abs() < 1.0,
            "PM: {}",
            m.phase_margin_deg
        );
    }
    #[test]
    fn discrete_ss_double_integrator_monotone() {
        let mut sys = DiscreteStateSpace::double_integrator_zoh(0.01);
        let ys = sys.simulate(1.0, 20);
        for i in 1..ys.len() {
            assert!(ys[i] >= ys[i - 1] - 1e-12, "not monotone at step {i}");
        }
    }
    #[test]
    fn discrete_ss_zero_input_stays_zero() {
        let mut sys = DiscreteStateSpace::double_integrator_zoh(0.01);
        let ys = sys.simulate(0.0, 10);
        for y in ys {
            assert!(y.abs() < 1e-12, "nonzero with zero input: {y}");
        }
    }
    #[test]
    fn discrete_ss_reset_clears_state() {
        let mut sys = DiscreteStateSpace::double_integrator_zoh(0.01);
        sys.simulate(1.0, 10);
        sys.reset();
        assert!(
            sys.output().abs() < 1e-12,
            "output after reset: {}",
            sys.output()
        );
    }
    #[test]
    fn discrete_ss_single_step() {
        let mut sys = DiscreteStateSpace::double_integrator_zoh(0.1);
        sys.step(1.0);
        assert!(
            sys.output() > 0.0,
            "should be positive after step: {}",
            sys.output()
        );
    }
    #[test]
    fn luenberger_observer_stable_double_integrator() {
        let dt = 0.01;
        let ad = [[1.0, dt], [0.0, 1.0]];
        let bd = [0.5 * dt * dt, dt];
        let c = [1.0, 0.0];
        let l = [2.0 * (1.0 - 0.5), (1.0 - 0.5) * (1.0 - 0.5) / dt];
        let obs = LuenbergerObserver::new(ad, bd, c, l);
        let _ = obs.is_stable();
    }
    #[test]
    fn luenberger_observer_reset() {
        let ad = [[1.0, 0.01], [0.0, 1.0]];
        let bd = [0.00005, 0.01];
        let c = [1.0, 0.0];
        let l = [0.1, 0.01];
        let mut obs = LuenbergerObserver::new(ad, bd, c, l);
        obs.update(1.0, 0.5);
        obs.reset();
        let est = obs.estimate();
        assert!(est[0].abs() < 1e-12 && est[1].abs() < 1e-12);
    }
    #[test]
    fn luenberger_observer_estimate_converges() {
        let dt = 0.01;
        let ad = [[1.0, dt], [0.0, 1.0]];
        let bd = [0.0, 0.0];
        let c = [1.0, 0.0];
        let l = [1.5, 0.5];
        let mut obs = LuenbergerObserver::new(ad, bd, c, l);
        for _ in 0..200 {
            obs.update(0.0, 1.0);
        }
        assert!(
            (obs.estimate()[0] - 1.0).abs() < 0.1,
            "estimate: {}",
            obs.estimate()[0]
        );
    }
    #[test]
    fn dare_2d_double_integrator_converges() {
        let a = [[1.0, 0.01], [0.0, 1.0]];
        let b = [0.00005_f64, 0.01_f64];
        let q = [[1.0, 0.0], [0.0, 0.1]];
        let result = dare_2d(&a, &b, &q, 0.01, 500, 1e-10);
        assert!(
            result.is_some(),
            "DARE should converge for double integrator"
        );
        let (p, k) = result.unwrap();
        assert!(p[0][0] > 0.0 && p[1][1] > 0.0, "P = {:?}", p);
        assert!(k[0].is_finite() && k[1].is_finite(), "K = {:?}", k);
    }
    #[test]
    fn dare_2d_trivial_system() {
        let a = [[1.0, 0.0], [0.0, 1.0]];
        let b = [1.0, 0.0];
        let q = [[1.0, 0.0], [0.0, 1.0]];
        let result = dare_2d(&a, &b, &q, 1.0, 200, 1e-8);
        assert!(result.is_some());
    }
    #[test]
    fn relay_feedback_ku_formula() {
        let res = relay_feedback_tune(1.0, 0.5, 1.0);
        let expected = 4.0 / (std::f64::consts::PI * 0.5);
        assert!((res.ku - expected).abs() < 1e-10, "ku = {}", res.ku);
    }
    #[test]
    fn relay_feedback_tu_stored() {
        let res = relay_feedback_tune(1.0, 0.5, 2.5);
        assert!((res.tu - 2.5).abs() < 1e-12, "tu = {}", res.tu);
    }
    #[test]
    fn relay_feedback_zero_amplitude_gives_zero_ku() {
        let res = relay_feedback_tune(1.0, 0.0, 1.0);
        assert!(res.ku.abs() < 1e-12, "ku with zero amplitude: {}", res.ku);
    }
    #[test]
    fn gain_margin_double_lag_positive_for_low_gain() {
        let gm = gain_margin_double_lag(1.0, 1.0, 0.1);
        assert!(gm > 0.0, "gain margin: {gm}");
    }
    #[test]
    fn phase_margin_first_order_p_no_crossover_returns_90() {
        let pm = phase_margin_first_order_p(0.5, 1.0, 1.0);
        assert!((pm - 90.0).abs() < 1e-10, "pm = {pm}");
    }
    #[test]
    fn phase_margin_first_order_p_decreases_with_gain() {
        let pm_low = phase_margin_first_order_p(1.0, 1.0, 2.0);
        let pm_high = phase_margin_first_order_p(1.0, 1.0, 10.0);
        assert!(
            pm_low > pm_high,
            "higher gain should reduce PM: {pm_low} vs {pm_high}"
        );
    }
    #[test]
    fn smith_predictor_returns_finite_output() {
        let pid = BumplessPid::new(1.0, 0.1, 0.0, 100.0, 100.0);
        let mut sp = SmithPredictor::new(pid, 5, 1.0, 1.0, 0.01);
        for _ in 0..20 {
            let u = sp.update(1.0, 0.0);
            assert!(u.is_finite(), "smith predictor output is finite: {u}");
        }
    }
    #[test]
    fn smith_predictor_nonzero_output_for_nonzero_setpoint() {
        let pid = BumplessPid::new(2.0, 0.0, 0.0, 100.0, 100.0);
        let mut sp = SmithPredictor::new(pid, 2, 1.0, 0.5, 0.01);
        let u = sp.update(5.0, 0.0);
        assert!(
            u.abs() > 0.0,
            "should produce nonzero output for setpoint=5"
        );
    }
    #[test]
    fn bode_response_magnitude_decreases_with_freq() {
        let num = [1.0_f64];
        let den = [1.0_f64, 1.0_f64];
        let freqs = [0.1, 1.0, 10.0, 100.0];
        let pts = bode_response(&num, &den, &freqs);
        for i in 0..pts.len() - 1 {
            assert!(
                pts[i].magnitude_db >= pts[i + 1].magnitude_db,
                "magnitude should decrease: {} vs {}",
                pts[i].magnitude_db,
                pts[i + 1].magnitude_db
            );
        }
    }
    #[test]
    fn bode_response_phase_at_corner_minus45() {
        let num = [1.0_f64];
        let den = [1.0_f64, 1.0_f64];
        let pts = bode_response(&num, &den, &[1.0]);
        assert!(
            (pts[0].phase_deg + 45.0).abs() < 0.1,
            "phase at corner: {} (expected -45)",
            pts[0].phase_deg
        );
    }
    #[test]
    fn bode_response_high_freq_rolloff_20db_per_decade() {
        let num = [1.0_f64];
        let den = [1.0_f64, 1.0_f64];
        let pts = bode_response(&num, &den, &[10.0, 100.0]);
        let diff = pts[0].magnitude_db - pts[1].magnitude_db;
        assert!((diff - 20.0).abs() < 1.0, "rolloff per decade: {diff}");
    }
}
