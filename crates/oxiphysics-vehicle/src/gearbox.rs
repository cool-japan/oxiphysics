// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gearbox and transmission modeling.
//!
//! Provides a complete transmission model including gear ratios, shift logic,
//! torque converter, clutch dynamics, friction losses, and dual-clutch shifting.

// ── Gearbox ───────────────────────────────────────────────────────────────────

/// A multi-speed gearbox with gear ratios and mechanical efficiency maps.
///
/// Gear indices start at 1 (first gear).  Gear 0 is neutral (ratio = 0).
#[derive(Debug, Clone)]
pub struct Gearbox {
    /// Gear ratios for gears 1..N (index 0 = gear 1).
    pub gear_ratios: Vec<f64>,
    /// Mechanical efficiency for each gear (parallel to `gear_ratios`).
    pub efficiency: Vec<f64>,
    /// Final-drive ratio (axle ratio).
    pub final_drive: f64,
    /// Currently selected gear (1-based; 0 = neutral).
    pub current_gear: usize,
    /// Minimum engine RPM for upshift consideration.
    pub upshift_rpm: f64,
    /// Engine RPM at which to consider a downshift.
    pub downshift_rpm: f64,
    /// Shift completion time \[s\] (used by clutch model).
    pub shift_time: f64,
    /// Torque converter lock-up speed threshold \[rpm\].
    pub lockup_threshold_rpm: f64,
    /// Clutch engagement state: 0.0 = fully disengaged, 1.0 = fully engaged.
    pub clutch_engagement: f64,
    /// Internal shift progress \[0, 1\] used by dual-clutch model.
    pub dct_shift_progress: f64,
}

impl Gearbox {
    /// Create a new gearbox with the given gear ratios and final-drive ratio.
    ///
    /// All efficiencies default to 0.97.
    ///
    /// # Panics
    /// Panics if `gear_ratios` is empty.
    pub fn new(gear_ratios: Vec<f64>, final_drive: f64) -> Self {
        assert!(!gear_ratios.is_empty(), "gear_ratios must not be empty");
        let n = gear_ratios.len();
        Self {
            efficiency: vec![0.97; n],
            gear_ratios,
            final_drive,
            current_gear: 1,
            upshift_rpm: 6000.0,
            downshift_rpm: 2000.0,
            shift_time: 0.2,
            lockup_threshold_rpm: 1500.0,
            clutch_engagement: 1.0,
            dct_shift_progress: 0.0,
        }
    }

    /// Create a typical 6-speed gearbox.
    pub fn typical_6speed() -> Self {
        Self::new(vec![3.82, 2.20, 1.52, 1.15, 0.89, 0.73], 3.73)
    }

    /// Number of forward gears.
    pub fn num_gears(&self) -> usize {
        self.gear_ratios.len()
    }
}

// ── gear_ratio ────────────────────────────────────────────────────────────────

/// Return the overall drive ratio for the current gear (gear ratio × final drive).
///
/// Returns 0.0 when in neutral (gear 0) or for an invalid gear number.
///
/// # Arguments
/// * `gearbox` – The gearbox.
pub fn gear_ratio(gearbox: &Gearbox) -> f64 {
    if gearbox.current_gear == 0 || gearbox.current_gear > gearbox.gear_ratios.len() {
        return 0.0;
    }
    gearbox.gear_ratios[gearbox.current_gear - 1] * gearbox.final_drive
}

/// Return the individual internal gear ratio (without final drive) for a given gear.
///
/// Returns `None` if the gear number is out of range.
///
/// # Arguments
/// * `gearbox` – The gearbox.
/// * `gear` – Gear number (1-based).
pub fn gear_ratio_for(gearbox: &Gearbox, gear: usize) -> Option<f64> {
    if gear == 0 || gear > gearbox.gear_ratios.len() {
        return None;
    }
    Some(gearbox.gear_ratios[gear - 1])
}

// ── shift_up / shift_down ─────────────────────────────────────────────────────

/// Attempt to shift up one gear.
///
/// No-ops if already in the highest gear.
///
/// # Arguments
/// * `gearbox` – Mutable reference to the gearbox.
///
/// Returns `true` if a shift occurred.
pub fn shift_up(gearbox: &mut Gearbox) -> bool {
    if gearbox.current_gear < gearbox.gear_ratios.len() {
        gearbox.current_gear += 1;
        gearbox.clutch_engagement = 0.0; // disengage during shift
        true
    } else {
        false
    }
}

/// Attempt to shift down one gear.
///
/// No-ops if already in first gear.
///
/// # Arguments
/// * `gearbox` – Mutable reference to the gearbox.
///
/// Returns `true` if a shift occurred.
pub fn shift_down(gearbox: &mut Gearbox) -> bool {
    if gearbox.current_gear > 1 {
        gearbox.current_gear -= 1;
        gearbox.clutch_engagement = 0.0;
        true
    } else {
        false
    }
}

// ── optimal_shift_point ───────────────────────────────────────────────────────

/// Determine the optimal shift point based on current RPM and gear.
///
/// Returns `ShiftCommand::Up`, `ShiftCommand::Down`, or `ShiftCommand::Hold`.
///
/// Simple strategy: shift up at `gearbox.upshift_rpm`, shift down at `gearbox.downshift_rpm`.
///
/// # Arguments
/// * `gearbox` – The gearbox.
/// * `engine_rpm` – Current engine speed \[rpm\].
pub fn optimal_shift_point(gearbox: &Gearbox, engine_rpm: f64) -> ShiftCommand {
    if engine_rpm >= gearbox.upshift_rpm && gearbox.current_gear < gearbox.num_gears() {
        ShiftCommand::Up
    } else if engine_rpm <= gearbox.downshift_rpm && gearbox.current_gear > 1 {
        ShiftCommand::Down
    } else {
        ShiftCommand::Hold
    }
}

/// Shift command returned by `optimal_shift_point`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShiftCommand {
    /// Shift up one gear.
    Up,
    /// Shift down one gear.
    Down,
    /// Maintain current gear.
    Hold,
}

// ── torque_converter ──────────────────────────────────────────────────────────

/// Compute the torque multiplication factor of a fluid torque converter.
///
/// At low speed ratios (stall), the converter multiplies engine torque.
/// Above the coupling point, torque ratio approaches 1.0.
///
/// Uses a simplified model: `TR = k_max / (1 + (sr / sr_coupling)^n) + 1`
/// clamped to \[1.0, k_max\].
///
/// # Arguments
/// * `speed_ratio` – Output speed / input speed of the converter \[0, 1\].
/// * `stall_torque_ratio` – Maximum torque multiplication at stall (typically 1.8–2.5).
/// * `coupling_speed_ratio` – Speed ratio at which coupling begins (typically 0.85).
pub fn torque_converter(
    speed_ratio: f64,
    stall_torque_ratio: f64,
    coupling_speed_ratio: f64,
) -> f64 {
    let sr = speed_ratio.clamp(0.0, 1.0);
    let k = stall_torque_ratio - 1.0;
    let ratio = sr / coupling_speed_ratio.max(1e-6);
    let tr = k / (1.0 + ratio.powi(4)) + 1.0;
    tr.clamp(1.0, stall_torque_ratio)
}

// ── clutch_model ──────────────────────────────────────────────────────────────

/// Update the clutch engagement state for one time step.
///
/// Engagement increases toward 1.0 at rate `1/shift_time` per second.
///
/// # Arguments
/// * `gearbox` – Mutable gearbox (updates `clutch_engagement`).
/// * `dt` – Time step \[s\].
pub fn clutch_model(gearbox: &mut Gearbox, dt: f64) {
    let rate = 1.0 / gearbox.shift_time.max(1e-6);
    gearbox.clutch_engagement = (gearbox.clutch_engagement + rate * dt).min(1.0);
}

/// Compute the transmitted torque through the clutch.
///
/// `T_out = T_in * engagement^2`  (quadratic engagement curve)
///
/// # Arguments
/// * `input_torque` – Engine torque at the clutch input \[N·m\].
/// * `engagement` – Clutch engagement \[0, 1\].
pub fn clutch_transmitted_torque(input_torque: f64, engagement: f64) -> f64 {
    let e = engagement.clamp(0.0, 1.0);
    input_torque * e * e
}

// ── transmission_loss ─────────────────────────────────────────────────────────

/// Compute the friction power loss in the current gear.
///
/// `P_loss = T_in * ω_in * (1 - η)`
///
/// # Arguments
/// * `gearbox` – The gearbox.
/// * `input_torque` – Input torque \[N·m\].
/// * `input_angular_velocity` – Input shaft angular velocity \[rad/s\].
pub fn transmission_loss(gearbox: &Gearbox, input_torque: f64, input_angular_velocity: f64) -> f64 {
    let eta = if gearbox.current_gear == 0 || gearbox.current_gear > gearbox.efficiency.len() {
        1.0
    } else {
        gearbox.efficiency[gearbox.current_gear - 1]
    };
    input_torque * input_angular_velocity * (1.0 - eta)
}

// ── dual_clutch_shift ─────────────────────────────────────────────────────────

/// Simulate a Dual Clutch Transmission (DCT) shift step.
///
/// In a DCT, the outgoing clutch disengages while the incoming clutch engages
/// simultaneously, providing near-instantaneous torque transfer.
///
/// This function advances the shift progress by `dt / shift_time` and applies
/// the resulting gear change when progress reaches 1.0.
///
/// # Arguments
/// * `gearbox` – Mutable gearbox.
/// * `target_gear` – The gear to shift to (1-based).
/// * `dt` – Time step \[s\].
///
/// Returns `true` when the shift is complete.
pub fn dual_clutch_shift(gearbox: &mut Gearbox, target_gear: usize, dt: f64) -> bool {
    if gearbox.current_gear == target_gear {
        gearbox.dct_shift_progress = 0.0;
        return true;
    }

    let rate = dt / gearbox.shift_time.max(1e-6);
    gearbox.dct_shift_progress = (gearbox.dct_shift_progress + rate).min(1.0);

    // Blend: outgoing clutch = 1 - progress, incoming = progress
    // We model the effective engagement as min(incoming, outgoing) transition
    gearbox.clutch_engagement = gearbox.dct_shift_progress;

    if gearbox.dct_shift_progress >= 1.0 {
        if target_gear >= 1 && target_gear <= gearbox.gear_ratios.len() {
            gearbox.current_gear = target_gear;
        }
        gearbox.clutch_engagement = 1.0;
        gearbox.dct_shift_progress = 0.0;
        true
    } else {
        false
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn six_speed() -> Gearbox {
        Gearbox::typical_6speed()
    }

    // ── gear_ratio ──────────────────────────────────────────────────────────

    #[test]
    fn gear_ratio_first_gear() {
        let gb = six_speed();
        let r = gear_ratio(&gb);
        let expected = 3.82 * 3.73;
        assert!((r - expected).abs() < 1e-6, "r={r}, expected={expected}");
    }

    #[test]
    fn gear_ratio_neutral_is_zero() {
        let mut gb = six_speed();
        gb.current_gear = 0;
        assert_eq!(gear_ratio(&gb), 0.0);
    }

    #[test]
    fn gear_ratio_decreases_with_higher_gear() {
        let mut gb = six_speed();
        gb.current_gear = 1;
        let r1 = gear_ratio(&gb);
        gb.current_gear = 6;
        let r6 = gear_ratio(&gb);
        assert!(r6 < r1, "r6={r6} should be < r1={r1}");
    }

    #[test]
    fn gear_ratio_for_valid_gear() {
        let gb = six_speed();
        let r = gear_ratio_for(&gb, 1);
        assert!(r.is_some());
        assert!((r.unwrap() - 3.82).abs() < 1e-6);
    }

    #[test]
    fn gear_ratio_for_invalid_gear() {
        let gb = six_speed();
        assert!(gear_ratio_for(&gb, 0).is_none());
        assert!(gear_ratio_for(&gb, 99).is_none());
    }

    // ── shift_up / shift_down ───────────────────────────────────────────────

    #[test]
    fn shift_up_increments_gear() {
        let mut gb = six_speed();
        gb.current_gear = 3;
        let shifted = shift_up(&mut gb);
        assert!(shifted);
        assert_eq!(gb.current_gear, 4);
    }

    #[test]
    fn shift_up_at_top_gear_does_nothing() {
        let mut gb = six_speed();
        gb.current_gear = 6;
        let shifted = shift_up(&mut gb);
        assert!(!shifted);
        assert_eq!(gb.current_gear, 6);
    }

    #[test]
    fn shift_down_decrements_gear() {
        let mut gb = six_speed();
        gb.current_gear = 4;
        let shifted = shift_down(&mut gb);
        assert!(shifted);
        assert_eq!(gb.current_gear, 3);
    }

    #[test]
    fn shift_down_at_first_gear_does_nothing() {
        let mut gb = six_speed();
        gb.current_gear = 1;
        let shifted = shift_down(&mut gb);
        assert!(!shifted);
        assert_eq!(gb.current_gear, 1);
    }

    #[test]
    fn shift_up_disengages_clutch() {
        let mut gb = six_speed();
        gb.current_gear = 2;
        gb.clutch_engagement = 1.0;
        shift_up(&mut gb);
        assert_eq!(gb.clutch_engagement, 0.0);
    }

    #[test]
    fn shift_down_disengages_clutch() {
        let mut gb = six_speed();
        gb.current_gear = 3;
        gb.clutch_engagement = 1.0;
        shift_down(&mut gb);
        assert_eq!(gb.clutch_engagement, 0.0);
    }

    // ── optimal_shift_point ────────────────────────────────────────────────

    #[test]
    fn optimal_shift_upshift_above_rpm() {
        let mut gb = six_speed();
        gb.current_gear = 2;
        let cmd = optimal_shift_point(&gb, 7000.0);
        assert_eq!(cmd, ShiftCommand::Up);
    }

    #[test]
    fn optimal_shift_downshift_below_rpm() {
        let mut gb = six_speed();
        gb.current_gear = 3;
        let cmd = optimal_shift_point(&gb, 1500.0);
        assert_eq!(cmd, ShiftCommand::Down);
    }

    #[test]
    fn optimal_shift_hold_in_range() {
        let gb = six_speed();
        let cmd = optimal_shift_point(&gb, 3500.0);
        assert_eq!(cmd, ShiftCommand::Hold);
    }

    #[test]
    fn optimal_shift_no_upshift_in_top_gear() {
        let mut gb = six_speed();
        gb.current_gear = 6;
        let cmd = optimal_shift_point(&gb, 8000.0);
        // top gear → cannot upshift
        assert_eq!(cmd, ShiftCommand::Hold);
    }

    #[test]
    fn optimal_shift_no_downshift_in_first_gear() {
        let mut gb = six_speed();
        gb.current_gear = 1;
        let cmd = optimal_shift_point(&gb, 500.0);
        assert_eq!(cmd, ShiftCommand::Hold);
    }

    // ── torque_converter ───────────────────────────────────────────────────

    #[test]
    fn torque_converter_stall_ratio() {
        // At speed_ratio = 0, should return stall_torque_ratio
        let tr = torque_converter(0.0, 2.0, 0.85);
        assert!((tr - 2.0).abs() < 0.05, "tr={tr}");
    }

    #[test]
    fn torque_converter_high_speed_ratio() {
        // At speed_ratio = 1, torque ratio should be >= 1 and <= stall_ratio
        let tr = torque_converter(1.0, 2.0, 0.85);
        assert!((1.0..=2.0).contains(&tr), "tr={tr}");
    }

    #[test]
    fn torque_converter_decreases_with_speed_ratio() {
        let tr0 = torque_converter(0.0, 2.0, 0.85);
        let tr1 = torque_converter(0.9, 2.0, 0.85);
        assert!(tr0 > tr1, "tr0={tr0} should be > tr1={tr1}");
    }

    #[test]
    fn torque_converter_minimum_one() {
        for sr in [0.0, 0.5, 1.0, 1.5] {
            let tr = torque_converter(sr, 2.0, 0.85);
            assert!(tr >= 1.0, "tr={tr} at sr={sr}");
        }
    }

    // ── clutch_model ───────────────────────────────────────────────────────

    #[test]
    fn clutch_model_engagement_increases() {
        let mut gb = six_speed();
        gb.clutch_engagement = 0.0;
        gb.shift_time = 0.2;
        clutch_model(&mut gb, 0.1);
        assert!(gb.clutch_engagement > 0.0 && gb.clutch_engagement < 1.0);
    }

    #[test]
    fn clutch_model_capped_at_one() {
        let mut gb = six_speed();
        gb.clutch_engagement = 0.9;
        gb.shift_time = 0.01;
        clutch_model(&mut gb, 1.0);
        assert_eq!(gb.clutch_engagement, 1.0);
    }

    #[test]
    fn clutch_transmitted_torque_zero_engagement() {
        assert_eq!(clutch_transmitted_torque(200.0, 0.0), 0.0);
    }

    #[test]
    fn clutch_transmitted_torque_full_engagement() {
        assert!((clutch_transmitted_torque(200.0, 1.0) - 200.0).abs() < 1e-10);
    }

    #[test]
    fn clutch_transmitted_torque_partial() {
        let t = clutch_transmitted_torque(100.0, 0.5);
        assert!((t - 25.0).abs() < 1e-10, "t={t}"); // 100 * 0.25
    }

    // ── transmission_loss ──────────────────────────────────────────────────

    #[test]
    fn transmission_loss_neutral_zero() {
        let mut gb = six_speed();
        gb.current_gear = 0;
        let loss = transmission_loss(&gb, 200.0, 100.0);
        assert_eq!(loss, 0.0);
    }

    #[test]
    fn transmission_loss_positive_in_gear() {
        let gb = six_speed(); // gear 1, η=0.97
        let loss = transmission_loss(&gb, 200.0, 100.0);
        let expected = 200.0 * 100.0 * 0.03;
        assert!((loss - expected).abs() < 1e-6, "loss={loss}");
    }

    #[test]
    fn transmission_loss_scales_with_torque() {
        let gb = six_speed();
        let l1 = transmission_loss(&gb, 100.0, 100.0);
        let l2 = transmission_loss(&gb, 200.0, 100.0);
        assert!((l2 - 2.0 * l1).abs() < 1e-10, "l1={l1}, l2={l2}");
    }

    // ── dual_clutch_shift ──────────────────────────────────────────────────

    #[test]
    fn dct_shift_same_gear_instant() {
        let mut gb = six_speed();
        let done = dual_clutch_shift(&mut gb, 1, 1.0);
        assert!(done);
        assert_eq!(gb.current_gear, 1);
    }

    #[test]
    fn dct_shift_completes_over_time() {
        let mut gb = six_speed();
        gb.current_gear = 2;
        gb.shift_time = 0.2;
        let mut done = false;
        for _ in 0..20 {
            done = dual_clutch_shift(&mut gb, 3, 0.02);
        }
        assert!(done);
        assert_eq!(gb.current_gear, 3);
    }

    #[test]
    fn dct_shift_engagement_increases() {
        let mut gb = six_speed();
        gb.current_gear = 1;
        gb.shift_time = 0.5;
        gb.clutch_engagement = 0.0;
        dual_clutch_shift(&mut gb, 2, 0.1);
        assert!(gb.clutch_engagement > 0.0);
    }

    #[test]
    fn dct_shift_without_torque_interruption() {
        let mut gb = six_speed();
        gb.current_gear = 3;
        gb.shift_time = 0.1;
        // Advance in small steps — engine torque path should remain continuous
        let mut prev_engagement = 0.0_f64;
        gb.clutch_engagement = 0.0;
        for _ in 0..5 {
            dual_clutch_shift(&mut gb, 4, 0.02);
            assert!(
                gb.clutch_engagement >= prev_engagement,
                "engagement should be non-decreasing: {} < {}",
                gb.clutch_engagement,
                prev_engagement
            );
            prev_engagement = gb.clutch_engagement;
        }
    }

    #[test]
    fn gearbox_num_gears() {
        let gb = six_speed();
        assert_eq!(gb.num_gears(), 6);
    }
}
