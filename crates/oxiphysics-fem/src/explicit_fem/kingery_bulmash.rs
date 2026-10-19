// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Kingery–Bulmash blast-wave parameters of a TNT surface burst (the ConWep model).
//!
//! Kingery and Bulmash (BRL report ARBRL-TR-02555, 1984) tabulated the
//! positive-phase parameters of the air blast from a hemispherical TNT charge
//! on the ground as functions of the Hopkinson–Cranz scaled distance
//! `Z = R / W^(1/3)` (`R` in m, `W` in kg of TNT).  Swisdak ("Simplified
//! Kingery Airblast Calculations", Proc. 26th DoD Explosives Safety Seminar,
//! 1994, Table 1, metric units) condensed the curves into piecewise polynomials
//! of `ln Z`, accurate to about 1 % on average:
//!
//! ```text
//! parameter = exp( c0 + c1·ln Z + c2·(ln Z)² + … )
//! ```
//!
//! This module implements four of them:
//!
//! | parameter                    | unit               | scaling        |
//! |------------------------------|--------------------|----------------|
//! | peak incident overpressure `P_s` | kPa            | none           |
//! | time of arrival `t_a`        | ms / kg^(1/3)      | × `W^(1/3)`    |
//! | positive-phase duration `t_d`| ms / kg^(1/3)      | × `W^(1/3)`    |
//! | positive incident impulse `i_s` | kPa·ms / kg^(1/3) | × `W^(1/3)`  |
//!
//! The coefficients were cross-checked against the independent English-unit
//! fits of the same table (agreement within about 1 % for `0.25 ≤ Z ≤ 30`) and
//! reproduce the pressure, impulse and arrival times quoted for 1 t, 10 t and
//! 100 t of TNT at 50 m (UN IATG 01.80, Table 5).
//!
//! # Range handling
//!
//! - The fits are valid for `0.2 ≤ Z ≤ 40` (pressure and impulse up to 198.5 and
//!   158.7).  Below `Z = 0.2` the scaled distance is clamped, so the load stays
//!   bounded close to the charge.
//! - Above the upper end of the last polynomial the curve continues as a power
//!   law with the local logarithmic slope of the fit, so `P_s → 0` as `Z → ∞`.
//! - Neighbouring polynomial pieces disagree by up to 2.4 % at their junctions.
//!   They are blended with a C¹ smooth-step over ±5 % of the junction's scaled
//!   distance, so every parameter is a continuous function of the range and the
//!   overpressure decreases monotonically.
//!
//! # Waveform
//!
//! The pressure at a point is the Friedlander pulse
//! `p(t) = P_s (1 − τ) e^(−b τ)`, `τ = (t − t_a)/t_d`, for `t_a ≤ t ≤ t_a + t_d`
//! and zero otherwise.  The decay coefficient `b` is chosen so that the pulse
//! carries exactly the tabulated impulse:
//! `i_s = P_s t_d (b − 1 + e^(−b)) / b²`.

/// Lower end of the scaled-distance range (m/kg^(1/3)) of the fits; closer
/// points are evaluated at this distance.
pub const MIN_SCALED_DISTANCE: f64 = 0.2;

/// Upper end of the scaled-distance range (m/kg^(1/3)) of the arrival-time and
/// duration fits; farther points continue the last polynomial as a power law.
pub const MAX_SCALED_DISTANCE: f64 = 40.0;

/// Half-width, in `ln Z`, of the blending band around a junction of two pieces.
const BLEND_HALF_WIDTH: f64 = 0.05;

/// One polynomial piece `Y(u) = Σ cᵢ uⁱ`, `u = ln Z`, valid up to the scaled
/// distance `upper` (the lower end is the upper end of the previous piece).
struct Piece {
    upper: f64,
    coefficients: &'static [f64],
}

impl Piece {
    fn value(&self, u: f64) -> f64 {
        self.coefficients
            .iter()
            .rev()
            .fold(0.0, |acc, &c| acc * u + c)
    }

    /// Derivative `dY/du`.
    fn slope(&self, u: f64) -> f64 {
        self.coefficients
            .iter()
            .enumerate()
            .skip(1)
            .rev()
            .fold(0.0, |acc, (i, &c)| acc * u + i as f64 * c)
    }
}

/// Time of arrival, ms per kg^(1/3).
const ARRIVAL_TIME: &[Piece] = &[
    Piece {
        upper: 1.50,
        coefficients: &[-0.7604, 1.8058, 0.1257, -0.0437, -0.0310, -0.00669],
    },
    Piece {
        upper: 40.0,
        coefficients: &[-0.7137, 1.5732, 0.5561, -0.4213, 0.1054, -0.00929],
    },
];

/// Peak incident overpressure, kPa.
const PEAK_OVERPRESSURE: &[Piece] = &[
    Piece {
        upper: 2.9,
        coefficients: &[7.2106, -2.1069, -0.3229, 0.1117, 0.0685],
    },
    Piece {
        upper: 23.8,
        coefficients: &[7.5938, -3.0523, 0.40977, 0.0261, -0.01267],
    },
    Piece {
        upper: 198.5,
        coefficients: &[6.0536, -1.4066],
    },
];

/// Positive-phase duration, ms per kg^(1/3).
const POSITIVE_DURATION: &[Piece] = &[
    Piece {
        upper: 1.02,
        coefficients: &[0.5426, 3.2299, -1.5931, -5.9667, -4.0815, -0.9149],
    },
    Piece {
        upper: 2.80,
        coefficients: &[0.5440, 2.7082, -9.7354, 14.3425, -9.7791, 2.8535],
    },
    Piece {
        upper: 40.0,
        coefficients: &[-2.4608, 7.1639, -5.6215, 2.2711, -0.44994, 0.03486],
    },
];

/// Positive incident impulse, kPa·ms per kg^(1/3).
const POSITIVE_IMPULSE: &[Piece] = &[
    Piece {
        upper: 0.96,
        coefficients: &[5.522, 1.117, 0.6, -0.292, -0.087],
    },
    Piece {
        upper: 2.38,
        coefficients: &[5.465, -0.308, -1.464, 1.362, -0.432],
    },
    Piece {
        upper: 33.7,
        coefficients: &[5.2749, -0.4677, -0.2499, 0.0588, -0.00554],
    },
    Piece {
        upper: 158.7,
        coefficients: &[5.9825, -1.062],
    },
];

/// C¹ smooth-step from 0 to 1 over `s ∈ [0, 1]`.
fn smooth_step(s: f64) -> f64 {
    let s = s.clamp(0.0, 1.0);
    s * s * (3.0 - 2.0 * s)
}

/// Evaluate a piecewise fit at the scaled distance `z` (see the module docs for
/// the clamping, blending and extension rules).
fn evaluate(pieces: &[Piece], z: f64) -> f64 {
    let z = z.max(MIN_SCALED_DISTANCE);
    let u = z.ln();
    for (k, piece) in pieces.iter().enumerate() {
        let Some(next) = pieces.get(k + 1) else {
            // Last piece: continue as a power law beyond its upper limit.
            let u_hi = piece.upper.ln();
            let y = if u > u_hi {
                piece.value(u_hi) + piece.slope(u_hi) * (u - u_hi)
            } else {
                piece.value(u)
            };
            return y.exp();
        };
        let u_junction = piece.upper.ln();
        if u <= u_junction - BLEND_HALF_WIDTH {
            return piece.value(u).exp();
        }
        if u < u_junction + BLEND_HALF_WIDTH {
            let weight =
                smooth_step((u - (u_junction - BLEND_HALF_WIDTH)) / (2.0 * BLEND_HALF_WIDTH));
            return ((1.0 - weight) * piece.value(u) + weight * next.value(u)).exp();
        }
    }
    f64::NAN
}

/// Kingery–Bulmash surface-burst parameters of a 1 kg TNT charge at one scaled
/// distance (see the module docs for the units).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaledSurfaceBurst {
    /// Peak incident overpressure (kPa).
    pub peak_overpressure_kpa: f64,
    /// Shock arrival time (ms per kg^(1/3)).
    pub arrival_time_ms: f64,
    /// Positive-phase duration (ms per kg^(1/3)).
    pub positive_duration_ms: f64,
    /// Positive incident impulse (kPa·ms per kg^(1/3)).
    pub positive_impulse_kpa_ms: f64,
}

impl ScaledSurfaceBurst {
    /// Parameters at the scaled distance `z` (m/kg^(1/3)); values below
    /// [`MIN_SCALED_DISTANCE`] are evaluated at that distance.
    pub fn at(z: f64) -> Self {
        Self {
            peak_overpressure_kpa: evaluate(PEAK_OVERPRESSURE, z),
            arrival_time_ms: evaluate(ARRIVAL_TIME, z),
            positive_duration_ms: evaluate(POSITIVE_DURATION, z),
            positive_impulse_kpa_ms: evaluate(POSITIVE_IMPULSE, z),
        }
    }
}

/// Dimensionless impulse `(b − 1 + e^(−b)) / b²` of the Friedlander pulse
/// `(1 − τ) e^(−bτ)` over `τ ∈ [0, 1]` (it equals `½` for `b = 0` and decreases
/// monotonically with `b`).
pub fn friedlander_impulse_factor(b: f64) -> f64 {
    if b.abs() < 1.0e-2 {
        // Taylor series: the closed form cancels catastrophically near b = 0.
        0.5 - b / 6.0 + b * b / 24.0 - b * b * b / 120.0 + b * b * b * b / 720.0
    } else {
        (b - 1.0 + (-b).exp()) / (b * b)
    }
}

/// Decay coefficient `b ≥ 0` of the Friedlander pulse whose impulse is the
/// fraction `ratio = i_s / (P_s t_d)` of the peak-pressure rectangle.
///
/// The pulse impulse never exceeds the triangular value `½`, so `ratio ≥ ½`
/// yields `b = 0`.  A vanishing `ratio` (infinitely sharp decay) returns
/// `f64::INFINITY`, and a non-finite `ratio` returns NaN.
pub fn friedlander_decay_coefficient(ratio: f64) -> f64 {
    if !ratio.is_finite() {
        return f64::NAN;
    }
    if ratio >= 0.5 {
        return 0.0;
    }
    if ratio <= 0.0 {
        return f64::INFINITY;
    }
    // g(b) < 1/b, so g(1/ratio) < ratio < g(0) = 1/2 brackets the root.
    let mut low = 0.0_f64;
    let mut high = 1.0 / ratio;
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if mid <= low || mid >= high {
            break;
        }
        if friedlander_impulse_factor(mid) > ratio {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

/// Positive-phase blast wave at a point: arrival, duration, peak overpressure
/// and the Friedlander decay coefficient that reproduces the positive impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlastWave {
    /// Peak incident overpressure `P_s` (Pa).
    pub peak_overpressure: f64,
    /// Shock arrival time `t_a` (s).
    pub arrival_time: f64,
    /// Positive-phase duration `t_d` (s).
    pub positive_duration: f64,
    /// Positive incident impulse `i_s` (Pa·s).
    pub positive_impulse: f64,
    /// Friedlander decay coefficient `b` matching `positive_impulse`.
    pub decay_coeff: f64,
}

impl BlastWave {
    /// Blast wave of a TNT-equivalent surface burst of `charge_mass` kg at the
    /// distance `standoff` m, from the Kingery–Bulmash fits of this module.
    ///
    /// Returns `None` for a non-positive or non-finite charge mass, or a
    /// negative or non-finite standoff.  A zero standoff is evaluated at the
    /// lower end of the fit range.
    pub fn surface_burst(standoff: f64, charge_mass: f64) -> Option<Self> {
        if !standoff.is_finite() || !charge_mass.is_finite() || standoff < 0.0 || charge_mass <= 0.0
        {
            return None;
        }
        let cbrt_mass = charge_mass.cbrt();
        let scaled = ScaledSurfaceBurst::at(standoff / cbrt_mass);
        let peak_overpressure = scaled.peak_overpressure_kpa * 1.0e3;
        let positive_duration = scaled.positive_duration_ms * 1.0e-3 * cbrt_mass;
        // kPa·ms = Pa·s
        let positive_impulse = scaled.positive_impulse_kpa_ms * cbrt_mass;
        let ratio = positive_impulse / (peak_overpressure * positive_duration);
        Some(Self {
            peak_overpressure,
            arrival_time: scaled.arrival_time_ms * 1.0e-3 * cbrt_mass,
            positive_duration,
            positive_impulse,
            decay_coeff: friedlander_decay_coefficient(ratio),
        })
    }

    /// Overpressure (Pa) at time `t` (s after detonation): zero before the
    /// shock arrives and after the positive phase, the Friedlander pulse in
    /// between.
    pub fn overpressure_at(&self, t: f64) -> f64 {
        let tau = (t - self.arrival_time) / self.positive_duration;
        if !(0.0..=1.0).contains(&tau) {
            return 0.0;
        }
        self.peak_overpressure * (1.0 - tau) * (-self.decay_coeff * tau).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_rel_eq(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let scale = expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            (actual - expected).abs() <= tolerance * scale,
            "{what}: {actual} vs {expected} (relative tolerance {tolerance})"
        );
    }

    /// Hemispherical TNT surface burst at 50 m as tabulated in the UN IATG 01.80
    /// (3rd ed., 2021) Table 5: `(mass kg, P_s bar, i_s bar·ms, t_a ms)`.
    const IATG_TABLE_5: [(f64, f64, f64, f64); 3] = [
        (1_000.0, 0.43, 5.9, 82.4),
        (10_000.0, 2.02, 25.2, 48.1),
        (100_000.0, 11.5, 106.0, 24.8),
    ];

    #[test]
    fn matches_published_surface_burst_values() {
        for (mass, p_bar, impulse_bar_ms, arrival_ms) in IATG_TABLE_5 {
            let wave = BlastWave::surface_burst(50.0, mass).expect("valid inputs");
            assert_rel_eq(wave.peak_overpressure, p_bar * 1.0e5, 0.015, "P_s");
            // 1 bar·ms = 100 Pa·s
            assert_rel_eq(wave.positive_impulse, impulse_bar_ms * 100.0, 0.015, "i_s");
            assert_rel_eq(wave.arrival_time, arrival_ms * 1.0e-3, 0.015, "t_a");
        }
    }

    #[test]
    fn hopkinson_cranz_cube_root_scaling() {
        // Lengths x2 and mass x8 keep Z; times and impulse double, pressure stays.
        let small = BlastWave::surface_burst(7.0, 3.0).expect("valid");
        let large = BlastWave::surface_burst(14.0, 24.0).expect("valid");
        assert_rel_eq(
            large.peak_overpressure,
            small.peak_overpressure,
            1e-12,
            "P_s",
        );
        assert_rel_eq(large.arrival_time, 2.0 * small.arrival_time, 1e-12, "t_a");
        assert_rel_eq(
            large.positive_duration,
            2.0 * small.positive_duration,
            1e-12,
            "t_d",
        );
        assert_rel_eq(
            large.positive_impulse,
            2.0 * small.positive_impulse,
            1e-12,
            "i_s",
        );
        assert_rel_eq(large.decay_coeff, small.decay_coeff, 1e-9, "b");
    }

    #[test]
    fn peak_overpressure_is_continuous_and_strictly_decreasing() {
        let mut z = MIN_SCALED_DISTANCE;
        let mut previous = ScaledSurfaceBurst::at(z).peak_overpressure_kpa;
        while z < 5_000.0 {
            z *= 1.001;
            let p = ScaledSurfaceBurst::at(z).peak_overpressure_kpa;
            assert!(
                p < previous,
                "P_s not decreasing at Z = {z}: {p} >= {previous}"
            );
            assert!(
                (previous - p) / previous < 0.005,
                "P_s jumps at Z = {z}: {previous} -> {p}"
            );
            previous = p;
        }
    }

    #[test]
    fn all_parameters_are_continuous_across_the_polynomial_junctions() {
        for junction in [0.96, 1.02, 1.5, 2.38, 2.8, 2.9, 23.8, 33.7] {
            let below = ScaledSurfaceBurst::at(junction * (1.0 - 1e-9));
            let above = ScaledSurfaceBurst::at(junction * (1.0 + 1e-9));
            assert_rel_eq(
                below.peak_overpressure_kpa,
                above.peak_overpressure_kpa,
                1e-6,
                "P_s",
            );
            assert_rel_eq(below.arrival_time_ms, above.arrival_time_ms, 1e-6, "t_a");
            assert_rel_eq(
                below.positive_duration_ms,
                above.positive_duration_ms,
                1e-6,
                "t_d",
            );
            assert_rel_eq(
                below.positive_impulse_kpa_ms,
                above.positive_impulse_kpa_ms,
                1e-6,
                "i_s",
            );
        }
        // Dense scan: no parameter ever jumps by more than a few percent between
        // neighbouring samples 0.1 % apart.
        let mut z = MIN_SCALED_DISTANCE;
        let mut previous = ScaledSurfaceBurst::at(z);
        while z < 100.0 {
            z *= 1.001;
            let current = ScaledSurfaceBurst::at(z);
            for (a, b, what) in [
                (previous.arrival_time_ms, current.arrival_time_ms, "t_a"),
                (
                    previous.positive_duration_ms,
                    current.positive_duration_ms,
                    "t_d",
                ),
                (
                    previous.positive_impulse_kpa_ms,
                    current.positive_impulse_kpa_ms,
                    "i_s",
                ),
            ] {
                assert!(
                    (a - b).abs() / a < 0.02,
                    "{what} jumps at Z = {z}: {a} -> {b}"
                );
            }
            previous = current;
        }
    }

    #[test]
    fn overpressure_vanishes_far_from_the_charge() {
        let at_40 = ScaledSurfaceBurst::at(MAX_SCALED_DISTANCE).peak_overpressure_kpa;
        assert!(at_40 < 5.0, "P_s(Z = 40) = {at_40} kPa");
        assert!(ScaledSurfaceBurst::at(1.0e4).peak_overpressure_kpa < 0.01);
        let far = BlastWave::surface_burst(1.0e9, 1.0).expect("valid");
        assert!(
            far.peak_overpressure < 1.0e-3,
            "{} Pa",
            far.peak_overpressure
        );
        assert!(far.peak_overpressure >= 0.0);
        // The old fit tended to 0.84 MPa (about 8 atm) at infinite standoff.
        assert!(ScaledSurfaceBurst::at(1.0e12).peak_overpressure_kpa < 1.0e-6);
    }

    #[test]
    fn scaled_distance_below_the_fit_range_is_clamped() {
        let at_limit = ScaledSurfaceBurst::at(MIN_SCALED_DISTANCE);
        assert!(
            at_limit.peak_overpressure_kpa.is_finite() && at_limit.peak_overpressure_kpa > 1.0e3
        );
        for z in [0.0, 1.0e-6, 0.05, 0.1999] {
            assert_eq!(ScaledSurfaceBurst::at(z), at_limit, "Z = {z}");
        }
        let wave = BlastWave::surface_burst(0.0, 2.0).expect("zero standoff is clamped");
        assert!(wave.peak_overpressure.is_finite() && wave.peak_overpressure > 0.0);
        assert!(wave.decay_coeff.is_finite());
    }

    #[test]
    fn fits_continue_smoothly_beyond_their_upper_limits() {
        // Arrival time keeps growing with the slope reached at Z = 40 (C¹).
        let h = 1.0e-4;
        let left = (ScaledSurfaceBurst::at(40.0).arrival_time_ms
            - ScaledSurfaceBurst::at(40.0 - h).arrival_time_ms)
            / h;
        let right = (ScaledSurfaceBurst::at(40.0 + h).arrival_time_ms
            - ScaledSurfaceBurst::at(40.0).arrival_time_ms)
            / h;
        assert_rel_eq(right, left, 2.0e-3, "dt_a/dZ across Z = 40");
        let mut previous = ScaledSurfaceBurst::at(40.0);
        for z in [41.0, 60.0, 100.0, 500.0, 5_000.0] {
            let current = ScaledSurfaceBurst::at(z);
            assert!(current.arrival_time_ms > previous.arrival_time_ms);
            assert!(current.positive_duration_ms > previous.positive_duration_ms);
            assert!(current.positive_impulse_kpa_ms < previous.positive_impulse_kpa_ms);
            previous = current;
        }
        // The shock never outruns the speed of sound by much: dZ/dt_a at Z = 5000
        // is the shock speed in (m per kg^(1/3)) / (ms per kg^(1/3)) = km/s.
        let speed = 5_000.0 / ScaledSurfaceBurst::at(5_000.0).arrival_time_ms * 1.0e3;
        assert!(
            speed > 250.0 && speed < 450.0,
            "far-field shock speed {speed} m/s"
        );
    }

    #[test]
    fn friedlander_impulse_factor_properties() {
        assert!((friedlander_impulse_factor(0.0) - 0.5).abs() < 1e-15);
        assert!((friedlander_impulse_factor(1.0) - (-1.0_f64).exp()).abs() < 1e-15);
        // Continuous across the switch to the Taylor series.
        let below = friedlander_impulse_factor(0.01 - 1e-12);
        let above = friedlander_impulse_factor(0.01 + 1e-12);
        assert!((below - above).abs() < 1e-10);
        // Strictly decreasing, bounded by 1/b.
        let mut previous = 0.5;
        for k in 1..400 {
            let b = 0.05 * f64::from(k);
            let g = friedlander_impulse_factor(b);
            assert!(g < previous && g > 0.0 && g < 1.0 / b, "b = {b}");
            previous = g;
        }
    }

    #[test]
    fn decay_coefficient_inverts_the_impulse_factor() {
        for ratio in [0.499, 0.45, 0.32, 0.2, 0.1, 0.05, 0.01, 0.001] {
            let b = friedlander_decay_coefficient(ratio);
            assert!(b > 0.0 && b.is_finite(), "ratio {ratio}: b = {b}");
            assert_rel_eq(friedlander_impulse_factor(b), ratio, 1e-12, "g(b)");
        }
        // g(2) = (1 + e^-2) / 4
        let g2 = (1.0 + (-2.0_f64).exp()) / 4.0;
        assert!((friedlander_decay_coefficient(g2) - 2.0).abs() < 1e-10);
        // At or above the triangular limit the pulse cannot decay any slower.
        assert_eq!(friedlander_decay_coefficient(0.5), 0.0);
        assert_eq!(friedlander_decay_coefficient(0.8), 0.0);
        assert_eq!(friedlander_decay_coefficient(0.0), f64::INFINITY);
        assert_eq!(friedlander_decay_coefficient(-1.0), f64::INFINITY);
        assert!(friedlander_decay_coefficient(f64::NAN).is_nan());
    }

    #[test]
    fn tabulated_parameters_admit_a_physical_friedlander_pulse() {
        // i_s/(P_s t_d) must lie in (0, 1/2) over the whole range, otherwise no
        // monotonically decaying Friedlander pulse can carry the impulse.
        let mut z = MIN_SCALED_DISTANCE;
        while z < 5_000.0 {
            let s = ScaledSurfaceBurst::at(z);
            let ratio =
                s.positive_impulse_kpa_ms / (s.peak_overpressure_kpa * s.positive_duration_ms);
            assert!(ratio > 0.0 && ratio < 0.5, "Z = {z}: impulse ratio {ratio}");
            z *= 1.05;
        }
    }

    #[test]
    fn waveform_timing_and_zero_outside_the_positive_phase() {
        let wave = BlastWave::surface_burst(10.0, 100.0).expect("valid");
        assert!(wave.arrival_time > 0.0 && wave.positive_duration > 0.0);
        assert_eq!(wave.overpressure_at(0.0), 0.0);
        assert_eq!(wave.overpressure_at(-1.0), 0.0);
        assert_eq!(wave.overpressure_at(0.999 * wave.arrival_time), 0.0);
        // The front is a jump to the peak overpressure.
        assert_eq!(
            wave.overpressure_at(wave.arrival_time),
            wave.peak_overpressure
        );
        let end = wave.arrival_time + wave.positive_duration;
        assert_eq!(wave.overpressure_at(end * 1.000_001), 0.0);
        assert_eq!(wave.overpressure_at(f64::NAN), 0.0);
        // Strictly decaying inside the positive phase.
        let mut previous = wave.peak_overpressure;
        for k in 1..100 {
            let t = wave.arrival_time + wave.positive_duration * f64::from(k) / 100.0;
            let p = wave.overpressure_at(t);
            assert!(p > 0.0 && p < previous, "k = {k}: {p} vs {previous}");
            previous = p;
        }
    }

    #[test]
    fn waveform_integrates_to_the_tabulated_impulse() {
        for (standoff, mass) in [
            (0.4, 1.0),
            (1.0, 1.0),
            (4.0, 1.0),
            (20.0, 1.0),
            (39.0, 1.0),
            (300.0, 100.0),
            (2.0, 0.001),
            (5_000.0, 10.0),
        ] {
            let wave = BlastWave::surface_burst(standoff, mass).expect("valid");
            // Composite Simpson rule over the positive phase.
            let n = 4_000usize;
            let h = wave.positive_duration / n as f64;
            let mut sum = wave.overpressure_at(wave.arrival_time);
            for k in 1..n {
                let weight = if k % 2 == 1 { 4.0 } else { 2.0 };
                sum += weight * wave.overpressure_at(wave.arrival_time + h * k as f64);
            }
            sum += wave.overpressure_at(wave.arrival_time + h * n as f64);
            let integral = sum * h / 3.0;
            assert_rel_eq(
                integral,
                wave.positive_impulse,
                1e-6,
                &format!("impulse at R = {standoff} m, W = {mass} kg"),
            );
        }
    }

    #[test]
    fn invalid_inputs_give_no_wave() {
        assert!(BlastWave::surface_burst(5.0, 0.0).is_none());
        assert!(BlastWave::surface_burst(5.0, -3.0).is_none());
        assert!(BlastWave::surface_burst(5.0, f64::NAN).is_none());
        assert!(BlastWave::surface_burst(5.0, f64::INFINITY).is_none());
        assert!(BlastWave::surface_burst(-1.0, 10.0).is_none());
        assert!(BlastWave::surface_burst(f64::NAN, 10.0).is_none());
        assert!(BlastWave::surface_burst(f64::INFINITY, 10.0).is_none());
        assert!(BlastWave::surface_burst(5.0, 10.0).is_some());
    }

    /// Kinney–Graham (1985) free-air peak overpressure ratio, an independent fit.
    fn kinney_graham_ratio(z: f64) -> f64 {
        let numerator = 808.0 * (1.0 + (z / 4.5).powi(2));
        let denominator = (1.0 + (z / 0.048).powi(2)).sqrt()
            * (1.0 + (z / 0.32).powi(2)).sqrt()
            * (1.0 + (z / 1.35).powi(2)).sqrt();
        numerator / denominator
    }

    #[test]
    fn consistent_with_the_kinney_graham_free_air_fit() {
        // A surface burst of W kg acts like a free-air burst of about 1.8 W kg
        // (ground reflection), i.e. at the scaled distance Z / 1.8^(1/3).
        for z in [0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0] {
            let surface = ScaledSurfaceBurst::at(z).peak_overpressure_kpa;
            let free_air = 101.325 * kinney_graham_ratio(z / 1.8_f64.cbrt());
            let ratio = surface / free_air;
            assert!(
                (0.8..1.2).contains(&ratio),
                "Z = {z}: surface {surface} kPa vs free-air {free_air} kPa"
            );
        }
    }
}
