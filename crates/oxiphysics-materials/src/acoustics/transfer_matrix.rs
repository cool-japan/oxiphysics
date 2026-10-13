// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Transfer-matrix method for the normal-incidence transmission of layered
//! acoustic media.

use std::f64::consts::PI;
use std::ops::{Add, Mul};

/// Minimal complex number for the transfer-matrix product (the layer matrix
/// of a lossless acoustic layer has purely imaginary off-diagonal entries).
#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };
    const ONE: Self = Self { re: 1.0, im: 0.0 };

    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn scale(self, factor: f64) -> Self {
        Self::new(self.re * factor, self.im * factor)
    }

    fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.re + other.re, self.im + other.im)
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )
    }
}

/// 2×2 complex transfer matrix `[[a, b], [c, d]]` relating (pressure, particle
/// velocity) on the two sides of a layer.
#[derive(Debug, Clone, Copy)]
struct TransferMatrix {
    a: Complex,
    b: Complex,
    c: Complex,
    d: Complex,
}

impl TransferMatrix {
    const IDENTITY: Self = Self {
        a: Complex::ONE,
        b: Complex::ZERO,
        c: Complex::ZERO,
        d: Complex::ONE,
    };

    /// Matrix of a lossless layer of characteristic impedance `z`,
    /// wavenumber `k` and thickness `d`:
    /// `[[cos kd, j·Z·sin kd], [j·sin kd / Z, cos kd]]`.
    fn layer(z: f64, k: f64, thickness: f64) -> Self {
        let (sin_kd, cos_kd) = (k * thickness).sin_cos();
        Self {
            a: Complex::new(cos_kd, 0.0),
            b: Complex::new(0.0, z * sin_kd),
            c: Complex::new(0.0, sin_kd / z),
            d: Complex::new(cos_kd, 0.0),
        }
    }
}

impl Mul for TransferMatrix {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }
}

/// Transfer-matrix method for a multilayer system at normal incidence.
///
/// The stack consists of N media: an incident half-space, N−2 finite layers
/// and a transmission half-space, listed in order from the incident side.
///
/// * `impedances` – characteristic impedances `Z_i = ρ_i·c_i` of all N media
///   (Pa·s/m).
/// * `sound_speeds` – sound speeds `c_i` of all N media (m/s); the entries of
///   the two half-spaces are not used but must be present so that the slice
///   is parallel to `impedances`.
/// * `thicknesses` – thicknesses of the N−2 interior layers (m).
/// * `freq` – frequency in Hz.
///
/// Each layer contributes `[[cos k·d, j·Z·sin k·d], [j·sin k·d / Z, cos k·d]]`
/// with its own wavenumber `k = 2π·f / c_i`; the layer matrices are multiplied
/// in propagation order.  With the total matrix `M` the pressure transmission
/// coefficient is
/// `T = 2 / (M₁₁ + M₁₂/Z_N + Z_1·M₂₁ + Z_1·M₂₂/Z_N)` and the returned value is
/// the *intensity* (power) transmission coefficient
/// `τ = (Z_1/Z_N)·|T|²`, which lies in `[0, 1]` for lossless media.
///
/// # Errors
///
/// * [`crate::Error::DimensionMismatch`] if fewer than two media are given,
///   `sound_speeds.len() != N` or `thicknesses.len() != N - 2`.
/// * [`crate::Error::ParameterOutOfRange`] for a non-positive or non-finite
///   impedance or layer sound speed, a negative or non-finite thickness, or a
///   negative or non-finite frequency.
pub fn multilayer_transmission(
    impedances: &[f64],
    sound_speeds: &[f64],
    thicknesses: &[f64],
    freq: f64,
) -> Result<f64, crate::Error> {
    let n = impedances.len();
    if n < 2 {
        return Err(crate::Error::DimensionMismatch {
            expected: "at least 2 impedances (incident medium, layers..., transmission medium)"
                .to_string(),
            actual: format!("{n} impedance(s)"),
        });
    }
    if sound_speeds.len() != n {
        return Err(crate::Error::DimensionMismatch {
            expected: format!("{n} sound speeds (one per impedance)"),
            actual: format!("{} sound speed(s)", sound_speeds.len()),
        });
    }
    let n_layers = n - 2;
    if thicknesses.len() != n_layers {
        return Err(crate::Error::DimensionMismatch {
            expected: format!("{n_layers} layer thicknesses (impedances.len() - 2)"),
            actual: format!("{} thickness(es)", thicknesses.len()),
        });
    }
    if !freq.is_finite() || freq < 0.0 {
        return Err(crate::Error::below_minimum("freq", freq, 0.0));
    }
    for &z in impedances {
        if !z.is_finite() || z <= 0.0 {
            return Err(crate::Error::below_minimum("impedance", z, 0.0));
        }
    }
    let interior_speeds = sound_speeds.iter().skip(1).take(n_layers);
    for (&c, &d) in interior_speeds.zip(thicknesses.iter()) {
        if !c.is_finite() || c <= 0.0 {
            return Err(crate::Error::below_minimum("sound_speed", c, 0.0));
        }
        if !d.is_finite() || d < 0.0 {
            return Err(crate::Error::below_minimum("thickness", d, 0.0));
        }
    }

    let omega = 2.0 * PI * freq;
    let mut total = TransferMatrix::IDENTITY;
    let layers = impedances
        .iter()
        .zip(sound_speeds.iter())
        .skip(1)
        .take(n_layers)
        .zip(thicknesses.iter());
    for ((&z, &c), &d) in layers {
        total = total * TransferMatrix::layer(z, omega / c, d);
    }

    let z_in = impedances[0];
    let z_out = impedances[n - 1];
    let denominator =
        total.a + total.b.scale(1.0 / z_out) + total.c.scale(z_in) + total.d.scale(z_in / z_out);
    // |T|² = 4 / |denominator|²; the denominator of a lossless stack never
    // vanishes (|denominator|² >= 4·Z_1/Z_N), so the division is safe.
    Ok(z_in / z_out * 4.0 / denominator.norm_sqr())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::acoustics::transmission_loss_db;

    /// Independent reference for the intensity transmission of a lossless
    /// stack: `τ = 1 - |R|²` from the recursive input impedance
    /// `Z_in = Z (Z_L + j Z tan kd) / (Z + j Z_L tan kd)`, evaluated with
    /// hand-written complex arithmetic.
    fn multilayer_reference(z: &[f64], c: &[f64], d: &[f64], freq: f64) -> f64 {
        let omega = 2.0 * PI * freq;
        let n = z.len();
        let mut load = (z[n - 1], 0.0_f64);
        for i in (1..n - 1).rev() {
            let t = (omega * d[i - 1] / c[i]).tan();
            let zi = z[i];
            let num = (zi * load.0, zi * (load.1 + zi * t));
            let den = (zi - load.1 * t, load.0 * t);
            let den2 = den.0 * den.0 + den.1 * den.1;
            load = (
                (num.0 * den.0 + num.1 * den.1) / den2,
                (num.1 * den.0 - num.0 * den.1) / den2,
            );
        }
        let r_num = (load.0 - z[0], load.1);
        let r_den = (load.0 + z[0], load.1);
        1.0 - (r_num.0 * r_num.0 + r_num.1 * r_num.1) / (r_den.0 * r_den.0 + r_den.1 * r_den.1)
    }
    #[test]
    fn test_multilayer_trivial_no_layers() {
        // Two half-spaces: τ = 4·Z1·Z2/(Z1+Z2)², identical to the interface
        // formula behind `transmission_loss_db`.
        let tau = multilayer_transmission(&[400.0, 1600.0], &[343.0, 1500.0], &[], 1000.0).unwrap();
        assert!((tau - 0.64).abs() < 1e-12, "τ = {tau}");
        let tl = transmission_loss_db(400.0, 1600.0);
        assert!((-10.0 * tau.log10() - tl).abs() < 1e-10);
        let same = multilayer_transmission(&[400.0, 400.0], &[343.0, 343.0], &[], 1000.0).unwrap();
        assert!((same - 1.0).abs() < 1e-12, "τ = {same}");
    }
    #[test]
    fn test_multilayer_matched_medium_is_transparent() {
        // Regression: a matched [400, 400, 400] stack used to give |T|² ≈ 2.43.
        for &freq in &[1.0, 100.0, 1000.0, 12_345.6] {
            let tau = multilayer_transmission(
                &[400.0, 400.0, 400.0],
                &[343.0, 343.0, 343.0],
                &[0.1],
                freq,
            )
            .unwrap();
            assert!((tau - 1.0).abs() < 1e-12, "f = {freq}: τ = {tau}");
        }
    }
    #[test]
    fn test_multilayer_quarter_wave_matching_layer() {
        // A quarter-wave layer with Z = sqrt(Z1·Z2) matches the two media
        // perfectly at its design frequency.
        let (z1, z2) = (400.0_f64, 1600.0_f64);
        let z_layer = (z1 * z2).sqrt();
        let c_layer = 1000.0;
        let f0 = 1000.0;
        let d = c_layer / (4.0 * f0);
        let z = [z1, z_layer, z2];
        let c = [343.0, c_layer, 1500.0];
        let at_design = multilayer_transmission(&z, &c, &[d], f0).unwrap();
        assert!((at_design - 1.0).abs() < 1e-12, "τ(f0) = {at_design}");
        // At 2·f0 the layer is a half wave thick and electrically invisible.
        let half_wave = multilayer_transmission(&z, &c, &[d], 2.0 * f0).unwrap();
        assert!((half_wave - 0.64).abs() < 1e-12, "τ(2·f0) = {half_wave}");
        // Away from the design frequency the match degrades.
        let off_design = multilayer_transmission(&z, &c, &[d], 0.5 * f0).unwrap();
        assert!(off_design < at_design - 1e-3, "τ(f0/2) = {off_design}");
    }
    #[test]
    fn test_multilayer_matches_input_impedance_recursion() {
        let z = [413.0, 2.0e5, 1.1e6, 3.0e5, 1.5e6];
        let c = [343.0, 800.0, 3000.0, 1200.0, 1500.0];
        let d = [0.02, 0.05, 0.01];
        for &freq in &[137.0, 777.0, 2500.0, 9100.0] {
            let got = multilayer_transmission(&z, &c, &d, freq).unwrap();
            let want = multilayer_reference(&z, &c, &d, freq);
            assert!(
                (got - want).abs() < 1e-9 * want.abs().max(1e-6),
                "f = {freq}: τ = {got}, reference {want}"
            );
        }
    }
    #[test]
    fn test_multilayer_layer_order_follows_propagation_direction() {
        // Asymmetric stack with different outer media: reversing the layer
        // order changes the result, and the reversed stack seen from the
        // other side (reciprocity) gives the original value back.
        let z = [413.0, 2.0e5, 1.1e6, 1.5e6];
        let c = [343.0, 800.0, 3000.0, 1500.0];
        let d = [0.02, 0.05];
        let freq = 1500.0;
        let forward = multilayer_transmission(&z, &c, &d, freq).unwrap();
        let reference = multilayer_reference(&z, &c, &d, freq);
        assert!((forward - reference).abs() < 1e-9 * reference);

        let z_rev = [z[3], z[2], z[1], z[0]];
        let c_rev = [c[3], c[2], c[1], c[0]];
        let d_rev = [d[1], d[0]];
        let reciprocal = multilayer_transmission(&z_rev, &c_rev, &d_rev, freq).unwrap();
        assert!(
            (forward - reciprocal).abs() < 1e-12,
            "{forward} vs {reciprocal}"
        );

        // Same outer media, layers swapped: physically a different stack.
        let z_swapped = [z[0], z[2], z[1], z[3]];
        let c_swapped = [c[0], c[2], c[1], c[3]];
        let d_swapped = [d[1], d[0]];
        let swapped = multilayer_transmission(&z_swapped, &c_swapped, &d_swapped, freq).unwrap();
        let swapped_reference = multilayer_reference(&z_swapped, &c_swapped, &d_swapped, freq);
        assert!((swapped - swapped_reference).abs() < 1e-9 * swapped_reference);
        assert!(
            (forward - swapped).abs() > 0.1 * forward,
            "{forward} vs {swapped}"
        );
    }
    #[test]
    fn test_multilayer_uses_the_supplied_sound_speeds() {
        // The wavenumber must come from the layer's own sound speed, not from
        // an assumed density.
        let z = [413.0, 1.0e6, 1.5e6];
        let slow = multilayer_transmission(&z, &[343.0, 500.0, 1500.0], &[0.05], 2000.0).unwrap();
        let fast = multilayer_transmission(&z, &[343.0, 5000.0, 1500.0], &[0.05], 2000.0).unwrap();
        assert!((slow - fast).abs() > 1e-3, "slow {slow}, fast {fast}");
    }
    #[test]
    fn test_multilayer_energy_is_conserved_over_a_sweep() {
        let z = [413.0, 3.0e5, 9.0e5, 1.5e6];
        let c = [343.0, 900.0, 2400.0, 1500.0];
        let d = [0.013, 0.007];
        let mut freq = 20.0;
        while freq < 20_000.0 {
            let tau = multilayer_transmission(&z, &c, &d, freq).unwrap();
            assert!(
                tau > 0.0 && tau <= 1.0 + 1e-12,
                "f = {freq}: τ = {tau} violates 0 < τ <= 1"
            );
            freq *= 1.37;
        }
    }
    #[test]
    fn test_multilayer_rejects_inconsistent_slice_lengths() {
        use crate::Error;
        let z = [400.0, 400.0, 400.0];
        let c = [343.0, 343.0, 343.0];
        // thicknesses must have N - 2 entries
        assert!(matches!(
            multilayer_transmission(&z, &c, &[], 1000.0),
            Err(Error::DimensionMismatch { .. })
        ));
        assert!(matches!(
            multilayer_transmission(&z, &c, &[0.1, 0.2], 1000.0),
            Err(Error::DimensionMismatch { .. })
        ));
        // sound speeds must have N entries
        assert!(matches!(
            multilayer_transmission(&z, &c[..2], &[0.1], 1000.0),
            Err(Error::DimensionMismatch { .. })
        ));
        // fewer than two media is not a stack
        assert!(matches!(
            multilayer_transmission(&[400.0], &[343.0], &[], 1000.0),
            Err(Error::DimensionMismatch { .. })
        ));
        assert!(matches!(
            multilayer_transmission(&[], &[], &[], 1000.0),
            Err(Error::DimensionMismatch { .. })
        ));
    }
    #[test]
    fn test_multilayer_rejects_non_physical_values() {
        use crate::Error;
        let z = [400.0, 400.0, 400.0];
        let c = [343.0, 343.0, 343.0];
        let out_of_range =
            |r: Result<f64, Error>| matches!(r, Err(Error::ParameterOutOfRange { .. }));
        assert!(out_of_range(multilayer_transmission(
            &[400.0, 0.0, 400.0],
            &c,
            &[0.1],
            1000.0
        )));
        assert!(out_of_range(multilayer_transmission(
            &[400.0, -5.0, 400.0],
            &c,
            &[0.1],
            1000.0
        )));
        assert!(out_of_range(multilayer_transmission(
            &z,
            &[343.0, 0.0, 343.0],
            &[0.1],
            1000.0
        )));
        assert!(out_of_range(multilayer_transmission(
            &z,
            &c,
            &[-0.1],
            1000.0
        )));
        assert!(out_of_range(multilayer_transmission(&z, &c, &[0.1], -1.0)));
        assert!(out_of_range(multilayer_transmission(
            &z,
            &c,
            &[0.1],
            f64::NAN
        )));
        // Unused half-space sound speeds are not validated.
        assert!(multilayer_transmission(&z, &[0.0, 343.0, 0.0], &[0.1], 1000.0).is_ok());
    }
}
