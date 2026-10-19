// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Powder X-ray diffraction pattern simulation.
//!
//! Every reflection `(h, k, l)` whose Bragg angle falls inside the requested
//! 2-theta window is weighted by its squared structure factor and the powder
//! Lorentz-polarization factor and drawn as a Gaussian whose width is given in
//! degrees of 2-theta.

use super::lattice_math::MIN_LATTICE_DET;
use super::{CrystalStructure, XrdPattern};
use crate::{Error, Result};
use std::f64::consts::PI;

/// Default full width at half maximum of a simulated Bragg peak, in degrees of
/// 2-theta.
pub const DEFAULT_XRD_FWHM_DEG: f64 = 0.1;

/// Ratio of a Gaussian's FWHM to its standard deviation, `2 sqrt(2 ln 2)`.
const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;

/// Peaks are evaluated out to this many standard deviations from their centre
/// (the Gaussian is below `4e-6` of its height there).
const PEAK_CUTOFF_SIGMAS: f64 = 5.0;

/// Reflections whose intensity is below this fraction of the strongest one are
/// dropped.  Systematically absent reflections come out of the structure
/// factor sum as rounding noise around `1e-28` of a normal reflection, which
/// would otherwise show up as spurious local maxima.
const NOISE_FLOOR: f64 = 1e-10;

/// Largest number of `(h, k, l)` triples the reflection search may visit.
const MAX_SEARCH_TRIPLES: f64 = 4.0e8;

/// Bragg angles whose cosine is below this value are skipped: the powder
/// Lorentz factor diverges at exact backscattering (2-theta = 180 degrees).
const MIN_COS_THETA: f64 = 1e-9;

/// Lorentz-polarization factor of an unpolarised powder diffractometer,
/// `(1 + cos^2 2θ) / (sin^2 θ cos θ)`, for the Bragg angle `theta` in radians.
///
/// The factor is relative; it equals `2 sqrt(2)` at `θ = 45°` and diverges
/// only towards forward (`θ → 0`) and backward (`θ → 90°`) scattering.
pub fn lorentz_polarization_factor(theta: f64) -> f64 {
    let two_theta = 2.0 * theta;
    (1.0 + two_theta.cos().powi(2)) / (theta.sin().powi(2) * theta.cos())
}

/// One Bragg reflection inside the requested angular window.
#[derive(Debug, Clone, Copy)]
struct Reflection {
    /// Peak centre in degrees of 2-theta.
    two_theta_deg: f64,
    /// `|F|^2 * LP`, the peak height.
    strength: f64,
}

/// Largest `|h|`, `|k|`, `|l|` that can reach spacings down to `d_min`.
///
/// A reflection with spacing `d` has `|h| = |G . a1| <= |a1| / d` for
/// `|G| = 1 / d`, and likewise for `k` and `l`, so each axis needs only
/// `ceil(|a_i| / d_min)` indices, however skewed the cell is.
fn hkl_limits(lattice: &[[f64; 3]; 3], d_min: f64) -> [f64; 3] {
    lattice.map(|a| {
        let length = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        (length / d_min).ceil()
    })
}

/// All reflections with a 2-theta in `[two_theta_lo, two_theta_hi]` degrees.
fn collect_reflections(
    crystal: &CrystalStructure,
    wavelength: f64,
    two_theta_lo: f64,
    two_theta_hi: f64,
) -> Result<Vec<Reflection>> {
    let theta_hi = (0.5 * two_theta_hi).to_radians();
    let sin_hi = theta_hi.sin();
    let sin_lo = (0.5 * two_theta_lo).to_radians().sin();
    if sin_hi.is_nan() || sin_hi <= 0.0 {
        return Ok(Vec::new());
    }
    let d_min = wavelength / (2.0 * sin_hi);
    let limits = hkl_limits(&crystal.lattice, d_min);
    let triples: f64 = limits.iter().map(|n| 2.0 * n + 1.0).product();
    if !triples.is_finite() || triples > MAX_SEARCH_TRIPLES {
        return Err(Error::General(format!(
            "XRD: reflection search over {triples:e} Miller index triples is too large \
             (limit {MAX_SEARCH_TRIPLES:e}); reduce the cell size or the 2-theta range"
        )));
    }
    let [h_max, k_max, l_max] = limits.map(|n| n as i32);
    let rec = crystal.reciprocal_lattice();
    let mut reflections = Vec::new();
    for h in -h_max..=h_max {
        for k in -k_max..=k_max {
            for l in -l_max..=l_max {
                let (hf, kf, lf) = (f64::from(h), f64::from(k), f64::from(l));
                let g = [
                    hf * rec[0][0] + kf * rec[1][0] + lf * rec[2][0],
                    hf * rec[0][1] + kf * rec[1][1] + lf * rec[2][1],
                    hf * rec[0][2] + kf * rec[1][2] + lf * rec[2][2],
                ];
                let g_len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if g_len < 1e-12 {
                    continue;
                }
                // Bragg's law with |G| = 2 pi / d: sin(theta) = lambda |G| / (4 pi).
                let sin_theta = wavelength * g_len / (4.0 * PI);
                if sin_theta < sin_lo || sin_theta > sin_hi {
                    continue;
                }
                let theta = sin_theta.asin();
                if theta.cos() < MIN_COS_THETA {
                    continue;
                }
                let (re, im) = crystal.structure_factor(h, k, l);
                reflections.push(Reflection {
                    two_theta_deg: 2.0 * theta.to_degrees(),
                    strength: (re * re + im * im) * lorentz_polarization_factor(theta),
                });
            }
        }
    }
    Ok(reflections)
}

/// Simulate a powder XRD pattern for the given crystal with the default peak
/// width [`DEFAULT_XRD_FWHM_DEG`].
///
/// `wavelength` in Å, `theta_range` in degrees (2θ_min, 2θ_max), `n_pts`
/// sample points.  Invalid arguments (non-positive wavelength, reversed or
/// non-finite range, a search that would be too large) give an all-zero
/// pattern; use [`simulate_xrd_with_fwhm`] to receive the error instead.
pub fn simulate_xrd(
    crystal: &CrystalStructure,
    wavelength: f64,
    theta_range: (f64, f64),
    n_pts: usize,
) -> XrdPattern {
    simulate_xrd_with_fwhm(
        crystal,
        wavelength,
        theta_range,
        n_pts,
        DEFAULT_XRD_FWHM_DEG,
    )
    .unwrap_or_else(|_| {
        let (two_theta, _) = sample_grid(theta_range, n_pts);
        XrdPattern::new(two_theta, vec![0.0; n_pts])
    })
}

/// Evenly spaced 2-theta sample points and their spacing.
fn sample_grid(theta_range: (f64, f64), n_pts: usize) -> (Vec<f64>, f64) {
    let (two_theta_min, two_theta_max) = theta_range;
    let step = if n_pts > 1 {
        (two_theta_max - two_theta_min) / (n_pts - 1) as f64
    } else {
        0.0
    };
    let grid = (0..n_pts)
        .map(|i| two_theta_min + i as f64 * step)
        .collect();
    (grid, step)
}

/// Simulate a powder XRD pattern with peaks of the given full width at half
/// maximum.
///
/// * `wavelength` in Å.
/// * `theta_range` is the 2θ window `(min, max)` in degrees; the pattern has
///   `n_pts` evenly spaced samples across it.
/// * `fwhm_deg` is the FWHM of every Gaussian peak in degrees of 2θ, so the
///   peak width does not depend on `n_pts`.
///
/// Every reflection up to the largest 2θ of the window is included (the Miller
/// index range follows from the cell geometry, not from a fixed cap), and the
/// tails of peaks centred just outside the window are drawn.  Peak heights are
/// `|F|² · LP` with [`lorentz_polarization_factor`].
///
/// # Errors
///
/// Returns an error when the wavelength or `fwhm_deg` is not positive and
/// finite, when the 2θ window is reversed or not finite, or when the
/// reflection search would visit more than `4e8` Miller index triples.
pub fn simulate_xrd_with_fwhm(
    crystal: &CrystalStructure,
    wavelength: f64,
    theta_range: (f64, f64),
    n_pts: usize,
    fwhm_deg: f64,
) -> Result<XrdPattern> {
    let (two_theta_min, two_theta_max) = theta_range;
    if !(wavelength.is_finite() && wavelength > 0.0) {
        return Err(Error::General(format!(
            "XRD: wavelength must be positive and finite, got {wavelength}"
        )));
    }
    if !(fwhm_deg.is_finite() && fwhm_deg > 0.0) {
        return Err(Error::General(format!(
            "XRD: peak FWHM must be positive and finite, got {fwhm_deg}"
        )));
    }
    if !(two_theta_min.is_finite() && two_theta_max.is_finite() && two_theta_min <= two_theta_max) {
        return Err(Error::General(format!(
            "XRD: 2-theta range must be finite with min <= max, got ({two_theta_min}, {two_theta_max})"
        )));
    }
    let (two_theta, step) = sample_grid(theta_range, n_pts);
    let mut intensity = vec![0.0_f64; n_pts];
    let volume = crystal.volume();
    if n_pts == 0 || volume.is_nan() || volume <= MIN_LATTICE_DET {
        return Ok(XrdPattern::new(two_theta, intensity));
    }

    let sigma = fwhm_deg / FWHM_PER_SIGMA;
    let half_width = PEAK_CUTOFF_SIGMAS * sigma;
    let window_lo = (two_theta_min - half_width).max(0.0);
    let window_hi = (two_theta_max + half_width).min(180.0);
    if window_lo > window_hi {
        return Ok(XrdPattern::new(two_theta, intensity));
    }
    let reflections = collect_reflections(crystal, wavelength, window_lo, window_hi)?;
    let strongest = reflections
        .iter()
        .map(|r| r.strength)
        .fold(0.0_f64, f64::max);
    let threshold = strongest * NOISE_FLOOR;

    for reflection in reflections.iter().filter(|r| r.strength > threshold) {
        let centre = reflection.two_theta_deg;
        let (first, last) = if step > 0.0 {
            let lo = ((centre - half_width - two_theta_min) / step)
                .ceil()
                .max(0.0);
            let hi = ((centre + half_width - two_theta_min) / step)
                .floor()
                .min((n_pts - 1) as f64);
            if lo > hi {
                continue;
            }
            (lo as usize, hi as usize)
        } else if (centre - two_theta_min).abs() <= half_width {
            (0, n_pts - 1)
        } else {
            continue;
        };
        let samples = two_theta[first..=last]
            .iter()
            .zip(intensity[first..=last].iter_mut());
        for (angle, value) in samples {
            let z = (angle - centre) / sigma;
            *value += reflection.strength * (-0.5 * z * z).exp();
        }
    }
    Ok(XrdPattern::new(two_theta, intensity))
}

#[cfg(test)]
mod tests {
    use super::super::CrystalAtom;
    use super::*;

    const CU_KA: f64 = 1.5406;

    fn atom(element: &str, pos: [f64; 3]) -> CrystalAtom {
        CrystalAtom {
            element: element.to_string(),
            fractional_pos: pos,
            occupancy: 1.0,
            b_factor: 0.0,
        }
    }

    /// Simple-cubic cell with one `Fe` atom (`f = 26`, so `F^2 = 676`).
    fn simple_cubic(a: f64) -> CrystalStructure {
        let mut crystal = CrystalStructure::new_cubic(a);
        crystal.atoms.push(atom("Fe", [0.0, 0.0, 0.0]));
        crystal
    }

    fn silicon() -> CrystalStructure {
        let mut crystal = CrystalStructure::new_cubic(5.431);
        let fcc = [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ];
        for base in fcc {
            crystal.atoms.push(atom("Si", base));
            crystal
                .atoms
                .push(atom("Si", [base[0] + 0.25, base[1] + 0.25, base[2] + 0.25]));
        }
        crystal
    }

    /// 2θ in degrees of the reflection with `h^2 + k^2 + l^2 = n` in a cubic cell.
    fn cubic_two_theta(a: f64, n: f64) -> f64 {
        2.0 * (CU_KA * n.sqrt() / (2.0 * a)).asin().to_degrees()
    }

    fn peak_height(pattern: &XrdPattern, centre: f64) -> f64 {
        pattern
            .two_theta
            .iter()
            .zip(&pattern.intensity)
            .filter(|(t, _)| (**t - centre).abs() < 0.01)
            .map(|(_, i)| *i)
            .fold(0.0_f64, f64::max)
    }

    /// Full width at half maximum of the highest peak, by linear interpolation.
    fn measured_fwhm(pattern: &XrdPattern) -> f64 {
        let (top, max) = pattern.intensity.iter().copied().enumerate().fold(
            (0usize, 0.0_f64),
            |best, (i, v)| if v > best.1 { (i, v) } else { best },
        );
        let half = 0.5 * max;
        let crossing = |from: usize, towards: isize| -> f64 {
            let mut i = from as isize;
            while i + towards >= 0
                && (i + towards) < pattern.intensity.len() as isize
                && pattern.intensity[(i + towards) as usize] > half
            {
                i += towards;
            }
            let inside = i as usize;
            let outside = (i + towards) as usize;
            let (t0, v0) = (pattern.two_theta[inside], pattern.intensity[inside]);
            let (t1, v1) = (pattern.two_theta[outside], pattern.intensity[outside]);
            t0 + (v0 - half) / (v0 - v1) * (t1 - t0)
        };
        crossing(top, 1) - crossing(top, -1)
    }

    #[test]
    fn lorentz_polarization_matches_the_textbook_value_at_45_degrees() {
        let lp = lorentz_polarization_factor(PI / 4.0);
        assert!(lp.is_finite());
        assert!((lp - 2.0 * 2.0_f64.sqrt()).abs() < 1e-12, "LP(45°) = {lp}");
    }

    #[test]
    fn lorentz_polarization_is_smooth_through_two_theta_90() {
        // The old denominator contained cos(2θ) and blew up at 2θ = 90°.
        let below = lorentz_polarization_factor(44.9_f64.to_radians());
        let at = lorentz_polarization_factor(45.0_f64.to_radians());
        let above = lorentz_polarization_factor(45.1_f64.to_radians());
        assert!((below - at).abs() < 0.02 && (above - at).abs() < 0.02);
        // Independent evaluation at θ = 14.22°.
        let theta = 14.22_f64.to_radians();
        let direct = (1.0 + (2.0 * theta).cos().powi(2)) / (theta.sin().powi(2) * theta.cos());
        assert_eq!(lorentz_polarization_factor(theta), direct);
    }

    #[test]
    fn fwhm_constant_is_two_root_two_ln_two() {
        let exact = 2.0 * (2.0 * 2.0_f64.ln()).sqrt();
        assert!((FWHM_PER_SIGMA - exact).abs() < 1e-14);
    }

    #[test]
    fn reflections_with_indices_above_five_are_simulated() {
        // 10 Å cubic cell, Cu Kα, 2θ ≤ 90°: h²+k²+l² = 76..84 needs an index
        // of at least 6, beyond the old fixed ±5 box.
        let crystal = simple_cubic(10.0);
        let pattern = simulate_xrd(&crystal, CU_KA, (84.0, 84.8), 801);
        let centre = cubic_two_theta(10.0, 76.0);
        // 76 = 6² + 6² + 2² only: 3 axis choices for the 2, times 8 sign choices.
        let theta = (0.5 * centre).to_radians();
        let expected = 24.0 * 676.0 * lorentz_polarization_factor(theta);
        let height = peak_height(&pattern, centre);
        assert!(height > 0.0, "no peak found at 2θ = {centre}");
        assert!(
            (height - expected).abs() / expected < 2e-4,
            "peak height {height}, expected {expected}"
        );
    }

    #[test]
    fn no_intensity_spike_at_two_theta_90() {
        // a = λ puts the {110} reflections exactly at 2θ = 90°.
        let crystal = simple_cubic(CU_KA);
        let pattern = simulate_xrd(&crystal, CU_KA, (80.0, 100.0), 2001);
        let expected = 12.0 * 676.0 * lorentz_polarization_factor(PI / 4.0);
        let height = peak_height(&pattern, 90.0);
        assert!(
            (height - expected).abs() / expected < 1e-9,
            "height {height}, expected {expected}"
        );
        let max = pattern.intensity.iter().copied().fold(0.0_f64, f64::max);
        assert!(max <= expected * (1.0 + 1e-9), "spike of {max}");
    }

    #[test]
    fn silicon_has_the_textbook_diamond_peaks() {
        let crystal = silicon();
        let pattern = simulate_xrd(&crystal, CU_KA, (20.0, 90.0), 7001);
        for (n, textbook) in [(3.0, 28.44), (8.0, 47.30), (11.0, 56.12)] {
            let centre = cubic_two_theta(5.431, n);
            assert!((centre - textbook).abs() < 0.05, "{centre} vs {textbook}");
            let (best_angle, best_value) = pattern
                .two_theta
                .iter()
                .zip(&pattern.intensity)
                .filter(|(t, _)| (**t - centre).abs() < 0.5)
                .fold(
                    (0.0, 0.0),
                    |best, (t, i)| if *i > best.1 { (*t, *i) } else { best },
                );
            assert!(best_value > 0.0, "no intensity near {textbook}");
            assert!(
                (best_angle - textbook).abs() <= 0.05,
                "strongest point near {textbook} is at {best_angle}"
            );
        }
        // (200) and (222) are forbidden in the diamond structure.
        for n in [4.0, 12.0] {
            let centre = cubic_two_theta(5.431, n);
            let height = pattern
                .two_theta
                .iter()
                .zip(&pattern.intensity)
                .filter(|(t, _)| (**t - centre).abs() < 0.3)
                .map(|(_, i)| *i)
                .fold(0.0_f64, f64::max);
            assert_eq!(height, 0.0, "forbidden reflection (n = {n}) has intensity");
        }
    }

    #[test]
    fn peak_width_is_in_degrees_whatever_the_sampling() {
        let crystal = simple_cubic(3.0);
        let window = (28.0, 32.0);
        for n_pts in [401usize, 1601, 3201] {
            let narrow = simulate_xrd_with_fwhm(&crystal, CU_KA, window, n_pts, 0.2).unwrap();
            let width = measured_fwhm(&narrow);
            assert!((width - 0.2).abs() < 0.004, "n_pts = {n_pts}: FWHM {width}");
            let default = simulate_xrd(&crystal, CU_KA, window, n_pts);
            let width = measured_fwhm(&default);
            assert!(
                (width - DEFAULT_XRD_FWHM_DEG).abs() < 0.004,
                "n_pts = {n_pts}: default FWHM {width}"
            );
        }
    }

    #[test]
    fn tails_of_peaks_just_outside_the_window_are_drawn() {
        // {100} sits at 2θ = 60° for a = λ; the window starts half a FWHM
        // above it, where a Gaussian is at exactly half its height.
        let crystal = simple_cubic(CU_KA);
        let pattern = simulate_xrd_with_fwhm(&crystal, CU_KA, (60.05, 62.0), 40, 0.1).unwrap();
        let full = 6.0 * 676.0 * lorentz_polarization_factor(30.0_f64.to_radians());
        assert!(
            (pattern.intensity[0] - 0.5 * full).abs() / full < 1e-6,
            "{} vs {}",
            pattern.intensity[0],
            0.5 * full
        );
    }

    #[test]
    fn single_sample_and_empty_grids_are_handled() {
        let crystal = simple_cubic(CU_KA);
        let single = simulate_xrd_with_fwhm(&crystal, CU_KA, (60.0, 60.0), 1, 0.1).unwrap();
        let full = 6.0 * 676.0 * lorentz_polarization_factor(30.0_f64.to_radians());
        assert_eq!(single.two_theta, vec![60.0]);
        assert!((single.intensity[0] - full).abs() / full < 1e-12);
        let far = simulate_xrd_with_fwhm(&crystal, CU_KA, (20.0, 20.0), 1, 0.1).unwrap();
        assert_eq!(far.intensity, vec![0.0]);
        let empty = simulate_xrd_with_fwhm(&crystal, CU_KA, (20.0, 80.0), 0, 0.1).unwrap();
        assert!(empty.two_theta.is_empty() && empty.intensity.is_empty());
    }

    #[test]
    fn invalid_arguments_are_errors_not_panics() {
        let crystal = simple_cubic(3.0);
        let window = (20.0, 80.0);
        for wavelength in [0.0, -1.5406, f64::NAN, f64::INFINITY] {
            assert!(simulate_xrd_with_fwhm(&crystal, wavelength, window, 100, 0.1).is_err());
        }
        for fwhm in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(simulate_xrd_with_fwhm(&crystal, CU_KA, window, 100, fwhm).is_err());
        }
        for range in [(80.0, 20.0), (f64::NAN, 80.0), (20.0, f64::INFINITY)] {
            assert!(simulate_xrd_with_fwhm(&crystal, CU_KA, range, 100, 0.1).is_err());
        }
        // The convenience wrapper degrades to an all-zero pattern.
        let zero = simulate_xrd(&crystal, f64::NAN, window, 50);
        assert_eq!(zero.intensity, vec![0.0; 50]);
        assert_eq!(zero.two_theta.len(), 50);
    }

    #[test]
    fn degenerate_lattices_give_a_zero_pattern() {
        let flat = CrystalStructure {
            lattice: [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            atoms: vec![atom("Fe", [0.0, 0.0, 0.0])],
            space_group: 1,
        };
        let pattern = simulate_xrd_with_fwhm(&flat, CU_KA, (20.0, 80.0), 30, 0.1).unwrap();
        assert_eq!(pattern.intensity, vec![0.0; 30]);
    }

    #[test]
    fn enormous_search_boxes_are_rejected() {
        let crystal = simple_cubic(1.0e5);
        let err = simulate_xrd_with_fwhm(&crystal, CU_KA, (20.0, 90.0), 100, 0.1).unwrap_err();
        assert!(err.to_string().contains("too large"), "{err}");
        // The convenience wrapper answers with zeros instead of hanging.
        let zero = simulate_xrd(&crystal, CU_KA, (20.0, 90.0), 10);
        assert_eq!(zero.intensity, vec![0.0; 10]);
    }

    #[test]
    fn hkl_limits_cover_a_skewed_cell_exactly() {
        // Triclinic cell: compare the per-axis search with a brute-force scan
        // over a much larger index box.
        let a = 4.2_f64;
        let b = 5.1_f64;
        let c = 6.3_f64;
        let (alpha, beta, gamma) = (83.0_f64, 95.0_f64, 101.0_f64);
        let (ca, cb, cg) = (
            alpha.to_radians().cos(),
            beta.to_radians().cos(),
            gamma.to_radians().cos(),
        );
        let sg = gamma.to_radians().sin();
        let cy = (ca - cb * cg) / sg;
        let cz = (1.0 - cb * cb - cy * cy).sqrt();
        let mut crystal = CrystalStructure {
            lattice: [
                [a, 0.0, 0.0],
                [b * cg, b * sg, 0.0],
                [c * cb, c * cy, c * cz],
            ],
            atoms: Vec::new(),
            space_group: 1,
        };
        crystal.atoms.push(atom("Fe", [0.1, 0.2, 0.3]));
        crystal.atoms.push(atom("O", [0.6, 0.7, 0.1]));

        let (lo, hi) = (10.0_f64, 70.0_f64);
        let mut found: Vec<f64> = collect_reflections(&crystal, CU_KA, lo, hi)
            .unwrap()
            .iter()
            .map(|r| r.two_theta_deg)
            .collect();
        found.sort_by(|x, y| x.total_cmp(y));

        let mut brute = Vec::new();
        for h in -14..=14 {
            for k in -14..=14 {
                for l in -14..=14 {
                    if h == 0 && k == 0 && l == 0 {
                        continue;
                    }
                    let d = crystal.d_spacing(h, k, l);
                    let sin_theta = CU_KA / (2.0 * d);
                    if sin_theta > 1.0 {
                        continue;
                    }
                    let two_theta = 2.0 * sin_theta.asin().to_degrees();
                    if (lo..=hi).contains(&two_theta) {
                        brute.push(two_theta);
                    }
                }
            }
        }
        brute.sort_by(|x, y| x.total_cmp(y));
        assert!(
            brute.len() > 100,
            "scan found only {} reflections",
            brute.len()
        );
        assert_eq!(found.len(), brute.len());
        for (x, y) in found.iter().zip(&brute) {
            assert!((x - y).abs() < 1e-9, "{x} vs {y}");
        }
    }
}
