//! Regression tests for the settling-velocity and SPH advection-diffusion
//! fixes in [`super::types`].

use super::types::{SettlingVelocity, SphSedimentAdvection};
use crate::kernel::{CubicSplineKernel, SphKernel};

const RHO_QUARTZ: f64 = 2650.0;
const RHO_WATER: f64 = 1000.0;
const NU_WATER: f64 = 1.0e-6;
const GRAVITY: f64 = 9.81;

fn quartz(diameter: f64) -> SettlingVelocity {
    SettlingVelocity::new(diameter, RHO_QUARTZ, RHO_WATER, NU_WATER)
}

#[test]
fn rubey_reduces_to_stokes_for_fine_silt() {
    let sv = quartz(10.0e-6);
    let ratio = sv.rubey() / sv.stokes();
    assert!(
        (ratio - 1.0).abs() < 0.02,
        "rubey/stokes = {ratio} in the Stokes limit"
    );
}

#[test]
fn rubey_for_one_tenth_millimetre_quartz() {
    // xi = 36 nu^2 / ((s - 1) g d^3) = 2.224, F = sqrt(2/3 + xi) - sqrt(xi) = 0.2089,
    // w = F sqrt((s - 1) g d) = 8.40e-3 m/s.
    let w = quartz(1.0e-4).rubey();
    assert!(
        (w - 8.40e-3).abs() < 0.03 * 8.40e-3,
        "rubey(0.1 mm) = {w} m/s"
    );
}

#[test]
fn rubey_matches_the_closed_form_for_a_range_of_grain_sizes() {
    for &d in &[2.0e-5_f64, 5.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 5.0e-3] {
        let s = RHO_QUARTZ / RHO_WATER;
        let xi = 36.0 * NU_WATER * NU_WATER / ((s - 1.0) * GRAVITY * d * d * d);
        let expected = ((2.0 / 3.0 + xi).sqrt() - xi.sqrt()) * ((s - 1.0) * GRAVITY * d).sqrt();
        let w = quartz(d).rubey();
        assert!(
            (w - expected).abs() <= 1e-9 * expected,
            "d = {d}: rubey = {w}, closed form = {expected}"
        );
    }
}

#[test]
fn rubey_tends_to_the_inviscid_limit_for_gravel() {
    let d = 0.05;
    let s = RHO_QUARTZ / RHO_WATER;
    let inviscid = (2.0_f64 / 3.0).sqrt() * ((s - 1.0) * GRAVITY * d).sqrt();
    let w = quartz(d).rubey();
    assert!(
        (w - inviscid).abs() < 1e-3 * inviscid,
        "rubey(50 mm) = {w}, inviscid limit = {inviscid}"
    );
}

#[test]
fn rubey_never_exceeds_stokes_and_grows_with_diameter() {
    let mut previous = 0.0;
    for k in 0..40 {
        let d = 1.0e-6 * 10.0_f64.powf(k as f64 * 0.1);
        let sv = quartz(d);
        let w = sv.rubey();
        assert!(w > previous, "d = {d}: rubey {w} must exceed {previous}");
        assert!(
            w <= sv.stokes() * (1.0 + 1e-12),
            "d = {d}: rubey {w} exceeds Stokes {}",
            sv.stokes()
        );
        previous = w;
    }
}

#[test]
fn rubey_is_zero_for_grains_that_do_not_sink() {
    assert_eq!(quartz(0.0).rubey(), 0.0);
    assert_eq!(quartz(-1.0e-4).rubey(), 0.0);
    assert_eq!(
        SettlingVelocity::new(1.0e-4, RHO_WATER, RHO_WATER, NU_WATER).rubey(),
        0.0
    );
    assert_eq!(
        SettlingVelocity::new(1.0e-4, 900.0, RHO_WATER, NU_WATER).rubey(),
        0.0
    );
}

fn advection() -> SphSedimentAdvection {
    SphSedimentAdvection::new(0.1, 1.0e-3, 1.0e-5)
}

/// Two neighbours at `(±0.05, 0, 0)` with `m = 1`, `rho = 1000`.
fn symmetric_neighbours(c_j: f64) -> [([f64; 3], f64, f64, f64); 2] {
    [
        ([0.05, 0.0, 0.0], c_j, 1.0, 1000.0),
        ([-0.05, 0.0, 0.0], c_j, 1.0, 1000.0),
    ]
}

#[test]
fn diffusion_drains_a_concentration_peak() {
    let rate = advection().diffusion_rate(1.0, &symmetric_neighbours(0.0));
    assert!(rate < 0.0, "peak must lose sediment, got {rate}");
}

#[test]
fn diffusion_fills_a_concentration_trough() {
    let rate = advection().diffusion_rate(0.0, &symmetric_neighbours(1.0));
    assert!(rate > 0.0, "trough must gain sediment, got {rate}");
}

#[test]
fn diffusion_of_a_uniform_field_vanishes() {
    let rate = advection().diffusion_rate(0.3, &symmetric_neighbours(0.3));
    assert!(rate.abs() < 1e-18, "uniform field rate = {rate}");
}

#[test]
fn diffusion_rate_matches_the_brookshaw_laplacian() {
    let adv = advection();
    let (r, m, rho) = (0.05, 1.0, 1000.0);
    let dw_dr = CubicSplineKernel.grad_w(r, adv.h);
    // Two neighbours, each contributing 2 (m/rho) (c_i - c_j) r W'(r) / r^2.
    let expected = adv.diffusivity * 2.0 * (2.0 * (m / rho) * (1.0 - 0.0) * r * dw_dr / (r * r));
    let rate = adv.diffusion_rate(1.0, &symmetric_neighbours(0.0));
    assert!(
        (rate - expected).abs() <= 1e-12 * expected.abs(),
        "rate = {rate}, expected = {expected}"
    );
}

#[test]
fn diffusion_on_a_particle_row_flattens_the_peak_and_conserves_mass() {
    let adv = advection();
    let spacing = 0.05;
    let n = 9;
    let xs: Vec<f64> = (0..n).map(|k| k as f64 * spacing).collect();
    let conc: Vec<f64> = xs
        .iter()
        .map(|&x| (-((x - 0.2) / 0.08).powi(2)).exp())
        .collect();
    let (m, rho) = (1.0, 1000.0);
    let rates: Vec<f64> = (0..n)
        .map(|i| {
            let neighbours: Vec<([f64; 3], f64, f64, f64)> = (0..n)
                .filter(|&j| j != i && (xs[i] - xs[j]).abs() < 2.0 * adv.h)
                .map(|j| ([xs[i] - xs[j], 0.0, 0.0], conc[j], m, rho))
                .collect();
            adv.diffusion_rate(conc[i], &neighbours)
        })
        .collect();
    let peak = n / 2;
    assert!(rates[peak] < 0.0, "peak rate {}", rates[peak]);
    assert!(rates[0] > 0.0 && rates[n - 1] > 0.0, "tail rates {rates:?}");
    let total: f64 = rates.iter().map(|r| r * (m / rho)).sum();
    let scale: f64 = rates.iter().map(|r| (r * (m / rho)).abs()).sum();
    assert!(
        total.abs() <= 1e-12 * scale,
        "diffusion must conserve sediment: sum = {total}, scale = {scale}"
    );
}
