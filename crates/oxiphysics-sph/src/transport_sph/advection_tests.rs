//! Tests of the flux-limited [`AdvectionSchemeKind::WenoSph`] scheme (van Leer
//! limiter on the gradient ratio) and of the semi-Lagrangian back-tracing of
//! [`AdvectionSPH`].

use super::functions::{cubic_kernel, cubic_kernel_grad};
use super::types::{
    AdvectionSPH, AdvectionScheme, AdvectionSchemeKind, TransportParticle, van_leer_limiter,
};

const RHO: f64 = 1000.0;
const DX: f64 = 0.05;
const H: f64 = 1.25 * DX;

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

/// Slab of particles with `ix ∈ [lo, hi)` and `|iy|, |iz| ≤ ny`, velocity `v`, and the
/// concentration `field(x)`.  Returns the particles and the indices of the centre line
/// (`iy = iz = 0`) ordered by `ix`.
fn slab(
    lo: i32,
    hi: i32,
    ny: i32,
    velocity: [f64; 3],
    field: impl Fn([f64; 3]) -> f64,
) -> (Vec<TransportParticle>, Vec<usize>) {
    let mut particles = Vec::new();
    let mut line = Vec::new();
    for ix in lo..hi {
        for iy in -ny..=ny {
            for iz in -ny..=ny {
                let pos = [ix as f64 * DX, iy as f64 * DX, iz as f64 * DX];
                if iy == 0 && iz == 0 {
                    line.push(particles.len());
                }
                particles.push(TransportParticle::new(
                    pos,
                    velocity,
                    RHO * DX * DX * DX,
                    RHO,
                    H,
                    field(pos),
                    300.0,
                ));
            }
        }
    }
    (particles, line)
}

fn rates(kind: AdvectionSchemeKind, particles: &[TransportParticle], line: &[usize]) -> Vec<f64> {
    let scheme = AdvectionScheme::new(kind, 0.0);
    line.iter()
        .map(|&i| scheme.compute_rate_i(particles, i))
        .collect()
}

// --- the limiter --------------------------------------------------------------------

#[test]
fn van_leer_limiter_has_the_tvd_properties() {
    let psi = van_leer_limiter;
    assert_eq!(psi(1.0), 1.0);
    // Local extrema (r ≤ 0) are limited to first order.
    for r in [0.0, -0.1, -1.0, -50.0] {
        assert_eq!(psi(r), 0.0, "ψ({r})");
    }
    // 0 ≤ ψ < 2, increasing, and the second-order symmetry ψ(r)/r = ψ(1/r).
    let mut previous = 0.0;
    for k in 1..=400 {
        let r = 0.05 * k as f64;
        let value = psi(r);
        assert!((0.0..2.0).contains(&value), "ψ({r}) = {value}");
        assert!(value >= previous, "ψ must not decrease");
        previous = value;
        assert!((value / r - psi(1.0 / r)).abs() < 1e-12, "r = {r}");
    }
    // Below the diagonal ψ(r) ≤ 2r for r ≤ 1 (second-order TVD region).
    for k in 1..=20 {
        let r = 0.05 * k as f64;
        assert!(psi(r) <= 2.0 * r + 1e-15);
    }
}

// --- WenoSph scheme -----------------------------------------------------------------

#[test]
fn weno_reduces_to_the_standard_scheme_on_a_linear_profile() {
    // Smooth monotone data has gradient ratio r ≈ 1 and limiter ≈ 1 (the SPH gradient is
    // within 1 % of the exact one), so both schemes give −v ∂c/∂x = −3.
    let (particles, line) = slab(-10, 10, 2, [1.0, 0.0, 0.0], |x| 3.0 * x[0]);
    let standard = rates(AdvectionSchemeKind::Standard, &particles, &line);
    let weno = rates(AdvectionSchemeKind::WenoSph, &particles, &line);
    for k in 4..line.len() - 4 {
        assert!(
            rel_err(standard[k], -3.0) < 0.03,
            "standard {}",
            standard[k]
        );
        assert!(
            rel_err(weno[k], standard[k]) < 0.05,
            "weno {} vs standard {}",
            weno[k],
            standard[k]
        );
    }
}

#[test]
fn weno_limiter_is_not_identically_one() {
    // A step profile: the unlimited central scheme and the limited one must differ
    // clearly next to the step (the former limiter was ≈ 1 for every profile).
    let (particles, line) = slab(-10, 10, 2, [1.0, 0.0, 0.0], |x| {
        if x[0] < -1.0e-9 { 1.0 } else { 0.0 }
    });
    let standard = rates(AdvectionSchemeKind::Standard, &particles, &line);
    let weno = rates(AdvectionSchemeKind::WenoSph, &particles, &line);
    // line[k] has ix = k − 10; the step sits between ix = −1 and ix = 0.
    for k in [9, 10] {
        assert!(
            standard[k].abs() > 5.0,
            "central scheme reacts strongly: {}",
            standard[k]
        );
        assert!(
            weno[k].abs() < 0.2 * standard[k].abs(),
            "limited rate {} vs central {}",
            weno[k],
            standard[k]
        );
    }
    // Away from the step both are zero.
    for k in [3, 4, 5, 15, 16] {
        assert!(standard[k].abs() < 1e-9 && weno[k].abs() < 1e-9);
    }
}

#[test]
fn weno_suppresses_the_overshoot_of_the_central_scheme_at_a_step() {
    // One explicit step of a unit step profile advected with v = 1 at CFL 0.25.  The
    // central scheme pushes the high side of the step above 1; the limited scheme
    // keeps the new values inside [0, 1] to within 1 %.
    let (particles, line) = slab(-10, 10, 2, [1.0, 0.0, 0.0], |x| {
        if x[0] < -1.0e-9 { 1.0 } else { 0.0 }
    });
    let dt = 0.25 * DX;
    let advance = |kind: AdvectionSchemeKind| -> (f64, f64) {
        let r = rates(kind, &particles, &line);
        let new: Vec<f64> = line
            .iter()
            .zip(r.iter())
            .map(|(&i, rate)| particles[i].conc() + rate * dt)
            .collect();
        let max = new.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min = new.iter().cloned().fold(f64::INFINITY, f64::min);
        (max, min)
    };
    let (standard_max, standard_min) = advance(AdvectionSchemeKind::Standard);
    let (weno_max, weno_min) = advance(AdvectionSchemeKind::WenoSph);
    assert!(
        standard_max > 1.05,
        "central overshoot expected: {standard_max}"
    );
    assert!(weno_max < 1.01, "limited overshoot {weno_max}");
    assert!(weno_max - 1.0 < 0.2 * (standard_max - 1.0));
    assert!(standard_min > -1e-9 && weno_min > -1e-9);
}

#[test]
fn weno_gives_zero_rate_for_a_uniform_field() {
    let (particles, line) = slab(-6, 6, 2, [1.0, 0.5, 0.0], |_| 0.7);
    for rate in rates(AdvectionSchemeKind::WenoSph, &particles, &line) {
        assert_eq!(rate, 0.0);
    }
}

// --- AdvectionSPH back-tracing ------------------------------------------------------

/// Tent kernel `W(r, h) = max(0, 1 − r/(2h))` (a value-only kernel for hand calculations).
fn tent(r: f64, h: f64) -> f64 {
    (1.0 - r / (2.0 * h)).max(0.0)
}

fn line_particle(x: f64, c: f64, v: f64) -> TransportParticle {
    TransportParticle::new([x, 0.0, 0.0], [v, 0.0, 0.0], 1.0, RHO, 0.1, c, 300.0)
}

#[test]
fn back_tracing_interpolates_with_the_supplied_kernel() {
    // Particles at x = 0, 0.05, 0.10 with c = (0, 1, 0); particle 1 moves with v = 1 and
    // Δt = 0.02, so x_dep = 0.03.  The tent weights are 1 − r/0.2 for r = 0.03, 0.02, 0.07
    // (equal volumes): c_dep = 0.9/(0.85 + 0.9 + 0.65) = 0.375.  At the particle itself the
    // weights are 0.75, 1.0, 0.75, so the interpolant there is ĉ = 1.0/2.5 = 0.4 and the
    // rate is (0.375 − 0.4)/0.02.
    let particles = vec![
        line_particle(0.0, 0.0, 1.0),
        line_particle(0.05, 1.0, 1.0),
        line_particle(0.10, 0.0, 1.0),
    ];
    let adv = AdvectionSPH::new(0.02, 0.0, true);
    let rate = adv.advect_scalar(&particles, 1, tent, cubic_kernel_grad);
    let expected = (0.9 / 2.4 - 1.0 / 2.5) / 0.02;
    assert!(
        rel_err(rate, expected) < 1e-12,
        "rate {rate} vs hand-computed {expected}"
    );
    // A different kernel gives a different answer: the argument is really used.  With
    // support h the weights at x_dep are 0.7, 0.8 and 0.3 (c_dep = 0.8/1.8) and at the
    // particle 0.5, 1.0, 0.5 (ĉ = 1.0/2.0).
    let narrow = adv.advect_scalar(
        &particles,
        1,
        |r, h| (1.0 - r / h).max(0.0),
        cubic_kernel_grad,
    );
    let expected_narrow = (0.8 / 1.8 - 1.0 / 2.0) / 0.02;
    assert!(
        rel_err(narrow, expected_narrow) < 1e-12,
        "{narrow} vs {expected_narrow}"
    );
    assert!((narrow - rate).abs() > 1.0);
}

#[test]
fn back_tracing_does_not_smooth_a_profile_at_rest() {
    // The Shepard interpolant is not nodally exact (it is 0.4 at the peak particle
    // below, whose value is 1), so a rate built from c_i would smooth the peak even
    // without any flow.  Comparing the departure value with the interpolant at the
    // particle itself gives exactly zero for v = 0 on any profile.
    let particles = vec![
        line_particle(0.0, 0.0, 0.0),
        line_particle(0.05, 1.0, 0.0),
        line_particle(0.10, 0.0, 0.0),
    ];
    let adv = AdvectionSPH::new(0.02, 0.0, true);
    for i in 0..particles.len() {
        let rate = adv.advect_scalar(&particles, i, tent, cubic_kernel_grad);
        assert_eq!(rate, 0.0, "particle {i} at rest must not change");
    }
}

#[test]
fn back_tracing_of_a_uniform_field_is_exactly_zero() {
    let (particles, line) = slab(-5, 5, 3, [1.0, 0.0, 0.0], |_| 0.42);
    let adv = AdvectionSPH::new(0.01, 0.0, true);
    for &i in &line {
        let rate = adv.advect_scalar(&particles, i, cubic_kernel, cubic_kernel_grad);
        assert!(rate.abs() < 1e-9, "uniform field: {rate}");
    }
}

#[test]
fn back_tracing_without_velocity_or_time_interval_has_no_convective_term() {
    let (mut particles, line) = slab(-5, 5, 3, [0.0; 3], |x| 2.0 * x[0]);
    let centre = line[5];
    let at_rest = AdvectionSPH::new(0.01, 0.0, true);
    assert!(
        at_rest
            .advect_scalar(&particles, centre, cubic_kernel, cubic_kernel_grad)
            .abs()
            < 1e-9
    );
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    for dt in [0.0, -0.01, f64::NAN] {
        let adv = AdvectionSPH::new(dt, 0.0, true);
        let rate = adv.advect_scalar(&particles, centre, cubic_kernel, cubic_kernel_grad);
        assert_eq!(rate, 0.0, "dt = {dt}");
    }
}

#[test]
fn back_tracing_is_zero_when_the_departure_point_has_no_neighbours() {
    // Fast particle, tiny kernel support: x_dep is outside every particle's support.
    let particles = vec![
        line_particle(0.0, 1.0, 100.0),
        line_particle(0.05, 0.0, 100.0),
    ];
    let adv = AdvectionSPH::new(0.1, 0.0, true);
    assert_eq!(
        adv.advect_scalar(&particles, 0, cubic_kernel, cubic_kernel_grad),
        0.0
    );
}

#[test]
fn artificial_diffusion_stabilises_the_back_traced_scheme() {
    // Two-particle step: the stabiliser pulls |c_i − c_j| together.
    let particles = vec![line_particle(0.0, 1.0, 1.0), line_particle(0.05, 0.0, 1.0)];
    let plain = AdvectionSPH::new(0.001, 0.0, false);
    let damped = AdvectionSPH::new(0.001, 0.5, false);
    let r0_plain = plain.advect_scalar(&particles, 0, cubic_kernel, cubic_kernel_grad);
    let r1_plain = plain.advect_scalar(&particles, 1, cubic_kernel, cubic_kernel_grad);
    assert_eq!((r0_plain, r1_plain), (0.0, 0.0));
    let r0 = damped.advect_scalar(&particles, 0, cubic_kernel, cubic_kernel_grad);
    let r1 = damped.advect_scalar(&particles, 1, cubic_kernel, cubic_kernel_grad);
    assert!(r0 < 0.0 && r1 > 0.0, "{r0}, {r1}");
    let gap_before = particles[0].conc() - particles[1].conc();
    let gap_after = (particles[0].conc() + r0 * 0.001) - (particles[1].conc() + r1 * 0.001);
    assert!(gap_after < gap_before && gap_after > 0.0);
}
