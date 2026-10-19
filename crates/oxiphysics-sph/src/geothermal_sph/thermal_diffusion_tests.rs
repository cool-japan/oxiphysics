// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Regression tests for the SPH thermal diffusion of the geothermal solver.
//
// The pair term must be the Cleary & Monaghan (1999) conduction operator: it
// carries `1/(ρ_i cₚ)` (so the result is a temperature rate in K/s), the
// factor 4 κ_i κ_j/(κ_i + κ_j) and the Brookshaw regularisation, and it must
// reproduce the diffusion equation `∂T/∂t = α ∇²T`, `α = κ/(ρ cₚ)`.

use std::f64::consts::PI;

use super::*;

/// Thermal conductivity of the synthetic fluid (W/(m K)).
const KAPPA: f64 = 0.6;
/// Density of the synthetic fluid (kg/m³).
const RHO: f64 = 1000.0;
/// Specific heat of the synthetic fluid (J/(kg K)).
const CP: f64 = 4184.0;
/// Thermal diffusivity `κ/(ρ cₚ)` of the synthetic fluid (m²/s).
const ALPHA: f64 = KAPPA / (RHO * CP);

/// Cubic-spline `dW/dr` (support `2h`) written independently of the crate's
/// kernel routines.
fn cubic_spline_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        -alpha / h * 0.75 * (2.0 - q) * (2.0 - q)
    } else {
        0.0
    }
}

/// Pair parameters of the synthetic fluid for two temperatures and positions.
fn synthetic_pair(pos_i: [f64; 3], pos_j: [f64; 3], t_i: f64, t_j: f64) -> ThermalPairParams {
    ThermalPairParams {
        pos_i,
        pos_j,
        rho_i: RHO,
        rho_j: RHO,
        t_i,
        t_j,
        kappa_i: KAPPA,
        kappa_j: KAPPA,
        cp_i: CP,
    }
}

/// Complementary error function (Abramowitz & Stegun 7.1.26, |error| < 1.5e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * z);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let tail = poly * (-z * z).exp();
    if x >= 0.0 { tail } else { 2.0 - tail }
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Positions of an `n³` cubic lattice with spacing `dx`.
fn lattice_positions(n: usize, dx: f64) -> Vec<[f64; 3]> {
    let mut positions = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                positions.push([ix as f64 * dx, iy as f64 * dx, iz as f64 * dx]);
            }
        }
    }
    positions
}

/// Closed-form pair value `(m/ρ²) · 2κ · ΔT · r dW/dr / (c_p (r² + η²))` for the
/// synthetic fluid.
#[test]
fn test_thermal_diffusion_pair_matches_the_cleary_monaghan_formula() {
    let (h, r, dx) = (0.05, 0.04, 0.02);
    let mass = RHO * dx * dx * dx;
    let pair = synthetic_pair([0.0; 3], [r, 0.0, 0.0], 310.0, 300.0);
    let got = sph_thermal_diffusion(pair, mass, h);
    let expected = mass / (RHO * RHO) * (2.0 * KAPPA) * 10.0 * (r * cubic_spline_dw_dr(r, h))
        / (CP * (r * r + 0.01 * h * h));
    assert!(expected < 0.0, "the hotter particle must cool");
    assert!(
        (got - expected).abs() < 1e-10 * expected.abs(),
        "pair rate {got} K/s, expected {expected} K/s"
    );
    // The colder neighbour is heated by exactly the same energy flux.
    let reverse = synthetic_pair([r, 0.0, 0.0], [0.0; 3], 300.0, 310.0);
    let heating = sph_thermal_diffusion(reverse, mass, h);
    assert!((heating + got).abs() < 1e-12 * got.abs());
    assert!(heating > 0.0);
}

/// Without a temperature difference, a coincident neighbour or a degenerate
/// material there is no heat flux (and no NaN).
#[test]
fn test_thermal_diffusion_vanishes_without_a_driving_force() {
    let h = 0.05;
    let isothermal = synthetic_pair([0.0; 3], [0.03, 0.0, 0.0], 300.0, 300.0);
    assert_eq!(sph_thermal_diffusion(isothermal, 1e-3, h), 0.0);
    let coincident = synthetic_pair([0.0; 3], [0.0; 3], 310.0, 300.0);
    assert_eq!(sph_thermal_diffusion(coincident, 1e-3, h), 0.0);
    let outside_support = synthetic_pair([0.0; 3], [0.2, 0.0, 0.0], 310.0, 300.0);
    assert_eq!(sph_thermal_diffusion(outside_support, 1e-3, h), 0.0);
    let mut no_conductor = synthetic_pair([0.0; 3], [0.03, 0.0, 0.0], 310.0, 300.0);
    no_conductor.kappa_i = 0.0;
    no_conductor.kappa_j = 0.0;
    assert_eq!(sph_thermal_diffusion(no_conductor, 1e-3, h), 0.0);
}

/// For `T = T₀ + c x²` the SPH sum is `α ∇²T = 2 α c` at an interior lattice
/// particle. The old expression lacked `1/ρ_i` (a factor ~1000 here).
#[test]
fn test_thermal_diffusion_of_a_quadratic_profile_is_alpha_times_the_laplacian() {
    let dx = 0.01;
    let h = 1.3 * dx;
    let n = 7;
    let curvature = 1.0e3; // K/m²
    let mass = RHO * dx * dx * dx;
    let positions = lattice_positions(n, dx);
    let temperature = |p: &[f64; 3]| 300.0 + curvature * p[0] * p[0];
    let centre = 3 * n * n + 3 * n + 3;
    let mut rate = 0.0;
    for (j, pos_j) in positions.iter().enumerate() {
        if j == centre {
            continue;
        }
        let pair = synthetic_pair(
            positions[centre],
            *pos_j,
            temperature(&positions[centre]),
            temperature(pos_j),
        );
        rate += sph_thermal_diffusion(pair, mass, h);
    }
    let expected = ALPHA * 2.0 * curvature;
    assert!(
        (rate - expected).abs() < 0.03 * expected,
        "dT/dt = {rate} K/s, expected {expected} K/s"
    );
    // ~ 3e-4 K/s: nothing like the ~1e-1 K/s of the dimensionally wrong sum.
    assert!(rate > 2.0e-4 && rate < 4.0e-4);
}

/// Total thermal energy is conserved (pair terms cancel in `m c_p dT/dt`) for
/// unequal masses, densities, conductivities and heat capacities, and every
/// pair moves heat down the temperature gradient.
#[test]
fn test_thermal_diffusion_conserves_energy_for_inhomogeneous_particles() {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let h = 0.05;
    let count = 30;
    let mut particles = Vec::new();
    for _ in 0..count {
        let pos = [
            0.12 * lcg_unit(&mut state),
            0.12 * lcg_unit(&mut state),
            0.12 * lcg_unit(&mut state),
        ];
        let mass = 5.0e-3 * (0.5 + lcg_unit(&mut state));
        let rho = 800.0 + 400.0 * lcg_unit(&mut state);
        let temperature = 300.0 + 80.0 * lcg_unit(&mut state);
        let kappa = 0.3 + 0.6 * lcg_unit(&mut state);
        let cp = 3000.0 + 2000.0 * lcg_unit(&mut state);
        particles.push((pos, mass, rho, temperature, kappa, cp));
    }
    let mut energy_rate = 0.0f64;
    let mut scale = 0.0f64;
    for (i, pi) in particles.iter().enumerate() {
        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let pair = ThermalPairParams {
                pos_i: pi.0,
                pos_j: pj.0,
                rho_i: pi.2,
                rho_j: pj.2,
                t_i: pi.3,
                t_j: pj.3,
                kappa_i: pi.4,
                kappa_j: pj.4,
                cp_i: pi.5,
            };
            let rate = sph_thermal_diffusion(pair, pj.1, h);
            if pi.3 > pj.3 {
                assert!(rate <= 0.0, "hot particle {i} heated by cold {j}: {rate}");
            } else if pi.3 < pj.3 {
                assert!(rate >= 0.0, "cold particle {i} cooled by hot {j}: {rate}");
            }
            let power = pi.1 * pi.5 * rate;
            energy_rate += power;
            scale += power.abs();
        }
    }
    assert!(
        scale > 1.0e-2,
        "the particles must exchange heat: {scale} W"
    );
    assert!(
        energy_rate.abs() < 1e-12 * scale,
        "net energy rate {energy_rate} (scale {scale})"
    );
}

/// A hot half-space in contact with a cold one: the SPH diffusion of a step
/// profile on a lattice (periodic in the two lateral directions) follows the
/// erf solution `T = T_c + ΔT/2 · erfc((x − x₀)/(2 √(α t)))` of the 1-D heat
/// equation with `α = κ/(ρ c_p)`.
#[test]
fn test_thermal_diffusion_of_a_step_profile_follows_the_erf_solution() {
    let dx = 0.01;
    let h = 1.3 * dx;
    let (nx, ny, nz) = (48usize, 6usize, 6usize);
    let (period_y, period_z) = (ny as f64 * dx, nz as f64 * dx);
    let mass = RHO * dx * dx * dx;
    let (t_hot, t_cold) = (400.0, 300.0);
    let x0 = 23.5 * dx;

    let mut positions = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                positions.push([ix as f64 * dx, iy as f64 * dx, iz as f64 * dx]);
            }
        }
    }
    // Static neighbour lists: (index, nearest periodic image of the neighbour).
    let wrap = |d: f64, period: f64| d - period * (d / period).round();
    let neighbours: Vec<Vec<(usize, [f64; 3])>> = positions
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            positions
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .filter_map(|(j, pj)| {
                    let image = [
                        pj[0],
                        pi[1] + wrap(pj[1] - pi[1], period_y),
                        pi[2] + wrap(pj[2] - pi[2], period_z),
                    ];
                    let d = sub3(image, *pi);
                    (dot3(d, d) < 4.0 * h * h).then_some((j, image))
                })
                .collect()
        })
        .collect();

    // Heating rates of all particles for the temperatures `temps`.
    let rates = |temps: &[f64]| -> Vec<f64> {
        neighbours
            .iter()
            .enumerate()
            .map(|(i, list)| {
                list.iter()
                    .map(|&(j, image)| {
                        let pair = synthetic_pair(positions[i], image, temps[i], temps[j]);
                        sph_thermal_diffusion(pair, mass, h)
                    })
                    .sum()
            })
            .collect()
    };

    // Stable explicit step: 0.4 / (largest sum of the pair weights).
    let max_weight = neighbours
        .iter()
        .enumerate()
        .map(|(i, list)| {
            list.iter()
                .map(|&(_, image)| {
                    let unit_step = synthetic_pair(positions[i], image, 1.0, 0.0);
                    -sph_thermal_diffusion(unit_step, mass, h)
                })
                .sum::<f64>()
        })
        .fold(0.0f64, f64::max);
    assert!(max_weight > 0.0);

    let diffusion_length = 8.0 * dx; // 2 √(α t)
    let t_end = (diffusion_length / 2.0).powi(2) / ALPHA;
    let steps = (t_end * max_weight / 0.4).ceil() as usize;
    let dt = t_end / steps as f64;
    let mut temps: Vec<f64> = positions
        .iter()
        .map(|p| if p[0] < x0 { t_hot } else { t_cold })
        .collect();
    for _ in 0..steps {
        let rate = rates(&temps);
        for (t, r) in temps.iter_mut().zip(&rate) {
            *t += r * dt;
        }
    }

    let spread = t_hot - t_cold;
    let mut worst = 0.0f64;
    for (p, t) in positions.iter().zip(&temps) {
        let exact = t_cold + 0.5 * spread * erfc((p[0] - x0) / (2.0 * (ALPHA * t_end).sqrt()));
        worst = worst.max((t - exact).abs() / spread);
    }
    assert!(
        worst < 0.03,
        "largest deviation from the erf solution: {:.3} % of the step ({steps} steps)",
        100.0 * worst
    );
    // The profile really diffused: the first cold-side particle (x = x₀ + dx/2)
    // has been heated well above the cold temperature.
    let first_cold_side = temps[24 * ny * nz];
    assert!(
        first_cold_side > t_cold + 0.3 * spread && first_cold_side < t_hot - 0.2 * spread,
        "temperature next to the interface: {first_cold_side}"
    );
}

/// One `step_forward` on a lattice with `T = T₀ + c x²`: the centre particle
/// warms at `α ∇²T dt` with the actual `κ(T)`, `c_p(T)` and summation density
/// of the particle (all dependences of the real solver), not 500-1000 times
/// that.
#[test]
fn test_step_forward_diffuses_heat_at_the_physical_rate() {
    let dx = 0.01;
    let n = 7;
    let h = 1.3 * dx;
    let dt = 0.1;
    let curvature = 1.0e3; // K/m²
    let config = GeothermalConfig {
        smoothing_length: h,
        dt,
        enable_chemistry: false,
        ..GeothermalConfig::default()
    };
    let mut sim = GeothermalSimulation::new(config);
    let mass = 1000.0 * dx * dx * dx;
    let positions = lattice_positions(n, dx);
    for pos in &positions {
        let temperature = 300.0 + curvature * pos[0] * pos[0];
        sim.add_particle(GeoParticle::new(*pos, mass, temperature, P_ATM));
    }
    let centre = 3 * n * n + 3 * n + 3;
    let t_before = sim.particles[centre].temperature;

    // Summation density of the centre particle exactly as the solver computes it.
    let neighbours: Vec<(f64, [f64; 3])> = positions
        .iter()
        .filter(|p| norm3(sub3(positions[centre], **p)) < 2.0 * h)
        .map(|p| (mass, *p))
        .collect();
    let rho_centre = sph_density_sum(positions[centre], &neighbours, h);
    assert!(
        (rho_centre - 1000.0).abs() < 50.0,
        "summation density {rho_centre}"
    );

    sim.step_forward();

    let alpha = water_thermal_conductivity(t_before) / (rho_centre * water_cp(t_before));
    let expected = alpha * 2.0 * curvature * dt;
    let got = sim.particles[centre].temperature - t_before;
    assert!(
        got > 0.0,
        "the centre particle sits at the minimum of the profile and must warm: {got}"
    );
    assert!(
        (got - expected).abs() < 0.05 * expected,
        "ΔT = {got:e} K, expected {expected:e} K"
    );
}
