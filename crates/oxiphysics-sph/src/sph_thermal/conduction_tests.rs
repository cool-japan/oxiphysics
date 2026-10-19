//! Regression tests for the conduction operator ([`compute_heat_flux`]) and for
//! the enthalpy method ([`thermal_conduction_step_with_phase`]).
//!
//! The conduction tests pin the operator against the closed-form Cleary–Monaghan
//! pair term, against `λ∇²T` on a quadratic field, and against the exact decay
//! `exp(−α k² t)` of a sinusoidal temperature mode.  The phase-change tests pin
//! the latent-heat plateau: while a particle is inside the mushy zone the heat it
//! receives melts material instead of raising its temperature.

use super::*;

/// Density of the lattice particles (kg/m³).
const RHO: f64 = 1000.0;

/// Particle of a cubic lattice with spacing `dx`: `m = ρ dx³`, `h = 1.3 dx`.
fn lattice_particle(
    position: [f64; 3],
    dx: f64,
    temperature: f64,
    conductivity: f64,
    specific_heat: f64,
) -> ThermalParticle {
    ThermalParticle::new(
        position,
        RHO * dx * dx * dx,
        RHO,
        1.3 * dx,
        temperature,
        conductivity,
        specific_heat,
    )
}

/// Relative difference `|a − b| / |b|`.
fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

/// Ice-like particle with `c_p = 2100 J/(kg K)` at temperature `t` (K).
fn ice_particle(t: f64) -> ThermalParticle {
    ThermalParticle::new([0.0; 3], 1.0e-3, 920.0, 0.01, t, 2.2, 2100.0)
}

// ---------------------------------------------------------------------------
// compute_heat_flux — magnitude and units (F4665)
// ---------------------------------------------------------------------------

#[test]
fn heat_flux_pair_matches_cleary_monaghan_formula() {
    // d T_i/dt = (1/(ρ_i c_p)) (m_j/ρ_j) [4 λ_i λ_j/(λ_i+λ_j)] (T_i − T_j) r ∂W/∂r / (r² + ε)
    let (mass, h, lam, cp) = (2.0e-3, 0.05, 0.6, 4186.0);
    let r = 0.04;
    let mut particles = vec![
        ThermalParticle::new([0.0; 3], mass, RHO, h, 400.0, lam, cp),
        ThermalParticle::new([r, 0.0, 0.0], mass, RHO, h, 300.0, lam, cp),
    ];
    compute_heat_flux(&mut particles);

    let dw_dr = cubic_kernel_deriv(r, h);
    assert!(dw_dr < 0.0);
    let eps = 1.0e-4 * h * h;
    // 4 λ²/(2 λ) = 2 λ for a uniform conductivity.
    let expected = (mass / RHO) * (2.0 * lam) * 100.0 * (r * dw_dr) / (r * r + eps) / (RHO * cp);
    assert!(expected < 0.0, "the hot particle must cool: {expected}");
    assert!(
        rel_err(particles[0].d_temperature, expected) < 1.0e-12,
        "{} vs {expected}",
        particles[0].d_temperature
    );
    assert!(
        rel_err(particles[1].d_temperature, -expected) < 1.0e-12,
        "{} vs {}",
        particles[1].d_temperature,
        -expected
    );
}

#[test]
fn heat_flux_carries_the_radial_factor_and_no_extra_two() {
    // The old operator divided by r (missing r in r·∇W) and doubled the result;
    // the pair rate must therefore scale like r ∂W/∂r / (r²+ε), not like ∂W/∂r / (r²+ε).
    let (mass, h, lam, cp) = (1.0e-3, 0.05, 1.0, 1000.0);
    let rate_at = |r: f64| {
        let mut p = vec![
            ThermalParticle::new([0.0; 3], mass, RHO, h, 310.0, lam, cp),
            ThermalParticle::new([r, 0.0, 0.0], mass, RHO, h, 300.0, lam, cp),
        ];
        compute_heat_flux(&mut p);
        p[0].d_temperature
    };
    let (r1, r2) = (0.03, 0.06);
    let eps = 1.0e-4 * h * h;
    let expected_ratio = (r2 * cubic_kernel_deriv(r2, h) / (r2 * r2 + eps))
        / (r1 * cubic_kernel_deriv(r1, h) / (r1 * r1 + eps));
    let measured_ratio = rate_at(r2) / rate_at(r1);
    assert!(
        rel_err(measured_ratio, expected_ratio) < 1.0e-12,
        "{measured_ratio} vs {expected_ratio}"
    );
}

#[test]
fn heat_flux_scales_with_conductivity_and_inverse_heat_capacity() {
    let rate = |lam: f64, cp: f64| {
        let mut p = vec![
            ThermalParticle::new([0.0; 3], 1.0e-3, RHO, 0.05, 350.0, lam, cp),
            ThermalParticle::new([0.04, 0.0, 0.0], 1.0e-3, RHO, 0.05, 300.0, lam, cp),
        ];
        compute_heat_flux(&mut p);
        p[0].d_temperature
    };
    let base = rate(1.0, 1000.0);
    assert!(rel_err(rate(3.0, 1000.0), 3.0 * base) < 1.0e-12);
    assert!(rel_err(rate(1.0, 4000.0), 0.25 * base) < 1.0e-12);
}

#[test]
fn heat_flux_uses_harmonic_mean_conductivity() {
    // 4 λ_i λ_j/(λ_i+λ_j) with (1, 3) is 3, the value 2 λ takes for the uniform λ = 1.5.
    let make = |lam_a: f64, lam_b: f64| {
        let mut p = vec![
            ThermalParticle::new([0.0; 3], 1.0e-3, RHO, 0.05, 340.0, lam_a, 1000.0),
            ThermalParticle::new([0.04, 0.0, 0.0], 1.0e-3, RHO, 0.05, 300.0, lam_b, 1000.0),
        ];
        compute_heat_flux(&mut p);
        (p[0].d_temperature, p[1].d_temperature)
    };
    let (mixed_a, mixed_b) = make(1.0, 3.0);
    let (uniform_a, uniform_b) = make(1.5, 1.5);
    assert!(rel_err(mixed_a, uniform_a) < 1.0e-12);
    assert!(rel_err(mixed_b, uniform_b) < 1.0e-12);
}

#[test]
fn heat_flux_ignores_neighbours_without_density() {
    let mut particles = vec![
        ThermalParticle::new([0.0; 3], 1.0e-3, RHO, 0.05, 400.0, 0.6, 4186.0),
        ThermalParticle::new([0.04, 0.0, 0.0], 1.0e-3, 0.0, 0.05, 300.0, 0.6, 4186.0),
    ];
    compute_heat_flux(&mut particles);
    assert_eq!(particles[0].d_temperature, 0.0);
    assert!(particles[1].d_temperature.is_finite());
}

/// Cubic lattice of `n³` particles (spacing `dx`) centred on the origin whose
/// temperature is `field(x)`; returns the particles and the central index.
fn centred_lattice(
    n: usize,
    dx: f64,
    lam: f64,
    cp: f64,
    field: impl Fn([f64; 3]) -> f64,
) -> (Vec<ThermalParticle>, usize) {
    let half = (n as f64 - 1.0) / 2.0;
    let mut particles = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [
                    (ix as f64 - half) * dx,
                    (iy as f64 - half) * dx,
                    (iz as f64 - half) * dx,
                ];
                particles.push(lattice_particle(pos, dx, field(pos), lam, cp));
            }
        }
    }
    let centre = (n / 2) * n * n + (n / 2) * n + n / 2;
    (particles, centre)
}

#[test]
fn heat_flux_of_quadratic_field_is_alpha_times_laplacian() {
    // T = 300 + x²  ⇒  ∇²T = 2  ⇒  dT/dt = 2 λ/(ρ c_p) (K/s).
    let (lam, cp) = (0.6, 4186.0);
    let alpha = lam / (RHO * cp);
    let (mut particles, centre) = centred_lattice(7, 0.05, lam, cp, |x| 300.0 + x[0] * x[0]);
    compute_heat_flux(&mut particles);
    let rate = particles[centre].d_temperature;
    assert!(
        rel_err(rate, 2.0 * alpha) < 0.03,
        "dT/dt = {rate}, expected 2 α = {}",
        2.0 * alpha
    );
}

#[test]
fn heat_flux_vanishes_for_linear_temperature_profile() {
    let (mut particles, centre) = centred_lattice(7, 0.05, 0.6, 4186.0, |x| 300.0 + 40.0 * x[1]);
    compute_heat_flux(&mut particles);
    let alpha = 0.6 / (RHO * 4186.0);
    // Reference scale: the rate a unit-curvature field would produce.
    assert!(
        particles[centre].d_temperature.abs() < 1.0e-6 * alpha,
        "∇²T = 0 must give no heating: {}",
        particles[centre].d_temperature
    );
}

#[test]
fn heat_flux_conserves_total_heat_on_a_lattice() {
    let (mut particles, _) = centred_lattice(5, 0.05, 0.6, 4186.0, |x| {
        300.0 + 25.0 * (7.0 * x[0]).sin() + 10.0 * x[2]
    });
    compute_heat_flux(&mut particles);
    let total: f64 = particles
        .iter()
        .map(|p| p.mass * p.specific_heat * p.d_temperature)
        .sum();
    let scale: f64 = particles
        .iter()
        .map(|p| (p.mass * p.specific_heat * p.d_temperature).abs())
        .sum();
    assert!(scale > 0.0);
    assert!(total.abs() < 1.0e-10 * scale, "Σ m c_p dT/dt = {total}");
}

#[test]
fn sinusoidal_temperature_mode_decays_at_the_diffusive_rate() {
    // T = T0 + A sin(k x) must decay as exp(−α k² t).  The lattice is surrounded by
    // ghost layers that are re-imposed from the analytic solution after every step,
    // so the core particles always see a complete kernel support.
    let dx = 0.05;
    let wavelength_cells = 12_i32;
    let halo = 3_i32;
    let k = 2.0 * PI / (wavelength_cells as f64 * dx);
    // α = λ/(ρ c_p) = 1e-3 m²/s;  α dt/dx² = 0.1.
    let (lam, cp) = (1000.0, 1000.0);
    let alpha = lam / (RHO * cp);
    let dt = 0.1 * dx * dx / alpha;
    let (t0, amp) = (300.0, 10.0);

    let mut particles = Vec::new();
    let mut is_core = Vec::new();
    for ix in -halo..(wavelength_cells + halo) {
        for iy in -halo..(1 + halo) {
            for iz in -halo..(1 + halo) {
                let pos = [ix as f64 * dx, iy as f64 * dx, iz as f64 * dx];
                particles.push(lattice_particle(
                    pos,
                    dx,
                    t0 + amp * (k * pos[0]).sin(),
                    lam,
                    cp,
                ));
                is_core.push((0..wavelength_cells).contains(&ix) && iy == 0 && iz == 0);
            }
        }
    }

    let mode_amplitude = |particles: &[ThermalParticle]| {
        let (mut num, mut den) = (0.0, 0.0);
        for (p, &core) in particles.iter().zip(is_core.iter()) {
            if core {
                let s = (k * p.position[0]).sin();
                num += (p.temperature - t0) * s;
                den += s * s;
            }
        }
        num / den
    };
    assert!((mode_amplitude(&particles) - amp).abs() < 1.0e-12);

    let steps = 8;
    for step in 1..=steps {
        thermal_conduction_step(&mut particles, dt);
        let decayed = amp * (-alpha * k * k * dt * step as f64).exp();
        for (p, &core) in particles.iter_mut().zip(is_core.iter()) {
            if !core {
                p.temperature = t0 + decayed * (k * p.position[0]).sin();
            }
        }
    }
    let measured = mode_amplitude(&particles);
    let t_final = dt * steps as f64;
    let measured_exponent = -(measured / amp).ln();
    let expected_exponent = alpha * k * k * t_final;
    assert!(
        measured < amp,
        "the mode must decay, not grow: {amp} -> {measured}"
    );
    assert!(
        rel_err(measured_exponent, expected_exponent) < 0.10,
        "decay exponent {measured_exponent} vs α k² t = {expected_exponent}"
    );
}

// ---------------------------------------------------------------------------
// phase_change and the enthalpy method (F4666)
// ---------------------------------------------------------------------------

#[test]
fn phase_change_synchronises_fraction_and_enthalpy_without_touching_temperature() {
    let material = PhaseChangeMaterial::water_ice();
    let mut p = ice_particle(273.0);
    phase_change(&mut p, &material);
    assert_eq!(p.temperature, 273.0);
    assert!((p.phase - 0.5).abs() < 1.0e-12, "phase = {}", p.phase);
    let expected = 2100.0 * 273.0 + 0.5 * material.latent_heat;
    assert!(
        rel_err(p.internal_energy, expected) < 1.0e-12,
        "{} vs {expected}",
        p.internal_energy
    );
}

#[test]
fn enthalpy_inversion_roundtrips_in_every_regime() {
    for material in [
        PhaseChangeMaterial::water_ice(),
        PhaseChangeMaterial::aluminium(),
    ] {
        let cp = 2100.0;
        let ts = material.solidus_temperature;
        let tl = material.liquidus_temperature;
        let (t_min, t_max) = (0.8 * ts, 1.2 * tl);
        for k in 0..=80 {
            let t = t_min + (t_max - t_min) * k as f64 / 80.0;
            let f = material.liquid_fraction(t);
            let e = material.specific_enthalpy(cp, t, f);
            let (t_back, f_back) = material.temperature_and_fraction(cp, e);
            assert!(
                (t_back - t).abs() < 1.0e-9 * t,
                "T = {t}: recovered {t_back}"
            );
            assert!((f_back - f).abs() < 1.0e-9, "T = {t}: f {f} vs {f_back}");
        }
    }
}

#[test]
fn enthalpy_is_strictly_increasing_across_the_mushy_zone() {
    let material = PhaseChangeMaterial::water_ice();
    let mut previous = f64::NEG_INFINITY;
    for k in 0..=400 {
        let t = 270.0 + 6.0 * k as f64 / 400.0;
        let e = material.specific_enthalpy(2100.0, t, material.liquid_fraction(t));
        assert!(e > previous, "e(T) must increase: T = {t}");
        previous = e;
    }
    // The mushy zone [272, 274] K holds c_p ΔT + L of enthalpy.
    let e_lo = material.specific_enthalpy(2100.0, 272.0, 0.0);
    let e_hi = material.specific_enthalpy(2100.0, 274.0, 1.0);
    assert!(rel_err(e_hi - e_lo, 2100.0 * 2.0 + 334_000.0) < 1.0e-12);
}

#[test]
fn latent_heat_holds_the_temperature_inside_the_mushy_zone() {
    // Constant heating of ice that starts at 271 K: c_p dT/dt = c_p · (1 K/s).
    let material = PhaseChangeMaterial::water_ice();
    let cp = 2100.0;
    let mut particles = vec![ice_particle(271.0)];
    phase_change(&mut particles[0], &material);
    particles[0].d_temperature = 1.0;
    let dt = 1.0;
    let mushy_capacity = cp * 2.0 + material.latent_heat; // per kg, between T_s and T_l

    // 1 s of heating brings the ice to the solidus and no further.
    integrate_enthalpy(&mut particles, dt, &material);
    assert!((particles[0].temperature - 272.0).abs() < 1.0e-9);
    assert_eq!(particles[0].phase, 0.0);

    // Every further second adds c_p J/kg, which melts material instead of heating it.
    for n in 1..=100 {
        integrate_enthalpy(&mut particles, dt, &material);
        let p = &particles[0];
        let fraction = n as f64 * cp / mushy_capacity;
        assert!(
            p.temperature > 272.0 && p.temperature < 274.0,
            "step {n}: T = {} left the mushy zone",
            p.temperature
        );
        assert!(
            (p.phase - fraction).abs() < 1.0e-9,
            "step {n}: f = {} vs {fraction}",
            p.phase
        );
        assert!(
            (p.temperature - (272.0 + 2.0 * fraction)).abs() < 1.0e-9,
            "step {n}: T = {}",
            p.temperature
        );
    }
    // After 101 s a sensible-heat-only integration would be at 372 K; the latent heat
    // has kept the ice at 273.24 K, 62 % molten.
    assert!(particles[0].temperature > 273.0 && particles[0].temperature < 273.5);
    assert!(particles[0].phase > 0.5 && particles[0].phase < 0.7);
}

#[test]
fn melting_completes_after_the_latent_heat_has_been_absorbed() {
    let material = PhaseChangeMaterial::water_ice();
    let cp = 2100.0;
    let mut particles = vec![ice_particle(272.0)];
    phase_change(&mut particles[0], &material);
    particles[0].d_temperature = 1.0;
    let e_start = particles[0].internal_energy;

    let mushy_capacity = cp * 2.0 + material.latent_heat;
    let steps_to_melt = (mushy_capacity / cp).ceil() as usize; // 162 one-second steps
    let mut steps_in_zone = 0;
    for _ in 0..steps_to_melt + 20 {
        integrate_enthalpy(&mut particles, 1.0, &material);
        if particles[0].phase < 1.0 {
            steps_in_zone += 1;
            assert!(particles[0].temperature < 274.0 + 1.0e-12);
        }
    }
    assert_eq!(steps_in_zone, steps_to_melt - 1);
    let p = &particles[0];
    assert_eq!(p.phase, 1.0);
    // Liquid: T = (e − L)/c_p with e = e_start + c_p · steps.
    let steps = (steps_to_melt + 20) as f64;
    let expected_t = (e_start + cp * steps - material.latent_heat) / cp;
    assert!(
        (p.temperature - expected_t).abs() < 1.0e-9,
        "{} vs {expected_t}",
        p.temperature
    );
    assert!(p.temperature > 274.0);
}

#[test]
fn freezing_holds_the_temperature_until_the_latent_heat_is_released() {
    let material = PhaseChangeMaterial::water_ice();
    let cp = 2100.0;
    let mut particles = vec![ice_particle(275.0)];
    phase_change(&mut particles[0], &material);
    assert_eq!(particles[0].phase, 1.0);
    particles[0].d_temperature = -1.0;

    // 1 s of cooling reaches the liquidus, then the zone is crossed with f falling.
    integrate_enthalpy(&mut particles, 1.0, &material);
    assert!((particles[0].temperature - 274.0).abs() < 1.0e-9);
    let mushy_capacity = cp * 2.0 + material.latent_heat;
    let mut previous_fraction = particles[0].phase;
    for n in 1..=120 {
        integrate_enthalpy(&mut particles, 1.0, &material);
        let p = &particles[0];
        assert!(p.temperature > 272.0 && p.temperature <= 274.0 + 1.0e-12);
        assert!(p.phase < previous_fraction, "f must fall while freezing");
        previous_fraction = p.phase;
        let expected_fraction = 1.0 - n as f64 * cp / mushy_capacity;
        assert!((p.phase - expected_fraction).abs() < 1.0e-9);
    }
}

#[test]
fn isothermal_melting_stays_at_the_melting_point_until_fully_liquid() {
    let material = PhaseChangeMaterial {
        melting_temperature: 273.15,
        latent_heat: 334_000.0,
        solidus_temperature: 273.15,
        liquidus_temperature: 273.15,
    };
    assert_eq!(material.mushy_zone_width(), 0.0);
    let cp = 2100.0;
    let mut particles = vec![ice_particle(273.15)];
    particles[0].phase = 0.0; // a solid exactly at its melting point
    particles[0].internal_energy = cp * 273.15;
    particles[0].d_temperature = 2.0; // c_p · 2 K/s

    let fraction_per_step = cp * 2.0 / material.latent_heat;
    for n in 1..=50 {
        integrate_enthalpy(&mut particles, 1.0, &material);
        let p = &particles[0];
        assert!(
            (p.temperature - 273.15).abs() < 1.0e-9,
            "step {n}: T = {}",
            p.temperature
        );
        assert!(
            (p.phase - n as f64 * fraction_per_step).abs() < 1.0e-9,
            "step {n}: f = {}",
            p.phase
        );
    }
    // Melt the remainder, then the temperature rises above the melting point.
    for _ in 0..200 {
        integrate_enthalpy(&mut particles, 1.0, &material);
    }
    assert_eq!(particles[0].phase, 1.0);
    assert!(particles[0].temperature > 273.15);
}

#[test]
fn integrate_enthalpy_adds_exactly_the_conducted_heat() {
    let material = PhaseChangeMaterial::water_ice();
    let cp = 2100.0;
    let mut particles = vec![ice_particle(272.5), ice_particle(300.0)];
    for p in particles.iter_mut() {
        phase_change(p, &material);
    }
    particles[0].d_temperature = 0.4;
    particles[1].d_temperature = -0.9;
    let before: Vec<f64> = particles.iter().map(|p| p.internal_energy).collect();
    integrate_enthalpy(&mut particles, 0.5, &material);
    assert!((particles[0].internal_energy - (before[0] + cp * 0.4 * 0.5)).abs() < 1.0e-6);
    assert!((particles[1].internal_energy - (before[1] - cp * 0.9 * 0.5)).abs() < 1.0e-6);
    // The liquid particle only changes its sensible heat.
    assert!((particles[1].temperature - (300.0 - 0.45)).abs() < 1.0e-9);
    assert_eq!(particles[1].phase, 1.0);
}

#[test]
fn integrate_enthalpy_skips_particles_without_a_usable_heat_capacity() {
    let material = PhaseChangeMaterial::water_ice();
    let mut zero_cp = ThermalParticle::new([0.0; 3], 1.0e-3, 920.0, 0.01, 272.5, 2.2, 0.0);
    zero_cp.d_temperature = 5.0;
    let mut nan_cp = ThermalParticle::new([0.0; 3], 1.0e-3, 920.0, 0.01, 272.5, 2.2, f64::NAN);
    nan_cp.d_temperature = 5.0;
    let mut particles = vec![zero_cp, nan_cp];
    integrate_enthalpy(&mut particles, 1.0, &material);
    assert_eq!(particles[0].temperature, 272.5);
    assert_eq!(particles[1].temperature, 272.5);
}

#[test]
fn integrate_enthalpy_clamps_at_absolute_zero() {
    let material = PhaseChangeMaterial::water_ice();
    let mut particles = vec![ice_particle(1.0)];
    particles[0].d_temperature = -1.0e6;
    integrate_enthalpy(&mut particles, 1.0, &material);
    assert_eq!(particles[0].temperature, 0.0);
    assert!(particles[0].internal_energy.is_finite());
    assert_eq!(particles[0].phase, 0.0);
}

#[test]
fn conduction_with_phase_change_delays_the_heating_of_ice() {
    // A cold ice particle at the solidus is heated by a hot liquid neighbour.  A
    // sensible-heat-only integration pushes it through the mushy zone within a few
    // seconds; with latent heat the heat melts ice and T stays below the liquidus.
    let material = PhaseChangeMaterial::water_ice();
    let build = || {
        let mut ice = ThermalParticle::new([0.0; 3], 1.0e-3, RHO, 0.05, 272.0, 6.0e3, 2100.0);
        let mut water =
            ThermalParticle::new([0.04, 0.0, 0.0], 1.0e-3, RHO, 0.05, 300.0, 6.0e3, 2100.0);
        phase_change(&mut ice, &material);
        phase_change(&mut water, &material);
        vec![ice, water]
    };
    let dt = 0.5;

    let mut sensible_only = build();
    let mut with_latent = build();
    let enthalpy_before: f64 = with_latent.iter().map(|p| p.mass * p.internal_energy).sum();
    for _ in 0..40 {
        thermal_conduction_step(&mut sensible_only, dt);
        thermal_conduction_step_with_phase(&mut with_latent, dt, &material);
    }
    assert!(
        sensible_only[0].temperature > 274.0,
        "without latent heat the ice warms straight through: {}",
        sensible_only[0].temperature
    );
    let ice = &with_latent[0];
    assert!(
        ice.temperature > 272.0 && ice.temperature < 274.0,
        "with latent heat the ice stays in the mushy zone: {}",
        ice.temperature
    );
    assert!(
        ice.phase > 0.0 && ice.phase < 0.5,
        "melting must have started: f = {}",
        ice.phase
    );
    assert!(with_latent[1].temperature < 300.0, "the water must cool");

    // The conducted heat is exchanged between the particles: Σ m e is conserved.
    let enthalpy_after: f64 = with_latent.iter().map(|p| p.mass * p.internal_energy).sum();
    assert!(
        rel_err(enthalpy_after, enthalpy_before) < 1.0e-12,
        "Σ m e: {enthalpy_before} -> {enthalpy_after}"
    );
}

#[test]
fn conduction_with_phase_matches_plain_step_when_no_latent_heat_is_involved() {
    // Both particles stay liquid: only the sensible heat changes, the two
    // integrators must agree.
    let material = PhaseChangeMaterial::water_ice();
    let build = || {
        let mut a = ThermalParticle::new([0.0; 3], 1.0e-3, RHO, 0.05, 320.0, 6.0e3, 4186.0);
        let mut b = ThermalParticle::new([0.04, 0.0, 0.0], 1.0e-3, RHO, 0.05, 300.0, 6.0e3, 4186.0);
        phase_change(&mut a, &material);
        phase_change(&mut b, &material);
        vec![a, b]
    };
    let mut plain = build();
    let mut phased = build();
    for _ in 0..10 {
        thermal_conduction_step(&mut plain, 0.5);
        thermal_conduction_step_with_phase(&mut phased, 0.5, &material);
    }
    for (p, q) in plain.iter().zip(phased.iter()) {
        assert!(
            (p.temperature - q.temperature).abs() < 1.0e-9,
            "{} vs {}",
            p.temperature,
            q.temperature
        );
        assert_eq!(q.phase, 1.0);
    }
}
