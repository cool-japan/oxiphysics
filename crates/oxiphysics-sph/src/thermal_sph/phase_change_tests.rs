//! Tests for the enthalpy-method phase change of `thermal_sph`: the
//! temperature <-> enthalpy mapping and its use by `ThermalSphSimulation::step`.

use super::functions::{enthalpy_from_temp, liquid_fraction_from_temp, temp_from_enthalpy};
use super::types::{PhaseChangeSph, ThermalSphSimulation};
use super::types_ext::ThermalParticle;

const T_MELT: f64 = 273.15;
const T_LIQUIDUS: f64 = 273.65;
const LATENT: f64 = 334_000.0;
const CP: f64 = 4182.0;
/// Particle mass: 1000 kg/m³ in a (0.01 m)³ cell.
const MASS: f64 = 1.0e-3;

fn mushy_model() -> PhaseChangeSph {
    PhaseChangeSph::new(T_MELT, T_LIQUIDUS, LATENT, CP)
}

fn isothermal_model() -> PhaseChangeSph {
    PhaseChangeSph::new(T_MELT, T_MELT, LATENT, CP)
}

/// Two water-like particles one smoothing length (0.01 m) apart.  With the
/// default conductivity 0.6 W/(m K) the pair relaxes with a rate constant of
/// about `7e-4 / s`, so `dt = 10 s` is a stable explicit step.
fn two_particles(t_a: f64, t_b: f64, dt: f64) -> ThermalSphSimulation {
    let mut sim = ThermalSphSimulation::new(0.01, dt);
    sim.add_particle(ThermalParticle::new(
        [0.0, 0.0, 0.0],
        [0.0; 3],
        MASS,
        1000.0,
        t_a,
    ));
    sim.add_particle(ThermalParticle::new(
        [0.01, 0.0, 0.0],
        [0.0; 3],
        MASS,
        1000.0,
        t_b,
    ));
    sim
}

// --- temperature <-> enthalpy mapping ----------------------------------------

#[test]
fn liquid_fraction_regimes_and_isothermal_convention() {
    assert_eq!(liquid_fraction_from_temp(200.0, T_MELT, T_LIQUIDUS), 0.0);
    assert_eq!(liquid_fraction_from_temp(T_MELT, T_MELT, T_LIQUIDUS), 0.0);
    assert_eq!(
        liquid_fraction_from_temp(T_LIQUIDUS, T_MELT, T_LIQUIDUS),
        1.0
    );
    assert_eq!(liquid_fraction_from_temp(400.0, T_MELT, T_LIQUIDUS), 1.0);
    let mid = liquid_fraction_from_temp(T_MELT + 0.125, T_MELT, T_LIQUIDUS);
    assert!((mid - 0.25).abs() < 1e-12, "{mid}");
    // Isothermal transition: the melting temperature itself is solid.
    assert_eq!(liquid_fraction_from_temp(T_MELT, T_MELT, T_MELT), 0.0);
    assert_eq!(
        liquid_fraction_from_temp(T_MELT + 1e-9, T_MELT, T_MELT),
        1.0
    );
}

#[test]
fn enthalpy_matches_the_closed_form_in_each_regime() {
    let pc = mushy_model();
    assert_eq!(pc.enthalpy(200.0), CP * 200.0);
    assert!((pc.enthalpy(T_MELT) - CP * T_MELT).abs() < 1e-9);
    assert!((pc.enthalpy(T_LIQUIDUS) - (CP * T_LIQUIDUS + LATENT)).abs() < 1e-6);
    let mid = 0.5 * (T_MELT + T_LIQUIDUS);
    assert!((pc.enthalpy(mid) - (CP * mid + 0.5 * LATENT)).abs() < 1e-5);
    assert!((pc.enthalpy(400.0) - (CP * 400.0 + LATENT)).abs() < 1e-9);
    assert_eq!(
        pc.enthalpy(310.0),
        enthalpy_from_temp(310.0, CP, T_MELT, T_LIQUIDUS, LATENT)
    );
}

#[test]
fn enthalpy_is_continuous_at_the_solidus_and_the_liquidus() {
    // The former H(T) jumped by c_p (T_liquidus - T_melt) = 2091 J/kg at the liquidus.
    let pc = mushy_model();
    let eps = 1.0e-9;
    for edge in [T_MELT, T_LIQUIDUS] {
        let below = pc.enthalpy(edge - eps);
        let at = pc.enthalpy(edge);
        let above = pc.enthalpy(edge + eps);
        assert!(
            (at - below).abs() < 1.0e-2 && (above - at).abs() < 1.0e-2,
            "H jumps at T = {edge}: {below} / {at} / {above}"
        );
    }
}

#[test]
fn temperature_enthalpy_round_trip_in_every_regime() {
    let pc = mushy_model();
    for t in [
        100.0, 273.0, T_MELT, 273.2, 273.4, 273.6, T_LIQUIDUS, 273.7, 300.0, 1000.0,
    ] {
        let back = pc.temperature_from_enthalpy(pc.enthalpy(t));
        assert!((back - t).abs() < 1e-9, "T = {t}: round trip gave {back}");
    }
}

#[test]
fn enthalpy_temperature_round_trip_is_continuous_and_monotone() {
    let pc = mushy_model();
    let h_start = pc.enthalpy(250.0);
    let h_end = pc.enthalpy(300.0);
    let count = 20_000_u32;
    let mut previous = pc.temperature_from_enthalpy(h_start);
    for k in 1..=count {
        let h = h_start + (h_end - h_start) * f64::from(k) / f64::from(count);
        let t = pc.temperature_from_enthalpy(h);
        assert!(
            t >= previous - 1e-12,
            "T(H) must be non-decreasing: {previous} -> {t} at H = {h}"
        );
        assert!(
            t - previous < 0.05,
            "T(H) must not jump: {previous} -> {t} at H = {h}"
        );
        let h_back = pc.enthalpy(t);
        assert!((h_back - h).abs() < 1e-5, "H = {h}: H(T(H)) = {h_back}");
        previous = t;
    }
}

#[test]
fn mushy_zone_has_the_apparent_heat_capacity_cp_plus_l_over_delta_t() {
    let pc = mushy_model();
    let h_solidus = CP * T_MELT;
    let apparent = CP + LATENT / (T_LIQUIDUS - T_MELT);
    let t = pc.temperature_from_enthalpy(h_solidus + 1000.0);
    assert!(
        (t - (T_MELT + 1000.0 / apparent)).abs() < 1e-12,
        "T = {t}, expected {}",
        T_MELT + 1000.0 / apparent
    );
}

#[test]
fn narrow_melting_range_round_trips_without_losing_digits() {
    let width = 1.0e-6;
    let pc = PhaseChangeSph::new(T_MELT, T_MELT + width, LATENT, CP);
    for frac in [0.0, 0.1, 0.5, 0.9, 1.0] {
        let t = T_MELT + frac * width;
        let back = pc.temperature_from_enthalpy(pc.enthalpy(t));
        assert!((back - t).abs() < 1e-11, "T = {t}: round trip gave {back}");
    }
    // Half of the latent heat: resolved far below the width of the melting range.
    let h = CP * T_MELT + 0.5 * (CP * width + LATENT);
    let t = pc.temperature_from_enthalpy(h);
    assert!(
        (t - (T_MELT + 0.5 * width)).abs() < 1e-9,
        "T = {t} at half the latent heat"
    );
    assert!((pc.liquid_fraction_from_enthalpy(h) - 0.5).abs() < 1e-9);
}

#[test]
fn isothermal_transition_has_an_exact_temperature_plateau() {
    let pc = isothermal_model();
    let h_solidus = CP * T_MELT;
    assert_eq!(pc.enthalpy(T_MELT), h_solidus);
    assert_eq!(pc.liquid_fraction(T_MELT), 0.0);
    let just_above = T_MELT + 1.0e-9;
    assert!((pc.enthalpy(just_above) - (CP * just_above + LATENT)).abs() < 1e-9);

    for frac in [1.0e-6, 0.25, 0.5, 0.75, 1.0 - 1.0e-6] {
        let h = h_solidus + frac * LATENT;
        assert_eq!(pc.temperature_from_enthalpy(h), T_MELT, "plateau at {frac}");
        assert!((pc.liquid_fraction_from_enthalpy(h) - frac).abs() < 1e-9);
    }
    let below = pc.temperature_from_enthalpy(h_solidus - 10.0 * CP);
    assert!((below - (T_MELT - 10.0)).abs() < 1e-9, "{below}");
    let above = pc.temperature_from_enthalpy(h_solidus + LATENT + 10.0 * CP);
    assert!((above - (T_MELT + 10.0)).abs() < 1e-9, "{above}");
    assert_eq!(pc.liquid_fraction_from_enthalpy(h_solidus - 1.0), 0.0);
    assert_eq!(
        pc.liquid_fraction_from_enthalpy(h_solidus + LATENT + 1.0),
        1.0
    );
}

#[test]
fn enthalpy_based_liquid_fraction_matches_the_temperature_based_one() {
    let pc = mushy_model();
    for k in 0..=50 {
        let t = 273.0 + 0.015 * f64::from(k);
        let by_temperature = pc.liquid_fraction(t);
        let by_enthalpy = pc.liquid_fraction_from_enthalpy(pc.enthalpy(t));
        assert!(
            (by_temperature - by_enthalpy).abs() < 1e-9,
            "T = {t}: {by_temperature} vs {by_enthalpy}"
        );
    }
}

#[test]
fn temp_from_enthalpy_without_latent_heat_is_plain_sensible_heating() {
    for h in [1.0e5, CP * T_MELT, CP * 273.4, CP * T_LIQUIDUS, 2.0e6] {
        let t = temp_from_enthalpy(h, CP, T_MELT, T_LIQUIDUS, 0.0);
        assert!((t - h / CP).abs() < 1e-9, "H = {h}: T = {t}");
    }
}

// --- time integration --------------------------------------------------------

#[test]
fn isothermal_melting_holds_the_melting_point_until_the_latent_heat_is_absorbed() {
    let t_hot = 373.15;
    let dt = 10.0;
    let mut sim = two_particles(T_MELT, t_hot, dt).with_phase_change(isothermal_model());
    let mut plain = two_particles(T_MELT, t_hot, dt);
    let enthalpy_before = sim.total_enthalpy();

    // Without the phase-change model the solid particle heats up immediately.
    plain.step();
    assert!(
        plain.particles[0].temperature > T_MELT + 0.1,
        "plain model: {}",
        plain.particles[0].temperature
    );

    let mut previous_fraction = 0.0;
    let mut steps_on_plateau = 0_usize;
    let mut finished = false;
    for _ in 0..20_000 {
        sim.step();
        let t_a = sim.particles[0].temperature;
        let fraction = sim.liquid_fraction(0).expect("phase change is active");
        if t_a > T_MELT {
            finished = true;
            assert!(fraction >= 1.0 - 1e-12, "melting is complete: {fraction}");
            break;
        }
        assert_eq!(t_a, T_MELT, "the solid/liquid particle must sit at T_melt");
        assert!(
            fraction >= previous_fraction && fraction < 1.0,
            "melt fraction must grow monotonically: {previous_fraction} -> {fraction}"
        );
        // Heat balance: what the hot particle lost went into the latent heat.
        let absorbed = CP * (t_hot - sim.particles[1].temperature);
        assert!(
            (absorbed - fraction * LATENT).abs() < 1.0e-6 * LATENT,
            "absorbed {absorbed} J/kg vs f_l L = {}",
            fraction * LATENT
        );
        previous_fraction = fraction;
        steps_on_plateau += 1;
    }
    assert!(finished, "the particle never finished melting");
    assert!(
        steps_on_plateau > 100,
        "plateau lasted only {steps_on_plateau} steps"
    );

    // Total enthalpy is conserved; the hot particle paid L / c_p (to within one step).
    let enthalpy_after = sim.total_enthalpy();
    assert!(
        (enthalpy_after - enthalpy_before).abs() < 1e-9 * enthalpy_before.abs(),
        "enthalpy drifted: {enthalpy_before} -> {enthalpy_after}"
    );
    let drop = t_hot - sim.particles[1].temperature;
    assert!(
        drop * CP >= LATENT * (1.0 - 1e-9) && drop * CP < LATENT + 1.0e4,
        "the hot particle gave up {} J/kg, expected about {LATENT}",
        drop * CP
    );
    assert!(sim.particles[0].temperature - T_MELT < 1.0);
}

#[test]
fn mushy_zone_dwell_is_much_longer_than_without_phase_change() {
    let t_hot = 373.15;
    let dt = 10.0;
    let mut with_pc = two_particles(T_MELT, t_hot, dt).with_phase_change(mushy_model());
    let mut plain = two_particles(T_MELT, t_hot, dt);

    let mut plain_steps = 0_usize;
    while plain.particles[0].temperature <= T_LIQUIDUS && plain_steps < 20_000 {
        plain.step();
        plain_steps += 1;
    }
    assert!((1..=3).contains(&plain_steps), "plain: {plain_steps}");

    let mut previous = T_MELT;
    let mut pc_steps = 0_usize;
    while with_pc.particles[0].temperature <= T_LIQUIDUS && pc_steps < 20_000 {
        with_pc.step();
        pc_steps += 1;
        let t_a = with_pc.particles[0].temperature;
        if t_a <= T_LIQUIDUS {
            assert!(
                t_a >= T_MELT && t_a >= previous - 1e-12,
                "the temperature must stay in the melting range and not fall: {previous} -> {t_a}"
            );
            previous = t_a;
        }
    }
    assert!(
        pc_steps < 20_000,
        "the particle never left the melting range"
    );
    assert!(
        pc_steps >= 100 * plain_steps,
        "with latent heat: {pc_steps} steps, without: {plain_steps}"
    );
}

#[test]
fn solidification_releases_latent_heat_and_dwells_in_the_melting_range() {
    let t_cold = 150.0;
    let dt = 10.0;
    let mut sim = two_particles(T_LIQUIDUS, t_cold, dt).with_phase_change(mushy_model());
    let mut plain = two_particles(T_LIQUIDUS, t_cold, dt);
    let enthalpy_before = sim.total_enthalpy();

    plain.step();
    assert!(
        plain.particles[0].temperature < T_MELT,
        "plain model: the particle drops below the melting point at once ({})",
        plain.particles[0].temperature
    );

    let mut previous_fraction = 1.0;
    for step in 0..50 {
        sim.step();
        let t_a = sim.particles[0].temperature;
        let fraction = sim.liquid_fraction(0).expect("phase change is active");
        assert!(
            (T_MELT..=T_LIQUIDUS).contains(&t_a),
            "step {step}: T = {t_a} left the melting range"
        );
        assert!(
            fraction < previous_fraction,
            "step {step}: the liquid fraction must fall: {previous_fraction} -> {fraction}"
        );
        previous_fraction = fraction;
    }
    assert!(
        previous_fraction > 0.5,
        "latent heat release is slow, f_l = {previous_fraction}"
    );
    let enthalpy_after = sim.total_enthalpy();
    assert!(
        (enthalpy_after - enthalpy_before).abs() < 1e-11 * enthalpy_before.abs(),
        "enthalpy drifted: {enthalpy_before} -> {enthalpy_after}"
    );
    assert!(sim.particles[1].temperature > t_cold);
}

#[test]
fn step_conserves_total_enthalpy_for_identical_particles() {
    let temperatures = [255.0, 268.0, 273.4, 279.0, 300.0];
    let mut sim = ThermalSphSimulation::new(0.01, 5.0).with_phase_change(mushy_model());
    for (k, &t) in temperatures.iter().enumerate() {
        sim.add_particle(ThermalParticle::new(
            [0.01 * k as f64, 0.0, 0.0],
            [0.0; 3],
            MASS,
            1000.0,
            t,
        ));
    }
    let before = sim.total_enthalpy();
    let spread_before = sim.max_temperature() - sim.min_temperature();
    for _ in 0..400 {
        sim.step();
    }
    let after = sim.total_enthalpy();
    assert!(
        (after - before).abs() < 1e-11 * before.abs(),
        "total enthalpy drifted: {before} -> {after}"
    );
    let spread_after = sim.max_temperature() - sim.min_temperature();
    assert!(
        spread_after < spread_before,
        "conduction must smooth the profile: {spread_before} -> {spread_after}"
    );
    assert!((sim.time - 400.0 * 5.0).abs() < 1e-9);
}

#[test]
fn externally_modified_temperatures_and_new_particles_are_resynchronised() {
    let mut sim = two_particles(T_MELT, 373.15, 10.0).with_phase_change(isothermal_model());
    for _ in 0..10 {
        sim.step();
    }
    assert_eq!(sim.particles[0].temperature, T_MELT);
    let melted = sim.liquid_fraction(0).expect("phase change is active");
    assert!(melted > 0.0 && melted < 1.0, "partially molten: {melted}");

    // A boundary condition resets the hot particle: the stale enthalpy must not be used.
    sim.particles[1].temperature = 500.0;
    sim.step();
    let t_b = sim.particles[1].temperature;
    assert!(t_b > 495.0 && t_b < 500.0, "T_B = {t_b}");

    // A particle appended directly to the public vector is picked up as well.
    sim.particles.push(ThermalParticle::new(
        [0.02, 0.0, 0.0],
        [0.0; 3],
        MASS,
        1000.0,
        290.0,
    ));
    sim.step();
    let t_c = sim.particles[2].temperature;
    assert!(t_c.is_finite() && (t_c - 290.0).abs() < 5.0, "T_C = {t_c}");

    // The plateau particle keeps its melt fraction through all of this.
    assert_eq!(sim.particles[0].temperature, T_MELT);
    assert!(sim.liquid_fraction(0).expect("active") > melted);
}

#[test]
fn disabled_phase_change_reproduces_the_plain_temperature_update() {
    let mut disabled = two_particles(T_MELT, 373.15, 10.0).with_phase_change(mushy_model());
    disabled.enable_phase_change = false;
    let mut plain = two_particles(T_MELT, 373.15, 10.0);
    for _ in 0..25 {
        disabled.step();
        plain.step();
    }
    for (a, b) in disabled.particles.iter().zip(plain.particles.iter()) {
        assert_eq!(a.temperature, b.temperature);
    }
    assert_eq!(disabled.liquid_fraction(0), None);
    assert_eq!(disabled.specific_enthalpy(0), None);
    assert_eq!(disabled.total_enthalpy(), disabled.total_thermal_energy());
}

#[test]
fn particles_away_from_the_melting_range_behave_like_the_plain_model() {
    for (t_a, t_b) in [(350.0, 380.0), (200.0, 250.0)] {
        let mut with_pc = two_particles(t_a, t_b, 10.0).with_phase_change(mushy_model());
        let mut plain = two_particles(t_a, t_b, 10.0);
        // Particle 0 uses a specific heat that differs from the model's.
        with_pc.particles[0].specific_heat = 900.0;
        plain.particles[0].specific_heat = 900.0;
        for _ in 0..100 {
            with_pc.step();
            plain.step();
        }
        for (a, b) in with_pc.particles.iter().zip(plain.particles.iter()) {
            assert!(
                (a.temperature - b.temperature).abs() < 1e-8,
                "T = {} with latent heat vs {} without",
                a.temperature,
                b.temperature
            );
        }
        assert!(with_pc.particles[0].temperature != t_a, "heat must flow");
    }
}

#[test]
fn enthalpy_accessors_report_sensible_plus_latent_heat() {
    let sim = two_particles(200.0, 400.0, 1.0).with_phase_change(mushy_model());
    assert_eq!(sim.liquid_fraction(0), Some(0.0));
    assert_eq!(sim.liquid_fraction(1), Some(1.0));
    assert_eq!(sim.liquid_fraction(2), None);
    assert_eq!(sim.specific_enthalpy(0), Some(CP * 200.0));
    let h_liquid = sim.specific_enthalpy(1).expect("particle 1 exists");
    assert!((h_liquid - (CP * 400.0 + LATENT)).abs() < 1e-9);
    assert_eq!(sim.specific_enthalpy(5), None);

    let expected = MASS * (CP * 200.0) + MASS * (CP * 400.0 + LATENT);
    assert!((sim.total_enthalpy() - expected).abs() < 1e-9 * expected);
    // The sensible-only energy deliberately excludes the latent heat.
    let sensible = sim.total_thermal_energy();
    assert!((sensible - MASS * CP * 600.0).abs() < 1e-9 * sensible);

    let without_model = two_particles(200.0, 400.0, 1.0);
    assert_eq!(without_model.liquid_fraction(0), None);
    assert_eq!(without_model.specific_enthalpy(0), None);
}
