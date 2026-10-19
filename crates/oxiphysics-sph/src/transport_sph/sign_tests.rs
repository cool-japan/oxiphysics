//! Direction, conservation and consistency tests for the SPH transport operators.
//!
//! All Brookshaw/Morris style Laplacians in this module follow the convention
//! `r_ij = x_i − x_j` with `∇_iW_ij = (dW/dr) r_ij/|r_ij|`, so `r_ij·∇_iW_ij < 0`
//! and the diffusive pair term is `2 V_j (φ_i − φ_j) (r_ij·∇_iW_ij) / (r²_ij + η²)`.
//! These tests pin the *physical* behaviour (scalar and heat flow down the
//! gradient, magnitude of `∇²` on a quadratic field, pairwise conservation)
//! rather than finiteness.

use super::{
    AdvectionSPH, AdvectionScheme, AdvectionSchemeKind, DiffusionSph, FickianDiffusion,
    GrayScottSph, HeatTransferSPH, ReactionDiffusion, SoluteTransport, ThermalConduction,
    TransportParticle, cubic_kernel_grad,
};

/// Particle on the x axis with `m = 1 kg`, `ρ = 1000 kg/m³` (`V = 1e-3 m³`) and `h = 0.1 m`.
fn particle(x: f64, c: f64, t: f64) -> TransportParticle {
    TransportParticle::new([x, 0.0, 0.0], [0.0; 3], 1.0, 1000.0, 0.1, c, t)
}

/// Two particles 0.05 m apart (`q = r/h = 0.5`).
fn pair(c0: f64, c1: f64, t0: f64, t1: f64) -> Vec<TransportParticle> {
    vec![particle(0.0, c0, t0), particle(0.05, c1, t1)]
}

/// Three particles at x = 0, 0.05, 0.1 with a concentration peak in the middle.
fn peak_triple() -> Vec<TransportParticle> {
    vec![
        particle(0.0, 0.0, 300.0),
        particle(0.05, 1.0, 300.0),
        particle(0.1, 0.0, 300.0),
    ]
}

/// Cubic lattice of `n³` particles (spacing `dx`, `h = 1.25 dx`, `m = ρ dx³`)
/// centred on the origin; `field(x)` supplies the (concentration, temperature)
/// sample.  Returns the particles and the index of the central one.
fn lattice(
    n: usize,
    dx: f64,
    rho: f64,
    field: impl Fn([f64; 3]) -> (f64, f64),
) -> (Vec<TransportParticle>, usize) {
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
                let (c, t) = field(pos);
                particles.push(TransportParticle::new(
                    pos,
                    [0.0; 3],
                    rho * dx * dx * dx,
                    rho,
                    1.25 * dx,
                    c,
                    t,
                ));
            }
        }
    }
    let centre = (n / 2) * n * n + (n / 2) * n + n / 2;
    (particles, centre)
}

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

// --- DiffusionSph -----------------------------------------------------------

#[test]
fn diffusion_pair_flows_down_gradient_and_conserves_scalar() {
    let particles = pair(0.0, 1.0, 300.0, 300.0);
    let rates = DiffusionSph::new(1.0e-3).compute_dphidt(&particles);
    assert!(rates[0] > 0.0, "c = 0 particle must gain: {}", rates[0]);
    assert!(rates[1] < 0.0, "c = 1 particle must lose: {}", rates[1]);
    let flux_0 = particles[0].volume() * rates[0];
    let flux_1 = particles[1].volume() * rates[1];
    assert!(
        (flux_0 + flux_1).abs() < 1e-12 * flux_0.abs(),
        "pair exchange must be antisymmetric: {flux_0} vs {flux_1}"
    );
}

#[test]
fn diffusion_rate_magnitude_matches_brookshaw_pair_formula() {
    // dφ_0/dt = 2 D V_1 (φ_0 − φ_1)(r_01·∇W_01) / (r² + 0.01 h²), with r_01 = (−0.05, 0, 0).
    let d = 2.0e-3;
    let particles = pair(0.0, 1.0, 300.0, 300.0);
    let rates = DiffusionSph::new(d).compute_dphidt(&particles);
    let r_01 = [-0.05, 0.0, 0.0];
    let grad_w = cubic_kernel_grad(r_01, 0.1);
    let r_dot_grad = r_01[0] * grad_w[0];
    assert!(
        r_dot_grad < 0.0,
        "r_ij·∇_iW_ij must be negative: {r_dot_grad}"
    );
    let expected = 2.0 * d * 1.0e-3 * (0.0 - 1.0) * r_dot_grad / (0.05 * 0.05 + 0.01 * 0.1 * 0.1);
    assert!(
        rel_err(rates[0], expected) < 1e-12,
        "rate {} vs Brookshaw pair formula {expected}",
        rates[0]
    );
}

#[test]
fn diffusion_laplacian_of_quadratic_field_is_positive_two() {
    // c = x²  =>  ∇²c = 2  =>  dc/dt = 2 D at an interior particle.
    let d = 1.0e-3;
    let (particles, centre) = lattice(7, 0.05, 1000.0, |x| (x[0] * x[0], 300.0));
    let rates = DiffusionSph::new(d).compute_dphidt(&particles);
    assert!(
        rel_err(rates[centre], 2.0 * d) < 0.08,
        "dc/dt = {} but D ∇²c = {}",
        rates[centre],
        2.0 * d
    );
}

#[test]
fn diffusion_step_smooths_gaussian_peak_and_conserves_scalar() {
    let sigma = 0.1_f64;
    let mut particles: Vec<TransportParticle> = (0..11)
        .map(|i| {
            let x = (i as f64 - 5.0) * 0.05;
            particle(x, (-x * x / (2.0 * sigma * sigma)).exp(), 300.0)
        })
        .collect();
    let total_before: f64 = particles.iter().map(|p| p.volume() * p.conc()).sum();
    let peak_before = particles[5].conc();
    let max_before = particles.iter().map(|p| p.conc()).fold(f64::MIN, f64::max);
    DiffusionSph::new(1.0e-2).step(&mut particles, 0.01);
    let total_after: f64 = particles.iter().map(|p| p.volume() * p.conc()).sum();
    let max_after = particles.iter().map(|p| p.conc()).fold(f64::MIN, f64::max);
    assert!(
        particles[5].conc() < peak_before,
        "the peak must decay: {peak_before} -> {}",
        particles[5].conc()
    );
    assert!(max_after < max_before, "maximum must not grow");
    assert!(
        rel_err(total_after, total_before) < 1e-12,
        "total scalar must be conserved: {total_before} -> {total_after}"
    );
}

// --- ThermalConduction ------------------------------------------------------

#[test]
fn thermal_conduction_pair_heats_cold_particle_and_conserves_energy() {
    let mut particles = pair(0.5, 0.5, 200.0, 400.0);
    let tc = ThermalConduction::new(1.0, 1000.0);
    let rates = tc.compute_dtdt(&particles);
    assert!(rates[0] > 0.0, "200 K particle must warm: {}", rates[0]);
    assert!(rates[1] < 0.0, "400 K particle must cool: {}", rates[1]);

    let energy_before = tc.total_thermal_energy(&particles);
    let gap_before = particles[1].temperature - particles[0].temperature;
    tc.step(&mut particles, 0.001);
    let energy_after = tc.total_thermal_energy(&particles);
    let gap_after = particles[1].temperature - particles[0].temperature;
    assert!(
        rel_err(energy_after, energy_before) < 1e-12,
        "thermal energy must be conserved: {energy_before} -> {energy_after}"
    );
    assert!(
        gap_after > 0.0 && gap_after < gap_before,
        "temperature gap must shrink without overshoot: {gap_before} -> {gap_after}"
    );
}

#[test]
fn thermal_conduction_laplacian_of_quadratic_field() {
    // T = x²  =>  dT/dt = k/(ρ c_v) ∇²T = 2 k/(ρ c_v).
    let k = 0.6;
    let rho_cv = 4.18e6;
    let (particles, centre) = lattice(7, 0.05, 1000.0, |x| (0.0, x[0] * x[0]));
    let rates = ThermalConduction::new(k, rho_cv).compute_dtdt(&particles);
    let expected = 2.0 * k / rho_cv;
    assert!(
        rel_err(rates[centre], expected) < 0.08,
        "dT/dt = {} but α ∇²T = {expected}",
        rates[centre]
    );
}

// --- FickianDiffusion -------------------------------------------------------

#[test]
fn fickian_laplacian_negative_at_concentration_peak() {
    let particles = peak_triple();
    let fd = FickianDiffusion::isotropic(1.0e-3);
    let centre = fd.sph_laplacian(&particles, 1, cubic_kernel_grad);
    let left = fd.sph_laplacian(&particles, 0, cubic_kernel_grad);
    let right = fd.sph_laplacian(&particles, 2, cubic_kernel_grad);
    assert!(
        centre < 0.0,
        "Laplacian at a peak must be negative: {centre}"
    );
    assert!(
        left > 0.0 && right > 0.0,
        "flanks must gain: {left}, {right}"
    );
    assert!(
        (left - right).abs() < 1e-12 * left.abs(),
        "symmetric configuration: {left} vs {right}"
    );
}

#[test]
fn fickian_laplacian_of_quadratic_field() {
    let d = 5.0e-4;
    let (particles, centre) = lattice(7, 0.05, 1000.0, |x| (x[1] * x[1], 300.0));
    let lap = FickianDiffusion::isotropic(d).sph_laplacian(&particles, centre, cubic_kernel_grad);
    assert!(
        rel_err(lap, 2.0 * d) < 0.08,
        "D ∇²c = {lap}, expected {}",
        2.0 * d
    );
}

#[test]
fn fickian_laplacian_respects_tortuosity() {
    let particles = peak_triple();
    let plain = FickianDiffusion::isotropic(1.0e-3).sph_laplacian(&particles, 1, cubic_kernel_grad);
    let mut tortuous = FickianDiffusion::isotropic(1.0e-3);
    tortuous.tortuosity = 2.0;
    let slowed = tortuous.sph_laplacian(&particles, 1, cubic_kernel_grad);
    assert!(
        rel_err(slowed, 0.5 * plain) < 1e-12,
        "tortuosity 2 must halve D_eff: {slowed} vs {plain}"
    );
}

// --- AdvectionSPH / AdvectionScheme ----------------------------------------

#[test]
fn advection_sph_artificial_diffusion_damps_peak() {
    let mut particles = peak_triple();
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let adv = AdvectionSPH::new(0.01, 0.5, false);
    let centre = adv.advect_scalar(&particles, 1, super::cubic_kernel, cubic_kernel_grad);
    let left = adv.advect_scalar(&particles, 0, super::cubic_kernel, cubic_kernel_grad);
    let right = adv.advect_scalar(&particles, 2, super::cubic_kernel, cubic_kernel_grad);
    assert!(
        centre < 0.0,
        "artificial diffusion must lower the peak: {centre}"
    );
    assert!(
        left > 0.0 && right > 0.0,
        "flanks must rise: {left}, {right}"
    );
}

#[test]
fn advection_sph_back_tracing_transports_a_ramp_downstream() {
    // c = 2x with the flow along +x and fixed particles: ∂c/∂t = −v ∂c/∂x = −2 per
    // second.  Back-tracing to x − v Δt (Shepard interpolation with the cubic kernel)
    // reproduces this to within a few per cent on an interior particle.
    let (mut particles, centre) = lattice(7, 0.05, 1000.0, |x| (2.0 * x[0], 300.0));
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let adv = AdvectionSPH::new(0.01, 0.0, true);
    let rate = adv.advect_scalar(&particles, centre, super::cubic_kernel, cubic_kernel_grad);
    assert!(
        rel_err(rate, -2.0) < 0.03,
        "back-traced rate {rate} vs −v ∂c/∂x = −2"
    );
    // The same ramp advected against the flow gains scalar.
    for p in &mut particles {
        p.velocity = [-1.0, 0.0, 0.0];
    }
    let reversed = adv.advect_scalar(&particles, centre, super::cubic_kernel, cubic_kernel_grad);
    assert!(
        rel_err(reversed, 2.0) < 0.03,
        "reversed flow rate {reversed} vs +2"
    );
}

#[test]
fn advection_sph_without_semi_lagrangian_flag_has_no_convective_term() {
    // Lagrangian particles carry the scalar with them: without artificial diffusion the
    // rate vanishes even though the velocity and the gradient are non-zero.
    let (mut particles, centre) = lattice(7, 0.05, 1000.0, |x| (2.0 * x[0], 300.0));
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let lagrangian = AdvectionSPH::new(0.01, 0.0, false);
    let rate = lagrangian.advect_scalar(&particles, centre, super::cubic_kernel, cubic_kernel_grad);
    assert_eq!(rate, 0.0);
    // The flag matters: the back-traced scheme sees the ramp.
    let traced = AdvectionSPH::new(0.01, 0.0, true);
    let traced_rate =
        traced.advect_scalar(&particles, centre, super::cubic_kernel, cubic_kernel_grad);
    assert!(traced_rate < -1.0);
    // With artificial diffusion a Lagrangian particle only feels the Laplacian of a
    // quadratic profile, which is positive (the profile is convex) for c = x².
    let (mut convex, centre) = lattice(7, 0.05, 1000.0, |x| (x[0] * x[0], 300.0));
    for p in &mut convex {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let diffusive = AdvectionSPH::new(0.01, 1.0, false);
    let smoothing =
        diffusive.advect_scalar(&convex, centre, super::cubic_kernel, cubic_kernel_grad);
    // ε_art |v| h ∇²c = 1 · 1 · 0.0625 · 2
    assert!(
        rel_err(smoothing, 0.125) < 0.08,
        "artificial diffusion rate {smoothing} vs 0.125"
    );
}

#[test]
fn upwind_scheme_adds_diffusion_relative_to_standard() {
    let mut particles = peak_triple();
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let standard = AdvectionScheme::new(AdvectionSchemeKind::Standard, 0.0);
    let upwind = AdvectionScheme::new(AdvectionSchemeKind::Upwind, 0.5);
    let no_eps = AdvectionScheme::new(AdvectionSchemeKind::Upwind, 0.0);
    let r_std = standard.compute_rate_i(&particles, 1);
    let r_up = upwind.compute_rate_i(&particles, 1);
    let r_zero = no_eps.compute_rate_i(&particles, 1);
    assert!(
        (r_zero - r_std).abs() < 1e-15 * r_std.abs().max(1.0),
        "ε = 0 upwind must equal the standard scheme: {r_zero} vs {r_std}"
    );
    assert!(
        r_up - r_std < -1.0e-6,
        "stabilising term must lower the peak: upwind {r_up}, standard {r_std}"
    );
}

#[test]
fn upwind_stabilisation_only_acts_on_downwind_neighbours() {
    // v = +x: only the neighbour with v·r_ij < 0 (the one at larger x) contributes.
    let mut particles = peak_triple();
    for p in &mut particles {
        p.velocity = [1.0, 0.0, 0.0];
    }
    let eps = 0.5;
    let standard = AdvectionScheme::new(AdvectionSchemeKind::Standard, 0.0);
    let upwind = AdvectionScheme::new(AdvectionSchemeKind::Upwind, eps);
    let extra = upwind.compute_rate_i(&particles, 1) - standard.compute_rate_i(&particles, 1);
    // r_12 = x_1 − x_2 = (−0.05, 0, 0) for the downwind neighbour (particle 2).
    let r_12 = [-0.05, 0.0, 0.0];
    let grad = cubic_kernel_grad(r_12, 0.1);
    let expected =
        eps * 1.0e-3 * 2.0 * (1.0 - 0.0) * (r_12[0] * grad[0]) / (0.05 * 0.05 + 0.01 * 0.1 * 0.1);
    assert!(
        rel_err(extra, expected) < 1e-12,
        "upwind increment {extra} vs single downwind-neighbour term {expected}"
    );
}

// --- HeatTransferSPH --------------------------------------------------------

#[test]
fn heat_transfer_pair_flows_hot_to_cold_and_conserves_energy() {
    let ht = HeatTransferSPH::new(0.6, 4182.0);
    let particles = pair(0.0, 0.0, 200.0, 400.0);
    let r_cold = ht.heat_conduction_rate(&particles, 0, cubic_kernel_grad);
    let r_hot = ht.heat_conduction_rate(&particles, 1, cubic_kernel_grad);
    assert!(r_cold > 0.0, "cold particle must warm: {r_cold}");
    assert!(r_hot < 0.0, "hot particle must cool: {r_hot}");
    let q_cold = particles[0].mass * r_cold;
    let q_hot = particles[1].mass * r_hot;
    assert!(
        (q_cold + q_hot).abs() < 1e-12 * q_cold.abs(),
        "m0 dT0/dt = -m1 dT1/dt expected: {q_cold} vs {q_hot}"
    );
}

#[test]
fn heat_transfer_matches_thermal_conduction_for_same_diffusivity() {
    let k = 0.6;
    let cp = 4182.0;
    let rho = 1000.0;
    let particles = pair(0.0, 0.0, 200.0, 400.0);
    let ht = HeatTransferSPH::new(k, cp);
    let tc = ThermalConduction::new(k, rho * cp);
    let dtdt = tc.compute_dtdt(&particles);
    for (i, expected) in dtdt.iter().enumerate() {
        let rate = ht.heat_conduction_rate(&particles, i, cubic_kernel_grad);
        assert!(
            rel_err(rate, *expected) < 1e-12,
            "particle {i}: HeatTransferSPH {rate} vs ThermalConduction {expected}"
        );
    }
}

#[test]
fn heat_transfer_laplacian_of_quadratic_field_has_correct_conductivity() {
    // Guards against a doubled conductivity: dT/dt = 2 k /(ρ c_p) for T = x².
    let k = 0.6;
    let cp = 4182.0;
    let (particles, centre) = lattice(7, 0.05, 1000.0, |x| (0.0, x[0] * x[0]));
    let rate =
        HeatTransferSPH::new(k, cp).heat_conduction_rate(&particles, centre, cubic_kernel_grad);
    let expected = 2.0 * k / (1000.0 * cp);
    assert!(
        rel_err(rate, expected) < 0.08,
        "dT/dt = {rate}, expected {expected}"
    );
}

#[test]
fn heat_transfer_per_particle_conductivity_uses_harmonic_mean() {
    let cp = 4182.0;
    let ht = HeatTransferSPH::new(1.0, cp);
    let particles = pair(0.0, 0.0, 200.0, 400.0);
    let ks = [1.0, 3.0];
    let r0 = ht.heat_conduction_rate_with_k(&particles, 0, &ks, cubic_kernel_grad);
    let r1 = ht.heat_conduction_rate_with_k(&particles, 1, &ks, cubic_kernel_grad);
    // Pair conductivity 4 k0 k1/(k0+k1) = 3 = 2 k_eff  =>  k_eff = 1.5.
    let reference = HeatTransferSPH::new(1.5, cp);
    let e0 = reference.heat_conduction_rate(&particles, 0, cubic_kernel_grad);
    let e1 = reference.heat_conduction_rate(&particles, 1, cubic_kernel_grad);
    assert!(rel_err(r0, e0) < 1e-12, "{r0} vs {e0}");
    assert!(rel_err(r1, e1) < 1e-12, "{r1} vs {e1}");
    // Energy exchange stays antisymmetric for a conductivity jump.
    assert!((particles[0].mass * r0 + particles[1].mass * r1).abs() < 1e-12 * r0.abs());
}

#[test]
fn heat_transfer_missing_conductivity_entries_fall_back_to_uniform_value() {
    let ht = HeatTransferSPH::new(0.8, 4182.0);
    let particles = pair(0.0, 0.0, 250.0, 350.0);
    let uniform = ht.heat_conduction_rate(&particles, 0, cubic_kernel_grad);
    let fallback = ht.heat_conduction_rate_with_k(&particles, 0, &[], cubic_kernel_grad);
    assert!(rel_err(fallback, uniform) < 1e-14);
}

#[test]
fn heat_transfer_heat_source_adds_uniform_heating() {
    let mut ht = HeatTransferSPH::new(0.6, 4182.0);
    ht.heat_source = 4.182e6; // W/m³  =>  1 K/s for ρ = 1000, c_p = 4182
    let particles = vec![particle(0.0, 0.0, 300.0)];
    let rate = ht.heat_conduction_rate(&particles, 0, cubic_kernel_grad);
    assert!((rate - 1.0).abs() < 1e-12, "rate = {rate}");
}

// --- Derived operators that inherit the Laplacian ---------------------------

#[test]
fn reaction_diffusion_spreads_scalar_down_gradient() {
    let rd = ReactionDiffusion::new(1.0e-3, 0.0);
    let particles = pair(0.2, 0.8, 300.0, 300.0);
    let rates = rd.compute_dphidt_fisher(&particles);
    assert!(rates[0] > 0.0 && rates[1] < 0.0, "{rates:?}");

    let mut stepped = pair(0.0, 1.0, 300.0, 300.0);
    rd.step_custom(&mut stepped, 0.01, |_| 0.0);
    assert!(stepped[0].conc() > 0.0 && stepped[1].conc() < 1.0);
}

#[test]
fn gray_scott_diffusion_moves_species_down_their_gradients() {
    let gs = GrayScottSph::new(1.0e-3, 1.0e-3, 0.0, 0.0);
    // u gradient with v = 0: pure diffusion of u.
    let mut particles = pair(0.0, 0.0, 300.0, 300.0);
    particles[0].concentration = vec![0.2, 0.0];
    particles[1].concentration = vec![0.8, 0.0];
    let (du, dv) = gs.compute_rates(&particles);
    assert!(du[0] > 0.0 && du[1] < 0.0, "u must diffuse: {du:?}");
    assert!(dv.iter().all(|x| x.abs() < 1e-15), "v untouched: {dv:?}");

    // v gradient with u = 0: reaction terms vanish, pure diffusion of v.
    particles[0].concentration = vec![0.0, 0.0];
    particles[1].concentration = vec![0.0, 0.6];
    let (du, dv) = gs.compute_rates(&particles);
    assert!(dv[0] > 0.0 && dv[1] < 0.0, "v must diffuse: {dv:?}");
    assert!(du.iter().all(|x| x.abs() < 1e-15), "u untouched: {du:?}");
}

#[test]
fn solute_transport_species_diffuse_independently_with_cross_terms() {
    let st = SoluteTransport::diagonal(&[1.0e-3, 1.0e-3]);
    let mut particles = pair(0.0, 0.0, 300.0, 300.0);
    particles[0].concentration = vec![0.0, 1.0];
    particles[1].concentration = vec![1.0, 0.0];
    let rates = st.compute_all_rates(&particles);
    assert!(
        rates[0][0] > 0.0 && rates[0][1] < 0.0,
        "species 0: {:?}",
        rates[0]
    );
    assert!(
        rates[1][0] < 0.0 && rates[1][1] > 0.0,
        "species 1: {:?}",
        rates[1]
    );

    // Cross diffusion: species 0 is driven by the gradient of species 1.
    let cross = SoluteTransport::new(2, vec![0.0, 1.0e-3, 0.0, 0.0]);
    let cross_rates = cross.compute_all_rates(&particles);
    assert!(
        cross_rates[0][0] < 0.0 && cross_rates[0][1] > 0.0,
        "species 0 follows the gradient of species 1: {:?}",
        cross_rates[0]
    );
    assert!(cross_rates[1].iter().all(|x| x.abs() < 1e-15));
}
