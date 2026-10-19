//! Tests of the Brookshaw/Morris SPH Laplacian that replaces the former "kernel value
//! only" operators: [`morris_laplacian`], [`sph_diffusion_laplacian`],
//! [`advection_diffusion_rhs`], [`energy_equation_rhs`] and the tuple-API models
//! [`SPHDiffusion`], [`SPHMultispeciesTransport`], [`SPHReactiveTransport`] and
//! [`SPHHeatTransfer`].
//!
//! Convention: `r_ij = x_i − x_j`, `∇_iW_ij = cubic_kernel_grad(r_ij, h)` and the
//! pair term `2 V_j (φ_i − φ_j)(r_ij·∇_iW_ij)/(r_ij² + 0.01 h²)`, which is negative
//! when `φ_i > φ_j`.

use super::functions::{
    advection_diffusion_rhs, cubic_kernel_grad, dot3, energy_equation_rhs, morris_laplacian,
    sph_diffusion_laplacian,
};
use super::types::{
    MultiTransportParticle, SPHDiffusion, SPHHeatTransfer, SPHMultispeciesTransport,
    SPHReactiveTransport,
};

const RHO: f64 = 1000.0;
const DX: f64 = 0.05;
/// Smoothing length of the lattice fixtures: the kernel support (`2h = 2.6 dx`) covers
/// the first four neighbour shells and the Laplacian of a quadratic field is
/// reproduced to better than 3 %.
const H: f64 = 1.3 * DX;

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

/// `n³` lattice of multi-species particles centred on the origin; species 0 holds
/// `c(x)` and the temperature is `t(x)`.  Returns the particles and the central index.
fn lattice(
    n: usize,
    c: impl Fn([f64; 3]) -> f64,
    t: impl Fn([f64; 3]) -> f64,
) -> (Vec<MultiTransportParticle>, usize) {
    let half = (n as f64 - 1.0) / 2.0;
    let mut particles = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [
                    (ix as f64 - half) * DX,
                    (iy as f64 - half) * DX,
                    (iz as f64 - half) * DX,
                ];
                let mut p = MultiTransportParticle::new(pos, [0.0; 3], 2);
                p.rho = RHO;
                p.mass = RHO * DX * DX * DX;
                p.h = H;
                p.concentration = vec![c(pos), 0.25];
                p.temperature = t(pos);
                particles.push(p);
            }
        }
    }
    let centre = (n / 2) * n * n + (n / 2) * n + n / 2;
    (particles, centre)
}

/// `(c_j, r_ij, m_j, ρ_j)` tuples of every other particle for the tuple-API models.
fn neighbour_tuples(
    particles: &[MultiTransportParticle],
    i: usize,
    value: impl Fn(&MultiTransportParticle) -> f64,
) -> Vec<(f64, [f64; 3], f64, f64)> {
    particles
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, pj)| {
            let r_ij = [
                particles[i].pos[0] - pj.pos[0],
                particles[i].pos[1] - pj.pos[1],
                particles[i].pos[2] - pj.pos[2],
            ];
            (value(pj), r_ij, pj.mass, pj.rho)
        })
        .collect()
}

// --- the pair term ------------------------------------------------------------------

#[test]
fn morris_pair_term_matches_the_closed_form_and_flows_down_the_gradient() {
    let r_ij = [-0.04, 0.0, 0.0];
    let (vol, h) = (2.0e-3, 0.05);
    let grad = cubic_kernel_grad(r_ij, h);
    let r_dot_grad = dot3(r_ij, grad);
    assert!(
        r_dot_grad < 0.0,
        "r·∇W = r dW/dr must be negative: {r_dot_grad}"
    );
    let expected = 2.0 * vol * (3.0 - 1.0) * r_dot_grad / (0.04 * 0.04 + 0.01 * h * h);
    let pair = morris_laplacian(3.0, 1.0, r_ij, vol, h);
    assert!(rel_err(pair, expected) < 1e-14, "{pair} vs {expected}");
    assert!(pair < 0.0, "φ_i > φ_j: the pair term lowers φ_i");
    assert!(morris_laplacian(1.0, 3.0, r_ij, vol, h) > 0.0);
    // Equal values do not diffuse, the term is linear in the difference and in V_j.
    assert_eq!(morris_laplacian(2.0, 2.0, r_ij, vol, h), 0.0);
    assert!(rel_err(morris_laplacian(5.0, 1.0, r_ij, vol, h), 2.0 * pair) < 1e-14);
    assert!(rel_err(morris_laplacian(3.0, 1.0, r_ij, 3.0 * vol, h), 3.0 * pair) < 1e-14);
}

#[test]
fn morris_pair_term_does_not_depend_on_the_orientation_of_the_separation() {
    let r = [0.02, -0.03, 0.01];
    let flipped = [-0.02, 0.03, -0.01];
    let a = morris_laplacian(2.0, 0.5, r, 1.0e-3, 0.06);
    let b = morris_laplacian(2.0, 0.5, flipped, 1.0e-3, 0.06);
    assert!(rel_err(a, b) < 1e-14);
}

#[test]
fn morris_pair_term_is_bounded_for_coincident_particles() {
    // The η² regulariser keeps the pair finite (and zero) as r → 0.
    for r in [1.0e-3, 1.0e-6, 1.0e-9, 0.0] {
        let pair = morris_laplacian(1.0, 0.0, [r, 0.0, 0.0], 1.0e-3, 0.05);
        assert!(pair.is_finite() && pair <= 0.0, "r = {r}: {pair}");
    }
    assert_eq!(morris_laplacian(1.0, 0.0, [0.0; 3], 1.0e-3, 0.05), 0.0);
}

#[test]
fn morris_pair_term_vanishes_outside_the_kernel_support() {
    assert_eq!(
        morris_laplacian(1.0, 0.0, [0.2, 0.0, 0.0], 1.0e-3, 0.05),
        0.0
    );
}

// --- sph_diffusion_laplacian / advection_diffusion_rhs ------------------------------

#[test]
fn diffusion_laplacian_of_a_quadratic_field_is_two_times_the_diffusivity() {
    // c = x²  ⇒  ∇²c = 2  ⇒  D ∇²c = 2 D at an interior particle (the former
    // kernel-value-only operator was not a consistent Laplacian).
    let d = 1.0e-3;
    let (particles, centre) = lattice(7, |x| x[0] * x[0], |_| 300.0);
    let value = sph_diffusion_laplacian(&particles, centre, 0, cubic_kernel_grad, d);
    assert!(
        rel_err(value, 2.0 * d) < 0.03,
        "D ∇²c = {value}, expected {}",
        2.0 * d
    );
    // The y² profile of species 0 gives the same answer in another direction.
    let (particles, centre) = lattice(7, |x| x[1] * x[1], |_| 300.0);
    let value = sph_diffusion_laplacian(&particles, centre, 0, cubic_kernel_grad, d);
    assert!(rel_err(value, 2.0 * d) < 0.03, "{value}");
}

#[test]
fn diffusion_laplacian_is_negative_at_a_peak_and_positive_in_a_trough() {
    let peak = |x: [f64; 3]| (-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / (2.0 * 0.1 * 0.1)).exp();
    let (particles, centre) = lattice(7, peak, |_| 300.0);
    assert!(sph_diffusion_laplacian(&particles, centre, 0, cubic_kernel_grad, 1.0e-3) < 0.0);
    let (particles, centre) = lattice(7, |x| -peak(x), |_| 300.0);
    assert!(sph_diffusion_laplacian(&particles, centre, 0, cubic_kernel_grad, 1.0e-3) > 0.0);
}

#[test]
fn diffusion_laplacian_pair_is_antisymmetric_and_conserves_the_scalar() {
    let mut a = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
    let mut b = MultiTransportParticle::new([0.05, 0.0, 0.0], [0.0; 3], 1);
    for p in [&mut a, &mut b] {
        p.h = 0.1;
        p.mass = 1.0;
        p.rho = RHO;
    }
    a.concentration = vec![1.0];
    b.concentration = vec![0.0];
    let particles = vec![a, b];
    let da = sph_diffusion_laplacian(&particles, 0, 0, cubic_kernel_grad, 1.0e-3);
    let db = sph_diffusion_laplacian(&particles, 1, 0, cubic_kernel_grad, 1.0e-3);
    assert!(da < 0.0 && db > 0.0, "{da}, {db}");
    assert!(
        (particles[0].mass * da + particles[1].mass * db).abs() < 1e-12 * da.abs(),
        "total scalar must be conserved: {da} vs {db}"
    );
}

#[test]
fn diffusion_laplacian_ignores_missing_species_and_invalid_indices() {
    let (particles, centre) = lattice(5, |x| x[0] * x[0], |_| 300.0);
    // A species the particles do not carry counts as concentration zero everywhere.
    assert_eq!(
        sph_diffusion_laplacian(&particles, centre, 7, cubic_kernel_grad, 1.0),
        0.0
    );
    assert_eq!(
        sph_diffusion_laplacian(&particles, particles.len(), 0, cubic_kernel_grad, 1.0),
        0.0
    );
}

#[test]
fn advection_diffusion_rhs_is_the_diffusion_operator() {
    let (particles, centre) = lattice(7, |x| x[0] * x[0], |_| 300.0);
    let rhs = advection_diffusion_rhs(&particles, centre, 0, cubic_kernel_grad, 2.0e-3);
    let lap = sph_diffusion_laplacian(&particles, centre, 0, cubic_kernel_grad, 2.0e-3);
    assert_eq!(rhs, lap);
    assert!(rel_err(rhs, 4.0e-3) < 0.03);
}

// --- energy_equation_rhs ------------------------------------------------------------

#[test]
fn energy_equation_rhs_is_alpha_times_the_temperature_laplacian() {
    // T = x²  ⇒  dT/dt = α ∇²T = 2 α.  The argument is the thermal diffusivity
    // α = λ/(ρ c_v), not the conductivity.
    let alpha = 1.4e-7;
    let (particles, centre) = lattice(7, |_| 0.0, |x| 300.0 + x[0] * x[0]);
    let rate = energy_equation_rhs(&particles, centre, cubic_kernel_grad, alpha);
    assert!(rel_err(rate, 2.0 * alpha) < 0.03, "dT/dt = {rate}");
}

#[test]
fn energy_equation_rhs_cools_a_hot_spot_and_heats_a_cold_one() {
    let spot = |x: [f64; 3]| (-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / (2.0 * 0.1 * 0.1)).exp();
    let (hot, centre) = lattice(7, |_| 0.0, |x| 300.0 + 100.0 * spot(x));
    assert!(energy_equation_rhs(&hot, centre, cubic_kernel_grad, 1.0e-6) < 0.0);
    let (cold, centre) = lattice(7, |_| 0.0, |x| 300.0 - 100.0 * spot(x));
    assert!(energy_equation_rhs(&cold, centre, cubic_kernel_grad, 1.0e-6) > 0.0);
    // Uniform temperature: no heat flow.
    let (flat, centre) = lattice(5, |_| 0.0, |_| 321.0);
    assert!(energy_equation_rhs(&flat, centre, cubic_kernel_grad, 1.0e-6).abs() < 1e-18);
}

#[test]
fn energy_equation_rhs_pair_conserves_heat() {
    let mut a = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
    let mut b = MultiTransportParticle::new([0.05, 0.0, 0.0], [0.0; 3], 1);
    for p in [&mut a, &mut b] {
        p.h = 0.1;
        p.mass = 1.0;
        p.rho = RHO;
    }
    a.temperature = 400.0;
    b.temperature = 300.0;
    let particles = vec![a, b];
    let ra = energy_equation_rhs(&particles, 0, cubic_kernel_grad, 1.0e-3);
    let rb = energy_equation_rhs(&particles, 1, cubic_kernel_grad, 1.0e-3);
    assert!(ra < 0.0 && rb > 0.0);
    assert!((ra + rb).abs() < 1e-12 * ra.abs());
}

// --- tuple-API models ---------------------------------------------------------------

#[test]
fn sph_diffusion_model_reproduces_the_laplacian_of_a_quadratic_field() {
    let (particles, centre) = lattice(7, |x| x[2] * x[2], |_| 300.0);
    let neighbours = neighbour_tuples(&particles, centre, |p| p.concentration[0]);
    let c_i = particles[centre].concentration[0];
    let model = SPHDiffusion::new(1.0e-3, H);
    let lap = model.concentration_laplacian(c_i, &neighbours);
    assert!(rel_err(lap, 2.0) < 0.03, "∇²c = {lap}");
    assert!(rel_err(model.dcdt(c_i, &neighbours), 2.0e-3) < 0.03);
}

#[test]
fn sph_diffusion_model_is_negative_at_a_concentration_peak() {
    let model = SPHDiffusion::new(1.0, H);
    let neighbours = [
        (0.0, [-0.05, 0.0, 0.0], 1.0, RHO),
        (0.0, [0.05, 0.0, 0.0], 1.0, RHO),
    ];
    assert!(model.concentration_laplacian(1.0, &neighbours) < 0.0);
    // The pair term depends on r_ij only through |r_ij| and r_ij·∇W, so reversing
    // every separation vector leaves the result unchanged.
    let reversed = [
        (0.0, [0.05, 0.0, 0.0], 1.0, RHO),
        (0.0, [-0.05, 0.0, 0.0], 1.0, RHO),
    ];
    assert_eq!(
        model.concentration_laplacian(1.0, &neighbours),
        model.concentration_laplacian(1.0, &reversed)
    );
}

#[test]
fn sph_diffusion_model_skips_neighbours_without_density() {
    let model = SPHDiffusion::new(1.0, H);
    let neighbours = [(0.0, [0.05, 0.0, 0.0], 1.0, 0.0)];
    assert_eq!(model.concentration_laplacian(1.0, &neighbours), 0.0);
}

#[test]
fn multispecies_transport_diffuses_each_species_with_its_own_coefficient() {
    let (particles, centre) = lattice(7, |x| x[0] * x[0], |_| 300.0);
    let neighbours = neighbour_tuples(&particles, centre, |p| p.concentration[0]);
    let c_i = particles[centre].concentration[0];
    let model =
        SPHMultispeciesTransport::new(vec![1.0e-3, 4.0e-3], vec![-1.0, 1.0], vec![1.0, 1.0]);
    let slow = model.dcdt(0, c_i, &neighbours, H, 0.0);
    let fast = model.dcdt(1, c_i, &neighbours, H, 0.0);
    assert!(rel_err(slow, 2.0e-3) < 0.03, "{slow}");
    assert!(rel_err(fast, 8.0e-3) < 0.03, "{fast}");
    // The reaction source is added on top of the diffusion.
    let with_source = model.dcdt(0, c_i, &neighbours, H, 0.5);
    assert!((with_source - slow - 0.5).abs() < 1e-12);
}

#[test]
fn reactive_transport_combines_diffusion_decay_and_source() {
    let (particles, centre) = lattice(7, |x| x[1] * x[1], |_| 300.0);
    let neighbours = neighbour_tuples(&particles, centre, |p| p.concentration[0]);
    let c_i = particles[centre].concentration[0];
    let model = SPHReactiveTransport::new(2.0e-3, 0.4, H);
    let diffusion_only = SPHReactiveTransport::new(2.0e-3, 0.0, H).dcdt(c_i, &neighbours, 0.0);
    assert!(rel_err(diffusion_only, 4.0e-3) < 0.03, "{diffusion_only}");
    let full = model.dcdt(c_i, &neighbours, 0.25);
    assert!((full - (diffusion_only - 0.4 * c_i + 0.25)).abs() < 1e-12);
}

#[test]
fn heat_transfer_model_conducts_down_the_temperature_gradient() {
    let (k, cp) = (0.6, 4182.0);
    let alpha = k / (RHO * cp);
    let model = SPHHeatTransfer::new(k, cp, RHO);
    // T = x²: dT/dt = α ∇²T = 2 α.
    let (particles, centre) = lattice(7, |_| 0.0, |x| 300.0 + x[0] * x[0]);
    let neighbours = neighbour_tuples(&particles, centre, |p| p.temperature);
    let t_i = particles[centre].temperature;
    let rate = model.dtemp_dt(t_i, &neighbours, H);
    assert!(rel_err(rate, 2.0 * alpha) < 0.03, "dT/dt = {rate}");

    // A single hot-cold pair: the hot particle loses heat, the cold one gains it.
    let r_ij = [0.05, 0.0, 0.0];
    let hot = model.heat_exchange(400.0, 300.0, 1.0, RHO, r_ij, 0.1);
    let cold = model.heat_exchange(300.0, 400.0, 1.0, RHO, [-0.05, 0.0, 0.0], 0.1);
    assert!(hot < 0.0 && cold > 0.0, "{hot}, {cold}");
    assert!((hot + cold).abs() < 1e-12 * hot.abs());
    // heat_exchange is the conduction power density k ∇²T of the pair: ρ c_p α times
    // the temperature rate of the same pair.
    let pair_rate = model.dtemp_dt(400.0, &[(300.0, r_ij, 1.0, RHO)], 0.1);
    assert!(
        rel_err(hot, RHO * cp * pair_rate) < 1e-12,
        "{hot} vs {}",
        RHO * cp * pair_rate
    );
    // Neighbours without density are skipped.
    assert_eq!(model.heat_exchange(400.0, 300.0, 1.0, 0.0, r_ij, 0.1), 0.0);
}

#[test]
fn heat_transfer_model_scales_with_conductivity() {
    let r_ij = [0.05, 0.0, 0.0];
    let one = SPHHeatTransfer::new(1.0, 1.0, 1.0).heat_exchange(400.0, 300.0, 1.0, 1.0, r_ij, 0.5);
    let three =
        SPHHeatTransfer::new(3.0, 1.0, 1.0).heat_exchange(400.0, 300.0, 1.0, 1.0, r_ij, 0.5);
    assert!(rel_err(three, 3.0 * one) < 1e-14);
}
