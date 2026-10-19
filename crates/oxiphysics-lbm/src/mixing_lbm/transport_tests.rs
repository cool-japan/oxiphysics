//! Transport regression tests of the mixing lattices: periodic streaming,
//! scalar conservation, lattice diffusivity, advection by the flow, momentum
//! conservation of the fluid and cross-diffusion between species.

use super::functions::{
    D2Q9_CX, D2Q9_CY, D2Q9_W, d2q9_equilibrium, d2q9_moments, laplacian_d2q9_periodic,
    lattice_diffusivity, stream_d2q9_periodic_flat, sync_populations,
};
use super::types::{DoubleDistributionFunction, MixingLBM, MultiComponentLbm, PassiveScalarField};
use crate::error::Error;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Zeroth, first and second spatial moments of a row-major scalar field.
struct FieldMoments {
    mass: f64,
    centroid_x: f64,
    centroid_y: f64,
    variance_x: f64,
    variance_y: f64,
}

fn field_moments(field: &[f64], nx: usize, ny: usize) -> FieldMoments {
    let mut mass = 0.0;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let c = field[j * nx + i];
            mass += c;
            sum_x += c * i as f64;
            sum_y += c * j as f64;
        }
    }
    let centroid_x = sum_x / mass;
    let centroid_y = sum_y / mass;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let c = field[j * nx + i];
            var_x += c * (i as f64 - centroid_x).powi(2);
            var_y += c * (j as f64 - centroid_y).powi(2);
        }
    }
    FieldMoments {
        mass,
        centroid_x,
        centroid_y,
        variance_x: var_x / mass,
        variance_y: var_y / mass,
    }
}

/// Gaussian blob `exp(-r² / (2σ²))` on an `n × n` grid.
fn gaussian_field(n: usize, cx: f64, cy: f64, sigma: f64) -> Vec<f64> {
    let mut field = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let r2 = (i as f64 - cx).powi(2) + (j as f64 - cy).powi(2);
            field[j * n + i] = (-r2 / (2.0 * sigma * sigma)).exp();
        }
    }
    field
}

/// Exact growth of the variance (per direction) of a blob that starts at
/// equilibrium in a D2Q9 BGK diffusion lattice at rest and evolves for `t`
/// steps: `2 D t − 2 c_s² τ (τ − 1) (1 − (1 − 1/τ)^t)`.  The transient term
/// comes from the first-moment (flux) relaxation; for `τ = 1` it vanishes.
fn exact_variance_growth(tau: f64, steps: usize) -> f64 {
    let d = lattice_diffusivity(tau);
    let r = 1.0 - 1.0 / tau;
    2.0 * d * steps as f64 - 2.0 * (1.0 / 3.0) * tau * (tau - 1.0) * (1.0 - r.powi(steps as i32))
}

/// Grid size, step count, blob position and width of the transport tests: the
/// blob stays far enough from the periodic seam for the wrap-around to change
/// its variance by less than 1e-5.
const GRID: usize = 48;
const STEPS: usize = 40;
const CENTRE: f64 = 24.0;
const SIGMA: f64 = 2.0;

/// Peclet number that gives the scalar relaxation time `tau` in [`MixingLBM`].
fn pe_for_tau(tau: f64) -> f64 {
    1.0 / (6.0 * (tau - 0.5))
}

// ---------------------------------------------------------------------------
// Streaming helper
// ---------------------------------------------------------------------------

/// Every population of every direction moves exactly one link and wraps around
/// at all four edges and corners, also on grids as small as 2 × 3.
#[test]
fn streaming_moves_each_population_one_link_with_wrap_around() {
    for &(nx, ny) in &[(5usize, 4usize), (2, 3), (1, 1), (3, 1)] {
        let sources = [
            (0, 0),
            (nx - 1, 0),
            (0, ny - 1),
            (nx - 1, ny - 1),
            (nx / 2, ny / 2),
        ];
        for &(i0, j0) in &sources {
            for q in 0..9 {
                let mut src = vec![0.0; 9 * nx * ny];
                src[9 * (j0 * nx + i0) + q] = 1.0;
                let mut dst = vec![f64::NAN; 9 * nx * ny];
                stream_d2q9_periodic_flat(&src, &mut dst, nx, ny).expect("lengths match");

                let ti = (i0 as i64 + D2Q9_CX[q] as i64).rem_euclid(nx as i64) as usize;
                let tj = (j0 as i64 + D2Q9_CY[q] as i64).rem_euclid(ny as i64) as usize;
                let target = 9 * (tj * nx + ti) + q;
                for (idx, &value) in dst.iter().enumerate() {
                    let expected = if idx == target { 1.0 } else { 0.0 };
                    assert_eq!(
                        value, expected,
                        "{nx}x{ny}, source ({i0},{j0}) q={q}: entry {idx}"
                    );
                }
            }
        }
    }
}

/// The two named cases of the audit: `q = 1` at `(0, 0)` lands on `(1, 0)`,
/// `q = 5` on `(1, 1)`, and both wrap at the far edge.
#[test]
fn streaming_east_and_north_east_examples() {
    let (nx, ny) = (4, 3);
    let node = |i: usize, j: usize, q: usize| 9 * (j * nx + i) + q;
    let mut src = vec![0.0; 9 * nx * ny];
    src[node(0, 0, 1)] = 1.0;
    src[node(0, 0, 5)] = 2.0;
    src[node(nx - 1, ny - 1, 1)] = 3.0;
    src[node(nx - 1, ny - 1, 5)] = 4.0;
    let mut dst = vec![0.0; src.len()];
    stream_d2q9_periodic_flat(&src, &mut dst, nx, ny).expect("lengths match");
    assert_eq!(dst[node(1, 0, 1)], 1.0);
    assert_eq!(dst[node(1, 1, 5)], 2.0);
    assert_eq!(dst[node(0, ny - 1, 1)], 3.0);
    assert_eq!(dst[node(0, 0, 5)], 4.0);
}

#[test]
fn streaming_is_a_permutation_that_conserves_every_direction() {
    let (nx, ny) = (7, 5);
    let src: Vec<f64> = (0..9 * nx * ny)
        .map(|k| (k as f64 * 0.37).sin() + 2.0)
        .collect();
    let mut dst = vec![0.0; src.len()];
    stream_d2q9_periodic_flat(&src, &mut dst, nx, ny).expect("lengths match");
    for q in 0..9 {
        let before: f64 = src.iter().skip(q).step_by(9).sum();
        let after: f64 = dst.iter().skip(q).step_by(9).sum();
        assert!((before - after).abs() < 1e-12, "direction {q}");
    }
    let mut sorted_src = src.clone();
    let mut sorted_dst = dst.clone();
    sorted_src.sort_by(f64::total_cmp);
    sorted_dst.sort_by(f64::total_cmp);
    assert_eq!(sorted_src, sorted_dst);
}

#[test]
fn streaming_rejects_mismatched_buffers_without_writing() {
    let src = vec![1.0; 9 * 6];
    let mut dst = vec![7.0; 9 * 6];
    assert!(matches!(
        stream_d2q9_periodic_flat(&src, &mut dst, 3, 3),
        Err(Error::IncompatibleParameters { .. })
    ));
    assert!(dst.iter().all(|&x| x == 7.0), "destination was modified");
    let mut short = vec![0.0; 9 * 6 - 1];
    assert!(stream_d2q9_periodic_flat(&src, &mut short, 3, 2).is_err());
    let mut empty: Vec<f64> = Vec::new();
    assert!(stream_d2q9_periodic_flat(&[], &mut empty, 0, 4).is_ok());
}

// ---------------------------------------------------------------------------
// Node helpers
// ---------------------------------------------------------------------------

#[test]
fn equilibrium_has_the_expected_moments() {
    let (rho, ux, uy) = (1.3, 0.04, -0.03);
    let feq: Vec<f64> = (0..9).map(|q| d2q9_equilibrium(q, rho, ux, uy)).collect();
    let (r, u, v) = d2q9_moments(&feq);
    assert!((r - rho).abs() < 1e-14);
    assert!((u - ux).abs() < 1e-14);
    assert!((v - uy).abs() < 1e-14);
    // Second moment: rho (c_s² + u²).
    let pxx: f64 = feq
        .iter()
        .enumerate()
        .map(|(q, f)| f * D2Q9_CX[q].powi(2))
        .sum();
    let pxy: f64 = feq
        .iter()
        .enumerate()
        .map(|(q, f)| f * D2Q9_CX[q] * D2Q9_CY[q])
        .sum();
    assert!((pxx - rho * (1.0 / 3.0 + ux * ux)).abs() < 1e-14);
    assert!((pxy - rho * ux * uy).abs() < 1e-14);
    // The weights themselves sum to one.
    assert!((D2Q9_W.iter().sum::<f64>() - 1.0).abs() < 1e-15);
}

#[test]
fn moments_of_a_vacuum_node_have_zero_velocity() {
    let (rho, ux, uy) = d2q9_moments(&[0.0; 9]);
    assert_eq!((rho, ux, uy), (0.0, 0.0, 0.0));
}

#[test]
fn sync_moves_only_the_missing_mass_onto_the_weights() {
    let mut pops: Vec<f64> = (0..9).map(|q| 0.1 * (q as f64 + 1.0)).collect();
    pops.extend((0..9).map(|q| D2Q9_W[q] * 2.0));
    let before = pops.clone();
    // Node 0 gets a new target, node 1 is already consistent.
    let consistent: f64 = pops[9..].iter().sum();
    sync_populations(&mut pops, &[5.0, consistent]);
    let node0: f64 = pops[..9].iter().sum();
    assert!((node0 - 5.0).abs() < 1e-14);
    let missing = 5.0 - before[..9].iter().sum::<f64>();
    for q in 0..9 {
        assert!((pops[q] - (before[q] + D2Q9_W[q] * missing)).abs() < 1e-14);
    }
    assert_eq!(
        &pops[9..],
        &before[9..],
        "a consistent node must not change"
    );
}

#[test]
fn laplacian_matches_second_derivatives() {
    let (nx, ny) = (16, 12);
    // Constant field: zero.
    let constant = [3.0_f64; 16 * 12];
    assert!(
        laplacian_d2q9_periodic(&constant, nx, ny)
            .iter()
            .all(|&l| l.abs() < 1e-14)
    );
    // sin(kx) is an eigenfunction with eigenvalue 2 (cos k − 1).
    let k = 2.0 * std::f64::consts::PI / nx as f64;
    let field: Vec<f64> = (0..nx * ny).map(|n| (k * (n % nx) as f64).sin()).collect();
    let lap = laplacian_d2q9_periodic(&field, nx, ny);
    for (f, l) in field.iter().zip(lap.iter()) {
        assert!((l - 2.0 * (k.cos() - 1.0) * f).abs() < 1e-13);
    }
    // x² + 2y² away from the periodic seam: ∇² = 6.
    let quad: Vec<f64> = (0..nx * ny)
        .map(|n| {
            let (x, y) = ((n % nx) as f64, (n / nx) as f64);
            x * x + 2.0 * y * y
        })
        .collect();
    let lap = laplacian_d2q9_periodic(&quad, nx, ny);
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            assert!((lap[j * nx + i] - 6.0).abs() < 1e-12, "({i},{j})");
        }
    }
    // The periodic Laplacian sums to zero; a wrong-sized field gives zeros.
    assert!(lap.iter().sum::<f64>().abs() < 1e-10);
    assert!(
        laplacian_d2q9_periodic(&[1.0, 2.0], nx, ny)
            .iter()
            .all(|&l| l == 0.0)
    );
}

// ---------------------------------------------------------------------------
// MixingLBM
// ---------------------------------------------------------------------------

#[test]
fn mixing_lbm_starts_as_a_unit_density_fluid_at_rest() {
    let m = MixingLBM::new(6, 5, 10.0, 2.0);
    for node in m.f_u.chunks_exact(9) {
        assert_eq!(node, &D2Q9_W[..], "f_u must be the rest equilibrium w_q");
        let (rho, ux, uy) = d2q9_moments(node);
        assert!((rho - 1.0).abs() < 1e-15);
        assert_eq!((ux, uy), (0.0, 0.0));
    }
    assert!(m.f_c.iter().all(|&p| p == 0.0));
    assert_eq!(m.velocity_at(0, 0), Some((0.0, 0.0)));
    assert_eq!(m.velocity_at(6, 0), None);
}

#[test]
fn mixing_lbm_relaxation_times_follow_the_dimensionless_numbers() {
    let m = MixingLBM::new(4, 4, 10.0, 2.0);
    assert!((m.tau_c() - (0.5 + 1.0 / 60.0)).abs() < 1e-15);
    assert!((m.tau_u() - (0.5 + 1.0 / 12.0)).abs() < 1e-15);
    assert!((m.scalar_diffusivity() - (1.0 / 3.0) * (m.tau_c() - 0.5)).abs() < 1e-15);
    assert!((m.flow_viscosity() - (1.0 / 3.0) * (m.tau_u() - 0.5)).abs() < 1e-15);
    assert!((pe_for_tau(0.8) - 1.0 / 1.8).abs() < 1e-12);
}

/// `init_*` writes the scalar populations `w_q c` — the previous code left
/// `f_c` at zero, so the first step turned `c` into `c / τ_c`.
#[test]
fn mixing_lbm_init_sets_the_scalar_populations() {
    let mut m = MixingLBM::new(12, 10, 10.0, 1.0);
    m.init_gaussian_blob(6.0, 5.0, 2.0, 1.5);
    for (c, node) in m.c_scalar.iter().zip(m.f_c.chunks_exact(9)) {
        for (q, &p) in node.iter().enumerate() {
            assert!((p - D2Q9_W[q] * c).abs() < 1e-15);
        }
    }
    m.init_stripe_concentration(3);
    for (c, node) in m.c_scalar.iter().zip(m.f_c.chunks_exact(9)) {
        assert!((node.iter().sum::<f64>() - c).abs() < 1e-15);
    }
}

#[test]
fn mixing_lbm_first_step_keeps_the_scalar() {
    let mut m = MixingLBM::new(16, 16, 10.0, 1.0);
    assert!((m.tau_c() - 1.0).abs() > 0.1, "τ_c must differ from 1");
    m.init_gaussian_blob(8.0, 8.0, 2.0, 1.0);
    let before = m.total_scalar();
    m.step();
    let after = m.total_scalar();
    assert!(
        (after - before).abs() < 1e-12 * before,
        "total scalar {before} -> {after}"
    );
}

#[test]
fn mixing_lbm_conserves_the_scalar_over_many_steps() {
    let mut m = MixingLBM::new(32, 32, 10.0, 1.0);
    m.init_gaussian_blob(16.0, 16.0, 2.5, 1.0);
    let before = m.total_scalar();
    for _ in 0..100 {
        m.step();
    }
    let after = m.total_scalar();
    assert!(
        (after - before).abs() < 1e-12 * before,
        "total scalar {before} -> {after}"
    );
    let from_populations: f64 = m.f_c.iter().sum();
    assert!((from_populations - after).abs() < 1e-10);
}

/// The fluid at rest stays at rest with unit density (the previous code
/// relaxed `f_u` toward weights that summed to 4/3).
#[test]
fn mixing_lbm_fluid_at_rest_stays_at_rest() {
    let mut m = MixingLBM::new(10, 8, 10.0, 5.0);
    m.init_stripe_concentration(2);
    for _ in 0..25 {
        m.step();
    }
    for node in m.f_u.chunks_exact(9) {
        let (rho, ux, uy) = d2q9_moments(node);
        assert!((rho - 1.0).abs() < 1e-13, "ρ = {rho}");
        assert!(ux.abs() < 1e-13 && uy.abs() < 1e-13, "u = ({ux}, {uy})");
    }
}

/// A uniform flow is a fixed point and carries the density along unchanged.
#[test]
fn mixing_lbm_uniform_flow_persists() {
    let mut m = MixingLBM::new(12, 9, 10.0, 5.0);
    m.set_uniform_velocity(0.05, -0.02);
    for _ in 0..40 {
        m.step();
    }
    for j in 0..9 {
        for i in 0..12 {
            let (ux, uy) = m.velocity_at(i, j).expect("inside the grid");
            assert!((ux - 0.05).abs() < 1e-12 && (uy + 0.02).abs() < 1e-12);
        }
    }
    let rho: f64 = m.f_u.iter().sum();
    assert!((rho - 12.0 * 9.0).abs() < 1e-10);
}

/// The scalar spreads with `D = c_s² (τ_c − ½)`: for a blob that starts at
/// equilibrium the variance growth is the exact lattice result.
#[test]
fn mixing_lbm_variance_grows_with_the_lattice_diffusivity() {
    let n = GRID;
    let steps = STEPS;
    for &tau in &[0.8, 1.0] {
        let mut m = MixingLBM::new(n, n, pe_for_tau(tau), 1.0);
        assert!((m.tau_c() - tau).abs() < 1e-12);
        m.init_gaussian_blob(CENTRE, CENTRE, SIGMA, 1.0);
        let start = field_moments(&m.c_scalar, n, n);
        for _ in 0..steps {
            m.step();
        }
        let end = field_moments(&m.c_scalar, n, n);
        let expected = exact_variance_growth(tau, steps);
        assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
        assert!((end.centroid_x - start.centroid_x).abs() < 1e-5);
        assert!((end.centroid_y - start.centroid_y).abs() < 1e-5);
        assert!(
            (end.variance_x - start.variance_x - expected).abs() < 1e-5,
            "τ={tau}: Δvar_x = {}, expected {expected}",
            end.variance_x - start.variance_x
        );
        assert!(
            (end.variance_y - start.variance_y - expected).abs() < 1e-5,
            "τ={tau}: Δvar_y = {}, expected {expected}",
            end.variance_y - start.variance_y
        );
        // 2 D t is the leading term.
        let two_dt = 2.0 * m.scalar_diffusivity() * steps as f64;
        assert!((expected - two_dt).abs() < 0.2);
    }
}

/// The scalar is carried by the flow: the centroid moves by `u t` and the
/// variance still grows with the lattice diffusivity.
#[test]
fn mixing_lbm_advects_the_scalar_with_the_flow() {
    let n = GRID;
    let steps = STEPS;
    let tau = 0.8;
    let mut m = MixingLBM::new(n, n, pe_for_tau(tau), 5.0);
    m.set_uniform_velocity(0.05, 0.02);
    m.init_gaussian_blob(CENTRE - 4.0, CENTRE - 2.0, SIGMA, 1.0);
    let start = field_moments(&m.c_scalar, n, n);
    for _ in 0..steps {
        m.step();
    }
    let end = field_moments(&m.c_scalar, n, n);
    assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
    assert!(
        (end.centroid_x - start.centroid_x - 0.05 * steps as f64).abs() < 1e-5,
        "x centroid moved by {}",
        end.centroid_x - start.centroid_x
    );
    assert!(
        (end.centroid_y - start.centroid_y - 0.02 * steps as f64).abs() < 1e-5,
        "y centroid moved by {}",
        end.centroid_y - start.centroid_y
    );
    let expected = exact_variance_growth(tau, steps);
    assert!((end.variance_x - start.variance_x - expected).abs() < 1e-4);
    assert!((end.variance_y - start.variance_y - expected).abs() < 1e-4);
}

/// Stripes diffuse into each other: the variance decays monotonically and the
/// total is conserved.
#[test]
fn mixing_lbm_stripes_mix() {
    let (nx, ny) = (32, 6);
    let mut m = MixingLBM::new(nx, ny, 1.0, 1.0);
    m.init_stripe_concentration(2);
    let total = m.total_scalar();
    let variance = |m: &MixingLBM| {
        let mean = m.total_scalar() / m.c_scalar.len() as f64;
        m.c_scalar.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / m.c_scalar.len() as f64
    };
    let mut previous = variance(&m);
    let initial = previous;
    for _ in 0..8 {
        for _ in 0..25 {
            m.step();
        }
        let now = variance(&m);
        assert!(
            now < previous,
            "variance did not decrease: {previous} -> {now}"
        );
        previous = now;
    }
    assert!(
        previous < 0.9 * initial,
        "mixing too slow: {initial} -> {previous}"
    );
    assert!((m.total_scalar() - total).abs() < 1e-11 * total);
}

/// Editing `c_scalar` directly is allowed: the edit is kept and conserved.
#[test]
fn mixing_lbm_direct_concentration_edits_are_conserved() {
    let mut m = MixingLBM::new(10, 10, 5.0, 1.0);
    let idx = m.index(4, 5);
    m.c_scalar[idx] = 2.0;
    for _ in 0..30 {
        m.step();
    }
    assert!((m.total_scalar() - 2.0).abs() < 1e-12);
    assert!(m.c_scalar[idx] < 2.0, "the point source must spread");
}

#[test]
fn mixing_lbm_set_concentration_validates_the_length() {
    let mut m = MixingLBM::new(4, 3, 5.0, 1.0);
    assert!(matches!(
        m.set_concentration(&[1.0; 5]),
        Err(Error::IncompatibleParameters { .. })
    ));
    let field: Vec<f64> = (0..12).map(f64::from).collect();
    m.set_concentration(&field).expect("length matches");
    assert_eq!(m.c_scalar, field);
    for (c, node) in m.c_scalar.iter().zip(m.f_c.chunks_exact(9)) {
        assert!((node.iter().sum::<f64>() - c).abs() < 1e-13);
    }
}

// ---------------------------------------------------------------------------
// PassiveScalarField
// ---------------------------------------------------------------------------

fn passive_blob(n: usize, tau: f64, centre: (f64, f64)) -> PassiveScalarField {
    let mut ps = PassiveScalarField::new(n, n, tau, 0.0);
    ps.set_scalar_field(&gaussian_field(n, centre.0, centre.1, SIGMA))
        .expect("field has nx * ny values");
    ps
}

#[test]
fn passive_scalar_variance_grows_with_the_lattice_diffusivity() {
    let n = GRID;
    let steps = STEPS;
    for &tau in &[0.8, 1.0] {
        let mut ps = passive_blob(n, tau, (CENTRE, CENTRE));
        assert!((ps.diffusivity() - lattice_diffusivity(tau)).abs() < 1e-15);
        let start = field_moments(&ps.phi, n, n);
        for _ in 0..steps {
            ps.step();
        }
        let end = field_moments(&ps.phi, n, n);
        let expected = exact_variance_growth(tau, steps);
        assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
        assert!(
            (end.variance_x - start.variance_x - expected).abs() < 1e-5,
            "τ={tau}: {} vs {expected}",
            end.variance_x - start.variance_x
        );
        assert!((end.variance_y - start.variance_y - expected).abs() < 1e-5);
    }
}

#[test]
fn passive_scalar_is_advected_by_the_velocity_field() {
    let n = GRID;
    let steps = STEPS;
    let tau = 0.8;
    let mut ps = PassiveScalarField::new(n, n, tau, 0.0);
    ps.set_velocity(vec![0.05; n * n], vec![-0.02; n * n]);
    ps.set_scalar_field(&gaussian_field(n, CENTRE - 4.0, CENTRE + 2.0, SIGMA))
        .expect("field has nx * ny values");
    let start = field_moments(&ps.phi, n, n);
    for _ in 0..steps {
        ps.step();
    }
    let end = field_moments(&ps.phi, n, n);
    assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
    assert!((end.centroid_x - start.centroid_x - 0.05 * steps as f64).abs() < 1e-5);
    assert!((end.centroid_y - start.centroid_y + 0.02 * steps as f64).abs() < 1e-5);
    let expected = exact_variance_growth(tau, steps);
    assert!((end.variance_x - start.variance_x - expected).abs() < 1e-4);
}

/// `collide` alone is local: a point value stays on its node, relaxes toward
/// `w_q φ` and keeps its mass.
#[test]
fn passive_scalar_collide_is_local_and_streaming_is_separate() {
    let mut ps = PassiveScalarField::new(8, 8, 0.7, 0.0);
    let mut field = vec![0.0; 64];
    field[3 * 8 + 4] = 1.0;
    ps.set_scalar_field(&field)
        .expect("field has nx * ny values");
    // Perturb the populations away from equilibrium without changing the mass.
    let node = 9 * (3 * 8 + 4);
    ps.g[node + 1] += 0.05;
    ps.g[node + 3] -= 0.05;
    ps.collide();
    assert!((ps.total_scalar() - 1.0).abs() < 1e-14);
    assert!(
        (ps.phi[3 * 8 + 4] - 1.0).abs() < 1e-14,
        "no transport during collision"
    );
    ps.stream();
    assert!((ps.total_scalar() - 1.0).abs() < 1e-14);
    assert!(
        ps.phi[3 * 8 + 4] < 1.0 - 1e-3,
        "streaming must move the scalar"
    );
}

#[test]
fn passive_scalar_direct_phi_edits_are_kept() {
    let mut ps = PassiveScalarField::new(6, 6, 1.0, 1.0);
    ps.phi[10] = 4.0;
    let expected = 36.0 + 3.0;
    for _ in 0..10 {
        ps.step();
    }
    assert!((ps.total_scalar() - expected).abs() < 1e-11);
}

#[test]
fn passive_scalar_set_scalar_field_validates_lengths() {
    let mut ps = PassiveScalarField::new(4, 4, 1.0, 0.0);
    assert!(matches!(
        ps.set_scalar_field(&[0.0; 3]),
        Err(Error::IncompatibleParameters { .. })
    ));
    ps.set_velocity(vec![0.0; 3], vec![0.0; 3]);
    assert!(ps.set_scalar_field(&[0.0; 16]).is_err());
}

// ---------------------------------------------------------------------------
// DoubleDistributionFunction
// ---------------------------------------------------------------------------

/// A density excess streams to the neighbours: after one step node `(x+1, y)`
/// has gained `w_1 δρ`, the diagonal neighbour `w_5 δρ`, and the centre keeps
/// only its rest population `w_0 δρ`.  The previous `step` never moved anything.
#[test]
fn ddf_density_pulse_propagates_one_link_per_step() {
    let n = 9;
    let delta = 0.01;
    let mut ddf = DoubleDistributionFunction::new(n, n, 0.8, 1.0, 1.0, 0.0);
    ddf.set_flow_at(4, 4, 1.0 + delta, 0.0, 0.0)
        .expect("inside");
    let mass0 = ddf.total_mass();
    ddf.step();
    let rho = |x: usize, y: usize| ddf.rho[y * n + x];
    // τ_f ≠ 1 does not matter: equilibrium populations are not changed by
    // the collision, only moved by the streaming.
    assert!((rho(4, 4) - (1.0 + 4.0 * delta / 9.0)).abs() < 1e-14);
    assert!((rho(5, 4) - (1.0 + delta / 9.0)).abs() < 1e-14);
    assert!((rho(3, 4) - (1.0 + delta / 9.0)).abs() < 1e-14);
    assert!((rho(4, 5) - (1.0 + delta / 9.0)).abs() < 1e-14);
    assert!((rho(5, 5) - (1.0 + delta / 36.0)).abs() < 1e-14);
    assert!((rho(3, 3) - (1.0 + delta / 36.0)).abs() < 1e-14);
    assert!((rho(6, 4) - 1.0).abs() < 1e-14, "only one link per step");
    assert!((ddf.total_mass() - mass0).abs() < 1e-12);
}

/// The fluid carries momentum: a uniform flow keeps its velocity (the previous
/// collision relaxed toward `w ρ` and erased the momentum in one step).
#[test]
fn ddf_conserves_and_keeps_momentum() {
    let (nx, ny) = (10, 6);
    let mut ddf = DoubleDistributionFunction::new(nx, ny, 0.9, 1.2, 1.0, 1.0);
    for y in 0..ny {
        for x in 0..nx {
            ddf.set_flow_at(x, y, 1.0, 0.04, -0.01).expect("inside");
        }
    }
    for _ in 0..50 {
        ddf.step();
    }
    for y in 0..ny {
        for x in 0..nx {
            let (ux, uy) = ddf.velocity_at(x, y).expect("inside");
            assert!((ux - 0.04).abs() < 1e-12 && (uy + 0.01).abs() < 1e-12);
        }
    }
    let jx: f64 = ddf
        .f
        .chunks_exact(9)
        .map(|node| {
            node.iter()
                .enumerate()
                .map(|(q, f)| f * D2Q9_CX[q])
                .sum::<f64>()
        })
        .sum();
    assert!((jx - 0.04 * (nx * ny) as f64).abs() < 1e-10);
}

/// A shear wave of the fluid decays with `ν = c_s² (τ_f − ½)`.
#[test]
fn ddf_shear_wave_decays_with_the_lattice_viscosity() {
    let (nx, ny) = (4, 32);
    let amplitude = 0.01;
    let steps = 100;
    let mut ddf = DoubleDistributionFunction::new(nx, ny, 0.8, 1.0, 1.0, 0.0);
    let phase = |y: usize| (2.0 * std::f64::consts::PI * y as f64 / ny as f64).sin();
    for y in 0..ny {
        for x in 0..nx {
            ddf.set_flow_at(x, y, 1.0, amplitude * phase(y), 0.0)
                .expect("inside");
        }
    }
    for _ in 0..steps {
        ddf.step();
    }
    let mut projection = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            let (ux, _) = ddf.velocity_at(x, y).expect("inside");
            projection += ux * phase(y);
        }
    }
    let measured = 2.0 * projection / (nx * ny) as f64 / amplitude;
    let k = 2.0 * std::f64::consts::PI / ny as f64;
    let expected = (-ddf.viscosity() * k * k * steps as f64).exp();
    assert!(
        (measured / expected - 1.0).abs() < 0.01,
        "measured {measured}, expected {expected}"
    );
}

#[test]
fn ddf_scalar_diffuses_with_its_own_relaxation_time() {
    let n = GRID;
    let steps = STEPS;
    for &tau_g in &[0.8, 1.0] {
        let mut ddf = DoubleDistributionFunction::new(n, n, 1.0, tau_g, 1.0, 0.0);
        let field = gaussian_field(n, CENTRE, CENTRE, SIGMA);
        for j in 0..n {
            for i in 0..n {
                ddf.set_scalar_at(i, j, field[j * n + i]).expect("inside");
            }
        }
        let start = field_moments(&ddf.phi, n, n);
        for _ in 0..steps {
            ddf.step();
        }
        let end = field_moments(&ddf.phi, n, n);
        let expected = exact_variance_growth(tau_g, steps);
        assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
        assert!(
            (end.variance_x - start.variance_x - expected).abs() < 1e-5,
            "τ_g={tau_g}: {} vs {expected}",
            end.variance_x - start.variance_x
        );
        assert!((ddf.total_scalar() - start.mass).abs() < 1e-10);
    }
}

/// The scalar rides on the fluid velocity.
#[test]
fn ddf_scalar_is_advected_by_the_fluid() {
    let n = GRID;
    let steps = STEPS;
    let mut ddf = DoubleDistributionFunction::new(n, n, 0.9, 0.8, 1.0, 0.0);
    for j in 0..n {
        for i in 0..n {
            ddf.set_flow_at(i, j, 1.0, 0.05, 0.02).expect("inside");
        }
    }
    let field = gaussian_field(n, CENTRE - 4.0, CENTRE - 2.0, SIGMA);
    for j in 0..n {
        for i in 0..n {
            ddf.set_scalar_at(i, j, field[j * n + i]).expect("inside");
        }
    }
    let start = field_moments(&ddf.phi, n, n);
    for _ in 0..steps {
        ddf.step();
    }
    let end = field_moments(&ddf.phi, n, n);
    assert!((end.centroid_x - start.centroid_x - 0.05 * steps as f64).abs() < 1e-5);
    assert!((end.centroid_y - start.centroid_y - 0.02 * steps as f64).abs() < 1e-5);
    assert!((end.mass - start.mass).abs() < 1e-12 * start.mass);
}

#[test]
fn ddf_direct_field_edits_are_kept_and_indices_are_checked() {
    let mut ddf = DoubleDistributionFunction::new(5, 5, 1.0, 1.0, 1.0, 1.0);
    ddf.rho[7] = 1.5;
    ddf.phi[8] = 3.0;
    let mass = ddf.total_mass();
    let scalar = ddf.total_scalar();
    assert!((mass - 25.5).abs() < 1e-14 && (scalar - 27.0).abs() < 1e-14);
    for _ in 0..12 {
        ddf.step();
    }
    assert!((ddf.total_mass() - mass).abs() < 1e-11);
    assert!((ddf.total_scalar() - scalar).abs() < 1e-11);
    assert!(matches!(
        ddf.set_flow_at(5, 0, 1.0, 0.0, 0.0),
        Err(Error::IndexOutOfBounds { .. })
    ));
    assert!(ddf.set_scalar_at(0, 5, 1.0).is_err());
    assert_eq!(ddf.velocity_at(5, 5), None);
}

// ---------------------------------------------------------------------------
// MultiComponentLbm
// ---------------------------------------------------------------------------

fn species_blob(n: usize, centre: (f64, f64)) -> Vec<f64> {
    gaussian_field(n, centre.0, centre.1, SIGMA)
}

/// Each species spreads with its own `D_α = c_s² (τ_α − ½)` and none of them
/// leaks mass.
#[test]
fn multi_component_species_diffuse_independently() {
    let n = GRID;
    let steps = STEPS;
    let taus = [0.8, 1.0];
    let mut mc = MultiComponentLbm::new(n, n, taus.to_vec(), vec![0.0, 0.0]);
    mc.set_species_field(0, &species_blob(n, (CENTRE, CENTRE)))
        .expect("valid");
    mc.set_species_field(1, &species_blob(n, (CENTRE - 1.0, CENTRE - 1.0)))
        .expect("valid");
    let starts: Vec<FieldMoments> = (0..2).map(|a| field_moments(&mc.phi[a], n, n)).collect();
    for _ in 0..steps {
        mc.step();
    }
    for (alpha, &tau) in taus.iter().enumerate() {
        let end = field_moments(&mc.phi[alpha], n, n);
        let start = &starts[alpha];
        let expected = exact_variance_growth(tau, steps);
        assert!(
            (end.mass - start.mass).abs() < 1e-12 * start.mass,
            "species {alpha}"
        );
        assert!(
            (end.variance_x - start.variance_x - expected).abs() < 1e-5,
            "species {alpha}: {} vs {expected}",
            end.variance_x - start.variance_x
        );
        assert!(
            (end.variance_y - start.variance_y - expected).abs() < 1e-5,
            "species {alpha}: Δvar_y = {}, expected {expected}",
            end.variance_y - start.variance_y
        );
        assert!((end.centroid_x - start.centroid_x).abs() < 1e-5);
        assert!((end.centroid_y - start.centroid_y).abs() < 1e-5);
        let d = mc.species_diffusivity(alpha).expect("species exists");
        assert!((d - lattice_diffusivity(tau)).abs() < 1e-15);
    }
    assert_eq!(mc.species_diffusivity(2), None);
}

/// A unit mass on one node is shared with the neighbours by the lattice
/// weights after one step.
#[test]
fn multi_component_step_streams_a_point_source() {
    let n = 7;
    let mut mc = MultiComponentLbm::new(n, n, vec![1.0], vec![0.0]);
    let mut field = vec![0.0; n * n];
    field[3 * n + 3] = 1.0;
    mc.set_species_field(0, &field).expect("valid");
    mc.step();
    let at = |x: usize, y: usize| mc.phi[0][y * n + x];
    assert!((at(3, 3) - 4.0 / 9.0).abs() < 1e-15, "w_0 stays");
    assert!((at(4, 3) - 1.0 / 9.0).abs() < 1e-15);
    assert!((at(3, 4) - 1.0 / 9.0).abs() < 1e-15);
    assert!((at(2, 3) - 1.0 / 9.0).abs() < 1e-15);
    assert!((at(3, 2) - 1.0 / 9.0).abs() < 1e-15);
    assert!((at(4, 4) - 1.0 / 36.0).abs() < 1e-15);
    assert!((at(2, 2) - 1.0 / 36.0).abs() < 1e-15);
    assert!((mc.total_concentration(0) - 1.0).abs() < 1e-15);
}

/// With cross-diffusion switched on every species keeps its total.
#[test]
fn multi_component_cross_diffusion_conserves_each_species() {
    let n = 32;
    let mut mc = MultiComponentLbm::new(n, n, vec![0.8, 1.0, 0.9], vec![0.5, 0.2, 0.0]);
    mc.set_species_field(2, &species_blob(n, (16.0, 16.0)))
        .expect("valid");
    mc.set_cross_diff(0, 2, 0.05);
    mc.set_cross_diff(1, 2, -0.03);
    mc.set_cross_diff(2, 0, 0.02);
    let totals: Vec<f64> = (0..3).map(|a| mc.total_concentration(a)).collect();
    for _ in 0..100 {
        mc.step();
    }
    for (alpha, &before) in totals.iter().enumerate() {
        let after = mc.total_concentration(alpha);
        assert!(
            (after - before).abs() < 1e-11 * before.max(1.0),
            "species {alpha}: {before} -> {after}"
        );
    }
}

/// One step of the source on a uniform species 0 and a sinusoidal species 1:
/// `φ₀ − φ̄ = D₀₁ λ_k T_k A sin(kx)` with the lattice Laplacian eigenvalue
/// `λ_k = 2 (cos k − 1)` and the streaming factor `T_k = 2/3 + cos(k)/3`.
#[test]
fn multi_component_cross_diffusion_first_step_matches_the_stencil() {
    let (nx, ny) = (32, 8);
    let k = 2.0 * std::f64::consts::PI / nx as f64;
    let amplitude = 0.1;
    let d01 = 0.05;
    let mut mc = MultiComponentLbm::new(nx, ny, vec![1.0, 1.0], vec![1.0, 1.0]);
    let wave: Vec<f64> = (0..nx * ny)
        .map(|n| 1.0 + amplitude * (k * (n % nx) as f64).sin())
        .collect();
    mc.set_species_field(1, &wave).expect("valid");
    mc.set_cross_diff(0, 1, d01);
    mc.step();
    let lambda = 2.0 * (k.cos() - 1.0);
    let factor = 2.0 / 3.0 + k.cos() / 3.0;
    for (n, &value) in mc.phi[0].iter().enumerate() {
        let expected = 1.0 + d01 * lambda * factor * amplitude * (k * (n % nx) as f64).sin();
        assert!(
            (value - expected).abs() < 1e-13,
            "node {n}: {value} vs {expected}"
        );
    }
    // Species 1 is not driven by species 0 (D₁₀ = 0): its amplitude only
    // shrinks by the streaming factor T_k.
    let proj: f64 = mc.phi[1]
        .iter()
        .enumerate()
        .map(|(n, v)| (v - 1.0) * (k * (n % nx) as f64).sin())
        .sum::<f64>()
        * 2.0
        / (nx * ny) as f64;
    assert!((proj / amplitude - factor).abs() < 1e-12);
}

/// Several steps follow the continuum solution `a₀(t) = −D₀₁ k² A t e^{−D k² t}`
/// (equal self-diffusivities) of the coupled diffusion equations.
#[test]
fn multi_component_cross_diffusion_follows_the_continuum_solution() {
    let (nx, ny) = (32, 8);
    let k = 2.0 * std::f64::consts::PI / nx as f64;
    let amplitude = 0.1;
    let d01 = 0.05;
    let steps = 8;
    let mut mc = MultiComponentLbm::new(nx, ny, vec![1.0, 1.0], vec![1.0, 1.0]);
    let wave: Vec<f64> = (0..nx * ny)
        .map(|n| 1.0 + amplitude * (k * (n % nx) as f64).sin())
        .collect();
    mc.set_species_field(1, &wave).expect("valid");
    mc.set_cross_diff(0, 1, d01);
    for _ in 0..steps {
        mc.step();
    }
    let projection: f64 = mc.phi[0]
        .iter()
        .enumerate()
        .map(|(n, v)| (v - 1.0) * (k * (n % nx) as f64).sin())
        .sum::<f64>()
        * 2.0
        / (nx * ny) as f64;
    let t = steps as f64;
    let d_self = lattice_diffusivity(1.0);
    let expected = -d01 * k * k * amplitude * t * (-d_self * k * k * t).exp();
    assert!(expected < 0.0 && projection < 0.0, "sign: {projection}");
    assert!(
        (projection / expected - 1.0).abs() < 0.02,
        "projection {projection}, continuum {expected}"
    );
}

/// Species α is depleted where species β has a maximum and `D_αβ > 0`; a
/// negative coefficient reverses the sign.
#[test]
fn multi_component_cross_diffusion_sign_follows_the_coefficient() {
    let n = 24;
    let centre = n / 2 * n + n / 2;
    for &(d_ab, sign) in &[(0.06, -1.0), (-0.06, 1.0)] {
        let mut mc = MultiComponentLbm::new(n, n, vec![1.0, 1.0], vec![1.0, 0.0]);
        mc.set_species_field(1, &species_blob(n, (12.0, 12.0)))
            .expect("valid");
        mc.set_cross_diff(0, 1, d_ab);
        for _ in 0..5 {
            mc.step();
        }
        let shift = mc.phi[0][centre] - 1.0;
        assert!(
            shift * sign > 1e-4,
            "D={d_ab}: φ₀ at the centre changed by {shift}"
        );
    }
}

/// Zero cross-diffusion and the ignored diagonal reproduce the uncoupled
/// collide-then-stream sequence bit for bit.
#[test]
fn multi_component_without_coupling_equals_collide_and_stream() {
    let n = 16;
    let build = || {
        let mut mc = MultiComponentLbm::new(n, n, vec![0.7, 1.1], vec![1.0, 0.0]);
        mc.set_species_field(1, &species_blob(n, (8.0, 8.0)))
            .expect("valid");
        mc
    };
    let mut coupled = build();
    coupled.set_cross_diff(0, 0, 0.2);
    coupled.set_cross_diff(1, 1, 0.1);
    let mut manual = build();
    for _ in 0..10 {
        coupled.step();
        manual.collide_uncoupled();
        manual.stream();
    }
    assert_eq!(coupled.g, manual.g);
    assert_eq!(coupled.phi, manual.phi);
}

#[test]
fn multi_component_direct_phi_edits_are_kept() {
    let mut mc = MultiComponentLbm::new(8, 8, vec![0.9], vec![1.0]);
    mc.phi[0][20] += 5.0;
    for _ in 0..10 {
        mc.step();
    }
    assert!((mc.total_concentration(0) - 69.0).abs() < 1e-11);
}

#[test]
fn multi_component_set_species_field_validates_its_arguments() {
    let mut mc = MultiComponentLbm::new(4, 4, vec![1.0, 1.0], vec![0.0, 0.0]);
    assert!(matches!(
        mc.set_species_field(2, &[0.0; 16]),
        Err(Error::InvalidSpeciesCount {
            min_species: 3,
            actual: 2
        })
    ));
    assert!(matches!(
        mc.set_species_field(0, &[0.0; 15]),
        Err(Error::IncompatibleParameters { .. })
    ));
    let field: Vec<f64> = (0..16).map(f64::from).collect();
    mc.set_species_field(1, &field).expect("valid");
    assert_eq!(mc.phi[1], field);
    for (c, node) in field.iter().zip(mc.g[1].chunks_exact(9)) {
        assert!((node.iter().sum::<f64>() - c).abs() < 1e-13);
    }
}
