// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the curvature / twist strain measures and the nodal
//! moments of [`InternalForces`].
//!
//! The Darboux vector `Ω = κ₁ d₁ + κ₂ d₂ + τ d₃` satisfies `d₃' = κ₂ d₁ − κ₁ d₂`:
//! a bend about **d**₁ tilts **d**₃ along **−d**₂.  The earlier implementation
//! projected `d₃' ` onto `d₁` for `κ₁` and onto `d₂` for `κ₂`, so the moment that
//! resists a bend about one axis was applied about the *other* axis, and the
//! nodal moment was divided by the segment length.  Every earlier test used
//! segments of length 1, where that division is invisible.

use super::*;

const X: Vec3 = [1.0, 0.0, 0.0];
const Y: Vec3 = [0.0, 1.0, 0.0];
const Z: Vec3 = [0.0, 0.0, 1.0];

/// Stiffness with only bending and twist: no axial and no shear coupling, so that
/// the nodal moments are exactly the joint moments.
fn bending_only(ei1: f64, ei2: f64, gj: f64) -> RodStiffness {
    RodStiffness {
        ei1,
        ei2,
        gj,
        ea: 0.0,
        ga: 0.0,
    }
}

/// Two-node rod of the given length along z.
fn two_node_rod(length: f64, stiffness: RodStiffness) -> CosseratRod {
    CosseratRod::new_straight(2, length, stiffness, 1.0, 1.0e-3)
}

/// Rotate the directors of node `node` by `angle` about the world axis `axis`.
fn rotate_node(rod: &mut CosseratRod, node: usize, axis: Vec3, angle: f64) {
    rod.nodes[node].rotate_directors(axis, angle);
}

fn assert_vec_close(actual: Vec3, expected: Vec3, tol: f64, what: &str) {
    for d in 0..3 {
        assert!(
            (actual[d] - expected[d]).abs() <= tol,
            "{what}: component {d}: {actual:?} vs {expected:?} (tol {tol:e})"
        );
    }
}

// A rotation of node b about d1 is a curvature about d1: κ₁ = sin φ / ℓ, κ₂ = τ = 0.
#[test]
fn rotation_about_d1_is_curvature_about_d1() {
    let length = 0.5;
    let phi = 0.01_f64;
    let mut rod = two_node_rod(length, bending_only(3.0, 5.0, 7.0));
    rotate_node(&mut rod, 1, X, phi);

    let kappa = InternalForces::darboux_vector(&rod, 0);
    assert_vec_close(
        kappa,
        [phi.sin() / length, 0.0, 0.0],
        1e-12,
        "Darboux vector",
    );

    let (_eps, _g1, _g2, dk1, dk2, dtau) = InternalForces::segment_strains(&rod, 0);
    assert!((dk1 - phi.sin() / length).abs() < 1e-12, "Δκ₁ = {dk1}");
    assert!(dk2.abs() < 1e-12, "Δκ₂ = {dk2}");
    assert!(dtau.abs() < 1e-12, "Δτ = {dtau}");
}

// The moment that resists the bend about d1 acts about d1: on node b it is
// anti-parallel to d1, on node a parallel, and its size is EI₁ κ₁ (not EI₁ κ₁ / ℓ).
#[test]
fn bend_about_d1_is_resisted_by_a_moment_about_d1() {
    let length = 0.5;
    let phi = 0.01_f64;
    let ei1 = 3.0;
    let mut rod = two_node_rod(length, bending_only(ei1, 5.0, 7.0));
    rotate_node(&mut rod, 1, X, phi);

    let (f_a, m_a, f_b, m_b) = InternalForces::nodal_forces(&rod, 0);
    let moment = ei1 * phi.sin() / length;
    assert_vec_close(m_b, [-moment, 0.0, 0.0], 1e-12 * moment, "moment on b");
    assert_vec_close(m_a, [moment, 0.0, 0.0], 1e-12 * moment, "moment on a");
    assert!(m_b[0] < 0.0, "moment on b is anti-parallel to d1: {m_b:?}");
    assert_vec_close(f_a, [0.0; 3], 1e-12, "no force on a");
    assert_vec_close(f_b, [0.0; 3], 1e-12, "no force on b");
}

// A rotation about d2 is a curvature about d2 and resisted about d2 with EI₂.
#[test]
fn bend_about_d2_is_resisted_by_a_moment_about_d2() {
    let length = 0.8;
    let phi = 0.02_f64;
    let ei2 = 5.0;
    let mut rod = two_node_rod(length, bending_only(3.0, ei2, 7.0));
    // A rotation about -d2 is a negative curvature about d2.
    rotate_node(&mut rod, 1, vec3_scale(Y, -1.0), phi);

    let kappa = InternalForces::darboux_vector(&rod, 0);
    assert_vec_close(
        kappa,
        [0.0, -phi.sin() / length, 0.0],
        1e-12,
        "Darboux vector",
    );

    let (_, m_a, _, m_b) = InternalForces::nodal_forces(&rod, 0);
    let moment = ei2 * phi.sin() / length;
    assert_vec_close(m_b, [0.0, moment, 0.0], 1e-12 * moment, "moment on b");
    assert_vec_close(m_a, [0.0, -moment, 0.0], 1e-12 * moment, "moment on a");
}

// A rotation about d3 is twist, resisted about d3 with GJ.
#[test]
fn rotation_about_d3_is_twist() {
    let length = 0.25;
    let phi = 0.03_f64;
    let gj = 7.0;
    let mut rod = two_node_rod(length, bending_only(3.0, 5.0, gj));
    rotate_node(&mut rod, 1, Z, phi);

    let kappa = InternalForces::darboux_vector(&rod, 0);
    assert_vec_close(
        kappa,
        [0.0, 0.0, phi.sin() / length],
        1e-12,
        "Darboux vector",
    );

    let (_, m_a, _, m_b) = InternalForces::nodal_forces(&rod, 0);
    let torque = gj * phi.sin() / length;
    assert_vec_close(m_b, [0.0, 0.0, -torque], 1e-12 * torque, "torque on b");
    assert_vec_close(m_a, [0.0, 0.0, torque], 1e-12 * torque, "torque on a");
}

// The sign convention follows the director derivatives d₃' = κ₂ d₁ − κ₁ d₂,
// d₁' = τ d₂ − κ₂ d₃: build a segment from a prescribed Darboux vector and read it
// back.
#[test]
fn darboux_vector_matches_the_director_derivative_convention() {
    let length = 1.0e-3;
    let omega = [0.7, -1.3, 0.4];
    let mut rod = two_node_rod(length, bending_only(1.0, 1.0, 1.0));
    // Node b = node a rotated by Ω ℓ expressed in the (identity) frame of a.
    let angle = vec3_norm(omega) * length;
    let axis = vec3_normalize(omega);
    rotate_node(&mut rod, 1, axis, angle);

    let kappa = InternalForces::darboux_vector(&rod, 0);
    assert_vec_close(kappa, omega, 1e-6, "recovered Darboux vector");

    // d₃' = κ₂ d₁ − κ₁ d₂ and d₁' = τ d₂ − κ₂ d₃ to first order.
    let a = &rod.nodes[0];
    let b = &rod.nodes[1];
    let dd3 = vec3_scale(vec3_sub(b.d3, a.d3), 1.0 / length);
    let expected_dd3 = vec3_sub(vec3_scale(a.d1, omega[1]), vec3_scale(a.d2, omega[0]));
    assert_vec_close(dd3, expected_dd3, 1e-2, "d₃'");
    let dd1 = vec3_scale(vec3_sub(b.d1, a.d1), 1.0 / length);
    let expected_dd1 = vec3_sub(vec3_scale(a.d2, omega[2]), vec3_scale(a.d3, omega[1]));
    assert_vec_close(dd1, expected_dd1, 1e-2, "d₁'");
}

// Nodes on a circular arc whose directors turn about d1 with the arc have no
// shear and no twist, and exactly the rest curvature about d1 in the Δκ₁ slot.
#[test]
fn circular_arc_about_d1_has_only_curvature_about_d1() {
    let kappa = 0.8_f64;
    let n_nodes = 6;
    let ds = 0.05_f64;
    let stiffness = bending_only(2.0, 3.0, 4.0);
    let mut rod =
        CosseratRod::new_straight(n_nodes, ds * (n_nodes - 1) as f64, stiffness, 1.0, 1e-3);
    for (i, node) in rod.nodes.iter_mut().enumerate() {
        let s = i as f64 * ds;
        let phi = kappa * s;
        node.position = [0.0, (phi.cos() - 1.0) / kappa, phi.sin() / kappa];
        node.d1 = X;
        node.d2 = [0.0, phi.cos(), phi.sin()];
        node.d3 = [0.0, -phi.sin(), phi.cos()];
    }
    let discrete_kappa = (kappa * ds).sin() / ds;
    for seg in 0..rod.n_segments() {
        rod.kappa1_0[seg] = discrete_kappa;
        let (eps, g1, g2, dk1, dk2, dtau) = InternalForces::segment_strains(&rod, seg);
        // The chord of an arc is shorter than the arc by (κ ds)² / 24.
        assert!(
            eps.abs() < (kappa * ds).powi(2) / 20.0,
            "axial strain {eps}"
        );
        assert!(g1.abs() < 1e-12 && g2.abs() < 1e-12, "shear {g1} {g2}");
        assert!(dk1.abs() < 1e-12, "seg {seg}: Δκ₁ = {dk1}");
        assert!(dk2.abs() < 1e-12, "seg {seg}: Δκ₂ = {dk2}");
        assert!(dtau.abs() < 1e-12, "seg {seg}: Δτ = {dtau}");
        let loads = InternalForces::compute(&rod, seg);
        assert!(loads.moment1.abs() < 1e-10 && loads.moment2.abs() < 1e-10);
    }
    // Without the rest curvature the same shape is bent about d1 only.
    rod.kappa1_0 = vec![0.0; rod.n_segments()];
    for seg in 0..rod.n_segments() {
        let (_, _, _, dk1, dk2, dtau) = InternalForces::segment_strains(&rod, seg);
        assert!(
            (dk1 - discrete_kappa).abs() < 1e-12,
            "seg {seg}: Δκ₁ = {dk1}"
        );
        assert!(dk2.abs() < 1e-12 && dtau.abs() < 1e-12);
    }
}

/// Deterministic pseudo-random numbers in `[-1, 1]` (xorshift).
fn noise(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 52) as f64 * 2.0 - 1.0
}

/// A five-node rod with arbitrarily displaced nodes and arbitrarily rotated frames.
fn scrambled_rod(seed: u64) -> CosseratRod {
    let cs = RodCrossSection::new(CrossSectionShape::Rectangular {
        width: 0.03,
        height: 0.02,
    });
    let stiffness = RodStiffness::from_material(2.0e5, 7.0e4, &cs);
    let mut rod = CosseratRod::new_straight(5, 0.8, stiffness, 1.0, 1.0e-3);
    let mut state = seed;
    for node in rod.nodes.iter_mut() {
        for p in node.position.iter_mut() {
            *p += 0.04 * noise(&mut state);
        }
        let axis = vec3_normalize([noise(&mut state), noise(&mut state), noise(&mut state)]);
        node.rotate_directors(axis, 0.3 * noise(&mut state).abs());
        node.orthonormalize();
    }
    rod.kappa1_0 = vec![0.2; 4];
    rod.kappa2_0 = vec![-0.1; 4];
    rod.tau_0 = vec![0.05; 4];
    rod
}

// The elastic loads of a free rod carry no net force and no net torque about the
// origin, whatever the configuration: this is the algebraic form of linear and
// angular momentum conservation.
#[test]
fn internal_loads_carry_no_net_force_or_torque() {
    for seed in [
        0x9E37_79B9_7F4A_7C15_u64,
        0xD1B5_4A32_D192_ED03,
        0x2545_F491_4F6C_DD1D,
    ] {
        let rod = scrambled_rod(seed);
        let mut force = [0.0_f64; 3];
        let mut torque = [0.0_f64; 3];
        let mut scale = 0.0_f64;
        for seg in 0..rod.n_segments() {
            let (f_a, m_a, f_b, m_b) = InternalForces::nodal_forces(&rod, seg);
            let x_a = rod.nodes[seg].position;
            let x_b = rod.nodes[seg + 1].position;
            for (f, m, x) in [(f_a, m_a, x_a), (f_b, m_b, x_b)] {
                let lever = vec3_cross(x, f);
                for d in 0..3 {
                    force[d] += f[d];
                    torque[d] += lever[d] + m[d];
                }
                scale = scale.max(vec3_norm(f)).max(vec3_norm(m));
            }
        }
        assert!(scale > 1.0, "the scrambled rod must be loaded: {scale}");
        for d in 0..3 {
            assert!(
                force[d].abs() < 1e-9 * scale,
                "net force {force:?} (scale {scale})"
            );
            assert!(
                torque[d].abs() < 1e-9 * scale,
                "net torque {torque:?} (scale {scale})"
            );
        }
    }
}

/// Cantilever along z with a rectangular section: width 0.04 along d1, height 0.02
/// along d2, so that EI₂ = E h b³ / 12 is four times EI₁ = E b h³ / 12.
fn rectangular_cantilever(
    n_nodes: usize,
    length: f64,
    young: f64,
    poisson: f64,
) -> (CosseratRod, f64, f64, f64) {
    let (b, h) = (0.04, 0.02);
    let cs = RodCrossSection::new(CrossSectionShape::Rectangular {
        width: b,
        height: h,
    });
    let shear_modulus = young / (2.0 * (1.0 + poisson));
    let stiffness = RodStiffness::from_material(young, shear_modulus, &cs);
    let ga = stiffness.ga;
    let rod = CosseratRod::new_straight(n_nodes, length, stiffness, 0.5, 2.0e-3);
    (rod, cs.i1 * young, cs.i2 * young, ga)
}

/// Relax the cantilever under the current external loads with damped explicit
/// dynamics and return the tip displacement.
fn relax_tip_displacement(rod: &mut CosseratRod, dt: f64, steps: usize, damping_rate: f64) -> Vec3 {
    let bcs = BoundaryConditions::cantilever([0.0; 3], X, Y, Z);
    let tip = rod.n_nodes() - 1;
    let factor = 1.0 - damping_rate * dt;
    for _ in 0..steps {
        CosseratDynamics::step_euler(rod, dt);
        bcs.apply(rod);
        for node in rod.nodes.iter_mut() {
            node.velocity = vec3_scale(node.velocity, factor);
            node.angular_velocity = vec3_scale(node.angular_velocity, factor);
        }
    }
    let rest_tip = [0.0, 0.0, rod.total_rest_length()];
    vec3_sub(rod.nodes[tip].position, rest_tip)
}

// A tip load along d1 bends the cantilever about d2 (stiffness EI₂), a tip load
// along d2 bends it about d1 (stiffness EI₁): the deflections follow the
// Timoshenko beam formula δ = P L³ / (3 EI) + P L / (G A_s) with the matching
// bending stiffness, so the two differ by the ratio EI₂ / EI₁ = 4.
#[test]
fn cantilever_tip_deflection_uses_the_bending_stiffness_of_the_load_axis() {
    let length = 1.0;
    let load = 2.0e-2;
    let dt = 5.0e-4;
    let mut deflections = Vec::new();
    for (axis, direction) in [(0usize, X), (1usize, Y)] {
        let (mut rod, ei1, ei2, ga) = rectangular_cantilever(11, length, 1.0e7, 0.3);
        let tip = rod.n_nodes() - 1;
        rod.nodes[tip].external_force = vec3_scale(direction, load);
        let ei = if axis == 0 { ei2 } else { ei1 };
        let expected = load * length.powi(3) / (3.0 * ei) + load * length / ga;
        // Damping close to critical for the fundamental mode (about 2.6 rad/s for
        // the softer axis).
        let tip_displacement = relax_tip_displacement(&mut rod, dt, 20_000, 5.0);
        let delta = tip_displacement[axis];
        assert!(
            (delta - expected).abs() < 0.01 * expected,
            "load along axis {axis}: tip deflection {delta:e}, Timoshenko {expected:e}"
        );
        let other = tip_displacement[1 - axis];
        assert!(other.abs() < 1e-3 * delta, "cross deflection {other:e}");
        deflections.push(delta);
    }
    let ratio = deflections[1] / deflections[0];
    assert!(
        (ratio - 4.0).abs() < 0.05,
        "δ(d2 load) / δ(d1 load) = {ratio}, expected EI₂/EI₁ = 4"
    );
}
