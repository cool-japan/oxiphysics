// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the edge fluxes of [`FluidFilmSimulation::step`].
//!
//! Gravity has to carry film fluid *downhill*, i.e. from node `i` to node `j`
//! when the surface-tangent gravity component along `x_j - x_i` is positive,
//! and the hydrostatic levelling term has to carry it from the thicker to the
//! thinner node.  The gravity flux used to be added to the wrong node, so the
//! film flowed uphill and gravity and levelling reinforced each other instead of
//! balancing at the equilibrium slope.

use super::*;

const H0: f64 = 1.0e-3;
const SPACING: f64 = 0.01;
const DT: f64 = 1.0e-3;

/// A viscous oil film: slow enough that the explicit scheme is far from its
/// diffusion limit (`h_min` never matters, volume stays conserved).
fn oil() -> FluidFilmProperties {
    FluidFilmProperties::new(1000.0, 1.0, 0.03, 293.15, 0.0)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn thickness(sim: &FluidFilmSimulation) -> Vec<f64> {
    sim.nodes.iter().map(|n| n.thickness).collect()
}

/// Index of the node that lies furthest along the tangential gravity of the
/// first node (the downhill end) and of the one furthest against it (uphill).
fn downhill_and_uphill(sim: &FluidFilmSimulation) -> (usize, usize) {
    let g_tangent = sim.nodes[0].tangent_gravity(sim.gravity);
    let depth = |k: usize| dot(sim.nodes[k].position, g_tangent);
    let mut downhill = 0;
    let mut uphill = 0;
    for k in 1..sim.nodes.len() {
        if depth(k) > depth(downhill) {
            downhill = k;
        }
        if depth(k) < depth(uphill) {
            uphill = k;
        }
    }
    (downhill, uphill)
}

fn assert_film_dynamics_are_sane(sim: &FluidFilmSimulation, volume_before: f64) {
    let volume_after = sim.total_volume();
    assert!(
        (volume_after - volume_before).abs() < 1.0e-12 * volume_before,
        "the film volume must be conserved: {volume_before} -> {volume_after}"
    );
    let thinnest = thickness(sim).into_iter().fold(f64::INFINITY, f64::min);
    assert!(
        thinnest > 0.5 * H0 && thinnest > 1.0e6 * sim.h_min,
        "the numerical floor h_min must never be active, thinnest film {thinnest}"
    );
}

#[test]
fn gravity_drains_the_uphill_node_of_an_inclined_plane() {
    // A plane inclined by 30 degrees: normal (sin, cos, 0), steepest descent
    // direction (cos, -sin, 0).  Nodes are laid out along the descent direction.
    let theta = 30.0_f64.to_radians();
    let normal = [theta.sin(), theta.cos(), 0.0];
    let descent = [theta.cos(), -theta.sin(), 0.0];
    let nodes: Vec<FilmNode> = (0..3)
        .map(|k| {
            let s = k as f64 * SPACING;
            FilmNode::new([s * descent[0], s * descent[1], s * descent[2]], normal, H0)
        })
        .collect();
    let gravity = [0.0, -G_ACCEL, 0.0];
    let mut sim = FluidFilmSimulation::new(nodes, vec![(0, 1), (1, 2)], oil(), gravity);
    let (downhill, uphill) = downhill_and_uphill(&sim);
    assert_eq!(
        (downhill, uphill),
        (2, 0),
        "test set-up: node 2 is the lowest"
    );
    let volume_before = sim.total_volume();

    let steps = 100;
    for _ in 0..steps {
        sim.step(DT);
    }

    let h = thickness(&sim);
    assert!(h[uphill] < H0, "the uphill node loses fluid: {h:?}");
    assert!(h[downhill] > H0, "the downhill node gains fluid: {h:?}");
    assert!(
        h[uphill] < h[1] && h[1] < h[downhill],
        "the film must thicken downhill: {h:?}"
    );
    // Magnitude: q = rho g sin(theta) h^3 / (3 mu) leaves the uphill node.
    let flux = 1000.0 * G_ACCEL * theta.sin() * H0.powi(3) / 3.0;
    let expected_loss = steps as f64 * flux * DT / SPACING;
    let loss = H0 - h[uphill];
    assert!(
        (loss - expected_loss).abs() < 0.03 * expected_loss,
        "uphill thickness loss {loss} vs analytic {expected_loss}"
    );
    assert_film_dynamics_are_sane(&sim, volume_before);
}

#[test]
fn tilted_normals_along_x_move_fluid_towards_the_gravity_side() {
    // Nodes along x with the surface normals tilted about z (both tilt signs):
    // the downhill end is whichever node the tangential gravity points to.
    for tilt_deg in [30.0_f64, -30.0] {
        let tilt = tilt_deg.to_radians();
        let normal = [tilt.sin(), tilt.cos(), 0.0];
        let nodes: Vec<FilmNode> = (0..3)
            .map(|k| FilmNode::new([k as f64 * SPACING, 0.0, 0.0], normal, H0))
            .collect();
        let mut sim =
            FluidFilmSimulation::new(nodes, vec![(0, 1), (1, 2)], oil(), [0.0, -9.81, 0.0]);
        let (downhill, uphill) = downhill_and_uphill(&sim);
        assert_ne!(downhill, uphill);
        let volume_before = sim.total_volume();

        for _ in 0..100 {
            sim.step(DT);
        }

        let h = thickness(&sim);
        let thickest = (0..3)
            .max_by(|&a, &b| h[a].partial_cmp(&h[b]).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(0);
        let thinnest = (0..3)
            .min_by(|&a, &b| h[a].partial_cmp(&h[b]).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(0);
        assert_eq!(thickest, downhill, "tilt {tilt_deg}: thicknesses {h:?}");
        assert_eq!(thinnest, uphill, "tilt {tilt_deg}: thicknesses {h:?}");
        assert_film_dynamics_are_sane(&sim, volume_before);
    }
}

#[test]
fn single_edge_gravity_flux_has_the_analytic_direction_and_magnitude() {
    let props = oil();
    let h: f64 = 2.0e-4;
    let flux_scale = props.density * h.powi(3) / (3.0 * props.viscosity);
    // Surface tangent gravity along +x and along -x.
    for g_x in [2.0_f64, -2.0] {
        let nodes = vec![
            FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], h),
            FilmNode::new([SPACING, 0.0, 0.0], [0.0, 1.0, 0.0], h),
        ];
        let mut sim = FluidFilmSimulation::new(nodes, vec![(0, 1)], props.clone(), [g_x, 0.0, 0.0]);
        sim.step(DT);

        // Fluid moves in the direction gravity points: node 1 gains for g_x > 0.
        let expected_gain = flux_scale * g_x * DT / SPACING;
        let h_after = thickness(&sim);
        assert!(
            (h_after[1] - h - expected_gain).abs() < 1.0e-9 * expected_gain.abs(),
            "g_x = {g_x}: node 1 changed by {} but {expected_gain} was expected",
            h_after[1] - h
        );
        assert!(
            (h_after[0] - h + expected_gain).abs() < 1.0e-9 * expected_gain.abs(),
            "g_x = {g_x}: node 0 changed by {} but {} was expected",
            h_after[0] - h,
            -expected_gain
        );
        assert_eq!(
            h_after[1] > h,
            g_x > 0.0,
            "the node gravity points to must gain fluid"
        );
    }
}

#[test]
fn levelling_balances_gravity_at_the_equilibrium_slope() {
    // Gravity along +x pushes fluid towards node 1; the film surface slope
    // dh/dx = g_along / G makes the hydrostatic levelling flux cancel it exactly.
    // With both fluxes added (the old sign error) the edge would not balance.
    let props = oil();
    let g_along = 0.01;
    let h_i = H0;
    let h_j = H0 + SPACING * g_along / G_ACCEL;
    let nodes = vec![
        FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], h_i),
        FilmNode::new([SPACING, 0.0, 0.0], [0.0, 1.0, 0.0], h_j),
    ];
    let mut sim = FluidFilmSimulation::new(nodes, vec![(0, 1)], props.clone(), [g_along, 0.0, 0.0]);
    sim.step(DT);

    let h_avg = 0.5 * (h_i + h_j);
    let single_flux =
        props.density * g_along * h_avg.powi(3) / (3.0 * props.viscosity) * DT / SPACING;
    let h_after = thickness(&sim);
    assert!(single_flux > 0.0);
    assert!(
        (h_after[0] - h_i).abs() < 1.0e-6 * single_flux,
        "node 0 must stay in equilibrium, changed by {}",
        h_after[0] - h_i
    );
    assert!(
        (h_after[1] - h_j).abs() < 1.0e-6 * single_flux,
        "node 1 must stay in equilibrium, changed by {}",
        h_after[1] - h_j
    );
}

#[test]
fn levelling_moves_fluid_from_the_thick_to_the_thin_node() {
    let props = oil();
    let (h_thick, h_thin) = (3.0e-4, 1.0e-4);
    let nodes = vec![
        FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], h_thick),
        FilmNode::new([SPACING, 0.0, 0.0], [0.0, 1.0, 0.0], h_thin),
    ];
    let mut sim = FluidFilmSimulation::new(nodes, vec![(0, 1)], props.clone(), [0.0; 3]);
    let volume_before = sim.total_volume();
    sim.step(DT);

    let h_avg = 0.5 * (h_thick + h_thin);
    let dh_dx = (h_thin - h_thick) / SPACING;
    // Flux from i to j is -rho G h^3 / (3 mu) dh/dx (positive: thick -> thin).
    let expected_loss =
        -props.density * G_ACCEL * h_avg.powi(3) * dh_dx / (3.0 * props.viscosity) * DT / SPACING;
    let h = thickness(&sim);
    assert!(expected_loss > 0.0);
    assert!((h_thick - h[0] - expected_loss).abs() < 1.0e-9 * expected_loss);
    assert!((h[1] - h_thin - expected_loss).abs() < 1.0e-9 * expected_loss);
    assert!((sim.total_volume() - volume_before).abs() < 1.0e-18);
}
