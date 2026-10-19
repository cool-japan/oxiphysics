// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the curved boundary kernels, lattice descriptor and set.

use super::set::symmetric_pinv_solve;
use super::*;
use crate::lattice::CS2;

const D2Q9_NE: usize = 5;
const D2Q9_SW: usize = 7;

/// Link populations of a uniform equilibrium state `f_eq(ρ, u)`.
fn equilibrium_link(lattice: &CurvedLattice, dir: usize, rho: f64, u: [f64; 3]) -> LinkPopulations {
    let opp = lattice.opposite(dir).expect("valid dir");
    let fi = lattice.equilibrium(dir, rho, u);
    let fo = lattice.equilibrium(opp, rho, u);
    LinkPopulations {
        f_i_xf: fi,
        f_i_xff: fi,
        f_i_xfff: fi,
        f_opp_xf: fo,
        f_opp_xff: fo,
    }
}

fn d2q9_link(dir: usize, q: f64, u_wall: [f64; 3], rho: f64) -> WallLink {
    let lattice = CurvedLattice::d2q9();
    WallLink::new(
        q,
        lattice.weight(dir).expect("weight"),
        lattice.velocity_f64(dir).expect("velocity"),
        u_wall,
        rho,
    )
    .expect("valid link")
}

// --- CurvedLattice ---------------------------------------------------------

#[test]
fn test_d2q9_opposites_are_lattice_links() {
    let l = CurvedLattice::d2q9();
    assert_eq!(l.opposite(1), Some(3), "opposite of E must be W");
    assert_eq!(l.opposite(2), Some(4), "opposite of N must be S");
    assert_eq!(
        l.opposite(D2Q9_NE),
        Some(D2Q9_SW),
        "opposite of NE must be SW"
    );
    for i in 0..l.num_directions() {
        let c = l.velocity_f64(i).expect("c");
        let co = l.velocity_f64(l.opposite(i).expect("opp")).expect("c_opp");
        for d in 0..3 {
            assert_eq!(c[d], -co[d]);
        }
    }
}

#[test]
fn test_lattice_new_rejects_inconsistent_opposites() {
    let l = CurvedLattice::d2q9();
    let mut opp = l.opposites().to_vec();
    opp.swap(1, 2);
    let res = CurvedLattice::new(l.velocities().to_vec(), l.weights().to_vec(), opp, 2);
    assert!(res.is_err());
    let ok = CurvedLattice::new(
        l.velocities().to_vec(),
        l.weights().to_vec(),
        l.opposites().to_vec(),
        2,
    );
    assert_eq!(ok, Ok(l));
}

#[test]
fn test_lattice_equilibrium_moments() {
    for l in [
        CurvedLattice::d2q9(),
        CurvedLattice::d3q19(),
        CurvedLattice::d3q27(),
    ] {
        let u = if l.dimensions() == 2 {
            [0.04, -0.02, 0.0]
        } else {
            [0.04, -0.02, 0.03]
        };
        let f: Vec<f64> = (0..l.num_directions())
            .map(|i| l.equilibrium(i, 1.3, u))
            .collect();
        let (rho, j) = l.moments(&f).expect("moments");
        assert!((rho - 1.3).abs() < 1e-13);
        for d in 0..3 {
            assert!((j[d] - 1.3 * u[d]).abs() < 1e-13);
        }
    }
}

// --- Moving-wall term -------------------------------------------------------

#[test]
fn test_moving_wall_term_value() {
    // 2 w ρ (c·u)/cs² = 2 · (1/36) · 1.2 · 0.15 / (1/3) = 0.03
    let v = bouzidi_wall_correction(1.2, 1.0 / 36.0, [1.0, 1.0, 0.0], [0.1, 0.05, 0.0]);
    assert!((v - 0.03).abs() < 1e-15, "v={v}");
    let link = d2q9_link(D2Q9_NE, 0.5, [0.1, 0.05, 0.0], 1.2);
    assert!((link.moving_wall_term() - 0.03).abs() < 1e-15);
    // Perpendicular wall motion contributes nothing.
    let w = bouzidi_wall_correction(1.0, 1.0 / 9.0, [1.0, 0.0, 0.0], [0.0, 0.3, 0.0]);
    assert_eq!(w, 0.0);
}

#[test]
fn test_wall_link_rejects_invalid_q() {
    assert!(WallLink::new(0.0, 1.0 / 9.0, [1.0, 0.0, 0.0], [0.0; 3], 1.0).is_err());
    assert!(WallLink::new(1.5, 1.0 / 9.0, [1.0, 0.0, 0.0], [0.0; 3], 1.0).is_err());
    assert!(WallLink::new(f64::NAN, 1.0 / 9.0, [1.0, 0.0, 0.0], [0.0; 3], 1.0).is_err());
}

// --- Bouzidi -----------------------------------------------------------------

#[test]
fn test_bouzidi_q_half_is_half_way_bounce_back() {
    let pops = LinkPopulations {
        f_i_xf: 0.31,
        f_i_xff: 0.17,
        f_i_xfff: 0.05,
        f_opp_xf: 0.23,
        f_opp_xff: 0.11,
    };
    let u_w = [0.1, 0.05, 0.0];
    let link = d2q9_link(D2Q9_NE, 0.5, u_w, 1.2);
    let expected = 0.31 - 0.03;
    assert!((bouzidi_bc(&pops, &link) - expected).abs() < 1e-14);
    assert!((bouzidi_quadratic_bc(&pops, &link) - expected).abs() < 1e-14);
    let fh = filippova_hanel(0.31, &link, [0.02, 0.0, 0.0], None, 1.2).expect("fh");
    assert!((fh - expected).abs() < 1e-14, "FH at q=1/2 is half-way BB");
}

#[test]
fn test_bouzidi_linear_branch_values() {
    let pops = LinkPopulations {
        f_i_xf: 0.3,
        f_i_xff: 0.1,
        f_i_xfff: 0.0,
        f_opp_xf: 0.2,
        f_opp_xff: 0.0,
    };
    // q < 1/2: 2q f_i(x_f) + (1−2q) f_i(x_ff)
    let link = d2q9_link(1, 0.25, [0.0; 3], 1.0);
    assert!((bouzidi_bc(&pops, &link) - (0.5 * 0.3 + 0.5 * 0.1)).abs() < 1e-15);
    // q ≥ 1/2: f_i/(2q) + (2q−1)/(2q) f_opp
    let link = d2q9_link(1, 0.8, [0.0; 3], 1.0);
    let expected = 0.3 / 1.6 + 0.6 / 1.6 * 0.2;
    assert!((bouzidi_bc(&pops, &link) - expected).abs() < 1e-15);
}

#[test]
fn test_bouzidi_moving_wall_term_scaled_by_inverse_2q() {
    let zero = LinkPopulations::default();
    let u_w = [0.1, 0.05, 0.0];
    let near = d2q9_link(D2Q9_NE, 0.3, u_w, 1.2);
    assert!((bouzidi_bc(&zero, &near) + 0.03).abs() < 1e-15);
    let far = d2q9_link(D2Q9_NE, 0.75, u_w, 1.2);
    assert!((bouzidi_bc(&zero, &far) + 0.03 / 1.5).abs() < 1e-15);
}

#[test]
fn test_all_schemes_exact_for_uniform_flow_with_comoving_wall() {
    let lattice = CurvedLattice::d2q9();
    let rho = 1.07;
    let u = [0.06, -0.035, 0.0];
    for dir in 1..9 {
        let opp = lattice.opposite(dir).expect("opp");
        let target = lattice.equilibrium(opp, rho, u);
        let pops = equilibrium_link(&lattice, dir, rho, u);
        for &q in &[0.05, 0.2, 0.45, 0.5, 0.55, 0.8, 1.0] {
            let link = d2q9_link(dir, q, u, rho);
            let tol = 1e-14;
            assert!(
                (bouzidi_bc(&pops, &link) - target).abs() < tol,
                "BFL dir={dir} q={q}"
            );
            assert!(
                (bouzidi_quadratic_bc(&pops, &link) - target).abs() < tol,
                "QBFL dir={dir} q={q}"
            );
            assert!(
                (interpolated_bounce_back(&pops, &link) - target).abs() < tol,
                "Yu dir={dir} q={q}"
            );
            let mls = filippova_hanel(pops.f_i_xf, &link, u, Some(u), 1.0 / 0.8).expect("mls");
            assert!((mls - target).abs() < tol, "MLS dir={dir} q={q}");
            let fh = filippova_hanel(pops.f_i_xf, &link, u, None, 1.0 / 0.8).expect("fh");
            assert!((fh - target).abs() < tol, "FH dir={dir} q={q}");
        }
    }
}

#[test]
fn test_quadratic_coefficients_sum_to_one() {
    let ones = LinkPopulations {
        f_i_xf: 1.0,
        f_i_xff: 1.0,
        f_i_xfff: 1.0,
        f_opp_xf: 1.0,
        f_opp_xff: 1.0,
    };
    for &q in &[0.1, 0.3, 0.49, 0.5, 0.6, 0.9, 1.0] {
        let link = d2q9_link(2, q, [0.0; 3], 1.0);
        assert!(
            (bouzidi_quadratic_bc(&ones, &link) - 1.0).abs() < 1e-14,
            "q={q}"
        );
        assert!(
            (interpolated_bounce_back(&ones, &link) - 1.0).abs() < 1e-14,
            "q={q}"
        );
        assert!((bouzidi_bc(&ones, &link) - 1.0).abs() < 1e-14, "q={q}");
    }
}

#[test]
fn test_quadratic_exact_for_quadratic_profile_q_small() {
    // Populations varying quadratically along the link line: g(s) = 1 + s + s².
    // For q < 1/2 the scheme samples g at s = 0, −1, −2 and evaluates at 2q − 1.
    let g = |s: f64| 1.0 + s + s * s;
    let pops = LinkPopulations {
        f_i_xf: g(0.0),
        f_i_xff: g(-1.0),
        f_i_xfff: g(-2.0),
        f_opp_xf: 0.0,
        f_opp_xff: 0.0,
    };
    let q = 0.3;
    let link = d2q9_link(2, q, [0.0; 3], 1.0);
    assert!((bouzidi_quadratic_bc(&pops, &link) - g(2.0 * q - 1.0)).abs() < 1e-14);
}

// --- Filippova–Hänel ---------------------------------------------------------

#[test]
fn test_mls_chi_values_and_singularities() {
    let omega = 1.0 / 0.8;
    // q ≥ 1/2: χ = (2q − 1)/(τ + 1/2)
    let chi = mei_luo_shyy_chi(0.75, omega, true).expect("chi");
    assert!((chi - 0.5 / (0.8 + 0.5)).abs() < 1e-15);
    let chi_no_up = mei_luo_shyy_chi(0.75, omega, false).expect("chi");
    assert_eq!(chi, chi_no_up, "q ≥ 1/2 does not use the upstream node");
    let chi = mei_luo_shyy_chi(0.25, omega, true).expect("chi");
    assert!((chi - omega * (-0.5) / (1.0 - 2.0 * omega)).abs() < 1e-15);
    let chi = mei_luo_shyy_chi(0.25, omega, false).expect("chi");
    assert!((chi - omega * (-0.5) / (1.0 - omega)).abs() < 1e-15);
    // τ = 2 is singular for MLS with q < 1/2, τ = 1 for the original FH.
    assert!(mei_luo_shyy_chi(0.25, 0.5, true).is_err());
    assert!(mei_luo_shyy_chi(0.25, 1.0, false).is_err());
    assert!(mei_luo_shyy_chi(0.75, 2.5, true).is_err());
}

#[test]
fn test_fictitious_equilibrium_matches_equilibrium_when_u_bf_equals_u_f() {
    let lattice = CurvedLattice::d2q9();
    let u = [0.03, 0.01, 0.0];
    let link = d2q9_link(D2Q9_NE, 0.4, [0.0; 3], 1.1);
    let f_star = fictitious_equilibrium(&link, u, u);
    assert!((f_star - lattice.equilibrium(D2Q9_NE, 1.1, u)).abs() < 1e-15);
}

// --- Quadratic helper, half-way BB, moving wall ------------------------------

#[test]
fn test_quadratic_interpolation_bc_nodes() {
    assert!((quadratic_interpolation_bc(1.0, 2.0, 3.0, 0.0) - 1.0).abs() < 1e-12);
    assert!((quadratic_interpolation_bc(1.0, 2.0, 3.0, 1.0) - 2.0).abs() < 1e-12);
    assert!((quadratic_interpolation_bc(1.0, 2.0, 3.0, 2.0) - 3.0).abs() < 1e-12);
    // Exact for a parabola: p(t) = t² → p(0.5) = 0.25
    assert!((quadratic_interpolation_bc(0.0, 1.0, 4.0, 0.5) - 0.25).abs() < 1e-12);
}

#[test]
fn test_half_way_bounce_back_returns_outgoing() {
    let f = vec![0.1, 0.2, 0.3];
    assert!((half_way_bounce_back(&f, 1) - 0.2).abs() < 1e-12);
    assert_eq!(half_way_bounce_back(&f, 5), 0.0);
}

#[test]
fn test_moving_wall_bc_uses_true_opposite_and_dot_product() {
    let lattice = CurvedLattice::d2q9();
    let f_post: Vec<f64> = (0..9).map(|i| 0.1 + 0.01 * i as f64).collect();
    let mut f = vec![0.0; 9];
    // Link E (dir 1) into a wall moving with u = (0.1, 0.2): c·u = 0.1.
    let v = moving_wall_bc(&f_post, &mut f, &lattice, 1, [0.1, 0.2, 0.0], 1.0).expect("bc");
    let expected = f_post[1] - 2.0 * (1.0 / 9.0) * 0.1 / CS2;
    assert!((v - expected).abs() < 1e-15);
    assert!(
        (f[3] - expected).abs() < 1e-15,
        "W (3) is the opposite of E (1)"
    );
    assert_eq!(f[7], 0.0, "direction 7 must not be touched");
    // Link N (dir 2) with the same wall: c·u = 0.2.
    let v = moving_wall_bc(&f_post, &mut f, &lattice, 2, [0.1, 0.2, 0.0], 1.0).expect("bc");
    assert!((v - (f_post[2] - 2.0 * (1.0 / 9.0) * 0.2 / CS2)).abs() < 1e-15);
    assert!((f[4] - v).abs() < 1e-15);
    assert!(moving_wall_bc(&f_post, &mut f, &lattice, 9, [0.0; 3], 1.0).is_err());
    assert!(moving_wall_bc(&f_post[..4], &mut f, &lattice, 1, [0.0; 3], 1.0).is_err());
}

// --- direct forcing ----------------------------------------------------------

#[test]
fn test_direct_forcing_blends() {
    let mut u = [0.0; 3];
    direct_forcing(&mut u, [1.0, 2.0, 3.0], 1.0);
    assert_eq!(u, [1.0, 2.0, 3.0]);
    let mut u = [1.0, 2.0, 3.0];
    direct_forcing(&mut u, [0.0; 3], 0.0);
    assert_eq!(u, [1.0, 2.0, 3.0]);
    let mut u = [0.0; 3];
    direct_forcing(&mut u, [2.0, 0.0, 0.0], 0.5);
    assert!((u[0] - 1.0).abs() < 1e-12);
}

// --- equilibrium and SDF helpers ---------------------------------------------

#[test]
fn test_d2q9_equilibrium_moments() {
    let feq = d2q9_equilibrium(1.0, 0.1, 0.05);
    let (rho, j) = CurvedLattice::d2q9().moments(&feq).expect("moments");
    assert!((rho - 1.0).abs() < 1e-12);
    assert!((j[0] - 0.1).abs() < 1e-12);
    assert!((j[1] - 0.05).abs() < 1e-12);
}

#[test]
fn test_sdf_helpers() {
    assert!((sphere_sdf([3.0, 0.0, 0.0], [0.0; 3], 1.0) - 2.0).abs() < 1e-12);
    assert!((sphere_sdf([0.0; 3], [0.0; 3], 1.0) + 1.0).abs() < 1e-12);
    assert!(box_sdf([0.0; 3], [0.0; 3], [1.0; 3]) <= 0.0);
    let d = box_sdf([2.0, 2.0, 2.0], [0.0; 3], [1.0; 3]);
    assert!((d - 3.0_f64.sqrt()).abs() < 1e-12);
}

// --- CurvedWallNode / LatticeDomain ------------------------------------------

#[test]
fn test_curved_wall_node_construction() {
    let n = CurvedWallNode::new(3, 7, 2, 0.4)
        .expect("node")
        .with_wall_velocity([0.1, 0.0, 0.0])
        .with_upstream(Some(1), None);
    assert_eq!((n.fluid_idx, n.wall_idx, n.dir), (3, 7, 2));
    assert!((n.q - 0.4).abs() < 1e-15);
    assert_eq!(n.u_wall, [0.1, 0.0, 0.0]);
    assert_eq!((n.upstream, n.upstream2), (Some(1), None));
    assert!(CurvedWallNode::new(0, 1, 1, 0.0).is_err());
    assert!(CurvedWallNode::new(0, 1, 1, 1.2).is_err());
}

#[test]
fn test_lattice_domain_indexing_and_periodicity() {
    let d = LatticeDomain::new(4, 3, 2, 0.5).expect("domain");
    let idx = d.linear_index(3, 2, 1);
    assert_eq!(d.coords(idx), [3, 2, 1]);
    assert_eq!(d.position([3, 2, 1]), [1.5, 1.0, 0.5]);
    assert_eq!(d.neighbour([3, 0, 0], [1, 0, 0]), None);
    let p = d.with_periodic([true, false, false]);
    assert_eq!(p.neighbour([3, 0, 0], [1, 1, 0]), Some([0, 1, 0]));
    assert_eq!(p.neighbour([0, 0, 0], [-1, -1, 0]), None);
    assert!(LatticeDomain::new(0, 3, 1, 1.0).is_err());
    assert!(LatticeDomain::new(2, 3, 1, 0.0).is_err());
}

// --- from_geometry -------------------------------------------------------------

#[test]
fn test_from_geometry_flat_wall_scans_diagonals_and_periodic_links() {
    let lattice = CurvedLattice::d2q9();
    let domain = LatticeDomain::new(4, 5, 1, 1.0)
        .expect("domain")
        .with_periodic([true, false, false]);
    // Wall at y = 0.7: row 0 solid, rows ≥ 1 fluid; every downward link from
    // row 1 is cut at q = 0.3 (also the diagonal ones).
    let set = CurvedBoundarySet::from_geometry(&|p: [f64; 3]| p[1] - 0.7, &domain, &lattice);
    assert_eq!(set.method, CurvedBcMethod::Bouzidi);
    assert_eq!(
        set.nodes.len(),
        4 * 3,
        "S, SW and SE links of each bottom node"
    );
    for n in &set.nodes {
        let [ix, iy, _] = domain.coords(n.fluid_idx);
        assert_eq!(iy, 1);
        assert!([4, 7, 8].contains(&n.dir), "dir={}", n.dir);
        assert!((n.q - 0.3).abs() < 1e-12, "q={}", n.q);
        assert_eq!(domain.coords(n.wall_idx)[1], 0);
        let c = lattice.velocities()[n.dir];
        let back = domain
            .neighbour([ix, iy, 0], [-c[0], -c[1], 0])
            .map(|m| domain.linear_index(m[0], m[1], m[2]));
        assert_eq!(n.upstream, back);
        assert!(n.upstream2.is_some());
    }
}

#[test]
fn test_from_geometry_circle_classifies_fluid_and_solid() {
    let lattice = CurvedLattice::d2q9();
    let domain = LatticeDomain::new(12, 12, 1, 1.0).expect("domain");
    let center = [5.3, 5.6, 0.0];
    let radius = 3.1;
    let set = CurvedBoundarySet::from_geometry(
        &|p: [f64; 3]| sphere_sdf(p, center, radius),
        &domain,
        &lattice,
    );
    assert!(!set.nodes.is_empty());
    let mut saw_diagonal = false;
    for n in &set.nodes {
        let pf = domain.position(domain.coords(n.fluid_idx));
        let pw = domain.position(domain.coords(n.wall_idx));
        assert!(sphere_sdf(pf, center, radius) > 0.0);
        assert!(sphere_sdf(pw, center, radius) <= 0.0);
        let c = lattice.velocity_f64(n.dir).expect("c");
        let hit = [pf[0] + n.q * c[0], pf[1] + n.q * c[1], 0.0];
        assert!(
            sphere_sdf(hit, center, radius).abs() < 1e-9,
            "q must hit the surface"
        );
        saw_diagonal |= n.dir >= 5;
    }
    assert!(saw_diagonal, "diagonal links must be detected");
}

// --- RotatingCylinder ----------------------------------------------------------

#[test]
fn test_rotating_cylinder_kinematics() {
    let cyl = RotatingCylinder::new([0.0, 0.0], 1.0, 2.0);
    assert_eq!(cyl.wall_velocity([1.0, 0.0]), [0.0, 2.0, 0.0]);
    assert_eq!(cyl.wall_velocity([0.0, 1.0]), [-2.0, 0.0, 0.0]);
    assert!(cyl.contains([0.5, 0.5]));
    assert!(!cyl.contains([5.0, 0.0]));
    assert!(cyl.signed_distance([3.0, 0.0]) > 0.0);
    assert!(cyl.signed_distance([0.2, 0.0]) < 0.0);
    let q = cyl.ray_q([-1.5, 0.0], [1.0, 0.0]).expect("hit");
    assert!((q - 0.5).abs() < 1e-12);
    assert_eq!(cyl.ray_q([-5.0, 0.0], [1.0, 0.0]), None);
}

#[test]
fn test_rotating_cylinder_boundary_set_orientation_and_wall_velocity() {
    let omega_rot = 0.01;
    let cyl = RotatingCylinder::new([10.2, 9.7], 4.3, omega_rot);
    let set = cyl.boundary_set(21, 20).expect("set");
    let lattice = CurvedLattice::d2q9();
    let domain = LatticeDomain::new(21, 20, 1, 1.0).expect("domain");
    assert!(!set.nodes.is_empty());
    for n in &set.nodes {
        let pf = domain.position(domain.coords(n.fluid_idx));
        let pw = domain.position(domain.coords(n.wall_idx));
        assert!(!cyl.contains([pf[0], pf[1]]), "fluid node inside cylinder");
        assert!(cyl.contains([pw[0], pw[1]]), "wall node outside cylinder");
        let c = lattice.velocity_f64(n.dir).expect("c");
        let q_exact = cyl.ray_q([pf[0], pf[1]], [c[0], c[1]]).expect("ray hit");
        assert!((n.q - q_exact).abs() < 1e-9, "q={} exact={q_exact}", n.q);
        let speed = n.u_wall[0].hypot(n.u_wall[1]);
        assert!(
            (speed - omega_rot * cyl.radius).abs() < 1e-9,
            "|u_w|={speed}"
        );
        // Tangential: u_w ⟂ (hit − centre)
        let hit = [pf[0] + n.q * c[0] - 10.2, pf[1] + n.q * c[1] - 9.7];
        assert!((n.u_wall[0] * hit[0] + n.u_wall[1] * hit[1]).abs() < 1e-9);
    }
}

// --- CurvedBoundarySet::apply ----------------------------------------------------

fn uniform_field(n: usize, rho: f64, u: [f64; 3]) -> Vec<Vec<f64>> {
    let l = CurvedLattice::d2q9();
    (0..n)
        .map(|_| (0..9).map(|i| l.equilibrium(i, rho, u)).collect())
        .collect()
}

#[test]
fn test_apply_rejects_bad_inputs() {
    let lattice = CurvedLattice::d2q9();
    let mut set = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
    set.add_node(CurvedWallNode::new(0, 1, 4, 0.4).expect("node"));
    let f_post = uniform_field(2, 1.0, [0.0; 3]);
    let mut f = uniform_field(3, 1.0, [0.0; 3]);
    assert!(set.apply(&f_post, &mut f, &lattice, 1.0, None).is_err());
    let mut f = uniform_field(2, 1.0, [0.0; 3]);
    assert!(set.apply(&f_post, &mut f, &lattice, 2.0, None).is_err());
    let mut bad = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
    bad.add_node(CurvedWallNode::new(0, 1, 9, 0.4).expect("node"));
    assert!(bad.apply(&f_post, &mut f, &lattice, 1.0, None).is_err());
    let mut rest = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
    rest.add_node(CurvedWallNode::new(0, 1, 0, 0.4).expect("node"));
    assert!(rest.apply(&f_post, &mut f, &lattice, 1.0, None).is_err());
    let mut oob = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
    oob.add_node(CurvedWallNode::new(5, 1, 4, 0.4).expect("node"));
    assert!(oob.apply(&f_post, &mut f, &lattice, 1.0, None).is_err());
}

#[test]
fn test_apply_writes_only_the_unknown_population() {
    let lattice = CurvedLattice::d2q9();
    let f_post: Vec<Vec<f64>> = (0..3)
        .map(|n| (0..9).map(|i| 0.05 + 0.01 * (i + 3 * n) as f64).collect())
        .collect();
    let original: Vec<Vec<f64>> = (0..3).map(|_| vec![0.111; 9]).collect();
    for method in [
        CurvedBcMethod::Bouzidi,
        CurvedBcMethod::QuadraticInterpolation,
        CurvedBcMethod::FilippovaHanel,
    ] {
        let mut set = CurvedBoundarySet::new(method, [0.0; 3]);
        // Fluid node 0, link S (dir 4) into the wall, upstream nodes 1 and 2.
        set.add_node(
            CurvedWallNode::new(0, 9, 4, 0.3)
                .expect("node")
                .with_wall_velocity([0.02, 0.0, 0.0])
                .with_upstream(Some(1), Some(2)),
        );
        let mut f = original.clone();
        set.apply(&f_post, &mut f, &lattice, 1.25, None)
            .expect("apply");
        for (n, row) in f.iter().enumerate() {
            for (i, v) in row.iter().enumerate() {
                if n == 0 && i == 2 {
                    assert!((v - 0.111).abs() > 1e-6, "{method:?}: N must be rebuilt");
                } else {
                    assert_eq!(*v, 0.111, "{method:?}: f[{n}][{i}] changed");
                }
            }
        }
        // Linear BFL value for q = 0.3: 2q f_S(x_f) + (1−2q) f_S(x_ff) − M.
        if method == CurvedBcMethod::Bouzidi {
            let rho: f64 = f_post[0].iter().sum();
            let m = bouzidi_wall_correction(rho, 1.0 / 9.0, [0.0, -1.0, 0.0], [0.02, 0.0, 0.0]);
            let expected = 0.6 * f_post[0][4] + 0.4 * f_post[1][4] - m;
            assert!((f[0][2] - expected).abs() < 1e-14);
        }
    }
}

#[test]
fn test_apply_moment_based_imposes_interpolated_velocity() {
    let lattice = CurvedLattice::d2q9();
    // Column of three fluid nodes (0 nearest the wall) above a bottom wall at
    // q = 0.4 below node 0, wall moving with u_w = (0.05, 0).
    let u_up = [0.08, 0.0, 0.0];
    let mut f_post = uniform_field(3, 1.0, [0.06, 0.0, 0.0]);
    f_post[1] = uniform_field(1, 1.0, u_up).remove(0);
    let q = 0.4;
    let u_w = [0.05, 0.0, 0.0];
    let mut set = CurvedBoundarySet::new(CurvedBcMethod::MomentBased, [0.0; 3]);
    // S, SW, SE links from node 0; their upstream nodes all carry u_up.
    for dir in [4, 7, 8] {
        set.add_node(
            CurvedWallNode::new(0, 9, dir, q)
                .expect("node")
                .with_wall_velocity(u_w)
                .with_upstream(Some(1), Some(2)),
        );
    }
    let mut f = uniform_field(3, 1.0, [0.03, 0.0, 0.0]);
    set.apply(&f_post, &mut f, &lattice, 1.0, None)
        .expect("apply");
    let (rho, j) = lattice.moments(&f[0]).expect("moments");
    let u_t = (u_w[0] + q * u_up[0]) / (1.0 + q);
    assert!(
        (j[0] / rho - u_t).abs() < 1e-13,
        "ux={} target={u_t}",
        j[0] / rho
    );
    assert!((j[1] / rho).abs() < 1e-13, "uy={}", j[1] / rho);
}

#[test]
fn test_symmetric_pinv_solve_full_and_rank_deficient() {
    let m = [[4.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 0.0]];
    let r = [1.0, 2.0, 0.0];
    let x = symmetric_pinv_solve(&m, &r, 2);
    assert!((4.0 * x[0] + x[1] - 1.0).abs() < 1e-12);
    assert!((x[0] + 3.0 * x[1] - 2.0).abs() < 1e-12);
    // Rank one: M = v vᵀ with v = (1, 1); pinv solution of M x = (2, 2)
    // is x = (1, 1) (minimal norm).
    let m = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0; 3]];
    let x = symmetric_pinv_solve(&m, &[2.0, 2.0, 0.0], 2);
    assert!(
        (x[0] - 1.0).abs() < 1e-12 && (x[1] - 1.0).abs() < 1e-12,
        "{x:?}"
    );
}

#[test]
fn test_apply_fh_without_upstream_falls_back_at_tau_one() {
    // ω = 1 makes the original Filippova–Hänel weighting singular; a q < 1/2
    // link without an upstream node must use the bounce-back fallback.
    let lattice = CurvedLattice::d2q9();
    let f_post = uniform_field(2, 1.0, [0.01, 0.0, 0.0]);
    let mut f = uniform_field(2, 1.0, [0.0; 3]);
    let u_w = [0.02, 0.0, 0.0];
    let mut set = CurvedBoundarySet::new(CurvedBcMethod::FilippovaHanel, [0.0; 3]);
    set.add_node(
        CurvedWallNode::new(0, 1, 4, 0.3)
            .expect("node")
            .with_wall_velocity(u_w),
    );
    set.apply(&f_post, &mut f, &lattice, 1.0, None)
        .expect("fallback must not error");
    let rho: f64 = f_post[0].iter().sum();
    // c_S · u_w = 0 → plain half-way bounce-back of f̃_S.
    let expected = f_post[0][4] - bouzidi_wall_correction(rho, 1.0 / 9.0, [0.0, -1.0, 0.0], u_w);
    assert!((f[0][2] - expected).abs() < 1e-15);
    // With an upstream node MLS is used and τ = 1 is regular.
    let mut f_post_up = f_post.clone();
    f_post_up[1] = uniform_field(1, 1.0, [0.01, -0.005, 0.0]).remove(0);
    let mut mls = CurvedBoundarySet::new(CurvedBcMethod::FilippovaHanel, [0.0; 3]);
    mls.add_node(
        CurvedWallNode::new(0, 1, 4, 0.3)
            .expect("node")
            .with_upstream(Some(1), None),
    );
    let mut g = uniform_field(2, 1.0, [0.0; 3]);
    mls.apply(&f_post_up, &mut g, &lattice, 1.0, None)
        .expect("mls");
    let link = WallLink::new(0.3, 1.0 / 9.0, [0.0, -1.0, 0.0], [0.0; 3], rho).expect("link");
    let expected_mls = filippova_hanel(
        f_post_up[0][4],
        &link,
        [0.01, 0.0, 0.0],
        Some([0.01, -0.005, 0.0]),
        1.0,
    )
    .expect("kernel");
    assert!((g[0][2] - expected_mls).abs() < 1e-14);
    assert!(
        (g[0][2] - f[0][2]).abs() > 1e-6,
        "MLS must differ from the fallback"
    );
}
