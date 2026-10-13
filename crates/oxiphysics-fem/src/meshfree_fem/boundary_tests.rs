// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the essential-boundary-condition methods for MLS approximations.
//!
//! The decisive properties: (1) all methods act on the *field* `u_h(x_b)`, i.e.
//! on the rows `φ_J(x_b)`, never on a single nodal parameter; (2) Nitsche's
//! method contains the consistency terms built from `∂φ/∂x` and the outward
//! normal, so it is genuinely different from (and symmetric unlike) a penalty.

use super::functions::{Point2, plane_stress_d};
use super::*;

const KERNEL: KernelType = KernelType::QuarticSpline;

fn square_cloud(n: usize, support_factor: f64) -> SupportDomain {
    SupportDomain::new(
        generate_regular_nodes(0.0, 0.0, 1.0, 1.0, n, n, support_factor),
        1.0,
    )
}

fn material() -> EfgMaterial {
    EfgMaterial {
        young: 1.0e3,
        poisson: 0.3,
        thickness: 0.5,
    }
}

fn zero_system(domain: &SupportDomain) -> (Vec<f64>, Vec<f64>, usize) {
    let ndof = 2 * domain.node_count();
    (vec![0.0; ndof * ndof], vec![0.0; ndof], ndof)
}

/// Nodal parameters `û_I = u_ex(x_I)` of a field.
fn nodal_values(domain: &SupportDomain, field: impl Fn(Point2) -> [f64; 2]) -> Vec<f64> {
    domain
        .nodes
        .iter()
        .flat_map(|n| field(n.pos))
        .collect::<Vec<f64>>()
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}

// ---------------------------------------------------------------------------
// Operator level
// ---------------------------------------------------------------------------

#[test]
fn penalty_couples_all_nodes_in_range_through_the_shape_functions() {
    let domain = square_cloud(5, 2.5);
    let (mut k, mut f, ndof) = zero_system(&domain);
    let bc = EssentialBc::new([0.37, 0.0], [0.0, -1.0], 0.1, 0, 0.2);
    apply_penalty_bc(&mut k, &mut f, &domain, KERNEL, &[bc], 1e6).expect("valid condition");

    let shape = mls_shape_2d(bc.pos, &domain, KERNEL).expect("shape functions exist");
    assert!(shape.neighbours.len() >= 6, "several nodes are in range");
    let mut touched = 0;
    for (a, &ia) in shape.neighbours.iter().enumerate() {
        assert!((f[2 * ia] - 1e6 * 0.1 * shape.phi[a] * 0.2).abs() < 1e-9);
        assert_eq!(f[2 * ia + 1], 0.0, "y component is not constrained");
        for (b, &ib) in shape.neighbours.iter().enumerate() {
            let expected = 1e6 * 0.1 * shape.phi[a] * shape.phi[b];
            let got = k[2 * ia * ndof + 2 * ib];
            assert!((got - expected).abs() <= 1e-12 * expected.abs().max(1.0));
            if ia != ib && expected.abs() > 0.0 {
                touched += 1;
            }
        }
    }
    // The old code only added `alpha` to one diagonal entry.
    assert!(touched > 10, "off-diagonal couplings: {touched}");
    // nothing in the y rows / columns
    for row in 0..ndof {
        for col in 0..ndof {
            if row % 2 == 1 || col % 2 == 1 {
                assert_eq!(k[row * ndof + col], 0.0);
            }
        }
    }
}

#[test]
fn lagrange_rows_hold_shape_function_values() {
    let domain = square_cloud(5, 2.5);
    let (k, f, ndof) = zero_system(&domain);
    let bcs = [
        EssentialBc::new([0.37, 0.0], [0.0, -1.0], 0.1, 0, 0.2),
        EssentialBc::new([0.0, 0.61], [-1.0, 0.0], 0.1, 1, -0.4),
    ];
    let (k_aug, f_aug) =
        apply_lagrange_bc(&k, &f, &domain, KERNEL, &bcs).expect("valid conditions");
    let n_aug = ndof + 2;
    assert_eq!(k_aug.len(), n_aug * n_aug);
    assert_eq!(f_aug.len(), n_aug);
    for (c, bc) in bcs.iter().enumerate() {
        let shape = mls_shape_2d(bc.pos, &domain, KERNEL).expect("shape functions exist");
        let row = &k_aug[(ndof + c) * n_aug..(ndof + c + 1) * n_aug];
        let sum: f64 = row.iter().sum();
        // the row of an MLS constraint sums to the partition of unity
        assert!((sum - 1.0).abs() < 1e-10, "row sum {sum}");
        for (a, &node) in shape.neighbours.iter().enumerate() {
            assert!((row[2 * node + bc.dof] - shape.phi[a]).abs() < 1e-14);
            // symmetric placement in the column block
            assert_eq!(
                k_aug[(2 * node + bc.dof) * n_aug + ndof + c],
                row[2 * node + bc.dof]
            );
        }
        assert!((f_aug[ndof + c] - bc.value).abs() < 1e-15);
    }
}

/// Independent evaluation of the Nitsche contribution with explicit N, D, B
/// matrices (no shared code with the production traction operator).
fn nitsche_reference(
    domain: &SupportDomain,
    mat: &EfgMaterial,
    bc: &EssentialBc,
    beta: f64,
) -> (Vec<f64>, Vec<f64>) {
    let shape = mls_shape_2d(bc.pos, domain, KERNEL).expect("shape functions exist");
    let d = plane_stress_d(mat.young, mat.poisson);
    let n_nodes = domain.node_count();
    let ndof = 2 * n_nodes;
    let mut k = vec![0.0; ndof * ndof];
    let mut f = vec![0.0; ndof];
    // N (2 x 3): traction from (sxx, syy, txy)
    let n_op = [
        [bc.normal[0], 0.0, bc.normal[1]],
        [0.0, bc.normal[1], bc.normal[0]],
    ];
    let spacing = local_node_spacing(domain, bc.pos);
    let gamma = beta * d[0][0] / spacing;
    // T_I = N D B_I
    let traction: Vec<[[f64; 2]; 2]> = (0..shape.neighbours.len())
        .map(|a| {
            let b = [
                [shape.dphi_dx[a], 0.0],
                [0.0, shape.dphi_dy[a]],
                [shape.dphi_dy[a], shape.dphi_dx[a]],
            ];
            let mut db = [[0.0; 2]; 3];
            for r in 0..3 {
                for c in 0..2 {
                    db[r][c] = (0..3).map(|s| d[r][s] * b[s][c]).sum();
                }
            }
            let mut t = [[0.0; 2]; 2];
            for r in 0..2 {
                for c in 0..2 {
                    t[r][c] = (0..3).map(|s| n_op[r][s] * db[s][c]).sum();
                }
            }
            t
        })
        .collect();
    let scale = bc.weight * mat.thickness;
    let c = bc.dof;
    for (a, &ia) in shape.neighbours.iter().enumerate() {
        for (b, &ib) in shape.neighbours.iter().enumerate() {
            for i in 0..2 {
                for j in 0..2 {
                    let mut value = 0.0;
                    if i == c {
                        value -= shape.phi[a] * traction[b][c][j];
                    }
                    if j == c {
                        value -= traction[a][c][i] * shape.phi[b];
                    }
                    if i == c && j == c {
                        value += gamma * shape.phi[a] * shape.phi[b];
                    }
                    k[(2 * ia + i) * ndof + 2 * ib + j] += scale * value;
                }
            }
        }
        for i in 0..2 {
            let mut value = -traction[a][c][i] * bc.value;
            if i == c {
                value += gamma * shape.phi[a] * bc.value;
            }
            f[2 * ia + i] += scale * value;
        }
    }
    (k, f)
}

#[test]
fn nitsche_matches_independent_reference_and_is_symmetric() {
    let domain = square_cloud(6, 2.5);
    let mat = material();
    let bcs = [
        EssentialBc::new([0.0, 0.33], [-1.0, 0.0], 0.05, 0, 0.01),
        EssentialBc::new([0.0, 0.33], [-1.0, 0.0], 0.05, 1, -0.02),
        EssentialBc::new([0.58, 0.0], [0.0, -1.0], 0.07, 0, 0.0),
        EssentialBc::new([1.0, 1.0], [0.6, 0.8], 0.03, 1, 0.005),
    ];
    let beta = 40.0;
    let (mut k, mut f, ndof) = zero_system(&domain);
    apply_nitsche_bc(&mut k, &mut f, &domain, KERNEL, &mat, &bcs, beta).expect("valid conditions");

    let mut k_ref = vec![0.0; ndof * ndof];
    let mut f_ref = vec![0.0; ndof];
    for bc in &bcs {
        let (kb, fb) = nitsche_reference(&domain, &mat, bc, beta);
        for (acc, v) in k_ref.iter_mut().zip(kb.iter()) {
            *acc += v;
        }
        for (acc, v) in f_ref.iter_mut().zip(fb.iter()) {
            *acc += v;
        }
    }
    let scale = max_abs(&k_ref);
    assert!(scale > 0.0);
    for (a, b) in k.iter().zip(k_ref.iter()) {
        assert!((a - b).abs() < 1e-11 * scale, "K mismatch {a} vs {b}");
    }
    let f_scale = max_abs(&f_ref).max(1e-30);
    for (a, b) in f.iter().zip(f_ref.iter()) {
        assert!((a - b).abs() < 1e-11 * f_scale, "f mismatch {a} vs {b}");
    }
    // Nitsche's method is symmetric
    for i in 0..ndof {
        for j in 0..ndof {
            assert!((k[i * ndof + j] - k[j * ndof + i]).abs() < 1e-11 * scale);
        }
    }
}

#[test]
fn nitsche_is_not_the_penalty_method() {
    // The old apply_nitsche_bc added `param` to one diagonal entry and `param * g`
    // to one force entry -- bit for bit the penalty method.  The consistency
    // terms are what make the difference.
    let domain = square_cloud(6, 2.5);
    let mat = material();
    let bc = EssentialBc::new([0.0, 0.33], [-1.0, 0.0], 0.05, 0, 0.01);
    let beta = 40.0;
    let gamma =
        beta * plane_stress_d(mat.young, mat.poisson)[0][0] / local_node_spacing(&domain, bc.pos);

    let (mut k_n, mut f_n, ndof) = zero_system(&domain);
    apply_nitsche_bc(&mut k_n, &mut f_n, &domain, KERNEL, &mat, &[bc], beta).expect("valid");
    let (mut k_p, mut f_p, _) = zero_system(&domain);
    apply_penalty_bc(
        &mut k_p,
        &mut f_p,
        &domain,
        KERNEL,
        &[bc],
        gamma * mat.thickness,
    )
    .expect("valid");

    // same stabilisation block ...
    let difference: Vec<f64> = k_n.iter().zip(k_p.iter()).map(|(a, b)| a - b).collect();
    // ... but Nitsche additionally couples the x and y components through the
    // traction operator, which a penalty on one component never does
    let mut cross_coupling = 0.0_f64;
    for row in 0..ndof {
        for col in 0..ndof {
            if row % 2 != col % 2 {
                cross_coupling = cross_coupling.max(k_n[row * ndof + col].abs());
            }
        }
    }
    assert!(cross_coupling > 0.0, "consistency terms couple u_x and u_y");
    assert!(max_abs(&difference) > 1e-3 * max_abs(&k_n));
    // and the right-hand side carries the adjoint-consistency term -T^T g on
    // the *other* component too
    let shape = mls_shape_2d(bc.pos, &domain, KERNEL).expect("shape functions exist");
    let other_component_load = shape
        .neighbours
        .iter()
        .map(|&node| f_n[2 * node + 1].abs())
        .fold(0.0, f64::max);
    assert!(other_component_load > 0.0);
    assert!(f_p.iter().skip(1).step_by(2).all(|v| *v == 0.0));
}

/// MLS reproduces linear fields, so the boundary term of an exact linear
/// solution is the pure consistency term -∮ φ_I (σ·n) dΓ: stabilisation and
/// adjoint-consistency vanish because u_h(x_b) = g exactly.
#[test]
fn nitsche_residual_of_a_linear_field_is_the_consistency_term() {
    let domain = square_cloud(7, 2.5);
    let mat = material();
    // u = A x + c
    let a = [[0.011, -0.004], [0.006, 0.009]];
    let c = [0.002, -0.001];
    let field = |p: Point2| {
        [
            a[0][0] * p[0] + a[0][1] * p[1] + c[0],
            a[1][0] * p[0] + a[1][1] * p[1] + c[1],
        ]
    };
    let u_hat = nodal_values(&domain, field);
    let d = plane_stress_d(mat.young, mat.poisson);
    let eps = [a[0][0], a[1][1], a[0][1] + a[1][0]];
    let sigma = [
        d[0][0] * eps[0] + d[0][1] * eps[1],
        d[1][0] * eps[0] + d[1][1] * eps[1],
        d[2][2] * eps[2],
    ];

    let mut bcs = Vec::new();
    for bc in essential_bcs_on_segment([0.0, 1.0], [0.0, 0.0], 6, 3, |p| {
        let u = field(p);
        [Some(u[0]), Some(u[1])]
    }) {
        bcs.push(bc);
    }
    for bc in essential_bcs_on_segment([1.0, 0.0], [1.0, 1.0], 6, 3, |p| {
        let u = field(p);
        [Some(u[0]), Some(u[1])]
    }) {
        bcs.push(bc);
    }
    let (mut k, mut f, ndof) = zero_system(&domain);
    apply_nitsche_bc(&mut k, &mut f, &domain, KERNEL, &mat, &bcs, 25.0).expect("valid");
    let residual: Vec<f64> = (0..ndof)
        .map(|i| (0..ndof).map(|j| k[i * ndof + j] * u_hat[j]).sum::<f64>() - f[i])
        .collect();

    let mut expected = vec![0.0; ndof];
    for bc in &bcs {
        // only count each point once: add the component of this condition
        let shape = mls_shape_2d(bc.pos, &domain, KERNEL).expect("shape functions exist");
        let n = bc.normal;
        let traction = [
            n[0] * sigma[0] + n[1] * sigma[2],
            n[0] * sigma[2] + n[1] * sigma[1],
        ];
        for (aidx, &node) in shape.neighbours.iter().enumerate() {
            // -scale * phi_I * t_c, on row (I, c)
            expected[2 * node + bc.dof] -=
                bc.weight * mat.thickness * shape.phi[aidx] * traction[bc.dof];
        }
    }
    let scale = max_abs(&expected);
    assert!(scale > 0.0);
    for (r, e) in residual.iter().zip(expected.iter()) {
        assert!((r - e).abs() < 1e-9 * scale, "residual {r} vs expected {e}");
    }
}

// ---------------------------------------------------------------------------
// Boundary point generation
// ---------------------------------------------------------------------------

#[test]
fn segment_points_carry_length_weights_and_outward_normals() {
    // counter-clockwise traversal of the unit square: the interior is on the left
    let edges = [
        ([0.0, 0.0], [1.0, 0.0], [0.0, -1.0]),
        ([1.0, 0.0], [1.0, 1.0], [1.0, 0.0]),
        ([1.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        ([0.0, 1.0], [0.0, 0.0], [-1.0, 0.0]),
    ];
    for (start, end, outward) in edges {
        let bcs = essential_bcs_on_segment(start, end, 4, 3, |_| [Some(1.0), None]);
        assert_eq!(bcs.len(), 4 * 3, "one component prescribed");
        let total: f64 = bcs.iter().map(|b| b.weight).sum();
        assert!(
            (total - 1.0).abs() < 1e-14,
            "weights sum to the edge length"
        );
        for bc in &bcs {
            assert!((bc.normal[0] - outward[0]).abs() < 1e-14);
            assert!((bc.normal[1] - outward[1]).abs() < 1e-14);
            assert_eq!(bc.dof, 0);
        }
    }
    // both components -> two conditions per point
    let both = essential_bcs_on_segment([0.0, 0.0], [2.0, 0.0], 2, 2, |_| [Some(0.0), Some(0.0)]);
    assert_eq!(both.len(), 2 * 2 * 2);
    // quadrature exactness: integral of x over [0, 2] is 2
    let integral: f64 = both
        .iter()
        .filter(|b| b.dof == 0)
        .map(|b| b.weight * b.pos[0])
        .sum();
    assert!((integral - 2.0).abs() < 1e-14);
    // degenerate input
    assert!(essential_bcs_on_segment([0.0; 2], [0.0; 2], 2, 2, |_| [Some(0.0), None]).is_empty());
    assert!(essential_bcs_on_segment([0.0; 2], [1.0, 0.0], 0, 2, |_| [Some(0.0), None]).is_empty());
    assert!(essential_bcs_on_segment([0.0; 2], [1.0, 0.0], 2, 0, |_| [Some(0.0), None]).is_empty());
}

#[test]
fn local_spacing_of_a_regular_cloud() {
    let domain = square_cloud(9, 2.5);
    let h = 1.0 / 8.0;
    // at an interior node the four nearest neighbours are at distance h
    let at_node = local_node_spacing(&domain, [0.5, 0.5]);
    assert!((at_node - h).abs() < 1e-12, "{at_node}");
    // between two nodes the estimate stays within a factor of two of h
    let between = local_node_spacing(&domain, [0.5 + 0.5 * h, 0.5]);
    assert!(between > 0.5 * h && between < 1.5 * h, "{between}");
    // empty cloud
    let empty = SupportDomain::new(Vec::new(), 1.0);
    assert_eq!(local_node_spacing(&empty, [0.0, 0.0]), 0.0);
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

#[test]
fn malformed_conditions_are_rejected_and_leave_the_system_untouched() {
    let domain = square_cloud(5, 2.5);
    let mat = material();
    let (mut k, mut f, ndof) = zero_system(&domain);
    let good = EssentialBc::new([0.5, 0.0], [0.0, -1.0], 0.1, 0, 0.0);
    let bad_dof = EssentialBc::new([0.5, 0.0], [0.0, -1.0], 0.1, 2, 0.0);
    let no_weight = EssentialBc::new([0.5, 0.0], [0.0, -1.0], 0.0, 0, 0.0);
    let nan_value = EssentialBc::new([0.5, 0.0], [0.0, -1.0], 0.1, 0, f64::NAN);
    let outside = EssentialBc::new([50.0, 50.0], [0.0, -1.0], 0.1, 0, 0.0);
    let bad_normal = EssentialBc::new([0.5, 0.0], [0.0, -3.0], 0.1, 0, 0.0);

    for bc in [bad_dof, no_weight, nan_value, outside] {
        assert!(apply_penalty_bc(&mut k, &mut f, &domain, KERNEL, &[good, bc], 1e3).is_err());
        assert!(
            apply_nitsche_bc(&mut k, &mut f, &domain, KERNEL, &mat, &[good, bc], 10.0).is_err()
        );
        assert!(apply_lagrange_bc(&k, &f, &domain, KERNEL, &[good, bc]).is_err());
    }
    // The valid condition listed *before* the malformed one must not have been
    // applied: a failed call leaves the system exactly as it was.
    assert!(k.iter().all(|v| *v == 0.0), "stiffness was modified");
    assert!(f.iter().all(|v| *v == 0.0), "force was modified");

    // a normal that is not of unit length is only an error for Nitsche
    assert!(apply_nitsche_bc(&mut k, &mut f, &domain, KERNEL, &mat, &[bad_normal], 10.0).is_err());
    assert!(k.iter().all(|v| *v == 0.0) && f.iter().all(|v| *v == 0.0));
    let (mut k_ok, mut f_ok, _) = zero_system(&domain);
    assert!(apply_penalty_bc(&mut k_ok, &mut f_ok, &domain, KERNEL, &[bad_normal], 1e3).is_ok());
    assert!(
        k_ok.iter().any(|v| *v != 0.0),
        "the penalty condition was applied"
    );

    // parameters
    assert!(apply_penalty_bc(&mut k, &mut f, &domain, KERNEL, &[good], 0.0).is_err());
    assert!(apply_penalty_bc(&mut k, &mut f, &domain, KERNEL, &[good], f64::NAN).is_err());
    assert!(apply_nitsche_bc(&mut k, &mut f, &domain, KERNEL, &mat, &[good], -1.0).is_err());
    assert!(k.iter().all(|v| *v == 0.0) && f.iter().all(|v| *v == 0.0));
    // array sizes
    let mut short = vec![0.0; ndof * ndof - 1];
    assert!(apply_penalty_bc(&mut short, &mut f, &domain, KERNEL, &[good], 1e3).is_err());
    assert!(apply_lagrange_bc(&short, &f, &domain, KERNEL, &[good]).is_err());
}

#[test]
fn bc_at_node_uses_the_node_position() {
    let domain = square_cloud(3, 2.5);
    let bc = EssentialBc::at_node(&domain, 4, [1.0, 0.0], 0.5, 1, 0.3).expect("node exists");
    assert_eq!(bc.pos, domain.nodes[4].pos);
    assert_eq!(bc.dof, 1);
    assert!(EssentialBc::at_node(&domain, 9, [1.0, 0.0], 0.5, 1, 0.3).is_none());
}

// ---------------------------------------------------------------------------
// End to end
// ---------------------------------------------------------------------------

/// Boundary points of the unit square along its four edges (counter-clockwise),
/// both components prescribed by `field`.
fn boundary_bcs_of_unit_square(
    n_pieces: usize,
    field: impl Fn(Point2) -> [f64; 2] + Copy,
) -> Vec<EssentialBc> {
    let corners = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut bcs = Vec::new();
    for e in 0..4 {
        bcs.extend(essential_bcs_on_segment(
            corners[e],
            corners[(e + 1) % 4],
            n_pieces,
            3,
            |p| {
                let u = field(p);
                [Some(u[0]), Some(u[1])]
            },
        ));
    }
    bcs
}

/// Point-collocated conditions (one per component) at the given nodes.
fn nodal_bcs(
    domain: &SupportDomain,
    nodes: &[usize],
    prescribed: impl Fn(Point2) -> [Option<f64>; 2],
) -> Vec<EssentialBc> {
    let mut bcs = Vec::new();
    for &node in nodes {
        let values = prescribed(domain.nodes[node].pos);
        for (dof, value) in values.iter().enumerate() {
            if let Some(value) = value
                && let Some(bc) = EssentialBc::at_node(domain, node, [0.0, 0.0], 1.0, dof, *value)
            {
                bcs.push(bc);
            }
        }
    }
    bcs
}

/// Indices of the nodes on the boundary of the unit-square cloud.
fn boundary_nodes(domain: &SupportDomain) -> Vec<usize> {
    domain
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| {
            n.pos
                .iter()
                .any(|c| c.abs() < 1e-12 || (c - 1.0).abs() < 1e-12)
        })
        .map(|(i, _)| i)
        .collect()
}

/// Largest error of `u_h` against `field` over a grid of interior sample points.
fn field_error(
    result: &EfgResult,
    domain: &SupportDomain,
    field: impl Fn(Point2) -> [f64; 2],
) -> f64 {
    let mut worst = 0.0_f64;
    for i in 1..8 {
        for j in 1..8 {
            let p = [i as f64 / 8.0 + 0.013, j as f64 / 8.0 - 0.011];
            let u = evaluate_displacement_2d(p, &result.displacements, domain, KERNEL)
                .expect("inside the cloud");
            let e = field(p);
            worst = worst.max((u[0] - e[0]).abs()).max((u[1] - e[1]).abs());
        }
    }
    worst
}

fn linear_field(p: Point2) -> [f64; 2] {
    [
        0.010 * p[0] + 0.004 * p[1] + 0.002,
        -0.003 * p[0] + 0.007 * p[1] - 0.001,
    ]
}

/// Largest amplitude of [`linear_field`] over the unit square.
const LINEAR_FIELD_AMPLITUDE: f64 = 0.016;

#[test]
fn patch_test_all_methods_recover_a_linear_field_up_to_quadrature_error() {
    // The boundary conditions reproduce the linear field exactly; what remains
    // is the error of the 1-point background integration (the discrete
    // integration constraint is not satisfied by the centroid rule).
    let domain = square_cloud(9, 2.5);
    let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 32, 32);
    let mat = material();
    let weak = boundary_bcs_of_unit_square(8, linear_field);
    let nodal = nodal_bcs(&domain, &boundary_nodes(&domain), |p| {
        let u = linear_field(p);
        [Some(u[0]), Some(u[1])]
    });
    // Point collocation at the boundary *nodes* only (Lagrange) constrains the
    // field less than the boundary integrals do, so it gets a looser bound.
    for (name, method, bcs, parameter, tolerance) in [
        ("penalty", BcEnforcementMethod::Penalty, &weak, 1e9, 2e-3),
        (
            "lagrange",
            BcEnforcementMethod::LagrangeMultiplier,
            &nodal,
            0.0,
            8e-3,
        ),
        ("nitsche", BcEnforcementMethod::Nitsche, &weak, 50.0, 2e-3),
    ] {
        let result = efg_analysis_2d(
            &domain,
            &cells,
            &mat,
            KERNEL,
            [0.0, 0.0],
            bcs,
            method,
            parameter,
        )
        .unwrap_or_else(|| panic!("{name}: analysis failed"));
        let error = field_error(&result, &domain, linear_field);
        assert!(
            error < tolerance * LINEAR_FIELD_AMPLITUDE,
            "{name}: max error {error:e}"
        );
    }
}

#[test]
fn patch_test_error_is_quadrature_limited_and_converges_with_the_cells() {
    let domain = square_cloud(9, 2.5);
    let mat = material();
    let weak = boundary_bcs_of_unit_square(8, linear_field);
    let mut previous = f64::INFINITY;
    for n_cells in [8, 16, 32] {
        let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, n_cells, n_cells);
        let result = efg_analysis_2d(
            &domain,
            &cells,
            &mat,
            KERNEL,
            [0.0, 0.0],
            &weak,
            BcEnforcementMethod::Nitsche,
            50.0,
        )
        .expect("analysis succeeds");
        let error = field_error(&result, &domain, linear_field);
        assert!(
            error < 0.7 * previous,
            "{n_cells} cells: {error:e} vs {previous:e}"
        );
        previous = error;
    }
}

/// Rigid motions are reproduced to rounding by every method, whatever the
/// quadrature: the strain is exactly zero, so only the boundary terms act.
#[test]
fn rigid_motion_is_reproduced_to_rounding() {
    let domain = square_cloud(7, 2.5);
    let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 12, 12);
    let mat = material();
    // translation + small rotation about (0.3, 0.6)
    let rigid = |p: Point2| {
        let (dx, dy) = (p[0] - 0.3, p[1] - 0.6);
        [0.004 - 0.02 * dy, -0.003 + 0.02 * dx]
    };
    let weak = boundary_bcs_of_unit_square(6, rigid);
    let nodal = nodal_bcs(&domain, &boundary_nodes(&domain), |p| {
        let u = rigid(p);
        [Some(u[0]), Some(u[1])]
    });
    for (name, method, bcs, parameter) in [
        ("penalty", BcEnforcementMethod::Penalty, &weak, 1e8),
        (
            "lagrange",
            BcEnforcementMethod::LagrangeMultiplier,
            &nodal,
            0.0,
        ),
        ("nitsche", BcEnforcementMethod::Nitsche, &weak, 30.0),
    ] {
        let result = efg_analysis_2d(
            &domain,
            &cells,
            &mat,
            KERNEL,
            [0.0, 0.0],
            bcs,
            method,
            parameter,
        )
        .unwrap_or_else(|| panic!("{name}: analysis failed"));
        let error = field_error(&result, &domain, rigid);
        // Lagrange and Nitsche are exact; the penalty leaves ~ reaction / alpha
        let tolerance = if name == "penalty" { 1e-7 } else { 1e-9 };
        assert!(error < tolerance, "{name}: max error {error:e}");
    }
}

/// Largest `|u_h(x_b) - g(x_b)|` over the boundary points of the conditions.
fn boundary_mismatch(result: &EfgResult, domain: &SupportDomain, bcs: &[EssentialBc]) -> f64 {
    bcs.iter().fold(0.0_f64, |worst, bc| {
        let u = evaluate_displacement_2d(bc.pos, &result.displacements, domain, KERNEL)
            .expect("boundary point is inside the cloud");
        worst.max((u[bc.dof] - bc.value).abs())
    })
}

/// Left edge pushed by a non-linear profile, the rest of the boundary free.
fn clamp_left_edge_with_profile(
    domain: &SupportDomain,
) -> (impl Fn(Point2) -> [Option<f64>; 2] + Copy, Vec<usize>) {
    let profile = |p: Point2| [Some(0.02 * (std::f64::consts::PI * p[1]).sin()), Some(0.0)];
    let left: Vec<usize> = domain
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.pos[0].abs() < 1e-12)
        .map(|(i, _)| i)
        .collect();
    (profile, left)
}

#[test]
fn lagrange_enforces_a_nonlinear_profile_exactly_at_the_constraint_points() {
    let domain = square_cloud(9, 2.5);
    let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 24, 24);
    let mat = material();
    let (profile, left) = clamp_left_edge_with_profile(&domain);
    let bcs = nodal_bcs(&domain, &left, profile);
    let result = efg_analysis_2d(
        &domain,
        &cells,
        &mat,
        KERNEL,
        [0.0, 0.0],
        &bcs,
        BcEnforcementMethod::LagrangeMultiplier,
        0.0,
    )
    .expect("analysis succeeds");
    let mismatch = boundary_mismatch(&result, &domain, &bcs);
    assert!(mismatch < 1e-9, "u_h(x_b) - g = {mismatch:e}");

    // Assigning g(x_I) to the nodal *parameters* (what the old code effectively
    // did) does not make the field take the value g at the node.
    let middle = left[left.len() / 2];
    let x_i = domain.nodes[middle].pos;
    let shape = mls_shape_2d(x_i, &domain, KERNEL).expect("shape functions exist");
    let naive: f64 = shape
        .neighbours
        .iter()
        .zip(shape.phi.iter())
        .map(|(&j, phi)| phi * profile(domain.nodes[j].pos)[0].unwrap_or(0.0))
        .sum();
    let wanted = profile(x_i)[0].unwrap_or(0.0);
    assert!((naive - wanted).abs() > 1e-5, "{naive} vs {wanted}");
}

/// Boundary mismatch of the left-edge profile for a given method and parameter.
fn profile_mismatch(method: BcEnforcementMethod, parameter: f64) -> f64 {
    let domain = square_cloud(9, 2.5);
    let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 24, 24);
    let bcs = essential_bcs_on_segment([0.0, 1.0], [0.0, 0.0], 8, 3, |p| {
        [Some(0.02 * (std::f64::consts::PI * p[1]).sin()), Some(0.0)]
    });
    let result = efg_analysis_2d(
        &domain,
        &cells,
        &material(),
        KERNEL,
        [0.0, 0.0],
        &bcs,
        method,
        parameter,
    )
    .expect("analysis succeeds");
    boundary_mismatch(&result, &domain, &bcs)
}

#[test]
fn penalty_boundary_mismatch_shrinks_with_the_penalty_parameter() {
    let coarse = profile_mismatch(BcEnforcementMethod::Penalty, 1e3);
    let medium = profile_mismatch(BcEnforcementMethod::Penalty, 1e5);
    let fine = profile_mismatch(BcEnforcementMethod::Penalty, 1e7);
    assert!(medium < 0.1 * coarse, "{medium:e} vs {coarse:e}");
    assert!(fine < 0.2 * medium, "{fine:e} vs {medium:e}");
    // what is left is the MLS approximation error of the sin profile
    assert!(fine < 1e-4, "{fine:e}");
}

#[test]
fn nitsche_boundary_mismatch_shrinks_with_beta() {
    let loose = profile_mismatch(BcEnforcementMethod::Nitsche, 10.0);
    let tight = profile_mismatch(BcEnforcementMethod::Nitsche, 100.0);
    let tighter = profile_mismatch(BcEnforcementMethod::Nitsche, 1000.0);
    assert!(tight < 0.25 * loose, "{tight:e} vs {loose:e}");
    // beyond beta ~ 100 only the MLS approximation error of the profile remains
    assert!(tighter < 1.5 * tight, "{tighter:e} vs {tight:e}");
    assert!(tighter < 1e-4, "{tighter:e}");
}

/// With the *same* stabilisation stiffness, the penalty method (what the old
/// `apply_nitsche_bc` amounted to) misses the linear patch field by far more
/// than Nitsche's method, whose consistency terms make it exact for any beta.
#[test]
fn nitsche_is_consistent_where_the_equivalent_penalty_is_not() {
    let domain = square_cloud(9, 2.5);
    let cells = generate_tri_cells(0.0, 0.0, 1.0, 1.0, 32, 32);
    let mat = material();
    let weak = boundary_bcs_of_unit_square(8, linear_field);
    let beta = 5.0;

    let nitsche = efg_analysis_2d(
        &domain,
        &cells,
        &mat,
        KERNEL,
        [0.0, 0.0],
        &weak,
        BcEnforcementMethod::Nitsche,
        beta,
    )
    .expect("analysis succeeds");
    let nitsche_error = field_error(&nitsche, &domain, linear_field);

    // penalty with the stiffness gamma = beta D11 / h of Nitsche's stabilisation
    let h = 1.0 / 8.0;
    let gamma = beta * plane_stress_d(mat.young, mat.poisson)[0][0] / h;
    let penalty = efg_analysis_2d(
        &domain,
        &cells,
        &mat,
        KERNEL,
        [0.0, 0.0],
        &weak,
        BcEnforcementMethod::Penalty,
        gamma * mat.thickness,
    )
    .expect("analysis succeeds");
    let penalty_error = field_error(&penalty, &domain, linear_field);
    assert!(nitsche_error < 4e-3 * LINEAR_FIELD_AMPLITUDE);
    assert!(
        penalty_error > 4.0 * nitsche_error,
        "penalty {penalty_error:e} vs nitsche {nitsche_error:e}"
    );
}

#[test]
fn bar_under_body_force_matches_the_analytical_solution() {
    // Bar [0,1] x [0,0.5], clamped at x = 0, loaded by a body force f in x with
    // nu = 0: sigma_xx = f (L - x), u_x = f (L x - x^2/2) / E, u_y = 0.
    let (length, height, young, force) = (1.0, 0.5, 1.0e3, 10.0);
    let domain = SupportDomain::new(
        generate_regular_nodes(0.0, 0.0, length, height, 11, 6, 2.5),
        1.0,
    );
    let cells = generate_tri_cells(0.0, 0.0, length, height, 20, 10);
    let mat = EfgMaterial {
        young,
        poisson: 0.0,
        thickness: 1.0,
    };
    let exact = |x: f64| force * (length * x - 0.5 * x * x) / young;

    let weak =
        essential_bcs_on_segment([0.0, height], [0.0, 0.0], 5, 3, |_| [Some(0.0), Some(0.0)]);
    let left: Vec<usize> = domain
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.pos[0].abs() < 1e-12)
        .map(|(i, _)| i)
        .collect();
    let nodal = nodal_bcs(&domain, &left, |_| [Some(0.0), Some(0.0)]);

    for (name, method, bcs, parameter) in [
        ("penalty", BcEnforcementMethod::Penalty, &weak, 1e8),
        (
            "lagrange",
            BcEnforcementMethod::LagrangeMultiplier,
            &nodal,
            0.0,
        ),
        ("nitsche", BcEnforcementMethod::Nitsche, &weak, 50.0),
    ] {
        let result = efg_analysis_2d(
            &domain,
            &cells,
            &mat,
            KERNEL,
            [force, 0.0],
            bcs,
            method,
            parameter,
        )
        .unwrap_or_else(|| panic!("{name}: analysis failed"));
        let mut worst = 0.0_f64;
        let mut worst_y = 0.0_f64;
        for x in [0.3, 0.5, 0.7, 0.9] {
            for y in [0.1, 0.25, 0.4] {
                let u = evaluate_displacement_2d([x, y], &result.displacements, &domain, KERNEL)
                    .expect("inside the cloud");
                worst = worst.max((u[0] - exact(x)).abs());
                worst_y = worst_y.max(u[1].abs());
            }
        }
        let tip = exact(length);
        assert!(
            worst < 0.02 * tip,
            "{name}: u_x error {worst:e} (tip {tip:e})"
        );
        assert!(worst_y < 0.03 * tip, "{name}: u_y {worst_y:e}");
    }
}
