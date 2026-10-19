// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the boundary element influence matrices.
//!
//! * `kelvin_laplace_3d_grad` (and with it `half_space_laplace_3d_grad`,
//!   `helmholtz_3d_grad`) returned `-dG/dn_y`, so the interior representation
//!   `u(x) = ∮ (G q − dG/dn u) dS` evaluated a constant potential as `−1` and
//!   `solve_laplace_bem` returned the negated flux.
//! * `BoundaryMesh::sphere` wound the triangles of the latitude strips
//!   *inward* (192 of 224 elements of `sphere(1, 8, 16)`) although the doc and
//!   the stored normals claim outward orientation.
//! * `assemble_elastostatic_bem` used a self `G` that is 3.6–4.6× too small,
//!   evaluated the traction kernel with swapped source / field arguments
//!   (negating `H`) and set only the diagonal of the rigid-body block of `H`.

use oxiphysics_fem::boundary_element::{
    BemBc, BemBcType, BoundaryElement, BoundaryMesh, ElementOrder, assemble_elastostatic_bem,
    assemble_laplace_bem, half_space_laplace_3d, half_space_laplace_3d_grad, helmholtz_3d,
    helmholtz_3d_grad, interior_potential, kelvin_displacement_3d, kelvin_laplace_3d,
    kelvin_laplace_3d_grad, solve_laplace_bem,
};
use std::f64::consts::PI;

fn unit(v: [f64; 3]) -> [f64; 3] {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / norm, v[1] / norm, v[2] / norm]
}

fn offset(y: &[f64; 3], n: &[f64; 3], h: f64) -> [f64; 3] {
    [y[0] + h * n[0], y[1] + h * n[1], y[2] + h * n[2]]
}

fn centroid(mesh: &BoundaryMesh, index: usize) -> [f64; 3] {
    mesh.element_centroid(index)
}

// ---------------------------------------------------------------------------
// Kernel signs
// ---------------------------------------------------------------------------

#[test]
fn laplace_double_layer_kernel_is_the_normal_derivative_of_the_fundamental_solution() {
    let x = [0.3, -0.2, 0.5];
    let y = [1.1, 0.4, -0.3];
    let n = unit([0.6, -0.3, 0.7]);
    let h = 1.0e-6;
    let numeric = (kelvin_laplace_3d(&x, &offset(&y, &n, h))
        - kelvin_laplace_3d(&x, &offset(&y, &n, -h)))
        / (2.0 * h);
    let analytic = kelvin_laplace_3d_grad(&x, &y, &n);
    assert!(
        (numeric - analytic).abs() <= 1e-7 * analytic.abs(),
        "FD {numeric} vs kernel {analytic}"
    );
    // A normal pointing away from the source lets G decay: dG/dn < 0.
    let away = unit([y[0] - x[0], y[1] - x[1], y[2] - x[2]]);
    assert!(kelvin_laplace_3d_grad(&x, &y, &away) < 0.0);
    // Pointing towards the source it grows.
    let towards = [-away[0], -away[1], -away[2]];
    assert!(kelvin_laplace_3d_grad(&x, &y, &towards) > 0.0);
    // Magnitude on the axis: |dG/dn| = 1 / (4 pi r^2).
    let r = ((y[0] - x[0]).powi(2) + (y[1] - x[1]).powi(2) + (y[2] - x[2]).powi(2)).sqrt();
    assert!((kelvin_laplace_3d_grad(&x, &y, &away) + 1.0 / (4.0 * PI * r * r)).abs() < 1e-14);
}

#[test]
fn helmholtz_double_layer_kernel_is_the_normal_derivative_of_the_fundamental_solution() {
    let x = [0.3, -0.2, 0.5];
    let y = [1.1, 0.4, -0.3];
    let n = unit([-0.2, 0.9, 0.4]);
    let h = 1.0e-6;
    for k in [0.7, 2.3, 6.1] {
        let (gp_re, gp_im) = helmholtz_3d(&x, &offset(&y, &n, h), k);
        let (gm_re, gm_im) = helmholtz_3d(&x, &offset(&y, &n, -h), k);
        let numeric = ((gp_re - gm_re) / (2.0 * h), (gp_im - gm_im) / (2.0 * h));
        let (re, im) = helmholtz_3d_grad(&x, &y, &n, k);
        let scale = re.abs().max(im.abs());
        assert!(
            (numeric.0 - re).abs() <= 1e-6 * scale,
            "k = {k}: re {} vs {re}",
            numeric.0
        );
        assert!(
            (numeric.1 - im).abs() <= 1e-6 * scale,
            "k = {k}: im {} vs {im}",
            numeric.1
        );
    }
    // k -> 0 recovers the Laplace kernel, with the same sign.
    let (re0, im0) = helmholtz_3d_grad(&x, &y, &n, 0.0);
    assert!((re0 - kelvin_laplace_3d_grad(&x, &y, &n)).abs() < 1e-15);
    assert!(im0.abs() < 1e-15);
}

/// The image term of the free-surface Green's function is differentiated along
/// the *mirrored* normal; on the surface itself the normal derivative of the
/// Neumann function vanishes.
#[test]
fn half_space_kernel_is_the_normal_derivative_of_the_image_green_function() {
    let x = [0.3, -0.2, 0.5];
    let y = [1.1, 0.4, 0.8];
    let h = 1.0e-6;
    for n in [
        unit([0.6, -0.3, 0.7]),
        [0.0, 0.0, 1.0],
        unit([1.0, 1.0, -1.0]),
    ] {
        let numeric = (half_space_laplace_3d(&x, &offset(&y, &n, h))
            - half_space_laplace_3d(&x, &offset(&y, &n, -h)))
            / (2.0 * h);
        let analytic = half_space_laplace_3d_grad(&x, &y, &n);
        assert!(
            (numeric - analytic).abs() <= 1e-6 * analytic.abs().max(1e-3),
            "n = {n:?}: FD {numeric} vs kernel {analytic}"
        );
    }
    // y on the free surface z = 0 with n = e_z: dG_hs/dn = 0.
    let surface = [0.7, -0.4, 0.0];
    let flux = half_space_laplace_3d_grad(&x, &surface, &[0.0, 0.0, 1.0]);
    assert!(flux.abs() < 1e-15, "{flux}");
}

// ---------------------------------------------------------------------------
// Sphere mesh orientation
// ---------------------------------------------------------------------------

#[test]
fn every_sphere_triangle_is_wound_outward_and_agrees_with_its_stored_normal() {
    for (radius, n_lat, n_lon) in [(1.0, 8, 16), (2.5, 6, 12), (0.4, 12, 20), (1.0, 3, 5)] {
        let mesh = BoundaryMesh::sphere(radius, n_lat, n_lon);
        for (index, element) in mesh.elements.iter().enumerate() {
            let [a, b, c] = [
                element.node_ids[0],
                element.node_ids[1],
                element.node_ids[2],
            ];
            let (pa, pb, pc) = (mesh.nodes[a], mesh.nodes[b], mesh.nodes[c]);
            let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
            let ac = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
            let wound = unit([
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ]);
            // The stored normal is the node-order normal ...
            for (wound_axis, stored_axis) in wound.iter().zip(&element.normal) {
                assert!(
                    (wound_axis - stored_axis).abs() < 1e-12,
                    "sphere({radius}, {n_lat}, {n_lon}) element {index}"
                );
            }
            // ... and it points away from the centre.
            let centre = centroid(&mesh, index);
            let outward = element.normal[0] * centre[0]
                + element.normal[1] * centre[1]
                + element.normal[2] * centre[2];
            assert!(
                outward > 0.0,
                "sphere({radius}, {n_lat}, {n_lon}) element {index} is wound inward"
            );
        }
    }
}

#[test]
fn sphere_surface_area_and_closure() {
    let mesh = BoundaryMesh::sphere(1.0, 12, 24);
    let area: f64 = mesh.elements.iter().map(|e| e.area).sum();
    // Inscribed polyhedron: slightly below 4 pi.
    assert!(area < 4.0 * PI && area > 0.97 * 4.0 * PI, "area {area}");
    // A closed surface: sum of n_i A_i = 0 (outward normals integrate to zero).
    let mut total = [0.0_f64; 3];
    for element in &mesh.elements {
        for (sum, component) in total.iter_mut().zip(&element.normal) {
            *sum += component * element.area;
        }
    }
    assert!(total.iter().all(|t| t.abs() < 1e-12), "{total:?}");
}

// ---------------------------------------------------------------------------
// Laplace BEM on the sphere
// ---------------------------------------------------------------------------

fn dirichlet_bcs(values: &[f64]) -> Vec<BemBc> {
    values
        .iter()
        .enumerate()
        .map(|(element_id, &value)| BemBc {
            element_id,
            bc_type: BemBcType::Dirichlet,
            values: vec![value],
        })
        .collect()
}

/// `u ≡ 1`, `q = 0` represents the constant potential 1: the interior formula
/// must return 1 (it returned −1 with the flipped kernel, and 0.87 on top of
/// that with the inward-wound strips).
#[test]
fn interior_representation_reproduces_a_constant_potential() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let n = mesh.num_elements();
    let u = vec![1.0; n];
    let q = vec![0.0; n];
    for point in [[0.0, 0.0, 0.0], [0.3, -0.2, 0.1], [-0.5, 0.1, 0.4]] {
        let value = interior_potential(&mesh, &u, &q, &point);
        assert!((value - 1.0).abs() < 0.05, "u({point:?}) = {value}");
    }
    // Refinement converges towards the exact value.
    let fine = BoundaryMesh::sphere(1.0, 16, 32);
    let value = interior_potential(
        &fine,
        &vec![1.0; fine.num_elements()],
        &vec![0.0; fine.num_elements()],
        &[0.0; 3],
    );
    assert!((value - 1.0).abs() < 0.01, "fine mesh: {value}");
}

/// The harmonic function `u = z` has flux `q = ∂u/∂n = n_z` on the unit sphere.
#[test]
fn dirichlet_problem_returns_the_outward_flux_with_the_correct_sign() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let n = mesh.num_elements();
    let z: Vec<f64> = (0..n).map(|i| centroid(&mesh, i)[2]).collect();
    let (h, g) = assemble_laplace_bem(&mesh);
    let (u, q) = solve_laplace_bem(&h, &g, &dirichlet_bcs(&z));
    assert_eq!(u, z);

    let exact: Vec<f64> = mesh.elements.iter().map(|e| e.normal[2]).collect();
    let error: f64 = q
        .iter()
        .zip(&exact)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    let norm: f64 = exact.iter().map(|b| b * b).sum::<f64>().sqrt();
    assert!(error < 0.10 * norm, "relative L2 error {}", error / norm);
    // Element by element the flux is within 0.2 of the exact n_z.
    for (i, (computed, expected)) in q.iter().zip(&exact).enumerate() {
        assert!(
            (computed - expected).abs() < 0.2,
            "element {i}: q = {computed}, n_z = {expected}"
        );
    }
    // Sign: outward flux on the northern cap, inward on the southern one.
    for (i, element) in mesh.elements.iter().enumerate() {
        if element.normal[2].abs() > 0.5 {
            assert!(q[i] * element.normal[2] > 0.0, "element {i}: q = {}", q[i]);
        }
    }
    // The solved flux integrates to zero over the closed surface (Gauss).
    let net: f64 = q
        .iter()
        .zip(&mesh.elements)
        .map(|(qi, e)| qi * e.area)
        .sum();
    assert!(net.abs() < 0.05, "net flux {net}");
}

/// Dirichlet data on the northern half and Neumann data (from `u = z`) on the
/// southern half: the unknown potential on the south must reproduce `z`.
#[test]
fn mixed_problem_recovers_the_harmonic_potential() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let n = mesh.num_elements();
    let (h, g) = assemble_laplace_bem(&mesh);
    let mut bcs = Vec::with_capacity(n);
    for (i, element) in mesh.elements.iter().enumerate() {
        let c = centroid(&mesh, i);
        if c[2] > 0.0 {
            bcs.push(BemBc {
                element_id: i,
                bc_type: BemBcType::Dirichlet,
                values: vec![c[2]],
            });
        } else {
            bcs.push(BemBc {
                element_id: i,
                bc_type: BemBcType::Neumann,
                values: vec![element.normal[2]],
            });
        }
    }
    let (u, _q) = solve_laplace_bem(&h, &g, &bcs);
    let mut error = 0.0_f64;
    let mut norm = 0.0_f64;
    for (i, &potential) in u.iter().enumerate() {
        let z = centroid(&mesh, i)[2];
        if z <= 0.0 {
            error += (potential - z).powi(2);
            norm += z * z;
        }
    }
    assert!(
        error.sqrt() < 0.03 * norm.sqrt(),
        "relative error {}",
        (error / norm).sqrt()
    );
}

#[test]
fn interior_representation_of_a_linear_harmonic_field() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let n = mesh.num_elements();
    let u: Vec<f64> = (0..n).map(|i| centroid(&mesh, i)[2]).collect();
    let q: Vec<f64> = mesh.elements.iter().map(|e| e.normal[2]).collect();
    let point = [0.2, 0.1, 0.3];
    let value = interior_potential(&mesh, &u, &q, &point);
    assert!(
        (value - 0.3).abs() < 0.02,
        "u({point:?}) = {value}, exact 0.3"
    );
}

// ---------------------------------------------------------------------------
// Elastostatic BEM
// ---------------------------------------------------------------------------

const SHEAR: f64 = 0.8;
const NU: f64 = 0.3;

/// A mesh holding a single element of the given area and normal.
fn single_element(area: f64, normal: [f64; 3]) -> BoundaryMesh {
    BoundaryMesh {
        nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        elements: vec![BoundaryElement {
            node_ids: vec![0, 1, 2],
            normal,
            area,
            order: ElementOrder::Constant,
        }],
    }
}

/// `∫∫ U_ij(0, y) dS(y)` over the disk of radius `radius` in the plane
/// `z = 0`, by composite Gauss-Legendre in `r` and the trapezoid rule in `θ`.
fn disk_integral_of_u(radius: f64, i: usize, j: usize) -> f64 {
    let nodes = [
        (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
        (-0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (0.0, 0.568_888_888_888_889),
        (0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_1),
    ];
    let panels = 8;
    let angles = 64;
    let origin = [0.0, 0.0, 0.0];
    let mut total = 0.0;
    for panel in 0..panels {
        let r0 = radius * panel as f64 / panels as f64;
        let width = radius / panels as f64;
        for &(node, weight) in &nodes {
            let r = r0 + 0.5 * width * (node + 1.0);
            for k in 0..angles {
                let theta = 2.0 * PI * (k as f64 + 0.5) / angles as f64;
                let y = [r * theta.cos(), r * theta.sin(), 0.0];
                let u = kelvin_displacement_3d(&origin, &y, i, j, SHEAR, NU);
                total += 0.5 * width * weight * r * u * (2.0 * PI / angles as f64);
            }
        }
    }
    total
}

/// The self block of `G` is the exact disk integral of the Kelvin kernel:
/// `R(7 − 8ν)/(16G(1−ν))` tangentially and `R(3 − 4ν)/(8G(1−ν))` normally, with
/// no coupling.  The old block was `R/(16G(1−ν))·I`.
#[test]
fn self_influence_block_is_the_disk_integral_of_the_kelvin_kernel() {
    let radius = 0.3;
    let area = PI * radius * radius;
    let mesh = single_element(area, [0.0, 0.0, 1.0]);
    let (h, g) = assemble_elastostatic_bem(&mesh, SHEAR, NU);
    assert_eq!((g.rows, g.cols), (3, 3));
    // A lone element has no off-diagonal influence, so H vanishes.
    assert!(h.data.iter().all(|v| v.abs() < 1e-15), "{:?}", h.data);
    for i in 0..3 {
        for j in 0..3 {
            let numeric = disk_integral_of_u(radius, i, j);
            let got = g.get(i, j);
            let scale = numeric.abs().max(g.get(0, 0));
            assert!(
                (got - numeric).abs() <= 1e-9 * scale,
                "G[{i}][{j}] = {got}, disk integral {numeric}"
            );
        }
    }
    let tangential = radius * (7.0 - 8.0 * NU) / (16.0 * SHEAR * (1.0 - NU));
    let normal = radius * (3.0 - 4.0 * NU) / (8.0 * SHEAR * (1.0 - NU));
    assert!((g.get(0, 0) - tangential).abs() < 1e-14);
    assert!((g.get(1, 1) - tangential).abs() < 1e-14);
    assert!((g.get(2, 2) - normal).abs() < 1e-14);
    // The old block R/(16 G (1-nu)) was 4.6x / 3.6x smaller.
    let old = radius / (16.0 * SHEAR * (1.0 - NU));
    assert!((tangential / old - (7.0 - 8.0 * NU)).abs() < 1e-12);
    assert!((normal / old - 2.0 * (3.0 - 4.0 * NU)).abs() < 1e-12);
}

#[test]
fn self_influence_block_follows_the_element_normal() {
    let radius = 0.2;
    let normal = unit([1.0, 2.0, -2.0]);
    let mesh = single_element(PI * radius * radius, normal);
    let (_, g) = assemble_elastostatic_bem(&mesh, SHEAR, NU);
    let tangential = radius * (7.0 - 8.0 * NU) / (16.0 * SHEAR * (1.0 - NU));
    let normal_value = radius * (3.0 - 4.0 * NU) / (8.0 * SHEAR * (1.0 - NU));
    // Symmetric, n is an eigenvector with the normal value, tangents with the other.
    for i in 0..3 {
        for j in 0..3 {
            assert!((g.get(i, j) - g.get(j, i)).abs() < 1e-15);
        }
    }
    let apply = |v: [f64; 3]| -> [f64; 3] {
        let mut out = [0.0; 3];
        for (i, o) in out.iter_mut().enumerate() {
            *o = (0..3).map(|j| g.get(i, j) * v[j]).sum();
        }
        out
    };
    let gn = apply(normal);
    for axis in 0..3 {
        assert!((gn[axis] - normal_value * normal[axis]).abs() < 1e-14);
    }
    let tangent = unit([2.0, 1.0, 2.0]); // orthogonal to (1, 2, -2)
    let gt = apply(tangent);
    for axis in 0..3 {
        assert!((gt[axis] - tangential * tangent[axis]).abs() < 1e-14);
    }
}

fn translation(n_elements: usize, direction: usize) -> Vec<f64> {
    let mut u = vec![0.0; 3 * n_elements];
    for i in 0..n_elements {
        u[3 * i + direction] = 1.0;
    }
    u
}

fn matvec(m: &oxiphysics_fem::boundary_element::DenseMatrix, x: &[f64]) -> Vec<f64> {
    m.matvec(x)
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// A rigid translation produces no traction: `H·(e_k, e_k, ...) = 0` for every
/// direction.  With only the diagonal of the self block set, the off-diagonal
/// entries `H[3i+a, 3i+b]` (a ≠ b) broke this at every element.
#[test]
fn rigid_translations_are_in_the_null_space_of_h() {
    let mesh = BoundaryMesh::sphere(1.0, 6, 12);
    let n = mesh.num_elements();
    let (h, _) = assemble_elastostatic_bem(&mesh, SHEAR, NU);
    let h_scale = h.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for direction in 0..3 {
        let residual = matvec(&h, &translation(n, direction));
        let worst = residual.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(worst < 1e-12 * h_scale, "direction {direction}: {worst:e}");
    }
    // The diagonal blocks are full 3x3 blocks, not diagonal ones.
    let mut coupled = false;
    for i in 0..n {
        for a in 0..3 {
            for b in 0..3 {
                if a != b && h.get(3 * i + a, 3 * i + b).abs() > 1e-6 * h_scale {
                    coupled = true;
                }
            }
        }
    }
    assert!(
        coupled,
        "H must have off-diagonal terms inside its diagonal blocks"
    );
}

fn elastic_modulus() -> f64 {
    2.0 * SHEAR * (1.0 + NU)
}

/// Relative residual `‖H u − G t‖ / ‖G t‖` of a known elastic field on the unit
/// sphere: displacement `u(x)` and traction `t = σ n` at the element centroids.
fn field_residual(
    mesh: &BoundaryMesh,
    displacement: impl Fn(&[f64; 3]) -> [f64; 3],
    stress: impl Fn(&[f64; 3]) -> [[f64; 3]; 3],
) -> f64 {
    let n = mesh.num_elements();
    let (h, g) = assemble_elastostatic_bem(mesh, SHEAR, NU);
    let mut u = vec![0.0; 3 * n];
    let mut t = vec![0.0; 3 * n];
    for (i, element) in mesh.elements.iter().enumerate() {
        let x = centroid(mesh, i);
        let ux = displacement(&x);
        let sigma = stress(&x);
        for a in 0..3 {
            u[3 * i + a] = ux[a];
            t[3 * i + a] = (0..3).map(|b| sigma[a][b] * element.normal[b]).sum();
        }
    }
    let lhs = matvec(&h, &u);
    let rhs = matvec(&g, &t);
    let diff: Vec<f64> = lhs.iter().zip(&rhs).map(|(a, b)| a - b).collect();
    norm(&diff) / norm(&rhs)
}

/// Hydrostatic compression `σ = −p·I`, `u = −p(1 − 2ν)/E · x`.  The wrong
/// self block left a 55 % residual, the swapped traction arguments 250 %.
#[test]
fn hydrostatic_state_satisfies_the_boundary_integral_equation() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let p = 1.0;
    let coefficient = -p * (1.0 - 2.0 * NU) / elastic_modulus();
    let residual = field_residual(
        &mesh,
        |x| [coefficient * x[0], coefficient * x[1], coefficient * x[2]],
        |_| [[-p, 0.0, 0.0], [0.0, -p, 0.0], [0.0, 0.0, -p]],
    );
    assert!(residual < 0.08, "relative residual {residual}");
}

/// Uniaxial tension `σ_zz = 1`.
#[test]
fn uniaxial_state_satisfies_the_boundary_integral_equation() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let e = elastic_modulus();
    let residual = field_residual(
        &mesh,
        |x| [-NU / e * x[0], -NU / e * x[1], x[2] / e],
        |_| [[0.0; 3], [0.0; 3], [0.0, 0.0, 1.0]],
    );
    assert!(residual < 0.08, "relative residual {residual}");
}

/// Simple shear `σ_xy = 1`: `u = (y, x, 0)/(2G)`.
#[test]
fn shear_state_satisfies_the_boundary_integral_equation() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let residual = field_residual(
        &mesh,
        |x| [x[1] / (2.0 * SHEAR), x[0] / (2.0 * SHEAR), 0.0],
        |_| [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    );
    assert!(residual < 0.08, "relative residual {residual}");
}

/// The residual shrinks (or at least does not grow) when the mesh is refined.
#[test]
fn elastostatic_residual_does_not_grow_under_refinement() {
    let p = 1.0;
    let coefficient = -p * (1.0 - 2.0 * NU) / elastic_modulus();
    let residual = |n_lat: usize, n_lon: usize| {
        field_residual(
            &BoundaryMesh::sphere(1.0, n_lat, n_lon),
            |x| [coefficient * x[0], coefficient * x[1], coefficient * x[2]],
            |_| [[-p, 0.0, 0.0], [0.0, -p, 0.0], [0.0, 0.0, -p]],
        )
    };
    let coarse = residual(8, 16);
    let fine = residual(12, 24);
    assert!(fine <= coarse * 1.05, "coarse {coarse}, fine {fine}");
}

/// A rigid rotation `u = ω × x` is traction free: `H u` must be small compared
/// with the entries it is built from.
#[test]
fn rigid_rotation_is_traction_free() {
    let mesh = BoundaryMesh::sphere(1.0, 8, 16);
    let n = mesh.num_elements();
    let (h, _) = assemble_elastostatic_bem(&mesh, SHEAR, NU);
    let mut u = vec![0.0; 3 * n];
    for i in 0..n {
        let x = centroid(&mesh, i);
        u[3 * i] = -x[1];
        u[3 * i + 1] = x[0];
    }
    let hu = matvec(&h, &u);
    let frobenius = h.data.iter().map(|v| v * v).sum::<f64>().sqrt();
    let ratio = norm(&hu) / (frobenius * norm(&u) / ((3 * n) as f64).sqrt());
    assert!(ratio < 0.06, "‖H u‖ / (‖H‖ ‖u‖ / √N) = {ratio}");
}
