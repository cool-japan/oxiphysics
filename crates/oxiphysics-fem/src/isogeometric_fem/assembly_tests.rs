// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the local sparse IGA assembly (`assembly.rs`).
//!
//! The central cross-check is a brute-force *dense* reference (every element
//! spans every control point, basis functions by Cox-de Boor recursion) that
//! shares no code with the production path.  The reference also exists in a
//! `legacy_inverse` flavour that reproduces the old, wrong physical-gradient
//! formula `∇ₓR = J⁻¹ ∇_ξR` (instead of `J⁻ᵀ ∇_ξR`): it agrees with the correct
//! one only when the Jacobian is diagonal, i.e. on axis-aligned geometries.

use super::*;
use std::f64::consts::FRAC_1_SQRT_2;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Tensor-product surface whose control point `(i, j)` sits at
/// `map(i / (nu - 1), j / (nv - 1))` with weight `weight(i, j)`.
fn mapped_surface(
    nu: usize,
    nv: usize,
    degree: (usize, usize),
    map: impl Fn(f64, f64) -> [f64; 2],
    weight: impl Fn(usize, usize) -> f64,
) -> NurbsSurface {
    let mut cps = Vec::with_capacity(nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let s = i as f64 / (nu - 1) as f64;
            let t = j as f64 / (nv - 1) as f64;
            let xy = map(s, t);
            cps.push(WeightedPoint::new([xy[0], xy[1], 0.0], weight(i, j)));
        }
    }
    NurbsSurface::new(
        BSplineBasis::uniform_open(nu, degree.0),
        BSplineBasis::uniform_open(nv, degree.1),
        cps,
        nu,
        nv,
    )
}

/// `width x height` parallelogram: the top edge is shifted by `shear * height`.
fn sheared_rectangle(
    width: f64,
    height: f64,
    shear: f64,
    nu: usize,
    nv: usize,
    degree: usize,
) -> NurbsSurface {
    mapped_surface(
        nu,
        nv,
        (degree, degree),
        |s, t| [width * s + shear * height * t, height * t],
        |_, _| 1.0,
    )
}

/// Curved, strongly non-affine, genuinely rational surface.
fn distorted_rational_surface() -> NurbsSurface {
    mapped_surface(
        6,
        5,
        (2, 2),
        |s, t| {
            [
                2.0 * s + 0.35 * (3.0 * t).sin() + 0.2 * s * t,
                1.5 * t + 0.25 * (2.5 * s).cos() * t + 0.3 * s * s,
            ]
        },
        |i, j| 1.0 + 0.3 * ((i + 2 * j) as f64).sin(),
    )
}

/// Dense brute-force element stiffness sum (all elements span all control points).
fn dense_reference_stiffness(
    surface: &NurbsSurface,
    material: &IgaMaterial,
    num_gauss: usize,
    legacy_inverse: bool,
) -> Vec<f64> {
    let (nu_cp, nv_cp) = (surface.num_u, surface.num_v);
    let n_total = nu_cp * nv_cp;
    let ndof = 2 * n_total;
    let mut k = vec![0.0; ndof * ndof];
    let d_mat = material.plane_stress_matrix();
    let (gp_ref, gw_ref) = gauss_legendre(num_gauss);
    let ku = &surface.basis_u.knots;
    let kv = &surface.basis_v.knots;
    for iu in 0..ku.len() - 1 {
        if ku[iu + 1] - ku[iu] < 1e-14 {
            continue;
        }
        for iv in 0..kv.len() - 1 {
            if kv[iv + 1] - kv[iv] < 1e-14 {
                continue;
            }
            let (gu, wu) = gauss_map(&gp_ref, &gw_ref, ku[iu], ku[iu + 1]);
            let (gv, wv) = gauss_map(&gp_ref, &gw_ref, kv[iv], kv[iv + 1]);
            for (u, w_u) in gu.iter().zip(wu.iter()) {
                for (v, w_v) in gv.iter().zip(wv.iter()) {
                    let nu: Vec<f64> = (0..nu_cp)
                        .map(|i| surface.basis_u.evaluate(i, *u))
                        .collect();
                    let dnu: Vec<f64> = (0..nu_cp)
                        .map(|i| surface.basis_u.derivative(i, *u))
                        .collect();
                    let nv: Vec<f64> = (0..nv_cp)
                        .map(|j| surface.basis_v.evaluate(j, *v))
                        .collect();
                    let dnv: Vec<f64> = (0..nv_cp)
                        .map(|j| surface.basis_v.derivative(j, *v))
                        .collect();
                    let mut nw = vec![0.0; n_total];
                    let mut dnw_du = vec![0.0; n_total];
                    let mut dnw_dv = vec![0.0; n_total];
                    let (mut w_sum, mut dw_du, mut dw_dv) = (0.0, 0.0, 0.0);
                    for j in 0..nv_cp {
                        for i in 0..nu_cp {
                            let idx = j * nu_cp + i;
                            let wi = surface.control_points[idx].weight;
                            nw[idx] = nu[i] * nv[j] * wi;
                            dnw_du[idx] = dnu[i] * nv[j] * wi;
                            dnw_dv[idx] = nu[i] * dnv[j] * wi;
                            w_sum += nw[idx];
                            dw_du += dnw_du[idx];
                            dw_dv += dnw_dv[idx];
                        }
                    }
                    let dxu = surface.deriv_u(*u, *v);
                    let dxv = surface.deriv_v(*u, *v);
                    let jac = [[dxu[0], dxv[0]], [dxu[1], dxv[1]]];
                    let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
                    let inv = [
                        [jac[1][1] / det, -jac[0][1] / det],
                        [-jac[1][0] / det, jac[0][0] / det],
                    ];
                    let mut b = vec![0.0; 3 * ndof];
                    for idx in 0..n_total {
                        let dr_du = (dnw_du[idx] * w_sum - nw[idx] * dw_du) / (w_sum * w_sum);
                        let dr_dv = (dnw_dv[idx] * w_sum - nw[idx] * dw_dv) / (w_sum * w_sum);
                        let (dx, dy) = if legacy_inverse {
                            // the old, wrong formula: rows of J^-1
                            (
                                inv[0][0] * dr_du + inv[0][1] * dr_dv,
                                inv[1][0] * dr_du + inv[1][1] * dr_dv,
                            )
                        } else {
                            // J^-T: columns of J^-1
                            (
                                inv[0][0] * dr_du + inv[1][0] * dr_dv,
                                inv[0][1] * dr_du + inv[1][1] * dr_dv,
                            )
                        };
                        b[2 * idx] = dx;
                        b[ndof + 2 * idx + 1] = dy;
                        b[2 * ndof + 2 * idx] = dy;
                        b[2 * ndof + 2 * idx + 1] = dx;
                    }
                    let weight = w_u * w_v * det.abs() * material.thickness;
                    for r in 0..ndof {
                        for c in 0..ndof {
                            let mut val = 0.0;
                            for s in 0..3 {
                                for t in 0..3 {
                                    val += b[s * ndof + r] * d_mat[s][t] * b[t * ndof + c];
                                }
                            }
                            k[r * ndof + c] += val * weight;
                        }
                    }
                }
            }
        }
    }
    k
}

/// `max |a - b| / max |b|`.
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0_f64, |m, x| m.max(x.abs())).max(1e-300);
    a.iter()
        .zip(b.iter())
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
        / scale
}

fn material() -> IgaMaterial {
    IgaMaterial::new(2.0e5, 0.3, 0.05)
}

fn assemble_dense(surface: &NurbsSurface, num_gauss: usize) -> Vec<f64> {
    let assembler = IgaAssembler::new(surface.num_u * surface.num_v, material(), num_gauss);
    assembler
        .assemble(surface)
        .expect("valid surface")
        .dense_stiffness()
}

/// `y = K x` for a dense row-major matrix.
fn dense_matvec(k: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| (0..n).map(|j| k[i * n + j] * x[j]).sum())
        .collect()
}

/// Control-point displacements `u_a = (g(x_a), h(x_a))` for a displacement field.
fn nodal_field(surface: &NurbsSurface, field: impl Fn(f64, f64) -> [f64; 2]) -> Vec<f64> {
    let mut u = Vec::with_capacity(2 * surface.control_points.len());
    for cp in &surface.control_points {
        let d = field(cp.pos[0], cp.pos[1]);
        u.push(d[0]);
        u.push(d[1]);
    }
    u
}

// ---------------------------------------------------------------------------
// Local basis
// ---------------------------------------------------------------------------

#[test]
fn elements_carry_local_connectivity() {
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 6, 5, 2);
    let assembler = IgaAssembler::new(30, material(), 3);
    let elements = assembler.create_elements(&surf);
    // uniform_open(6, 2) has 4 spans, uniform_open(5, 2) has 3 spans.
    assert_eq!(elements.len(), 4 * 3);
    let mut covered = [false; 30];
    for elem in &elements {
        assert_eq!(
            elem.basis_indices.len(),
            9,
            "(p + 1)(q + 1) active functions"
        );
        let mut sorted = elem.basis_indices.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 9, "no duplicate indices");
        assert!(sorted.iter().all(|&g| g < 30));
        for &g in &sorted {
            covered[g] = true;
        }
        // the 3 x 3 window is contiguous in the control net
        let (i0, j0) = (sorted[0] % 6, sorted[0] / 6);
        for b in 0..3 {
            for a in 0..3 {
                assert_eq!(elem.basis_indices[b * 3 + a], (j0 + b) * 6 + (i0 + a));
            }
        }
    }
    assert!(
        covered.iter().all(|&c| c),
        "every control point is active somewhere"
    );
}

#[test]
fn local_basis_is_partition_of_unity_with_vanishing_gradient_sum() {
    let surf = distorted_rational_surface();
    let assembler = IgaAssembler::new(30, material(), 3);
    for elem in assembler.create_elements(&surf) {
        for fu in [0.0, 0.21, 0.5, 0.93, 1.0] {
            for fv in [0.0, 0.37, 0.5, 0.8, 1.0] {
                let u = elem.param_domain[0] + fu * (elem.param_domain[1] - elem.param_domain[0]);
                let v = elem.param_domain[2] + fv * (elem.param_domain[3] - elem.param_domain[2]);
                let basis = elem.local_basis(&surf, u, v).expect("valid element");
                let sum: f64 = basis.values.iter().sum();
                let sum_du: f64 = basis.du.iter().sum();
                let sum_dv: f64 = basis.dv.iter().sum();
                assert!((sum - 1.0).abs() < 1e-12, "sum R = {sum}");
                assert!(sum_du.abs() < 1e-10, "sum dR/du = {sum_du}");
                assert!(sum_dv.abs() < 1e-10, "sum dR/dv = {sum_dv}");
            }
        }
    }
}

#[test]
fn local_basis_matches_dense_cox_de_boor() {
    let surf = distorted_rational_surface();
    let (nu_cp, nv_cp) = (surf.num_u, surf.num_v);
    let assembler = IgaAssembler::new(nu_cp * nv_cp, material(), 3);
    for elem in assembler.create_elements(&surf) {
        for fu in [0.11, 0.5, 0.83] {
            for fv in [0.2, 0.5, 0.97] {
                let u = elem.param_domain[0] + fu * (elem.param_domain[1] - elem.param_domain[0]);
                let v = elem.param_domain[2] + fv * (elem.param_domain[3] - elem.param_domain[2]);
                let local = elem.local_basis(&surf, u, v).unwrap();

                let nu: Vec<f64> = (0..nu_cp).map(|i| surf.basis_u.evaluate(i, u)).collect();
                let dnu: Vec<f64> = (0..nu_cp).map(|i| surf.basis_u.derivative(i, u)).collect();
                let nv: Vec<f64> = (0..nv_cp).map(|j| surf.basis_v.evaluate(j, v)).collect();
                let dnv: Vec<f64> = (0..nv_cp).map(|j| surf.basis_v.derivative(j, v)).collect();
                let (mut w, mut dw_du, mut dw_dv) = (0.0, 0.0, 0.0);
                for j in 0..nv_cp {
                    for i in 0..nu_cp {
                        let wi = surf.control_points[j * nu_cp + i].weight;
                        w += wi * nu[i] * nv[j];
                        dw_du += wi * dnu[i] * nv[j];
                        dw_dv += wi * nu[i] * dnv[j];
                    }
                }
                let mut covered = 0.0;
                for (loc, &g) in elem.basis_indices.iter().enumerate() {
                    let (i, j) = (g % nu_cp, g / nu_cp);
                    let wi = surf.control_points[g].weight;
                    let r = wi * nu[i] * nv[j] / w;
                    let dr_du = wi * (dnu[i] * nv[j] * w - nu[i] * nv[j] * dw_du) / (w * w);
                    let dr_dv = wi * (nu[i] * dnv[j] * w - nu[i] * nv[j] * dw_dv) / (w * w);
                    assert!((local.values[loc] - r).abs() < 1e-12, "R mismatch");
                    assert!((local.du[loc] - dr_du).abs() < 1e-9, "dR/du mismatch");
                    assert!((local.dv[loc] - dr_dv).abs() < 1e-9, "dR/dv mismatch");
                    covered += r;
                }
                // every other control point contributes nothing on this element
                assert!(
                    (covered - 1.0).abs() < 1e-12,
                    "inactive functions are not zero"
                );
            }
        }
    }
}

#[test]
fn inconsistent_elements_are_rejected() {
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 4, 4, 2);
    let assembler = IgaAssembler::new(16, material(), 2);
    let good = assembler.create_elements(&surf).remove(0);
    assert!(good.local_basis(&surf, 0.1, 0.1).is_ok());

    // empty connectivity (what the old dense code accepted)
    let empty = IgaElement::new(vec![], good.param_domain, 2);
    assert!(empty.local_basis(&surf, 0.1, 0.1).is_err());
    assert!(empty.stiffness_matrix(&surf, &material()).is_err());

    // right length but the wrong control points for the knot span
    let mut wrong = good.basis_indices.clone();
    wrong.rotate_left(1);
    let wrong = IgaElement::new(wrong, good.param_domain, 2);
    assert!(wrong.local_basis(&surf, 0.1, 0.1).is_err());

    // parametric domain outside the knot range
    let outside = IgaElement::new(good.basis_indices.clone(), [2.0, 3.0, 0.0, 0.5], 2);
    assert!(outside.local_basis(&surf, 2.5, 0.1).is_err());
}

// ---------------------------------------------------------------------------
// Stiffness: equivalence with the dense reference, J^-T
// ---------------------------------------------------------------------------

#[test]
fn matches_dense_reference_on_axis_aligned_4x4_rectangle() {
    // On an axis-aligned rectangle J is diagonal, so J^-1 == J^-T and the old
    // (wrong) and the correct physical-gradient formulas coincide.
    let surf = NurbsSurface::flat_rectangle(2.0, 1.0, 4, 4, 2);
    let k = assemble_dense(&surf, 3);
    let reference = dense_reference_stiffness(&surf, &material(), 3, false);
    let legacy = dense_reference_stiffness(&surf, &material(), 3, true);
    assert!(
        relative_difference(&k, &reference) < 1e-10,
        "local vs dense: {}",
        relative_difference(&k, &reference)
    );
    assert!(
        relative_difference(&legacy, &reference) < 1e-10,
        "J^-1 == J^-T for a diagonal Jacobian"
    );
}

#[test]
fn matches_correct_dense_reference_on_sheared_geometry_but_not_the_legacy_one() {
    let surf = sheared_rectangle(2.0, 1.0, 0.6, 5, 4, 2);
    let k = assemble_dense(&surf, 3);
    let reference = dense_reference_stiffness(&surf, &material(), 3, false);
    let legacy = dense_reference_stiffness(&surf, &material(), 3, true);
    assert!(
        relative_difference(&k, &reference) < 1e-10,
        "local vs dense J^-T: {}",
        relative_difference(&k, &reference)
    );
    assert!(
        relative_difference(&legacy, &reference) > 1e-2,
        "the transposed-inverse formula must differ once J is not diagonal"
    );
    assert!(relative_difference(&k, &legacy) > 1e-2);
}

#[test]
fn matches_dense_reference_on_distorted_rational_surface() {
    let surf = distorted_rational_surface();
    let k = assemble_dense(&surf, 4);
    let reference = dense_reference_stiffness(&surf, &material(), 4, false);
    assert!(
        relative_difference(&k, &reference) < 1e-9,
        "local vs dense: {}",
        relative_difference(&k, &reference)
    );
}

#[test]
fn assembled_stiffness_is_symmetric() {
    let surf = distorted_rational_surface();
    let k = assemble_dense(&surf, 3);
    let n = 2 * surf.control_points.len();
    let scale = k.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    for i in 0..n {
        for j in 0..n {
            assert!((k[i * n + j] - k[j * n + i]).abs() < 1e-10 * scale);
        }
    }
}

// ---------------------------------------------------------------------------
// Stiffness: physics (rigid-body modes, uniform strain)
// ---------------------------------------------------------------------------

#[test]
fn rigid_body_modes_are_in_the_null_space() {
    // Needs J^-T: with J^-1 a rotation produces spurious strain as soon as the
    // Jacobian is not diagonal.
    for surf in [
        sheared_rectangle(2.0, 1.0, 0.6, 5, 4, 2),
        distorted_rational_surface(),
    ] {
        let k = assemble_dense(&surf, 4);
        let scale = k.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let translation_x = nodal_field(&surf, |_, _| [1.0, 0.0]);
        let translation_y = nodal_field(&surf, |_, _| [0.0, 1.0]);
        let rotation = nodal_field(&surf, |x, y| [-y, x]);
        for (name, mode) in [
            ("x translation", translation_x),
            ("y translation", translation_y),
            ("rotation", rotation),
        ] {
            let force = dense_matvec(&k, &mode);
            let worst = force.iter().fold(0.0_f64, |m, f| m.max(f.abs()));
            assert!(
                worst < 1e-9 * scale,
                "{name}: K u = {worst:e} (|K| = {scale:e})"
            );
        }
    }
}

#[test]
fn uniform_strain_energy_is_exact_on_sheared_rectangle() {
    let (width, height, shear) = (2.0, 1.5, 0.7);
    let surf = sheared_rectangle(width, height, shear, 6, 5, 2);
    let mat = material();
    let k = assemble_dense(&surf, 3);

    // u = G x with a general (non-symmetric) displacement gradient.
    let g = [[0.012, -0.007], [0.004, 0.009]];
    let u = nodal_field(&surf, |x, y| {
        [g[0][0] * x + g[0][1] * y, g[1][0] * x + g[1][1] * y]
    });
    let ku = dense_matvec(&k, &u);
    let energy2: f64 = u.iter().zip(ku.iter()).map(|(a, b)| a * b).sum();

    let eps = [g[0][0], g[1][1], g[0][1] + g[1][0]];
    let d = mat.plane_stress_matrix();
    let density: f64 = d
        .iter()
        .zip(eps.iter())
        .map(|(row, e_r)| {
            e_r * row
                .iter()
                .zip(eps.iter())
                .map(|(d_rc, e_c)| d_rc * e_c)
                .sum::<f64>()
        })
        .sum();
    let expected = density * (width * height) * mat.thickness;
    assert!(
        (energy2 - expected).abs() < 1e-9 * expected,
        "u^T K u = {energy2}, expected {expected}"
    );
}

#[test]
fn body_force_resultant_equals_force_times_area_times_thickness() {
    let (width, height, shear) = (2.0, 1.5, 0.7);
    let surf = sheared_rectangle(width, height, shear, 6, 5, 2);
    let mat = material();
    let assembler = IgaAssembler::new(30, mat, 3);
    let mut system = assembler.assemble(&surf).unwrap();
    assembler
        .apply_body_force(&surf, &mut system, [3.0, -1.5])
        .unwrap();
    let (mut fx, mut fy) = (0.0, 0.0);
    for pair in system.load.chunks(2) {
        fx += pair[0];
        fy += pair[1];
    }
    let volume = width * height * mat.thickness;
    assert!((fx - 3.0 * volume).abs() < 1e-12 * volume, "Fx = {fx}");
    assert!((fy + 1.5 * volume).abs() < 1e-12 * volume, "Fy = {fy}");
}

#[test]
fn body_force_on_rational_surface_uses_local_shape_functions() {
    let surf = distorted_rational_surface();
    let mat = material();
    let assembler = IgaAssembler::new(30, mat, 4);
    let mut system = assembler.assemble(&surf).unwrap();
    assembler
        .apply_body_force(&surf, &mut system, [1.0, 0.0])
        .unwrap();
    // Resultant = integral of the dense Jacobian determinant over the surface.
    let (gp, gw) = gauss_legendre(4);
    let mut area = 0.0;
    for elem in assembler.create_elements(&surf) {
        let (gu, wu) = gauss_map(&gp, &gw, elem.param_domain[0], elem.param_domain[1]);
        let (gv, wv) = gauss_map(&gp, &gw, elem.param_domain[2], elem.param_domain[3]);
        for (u, w_u) in gu.iter().zip(wu.iter()) {
            for (v, w_v) in gv.iter().zip(wv.iter()) {
                area += w_u * w_v * elem.jacobian_det(&surf, *u, *v).abs();
            }
        }
    }
    let fx: f64 = system.load.chunks(2).map(|c| c[0]).sum();
    let fy: f64 = system.load.chunks(2).map(|c| c[1]).sum();
    assert!(
        (fx - area * mat.thickness).abs() < 1e-12 * area,
        "Fx = {fx}"
    );
    assert!(fy.abs() < 1e-14);
}

// ---------------------------------------------------------------------------
// Sparse structure, scalability, validation
// ---------------------------------------------------------------------------

#[test]
fn assembly_is_genuinely_sparse() {
    let n = 20;
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, n, n, 2);
    let assembler = IgaAssembler::new(n * n, material(), 3);
    let system = assembler.assemble(&surf).unwrap();
    let ndof = 2 * n * n;
    assert_eq!(system.ndof, ndof);
    assert_eq!(system.stiffness.nrows, ndof);
    // A control point couples with at most (2p + 1)^2 = 25 control points, i.e.
    // 50 columns per DOF row.
    let max_row = (0..ndof)
        .map(|r| system.stiffness.row_ptr[r + 1] - system.stiffness.row_ptr[r])
        .max()
        .unwrap();
    assert!(max_row <= 50, "row with {max_row} entries");
    assert!(system.stiffness.nnz() <= 50 * ndof);
    assert!(
        (system.stiffness.nnz() as f64) < 0.1 * (ndof * ndof) as f64,
        "nnz = {} of {}",
        system.stiffness.nnz(),
        ndof * ndof
    );
}

#[test]
fn large_net_assembles_without_dense_storage() {
    // 40 x 40 control points = 3200 DOFs.  The former implementation allocated
    // a dense 3200 x 3200 matrix (82 MB) for *every* one of the 38 x 38
    // elements and evaluated every basis function by recursion.
    let n = 40;
    let surf = NurbsSurface::flat_rectangle(4.0, 4.0, n, n, 2);
    let assembler = IgaAssembler::new(n * n, material(), 3);
    let system = assembler.assemble(&surf).unwrap();
    assert_eq!(system.ndof, 2 * n * n);
    assert!(system.stiffness.nnz() < 50 * system.ndof);
    // K u_rigid = 0 for a translation
    let ux = vec![1.0; system.ndof / 2]
        .into_iter()
        .flat_map(|v| [v, 0.0])
        .collect::<Vec<f64>>();
    let f = system.stiffness.mul_vec(&ux);
    let scale = system
        .stiffness
        .values
        .iter()
        .fold(0.0_f64, |m, x| m.max(x.abs()));
    assert!(f.iter().all(|x| x.abs() < 1e-9 * scale));
}

#[test]
fn assemble_rejects_wrong_control_point_count() {
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 4, 4, 2);
    // too few: the old code panicked with an out-of-bounds index here
    let too_few = IgaAssembler::new(10, material(), 2);
    let err = too_few.assemble(&surf).unwrap_err();
    assert!(matches!(err, crate::Error::InvalidInput(_)), "{err}");
    // too many: the old code silently produced a larger, partly empty matrix
    let too_many = IgaAssembler::new(20, material(), 2);
    assert!(too_many.assemble(&surf).is_err());
    let ok = IgaAssembler::new(16, material(), 2);
    assert!(ok.assemble(&surf).is_ok());

    let mut system = ok.assemble(&surf).unwrap();
    assert!(
        too_few
            .apply_body_force(&surf, &mut system, [1.0, 0.0])
            .is_err()
    );
    system.load.pop();
    assert!(ok.apply_body_force(&surf, &mut system, [1.0, 0.0]).is_err());
}

// ---------------------------------------------------------------------------
// End to end: patch test and solve
// ---------------------------------------------------------------------------

/// Control points on the outer ring of the net.
fn boundary_control_points(surf: &NurbsSurface) -> Vec<usize> {
    let (nu, nv) = (surf.num_u, surf.num_v);
    (0..nu * nv)
        .filter(|&g| {
            let (i, j) = (g % nu, g / nu);
            i == 0 || j == 0 || i == nu - 1 || j == nv - 1
        })
        .collect()
}

#[test]
fn patch_test_recovers_linear_displacement_on_sheared_geometry() {
    // Constant-strain patch test: prescribe the exact linear displacement on the
    // boundary control points; the interior must reproduce it.  A wrong
    // physical gradient (J^-1) breaks this as soon as the mesh is sheared.
    let surf = sheared_rectangle(2.0, 1.5, 0.7, 7, 6, 2);
    let n = surf.control_points.len();
    let assembler = IgaAssembler::new(n, material(), 3);
    let mut system = assembler.assemble(&surf).unwrap();

    let exact = |x: f64, y: f64| {
        [
            0.010 * x + 0.004 * y + 0.002,
            -0.003 * x + 0.007 * y - 0.001,
        ]
    };
    let u_exact = nodal_field(&surf, exact);
    let mut dofs = Vec::new();
    let mut values = Vec::new();
    for g in boundary_control_points(&surf) {
        for c in 0..2 {
            dofs.push(2 * g + c);
            values.push(u_exact[2 * g + c]);
        }
    }
    // The penalty error is (boundary reaction) / penalty ~ 1e-9 and the CG error
    // (relative residual 1e-14 of a right-hand side dominated by the penalty
    // terms) ~ 1e-8; both are far below the O(1e-2) error that a wrong physical
    // gradient produces on this sheared mesh.
    let penalty = 1e10;
    IgaBoundaryCondition::new(dofs, values, IgaBcType::Dirichlet)
        .apply(&mut system, penalty)
        .unwrap();

    let u = system
        .solve(20_000, 1e-14)
        .expect("patch problem converges");
    let scale = u_exact.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let worst = u
        .iter()
        .zip(u_exact.iter())
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(
        worst < 1e-5 * scale,
        "max error {worst:e} (scale {scale:e})"
    );
}

#[test]
fn patch_test_on_curved_rational_surface() {
    // A linear field is reproduced exactly even by a curved rational geometry
    // (sum R_a = 1, sum R_a x_a = x), so the patch test must also pass there.
    let surf = distorted_rational_surface();
    let n = surf.control_points.len();
    let assembler = IgaAssembler::new(n, material(), 4);
    let mut system = assembler.assemble(&surf).unwrap();
    let exact = |x: f64, y: f64| [0.01 * x - 0.004 * y, 0.006 * x + 0.002 * y + 0.001];
    let u_exact = nodal_field(&surf, exact);
    let mut dofs = Vec::new();
    let mut values = Vec::new();
    for g in boundary_control_points(&surf) {
        for c in 0..2 {
            dofs.push(2 * g + c);
            values.push(u_exact[2 * g + c]);
        }
    }
    IgaBoundaryCondition::new(dofs, values, IgaBcType::Dirichlet)
        .apply(&mut system, 1e10)
        .unwrap();
    let u = system.solve(20_000, 1e-14).expect("converges");
    let scale = u_exact.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let worst = u
        .iter()
        .zip(u_exact.iter())
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    // the interior stiffness rows only see strains from the quadrature of
    // rational functions, so allow a modest quadrature error
    assert!(
        worst < 5e-4 * scale,
        "max error {worst:e} (scale {scale:e})"
    );
}

#[test]
fn solve_reports_non_convergence_for_an_unconstrained_structure() {
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 4, 4, 2);
    let assembler = IgaAssembler::new(16, material(), 3);
    let mut system = assembler.assemble(&surf).unwrap();
    assembler
        .apply_body_force(&surf, &mut system, [1.0, 0.0])
        .unwrap();
    // free-floating body under a net force: no equilibrium exists
    let err = system.solve(200, 1e-10).unwrap_err();
    assert!(
        matches!(
            err,
            crate::Error::NotConverged { .. } | crate::Error::SingularSystem
        ),
        "{err}"
    );
}

#[test]
fn solve_validates_system_dimensions() {
    let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 3, 3, 1);
    let assembler = IgaAssembler::new(9, material(), 2);
    let mut system = assembler.assemble(&surf).unwrap();
    system.load.pop();
    assert!(matches!(
        system.solve(10, 1e-8),
        Err(crate::Error::InvalidInput(_))
    ));
}

#[test]
fn quarter_annulus_with_rational_arc_has_expected_area() {
    // Exact rational quarter annulus (quadratic arc x linear radial): the
    // integral of |det J| over the single element is pi (r_out^2 - r_in^2) / 4.
    let (r_in, r_out) = (1.0, 2.0);
    let arc = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let arc_w = [1.0, FRAC_1_SQRT_2, 1.0];
    let mut cps = Vec::new();
    for r in [r_in, r_out] {
        for (a, w) in arc.iter().zip(arc_w.iter()) {
            cps.push(WeightedPoint::new([r * a[0], r * a[1], 0.0], *w));
        }
    }
    let surf = NurbsSurface::new(
        BSplineBasis::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2),
        BSplineBasis::new(vec![0.0, 0.0, 1.0, 1.0], 1),
        cps,
        3,
        2,
    );
    let assembler = IgaAssembler::new(6, material(), 5);
    let elements = assembler.create_elements(&surf);
    assert_eq!(elements.len(), 1);
    let elem = &elements[0];
    let (gp, gw) = gauss_legendre(5);
    let (gu, wu) = gauss_map(&gp, &gw, 0.0, 1.0);
    let (gv, wv) = gauss_map(&gp, &gw, 0.0, 1.0);
    let mut area_local = 0.0;
    for (u, w_u) in gu.iter().zip(wu.iter()) {
        for (v, w_v) in gv.iter().zip(wv.iter()) {
            let basis = elem.local_basis(&surf, *u, *v).unwrap();
            let jac = elem.local_jacobian(&surf, &basis);
            area_local += w_u * w_v * (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]).abs();
        }
    }
    let expected = std::f64::consts::PI * (r_out * r_out - r_in * r_in) / 4.0;
    assert!(
        (area_local - expected).abs() < 1e-6 * expected,
        "area {area_local}, expected {expected}"
    );
}

// ---------------------------------------------------------------------------
// B-matrix: physical gradients use J^-T
// ---------------------------------------------------------------------------

/// Strain `(eps_xx, eps_yy, gamma_xy) = B u_e` of the nodal displacement field
/// `nodal` (layout `2 * control_point + component`) at `(u, v)` of `elem`.
fn element_strain(
    surface: &NurbsSurface,
    elem: &IgaElement,
    u: f64,
    v: f64,
    nodal: &[f64],
) -> ([f64; 3], f64) {
    let (b, det_j) = elem.b_matrix(surface, u, v).expect("valid element");
    let ndof = 2 * elem.basis_indices.len();
    assert_eq!(b.len(), 3 * ndof);
    let mut local = Vec::with_capacity(ndof);
    for &g in &elem.basis_indices {
        local.push(nodal[2 * g]);
        local.push(nodal[2 * g + 1]);
    }
    let mut strain = [0.0; 3];
    for (row, e) in strain.iter_mut().enumerate() {
        *e = (0..ndof).map(|c| b[row * ndof + c] * local[c]).sum();
    }
    (strain, det_j)
}

/// Displacement fields `u = G x` with their exact constant strains.
fn linear_fields(surface: &NurbsSurface) -> Vec<(&'static str, Vec<f64>, [f64; 3])> {
    let affine = |g: [[f64; 2]; 2]| {
        nodal_field(surface, move |x, y| {
            [g[0][0] * x + g[0][1] * y, g[1][0] * x + g[1][1] * y]
        })
    };
    vec![
        // u_x = x: eps = [1, 0, 0]
        ("stretch", affine([[1.0, 0.0], [0.0, 0.0]]), [1.0, 0.0, 0.0]),
        // u_x = y: gamma_xy = 1
        ("shear", affine([[0.0, 1.0], [0.0, 0.0]]), [0.0, 0.0, 1.0]),
        // rigid rotation u = (-y, x): no strain at all
        ("rotation", affine([[0.0, -1.0], [1.0, 0.0]]), [0.0; 3]),
        // general affine field
        (
            "affine",
            affine([[0.3, -0.2], [0.1, 0.4]]),
            [0.3, 0.4, -0.2 + 0.1],
        ),
    ]
}

/// Parallelogram `x = s + 0.5 t, y = t`: `x_v = 0.5 != y_u = 0`, so `J^-1` and
/// `J^-T` differ.  The stretch `u_x = x` must give `eps = [1, 0, 0]` and a rigid
/// rotation no strain at all; with the former `J^-1` formula both are wrong.
#[test]
fn b_matrix_reproduces_constant_strain_on_a_sheared_parallelogram() {
    let surf = mapped_surface(5, 4, (2, 2), |s, t| [s + 0.5 * t, t], |_, _| 1.0);
    let assembler = IgaAssembler::new(20, material(), 3);
    let fields = linear_fields(&surf);
    let (gp, gw) = gauss_legendre(3);
    let mut area = 0.0;
    let mut points = 0;
    for elem in assembler.create_elements(&surf) {
        let (gu, wu) = gauss_map(&gp, &gw, elem.param_domain[0], elem.param_domain[1]);
        let (gv, wv) = gauss_map(&gp, &gw, elem.param_domain[2], elem.param_domain[3]);
        for (u, w_u) in gu.iter().zip(wu.iter()) {
            for (v, w_v) in gv.iter().zip(wv.iter()) {
                let mut det_at_point = 0.0;
                for (name, nodal, exact) in &fields {
                    let (strain, det_j) = element_strain(&surf, &elem, *u, *v, nodal);
                    assert!(det_j > 0.0, "det J = {det_j}");
                    det_at_point = det_j;
                    for (e, x) in strain.iter().zip(exact.iter()) {
                        assert!(
                            (e - x).abs() < 1e-10,
                            "{name}: strain {strain:?}, expected {exact:?} at ({u}, {v})"
                        );
                    }
                }
                area += w_u * w_v * det_at_point;
                points += 1;
            }
        }
    }
    // uniform_open(5, 2) has 3 knot spans, uniform_open(4, 2) has 2; 3 x 3 Gauss points
    assert_eq!(points, 3 * 2 * 9);
    // shearing a unit square does not change its area
    assert!((area - 1.0).abs() < 1e-12, "area {area}");
}

#[test]
fn b_matrix_reproduces_constant_strain_on_a_distorted_rational_surface() {
    let surf = distorted_rational_surface();
    let assembler = IgaAssembler::new(30, material(), 4);
    let fields = linear_fields(&surf);
    let (gp, gw) = gauss_legendre(4);
    for elem in assembler.create_elements(&surf) {
        let (gu, _) = gauss_map(&gp, &gw, elem.param_domain[0], elem.param_domain[1]);
        let (gv, _) = gauss_map(&gp, &gw, elem.param_domain[2], elem.param_domain[3]);
        for u in &gu {
            for v in &gv {
                for (name, nodal, exact) in &fields {
                    let (strain, _) = element_strain(&surf, &elem, *u, *v, nodal);
                    for (e, x) in strain.iter().zip(exact.iter()) {
                        assert!(
                            (e - x).abs() < 1e-9,
                            "{name}: strain {strain:?}, expected {exact:?} at ({u}, {v})"
                        );
                    }
                }
            }
        }
    }
}

/// The B-matrix of a degenerate (collapsed) Jacobian is zero instead of NaN.
#[test]
fn b_matrix_of_a_collapsed_element_is_zero() {
    // all control points on one line: det J = 0 everywhere
    let surf = mapped_surface(4, 4, (2, 2), |s, t| [s + t, s + t], |_, _| 1.0);
    let assembler = IgaAssembler::new(16, material(), 2);
    let elem = assembler.create_elements(&surf).remove(0);
    let (b, det_j) = elem.b_matrix(&surf, 0.2, 0.3).expect("valid element");
    assert!(det_j.abs() < 1e-12, "det J = {det_j}");
    assert_eq!(b.len(), 3 * 2 * 9);
    assert!(b.iter().all(|x| *x == 0.0));
    let k = elem.stiffness_matrix(&surf, &material()).expect("valid");
    assert!(k.iter().all(|x| *x == 0.0));
}
