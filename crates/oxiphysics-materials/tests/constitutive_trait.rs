// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the unified UMAT-style `ConstitutiveModel` trait:
//! parity of the trait path against each model's concrete arithmetic, virgin
//! response, state-variable counts, and a generic single-step driver.

use oxiphysics_materials::constitutive::{
    ConstitutiveModel, ConstitutiveResponse, J2State, ViscoelasticState,
};
use oxiphysics_materials::elastic::{
    IsotropicElastic, LinearElastic, OrthotropicElastic, TransverselyIsotropicElastic,
};
use oxiphysics_materials::plasticity::{J2ConsistentTangent, J2ReturnMapping};
use oxiphysics_materials::viscoelastic::types::{GeneralizedMaxwell, MaxwellModel};

const EPS: [f64; 6] = [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.001];

fn close(a: f64, b: f64) {
    let tol = a.abs().max(1.0) * 1e-12;
    assert!(
        (a - b).abs() <= tol,
        "parity failure: a={a}, b={b}, |a-b|={}, tol={tol}",
        (a - b).abs()
    );
}

fn assert_stress_parity(trait_stress: &[f64; 6], concrete: &[f64; 6]) {
    for (a, b) in trait_stress.iter().zip(concrete.iter()) {
        close(*a, *b);
    }
}

fn assert_tangent_parity(trait_tan: &[[f64; 6]; 6], concrete: &[[f64; 6]; 6]) {
    for (ra, rb) in trait_tan.iter().zip(concrete.iter()) {
        for (a, b) in ra.iter().zip(rb.iter()) {
            close(*a, *b);
        }
    }
}

fn unflatten(flat: &[f64; 36]) -> [[f64; 6]; 6] {
    let mut out = [[0.0_f64; 6]; 6];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = flat[i * 6 + j];
        }
    }
    out
}

fn mat6_vec6(c: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, row) in c.iter().enumerate() {
        for (j, &cij) in row.iter().enumerate() {
            out[i] += cij * v[j];
        }
    }
    out
}

fn flat36_vec6(c: &[f64; 36], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, o) in out.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (j, &vj) in v.iter().enumerate() {
            acc += c[i * 6 + j] * vj;
        }
        *o = acc;
    }
    out
}

fn cfrp_orthotropic() -> OrthotropicElastic {
    OrthotropicElastic {
        e1: 140e9,
        e2: 10e9,
        e3: 10e9,
        g12: 5e9,
        g23: 3.5e9,
        g13: 5e9,
        nu12: 0.3,
        nu23: 0.4,
        nu13: 0.3,
    }
}

#[test]
fn parity_linear_elastic() {
    let m = LinearElastic::new(200e9, 0.3);
    let c = m.stress_strain_matrix_3d();
    let sigma_concrete = mat6_vec6(&c, &EPS);
    let r = m.stress_update(&EPS, &(), 1e-3);
    assert_stress_parity(&r.stress, &sigma_concrete);
    assert_tangent_parity(&r.tangent, &c);
}

#[test]
fn parity_isotropic_elastic() {
    let m = IsotropicElastic::new(200e9, 0.3);
    let c = LinearElastic::new(200e9, 0.3).stress_strain_matrix_3d();
    let sigma_concrete = mat6_vec6(&c, &EPS);
    let r = m.stress_update(&EPS, &(), 1e-3);
    assert_stress_parity(&r.stress, &sigma_concrete);
    assert_tangent_parity(&r.tangent, &c);
}

#[test]
fn parity_orthotropic_elastic() {
    let m = cfrp_orthotropic();
    let c = unflatten(&m.stiffness_voigt());
    let sigma_concrete = mat6_vec6(&c, &EPS);
    let r = m.stress_update(&EPS, &(), 1e-3);
    assert_stress_parity(&r.stress, &sigma_concrete);
    assert_tangent_parity(&r.tangent, &c);
}

#[test]
fn parity_transversely_isotropic_elastic() {
    let m = TransverselyIsotropicElastic::new(140e9, 10e9, 5e9, 0.3, 0.3);
    let c = unflatten(&m.stiffness_voigt());
    let sigma_concrete = mat6_vec6(&c, &EPS);
    let r = m.stress_update(&EPS, &(), 1e-3);
    assert_stress_parity(&r.stress, &sigma_concrete);
    assert_tangent_parity(&r.tangent, &c);
}

#[test]
fn parity_j2_yielding() {
    let eps_y = [0.01, -0.003, -0.003, 0.0, 0.0, 0.005];
    let rm = J2ReturnMapping::from_young_poisson(200e9, 0.3, 250e6, 2e9);
    let tan = J2ConsistentTangent::from_young_poisson(200e9, 0.3, 2e9);
    let c_e = tan.elastic_stiffness();
    let trial = flat36_vec6(&c_e, &eps_y);
    let (s_concrete, dep_concrete, dg) = rm.return_map(&trial, 0.0);
    let tangent_concrete = unflatten(&tan.compute_consistent_tangent(&trial, dg));

    let r = rm.stress_update(&eps_y, &J2State::default(), 1e-3);
    assert_stress_parity(&r.stress, &s_concrete);
    assert_tangent_parity(&r.tangent, &tangent_concrete);
    close(r.state.equiv_plastic_strain, dg);
    for (a, b) in r.state.plastic_strain.iter().zip(dep_concrete.iter()) {
        close(*a, *b);
    }
    // The yielding branch must actually be exercised.
    assert!(dg > 0.0, "test strain path did not yield (dg={dg})");
}

fn two_branch_maxwell() -> GeneralizedMaxwell {
    let mut gm = GeneralizedMaxwell::new(1.0e9);
    gm.add_element(MaxwellModel::new(2.0e9, 4.0e9), 1.0);
    gm.add_element(MaxwellModel::new(0.5e9, 5.0e8), 0.5);
    gm
}

#[test]
fn parity_generalized_maxwell() {
    let gm = two_branch_maxwell();
    let dt = 1e-3;
    // Advance one step from virgin to obtain a non-trivial history.
    let step1 = gm.stress_update(&EPS, &ViscoelasticState::default(), dt);
    let state1 = step1.state.clone();

    // Second step: trait vs concrete.
    let (s_concrete, tan_concrete, hist_concrete) =
        gm.voigt_stress_update(&EPS, &state1.branch_stress, dt);
    let r = gm.stress_update(&EPS, &state1, dt);
    assert_stress_parity(&r.stress, &s_concrete);
    assert_tangent_parity(&r.tangent, &tan_concrete);
    assert_eq!(r.state.branch_stress.len(), hist_concrete.len());
    for (rb, cb) in r.state.branch_stress.iter().zip(hist_concrete.iter()) {
        for (a, b) in rb.iter().zip(cb.iter()) {
            close(*a, *b);
        }
    }
}

#[test]
fn virgin_step_elastic_zero_stress() {
    let zero = [0.0_f64; 6];
    let lin = LinearElastic::new(200e9, 0.3);
    let iso = IsotropicElastic::new(200e9, 0.3);
    let ortho = cfrp_orthotropic();
    let trans = TransverselyIsotropicElastic::new(140e9, 10e9, 5e9, 0.3, 0.3);
    for s in lin.stress_update(&zero, &(), 1e-3).stress {
        assert!(s.abs() < 1e-12);
    }
    for s in iso.stress_update(&zero, &(), 1e-3).stress {
        assert!(s.abs() < 1e-12);
    }
    for s in ortho.stress_update(&zero, &(), 1e-3).stress {
        assert!(s.abs() < 1e-12);
    }
    for s in trans.stress_update(&zero, &(), 1e-3).stress {
        assert!(s.abs() < 1e-12);
    }
}

#[test]
fn virgin_step_j2_zero_stress() {
    let zero = [0.0_f64; 6];
    let rm = J2ReturnMapping::from_young_poisson(200e9, 0.3, 250e6, 2e9);
    let r = rm.stress_update(&zero, &J2State::default(), 1e-3);
    for s in r.stress {
        assert!(s.abs() < 1e-12);
    }
    assert_eq!(r.state.equiv_plastic_strain, 0.0);
}

#[test]
fn virgin_step_generalized_maxwell_zero_stress() {
    let zero = [0.0_f64; 6];
    let gm = two_branch_maxwell();
    let r = gm.stress_update(&zero, &ViscoelasticState::default(), 1e-3);
    for s in r.stress {
        assert!(s.abs() < 1e-12);
    }
}

#[test]
fn n_state_vars_counts() {
    assert_eq!(LinearElastic::new(200e9, 0.3).n_state_vars(), 0);
    assert_eq!(IsotropicElastic::new(200e9, 0.3).n_state_vars(), 0);
    assert_eq!(cfrp_orthotropic().n_state_vars(), 0);
    assert_eq!(
        TransverselyIsotropicElastic::new(140e9, 10e9, 5e9, 0.3, 0.3).n_state_vars(),
        0
    );
    assert_eq!(
        J2ReturnMapping::from_young_poisson(200e9, 0.3, 250e6, 2e9).n_state_vars(),
        7
    );
    assert_eq!(two_branch_maxwell().n_state_vars(), 12);
}

fn drive_one_step<M: ConstitutiveModel>(
    model: &M,
    state: &M::State,
    strain: &[f64; 6],
) -> ConstitutiveResponse<M::State> {
    model.stress_update(strain, state, 1e-3)
}

#[test]
fn generic_driver_runs_all_models() {
    let lin = LinearElastic::new(200e9, 0.3);
    let iso = IsotropicElastic::new(200e9, 0.3);
    let ortho = cfrp_orthotropic();
    let trans = TransverselyIsotropicElastic::new(140e9, 10e9, 5e9, 0.3, 0.3);
    let rm = J2ReturnMapping::from_young_poisson(200e9, 0.3, 250e6, 2e9);
    let gm = two_branch_maxwell();

    assert!(drive_one_step(&lin, &(), &EPS).stress[0].is_finite());
    assert!(drive_one_step(&iso, &(), &EPS).stress[0].is_finite());
    assert!(drive_one_step(&ortho, &(), &EPS).stress[0].is_finite());
    assert!(drive_one_step(&trans, &(), &EPS).stress[0].is_finite());
    assert!(drive_one_step(&rm, &J2State::default(), &EPS).stress[0].is_finite());
    assert!(drive_one_step(&gm, &ViscoelasticState::default(), &EPS).stress[0].is_finite());
}
