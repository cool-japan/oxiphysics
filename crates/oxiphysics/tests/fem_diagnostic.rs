// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Diagnostic tests for the linear-tet FEM cantilever locking / stiffness bug.
// These are all #[ignore]'d and must be run explicitly with
//   cargo test -p oxiphysics --test fem_diagnostic -- --ignored --nocapture
// They do NOT assert pass/fail; they print numbers that distinguish pure
// CST shear locking from solver/mesh bugs.

use oxiphysics_core::math::Vec3;
use oxiphysics_fem::analysis::LinearStaticAnalysis;
use oxiphysics_fem::boundary::{DirichletBc, NeumannBc};
use oxiphysics_fem::constitutive::LinearElasticMaterial;
use oxiphysics_fem::mesh::TetrahedralMesh;

const E_MOD: f64 = 210.0e9;
const NU: f64 = 0.3;

// ---------------------------------------------------------------------------
// Helper: count signed-volume-negative tets.
// ---------------------------------------------------------------------------

fn count_inverted_tets(mesh: &TetrahedralMesh) -> (usize, Vec<usize>) {
    let mut negatives = 0usize;
    let mut first_bad = Vec::new();
    for (idx, elem) in mesh.elements.iter().enumerate() {
        let n = [
            mesh.nodes[elem[0]],
            mesh.nodes[elem[1]],
            mesh.nodes[elem[2]],
            mesh.nodes[elem[3]],
        ];
        let x10 = n[1] - n[0];
        let x20 = n[2] - n[0];
        let x30 = n[3] - n[0];
        let det = x10.x * (x20.y * x30.z - x20.z * x30.y) - x10.y * (x20.x * x30.z - x20.z * x30.x)
            + x10.z * (x20.x * x30.y - x20.y * x30.x);
        let v = det / 6.0;
        if v < 0.0 {
            negatives += 1;
            if first_bad.len() < 6 {
                first_bad.push(idx);
            }
        }
    }
    (negatives, first_bad)
}

// ---------------------------------------------------------------------------
// Diagnostic 1: patch test with prescribed linear displacement field.
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn fem_diagnostic_patch_test() {
    // Unit cube [0,1]^3 meshed with generate_beam.
    // generate_beam(L,W,H,nx,ny,nz) -> beam along x from 0..L.
    // We use L=W=H=1, nx=ny=nz=1 to get one hex split into 5 tets (no interior nodes).
    // To have any interior nodes we need a finer split, so use (2,2,2).
    let mesh = TetrahedralMesh::generate_beam(1.0, 1.0, 1.0, 2, 2, 2);
    let material = LinearElasticMaterial::new(E_MOD, NU);

    // Linear displacement field consistent with uniaxial tension under ν = 0.3:
    // stretch alpha in x, contraction -ν*alpha in y and z.
    let alpha = 1.0e-2_f64;
    let beta = -NU * alpha;
    let gamma = -NU * alpha;

    let field = |p: Vec3| -> [f64; 3] { [alpha * p.x, beta * p.y, gamma * p.z] };

    // Identify boundary nodes: anything on the surface of the unit cube.
    let is_boundary = |p: Vec3| -> bool {
        let eps = 1e-9;
        p.x.abs() < eps
            || (p.x - 1.0).abs() < eps
            || p.y.abs() < eps
            || (p.y - 1.0).abs() < eps
            || p.z.abs() < eps
            || (p.z - 1.0).abs() < eps
    };

    // Prescribe the linear field on every boundary node.
    let mut dirichlet = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if is_boundary(*n) {
            let u = field(*n);
            dirichlet.push(DirichletBc::new(i, 0, u[0]));
            dirichlet.push(DirichletBc::new(i, 1, u[1]));
            dirichlet.push(DirichletBc::new(i, 2, u[2]));
        }
    }

    let analysis = LinearStaticAnalysis {
        max_iter: 200_000,
        tolerance: 1e-14,
    };
    let result = analysis.solve(&mesh, &material, &dirichlet, &[], &Vec3::new(0.0, 0.0, 0.0));

    // Check every interior node matches the field.
    let mut worst = 0.0_f64;
    let mut n_interior = 0usize;
    let mut worst_node = 0usize;
    for (i, n) in mesh.nodes.iter().enumerate() {
        if !is_boundary(*n) {
            n_interior += 1;
            let expected = field(*n);
            let u = result.displacements[i];
            let err = ((u.x - expected[0]).powi(2)
                + (u.y - expected[1]).powi(2)
                + (u.z - expected[2]).powi(2))
            .sqrt();
            if err > worst {
                worst = err;
                worst_node = i;
            }
        }
    }

    let pass = worst < 1e-9 || n_interior == 0;
    println!("\n=== DIAGNOSTIC 1: PATCH TEST ===");
    println!("n_nodes_total = {}", mesh.nodes.len());
    println!("n_interior = {}", n_interior);
    println!(
        "prescribed field: u(x,y,z) = ({:.3e}*x, {:.3e}*y, {:.3e}*z)",
        alpha, beta, gamma
    );
    println!(
        "worst interior-node error = {:.6e} m  (at node {})",
        worst, worst_node
    );
    println!(
        "tolerance = 1e-9 m -> patch_test: {}",
        if pass { "pass" } else { "fail" }
    );
    println!(
        "YAML_OUT: patch_test: {}",
        if pass { "pass" } else { "fail" }
    );
    println!("YAML_OUT: patch_test_max_err: {:.6e}", worst);
}

// ---------------------------------------------------------------------------
// Diagnostic 2: uniaxial tension — CST should be exact in 1D stretch.
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn fem_diagnostic_tension_test() {
    // Beam matches the failing cantilever setup: L=1, cs 0.1x0.1.
    let length = 1.0_f64;
    let width = 0.1_f64;
    let height = 0.1_f64;
    let nx = 10usize;
    let ny = 2usize;
    let nz = 2usize;

    let mesh = TetrahedralMesh::generate_beam(length, width, height, nx, ny, nz);
    let material = LinearElasticMaterial::new(E_MOD, NU);

    // Clamp x=0 face (all DOFs = 0).
    let mut dirichlet = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if n.x.abs() < 1e-10 {
            dirichlet.push(DirichletBc::new(i, 0, 0.0));
            dirichlet.push(DirichletBc::new(i, 1, 0.0));
            dirichlet.push(DirichletBc::new(i, 2, 0.0));
        }
    }

    // Uniform axial traction on x=L face, total Fx = +1000 N evenly distributed.
    let tip_nodes: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - length).abs() < 1e-10)
        .map(|(i, _)| i)
        .collect();
    let total_force = 1000.0_f64;
    let per_node = total_force / tip_nodes.len() as f64;
    let neumann: Vec<NeumannBc> = tip_nodes
        .iter()
        .map(|&i| NeumannBc::new(i, [per_node, 0.0, 0.0]))
        .collect();

    let analysis = LinearStaticAnalysis {
        max_iter: 200_000,
        tolerance: 1e-14,
    };
    let result = analysis.solve(
        &mesh,
        &material,
        &dirichlet,
        &neumann,
        &Vec3::new(0.0, 0.0, 0.0),
    );

    // Analytical: u_x(L) = (F/A) L / E.
    let area = width * height;
    let sigma = total_force / area;
    let analytical = sigma * length / E_MOD;

    let measured = tip_nodes
        .iter()
        .map(|&i| result.displacements[i].x)
        .fold(f64::NEG_INFINITY, f64::max);

    let rel_err = ((measured - analytical) / analytical).abs();
    let pass = rel_err < 0.01;

    println!("\n=== DIAGNOSTIC 2: UNIAXIAL TENSION ===");
    println!("mesh = {}x{}x{} ({} tets)", nx, ny, nz, mesh.elements.len());
    println!("analytical u_x(L) = {:.6e} m", analytical);
    println!("measured max u_x at tip = {:.6e} m", measured);
    println!(
        "rel_err = {:.3}%  -> tension_test: {}",
        rel_err * 100.0,
        if pass { "pass" } else { "fail" }
    );
    println!(
        "YAML_OUT: tension_test: {}",
        if pass { "pass" } else { "fail" }
    );
    println!("YAML_OUT: tension_test_measured: {:.6e}", measured);
    println!("YAML_OUT: tension_test_analytical: {:.6e}", analytical);
    println!("YAML_OUT: tension_test_rel_err: {:.6e}", rel_err);
}

// ---------------------------------------------------------------------------
// Diagnostic 3 + 4: h-convergence scan + inverted-tet counts.
// ---------------------------------------------------------------------------

fn solve_cantilever_tip(nx: usize, ny: usize, nz: usize) -> (f64, usize) {
    let length = 1.0_f64;
    let width = 0.1_f64;
    let height = 0.1_f64;
    let tip_load = -1000.0_f64;

    let mesh = TetrahedralMesh::generate_beam(length, width, height, nx, ny, nz);
    let n_tets = mesh.elements.len();

    let material = LinearElasticMaterial::new(E_MOD, NU);

    let mut dirichlet = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if n.x.abs() < 1e-10 {
            dirichlet.push(DirichletBc::new(i, 0, 0.0));
            dirichlet.push(DirichletBc::new(i, 1, 0.0));
            dirichlet.push(DirichletBc::new(i, 2, 0.0));
        }
    }

    let tip_nodes: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - length).abs() < 1e-10)
        .map(|(i, _)| i)
        .collect();
    let n_tip = tip_nodes.len() as f64;
    let neumann: Vec<NeumannBc> = tip_nodes
        .iter()
        .map(|&i| NeumannBc::new(i, [0.0, tip_load / n_tip, 0.0]))
        .collect();

    let analysis = LinearStaticAnalysis {
        max_iter: 200_000,
        tolerance: 1e-12,
    };
    let result = analysis.solve(
        &mesh,
        &material,
        &dirichlet,
        &neumann,
        &Vec3::new(0.0, 0.0, 0.0),
    );

    let min_y = tip_nodes
        .iter()
        .map(|&i| result.displacements[i].y)
        .fold(f64::INFINITY, f64::min);

    (min_y.abs(), n_tets)
}

#[test]
#[ignore]
fn fem_diagnostic_h_convergence_and_inversions() {
    let length = 1.0_f64;
    let width = 0.1_f64;
    let height = 0.1_f64;

    // I = b h^3 / 12 for the square cross-section.
    let inertia = width * height.powi(3) / 12.0;
    let analytical = 1000.0_f64 * length.powi(3) / (3.0 * E_MOD * inertia);

    let resolutions = [(10, 2, 2), (20, 4, 4), (30, 6, 6), (40, 8, 8)];
    let mut prev_err = f64::INFINITY;
    let mut monotone = true;
    let mut rows = Vec::new();
    let mut invert_rows = Vec::new();

    for &(nx, ny, nz) in &resolutions {
        let mesh = TetrahedralMesh::generate_beam(length, width, height, nx, ny, nz);
        let (inv_count, first_bad) = count_inverted_tets(&mesh);
        let total = mesh.elements.len();
        invert_rows.push((nx, ny, nz, inv_count, total, first_bad.clone()));

        let (tip_defl, n_tets) = solve_cantilever_tip(nx, ny, nz);
        let rel_err = ((tip_defl - analytical) / analytical).abs();
        rows.push((nx, ny, nz, n_tets, tip_defl, rel_err));

        if rel_err >= prev_err {
            monotone = false;
        }
        prev_err = rel_err;
    }

    println!("\n=== DIAGNOSTIC 3: H-CONVERGENCE ===");
    println!("analytical tip deflection = {:.6e} m", analytical);
    for (nx, ny, nz, n_tets, tip_defl, rel_err) in &rows {
        println!(
            "  {}x{}x{:>2} tets={:>6}  tip={:.6e}  rel_err={:.4}% ",
            nx,
            ny,
            nz,
            n_tets,
            tip_defl,
            rel_err * 100.0
        );
    }
    println!("convergence_monotone = {}", monotone);

    println!("\n=== DIAGNOSTIC 4: INVERTED-TET COUNT ===");
    for (nx, ny, nz, inv, total, first) in &invert_rows {
        println!(
            "  {}x{}x{:>2}  inverted = {}/{}  first_bad_ids = {:?}",
            nx, ny, nz, inv, total, first
        );
    }

    // YAML-friendly emission
    println!("\n=== YAML SUMMARY ===");
    for (nx, ny, nz, n_tets, tip_defl, rel_err) in &rows {
        println!(
            "YAML_OUT: h_convergence: {{ n_cells: \"{}x{}x{}\", n_tets: {}, tip_defl: {:.6e}, rel_err: {:.6e} }}",
            nx, ny, nz, n_tets, tip_defl, rel_err
        );
    }
    println!("YAML_OUT: convergence_monotone: {}", monotone);
    for (nx, ny, nz, inv, total, _) in &invert_rows {
        println!(
            "YAML_OUT: inverted_tets: {{ mesh: \"{}x{}x{}\", count: {}, total: {} }}",
            nx, ny, nz, inv, total
        );
    }
}
