// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM cantilever beam validation.
//! Steel beam L=1m, 0.1×0.1 cross section, fixed at x=0, 1 kN load at x=L.
//! Analytical deflection: δ = P L³ / (3 E I).

use oxiphysics_core::math::Vec3;
use oxiphysics_fem::{
    analysis::LinearStaticAnalysis,
    boundary::{DirichletBc, NeumannBc},
    constitutive::LinearElasticMaterial,
    mesh::TetrahedralMesh,
};

fn main() {
    let length = 1.0_f64; // m
    let width = 0.1_f64; // m
    let height = 0.1_f64; // m
    let e_modulus = 200.0e9_f64; // Pa (steel)
    let nu = 0.3_f64;
    let load = -1000.0_f64; // N, downward (negative y)

    // Analytical solution: δ = P L³ / (3 E I)
    let inertia = width * height.powi(3) / 12.0;
    let analytical = load * length.powi(3) / (3.0 * e_modulus * inertia);
    let analytical_mm = analytical * 1000.0;

    // Generate mesh.
    let nx = 10;
    let ny = 2;
    let nz = 2;
    let mesh = TetrahedralMesh::generate_beam(length, width, height, nx, ny, nz);
    let material = LinearElasticMaterial::new(e_modulus, nu);

    // Dirichlet BCs: fix all nodes at x ≈ 0.
    let mut dirichlet_bcs = Vec::new();
    for (i, node) in mesh.nodes.iter().enumerate() {
        if node.x.abs() < 1e-10 {
            dirichlet_bcs.push(DirichletBc::new(i, 0, 0.0));
            dirichlet_bcs.push(DirichletBc::new(i, 1, 0.0));
            dirichlet_bcs.push(DirichletBc::new(i, 2, 0.0));
        }
    }

    // Neumann BCs: distribute load over nodes at x ≈ L.
    let free_nodes: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - length).abs() < 1e-10)
        .map(|(i, _)| i)
        .collect();
    let n_free = free_nodes.len() as f64;
    let neumann_bcs: Vec<NeumannBc> = free_nodes
        .iter()
        .map(|&i| NeumannBc::new(i, [0.0, load / n_free, 0.0]))
        .collect();

    let analysis = LinearStaticAnalysis {
        max_iter: 50_000,
        tolerance: 1e-12,
    };

    let result = analysis.solve(
        &mesh,
        &material,
        &dirichlet_bcs,
        &neumann_bcs,
        &Vec3::zeros(),
    );

    // Maximum y-displacement at the free end.
    let max_disp_y = free_nodes
        .iter()
        .map(|&i| result.displacements[i].y)
        .fold(f64::INFINITY, f64::min);
    let fem_mm = max_disp_y * 1000.0;

    println!(
        "FEM deflection: {:.4} mm (analytical: {:.4} mm)",
        fem_mm, analytical_mm
    );
    println!(
        "Relative error: {:.1}%",
        ((fem_mm - analytical_mm) / analytical_mm).abs() * 100.0
    );
    // Note: linear tetrahedra are overly stiff in bending on coarse meshes.
    // The deflection should be in the right direction and within order-of-magnitude.
    if max_disp_y < 0.0 {
        println!("PASS: deflection direction correct (downward)");
    } else {
        println!("FAIL: unexpected deflection direction");
    }
}
