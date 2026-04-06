// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cloth simulation demo using XPBD position-based dynamics.
//!
//! Creates a 10x10 cloth grid, fixes the top row, applies gravity,
//! and runs the PBD solver for 300 iterations. Prints the final
//! particle positions showing the draping shape.

use oxiphysics_core::math::Vec3;
use oxiphysics_softbody::{DistanceConstraint, SoftBody, SoftConstraint, SoftParticle, XpbdSolver};

fn main() {
    let nx = 10_usize;
    let ny = 10_usize;
    let width = 2.0_f64; // metres
    let height = 2.0_f64;
    let mass_per_particle = 0.05_f64; // 50 g per node
    let compliance = 0.0001_f64; // fairly stiff cloth

    let dx = width / (nx - 1) as f64;
    let dy = height / (ny - 1) as f64;

    // Create particles in the XZ plane at y=2 (will drape downward).
    let mut particles = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let pos = Vec3::new(i as f64 * dx, 2.0, j as f64 * dy);
            // Top row (j == 0) is fixed.
            if j == 0 {
                particles.push(SoftParticle::new_static(pos));
            } else {
                particles.push(SoftParticle::new(pos, mass_per_particle));
            }
        }
    }

    let idx = |ix: usize, iy: usize| -> usize { iy * nx + ix };

    // Build distance constraints: structural (horizontal + vertical) + shear.
    let mut constraints: Vec<Box<dyn SoftConstraint>> = Vec::new();

    for j in 0..ny {
        for i in 0..nx {
            // Horizontal neighbour.
            if i + 1 < nx {
                constraints.push(Box::new(DistanceConstraint::from_particles(
                    idx(i, j),
                    idx(i + 1, j),
                    &particles,
                    compliance,
                )));
            }
            // Vertical neighbour.
            if j + 1 < ny {
                constraints.push(Box::new(DistanceConstraint::from_particles(
                    idx(i, j),
                    idx(i, j + 1),
                    &particles,
                    compliance,
                )));
            }
            // Diagonal shear.
            if i + 1 < nx && j + 1 < ny {
                constraints.push(Box::new(DistanceConstraint::from_particles(
                    idx(i, j),
                    idx(i + 1, j + 1),
                    &particles,
                    compliance,
                )));
            }
        }
    }

    let mut body = SoftBody::from_particles(particles);
    body.damping = 0.02;

    // Gravity force per particle: F = m * g.
    let gravity = Vec3::new(0.0, -9.81 * mass_per_particle, 0.0);

    let mut solver = XpbdSolver::with_iterations(4, 10);

    let dt = 1.0 / 60.0;
    let n_steps = 300_usize;

    println!("=== cloth_simulation: {nx}x{ny} XPBD cloth ===");
    println!("Constraints: {}", constraints.len());
    println!("Running {} steps at dt={:.4} s ...", n_steps, dt);

    for _step in 0..n_steps {
        // Apply gravity to dynamic particles.
        body.clear_forces();
        body.apply_force(&gravity);

        solver.solve(&mut body, &mut constraints, dt);
    }

    // Print final positions of a subset of particles.
    println!("\nFinal particle positions (selected rows):");
    println!(
        "{:>4}  {:>4}  {:>10}  {:>10}  {:>10}",
        "i", "j", "x", "y", "z"
    );
    for j in (0..ny).step_by(3) {
        for i in (0..nx).step_by(3) {
            let p = &body.particles[idx(i, j)];
            println!(
                "{:>4}  {:>4}  {:>10.4}  {:>10.4}  {:>10.4}",
                i, j, p.position.x, p.position.y, p.position.z
            );
        }
    }

    // Verify: top row should remain fixed, bottom row should have dropped.
    let top_y = body.particles[idx(0, 0)].position.y;
    let bottom_y = body.particles[idx(nx / 2, ny - 1)].position.y;
    let center_y = body.particles[idx(nx / 2, ny / 2)].position.y;

    println!("\nTop-left y:    {:.4}", top_y);
    println!("Center y:      {:.4}", center_y);
    println!("Bottom-mid y:  {:.4}", bottom_y);

    if top_y > center_y && center_y > bottom_y {
        println!("PASS: cloth drapes correctly (top > center > bottom)");
    } else if (top_y - center_y).abs() < 1e-6 {
        println!("INFO: cloth barely moved (increase steps or reduce compliance)");
    } else {
        println!("WARN: unexpected draping pattern");
    }
}
