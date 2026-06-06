// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Stress computation utilities and contact handling.

use super::math_helpers::{det3x3, mul3x3, transpose3x3};

// ---------------------------------------------------------------------------
// Stress computation utilities
// ---------------------------------------------------------------------------

/// Compute the Cauchy stress tensor sigma = (1/J) P F^T from first Piola-Kirchhoff stress P.
pub fn cauchy_stress(piola: [[f64; 3]; 3], f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = det3x3(f);
    if j.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let ft = transpose3x3(f);
    let pft = mul3x3(piola, ft);
    let inv_j = 1.0 / j;
    let mut sigma = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            sigma[i][k] = inv_j * pft[i][k];
        }
    }
    sigma
}

/// Compute the von Mises stress from a Cauchy stress tensor.
pub fn von_mises_stress(sigma: [[f64; 3]; 3]) -> f64 {
    let s11 = sigma[0][0];
    let s22 = sigma[1][1];
    let s33 = sigma[2][2];
    let s12 = sigma[0][1];
    let s23 = sigma[1][2];
    let s13 = sigma[0][2];

    let term1 = (s11 - s22) * (s11 - s22) + (s22 - s33) * (s22 - s33) + (s33 - s11) * (s33 - s11);
    let term2 = 6.0 * (s12 * s12 + s23 * s23 + s13 * s13);

    ((term1 + term2) / 2.0).sqrt()
}

/// Compute the Green-Lagrange strain tensor E = 0.5 (F^T F - I).
pub fn green_lagrange_strain(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ft = transpose3x3(f);
    let ftf = mul3x3(ft, f);
    let mut e = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            e[i][j] = 0.5 * (ftf[i][j] - delta);
        }
    }
    e
}

// ---------------------------------------------------------------------------
// Contact handling (simple penalty-based ground plane contact)
// ---------------------------------------------------------------------------

/// Simple penalty-based contact response against a ground plane.
pub fn apply_ground_contact(
    nodes: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    ground_y: f64,
    stiffness: f64,
    damping: f64,
    dt: f64,
) {
    for i in 0..nodes.len() {
        let penetration = ground_y - nodes[i][1];
        if penetration > 0.0 {
            // Penalty force
            let force_y = stiffness * penetration - damping * velocities[i][1];
            velocities[i][1] += force_y * dt;
            nodes[i][1] = ground_y;
        }
    }
}
