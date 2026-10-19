//! Shared helpers of the collision-operator regression tests.

use std::f64::consts::PI;

use super::functions::{compute_equilibrium_d3q19, compute_macroscopic};
use crate::lattice::{D3Q19_VELOCITIES, D3Q19_WEIGHTS};

/// Deterministic pseudo-random populations: the equilibrium of a moving fluid multiplied by
/// `1 + amplitude * r_i` with `r_i` uniform in `[-1, 1)`.
pub(super) fn perturbed_populations(seed: u64, amplitude: f64) -> [f64; 19] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let feq = compute_equilibrium_d3q19(1.15, [0.03, -0.02, 0.025]);
    std::array::from_fn(|i| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        feq[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
    })
}

/// Deterministic pseudo-random populations around the rest state: `w_i (1 + amplitude * r_i)`
/// with `r_i` uniform in `[-1, 1)`.
pub(super) fn perturbed_rest_populations(seed: u64, amplitude: f64) -> [f64; 19] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    std::array::from_fn(|i| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        D3Q19_WEIGHTS[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
    })
}

/// Equilibrium of a moving fluid whose six face populations (directions 1..=6) are scaled by 1.01.
pub(super) fn face_scaled_equilibrium() -> [f64; 19] {
    let mut f = compute_equilibrium_d3q19(1.0, [0.03, -0.02, 0.01]);
    for f_i in f.iter_mut().skip(1).take(6) {
        *f_i *= 1.01;
    }
    f
}

/// Mass and momentum `[rho, jx, jy, jz]` of a population set.
pub(super) fn mass_and_momentum(f: &[f64; 19]) -> [f64; 4] {
    let mut m = [0.0_f64; 4];
    for (&fi, c) in f.iter().zip(D3Q19_VELOCITIES.iter()) {
        m[0] += fi;
        m[1] += fi * f64::from(c[0]);
        m[2] += fi * f64::from(c[1]);
        m[3] += fi * f64::from(c[2]);
    }
    m
}

/// Assert that mass and momentum of `after` equal those of `before` to `tolerance`.
pub(super) fn assert_conserved(before: &[f64; 19], after: &[f64; 19], tolerance: f64, ctx: &str) {
    let a = mass_and_momentum(before);
    let b = mass_and_momentum(after);
    for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
        assert!(
            (a[c] - b[c]).abs() < tolerance,
            "{ctx}: {name} changed from {} to {}",
            a[c],
            b[c]
        );
    }
}

/// Amplitude ratio `A(steps) / A(0)` of a periodic shear wave `u_x = U sin(2 pi y / ny)` evolved
/// with the given collision operator on a quasi-one-dimensional D3Q19 domain (one node in x and
/// z, `ny` nodes in y, periodic streaming).
///
/// The collision closure receives `(f, rho, u)` of a node and returns its post-collision
/// populations. A collision operator with kinematic viscosity `nu` must reproduce
/// `exp(-nu (2 pi / ny)^2 steps)` up to the small discretisation error of BGK itself.
pub(super) fn shear_wave_amplitude_ratio<F>(
    collide: F,
    ny: usize,
    amplitude: f64,
    steps: usize,
) -> f64
where
    F: Fn(&[f64; 19], f64, [f64; 3]) -> [f64; 19],
{
    let k = 2.0 * PI / ny as f64;
    let mut cells: Vec<[f64; 19]> = (0..ny)
        .map(|j| compute_equilibrium_d3q19(1.0, [amplitude * (k * j as f64).sin(), 0.0, 0.0]))
        .collect();
    let measure = |cells: &[[f64; 19]]| -> f64 {
        let mut acc = 0.0;
        for (j, cell) in cells.iter().enumerate() {
            let (_, u) = compute_macroscopic(cell);
            acc += u[0] * (k * j as f64).sin();
        }
        2.0 * acc / ny as f64
    };
    let initial = measure(&cells);
    for _ in 0..steps {
        let mut next = vec![[0.0_f64; 19]; ny];
        for (j, cell) in cells.iter().enumerate() {
            let (rho, u) = compute_macroscopic(cell);
            let post = collide(cell, rho, u);
            for (i, &value) in post.iter().enumerate() {
                let dest = (j as isize + D3Q19_VELOCITIES[i][1] as isize).rem_euclid(ny as isize);
                next[dest as usize][i] = value;
            }
        }
        cells = next;
    }
    measure(&cells) / initial
}

/// Outcome of [`double_shear_layer_survival`].
pub(super) struct ShearLayerOutcome {
    /// Number of steps completed before a population became non-positive (`steps` if none did).
    pub completed_steps: usize,
    /// Change of the total mass over the completed steps.
    pub mass_drift: f64,
    /// Change of the total momentum over the completed steps.
    pub momentum_drift: [f64; 3],
}

/// Under-resolved double shear layer on an `n x n` periodic lattice (one node in z):
/// `u_x = u0 tanh(80 (y - 1/4))` for `y <= 1/2` and `u0 tanh(80 (3/4 - y))` above, with a
/// sinusoidal `u_y` perturbation of relative size 5 %. The layers are far thinner than a lattice
/// spacing, which drives BGK at low viscosity to negative populations within a few hundred
/// steps.
pub(super) fn double_shear_layer_survival<F>(
    collide: F,
    n: usize,
    u0: f64,
    steps: usize,
) -> ShearLayerOutcome
where
    F: Fn(&[f64; 19], f64, [f64; 3]) -> [f64; 19],
{
    let totals = |cells: &[[f64; 19]]| -> (f64, [f64; 3]) {
        let mut mass = 0.0;
        let mut momentum = [0.0_f64; 3];
        for cell in cells {
            let m = mass_and_momentum(cell);
            mass += m[0];
            for (a, b) in momentum.iter_mut().zip(m[1..].iter()) {
                *a += b;
            }
        }
        (mass, momentum)
    };
    let mut cells: Vec<[f64; 19]> = Vec::with_capacity(n * n);
    for iy in 0..n {
        for ix in 0..n {
            let x = ix as f64 / n as f64;
            let y = iy as f64 / n as f64;
            let ux = if y <= 0.5 {
                u0 * (80.0 * (y - 0.25)).tanh()
            } else {
                u0 * (80.0 * (0.75 - y)).tanh()
            };
            let uy = 0.05 * u0 * (2.0 * PI * (x + 0.25)).sin();
            cells.push(compute_equilibrium_d3q19(1.0, [ux, uy, 0.0]));
        }
    }
    let (mass0, momentum0) = totals(&cells);
    let mut completed = 0;
    for _ in 0..steps {
        if cells
            .iter()
            .any(|cell| cell.iter().any(|&v| v <= 0.0 || !v.is_finite()))
        {
            break;
        }
        let mut next = vec![[0.0_f64; 19]; n * n];
        for iy in 0..n {
            for ix in 0..n {
                let cell = &cells[ix + iy * n];
                let (rho, u) = compute_macroscopic(cell);
                let post = collide(cell, rho, u);
                for (i, &value) in post.iter().enumerate() {
                    let c = D3Q19_VELOCITIES[i];
                    let dx = (ix as isize + c[0] as isize).rem_euclid(n as isize) as usize;
                    let dy = (iy as isize + c[1] as isize).rem_euclid(n as isize) as usize;
                    next[dx + dy * n][i] = value;
                }
            }
        }
        cells = next;
        completed += 1;
    }
    let (mass, momentum) = totals(&cells);
    ShearLayerOutcome {
        completed_steps: completed,
        mass_drift: mass - mass0,
        momentum_drift: std::array::from_fn(|k| momentum[k] - momentum0[k]),
    }
}

/// Analytic decay `exp(-nu k^2 steps)` of a shear wave with wave number `2 pi / ny`.
pub(super) fn shear_wave_decay(nu: f64, ny: usize, steps: usize) -> f64 {
    let k = 2.0 * PI / ny as f64;
    (-nu * k * k * steps as f64).exp()
}
