//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use std::f64::consts::PI;

use super::types::{Atom, HmmConfig, HomogenizationResult, NebImage, PhaseFieldParams, UnitCell};

/// Dot product of two 3D vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Element-wise addition of two 3D vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Element-wise subtraction of two 3D vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3D vector by a scalar.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Euclidean norm of a 3D vector.
#[inline]
pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3D vector to unit length.
#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}
/// Multiply a 3×3 row-major matrix by a 3-vector.
#[inline]
pub fn mat3_vec(m: &[f64; 9], v: [f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}
/// Add two 3×3 matrices stored row-major.
#[inline]
pub fn mat3_add(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    let mut r = [0.0f64; 9];
    for i in 0..9 {
        r[i] = a[i] + b[i];
    }
    r
}
/// Scale a 3×3 matrix.
#[inline]
pub fn mat3_scale(m: &[f64; 9], s: f64) -> [f64; 9] {
    let mut r = [0.0f64; 9];
    for i in 0..9 {
        r[i] = m[i] * s;
    }
    r
}
/// Identity 3×3 matrix.
#[inline]
pub fn mat3_identity() -> [f64; 9] {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}
/// Frobenius norm of a 3×3 matrix.
#[inline]
pub fn mat3_frob(m: &[f64; 9]) -> f64 {
    m.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Compute homogenized (effective) properties of a two-phase unit cell.
///
/// Uses arithmetic/harmonic averaging for Voigt/Reuss bounds and the
/// classical Hashin-Shtrikman estimates for a two-phase composite.
///
/// # Arguments
/// * `cell`  – The periodic unit cell.
/// * `prop1` – Property of phase 1 (larger value).
/// * `prop2` – Property of phase 2 (smaller value).
pub fn homogenize_two_phase(cell: &UnitCell, prop1: f64, prop2: f64) -> HomogenizationResult {
    let n = cell.property.len() as f64;
    let n1 = cell
        .property
        .iter()
        .filter(|&&p| (p - prop1).abs() < 1e-10)
        .count();
    let f1 = n1 as f64 / n;
    let f2 = 1.0 - f1;
    let voigt = f1 * prop1 + f2 * prop2;
    let reuss = 1.0 / (f1 / prop1 + f2 / prop2);
    let hs_upper = prop1 + f2 / (1.0 / (prop2 - prop1) + f1 / (3.0 * prop1));
    let hs_lower = prop2 + f1 / (1.0 / (prop1 - prop2) + f2 / (3.0 * prop2));
    let effective = 0.5 * (hs_upper + hs_lower);
    HomogenizationResult {
        effective_property: effective,
        voigt_bound: voigt,
        reuss_bound: reuss,
        hs_upper,
        hs_lower,
        volume_fraction: f1,
    }
}
/// Solve the scalar periodic cell problem via finite differences to obtain
/// the effective conductivity in the x-direction.
///
/// The cell problem: -∇·(k ∇χ) = ∂k/∂x with periodic BCs.
/// Here we use a simple 1-D iterative solve across a column of voxels.
///
/// # Arguments
/// * `cell`       – Unit cell with nx × ny voxels (nz=1 assumed for 2-D).
/// * `max_iter`   – Maximum Gauss-Seidel iterations.
/// * `tol`        – Convergence tolerance on relative residual.
pub fn cell_problem_1d(cell: &UnitCell, max_iter: usize, tol: f64) -> Vec<f64> {
    let n = cell.nx;
    let mut chi = vec![0.0f64; n];
    for _iter in 0..max_iter {
        let mut max_change = 0.0f64;
        for i in 0..n {
            let im = (i + n - 1) % n;
            let ip = (i + 1) % n;
            let km = 0.5 * (cell.property[i] + cell.property[im]);
            let kp = 0.5 * (cell.property[i] + cell.property[ip]);
            let rhs = (kp - km) / cell.h;
            let new_val = (km * chi[im] + kp * chi[ip] - rhs * cell.h * cell.h) / (km + kp);
            max_change = max_change.max((new_val - chi[i]).abs());
            chi[i] = new_val;
        }
        if max_change < tol {
            break;
        }
    }
    chi
}
/// Compress micro-scale simulation data into a macro-scale flux estimate.
///
/// This implements the data-estimation step of HMM: run the micro solver
/// for `config.micro_steps` steps, then average the micro fluxes over the
/// steady-state window.
///
/// # Arguments
/// * `config`       – HMM configuration parameters.
/// * `macro_strain` – Applied strain at the macro quadrature point.
/// * `property_fn`  – Function giving the local micro property at position x ∈ \[0, L\].
pub fn hmm_estimate_flux<F>(config: &HmmConfig, macro_strain: f64, property_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let n = 16usize;
    let h = config.micro_domain_size / n as f64;
    let mut t = vec![0.0f64; n];
    for i in 0..n {
        let x = (i as f64 + 0.5) * h;
        t[i] = macro_strain * x;
    }
    for _step in 0..config.micro_steps {
        let mut t_new = t.clone();
        for i in 1..(n - 1) {
            let x = (i as f64) * h;
            let km = property_fn(x - 0.5 * h);
            let kp = property_fn(x + 0.5 * h);
            t_new[i] = (km * t[i - 1] + kp * t[i + 1]) / (km + kp);
        }
        t = t_new;
    }
    let flux_sum: f64 = (0..n - 1)
        .map(|i| {
            let x = (i as f64 + 0.5) * h;
            let k = property_fn(x);
            -k * (t[i + 1] - t[i]) / h
        })
        .sum();
    flux_sum / (n - 1) as f64
}
/// Lennard-Jones pair potential energy between two atoms.
///
/// `U(r) = 4ε[(σ/r)^12 - (σ/r)^6]`
pub fn lennard_jones_potential(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let s_over_r = sigma / r.max(1e-15);
    let sr6 = s_over_r.powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}
/// Lennard-Jones force magnitude (scalar, positive = repulsive).
///
/// `F(r) = 24ε/r [2(σ/r)^12 - (σ/r)^6]`
pub fn lennard_jones_force(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let s_over_r = sigma / r.max(1e-15);
    let sr6 = s_over_r.powi(6);
    24.0 * epsilon / r * (2.0 * sr6 * sr6 - sr6)
}
/// Cauchy-Born rule: compute the first Piola-Kirchhoff stress from a
/// deformation gradient `F` using a linearised atomistic model.
///
/// For a simple cubic lattice with LJ bonds, the effective stress in the
/// x-direction (1D simplification) is:
///
/// `P ≈ C * (F - 1)` where `C` is the tangent modulus at the reference bond length.
pub fn cauchy_born_stress_1d(
    f_deformation: f64,
    epsilon: f64,
    sigma: f64,
    lattice_const: f64,
) -> f64 {
    let r0 = lattice_const;
    let r = f_deformation * r0;
    let force = lennard_jones_force(r, epsilon, sigma);
    force * f_deformation
}
/// Hardy stress tensor estimator for the atomistic region.
///
/// Computes the virial stress tensor (3×3) by summing pair contributions
/// from a list of atom pairs.
///
/// # Arguments
/// * `atoms`   – Slice of all atoms.
/// * `pairs`   – Pairs `(i, j)` of interacting atoms.
/// * `epsilon` – LJ well depth.
/// * `sigma`   – LJ length scale.
/// * `volume`  – Reference volume of the atomistic region.
pub fn virial_stress_tensor(
    atoms: &[Atom],
    pairs: &[(usize, usize)],
    epsilon: f64,
    sigma: f64,
    volume: f64,
) -> [f64; 9] {
    let mut sigma_v = [0.0f64; 9];
    for &(i, j) in pairs {
        if i >= atoms.len() || j >= atoms.len() {
            continue;
        }
        let r_vec = sub3(atoms[j].pos, atoms[i].pos);
        let r = norm3(r_vec);
        if r < 1e-12 {
            continue;
        }
        let f_mag = lennard_jones_force(r, epsilon, sigma);
        let f_vec = scale3(normalize3(r_vec), f_mag);
        for a in 0..3 {
            for b in 0..3 {
                sigma_v[a * 3 + b] -= r_vec[a] * f_vec[b];
            }
        }
    }
    mat3_scale(&sigma_v, 1.0 / volume.max(1e-300))
}
/// Solve 1-D steady-state diffusion with heterogeneous effective conductivity
/// using a finite-difference discretisation.
///
/// `-d/dx [K_eff(x) dT/dx] = f(x)` with T(0) = T_left, T(L) = T_right.
///
/// # Arguments
/// * `n`          – Number of internal nodes (total nodes = n+2 including BCs).
/// * `length`     – Domain length L.
/// * `k_eff`      – Effective conductivity at each of the n internal nodes.
/// * `source`     – Source term at each internal node.
/// * `t_left`     – Left Dirichlet BC.
/// * `t_right`    – Right Dirichlet BC.
///
/// Returns the temperature field at all n+2 nodes (including boundaries).
pub fn solve_macro_diffusion_1d(
    n: usize,
    length: f64,
    k_eff: &[f64],
    source: &[f64],
    t_left: f64,
    t_right: f64,
) -> Vec<f64> {
    let h = length / (n + 1) as f64;
    let mut t = vec![0.0f64; n + 2];
    t[0] = t_left;
    t[n + 1] = t_right;
    for _iter in 0..10_000 {
        let mut max_change = 0.0f64;
        for i in 1..=n {
            let ii = i - 1;
            let km = if ii == 0 {
                k_eff[0]
            } else {
                0.5 * (k_eff[ii] + k_eff[ii - 1])
            };
            let kp = if ii + 1 >= k_eff.len() {
                k_eff[k_eff.len() - 1]
            } else {
                0.5 * (k_eff[ii] + k_eff[ii + 1])
            };
            let rhs = source[ii] * h * h;
            let t_new = (km * t[i - 1] + kp * t[i + 1] + rhs) / (km + kp);
            max_change = max_change.max((t_new - t[i]).abs());
            t[i] = t_new;
        }
        if max_change < 1e-12 {
            break;
        }
    }
    t
}
/// Compute the centre-of-mass velocity of a group of atoms.
pub fn centre_of_mass_velocity(atoms: &[Atom]) -> [f64; 3] {
    let mut mom = [0.0f64; 3];
    let mut total_mass = 0.0f64;
    for a in atoms {
        mom = add3(mom, scale3(a.vel, a.mass));
        total_mass += a.mass;
    }
    if total_mass < 1e-300 {
        [0.0; 3]
    } else {
        scale3(mom, 1.0 / total_mass)
    }
}
/// Compute the instantaneous temperature of a group of atoms via the
/// equipartition theorem: `T = 2 KE / (3 N k_B)`.
pub fn atomistic_temperature(atoms: &[Atom], k_boltzmann: f64) -> f64 {
    let ke: f64 = atoms.iter().map(|a| a.kinetic_energy()).sum();
    let n = atoms.len() as f64;
    if n < 1.0 || k_boltzmann < 1e-300 {
        0.0
    } else {
        2.0 * ke / (3.0 * n * k_boltzmann)
    }
}
/// Compute the mean squared displacement (MSD) of atoms relative to their
/// reference positions.
pub fn mean_squared_displacement(atoms: &[Atom], reference: &[[f64; 3]]) -> f64 {
    let n = atoms.len().min(reference.len());
    if n == 0 {
        return 0.0;
    }
    let msd: f64 = atoms
        .iter()
        .take(n)
        .zip(reference.iter())
        .map(|(a, r)| {
            let d = sub3(a.pos, *r);
            dot3(d, d)
        })
        .sum();
    msd / n as f64
}
/// Restriction operator: project a fine-grid vector onto a coarse grid
/// by L2-projection (average over blocks).
pub fn restriction_l2(fine: &[f64], ratio: usize) -> Vec<f64> {
    fine.chunks(ratio)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect()
}
/// Prolongation operator: interpolate a coarse-grid vector to the fine grid
/// using piecewise-linear interpolation.
pub fn prolongation_linear(coarse: &[f64], ratio: usize) -> Vec<f64> {
    let n_fine = coarse.len() * ratio;
    let mut fine = vec![0.0f64; n_fine];
    for fi in 0..n_fine {
        let ci_f = fi as f64 / ratio as f64;
        let ci0 = (ci_f as usize).min(coarse.len() - 1);
        let ci1 = (ci0 + 1).min(coarse.len() - 1);
        let alpha = ci_f - ci0 as f64;
        fine[fi] = (1.0 - alpha) * coarse[ci0] + alpha * coarse[ci1];
    }
    fine
}
/// Multigrid V-cycle for a 1-D Poisson problem.
///
/// Recursively applies pre-smoothing, coarse-grid correction, and
/// post-smoothing via Gauss-Seidel.
///
/// # Arguments
/// * `u`         – Current solution approximation (modified in place).
/// * `f`         – Right-hand side.
/// * `h`         – Grid spacing.
/// * `n_smooth`  – Number of Gauss-Seidel smoothing steps.
/// * `depth`     – Remaining recursion depth (0 = direct solve).
pub fn multigrid_vcycle(u: &mut Vec<f64>, f: &[f64], h: f64, n_smooth: usize, depth: usize) {
    let n = u.len();
    for _ in 0..n_smooth {
        for i in 1..(n - 1) {
            u[i] = 0.5 * (u[i - 1] + u[i + 1] - h * h * f[i]);
        }
    }
    if depth == 0 || n <= 3 {
        return;
    }
    let mut res = vec![0.0f64; n];
    for i in 1..(n - 1) {
        res[i] = f[i] - (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (h * h);
    }
    let coarse_res = restriction_l2(&res, 2);
    let nc = coarse_res.len();
    let mut e_coarse = vec![0.0f64; nc];
    multigrid_vcycle(&mut e_coarse, &coarse_res, 2.0 * h, n_smooth, depth - 1);
    let e_fine = prolongation_linear(&e_coarse, 2);
    for i in 0..n.min(e_fine.len()) {
        u[i] += e_fine[i];
    }
    for _ in 0..n_smooth {
        for i in 1..(n - 1) {
            u[i] = 0.5 * (u[i - 1] + u[i + 1] - h * h * f[i]);
        }
    }
}
/// Advance the Allen-Cahn phase-field `φ` by one explicit time step.
///
/// `∂φ/∂t = M [ε² ∇²φ - W φ(φ-1)(φ-0.5) × 2]`
///
/// Uses periodic boundary conditions.
///
/// # Arguments
/// * `phi`    – Phase-field (in-place update).
/// * `params` – Phase-field parameters.
/// * `dx`     – Grid spacing.
/// * `dt`     – Time step.
pub fn allen_cahn_step(phi: &mut Vec<f64>, params: &PhaseFieldParams, dx: f64, dt: f64) {
    let n = phi.len();
    let mut dphi = vec![0.0f64; n];
    for i in 0..n {
        let im = (i + n - 1) % n;
        let ip = (i + 1) % n;
        let laplacian = (phi[im] - 2.0 * phi[i] + phi[ip]) / (dx * dx);
        let p = phi[i];
        let dw = 2.0 * params.well_height * p * (p - 1.0) * (2.0 * p - 1.0);
        dphi[i] = params.mobility * (params.epsilon_sq * laplacian - dw);
    }
    for i in 0..n {
        phi[i] += dt * dphi[i];
        phi[i] = phi[i].clamp(0.0, 1.0);
    }
}
/// Compute the coarse-grained potential energy surface along a reaction coordinate.
///
/// Uses umbrella sampling with a harmonic bias potential to reconstruct
/// the free energy profile via the weighted histogram analysis.
///
/// # Arguments
/// * `samples` – Sampled values of the reaction coordinate ξ.
/// * `xi_0`    – Centre of the harmonic bias window.
/// * `k_bias`  – Spring constant of the bias potential.
/// * `beta`    – Inverse temperature 1/(kT).
pub fn umbrella_free_energy(samples: &[f64], xi_0: f64, k_bias: f64, beta: f64) -> f64 {
    let n = samples.len() as f64;
    if n < 1.0 {
        return 0.0;
    }
    let mean_bias: f64 = samples
        .iter()
        .map(|&xi| 0.5 * k_bias * (xi - xi_0).powi(2))
        .sum::<f64>()
        / n;
    let mean_xi: f64 = samples.iter().sum::<f64>() / n;
    let var_xi: f64 = samples
        .iter()
        .map(|&xi| (xi - mean_xi).powi(2))
        .sum::<f64>()
        / n;
    let entropy_term = if var_xi > 1e-300 {
        -0.5 * (2.0 * PI * var_xi).ln() / beta
    } else {
        0.0
    };
    entropy_term + mean_bias
}
/// Check scale separation: verifies that the micro length scale is at
/// least `factor` times smaller than the macro length scale.
pub fn check_scale_separation(micro_length: f64, macro_length: f64, factor: f64) -> bool {
    macro_length >= factor * micro_length
}
/// Compute the Hill-Mandel macrohomogeneity condition error.
///
/// For a periodic unit cell under macroscopic strain ε*, the condition
/// requires `<σ : δε> = σ* : δε*`. Returns the relative violation.
pub fn hill_mandel_error(
    micro_stress: &[f64],
    micro_strain: &[f64],
    macro_stress: [f64; 9],
    macro_strain: [f64; 9],
) -> f64 {
    let n = micro_stress.len().min(micro_strain.len()) as f64;
    if n < 1.0 {
        return 0.0;
    }
    let micro_power: f64 = micro_stress
        .iter()
        .zip(micro_strain.iter())
        .map(|(s, e)| s * e)
        .sum::<f64>()
        / n;
    let macro_power: f64 = macro_stress
        .iter()
        .zip(macro_strain.iter())
        .map(|(s, e)| s * e)
        .sum();
    (micro_power - macro_power).abs() / macro_power.abs().max(1e-300)
}
/// Compute harmonic bond force between two CG beads.
///
/// `F = -k (r - r0) * r_hat`
///
/// Returns force on bead `a` (force on `b` is equal and opposite).
pub fn cg_bond_force(pos_a: [f64; 3], pos_b: [f64; 3], k: f64, r0: f64) -> [f64; 3] {
    let r_vec = sub3(pos_b, pos_a);
    let r = norm3(r_vec);
    if r < 1e-12 {
        return [0.0; 3];
    }
    let f_mag = k * (r - r0);
    scale3(normalize3(r_vec), f_mag)
}
/// Linear elasticity tangent modulus from a set of atom-pair bond stiffnesses.
///
/// Under the Cauchy-Born assumption on a simple cubic lattice, the elastic
/// modulus (Young's modulus in 1D) is:
///
/// `E = (1/V₀) Σ_bonds k_bond * (r_hat ⊗ r_hat) : I`
pub fn cauchy_born_modulus(bond_stiffnesses: &[f64], bond_lengths: &[f64], volume: f64) -> f64 {
    let n = bond_stiffnesses.len().min(bond_lengths.len());
    let sum: f64 = (0..n)
        .map(|i| bond_stiffnesses[i] * bond_lengths[i].powi(2))
        .sum();
    sum / volume.max(1e-300)
}
/// Handshake zone energy correction to remove double-counting in
/// concurrent multiscale coupling.
///
/// The Arlequin method blends atomistic and continuum energies in an
/// overlap zone using a weight function `w(x) ∈ [0,1]`.
///
/// Returns blended energy density.
pub fn arlequin_blend_energy(e_atomistic: f64, e_continuum: f64, weight: f64) -> f64 {
    let w = weight.clamp(0.0, 1.0);
    w * e_atomistic + (1.0 - w) * e_continuum
}
/// Perform one NEB step with spring forces.
///
/// # Arguments
/// * `images`          – Mutable slice of NEB images.
/// * `spring_constant` – Spring constant k.
/// * `grad_fn`         – Returns -(gradient of energy) at a configuration.
/// * `dt`              – Step size.
#[allow(clippy::too_many_arguments)]
pub fn neb_step<F>(images: &mut Vec<NebImage>, spring_constant: f64, grad_fn: &F, dt: f64)
where
    F: Fn(&[f64]) -> (f64, Vec<f64>),
{
    let n = images.len();
    if n < 3 {
        return;
    }
    for img in images.iter_mut() {
        let (e, g) = grad_fn(&img.coords);
        img.energy = e;
        img.force = g;
    }
    for i in 1..(n - 1) {
        let prev = images[i - 1].coords.clone();
        let next = images[i + 1].coords.clone();
        let curr = images[i].coords.clone();
        let tau: Vec<f64> = curr
            .iter()
            .zip(prev.iter())
            .zip(next.iter())
            .map(|((&c, &p), &nx)| {
                let dp = c - p;
                let dn = nx - c;
                0.5 * (dn + dp)
            })
            .collect();
        let tau_norm: f64 = tau.iter().map(|t| t * t).sum::<f64>().sqrt().max(1e-300);
        let d_prev: f64 = curr
            .iter()
            .zip(prev.iter())
            .map(|(&c, &p)| (c - p).powi(2))
            .sum::<f64>()
            .sqrt();
        let d_next: f64 = next
            .iter()
            .zip(curr.iter())
            .map(|(&nx, &c)| (nx - c).powi(2))
            .sum::<f64>()
            .sqrt();
        let f_spring = spring_constant * (d_next - d_prev);
        let f_real = images[i].force.clone();
        let f_dot_tau: f64 = f_real
            .iter()
            .zip(tau.iter())
            .map(|(&f, &t)| f * t / tau_norm)
            .sum();
        let f_perp: Vec<f64> = f_real
            .iter()
            .zip(tau.iter())
            .map(|(&f, &t)| f - f_dot_tau * t / tau_norm)
            .collect();
        let f_spring_vec: Vec<f64> = tau.iter().map(|&t| f_spring * t / tau_norm).collect();
        let f_neb: Vec<f64> = f_perp
            .iter()
            .zip(f_spring_vec.iter())
            .map(|(fp, fs)| fp + fs)
            .collect();
        for (c, f) in images[i].coords.iter_mut().zip(f_neb.iter()) {
            *c += dt * f;
        }
    }
}
