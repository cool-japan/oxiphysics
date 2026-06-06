//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BondOrder, OverCoordParams, OxidationState, ReactiveSite, ReaxFFBond};

/// Atom index pair representing a bond.
type BondPair = (usize, usize);

/// Morse pair potential.
///
/// ```text
/// V(r) = D_e * (1 - exp(-beta*(r - r0)))^2 - D_e
/// ```
///
/// Minimum value -D_e is at r = r0.
pub fn morse_potential(r: f64, params: &BondOrder) -> f64 {
    let x = (-params.beta * (r - params.r0)).exp();
    params.d_e * (1.0 - x) * (1.0 - x) - params.d_e
}
/// Derivative of the Morse potential with respect to r (i.e. the scalar force
/// along the bond axis, positive = repulsive).
///
/// ```text
/// dV/dr = 2 * D_e * beta * exp(-beta*(r-r0)) * (1 - exp(-beta*(r-r0)))
/// ```
pub fn morse_force(r: f64, params: &BondOrder) -> f64 {
    let x = (-params.beta * (r - params.r0)).exp();
    2.0 * params.d_e * params.beta * x * (1.0 - x)
}
/// Morse potential weighted by bond order: E = BO * V_morse(r).
pub fn morse_potential_weighted(r: f64, params: &BondOrder, bond_order: f64) -> f64 {
    bond_order * morse_potential(r, params)
}
/// Morse force weighted by bond order: F = BO * dV/dr.
pub fn morse_force_weighted(r: f64, params: &BondOrder, bond_order: f64) -> f64 {
    bond_order * morse_force(r, params)
}
/// Three-body harmonic angle energy.
///
/// ```text
/// E_angle = 0.5 * k_theta * (theta - theta_0)^2
/// ```
pub fn three_body_angle_energy(
    _r_ij: f64,
    _r_jk: f64,
    theta: f64,
    k_theta: f64,
    theta_0: f64,
) -> f64 {
    let d = theta - theta_0;
    0.5 * k_theta * d * d
}
/// Cosine-based angle energy: E = k * (1 - cos(theta - theta_0)).
///
/// Smoother than harmonic for large deviations.
pub fn cosine_angle_energy(theta: f64, k_theta: f64, theta_0: f64) -> f64 {
    k_theta * (1.0 - (theta - theta_0).cos())
}
/// Compute the ReaxFF bond order at distance `r`.
///
/// ```text
/// BO(r) = exp(p_bo1 * (r/r0)^p_bo2) + exp(p_bo3 * (r/r0)^p_bo4)
/// ```
pub fn compute_bond_order(r: f64, bond: &ReaxFFBond) -> f64 {
    let rho = r / bond.r0;
    let sigma = (bond.p_bo1 * rho.powf(bond.p_bo2)).exp();
    let pi = (bond.p_bo3 * rho.powf(bond.p_bo4)).exp();
    sigma + pi
}
/// Compute the derivative of bond order with respect to r.
///
/// ```text
/// dBO/dr = (p_bo1*p_bo2/r0) * (r/r0)^(p_bo2-1) * exp(p_bo1*(r/r0)^p_bo2)
///        + (p_bo3*p_bo4/r0) * (r/r0)^(p_bo4-1) * exp(p_bo3*(r/r0)^p_bo4)
/// ```
pub fn bond_order_derivative(r: f64, bond: &ReaxFFBond) -> f64 {
    let rho = r / bond.r0;
    let dsigma = bond.p_bo1 * bond.p_bo2 / bond.r0
        * rho.powf(bond.p_bo2 - 1.0)
        * (bond.p_bo1 * rho.powf(bond.p_bo2)).exp();
    let dpi = bond.p_bo3 * bond.p_bo4 / bond.r0
        * rho.powf(bond.p_bo4 - 1.0)
        * (bond.p_bo3 * rho.powf(bond.p_bo4)).exp();
    dsigma + dpi
}
/// Compute the ReaxFF-style bond energy at distance `r`.
///
/// Uses a Morse-like expression:
/// ```text
/// E(r) = D_e * (1 - exp(-beta*(r - r0)))^2 - D_e
/// ```
pub fn compute_reaxff_bond_energy(r: f64, bond: &ReaxFFBond) -> f64 {
    let x = (-bond.beta * (r - bond.r0)).exp();
    bond.d_e * (1.0 - x) * (1.0 - x) - bond.d_e
}
/// Smooth the raw bond order using a cutoff function.
///
/// ```text
/// BO_smooth = BO * (1 - (BO / bo_max)^n)  if BO < bo_max
///           = 0                             otherwise
/// ```
///
/// This ensures bond orders smoothly go to zero above bo_max.
pub fn smooth_bond_order(bo: f64, bo_max: f64, n: f64) -> f64 {
    if bo <= 0.0 || bo >= bo_max {
        return 0.0;
    }
    bo * (1.0 - (bo / bo_max).powf(n))
}
/// Build a bond list from atom positions: pairs within `max_bond_dist` are bonded.
pub fn topology_from_distances(positions: &[[f64; 3]], max_bond_dist: f64) -> Vec<(usize, usize)> {
    let n = positions.len();
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if dist(positions[i], positions[j]) < max_bond_dist {
                bonds.push((i, j));
            }
        }
    }
    bonds
}
/// Redistribute charges based on bond orders using electronegativity equalization.
///
/// This simplified model assigns charge to each atom proportional to its
/// total bond order relative to the system average.
///
/// ```text
/// q_i = q_total * (TBO_i / sum_j TBO_j) - chi_i * sum_j (TBO_ij / r_ij)
/// ```
///
/// For simplicity, we use: q_i = q_ref * (1 - TBO_i / TBO_max)
/// where TBO_max is the maximum expected total bond order.
pub fn redistribute_charges(
    sites: &[ReactiveSite],
    base_charges: &[f64],
    bo_scale: f64,
) -> Vec<f64> {
    let n = base_charges.len();
    let mut charges = base_charges.to_vec();
    for site in sites {
        let idx = site.central_atom;
        if idx < n {
            let tbo = site.total_bond_order();
            charges[idx] = base_charges[idx] * (1.0 - bo_scale * tbo).max(0.0);
        }
    }
    charges
}
/// Update oxidation states based on bond orders and electronegativities.
///
/// For each bond, the more electronegative atom gets a negative contribution,
/// the less electronegative atom gets a positive contribution, scaled by BO.
pub fn update_oxidation_states(
    states: &mut [OxidationState],
    bonds: &[(usize, usize)],
    bond_orders: &[f64],
) {
    for s in states.iter_mut() {
        s.state = 0;
    }
    for (bond_idx, &(i, j)) in bonds.iter().enumerate() {
        if bond_idx >= bond_orders.len() {
            break;
        }
        let bo = bond_orders[bond_idx];
        let bo_int = bo.round() as i32;
        if i < states.len() && j < states.len() {
            if states[i].electronegativity > states[j].electronegativity {
                states[i].state -= bo_int;
                states[j].state += bo_int;
            } else {
                states[j].state -= bo_int;
                states[i].state += bo_int;
            }
        }
    }
}
pub(super) fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Four-body (dihedral) torsion potential using a cosine Fourier series.
///
/// ```text
/// E_torsion = Σ_n V_n / 2 * (1 + cos(n * phi - gamma_n))
/// ```
///
/// where `coefficients[n] = (V_n, gamma_n)`.
pub fn torsion_energy(phi: f64, coefficients: &[(f64, f64)]) -> f64 {
    coefficients
        .iter()
        .enumerate()
        .map(|(n, &(v_n, gamma))| v_n / 2.0 * (1.0 + ((n as f64 + 1.0) * phi - gamma).cos()))
        .sum()
}
/// Derivative of the torsion energy with respect to the dihedral angle phi.
///
/// `dE/dphi = Σ_n -V_n/2 * n * sin(n * phi - gamma_n)`
pub fn torsion_energy_derivative(phi: f64, coefficients: &[(f64, f64)]) -> f64 {
    coefficients
        .iter()
        .enumerate()
        .map(|(n, &(v_n, gamma))| {
            let m = (n as f64) + 1.0;
            -v_n / 2.0 * m * (m * phi - gamma).sin()
        })
        .sum()
}
/// Lennard-Jones non-bonded energy with 1-4 scaling (for reactive MD).
///
/// ```text
/// E_LJ = 4ε * [(σ/r)^12 - (σ/r)^6]
/// ```
///
/// A pair-specific scaling factor `scale` (typically 0.5 for 1-4 pairs) is applied.
pub fn lj_nonbonded(r: f64, epsilon: f64, sigma: f64, scale: f64) -> f64 {
    if r < 1e-14 {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    scale * 4.0 * epsilon * (sr6 * sr6 - sr6)
}
/// Lennard-Jones force magnitude along the bond axis (∂E/∂r).
///
/// Positive values are repulsive.
pub fn lj_nonbonded_force(r: f64, epsilon: f64, sigma: f64, scale: f64) -> f64 {
    if r < 1e-14 {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    scale * 4.0 * epsilon * (-12.0 * sr6 * sr6 + 6.0 * sr6) / r
}
/// ReaxFF over/under-coordination energy penalty for atom i.
///
/// Penalises deviation of the total bond order `total_bo` from the
/// expected valence `valence`:
///
/// ```text
/// E_penalty = p_ovun1 * (total_bo - valence)^2 / (total_bo + eps)
/// ```
pub fn reaxff_valence_penalty(total_bo: f64, valence: f64, p_ovun1: f64) -> f64 {
    let eps = 1e-12;
    let delta = total_bo - valence;
    p_ovun1 * delta * delta / (total_bo + eps)
}
/// Derivative of the over/under-coordination penalty with respect to `total_bo`.
pub fn reaxff_valence_penalty_deriv(total_bo: f64, valence: f64, p_ovun1: f64) -> f64 {
    let eps = 1e-12;
    let delta = total_bo - valence;
    let denom = total_bo + eps;
    p_ovun1 * (2.0 * delta / denom - delta * delta / (denom * denom))
}
/// Fermi-Dirac-style bond-order switching function.
///
/// Provides a smooth cutoff between `r_inner` and `r_outer`:
/// ```text
/// S(r) = 1 / (1 + exp(a * (r - r_mid)))
/// ```
/// where `r_mid = (r_inner + r_outer) / 2` and `a` is chosen so that
/// `S(r_inner) ≈ 1` and `S(r_outer) ≈ 0`.
pub fn rebo_switching(r: f64, r_inner: f64, r_outer: f64) -> f64 {
    if r <= r_inner {
        return 1.0;
    }
    if r >= r_outer {
        return 0.0;
    }
    let r_mid = 0.5 * (r_inner + r_outer);
    let width = (r_outer - r_inner).max(1e-14);
    let a = 10.0 / width;
    1.0 / (1.0 + (a * (r - r_mid)).exp())
}
/// Derivative of the REBO switching function with respect to r.
pub fn rebo_switching_derivative(r: f64, r_inner: f64, r_outer: f64) -> f64 {
    if r <= r_inner || r >= r_outer {
        return 0.0;
    }
    let r_mid = 0.5 * (r_inner + r_outer);
    let width = (r_outer - r_inner).max(1e-14);
    let a = 10.0 / width;
    let e = (a * (r - r_mid)).exp();
    let s = 1.0 / (1.0 + e);
    -a * e * s * s
}
/// Conjugation correction for pi bonds (simplified Brenner model).
///
/// Returns an additive correction to the bond energy when two adjacent bonds
/// are both partial double bonds (each with pi-BO near 0.5):
///
/// ```text
/// E_conj = -F_conj * exp(-k * (BO_pi_1 - 0.5)^2) * exp(-k * (BO_pi_2 - 0.5)^2)
/// ```
pub fn conjugation_correction(bo_pi_1: f64, bo_pi_2: f64, f_conj: f64, k: f64) -> f64 {
    -f_conj * (-(k * (bo_pi_1 - 0.5).powi(2))).exp() * (-(k * (bo_pi_2 - 0.5).powi(2))).exp()
}
/// Electronegativity equalization for charges in a reactive system.
///
/// Iteratively updates charges on each atom based on the EEM model:
/// ```text
/// chi_i + η_i * q_i + Σ_{j≠i} J_ij * q_j = λ  (for all i)
/// Σ_i q_i = Q_total
/// ```
///
/// Uses a simple Jacobi-style iteration for the equalisation.
/// Returns updated charges; convergence to `tol` or `max_iter` iterations.
pub fn eem_charges(
    positions: &[[f64; 3]],
    electronegativities: &[f64],
    hardnesses: &[f64],
    total_charge: f64,
    tol: f64,
    max_iter: usize,
) -> Vec<f64> {
    let n = positions.len();
    if n == 0 {
        return Vec::new();
    }
    let ke = 14.3996_f64;
    let mut q = vec![total_charge / n as f64; n];
    for _ in 0..max_iter {
        let lambda = {
            let sum_chi_eta: f64 = (0..n)
                .map(|i| electronegativities[i] + hardnesses[i] * q[i])
                .sum::<f64>();
            let coulomb_sum: f64 = (0..n)
                .map(|i| {
                    (0..n)
                        .filter(|&j| j != i)
                        .map(|j| {
                            let r = dist(positions[i], positions[j]);
                            if r > 1e-10 { ke * q[j] / r } else { 0.0 }
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            (sum_chi_eta + coulomb_sum) / n as f64
        };
        let mut q_new = vec![0.0; n];
        for i in 0..n {
            let coulomb_i: f64 = (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    let r = dist(positions[i], positions[j]);
                    if r > 1e-10 { ke * q[j] / r } else { 0.0 }
                })
                .sum();
            q_new[i] = (lambda - electronegativities[i] - coulomb_i) / hardnesses[i].max(1e-12);
        }
        let q_sum: f64 = q_new.iter().sum();
        if q_sum.abs() > 1e-12 {
            let scale = total_charge / q_sum;
            for v in q_new.iter_mut() {
                *v *= scale;
            }
        }
        let delta: f64 = q
            .iter()
            .zip(q_new.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        q = q_new;
        if delta < tol {
            break;
        }
    }
    q
}
/// ReaxFF-style valence angle energy.
///
/// ```text
/// E_val = f_ang * exp(-p_val * (BO_ij - 1)^2)
///       * [theta_0(BO) - theta]^2
/// ```
///
/// A simplified model: `theta_0` depends on bond orders, approximated as
/// `theta_0 = pi * (1 - 0.5 * (BO_ij + BO_jk - 2))` (capped to \[0, pi\]).
pub fn reaxff_valence_angle_energy(
    theta: f64,
    bo_ij: f64,
    bo_jk: f64,
    f_ang: f64,
    p_val: f64,
) -> f64 {
    use std::f64::consts::PI;
    let bo_factor =
        (-(p_val * (bo_ij - 1.0).powi(2))).exp() * (-(p_val * (bo_jk - 1.0).powi(2))).exp();
    let theta_0 = (PI - 0.5 * (bo_ij + bo_jk - 2.0).abs()).clamp(0.0, PI);
    let d_theta = theta_0 - theta;
    f_ang * bo_factor * d_theta * d_theta
}
/// Cosine representation of the valence angle energy (smoother form).
///
/// ```text
/// E_val_cos = k_ang * (1 - cos(theta - theta_0))
/// ```
pub fn valence_angle_cosine_energy(theta: f64, theta_0: f64, k_ang: f64) -> f64 {
    k_ang * (1.0 - (theta - theta_0).cos())
}
/// Named-parameter wrapper for torsion energy.
///
/// Computes the OPLS-style torsion:
/// ```text
/// E = V1/2*(1+cos(phi)) + V2/2*(1-cos(2*phi)) + V3/2*(1+cos(3*phi))
/// ```
pub fn opls_torsion_energy(phi: f64, v1: f64, v2: f64, v3: f64) -> f64 {
    v1 / 2.0 * (1.0 + phi.cos())
        + v2 / 2.0 * (1.0 - (2.0 * phi).cos())
        + v3 / 2.0 * (1.0 + (3.0 * phi).cos())
}
/// Simplified hydrogen-bond energy term (Dreiding model).
///
/// ```text
/// E_HB = D_hb * cos^n(theta) * [5*(r0/r)^12 - 6*(r0/r)^10]
/// ```
///
/// where `theta` is the donor-H…acceptor angle and `r` is the H…acceptor distance.
pub fn hydrogen_bond_energy(r: f64, theta: f64, d_hb: f64, r0: f64, n: u32) -> f64 {
    if r < 1e-14 {
        return 0.0;
    }
    let rho = r0 / r;
    let radial = 5.0 * rho.powi(12) - 6.0 * rho.powi(10);
    let angular = theta.cos().powi(n as i32);
    d_hb * angular * radial
}
/// Returns true if a hydrogen bond is detected based on distance and angle criteria.
///
/// Standard IUPAC criteria: r(H…A) < 2.5 Å, theta(D-H…A) > 110°.
pub fn is_hydrogen_bond(r_ha: f64, theta_dha: f64) -> bool {
    r_ha < 2.5 && theta_dha > 110.0_f64.to_radians()
}
/// Extended charge equilibration (QEq/EEM) iteration step.
///
/// Performs one Jacobi update of all atomic charges based on the EEM equations
/// and returns the updated charges together with the residual (RMS change).
pub fn qeq_iteration_step(
    positions: &[[f64; 3]],
    charges: &[f64],
    electronegativities: &[f64],
    hardnesses: &[f64],
    total_charge: f64,
) -> (Vec<f64>, f64) {
    let n = charges.len();
    if n == 0 {
        return (Vec::new(), 0.0);
    }
    let ke = 14.3996_f64;
    let sum_chi: f64 = (0..n)
        .map(|i| electronegativities[i] + hardnesses[i] * charges[i])
        .sum();
    let sum_j: f64 = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    let r = dist(positions[i], positions[j]);
                    if r > 1e-10 { ke * charges[j] / r } else { 0.0 }
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    let lambda = (sum_chi + sum_j) / n as f64;
    let mut q_new = vec![0.0; n];
    for i in 0..n {
        let j_sum: f64 = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let r = dist(positions[i], positions[j]);
                if r > 1e-10 { ke * charges[j] / r } else { 0.0 }
            })
            .sum();
        q_new[i] = (lambda - electronegativities[i] - j_sum) / hardnesses[i].max(1e-12);
    }
    let q_sum: f64 = q_new.iter().sum();
    if q_sum.abs() > 1e-12 {
        let scale = total_charge / q_sum;
        for v in q_new.iter_mut() {
            *v *= scale;
        }
    }
    let residual = charges
        .iter()
        .zip(q_new.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    (q_new, residual)
}
/// Determine which bonds have broken or formed by comparing two topology lists.
///
/// Returns `(broken, formed)` where each is a vector of `(atom_i, atom_j)` pairs.
pub fn detect_topology_changes(
    prev_bonds: &[BondPair],
    curr_bonds: &[BondPair],
) -> (Vec<BondPair>, Vec<BondPair>) {
    let broken: Vec<(usize, usize)> = prev_bonds
        .iter()
        .filter(|b| !curr_bonds.contains(b))
        .cloned()
        .collect();
    let formed: Vec<(usize, usize)> = curr_bonds
        .iter()
        .filter(|b| !prev_bonds.contains(b))
        .cloned()
        .collect();
    (broken, formed)
}
/// Apply a bond-order cutoff to decide if each pair is bonded.
///
/// Returns a list of bonded pairs `(i, j)` with `i < j` whose bond order
/// exceeds `threshold`.
pub fn bond_list_from_bond_orders(
    positions: &[[f64; 3]],
    bond: &ReaxFFBond,
    threshold: f64,
) -> Vec<(usize, usize)> {
    let n = positions.len();
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let r = dist(positions[i], positions[j]);
            if compute_bond_order(r, bond) >= threshold {
                bonds.push((i, j));
            }
        }
    }
    bonds
}
/// Identifies distinct molecular fragments (connected components) in a bond list.
///
/// Returns a vector of length `n_atoms` where each element is the fragment index
/// of that atom (0-based).
pub fn find_molecules(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<usize> {
    let mut label = (0..n_atoms).collect::<Vec<usize>>();
    fn find(label: &mut Vec<usize>, x: usize) -> usize {
        if label[x] != x {
            label[x] = find(label, label[x]);
        }
        label[x]
    }
    for &(i, j) in bonds {
        let ri = find(&mut label, i);
        let rj = find(&mut label, j);
        if ri != rj {
            label[rj] = ri;
        }
    }
    // Path-compression pass: collect all roots first, then flatten labels.
    let roots: Vec<usize> = (0..n_atoms).map(|k| find(&mut label, k)).collect();
    for (lbl, &root) in label.iter_mut().zip(roots.iter()) {
        *lbl = root;
    }
    let mut map = std::collections::HashMap::new();
    let mut next_id = 0_usize;
    let mut result = vec![0_usize; n_atoms];
    for (k, &l) in label.iter().enumerate() {
        let frag = *map.entry(l).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
        result[k] = frag;
    }
    result
}
/// Count the number of distinct molecular fragments.
pub fn count_molecules(n_atoms: usize, bonds: &[(usize, usize)]) -> usize {
    let frags = find_molecules(n_atoms, bonds);
    frags
        .iter()
        .cloned()
        .collect::<std::collections::HashSet<_>>()
        .len()
}
/// Compute the composition (atom count per fragment) as a sorted vector of fragment sizes.
pub fn fragment_sizes(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<usize> {
    let frags = find_molecules(n_atoms, bonds);
    let mut counts = std::collections::HashMap::<usize, usize>::new();
    for &f in &frags {
        *counts.entry(f).or_insert(0) += 1;
    }
    let mut sizes: Vec<usize> = counts.values().cloned().collect();
    sizes.sort_unstable();
    sizes
}
/// ReaxFF over-coordination penalty energy for a single atom.
///
/// ```text
/// ΔBO_i = sum_j(BO_ij) − val_i
/// E_over = p_ovun2 * ΔBO_i / (1 + exp(p_ovun1 * ΔBO_i))   if ΔBO_i > 0
///        = 0                                                  otherwise
/// ```
///
/// # Arguments
/// * `bo_sum`  — sum of bond orders on the atom (Σ_j BO_ij)
/// * `params`  — over-coordination parameters
///
/// Returns the penalty energy (kJ/mol-like units, same as `p_ovun2`).
pub fn reaxff_over_coordination_penalty(bo_sum: f64, params: &OverCoordParams) -> f64 {
    let delta = bo_sum - params.val_i;
    if delta <= 0.0 {
        return 0.0;
    }
    params.p_ovun2 * delta / (1.0 + (params.p_ovun1 * delta).exp())
}
/// Derivative of the over-coordination penalty with respect to Δ BO.
pub fn reaxff_over_coordination_deriv(bo_sum: f64, params: &OverCoordParams) -> f64 {
    let delta = bo_sum - params.val_i;
    if delta <= 0.0 {
        return 0.0;
    }
    let ex = (params.p_ovun1 * delta).exp();
    let denom = 1.0 + ex;
    params.p_ovun2 * (1.0 / denom - params.p_ovun1 * delta * ex / (denom * denom))
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::reactive::BondEnergyTable;
    use crate::reactive::BondOrderPotential;
    use crate::reactive::ReactiveAtomParams;
    use crate::reactive::ReactiveMDStats;
    use crate::reactive::ReactiveSimulation;
    use std::f64::consts::PI;
    fn default_params() -> BondOrder {
        BondOrder {
            r0: 1.5,
            d_e: 4.0,
            beta: 2.0,
            s: 1.0,
        }
    }
    fn default_reaxff() -> ReaxFFBond {
        ReaxFFBond::carbon_carbon()
    }
    #[test]
    fn test_morse_minimum_at_r0() {
        let p = default_params();
        let v_min = morse_potential(p.r0, &p);
        assert!(
            (v_min - (-p.d_e)).abs() < 1e-12,
            "Morse minimum should be -D_e = {}, got {v_min}",
            -p.d_e
        );
    }
    #[test]
    fn test_morse_force_zero_at_r0() {
        let p = default_params();
        let f = morse_force(p.r0, &p);
        assert!(f.abs() < 1e-12, "Morse force at r0 should be 0, got {f}");
    }
    #[test]
    fn test_morse_repulsive_below_r0() {
        let p = default_params();
        let v_short = morse_potential(p.r0 - 0.3, &p);
        let v_min = morse_potential(p.r0, &p);
        assert!(v_short > v_min, "Energy below r0 should exceed minimum");
    }
    #[test]
    fn test_morse_dissociation_limit() {
        let p = default_params();
        let v_inf = morse_potential(1000.0, &p);
        assert!(
            v_inf.abs() < 1e-6,
            "Morse at r->inf should -> 0, got {v_inf}"
        );
    }
    #[test]
    fn test_bond_energy_at_equilibrium() {
        let p = default_params();
        let r0 = p.r0;
        let system = ReaxFFSimple {
            atoms: vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]],
            box_lengths: [100.0, 100.0, 100.0],
            bond_params: vec![p.clone()],
            pairs: vec![(0, 1)],
        };
        let e = system.compute_bond_energy();
        assert!(
            (e - (-p.d_e)).abs() < 1e-12,
            "bond energy at r0 should be -D_e, got {e}"
        );
    }
    #[test]
    fn test_bond_forces_at_equilibrium_zero() {
        let p = default_params();
        let r0 = p.r0;
        let system = ReaxFFSimple {
            atoms: vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]],
            box_lengths: [100.0, 100.0, 100.0],
            bond_params: vec![p],
            pairs: vec![(0, 1)],
        };
        let forces = system.compute_bond_forces();
        for &f in forces.iter().flatten() {
            assert!(f.abs() < 1e-12, "force at equilibrium should be 0, got {f}");
        }
    }
    #[test]
    fn test_newton_third_law() {
        let p = default_params();
        let system = ReaxFFSimple {
            atoms: vec![[0.0, 0.0, 0.0], [p.r0 + 0.5, 0.0, 0.0]],
            box_lengths: [100.0, 100.0, 100.0],
            bond_params: vec![p],
            pairs: vec![(0, 1)],
        };
        let forces = system.compute_bond_forces();
        for (k, (&f0, &f1)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            assert!(
                (f0 + f1).abs() < 1e-12,
                "Newton III violated for component {k}"
            );
        }
    }
    #[test]
    fn test_three_body_angle_at_equilibrium() {
        let e = three_body_angle_energy(1.5, 1.5, 1.9, 100.0, 1.9);
        assert!(
            e.abs() < 1e-12,
            "angle energy at theta_0 should be 0, got {e}"
        );
    }
    #[test]
    fn test_three_body_angle_positive() {
        let e = three_body_angle_energy(1.5, 1.5, 2.0, 100.0, 1.9);
        assert!(e > 0.0, "angle energy away from theta_0 should be positive");
    }
    #[test]
    fn test_bond_order_at_equilibrium_approx_one() {
        let bond = default_reaxff();
        let bo = compute_bond_order(bond.r0, &bond);
        assert!(bo > 0.5, "Bond order at r0 should be > 0.5, got {bo}");
        assert!(bo < 5.0, "Bond order at r0 should be < 5.0, got {bo}");
    }
    #[test]
    fn test_bond_order_decays_at_large_r() {
        let bond = default_reaxff();
        let bo_near = compute_bond_order(bond.r0, &bond);
        let bo_far = compute_bond_order(bond.r0 * 5.0, &bond);
        assert!(
            bo_far < bo_near,
            "Bond order should decay at large r; near={bo_near}, far={bo_far}"
        );
        assert!(
            bo_far < 0.01,
            "Bond order at 5*r0 should be ~0, got {bo_far}"
        );
    }
    #[test]
    fn test_bond_order_always_nonneg() {
        let bond = default_reaxff();
        for i in 0..20 {
            let r = 0.5 + i as f64 * 0.3;
            let bo = compute_bond_order(r, &bond);
            assert!(
                bo >= 0.0,
                "Bond order must be non-negative at r={r}, got {bo}"
            );
        }
    }
    #[test]
    fn test_reaxff_energy_minimum_at_r0() {
        let bond = default_reaxff();
        let e_eq = compute_reaxff_bond_energy(bond.r0, &bond);
        let e_near = compute_reaxff_bond_energy(bond.r0 + 0.2, &bond);
        let e_far = compute_reaxff_bond_energy(bond.r0 + 2.0, &bond);
        assert!(e_eq < e_near, "Energy at r0 should be < energy at r0+0.2");
        assert!(e_eq < e_far, "Energy at r0 should be < energy at large r");
    }
    #[test]
    fn test_break_event_construction() {
        let evt = BondBreakingEvent::new(2, 5, 10.5, BreakReason::Stretched);
        assert_eq!(evt.atom_i, 2);
        assert_eq!(evt.atom_j, 5);
        assert!((evt.time - 10.5).abs() < 1e-12);
        assert_eq!(evt.reason, BreakReason::Stretched);
    }
    #[test]
    fn test_topology_from_distances() {
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let bonds = topology_from_distances(&positions, 2.0);
        assert_eq!(bonds.len(), 1, "Should find 1 bond, found {}", bonds.len());
        assert_eq!(bonds[0], (0, 1));
    }
    #[test]
    fn test_topology_manager_detects_break() {
        let mut mgr = TopologyManager::new(2.0, 1.5);
        mgr.current_bonds.push((0, 1));
        let positions = vec![[0.0f64, 0.0, 0.0], [3.0, 0.0, 0.0]];
        mgr.update(&positions, 1.0);
        assert_eq!(mgr.break_events.len(), 1, "Should have 1 break event");
        assert_eq!(mgr.break_events[0].atom_i, 0);
        assert_eq!(mgr.break_events[0].atom_j, 1);
    }
    #[test]
    fn test_topology_manager_forms_bond() {
        let mut mgr = TopologyManager::new(2.5, 1.5);
        let positions = vec![[0.0f64, 0.0, 0.0], [1.2, 0.0, 0.0]];
        mgr.update(&positions, 0.0);
        assert_eq!(mgr.current_bonds.len(), 1, "Should have formed 1 bond");
    }
    #[test]
    fn test_reactive_site_total_bond_order() {
        let mut site = ReactiveSite::new(0);
        site.add_neighbor(1, 1.0);
        site.add_neighbor(2, 2.0);
        site.add_neighbor(3, 1.5);
        let total = site.total_bond_order();
        assert!(
            (total - 4.5).abs() < 1e-12,
            "Total bond order should be 4.5, got {total}"
        );
    }
    #[test]
    fn test_morse_weighted_scales_with_bo() {
        let p = default_params();
        let r = p.r0 + 0.3;
        let v_full = morse_potential(r, &p);
        let v_half = morse_potential_weighted(r, &p, 0.5);
        assert!(
            (v_half - 0.5 * v_full).abs() < 1e-12,
            "Weighted potential should scale with BO"
        );
    }
    #[test]
    fn test_morse_force_weighted_scales() {
        let p = default_params();
        let r = p.r0 + 0.3;
        let f_full = morse_force(r, &p);
        let f_half = morse_force_weighted(r, &p, 0.5);
        assert!(
            (f_half - 0.5 * f_full).abs() < 1e-12,
            "Weighted force should scale with BO"
        );
    }
    #[test]
    fn test_cosine_angle_energy_at_equilibrium() {
        let e = cosine_angle_energy(1.9, 100.0, 1.9);
        assert!(
            e.abs() < 1e-12,
            "Cosine angle energy at theta_0 should be 0"
        );
    }
    #[test]
    fn test_cosine_angle_energy_positive() {
        let e = cosine_angle_energy(2.0, 100.0, 1.9);
        assert!(
            e > 0.0,
            "Cosine angle energy away from theta_0 should be positive"
        );
    }
    #[test]
    fn test_bond_order_derivative_at_r0() {
        let bond = default_reaxff();
        let eps = 1e-6;
        let bo_plus = compute_bond_order(bond.r0 + eps, &bond);
        let bo_minus = compute_bond_order(bond.r0 - eps, &bond);
        let numerical = (bo_plus - bo_minus) / (2.0 * eps);
        let analytical = bond_order_derivative(bond.r0, &bond);
        assert!(
            (analytical - numerical).abs() < 1e-3,
            "BO derivative: analytical={analytical}, numerical={numerical}"
        );
    }
    #[test]
    fn test_smooth_bond_order() {
        let bo = smooth_bond_order(0.5, 2.0, 2.0);
        assert!(bo > 0.0, "Smoothed BO should be positive for BO < bo_max");
        let bo_max = smooth_bond_order(2.0, 2.0, 2.0);
        assert!(bo_max.abs() < 1e-14, "Smoothed BO at bo_max should be 0");
        let bo_neg = smooth_bond_order(-0.1, 2.0, 2.0);
        assert!(
            bo_neg.abs() < 1e-14,
            "Smoothed BO for negative input should be 0"
        );
    }
    #[test]
    fn test_break_reason_bond_order_low() {
        let evt = BondBreakingEvent::new(0, 1, 5.0, BreakReason::BondOrderLow);
        assert_eq!(evt.reason, BreakReason::BondOrderLow);
    }
    #[test]
    fn test_reactive_site_prune_weak_bonds() {
        let mut site = ReactiveSite::new(0);
        site.add_neighbor(1, 0.1);
        site.add_neighbor(2, 1.5);
        site.add_neighbor(3, 0.05);
        site.prune_weak_bonds(0.2);
        assert_eq!(
            site.coordination_number(),
            1,
            "Only one bond above threshold"
        );
    }
    #[test]
    fn test_reactive_site_max_bond_order() {
        let mut site = ReactiveSite::new(0);
        site.add_neighbor(1, 0.5);
        site.add_neighbor(2, 2.0);
        site.add_neighbor(3, 1.0);
        assert!((site.max_bond_order() - 2.0).abs() < 1e-14);
    }
    #[test]
    fn test_topology_manager_bond_order_break() {
        let mut mgr = TopologyManager::new(10.0, 1.5);
        mgr.current_bonds.push((0, 1));
        let positions = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let bond = default_reaxff();
        mgr.update_by_bond_order(&positions, &bond, 0.01, 1.0);
        assert_eq!(mgr.num_bonds(), 0, "Bond should be broken by low BO");
        assert_eq!(mgr.num_break_events(), 1);
        assert_eq!(mgr.break_events[0].reason, BreakReason::BondOrderLow);
    }
    #[test]
    fn test_charge_redistribution() {
        let base_charges = vec![0.5, -0.3, 0.1];
        let mut site = ReactiveSite::new(0);
        site.add_neighbor(1, 1.0);
        let sites = vec![site];
        let new_charges = redistribute_charges(&sites, &base_charges, 0.5);
        assert!(
            (new_charges[0] - 0.25).abs() < 1e-12,
            "Atom 0 charge should be 0.25, got {}",
            new_charges[0]
        );
        assert!(
            (new_charges[1] - (-0.3)).abs() < 1e-12,
            "Atom 1 charge should be unchanged"
        );
    }
    #[test]
    fn test_oxidation_state_tracking() {
        let mut states = vec![
            OxidationState::new(0, 0, 3.44),
            OxidationState::new(1, 0, 2.20),
        ];
        let bonds = vec![(0, 1)];
        let bond_orders = vec![2.0];
        update_oxidation_states(&mut states, &bonds, &bond_orders);
        assert_eq!(states[0].state, -2, "O should have oxidation -2");
        assert_eq!(states[1].state, 2, "C should have oxidation +2");
    }
    #[test]
    fn test_carbon_hydrogen_bond_params() {
        let ch = ReaxFFBond::carbon_hydrogen();
        assert!((ch.r0 - 1.09).abs() < 1e-12);
        let bo = compute_bond_order(ch.r0, &ch);
        assert!(bo > 0.5, "CH bond order at r0 should be > 0.5");
    }
    #[test]
    fn test_compute_bond_orders_system() {
        let p = default_params();
        let system = ReaxFFSimple {
            atoms: vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]],
            box_lengths: [100.0, 100.0, 100.0],
            bond_params: vec![p],
            pairs: vec![(0, 1)],
        };
        let bond = default_reaxff();
        let bos = system.compute_bond_orders(&bond);
        assert_eq!(bos.len(), 1);
        assert!(bos[0] > 0.5, "BO at CC equilibrium should be > 0.5");
    }
    #[test]
    fn test_torsion_energy_zero_at_equilibrium() {
        let coeffs = vec![(2.0, 0.0)];
        let e = torsion_energy(std::f64::consts::PI, &coeffs);
        assert!(
            e.abs() < 1e-12,
            "torsion at phi=pi with gamma=0 should be 0, got {e}"
        );
    }
    #[test]
    fn test_torsion_energy_maximum() {
        let v = 3.0_f64;
        let coeffs = vec![(v, 0.0)];
        let e = torsion_energy(0.0, &coeffs);
        assert!(
            (e - v).abs() < 1e-12,
            "torsion max should be V={v}, got {e}"
        );
    }
    #[test]
    fn test_torsion_derivative_numerical() {
        let coeffs = vec![(2.0, 0.3), (1.0, 0.1)];
        let phi = 1.2_f64;
        let h = 1e-6;
        let numerical =
            (torsion_energy(phi + h, &coeffs) - torsion_energy(phi - h, &coeffs)) / (2.0 * h);
        let analytical = torsion_energy_derivative(phi, &coeffs);
        assert!(
            (numerical - analytical).abs() < 1e-6,
            "dE/dphi mismatch: {analytical} vs {numerical}"
        );
    }
    #[test]
    fn test_lj_nonbonded_energy_positive_at_short_range() {
        let e = lj_nonbonded(0.5, 1.0, 1.0, 1.0);
        assert!(e > 0.0, "LJ at r < sigma should be repulsive, got {e}");
    }
    #[test]
    fn test_lj_nonbonded_minimum() {
        let sigma = 1.0_f64;
        let eps = 1.0_f64;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let e_min = lj_nonbonded(r_min, eps, sigma, 1.0);
        assert!(
            (e_min - (-eps)).abs() < 1e-10,
            "LJ minimum should be -eps, got {e_min}"
        );
    }
    #[test]
    fn test_lj_scaling_factor() {
        let r = 2.0_f64;
        let e_full = lj_nonbonded(r, 1.0, 1.0, 1.0);
        let e_half = lj_nonbonded(r, 1.0, 1.0, 0.5);
        assert!(
            (e_half - 0.5 * e_full).abs() < 1e-12,
            "LJ scale=0.5 should halve energy"
        );
    }
    #[test]
    fn test_reaxff_valence_penalty_zero_at_equilibrium() {
        let e = reaxff_valence_penalty(4.0, 4.0, 0.1);
        assert!(
            e.abs() < 1e-10,
            "Penalty at total_bo=valence should be 0, got {e}"
        );
    }
    #[test]
    fn test_reaxff_valence_penalty_positive() {
        let e = reaxff_valence_penalty(5.0, 4.0, 1.0);
        assert!(
            e > 0.0,
            "Over-coordination penalty should be positive, got {e}"
        );
    }
    #[test]
    fn test_reaxff_valence_penalty_deriv_numerical() {
        let bo = 3.5_f64;
        let valence = 4.0_f64;
        let p = 1.0_f64;
        let h = 1e-6;
        let numerical = (reaxff_valence_penalty(bo + h, valence, p)
            - reaxff_valence_penalty(bo - h, valence, p))
            / (2.0 * h);
        let analytical = reaxff_valence_penalty_deriv(bo, valence, p);
        assert!(
            (numerical - analytical).abs() < 1e-5,
            "penalty deriv mismatch: {analytical} vs {numerical}"
        );
    }
    #[test]
    fn test_rebo_switching_endpoints() {
        let r_i = 1.2_f64;
        let r_o = 2.0_f64;
        assert!(
            rebo_switching(r_i - 0.1, r_i, r_o) > 0.99,
            "Switch should be ≈1 inside inner cutoff"
        );
        assert!(
            rebo_switching(r_o + 0.1, r_i, r_o) < 0.01,
            "Switch should be ≈0 outside outer cutoff"
        );
    }
    #[test]
    fn test_rebo_switching_monotone() {
        let r_i = 1.2_f64;
        let r_o = 2.0_f64;
        let mut prev = 1.0_f64;
        for k in 0..20 {
            let r = r_i + (r_o - r_i) * (k as f64 / 19.0);
            let s = rebo_switching(r, r_i, r_o);
            assert!(
                s <= prev + 1e-12,
                "Switch function must be monotone decreasing at r={r}"
            );
            prev = s;
        }
    }
    #[test]
    fn test_rebo_switching_derivative_numerical() {
        let r_i = 1.0_f64;
        let r_o = 2.0_f64;
        let r = 1.5_f64;
        let h = 1e-6;
        let numerical =
            (rebo_switching(r + h, r_i, r_o) - rebo_switching(r - h, r_i, r_o)) / (2.0 * h);
        let analytical = rebo_switching_derivative(r, r_i, r_o);
        assert!(
            (numerical - analytical).abs() < 1e-4,
            "Switch deriv mismatch: {analytical} vs {numerical}"
        );
    }
    #[test]
    fn test_conjugation_correction_maximum_at_half_bo() {
        let e_max = conjugation_correction(0.5, 0.5, 1.0, 10.0);
        let e_off = conjugation_correction(0.0, 0.0, 1.0, 10.0);
        assert!(
            e_max < e_off,
            "Conjugation correction should be most negative at BO_pi=0.5"
        );
    }
    #[test]
    fn test_reactive_md_stats_recording() {
        let mut stats = ReactiveMDStats::new();
        assert_eq!(stats.n_steps(), 0);
        stats.record(10.0, 300.0, 0, 1);
        stats.record(9.5, 305.0, 1, 1);
        assert_eq!(stats.n_steps(), 2);
        let mean_e = stats.mean_energy();
        assert!(
            (mean_e - 9.75).abs() < 1e-10,
            "mean energy should be 9.75, got {mean_e}"
        );
        assert_eq!(stats.total_breaks(), 1);
    }
    #[test]
    fn test_reactive_atom_params_carbon() {
        let c = ReactiveAtomParams::carbon();
        assert_eq!(c.element, "C");
        assert!((c.valence - 4.0).abs() < 1e-12);
        assert!(c.mass > 12.0 && c.mass < 13.0);
    }
    #[test]
    fn test_eem_charges_sum_to_total() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let chi = vec![4.53, 6.93, 4.53];
        let eta = vec![5.0, 8.0, 5.0];
        let total_q = 0.0_f64;
        let charges = eem_charges(&positions, &chi, &eta, total_q, 1e-6, 100);
        let sum: f64 = charges.iter().sum();
        assert!(
            (sum - total_q).abs() < 1e-4,
            "EEM charges should sum to {total_q}, got {sum}"
        );
    }
    #[test]
    fn test_lj_nonbonded_force_at_minimum_zero() {
        let sigma = 1.0_f64;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let f = lj_nonbonded_force(r_min, 1.0, sigma, 1.0);
        assert!(f.abs() < 1e-8, "LJ force at r_min should be 0, got {f}");
    }
    #[test]
    fn test_bond_energy_table_interpolation() {
        let table = BondEnergyTable::new(vec![(0.0, 0.0), (1.0, -4.0), (2.0, -6.0), (3.0, -7.5)]);
        let e = table.lookup(1.5);
        assert!(
            (e - (-5.0)).abs() < 1e-12,
            "Interpolated energy should be -5.0, got {e}"
        );
    }
    #[test]
    fn test_bond_energy_table_clamp_low() {
        let table = BondEnergyTable::new(vec![(1.0, -2.0), (2.0, -4.0)]);
        let e = table.lookup(0.0);
        assert!(
            (e - (-2.0)).abs() < 1e-12,
            "Should clamp to first entry at BO < min, got {e}"
        );
    }
    #[test]
    fn test_bond_energy_table_clamp_high() {
        let table = BondEnergyTable::new(vec![(0.0, 0.0), (2.0, -6.0)]);
        let e = table.lookup(5.0);
        assert!(
            (e - (-6.0)).abs() < 1e-12,
            "Should clamp to last entry at BO > max, got {e}"
        );
    }
    #[test]
    fn test_bond_energy_table_unsorted_input() {
        let table = BondEnergyTable::new(vec![(2.0, -6.0), (0.0, 0.0), (1.0, -3.0)]);
        let e = table.lookup(1.0);
        assert!(
            (e - (-3.0)).abs() < 1e-12,
            "Entry at BO=1 should be -3.0, got {e}"
        );
    }
    #[test]
    fn test_reaxff_valence_angle_positive() {
        let theta = PI / 3.0;
        let e = reaxff_valence_angle_energy(theta, 1.0, 1.0, 1.0, 0.5);
        assert!(
            e >= 0.0,
            "Valence angle energy should be non-negative, got {e}"
        );
    }
    #[test]
    fn test_valence_angle_cosine_at_equilibrium() {
        let theta_0 = std::f64::consts::PI / 2.0;
        let e = valence_angle_cosine_energy(theta_0, theta_0, 100.0);
        assert!(
            e.abs() < 1e-12,
            "Cosine angle energy at theta_0 should be 0, got {e}"
        );
    }
    #[test]
    fn test_valence_angle_cosine_positive_away_from_eq() {
        let theta_0 = std::f64::consts::PI / 2.0;
        let e = valence_angle_cosine_energy(theta_0 + 0.3, theta_0, 1.0);
        assert!(
            e > 0.0,
            "Cosine angle energy away from theta_0 should be positive, got {e}"
        );
    }
    #[test]
    fn test_opls_torsion_at_zero() {
        let e = opls_torsion_energy(0.0, 1.0, 2.0, 3.0);
        assert!(
            (e - (1.0 + 3.0)).abs() < 1e-12,
            "OPLS torsion at phi=0 should be V1+V3, got {e}"
        );
    }
    #[test]
    fn test_opls_torsion_finite() {
        let e = opls_torsion_energy(1.2, 1.0, 0.5, 0.3);
        assert!(e.is_finite(), "OPLS torsion should be finite, got {e}");
    }
    #[test]
    fn test_hydrogen_bond_energy_finite() {
        let e = hydrogen_bond_energy(2.0, 0.3, 5.0, 1.8, 4);
        assert!(e.is_finite(), "H-bond energy should be finite, got {e}");
    }
    #[test]
    fn test_hydrogen_bond_energy_zero_at_tiny_r() {
        let e = hydrogen_bond_energy(0.0, 0.5, 5.0, 1.8, 4);
        assert!(e.abs() < 1e-14, "H-bond energy at r=0 should be 0, got {e}");
    }
    #[test]
    fn test_is_hydrogen_bond_positive() {
        assert!(
            is_hydrogen_bond(2.0, 150.0_f64.to_radians()),
            "Should be H-bond"
        );
    }
    #[test]
    fn test_is_hydrogen_bond_negative_distance() {
        assert!(
            !is_hydrogen_bond(3.0, 150.0_f64.to_radians()),
            "Too far for H-bond"
        );
    }
    #[test]
    fn test_is_hydrogen_bond_negative_angle() {
        assert!(
            !is_hydrogen_bond(2.0, 60.0_f64.to_radians()),
            "Too small angle for H-bond"
        );
    }
    #[test]
    fn test_qeq_step_conserves_charge() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let charges = vec![0.1, -0.1];
        let chi = vec![4.53, 6.93];
        let eta = vec![5.0, 8.0];
        let total_q = 0.0_f64;
        let (q_new, _) = qeq_iteration_step(&positions, &charges, &chi, &eta, total_q);
        let sum: f64 = q_new.iter().sum();
        assert!(
            (sum - total_q).abs() < 1e-10,
            "QEq charges must sum to total: {sum}"
        );
    }
    #[test]
    fn test_qeq_step_returns_residual() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let charges = vec![0.2, -0.2];
        let chi = vec![4.53, 6.93];
        let eta = vec![5.0, 8.0];
        let (_, res) = qeq_iteration_step(&positions, &charges, &chi, &eta, 0.0);
        assert!(
            res.is_finite() && res >= 0.0,
            "Residual must be non-negative finite, got {res}"
        );
    }
    #[test]
    fn test_detect_topology_changes_break() {
        let prev = vec![(0, 1), (1, 2)];
        let curr = vec![(1, 2)];
        let (broken, formed) = detect_topology_changes(&prev, &curr);
        assert_eq!(broken.len(), 1, "Should detect 1 broken bond");
        assert_eq!(broken[0], (0, 1));
        assert!(formed.is_empty(), "No bonds should have formed");
    }
    #[test]
    fn test_detect_topology_changes_form() {
        let prev = vec![(0, 1)];
        let curr = vec![(0, 1), (0, 2)];
        let (broken, formed) = detect_topology_changes(&prev, &curr);
        assert!(broken.is_empty(), "No bonds broken");
        assert_eq!(formed.len(), 1, "Should detect 1 new bond");
        assert_eq!(formed[0], (0, 2));
    }
    #[test]
    fn test_bond_list_from_bond_orders_threshold() {
        let bond = ReaxFFBond::carbon_carbon();
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.54, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let bonds = bond_list_from_bond_orders(&positions, &bond, 0.5);
        assert!(bonds.contains(&(0, 1)), "Pair (0,1) at r0 should be bonded");
        assert!(
            !bonds.contains(&(0, 2)),
            "Pair (0,2) at large r should not be bonded"
        );
        assert!(
            !bonds.contains(&(1, 2)),
            "Pair (1,2) at large r should not be bonded"
        );
    }
    #[test]
    fn test_find_molecules_single_molecule() {
        let bonds = vec![(0, 1), (1, 2)];
        let frags = find_molecules(3, &bonds);
        assert_eq!(
            frags[0], frags[1],
            "Atoms 0 and 1 should be in same fragment"
        );
        assert_eq!(
            frags[1], frags[2],
            "Atoms 1 and 2 should be in same fragment"
        );
    }
    #[test]
    fn test_find_molecules_two_fragments() {
        let bonds = vec![(0, 1)];
        let frags = find_molecules(3, &bonds);
        assert_eq!(frags[0], frags[1], "Atoms 0 and 1 must be in same fragment");
        assert_ne!(frags[0], frags[2], "Atom 2 must be in a different fragment");
    }
    #[test]
    fn test_count_molecules() {
        let bonds = vec![(0, 1), (2, 3)];
        let n = count_molecules(4, &bonds);
        assert_eq!(n, 2, "Should have 2 fragments: (0,1) and (2,3)");
    }
    #[test]
    fn test_count_molecules_all_isolated() {
        let n = count_molecules(5, &[]);
        assert_eq!(n, 5, "5 isolated atoms → 5 fragments");
    }
    #[test]
    fn test_fragment_sizes_sorted() {
        let bonds = vec![(0, 1), (1, 2), (3, 4)];
        let sizes = fragment_sizes(6, &bonds);
        assert_eq!(
            sizes,
            vec![1, 2, 3],
            "Fragment sizes should be sorted ascending"
        );
    }
    #[test]
    fn test_over_coord_penalty_zero_at_equilibrium() {
        let params = OverCoordParams::carbon();
        let e = reaxff_over_coordination_penalty(4.0, &params);
        assert_eq!(e, 0.0, "No over-coordination: penalty should be 0");
    }
    #[test]
    fn test_over_coord_penalty_zero_under_coordinated() {
        let params = OverCoordParams::carbon();
        let e = reaxff_over_coordination_penalty(2.0, &params);
        assert_eq!(e, 0.0, "Under-coordinated: penalty should be 0");
    }
    #[test]
    fn test_over_coord_penalty_positive_when_over() {
        let params = OverCoordParams::carbon();
        let e = reaxff_over_coordination_penalty(5.5, &params);
        assert!(
            e > 0.0,
            "Over-coordinated: penalty should be positive, got {e}"
        );
    }
    #[test]
    fn test_over_coord_penalty_increases_with_over_coord() {
        let params = OverCoordParams::carbon();
        let e1 = reaxff_over_coordination_penalty(5.0, &params);
        let e2 = reaxff_over_coordination_penalty(6.0, &params);
        assert!(e2 > e1, "More over-coordinated should give higher penalty");
    }
    #[test]
    fn test_pi_bond_zero_when_bo_pi_zero() {
        let bop = BondOrderPotential::carbon_double();
        let e = bop.compute_pi_bond_contribution(0.0);
        assert_eq!(e, 0.0, "Zero pi BO should give zero pi energy");
    }
    #[test]
    fn test_pi_bond_positive_for_nonzero_bo_pi() {
        let bop = BondOrderPotential::carbon_double();
        let e = bop.compute_pi_bond_contribution(0.5);
        assert!(
            e > 0.0,
            "Pi bond energy should be positive for BO_pi=0.5, got {e}"
        );
    }
    #[test]
    fn test_pi_bond_scales_with_d_e_pi() {
        let bop1 = BondOrderPotential::carbon_double();
        let mut bop2 = bop1.clone();
        bop2.pi.d_e_pi *= 2.0;
        let e1 = bop1.compute_pi_bond_contribution(0.4);
        let e2 = bop2.compute_pi_bond_contribution(0.4);
        assert!(
            (e2 - 2.0 * e1).abs() < 1e-10,
            "Pi energy should scale with D_e_pi"
        );
    }
    #[test]
    fn test_reactive_sim_no_events_when_positions_unchanged() {
        let bond = ReaxFFBond::carbon_carbon();
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.54, 0.0, 0.0]];
        let mut sim = ReactiveSimulation::new(positions.clone(), bond, 0.5);
        let n = sim.detect_reaction_events(&positions, 1);
        assert_eq!(n, 0, "No position change → no reaction events");
    }
    #[test]
    fn test_reactive_sim_detects_bond_break() {
        let bond = ReaxFFBond::carbon_carbon();
        let positions_init = vec![[0.0_f64, 0.0, 0.0], [1.54, 0.0, 0.0]];
        let mut sim = ReactiveSimulation::new(positions_init, bond, 0.5);
        let positions_far = vec![[0.0_f64, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let n = sim.detect_reaction_events(&positions_far, 2);
        assert!(
            n > 0,
            "Moving atom far away should trigger bond-break event"
        );
        assert!(
            sim.count_breaks() > 0,
            "Should have a break in the event log"
        );
    }
    #[test]
    fn test_reactive_sim_event_log_step_recorded() {
        let bond = ReaxFFBond::carbon_carbon();
        let positions_init = vec![[0.0_f64, 0.0, 0.0], [1.54, 0.0, 0.0]];
        let mut sim = ReactiveSimulation::new(positions_init, bond, 0.5);
        let positions_far = vec![[0.0_f64, 0.0, 0.0], [100.0, 0.0, 0.0]];
        sim.detect_reaction_events(&positions_far, 42);
        if let Some(ev) = sim.event_log.first() {
            assert_eq!(ev.step, 42, "Event step should be 42, got {}", ev.step);
        }
    }
}
