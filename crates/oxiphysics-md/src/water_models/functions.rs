//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{WaterGeometry, WaterModelSummary, WaterModelType, WaterMolecule, WaterParams, WaterVelocities};
use crate::MdError;

pub(super) fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}
/// Lennard-Jones 12-6 potential for oxygen-oxygen interactions.
///
/// V(r) = 4*epsilon * \[ (sigma/r)^12 - (sigma/r)^6 \]
pub fn water_lj_energy(r: f64, params: &WaterParams) -> f64 {
    let sr = params.sigma_o / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    4.0 * params.epsilon_o * (sr12 - sr6)
}
/// Lennard-Jones force magnitude (radial): F(r) = 24*epsilon/r * \[2*(sigma/r)^12 - (sigma/r)^6\].
///
/// Positive = repulsive, negative = attractive.
pub fn water_lj_force(r: f64, params: &WaterParams) -> f64 {
    let sr = params.sigma_o / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    24.0 * params.epsilon_o / r * (2.0 * sr12 - sr6)
}
/// Coulomb electrostatic interaction between two point charges.
///
/// V = k * q1 * q2 / r
pub fn water_coulomb_energy(q1: f64, q2: f64, r: f64, coulomb_k: f64) -> f64 {
    coulomb_k * q1 * q2 / r
}
/// Total pair interaction energy between two 3-site water molecules.
///
/// Sums all 9 site-site interactions (O-O has LJ + Coulomb, rest Coulomb only).
pub fn water_total_pair_energy(
    mol_a: &WaterMolecule,
    mol_b: &WaterMolecule,
    params: &WaterParams,
    coulomb_k: f64,
) -> f64 {
    let qo = params.q_o;
    let qh = params.q_h;
    let r_oo = dist(&mol_a.oxygen, &mol_b.oxygen);
    let e_oo = water_lj_energy(r_oo, params) + water_coulomb_energy(qo, qo, r_oo, coulomb_k);
    let r_oh1b = dist(&mol_a.oxygen, &mol_b.hydrogen1);
    let e_oh1b = water_coulomb_energy(qo, qh, r_oh1b, coulomb_k);
    let r_oh2b = dist(&mol_a.oxygen, &mol_b.hydrogen2);
    let e_oh2b = water_coulomb_energy(qo, qh, r_oh2b, coulomb_k);
    let r_h1ao = dist(&mol_a.hydrogen1, &mol_b.oxygen);
    let e_h1ao = water_coulomb_energy(qh, qo, r_h1ao, coulomb_k);
    let r_h1ah1b = dist(&mol_a.hydrogen1, &mol_b.hydrogen1);
    let e_h1ah1b = water_coulomb_energy(qh, qh, r_h1ah1b, coulomb_k);
    let r_h1ah2b = dist(&mol_a.hydrogen1, &mol_b.hydrogen2);
    let e_h1ah2b = water_coulomb_energy(qh, qh, r_h1ah2b, coulomb_k);
    let r_h2ao = dist(&mol_a.hydrogen2, &mol_b.oxygen);
    let e_h2ao = water_coulomb_energy(qh, qo, r_h2ao, coulomb_k);
    let r_h2ah1b = dist(&mol_a.hydrogen2, &mol_b.hydrogen1);
    let e_h2ah1b = water_coulomb_energy(qh, qh, r_h2ah1b, coulomb_k);
    let r_h2ah2b = dist(&mol_a.hydrogen2, &mol_b.hydrogen2);
    let e_h2ah2b = water_coulomb_energy(qh, qh, r_h2ah2b, coulomb_k);
    e_oo + e_oh1b + e_oh2b + e_h1ao + e_h1ah1b + e_h1ah2b + e_h2ao + e_h2ah1b + e_h2ah2b
}
/// Total pair interaction energy between two TIP4P water molecules.
///
/// LJ is between oxygen sites; Coulomb is between H and M sites
/// (oxygen carries no charge in TIP4P).
pub fn water_tip4p_pair_energy(
    mol_a: &WaterMolecule,
    mol_b: &WaterMolecule,
    params: &WaterParams,
    coulomb_k: f64,
) -> f64 {
    let r_oo = dist(&mol_a.oxygen, &mol_b.oxygen);
    let e_lj = water_lj_energy(r_oo, params);
    let qh = params.q_h;
    let qm = params.q_m;
    let ma = mol_a.m_site.unwrap_or(mol_a.oxygen);
    let mb = mol_b.m_site.unwrap_or(mol_b.oxygen);
    let sites_a: [([f64; 3], f64); 3] = [(mol_a.hydrogen1, qh), (mol_a.hydrogen2, qh), (ma, qm)];
    let sites_b: [([f64; 3], f64); 3] = [(mol_b.hydrogen1, qh), (mol_b.hydrogen2, qh), (mb, qm)];
    let mut e_coul = 0.0;
    for &(pos_a, qa) in &sites_a {
        for &(pos_b, qb) in &sites_b {
            let r = dist(&pos_a, &pos_b);
            if r > 1e-12 {
                e_coul += water_coulomb_energy(qa, qb, r, coulomb_k);
            }
        }
    }
    e_lj + e_coul
}
/// Self-polarisation correction energy for SPC/E.
///
/// E_pol = (1/2*alpha) * (mu^2 - mu_gas^2)
/// where mu_gas = 1.85 D is the gas-phase dipole and alpha = 1.608 Å³.
/// This is approximately 5.22 kJ/mol per molecule for SPC/E.
pub fn spce_self_polarisation_energy() -> f64 {
    5.22
}
/// Backward-compatible SETTLE entry point.
///
/// Real SETTLE ([`settle_positions`]) needs BOTH the start-of-step reference
/// positions and the post-move unconstrained positions. This compat wrapper is
/// used when only a single molecule is available: it synthesizes an *idealized*
/// reference triangle via [`WaterMolecule::with_geometry`] (rigid TIP3P geometry
/// centered on the current oxygen) and then applies the analytic
/// [`settle_positions`] to constrain `mol` onto the rigid reference. For true MD
/// stepping, call [`settle_positions`] directly with the genuine reference.
pub fn rigid_water_settle(mol: &mut WaterMolecule, _dt: f64) -> Result<(), MdError> {
    let geom = WaterGeometry::tip3p();
    let reference = WaterMolecule::with_geometry(mol.oxygen, &geom);
    let params = WaterParams::tip3p();
    settle_positions(&reference, mol, &geom, params.mass_o, params.mass_h)
}
/// Iterative distance-projection (SHAKE-like) constraint with explicit geometry.
///
/// This is NOT the analytic SETTLE: it repeatedly rescales each bond vector to
/// its target length (up to 50 sweeps) until the O-H and H-H distances converge.
/// It does not preserve the molecular centre of mass and does not use a reference
/// orientation. For the closed-form analytic algorithm use [`settle_positions`].
pub fn rigid_water_project_geom(mol: &mut WaterMolecule, geom: &WaterGeometry) {
    let target_oh = geom.r_oh;
    let target_hh = geom.r_hh();
    for _ in 0..50 {
        project_to_distance(&mol.oxygen, &mut mol.hydrogen1, target_oh);
        project_to_distance(&mol.oxygen, &mut mol.hydrogen2, target_oh);
        project_to_distance(&mol.hydrogen1, &mut mol.hydrogen2, target_hh);
        project_to_distance(&mol.oxygen, &mut mol.hydrogen1, target_oh);
        project_to_distance(&mol.oxygen, &mut mol.hydrogen2, target_oh);
        if mol.check_constraints(geom, 1e-10) {
            break;
        }
    }
    if geom.has_m_site() {
        mol.update_m_site(geom);
    }
}
/// Iterative distance-projection constraint using the default TIP3P geometry.
///
/// Thin public alias of [`rigid_water_project_geom`] preserving the original
/// iterative entry point. NOT analytic SETTLE — see [`settle_positions`].
pub fn rigid_water_project(mol: &mut WaterMolecule, geom: &WaterGeometry) {
    rigid_water_project_geom(mol, geom);
}
/// Move `b` along the b-a axis so that |b - a| = `target_dist`.
pub(super) fn project_to_distance(a: &[f64; 3], b: &mut [f64; 3], target_dist: f64) {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    if r < 1e-12 {
        return;
    }
    let scale = target_dist / r;
    b[0] = a[0] + dx * scale;
    b[1] = a[1] + dy * scale;
    b[2] = a[2] + dz * scale;
}
/// Place `n_molecules` water molecules on a 3D cubic grid with given spacing.
///
/// Returns a vector of water molecules placed at grid points.
pub fn create_water_box(
    n_molecules: usize,
    spacing: f64,
    geom: &WaterGeometry,
) -> Vec<WaterMolecule> {
    let n_side = (n_molecules as f64).cbrt().ceil() as usize;
    let mut molecules = Vec::with_capacity(n_molecules);
    let mut count = 0;
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                if count >= n_molecules {
                    break;
                }
                let o = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                molecules.push(WaterMolecule::with_geometry(o, geom));
                count += 1;
            }
        }
    }
    molecules
}
/// Generate comparison summaries for all standard 3-site and 4-site models.
pub fn compare_water_models() -> Vec<WaterModelSummary> {
    let models = [
        WaterModelType::Tip3p,
        WaterModelType::Spc,
        WaterModelType::Spce,
        WaterModelType::Tip4p,
    ];
    models
        .iter()
        .map(|&mt| {
            let params = WaterParams::from_model_type(mt);
            let mol = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &params.geometry);
            let dipole = mol.dipole_magnitude_debye(&params);
            WaterModelSummary {
                name: params.name.clone(),
                dipole_debye: dipole,
                r_hh: params.geometry.r_hh(),
                sigma: params.sigma_o,
                epsilon: params.epsilon_o,
                has_virtual_site: params.geometry.has_m_site(),
            }
        })
        .collect()
}
/// Radial distribution function (RDF) for oxygen-oxygen distances.
///
/// Computes g_OO(r) from a list of oxygen positions in a cubic box.
///
/// # Arguments
/// * `oxygens`   – oxygen positions (Å)
/// * `box_l`     – cubic box length (Å)
/// * `r_max`     – maximum distance (Å)
/// * `n_bins`    – number of histogram bins
///
/// # Returns
/// `(r_centers, g_oo)` vectors of length `n_bins`.
pub fn oxygen_rdf(
    oxygens: &[[f64; 3]],
    box_l: f64,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = oxygens.len();
    if n < 2 || n_bins == 0 {
        return (vec![], vec![]);
    }
    let bw = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    for i in 0..n {
        for j in i + 1..n {
            let mut dx = oxygens[j][0] - oxygens[i][0];
            let mut dy = oxygens[j][1] - oxygens[i][1];
            let mut dz = oxygens[j][2] - oxygens[i][2];
            dx -= box_l * (dx / box_l).round();
            dy -= box_l * (dy / box_l).round();
            dz -= box_l * (dz / box_l).round();
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let b = (r / bw) as usize;
                if b < n_bins {
                    hist[b] += 2;
                }
            }
        }
    }
    let rho = n as f64 / (box_l * box_l * box_l);
    let r_centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * bw).collect();
    let g_oo: Vec<f64> = (0..n_bins)
        .map(|b| {
            let r = r_centers[b];
            let vol_shell = 4.0 * PI * r * r * bw;
            let n_ideal = rho * vol_shell * n as f64;
            if n_ideal < 1e-20 {
                0.0
            } else {
                hist[b] as f64 / n_ideal
            }
        })
        .collect();
    (r_centers, g_oo)
}
/// Compute the average number of hydrogen bonds per water molecule.
///
/// Uses a geometric criterion:
/// - O-O distance < `r_oo_max` (Å)
/// - O-H···O angle > `angle_min` (degrees)
///
/// Simple implementation: counts OO pairs within cutoff (not full angle check).
pub fn count_hydrogen_bonds_simple(molecules: &[WaterMolecule], box_l: f64, r_oo_max: f64) -> f64 {
    let n = molecules.len();
    if n == 0 {
        return 0.0;
    }
    let mut n_hbonds = 0u64;
    for i in 0..n {
        for j in i + 1..n {
            let mut dx = molecules[j].oxygen[0] - molecules[i].oxygen[0];
            let mut dy = molecules[j].oxygen[1] - molecules[i].oxygen[1];
            let mut dz = molecules[j].oxygen[2] - molecules[i].oxygen[2];
            dx -= box_l * (dx / box_l).round();
            dy -= box_l * (dy / box_l).round();
            dz -= box_l * (dz / box_l).round();
            let r_oo = (dx * dx + dy * dy + dz * dz).sqrt();
            if r_oo < r_oo_max {
                n_hbonds += 2;
            }
        }
    }
    n_hbonds as f64 / n as f64
}
/// Compute the tetrahedral order parameter q for each water molecule.
///
/// q = 1 − (3/8) Σ_{j<k} (cos θ_jik + 1/3)²
///
/// where the sum is over the 4 nearest neighbours j, k.
/// Returns values close to 1 for perfect tetrahedral order, ~0 for random.
pub fn tetrahedral_order_parameter(oxygens: &[[f64; 3]], box_l: f64) -> Vec<f64> {
    let n = oxygens.len();
    if n < 5 {
        return vec![0.0; n];
    }
    let mut q_values = Vec::with_capacity(n);
    for i in 0..n {
        let mut dists: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let mut dx = oxygens[j][0] - oxygens[i][0];
                let mut dy = oxygens[j][1] - oxygens[i][1];
                let mut dz = oxygens[j][2] - oxygens[i][2];
                dx -= box_l * (dx / box_l).round();
                dy -= box_l * (dy / box_l).round();
                dz -= box_l * (dz / box_l).round();
                (dx * dx + dy * dy + dz * dz, j)
            })
            .collect();
        dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let nn: Vec<usize> = dists[..4.min(dists.len())]
            .iter()
            .map(|&(_, j)| j)
            .collect();
        if nn.len() < 4 {
            q_values.push(0.0);
            continue;
        }
        let mut sum = 0.0_f64;
        for j_idx in 0..4 {
            for k_idx in j_idx + 1..4 {
                let j = nn[j_idx];
                let k = nn[k_idx];
                let mk_vec = |a: usize, b: usize| -> [f64; 3] {
                    let mut dx = oxygens[b][0] - oxygens[a][0];
                    let mut dy = oxygens[b][1] - oxygens[a][1];
                    let mut dz = oxygens[b][2] - oxygens[a][2];
                    dx -= box_l * (dx / box_l).round();
                    dy -= box_l * (dy / box_l).round();
                    dz -= box_l * (dz / box_l).round();
                    [dx, dy, dz]
                };
                let rij = mk_vec(i, j);
                let rik = mk_vec(i, k);
                let n_ij = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
                let n_ik = (rik[0] * rik[0] + rik[1] * rik[1] + rik[2] * rik[2]).sqrt();
                if n_ij < 1e-15 || n_ik < 1e-15 {
                    continue;
                }
                let cos_theta =
                    (rij[0] * rik[0] + rij[1] * rik[1] + rij[2] * rik[2]) / (n_ij * n_ik);
                let diff = cos_theta + 1.0 / 3.0;
                sum += diff * diff;
            }
        }
        q_values.push(1.0 - (3.0 / 8.0) * sum);
    }
    q_values
}
/// Compute the water dimer binding energy.
///
/// `E_bind = E_pair(A,B) - E_intramol(A) - E_intramol(B)`
///
/// For rigid models (no intramolecular energy), this is simply the pair
/// interaction energy.  Returns the binding energy in kJ/mol (negative = bound).
///
/// `coulomb_k` in kJ mol⁻¹ Å e⁻² (e.g. 1389.35 for nm-based units, or
/// 1389.35 * 0.01 for Å-based units).
pub fn water_dimer_binding_energy(
    mol_a: &WaterMolecule,
    mol_b: &WaterMolecule,
    params: &WaterParams,
    coulomb_k: f64,
) -> f64 {
    match params.model_type {
        WaterModelType::Tip4p => water_tip4p_pair_energy(mol_a, mol_b, params, coulomb_k),
        _ => water_total_pair_energy(mol_a, mol_b, params, coulomb_k),
    }
}
/// SPC/E dimer binding energy at the typical hydrogen-bond geometry.
///
/// Places molecule B at distance `r_oo` Å from molecule A along the x-axis,
/// oriented for an ideal linear hydrogen bond (H-bond donor A to acceptor B).
pub fn spce_dimer_energy_at_separation(r_oo: f64) -> f64 {
    let coulomb_k = 1389.35_f64;
    let params = WaterParams::spce();
    let mol_a = WaterMolecule::spce([0.0, 0.0, 0.0]);
    let mol_b = WaterMolecule::spce([r_oo, 0.0, 0.0]);
    water_dimer_binding_energy(&mol_a, &mol_b, &params, coulomb_k)
}
/// Estimate the dielectric constant from the Clausius-Mossotti relation.
///
/// For a liquid of molecules with dipole moment `mu` (Debye) at number density
/// `rho_n` (Å⁻³) and temperature T (K), the Clausius-Mossotti / Kirkwood
/// equation gives an estimate:
///
/// `(eps - 1) / (eps + 2) = rho_n * alpha / (3 * eps_0)`
///
/// For a polar liquid, we use the simplified Onsager-like relation:
///
/// `eps ≈ 1 + (n_mol * mu² ) / (3 * eps_0 * V * kBT)`
///
/// Here we use a simplified version with `mu` in Debye, `V` in Å³, `T` in K.
///
/// # Arguments
/// * `mu_debye`   – molecular dipole moment (Debye).
/// * `n_mol`      – number of molecules in the simulation box.
/// * `box_vol_A3` – simulation box volume (Å³).
/// * `temperature`– temperature (K).
///
/// Returns an approximate dielectric constant (dimensionless).
pub fn clausius_mossotti_dielectric(
    mu_debye: f64,
    n_mol: f64,
    box_vol_a3: f64,
    temperature: f64,
) -> f64 {
    if box_vol_a3 <= 0.0 || temperature <= 0.0 || n_mol <= 0.0 {
        return 1.0;
    }
    let mu_cm = mu_debye * 3.336e-30_f64;
    let vol_m3 = box_vol_a3 * 1e-30_f64;
    let kb = 1.380_649e-23_f64;
    let eps0 = 8.854_188e-12_f64;
    let y = n_mol * mu_cm * mu_cm / (3.0 * eps0 * vol_m3 * kb * temperature);
    1.0 + y
}
/// Kirkwood dipole correlation factor estimate.
///
/// Computes the Kirkwood g_K factor from a list of molecular dipole vectors
/// (in Debye or any consistent units).  g_K > 1 indicates parallel alignment.
pub fn kirkwood_g_factor(dipoles: &[[f64; 3]]) -> f64 {
    let n = dipoles.len();
    if n == 0 {
        return 0.0;
    }
    let mag2: Vec<f64> = dipoles
        .iter()
        .map(|d| d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        .collect();
    let mean_mag2 = mag2.iter().sum::<f64>() / n as f64;
    if mean_mag2 < 1e-30 {
        return 1.0;
    }
    let mut cross_sum = 0.0_f64;
    for i in 0..n {
        for j in i + 1..n {
            let dot = dipoles[i][0] * dipoles[j][0]
                + dipoles[i][1] * dipoles[j][1]
                + dipoles[i][2] * dipoles[j][2];
            cross_sum += 2.0 * dot;
        }
    }
    1.0 + cross_sum / (n as f64 * mean_mag2)
}
/// Radial distribution function for O-H pairs.
///
/// Computes g_OH(r) from oxygen and hydrogen positions in a cubic box.
///
/// # Arguments
/// * `oxygens`   – oxygen positions (Å)
/// * `hydrogens` – hydrogen positions (Å)
/// * `box_l`     – cubic box length (Å)
/// * `r_max`     – maximum distance (Å)
/// * `n_bins`    – number of histogram bins
///
/// # Returns
/// `(r_centers, g_oh)` vectors of length `n_bins`.
pub fn oh_rdf(
    oxygens: &[[f64; 3]],
    hydrogens: &[[f64; 3]],
    box_l: f64,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n_o = oxygens.len();
    let n_h = hydrogens.len();
    if n_o == 0 || n_h == 0 || n_bins == 0 {
        return (vec![], vec![]);
    }
    let bw = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    for o in oxygens.iter().take(n_o) {
        for h in hydrogens.iter().take(n_h) {
            let mut dx = h[0] - o[0];
            let mut dy = h[1] - o[1];
            let mut dz = h[2] - o[2];
            dx -= box_l * (dx / box_l).round();
            dy -= box_l * (dy / box_l).round();
            dz -= box_l * (dz / box_l).round();
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let b = (r / bw) as usize;
                if b < n_bins {
                    hist[b] += 1;
                }
            }
        }
    }
    let rho_h = n_h as f64 / (box_l * box_l * box_l);
    let r_centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * bw).collect();
    let g_oh: Vec<f64> = (0..n_bins)
        .map(|b| {
            let r = r_centers[b];
            let vol_shell = 4.0 * PI * r * r * bw;
            let n_ideal = rho_h * vol_shell * n_o as f64;
            if n_ideal < 1e-20 {
                0.0
            } else {
                hist[b] as f64 / n_ideal
            }
        })
        .collect();
    (r_centers, g_oh)
}
/// Radial distribution function for H-H pairs.
///
/// # Arguments
/// * `hydrogens` – all hydrogen positions (Å).
/// * `box_l`     – cubic box length (Å).
/// * `r_max`     – maximum r (Å).
/// * `n_bins`    – histogram bins.
pub fn hh_rdf(
    hydrogens: &[[f64; 3]],
    box_l: f64,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = hydrogens.len();
    if n < 2 || n_bins == 0 {
        return (vec![], vec![]);
    }
    let bw = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    for i in 0..n {
        for j in i + 1..n {
            let mut dx = hydrogens[j][0] - hydrogens[i][0];
            let mut dy = hydrogens[j][1] - hydrogens[i][1];
            let mut dz = hydrogens[j][2] - hydrogens[i][2];
            dx -= box_l * (dx / box_l).round();
            dy -= box_l * (dy / box_l).round();
            dz -= box_l * (dz / box_l).round();
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let b = (r / bw) as usize;
                if b < n_bins {
                    hist[b] += 2;
                }
            }
        }
    }
    let rho = n as f64 / (box_l * box_l * box_l);
    let r_centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * bw).collect();
    let g_hh: Vec<f64> = (0..n_bins)
        .map(|b| {
            let r = r_centers[b];
            let vol_shell = 4.0 * PI * r * r * bw;
            let n_ideal = rho * vol_shell * n as f64;
            if n_ideal < 1e-20 {
                0.0
            } else {
                hist[b] as f64 / n_ideal
            }
        })
        .collect();
    (r_centers, g_hh)
}
/// Estimate the expected density of liquid water from a water model by
/// computing the effective molecular volume from LJ sigma.
///
/// The packing fraction of a liquid is approximately 0.64 (random close packing).
/// `rho ≈ 0.64 * (M_water) / (N_A * V_mol)` where `V_mol = (4/3) π (sigma/2)³`.
///
/// Returns density in g/cm³.
pub fn water_model_estimated_density(params: &WaterParams, packing_fraction: f64) -> f64 {
    let sigma_cm = params.sigma_o * 1e-8_f64;
    let r = sigma_cm / 2.0;
    let v_mol_cm3 = (4.0 / 3.0) * PI * r * r * r;
    let m_water_g = 18.015 / 6.022e23_f64;
    packing_fraction * m_water_g / v_mol_cm3
}
/// Full SPC/E polarisation correction including environment dipole.
///
/// The self-polarisation energy per molecule is:
/// `E_pol = (mu - mu_gas)² / (2 * alpha)`
///
/// where `mu_gas` = 1.85 D is the gas-phase dipole, `alpha` = 1.608 Å³ is
/// the molecular polarisability, and `mu` is the model dipole moment.
///
/// Returns energy in kJ/mol.
pub fn spce_polarisation_correction(mu_model_debye: f64) -> f64 {
    let mu_gas = 1.85_f64;
    let alpha_a3 = 1.608_f64;
    let conv = 60.3_f64;
    let delta_mu = mu_model_debye - mu_gas;
    delta_mu * delta_mu / (2.0 * alpha_a3) * conv
}
// --- Small private vector helpers for the analytic SETTLE solver -------------
// (named with a `vec3_` prefix to avoid clashing with the existing `dist`/`dist3`.)

/// Vector subtraction `a - b`.
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Vector addition `a + b`.
#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scalar multiple `s * a`.
#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product `a · b`.
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product `a × b`.
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Euclidean norm `|a|`.
#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
/// Normalize `a`; returns `Err(SettleDegenerate)` if the norm is below `eps`
/// or the result is non-finite.
#[inline]
fn vec3_normalize(a: [f64; 3], eps: f64, what: &str) -> Result<[f64; 3], MdError> {
    let n = vec3_norm(a);
    if !n.is_finite() || n < eps {
        return Err(MdError::SettleDegenerate(format!(
            "cannot normalize {what}: |v| = {n}"
        )));
    }
    let inv = 1.0 / n;
    let v = [a[0] * inv, a[1] * inv, a[2] * inv];
    if !(v[0].is_finite() && v[1].is_finite() && v[2].is_finite()) {
        return Err(MdError::SettleDegenerate(format!(
            "non-finite normalized {what}"
        )));
    }
    Ok(v)
}
/// True if all three components are finite.
#[inline]
fn vec3_finite(a: [f64; 3]) -> bool {
    a[0].is_finite() && a[1].is_finite() && a[2].is_finite()
}
/// Rotate `(x, y, z)` about the body x' axis by angle with the given sin/cos.
#[inline]
fn rot_x(p: [f64; 3], s: f64, c: f64) -> [f64; 3] {
    [p[0], p[1] * c - p[2] * s, p[1] * s + p[2] * c]
}
/// Rotate `(x, y, z)` about the body y' axis by angle with the given sin/cos.
#[inline]
fn rot_y(p: [f64; 3], s: f64, c: f64) -> [f64; 3] {
    [p[0] * c + p[2] * s, p[1], -p[0] * s + p[2] * c]
}
/// Rotate `(x, y, z)` about the body z' axis by angle with the given sin/cos.
#[inline]
fn rot_z(p: [f64; 3], s: f64, c: f64) -> [f64; 3] {
    [p[0] * c - p[1] * s, p[0] * s + p[1] * c, p[2]]
}

/// Analytic Miyamoto-Kollman SETTLE for a single rigid water molecule.
///
/// Closed-form (NON-iterative) constraint solver from
/// Miyamoto & Kollman, *J. Comput. Chem.* **13**(8):952-962 (1992), the same
/// algorithm used by GROMACS and OpenMM. Given the start-of-step `reference`
/// triangle (which fixes the rigid bond/angle geometry and supplies the body
/// frame) and the post-move `unconstrained` positions, it computes the rigid
/// placement of the molecule that (a) satisfies all three distance constraints
/// (O-H1, O-H2, H1-H2) to machine precision *by construction*, (b) preserves the
/// centre of mass of `unconstrained` exactly, and (c) is the least-displacement
/// rigid fit to `unconstrained`. The constrained coordinates are written back
/// into `unconstrained`.
///
/// # Arguments
/// * `reference` — start-of-step positions (define the rigid geometry & frame).
/// * `unconstrained` — post-move positions; overwritten with constrained result.
/// * `geom` — water geometry (provides r_OH, the H-O-H angle, and the M-site).
/// * `m_o` — oxygen mass (must be > 0).
/// * `m_h` — hydrogen mass (must be > 0).
///
/// # Errors
/// Returns [`MdError::SettleDegenerate`] if masses or `r_oh` are non-positive,
/// any input coordinate is non-finite, the reference triangle is collinear /
/// zero-area, or an intermediate quantity becomes degenerate (e.g. the molecule
/// would have to rotate ~90° out of the reference plane).
pub fn settle_positions(
    reference: &WaterMolecule,
    unconstrained: &mut WaterMolecule,
    geom: &WaterGeometry,
    m_o: f64,
    m_h: f64,
) -> Result<(), MdError> {
    // ---- guards ----------------------------------------------------------
    if m_o <= 0.0 || m_o.is_nan() || m_h <= 0.0 || m_h.is_nan() {
        return Err(MdError::SettleDegenerate(format!(
            "non-positive mass: m_o = {m_o}, m_h = {m_h}"
        )));
    }
    if geom.r_oh <= 0.0 || geom.r_oh.is_nan() {
        return Err(MdError::SettleDegenerate(format!(
            "non-positive r_oh = {}",
            geom.r_oh
        )));
    }
    for v in [
        reference.oxygen,
        reference.hydrogen1,
        reference.hydrogen2,
        unconstrained.oxygen,
        unconstrained.hydrogen1,
        unconstrained.hydrogen2,
    ] {
        if !vec3_finite(v) {
            return Err(MdError::SettleDegenerate(
                "non-finite input coordinate".to_string(),
            ));
        }
    }

    let m_t = m_o + 2.0 * m_h;

    // ---- (1) canonical rigid-triangle half-dimensions in the body frame ---
    // O is on +y', the two H atoms at (±rc, -rb, 0); COM at the origin.
    let rc = 0.5 * geom.r_hh(); // half the H-H distance
    let doh_proj = geom.r_oh * geom.half_angle_rad().cos(); // O -> HH-midpoint distance
    let ra = (2.0 * m_h / m_t) * doh_proj; // O distance from COM along +y'
    let rb = (m_o / m_t) * doh_proj; // H distance from COM along -y'
    if !(ra.is_finite() && rb.is_finite() && rc.is_finite()) || rc <= 0.0 || ra <= 0.0 {
        return Err(MdError::SettleDegenerate(
            "degenerate canonical triangle".to_string(),
        ));
    }
    // Canonical body-frame points (COM at origin); this exact rigid triangle is
    // rotated below, so all constraints hold to machine precision regardless of
    // how well the orientation matches the unconstrained points.
    let o_canon = [0.0, ra, 0.0];
    let h1_canon = [-rc, -rb, 0.0];
    let h2_canon = [rc, -rb, 0.0];

    // ---- (2) centres of mass (preserved by SETTLE) ------------------------
    let com = vec3_scale(
        vec3_add(
            vec3_scale(unconstrained.oxygen, m_o),
            vec3_add(
                vec3_scale(unconstrained.hydrogen1, m_h),
                vec3_scale(unconstrained.hydrogen2, m_h),
            ),
        ),
        1.0 / m_t,
    );
    let com_ref = vec3_scale(
        vec3_add(
            vec3_scale(reference.oxygen, m_o),
            vec3_add(
                vec3_scale(reference.hydrogen1, m_h),
                vec3_scale(reference.hydrogen2, m_h),
            ),
        ),
        1.0 / m_t,
    );

    // ---- (3) orthonormal body frame from the REFERENCE triangle -----------
    let a0 = vec3_sub(reference.oxygen, com_ref);
    let b0 = vec3_sub(reference.hydrogen1, com_ref);
    let c0 = vec3_sub(reference.hydrogen2, com_ref);
    // z' = unit normal of the reference plane.
    let normal = vec3_cross(vec3_sub(b0, a0), vec3_sub(c0, a0));
    let axis_z = vec3_normalize(normal, 1e-10, "plane normal")?;
    // y' = in-plane component of the O direction (a0), normalized.
    let a0_perp = vec3_sub(a0, vec3_scale(axis_z, vec3_dot(a0, axis_z)));
    let axis_y = vec3_normalize(a0_perp, 1e-10, "in-plane O axis")?;
    // x' = y' × z'  (right-handed: x = y × z).
    let axis_x = vec3_cross(axis_y, axis_z);
    if !(vec3_finite(axis_x) && vec3_finite(axis_y) && vec3_finite(axis_z)) {
        return Err(MdError::SettleDegenerate(
            "non-finite body axis".to_string(),
        ));
    }

    // Project a lab vector into body coordinates (x', y', z').
    let to_body = |v: [f64; 3]| -> [f64; 3] {
        [
            vec3_dot(v, axis_x),
            vec3_dot(v, axis_y),
            vec3_dot(v, axis_z),
        ]
    };

    // ---- (4) unconstrained positions relative to COM, in body frame -------
    let a1 = to_body(vec3_sub(unconstrained.oxygen, com));
    let b1 = to_body(vec3_sub(unconstrained.hydrogen1, com));
    let c1 = to_body(vec3_sub(unconstrained.hydrogen2, com));

    // ---- (5) Miyamoto-Kollman analytic angles -----------------------------
    // phi: rotation about body x'. O's out-of-plane (z') coordinate is ra*sin(phi).
    let sinphi = (a1[2] / ra).clamp(-1.0, 1.0);
    let cosphi = (1.0 - sinphi * sinphi).sqrt();
    if !cosphi.is_finite() || cosphi < 1e-9 {
        return Err(MdError::SettleDegenerate(format!(
            "cosphi degenerate: cosphi = {cosphi}"
        )));
    }
    // psi: rotation about body y'. (zb - zc) = 2*rc*cosphi*sin(psi).
    let sinpsi = ((b1[2] - c1[2]) / (2.0 * rc * cosphi)).clamp(-1.0, 1.0);
    let cospsi = (1.0 - sinpsi * sinpsi).sqrt();
    if !cospsi.is_finite() {
        return Err(MdError::SettleDegenerate(
            "cospsi non-finite".to_string(),
        ));
    }

    // Apply the out-of-plane rotations R = Ry(psi) * Rx(phi) to the canonical
    // points. Because these are exact rotation matrices, the rigid triangle stays
    // rigid (all pairwise distances preserved to machine precision).
    let o_pp = rot_y(rot_x(o_canon, sinphi, cosphi), sinpsi, cospsi);
    let h1_pp = rot_y(rot_x(h1_canon, sinphi, cosphi), sinpsi, cospsi);
    let h2_pp = rot_y(rot_x(h2_canon, sinphi, cosphi), sinpsi, cospsi);

    // theta: in-plane rotation about z'. Solve the exact 2-D Procrustes / Kabsch
    // problem aligning the (already phi,psi-rotated) rigid triangle to the
    // unconstrained points projected on the body x'y' plane. Because both sets
    // share the COM (origin), the optimal angle is
    //   theta = atan2( Σ m (x_r y_u − y_r x_u),  Σ m (x_r x_u + y_r y_u) ).
    // This closed form equals the Miyamoto-Kollman θ.
    let rigid_xy = [(o_pp[0], o_pp[1]), (h1_pp[0], h1_pp[1]), (h2_pp[0], h2_pp[1])];
    let uncon_xy = [(a1[0], a1[1]), (b1[0], b1[1]), (c1[0], c1[1])];
    let masses = [m_o, m_h, m_h];
    let mut sum_sin = 0.0_f64;
    let mut sum_cos = 0.0_f64;
    for ((&(xr, yr), &(xu, yu)), &mass) in rigid_xy.iter().zip(uncon_xy.iter()).zip(masses.iter()) {
        sum_sin += mass * (xr * yu - yr * xu);
        sum_cos += mass * (xr * xu + yr * yu);
    }
    let theta = sum_sin.atan2(sum_cos);
    let costh = theta.cos();
    let sinth = theta.sin();

    // Apply the in-plane theta rotation about z' to complete R = Rz * Ry * Rx.
    let o_body = rot_z(o_pp, sinth, costh);
    let h1_body = rot_z(h1_pp, sinth, costh);
    let h2_body = rot_z(h2_pp, sinth, costh);

    // Constraints are satisfied by construction (rigid canonical triangle);
    // verify only in debug builds.
    debug_assert!(
        ((vec3_norm(vec3_sub(o_body, h1_body)) - geom.r_oh).abs() < 1e-9)
            && ((vec3_norm(vec3_sub(o_body, h2_body)) - geom.r_oh).abs() < 1e-9)
            && ((vec3_norm(vec3_sub(h1_body, h2_body)) - geom.r_hh()).abs() < 1e-9),
        "SETTLE body-frame triangle violates constraints"
    );

    // ---- (6) transform back to the lab frame and add the COM --------------
    let to_lab = |p: [f64; 3]| -> [f64; 3] {
        vec3_add(
            com,
            vec3_add(
                vec3_scale(axis_x, p[0]),
                vec3_add(vec3_scale(axis_y, p[1]), vec3_scale(axis_z, p[2])),
            ),
        )
    };
    let new_o = to_lab(o_body);
    let new_h1 = to_lab(h1_body);
    let new_h2 = to_lab(h2_body);

    // ---- (8) final non-finite guard ---------------------------------------
    if !(vec3_finite(new_o) && vec3_finite(new_h1) && vec3_finite(new_h2)) {
        return Err(MdError::SettleDegenerate(
            "non-finite output coordinate".to_string(),
        ));
    }
    unconstrained.oxygen = new_o;
    unconstrained.hydrogen1 = new_h1;
    unconstrained.hydrogen2 = new_h2;

    // ---- (7) reconstruct the M-site if present ----------------------------
    if geom.has_m_site() || unconstrained.m_site.is_some() {
        unconstrained.update_m_site(geom);
    }
    Ok(())
}

/// Analytic SETTLE velocity (RATTLE) step for a single rigid water molecule.
///
/// Removes the velocity components along the three bonds (O-H1, H1-H2, H2-O) so
/// that the time derivative of every distance constraint vanishes, i.e. for each
/// bonded pair (i, j): `(v_i − v_j) · (r_i − r_j) = 0` after correction. This is
/// the velocity half of SETTLE (Miyamoto & Kollman, 1992): three Lagrange
/// multipliers `tau` for the three bonds are found by solving a 3×3 linear system
/// (via Cramer's rule), then applied to the velocities. The correction conserves
/// linear (centre-of-mass) momentum exactly.
///
/// # Arguments
/// * `positions` — current (constrained) positions; define the bond directions.
/// * `velocities` — velocities to correct in place.
/// * `m_o` — oxygen mass (> 0).
/// * `m_h` — hydrogen mass (> 0).
///
/// # Errors
/// Returns [`MdError::SettleDegenerate`] for non-positive masses, non-finite
/// inputs, or a singular 3×3 system (degenerate geometry).
pub fn settle_velocities(
    positions: &WaterMolecule,
    velocities: &mut WaterVelocities,
    m_o: f64,
    m_h: f64,
) -> Result<(), MdError> {
    if m_o <= 0.0 || m_o.is_nan() || m_h <= 0.0 || m_h.is_nan() {
        return Err(MdError::SettleDegenerate(format!(
            "non-positive mass: m_o = {m_o}, m_h = {m_h}"
        )));
    }
    let r_a = positions.oxygen;
    let r_b = positions.hydrogen1;
    let r_c = positions.hydrogen2;
    let v_a = velocities.v_oxygen;
    let v_b = velocities.v_hydrogen1;
    let v_c = velocities.v_hydrogen2;
    for v in [r_a, r_b, r_c, v_a, v_b, v_c] {
        if !vec3_finite(v) {
            return Err(MdError::SettleDegenerate(
                "non-finite velocity-step input".to_string(),
            ));
        }
    }

    // Bond vectors (a=O, b=H1, c=H2):
    //   vab = r_a - r_b (O-H1), vbc = r_b - r_c (H1-H2), vca = r_c - r_a (H2-O).
    let vab = vec3_sub(r_a, r_b);
    let vbc = vec3_sub(r_b, r_c);
    let vca = vec3_sub(r_c, r_a);

    // Inverse masses.
    let w_a = 1.0 / m_o;
    let w_b = 1.0 / m_h;
    let w_c = 1.0 / m_h;

    // Corrected velocities are
    //   v_a' = v_a + w_a ( tau_ab vab - tau_ca vca )
    //   v_b' = v_b + w_b ( tau_bc vbc - tau_ab vab )
    //   v_c' = v_c + w_c ( tau_ca vca - tau_bc vbc )
    // Imposing (v_i' - v_j')·v_ij = 0 for each bond gives M·tau = rhs, where:
    //
    //   v_a' - v_b' = (v_a - v_b) + (w_a+w_b) tau_ab vab
    //                 - w_a tau_ca vca - w_b tau_bc vbc
    // dotted with vab ->
    //   row AB:  (w_a+w_b)|vab|^2 tau_ab - w_b (vbc·vab) tau_bc - w_a (vca·vab) tau_ca
    //            = -(v_a - v_b)·vab
    //   row BC: -w_b (vab·vbc) tau_ab + (w_b+w_c)|vbc|^2 tau_bc - w_c (vca·vbc) tau_ca
    //            = -(v_b - v_c)·vbc
    //   row CA: -w_a (vab·vca) tau_ab - w_c (vbc·vca) tau_bc + (w_c+w_a)|vca|^2 tau_ca
    //            = -(v_c - v_a)·vca
    let d_ab_ab = vec3_dot(vab, vab);
    let d_bc_bc = vec3_dot(vbc, vbc);
    let d_ca_ca = vec3_dot(vca, vca);
    let d_bc_ab = vec3_dot(vbc, vab);
    let d_ca_ab = vec3_dot(vca, vab);
    let d_ab_bc = d_bc_ab; // symmetric
    let d_ca_bc = vec3_dot(vca, vbc);
    let d_ab_ca = d_ca_ab; // symmetric
    let d_bc_ca = d_ca_bc; // symmetric

    let m11 = (w_a + w_b) * d_ab_ab;
    let m12 = -w_b * d_bc_ab;
    let m13 = -w_a * d_ca_ab;
    let m21 = -w_b * d_ab_bc;
    let m22 = (w_b + w_c) * d_bc_bc;
    let m23 = -w_c * d_ca_bc;
    let m31 = -w_a * d_ab_ca;
    let m32 = -w_c * d_bc_ca;
    let m33 = (w_c + w_a) * d_ca_ca;

    let rhs1 = -vec3_dot(vec3_sub(v_a, v_b), vab);
    let rhs2 = -vec3_dot(vec3_sub(v_b, v_c), vbc);
    let rhs3 = -vec3_dot(vec3_sub(v_c, v_a), vca);

    // Cramer's rule.
    let det = m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31)
        + m13 * (m21 * m32 - m22 * m31);
    if !det.is_finite() || det.abs() < 1e-14 {
        return Err(MdError::SettleDegenerate(format!(
            "singular velocity-constraint matrix: det = {det}"
        )));
    }
    let det_ab = rhs1 * (m22 * m33 - m23 * m32) - m12 * (rhs2 * m33 - m23 * rhs3)
        + m13 * (rhs2 * m32 - m22 * rhs3);
    let det_bc = m11 * (rhs2 * m33 - m23 * rhs3) - rhs1 * (m21 * m33 - m23 * m31)
        + m13 * (m21 * rhs3 - rhs2 * m31);
    let det_ca = m11 * (m22 * rhs3 - rhs2 * m32) - m12 * (m21 * rhs3 - rhs2 * m31)
        + rhs1 * (m21 * m32 - m22 * m31);
    let tau_ab = det_ab / det;
    let tau_bc = det_bc / det;
    let tau_ca = det_ca / det;
    if !(tau_ab.is_finite() && tau_bc.is_finite() && tau_ca.is_finite()) {
        return Err(MdError::SettleDegenerate(
            "non-finite velocity multiplier".to_string(),
        ));
    }

    let new_va = vec3_add(
        v_a,
        vec3_scale(vec3_sub(vec3_scale(vab, tau_ab), vec3_scale(vca, tau_ca)), w_a),
    );
    let new_vb = vec3_add(
        v_b,
        vec3_scale(vec3_sub(vec3_scale(vbc, tau_bc), vec3_scale(vab, tau_ab)), w_b),
    );
    let new_vc = vec3_add(
        v_c,
        vec3_scale(vec3_sub(vec3_scale(vca, tau_ca), vec3_scale(vbc, tau_bc)), w_c),
    );
    if !(vec3_finite(new_va) && vec3_finite(new_vb) && vec3_finite(new_vc)) {
        return Err(MdError::SettleDegenerate(
            "non-finite corrected velocity".to_string(),
        ));
    }
    velocities.v_oxygen = new_va;
    velocities.v_hydrogen1 = new_vb;
    velocities.v_hydrogen2 = new_vc;
    Ok(())
}

/// Euclidean distance between two 3-D points (Å).
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
