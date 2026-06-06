//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{WaterGeometry, WaterModelSummary, WaterModelType, WaterMolecule, WaterParams};

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
/// SETTLE rigid-body constraint algorithm stub.
///
/// Projects the hydrogen positions back so that:
/// - Each O-H bond length equals the target from `geom`.
/// - The H-H distance is consistent with the bond angle.
pub fn rigid_water_settle(mol: &mut WaterMolecule, _dt: f64) {
    rigid_water_settle_geom(mol, &WaterGeometry::tip3p());
}
/// SETTLE with explicit geometry parameters.
///
/// Iteratively projects bond distances until convergence.
pub fn rigid_water_settle_geom(mol: &mut WaterMolecule, geom: &WaterGeometry) {
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
/// Euclidean distance between two 3-D points (Å).
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
