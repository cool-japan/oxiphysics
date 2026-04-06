#![allow(clippy::too_many_arguments)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Microstructure, materials science, and computational chemistry I/O:
//! cluster expansion, EOS, convex hull, solvation, excitations, NMR,
//! phase field, percolation, heat equation, nucleation, molecular graph,
//! fragmentation, ionic conductivity, and final utility helpers.

#![allow(dead_code)]

use super::convenience::{write_f64_dataset, write_vlen_strings};
use super::file::Hdf5File;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Error, Hdf5Result};

// ── Cluster expansion ─────────────────────────────────────────────────────────

/// Write cluster expansion ECIs (effective cluster interactions).
#[allow(dead_code)]
pub fn write_cluster_expansion_eci(
    file: &mut Hdf5File,
    group: &str,
    cluster_sizes: &[usize],
    eci_values: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(cluster_sizes.len(), eci_values.len());
    let cs_f: Vec<f64> = cluster_sizes.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "cluster_sizes", &cs_f)?;
    write_f64_dataset(file, group, "eci_values", eci_values)
}

// ── Structure enumeration ─────────────────────────────────────────────────────

/// Write enumerated structures and their formation energies.
#[allow(dead_code)]
pub fn write_enumerated_structures(
    file: &mut Hdf5File,
    group: &str,
    struct_ids: &[usize],
    concentrations: &[f64],
    formation_energies: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(struct_ids.len(), concentrations.len());
    assert_eq!(struct_ids.len(), formation_energies.len());
    let id_f: Vec<f64> = struct_ids.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "struct_ids", &id_f)?;
    write_f64_dataset(file, group, "concentrations", concentrations)?;
    write_f64_dataset(file, group, "formation_energies", formation_energies)
}

// ── Equation of state ─────────────────────────────────────────────────────────

/// Write Birch-Murnaghan EOS fit results.
#[allow(dead_code)]
pub fn write_eos_fit(
    file: &mut Hdf5File,
    group: &str,
    volumes: &[f64],
    energies: &[f64],
    v0: f64,
    e0: f64,
    b0: f64,
    b0_prime: f64,
) -> Hdf5Result<()> {
    assert_eq!(volumes.len(), energies.len());
    write_f64_dataset(file, group, "volumes", volumes)?;
    write_f64_dataset(file, group, "energies", energies)?;
    file.open_group_mut(group)?
        .set_attr("V0", AttrValue::Float64(v0));
    file.open_group_mut(group)?
        .set_attr("E0", AttrValue::Float64(e0));
    file.open_group_mut(group)?
        .set_attr("B0_GPa", AttrValue::Float64(b0));
    file.open_group_mut(group)?
        .set_attr("B0_prime", AttrValue::Float64(b0_prime));
    Ok(())
}

// ── Convex hull ───────────────────────────────────────────────────────────────

/// Write convex hull stability data.
#[allow(dead_code)]
pub fn write_convex_hull(
    file: &mut Hdf5File,
    group: &str,
    compositions: &[f64],
    hull_energies: &[f64],
    stable_mask: &[bool],
) -> Hdf5Result<()> {
    assert_eq!(compositions.len(), hull_energies.len());
    assert_eq!(compositions.len(), stable_mask.len());
    write_f64_dataset(file, group, "compositions", compositions)?;
    write_f64_dataset(file, group, "hull_energies", hull_energies)?;
    let sm_f: Vec<f64> = stable_mask
        .iter()
        .map(|&b| if b { 1.0 } else { 0.0 })
        .collect();
    write_f64_dataset(file, group, "stable_mask", &sm_f)
}

// ── Solvation energy ──────────────────────────────────────────────────────────

/// Write implicit solvation energies for a set of solvents.
#[allow(dead_code)]
pub fn write_solvation_energies(
    file: &mut Hdf5File,
    group: &str,
    solvents: &[String],
    dg_solv: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(solvents.len(), dg_solv.len());
    write_vlen_strings(file, group, "solvents", solvents)?;
    write_f64_dataset(file, group, "dG_solv_kcal_mol", dg_solv)
}

// ── Excited states ────────────────────────────────────────────────────────────

/// Write TDDFT excitation energies and oscillator strengths.
#[allow(dead_code)]
pub fn write_excitations(
    file: &mut Hdf5File,
    group: &str,
    excitation_energies_ev: &[f64],
    oscillator_strengths: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(excitation_energies_ev.len(), oscillator_strengths.len());
    write_f64_dataset(
        file,
        group,
        "excitation_energies_eV",
        excitation_energies_ev,
    )?;
    write_f64_dataset(file, group, "oscillator_strengths", oscillator_strengths)
}

// ── NMR chemical shifts ───────────────────────────────────────────────────────

/// Write NMR shielding tensors (isotropic + anisotropy).
#[allow(dead_code)]
pub fn write_nmr_shielding(
    file: &mut Hdf5File,
    group: &str,
    sigma_iso: &[f64],
    sigma_aniso: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(sigma_iso.len(), sigma_aniso.len());
    write_f64_dataset(file, group, "sigma_iso", sigma_iso)?;
    write_f64_dataset(file, group, "sigma_aniso", sigma_aniso)
}

// ── Hyperfine coupling ────────────────────────────────────────────────────────

/// Write isotropic hyperfine coupling constants.
#[allow(dead_code)]
pub fn write_hyperfine_coupling(
    file: &mut Hdf5File,
    group: &str,
    a_iso_mhz: &[f64],
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "A_iso_MHz", a_iso_mhz)
}

// ── NICS aromaticity ──────────────────────────────────────────────────────────

/// Write NICS (nucleus-independent chemical shift) values.
#[allow(dead_code)]
pub fn write_nics(
    file: &mut Hdf5File,
    group: &str,
    heights_angstrom: &[f64],
    nics_values: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(heights_angstrom.len(), nics_values.len());
    write_f64_dataset(file, group, "heights_Angstrom", heights_angstrom)?;
    write_f64_dataset(file, group, "NICS_values", nics_values)
}

// ── Charge flow matrix ────────────────────────────────────────────────────────

/// Write atomic charge fluctuation matrix `δQ[i,j]`.
#[allow(dead_code)]
pub fn write_charge_flow_matrix(
    file: &mut Hdf5File,
    group: &str,
    n: usize,
    dq: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(dq.len(), n * n);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "charge_flow", vec![n, n], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "charge_flow")?.write_f64(dq)
}

// ── Fukui functions ───────────────────────────────────────────────────────────

/// Write Fukui reactivity indicators.
#[allow(dead_code)]
pub fn write_fukui_functions(
    file: &mut Hdf5File,
    group: &str,
    f_plus: &[f64],
    f_minus: &[f64],
    f_zero: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(f_plus.len(), f_minus.len());
    assert_eq!(f_plus.len(), f_zero.len());
    write_f64_dataset(file, group, "f_plus", f_plus)?;
    write_f64_dataset(file, group, "f_minus", f_minus)?;
    write_f64_dataset(file, group, "f_zero", f_zero)
}

// ── Final computational chemistry tests ───────────────────────────────────────

// ── Microstructure simulation ─────────────────────────────────────────────────

/// Write phase field order parameter on a 2-D grid.
#[allow(dead_code)]
pub fn write_phase_field_2d(
    file: &mut Hdf5File,
    group: &str,
    step: usize,
    ny: usize,
    nx: usize,
    phi: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(phi.len(), ny * nx);
    let sub = format!("{group}/step_{step:010}");
    file.create_group(&sub)?;
    let _ = file.create_dataset(&sub, "phi", vec![ny, nx], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "phi")?.write_f64(phi)
}

/// Write Cahn-Hilliard chemical potential field.
#[allow(dead_code)]
pub fn write_chemical_potential_field(
    file: &mut Hdf5File,
    group: &str,
    ny: usize,
    nx: usize,
    mu: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(mu.len(), ny * nx);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "chemical_potential",
        vec![ny, nx],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "chemical_potential")?
        .write_f64(mu)
}

// ── Percolation analysis ──────────────────────────────────────────────────────

/// Write percolation cluster statistics.
#[allow(dead_code)]
pub fn write_percolation_stats(
    file: &mut Hdf5File,
    group: &str,
    occupation_probs: &[f64],
    cluster_sizes: &[f64],
    percolates: &[bool],
) -> Hdf5Result<()> {
    let n = occupation_probs.len();
    assert_eq!(cluster_sizes.len(), n);
    assert_eq!(percolates.len(), n);
    write_f64_dataset(file, group, "occupation_probs", occupation_probs)?;
    write_f64_dataset(file, group, "mean_cluster_size", cluster_sizes)?;
    let pf: Vec<f64> = percolates
        .iter()
        .map(|&b| if b { 1.0 } else { 0.0 })
        .collect();
    write_f64_dataset(file, group, "percolates", &pf)
}

// ── Heat equation finite difference ──────────────────────────────────────────

/// Write temperature field from a heat diffusion simulation.
#[allow(dead_code)]
pub fn write_heat_field(
    file: &mut Hdf5File,
    group: &str,
    step: usize,
    ny: usize,
    nx: usize,
    temp: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(temp.len(), ny * nx);
    let sub = format!("{group}/step_{step:010}");
    file.create_group(&sub)?;
    let _ = file.create_dataset(&sub, "temperature", vec![ny, nx], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "temperature")?.write_f64(temp)
}

// ── Reaction-diffusion ────────────────────────────────────────────────────────

/// Write a 2-D concentration field from reaction-diffusion simulation.
#[allow(dead_code)]
pub fn write_concentration_field(
    file: &mut Hdf5File,
    group: &str,
    species: &str,
    step: usize,
    ny: usize,
    nx: usize,
    c: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(c.len(), ny * nx);
    let sub = format!("{group}/{species}/step_{step:010}");
    file.create_group(&sub)?;
    let _ = file.create_dataset(&sub, "concentration", vec![ny, nx], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "concentration")?.write_f64(c)
}

// ── Nucleation theory ─────────────────────────────────────────────────────────

/// Write classical nucleation theory data.
#[allow(dead_code)]
pub fn write_nucleation_data(
    file: &mut Hdf5File,
    group: &str,
    cluster_sizes: &[usize],
    free_energies: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(cluster_sizes.len(), free_energies.len());
    let cs: Vec<f64> = cluster_sizes.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "cluster_sizes", &cs)?;
    write_f64_dataset(file, group, "free_energies", free_energies)
}

// ── Tests: microstructure, percolation ───────────────────────────────────────

// ── Extra roundtrip tests ─────────────────────────────────────────────────────

// ── Parallel coordinates dataset ─────────────────────────────────────────────

/// Write a multi-dimensional dataset suitable for parallel-coordinates plots.
///
/// `data` shape `[n_samples, n_vars]`.
#[allow(dead_code)]
pub fn write_parallel_coordinates(
    file: &mut Hdf5File,
    group: &str,
    var_names: &[String],
    data: &[f64],
    n_samples: usize,
) -> Hdf5Result<()> {
    let n_vars = var_names.len();
    assert_eq!(data.len(), n_samples * n_vars);
    write_vlen_strings(file, group, "variable_names", var_names)?;
    file.create_group(group)?;
    let _ = file.create_dataset(group, "data", vec![n_samples, n_vars], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "data")?.write_f64(data)
}

// ── Simulation convergence log ────────────────────────────────────────────────

/// Write a convergence log (e.g. SCF iterations).
#[allow(dead_code)]
pub fn write_convergence_log(
    file: &mut Hdf5File,
    group: &str,
    iterations: &[usize],
    residuals: &[f64],
    converged: bool,
) -> Hdf5Result<()> {
    assert_eq!(iterations.len(), residuals.len());
    let it_f: Vec<f64> = iterations.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "iterations", &it_f)?;
    write_f64_dataset(file, group, "residuals", residuals)?;
    file.open_group_mut(group)?
        .set_attr("converged", AttrValue::Int32(if converged { 1 } else { 0 }));
    Ok(())
}

// ── Optimization trajectory ───────────────────────────────────────────────────

/// Write geometry optimization step data.
#[allow(dead_code)]
pub fn write_geopt_step(
    file: &mut Hdf5File,
    group: &str,
    step: usize,
    positions: &[[f64; 3]],
    forces: &[[f64; 3]],
    energy: f64,
    max_force: f64,
) -> Hdf5Result<()> {
    let sub = format!("{group}/geopt_{step:06}");
    let pos_f: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    let frc_f: Vec<f64> = forces.iter().flat_map(|f| f.iter().copied()).collect();
    file.create_group(&sub)?;
    let _ = file.create_dataset(
        &sub,
        "positions",
        vec![positions.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(&sub, "positions")?
        .write_f64(&pos_f)?;
    let _ = file.create_dataset(&sub, "forces", vec![forces.len(), 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "forces")?.write_f64(&frc_f)?;
    let g = file.open_group_mut(&sub)?;
    g.set_attr("energy", AttrValue::Float64(energy));
    g.set_attr("max_force", AttrValue::Float64(max_force));
    Ok(())
}

// ── Population analysis ───────────────────────────────────────────────────────

/// Write Mulliken population analysis.
#[allow(dead_code)]
pub fn write_mulliken_charges(
    file: &mut Hdf5File,
    group: &str,
    charges: &[f64],
    spin_densities: Option<&[f64]>,
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "mulliken_charges", charges)?;
    if let Some(sd) = spin_densities {
        assert_eq!(sd.len(), charges.len());
        write_f64_dataset(file, group, "spin_densities", sd)?;
    }
    Ok(())
}

// ── Final integration tests ───────────────────────────────────────────────────

// ── Molecular graph ───────────────────────────────────────────────────────────

/// Write a molecular graph (adjacency list) with bond orders.
#[allow(dead_code)]
pub fn write_molecular_graph(
    file: &mut Hdf5File,
    group: &str,
    bonds: &[(usize, usize, f64)],
) -> Hdf5Result<()> {
    let i_vec: Vec<f64> = bonds.iter().map(|&(i, _, _)| i as f64).collect();
    let j_vec: Vec<f64> = bonds.iter().map(|&(_, j, _)| j as f64).collect();
    let bo_vec: Vec<f64> = bonds.iter().map(|&(_, _, bo)| bo).collect();
    write_f64_dataset(file, group, "bond_i", &i_vec)?;
    write_f64_dataset(file, group, "bond_j", &j_vec)?;
    write_f64_dataset(file, group, "bond_order", &bo_vec)
}

// ── Fragmentation map ─────────────────────────────────────────────────────────

/// Write a fragmentation map: atom-to-fragment assignments.
#[allow(dead_code)]
pub fn write_fragmentation_map(
    file: &mut Hdf5File,
    group: &str,
    atom_fragment: &[usize],
) -> Hdf5Result<()> {
    let af_f: Vec<f64> = atom_fragment.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "atom_fragment", &af_f)
}

/// Write fragment masses.
#[allow(dead_code)]
pub fn write_fragment_masses(file: &mut Hdf5File, group: &str, masses: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "fragment_masses", masses)
}

// ── Final graph / fragmentation tests ────────────────────────────────────────

// ── Bulk transport properties ─────────────────────────────────────────────────

/// Write ionic conductivity vs temperature.
#[allow(dead_code)]
pub fn write_ionic_conductivity(
    file: &mut Hdf5File,
    group: &str,
    temps: &[f64],
    sigma: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(temps.len(), sigma.len());
    write_f64_dataset(file, group, "temperatures", temps)?;
    write_f64_dataset(file, group, "conductivity_S_m", sigma)
}

/// Write heat capacity data.
#[allow(dead_code)]
pub fn write_heat_capacity(
    file: &mut Hdf5File,
    group: &str,
    temps: &[f64],
    cp: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(temps.len(), cp.len());
    write_f64_dataset(file, group, "temperatures", temps)?;
    write_f64_dataset(file, group, "Cp_J_mol_K", cp)
}

// ── Bulk transport tests ──────────────────────────────────────────────────────

// ── Miscellaneous final helpers ───────────────────────────────────────────────

/// Write a named scalar value with units.
#[allow(dead_code)]
pub fn write_scalar_with_units(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    value: f64,
    units: &str,
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, name, &[value])?;
    file.set_dataset_attr(group, name, "units", AttrValue::String(units.to_string()))
}

/// Read a named scalar value.
#[allow(dead_code)]
pub fn read_scalar(file: &Hdf5File, group: &str, name: &str) -> Hdf5Result<f64> {
    let v = file.open_dataset(group, name)?.read_f64()?;
    v.into_iter()
        .next()
        .ok_or_else(|| Hdf5Error::Generic(format!("scalar '{name}' is empty")))
}

/// Return the number of groups directly under the file root.
#[allow(dead_code)]
pub fn root_group_count(file: &Hdf5File) -> usize {
    file.root.groups.len()
}

/// Return `true` if the file root has no datasets or subgroups.
#[allow(dead_code)]
pub fn is_empty_file(file: &Hdf5File) -> bool {
    file.root.groups.is_empty() && file.root.datasets.is_empty()
}
