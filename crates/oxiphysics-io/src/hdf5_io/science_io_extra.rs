// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended simulation I/O: phonon bands, thermal expansion, grain orientation,
//! electron density, magnetisation, spin-orbit coupling, Fermi surface, Wannier,
//! tight-binding, NEB, transition state, GW, optics, thermochemistry, force
//! constants, dipole, polarisability, charge transfer, adsorption, catalysis,
//! radiation shielding.

use super::convenience::{write_f64_dataset, write_vlen_strings};
use super::file::Hdf5File;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Result};

// ── Phonon band structure ─────────────────────────────────────────────────────

/// Write phonon frequencies along a q-path.
///
/// `q_points` shape `[nq, 3]`; `freq_meV` shape `[nq, n_modes]`.
pub fn write_phonon_band_structure(
    file: &mut Hdf5File,
    group: &str,
    nq: usize,
    n_modes: usize,
    q_points: &[f64],
    freq_mev: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(q_points.len(), nq * 3);
    assert_eq!(freq_mev.len(), nq * n_modes);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "q_points", vec![nq, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "q_points")?
        .write_f64(q_points)?;
    let _ = file.create_dataset(
        group,
        "frequencies_meV",
        vec![nq, n_modes],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "frequencies_meV")?
        .write_f64(freq_mev)
}

/// Write phonon group velocities.
pub fn write_phonon_group_velocities(
    file: &mut Hdf5File,
    group: &str,
    nq: usize,
    n_modes: usize,
    vg: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(vg.len(), nq * n_modes * 3);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "group_velocities",
        vec![nq, n_modes, 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "group_velocities")?
        .write_f64(vg)
}

// ── Thermal expansion ─────────────────────────────────────────────────────────

/// Write lattice parameters vs temperature.
pub fn write_lattice_parameters(
    file: &mut Hdf5File,
    group: &str,
    temperatures: &[f64],
    lattice_a: &[f64],
    lattice_b: &[f64],
    lattice_c: &[f64],
) -> Hdf5Result<()> {
    let n = temperatures.len();
    assert_eq!(lattice_a.len(), n);
    assert_eq!(lattice_b.len(), n);
    assert_eq!(lattice_c.len(), n);
    write_f64_dataset(file, group, "temperatures", temperatures)?;
    write_f64_dataset(file, group, "a", lattice_a)?;
    write_f64_dataset(file, group, "b", lattice_b)?;
    write_f64_dataset(file, group, "c", lattice_c)
}

// ── Grain orientation ─────────────────────────────────────────────────────────

/// Write Euler angles (φ1, Φ, φ2) for each grain.
pub fn write_euler_angles(file: &mut Hdf5File, group: &str, euler: &[[f64; 3]]) -> Hdf5Result<()> {
    let flat: Vec<f64> = euler.iter().flat_map(|e| e.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "euler_angles",
        vec![euler.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "euler_angles")?
        .write_f64(&flat)
}

// ── Electron density ──────────────────────────────────────────────────────────

/// Write electronic charge density.
pub fn write_electron_density(
    file: &mut Hdf5File,
    group: &str,
    nx: usize,
    ny: usize,
    nz: usize,
    rho: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(rho.len(), nx * ny * nz);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "rho_e", vec![nx, ny, nz], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "rho_e")?.write_f64(rho)
}

// ── Magnetisation ─────────────────────────────────────────────────────────────

/// Write per-site magnetic moments.
pub fn write_magnetic_moments(
    file: &mut Hdf5File,
    group: &str,
    moments: &[[f64; 3]],
) -> Hdf5Result<()> {
    let flat: Vec<f64> = moments.iter().flat_map(|m| m.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "magnetic_moments",
        vec![moments.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "magnetic_moments")?
        .write_f64(&flat)
}

// ── Spin-orbit coupling ───────────────────────────────────────────────────────

/// Write effective SOC matrix elements (2n×2n complex, stored interleaved real/imag).
pub fn write_soc_matrix(file: &mut Hdf5File, group: &str, n: usize, mat: &[f64]) -> Hdf5Result<()> {
    assert_eq!(mat.len(), (2 * n) * (2 * n) * 2);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "soc_matrix",
        vec![2 * n, 2 * n, 2],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "soc_matrix")?.write_f64(mat)
}

// ── Fermi surface ─────────────────────────────────────────────────────────────

/// Write Fermi surface k-points and weights.
pub fn write_fermi_surface(
    file: &mut Hdf5File,
    group: &str,
    kpoints: &[[f64; 3]],
    weights: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(kpoints.len(), weights.len());
    let flat: Vec<f64> = kpoints.iter().flat_map(|k| k.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "fs_kpoints",
        vec![kpoints.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "fs_kpoints")?
        .write_f64(&flat)?;
    write_f64_dataset(file, group, "fs_weights", weights)
}

// ── Wannier functions ─────────────────────────────────────────────────────────

/// Write Wannier function centres and spreads.
pub fn write_wannier_centres(
    file: &mut Hdf5File,
    group: &str,
    centres: &[[f64; 3]],
    spreads: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(centres.len(), spreads.len());
    let flat: Vec<f64> = centres.iter().flat_map(|c| c.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "wannier_centres",
        vec![centres.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "wannier_centres")?
        .write_f64(&flat)?;
    write_f64_dataset(file, group, "wannier_spreads", spreads)
}

// ── Tight-binding Hamiltonian ─────────────────────────────────────────────────

/// Write a tight-binding Hamiltonian as real and imaginary parts.
pub fn write_tight_binding_hamiltonian(
    file: &mut Hdf5File,
    group: &str,
    n_orb: usize,
    h_real: &[f64],
    h_imag: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(h_real.len(), n_orb * n_orb);
    assert_eq!(h_imag.len(), n_orb * n_orb);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "H_real", vec![n_orb, n_orb], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "H_real")?.write_f64(h_real)?;
    let _ = file.create_dataset(group, "H_imag", vec![n_orb, n_orb], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "H_imag")?.write_f64(h_imag)
}

// ── Nudged elastic band ───────────────────────────────────────────────────────

/// Write NEB path: one image per "image" group.
pub fn write_neb_path(
    file: &mut Hdf5File,
    group: &str,
    images: &[Vec<[f64; 3]>],
    energies: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(images.len(), energies.len());
    for (i, (pos, &e)) in images.iter().zip(energies.iter()).enumerate() {
        let sub = format!("{group}/image_{i:04}");
        let flat: Vec<f64> = pos.iter().flat_map(|p| p.iter().copied()).collect();
        file.create_group(&sub)?;
        let _ = file.create_dataset(&sub, "positions", vec![pos.len(), 3], Hdf5Dtype::Float64);
        file.open_dataset_mut(&sub, "positions")?.write_f64(&flat)?;
        file.open_group_mut(&sub)?
            .set_attr("energy", AttrValue::Float64(e));
    }
    Ok(())
}

/// Read energies from a NEB path.
pub fn read_neb_energies(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<f64>> {
    let g = file.open_group(group)?;
    let mut images: Vec<(String, f64)> = g
        .groups
        .iter()
        .filter(|(k, _)| k.starts_with("image_"))
        .map(|(k, sg)| {
            let e = match sg.attributes.get("energy") {
                Some(AttrValue::Float64(v)) => *v,
                _ => 0.0,
            };
            (k.clone(), e)
        })
        .collect();
    images.sort_by_key(|(k, _)| k.clone());
    Ok(images.into_iter().map(|(_, e)| e).collect())
}

// ── Transition state theory ───────────────────────────────────────────────────

/// Write transition state data: reactant/TS/product energies.
pub fn write_transition_state(
    file: &mut Hdf5File,
    group: &str,
    e_reactant: f64,
    e_ts: f64,
    e_product: f64,
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let g = file.open_group_mut(group)?;
    g.set_attr("E_reactant", AttrValue::Float64(e_reactant));
    g.set_attr("E_ts", AttrValue::Float64(e_ts));
    g.set_attr("E_product", AttrValue::Float64(e_product));
    g.set_attr("barrier", AttrValue::Float64(e_ts - e_reactant));
    Ok(())
}

// ── G0W0 quasiparticle energies ───────────────────────────────────────────────

/// Write G0W0 quasiparticle corrections.
pub fn write_gw_corrections(
    file: &mut Hdf5File,
    group: &str,
    dft_energies: &[f64],
    gw_energies: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(dft_energies.len(), gw_energies.len());
    write_f64_dataset(file, group, "dft_energies", dft_energies)?;
    write_f64_dataset(file, group, "gw_energies", gw_energies)
}

// ── Optical spectrum ──────────────────────────────────────────────────────────

/// Write optical absorption spectrum.
pub fn write_optical_spectrum(
    file: &mut Hdf5File,
    group: &str,
    energies_ev: &[f64],
    epsilon2: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(energies_ev.len(), epsilon2.len());
    write_f64_dataset(file, group, "energies_eV", energies_ev)?;
    write_f64_dataset(file, group, "epsilon2", epsilon2)
}

// ── Thermochemistry ───────────────────────────────────────────────────────────

/// Write thermochemical data: enthalpy, entropy, free energy vs temperature.
pub fn write_thermochemistry(
    file: &mut Hdf5File,
    group: &str,
    temps: &[f64],
    h: &[f64],
    s: &[f64],
    g: &[f64],
) -> Hdf5Result<()> {
    let n = temps.len();
    assert_eq!(h.len(), n);
    assert_eq!(s.len(), n);
    assert_eq!(g.len(), n);
    write_f64_dataset(file, group, "temperatures", temps)?;
    write_f64_dataset(file, group, "enthalpy", h)?;
    write_f64_dataset(file, group, "entropy", s)?;
    write_f64_dataset(file, group, "free_energy", g)
}

// ── Force constants ───────────────────────────────────────────────────────────

/// Write interatomic force constants (IFC) matrix.
///
/// `ifc` shape `[n_atoms, n_atoms, 3, 3]` stored flat.
pub fn write_force_constants(
    file: &mut Hdf5File,
    group: &str,
    n_atoms: usize,
    ifc: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(ifc.len(), n_atoms * n_atoms * 9);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "force_constants",
        vec![n_atoms, n_atoms, 3, 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "force_constants")?
        .write_f64(ifc)
}

// ── Dipole moment ─────────────────────────────────────────────────────────────

/// Write dipole moment trajectory.
pub fn write_dipole_trajectory(
    file: &mut Hdf5File,
    group: &str,
    dipoles: &[[f64; 3]],
) -> Hdf5Result<()> {
    let flat: Vec<f64> = dipoles.iter().flat_map(|d| d.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "dipole_moments",
        vec![dipoles.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "dipole_moments")?
        .write_f64(&flat)
}

// ── Polarisability ────────────────────────────────────────────────────────────

/// Write polarisability tensor trajectory (per step: 3×3).
pub fn write_polarisability(file: &mut Hdf5File, group: &str, alpha: &[f64]) -> Hdf5Result<()> {
    let n_steps = alpha.len() / 9;
    assert_eq!(n_steps * 9, alpha.len());
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "polarisability",
        vec![n_steps, 3, 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "polarisability")?
        .write_f64(alpha)
}

// ── Charge transfer ───────────────────────────────────────────────────────────

/// Write Bader charge decomposition.
pub fn write_bader_charges(file: &mut Hdf5File, group: &str, bader_q: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "bader_charges", bader_q)
}

// ── Adsorption data ───────────────────────────────────────────────────────────

/// Write adsorption energy data for different configurations.
pub fn write_adsorption_energies(
    file: &mut Hdf5File,
    group: &str,
    site_labels: &[&str],
    e_ads: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(site_labels.len(), e_ads.len());
    file.create_group(group)?;
    let g = file.open_group_mut(group)?;
    for (&lbl, &e) in site_labels.iter().zip(e_ads.iter()) {
        g.set_attr(lbl, AttrValue::Float64(e));
    }
    Ok(())
}

// ── Catalysis microkinetics ───────────────────────────────────────────────────

/// Write microkinetic model rate constants.
pub fn write_rate_constants(
    file: &mut Hdf5File,
    group: &str,
    rxn_labels: &[String],
    k_forward: &[f64],
    k_reverse: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(rxn_labels.len(), k_forward.len());
    assert_eq!(rxn_labels.len(), k_reverse.len());
    write_vlen_strings(file, group, "reaction_labels", rxn_labels)?;
    write_f64_dataset(file, group, "k_forward", k_forward)?;
    write_f64_dataset(file, group, "k_reverse", k_reverse)
}

// ── Coverage profile ──────────────────────────────────────────────────────────

/// Write surface coverage vs time.
pub fn write_coverage_profile(
    file: &mut Hdf5File,
    group: &str,
    time: &[f64],
    coverage: &[f64],
    species: &str,
) -> Hdf5Result<()> {
    assert_eq!(time.len(), coverage.len());
    let sub = format!("{group}/{species}");
    write_f64_dataset(file, &sub, "time", time)?;
    write_f64_dataset(file, &sub, "coverage", coverage)
}

// ── Additional tests (quantum chemistry, materials) ───────────────────────────

// ── Radiation shielding data ──────────────────────────────────────────────────

/// Write dose rate profile as a function of depth.
pub fn write_dose_rate_profile(
    file: &mut Hdf5File,
    group: &str,
    depth_cm: &[f64],
    dose_gy_hr: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(depth_cm.len(), dose_gy_hr.len());
    write_f64_dataset(file, group, "depth_cm", depth_cm)?;
    write_f64_dataset(file, group, "dose_Gy_hr", dose_gy_hr)
}

/// Write stopping power data (Bethe formula simulation).
pub fn write_stopping_power(
    file: &mut Hdf5File,
    group: &str,
    energy_mev: &[f64],
    sp: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(energy_mev.len(), sp.len());
    write_f64_dataset(file, group, "energy_MeV", energy_mev)?;
    write_f64_dataset(file, group, "stopping_power", sp)
}

// ── Final misc tests ──────────────────────────────────────────────────────────

// ── Cluster expansion ─────────────────────────────────────────────────────────
