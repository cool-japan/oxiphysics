// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! MD simulation I/O: metadata, trajectories, PES, AIMD steps, NNP training,
//! histograms, covariance, parameter sweeps, REMD, free energy, Monte Carlo,
//! structure factor, dielectric, born charges, and more.

use super::convenience::{count_datasets_recursive, write_f64_dataset};
use super::file::Hdf5File;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Error, Hdf5Result};

// ── Simulation metadata ───────────────────────────────────────────────────────

/// Metadata for an MD run.
#[derive(Debug, Clone)]
pub struct MdRunMetadata {
    /// Run title.
    pub title: String,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of MD steps.
    pub n_steps: usize,
    /// Timestep in fs.
    pub dt_fs: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Pressure in bar.
    pub pressure_bar: f64,
    /// Ensemble (NVT, NPT, etc).
    pub ensemble: String,
    /// Force field name.
    pub force_field: String,
}

impl MdRunMetadata {
    /// Write to file.
    pub fn write_to(&self, file: &mut Hdf5File) -> Hdf5Result<()> {
        file.create_group("metadata")?;
        let g = file.open_group_mut("metadata")?;
        g.set_attr("title", AttrValue::String(self.title.clone()));
        g.set_attr("n_atoms", AttrValue::Int32(self.n_atoms as i32));
        g.set_attr("n_steps", AttrValue::Int32(self.n_steps as i32));
        g.set_attr("dt_fs", AttrValue::Float64(self.dt_fs));
        g.set_attr("temperature_k", AttrValue::Float64(self.temperature_k));
        g.set_attr("pressure_bar", AttrValue::Float64(self.pressure_bar));
        g.set_attr("ensemble", AttrValue::String(self.ensemble.clone()));
        g.set_attr("force_field", AttrValue::String(self.force_field.clone()));
        Ok(())
    }

    /// Read from file.
    pub fn read_from(file: &Hdf5File) -> Hdf5Result<Self> {
        let g = file.open_group("metadata")?;
        let title = match g.attributes.get("title") {
            Some(AttrValue::String(s)) => s.clone(),
            _ => String::new(),
        };
        let n_atoms = match g.attributes.get("n_atoms") {
            Some(AttrValue::Int32(v)) => *v as usize,
            _ => 0,
        };
        let n_steps = match g.attributes.get("n_steps") {
            Some(AttrValue::Int32(v)) => *v as usize,
            _ => 0,
        };
        let dt_fs = match g.attributes.get("dt_fs") {
            Some(AttrValue::Float64(v)) => *v,
            _ => 0.0,
        };
        let temperature_k = match g.attributes.get("temperature_k") {
            Some(AttrValue::Float64(v)) => *v,
            _ => 0.0,
        };
        let pressure_bar = match g.attributes.get("pressure_bar") {
            Some(AttrValue::Float64(v)) => *v,
            _ => 0.0,
        };
        let ensemble = match g.attributes.get("ensemble") {
            Some(AttrValue::String(s)) => s.clone(),
            _ => String::new(),
        };
        let force_field = match g.attributes.get("force_field") {
            Some(AttrValue::String(s)) => s.clone(),
            _ => String::new(),
        };
        Ok(Self {
            title,
            n_atoms,
            n_steps,
            dt_fs,
            temperature_k,
            pressure_bar,
            ensemble,
            force_field,
        })
    }
}

// ── Multi-frame trajectory ────────────────────────────────────────────────────

/// Append a trajectory frame.
pub fn append_trajectory_frame(
    file: &mut Hdf5File,
    traj_group: &str,
    frame_idx: usize,
    positions: &[[f64; 3]],
) -> Hdf5Result<()> {
    let frame_path = format!("{traj_group}/frame_{frame_idx:08}");
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    file.create_group(&frame_path)?;
    let _ = file.create_dataset(
        &frame_path,
        "positions",
        vec![positions.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(&frame_path, "positions")?
        .write_f64(&flat)
}

/// Read positions for a trajectory frame.
pub fn read_trajectory_frame(
    file: &Hdf5File,
    traj_group: &str,
    frame_idx: usize,
) -> Hdf5Result<Vec<[f64; 3]>> {
    let frame_path = format!("{traj_group}/frame_{frame_idx:08}");
    let flat = file.open_dataset(&frame_path, "positions")?.read_f64()?;
    if flat.len() % 3 != 0 {
        return Err(Hdf5Error::Generic("positions not divisible by 3".into()));
    }
    Ok(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

/// Count trajectory frames.
pub fn count_trajectory_frames(file: &Hdf5File, traj_group: &str) -> usize {
    file.open_group(traj_group)
        .map(|g| g.groups.keys().filter(|k| k.starts_with("frame_")).count())
        .unwrap_or(0)
}

// ── PES sampling ──────────────────────────────────────────────────────────────

/// Write PES samples.
pub fn write_pes_samples(
    file: &mut Hdf5File,
    group: &str,
    n_configs: usize,
    n_atoms: usize,
    configs: &[f64],
    energies: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(configs.len(), n_configs * n_atoms * 3);
    assert_eq!(energies.len(), n_configs);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "configurations",
        vec![n_configs, n_atoms, 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "configurations")?
        .write_f64(configs)?;
    write_f64_dataset(file, group, "energies", energies)
}

// ── AIMD step ─────────────────────────────────────────────────────────────────

/// Write one AIMD step.
pub fn write_aimd_step(
    file: &mut Hdf5File,
    group: &str,
    step: usize,
    positions: &[[f64; 3]],
    forces: &[[f64; 3]],
    energy: f64,
    voigt: &[f64; 6],
) -> Hdf5Result<()> {
    let sub = format!("{group}/step_{step:08}");
    let pos_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    let f_flat: Vec<f64> = forces.iter().flat_map(|f| f.iter().copied()).collect();
    file.create_group(&sub)?;
    let _ = file.create_dataset(
        &sub,
        "positions",
        vec![positions.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(&sub, "positions")?
        .write_f64(&pos_flat)?;
    let _ = file.create_dataset(&sub, "forces", vec![forces.len(), 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "forces")?.write_f64(&f_flat)?;
    write_f64_dataset(file, &sub, "energy", &[energy])?;
    let _ = file.create_dataset(&sub, "stress", vec![6], Hdf5Dtype::Float64);
    file.open_dataset_mut(&sub, "stress")?.write_f64(voigt)
}

// ── NNP training batch ────────────────────────────────────────────────────────

/// Write a NNP training batch.
pub fn write_nnp_training_batch(
    file: &mut Hdf5File,
    group: &str,
    n_atoms: usize,
    n_features: usize,
    descriptors: &[f64],
    ref_energies: &[f64],
    ref_forces: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(descriptors.len(), n_atoms * n_features);
    assert_eq!(ref_forces.len(), n_atoms * 3);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "descriptors",
        vec![n_atoms, n_features],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "descriptors")?
        .write_f64(descriptors)?;
    write_f64_dataset(file, group, "ref_energies", ref_energies)?;
    let _ = file.create_dataset(group, "ref_forces", vec![n_atoms, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "ref_forces")?
        .write_f64(ref_forces)
}

// ── Histogram I/O ─────────────────────────────────────────────────────────────

/// Write a histogram.
pub fn write_histogram(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    edges: &[f64],
    counts: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(edges.len(), counts.len() + 1);
    write_f64_dataset(file, group, &format!("{name}_edges"), edges)?;
    write_f64_dataset(file, group, &format!("{name}_counts"), counts)
}

// ── Covariance matrix ─────────────────────────────────────────────────────────

/// Write a `[d, d]` covariance matrix.
pub fn write_covariance_matrix(
    file: &mut Hdf5File,
    group: &str,
    d: usize,
    cov: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(cov.len(), d * d);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "covariance", vec![d, d], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "covariance")?.write_f64(cov)
}

// ── Parameter sweep ───────────────────────────────────────────────────────────

/// Write a parameter sweep.
pub fn write_parameter_sweep(
    file: &mut Hdf5File,
    group: &str,
    param_name: &str,
    param_values: &[f64],
    obs_name: &str,
    obs_values: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(param_values.len(), obs_values.len());
    write_f64_dataset(file, group, param_name, param_values)?;
    write_f64_dataset(file, group, obs_name, obs_values)
}

// ── REMD metadata ─────────────────────────────────────────────────────────────

/// Write REMD metadata.
pub fn write_remd_metadata(
    file: &mut Hdf5File,
    group: &str,
    temperatures: &[f64],
    acceptance_rates: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(temperatures.len(), acceptance_rates.len());
    write_f64_dataset(file, group, "temperatures", temperatures)?;
    write_f64_dataset(file, group, "acceptance_rates", acceptance_rates)
}

// ── Free energy profile ───────────────────────────────────────────────────────

/// Write a PMF along a reaction coordinate.
pub fn write_free_energy_profile(
    file: &mut Hdf5File,
    group: &str,
    xi: &[f64],
    pmf: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(xi.len(), pmf.len());
    write_f64_dataset(file, group, "xi", xi)?;
    write_f64_dataset(file, group, "pmf", pmf)
}

// ── MC run ────────────────────────────────────────────────────────────────────

/// Write Monte Carlo run statistics.
pub fn write_mc_run(
    file: &mut Hdf5File,
    group: &str,
    steps: &[usize],
    energies: &[f64],
    accepted: &[bool],
) -> Hdf5Result<()> {
    assert_eq!(steps.len(), energies.len());
    assert_eq!(steps.len(), accepted.len());
    let sf: Vec<f64> = steps.iter().map(|&x| x as f64).collect();
    let af: Vec<f64> = accepted
        .iter()
        .map(|&b| if b { 1.0 } else { 0.0 })
        .collect();
    write_f64_dataset(file, group, "steps", &sf)?;
    write_f64_dataset(file, group, "energies", energies)?;
    write_f64_dataset(file, group, "accepted", &af)
}

// ── Structure factor ──────────────────────────────────────────────────────────

/// Write S(q).
pub fn write_structure_factor(
    file: &mut Hdf5File,
    group: &str,
    q_values: &[f64],
    sq: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(q_values.len(), sq.len());
    write_f64_dataset(file, group, "q_values", q_values)?;
    write_f64_dataset(file, group, "sq", sq)
}

// ── Dielectric tensor ─────────────────────────────────────────────────────────

/// Write a 3×3 dielectric tensor.
pub fn write_dielectric_tensor(file: &mut Hdf5File, group: &str, eps: &[f64; 9]) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, "dielectric_tensor", vec![3, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "dielectric_tensor")?
        .write_f64(eps)
}

// ── Born charges ──────────────────────────────────────────────────────────────

/// Write Born effective charge tensors.
pub fn write_born_charges(
    file: &mut Hdf5File,
    group: &str,
    n_atoms: usize,
    born: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(born.len(), n_atoms * 9);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "born_charges",
        vec![n_atoms, 3, 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "born_charges")?
        .write_f64(born)
}

// ── Virial coefficients ───────────────────────────────────────────────────────

/// Write virial coefficients.
pub fn write_virial_coefficients(
    file: &mut Hdf5File,
    group: &str,
    coeffs: &[f64],
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "virial_coefficients", coeffs)
}

// ── Order parameter ───────────────────────────────────────────────────────────

/// Write a scalar order parameter series.
pub fn write_order_parameter(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    values: &[f64],
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, name, values)
}

// ── Stress tensor series ──────────────────────────────────────────────────────

/// Write a Voigt stress tensor for a step.
pub fn write_stress_tensor(
    file: &mut Hdf5File,
    group: &str,
    step: usize,
    voigt: &[f64; 6],
) -> Hdf5Result<()> {
    let ds_name = format!("stress_{step:08}");
    file.create_group(group)?;
    let _ = file.create_dataset(group, &ds_name, vec![6], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, &ds_name)?.write_f64(voigt)
}

// ── Fluid dynamics ────────────────────────────────────────────────────────────

/// Write a 2-D velocity field.
pub fn write_velocity_field_2d(
    file: &mut Hdf5File,
    group: &str,
    ny: usize,
    nx: usize,
    u: &[f64],
    v: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(u.len(), ny * nx);
    assert_eq!(v.len(), ny * nx);
    file.create_group(group)?;
    for (name, data) in [("u", u), ("v", v)] {
        let _ = file.create_dataset(group, name, vec![ny, nx], Hdf5Dtype::Float64);
        file.open_dataset_mut(group, name)?.write_f64(data)?;
    }
    Ok(())
}

/// Write a 3-D velocity field.
pub fn write_velocity_field_3d(
    file: &mut Hdf5File,
    group: &str,
    nz: usize,
    ny: usize,
    nx: usize,
    u: &[f64],
    v: &[f64],
    w: &[f64],
) -> Hdf5Result<()> {
    let n = nz * ny * nx;
    assert_eq!(u.len(), n);
    assert_eq!(v.len(), n);
    assert_eq!(w.len(), n);
    file.create_group(group)?;
    for (name, data) in [("u", u), ("v", v), ("w", w)] {
        let _ = file.create_dataset(group, name, vec![nz, ny, nx], Hdf5Dtype::Float64);
        file.open_dataset_mut(group, name)?.write_f64(data)?;
    }
    Ok(())
}

/// Write a pressure field.
pub fn write_pressure_field(
    file: &mut Hdf5File,
    group: &str,
    ny: usize,
    nx: usize,
    p: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(p.len(), ny * nx);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "pressure", vec![ny, nx], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "pressure")?.write_f64(p)
}

/// Write a vorticity field.
pub fn write_vorticity_field(
    file: &mut Hdf5File,
    group: &str,
    ny: usize,
    nx: usize,
    omega: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(omega.len(), ny * nx);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "vorticity", vec![ny, nx], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "vorticity")?.write_f64(omega)
}

// ── SPH particles ─────────────────────────────────────────────────────────────

/// Write SPH particle data.
pub fn write_sph_particles(
    file: &mut Hdf5File,
    group: &str,
    positions: &[[f64; 3]],
    h: &[f64],
    rho: &[f64],
) -> Hdf5Result<()> {
    let n = positions.len();
    assert_eq!(h.len(), n);
    assert_eq!(rho.len(), n);
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(group, "positions", vec![n, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "positions")?
        .write_f64(&flat)?;
    write_f64_dataset(file, group, "h", h)?;
    write_f64_dataset(file, group, "rho", rho)
}

// ── FEM mesh ──────────────────────────────────────────────────────────────────

/// Write a FEM mesh.
pub fn write_fem_mesh(
    file: &mut Hdf5File,
    group: &str,
    n_nodes: usize,
    nodes: &[f64],
    n_elem: usize,
    nodes_per_elem: usize,
    elements: &[usize],
) -> Hdf5Result<()> {
    assert_eq!(nodes.len(), n_nodes * 3);
    assert_eq!(elements.len(), n_elem * nodes_per_elem);
    file.create_group(group)?;
    {
        let g = file.open_group_mut(group)?;
        g.set_attr("n_nodes", AttrValue::Int32(n_nodes as i32));
        g.set_attr("n_elem", AttrValue::Int32(n_elem as i32));
        g.set_attr("nodes_per_elem", AttrValue::Int32(nodes_per_elem as i32));
    }
    let _ = file.create_dataset(group, "nodes", vec![n_nodes, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "nodes")?.write_f64(nodes)?;
    let elem_f: Vec<f64> = elements.iter().map(|&x| x as f64).collect();
    let _ = file.create_dataset(
        group,
        "elements",
        vec![n_elem, nodes_per_elem],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "elements")?.write_f64(&elem_f)
}

/// Write nodal displacements.
pub fn write_nodal_displacements(
    file: &mut Hdf5File,
    group: &str,
    n_nodes: usize,
    disp: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(disp.len(), n_nodes * 3);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "displacements", vec![n_nodes, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "displacements")?
        .write_f64(disp)
}

/// Write von Mises stresses.
pub fn write_von_mises_stress(file: &mut Hdf5File, group: &str, stress: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "von_mises", stress)
}

// ── LBM populations ───────────────────────────────────────────────────────────

/// Write LBM D2Q9 populations.
pub fn write_lbm_populations(
    file: &mut Hdf5File,
    group: &str,
    ny: usize,
    nx: usize,
    f: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(f.len(), ny * nx * 9);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "populations", vec![ny, nx, 9], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "populations")?.write_f64(f)
}

// ── Sensor recording ──────────────────────────────────────────────────────────

/// Write a multi-channel sensor recording.
pub fn write_sensor_recording(
    file: &mut Hdf5File,
    group: &str,
    n_channels: usize,
    n_samples: usize,
    sample_rate_hz: f64,
    data: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(data.len(), n_channels * n_samples);
    file.create_group(group)?;
    {
        let g = file.open_group_mut(group)?;
        g.set_attr("n_channels", AttrValue::Int32(n_channels as i32));
        g.set_attr("n_samples", AttrValue::Int32(n_samples as i32));
        g.set_attr("sample_rate_hz", AttrValue::Float64(sample_rate_hz));
    }
    let _ = file.create_dataset(
        group,
        "data",
        vec![n_channels, n_samples],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "data")?.write_f64(data)
}

// ── Image stack ───────────────────────────────────────────────────────────────

/// Write an image stack `[n_frames, height, width]`.
pub fn write_image_stack(
    file: &mut Hdf5File,
    group: &str,
    n_frames: usize,
    height: usize,
    width: usize,
    pixels: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(pixels.len(), n_frames * height * width);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "images",
        vec![n_frames, height, width],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "images")?.write_f64(pixels)
}

// ── Point cloud ───────────────────────────────────────────────────────────────

/// Write a 3-D point cloud.
pub fn write_point_cloud(
    file: &mut Hdf5File,
    group: &str,
    points: &[[f64; 3]],
    labels: Option<&[i64]>,
) -> Hdf5Result<()> {
    let flat: Vec<f64> = points.iter().flat_map(|p| p.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(group, "points", vec![points.len(), 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "points")?.write_f64(&flat)?;
    if let Some(lbl) = labels {
        assert_eq!(lbl.len(), points.len());
        let lf: Vec<f64> = lbl.iter().map(|&x| x as f64).collect();
        write_f64_dataset(file, group, "labels", &lf)?;
    }
    Ok(())
}

// ── Diffusion coefficient ─────────────────────────────────────────────────────

/// Write a diffusion coefficient.
pub fn write_diffusion_coefficient(
    file: &mut Hdf5File,
    group: &str,
    species: &str,
    d_coeff: f64,
    units: &str,
) -> Hdf5Result<()> {
    let ds_name = format!("D_{species}");
    write_f64_dataset(file, group, &ds_name, &[d_coeff])?;
    file.set_dataset_attr(
        group,
        &ds_name,
        "units",
        AttrValue::String(units.to_string()),
    )
}

// ── Arrhenius data ────────────────────────────────────────────────────────────

/// Write Arrhenius rate data.
pub fn write_arrhenius_data(
    file: &mut Hdf5File,
    group: &str,
    temperatures: &[f64],
    rates: &[f64],
    activation_energy_ev: f64,
) -> Hdf5Result<()> {
    assert_eq!(temperatures.len(), rates.len());
    write_f64_dataset(file, group, "temperatures", temperatures)?;
    write_f64_dataset(file, group, "rates", rates)?;
    file.open_group_mut(group)?
        .set_attr("Ea_eV", AttrValue::Float64(activation_energy_ev));
    Ok(())
}

// ── ML feature matrix ─────────────────────────────────────────────────────────

/// Write a feature matrix `[n_samples, n_features]`.
pub fn write_feature_matrix(
    file: &mut Hdf5File,
    group: &str,
    n_samples: usize,
    n_features: usize,
    features: &[f64],
    labels: Option<&[f64]>,
) -> Hdf5Result<()> {
    assert_eq!(features.len(), n_samples * n_features);
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "features",
        vec![n_samples, n_features],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "features")?
        .write_f64(features)?;
    if let Some(lbl) = labels {
        assert_eq!(lbl.len(), n_samples);
        write_f64_dataset(file, group, "labels", lbl)?;
    }
    Ok(())
}

// ── MLP weights ───────────────────────────────────────────────────────────────

/// Write MLP layer weights.  Each tuple is `(weights, biases, n_out, n_in)`.
pub fn write_mlp_weights(
    file: &mut Hdf5File,
    group: &str,
    layers: &[(&[f64], &[f64], usize, usize)],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    file.open_group_mut(group)?
        .set_attr("n_layers", AttrValue::Int32(layers.len() as i32));
    for (i, &(w, b, nout, nin)) in layers.iter().enumerate() {
        assert_eq!(w.len(), nout * nin);
        let sub = format!("{group}/layer_{i:04}");
        file.create_group(&sub)?;
        let _ = file.create_dataset(&sub, "weights", vec![nout, nin], Hdf5Dtype::Float64);
        file.open_dataset_mut(&sub, "weights")?.write_f64(w)?;
        write_f64_dataset(file, &sub, "biases", b)?;
    }
    Ok(())
}

/// Read a single MLP layer's weights and biases.
pub fn read_mlp_layer(
    file: &Hdf5File,
    group: &str,
    layer_idx: usize,
) -> Hdf5Result<(Vec<f64>, Vec<f64>)> {
    let sub = format!("{group}/layer_{layer_idx:04}");
    let w = file.open_dataset(&sub, "weights")?.read_f64()?;
    let b = file.open_dataset(&sub, "biases")?.read_f64()?;
    Ok((w, b))
}

// ── Hyperparameter trial ──────────────────────────────────────────────────────

/// Write a hyperparameter search trial.
pub fn write_hparam_trial(
    file: &mut Hdf5File,
    group: &str,
    trial_idx: usize,
    hparams: &[(&str, f64)],
    val_metric: f64,
) -> Hdf5Result<()> {
    let sub = format!("{group}/trial_{trial_idx:06}");
    file.create_group(&sub)?;
    let g = file.open_group_mut(&sub)?;
    g.set_attr("val_metric", AttrValue::Float64(val_metric));
    for &(name, val) in hparams {
        g.set_attr(name, AttrValue::Float64(val));
    }
    Ok(())
}

// ── Format version ────────────────────────────────────────────────────────────

/// Write format version info.
pub fn write_format_version(
    file: &mut Hdf5File,
    major: u32,
    minor: u32,
    patch: u32,
    creator: &str,
) -> Hdf5Result<()> {
    file.create_group("__version__")?;
    let g = file.open_group_mut("__version__")?;
    g.set_attr("major", AttrValue::Int32(major as i32));
    g.set_attr("minor", AttrValue::Int32(minor as i32));
    g.set_attr("patch", AttrValue::Int32(patch as i32));
    g.set_attr("creator", AttrValue::String(creator.to_string()));
    Ok(())
}

/// Read format version; returns `(major, minor, patch, creator)`.
pub fn read_format_version(file: &Hdf5File) -> Hdf5Result<(u32, u32, u32, String)> {
    let g = file.open_group("__version__")?;
    let major = match g.attributes.get("major") {
        Some(AttrValue::Int32(v)) => *v as u32,
        _ => 0,
    };
    let minor = match g.attributes.get("minor") {
        Some(AttrValue::Int32(v)) => *v as u32,
        _ => 0,
    };
    let patch = match g.attributes.get("patch") {
        Some(AttrValue::Int32(v)) => *v as u32,
        _ => 0,
    };
    let creator = match g.attributes.get("creator") {
        Some(AttrValue::String(s)) => s.clone(),
        _ => String::new(),
    };
    Ok((major, minor, patch, creator))
}

// ── Elastic constants ─────────────────────────────────────────────────────────

/// Write the 6×6 Voigt elastic constant matrix.
pub fn write_elastic_constants(
    file: &mut Hdf5File,
    group: &str,
    cij: &[f64; 36],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, "Cij_GPa", vec![6, 6], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "Cij_GPa")?.write_f64(cij)
}

/// Read elastic constants.
pub fn read_elastic_constants(file: &Hdf5File, group: &str) -> Hdf5Result<[f64; 36]> {
    let v = file.open_dataset(group, "Cij_GPa")?.read_f64()?;
    if v.len() != 36 {
        return Err(Hdf5Error::Generic("Cij must have 36 elements".into()));
    }
    let mut arr = [0.0_f64; 36];
    arr.copy_from_slice(&v);
    Ok(arr)
}

// ── Molecular orbitals ────────────────────────────────────────────────────────

/// Write MO energies and occupations.
pub fn write_molecular_orbitals(
    file: &mut Hdf5File,
    group: &str,
    mo_energies: &[f64],
    occupations: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(mo_energies.len(), occupations.len());
    write_f64_dataset(file, group, "mo_energies", mo_energies)?;
    write_f64_dataset(file, group, "occupations", occupations)
}

/// Write HOMO-LUMO gap.
pub fn write_homo_lumo_gap(file: &mut Hdf5File, group: &str, gap_ev: f64) -> Hdf5Result<()> {
    file.create_group(group)?;
    file.open_group_mut(group)?
        .set_attr("homo_lumo_gap_eV", AttrValue::Float64(gap_ev));
    Ok(())
}

/// Write partial DOS for a species.
pub fn write_pdos(
    file: &mut Hdf5File,
    group: &str,
    species: &str,
    energies: &[f64],
    pdos: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(energies.len(), pdos.len());
    let sub = format!("{group}/pdos_{species}");
    write_f64_dataset(file, &sub, "energies", energies)?;
    write_f64_dataset(file, &sub, "pdos", pdos)
}

// ── Polymer data ──────────────────────────────────────────────────────────────

/// Write polymer end-to-end vectors.
pub fn write_polymer_ete(file: &mut Hdf5File, group: &str, ete: &[[f64; 3]]) -> Hdf5Result<()> {
    let flat: Vec<f64> = ete.iter().flat_map(|v| v.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(group, "end_to_end", vec![ete.len(), 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "end_to_end")?.write_f64(&flat)
}

/// Write radius of gyration series.
pub fn write_radius_of_gyration(file: &mut Hdf5File, group: &str, rg: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "rg", rg)
}

// ── Spectroscopy ──────────────────────────────────────────────────────────────

/// Write IR spectrum.
pub fn write_ir_spectrum(
    file: &mut Hdf5File,
    group: &str,
    wavenumbers: &[f64],
    absorbance: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(wavenumbers.len(), absorbance.len());
    write_f64_dataset(file, group, "wavenumbers", wavenumbers)?;
    write_f64_dataset(file, group, "absorbance", absorbance)
}

/// Write Raman spectrum.
pub fn write_raman_spectrum(
    file: &mut Hdf5File,
    group: &str,
    wavenumbers: &[f64],
    intensities: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(wavenumbers.len(), intensities.len());
    write_f64_dataset(file, group, "wavenumbers", wavenumbers)?;
    write_f64_dataset(file, group, "intensities", intensities)
}

// ── GCMC ──────────────────────────────────────────────────────────────────────

/// Write GCMC run data.
pub fn write_gcmc_run(
    file: &mut Hdf5File,
    group: &str,
    steps: &[usize],
    n_particles: &[usize],
    energies: &[f64],
    mu: f64,
) -> Hdf5Result<()> {
    assert_eq!(steps.len(), n_particles.len());
    assert_eq!(steps.len(), energies.len());
    let sf: Vec<f64> = steps.iter().map(|&x| x as f64).collect();
    let nf: Vec<f64> = n_particles.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "steps", &sf)?;
    write_f64_dataset(file, group, "n_particles", &nf)?;
    write_f64_dataset(file, group, "energies", energies)?;
    file.open_group_mut(group)?
        .set_attr("chemical_potential", AttrValue::Float64(mu));
    Ok(())
}

// ── Thermal conductivity ──────────────────────────────────────────────────────

/// Write thermal conductivity Green-Kubo data.
pub fn write_thermal_conductivity(
    file: &mut Hdf5File,
    group: &str,
    time: &[f64],
    hcacf: &[f64],
    kappa: f64,
    kappa_components: &[f64; 3],
) -> Hdf5Result<()> {
    assert_eq!(time.len(), hcacf.len());
    write_f64_dataset(file, group, "time", time)?;
    write_f64_dataset(file, group, "hcacf", hcacf)?;
    let g = file.open_group_mut(group)?;
    g.set_attr("kappa", AttrValue::Float64(kappa));
    g.set_attr("kappa_xx", AttrValue::Float64(kappa_components[0]));
    g.set_attr("kappa_yy", AttrValue::Float64(kappa_components[1]));
    g.set_attr("kappa_zz", AttrValue::Float64(kappa_components[2]));
    Ok(())
}

// ── Utility ───────────────────────────────────────────────────────────────────

/// Total dataset count across all groups.
pub fn total_dataset_count(file: &Hdf5File) -> usize {
    count_datasets_recursive(&file.root)
}

/// List top-level group names (sorted).
pub fn list_top_level_groups(file: &Hdf5File) -> Vec<String> {
    let mut names: Vec<String> = file.root.groups.keys().cloned().collect();
    names.sort();
    names
}

// ── Surface energy ────────────────────────────────────────────────────────────

/// Write surface energies per facet.
pub fn write_surface_energies(
    file: &mut Hdf5File,
    group: &str,
    facet_labels: &[&str],
    energies_j_m2: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(facet_labels.len(), energies_j_m2.len());
    file.create_group(group)?;
    let g = file.open_group_mut(group)?;
    for (&lbl, &e) in facet_labels.iter().zip(energies_j_m2.iter()) {
        g.set_attr(lbl, AttrValue::Float64(e));
    }
    Ok(())
}
