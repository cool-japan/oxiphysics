// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Domain-specific simulation I/O: phase-space, thermodynamics,
//! neighbours, RDF, MSD, bonds, charge density, bands, sparse matrices.

use super::convenience::write_f64_dataset;
use super::file::Hdf5File;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Error, Hdf5Result};

// ── Phase-space I/O helpers ───────────────────────────────────────────────────

/// Write per-atom velocities (one `[f64;3]` per atom) to a group.
///
/// The dataset `"velocities"` is created with shape `[n_atoms, 3]`.
pub fn write_velocities(
    file: &mut Hdf5File,
    group: &str,
    velocities: &[[f64; 3]],
) -> Hdf5Result<()> {
    let flat: Vec<f64> = velocities.iter().flat_map(|v| v.iter().copied()).collect();
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "velocities",
        vec![velocities.len(), 3],
        Hdf5Dtype::Float64,
    );
    file.open_dataset_mut(group, "velocities")?.write_f64(&flat)
}

/// Read velocities previously stored by [`write_velocities`].
pub fn read_velocities(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<[f64; 3]>> {
    let flat = file.open_dataset(group, "velocities")?.read_f64()?;
    if flat.len() % 3 != 0 {
        return Err(Hdf5Error::Generic(
            "velocity data not divisible by 3".into(),
        ));
    }
    Ok(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

/// Write per-atom masses as a 1-D dataset.
pub fn write_masses(file: &mut Hdf5File, group: &str, masses: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "masses", masses)
}

/// Read per-atom masses.
pub fn read_masses(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<f64>> {
    file.open_dataset(group, "masses")?.read_f64()
}

/// Write box vectors (3×3 matrix, row-major).
pub fn write_box_vectors(file: &mut Hdf5File, group: &str, box_vecs: &[f64; 9]) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, "box_vectors", vec![3, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "box_vectors")?
        .write_f64(box_vecs)
}

/// Read box vectors.
pub fn read_box_vectors(file: &Hdf5File, group: &str) -> Hdf5Result<[f64; 9]> {
    let v = file.open_dataset(group, "box_vectors")?.read_f64()?;
    if v.len() != 9 {
        return Err(Hdf5Error::Generic(
            "box_vectors must have 9 elements".into(),
        ));
    }
    let mut arr = [0.0_f64; 9];
    arr.copy_from_slice(&v);
    Ok(arr)
}

/// Write per-atom charges.
pub fn write_charges(file: &mut Hdf5File, group: &str, charges: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "charges", charges)
}

/// Read per-atom charges.
pub fn read_charges(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<f64>> {
    file.open_dataset(group, "charges")?.read_f64()
}

// ── Thermodynamic observables ─────────────────────────────────────────────────

/// Write a potential energy time series.
pub fn write_potential_energy_series(
    file: &mut Hdf5File,
    group: &str,
    energies: &[f64],
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "potential_energy", energies)
}

/// Write a temperature time series.
pub fn write_temperature_series(file: &mut Hdf5File, group: &str, temps: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "temperature", temps)
}

/// Write a pressure time series.
pub fn write_pressure_series(
    file: &mut Hdf5File,
    group: &str,
    pressures: &[f64],
) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "pressure", pressures)
}

/// Read a named scalar time series from a group.
pub fn read_scalar_series(file: &Hdf5File, group: &str, name: &str) -> Hdf5Result<Vec<f64>> {
    file.open_dataset(group, name)?.read_f64()
}

// ── Neighbour-list storage ────────────────────────────────────────────────────

/// Sparse neighbour list in CSR format.
#[derive(Debug, Clone)]
pub struct NeighbourList {
    /// Number of atoms.
    pub n_atoms: usize,
    /// CSR row pointer (length `n_atoms + 1`).
    pub row_ptr: Vec<usize>,
    /// Flat column indices.
    pub col_idx: Vec<usize>,
    /// Interatomic distances.
    pub distances: Vec<f64>,
}

impl NeighbourList {
    /// Create an empty neighbour list.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            row_ptr: vec![0; n_atoms + 1],
            col_idx: Vec::new(),
            distances: Vec::new(),
        }
    }

    /// Neighbour indices for atom `i`.
    pub fn neighbours(&self, i: usize) -> &[usize] {
        &self.col_idx[self.row_ptr[i]..self.row_ptr[i + 1]]
    }

    /// Distances for atom `i`'s neighbours.
    pub fn neighbour_distances(&self, i: usize) -> &[f64] {
        &self.distances[self.row_ptr[i]..self.row_ptr[i + 1]]
    }
}

/// Write a [`NeighbourList`].
pub fn write_neighbour_list(
    file: &mut Hdf5File,
    group: &str,
    nl: &NeighbourList,
) -> Hdf5Result<()> {
    let row_ptr_f: Vec<f64> = nl.row_ptr.iter().map(|&x| x as f64).collect();
    let col_idx_f: Vec<f64> = nl.col_idx.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "row_ptr", &row_ptr_f)?;
    write_f64_dataset(file, group, "col_idx", &col_idx_f)?;
    write_f64_dataset(file, group, "distances", &nl.distances)
}

/// Read a [`NeighbourList`].
pub fn read_neighbour_list(file: &Hdf5File, group: &str) -> Hdf5Result<NeighbourList> {
    let row_ptr: Vec<usize> = file
        .open_dataset(group, "row_ptr")?
        .read_f64()?
        .iter()
        .map(|&x| x as usize)
        .collect();
    let col_idx: Vec<usize> = file
        .open_dataset(group, "col_idx")?
        .read_f64()?
        .iter()
        .map(|&x| x as usize)
        .collect();
    let distances = file.open_dataset(group, "distances")?.read_f64()?;
    let n_atoms = row_ptr.len().saturating_sub(1);
    Ok(NeighbourList {
        n_atoms,
        row_ptr,
        col_idx,
        distances,
    })
}

// ── Radial distribution function ──────────────────────────────────────────────

/// Write a radial distribution function.
pub fn write_rdf(file: &mut Hdf5File, group: &str, r_bins: &[f64], gr: &[f64]) -> Hdf5Result<()> {
    assert_eq!(r_bins.len(), gr.len());
    write_f64_dataset(file, group, "r_bins", r_bins)?;
    write_f64_dataset(file, group, "gr", gr)
}

/// Read a radial distribution function; returns `(r_bins, gr)`.
pub fn read_rdf(file: &Hdf5File, group: &str) -> Hdf5Result<(Vec<f64>, Vec<f64>)> {
    let r = file.open_dataset(group, "r_bins")?.read_f64()?;
    let g = file.open_dataset(group, "gr")?.read_f64()?;
    Ok((r, g))
}

// ── Mean squared displacement ─────────────────────────────────────────────────

/// Write mean-squared displacement data.
pub fn write_msd(
    file: &mut Hdf5File,
    group: &str,
    time_lags: &[f64],
    msd: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(time_lags.len(), msd.len());
    write_f64_dataset(file, group, "time_lags", time_lags)?;
    write_f64_dataset(file, group, "msd", msd)
}

/// Read MSD data; returns `(time_lags, msd)`.
pub fn read_msd(file: &Hdf5File, group: &str) -> Hdf5Result<(Vec<f64>, Vec<f64>)> {
    let t = file.open_dataset(group, "time_lags")?.read_f64()?;
    let m = file.open_dataset(group, "msd")?.read_f64()?;
    Ok((t, m))
}

// ── Velocity autocorrelation ──────────────────────────────────────────────────

/// Write velocity autocorrelation function.
pub fn write_vacf(
    file: &mut Hdf5File,
    group: &str,
    time_lags: &[f64],
    vacf: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(time_lags.len(), vacf.len());
    write_f64_dataset(file, group, "time_lags", time_lags)?;
    write_f64_dataset(file, group, "vacf", vacf)
}

// ── Power spectrum ────────────────────────────────────────────────────────────

/// Write a phonon density-of-states.
pub fn write_power_spectrum(
    file: &mut Hdf5File,
    group: &str,
    freqs: &[f64],
    dos: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(freqs.len(), dos.len());
    write_f64_dataset(file, group, "frequencies", freqs)?;
    write_f64_dataset(file, group, "dos", dos)
}

// ── Bond topology ─────────────────────────────────────────────────────────────

/// Bond pair `(atom_i, atom_j)`.
pub type BondPair = (usize, usize);

/// Write a bond list.
pub fn write_bonds(file: &mut Hdf5File, group: &str, bonds: &[BondPair]) -> Hdf5Result<()> {
    let flat: Vec<f64> = bonds
        .iter()
        .flat_map(|&(a, b)| [a as f64, b as f64])
        .collect();
    file.create_group(group)?;
    let _ = file.create_dataset(group, "bonds", vec![bonds.len(), 2], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "bonds")?.write_f64(&flat)
}

/// Read a bond list.
pub fn read_bonds(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<BondPair>> {
    let flat = file.open_dataset(group, "bonds")?.read_f64()?;
    if flat.len() % 2 != 0 {
        return Err(Hdf5Error::Generic("bonds not divisible by 2".into()));
    }
    Ok(flat
        .chunks_exact(2)
        .map(|c| (c[0] as usize, c[1] as usize))
        .collect())
}

// ── Charge density grid ───────────────────────────────────────────────────────

/// Write a 3-D charge density grid.
pub fn write_charge_density(
    file: &mut Hdf5File,
    group: &str,
    dims: [usize; 3],
    data: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(data.len(), dims[0] * dims[1] * dims[2]);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "charge_density", dims.to_vec(), Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "charge_density")?
        .write_f64(data)
}

/// Read a charge density grid; returns `(dims, data)`.
pub fn read_charge_density(file: &Hdf5File, group: &str) -> Hdf5Result<(Vec<usize>, Vec<f64>)> {
    let ds = file.open_dataset(group, "charge_density")?;
    let dims = ds.shape.clone();
    let data = ds.read_f64()?;
    Ok((dims, data))
}

// ── Electronic band structure ─────────────────────────────────────────────────

/// Write electronic eigenvalues: `kpoints` shape `[nk, 3]`, `eigenvalues` shape `[nk, nbands]`.
pub fn write_band_structure(
    file: &mut Hdf5File,
    group: &str,
    nk: usize,
    nbands: usize,
    kpoints: &[f64],
    eigenvalues: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(kpoints.len(), nk * 3);
    assert_eq!(eigenvalues.len(), nk * nbands);
    file.create_group(group)?;
    let _ = file.create_dataset(group, "kpoints", vec![nk, 3], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "kpoints")?
        .write_f64(kpoints)?;
    let _ = file.create_dataset(group, "eigenvalues", vec![nk, nbands], Hdf5Dtype::Float64);
    file.open_dataset_mut(group, "eigenvalues")?
        .write_f64(eigenvalues)
}

// ── Sparse matrix I/O ─────────────────────────────────────────────────────────

/// Sparse matrix in COO format.
#[derive(Debug, Clone)]
pub struct SparseCoo {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row indices.
    pub row: Vec<usize>,
    /// Column indices.
    pub col: Vec<usize>,
    /// Non-zero values.
    pub data: Vec<f64>,
}

impl SparseCoo {
    /// Create an empty sparse matrix.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            row: Vec::new(),
            col: Vec::new(),
            data: Vec::new(),
        }
    }

    /// Push a non-zero entry.
    pub fn push(&mut self, r: usize, c: usize, v: f64) {
        self.row.push(r);
        self.col.push(c);
        self.data.push(v);
    }

    /// Number of non-zeros.
    pub fn nnz(&self) -> usize {
        self.data.len()
    }
}

/// Write a [`SparseCoo`] matrix.
pub fn write_sparse_coo(file: &mut Hdf5File, group: &str, mat: &SparseCoo) -> Hdf5Result<()> {
    let row_f: Vec<f64> = mat.row.iter().map(|&x| x as f64).collect();
    let col_f: Vec<f64> = mat.col.iter().map(|&x| x as f64).collect();
    write_f64_dataset(file, group, "row", &row_f)?;
    write_f64_dataset(file, group, "col", &col_f)?;
    write_f64_dataset(file, group, "data", &mat.data)?;
    file.open_group_mut(group)?
        .set_attr("nrows", AttrValue::Int32(mat.nrows as i32));
    file.open_group_mut(group)?
        .set_attr("ncols", AttrValue::Int32(mat.ncols as i32));
    Ok(())
}

/// Read a [`SparseCoo`] matrix.
pub fn read_sparse_coo(file: &Hdf5File, group: &str) -> Hdf5Result<SparseCoo> {
    let g = file.open_group(group)?;
    let nrows = match g.attributes.get("nrows") {
        Some(AttrValue::Int32(v)) => *v as usize,
        _ => 0,
    };
    let ncols = match g.attributes.get("ncols") {
        Some(AttrValue::Int32(v)) => *v as usize,
        _ => 0,
    };
    let row: Vec<usize> = file
        .open_dataset(group, "row")?
        .read_f64()?
        .iter()
        .map(|&x| x as usize)
        .collect();
    let col: Vec<usize> = file
        .open_dataset(group, "col")?
        .read_f64()?
        .iter()
        .map(|&x| x as usize)
        .collect();
    let data = file.open_dataset(group, "data")?.read_f64()?;
    Ok(SparseCoo {
        nrows,
        ncols,
        row,
        col,
        data,
    })
}
