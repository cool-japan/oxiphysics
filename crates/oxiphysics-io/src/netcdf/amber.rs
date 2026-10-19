// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The AMBER NetCDF trajectory convention
//! (<https://ambermd.org/netcdf/nctraj.xhtml>, version 1.0).
//!
//! Files built here and serialised with
//! [`write_cdf`](super::classic::write_cdf) /
//! [`CdfVersion::Offset64`] (the variant the convention requires) are read
//! by AMBER, cpptraj, VMD and MDAnalysis.

use super::classic::{
    CdfAttribute, CdfDataset, CdfVariable, CdfVersion, NcType, NcValues, read_cdf_dataset,
    write_cdf_dataset,
};
use super::types::{NetCdfFile, NetCdfVariable, NetcdfTrajectoryBuilder, TrajectoryFrame};
use crate::error::{Error, Result};

/// `scale_factor` of the AMBER `velocities` variable: stored values are in
/// Å per AKMA time unit and must be multiplied by this to give Å/ps.
pub const AMBER_VELOCITY_SCALE_FACTOR: f64 = 20.455;
/// Length of the `label` dimension (the longest cell-angle label, `alpha`).
const LABEL_LEN: usize = 5;
/// Contents of the `cell_angular(cell_angular, label)` variable.
const CELL_ANGULAR_LABELS: &str = "alphabeta gamma";
const ANGSTROM_PER_NM: f64 = 10.0;

fn dims(names: &[&str]) -> Vec<String> {
    names.iter().map(|n| (*n).to_string()).collect()
}

/// Create an AMBER-convention NetCDF trajectory as a [`NetCdfFile`].
///
/// The layout follows the AMBER trajectory convention 1.0:
/// - dimensions `frame` (unlimited, `n_frames` records), `spatial` (3) and
///   `atom`; with cell data also `cell_spatial` (3), `cell_angular` (3) and
///   `label` (5);
/// - `char spatial(spatial)` = `"xyz"`, and with cell data
///   `char cell_spatial(cell_spatial)` = `"abc"` and
///   `char cell_angular(cell_angular, label)` = `"alpha"`, `"beta "`,
///   `"gamma"`;
/// - `float coordinates(frame, atom, spatial)` in Å, `float time(frame)` in
///   ps, `float cell_lengths(frame, cell_spatial)` in Å and
///   `float cell_angles(frame, cell_angular)` in degrees;
/// - global attributes `Conventions`, `ConventionVersion`, `program` and
///   `programVersion`.
///
/// The convention includes `cell_angles` if and only if `cell_lengths` is
/// present; when only one of the two is supplied the other is written as
/// zeros, the convention's value for an undefined length or angle.
///
/// Write it with [`write_cdf`](super::classic::write_cdf) and
/// [`CdfVersion::Offset64`] to get a file other MD tools can read; the
/// data lengths are checked at that point.
pub fn create_amber_trajectory(
    n_atoms: usize,
    n_frames: usize,
    coordinates: Vec<f64>,
    time: Vec<f64>,
    cell_lengths: Option<Vec<f64>>,
    cell_angles: Option<Vec<f64>>,
) -> NetCdfFile {
    let mut file = NetCdfFile::new();
    file.add_unlimited_dimension("frame", n_frames);
    file.add_dimension("spatial", 3);
    file.add_dimension("atom", n_atoms);
    file.add_attribute("Conventions", "AMBER");
    file.add_attribute("ConventionVersion", "1.0");
    file.add_attribute("program", "OxiPhysics");
    file.add_attribute("programVersion", env!("CARGO_PKG_VERSION"));
    file.add_variable(NetCdfVariable::text("spatial", dims(&["spatial"]), "xyz"));
    let mut time_var =
        NetCdfVariable::new("time", dims(&["frame"]), time).with_nc_type(NcType::Float);
    time_var.add_attribute("units", "picosecond");
    file.add_variable(time_var);
    let mut coords_var = NetCdfVariable::new(
        "coordinates",
        dims(&["frame", "atom", "spatial"]),
        coordinates,
    )
    .with_nc_type(NcType::Float);
    coords_var.add_attribute("units", "angstrom");
    file.add_variable(coords_var);
    if cell_lengths.is_some() || cell_angles.is_some() {
        let zeros = || vec![0.0; n_frames * 3];
        file.add_dimension("cell_spatial", 3);
        file.add_dimension("cell_angular", 3);
        file.add_dimension("label", LABEL_LEN);
        file.add_variable(NetCdfVariable::text(
            "cell_spatial",
            dims(&["cell_spatial"]),
            "abc",
        ));
        file.add_variable(NetCdfVariable::text(
            "cell_angular",
            dims(&["cell_angular", "label"]),
            CELL_ANGULAR_LABELS,
        ));
        let mut lengths = NetCdfVariable::new(
            "cell_lengths",
            dims(&["frame", "cell_spatial"]),
            cell_lengths.unwrap_or_else(zeros),
        )
        .with_nc_type(NcType::Float);
        lengths.add_attribute("units", "angstrom");
        file.add_variable(lengths);
        let mut angles = NetCdfVariable::new(
            "cell_angles",
            dims(&["frame", "cell_angular"]),
            cell_angles.unwrap_or_else(zeros),
        )
        .with_nc_type(NcType::Float);
        angles.add_attribute("units", "degree");
        file.add_variable(angles);
    }
    file
}

fn invalid(msg: String) -> Error {
    Error::General(msg)
}

impl NetcdfTrajectoryBuilder {
    /// Build the AMBER-convention dataset of this trajectory.
    ///
    /// Positions (and box lengths) are converted to Å when the builder works
    /// in nm.  `cell_lengths` / `cell_angles` (90°, the boxes being
    /// orthorhombic) are written only when every frame has a box, and
    /// `velocities` only when every frame has velocities; these are converted
    /// from nm/ps to Å/ps and stored in AMBER units with the convention's
    /// `scale_factor` of 20.455.
    ///
    /// # Errors
    ///
    /// Fails when the frames do not all hold the same number of atoms, or a
    /// frame's velocities do not match its positions.
    pub fn to_amber_dataset(&self) -> Result<CdfDataset> {
        let n_atoms = self.n_atoms();
        let n_frames = self.frame_count();
        for (idx, frame) in self.frames.iter().enumerate() {
            if frame.n_atoms() != n_atoms {
                return Err(invalid(format!(
                    "trajectory frame {idx} holds {} atoms but frame 0 holds {n_atoms}",
                    frame.n_atoms()
                )));
            }
            if let Some(vels) = &frame.velocities
                && vels.len() != n_atoms
            {
                return Err(invalid(format!(
                    "trajectory frame {idx} holds {} velocities for {n_atoms} atoms",
                    vels.len()
                )));
            }
        }
        if n_atoms == 0 && n_frames > 0 {
            return Err(invalid(
                "an AMBER trajectory needs at least one atom".to_string(),
            ));
        }
        let length_scale = if self.use_angstroms {
            1.0
        } else {
            ANGSTROM_PER_NM
        };
        let has_cell = n_frames > 0 && self.frames.iter().all(|f| f.box_lengths.is_some());
        let has_velocities = n_frames > 0 && self.frames.iter().all(|f| f.velocities.is_some());

        let mut ds = CdfDataset::new();
        ds.add_unlimited_dimension("frame", n_frames)?;
        ds.add_dimension("spatial", 3)?;
        // A zero-atom (and therefore zero-frame) trajectory still needs a
        // non-zero fixed dimension length; the classic format reserves 0 for
        // the unlimited dimension.
        ds.add_dimension("atom", n_atoms.max(1))?;
        ds.set_text_attribute("Conventions", "AMBER");
        ds.set_text_attribute("ConventionVersion", "1.0");
        ds.set_text_attribute("program", &self.application);
        ds.set_text_attribute("programVersion", env!("CARGO_PKG_VERSION"));
        ds.set_text_attribute("application", &self.application);
        if !self.title.is_empty() {
            ds.set_text_attribute("title", &self.title);
        }
        ds.add_variable("spatial", &["spatial"], NcValues::text("xyz"))?;
        let time: Vec<f32> = self.frames.iter().map(|f| f.time_ps as f32).collect();
        let t = ds.add_variable("time", &["frame"], NcValues::Float(time))?;
        ds.variables[t].set_attribute(CdfAttribute::text("units", "picosecond"));
        let coords: Vec<f32> = self
            .frames
            .iter()
            .flat_map(|f| f.positions.iter().flatten())
            .map(|&v| (v * length_scale) as f32)
            .collect();
        let c = ds.add_variable(
            "coordinates",
            &["frame", "atom", "spatial"],
            NcValues::Float(coords),
        )?;
        ds.variables[c].set_attribute(CdfAttribute::text("units", "angstrom"));
        if has_cell {
            ds.add_dimension("cell_spatial", 3)?;
            ds.add_dimension("cell_angular", 3)?;
            ds.add_dimension("label", LABEL_LEN)?;
            ds.add_variable("cell_spatial", &["cell_spatial"], NcValues::text("abc"))?;
            ds.add_variable(
                "cell_angular",
                &["cell_angular", "label"],
                NcValues::text(CELL_ANGULAR_LABELS),
            )?;
            let lengths: Vec<f32> = self
                .frames
                .iter()
                .filter_map(|f| f.box_lengths)
                .flatten()
                .map(|v| (v * length_scale) as f32)
                .collect();
            let l = ds.add_variable(
                "cell_lengths",
                &["frame", "cell_spatial"],
                NcValues::Float(lengths),
            )?;
            ds.variables[l].set_attribute(CdfAttribute::text("units", "angstrom"));
            let a = ds.add_variable(
                "cell_angles",
                &["frame", "cell_angular"],
                NcValues::Float(vec![90.0; n_frames * 3]),
            )?;
            ds.variables[a].set_attribute(CdfAttribute::text("units", "degree"));
        }
        if has_velocities {
            let vels: Vec<f32> = self
                .frames
                .iter()
                .filter_map(|f| f.velocities.as_deref())
                .flat_map(|v| v.iter().flatten())
                .map(|&v| (v * ANGSTROM_PER_NM / AMBER_VELOCITY_SCALE_FACTOR) as f32)
                .collect();
            let v = ds.add_variable(
                "velocities",
                &["frame", "atom", "spatial"],
                NcValues::Float(vels),
            )?;
            ds.variables[v].set_attribute(CdfAttribute::text("units", "angstrom/picosecond"));
            ds.variables[v].set_attribute(CdfAttribute::new(
                "scale_factor",
                NcValues::Double(vec![AMBER_VELOCITY_SCALE_FACTOR]),
            ));
        }
        Ok(ds)
    }

    /// Write the trajectory as an AMBER NetCDF file (64-bit offset format,
    /// as the convention requires) — the binary counterpart of
    /// [`write_cdl`](Self::write_cdl).
    ///
    /// # Errors
    ///
    /// As [`to_amber_dataset`](Self::to_amber_dataset), plus encoding and
    /// I/O errors.
    pub fn write_amber_netcdf<W: std::io::Write>(&self, writer: &mut W) -> Result<()> {
        let bytes = write_cdf_dataset(&self.to_amber_dataset()?, CdfVersion::Offset64)?;
        writer.write_all(&bytes)?;
        Ok(())
    }
}

/// A trajectory read from an AMBER NetCDF file.
#[derive(Debug, Clone, PartialEq)]
pub struct AmberTrajectory {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Per-frame positions in Å.
    pub coordinates: Vec<Vec<[f64; 3]>>,
    /// Per-frame time in ps (`None` when the file has no `time` variable).
    pub time_ps: Option<Vec<f64>>,
    /// Per-frame cell lengths `a, b, c` in Å.
    pub cell_lengths: Option<Vec<[f64; 3]>>,
    /// Per-frame cell angles `alpha, beta, gamma` in degrees.
    pub cell_angles: Option<Vec<[f64; 3]>>,
    /// Per-frame velocities in Å/ps (`scale_factor` already applied).
    pub velocities: Option<Vec<Vec<[f64; 3]>>>,
    /// The `title` global attribute.
    pub title: Option<String>,
    /// The `program` global attribute.
    pub program: Option<String>,
}

impl AmberTrajectory {
    /// Number of frames.
    pub fn n_frames(&self) -> usize {
        self.coordinates.len()
    }

    /// Convert to [`TrajectoryFrame`]s (positions and box lengths in Å,
    /// velocities in nm/ps, time 0 when the file stores none).
    ///
    /// `box_lengths` is set only for orthorhombic cells (all angles 90° or
    /// all zero, i.e. undefined), since [`TrajectoryFrame`] cannot describe
    /// a triclinic cell; the full cell stays available in
    /// [`cell_lengths`](Self::cell_lengths) / [`cell_angles`](Self::cell_angles).
    pub fn to_frames(&self) -> Vec<TrajectoryFrame> {
        (0..self.n_frames())
            .map(|i| {
                let orthorhombic =
                    self.cell_angles
                        .as_ref()
                        .and_then(|a| a.get(i))
                        .is_none_or(|a| {
                            a.iter().all(|&x| (x - 90.0).abs() < 1e-3)
                                || a.iter().all(|&x| x == 0.0)
                        });
                TrajectoryFrame {
                    time_ps: self
                        .time_ps
                        .as_ref()
                        .and_then(|t| t.get(i).copied())
                        .unwrap_or(0.0),
                    positions: self.coordinates.get(i).cloned().unwrap_or_default(),
                    velocities: self.velocities.as_ref().and_then(|v| {
                        v.get(i).map(|frame| {
                            frame
                                .iter()
                                .map(|p| p.map(|x| x / ANGSTROM_PER_NM))
                                .collect()
                        })
                    }),
                    box_lengths: if orthorhombic {
                        self.cell_lengths.as_ref().and_then(|l| l.get(i).copied())
                    } else {
                        None
                    },
                }
            })
            .collect()
    }

    /// A builder (in Å) holding [`to_frames`](Self::to_frames).
    pub fn to_builder(&self) -> NetcdfTrajectoryBuilder {
        let mut builder = NetcdfTrajectoryBuilder::new().in_angstroms();
        if let Some(title) = &self.title {
            builder = builder.with_title(title);
        }
        if let Some(program) = &self.program {
            builder = builder.with_application(program);
        }
        for frame in self.to_frames() {
            builder.add_frame(frame);
        }
        builder
    }
}

fn parse_err(msg: String) -> Error {
    Error::Parse(format!("AMBER NetCDF: {msg}"))
}

/// Values of a numeric variable as `f64`, multiplied by its `scale_factor`
/// attribute when present.
fn scaled_values(var: &CdfVariable) -> Result<Vec<f64>> {
    if var.nc_type() == NcType::Char {
        return Err(parse_err(format!("variable '{}' is not numeric", var.name)));
    }
    let mut values = var.values.to_f64_vec();
    if let Some(attr) = var.attribute("scale_factor") {
        match attr.values.to_f64_vec().as_slice() {
            [factor] if attr.values.nc_type() != NcType::Char => {
                values.iter_mut().for_each(|v| *v *= factor);
            }
            _ => {
                return Err(parse_err(format!(
                    "variable '{}' has a malformed scale_factor",
                    var.name
                )));
            }
        }
    }
    Ok(values)
}

/// Read a per-frame variable of shape `(frame, <dims of length 3>...)`.
fn frame_triplets(
    ds: &CdfDataset,
    name: &str,
    expected_shape: &[usize],
) -> Result<Option<Vec<[f64; 3]>>> {
    let Some(var) = ds.variable(name) else {
        return Ok(None);
    };
    let shape = ds.variable_shape(var);
    if !ds.is_record_variable(var) || shape != expected_shape {
        return Err(parse_err(format!(
            "variable '{name}' has shape {shape:?}, expected {expected_shape:?} over the frame dimension"
        )));
    }
    let values = scaled_values(var)?;
    Ok(Some(
        values.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect(),
    ))
}

fn per_frame_atoms(triplets: Vec<[f64; 3]>, n_atoms: usize) -> Vec<Vec<[f64; 3]>> {
    triplets
        .chunks(n_atoms.max(1))
        .map(<[[f64; 3]]>::to_vec)
        .collect()
}

/// Parse an AMBER NetCDF trajectory (any classic variant; AMBER writes the
/// 64-bit offset one).
///
/// Requires `Conventions` to list `AMBER`, an unlimited `frame` dimension,
/// `spatial` of length 3 and a `coordinates(frame, atom, spatial)` variable.
/// `time`, `cell_lengths`, `cell_angles` and `velocities` are read when
/// present; `scale_factor` attributes are applied.
///
/// # Errors
///
/// [`Error::Parse`] when the bytes are not a valid NetCDF classic file or do
/// not follow the convention.
pub fn read_amber_trajectory(bytes: &[u8]) -> Result<AmberTrajectory> {
    let ds = read_cdf_dataset(bytes)?;
    let conventions = ds.text_attribute("Conventions").unwrap_or_default();
    if !conventions.split([',', ' ']).any(|c| c.trim() == "AMBER") {
        return Err(parse_err(format!(
            "Conventions is {conventions:?}, expected it to include \"AMBER\""
        )));
    }
    let frame_dim = ds
        .dimension("frame")
        .filter(|d| d.unlimited)
        .ok_or_else(|| parse_err("missing unlimited 'frame' dimension".to_string()))?;
    let n_frames = frame_dim.len;
    if ds.dimension("spatial").map(|d| d.len) != Some(3) {
        return Err(parse_err(
            "'spatial' dimension must have length 3".to_string(),
        ));
    }
    let n_atoms = ds
        .dimension("atom")
        .map(|d| d.len)
        .ok_or_else(|| parse_err("missing 'atom' dimension".to_string()))?;
    let coordinates = frame_triplets(&ds, "coordinates", &[n_frames, n_atoms, 3])?
        .ok_or_else(|| parse_err("missing 'coordinates' variable".to_string()))?;
    let time_ps = match ds.variable("time") {
        Some(var) => {
            if !ds.is_record_variable(var) || ds.variable_shape(var) != [n_frames] {
                return Err(parse_err("'time' must have shape (frame)".to_string()));
            }
            Some(scaled_values(var)?)
        }
        None => None,
    };
    let cell_lengths = frame_triplets(&ds, "cell_lengths", &[n_frames, 3])?;
    let cell_angles = frame_triplets(&ds, "cell_angles", &[n_frames, 3])?;
    let velocities = frame_triplets(&ds, "velocities", &[n_frames, n_atoms, 3])?
        .map(|v| per_frame_atoms(v, n_atoms));
    Ok(AmberTrajectory {
        n_atoms,
        coordinates: per_frame_atoms(coordinates, n_atoms),
        time_ps,
        cell_lengths,
        cell_angles,
        velocities,
        title: ds.text_attribute("title"),
        program: ds.text_attribute("program"),
    })
}
