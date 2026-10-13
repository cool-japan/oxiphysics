// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Column layout of an `ITEM: ATOMS` table and the simulation-box geometry used
//! to turn scaled coordinates into Cartesian ones.
//!
//! LAMMPS lets the user pick the dump columns and their order (`dump custom`),
//! and the common `dump atom` style writes *scaled* coordinates (`xs ys zs`).
//! The reader therefore maps columns by the names in the header instead of by
//! position.

use super::LammpsDumpFrame;
use std::collections::HashMap;

/// Columns that hold text instead of numbers.
///
/// [`LammpsDumpFrame`] stores per-atom values as `f64`, so these columns are
/// skipped (the `type` column still identifies the species).
const TEXT_COLUMNS: [&str; 1] = ["element"];

/// Position column triples in order of preference, with whether they are scaled.
const POSITION_FAMILIES: [([&str; 3], bool); 4] = [
    (["x", "y", "z"], false),
    (["xu", "yu", "zu"], false),
    (["xs", "ys", "zs"], true),
    (["xsu", "ysu", "zsu"], true),
];

/// Where each per-atom quantity lives in an `ITEM: ATOMS` row.
#[derive(Debug, Clone)]
pub(super) struct DumpLayout {
    /// Number of columns named in the header.
    n_columns: usize,
    id: usize,
    atom_type: usize,
    position: [usize; 3],
    /// The position columns are fractions of the box edges (`xs`, `xsu`).
    scaled: bool,
    velocity: Option<[usize; 3]>,
    /// `(column, name)` of every numeric column not consumed above.
    custom: Vec<(usize, String)>,
}

impl DumpLayout {
    /// Build the layout from the column names of an `ITEM: ATOMS` header.
    ///
    /// Requires `id`, `type` and a complete position triple (`x y z`,
    /// `xu yu zu`, `xs ys zs` or `xsu ysu zsu`).  `vx vy vz` become velocities
    /// when all three are present.  Every other column except the text
    /// columns ([`TEXT_COLUMNS`]) is kept as a custom property.
    pub(super) fn from_header(columns: &[&str]) -> Result<Self, String> {
        let mut index: HashMap<&str, usize> = HashMap::new();
        for (i, name) in columns.iter().enumerate() {
            index.entry(*name).or_insert(i);
        }
        let id = *index
            .get("id")
            .ok_or_else(|| missing_column("id", columns))?;
        let atom_type = *index
            .get("type")
            .ok_or_else(|| missing_column("type", columns))?;
        let (position, scaled) = POSITION_FAMILIES
            .iter()
            .find_map(|(names, scaled)| {
                let cols = [
                    *index.get(names[0])?,
                    *index.get(names[1])?,
                    *index.get(names[2])?,
                ];
                Some((cols, *scaled))
            })
            .ok_or_else(|| {
                format!(
                    "LAMMPS dump: ATOMS header {columns:?} has no complete position columns \
                     (x y z, xu yu zu, xs ys zs or xsu ysu zsu)"
                )
            })?;
        let velocity = match (index.get("vx"), index.get("vy"), index.get("vz")) {
            (Some(&vx), Some(&vy), Some(&vz)) => Some([vx, vy, vz]),
            _ => None,
        };
        let mut consumed = vec![id, atom_type];
        consumed.extend(position);
        if let Some(cols) = velocity {
            consumed.extend(cols);
        }
        let custom = columns
            .iter()
            .enumerate()
            .filter(|(i, name)| !consumed.contains(i) && !TEXT_COLUMNS.contains(name))
            .map(|(i, name)| (i, (*name).to_string()))
            .collect();
        Ok(Self {
            n_columns: columns.len(),
            id,
            atom_type,
            position,
            scaled,
            velocity,
            custom,
        })
    }

    /// Whether the table carries `vx vy vz`.
    pub(super) fn has_velocities(&self) -> bool {
        self.velocity.is_some()
    }

    /// Names of the custom property columns, in header order.
    pub(super) fn custom_names(&self) -> impl Iterator<Item = &str> {
        self.custom.iter().map(|(_, name)| name.as_str())
    }
}

fn missing_column(name: &str, columns: &[&str]) -> String {
    format!("LAMMPS dump: ATOMS header {columns:?} has no '{name}' column")
}

/// Per-atom values collected while reading the rows of one frame.
pub(super) struct AtomTable {
    pub(super) ids: Vec<u32>,
    pub(super) types: Vec<u32>,
    pub(super) positions: Vec<[f64; 3]>,
    pub(super) velocities: Option<Vec<[f64; 3]>>,
    pub(super) custom: Vec<Vec<f64>>,
}

impl AtomTable {
    /// Allocate room for at most `expected_atoms` rows.
    ///
    /// The count comes from the file, so the reservation is capped; a bogus
    /// header cannot make the reader allocate gigabytes up front.
    pub(super) fn new(layout: &DumpLayout, expected_atoms: usize) -> Self {
        let capacity = expected_atoms.min(1 << 20);
        Self {
            ids: Vec::with_capacity(capacity),
            types: Vec::with_capacity(capacity),
            positions: Vec::with_capacity(capacity),
            velocities: layout
                .has_velocities()
                .then(|| Vec::with_capacity(capacity)),
            custom: layout
                .custom
                .iter()
                .map(|_| Vec::with_capacity(capacity))
                .collect(),
        }
    }

    /// Parse one atom row; `atom` is the zero-based row number for messages.
    pub(super) fn push_row(
        &mut self,
        layout: &DumpLayout,
        dump_box: &DumpBox,
        line: &str,
        atom: usize,
    ) -> Result<(), String> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < layout.n_columns {
            return Err(format!(
                "LAMMPS dump: atom line {atom} has {} columns, the header names {}",
                parts.len(),
                layout.n_columns
            ));
        }
        let id = parse_cell::<u32>(&parts, layout.id, "id", atom)?;
        let atom_type = parse_cell::<u32>(&parts, layout.atom_type, "type", atom)?;
        let raw = [
            parse_cell::<f64>(&parts, layout.position[0], "position", atom)?,
            parse_cell::<f64>(&parts, layout.position[1], "position", atom)?,
            parse_cell::<f64>(&parts, layout.position[2], "position", atom)?,
        ];
        let position = if layout.scaled {
            dump_box.to_cartesian(raw)
        } else {
            raw
        };
        let velocity = match layout.velocity {
            Some(cols) => Some([
                parse_cell::<f64>(&parts, cols[0], "vx", atom)?,
                parse_cell::<f64>(&parts, cols[1], "vy", atom)?,
                parse_cell::<f64>(&parts, cols[2], "vz", atom)?,
            ]),
            None => None,
        };
        // Parse every custom value before storing anything, so a bad row leaves the
        // table untouched.
        let mut custom_values = Vec::with_capacity(layout.custom.len());
        for (col, name) in &layout.custom {
            custom_values.push(parse_cell::<f64>(&parts, *col, name, atom)?);
        }
        self.ids.push(id);
        self.types.push(atom_type);
        self.positions.push(position);
        if let (Some(store), Some(v)) = (self.velocities.as_mut(), velocity) {
            store.push(v);
        }
        for (store, value) in self.custom.iter_mut().zip(custom_values) {
            store.push(value);
        }
        Ok(())
    }
}

fn parse_cell<T: std::str::FromStr>(
    parts: &[&str],
    col: usize,
    label: &str,
    atom: usize,
) -> Result<T, String> {
    let token = parts
        .get(col)
        .ok_or_else(|| format!("LAMMPS dump: atom {atom} has no column {col} ({label})"))?;
    token.parse::<T>().map_err(|_| {
        format!("LAMMPS dump: atom {atom} column {col} ({label}): cannot parse '{token}'")
    })
}

/// The simulation box as a cell matrix, with the origin and tilts of a LAMMPS
/// triclinic box: `a = (lx, 0, 0)`, `b = (xy, ly, 0)`, `c = (xz, yz, lz)`.
#[derive(Debug, Clone, Copy)]
pub(super) struct DumpBox {
    origin: [f64; 3],
    lengths: [f64; 3],
    /// `[xy, xz, yz]`; all zero for an orthogonal box.
    tilt: [f64; 3],
}

impl DumpBox {
    /// Build the box from the values of a `BOX BOUNDS` item.
    ///
    /// For a triclinic dump the bounds are those of the *bounding box* of the
    /// cell (`xlo_bound`, `xhi_bound`, ...); the true edges follow from the
    /// tilts:
    ///
    /// ```text
    /// xlo = xlo_bound - min(0, xy, xz, xy + xz)   xhi = xhi_bound - max(0, xy, xz, xy + xz)
    /// ylo = ylo_bound - min(0, yz)                yhi = yhi_bound - max(0, yz)
    /// ```
    pub(super) fn from_bounds(bounds: [[f64; 2]; 3], tilt_factors: Option<[f64; 3]>) -> Self {
        let [xy, xz, yz] = tilt_factors.unwrap_or([0.0; 3]);
        let shift_lo = [0.0_f64.min(xy).min(xz).min(xy + xz), 0.0_f64.min(yz), 0.0];
        let shift_hi = [0.0_f64.max(xy).max(xz).max(xy + xz), 0.0_f64.max(yz), 0.0];
        let lo = [
            bounds[0][0] - shift_lo[0],
            bounds[1][0] - shift_lo[1],
            bounds[2][0] - shift_lo[2],
        ];
        let hi = [
            bounds[0][1] - shift_hi[0],
            bounds[1][1] - shift_hi[1],
            bounds[2][1] - shift_hi[2],
        ];
        Self {
            origin: lo,
            lengths: [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]],
            tilt: [xy, xz, yz],
        }
    }

    /// Box of a parsed frame.
    pub(super) fn from_frame(frame: &LammpsDumpFrame) -> Self {
        Self::from_bounds(frame.box_bounds, frame.tilt_factors)
    }

    /// Cell volume `lx * ly * lz` (tilts do not change it).
    pub(super) fn volume(&self) -> f64 {
        self.lengths[0] * self.lengths[1] * self.lengths[2]
    }

    /// Lowest corner of the box.
    pub(super) fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Whether every edge has a positive, finite length.
    pub(super) fn is_valid(&self) -> bool {
        self.lengths.iter().all(|l| l.is_finite() && *l > 0.0)
            && self.tilt.iter().all(|t| t.is_finite())
    }

    /// Cartesian position of the scaled coordinates `s`:
    /// `x = xlo + xs*lx + ys*xy + zs*xz`, `y = ylo + ys*ly + zs*yz`,
    /// `z = zlo + zs*lz`.
    pub(super) fn to_cartesian(self, s: [f64; 3]) -> [f64; 3] {
        let [lx, ly, lz] = self.lengths;
        let [xy, xz, yz] = self.tilt;
        [
            self.origin[0] + s[0] * lx + s[1] * xy + s[2] * xz,
            self.origin[1] + s[1] * ly + s[2] * yz,
            self.origin[2] + s[2] * lz,
        ]
    }

    /// Fractional coordinates of the displacement `d` (no origin shift).
    pub(super) fn to_fractional_displacement(self, d: [f64; 3]) -> [f64; 3] {
        let [lx, ly, lz] = self.lengths;
        let [xy, xz, yz] = self.tilt;
        let s2 = d[2] / lz;
        let s1 = (d[1] - yz * s2) / ly;
        let s0 = (d[0] - xy * s1 - xz * s2) / lx;
        [s0, s1, s2]
    }

    /// Displacement `d` replaced by its minimum image under the box's periodic
    /// boundaries.
    ///
    /// Each fractional component is reduced to `[-0.5, 0.5]`.  For a distance
    /// below half the smallest perpendicular width of the cell
    /// ([`DumpBox::max_periodic_distance`]) this is exactly the minimum image.
    pub(super) fn minimum_image(&self, d: [f64; 3]) -> [f64; 3] {
        if self.tilt == [0.0; 3] {
            return [
                d[0] - self.lengths[0] * (d[0] / self.lengths[0]).round(),
                d[1] - self.lengths[1] * (d[1] / self.lengths[1]).round(),
                d[2] - self.lengths[2] * (d[2] / self.lengths[2]).round(),
            ];
        }
        let s = self.to_fractional_displacement(d);
        let wrapped = [
            s[0] - s[0].round(),
            s[1] - s[1].round(),
            s[2] - s[2].round(),
        ];
        let [lx, ly, lz] = self.lengths;
        let [xy, xz, yz] = self.tilt;
        [
            wrapped[0] * lx + wrapped[1] * xy + wrapped[2] * xz,
            wrapped[1] * ly + wrapped[2] * yz,
            wrapped[2] * lz,
        ]
    }

    /// Distances between opposite cell faces: the widths of the box measured
    /// perpendicular to the `bc`, `ca` and `ab` faces.
    pub(super) fn perpendicular_widths(&self) -> [f64; 3] {
        let [_, ly, lz] = self.lengths;
        let [xy, xz, yz] = self.tilt;
        // |b x c| = |(ly*lz, -xy*lz, xy*yz - ly*xz)|, |c x a| = lx*|(0, lz, -yz)|.
        let bc = (ly * lz).hypot(xy * lz).hypot(xy * yz - ly * xz);
        [self.volume() / bc, ly * lz / lz.hypot(yz), lz]
    }

    /// Largest distance for which [`DumpBox::minimum_image`] is exact: half the
    /// smallest perpendicular width.
    pub(super) fn max_periodic_distance(&self) -> f64 {
        let widths = self.perpendicular_widths();
        0.5 * widths[0].min(widths[1]).min(widths[2])
    }
}
