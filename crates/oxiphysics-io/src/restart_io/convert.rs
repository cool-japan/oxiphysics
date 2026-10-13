// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! One-way exports of a restart snapshot to other text formats.
//!
//! These are lossy human-readable conversions (six decimals); use
//! [`RestartWriter`](super::RestartWriter) for lossless restart files.

use super::RestartData;

// ── Conversion utilities ──────────────────────────────────────────────────────

/// Convert a `RestartData` snapshot to an XYZ-format string.
///
/// The XYZ file contains a single frame. Particle types 0, 1, 2 and 3 are
/// written as the element symbols `H`, `C`, `N` and `O`; any other type, and
/// particles without a type entry, are written as `X`.
pub fn restart_to_xyz(data: &RestartData) -> String {
    let n = data.n_particles();
    let mut s = String::new();
    s.push_str(&format!("{n}\n"));
    s.push_str(&format!(
        "Restart step={} time={:.6}\n",
        data.metadata.step, data.metadata.time
    ));
    for (i, p) in data.positions.iter().enumerate() {
        let sym = if let Some(&t) = data.types.get(i) {
            match t {
                0 => "H",
                1 => "C",
                2 => "N",
                3 => "O",
                _ => "X",
            }
        } else {
            "X"
        };
        s.push_str(&format!("{sym} {:.6} {:.6} {:.6}\n", p[0], p[1], p[2]));
    }
    s
}

/// Convert a `RestartData` snapshot to a LAMMPS dump-format string.
///
/// Fields: `id type x y z vx vy vz fx fy fz mass`
pub fn restart_to_lammps_dump(data: &RestartData) -> String {
    let n = data.n_particles();
    let step = data.metadata.step;
    let mut s = String::new();
    s.push_str(&format!("ITEM: TIMESTEP\n{step}\n"));
    s.push_str(&format!("ITEM: NUMBER OF ATOMS\n{n}\n"));
    let bm = &data.box_matrix;
    s.push_str("ITEM: BOX BOUNDS pp pp pp\n");
    s.push_str(&format!(
        "0.0 {:.6}\n0.0 {:.6}\n0.0 {:.6}\n",
        bm[0][0], bm[1][1], bm[2][2]
    ));
    s.push_str("ITEM: ATOMS id type x y z vx vy vz fx fy fz mass\n");
    for i in 0..n {
        let p = data.positions.get(i).copied().unwrap_or([0.0; 3]);
        let v = data.velocities.get(i).copied().unwrap_or([0.0; 3]);
        let f = data.forces.get(i).copied().unwrap_or([0.0; 3]);
        let m = data.masses.get(i).copied().unwrap_or(1.0);
        let tp = data.types.get(i).copied().unwrap_or(0);
        s.push_str(&format!(
            "{} {} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\n",
            i + 1,
            tp,
            p[0],
            p[1],
            p[2],
            v[0],
            v[1],
            v[2],
            f[0],
            f[1],
            f[2],
            m
        ));
    }
    s
}
