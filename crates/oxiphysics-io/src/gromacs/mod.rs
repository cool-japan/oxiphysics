// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GROMACS GRO format reader and writer.
//!
//! GRO file format:
//! ```text
//! Title line
//! N_atoms
//! %5d%5s%5s%5d%8.3f%8.3f%8.3f[%8.4f%8.4f%8.4f]  (resid, resname, atomname, atomid, x, y, z, [vx, vy, vz])
//! box_x box_y box_z
//! ```
//! All coordinates are in nm; velocities are in nm/ps.

pub mod formats;
pub mod gro_file;
pub mod topology;

// Re-export everything for backwards compatibility
pub use formats::*;
pub use gro_file::*;
pub use topology::*;
