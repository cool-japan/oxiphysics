// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Frame and atom counts in a DCD header are untrusted: they must never size
//! an allocation, and a stream shorter than the header declares is an error.

use super::{SimpleDcdReader, SimpleDcdWriter};

const N_FRAMES_OFFSET: usize = 8;
const N_ATOMS_OFFSET: usize = 12;

fn trajectory(n_atoms: u32, n_frames: usize) -> Vec<u8> {
    let mut writer = SimpleDcdWriter::new(n_atoms, 0.002);
    for frame in 0..n_frames {
        let positions: Vec<[f32; 3]> = (0..n_atoms)
            .map(|a| [frame as f32, a as f32, 0.5])
            .collect();
        writer.add_frame(&positions);
    }
    writer.to_bytes()
}

/// The error message of a read that must fail.  (The reader's `Ok` type does
/// not implement `Debug`, so `Result::unwrap_err` is not available.)
fn read_error(bytes: &[u8]) -> String {
    match SimpleDcdReader::from_bytes(bytes) {
        Ok((header, frames)) => panic!(
            "expected an error, but read {} frames of {} atoms",
            frames.len(),
            header.n_atoms
        ),
        Err(message) => message,
    }
}

fn with_counts(mut bytes: Vec<u8>, n_frames: u32, n_atoms: u32) -> Vec<u8> {
    bytes[N_FRAMES_OFFSET..N_FRAMES_OFFSET + 4].copy_from_slice(&n_frames.to_le_bytes());
    bytes[N_ATOMS_OFFSET..N_ATOMS_OFFSET + 4].copy_from_slice(&n_atoms.to_le_bytes());
    bytes
}

#[test]
fn honest_trajectory_still_reads() {
    let bytes = trajectory(3, 4);
    let (header, frames) = SimpleDcdReader::from_bytes(&bytes).unwrap();
    assert_eq!((header.n_atoms, header.n_frames), (3, 4));
    assert_eq!(frames.len(), 4);
    assert_eq!(frames[2].x, vec![2.0; 3]);
    assert_eq!(frames[2].y, vec![0.0, 1.0, 2.0]);
}

#[test]
fn absurd_frame_count_is_an_error_not_an_abort() {
    let bytes = with_counts(trajectory(2, 1), u32::MAX, 2);
    let err = read_error(&bytes);
    assert!(err.contains("declares"), "{err}");
}

#[test]
fn absurd_atom_count_is_an_error_not_an_abort() {
    let bytes = with_counts(trajectory(2, 1), 1, u32::MAX);
    assert!(SimpleDcdReader::from_bytes(&bytes).is_err());
}

#[test]
fn both_counts_at_the_maximum_are_an_error() {
    let bytes = with_counts(trajectory(2, 1), u32::MAX, u32::MAX);
    assert!(SimpleDcdReader::from_bytes(&bytes).is_err());
}

#[test]
fn stream_shorter_than_declared_is_an_error() {
    let bytes = with_counts(trajectory(2, 3), 4, 2);
    let err = read_error(&bytes);
    assert!(err.contains("4 frames"), "{err}");
    // A stream cut in the middle of the last frame is rejected as well.
    let mut cut = trajectory(2, 3);
    cut.truncate(cut.len() - 3);
    assert!(SimpleDcdReader::from_bytes(&cut).is_err());
}

#[test]
fn zero_frames_need_no_coordinate_data() {
    let bytes = with_counts(trajectory(2, 0), 0, u32::MAX);
    let (header, frames) = SimpleDcdReader::from_bytes(&bytes).unwrap();
    assert_eq!(header.n_atoms, u32::MAX);
    assert!(frames.is_empty());
}
