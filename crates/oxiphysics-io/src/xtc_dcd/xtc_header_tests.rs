// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The frame and atom counts of an XTC stream are untrusted: they must never
//! size an allocation before they have been checked against the bytes that
//! remain, and a stream shorter than it declares is an error.

use super::{SimpleXtcReader, SimpleXtcWriter, XtcVelocityWriter, xtc_to_dcd};

/// Byte offset of the frame count in the stream header (after the magic).
const N_FRAMES_OFFSET: usize = 4;
/// Byte offset of the first frame's atom count: header (8) + step (4) + time (4).
const FIRST_FRAME_N_ATOMS_OFFSET: usize = 16;

fn stream(n_frames: usize, n_atoms: usize) -> Vec<u8> {
    let mut writer = SimpleXtcWriter::new();
    for frame in 0..n_frames {
        let positions: Vec<[f32; 3]> = (0..n_atoms)
            .map(|a| [frame as f32, a as f32, 0.25])
            .collect();
        writer.add_frame(frame as i32, frame as f32 * 0.5, positions);
    }
    writer.to_bytes()
}

fn with_u32_at(mut bytes: Vec<u8>, offset: usize, value: u32) -> Vec<u8> {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    bytes
}

/// The error message of a read that must fail. (The reader's `Ok` type does
/// not implement `Debug`, so `Result::unwrap_err` is not available.)
fn read_error(bytes: &[u8]) -> String {
    match SimpleXtcReader::from_bytes(bytes) {
        Ok(frames) => panic!("expected an error, but read {} frames", frames.len()),
        Err(message) => message,
    }
}

#[test]
fn honest_stream_still_reads() {
    let frames = SimpleXtcReader::from_bytes(&stream(3, 4)).unwrap();
    assert_eq!(frames.len(), 3);
    assert_eq!(frames[2].step, 2);
    assert_eq!(frames[2].positions.len(), 4);
    assert_eq!(frames[2].positions[3], [2.0, 3.0, 0.25]);
}

#[test]
fn magic_plus_maximum_frame_count_is_an_error_not_an_abort() {
    // 8 bytes: the magic and n_frames = 0xFFFFFFFF. Reserving 4 294 967 295
    // frames used to abort the process before any EOF check ran.
    let mut bytes = stream(0, 0);
    assert_eq!(bytes.len(), 8);
    bytes[N_FRAMES_OFFSET..].copy_from_slice(&u32::MAX.to_le_bytes());
    let err = read_error(&bytes);
    assert!(err.contains("declares 4294967295 frames"), "{err}");
}

#[test]
fn frame_count_beyond_the_stream_is_an_error() {
    // Two empty frames take 2 * 48 bytes; three cannot fit.
    let bytes = with_u32_at(stream(2, 0), N_FRAMES_OFFSET, 3);
    let err = read_error(&bytes);
    assert!(err.contains("declares 3 frames"), "{err}");
    // With atoms the early size check cannot tell, but the missing third
    // frame is still an error (end of stream), never a panic.
    let bytes = with_u32_at(stream(2, 3), N_FRAMES_OFFSET, 3);
    let err = read_error(&bytes);
    assert!(err.contains("unexpected EOF"), "{err}");
}

#[test]
fn absurd_atom_count_is_an_error_not_an_abort() {
    let bytes = with_u32_at(stream(1, 2), FIRST_FRAME_N_ATOMS_OFFSET, u32::MAX);
    let err = read_error(&bytes);
    assert!(err.contains("frame 0 declares 4294967295 atoms"), "{err}");
}

#[test]
fn atom_count_beyond_the_stream_is_an_error() {
    // Two atoms are present; the header claims three.
    let bytes = with_u32_at(stream(1, 2), FIRST_FRAME_N_ATOMS_OFFSET, 3);
    let err = read_error(&bytes);
    assert!(err.contains("frame 0 declares 3 atoms"), "{err}");
}

#[test]
fn hostile_atom_count_in_a_later_frame_is_an_error() {
    // Stream header (8) + first frame (step, time, count, box, 2 atoms) + the
    // second frame's step and time.
    let first_frame_bytes = 48 + 12 * 2;
    let second_frame_n_atoms = 8 + first_frame_bytes + 8;
    let bytes = with_u32_at(stream(2, 2), second_frame_n_atoms, u32::MAX);
    let err = read_error(&bytes);
    assert!(err.contains("frame 1 declares"), "{err}");
}

#[test]
fn truncated_stream_is_an_error() {
    let mut cut = stream(2, 3);
    cut.truncate(cut.len() - 5);
    assert!(SimpleXtcReader::from_bytes(&cut).is_err());
    // Cut inside the header.
    assert!(SimpleXtcReader::from_bytes(&stream(1, 1)[..6]).is_err());
    assert!(SimpleXtcReader::from_bytes(&[]).is_err());
}

#[test]
fn zero_frames_read_as_an_empty_trajectory() {
    let frames = SimpleXtcReader::from_bytes(&stream(0, 0)).unwrap();
    assert!(frames.is_empty());
}

#[test]
fn bytes_after_the_declared_frames_are_ignored() {
    let mut velocity_writer = XtcVelocityWriter::new();
    velocity_writer.add_frame(
        0,
        0.0,
        vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        vec![[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]],
    );
    let frames = SimpleXtcReader::from_bytes(&velocity_writer.to_bytes()).unwrap();
    assert_eq!(frames.len(), 1);
    assert_eq!(frames[0].positions, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
}

#[test]
fn consumers_of_the_reader_propagate_the_error() {
    let hostile = with_u32_at(stream(1, 1), N_FRAMES_OFFSET, u32::MAX);
    assert!(SimpleXtcReader::read_frame_at(&hostile, 0).is_err());
    assert!(xtc_to_dcd(&hostile).is_err());
}
