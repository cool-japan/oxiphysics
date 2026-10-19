// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the frame readers of the `lammps` module: the text dump reader,
//! the binary dump reader and the restart reader.  The binary readers must
//! survive corrupt counts without allocating from them.

use super::types::{
    LammpsAtom, LammpsBinaryDumpReader, LammpsBinaryDumpWriter, LammpsDumpReader, LammpsDumpWriter,
    LammpsRestartReader, LammpsRestartWriter,
};
use crate::Error;
use oxiphysics_core::math::Vec3;
use std::io::{Cursor, Read};

/// Small deterministic xorshift generator (the crate has no `rand` dependency).
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform index in `0..n` (`n` must not be zero).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

fn atom(id: usize, type_id: usize, x: f64) -> LammpsAtom {
    LammpsAtom {
        id,
        type_id,
        position: Vec3::new(x, x + 0.5, x + 1.0),
        velocity: Vec3::new(-x, 0.25 * x, 2.0),
    }
}

fn assert_vec(got: Vec3, want: [f64; 3]) {
    assert!(
        (got.x - want[0]).abs() < 1e-12
            && (got.y - want[1]).abs() < 1e-12
            && (got.z - want[2]).abs() < 1e-12,
        "{got:?} vs {want:?}"
    );
}

/// The message of an `Error::Parse`; anything else fails the test.
fn parse_message<T: std::fmt::Debug>(result: crate::Result<T>) -> String {
    match result {
        Err(Error::Parse(message)) => message,
        other => panic!("expected Error::Parse, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// text dump reader
// ---------------------------------------------------------------------------

/// One dump frame with a 10 x 20 x 30 box.
fn dump_frame(timestep: u64, n_declared: usize, columns: &str, rows: &[&str]) -> String {
    let mut text = format!(
        "ITEM: TIMESTEP\n{timestep}\nITEM: NUMBER OF ATOMS\n{n_declared}\n\
         ITEM: BOX BOUNDS pp pp pp\n0 10\n0 20\n0 30\nITEM: ATOMS {columns}\n"
    );
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

#[test]
fn five_column_dump_returns_every_atom_and_the_timestep_of_frame_one() {
    let first = dump_frame(
        100,
        3,
        "id type x y z",
        &["1 1 0.5 1.5 2.5", "2 2 3.0 4.0 5.0", "3 1 6.0 7.0 8.0"],
    );
    let second = dump_frame(
        200,
        3,
        "id type x y z",
        &["1 1 9 9 9", "2 2 8 8 8", "3 1 7 7 7"],
    );
    let stream = format!("{first}{second}");
    let (timestep, atoms) = LammpsDumpReader::read_frame(Cursor::new(stream)).unwrap();
    assert_eq!(timestep, 100, "the timestep of the first frame");
    assert_eq!(atoms.len(), 3);
    let want = [
        (1, 1, [0.5, 1.5, 2.5]),
        (2, 2, [3.0, 4.0, 5.0]),
        (3, 1, [6.0, 7.0, 8.0]),
    ];
    for (a, (id, type_id, position)) in atoms.iter().zip(want) {
        assert_eq!((a.id, a.type_id), (id, type_id));
        assert_vec(a.position, position);
        assert_vec(a.velocity, [0.0, 0.0, 0.0]);
    }
}

#[test]
fn eight_column_dump_keeps_velocities_and_ignores_later_frames() {
    let atoms = vec![atom(1, 1, 1.0), atom(2, 3, 2.0)];
    let mut stream = Vec::new();
    LammpsDumpWriter::write_frame(&mut stream, 7, [[0.0, 10.0]; 3], &atoms).unwrap();
    LammpsDumpWriter::write_frame(&mut stream, 8, [[0.0, 10.0]; 3], &[atom(9, 9, 9.0)]).unwrap();
    let (timestep, read) = LammpsDumpReader::read_frame(Cursor::new(stream)).unwrap();
    assert_eq!(timestep, 7);
    assert_eq!(read.len(), 2);
    for (got, want) in read.iter().zip(&atoms) {
        assert_eq!((got.id, got.type_id), (want.id, want.type_id));
        assert_vec(
            got.position,
            [want.position.x, want.position.y, want.position.z],
        );
        assert_vec(
            got.velocity,
            [want.velocity.x, want.velocity.y, want.velocity.z],
        );
    }
}

#[test]
fn columns_are_found_by_name_not_by_position() {
    let text = dump_frame(
        5,
        1,
        "id mol type vx vy vz x y z",
        &["7 3 2 0.1 0.2 0.3 1.0 2.0 3.0"],
    );
    let (_, atoms) = LammpsDumpReader::read_frame(Cursor::new(text)).unwrap();
    assert_eq!(atoms.len(), 1);
    assert_eq!((atoms[0].id, atoms[0].type_id), (7, 2));
    assert_vec(atoms[0].position, [1.0, 2.0, 3.0]);
    assert_vec(atoms[0].velocity, [0.1, 0.2, 0.3]);
}

#[test]
fn units_and_time_blocks_in_front_of_the_frame_are_skipped() {
    // `dump_modify units yes time yes` puts these two blocks before the frame; the
    // reader this one replaced ignored them, so they must not turn into an error.
    let frame = dump_frame(
        40,
        1,
        "id type x y z vx vy vz",
        &["3 2 1.0 2.0 3.0 0.5 0.25 0.125"],
    );
    let text = format!("ITEM: UNITS\nlj\nITEM: TIME\n0.2\n{frame}");
    let (timestep, atoms) = LammpsDumpReader::read_frame(Cursor::new(text)).unwrap();
    assert_eq!(timestep, 40);
    assert_eq!(atoms.len(), 1);
    assert_eq!((atoms[0].id, atoms[0].type_id), (3, 2));
    assert_vec(atoms[0].position, [1.0, 2.0, 3.0]);
    assert_vec(atoms[0].velocity, [0.5, 0.25, 0.125]);
}

#[test]
fn scaled_coordinates_become_cartesian() {
    let text = dump_frame(0, 1, "id type xs ys zs", &["1 1 0.5 0.25 0.1"]);
    let (_, atoms) = LammpsDumpReader::read_frame(Cursor::new(text)).unwrap();
    // The box is 10 x 20 x 30 starting at the origin.
    assert_vec(atoms[0].position, [5.0, 5.0, 3.0]);
}

#[test]
fn stream_without_a_frame_is_a_parse_error() {
    for text in ["", "\n\n  \n"] {
        let message = parse_message(LammpsDumpReader::read_frame(Cursor::new(text)));
        assert!(message.contains("no frame"), "{message}");
    }
}

#[test]
fn malformed_frames_are_parse_errors() {
    // Two rows for three declared atoms.
    let truncated = dump_frame(1, 3, "id type x y z", &["1 1 0 0 0", "2 1 1 1 1"]);
    parse_message(LammpsDumpReader::read_frame(Cursor::new(truncated)));
    // A value that is not a number.
    let garbage = dump_frame(1, 1, "id type x y z", &["1 1 oops 0 0"]);
    let message = parse_message(LammpsDumpReader::read_frame(Cursor::new(garbage)));
    assert!(message.contains("oops"), "{message}");
    // No position columns at all.
    let no_positions = dump_frame(1, 1, "id type q", &["1 1 0.5"]);
    parse_message(LammpsDumpReader::read_frame(Cursor::new(no_positions)));
    // A row shorter than the header.
    let short_row = dump_frame(1, 1, "id type x y z", &["1 1 0.0"]);
    parse_message(LammpsDumpReader::read_frame(Cursor::new(short_row)));
}

#[test]
fn a_huge_declared_atom_count_does_not_allocate() {
    let text = dump_frame(1, 999_999_999_999, "id type x y z", &["1 1 0 0 0"]);
    parse_message(LammpsDumpReader::read_frame(Cursor::new(text)));
}

/// Serves `good` and then fails with a custom I/O error.
struct FailAfter {
    good: Cursor<Vec<u8>>,
}

impl Read for FailAfter {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self.good.read(buf)? {
            0 => Err(std::io::Error::other("device vanished")),
            n => Ok(n),
        }
    }
}

#[test]
fn read_failures_are_reported_as_io_errors() {
    let intact = dump_frame(1, 2, "id type x y z", &["1 1 0 0 0", "2 1 1 1 1"]);
    // Cut inside the atom rows: the stream breaks before the frame is complete.
    let cut = intact.len() - 6;
    let reader = FailAfter {
        good: Cursor::new(intact.as_bytes()[..cut].to_vec()),
    };
    match LammpsDumpReader::read_frame(reader) {
        Err(Error::Io(error)) => {
            assert_eq!(error.kind(), std::io::ErrorKind::Other);
            assert!(error.to_string().contains("device vanished"), "{error}");
        }
        other => panic!("expected Error::Io, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// binary dump reader
// ---------------------------------------------------------------------------

/// Bytes of a frame with two atoms.
fn binary_frame() -> Vec<u8> {
    let mut buf = Vec::new();
    let atoms = [atom(1, 1, 1.0), atom(2, 2, 2.0)];
    LammpsBinaryDumpWriter::write_frame(
        &mut buf,
        42,
        [[0.0, 10.0], [0.0, 20.0], [0.0, 30.0]],
        &atoms,
    )
    .unwrap();
    buf
}

/// Offset of the atom count in a binary frame (after the timestep).
const BINARY_COUNT: std::ops::Range<usize> = 8..16;
/// Offset of the first atom record in a binary frame.
const BINARY_FIRST_ATOM: usize = 8 + 8 + 6 * 8;

#[test]
fn binary_frame_round_trips_every_field() {
    let buf = binary_frame();
    let (timestep, bounds, atoms, used) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
    assert_eq!(timestep, 42);
    assert_eq!(bounds, [[0.0, 10.0], [0.0, 20.0], [0.0, 30.0]]);
    assert_eq!(used, buf.len());
    assert_eq!(atoms.len(), 2);
    assert_eq!((atoms[1].id, atoms[1].type_id), (2, 2));
    assert_vec(atoms[1].position, [2.0, 2.5, 3.0]);
    assert_vec(atoms[1].velocity, [-2.0, 0.5, 2.0]);
}

#[test]
fn binary_negative_atom_count_is_an_error() {
    for count in [-1_i64, i64::MIN] {
        let mut buf = binary_frame();
        buf[BINARY_COUNT].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsBinaryDumpReader::read_frame(&buf, 0));
        assert!(message.contains("negative"), "{message}");
    }
}

#[test]
fn binary_atom_count_beyond_the_data_is_an_error_without_allocating() {
    for count in [3_i64, 1_000_000, 1_000_000_000_000_000, i64::MAX] {
        let mut buf = binary_frame();
        buf[BINARY_COUNT].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsBinaryDumpReader::read_frame(&buf, 0));
        assert!(message.contains("exceeds"), "{count}: {message}");
    }
}

#[test]
fn binary_smaller_count_reads_that_many_atoms_and_reports_the_bytes_used() {
    let mut buf = binary_frame();
    buf[BINARY_COUNT].copy_from_slice(&1_i64.to_le_bytes());
    let (_, _, atoms, used) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
    assert_eq!(atoms.len(), 1);
    assert_eq!(used, BINARY_FIRST_ATOM + 56);
    buf[BINARY_COUNT].copy_from_slice(&0_i64.to_le_bytes());
    let (_, _, atoms, used) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
    assert!(atoms.is_empty());
    assert_eq!(used, BINARY_FIRST_ATOM);
}

#[test]
fn binary_offsets_past_or_near_the_end_are_errors_not_overflows() {
    let buf = binary_frame();
    for offset in [
        buf.len(),
        buf.len() + 1,
        usize::MAX - 7,
        usize::MAX - 3,
        usize::MAX,
    ] {
        parse_message(LammpsBinaryDumpReader::read_frame(&buf, offset));
    }
}

#[test]
fn binary_truncated_frames_are_errors() {
    let buf = binary_frame();
    for cut in [
        0,
        7,
        15,
        BINARY_FIRST_ATOM - 1,
        BINARY_FIRST_ATOM + 55,
        buf.len() - 1,
    ] {
        parse_message(LammpsBinaryDumpReader::read_frame(&buf[..cut], 0));
    }
}

#[test]
fn binary_negative_ids_and_types_are_errors() {
    let mut buf = binary_frame();
    buf[BINARY_FIRST_ATOM..BINARY_FIRST_ATOM + 4].copy_from_slice(&(-5_i32).to_le_bytes());
    let message = parse_message(LammpsBinaryDumpReader::read_frame(&buf, 0));
    assert!(message.contains("negative atom id"), "{message}");
    let mut buf = binary_frame();
    buf[BINARY_FIRST_ATOM + 4..BINARY_FIRST_ATOM + 8].copy_from_slice(&i32::MIN.to_le_bytes());
    let message = parse_message(LammpsBinaryDumpReader::read_frame(&buf, 0));
    assert!(message.contains("negative atom type"), "{message}");
}

#[test]
fn binary_frames_are_read_one_after_another_through_the_offset() {
    let mut stream = binary_frame();
    let mut second = Vec::new();
    LammpsBinaryDumpWriter::write_frame(&mut second, 43, [[0.0, 1.0]; 3], &[atom(9, 4, 9.0)])
        .unwrap();
    stream.extend_from_slice(&second);
    let (first_step, _, first_atoms, used) =
        LammpsBinaryDumpReader::read_frame(&stream, 0).unwrap();
    assert_eq!((first_step, first_atoms.len()), (42, 2));
    let (second_step, _, second_atoms, used_second) =
        LammpsBinaryDumpReader::read_frame(&stream, used).unwrap();
    assert_eq!((second_step, second_atoms.len()), (43, 1));
    assert_eq!(used + used_second, stream.len());
    assert_eq!(second_atoms[0].id, 9);
}

#[test]
fn binary_writer_rejects_ids_that_do_not_fit_the_format() {
    let too_big = i32::MAX as usize + 1;
    for bad in [atom(too_big, 1, 0.0), atom(1, too_big, 0.0)] {
        let mut buf = Vec::new();
        let result = LammpsBinaryDumpWriter::write_frame(&mut buf, 0, [[0.0, 1.0]; 3], &[bad]);
        assert!(matches!(result, Err(Error::General(_))), "{result:?}");
        assert!(buf.is_empty(), "nothing is written before the check fails");
    }
    // The largest representable id still round-trips.
    let mut buf = Vec::new();
    let max_id = i32::MAX as usize;
    LammpsBinaryDumpWriter::write_frame(&mut buf, 0, [[0.0, 1.0]; 3], &[atom(max_id, 1, 0.0)])
        .unwrap();
    let (_, _, atoms, _) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
    assert_eq!(atoms[0].id, max_id);
}

// ---------------------------------------------------------------------------
// restart reader
// ---------------------------------------------------------------------------

/// A restart file with two masses and two atoms.
fn restart_bytes() -> Vec<u8> {
    let mut buf = Vec::new();
    LammpsRestartWriter::write(
        &mut buf,
        1000,
        &[12.0, 16.0],
        &[atom(1, 1, 1.0), atom(2, 2, 2.0)],
    )
    .unwrap();
    buf
}

/// Where the counts sit in a restart file: magic, version, timestep, then
/// `n_atoms: i64` and `n_types: i32`.
const RESTART_ATOMS: std::ops::Range<usize> = 16..24;
const RESTART_TYPES: std::ops::Range<usize> = 24..28;
const RESTART_MASSES: usize = 28;

#[test]
fn restart_round_trips_masses_atoms_and_velocities() {
    let buf = restart_bytes();
    let (timestep, masses, atoms) = LammpsRestartReader::read(&buf).unwrap();
    assert_eq!(timestep, 1000);
    assert_eq!(masses, vec![12.0, 16.0]);
    assert_eq!(atoms.len(), 2);
    assert_vec(atoms[0].position, [1.0, 1.5, 2.0]);
    assert_vec(atoms[0].velocity, [-1.0, 0.25, 2.0]);
    assert_eq!((atoms[1].id, atoms[1].type_id), (2, 2));
}

#[test]
fn restart_negative_counts_are_errors() {
    for count in [-1_i64, i64::MIN] {
        let mut buf = restart_bytes();
        buf[RESTART_ATOMS].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsRestartReader::read(&buf));
        assert!(message.contains("negative atom count"), "{message}");
    }
    for count in [-1_i32, i32::MIN] {
        let mut buf = restart_bytes();
        buf[RESTART_TYPES].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsRestartReader::read(&buf));
        assert!(message.contains("negative atom type count"), "{message}");
    }
}

#[test]
fn restart_counts_beyond_the_data_are_errors_without_allocating() {
    for count in [3_i64, 1_000_000_000_000_000, i64::MAX] {
        let mut buf = restart_bytes();
        buf[RESTART_ATOMS].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsRestartReader::read(&buf));
        assert!(message.contains("exceeds"), "{count}: {message}");
    }
    for count in [3_i32, 1_000_000, i32::MAX] {
        let mut buf = restart_bytes();
        buf[RESTART_TYPES].copy_from_slice(&count.to_le_bytes());
        let message = parse_message(LammpsRestartReader::read(&buf));
        assert!(message.contains("exceeds"), "{count}: {message}");
    }
}

#[test]
fn restart_mass_table_and_atoms_must_fit_together() {
    // Each count alone fits the 2 * 8 + 2 * 56 bytes after the header, but five
    // masses (40 bytes) leave only 88 bytes, room for one atom, not two.
    let mut buf = restart_bytes();
    buf[RESTART_TYPES].copy_from_slice(&5_i32.to_le_bytes());
    let message = parse_message(LammpsRestartReader::read(&buf));
    assert!(message.contains("atom count 2 exceeds"), "{message}");
    // The same table with one atom fewer is consistent.
    buf[RESTART_ATOMS].copy_from_slice(&1_i64.to_le_bytes());
    let (_, masses, atoms) = LammpsRestartReader::read(&buf).unwrap();
    assert_eq!((masses.len(), atoms.len()), (5, 1));
}

#[test]
fn restart_truncated_files_are_errors() {
    let buf = restart_bytes();
    for cut in [3, 4, 20, 27, RESTART_MASSES + 15, buf.len() - 1] {
        parse_message(LammpsRestartReader::read(&buf[..cut]));
    }
}

#[test]
fn restart_negative_ids_are_errors() {
    let mut buf = restart_bytes();
    let first_atom = RESTART_MASSES + 2 * 8;
    buf[first_atom..first_atom + 4].copy_from_slice(&(-1_i32).to_le_bytes());
    let message = parse_message(LammpsRestartReader::read(&buf));
    assert!(message.contains("negative atom id"), "{message}");
}

#[test]
fn restart_writer_rejects_ids_that_do_not_fit_the_format() {
    let mut buf = Vec::new();
    let result =
        LammpsRestartWriter::write(&mut buf, 0, &[1.0], &[atom(i32::MAX as usize + 1, 1, 0.0)]);
    assert!(matches!(result, Err(Error::General(_))), "{result:?}");
    assert!(buf.is_empty(), "nothing is written before the check fails");
}

// ---------------------------------------------------------------------------
// corrupted data
// ---------------------------------------------------------------------------

#[test]
fn corrupted_binary_data_never_panics_and_never_reads_more_than_it_holds() {
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    let (mut accepted, mut rejected) = (0_usize, 0_usize);
    for _ in 0..4000 {
        let mut dump = binary_frame();
        let mut restart = restart_bytes();
        for bytes in [&mut dump, &mut restart] {
            for _ in 0..1 + rng.below(3) {
                let at = rng.below(bytes.len());
                bytes[at] = rng.next_u64() as u8;
            }
        }
        // Damage the counts outright as well: they are what sizes the allocations.
        if rng.below(2) == 0 {
            dump[BINARY_COUNT].copy_from_slice(&rng.next_u64().to_le_bytes());
            restart[RESTART_ATOMS].copy_from_slice(&rng.next_u64().to_le_bytes());
            restart[RESTART_TYPES].copy_from_slice(&(rng.next_u64() as u32).to_le_bytes());
        }
        let offset = if rng.below(4) == 0 {
            rng.below(dump.len() + 8)
        } else {
            0
        };
        match LammpsBinaryDumpReader::read_frame(&dump, offset) {
            Ok((_, _, atoms, used)) => {
                accepted += 1;
                assert!(atoms.len() * 56 <= dump.len(), "{} atoms", atoms.len());
                assert!(offset + used <= dump.len());
            }
            Err(_) => rejected += 1,
        }
        match LammpsRestartReader::read(&restart) {
            Ok((_, masses, atoms)) => {
                accepted += 1;
                assert!(masses.len() * 8 + atoms.len() * 56 <= restart.len());
            }
            Err(_) => rejected += 1,
        }
        let cut = rng.below(restart.len() + 1);
        let _ = LammpsRestartReader::read(&restart[..cut]);
    }
    // The corruption reaches both outcomes, so the checks above are not vacuous.
    assert!(accepted > 1000, "{accepted} reads were accepted");
    assert!(rejected > 1000, "{rejected} reads were rejected");
}
