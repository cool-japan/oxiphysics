//! Tests of the DCD reader and writer.
//!
//! The reader fixtures are laid out byte by byte from the CHARMM/NAMD DCD
//! specification (absolute file offsets below), independently of the code
//! under test.  A writer-to-reader round trip alone cannot catch a layout
//! mistake that both sides share, which is how the original header bugs
//! (the time step in the fixed-atom slot, the unit-cell record read as `x`)
//! went unnoticed.

use std::io::{self, Cursor};

use super::ParticleFrame;
use super::dcd::{
    DcdHeader, DcdReader, DcdTrajectory, DcdUnitCell, DcdWriter, read_dcd_trajectory,
};

// Absolute offsets in a file with a single 80-byte title; `ICNTRL[k]` is the
// 32-bit word at `8 + 4 k`.
const OFF_MAGIC: usize = 4;
const OFF_NSET: usize = 8;
const OFF_ISTART: usize = 12;
const OFF_NSAVC: usize = 16;
const OFF_NSTEP: usize = 20;
const OFF_NAMNF: usize = 40;
const OFF_DELTA: usize = 44;
const OFF_CELL_FLAG: usize = 48;
const OFF_4D_FLAG: usize = 52;
const OFF_VERSION: usize = 84;
const OFF_HEADER_CLOSE: usize = 88;
const OFF_TITLES: usize = 92;
/// Title record (4 + 4 + 80 + 4 bytes) and atom count record (12 bytes).
const OFF_DATA: usize = OFF_TITLES + 92 + 12;
/// Size of the unit-cell record including its markers.
const CELL_RECORD: usize = 8 + 48;

// ── fixture builder ─────────────────────────────────────────────────────────

#[derive(Clone)]
struct FrameSpec {
    /// Raw unit-cell record `[A, cos(gamma), B, cos(beta), cos(alpha), C]`.
    cell: Option<[f64; 6]>,
    xyz: Vec<[f32; 3]>,
}

#[derive(Clone)]
struct Spec {
    big_endian: bool,
    magic: [u8; 4],
    nset: u32,
    istart: u32,
    nsavc: u32,
    nstep: u32,
    namnf: u32,
    /// Zero selects the X-PLOR layout (`DELTA` is a double over words 9, 10).
    charmm_version: u32,
    delta: f64,
    cell_flag: u32,
    four_d_flag: u32,
    titles: Vec<String>,
    n_atoms: u32,
    frames: Vec<FrameSpec>,
}

impl Spec {
    fn charmm(n_atoms: u32) -> Self {
        Self {
            big_endian: false,
            magic: *b"CORD",
            nset: 0,
            istart: 0,
            nsavc: 1,
            nstep: 0,
            namnf: 0,
            charmm_version: 24,
            delta: 0.5,
            cell_flag: 0,
            four_d_flag: 0,
            titles: vec!["fixture".to_string()],
            n_atoms,
            frames: Vec::new(),
        }
    }

    fn xplor(n_atoms: u32) -> Self {
        Self {
            charmm_version: 0,
            delta: 0.001,
            ..Self::charmm(n_atoms)
        }
    }

    /// The same file with `n` frames and `NSET` set to match.
    fn with_frames(mut self, frames: Vec<FrameSpec>) -> Self {
        self.nset = frames.len() as u32;
        self.frames = frames;
        self
    }

    fn word(&self, value: u32) -> [u8; 4] {
        if self.big_endian {
            value.to_be_bytes()
        } else {
            value.to_le_bytes()
        }
    }

    fn double(&self, value: f64) -> [u8; 8] {
        if self.big_endian {
            value.to_be_bytes()
        } else {
            value.to_le_bytes()
        }
    }

    fn put_word(&self, block: &mut [u8], index: usize, value: u32) {
        block[4 * index..4 * index + 4].copy_from_slice(&self.word(value));
    }

    /// One coordinate record: marker, values, marker.
    fn record(&self, values: impl Iterator<Item = f32>, n: usize) -> Vec<u8> {
        let mut out = self.word(4 * n as u32).to_vec();
        for v in values {
            out.extend_from_slice(&if self.big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            });
        }
        out.extend_from_slice(&self.word(4 * n as u32));
        out
    }

    fn build(&self) -> Vec<u8> {
        let mut out = Vec::new();
        // Record 1: 84 | CORD | ICNTRL[20] | 84
        out.extend_from_slice(&self.word(84));
        out.extend_from_slice(&self.magic);
        let mut icntrl = [0u8; 80];
        self.put_word(&mut icntrl, 0, self.nset);
        self.put_word(&mut icntrl, 1, self.istart);
        self.put_word(&mut icntrl, 2, self.nsavc);
        self.put_word(&mut icntrl, 3, self.nstep);
        self.put_word(&mut icntrl, 8, self.namnf);
        if self.charmm_version != 0 {
            self.put_word(&mut icntrl, 9, (self.delta as f32).to_bits());
            self.put_word(&mut icntrl, 10, self.cell_flag);
            self.put_word(&mut icntrl, 11, self.four_d_flag);
        } else {
            icntrl[36..44].copy_from_slice(&self.double(self.delta));
        }
        self.put_word(&mut icntrl, 19, self.charmm_version);
        out.extend_from_slice(&icntrl);
        out.extend_from_slice(&self.word(84));
        // Record 2: titles
        let title_len = 4 + 80 * self.titles.len() as u32;
        out.extend_from_slice(&self.word(title_len));
        out.extend_from_slice(&self.word(self.titles.len() as u32));
        for title in &self.titles {
            let mut line = [b' '; 80];
            line[..title.len()].copy_from_slice(title.as_bytes());
            out.extend_from_slice(&line);
        }
        out.extend_from_slice(&self.word(title_len));
        // Record 3: number of atoms
        out.extend_from_slice(&self.word(4));
        out.extend_from_slice(&self.word(self.n_atoms));
        out.extend_from_slice(&self.word(4));
        // Frames
        let n = self.n_atoms as usize;
        for frame in &self.frames {
            if let Some(cell) = frame.cell {
                out.extend_from_slice(&self.word(48));
                for value in cell {
                    out.extend_from_slice(&self.double(value));
                }
                out.extend_from_slice(&self.word(48));
            }
            for axis in 0..3 {
                out.extend(self.record(frame.xyz.iter().map(|p| p[axis]), n));
            }
            if self.charmm_version != 0 && self.four_d_flag == 1 {
                // The fourth dimension carries values that must never leak
                // into the positions.
                out.extend(self.record((0..n).map(|i| 9000.0 + i as f32), n));
            }
        }
        out
    }
}

fn xyz(n: usize, base: f32) -> Vec<[f32; 3]> {
    (0..n)
        .map(|i| {
            let i = i as f32;
            [base + i, base + 10.0 + i, base + 20.0 + i]
        })
        .collect()
}

fn parse(bytes: &[u8]) -> io::Result<DcdTrajectory> {
    read_dcd_trajectory(Cursor::new(bytes.to_vec()), bytes.len() as u64)
}

fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([
        bytes[offset],
        bytes[offset + 1],
        bytes[offset + 2],
        bytes[offset + 3],
    ])
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn temp_path(tag: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("oxiphysics_dcd_{}_{tag}.dcd", std::process::id()))
}

fn make_frame(n: usize) -> ParticleFrame {
    let positions: Vec<[f32; 3]> = (0..n).map(|i| [i as f32, 0.0, 0.0]).collect();
    ParticleFrame::new(0, positions, vec![1u8; n])
}

/// A file with a unit cell: 3 atoms, 2 frames.
fn cell_spec() -> Spec {
    Spec {
        istart: 100,
        nsavc: 50,
        nstep: 200,
        cell_flag: 1,
        ..Spec::charmm(3)
    }
    .with_frames(vec![
        FrameSpec {
            cell: Some([10.0, 0.0, 20.0, 0.0, 0.0, 30.0]),
            xyz: xyz(3, 1.0),
        },
        FrameSpec {
            cell: Some([11.0, 0.5, 21.0, 0.0, 0.0, 31.0]),
            xyz: xyz(3, 100.5),
        },
    ])
}

fn assert_cell_spec_decoded(t: &DcdTrajectory, big_endian: bool) {
    let h = &t.header;
    assert_eq!(h.n_frames, 2);
    assert_eq!(h.n_atoms, 3);
    assert_eq!(h.first_step, 100);
    assert_eq!(h.step_interval, 50);
    assert_eq!(h.n_steps, 200);
    assert_eq!(h.delta, 0.5);
    assert!(h.has_unit_cell);
    assert_eq!(h.charmm_version, 24);
    assert_eq!(h.big_endian, big_endian);
    assert_eq!(h.titles, vec!["fixture".to_string()]);

    assert_eq!(t.frames.len(), 2);
    assert_eq!(t.frames[0].positions, xyz(3, 1.0));
    assert_eq!(t.frames[1].positions, xyz(3, 100.5));
    assert_eq!(t.frames[0].n_particles, 3);
    assert_eq!(t.frames[1].timestep, 1);

    assert_eq!(t.unit_cells.len(), 2);
    assert_eq!(
        t.unit_cells[0],
        DcdUnitCell::orthorhombic([10.0, 20.0, 30.0])
    );
    assert_eq!(t.unit_cells[1].lengths, [11.0, 21.0, 31.0]);
    assert_eq!(t.unit_cells[1].angles[0], 90.0);
    assert_eq!(t.unit_cells[1].angles[1], 90.0);
    assert!((t.unit_cells[1].angles[2] - 60.0).abs() < 1e-12);
}

// ── reader: layout (F2258 / F2259) ──────────────────────────────────────────

#[test]
fn charmm_file_with_unit_cell_decodes_cell_and_coordinates() {
    // The 48-byte cell record sits in front of x.  The original reader took
    // it for the x coordinates, shifted y and z and misread the cell.
    let t = parse(&cell_spec().build()).unwrap();
    assert_cell_spec_decoded(&t, false);
}

#[test]
fn big_endian_file_is_read_with_swapped_markers_and_values() {
    let spec = Spec {
        big_endian: true,
        ..cell_spec()
    };
    let bytes = spec.build();
    // 84 swapped is 0x54000000.
    assert_eq!(&bytes[..4], &[0, 0, 0, 84]);
    let t = parse(&bytes).unwrap();
    assert_cell_spec_decoded(&t, true);
}

#[test]
fn read_all_and_read_header_on_a_file_with_unit_cells() {
    let path = temp_path("cell_file");
    std::fs::write(&path, cell_spec().build()).unwrap();
    let reader = DcdReader::new(&path);
    let (header, frames) = reader.read_all().unwrap();
    assert_eq!(header.n_frames, 2);
    assert_eq!(frames.len(), 2);
    assert_eq!(frames[0].positions, xyz(3, 1.0));
    assert_eq!(frames[1].positions, xyz(3, 100.5));
    let header_only = reader.read_header().unwrap();
    assert_eq!(header_only.n_frames, 2);
    assert!(header_only.has_unit_cell);
    assert_eq!(header_only.first_step, 100);
    let trajectory = reader.read_trajectory().unwrap();
    assert_eq!(trajectory.unit_cells.len(), 2);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn xplor_file_double_delta_is_not_mistaken_for_a_cell_flag() {
    // X-PLOR files (version word 0) keep DELTA as a double over words 9 and
    // 10, so word 10 is the upper half of a number, not the unit-cell flag.
    let spec = Spec::xplor(2).with_frames(vec![
        FrameSpec {
            cell: None,
            xyz: xyz(2, 1.0),
        },
        FrameSpec {
            cell: None,
            xyz: xyz(2, 5.0),
        },
    ]);
    let bytes = spec.build();
    assert_ne!(u32_at(&bytes, OFF_CELL_FLAG), 0, "fixture must set word 10");
    assert_eq!(u32_at(&bytes, OFF_VERSION), 0);
    let t = parse(&bytes).unwrap();
    assert!(!t.header.has_unit_cell);
    assert_eq!(t.header.charmm_version, 0);
    assert_eq!(t.header.delta, 0.001);
    assert!(t.unit_cells.is_empty());
    assert_eq!(t.frames.len(), 2);
    assert_eq!(t.frames[0].positions, xyz(2, 1.0));
    assert_eq!(t.frames[1].positions, xyz(2, 5.0));
}

#[test]
fn four_dimensional_frames_stay_aligned() {
    let spec = Spec {
        four_d_flag: 1,
        ..Spec::charmm(4)
    }
    .with_frames(vec![
        FrameSpec {
            cell: None,
            xyz: xyz(4, 1.0),
        },
        FrameSpec {
            cell: None,
            xyz: xyz(4, 50.0),
        },
    ]);
    let t = parse(&spec.build()).unwrap();
    assert_eq!(t.frames.len(), 2);
    assert_eq!(t.frames[0].positions, xyz(4, 1.0));
    assert_eq!(t.frames[1].positions, xyz(4, 50.0));
}

#[test]
fn titles_are_trimmed_and_all_returned() {
    let spec = Spec {
        titles: vec!["REMARKS FILENAME=x.dcd".to_string(), "second".to_string()],
        ..Spec::charmm(1)
    }
    .with_frames(vec![FrameSpec {
        cell: None,
        xyz: xyz(1, 0.0),
    }]);
    let t = parse(&spec.build()).unwrap();
    assert_eq!(
        t.header.titles,
        vec!["REMARKS FILENAME=x.dcd".to_string(), "second".to_string()]
    );
    assert_eq!(t.frames[0].positions, xyz(1, 0.0));
}

// ── reader: frame count policy ──────────────────────────────────────────────

fn three_frame_spec() -> Spec {
    Spec::charmm(2).with_frames(vec![
        FrameSpec {
            cell: None,
            xyz: xyz(2, 1.0),
        },
        FrameSpec {
            cell: None,
            xyz: xyz(2, 2.0),
        },
        FrameSpec {
            cell: None,
            xyz: xyz(2, 3.0),
        },
    ])
}

#[test]
fn header_that_was_never_finalised_does_not_hide_the_frames() {
    let spec = Spec {
        nset: 0,
        ..three_frame_spec()
    };
    let t = parse(&spec.build()).unwrap();
    assert_eq!(t.header.n_frames, 3);
    assert_eq!(t.frames.len(), 3);
    assert_eq!(t.frames[2].positions, xyz(2, 3.0));
}

#[test]
fn out_of_date_smaller_frame_count_is_corrected_by_the_file_size() {
    let spec = Spec {
        nset: 1,
        ..three_frame_spec()
    };
    let t = parse(&spec.build()).unwrap();
    assert_eq!(t.header.n_frames, 3);
    assert_eq!(t.frames.len(), 3);
}

#[test]
fn every_truncation_of_the_file_is_an_error_not_a_short_read() {
    let spec = cell_spec();
    let bytes = spec.build();
    assert!(parse(&bytes).is_ok());
    for cut in 0..bytes.len() {
        let err = parse(&bytes[..cut])
            .map(|t| t.frames.len())
            .expect_err(&format!("{cut} of {} bytes must not parse", bytes.len()));
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "cut {cut}: {err}");
    }
}

#[test]
fn file_ending_inside_the_second_frame_names_the_partial_frame() {
    let spec = Spec {
        nset: 0,
        ..three_frame_spec()
    };
    let bytes = spec.build();
    let frame = (bytes.len() - OFF_DATA) / 3;
    let err = parse(&bytes[..OFF_DATA + frame + 7]).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    assert!(err.to_string().contains("frame 1"), "{err}");
}

#[test]
fn header_announcing_more_frames_than_present_is_an_error() {
    let spec = Spec {
        nset: 9,
        ..three_frame_spec()
    };
    let err = parse(&spec.build()).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    assert!(err.to_string().contains("announces 9"), "{err}");
}

// ── reader: malformed and unsupported files ─────────────────────────────────

#[test]
fn header_record_of_the_wrong_length_is_invalid_data_and_never_panics() {
    // 4..15 used to underflow `rec_len - 4 - 12` and panic.
    for rec_len in [
        0u32,
        1,
        4,
        8,
        12,
        15,
        16,
        83,
        85,
        96,
        1000,
        0x5400_0001,
        u32::MAX,
    ] {
        let mut bytes = rec_len.to_le_bytes().to_vec();
        bytes.extend(vec![0u8; 400]);
        let err = parse(&bytes).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "rec_len {rec_len}");
        assert!(err.to_string().contains("84"), "{err}");
    }
}

#[test]
fn empty_file_and_a_lone_record_marker_are_unexpected_eof() {
    let err = parse(&[]).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    // A valid first marker, then nothing: the magic is missing.
    let err = parse(&84u32.to_le_bytes()).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    assert!(err.to_string().contains("the magic"), "{err}");
}

#[test]
fn wrong_magic_is_invalid_data() {
    let spec = Spec {
        magic: *b"XXXX",
        ..three_frame_spec()
    };
    let bytes = spec.build();
    assert_eq!(&bytes[OFF_MAGIC..OFF_MAGIC + 4], b"XXXX");
    let err = parse(&bytes).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("magic"), "{err}");
}

#[test]
fn velocity_trajectories_are_reported_as_unsupported() {
    let spec = Spec {
        magic: *b"VELO",
        ..three_frame_spec()
    };
    let err = parse(&spec.build()).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::Unsupported);
    assert!(err.to_string().contains("VELO"), "{err}");
}

#[test]
fn fixed_atoms_are_reported_as_unsupported() {
    let spec = Spec {
        namnf: 2,
        ..three_frame_spec()
    };
    let bytes = spec.build();
    assert_eq!(u32_at(&bytes, OFF_NAMNF), 2);
    let err = parse(&bytes).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::Unsupported);
    assert!(err.to_string().contains("fixed atoms"), "{err}");
}

#[test]
fn record_markers_that_disagree_with_the_layout_are_invalid_data() {
    let good = cell_spec().build();
    let x_record = OFF_DATA + CELL_RECORD;
    let y_record = x_record + 8 + 12;
    let tampered = |offset: usize, value: u32| {
        let mut bytes = good.clone();
        put_u32(&mut bytes, offset, value);
        parse(&bytes).unwrap_err()
    };
    // Header closing marker, atom count markers, unit cell markers, and the
    // leading and trailing marker of a coordinate record.
    for (offset, what) in [
        (OFF_HEADER_CLOSE, "header closing marker"),
        (OFF_DATA - 12, "atom count marker"),
        (OFF_DATA - 4, "atom count closing marker"),
        (OFF_DATA, "unit cell marker"),
        (OFF_DATA + CELL_RECORD - 4, "unit cell closing marker"),
        (y_record, "y leading marker"),
        (y_record + 4 + 12, "y trailing marker"),
    ] {
        let err = tampered(offset, 40);
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{what}: {err}");
    }
    // A coordinate record for the wrong number of atoms.
    let err = tampered(y_record, 16);
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("expected 12"), "{err}");
}

#[test]
fn inconsistent_title_record_is_invalid_data() {
    let good = three_frame_spec().build();
    let tampered = |offset: usize, value: u32| {
        let mut bytes = good.clone();
        put_u32(&mut bytes, offset, value);
        parse(&bytes).unwrap_err()
    };
    // Length that is not 4 + 80 n, title count that disagrees with the
    // length, and a closing marker that differs from the opening one.
    for (offset, value) in [
        (OFF_TITLES, 85),
        (OFF_TITLES, 3),
        (OFF_TITLES + 4, 5),
        (OFF_TITLES + 4 + 4 + 80, 164),
    ] {
        let err = tampered(offset, value);
        assert_eq!(
            err.kind(),
            io::ErrorKind::InvalidData,
            "offset {offset}: {err}"
        );
    }
    // A title record longer than the file.
    let err = tampered(OFF_TITLES, 4 + 80 * 1000);
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
}

#[test]
fn absurd_atom_counts_are_rejected_before_anything_is_allocated() {
    // Does not fit the 32-bit record length.
    let err = parse(&Spec::charmm(u32::MAX).build()).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("atoms"), "{err}");

    // A frame of 10^9 atoms is 12 GB: the header announces frames the file
    // cannot hold.
    let big = Spec {
        nset: 3,
        ..Spec::charmm(1_000_000_000)
    };
    let err = parse(&big.build()).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);

    // Without announced frames and without data it is simply an empty file.
    let empty = Spec::charmm(1_000_000_000);
    let t = parse(&empty.build()).unwrap();
    assert_eq!(t.frames.len(), 0);
    assert_eq!(t.header.n_atoms, 1_000_000_000);

    // A few stray bytes after the header are a partial frame.
    let mut stray = empty.build();
    stray.extend_from_slice(&[0u8; 9]);
    let err = parse(&stray).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
}

#[test]
fn zero_atom_frames_are_supported() {
    let spec = Spec::charmm(0).with_frames(vec![
        FrameSpec {
            cell: None,
            xyz: Vec::new(),
        },
        FrameSpec {
            cell: None,
            xyz: Vec::new(),
        },
    ]);
    let t = parse(&spec.build()).unwrap();
    assert_eq!(t.frames.len(), 2);
    assert!(t.frames.iter().all(|f| f.positions.is_empty()));
}

// ── unit cell conversion ────────────────────────────────────────────────────

#[test]
fn unit_cell_record_in_cosines_and_in_degrees() {
    let cosines = DcdUnitCell::from_record([4.0, -0.5, 5.0, 0.0, 0.0, 6.0]);
    assert_eq!(cosines.lengths, [4.0, 5.0, 6.0]);
    assert_eq!(cosines.angles[0], 90.0);
    assert_eq!(cosines.angles[1], 90.0);
    assert!((cosines.angles[2] - 120.0).abs() < 1e-12);

    // Older NAMD versions store the angles themselves.
    let degrees = DcdUnitCell::from_record([4.0, 120.0, 5.0, 100.0, 90.0, 6.0]);
    assert_eq!(degrees.lengths, [4.0, 5.0, 6.0]);
    assert_eq!(degrees.angles, [90.0, 100.0, 120.0]);
}

#[test]
fn unit_cell_record_roundtrip() {
    let cell = DcdUnitCell {
        lengths: [3.0, 4.0, 5.0],
        angles: [80.0, 90.0, 110.0],
    };
    let record = cell.to_record();
    assert_eq!(record[0], 3.0);
    assert_eq!(record[2], 4.0);
    assert_eq!(record[5], 5.0);
    assert_eq!(record[3], 0.0, "cos(90 degrees) is stored as exactly zero");
    let back = DcdUnitCell::from_record(record);
    assert_eq!(back.lengths, cell.lengths);
    for (a, b) in back.angles.iter().zip(cell.angles) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }
}

// ── writer: raw bytes (F2260) ───────────────────────────────────────────────

#[test]
fn writer_header_follows_the_charmm_layout() {
    let path = temp_path("raw_header");
    let mut w = DcdWriter::new(&path, 3, 0.002_f32).unwrap();
    w.write_frame(&make_frame(3)).unwrap();
    w.write_frame(&make_frame(3)).unwrap();
    w.finalize().unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let _ = std::fs::remove_file(&path);

    assert_eq!(u32_at(&bytes, 0), 84);
    assert_eq!(&bytes[OFF_MAGIC..OFF_MAGIC + 4], b"CORD");
    assert_eq!(u32_at(&bytes, OFF_NSET), 2, "NSET");
    assert_eq!(u32_at(&bytes, OFF_ISTART), 0, "ISTART");
    assert_eq!(u32_at(&bytes, OFF_NSAVC), 1, "NSAVC");
    assert_eq!(
        u32_at(&bytes, OFF_NSTEP),
        2,
        "NSTEP is ISTART + NSET * NSAVC"
    );
    // Words 4..8 are unused and word 8 is NAMNF: no fixed atoms.  The old
    // writer put DELTA here, which made VMD and NAMD expect a fixed-atom block.
    for word in 4..=8 {
        assert_eq!(u32_at(&bytes, 8 + 4 * word), 0, "ICNTRL[{word}]");
    }
    assert_eq!(u32_at(&bytes, OFF_NAMNF), 0, "NAMNF");
    assert_eq!(u32_at(&bytes, OFF_DELTA), 0.002_f32.to_bits(), "DELTA");
    assert_eq!(u32_at(&bytes, OFF_CELL_FLAG), 0, "unit cell flag");
    assert_eq!(u32_at(&bytes, OFF_4D_FLAG), 0, "4D flag");
    for word in 12..19 {
        assert_eq!(u32_at(&bytes, 8 + 4 * word), 0, "ICNTRL[{word}]");
    }
    assert_eq!(u32_at(&bytes, OFF_VERSION), 24, "CHARMM version");
    assert_eq!(u32_at(&bytes, OFF_HEADER_CLOSE), 84);

    // Title record, atom count record, then two frames of three 12-byte
    // coordinate records (3 atoms, 4 bytes each, 8 bytes of markers).
    assert_eq!(u32_at(&bytes, OFF_TITLES), 84);
    assert_eq!(u32_at(&bytes, OFF_TITLES + 4), 1);
    assert_eq!(u32_at(&bytes, OFF_TITLES + 88), 84);
    assert_eq!(u32_at(&bytes, OFF_DATA - 12), 4);
    assert_eq!(u32_at(&bytes, OFF_DATA - 8), 3);
    assert_eq!(u32_at(&bytes, OFF_DATA - 4), 4);
    assert_eq!(bytes.len(), OFF_DATA + 2 * 3 * (8 + 12));
    assert_eq!(u32_at(&bytes, OFF_DATA), 12);
    // x record of the first frame: atom i has x = i.
    for atom in 0..3usize {
        let bits = u32_at(&bytes, OFF_DATA + 4 + 4 * atom);
        assert_eq!(f32::from_bits(bits), atom as f32);
    }
}

#[test]
fn writer_finalize_patches_nset_and_nstep() {
    let path = temp_path("finalize");
    let mut w = DcdWriter::new(&path, 2, 1.0e-3_f32).unwrap();
    let before = std::fs::read(&path).unwrap();
    assert_eq!(u32_at(&before, OFF_NSET), 0);
    assert_eq!(u32_at(&before, OFF_NSTEP), 0);
    for _ in 0..5 {
        w.write_frame(&make_frame(2)).unwrap();
    }
    // Before finalize the header still says zero frames.
    let unfinished = std::fs::read(&path).unwrap();
    assert_eq!(u32_at(&unfinished, OFF_NSET), 0);
    // ... which the reader copes with.
    let (header, frames) = DcdReader::new(&path).read_all().unwrap();
    assert_eq!((header.n_frames, frames.len()), (5, 5));
    w.finalize().unwrap();
    let after = std::fs::read(&path).unwrap();
    let _ = std::fs::remove_file(&path);
    assert_eq!(u32_at(&after, OFF_NSET), 5);
    assert_eq!(u32_at(&after, OFF_NSTEP), 5);
    assert_eq!(after.len(), before.len() + 5 * 3 * (8 + 8));
}

#[test]
fn written_file_is_read_back_with_the_header_fields() {
    let path = temp_path("fields");
    let mut w = DcdWriter::new(&path, 4, 2.0e-3_f32).unwrap();
    w.write_frame(&make_frame(4)).unwrap();
    w.finalize().unwrap();
    let header = DcdReader::new(&path).read_header().unwrap();
    let _ = std::fs::remove_file(&path);
    assert_eq!(header.n_frames, 1);
    assert_eq!(header.n_atoms, 4);
    assert_eq!(header.first_step, 0);
    assert_eq!(header.step_interval, 1);
    assert_eq!(header.n_steps, 1);
    assert_eq!(header.delta, f64::from(2.0e-3_f32));
    assert_eq!(header.charmm_version, 24);
    assert!(!header.has_unit_cell);
    assert!(!header.big_endian);
    assert_eq!(header.titles, vec!["Created by oxiphysics DcdWriter"]);
}

// ── writer: unit cells ──────────────────────────────────────────────────────

#[test]
fn unit_cell_writer_roundtrip() {
    let path = temp_path("cell_rt");
    let mut w = DcdWriter::with_unit_cell(&path, 3, 1.0e-3_f32).unwrap();
    let cells = [
        DcdUnitCell::orthorhombic([10.0, 20.0, 30.0]),
        DcdUnitCell {
            lengths: [11.0, 12.0, 13.0],
            angles: [80.0, 90.0, 120.0],
        },
    ];
    for (i, cell) in cells.iter().enumerate() {
        let frame = ParticleFrame::from_positions(i as u64, xyz(3, i as f32 * 10.0));
        w.write_frame_with_cell(&frame, cell).unwrap();
    }
    w.finalize().unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let trajectory = DcdReader::new(&path).read_trajectory().unwrap();
    let _ = std::fs::remove_file(&path);

    assert_eq!(u32_at(&bytes, OFF_CELL_FLAG), 1, "unit cell flag");
    assert_eq!(u32_at(&bytes, OFF_DATA), 48);
    assert_eq!(bytes.len(), OFF_DATA + 2 * (CELL_RECORD + 3 * (8 + 12)));
    assert!(trajectory.header.has_unit_cell);
    assert_eq!(trajectory.frames.len(), 2);
    assert_eq!(trajectory.frames[1].positions, xyz(3, 10.0));
    assert_eq!(trajectory.unit_cells.len(), 2);
    assert_eq!(trajectory.unit_cells[0], cells[0]);
    assert_eq!(trajectory.unit_cells[1].lengths, cells[1].lengths);
    for (a, b) in trajectory.unit_cells[1].angles.iter().zip(cells[1].angles) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }
}

#[test]
fn writer_modes_reject_the_wrong_kind_of_frame() {
    let plain_path = temp_path("plain_mode");
    let cell_path = temp_path("cell_mode");
    let mut plain = DcdWriter::new(&plain_path, 2, 1.0e-3_f32).unwrap();
    let mut with_cell = DcdWriter::with_unit_cell(&cell_path, 2, 1.0e-3_f32).unwrap();
    let cell = DcdUnitCell::orthorhombic([1.0, 1.0, 1.0]);
    let frame = make_frame(2);
    let err = plain.write_frame_with_cell(&frame, &cell).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    let err = with_cell.write_frame(&frame).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    // Nothing was appended by the rejected calls.
    plain.finalize().unwrap();
    with_cell.finalize().unwrap();
    assert_eq!(
        DcdReader::new(&plain_path).read_header().unwrap().n_frames,
        0
    );
    assert_eq!(
        DcdReader::new(&cell_path).read_header().unwrap().n_frames,
        0
    );
    let _ = std::fs::remove_file(&plain_path);
    let _ = std::fs::remove_file(&cell_path);
}

#[test]
fn writer_rejects_atom_counts_that_overflow_the_record_length() {
    let path = temp_path("too_many_atoms");
    let err = DcdWriter::new(&path, u32::MAX / 2, 1.0e-3_f32).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    assert!(!path.exists());
}

// ── writer/reader round trips (moved from the parent module) ────────────────

#[test]
fn test_dcd_writer_creates_file() {
    let path = std::env::temp_dir().join("oxiphysics_pf_dcd_create.dcd");
    let _ = std::fs::remove_file(&path);
    let mut w = DcdWriter::new(&path, 3, 2.0e-3_f32).unwrap();
    let frame = make_frame(3);
    w.write_frame(&frame).unwrap();
    w.finalize().unwrap();
    assert!(path.exists());
    let _ = std::fs::remove_file(&path);
}

#[test]
fn test_dcd_roundtrip_single_frame() {
    let path = std::env::temp_dir().join("oxiphysics_pf_dcd_rt1.dcd");
    let _ = std::fs::remove_file(&path);
    let frame = make_frame(4);
    let mut w = DcdWriter::new(&path, 4, 1.0e-3_f32).unwrap();
    w.write_frame(&frame).unwrap();
    w.finalize().unwrap();

    let (header, frames) = DcdReader::new(&path).read_all().unwrap();
    assert_eq!(header.n_frames, 1);
    assert_eq!(header.n_atoms, 4);
    assert_eq!(frames.len(), 1);
    for i in 0..4 {
        assert!(
            (frames[0].positions[i][0] - i as f32).abs() < 1e-5,
            "position mismatch at atom {i}"
        );
    }
    let _ = std::fs::remove_file(&path);
}

#[test]
fn test_dcd_roundtrip_multi_frame() {
    let path = std::env::temp_dir().join("oxiphysics_pf_dcd_rt3.dcd");
    let _ = std::fs::remove_file(&path);
    let n_atoms = 5u32;
    let mut w = DcdWriter::new(&path, n_atoms, 2.0e-3_f32).unwrap();
    for ts in 0..3u64 {
        let pos: Vec<[f32; 3]> = (0..n_atoms as usize)
            .map(|i| [ts as f32 + i as f32, 0.0, 0.0])
            .collect();
        let types = vec![0u8; n_atoms as usize];
        let frame = ParticleFrame::new(ts, pos, types);
        w.write_frame(&frame).unwrap();
    }
    w.finalize().unwrap();

    let (header, frames) = DcdReader::new(&path).read_all().unwrap();
    assert_eq!(header.n_frames, 3);
    assert_eq!(frames.len(), 3);
    // Verify frame 2, atom 0: x == 2.0
    assert!((frames[2].positions[0][0] - 2.0).abs() < 1e-5);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn dcd_header_default_is_empty() {
    let h = DcdHeader::default();
    assert_eq!(h.n_frames, 0);
    assert!(h.titles.is_empty());
    assert!(!h.has_unit_cell);
}
