// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of reading frames from a stream: the optional items in front of a
//! frame, consecutive frames, line endings, and what may stand between frames.

use super::{LammpsDumpFrame, LammpsDumpReader};
use std::io::Cursor;

/// The `ITEM: UNITS` and `ITEM: TIME` blocks of `dump_modify units yes time yes`.
const UNITS_AND_TIME: &str = "ITEM: UNITS\nlj\nITEM: TIME\n0.005\n";

/// One frame in a 10 x 10 x 10 box with the columns `id type x y z`.
fn frame_text(timestep: u64, rows: &[&str]) -> String {
    let mut text = format!(
        "ITEM: TIMESTEP\n{timestep}\nITEM: NUMBER OF ATOMS\n{}\n\
         ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\nITEM: ATOMS id type x y z\n",
        rows.len()
    );
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

fn read_one(text: &str) -> LammpsDumpFrame {
    LammpsDumpReader::read_frame(Cursor::new(text))
        .expect("the text holds a frame")
        .expect("the frame parses")
}

fn timesteps(frames: &[LammpsDumpFrame]) -> Vec<u64> {
    frames.iter().map(|f| f.timestep).collect()
}

#[test]
fn units_and_time_blocks_in_front_of_a_frame_are_skipped() {
    let frame = frame_text(5, &["1 1 0.5 1.5 2.5", "2 2 3.0 4.0 5.0"]);
    let variants = [
        format!("{UNITS_AND_TIME}{frame}"),
        format!("ITEM: TIME\n0.005\nITEM: UNITS\nmetal\n{frame}"),
        format!("ITEM: UNITS\nreal\n{frame}"),
        format!("ITEM: TIME\n12.5\n{frame}"),
        format!("\n\n{UNITS_AND_TIME}{frame}"),
    ];
    for text in variants {
        let parsed = read_one(&text);
        assert_eq!(parsed.timestep, 5, "{text}");
        assert_eq!(parsed.atom_ids, vec![1, 2]);
        assert_eq!(parsed.positions, vec![[0.5, 1.5, 2.5], [3.0, 4.0, 5.0]]);
    }
}

#[test]
fn a_block_without_its_value_line_is_an_error() {
    for text in [
        "ITEM: UNITS\n",
        "ITEM: TIME",
        "ITEM: UNITS\nlj\nITEM: TIME\n",
    ] {
        let message = LammpsDumpReader::read_frame(Cursor::new(text))
            .expect("there is something to read")
            .expect_err(text);
        assert!(message.contains("unexpected EOF"), "{text}: {message}");
    }
}

#[test]
fn an_unknown_item_in_front_of_a_frame_is_an_error() {
    let text = format!("ITEM: FOO\n1\n{}", frame_text(0, &["1 1 0 0 0"]));
    let message = LammpsDumpReader::read_frame(Cursor::new(text))
        .expect("there is something to read")
        .unwrap_err();
    assert!(
        message.contains("expected 'ITEM: TIMESTEP', got 'ITEM: FOO'"),
        "{message}"
    );
}

#[test]
fn read_frame_stops_after_the_frame_so_the_stream_can_be_read_on() {
    let text = [
        frame_text(0, &["1 1 0 0 0"]),
        format!(
            "{UNITS_AND_TIME}{}",
            frame_text(10, &["1 1 1 1 1", "2 1 2 2 2"])
        ),
        "\n".to_string(),
    ]
    .concat();
    let mut stream = Cursor::new(text.as_bytes());
    let first = LammpsDumpReader::read_frame(&mut stream).unwrap().unwrap();
    assert_eq!((first.timestep, first.n_atoms), (0, 1));
    let second = LammpsDumpReader::read_frame(&mut stream).unwrap().unwrap();
    assert_eq!((second.timestep, second.n_atoms), (10, 2));
    assert!(LammpsDumpReader::read_frame(&mut stream).is_none());
}

#[test]
fn read_all_frames_reads_every_frame_whatever_stands_in_front_of_it() {
    let text = [
        UNITS_AND_TIME.to_string(),
        frame_text(0, &["1 1 0 0 0", "2 1 1 1 1"]),
        "\n".to_string(),
        "ITEM: TIME\n0.5\n".to_string(),
        frame_text(100, &["1 1 0.1 0.1 0.1", "2 1 1.1 1.1 1.1"]),
        frame_text(200, &["1 1 0.2 0.2 0.2", "2 1 1.2 1.2 1.2"]),
        "\n\n".to_string(),
    ]
    .concat();
    let frames = LammpsDumpReader::read_all_frames(&text).unwrap();
    assert_eq!(timesteps(&frames), vec![0, 100, 200]);
    assert_eq!(frames[2].positions, vec![[0.2, 0.2, 0.2], [1.2, 1.2, 1.2]]);
    assert_eq!(
        LammpsDumpReader::extract_timesteps(&text),
        vec![0, 100, 200]
    );
}

#[test]
fn read_all_frames_reads_windows_line_endings() {
    let text = [
        frame_text(1, &["1 1 0 0 0"]),
        frame_text(2, &["1 1 1 1 1"]),
        frame_text(3, &["1 1 2 2 2"]),
    ]
    .concat()
    .replace('\n', "\r\n");
    let frames = LammpsDumpReader::read_all_frames(&text).unwrap();
    assert_eq!(timesteps(&frames), vec![1, 2, 3]);
    assert_eq!(frames[1].positions, vec![[1.0, 1.0, 1.0]]);
}

#[test]
fn read_all_frames_without_content_has_no_frames() {
    for text in ["", "\n", "  \n\t\n"] {
        assert!(LammpsDumpReader::read_all_frames(text).unwrap().is_empty());
    }
}

#[test]
fn read_all_frames_rejects_garbage_between_frames() {
    let text = [
        frame_text(0, &["1 1 0 0 0"]),
        "stray line\n".to_string(),
        frame_text(1, &["1 1 1 1 1"]),
    ]
    .concat();
    let message = LammpsDumpReader::read_all_frames(&text).unwrap_err();
    assert!(
        message.contains("expected 'ITEM: TIMESTEP', got 'stray line'"),
        "{message}"
    );
    // So are atom rows beyond the declared count.
    let extra = [
        frame_text(0, &["1 1 0 0 0"]),
        "2 1 1 1 1\n".to_string(),
        frame_text(1, &["1 1 1 1 1"]),
    ]
    .concat();
    assert!(LammpsDumpReader::read_all_frames(&extra).is_err());
}

#[test]
fn read_all_frames_rejects_a_frame_that_is_cut_short() {
    let last_row = "2 1 1 1 1\n";
    let whole = [
        frame_text(0, &["1 1 0 0 0"]),
        frame_text(7, &["1 1 0 0 0", "2 1 1 1 1"]),
    ]
    .concat();
    assert_eq!(LammpsDumpReader::read_all_frames(&whole).unwrap().len(), 2);
    // The second frame declares two atoms but only one row is left.
    let cut = &whole[..whole.len() - last_row.len()];
    let message = LammpsDumpReader::read_all_frames(cut).unwrap_err();
    assert!(message.contains("unexpected EOF"), "{message}");
    // The same file cut inside the last row has a row that is too short.
    let torn = &whole[..whole.len() - 5];
    assert!(LammpsDumpReader::read_all_frames(torn).is_err());
}
