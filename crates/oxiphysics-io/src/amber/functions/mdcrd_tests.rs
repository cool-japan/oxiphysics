//! Tests of [`parse_mdcrd_checked`], the MDCRD parser that reports defects.

use super::mdcrd::{parse_mdcrd, parse_mdcrd_checked};
use crate::fortran_fields::test_support::{SplitMix64, mutate_text};

/// One frame as `10F8.3` lines followed by an optional `3F8.3` box line.
fn frame_text(coords: &[f64], box_dims: Option<[f64; 3]>) -> String {
    let mut text = String::new();
    for row in coords.chunks(10) {
        for value in row {
            text.push_str(&format!("{value:8.3}"));
        }
        text.push('\n');
    }
    if let Some(dims) = box_dims {
        for value in dims {
            text.push_str(&format!("{value:8.3}"));
        }
        text.push('\n');
    }
    text
}

#[test]
fn touching_fields_parse_exactly_and_agree_with_the_lenient_parser() {
    let coords = [-123.456, 999.999, 1000.000, -100.000, 1.0, -999.999];
    let text = format!(
        "title\n{}{}",
        frame_text(&coords, Some([1000.0, 999.5, 10.0])),
        frame_text(&coords, Some([1000.0, 1000.0, 1000.0]))
    );
    assert!(text.contains("999.9991000.000-100.000"), "{text}");
    let checked = parse_mdcrd_checked(&text, 2, true).unwrap();
    assert_eq!(checked.len(), 2);
    assert_eq!(checked[0].coords, coords.to_vec());
    assert_eq!(checked[0].box_dims, Some([1000.0, 999.5, 10.0]));
    assert_eq!(checked[1].box_dims, Some([1000.0, 1000.0, 1000.0]));
    let lenient = parse_mdcrd(&text, 2, true);
    assert_eq!(lenient.len(), checked.len());
    for (a, b) in lenient.iter().zip(&checked) {
        assert_eq!(a.coords, b.coords);
        assert_eq!(a.box_dims, b.box_dims);
    }
}

#[test]
fn empty_file_and_blank_trailing_lines_are_fine() {
    assert!(parse_mdcrd_checked("", 3, false).unwrap().is_empty());
    assert!(parse_mdcrd_checked("title\n", 3, true).unwrap().is_empty());
    let text = format!("t\n{}\n  \n\n", frame_text(&[1.0, 2.0, 3.0], None));
    assert_eq!(parse_mdcrd_checked(&text, 1, false).unwrap().len(), 1);
}

#[test]
fn non_numeric_field_names_the_frame_and_the_line() {
    let text = "t\n   1.000   2.000   3.000\n   1.000   x.000   3.000\n";
    let message = parse_mdcrd_checked(text, 1, false).unwrap_err();
    assert!(message.contains("frame 1"), "{message}");
    assert!(message.contains("line 3"), "{message}");
    assert!(message.contains("x.000"), "{message}");
    // The lenient parser returns the frame before the defect and says nothing.
    assert_eq!(parse_mdcrd(text, 1, false).len(), 1);
}

#[test]
fn wrong_number_of_values_is_reported() {
    let text = "t\n   1.000   2.000   3.000   4.000\n";
    let message = parse_mdcrd_checked(text, 2, false).unwrap_err();
    assert!(message.contains("frame 0"), "{message}");
    assert!(message.contains("4 values, expected 6"), "{message}");
    assert!(parse_mdcrd(text, 2, false).is_empty());
}

#[test]
fn truncated_trailing_frame_is_reported() {
    // 5 atoms = 15 values = two lines per frame; the second frame has one.
    let coords: Vec<f64> = (0..15).map(f64::from).collect();
    let text = format!(
        "t\n{}{}",
        frame_text(&coords, None),
        frame_text(&coords[..10], None)
    );
    let message = parse_mdcrd_checked(&text, 5, false).unwrap_err();
    assert!(message.contains("frame 1"), "{message}");
    assert!(message.contains("1 trailing lines"), "{message}");
    assert_eq!(parse_mdcrd(&text, 5, false).len(), 1);
}

#[test]
fn missing_or_short_box_record_is_reported() {
    let coords = [1.0, 2.0, 3.0];
    let with_box = frame_text(&coords, Some([50.0, 60.0, 70.0]));
    // The last frame lacks its box line.
    let text = format!("t\n{with_box}{}", frame_text(&coords, None));
    let message = parse_mdcrd_checked(&text, 1, true).unwrap_err();
    assert!(
        message.contains("frame 1") && message.contains("line 5"),
        "{message}"
    );
    assert!(message.contains("box record"), "{message}");
    assert_eq!(parse_mdcrd(&text, 1, true)[1].box_dims, None);
    // A box line with two numbers.
    let text = "t\n   1.000   2.000   3.000\n  50.000  60.000\n";
    let message = parse_mdcrd_checked(text, 1, true).unwrap_err();
    assert!(message.contains("line 3"), "{message}");
    assert_eq!(parse_mdcrd(text, 1, true)[0].box_dims, None);
}

#[test]
fn impossible_atom_counts_are_errors() {
    let text = "t\n   1.000   2.000   3.000\n";
    assert!(
        parse_mdcrd_checked(text, 0, false)
            .unwrap_err()
            .contains("0 atoms")
    );
    assert!(parse_mdcrd_checked(text, usize::MAX, false).is_err());
    assert!(parse_mdcrd_checked(text, usize::MAX / 3 + 1, false).is_err());
}

/// Thousandths in the range an `F8.3` field can hold, as the text a
/// writer produces for them.
fn random_field(rng: &mut SplitMix64) -> (f64, String) {
    // -999.999 ..= 9999.999 in steps of 0.001.
    let thousandths = i64::try_from(rng.below(10_999_999)).unwrap_or(0) - 999_999;
    let text = format!("{:.3}", thousandths as f64 / 1000.0);
    let value = text.parse::<f64>().unwrap_or(f64::NAN);
    (value, text)
}

#[test]
fn random_frames_with_touching_fields_round_trip_exactly() {
    let mut rng = SplitMix64::new(0x6D64_5201);
    // 7 atoms = 21 values = lines of 10, 10 and 1 value.
    let n_atoms = 7;
    for has_box in [false, true] {
        let mut text = String::from("random frames\n");
        let mut expected: Vec<(Vec<f64>, Option<[f64; 3]>)> = Vec::new();
        let mut touching_lines = 0usize;
        for _ in 0..200 {
            let mut coords = Vec::new();
            let mut fields = Vec::new();
            for _ in 0..3 * n_atoms {
                let (value, field) = random_field(&mut rng);
                coords.push(value);
                fields.push(field);
            }
            for row in fields.chunks(10) {
                let line: String = row.iter().map(|f| format!("{f:>8}")).collect();
                // Two full-width fields (or a negative one next to a wide one)
                // leave no blank between them.
                if line.split_whitespace().count() < row.len() {
                    touching_lines += 1;
                }
                text.push_str(&line);
                text.push('\n');
            }
            let mut box_dims = None;
            if has_box {
                let mut dims = [0.0; 3];
                let mut line = String::new();
                for slot in &mut dims {
                    let (value, field) = random_field(&mut rng);
                    *slot = value;
                    line.push_str(&format!("{field:>8}"));
                }
                text.push_str(&line);
                text.push('\n');
                box_dims = Some(dims);
            }
            expected.push((coords, box_dims));
        }
        assert!(touching_lines > 50, "only {touching_lines} touching lines");
        let lenient = parse_mdcrd(&text, n_atoms, has_box);
        let checked = parse_mdcrd_checked(&text, n_atoms, has_box).unwrap();
        assert_eq!(lenient.len(), expected.len());
        assert_eq!(checked.len(), expected.len());
        for (index, (coords, box_dims)) in expected.iter().enumerate() {
            assert_eq!(&lenient[index].coords, coords, "frame {index}");
            assert_eq!(&checked[index].coords, coords, "frame {index}");
            assert_eq!(lenient[index].box_dims, *box_dims, "frame {index}");
            assert_eq!(checked[index].box_dims, *box_dims, "frame {index}");
        }
    }
}

#[test]
fn mutated_trajectories_never_panic_and_stay_consistent() {
    let coords = [
        -123.456, 999.999, 1000.000, -100.000, 1.0, -999.999, 5.5, -5.5, 0.25,
    ];
    let with_box = format!(
        "title\n{}{}",
        frame_text(&coords, Some([1000.0, 999.5, 10.0])),
        frame_text(&coords, Some([50.0, 60.0, 70.0]))
    );
    let plain = format!(
        "t\n{}{}",
        frame_text(&coords, None),
        frame_text(&coords, None)
    );
    let mut rng = SplitMix64::new(0x6D64_0001);
    let mut accepted = 0usize;
    let mut rejected = 0usize;
    for (seed, has_box) in [(with_box, true), (plain, false)] {
        for round in 0..2000 {
            let text = mutate_text(&seed, &mut rng, 1 + round % 4);
            let lenient = parse_mdcrd(&text, 3, has_box);
            for frame in &lenient {
                assert_eq!(frame.coords.len(), 9, "{text}");
            }
            match parse_mdcrd_checked(&text, 3, has_box) {
                Ok(checked) => {
                    accepted += 1;
                    assert_eq!(checked.len(), lenient.len(), "{text}");
                    for (a, b) in checked.iter().zip(&lenient) {
                        assert_eq!(a.coords, b.coords, "{text}");
                        assert_eq!(a.box_dims, b.box_dims, "{text}");
                    }
                }
                Err(_) => rejected += 1,
            }
        }
    }
    assert!(accepted > 0, "no damaged trajectory was accepted");
    assert!(rejected > 0, "no damaged trajectory was rejected");
}
