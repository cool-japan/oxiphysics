//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::io::{BufRead, BufReader, Read};

use super::types::SpectrumRecord;

/// Magic bytes for the ANIF binary spectral format.
pub(super) const ANIF_MAGIC: &[u8; 4] = b"ANIF";
/// Read a simple two-column CSV file into a `SpectrumRecord`.
///
/// Lines starting with `#` are treated as comments (the first comment line
/// becomes the title). The first non-comment non-header line begins data.
pub fn read_csv_spectrum<R: Read>(reader: R) -> crate::Result<SpectrumRecord> {
    let buf = BufReader::new(reader);
    let mut record = SpectrumRecord::new();
    let mut saw_header = false;
    for line_res in buf.lines() {
        let line = line_res.map_err(crate::Error::Io)?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('#') {
            if record.metadata.title.is_empty() {
                record.metadata.title = rest.trim().to_string();
            }
            continue;
        }
        if !saw_header && trimmed.to_lowercase().starts_with('x') {
            saw_header = true;
            continue;
        }
        saw_header = true;
        let parts: Vec<&str> = trimmed.splitn(2, ',').collect();
        if parts.len() == 2
            && let (Ok(x), Ok(y)) = (
                parts[0].trim().parse::<f64>(),
                parts[1].trim().parse::<f64>(),
            )
        {
            record.x.push(x);
            record.y.push(y);
        }
    }
    Ok(record)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectroscopy_io::AnifEndian;
    use crate::spectroscopy_io::AnifReader;
    use crate::spectroscopy_io::JcampDxReader;
    use crate::spectroscopy_io::SgWindow;
    use crate::spectroscopy_io::SimilarityMetric;
    use crate::spectroscopy_io::SpectralAnalysis;
    use crate::spectroscopy_io::SpectralDatabase;
    use crate::spectroscopy_io::SpectralExport;
    use crate::spectroscopy_io::SpectrumMetadata;
    use crate::spectroscopy_io::SpectrumSource;
    #[test]
    fn test_spectrum_record_empty() {
        let r = SpectrumRecord::new();
        assert!(r.is_empty());
        assert_eq!(r.len(), 0);
    }
    #[test]
    fn test_spectrum_record_len() {
        let r = SpectrumRecord::from_arrays(
            vec![1.0, 2.0, 3.0],
            vec![0.1, 0.5, 0.9],
            SpectrumMetadata::default(),
        );
        assert_eq!(r.len(), 3);
        assert!(!r.is_empty());
    }
    #[test]
    fn test_spectrum_record_x_min_max() {
        let r = SpectrumRecord::from_arrays(
            vec![400.0, 2000.0, 4000.0],
            vec![0.0, 1.0, 0.5],
            SpectrumMetadata::default(),
        );
        assert!((r.x_min() - 400.0).abs() < 1e-10);
        assert!((r.x_max() - 4000.0).abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_record_y_min_max() {
        let r = SpectrumRecord::from_arrays(
            vec![1.0, 2.0, 3.0],
            vec![0.1, 5.0, 2.0],
            SpectrumMetadata::default(),
        );
        assert!((r.y_max() - 5.0).abs() < 1e-10);
        assert!((r.y_min() - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_record_normalize() {
        let mut r = SpectrumRecord::from_arrays(
            vec![1.0, 2.0, 3.0],
            vec![2.0, 4.0, 8.0],
            SpectrumMetadata::default(),
        );
        r.normalize();
        assert!((r.y[2] - 1.0).abs() < 1e-10);
        assert!((r.y[1] - 0.5).abs() < 1e-10);
        assert!((r.y[0] - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_record_interpolate_at_exact() {
        let r = SpectrumRecord::from_arrays(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 0.0],
            SpectrumMetadata::default(),
        );
        assert!((r.interpolate_at(1.0).unwrap() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_record_interpolate_at_midpoint() {
        let r = SpectrumRecord::from_arrays(
            vec![0.0, 2.0],
            vec![0.0, 4.0],
            SpectrumMetadata::default(),
        );
        let v = r.interpolate_at(1.0).unwrap();
        assert!((v - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_record_interpolate_out_of_range() {
        let r = SpectrumRecord::from_arrays(
            vec![1.0, 2.0],
            vec![1.0, 2.0],
            SpectrumMetadata::default(),
        );
        assert!(r.interpolate_at(0.0).is_none());
        assert!(r.interpolate_at(3.0).is_none());
    }
    #[test]
    fn test_spectrum_source_display() {
        assert_eq!(format!("{}", SpectrumSource::Ftir), "FTIR");
        assert_eq!(format!("{}", SpectrumSource::Raman), "Raman");
    }
    fn make_jdx_content(title: &str, x: &[f64], y: &[f64]) -> String {
        let mut s = String::new();
        s.push_str(&format!("##TITLE={}\n", title));
        s.push_str("##JCAMP-DX=4.24\n");
        s.push_str("##DATA TYPE=INFRARED SPECTRUM\n");
        s.push_str("##XUNITS=1/CM\n");
        s.push_str("##YUNITS=ABSORBANCE\n");
        s.push_str(&format!("##NPOINTS={}\n", x.len()));
        s.push_str(&format!("##FIRSTX={}\n", x.first().copied().unwrap_or(0.0)));
        s.push_str(&format!(
            "##DELTAX={}\n",
            if x.len() >= 2 { x[1] - x[0] } else { 1.0 }
        ));
        s.push_str("##XFACTOR=1.0\n");
        s.push_str("##YFACTOR=1.0\n");
        s.push_str("##XYDATA=(X++(Y..Y))\n");
        for (xi, yi) in x.iter().zip(y.iter()) {
            s.push_str(&format!("{} {}\n", xi, yi));
        }
        s.push_str("##END=\n");
        s
    }
    #[test]
    fn test_jdx_parse_simple() {
        let content = make_jdx_content("Test Compound", &[400.0, 500.0, 600.0], &[0.1, 0.5, 0.3]);
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 1);
        let rec = &reader.records[0];
        assert_eq!(rec.metadata.title, "Test Compound");
        assert_eq!(rec.metadata.source, SpectrumSource::Ftir);
    }
    #[test]
    fn test_jdx_parse_metadata_fields() {
        let _unused = make_jdx_content("Test", &[1.0], &[1.0]);
        let content = "##TITLE=MySpectrum\n##JCAMP-DX=4.24\n##DATA TYPE=RAMAN\n\
             ##DATE=2026-01-01\n##XUNITS=CM-1\n##YUNITS=COUNTS\n\
             ##NPOINTS=2\n##FIRSTX=100.0\n##DELTAX=50.0\n\
             ##XFACTOR=1.0\n##YFACTOR=1.0\n\
             ##XYDATA=(X++(Y..Y))\n100.0 10.0\n150.0 20.0\n##END=\n"
            .to_string();
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 1);
        let rec = &reader.records[0];
        assert_eq!(rec.metadata.title, "MySpectrum");
        assert_eq!(rec.metadata.date, "2026-01-01");
        assert_eq!(rec.metadata.source, SpectrumSource::Raman);
    }
    #[test]
    fn test_jdx_parse_xfactor_yfactor() {
        let content = "##TITLE=Scaled\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n\
                       ##XUNITS=1/CM\n##YUNITS=ABSORBANCE\n\
                       ##NPOINTS=2\n##FIRSTX=1000.0\n##DELTAX=100.0\n\
                       ##XFACTOR=2.0\n##YFACTOR=0.5\n\
                       ##XYDATA=(X++(Y..Y))\n1000.0 10.0\n1100.0 20.0\n##END=\n";
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 1);
        let rec = &reader.records[0];
        // FIRSTX and DELTAX are actual values: XFACTOR must not scale them, and
        // YFACTOR must scale each listed ordinate exactly once.
        assert_eq!(rec.x.len(), 2);
        assert_eq!(rec.y.len(), 2);
        assert!((rec.x[0] - 1000.0).abs() < 1e-9, "x[0]={}", rec.x[0]);
        assert!((rec.x[1] - 1100.0).abs() < 1e-9, "x[1]={}", rec.x[1]);
        assert!((rec.y[0] - 5.0).abs() < 1e-12, "y[0]={}", rec.y[0]);
        assert!((rec.y[1] - 10.0).abs() < 1e-12, "y[1]={}", rec.y[1]);
    }
    #[test]
    fn test_jdx_parse_multi_block() {
        let content = "\
##TITLE=Block1\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n\
##NPOINTS=2\n##FIRSTX=100.0\n##DELTAX=50.0\n\
##XFACTOR=1.0\n##YFACTOR=1.0\n##XYDATA=(X++(Y..Y))\n100.0 1.0\n150.0 2.0\n##END=\n\
##TITLE=Block2\n##JCAMP-DX=4.24\n##DATA TYPE=RAMAN\n\
##NPOINTS=2\n##FIRSTX=200.0\n##DELTAX=50.0\n\
##XFACTOR=1.0\n##YFACTOR=1.0\n##XYDATA=(X++(Y..Y))\n200.0 3.0\n250.0 4.0\n##END=\n";
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 2);
        assert_eq!(reader.records[0].metadata.title, "Block1");
        assert_eq!(reader.records[1].metadata.title, "Block2");
    }
    #[test]
    fn test_jdx_parse_comment_lines_skipped() {
        let content = "$$This is a comment\n\
                       ##TITLE=Test\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n\
                       ##NPOINTS=1\n##FIRSTX=500.0\n##DELTAX=1.0\n\
                       ##XFACTOR=1.0\n##YFACTOR=1.0\n##XYDATA=(X++(Y..Y))\n\
                       $$Another comment\n500.0 0.3\n##END=\n";
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 1);
    }
    #[test]
    fn test_jdx_empty_input() {
        let mut reader = JcampDxReader::new();
        reader.parse("".as_bytes()).unwrap();
        assert_eq!(reader.count(), 0);
    }

    fn parse_jdx(content: &str) -> crate::Result<JcampDxReader> {
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes())?;
        Ok(reader)
    }

    fn assert_close_slices(got: &[f64], want: &[f64], tol: f64) {
        assert_eq!(got.len(), want.len());
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            assert!((g - w).abs() <= tol, "index {i}: {g} vs {w}");
        }
    }

    fn roundtrip_record(original: &SpectrumRecord) -> (String, SpectrumRecord) {
        let mut buf = Vec::new();
        SpectralExport::write_jcamp_dx(original, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let reader = parse_jdx(&text).unwrap();
        assert_eq!(reader.count(), 1, "file was:\n{text}");
        let record = reader.records.into_iter().next().unwrap();
        (text, record)
    }

    fn sample_metadata(title: &str) -> SpectrumMetadata {
        SpectrumMetadata {
            title: title.to_string(),
            source: SpectrumSource::Ftir,
            x_label: "1/CM".to_string(),
            y_label: "ABSORBANCE".to_string(),
            ..Default::default()
        }
    }

    #[test]
    fn test_jdx_write_read_roundtrip_uniform_axis() {
        // 21 points: two full lines of eight ordinates plus a line of five.
        let x: Vec<f64> = (0..21).map(|i| 400.0 + 15.5 * i as f64).collect();
        let y: Vec<f64> = (0..21)
            .map(|i| (0.37 * i as f64).sin() * 1.0e-3 + 0.125 * i as f64)
            .collect();
        let original =
            SpectrumRecord::from_arrays(x.clone(), y.clone(), sample_metadata("Uniform"));
        let (text, back) = roundtrip_record(&original);
        assert!(text.contains("##XYDATA=(X++(Y..Y))"), "{text}");
        assert!(text.contains("##DELTAX="), "{text}");
        assert!(text.contains("##LASTX="), "{text}");
        assert!(!text.contains("##XYPOINTS"), "{text}");
        assert_eq!(back.x.len(), back.y.len());
        assert_close_slices(&back.x, &x, 1e-9);
        assert_close_slices(&back.y, &y, 1e-9);
        assert_eq!(back.metadata.title, "Uniform");
        assert_eq!(back.metadata.source, SpectrumSource::Ftir);
        assert_eq!(back.metadata.x_label, "1/CM");
        assert_eq!(back.metadata.y_label, "ABSORBANCE");
    }

    #[test]
    fn test_jdx_write_read_roundtrip_descending_axis() {
        // Infrared spectra are usually stored from high to low wavenumber.
        let x: Vec<f64> = (0..10).map(|i| 4000.0 - 400.0 * i as f64).collect();
        let y: Vec<f64> = (0..10).map(|i| 0.1 * (i * i) as f64).collect();
        let original = SpectrumRecord::from_arrays(x.clone(), y.clone(), sample_metadata("Desc"));
        let (text, back) = roundtrip_record(&original);
        assert!(text.contains("##FIRSTX=4e3"), "{text}");
        assert!(text.contains("##LASTX=4e2"), "{text}");
        assert_close_slices(&back.x, &x, 1e-9);
        assert_close_slices(&back.y, &y, 1e-9);
    }

    #[test]
    fn test_jdx_write_read_roundtrip_non_uniform_axis_uses_xypoints() {
        // An X++(Y..Y) table cannot describe this grid; it must not be forced
        // into one (DELTAX would silently move every point).
        let x = vec![400.0, 450.0, 500.0, 640.0, 1000.0, 1000.5, 2500.25];
        let y = vec![0.1, -0.2, 0.3, 1.0e-7, 5.5e3, 0.0, -1.25];
        let original = SpectrumRecord::from_arrays(x.clone(), y.clone(), sample_metadata("Pairs"));
        let (text, back) = roundtrip_record(&original);
        assert!(text.contains("##XYPOINTS=(XY..XY)"), "{text}");
        assert!(!text.contains("##DELTAX"), "{text}");
        assert_eq!(back.x, x);
        assert_eq!(back.y, y);
    }

    #[test]
    fn test_jdx_write_read_roundtrip_single_and_empty_records() {
        let one = SpectrumRecord::from_arrays(vec![633.0], vec![0.75], sample_metadata("One"));
        let (_, back) = roundtrip_record(&one);
        assert_eq!(back.x, vec![633.0]);
        assert_eq!(back.y, vec![0.75]);

        let empty = SpectrumRecord::from_arrays(vec![], vec![], sample_metadata("None"));
        let mut buf = Vec::new();
        SpectralExport::write_jcamp_dx(&empty, &mut buf).unwrap();
        let reader = parse_jdx(&String::from_utf8(buf).unwrap()).unwrap();
        assert_eq!(reader.count(), 0);
    }

    #[test]
    fn test_jdx_write_rejects_mismatched_lengths() {
        let bad =
            SpectrumRecord::from_arrays(vec![1.0, 2.0, 3.0], vec![1.0], sample_metadata("Bad"));
        let mut buf = Vec::new();
        assert!(SpectralExport::write_jcamp_dx(&bad, &mut buf).is_err());
    }

    #[test]
    fn test_jdx_missing_deltax_is_derived_from_lastx() {
        // Written by other tools: FIRSTX/LASTX/NPOINTS without DELTAX, three
        // ordinates per line, only the first abscissa of each line listed.
        let content = "##TITLE=NoDelta\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n\
                       ##NPOINTS=9\n##FIRSTX=100.0\n##LASTX=180.0\n\
                       ##XFACTOR=1.0\n##YFACTOR=2.0\n##XYDATA=(X++(Y..Y))\n\
                       100.0 1 2 3\n130.0 4 5 6\n160.0 7 8 9\n##END=\n";
        let reader = parse_jdx(content).unwrap();
        let rec = &reader.records[0];
        let x: Vec<f64> = (0..9).map(|i| 100.0 + 10.0 * i as f64).collect();
        let y: Vec<f64> = (1..=9).map(|v| 2.0 * v as f64).collect();
        assert_close_slices(&rec.x, &x, 1e-12);
        assert_close_slices(&rec.y, &y, 1e-12);
    }

    #[test]
    fn test_jdx_axis_falls_back_to_the_abscissa_check_values() {
        let content = "##TITLE=Checks\n##XFACTOR=0.5\n##YFACTOR=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n1000.0 1.0 2.0\n1008.0 3.0 4.0\n##END=\n";
        let reader = parse_jdx(content).unwrap();
        let rec = &reader.records[0];
        assert_close_slices(&rec.x, &[500.0, 502.0, 504.0, 506.0], 1e-12);
        assert_close_slices(&rec.y, &[1.0, 2.0, 3.0, 4.0], 1e-12);
    }

    #[test]
    fn test_jdx_axis_that_cannot_be_determined_is_an_error() {
        let content = "##TITLE=NoAxis\n##XYDATA=(X++(Y..Y))\n100.0 1.0 2.0 3.0\n##END=\n";
        let err = parse_jdx(content).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "{err}");
    }

    #[test]
    fn test_jdx_npoints_mismatch_is_an_error() {
        let content = "##TITLE=Short\n##NPOINTS=5\n##FIRSTX=1.0\n##DELTAX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n1.0 1.0 2.0 3.0\n##END=\n";
        let err = parse_jdx(content).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "{err}");
    }

    #[test]
    fn test_jdx_non_numeric_data_token_is_an_error() {
        // SQZ-compressed ordinates are not supported; they must not be skipped.
        let content = "##TITLE=Squeezed\n##NPOINTS=3\n##FIRSTX=1.0\n##DELTAX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n1.0 1 2 3A\n##END=\n";
        let err = parse_jdx(content).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "{err}");
    }

    #[test]
    fn test_jdx_blocks_without_end_markers_are_not_lost() {
        let content = "##TITLE=First\n##NPOINTS=2\n##FIRSTX=1.0\n##DELTAX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n1.0 10.0 20.0\n\
                       ##TITLE=Second\n##NPOINTS=2\n##FIRSTX=5.0\n##DELTAX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n5.0 30.0 40.0\n";
        let reader = parse_jdx(content).unwrap();
        assert_eq!(reader.count(), 2);
        assert_eq!(reader.records[0].metadata.title, "First");
        assert_eq!(reader.records[0].y, vec![10.0, 20.0]);
        assert_eq!(reader.records[1].metadata.title, "Second");
        assert_eq!(reader.records[1].x, vec![5.0, 6.0]);
        assert_eq!(reader.records[1].y, vec![30.0, 40.0]);
    }

    #[test]
    fn test_jdx_block_state_does_not_leak_between_blocks() {
        // The second block has no FIRSTX/LASTX/DELTAX of its own: it must use
        // its check values, not the first block's axis or scale factors.
        let content = "##TITLE=A\n##XFACTOR=10.0\n##YFACTOR=3.0\n##FIRSTX=0.0\n##LASTX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n0.0 1.0 1.0\n##END=\n\
                       ##TITLE=B\n##XYDATA=(X++(Y..Y))\n7.0 1.0 1.0\n8.0 2.0 2.0\n##END=\n";
        let reader = parse_jdx(content).unwrap();
        assert_eq!(reader.count(), 2);
        assert_close_slices(&reader.records[0].y, &[3.0, 3.0], 1e-12);
        assert_close_slices(&reader.records[1].x, &[7.0, 7.5, 8.0, 8.5], 1e-12);
        assert_close_slices(&reader.records[1].y, &[1.0, 1.0, 2.0, 2.0], 1e-12);
    }

    #[test]
    fn test_jdx_xypoints_pairs_are_scaled_once() {
        let content = "##TITLE=Pairs\n##XFACTOR=2.0\n##YFACTOR=0.5\n##NPOINTS=3\n\
                       ##XYPOINTS=(XY..XY)\n1,2; 3,4\n5 6\n##END=\n";
        let reader = parse_jdx(content).unwrap();
        let rec = &reader.records[0];
        assert_close_slices(&rec.x, &[2.0, 6.0, 10.0], 1e-12);
        assert_close_slices(&rec.y, &[1.0, 2.0, 3.0], 1e-12);
    }

    #[test]
    fn test_jdx_xypoints_with_an_odd_number_count_is_an_error() {
        let content = "##TITLE=Odd\n##XYPOINTS=(XY..XY)\n1,2; 3\n##END=\n";
        assert!(parse_jdx(content).is_err());
    }

    #[test]
    fn test_jdx_inline_comments_are_ignored_in_data_lines() {
        let content = "##TITLE=Inline\n##NPOINTS=2\n##FIRSTX=1.0\n##DELTAX=1.0\n\
                       ##XYDATA=(X++(Y..Y))\n1.0 0.5 0.25 $$ checksum 17\n##END=\n";
        let reader = parse_jdx(content).unwrap();
        assert_eq!(reader.records[0].y, vec![0.5, 0.25]);
    }

    #[test]
    fn test_jdx_every_record_has_matching_axis_lengths() {
        let inputs = [
            "##TITLE=A\n##FIRSTX=0\n##DELTAX=1\n##XYDATA=(X++(Y..Y))\n0 1 2 3 4\n##END=\n",
            "##TITLE=B\n##XYPOINTS=(XY..XY)\n1,2; 3,4\n##END=\n",
            "##TITLE=C\n##NTUPLES=NMR\n##PAGE=N=1\n1 2\n3 4\n##END NTUPLES=NMR\n##END=\n",
            "##TITLE=D\n##FIRSTX=1\n##LASTX=3\n##XYDATA=(X++(Y..Y))\n1 7 8 9\n",
        ];
        for content in inputs {
            let reader = parse_jdx(content).unwrap();
            assert_eq!(reader.count(), 1, "{content}");
            let rec = &reader.records[0];
            assert!(!rec.x.is_empty(), "{content}");
            assert_eq!(rec.x.len(), rec.y.len(), "{content}");
        }
    }
    #[test]
    fn test_anif_roundtrip() {
        let meta = SpectrumMetadata {
            title: "TestSpectrum".to_string(),
            instrument: "FT-1000".to_string(),
            date: "2026-03-24".to_string(),
            ..Default::default()
        };
        let original = SpectrumRecord::from_arrays(
            vec![400.0, 500.0, 600.0, 700.0],
            vec![0.1, 0.5, 0.8, 0.3],
            meta,
        );
        let bytes = AnifReader::encode(&original);
        let mut reader = AnifReader::new();
        reader.parse_bytes(&bytes).unwrap();
        let decoded = &reader.record;
        assert_eq!(decoded.len(), original.len());
        for (a, b) in decoded.y.iter().zip(original.y.iter()) {
            assert!((a - b).abs() < 1e-10, "Y mismatch: {} vs {}", a, b);
        }
    }
    #[test]
    fn test_anif_bad_magic() {
        let bytes = b"NOTANIF\x00".to_vec();
        let mut reader = AnifReader::new();
        assert!(reader.parse_bytes(&bytes).is_err());
    }
    #[test]
    fn test_anif_header_too_short() {
        let bytes = b"ANIF\x01\x00\x00\x00".to_vec();
        let mut reader = AnifReader::new();
        assert!(reader.parse_bytes(&bytes).is_err());
    }
    #[test]
    fn test_anif_encode_empty_spectrum() {
        let record = SpectrumRecord::new();
        let bytes = AnifReader::encode(&record);
        let mut reader = AnifReader::new();
        reader.parse_bytes(&bytes).unwrap();
        assert_eq!(reader.record.len(), 0);
    }
    #[test]
    fn test_anif_header_fields() {
        let meta = SpectrumMetadata {
            title: "IR Spectrum".to_string(),
            instrument: "Bruker Tensor 27".to_string(),
            date: "2026-01-01".to_string(),
            ..Default::default()
        };
        let record = SpectrumRecord::from_arrays(vec![1.0], vec![1.0], meta);
        let bytes = AnifReader::encode(&record);
        let mut reader = AnifReader::new();
        reader.parse_bytes(&bytes).unwrap();
        assert_eq!(reader.header.version, (1, 0));
        assert_eq!(reader.header.num_points, 1);
        assert_eq!(reader.header.endian, AnifEndian::Little);
    }
    fn make_gaussian_spectrum(center: f64, sigma: f64, n: usize) -> SpectrumRecord {
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let y: Vec<f64> = x
            .iter()
            .map(|&xi| (-((xi - center) / sigma).powi(2) / 2.0).exp())
            .collect();
        SpectrumRecord::from_arrays(x, y, SpectrumMetadata::default())
    }
    #[test]
    fn test_find_peaks_single_gaussian() {
        let rec = make_gaussian_spectrum(50.0, 5.0, 100);
        let peaks = SpectralAnalysis::find_peaks(&rec, 3, 0.1);
        assert!(!peaks.is_empty(), "should find at least one peak");
        let best = peaks
            .iter()
            .max_by(|a, b| a.y.partial_cmp(&b.y).unwrap())
            .unwrap();
        assert!((best.x - 50.0).abs() < 2.0, "peak x ≈ 50, got {}", best.x);
    }
    #[test]
    fn test_find_peaks_no_peaks_flat() {
        let rec = SpectrumRecord::from_arrays(
            (0..50).map(|i| i as f64).collect(),
            vec![0.5; 50],
            SpectrumMetadata::default(),
        );
        let peaks = SpectralAnalysis::find_peaks(&rec, 2, 0.3);
        assert!(peaks.is_empty(), "flat spectrum should have no peaks");
    }
    #[test]
    fn test_als_baseline_flat_signal() {
        let y = vec![1.0f64; 100];
        let corrected = SpectralAnalysis::als_baseline(&y, 1e5, 0.01, 10);
        for (i, &v) in corrected.iter().enumerate() {
            assert!(v.abs() < 0.5, "corrected[{}] = {} (too large)", i, v);
        }
    }
    #[test]
    fn test_als_baseline_returns_same_length() {
        let y: Vec<f64> = (0..200).map(|i| (i as f64 * 0.1).sin()).collect();
        let corrected = SpectralAnalysis::als_baseline(&y, 1e4, 0.01, 5);
        assert_eq!(corrected.len(), y.len());
    }
    #[test]
    fn test_sg_smoothing_preserves_length() {
        let rec = make_gaussian_spectrum(50.0, 5.0, 100);
        let smoothed = SpectralAnalysis::savitzky_golay(&rec.y, SgWindow::W7);
        assert_eq!(smoothed.len(), rec.y.len());
    }
    #[test]
    fn test_sg_smoothing_gaussian_preserved() {
        let rec = make_gaussian_spectrum(50.0, 8.0, 100);
        let smoothed = SpectralAnalysis::savitzky_golay(&rec.y, SgWindow::W5);
        let orig_sum: f64 = rec.y.iter().sum();
        let smooth_sum: f64 = smoothed.iter().sum();
        assert!((orig_sum - smooth_sum).abs() / orig_sum < 0.1);
    }
    #[test]
    fn test_moving_average_preserves_length() {
        let y: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let out = SpectralAnalysis::moving_average(&y, 3);
        assert_eq!(out.len(), y.len());
    }
    #[test]
    fn test_derivative_linear_signal() {
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 * xi).collect();
        let rec = SpectrumRecord::from_arrays(x, y, SpectrumMetadata::default());
        let deriv = SpectralAnalysis::derivative(&rec);
        for (i, &d) in deriv.iter().enumerate().skip(1).take(8) {
            assert!((d - 2.0).abs() < 1e-6, "deriv[{}] = {} expected 2.0", i, d);
        }
    }
    #[test]
    fn test_integrate_gaussian_area() {
        let rec = make_gaussian_spectrum(50.0, 5.0, 200);
        let area = SpectralAnalysis::integrate(&rec, 30.0, 70.0);
        assert!(area > 10.0, "integral should be ≈ 12.5, got {}", area);
    }
    #[test]
    fn test_fwhm_gaussian() {
        let rec = make_gaussian_spectrum(100.0, 10.0, 300);
        let peaks = SpectralAnalysis::find_peaks(&rec, 5, 0.5);
        assert!(!peaks.is_empty(), "should find a peak");
        let p = &peaks[0];
        assert!(
            p.fwhm > 15.0 && p.fwhm < 35.0,
            "FWHM={} expected ~23.5",
            p.fwhm
        );
    }
    fn make_test_record(title: &str, y_values: Vec<f64>) -> SpectrumRecord {
        let x: Vec<f64> = (0..y_values.len()).map(|i| i as f64).collect();
        let meta = SpectrumMetadata {
            title: title.to_string(),
            ..Default::default()
        };
        SpectrumRecord::from_arrays(x, y_values, meta)
    }
    #[test]
    fn test_database_insert_and_get() {
        let mut db = SpectralDatabase::new();
        let rec = make_test_record("Acetone", vec![0.1, 0.5, 1.0, 0.3]);
        let id = db.insert(rec.clone());
        assert!(db.get(id).is_some());
        assert_eq!(db.get(id).unwrap().metadata.title, "Acetone");
    }
    #[test]
    fn test_database_len_and_is_empty() {
        let mut db = SpectralDatabase::new();
        assert!(db.is_empty());
        db.insert(make_test_record("A", vec![1.0]));
        assert_eq!(db.len(), 1);
        assert!(!db.is_empty());
    }
    #[test]
    fn test_database_remove() {
        let mut db = SpectralDatabase::new();
        let id = db.insert(make_test_record("Ethanol", vec![0.2, 0.8, 0.4]));
        assert!(db.remove(id));
        assert!(db.get(id).is_none());
        assert!(!db.remove(id));
    }
    #[test]
    fn test_database_search_dot_product_self_similarity() {
        let mut db = SpectralDatabase::new();
        let rec = make_test_record("Caffeine", vec![0.1, 0.3, 0.9, 0.5, 0.2]);
        let id = db.insert(rec.clone());
        let results = db.search(&rec, SimilarityMetric::DotProduct, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].id, id);
        assert!(
            (results[0].score - 1.0).abs() < 1e-6,
            "self-similarity should be 1.0"
        );
    }
    #[test]
    fn test_database_search_euclidean_self_similarity() {
        let mut db = SpectralDatabase::new();
        let rec = make_test_record("Ethanol", vec![0.1, 0.9, 0.5]);
        let id = db.insert(rec.clone());
        let results = db.search(&rec, SimilarityMetric::Euclidean, 1);
        assert_eq!(results[0].id, id);
        assert!(
            results[0].score < 1e-10,
            "Euclidean distance to itself should be 0"
        );
    }
    #[test]
    fn test_database_search_pearson_self_similarity() {
        let mut db = SpectralDatabase::new();
        let rec = make_test_record("Water", vec![0.5, 1.0, 0.7, 0.3]);
        let id = db.insert(rec.clone());
        let results = db.search(&rec, SimilarityMetric::Pearson, 1);
        assert_eq!(results[0].id, id);
        assert!(
            (results[0].score - 1.0).abs() < 1e-6,
            "Pearson self-correlation ≈ 1"
        );
    }
    #[test]
    fn test_database_search_top_k() {
        let mut db = SpectralDatabase::new();
        for i in 0..5 {
            db.insert(make_test_record(
                &format!("Compound{}", i),
                vec![i as f64; 4],
            ));
        }
        let query = make_test_record("Query", vec![3.0; 4]);
        let results = db.search(&query, SimilarityMetric::DotProduct, 3);
        assert_eq!(results.len(), 3);
    }
    #[test]
    fn test_database_search_sam() {
        let mut db = SpectralDatabase::new();
        let a = make_test_record("A", vec![1.0, 0.0, 0.0]);
        let b = make_test_record("B", vec![0.0, 1.0, 0.0]);
        db.insert(a.clone());
        db.insert(b.clone());
        let results = db.search(&a, SimilarityMetric::SpectralAngleMapper, 2);
        assert_eq!(results.len(), 2);
    }
    #[test]
    fn test_database_iter() {
        let mut db = SpectralDatabase::new();
        db.insert(make_test_record("X", vec![1.0]));
        db.insert(make_test_record("Y", vec![2.0]));
        let titles: Vec<&str> = db
            .iter()
            .map(|e| e.record.metadata.title.as_str())
            .collect();
        assert!(titles.contains(&"X"));
        assert!(titles.contains(&"Y"));
    }
    #[test]
    fn test_export_csv_roundtrip() {
        let original = make_test_record("Test", vec![0.1, 0.5, 1.0, 0.5, 0.1]);
        let mut buf = Vec::new();
        SpectralExport::write_csv(&original, &mut buf).unwrap();
        let content = String::from_utf8(buf).unwrap();
        assert!(content.contains("x,y"));
        assert!(content.contains("0,0.1") || content.contains("0,"));
    }
    #[test]
    fn test_export_jcamp_dx_contains_title() {
        let rec = make_test_record("Aspirin", vec![0.3, 0.7, 1.0]);
        let mut buf = Vec::new();
        SpectralExport::write_jcamp_dx(&rec, &mut buf).unwrap();
        let content = String::from_utf8(buf).unwrap();
        assert!(content.contains("##TITLE=Aspirin"));
        assert!(content.contains("##END="));
    }
    #[test]
    fn test_export_plain_text_contains_header_comment() {
        let mut rec = make_test_record("Benzene", vec![0.2, 0.4]);
        rec.metadata.instrument = "FT-IR-9000".to_string();
        let mut buf = Vec::new();
        SpectralExport::write_plain_text(&rec, &mut buf).unwrap();
        let content = String::from_utf8(buf).unwrap();
        assert!(content.contains("# Spectrum: Benzene"));
        assert!(content.contains("FT-IR-9000"));
    }
    #[test]
    fn test_export_multi_csv() {
        let r1 = make_test_record("A", vec![1.0, 2.0]);
        let r2 = make_test_record("B", vec![3.0, 4.0]);
        let mut buf = Vec::new();
        SpectralExport::write_multi_csv(&[r1, r2], &mut buf).unwrap();
        let content = String::from_utf8(buf).unwrap();
        assert!(content.contains("A"));
        assert!(content.contains("B"));
    }
    #[test]
    fn test_export_jcamp_dx_parseable() {
        let original = SpectrumRecord::from_arrays(
            vec![400.0, 600.0, 800.0, 1000.0],
            vec![0.1, 0.5, 0.8, 0.3],
            SpectrumMetadata {
                title: "RoundTrip".to_string(),
                source: SpectrumSource::Ftir,
                x_label: "1/CM".to_string(),
                y_label: "ABSORBANCE".to_string(),
                ..Default::default()
            },
        );
        let mut buf = Vec::new();
        SpectralExport::write_jcamp_dx(&original, &mut buf).unwrap();
        let content = String::from_utf8(buf).unwrap();
        let mut reader = JcampDxReader::new();
        reader.parse(content.as_bytes()).unwrap();
        assert_eq!(reader.count(), 1);
        assert_eq!(reader.records[0].metadata.title, "RoundTrip");
        // The data survive the trip, not just the title.
        assert_close_slices(&reader.records[0].x, &original.x, 1e-9);
        assert_close_slices(&reader.records[0].y, &original.y, 1e-9);
    }
    #[test]
    fn test_read_csv_spectrum() {
        let csv = "# My IR Spectrum\nx,y\n400.0,0.1\n500.0,0.5\n600.0,0.9\n";
        let rec = read_csv_spectrum(csv.as_bytes()).unwrap();
        assert_eq!(rec.len(), 3);
        assert_eq!(rec.metadata.title, "My IR Spectrum");
        assert!((rec.x[0] - 400.0).abs() < 1e-10);
        assert!((rec.y[2] - 0.9).abs() < 1e-10);
    }
    #[test]
    fn test_read_csv_spectrum_empty() {
        let csv = "";
        let rec = read_csv_spectrum(csv.as_bytes()).unwrap();
        assert!(rec.is_empty());
    }
    #[test]
    fn test_sg_window_half_sizes() {
        assert_eq!(SgWindow::W5.half(), 2);
        assert_eq!(SgWindow::W11.half(), 5);
        assert_eq!(SgWindow::W25.half(), 12);
    }
    #[test]
    fn test_moving_average_constant() {
        let y = vec![3.0f64; 50];
        let out = SpectralAnalysis::moving_average(&y, 5);
        for (i, &v) in out.iter().enumerate().skip(5).take(40) {
            assert!(
                (v - 3.0).abs() < 1e-10,
                "moving_average[{}] = {} expected 3.0",
                i,
                v
            );
        }
    }
    #[test]
    fn test_derivative_constant_zero() {
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y = vec![5.0f64; 20];
        let rec = SpectrumRecord::from_arrays(x, y, SpectrumMetadata::default());
        let deriv = SpectralAnalysis::derivative(&rec);
        for (i, &d) in deriv.iter().enumerate().skip(1).take(18) {
            assert!(d.abs() < 1e-10, "deriv[{}] = {} should be 0", i, d);
        }
    }
    #[test]
    fn test_integrate_zero_range() {
        let rec = make_test_record("Z", vec![1.0, 2.0, 3.0]);
        let area = SpectralAnalysis::integrate(&rec, 5.0, 10.0);
        assert!(area.abs() < 1e-10);
    }
    #[test]
    fn test_spectrum_metadata_default() {
        let meta = SpectrumMetadata::default();
        assert!(meta.title.is_empty());
        assert!(meta.date.is_empty());
    }
}
