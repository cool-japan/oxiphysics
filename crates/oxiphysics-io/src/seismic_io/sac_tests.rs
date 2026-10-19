//! Tests for the binary SAC header layout and the `SacFile` reader/writer.

use super::sac_layout as at;
use super::*;

fn le_i32(buf: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes([
        buf[offset],
        buf[offset + 1],
        buf[offset + 2],
        buf[offset + 3],
    ])
}

fn le_f32(buf: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes([
        buf[offset],
        buf[offset + 1],
        buf[offset + 2],
        buf[offset + 3],
    ])
}

/// A header with a distinctive value in every field.
fn full_header(npts: i32) -> SacHeader {
    SacHeader {
        delta: 0.025,
        depmin: -3.5,
        depmax: 7.25,
        b: -1.5,
        e: 98.5,
        o: 0.125,
        a: 12.75,
        dist: 1234.5,
        az: 45.0,
        baz: 225.0,
        cmpaz: 90.0,
        cmpinc: 89.5,
        stla: 34.945_91,
        stlo: -106.457_1,
        stel: 1820.0,
        evla: -12.5,
        evlo: 100.25,
        evdp: 33.0,
        mag: 6.5,
        npts,
        nzyear: 2024,
        nzjday: 183,
        nzhour: 23,
        nzmin: 59,
        nzsec: 58,
        nzmsec: 999,
        nvhdr: 6,
        kstnm: "ANMO".to_string(),
        kcmpnm: "BHZ".to_string(),
        knetwk: "IU".to_string(),
    }
}

fn sample_file() -> SacFile {
    let data: Vec<f32> = (0..50).map(|i| (i as f32 * 0.37).sin() * 100.0).collect();
    SacFile {
        header: full_header(data.len() as i32),
        data,
    }
}

fn assert_headers_equal(got: &SacHeader, want: &SacHeader) {
    let floats = [
        (got.delta, want.delta),
        (got.depmin, want.depmin),
        (got.depmax, want.depmax),
        (got.b, want.b),
        (got.e, want.e),
        (got.o, want.o),
        (got.a, want.a),
        (got.dist, want.dist),
        (got.az, want.az),
        (got.baz, want.baz),
        (got.cmpaz, want.cmpaz),
        (got.cmpinc, want.cmpinc),
        (got.stla, want.stla),
        (got.stlo, want.stlo),
        (got.stel, want.stel),
        (got.evla, want.evla),
        (got.evlo, want.evlo),
        (got.evdp, want.evdp),
        (got.mag, want.mag),
    ];
    for (index, (g, w)) in floats.iter().enumerate() {
        assert_eq!(g.to_bits(), w.to_bits(), "float field #{index}");
    }
    assert_eq!(got.npts, want.npts);
    assert_eq!(got.nzyear, want.nzyear);
    assert_eq!(got.nzjday, want.nzjday);
    assert_eq!(got.nzhour, want.nzhour);
    assert_eq!(got.nzmin, want.nzmin);
    assert_eq!(got.nzsec, want.nzsec);
    assert_eq!(got.nzmsec, want.nzmsec);
    assert_eq!(got.nvhdr, want.nvhdr);
    assert_eq!(got.kstnm, want.kstnm);
    assert_eq!(got.kcmpnm, want.kcmpnm);
    assert_eq!(got.knetwk, want.knetwk);
}

/// Convert a little-endian SAC file to big-endian: swap every four-byte word
/// of the numeric header and of the data (character fields are bytes).
fn to_big_endian(le: &[u8]) -> Vec<u8> {
    let mut out = le.to_vec();
    for (index, word) in out.chunks_exact_mut(4).enumerate() {
        let offset = index * 4;
        if !(at::CHAR_START..SAC_HEADER_BYTES).contains(&offset) {
            word.reverse();
        }
    }
    out
}

#[test]
fn sac_header_has_nvhdr_at_byte_304_and_nzjday_at_284() {
    // Regression: NVHDR was written at byte 284, the NZJDAY word, so readers
    // rejected the file.
    let buf = sample_file().serialize_header();
    assert_eq!(buf.len(), 632);
    assert_eq!(&buf[304..308], &6i32.to_le_bytes());
    assert_eq!(le_i32(&buf, 284), 183, "NZJDAY");
    assert_eq!(le_i32(&buf, 280), 2024, "NZYEAR");
    assert_eq!(le_i32(&buf, 288), 23, "NZHOUR");
    assert_eq!(le_i32(&buf, 292), 59, "NZMIN");
    assert_eq!(le_i32(&buf, 296), 58, "NZSEC");
    assert_eq!(le_i32(&buf, 300), 999, "NZMSEC");
    assert_eq!(le_i32(&buf, 316), 50, "NPTS");
}

#[test]
fn sac_header_is_a_time_series_with_even_sampling() {
    let buf = sample_file().serialize_header();
    assert_eq!(le_i32(&buf, 340), 1, "IFTYPE must be ITIME");
    assert_eq!(le_i32(&buf, 420), 1, "LEVEN must be true");
}

#[test]
fn sac_header_writes_float_fields_at_their_word_offsets() {
    let buf = sample_file().serialize_header();
    let expected = [
        (0, 0.025f32),     // DELTA
        (4, -3.5),         // DEPMIN
        (8, 7.25),         // DEPMAX
        (20, -1.5),        // B
        (24, 98.5),        // E
        (28, 0.125),       // O
        (32, 12.75),       // A
        (124, 34.945_91),  // STLA
        (128, -106.457_1), // STLO
        (132, 1820.0),     // STEL
        (140, -12.5),      // EVLA
        (144, 100.25),     // EVLO
        (152, 33.0),       // EVDP
        (156, 6.5),        // MAG
        (200, 1234.5),     // DIST
        (204, 45.0),       // AZ
        (208, 225.0),      // BAZ
        (228, 90.0),       // CMPAZ
        (232, 89.5),       // CMPINC
    ];
    for (offset, value) in expected {
        assert_eq!(
            le_f32(&buf, offset).to_bits(),
            value.to_bits(),
            "offset {offset}"
        );
    }
}

#[test]
fn sac_header_marks_every_unset_word_as_undefined() {
    let buf = sample_file().serialize_header();
    // Float words that SacHeader has no field for: SCALE, ODELTA, T0, F, USER0.
    for offset in [12usize, 16, 40, 80, 160] {
        assert_eq!(le_f32(&buf, offset), -12345.0, "float word at {offset}");
    }
    // Integer words: NORID, NEVID, IDEP, IZTYPE, IEVTYP.
    for offset in [308usize, 312, 344, 348, 368] {
        assert_eq!(le_i32(&buf, offset), -12345, "integer word at {offset}");
    }
    // KEVNM (16 bytes), KHOLE, KO.
    assert_eq!(&buf[448..464], b"-12345  -12345  ");
    assert_eq!(&buf[464..472], b"-12345  ");
    assert_eq!(&buf[472..480], b"-12345  ");
    // The words we do set are not undefined.
    assert_ne!(le_f32(&buf, 0), -12345.0);
}

#[test]
fn sac_header_character_fields_are_blank_padded_and_truncated() {
    let mut file = sample_file();
    file.header.kstnm = "LONGSTATION".to_string();
    file.header.kcmpnm = "Z".to_string();
    file.header.knetwk = "ü-net".to_string();
    let buf = file.serialize_header();
    assert_eq!(&buf[440..448], b"LONGSTAT");
    assert_eq!(&buf[600..608], b"Z       ");
    assert_eq!(&buf[608..616], b"?-net   ");
}

#[test]
fn sac_file_roundtrips_every_header_field_and_sample() {
    let original = sample_file();
    let bytes = original.to_bytes().unwrap();
    assert_eq!(bytes.len(), 632 + 4 * original.data.len());
    let back = SacFile::from_bytes(&bytes).unwrap();
    assert_headers_equal(&back.header, &original.header);
    assert_eq!(back.data.len(), original.data.len());
    for (g, w) in back.data.iter().zip(&original.data) {
        assert_eq!(g.to_bits(), w.to_bits());
    }
}

#[test]
fn sac_default_header_roundtrips() {
    let original = SacFile::from_data(vec![1.0, -2.0, 3.5], 0.01);
    let back = SacFile::from_bytes(&original.to_bytes().unwrap()).unwrap();
    assert_eq!(back.header.npts, 3);
    assert_eq!(back.header.nvhdr, 6);
    assert_eq!(back.header.delta, 0.01);
    assert_eq!(back.header.e, original.header.e);
    assert_eq!(back.header.depmin, -2.0);
    assert_eq!(back.header.depmax, 3.5);
    assert_eq!(back.header.o, -12345.0);
    // "-12345  " is read back without its blank padding.
    assert_eq!(back.header.kstnm, "-12345");
    assert_eq!(back.data, vec![1.0, -2.0, 3.5]);
}

#[test]
fn sac_from_bytes_detects_big_endian_files() {
    let original = sample_file();
    let be = to_big_endian(&original.to_bytes().unwrap());
    assert_eq!(&be[304..308], &6i32.to_be_bytes());
    let back = SacFile::from_bytes(&be).unwrap();
    assert_headers_equal(&back.header, &original.header);
    for (g, w) in back.data.iter().zip(&original.data) {
        assert_eq!(g.to_bits(), w.to_bits());
    }
}

#[test]
fn sac_writer_always_emits_header_version_6() {
    let mut file = sample_file();
    file.header.nvhdr = 7;
    let bytes = file.to_bytes().unwrap();
    assert_eq!(le_i32(&bytes, 304), 6);
    assert_eq!(SacFile::from_bytes(&bytes).unwrap().header.nvhdr, 6);
}

#[test]
fn sac_from_bytes_reads_version_7_and_ignores_its_footer() {
    let file = sample_file();
    let mut bytes = file.to_bytes().unwrap();
    bytes[304..308].copy_from_slice(&7i32.to_le_bytes());
    bytes.extend_from_slice(&[0xAB; 176]);
    let back = SacFile::from_bytes(&bytes).unwrap();
    assert_eq!(back.header.nvhdr, 7);
    assert_eq!(back.data, file.data);
}

#[test]
fn sac_file_without_samples_roundtrips() {
    let file = SacFile {
        header: full_header(0),
        data: Vec::new(),
    };
    let bytes = file.to_bytes().unwrap();
    assert_eq!(bytes.len(), 632);
    let back = SacFile::from_bytes(&bytes).unwrap();
    assert!(back.data.is_empty());
    assert_eq!(back.header.npts, 0);
}

#[test]
fn sac_from_bytes_rejects_short_and_foreign_buffers() {
    for len in [0usize, 1, 100, 631] {
        let err = SacFile::from_bytes(&vec![0u8; len]).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "len {len}: {err}");
    }
    // 632 zero bytes: NVHDR is 0 in either byte order.
    assert!(SacFile::from_bytes(&vec![0u8; 632]).is_err());
    // A SEG-Y style buffer is not a SAC file either.
    assert!(SacFile::from_bytes(&[0x40u8; 4000]).is_err());
}

#[test]
fn sac_from_bytes_rejects_unsupported_header_versions() {
    for nvhdr in [1i32, 5, 8, 100, -6] {
        let mut bytes = sample_file().to_bytes().unwrap();
        bytes[304..308].copy_from_slice(&nvhdr.to_le_bytes());
        let err = SacFile::from_bytes(&bytes).unwrap_err();
        assert!(err.to_string().contains("NVHDR"), "nvhdr {nvhdr}: {err}");
    }
}

#[test]
fn sac_from_bytes_rejects_negative_and_oversized_sample_counts() {
    for npts in [-1i32, i32::MIN, 51, 1_000_000, i32::MAX] {
        let mut bytes = sample_file().to_bytes().unwrap();
        bytes[316..320].copy_from_slice(&npts.to_le_bytes());
        assert!(SacFile::from_bytes(&bytes).is_err(), "npts {npts}");
    }
}

#[test]
fn sac_from_bytes_rejects_truncated_data() {
    let bytes = sample_file().to_bytes().unwrap();
    for cut in [1usize, 2, 4, 100, 200] {
        assert!(
            SacFile::from_bytes(&bytes[..bytes.len() - cut]).is_err(),
            "cut {cut}"
        );
    }
}

#[test]
fn sac_from_bytes_rejects_non_time_series_files() {
    let mut spectral = sample_file().to_bytes().unwrap();
    spectral[340..344].copy_from_slice(&3i32.to_le_bytes()); // IAMPH
    let err = SacFile::from_bytes(&spectral).unwrap_err();
    assert!(err.to_string().contains("IFTYPE"), "{err}");

    let mut uneven = sample_file().to_bytes().unwrap();
    uneven[420..424].copy_from_slice(&0i32.to_le_bytes());
    let err = SacFile::from_bytes(&uneven).unwrap_err();
    assert!(err.to_string().contains("unevenly"), "{err}");
}

#[test]
fn sac_to_bytes_rejects_inconsistent_npts() {
    let mut file = sample_file();
    file.data.pop();
    let err = file.to_bytes().unwrap_err();
    assert!(matches!(err, crate::Error::General(_)), "{err}");
    let mut negative = sample_file();
    negative.header.npts = -4;
    assert!(negative.to_bytes().is_err());
}

#[test]
fn sac_npts_helper_agrees_with_the_header_layout() {
    let buf = sample_file().serialize_header();
    assert_eq!(SacFile::npts_from_header(&buf), 50);
}
