//! Regression tests for the streaming codecs and the chunk container.
//!
//! These cover byte-exact round trips of arbitrary binary data (literal `0xFF`
//! bytes used to be mistaken for run markers) and decoding of corrupt chunk
//! files (wrapping size arithmetic used to hang or abort the process).

use super::*;

/// Small deterministic xorshift generator so no dev-dependency is needed.
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

    /// Uniform integer in `0..bound`.
    fn below(&mut self, bound: u64) -> u64 {
        self.next_u64() % bound
    }
}

fn temp_path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "oxiphysics_stream_{name}_{}.bin",
        std::process::id()
    ))
}

fn make_chunk(frame_id: usize, n_particles: usize) -> SimulationChunk {
    let mut chunk = SimulationChunk::new(frame_id, frame_id as f64 * 0.25);
    for i in 0..n_particles {
        chunk.add_particle([i as f32, 1.5, -2.0], [0.5, i as f32 * 2.0, 0.0]);
    }
    chunk
}

/// Header-only chunk bytes: frame id, particle count, time (no payload).
fn header_bytes(frame_id: u64, particles: u64) -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&frame_id.to_le_bytes());
    bytes.extend_from_slice(&particles.to_le_bytes());
    bytes.extend_from_slice(&1.0f64.to_le_bytes());
    bytes
}

// ── RleEncoder ──────────────────────────────────────────────────────────────

#[test]
fn rle_escapes_a_literal_marker_byte() {
    // Regression: encode([0xFF, 0x05, 0x07]) returned the input unchanged and
    // decode then read it as "five times 0x07".
    let data = [0xFF, 0x05, 0x07];
    let encoded = RleEncoder::encode(&data);
    assert_eq!(encoded, vec![0xFF, 1, 0xFF, 0x05, 0x07]);
    assert_eq!(RleEncoder::decode(&encoded).unwrap(), data);
}

#[test]
fn rle_runs_of_marker_bytes_roundtrip() {
    for len in [1usize, 2, 3, 253, 254, 255, 508, 1000] {
        let data = vec![0xFFu8; len];
        let encoded = RleEncoder::encode(&data);
        assert_eq!(RleEncoder::decode(&encoded).unwrap(), data, "len = {len}");
    }
}

#[test]
fn rle_every_short_input_over_marker_alphabet_roundtrips() {
    let alphabet = [0x00u8, 0x05, 0xFF];
    for len in 1..=6usize {
        let combinations = alphabet.len().pow(len as u32);
        for code in 0..combinations {
            let mut rest = code;
            let mut data = Vec::with_capacity(len);
            for _ in 0..len {
                data.push(alphabet[rest % alphabet.len()]);
                rest /= alphabet.len();
            }
            let encoded = RleEncoder::encode(&data);
            let decoded = RleEncoder::decode(&encoded).unwrap();
            assert_eq!(decoded, data, "input {data:?} encoded as {encoded:?}");
        }
    }
}

#[test]
fn rle_random_vectors_with_many_marker_bytes_roundtrip() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    for _ in 0..10_000 {
        let len = rng.below(96) as usize;
        let mut data: Vec<u8> = Vec::with_capacity(len);
        for _ in 0..len {
            let roll = rng.below(100);
            let byte = if roll < 30 {
                0xFF
            } else if roll < 50 {
                data.last().copied().unwrap_or(0)
            } else {
                (rng.next_u64() & 0xFF) as u8
            };
            data.push(byte);
        }
        let encoded = RleEncoder::encode(&data);
        assert_eq!(RleEncoder::decode(&encoded).unwrap(), data);
    }
}

#[test]
fn rle_long_runs_are_split_into_tokens_of_at_most_254() {
    let mut data = vec![0x07u8; 1000];
    data.extend(vec![0xFFu8; 600]);
    let encoded = RleEncoder::encode(&data);
    for token in encoded.chunks(3) {
        assert_eq!(token[0], 0xFF);
        assert!((1..=254).contains(&token[1]), "token {token:?}");
    }
    assert_eq!(RleEncoder::decode(&encoded).unwrap(), data);
}

#[test]
fn rle_decode_rejects_truncated_and_zero_length_tokens() {
    assert!(RleEncoder::decode(&[0x01, 0xFF]).is_err());
    assert!(RleEncoder::decode(&[0xFF, 0x03]).is_err());
    assert!(RleEncoder::decode(&[0x01, 0x02, 0xFF, 0x03]).is_err());
    let zero_run = RleEncoder::decode(&[0xFF, 0x00, 0x07]);
    assert!(zero_run.is_err());
}

// ── CompressionStream ───────────────────────────────────────────────────────

fn assert_bits_roundtrip(data: &[f32]) {
    let compressed = CompressionStream::compress_f32(data);
    let decoded = CompressionStream::decompress_f32(&compressed).unwrap();
    let want: Vec<u32> = data.iter().map(|v| v.to_bits()).collect();
    let got: Vec<u32> = decoded.iter().map(|v| v.to_bits()).collect();
    assert_eq!(got, want);
}

#[test]
fn compress_f32_literals_whose_low_byte_is_ff_roundtrip() {
    // Regression: 0x3F80_00FF is a literal whose first byte equals the run
    // marker; the decoder read [0xFF, 0x00, ...] as a zero-count run and
    // swallowed six bytes, desynchronising every later value.
    assert_bits_roundtrip(&[
        f32::from_bits(0x3F80_00FF),
        1.0,
        1.0,
        1.0,
        f32::from_bits(0xFF00_00FF),
    ]);
}

#[test]
fn compress_f32_random_bit_patterns_roundtrip_bit_exactly() {
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    for round in 0..200 {
        let len = rng.below(64) as usize;
        let mut data: Vec<f32> = Vec::with_capacity(len);
        for _ in 0..len {
            let roll = rng.below(100);
            let value = if roll < 30 {
                // Low byte (the first byte on the wire) is the marker.
                f32::from_bits(((rng.next_u64() as u32) & 0xFFFF_FF00) | 0xFF)
            } else if roll < 55 {
                data.last().copied().unwrap_or(0.0)
            } else {
                f32::from_bits(rng.next_u64() as u32)
            };
            data.push(value);
        }
        let compressed = CompressionStream::compress_f32(&data);
        let decoded = CompressionStream::decompress_f32(&compressed).unwrap();
        let want: Vec<u32> = data.iter().map(|v| v.to_bits()).collect();
        let got: Vec<u32> = decoded.iter().map(|v| v.to_bits()).collect();
        assert_eq!(got, want, "round {round}");
    }
}

#[test]
fn compress_f32_keeps_signed_zeros_and_nan_payloads_apart() {
    // `==` would merge 0.0 with -0.0; the codec compares bit patterns.
    assert_bits_roundtrip(&[0.0, -0.0, 0.0, -0.0, -0.0, -0.0]);
    let nan_a = f32::from_bits(0x7FC0_1234);
    let nan_b = f32::from_bits(0x7FC0_1235);
    assert_bits_roundtrip(&[nan_a, nan_a, nan_a, nan_b, nan_b, nan_b]);
}

#[test]
fn compress_f32_long_runs_are_split() {
    let data = vec![3.5f32; 1000];
    let compressed = CompressionStream::compress_f32(&data);
    assert!(
        compressed.len() < 40,
        "compressed to {} bytes",
        compressed.len()
    );
    assert_bits_roundtrip(&data);
}

#[test]
fn decompress_f32_rejects_truncated_and_corrupt_streams() {
    // Run token with only four of its six bytes.
    assert!(CompressionStream::decompress_f32(&[0xFF, 0x03, 0x00, 0x00]).is_err());
    // Zero-count run.
    assert!(CompressionStream::decompress_f32(&[0xFF, 0x00, 0x00, 0x00, 0x80, 0x3F]).is_err());
    // Trailing partial literal.
    assert!(CompressionStream::decompress_f32(&[0x01, 0x02, 0x03]).is_err());
    let one = 1.0f32.to_le_bytes();
    let mut good_then_partial = one.to_vec();
    good_then_partial.extend_from_slice(&[0x00, 0x00]);
    assert!(CompressionStream::decompress_f32(&good_then_partial).is_err());
    assert!(CompressionStream::decompress_f32(&[]).unwrap().is_empty());
}

// ── chunk_from_bytes ────────────────────────────────────────────────────────

#[test]
fn chunk_from_bytes_returns_none_for_overflowing_particle_counts() {
    // Regression: `24 + n * 2 * 3 * 4` wrapped (or panicked in debug builds)
    // and `Vec::with_capacity(n)` then aborted.
    for particles in [u64::MAX, u64::MAX / 24 + 1, 1u64 << 62, 1u64 << 40] {
        assert!(
            chunk_from_bytes(&header_bytes(0, particles)).is_none(),
            "particles = {particles}"
        );
    }
}

#[test]
fn chunk_from_bytes_returns_none_for_short_payload() {
    let mut bytes = chunk_to_bytes(&make_chunk(3, 4));
    bytes.truncate(bytes.len() - 1);
    assert!(chunk_from_bytes(&bytes).is_none());
}

#[test]
fn chunk_from_bytes_roundtrips_exactly() {
    let chunk = make_chunk(9, 5);
    let back = chunk_from_bytes(&chunk_to_bytes(&chunk)).unwrap();
    assert_eq!(back.frame_id, 9);
    assert_eq!(back.particle_count, 5);
    assert_eq!(back.positions, chunk.positions);
    assert_eq!(back.velocities, chunk.velocities);
    assert_eq!(back.time, chunk.time);
}

// ── ChunkedReader ───────────────────────────────────────────────────────────

fn write_frames(path: &std::path::Path, frames: usize) {
    let _ = std::fs::remove_file(path);
    let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 16);
    for i in 0..frames {
        writer.write_chunk(&make_chunk(i, 3)).unwrap();
    }
}

#[test]
fn chunked_reader_wrapping_len_returns_err() {
    // Regression: frame_len = u64::MAX - 7 made `pos += 8 + frame_len` wrap to
    // a zero step, so `open` pushed offsets forever until the process ran out
    // of memory.
    let path = temp_path("wrapping_len");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(u64::MAX - 7).to_le_bytes());
    bytes.extend_from_slice(&[0u8; 32]);
    std::fs::write(&path, &bytes).unwrap();
    let result = ChunkedReader::open(path.to_str().unwrap_or(""));
    let _ = std::fs::remove_file(&path);
    assert!(result.is_err());
}

#[test]
fn chunked_reader_huge_len_returns_err() {
    let path = temp_path("huge_len");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&u64::MAX.to_le_bytes());
    std::fs::write(&path, &bytes).unwrap();
    let result = ChunkedReader::open(path.to_str().unwrap_or(""));
    let _ = std::fs::remove_file(&path);
    assert!(result.is_err());
}

#[test]
fn chunked_reader_truncated_last_frame_returns_err() {
    let path = temp_path("truncated_last");
    write_frames(&path, 3);
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&path, &bytes[..bytes.len() - 5]).unwrap();
    let result = ChunkedReader::open(path.to_str().unwrap_or(""));
    let _ = std::fs::remove_file(&path);
    let message = result.unwrap_err();
    assert!(message.contains("truncated"), "message = {message}");
}

#[test]
fn chunked_reader_partial_trailing_header_returns_err() {
    let path = temp_path("partial_header");
    write_frames(&path, 2);
    let mut bytes = std::fs::read(&path).unwrap();
    bytes.extend_from_slice(&[1, 2, 3]);
    std::fs::write(&path, &bytes).unwrap();
    let result = ChunkedReader::open(path.to_str().unwrap_or(""));
    let _ = std::fs::remove_file(&path);
    assert!(result.is_err());
}

#[test]
fn chunked_reader_empty_file_has_no_frames() {
    let path = temp_path("empty");
    std::fs::write(&path, []).unwrap();
    let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
    let _ = std::fs::remove_file(&path);
    assert_eq!(reader.total_frames(), 0);
}

#[test]
fn chunked_reader_indexes_every_complete_frame() {
    let path = temp_path("indexes");
    write_frames(&path, 4);
    let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(reader.total_frames(), 4);
    for i in 0..4 {
        let chunk = reader.read_chunk(i).unwrap();
        assert_eq!(chunk.frame_id, i);
        assert_eq!(chunk.particle_count, 3);
    }
    let _ = std::fs::remove_file(&path);
}

#[test]
fn chunked_reader_read_chunk_rejects_a_length_that_exceeds_the_file() {
    // The file is valid when indexed and corrupted afterwards, so `read_chunk`
    // has to defend itself: the stored length would otherwise be handed to
    // `vec![0u8; len]`.
    let path = temp_path("grown_len");
    write_frames(&path, 2);
    let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
    let mut bytes = std::fs::read(&path).unwrap();
    let second = usize::try_from(reader.offsets[1]).unwrap();
    bytes[second..second + 8].copy_from_slice(&(u64::MAX / 2).to_le_bytes());
    std::fs::write(&path, &bytes).unwrap();
    let result = reader.read_chunk(1);
    let _ = std::fs::remove_file(&path);
    let message = result.unwrap_err();
    assert!(message.contains("remain"), "message = {message}");
}

#[test]
fn chunked_reader_read_chunk_rejects_a_file_shrunk_after_open() {
    let path = temp_path("shrunk");
    write_frames(&path, 3);
    let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&path, &bytes[..bytes.len() / 2]).unwrap();
    let result = reader.read_chunk(2);
    let _ = std::fs::remove_file(&path);
    assert!(result.is_err());
}
