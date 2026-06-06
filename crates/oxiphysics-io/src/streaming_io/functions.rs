//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    RingBuffer, SimulationChunk, StreamingXYZReader, StreamingXYZWriter, TrajectorySampler,
};

#[cfg(test)]
mod tests {

    use crate::streaming_io::BinaryProtocol;

    use crate::streaming_io::CompressionStream;

    use crate::streaming_io::FrameInterpolator;
    use crate::streaming_io::InterpolationMethod;

    use crate::streaming_io::NetworkTransport;

    use crate::streaming_io::PhysicsStreamClient;
    use crate::streaming_io::PhysicsStreamServer;
    use crate::streaming_io::ProtocolFrame;

    use crate::streaming_io::StreamReader;
    use crate::streaming_io::StreamWriter;

    use crate::streaming_io::WebSocketBridge;
    use crate::streaming_io::WsFrame;
    use crate::streaming_io::WsFrameType;
    use crate::streaming_io::types::*;
    #[test]
    fn test_ring_buffer_push_pop() {
        let mut rb: RingBuffer<i32> = RingBuffer::new(3);
        assert!(rb.push(1));
        assert!(rb.push(2));
        assert!(rb.push(3));
        assert!(!rb.push(4));
        assert_eq!(rb.pop(), Some(1));
        assert_eq!(rb.len(), 2);
    }
    #[test]
    fn test_ring_buffer_empty() {
        let mut rb: RingBuffer<i32> = RingBuffer::new(3);
        assert!(rb.is_empty());
        assert_eq!(rb.pop(), None);
    }
    #[test]
    fn test_ring_buffer_capacity() {
        let rb: RingBuffer<i32> = RingBuffer::new(5);
        assert_eq!(rb.capacity(), 5);
    }
    #[test]
    fn test_ring_buffer_peek() {
        let mut rb: RingBuffer<i32> = RingBuffer::new(3);
        rb.push(42);
        assert_eq!(rb.peek(), Some(&42));
        assert_eq!(rb.len(), 1);
    }
    #[test]
    fn test_stream_writer_write() {
        let mut sw = StreamWriter::new(64, 1.0);
        sw.write_f64(3.125);
        assert!(sw.total_written > 0);
    }
    #[test]
    fn test_stream_writer_auto_flush() {
        let mut sw = StreamWriter::new(8, 1.0);
        sw.write_f64(1.0);
        sw.write_f64(2.0);
        assert!(sw.flushed_bytes() > 0);
    }
    #[test]
    fn test_stream_reader_frame() {
        let data: Vec<u8> = (0..32).map(|i| i as u8).collect();
        let mut sr = StreamReader::new(data, 8);
        let frame = sr.read_frame().unwrap();
        assert_eq!(frame.len(), 8);
        assert_eq!(sr.frames_read, 1);
    }
    #[test]
    fn test_stream_reader_all() {
        let data: Vec<u8> = vec![0u8; 40];
        let mut sr = StreamReader::new(data, 8);
        let mut count = 0;
        sr.read_all(|_| count += 1);
        assert_eq!(count, 5);
    }
    #[test]
    fn test_binary_protocol_encode_decode() {
        let proto = BinaryProtocol::oxiphysics();
        let frame = ProtocolFrame {
            timestamp: 1.5,
            seq: 42,
            payload: vec![1u8, 2, 3],
        };
        let encoded = proto.encode(&frame);
        let decoded = proto.decode(&encoded).unwrap();
        assert!((decoded.timestamp - 1.5).abs() < 1e-10);
        assert_eq!(decoded.seq, 42);
        assert_eq!(decoded.payload, vec![1u8, 2, 3]);
    }
    #[test]
    fn test_binary_protocol_wrong_magic() {
        let proto = BinaryProtocol::oxiphysics();
        let mut data = vec![0u8; 30];
        data[0] = 0xFF;
        assert!(proto.decode(&data).is_none());
    }
    #[test]
    fn test_compression_round_trip() {
        let data = vec![1.0f32, 1.0, 1.0, 1.0, 2.0, 2.0, 3.125];
        let compressed = CompressionStream::compress_f32(&data);
        let decompressed = CompressionStream::decompress_f32(&compressed);
        assert_eq!(decompressed.len(), data.len());
        for (a, b) in decompressed.iter().zip(data.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }
    #[test]
    fn test_compression_ratio_repeated() {
        let data: Vec<f32> = vec![1.0f32; 100];
        let compressed = CompressionStream::compress_f32(&data);
        let ratio = CompressionStream::compression_ratio(&data, &compressed);
        assert!(ratio < 1.0);
    }
    #[test]
    fn test_network_transport_send_receive() {
        let mut net = NetworkTransport::new(1.0, 0.0);
        net.send_packet("127.0.0.1:9000", vec![1, 2, 3]);
        net.deliver();
        let pkt = net.receive_packet().unwrap();
        assert_eq!(pkt.data, vec![1, 2, 3]);
    }
    #[test]
    fn test_network_transport_queue_depth() {
        let mut net = NetworkTransport::new(1.0, 0.0);
        net.send_packet("dst", vec![0u8; 10]);
        assert_eq!(net.queue_depth(), 1);
    }
    #[test]
    fn test_websocket_send_text() {
        let mut ws = WebSocketBridge::new(65536);
        ws.connect();
        ws.send_text("hello");
        assert_eq!(ws.pending_out(), 1);
    }
    #[test]
    fn test_websocket_ping_pong() {
        let mut ws = WebSocketBridge::new(65536);
        ws.connect();
        let ping = WsFrame {
            frame_type: WsFrameType::Ping,
            payload: b"ping".to_vec(),
            fin: true,
        };
        ws.inject_frame(ping);
        let pong = ws.outgoing.back().unwrap();
        assert_eq!(pong.frame_type, WsFrameType::Pong);
    }
    #[test]
    fn test_physics_stream_server_broadcast() {
        let mut server = PhysicsStreamServer::new(9000);
        server.start();
        server.broadcast(0.1, vec![0u8; 32]);
        assert_eq!(server.queue_depth(), 1);
    }
    #[test]
    fn test_physics_stream_server_clients() {
        let mut server = PhysicsStreamServer::new(9000);
        server.connect_client();
        server.connect_client();
        assert_eq!(server.n_clients, 2);
        server.disconnect_client();
        assert_eq!(server.n_clients, 1);
    }
    #[test]
    fn test_physics_stream_client_interpolate() {
        let mut client = PhysicsStreamClient::new(9000);
        client.connect();
        let payload: Vec<u8> = [1.0f64, 2.0f64]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        let frame = ProtocolFrame {
            timestamp: 0.1,
            seq: 1,
            payload,
        };
        client.receive_frame(frame);
        client.decode_latest();
        let interp = client.interpolated_state(0.5);
        assert!(!interp.is_empty());
    }
    #[test]
    fn test_frame_interpolator_linear() {
        let mut fi = FrameInterpolator::new(InterpolationMethod::Linear);
        fi.update(vec![0.0, 0.0], 0.0, vec![2.0, 4.0], 1.0);
        let interp = fi.interpolate(0.5);
        assert!((interp[0] - 1.0).abs() < 1e-10);
        assert!((interp[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_frame_interpolator_nearest() {
        let mut fi = FrameInterpolator::new(InterpolationMethod::Nearest);
        fi.update(vec![0.0], 0.0, vec![1.0], 1.0);
        let a = fi.interpolate(0.3);
        let b = fi.interpolate(0.7);
        assert_eq!(a, vec![0.0]);
        assert_eq!(b, vec![1.0]);
    }
    #[test]
    fn test_frame_interpolator_cubic() {
        let mut fi = FrameInterpolator::new(InterpolationMethod::CubicHermite);
        fi.update(vec![0.0], 0.0, vec![1.0], 1.0);
        let mid = fi.interpolate(0.5);
        assert!(mid[0] > 0.0 && mid[0] < 1.0);
    }
    #[test]
    fn test_stream_writer_maybe_flush() {
        let mut sw = StreamWriter::new(1024, 0.1);
        sw.write(&[0u8; 8]);
        sw.maybe_flush(0.2);
        assert!(sw.flushed_bytes() > 0);
    }
    #[test]
    fn test_stream_reader_seek() {
        let data: Vec<u8> = (0..32u8).collect();
        let mut sr = StreamReader::new(data, 8);
        sr.seek(16);
        let frame = sr.read_frame().unwrap();
        assert_eq!(frame[0], 16);
    }
}
/// Delta-encode `new` relative to `base`: output\[i\] = new\[i\] wrapping_sub base\[i\].
///
/// Suitable for trajectory data where successive frames are similar.
pub fn delta_encode(base: &[u8], new: &[u8]) -> Vec<u8> {
    let len = base.len().min(new.len());
    let mut out: Vec<u8> = (0..len).map(|i| new[i].wrapping_sub(base[i])).collect();
    if new.len() > len {
        out.extend_from_slice(&new[len..]);
    }
    out
}
/// Reconstruct `new` from `base` and a delta produced by `delta_encode`.
pub fn delta_decode(base: &[u8], delta: &[u8]) -> Vec<u8> {
    let len = base.len().min(delta.len());
    let mut out: Vec<u8> = (0..len).map(|i| base[i].wrapping_add(delta[i])).collect();
    if delta.len() > len {
        out.extend_from_slice(&delta[len..]);
    }
    out
}
/// FNV-1a 64-bit hash function.
pub fn fnv1a_hash(data: &[u8]) -> u64 {
    let mut hash: u64 = 14695981039346656037;
    for &byte in data {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(1099511628211);
    }
    hash
}
#[cfg(test)]
mod streaming_extension_tests {
    use super::*;
    use crate::streaming_io::AsyncIoStub;
    use crate::streaming_io::BackpressureBuffer;
    use crate::streaming_io::ChunkedFileBuffer;

    use crate::streaming_io::DataDeduplicator;
    use crate::streaming_io::IoStatus;
    use crate::streaming_io::ParallelIoCoordinator;
    use crate::streaming_io::RleEncoder;
    use crate::streaming_io::TrajectoryCache;
    use crate::streaming_io::VtkStreamWriter;
    use crate::streaming_io::types::*;
    #[test]
    fn test_chunked_buffer_write_read() {
        let mut cb = ChunkedFileBuffer::new(4);
        cb.write_chunk(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let c1 = cb.read_chunk().unwrap();
        assert_eq!(c1, vec![1, 2, 3, 4]);
        let c2 = cb.read_chunk().unwrap();
        assert_eq!(c2, vec![5, 6, 7, 8]);
        assert!(cb.read_chunk().is_none());
    }
    #[test]
    fn test_chunked_buffer_available_chunks() {
        let mut cb = ChunkedFileBuffer::new(8);
        cb.write_chunk(&[0u8; 24]);
        assert_eq!(cb.available_chunks(), 3);
    }
    #[test]
    fn test_chunked_buffer_reset_read() {
        let mut cb = ChunkedFileBuffer::new(4);
        cb.write_chunk(&[10, 20, 30, 40]);
        cb.read_chunk();
        assert!(cb.read_chunk().is_none());
        cb.reset_read();
        assert!(cb.read_chunk().is_some());
    }
    #[test]
    fn test_chunked_buffer_bytes_written() {
        let mut cb = ChunkedFileBuffer::new(8);
        cb.write_chunk(&[0u8; 16]);
        cb.write_chunk(&[0u8; 8]);
        assert_eq!(cb.bytes_written, 24);
    }
    #[test]
    fn test_backpressure_no_pressure_below_watermark() {
        let mut bp: BackpressureBuffer<i32> = BackpressureBuffer::new(10, 0.8, 0.2);
        for i in 0..7 {
            assert!(bp.push(i), "push {i} should succeed");
        }
        assert!(!bp.is_backpressure());
    }
    #[test]
    fn test_backpressure_activates_at_high_watermark() {
        let mut bp: BackpressureBuffer<i32> = BackpressureBuffer::new(10, 0.8, 0.2);
        for i in 0..9 {
            bp.push(i);
        }
        assert!(bp.is_backpressure() || bp.dropped > 0 || bp.len() >= 8);
    }
    #[test]
    fn test_backpressure_drops_when_active() {
        let mut bp: BackpressureBuffer<i32> = BackpressureBuffer::new(4, 0.5, 0.1);
        for i in 0..10 {
            bp.push(i);
        }
        assert!(bp.dropped > 0, "some pushes should be dropped");
    }
    #[test]
    fn test_backpressure_fill_fraction() {
        let mut bp: BackpressureBuffer<i32> = BackpressureBuffer::new(10, 0.9, 0.1);
        for i in 0..5 {
            bp.push(i);
        }
        assert!((bp.fill_fraction() - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_checkpoint_delta_encode_decode_roundtrip() {
        let base = vec![10u8, 20, 30, 40];
        let new = vec![12u8, 18, 35, 40];
        let delta = delta_encode(&base, &new);
        let recovered = delta_decode(&base, &delta);
        assert_eq!(recovered, new);
    }
    #[test]
    fn test_checkpoint_manager_first_snapshot_is_full() {
        let mut cm = CheckpointManager::full_interval(5);
        cm.record(0.1, &[1, 2, 3, 4], &[]);
        assert!(cm.snapshots[0].is_full);
    }
    #[test]
    fn test_checkpoint_manager_restore_latest() {
        let mut cm = CheckpointManager::full_interval(10);
        let state1 = vec![10u8, 20, 30];
        let state2 = vec![11u8, 19, 31];
        cm.record(0.0, &state1, &[]);
        cm.record(1.0, &state2, &state1);
        let restored = cm.restore_latest().unwrap();
        assert_eq!(restored, state2);
    }
    #[test]
    fn test_checkpoint_manager_snapshot_count() {
        let mut cm = CheckpointManager::full_interval(5);
        for i in 0..7u8 {
            let s = vec![i];
            let prev = if i > 0 { vec![i - 1] } else { vec![] };
            cm.record(i as f64, &s, &prev);
        }
        assert_eq!(cm.snapshot_count(), 7);
    }
    #[test]
    fn test_rle_encode_decode_roundtrip() {
        let data = vec![0u8, 0, 0, 0, 1, 2, 2, 2, 3];
        let enc = RleEncoder::encode(&data);
        let dec = RleEncoder::decode(&enc);
        assert_eq!(dec, data);
    }
    #[test]
    fn test_rle_encode_compresses_run() {
        let data = vec![7u8; 100];
        let enc = RleEncoder::encode(&data);
        assert!(
            enc.len() < data.len(),
            "RLE should compress repetitive data"
        );
    }
    #[test]
    fn test_rle_empty() {
        assert!(RleEncoder::encode(&[]).is_empty());
        assert!(RleEncoder::decode(&[]).is_empty());
    }
    #[test]
    fn test_csv_parser_basic() {
        let mut p = CsvStreamParser::new(',', true);
        p.feed_text("x,y,z\n1.0,2.0,3.0\n4.0,5.0,6.0");
        assert_eq!(p.row_count(), 2);
        assert_eq!(p.headers, vec!["x", "y", "z"]);
    }
    #[test]
    fn test_csv_parser_column_f64() {
        let mut p = CsvStreamParser::new(',', true);
        p.feed_text("time,energy\n0.1,10.5\n0.2,11.0\n0.3,9.8");
        let times = p.column_f64("time");
        assert_eq!(times.len(), 3);
        assert!((times[0] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_csv_parser_no_header() {
        let mut p = CsvStreamParser::new(',', false);
        p.feed_text("1,2,3\n4,5,6");
        assert_eq!(p.row_count(), 2);
        let col0 = p.column_f64_by_idx(0);
        assert!((col0[0] - 1.0).abs() < 1e-10);
        assert!((col0[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_csv_parser_missing_column_returns_empty() {
        let mut p = CsvStreamParser::new(',', true);
        p.feed_text("a,b\n1,2");
        assert!(p.column_f64("c").is_empty());
    }
    #[test]
    fn test_vtk_writer_header_contains_dimensions() {
        let mut w = VtkStreamWriter::new((5, 5, 1), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0), "test");
        w.write_header();
        let s = w.as_str();
        assert!(s.contains("DIMENSIONS 5 5 1"), "header should contain dims");
    }
    #[test]
    fn test_vtk_writer_scalar_field_increments_count() {
        let mut w = VtkStreamWriter::new((2, 2, 1), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0), "grid");
        w.write_header();
        w.write_scalar_field("pressure", &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(w.fields_written, 1);
    }
    #[test]
    fn test_vtk_writer_bytes_written_positive() {
        let mut w = VtkStreamWriter::new((2, 2, 1), (0.5, 0.5, 1.0), (0.0, 0.0, 0.0), "d");
        w.write_header();
        w.write_scalar_field("vel", &[0.0, 1.0, 2.0, 3.0]);
        assert!(w.bytes_written() > 0);
    }
    #[test]
    fn test_parallel_io_register_and_submit() {
        let mut coord = ParallelIoCoordinator::new(16);
        let w0 = coord.register_writer(0, 8);
        let w1 = coord.register_writer(8, 16);
        coord.submit(w0, &[1u8; 8]);
        coord.submit(w1, &[2u8; 8]);
        assert!(coord.all_done());
        assert_eq!(coord.total_written(), 16);
    }
    #[test]
    fn test_parallel_io_not_done_until_all_submit() {
        let mut coord = ParallelIoCoordinator::new(8);
        let _w0 = coord.register_writer(0, 4);
        let w1 = coord.register_writer(4, 8);
        assert!(!coord.all_done());
        coord.submit(w1, &[0u8; 4]);
        assert!(!coord.all_done());
    }
    #[test]
    fn test_parallel_io_data_written_correctly() {
        let mut coord = ParallelIoCoordinator::new(4);
        let w = coord.register_writer(0, 4);
        coord.submit(w, &[10, 20, 30, 40]);
        assert_eq!(coord.output, vec![10, 20, 30, 40]);
    }
    #[test]
    fn test_dedup_unique_blocks_stored_once() {
        let mut d = DataDeduplicator::new();
        let h1 = d.ingest(b"hello");
        let h2 = d.ingest(b"hello");
        assert_eq!(h1, h2, "same data → same hash");
        assert_eq!(d.unique_blocks(), 1);
        assert_eq!(d.hits, 1);
    }
    #[test]
    fn test_dedup_different_blocks() {
        let mut d = DataDeduplicator::new();
        d.ingest(b"foo");
        d.ingest(b"bar");
        assert_eq!(d.unique_blocks(), 2);
    }
    #[test]
    fn test_dedup_retrieve() {
        let mut d = DataDeduplicator::new();
        let h = d.ingest(b"physics");
        let retrieved = d.retrieve(h).unwrap();
        assert_eq!(retrieved, b"physics");
    }
    #[test]
    fn test_dedup_ratio_positive() {
        let mut d = DataDeduplicator::new();
        for _ in 0..10 {
            d.ingest(b"same block content here");
        }
        assert!(d.dedup_ratio() > 0.0);
    }
    #[test]
    fn test_traj_cache_insert_and_get() {
        let mut cache = TrajectoryCache::new(4);
        cache.insert(0, 0.0, vec![1.0, 2.0, 3.0]);
        let frame = cache.get(0).unwrap();
        assert!((frame.time - 0.0).abs() < 1e-10);
        assert_eq!(frame.positions.len(), 3);
    }
    #[test]
    fn test_traj_cache_miss() {
        let mut cache = TrajectoryCache::new(4);
        assert!(cache.get(99).is_none());
        assert_eq!(cache.misses, 1);
    }
    #[test]
    fn test_traj_cache_lru_eviction() {
        let mut cache = TrajectoryCache::new(2);
        cache.insert(0, 0.0, vec![0.0]);
        cache.insert(1, 1.0, vec![1.0]);
        cache.get(0);
        cache.insert(2, 2.0, vec![2.0]);
        assert!(cache.frames.contains_key(&0), "frame 0 should survive");
        assert!(cache.frames.contains_key(&2), "frame 2 should be present");
        assert_eq!(cache.evictions, 1);
    }
    #[test]
    fn test_traj_cache_hit_ratio() {
        let mut cache = TrajectoryCache::new(10);
        cache.insert(0, 0.0, vec![1.0]);
        cache.get(0);
        cache.get(1);
        let ratio = cache.hit_ratio();
        assert!((ratio - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_async_io_stub_pending_then_ready() {
        let mut stub = AsyncIoStub::new(vec![1, 2, 3, 4], 2);
        let (s1, _) = stub.poll_read(4);
        assert_eq!(s1, IoStatus::Pending);
        let (s2, _) = stub.poll_read(4);
        assert_eq!(s2, IoStatus::Pending);
        let (s3, data) = stub.poll_read(4);
        assert_eq!(s3, IoStatus::Ready);
        assert_eq!(data, vec![1, 2, 3, 4]);
    }
    #[test]
    fn test_async_io_stub_write() {
        let mut stub = AsyncIoStub::new(vec![], 0);
        let s = stub.poll_write(&[10, 20, 30]);
        assert_eq!(s, IoStatus::Ready);
        let (_, data) = stub.poll_read(3);
        assert_eq!(data, vec![10, 20, 30]);
    }
    #[test]
    fn test_async_io_stub_blocking_read_all() {
        let mut stub = AsyncIoStub::new(vec![5, 6, 7, 8], 3);
        let all = stub.blocking_read_all(2);
        assert_eq!(all, vec![5, 6, 7, 8]);
    }
}
/// Serialise a `SimulationChunk` to bytes (little-endian).
pub fn chunk_to_bytes(chunk: &SimulationChunk) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(&(chunk.frame_id as u64).to_le_bytes());
    buf.extend_from_slice(&(chunk.particle_count as u64).to_le_bytes());
    buf.extend_from_slice(&chunk.time.to_le_bytes());
    for p in &chunk.positions {
        for &v in p {
            buf.extend_from_slice(&v.to_le_bytes());
        }
    }
    for vel in &chunk.velocities {
        for &v in vel {
            buf.extend_from_slice(&v.to_le_bytes());
        }
    }
    buf
}
/// Deserialise a `SimulationChunk` from bytes.
pub fn chunk_from_bytes(data: &[u8]) -> Option<SimulationChunk> {
    if data.len() < 24 {
        return None;
    }
    let frame_id = u64::from_le_bytes(data[0..8].try_into().ok()?) as usize;
    let n = u64::from_le_bytes(data[8..16].try_into().ok()?) as usize;
    let time = f64::from_le_bytes(data[16..24].try_into().ok()?);
    let needed = 24 + n * 2 * 3 * 4;
    if data.len() < needed {
        return None;
    }
    let mut pos_off = 24usize;
    let mut positions = Vec::with_capacity(n);
    for _ in 0..n {
        let x = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        let y = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        let z = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        positions.push([x, y, z]);
    }
    let mut velocities = Vec::with_capacity(n);
    for _ in 0..n {
        let x = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        let y = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        let z = f32::from_le_bytes(data[pos_off..pos_off + 4].try_into().ok()?);
        pos_off += 4;
        velocities.push([x, y, z]);
    }
    Some(SimulationChunk {
        frame_id,
        particle_count: n,
        positions,
        velocities,
        time,
    })
}
/// Export particle positions from a slice of chunks to a CSV file.
///
/// Columns: `frame_id,time,particle_idx,x,y,z`.
pub fn stream_positions_to_csv(chunks: &[SimulationChunk], path: &str) -> Result<(), String> {
    use std::io::Write as IoWrite;
    let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "frame_id,time,particle_idx,x,y,z").map_err(|e| e.to_string())?;
    for chunk in chunks {
        for (i, p) in chunk.positions.iter().enumerate() {
            writeln!(
                w,
                "{},{},{},{},{},{}",
                chunk.frame_id, chunk.time, i, p[0], p[1], p[2]
            )
            .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
/// Compute streaming statistics across a slice of chunks.
///
/// Returns `(mean_kinetic_energy, max_speed, total_time)`.
///
/// `mean_kinetic_energy` = mean of ½|v|² over all particles.
/// `total_time` = time of last frame minus time of first frame.
pub fn compute_streaming_stats(chunks: &[SimulationChunk]) -> (f64, f64, f64) {
    if chunks.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut ke_sum = 0.0f64;
    let mut ke_count = 0usize;
    let mut max_speed = 0.0f64;
    for chunk in chunks {
        for vel in &chunk.velocities {
            let speed2 = vel[0] as f64 * vel[0] as f64
                + vel[1] as f64 * vel[1] as f64
                + vel[2] as f64 * vel[2] as f64;
            ke_sum += 0.5 * speed2;
            ke_count += 1;
            let speed = speed2.sqrt();
            if speed > max_speed {
                max_speed = speed;
            }
        }
    }
    let mean_ke = if ke_count > 0 {
        ke_sum / ke_count as f64
    } else {
        0.0
    };
    let total_time = chunks.last().map(|c| c.time).unwrap_or(0.0)
        - chunks.first().map(|c| c.time).unwrap_or(0.0);
    (mean_ke, max_speed, total_time)
}
/// Run-length encode a slice of `f32` values.
///
/// Consecutive equal values (bitwise) are collapsed to `(value, count)` pairs.
pub fn compress_chunk_rle(data: &[f32]) -> Vec<(f32, u32)> {
    if data.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut cur = data[0];
    let mut count = 1u32;
    for &v in &data[1..] {
        if v.to_bits() == cur.to_bits() {
            count += 1;
        } else {
            out.push((cur, count));
            cur = v;
            count = 1;
        }
    }
    out.push((cur, count));
    out
}
/// Decode a run-length encoded slice produced by [`compress_chunk_rle`].
pub fn decompress_chunk_rle(encoded: &[(f32, u32)]) -> Vec<f32> {
    let mut out = Vec::new();
    for &(v, c) in encoded {
        for _ in 0..c {
            out.push(v);
        }
    }
    out
}
#[cfg(test)]
mod chunk_tests {
    use super::*;

    use crate::streaming_io::ChunkedReader;
    use crate::streaming_io::ChunkedWriter;

    use crate::streaming_io::types::*;
    fn make_chunk(frame_id: usize, n_particles: usize) -> SimulationChunk {
        let mut c = SimulationChunk::new(frame_id, frame_id as f64 * 0.01);
        for i in 0..n_particles {
            c.add_particle([i as f32, 0.0, 0.0], [1.0, 0.0, 0.0]);
        }
        c
    }
    #[test]
    fn test_add_particle_increases_count() {
        let mut c = SimulationChunk::new(0, 0.0);
        assert_eq!(c.particle_count(), 0);
        c.add_particle([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        assert_eq!(c.particle_count(), 1);
        c.add_particle([4.0, 5.0, 6.0], [0.0, 0.0, 0.0]);
        assert_eq!(c.particle_count(), 2);
    }
    #[test]
    fn test_byte_size_positive_with_particles() {
        let c = make_chunk(0, 5);
        assert!(c.byte_size() > 0);
    }
    #[test]
    fn test_byte_size_scales_with_particle_count() {
        let c1 = make_chunk(0, 10);
        let c2 = make_chunk(0, 20);
        assert!(c2.byte_size() > c1.byte_size());
    }
    #[test]
    fn test_byte_size_empty_chunk() {
        let c = SimulationChunk::new(0, 0.0);
        assert_eq!(c.byte_size(), 24);
    }
    #[test]
    fn test_chunk_serialize_roundtrip() {
        let c = make_chunk(7, 4);
        let bytes = chunk_to_bytes(&c);
        let c2 = chunk_from_bytes(&bytes).unwrap();
        assert_eq!(c2.frame_id, 7);
        assert_eq!(c2.particle_count, 4);
        assert!((c2.time - c.time).abs() < 1e-10);
        assert!((c2.positions[0][0] - 0.0).abs() < 1e-6);
    }
    #[test]
    fn test_chunk_from_bytes_too_short() {
        assert!(chunk_from_bytes(&[0u8; 10]).is_none());
    }
    #[test]
    fn test_write_read_roundtrip() {
        let path = std::env::temp_dir().join("oxiphysics_chunk_test.bin");
        let _ = std::fs::remove_file(&path);
        let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 100);
        for i in 0..3 {
            let c = make_chunk(i, 5);
            writer.write_chunk(&c).unwrap();
        }
        assert_eq!(writer.frame_count(), 3);
        let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(reader.total_frames(), 3);
        let c = reader.read_chunk(1).unwrap();
        assert_eq!(c.frame_id, 1);
        assert_eq!(c.particle_count, 5);
    }
    #[test]
    fn test_frame_count_increases_with_writes() {
        let path = std::env::temp_dir().join("oxiphysics_chunk_count.bin");
        let _ = std::fs::remove_file(&path);
        let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 10);
        assert_eq!(writer.frame_count(), 0);
        writer.write_chunk(&make_chunk(0, 2)).unwrap();
        assert_eq!(writer.frame_count(), 1);
        writer.write_chunk(&make_chunk(1, 2)).unwrap();
        assert_eq!(writer.frame_count(), 2);
    }
    #[test]
    fn test_reader_total_frames_correct() {
        let path = std::env::temp_dir().join("oxiphysics_chunk_total.bin");
        let _ = std::fs::remove_file(&path);
        let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 10);
        for i in 0..5 {
            writer.write_chunk(&make_chunk(i, 1)).unwrap();
        }
        let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(reader.total_frames(), 5);
    }
    #[test]
    fn test_reader_out_of_range() {
        let path = std::env::temp_dir().join("oxiphysics_chunk_oor.bin");
        let _ = std::fs::remove_file(&path);
        let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 10);
        writer.write_chunk(&make_chunk(0, 1)).unwrap();
        let reader = ChunkedReader::open(path.to_str().unwrap_or("")).unwrap();
        assert!(reader.read_chunk(99).is_err());
    }
    #[test]
    fn test_finalize_ok() {
        let path = std::env::temp_dir().join("oxiphysics_finalize.bin");
        let _ = std::fs::remove_file(&path);
        let mut writer = ChunkedWriter::new(path.to_str().unwrap_or(""), 10);
        writer.write_chunk(&make_chunk(0, 1)).unwrap();
        assert!(writer.finalize().is_ok());
    }
    #[test]
    fn test_stream_positions_to_csv_creates_file() {
        let chunks = vec![make_chunk(0, 3), make_chunk(1, 3)];
        let path = std::env::temp_dir().join("oxiphysics_pos_stream.csv");
        assert!(stream_positions_to_csv(&chunks, path.to_str().unwrap_or("")).is_ok());
    }
    #[test]
    fn test_stream_positions_to_csv_row_count() {
        let chunks = vec![make_chunk(0, 4), make_chunk(1, 4)];
        let path = std::env::temp_dir().join("oxiphysics_pos_rows.csv");
        stream_positions_to_csv(&chunks, path.to_str().unwrap_or("")).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert_eq!(content.lines().count(), 9);
    }
    #[test]
    fn test_streaming_stats_empty_zero() {
        let (ke, ms, tt) = compute_streaming_stats(&[]);
        assert!((ke).abs() < 1e-15);
        assert!((ms).abs() < 1e-15);
        assert!((tt).abs() < 1e-15);
    }
    #[test]
    fn test_streaming_stats_mean_ke_positive() {
        let chunks = vec![make_chunk(0, 5)];
        let (ke, _ms, _tt) = compute_streaming_stats(&chunks);
        assert!((ke - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_streaming_stats_max_speed() {
        let mut c = SimulationChunk::new(0, 0.0);
        c.add_particle([0.0; 3], [3.0, 4.0, 0.0]);
        let (_ke, ms, _tt) = compute_streaming_stats(&[c]);
        assert!((ms - 5.0).abs() < 1e-5);
    }
    #[test]
    fn test_streaming_stats_total_time() {
        let c1 = SimulationChunk::new(0, 0.0);
        let c2 = SimulationChunk::new(1, 2.5);
        let (_ke, _ms, tt) = compute_streaming_stats(&[c1, c2]);
        assert!((tt - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_rle_compress_decompress_roundtrip() {
        let data = vec![1.0f32, 1.0, 1.0, 2.0, 3.0, 3.0];
        let enc = compress_chunk_rle(&data);
        let dec = decompress_chunk_rle(&enc);
        assert_eq!(dec, data);
    }
    #[test]
    fn test_rle_all_same() {
        let data = vec![5.0f32; 100];
        let enc = compress_chunk_rle(&data);
        assert_eq!(enc.len(), 1);
        assert_eq!(enc[0], (5.0f32, 100));
    }
    #[test]
    fn test_rle_all_different() {
        let data = vec![1.0f32, 2.0, 3.0, 4.0];
        let enc = compress_chunk_rle(&data);
        assert_eq!(enc.len(), 4);
    }
    #[test]
    fn test_rle_empty() {
        let enc = compress_chunk_rle(&[]);
        assert!(enc.is_empty());
        let dec = decompress_chunk_rle(&[]);
        assert!(dec.is_empty());
    }
    #[test]
    fn test_rle_roundtrip_mixed() {
        let data: Vec<f32> = (0..20).map(|i| (i / 3) as f32).collect();
        let enc = compress_chunk_rle(&data);
        let dec = decompress_chunk_rle(&enc);
        assert_eq!(dec, data);
    }
    #[test]
    fn test_rle_single_element() {
        let data = vec![42.0f32];
        let enc = compress_chunk_rle(&data);
        assert_eq!(enc, vec![(42.0f32, 1)]);
        assert_eq!(decompress_chunk_rle(&enc), data);
    }
}
/// Open a streaming XYZ reader for the given file path.
///
/// The file is not fully loaded; frames are read one at a time via
/// [`next_frame_xyz`].
pub fn open_streaming_xyz(path: &str) -> Result<StreamingXYZReader, std::io::Error> {
    std::fs::metadata(path)?;
    Ok(StreamingXYZReader {
        path: path.to_string(),
        current_frame: 0,
        n_atoms: 0,
        buffer: Vec::new(),
    })
}
/// Read the next frame from a streaming XYZ reader.
///
/// Returns `Some((positions, symbols))` when a frame is available, or
/// `None` when the end of file is reached.
pub fn next_frame_xyz(reader: &mut StreamingXYZReader) -> Option<(Vec<[f64; 3]>, Vec<String>)> {
    use std::io::BufRead;
    let file = std::fs::File::open(&reader.path).ok()?;
    let buf_reader = std::io::BufReader::new(file);
    let lines: Vec<String> = buf_reader.lines().map_while(Result::ok).collect();
    let mut idx = 0;
    let mut frame_count = 0;
    while idx < lines.len() {
        let n_atoms: usize = lines[idx].trim().parse().ok()?;
        idx += 1;
        if idx >= lines.len() {
            return None;
        }
        idx += 1;
        if frame_count == reader.current_frame {
            let mut positions = Vec::with_capacity(n_atoms);
            let mut symbols = Vec::with_capacity(n_atoms);
            for _ in 0..n_atoms {
                if idx >= lines.len() {
                    return None;
                }
                let parts: Vec<&str> = lines[idx].split_whitespace().collect();
                if parts.len() >= 4 {
                    symbols.push(parts[0].to_string());
                    let x: f64 = parts[1].parse().unwrap_or(0.0);
                    let y: f64 = parts[2].parse().unwrap_or(0.0);
                    let z: f64 = parts[3].parse().unwrap_or(0.0);
                    positions.push([x, y, z]);
                }
                idx += 1;
            }
            reader.n_atoms = n_atoms;
            reader.current_frame += 1;
            return Some((positions, symbols));
        } else {
            idx += n_atoms;
        }
        frame_count += 1;
    }
    None
}
/// Open a streaming XYZ writer (creates or truncates the file).
pub fn open_streaming_writer(path: &str) -> Result<StreamingXYZWriter, std::io::Error> {
    std::fs::File::create(path)?;
    Ok(StreamingXYZWriter {
        path: path.to_string(),
        frames_written: 0,
    })
}
/// Append one frame to a streaming XYZ writer.
pub fn write_frame_streaming(
    writer: &mut StreamingXYZWriter,
    positions: &[[f64; 3]],
    symbols: &[&str],
) -> Result<(), std::io::Error> {
    use std::io::Write;
    let mut file = std::fs::OpenOptions::new()
        .append(true)
        .open(&writer.path)?;
    writeln!(file, "{}", positions.len())?;
    writeln!(file, "Frame {}", writer.frames_written)?;
    let n = symbols.len().min(positions.len());
    for (sym, pos) in symbols.iter().zip(positions.iter()) {
        writeln!(file, "{} {:.6} {:.6} {:.6}", sym, pos[0], pos[1], pos[2])?;
    }
    for pos in positions.iter().skip(n) {
        writeln!(file, "X {:.6} {:.6} {:.6}", pos[0], pos[1], pos[2])?;
    }
    writer.frames_written += 1;
    Ok(())
}
/// Determine whether the current frame should be sampled.
///
/// Advances the internal counter.  Returns `true` when the frame should be
/// kept; `false` when it should be skipped.
pub fn should_sample(sampler: &mut TrajectorySampler) -> bool {
    let keep = sampler.current == 0;
    sampler.current += 1;
    if sampler.current > sampler.skip_frames {
        sampler.current = 0;
    }
    keep
}
/// Push an item onto a [`RingBuffer`], dropping it silently if the buffer is full.
pub fn ring_push<T: Clone>(rb: &mut RingBuffer<T>, item: T) {
    rb.push(item);
}
/// Pop the front item from a [`RingBuffer`], or return `None` if empty.
pub fn ring_pop<T: Clone>(rb: &mut RingBuffer<T>) -> Option<T> {
    rb.pop()
}
/// Compute statistics over all frames in a streaming XYZ file.
///
/// Returns `(frame_count, mean_position, std_position)` where each component
/// is averaged over all atoms and all frames.
pub fn streaming_statistics(reader: &mut StreamingXYZReader) -> (usize, [f64; 3], [f64; 3]) {
    reader.current_frame = 0;
    let mut count = 0usize;
    let mut sum = [0.0_f64; 3];
    let mut sum_sq = [0.0_f64; 3];
    let mut n_values = 0usize;
    while let Some((positions, _)) = next_frame_xyz(reader) {
        count += 1;
        for pos in &positions {
            for k in 0..3 {
                sum[k] += pos[k];
                sum_sq[k] += pos[k] * pos[k];
                n_values += 1;
            }
        }
        n_values -= n_values % 3;
    }
    if n_values == 0 {
        return (count, [0.0; 3], [0.0; 3]);
    }
    let n_per_component = (n_values / 3).max(1) as f64;
    let mut mean = [0.0_f64; 3];
    let mut std_dev = [0.0_f64; 3];
    for k in 0..3 {
        mean[k] = sum[k] / n_per_component;
        let variance = (sum_sq[k] / n_per_component - mean[k] * mean[k]).max(0.0);
        std_dev[k] = variance.sqrt();
    }
    (count, mean, std_dev)
}
#[cfg(test)]
mod streaming_xyz_tests {
    use super::*;

    use crate::streaming_io::types::*;
    use std::fs;
    use std::io::Write;
    fn write_test_xyz(path: &str, n_frames: usize, n_atoms: usize) {
        let mut f = fs::File::create(path).expect("create");
        for frame in 0..n_frames {
            writeln!(f, "{n_atoms}").unwrap();
            writeln!(f, "Frame {frame}").unwrap();
            for atom in 0..n_atoms {
                writeln!(
                    f,
                    "C {:.3} {:.3} {:.3}",
                    atom as f64 * 0.1,
                    frame as f64 * 0.1,
                    0.0
                )
                .unwrap();
            }
        }
    }
    #[test]
    fn test_streaming_xyz_reader_struct() {
        let tmpdir = std::env::temp_dir();
        let r = StreamingXYZReader {
            path: tmpdir.join("dummy.xyz").to_str().unwrap_or("").to_string(),
            current_frame: 0,
            n_atoms: 0,
            buffer: Vec::new(),
        };
        assert_eq!(r.current_frame, 0);
    }
    #[test]
    fn test_open_streaming_xyz_missing_file() {
        let path = std::env::temp_dir().join("nonexistent_xyz_file_12345.xyz");
        let result = open_streaming_xyz(path.to_str().unwrap_or(""));
        assert!(result.is_err());
    }
    #[test]
    fn test_open_streaming_xyz_existing() {
        let path = std::env::temp_dir().join("test_streaming_open.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 2, 3);
        let r = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        assert_eq!(r.path, path.to_str().unwrap_or(""));
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_next_frame_xyz_reads_first_frame() {
        let path = std::env::temp_dir().join("test_next_frame.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 3, 2);
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        let frame = next_frame_xyz(&mut reader);
        assert!(frame.is_some());
        let (pos, sym) = frame.unwrap();
        assert_eq!(pos.len(), 2);
        assert_eq!(sym.len(), 2);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_next_frame_xyz_advances_counter() {
        let path = std::env::temp_dir().join("test_next_frame_adv.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 3, 2);
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        next_frame_xyz(&mut reader);
        assert_eq!(reader.current_frame, 1);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_next_frame_xyz_all_frames() {
        let path = std::env::temp_dir().join("test_all_frames.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 4, 3);
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        let mut count = 0;
        while next_frame_xyz(&mut reader).is_some() {
            count += 1;
        }
        assert_eq!(count, 4);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_next_frame_xyz_none_at_eof() {
        let path = std::env::temp_dir().join("test_eof.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 2, 2);
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        next_frame_xyz(&mut reader);
        next_frame_xyz(&mut reader);
        let result = next_frame_xyz(&mut reader);
        assert!(result.is_none());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_streaming_xyz_writer_struct() {
        let tmpdir = std::env::temp_dir();
        let w = StreamingXYZWriter {
            path: tmpdir.join("dummy.xyz").to_str().unwrap_or("").to_string(),
            frames_written: 0,
        };
        assert_eq!(w.frames_written, 0);
    }
    #[test]
    fn test_open_streaming_writer_creates_file() {
        let path = std::env::temp_dir().join("test_stream_write.xyz");
        let w = open_streaming_writer(path.to_str().unwrap_or("")).expect("open writer");
        assert_eq!(w.frames_written, 0);
        assert!(path.exists());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_write_frame_streaming_increments() {
        let path = std::env::temp_dir().join("test_write_frame.xyz");
        let mut w = open_streaming_writer(path.to_str().unwrap_or("")).expect("open writer");
        let pos = vec![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let sym = vec!["C", "H"];
        write_frame_streaming(&mut w, &pos, &sym).expect("write frame");
        assert_eq!(w.frames_written, 1);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_write_read_streaming_roundtrip() {
        let path = std::env::temp_dir().join("test_stream_roundtrip.xyz");
        let mut w = open_streaming_writer(path.to_str().unwrap_or("")).expect("writer");
        let pos1 = vec![[1.0_f64, 2.0, 3.0]];
        let pos2 = vec![[4.0_f64, 5.0, 6.0]];
        write_frame_streaming(&mut w, &pos1, &["C"]).expect("frame 1");
        write_frame_streaming(&mut w, &pos2, &["N"]).expect("frame 2");
        let mut r = open_streaming_xyz(path.to_str().unwrap_or("")).expect("reader");
        let (p1, s1) = next_frame_xyz(&mut r).expect("frame 1");
        let (p2, s2) = next_frame_xyz(&mut r).expect("frame 2");
        assert!((p1[0][0] - 1.0).abs() < 1e-4);
        assert!((p2[0][0] - 4.0).abs() < 1e-4);
        assert_eq!(s1[0], "C");
        assert_eq!(s2[0], "N");
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_trajectory_sampler_new() {
        let s = TrajectorySampler::new(2);
        assert_eq!(s.skip_frames, 2);
        assert_eq!(s.current, 0);
    }
    #[test]
    fn test_should_sample_every_frame() {
        let mut s = TrajectorySampler::new(0);
        assert!(should_sample(&mut s));
        assert!(should_sample(&mut s));
        assert!(should_sample(&mut s));
    }
    #[test]
    fn test_should_sample_skip_2() {
        let mut s = TrajectorySampler::new(2);
        assert!(should_sample(&mut s));
        assert!(!should_sample(&mut s));
        assert!(!should_sample(&mut s));
        assert!(should_sample(&mut s));
    }
    #[test]
    fn test_ring_push_pop() {
        let mut rb = RingBuffer::new(3);
        ring_push(&mut rb, 10);
        ring_push(&mut rb, 20);
        assert_eq!(ring_pop(&mut rb), Some(10));
        assert_eq!(ring_pop(&mut rb), Some(20));
        assert_eq!(ring_pop::<i32>(&mut rb), None);
    }
    #[test]
    fn test_ring_push_full() {
        let mut rb = RingBuffer::new(2);
        ring_push(&mut rb, 1);
        ring_push(&mut rb, 2);
        ring_push(&mut rb, 3);
        assert_eq!(rb.len(), 2);
    }
    #[test]
    fn test_ring_pop_empty() {
        let mut rb: RingBuffer<i32> = RingBuffer::new(5);
        assert_eq!(ring_pop(&mut rb), None);
    }
    #[test]
    fn test_streaming_statistics_basic() {
        let path = std::env::temp_dir().join("test_streaming_stats.xyz");
        write_test_xyz(path.to_str().unwrap_or(""), 3, 2);
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        let (count, mean, _std) = streaming_statistics(&mut reader);
        assert_eq!(count, 3);
        assert!(mean[0] >= 0.0);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn test_streaming_statistics_empty_file() {
        let path = std::env::temp_dir().join("test_stats_empty.xyz");
        fs::write(&path, "").expect("write empty");
        let mut reader = open_streaming_xyz(path.to_str().unwrap_or("")).expect("open");
        let (count, _mean, _std) = streaming_statistics(&mut reader);
        assert_eq!(count, 0);
        let _ = fs::remove_file(&path);
    }
}
