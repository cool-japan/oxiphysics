//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, VecDeque};

use super::functions::{chunk_from_bytes, chunk_to_bytes, delta_decode, delta_encode, fnv1a_hash};

/// Run-length encoder/decoder for byte streams (trajectory data).
///
/// Format: repeated runs are encoded as `[0xFF, count, value]`; literals are
/// passed through.  This mirrors the f32 variant in `CompressionStream` but
/// operates on raw bytes.
#[derive(Debug, Clone, Default)]
pub struct RleEncoder;
impl RleEncoder {
    /// Encode a byte slice using RLE.
    pub fn encode(data: &[u8]) -> Vec<u8> {
        if data.is_empty() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut i = 0;
        while i < data.len() {
            let val = data[i];
            let mut run = 1usize;
            while i + run < data.len() && data[i + run] == val && run < 254 {
                run += 1;
            }
            if run >= 3 {
                out.push(0xFF);
                out.push(run as u8);
                out.push(val);
            } else {
                for j in 0..run {
                    out.push(data[i + j]);
                }
            }
            i += run;
        }
        out
    }
    /// Decode a byte slice produced by `encode`.
    pub fn decode(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = 0;
        while i < data.len() {
            if data[i] == 0xFF && i + 2 < data.len() {
                let count = data[i + 1] as usize;
                let val = data[i + 2];
                for _ in 0..count {
                    out.push(val);
                }
                i += 3;
            } else {
                out.push(data[i]);
                i += 1;
            }
        }
        out
    }
}
/// A cached trajectory frame.
#[derive(Debug, Clone)]
pub struct TrajFrame {
    /// Frame index.
    pub idx: usize,
    /// Simulation time.
    pub time: f64,
    /// Position data (flat: N particles × 3).
    pub positions: Vec<f64>,
    /// Last access count (used for LRU).
    pub last_access: u64,
}
/// WebSocket frame.
#[derive(Debug, Clone)]
pub struct WsFrame {
    /// Frame type.
    pub frame_type: WsFrameType,
    /// Payload data.
    pub payload: Vec<u8>,
    /// Final fragment flag.
    pub fin: bool,
}
/// Coordinates multiple parallel I/O writers that write non-overlapping data
/// regions to a shared output (e.g. a parallel file system or collective MPI
/// write).
#[derive(Debug)]
pub struct ParallelIoCoordinator {
    /// All registered writers.
    pub writers: Vec<WriterMeta>,
    /// Total output buffer (simulates shared storage).
    pub output: Vec<u8>,
    /// Total size of the output file.
    pub total_size: usize,
}
impl ParallelIoCoordinator {
    /// Create a new parallel I/O coordinator.
    pub fn new(total_size: usize) -> Self {
        Self {
            writers: Vec::new(),
            output: vec![0u8; total_size],
            total_size,
        }
    }
    /// Register a writer responsible for bytes `[start, end)`.
    ///
    /// Returns the writer ID.
    pub fn register_writer(&mut self, start: usize, end: usize) -> usize {
        let id = self.writers.len();
        self.writers.push(WriterMeta {
            id,
            byte_range: (start, end),
            bytes_written: 0,
            done: false,
        });
        id
    }
    /// Submit data from writer `id`.
    ///
    /// Writes `data` at the registered offset.  Marks the writer as done.
    pub fn submit(&mut self, id: usize, data: &[u8]) {
        if id >= self.writers.len() {
            return;
        }
        let (start, end) = self.writers[id].byte_range;
        let len = data.len().min(end.saturating_sub(start));
        if start + len <= self.output.len() {
            self.output[start..start + len].copy_from_slice(&data[..len]);
        }
        self.writers[id].bytes_written = len;
        self.writers[id].done = true;
    }
    /// Returns `true` when all writers have submitted their data.
    pub fn all_done(&self) -> bool {
        self.writers.iter().all(|w| w.done)
    }
    /// Total bytes written across all writers.
    pub fn total_written(&self) -> usize {
        self.writers.iter().map(|w| w.bytes_written).sum()
    }
}
/// Controls frame sampling from a trajectory — saves every `skip_frames+1`-th frame.
#[derive(Debug, Clone)]
pub struct TrajectorySampler {
    /// Number of frames to skip between sampled frames (0 = keep every frame).
    pub skip_frames: usize,
    /// Internal counter tracking how many frames have been seen since last sample.
    pub current: usize,
}
impl TrajectorySampler {
    /// Create a new sampler that retains every `(skip_frames+1)`-th frame.
    pub fn new(skip_frames: usize) -> Self {
        Self {
            skip_frames,
            current: 0,
        }
    }
}
/// Metadata for a parallel I/O writer.
#[derive(Debug, Clone)]
pub struct WriterMeta {
    /// Writer identifier.
    pub id: usize,
    /// Data range this writer is responsible for (byte offsets).
    pub byte_range: (usize, usize),
    /// Bytes written by this writer.
    pub bytes_written: usize,
    /// Whether this writer has finished.
    pub done: bool,
}
/// Binary protocol frame.
#[derive(Debug, Clone)]
pub struct ProtocolFrame {
    /// Frame timestamp.
    pub timestamp: f64,
    /// Frame sequence number.
    pub seq: u64,
    /// Payload data.
    pub payload: Vec<u8>,
}
/// Status of a non-blocking I/O operation.
#[derive(Debug, Clone, PartialEq)]
pub enum IoStatus {
    /// Operation completed successfully.
    Ready,
    /// Operation is pending (would block).
    Pending,
    /// An error occurred.
    Error(String),
}
/// Frame interpolator for smooth rendering between physics states.
#[derive(Debug, Clone)]
pub struct FrameInterpolator {
    /// Previous frame state.
    pub prev: Vec<f64>,
    /// Next frame state.
    pub next: Vec<f64>,
    /// Previous frame timestamp.
    pub t_prev: f64,
    /// Next frame timestamp.
    pub t_next: f64,
    /// Interpolation method.
    pub method: InterpolationMethod,
}
impl FrameInterpolator {
    /// Create a new frame interpolator.
    pub fn new(method: InterpolationMethod) -> Self {
        Self {
            prev: Vec::new(),
            next: Vec::new(),
            t_prev: 0.0,
            t_next: 1.0,
            method,
        }
    }
    /// Update with new frame pair.
    pub fn update(&mut self, prev: Vec<f64>, t_prev: f64, next: Vec<f64>, t_next: f64) {
        self.prev = prev;
        self.t_prev = t_prev;
        self.next = next;
        self.t_next = t_next;
    }
    /// Interpolate at time t.
    pub fn interpolate(&self, t: f64) -> Vec<f64> {
        let dt = self.t_next - self.t_prev;
        if dt.abs() < 1e-30 {
            return self.next.clone();
        }
        let alpha = ((t - self.t_prev) / dt).clamp(0.0, 1.0);
        match self.method {
            InterpolationMethod::Linear => self
                .prev
                .iter()
                .zip(self.next.iter())
                .map(|(&p, &n)| p + alpha * (n - p))
                .collect(),
            InterpolationMethod::CubicHermite => {
                let h00 = 2.0 * alpha * alpha * alpha - 3.0 * alpha * alpha + 1.0;
                let h01 = -2.0 * alpha * alpha * alpha + 3.0 * alpha * alpha;
                self.prev
                    .iter()
                    .zip(self.next.iter())
                    .map(|(&p, &n)| h00 * p + h01 * n)
                    .collect()
            }
            InterpolationMethod::Nearest => {
                if alpha < 0.5 {
                    self.prev.clone()
                } else {
                    self.next.clone()
                }
            }
        }
    }
    /// Error estimate between interpolated and ground truth.
    pub fn interpolation_error(&self, gt: &[f64], t: f64) -> f64 {
        let interp = self.interpolate(t);
        if gt.len() != interp.len() {
            return f64::INFINITY;
        }
        let sum: f64 = interp
            .iter()
            .zip(gt.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        (sum / interp.len().max(1) as f64).sqrt()
    }
}
/// Mock packet for network transport.
#[derive(Debug, Clone)]
pub struct Packet {
    /// Destination address (mock).
    pub dest: String,
    /// Data payload.
    pub data: Vec<u8>,
    /// Sequence number.
    pub seq: u64,
    /// Simulated latency in ms.
    pub latency_ms: f64,
}
/// Mock physics stream server.
///
/// Serves physics state frames to connected mock clients.
#[derive(Debug)]
pub struct PhysicsStreamServer {
    /// Port (mock).
    pub port: u16,
    /// Number of connected clients.
    pub n_clients: usize,
    /// Frame queue for broadcasting.
    pub frame_queue: VecDeque<ProtocolFrame>,
    /// Protocol encoder.
    pub protocol: BinaryProtocol,
    /// Sequence counter.
    pub(super) seq: u64,
    /// Running flag.
    pub running: bool,
}
impl PhysicsStreamServer {
    /// Create a new physics stream server.
    pub fn new(port: u16) -> Self {
        Self {
            port,
            n_clients: 0,
            frame_queue: VecDeque::new(),
            protocol: BinaryProtocol::oxiphysics(),
            seq: 0,
            running: false,
        }
    }
    /// Start serving.
    pub fn start(&mut self) {
        self.running = true;
    }
    /// Stop serving.
    pub fn stop(&mut self) {
        self.running = false;
    }
    /// Connect a mock client.
    pub fn connect_client(&mut self) {
        self.n_clients += 1;
    }
    /// Disconnect a mock client.
    pub fn disconnect_client(&mut self) {
        if self.n_clients > 0 {
            self.n_clients -= 1;
        }
    }
    /// Broadcast physics state to all clients.
    pub fn broadcast(&mut self, timestamp: f64, state: Vec<u8>) {
        self.seq += 1;
        let frame = ProtocolFrame {
            timestamp,
            seq: self.seq,
            payload: state,
        };
        self.frame_queue.push_back(frame);
    }
    /// Get next frame from queue.
    pub fn pop_frame(&mut self) -> Option<ProtocolFrame> {
        self.frame_queue.pop_front()
    }
    /// Queue depth.
    pub fn queue_depth(&self) -> usize {
        self.frame_queue.len()
    }
}
/// Physics state binary protocol.
///
/// Header: magic (4 bytes), version (1 byte), timestamp (8 bytes f64).
/// Body: data payload.
#[derive(Debug, Clone)]
pub struct BinaryProtocol {
    /// Magic bytes.
    pub magic: [u8; 4],
    /// Protocol version.
    pub version: u8,
}
impl BinaryProtocol {
    /// Create a new binary protocol.
    pub fn new(magic: [u8; 4], version: u8) -> Self {
        Self { magic, version }
    }
    /// Default physics protocol (magic: OXIP).
    pub fn oxiphysics() -> Self {
        Self::new(*b"OXIP", 1)
    }
    /// Encode a frame.
    pub fn encode(&self, frame: &ProtocolFrame) -> Vec<u8> {
        let payload_len = frame.payload.len() as u32;
        let mut buf = Vec::new();
        buf.extend_from_slice(&self.magic);
        buf.push(self.version);
        buf.extend_from_slice(&frame.timestamp.to_le_bytes());
        buf.extend_from_slice(&frame.seq.to_le_bytes());
        buf.extend_from_slice(&payload_len.to_le_bytes());
        buf.extend_from_slice(&frame.payload);
        buf
    }
    /// Decode a frame.
    pub fn decode(&self, data: &[u8]) -> Option<ProtocolFrame> {
        if data.len() < 25 {
            return None;
        }
        if data[0..4] != self.magic {
            return None;
        }
        if data[4] != self.version {
            return None;
        }
        let timestamp = f64::from_le_bytes(data[5..13].try_into().ok()?);
        let seq = u64::from_le_bytes(data[13..21].try_into().ok()?);
        let payload_len = u32::from_le_bytes(data[21..25].try_into().ok()?) as usize;
        if data.len() < 25 + payload_len {
            return None;
        }
        let payload = data[25..25 + payload_len].to_vec();
        Some(ProtocolFrame {
            timestamp,
            seq,
            payload,
        })
    }
    /// Frame overhead size (header bytes).
    pub fn header_size(&self) -> usize {
        25
    }
}
/// Mock TCP/UDP transport with latency simulation.
#[derive(Debug)]
pub struct NetworkTransport {
    /// Simulated latency in ms.
    pub latency_ms: f64,
    /// Packet loss rate (0.0..1.0).
    pub packet_loss: f64,
    /// Outgoing packet queue.
    pub send_queue: VecDeque<Packet>,
    /// Incoming packet buffer.
    pub recv_buffer: VecDeque<Packet>,
    /// Total packets sent.
    pub packets_sent: u64,
    /// Total packets dropped.
    pub packets_dropped: u64,
    /// Sequence counter.
    pub(super) seq: u64,
}
impl NetworkTransport {
    /// Create a new network transport.
    pub fn new(latency_ms: f64, packet_loss: f64) -> Self {
        Self {
            latency_ms,
            packet_loss,
            send_queue: VecDeque::new(),
            recv_buffer: VecDeque::new(),
            packets_sent: 0,
            packets_dropped: 0,
            seq: 0,
        }
    }
    /// Send a packet (enqueue).
    pub fn send_packet(&mut self, dest: &str, data: Vec<u8>) {
        self.seq += 1;
        let packet = Packet {
            dest: dest.to_string(),
            data,
            seq: self.seq,
            latency_ms: self.latency_ms,
        };
        if self.seq % 100 < (self.packet_loss * 100.0) as u64 {
            self.packets_dropped += 1;
            return;
        }
        self.send_queue.push_back(packet);
        self.packets_sent += 1;
    }
    /// Deliver all queued packets to recv_buffer (simulates delivery).
    pub fn deliver(&mut self) {
        while let Some(pkt) = self.send_queue.pop_front() {
            self.recv_buffer.push_back(pkt);
        }
    }
    /// Receive next packet.
    pub fn receive_packet(&mut self) -> Option<Packet> {
        self.recv_buffer.pop_front()
    }
    /// Queue depth.
    pub fn queue_depth(&self) -> usize {
        self.send_queue.len() + self.recv_buffer.len()
    }
}
/// Mock physics stream client.
///
/// Connects to a mock server, receives frames, decodes, and interpolates.
#[derive(Debug)]
pub struct PhysicsStreamClient {
    /// Server port.
    pub port: u16,
    /// Connected flag.
    pub connected: bool,
    /// Received frames buffer.
    pub frames: VecDeque<ProtocolFrame>,
    /// Protocol decoder.
    pub protocol: BinaryProtocol,
    /// Latest decoded state.
    pub latest_state: Vec<f64>,
    /// Previous state for interpolation.
    pub prev_state: Vec<f64>,
}
impl PhysicsStreamClient {
    /// Create a new physics stream client.
    pub fn new(port: u16) -> Self {
        Self {
            port,
            connected: false,
            frames: VecDeque::new(),
            protocol: BinaryProtocol::oxiphysics(),
            latest_state: Vec::new(),
            prev_state: Vec::new(),
        }
    }
    /// Connect to server.
    pub fn connect(&mut self) {
        self.connected = true;
    }
    /// Disconnect.
    pub fn disconnect(&mut self) {
        self.connected = false;
    }
    /// Inject a frame (from server).
    pub fn receive_frame(&mut self, frame: ProtocolFrame) {
        self.frames.push_back(frame);
    }
    /// Decode the latest frame.
    pub fn decode_latest(&mut self) {
        if let Some(frame) = self.frames.back() {
            self.prev_state = self.latest_state.clone();
            let n = frame.payload.len() / 8;
            self.latest_state = (0..n)
                .map(|i| {
                    f64::from_le_bytes(
                        frame.payload[i * 8..(i + 1) * 8]
                            .try_into()
                            .unwrap_or([0u8; 8]),
                    )
                })
                .collect();
        }
    }
    /// Get interpolated state at t in \[0, 1\].
    pub fn interpolated_state(&self, t: f64) -> Vec<f64> {
        if self.prev_state.len() != self.latest_state.len() {
            return self.latest_state.clone();
        }
        self.prev_state
            .iter()
            .zip(self.latest_state.iter())
            .map(|(&p, &l)| p + t * (l - p))
            .collect()
    }
    /// Number of buffered frames.
    pub fn buffered_frames(&self) -> usize {
        self.frames.len()
    }
}
/// Chunked file reader/writer that processes data in fixed-size blocks.
///
/// Suitable for large datasets that exceed available memory.  Both reading
/// and writing are abstracted as in-memory operations (no actual file I/O in
/// this mock implementation).
#[derive(Debug)]
pub struct ChunkedFileBuffer {
    /// Underlying byte storage (simulates a file).
    pub storage: Vec<u8>,
    /// Chunk size in bytes.
    pub chunk_size: usize,
    /// Current read position.
    pub read_pos: usize,
    /// Current write position.
    pub write_pos: usize,
    /// Total bytes written.
    pub bytes_written: usize,
}
impl ChunkedFileBuffer {
    /// Create a new chunked buffer with the given chunk size.
    pub fn new(chunk_size: usize) -> Self {
        Self {
            storage: Vec::new(),
            chunk_size,
            read_pos: 0,
            write_pos: 0,
            bytes_written: 0,
        }
    }
    /// Write a chunk of data to the buffer.
    pub fn write_chunk(&mut self, data: &[u8]) {
        self.storage.extend_from_slice(data);
        self.write_pos += data.len();
        self.bytes_written += data.len();
    }
    /// Read the next chunk. Returns `None` if no more data.
    pub fn read_chunk(&mut self) -> Option<Vec<u8>> {
        if self.read_pos >= self.storage.len() {
            return None;
        }
        let end = (self.read_pos + self.chunk_size).min(self.storage.len());
        let chunk = self.storage[self.read_pos..end].to_vec();
        self.read_pos = end;
        Some(chunk)
    }
    /// Number of full chunks available for reading.
    pub fn available_chunks(&self) -> usize {
        let remaining = self.storage.len().saturating_sub(self.read_pos);
        (remaining + self.chunk_size - 1) / self.chunk_size.max(1)
    }
    /// Reset read cursor to beginning.
    pub fn reset_read(&mut self) {
        self.read_pos = 0;
    }
    /// Total bytes stored.
    pub fn stored_bytes(&self) -> usize {
        self.storage.len()
    }
}
/// Streaming binary reader for files written by [`ChunkedWriter`].
#[derive(Debug)]
pub struct ChunkedReader {
    /// Path to the data file.
    pub path: String,
    /// Total number of frames discovered during [`open`](ChunkedReader::open).
    pub total_frames: usize,
    /// Byte offsets for each frame.
    pub(super) offsets: Vec<u64>,
}
impl ChunkedReader {
    /// Open the file and index all frames.
    pub fn open(path: &str) -> Result<Self, String> {
        let data = std::fs::read(path).map_err(|e| e.to_string())?;
        let mut offsets = Vec::new();
        let mut pos = 0usize;
        while pos + 8 <= data.len() {
            let frame_len = u64::from_le_bytes(
                data[pos..pos + 8]
                    .try_into()
                    .map_err(|_| "bad length field")?,
            ) as usize;
            offsets.push(pos as u64);
            pos += 8 + frame_len;
        }
        let total = offsets.len();
        Ok(Self {
            path: path.into(),
            total_frames: total,
            offsets,
        })
    }
    /// Total number of frames.
    pub fn total_frames(&self) -> usize {
        self.total_frames
    }
    /// Read the frame with the given sequential index.
    pub fn read_chunk(&self, frame_idx: usize) -> Result<SimulationChunk, String> {
        if frame_idx >= self.offsets.len() {
            return Err(format!(
                "frame_idx {frame_idx} out of range ({})",
                self.offsets.len()
            ));
        }
        use std::io::{Read, Seek, SeekFrom};
        let mut file = std::fs::File::open(&self.path).map_err(|e| e.to_string())?;
        file.seek(SeekFrom::Start(self.offsets[frame_idx]))
            .map_err(|e| e.to_string())?;
        let mut len_buf = [0u8; 8];
        file.read_exact(&mut len_buf).map_err(|e| e.to_string())?;
        let frame_len = u64::from_le_bytes(len_buf) as usize;
        let mut frame_data = vec![0u8; frame_len];
        file.read_exact(&mut frame_data)
            .map_err(|e| e.to_string())?;
        chunk_from_bytes(&frame_data).ok_or_else(|| "corrupt frame data".into())
    }
}
/// A non-blocking (async-compatible) read/write stub that simulates
/// waker-based I/O without requiring an async runtime.
///
/// In a real implementation this would wrap `tokio::fs::File` or similar.
#[derive(Debug)]
pub struct AsyncIoStub {
    /// Internal buffer.
    pub buffer: Vec<u8>,
    /// Current read position.
    pub pos: usize,
    /// Simulated readiness delay: ops will return `Pending` this many times
    /// before becoming `Ready`.
    pub pending_ticks: usize,
    /// Count of operations completed.
    pub ops_completed: usize,
}
impl AsyncIoStub {
    /// Create a new async I/O stub with an initial payload.
    pub fn new(data: Vec<u8>, pending_ticks: usize) -> Self {
        Self {
            buffer: data,
            pos: 0,
            pending_ticks,
            ops_completed: 0,
        }
    }
    /// Non-blocking read of up to `len` bytes.
    ///
    /// Returns `Pending` while `pending_ticks > 0`, `Ready` once drained.
    pub fn poll_read(&mut self, len: usize) -> (IoStatus, Vec<u8>) {
        if self.pending_ticks > 0 {
            self.pending_ticks -= 1;
            return (IoStatus::Pending, Vec::new());
        }
        let end = (self.pos + len).min(self.buffer.len());
        let data = self.buffer[self.pos..end].to_vec();
        self.pos = end;
        self.ops_completed += 1;
        (IoStatus::Ready, data)
    }
    /// Non-blocking write.
    ///
    /// Returns `Pending` while `pending_ticks > 0`.
    pub fn poll_write(&mut self, data: &[u8]) -> IoStatus {
        if self.pending_ticks > 0 {
            self.pending_ticks -= 1;
            return IoStatus::Pending;
        }
        self.buffer.extend_from_slice(data);
        self.ops_completed += 1;
        IoStatus::Ready
    }
    /// Drive the stub until `poll_read` returns `Ready`, collecting all bytes.
    pub fn blocking_read_all(&mut self, chunk_size: usize) -> Vec<u8> {
        let mut out = Vec::new();
        loop {
            let (status, chunk) = self.poll_read(chunk_size);
            if status == IoStatus::Pending {
                continue;
            }
            if chunk.is_empty() {
                break;
            }
            out.extend_from_slice(&chunk);
        }
        out
    }
}
/// A streaming CSV parser that processes data line by line.
///
/// Supports optional header rows, delimiter configuration, and type-safe
/// extraction of f64 columns.
#[derive(Debug)]
pub struct CsvStreamParser {
    /// Column delimiter character.
    pub delimiter: char,
    /// Column headers (populated if the first row is a header).
    pub headers: Vec<String>,
    /// Parsed rows (each row is a Vec of raw string values).
    pub rows: Vec<Vec<String>>,
    /// Whether the first non-empty line should be treated as a header.
    pub has_header: bool,
    /// Number of parse errors encountered.
    pub error_count: usize,
}
impl CsvStreamParser {
    /// Create a new CSV parser.
    pub fn new(delimiter: char, has_header: bool) -> Self {
        Self {
            delimiter,
            headers: Vec::new(),
            rows: Vec::new(),
            has_header,
            error_count: 0,
        }
    }
    /// Feed a single line of text to the parser.
    pub fn feed_line(&mut self, line: &str) {
        let line = line.trim();
        if line.is_empty() {
            return;
        }
        let fields: Vec<String> = line
            .split(self.delimiter)
            .map(|s| s.trim().to_string())
            .collect();
        if self.has_header && self.headers.is_empty() {
            self.headers = fields;
        } else {
            self.rows.push(fields);
        }
    }
    /// Feed multiple lines of text.
    pub fn feed_text(&mut self, text: &str) {
        for line in text.lines() {
            self.feed_line(line);
        }
    }
    /// Extract a column by name as f64 values.  Returns empty if header not found.
    pub fn column_f64(&self, name: &str) -> Vec<f64> {
        let col_idx = self.headers.iter().position(|h| h == name);
        match col_idx {
            None => Vec::new(),
            Some(idx) => self
                .rows
                .iter()
                .filter_map(|row| row.get(idx).and_then(|v| v.parse::<f64>().ok()))
                .collect(),
        }
    }
    /// Extract a column by index as f64 values.
    pub fn column_f64_by_idx(&self, idx: usize) -> Vec<f64> {
        self.rows
            .iter()
            .filter_map(|row| row.get(idx).and_then(|v| v.parse::<f64>().ok()))
            .collect()
    }
    /// Number of parsed data rows (excluding header).
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }
}
/// WebSocket message type.
#[derive(Debug, Clone, PartialEq)]
pub enum WsFrameType {
    /// Text frame.
    Text,
    /// Binary frame.
    Binary,
    /// Ping frame.
    Ping,
    /// Pong frame.
    Pong,
    /// Close frame.
    Close,
}
/// WebSocket bridge for text/binary frame framing.
#[derive(Debug)]
pub struct WebSocketBridge {
    /// Buffered outgoing frames.
    pub outgoing: VecDeque<WsFrame>,
    /// Buffered incoming frames.
    pub incoming: VecDeque<WsFrame>,
    /// Connection state.
    pub connected: bool,
    /// Max frame size.
    pub max_frame_size: usize,
}
impl WebSocketBridge {
    /// Create a new WebSocket bridge.
    pub fn new(max_frame_size: usize) -> Self {
        Self {
            outgoing: VecDeque::new(),
            incoming: VecDeque::new(),
            connected: false,
            max_frame_size,
        }
    }
    /// Connect.
    pub fn connect(&mut self) {
        self.connected = true;
    }
    /// Disconnect.
    pub fn disconnect(&mut self) {
        self.connected = false;
        let close = WsFrame {
            frame_type: WsFrameType::Close,
            payload: Vec::new(),
            fin: true,
        };
        self.outgoing.push_back(close);
    }
    /// Send text frame.
    pub fn send_text(&mut self, text: &str) {
        let frame = WsFrame {
            frame_type: WsFrameType::Text,
            payload: text.as_bytes().to_vec(),
            fin: true,
        };
        self.outgoing.push_back(frame);
    }
    /// Send binary frame.
    pub fn send_binary(&mut self, data: Vec<u8>) {
        let frame = WsFrame {
            frame_type: WsFrameType::Binary,
            payload: data,
            fin: true,
        };
        self.outgoing.push_back(frame);
    }
    /// Send ping.
    pub fn ping(&mut self, data: Vec<u8>) {
        let frame = WsFrame {
            frame_type: WsFrameType::Ping,
            payload: data,
            fin: true,
        };
        self.outgoing.push_back(frame);
    }
    /// Receive next frame.
    pub fn recv_frame(&mut self) -> Option<WsFrame> {
        self.incoming.pop_front()
    }
    /// Inject an incoming frame (for testing).
    pub fn inject_frame(&mut self, frame: WsFrame) {
        if frame.frame_type == WsFrameType::Ping {
            let pong = WsFrame {
                frame_type: WsFrameType::Pong,
                payload: frame.payload.clone(),
                fin: true,
            };
            self.outgoing.push_back(pong);
        }
        self.incoming.push_back(frame);
    }
    /// Number of pending outgoing frames.
    pub fn pending_out(&self) -> usize {
        self.outgoing.len()
    }
}
/// A simple content-addressed deduplication store.
///
/// Stores unique blocks identified by a 64-bit hash (FNV-1a). Duplicate
/// blocks are not stored twice; the caller receives back the canonical hash.
#[derive(Debug)]
pub struct DataDeduplicator {
    /// Canonical block store: hash → data.
    pub store: HashMap<u64, Vec<u8>>,
    /// Total bytes received (pre-dedup).
    pub bytes_received: usize,
    /// Total bytes actually stored (post-dedup).
    pub bytes_stored: usize,
    /// Dedup hit count.
    pub hits: usize,
}
impl DataDeduplicator {
    /// Create a new deduplicator.
    pub fn new() -> Self {
        Self {
            store: HashMap::new(),
            bytes_received: 0,
            bytes_stored: 0,
            hits: 0,
        }
    }
    /// Ingest a block of data.  Returns the content hash.
    ///
    /// If the block is new, it is stored.  If it already exists, the existing
    /// hash is returned and `hits` is incremented.
    pub fn ingest(&mut self, data: &[u8]) -> u64 {
        self.bytes_received += data.len();
        let hash = fnv1a_hash(data);
        if let std::collections::hash_map::Entry::Vacant(e) = self.store.entry(hash) {
            e.insert(data.to_vec());
            self.bytes_stored += data.len();
        } else {
            self.hits += 1;
        }
        hash
    }
    /// Retrieve a block by hash.
    pub fn retrieve(&self, hash: u64) -> Option<&[u8]> {
        self.store.get(&hash).map(|v| v.as_slice())
    }
    /// Deduplication ratio (1 - stored/received).
    pub fn dedup_ratio(&self) -> f64 {
        if self.bytes_received == 0 {
            return 0.0;
        }
        1.0 - self.bytes_stored as f64 / self.bytes_received as f64
    }
    /// Number of unique blocks stored.
    pub fn unique_blocks(&self) -> usize {
        self.store.len()
    }
}
/// A streaming writer for VTK legacy format (ASCII structured points).
///
/// Builds VTK output incrementally as physics data is produced, without
/// buffering the entire dataset in memory.
#[derive(Debug)]
pub struct VtkStreamWriter {
    /// Output buffer (simulates file I/O).
    pub output: Vec<u8>,
    /// Grid dimensions (nx, ny, nz).
    pub dims: (usize, usize, usize),
    /// Grid spacing (dx, dy, dz).
    pub spacing: (f64, f64, f64),
    /// Origin (ox, oy, oz).
    pub origin: (f64, f64, f64),
    /// Dataset name.
    pub dataset_name: String,
    /// Whether the header has been written.
    pub header_written: bool,
    /// Number of scalar fields written.
    pub fields_written: usize,
}
impl VtkStreamWriter {
    /// Create a new VTK stream writer.
    pub fn new(
        dims: (usize, usize, usize),
        spacing: (f64, f64, f64),
        origin: (f64, f64, f64),
        dataset_name: &str,
    ) -> Self {
        Self {
            output: Vec::new(),
            dims,
            spacing,
            origin,
            dataset_name: dataset_name.to_string(),
            header_written: false,
            fields_written: 0,
        }
    }
    /// Write the VTK file header (must be called before any field).
    pub fn write_header(&mut self) {
        let header = format!(
            "# vtk DataFile Version 3.0\n{}\nASCII\nDATASET STRUCTURED_POINTS\n\
             DIMENSIONS {} {} {}\nSPACING {} {} {}\nORIGIN {} {} {}\nPOINT_DATA {}\n",
            self.dataset_name,
            self.dims.0,
            self.dims.1,
            self.dims.2,
            self.spacing.0,
            self.spacing.1,
            self.spacing.2,
            self.origin.0,
            self.origin.1,
            self.origin.2,
            self.dims.0 * self.dims.1 * self.dims.2,
        );
        self.output.extend_from_slice(header.as_bytes());
        self.header_written = true;
    }
    /// Write a scalar field with the given name and data.
    ///
    /// Must call `write_header` first.
    pub fn write_scalar_field(&mut self, field_name: &str, data: &[f64]) {
        if !self.header_written {
            self.write_header();
        }
        let mut section = format!("SCALARS {} float 1\nLOOKUP_TABLE default\n", field_name);
        for &v in data {
            section.push_str(&format!("{:.6e} ", v));
        }
        section.push('\n');
        self.output.extend_from_slice(section.as_bytes());
        self.fields_written += 1;
    }
    /// Return a string view of the accumulated output.
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.output).unwrap_or("")
    }
    /// Total bytes written.
    pub fn bytes_written(&self) -> usize {
        self.output.len()
    }
}
/// Streaming XYZ writer that appends frames one at a time.
#[derive(Debug)]
pub struct StreamingXYZWriter {
    /// Path to the output XYZ file.
    pub path: String,
    /// Number of frames written so far.
    pub frames_written: usize,
}
/// Interpolation method for frame rendering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InterpolationMethod {
    /// Linear interpolation.
    Linear,
    /// Cubic Hermite spline.
    CubicHermite,
    /// Nearest neighbor.
    Nearest,
}
/// Lock-free single-producer single-consumer ring buffer.
///
/// Supports push, pop, capacity, and len operations.
#[derive(Debug)]
pub struct RingBuffer<T: Clone> {
    /// Internal storage.
    pub(super) data: VecDeque<T>,
    /// Maximum capacity.
    pub(super) capacity: usize,
}
impl<T: Clone> RingBuffer<T> {
    /// Create a new ring buffer with given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            data: VecDeque::with_capacity(capacity),
            capacity,
        }
    }
    /// Push an item. Returns true on success, false if full.
    pub fn push(&mut self, item: T) -> bool {
        if self.data.len() >= self.capacity {
            return false;
        }
        self.data.push_back(item);
        true
    }
    /// Pop an item. Returns None if empty.
    pub fn pop(&mut self) -> Option<T> {
        self.data.pop_front()
    }
    /// Current number of items.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Is the buffer empty?
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    /// Is the buffer full?
    pub fn is_full(&self) -> bool {
        self.data.len() >= self.capacity
    }
    /// Buffer capacity.
    pub fn capacity(&self) -> usize {
        self.capacity
    }
    /// Peek at front item without removing.
    pub fn peek(&self) -> Option<&T> {
        self.data.front()
    }
}
/// A bounded ring buffer that signals backpressure when it becomes too full.
///
/// The `high_watermark` fraction (e.g. 0.8) determines when backpressure is
/// activated; the `low_watermark` (e.g. 0.2) determines when it is released.
#[derive(Debug)]
pub struct BackpressureBuffer<T: Clone> {
    /// Internal ring buffer.
    pub(super) inner: RingBuffer<T>,
    /// Fraction of capacity at which backpressure is asserted (0.0–1.0).
    pub high_watermark: f64,
    /// Fraction of capacity at which backpressure is released (0.0–1.0).
    pub low_watermark: f64,
    /// Current backpressure state.
    pub backpressure: bool,
    /// Total items dropped due to backpressure.
    pub dropped: usize,
}
impl<T: Clone> BackpressureBuffer<T> {
    /// Create a new backpressure-aware buffer.
    pub fn new(capacity: usize, high_watermark: f64, low_watermark: f64) -> Self {
        Self {
            inner: RingBuffer::new(capacity),
            high_watermark,
            low_watermark,
            backpressure: false,
            dropped: 0,
        }
    }
    /// Try to push an item.  Returns `false` and increments `dropped` if
    /// backpressure is active or the buffer is full.
    pub fn push(&mut self, item: T) -> bool {
        let fill = self.inner.len() as f64 / self.inner.capacity() as f64;
        if fill >= self.high_watermark {
            self.backpressure = true;
        }
        if self.backpressure {
            self.dropped += 1;
            return false;
        }
        self.inner.push(item)
    }
    /// Pop an item.  May release backpressure if below `low_watermark`.
    pub fn pop(&mut self) -> Option<T> {
        let item = self.inner.pop();
        let fill = self.inner.len() as f64 / self.inner.capacity().max(1) as f64;
        if fill <= self.low_watermark {
            self.backpressure = false;
        }
        item
    }
    /// Current fill fraction (0.0–1.0).
    pub fn fill_fraction(&self) -> f64 {
        self.inner.len() as f64 / self.inner.capacity().max(1) as f64
    }
    /// Whether backpressure is currently active.
    pub fn is_backpressure(&self) -> bool {
        self.backpressure
    }
    /// Number of items in the buffer.
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    /// Is the buffer empty?
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}
/// Streaming binary writer for large simulation data sets.
///
/// Each call to [`write_chunk`](ChunkedWriter::write_chunk) appends a serialised
/// `SimulationChunk` to the output file.
#[derive(Debug)]
pub struct ChunkedWriter {
    /// Path to the output file.
    pub path: String,
    /// Number of frames written so far.
    pub frame_count: usize,
    /// Maximum particles per chunk (advisory; not enforced).
    pub chunk_size: usize,
}
impl ChunkedWriter {
    /// Create a new writer targeting `path` with advisory `chunk_size`.
    pub fn new(path: impl Into<String>, chunk_size: usize) -> Self {
        Self {
            path: path.into(),
            frame_count: 0,
            chunk_size,
        }
    }
    /// Append `chunk` to the file.
    ///
    /// Returns an error string on I/O failure.
    pub fn write_chunk(&mut self, chunk: &SimulationChunk) -> Result<(), String> {
        use std::io::Write as IoWrite;
        let file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)
            .map_err(|e| e.to_string())?;
        let mut writer = std::io::BufWriter::new(file);
        let data = chunk_to_bytes(chunk);
        let len = (data.len() as u64).to_le_bytes();
        writer.write_all(&len).map_err(|e| e.to_string())?;
        writer.write_all(&data).map_err(|e| e.to_string())?;
        self.frame_count += 1;
        Ok(())
    }
    /// Number of frames written.
    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
    /// Finalise the file (no-op in this implementation — flush is automatic).
    pub fn finalize(&mut self) -> Result<(), String> {
        Ok(())
    }
}
/// Lightweight LZ4-style run-length encoding for float arrays.
///
/// Encodes repeated identical values efficiently.
#[derive(Debug, Clone, Default)]
pub struct CompressionStream;
impl CompressionStream {
    /// Compress a slice of f32 values using run-length encoding.
    pub fn compress_f32(data: &[f32]) -> Vec<u8> {
        if data.is_empty() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut i = 0;
        while i < data.len() {
            let val = data[i];
            let mut run = 1usize;
            while i + run < data.len() && data[i + run] == val && run < 255 {
                run += 1;
            }
            if run > 2 {
                out.push(0xFF);
                out.push(run as u8);
                out.extend_from_slice(&val.to_le_bytes());
            } else {
                for j in 0..run {
                    out.extend_from_slice(&data[i + j].to_le_bytes());
                }
            }
            i += run;
        }
        out
    }
    /// Decompress f32 values.
    pub fn decompress_f32(data: &[u8]) -> Vec<f32> {
        let mut out = Vec::new();
        let mut i = 0;
        while i < data.len() {
            if i + 4 <= data.len() && data[i] == 0xFF && i + 1 < data.len() {
                let count = data[i + 1] as usize;
                if i + 6 <= data.len() {
                    let val =
                        f32::from_le_bytes([data[i + 2], data[i + 3], data[i + 4], data[i + 5]]);
                    for _ in 0..count {
                        out.push(val);
                    }
                    i += 6;
                    continue;
                }
            }
            if i + 4 <= data.len() {
                let val = f32::from_le_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]);
                out.push(val);
                i += 4;
            } else {
                break;
            }
        }
        out
    }
    /// Compression ratio (compressed / original).
    pub fn compression_ratio(original: &[f32], compressed: &[u8]) -> f64 {
        let orig_bytes = original.len() * 4;
        if orig_bytes == 0 {
            return 1.0;
        }
        compressed.len() as f64 / orig_bytes as f64
    }
}
/// Frame-by-frame streaming reader for XYZ trajectory files.
///
/// Reads lazily so that large files are not loaded into memory at once.
#[derive(Debug)]
pub struct StreamingXYZReader {
    /// Path to the XYZ file.
    pub path: String,
    /// Index of the next frame to be read (0-based).
    pub current_frame: usize,
    /// Number of atoms per frame (detected on first read).
    pub n_atoms: usize,
    /// Internal line buffer.
    pub buffer: Vec<u8>,
}
/// One frame of simulation data: particle positions, velocities, and metadata.
#[derive(Debug, Clone)]
pub struct SimulationChunk {
    /// Sequential frame identifier.
    pub frame_id: usize,
    /// Number of particles in this frame.
    pub particle_count: usize,
    /// Particle positions `[x, y, z]` in metres.
    pub positions: Vec<[f32; 3]>,
    /// Particle velocities `[vx, vy, vz]` in m/s.
    pub velocities: Vec<[f32; 3]>,
    /// Simulation time of this frame in seconds.
    pub time: f64,
}
impl SimulationChunk {
    /// Create a new empty chunk for `frame_id` at simulation time `time`.
    pub fn new(frame_id: usize, time: f64) -> Self {
        Self {
            frame_id,
            particle_count: 0,
            positions: Vec::new(),
            velocities: Vec::new(),
            time,
        }
    }
    /// Add a particle with position `pos` and velocity `vel`.
    pub fn add_particle(&mut self, pos: [f32; 3], vel: [f32; 3]) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.particle_count += 1;
    }
    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.particle_count
    }
    /// Approximate byte size of this chunk (positions + velocities + header).
    pub fn byte_size(&self) -> usize {
        self.particle_count * 2 * 3 * 4 + 24
    }
}
/// An incremental snapshot record for checkpoint/restart.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Snapshot sequence number.
    pub seq: u64,
    /// Simulation time at snapshot.
    pub sim_time: f64,
    /// Delta-encoded payload (diff from previous snapshot).
    pub delta_payload: Vec<u8>,
    /// Whether this is a full (base) snapshot.
    pub is_full: bool,
}
/// Manages incremental checkpoint/restart for a streaming simulation.
///
/// Full snapshots are stored periodically; incremental deltas reduce I/O cost.
#[derive(Debug)]
pub struct CheckpointManager {
    /// Stored snapshots (base + incremental).
    pub snapshots: Vec<Snapshot>,
    /// Interval (in steps) between full snapshots.
    pub full_interval: u64,
    /// Last full snapshot index in `snapshots`.
    pub last_full_idx: Option<usize>,
    /// Current step counter.
    pub step: u64,
}
impl CheckpointManager {
    /// Create a new checkpoint manager.
    pub fn full_interval(full_interval: u64) -> Self {
        Self {
            snapshots: Vec::new(),
            full_interval,
            last_full_idx: None,
            step: 0,
        }
    }
    /// Record a new snapshot.  Decides automatically whether to store a full
    /// or incremental (delta) snapshot.
    ///
    /// `state` is the raw state bytes; `prev_state` is the previous state for
    /// delta encoding (ignored for full snapshots).
    pub fn record(&mut self, sim_time: f64, state: &[u8], prev_state: &[u8]) {
        self.step += 1;
        let is_full = self.last_full_idx.is_none() || self.step.is_multiple_of(self.full_interval);
        let delta_payload = if is_full {
            state.to_vec()
        } else {
            delta_encode(prev_state, state)
        };
        let seq = self.step;
        let snap = Snapshot {
            seq,
            sim_time,
            delta_payload,
            is_full,
        };
        if is_full {
            self.last_full_idx = Some(self.snapshots.len());
        }
        self.snapshots.push(snap);
    }
    /// Restore the state at the latest checkpoint by replaying deltas.
    ///
    /// Returns `None` if no snapshots exist.
    pub fn restore_latest(&self) -> Option<Vec<u8>> {
        if self.snapshots.is_empty() {
            return None;
        }
        let base_idx = self.last_full_idx?;
        let mut state = self.snapshots[base_idx].delta_payload.clone();
        for snap in &self.snapshots[base_idx + 1..] {
            state = delta_decode(&state, &snap.delta_payload);
        }
        Some(state)
    }
    /// Number of stored snapshots.
    pub fn snapshot_count(&self) -> usize {
        self.snapshots.len()
    }
}
/// Async-style streaming writer with buffered writes and flush interval.
#[derive(Debug)]
pub struct StreamWriter {
    /// Output buffer.
    pub buffer: Vec<u8>,
    /// Maximum buffer size before auto-flush.
    pub buffer_size: usize,
    /// Flush interval in seconds.
    pub flush_interval: f64,
    /// Last flush time.
    pub last_flush: f64,
    /// Total bytes written.
    pub total_written: usize,
    /// Flushed data (in-memory output).
    pub flushed: Vec<u8>,
}
impl StreamWriter {
    /// Create a new stream writer.
    pub fn new(buffer_size: usize, flush_interval: f64) -> Self {
        Self {
            buffer: Vec::with_capacity(buffer_size),
            buffer_size,
            flush_interval,
            last_flush: 0.0,
            total_written: 0,
            flushed: Vec::new(),
        }
    }
    /// Write bytes to buffer.
    pub fn write(&mut self, data: &[u8]) {
        self.buffer.extend_from_slice(data);
        self.total_written += data.len();
        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }
    /// Write a f64 value.
    pub fn write_f64(&mut self, v: f64) {
        self.write(&v.to_le_bytes());
    }
    /// Write a u32 value.
    pub fn write_u32(&mut self, v: u32) {
        self.write(&v.to_le_bytes());
    }
    /// Flush buffer to flushed storage.
    pub fn flush(&mut self) {
        self.flushed.extend_from_slice(&self.buffer);
        self.buffer.clear();
    }
    /// Check and flush if interval elapsed.
    pub fn maybe_flush(&mut self, current_time: f64) {
        if current_time - self.last_flush >= self.flush_interval {
            self.flush();
            self.last_flush = current_time;
        }
    }
    /// Total flushed bytes.
    pub fn flushed_bytes(&self) -> usize {
        self.flushed.len()
    }
}
/// In-memory trajectory cache with LRU (least-recently-used) eviction.
///
/// Stores up to `capacity` frames.  When full, the least-recently-accessed
/// frame is evicted to make room for the new entry.
#[derive(Debug)]
pub struct TrajectoryCache {
    /// Cached frames indexed by frame index.
    pub frames: HashMap<usize, TrajFrame>,
    /// Maximum number of frames to cache.
    pub capacity: usize,
    /// Global access counter for LRU ordering.
    pub(super) access_counter: u64,
    /// Total cache misses.
    pub misses: usize,
    /// Total cache hits.
    pub hits: usize,
    /// Total evictions.
    pub evictions: usize,
}
impl TrajectoryCache {
    /// Create a new trajectory cache.
    pub fn new(capacity: usize) -> Self {
        Self {
            frames: HashMap::new(),
            capacity,
            access_counter: 0,
            misses: 0,
            hits: 0,
            evictions: 0,
        }
    }
    /// Insert a frame into the cache.  Evicts the LRU frame if at capacity.
    pub fn insert(&mut self, idx: usize, time: f64, positions: Vec<f64>) {
        if self.frames.contains_key(&idx) {
            self.access_counter += 1;
            if let Some(f) = self.frames.get_mut(&idx) {
                f.last_access = self.access_counter;
                f.positions = positions;
                f.time = time;
            }
            return;
        }
        if self.frames.len() >= self.capacity {
            self.evict_lru();
        }
        self.access_counter += 1;
        self.frames.insert(
            idx,
            TrajFrame {
                idx,
                time,
                positions,
                last_access: self.access_counter,
            },
        );
    }
    /// Look up a frame. Updates `hits`/`misses` accordingly.
    pub fn get(&mut self, idx: usize) -> Option<&TrajFrame> {
        if self.frames.contains_key(&idx) {
            self.access_counter += 1;
            let ac = self.access_counter;
            self.frames
                .get_mut(&idx)
                .expect("key must exist in map")
                .last_access = ac;
            self.hits += 1;
            self.frames.get(&idx)
        } else {
            self.misses += 1;
            None
        }
    }
    /// Evict the least-recently-used frame.
    fn evict_lru(&mut self) {
        if let Some(&evict_idx) = self
            .frames
            .iter()
            .min_by_key(|(_, f)| f.last_access)
            .map(|(k, _)| k)
        {
            self.frames.remove(&evict_idx);
            self.evictions += 1;
        }
    }
    /// Number of frames currently cached.
    pub fn len(&self) -> usize {
        self.frames.len()
    }
    /// Is the cache empty?
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
    /// Cache hit ratio.
    pub fn hit_ratio(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            return 0.0;
        }
        self.hits as f64 / total as f64
    }
}
/// Async-style streaming reader with frame callbacks.
#[derive(Debug)]
pub struct StreamReader {
    /// Input data.
    pub data: Vec<u8>,
    /// Read position.
    pub pos: usize,
    /// Frame size in bytes (fixed-size framing).
    pub frame_size: usize,
    /// Total frames read.
    pub frames_read: usize,
}
impl StreamReader {
    /// Create a new stream reader.
    pub fn new(data: Vec<u8>, frame_size: usize) -> Self {
        Self {
            data,
            pos: 0,
            frame_size,
            frames_read: 0,
        }
    }
    /// Read next frame. Returns None if insufficient data.
    pub fn read_frame(&mut self) -> Option<Vec<u8>> {
        if self.pos + self.frame_size > self.data.len() {
            return None;
        }
        let frame = self.data[self.pos..self.pos + self.frame_size].to_vec();
        self.pos += self.frame_size;
        self.frames_read += 1;
        Some(frame)
    }
    /// Read all frames with a callback.
    pub fn read_all<F: FnMut(&[u8])>(&mut self, mut callback: F) {
        while let Some(frame) = self.read_frame() {
            callback(&frame);
        }
    }
    /// Seek to byte offset.
    pub fn seek(&mut self, pos: usize) {
        self.pos = pos.min(self.data.len());
    }
    /// Remaining bytes.
    pub fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }
    /// Append more data (for streaming use).
    pub fn append_data(&mut self, more: &[u8]) {
        self.data.extend_from_slice(more);
    }
}
