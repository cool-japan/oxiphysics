#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU memory management abstractions.
//!
//! Provides CPU-side mock implementations of GPU buffer types, a pool allocator,
//! transfer queue, and memory statistics — all without requiring an actual GPU.

#![allow(dead_code)]
#![allow(missing_docs)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Buffer type tag
// ---------------------------------------------------------------------------

/// Identifies the logical type of a GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GpuBufferType {
    /// Uniform / constant buffer.
    Uniform,
    /// Large read-write storage buffer.
    Storage,
    /// Vertex data buffer.
    Vertex,
    /// Index data buffer.
    Index,
    /// Texture / image buffer.
    Texture,
    /// Generic staging (upload / download) buffer.
    Staging,
}

// ---------------------------------------------------------------------------
// GpuBuffer — typed CPU-side buffer
// ---------------------------------------------------------------------------

/// A CPU-side mock for a typed GPU buffer backed by a `Vec`u8`.
#[derive(Debug, Clone)]
pub struct GpuBuffer {
    /// Logical type.
    pub buffer_type: GpuBufferType,
    /// Raw byte data.
    data: Vec<u8>,
    /// Whether the buffer is mapped (accessible for reads/writes).
    mapped: bool,
}

impl GpuBuffer {
    /// Allocate a new zero-initialised buffer of `size` bytes.
    pub fn new(buffer_type: GpuBufferType, size: usize) -> Self {
        Self {
            buffer_type,
            data: vec![0u8; size],
            mapped: false,
        }
    }

    /// Current size in bytes.
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// Read bytes at `offset` into `dst`.  Returns the number of bytes copied.
    pub fn read(&self, offset: usize, dst: &mut [u8]) -> usize {
        if offset >= self.data.len() {
            return 0;
        }
        let len = dst.len().min(self.data.len() - offset);
        dst[..len].copy_from_slice(&self.data[offset..offset + len]);
        len
    }

    /// Write bytes from `src` at `offset`.  Returns bytes written.
    pub fn write(&mut self, offset: usize, src: &[u8]) -> usize {
        if offset >= self.data.len() {
            return 0;
        }
        let len = src.len().min(self.data.len() - offset);
        self.data[offset..offset + len].copy_from_slice(&src[..len]);
        len
    }

    /// Resize the buffer (data is zero-initialised for new regions).
    pub fn resize(&mut self, new_size: usize) {
        self.data.resize(new_size, 0);
    }

    /// Zero-fill the entire buffer.
    pub fn clear(&mut self) {
        self.data.fill(0);
    }

    /// Map the buffer for direct slice access (returns `None` if already mapped).
    pub fn map(&mut self) -> Option<&mut [u8]> {
        if self.mapped {
            return None;
        }
        self.mapped = true;
        Some(&mut self.data)
    }

    /// Unmap the buffer.
    pub fn unmap(&mut self) {
        self.mapped = false;
    }

    /// True when the buffer is currently mapped.
    pub fn is_mapped(&self) -> bool {
        self.mapped
    }
}

// ---------------------------------------------------------------------------
// GpuBufferHandle
// ---------------------------------------------------------------------------

/// An opaque handle to an allocated GPU buffer region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GpuBufferHandle {
    /// Unique 64-bit handle value.
    pub id: u64,
    /// Size of the allocation in bytes.
    pub size: usize,
    /// Byte offset within the backing pool buffer.
    pub offset: usize,
    /// Logical buffer type.
    pub buffer_type: GpuBufferType,
}

impl GpuBufferHandle {
    /// Null/invalid handle sentinel.
    pub fn null() -> Self {
        Self {
            id: 0,
            size: 0,
            offset: 0,
            buffer_type: GpuBufferType::Staging,
        }
    }

    /// True when this is the null handle.
    pub fn is_null(&self) -> bool {
        self.id == 0
    }
}

// ---------------------------------------------------------------------------
// GpuBufferPool — pool allocator
// ---------------------------------------------------------------------------

/// A pool allocator that sub-allocates from a large backing buffer.
#[derive(Debug)]
pub struct GpuBufferPool {
    /// Backing byte store.
    backing: Vec<u8>,
    /// Next free byte offset.
    next_offset: usize,
    /// All live allocations.
    allocations: HashMap<u64, GpuBufferHandle>,
    /// Next handle ID counter.
    next_id: u64,
    /// Total capacity.
    capacity: usize,
}

impl GpuBufferPool {
    /// Create a pool with `capacity` bytes.
    pub fn new(capacity: usize) -> Self {
        Self {
            backing: vec![0u8; capacity],
            next_offset: 0,
            allocations: HashMap::new(),
            next_id: 1,
            capacity,
        }
    }

    /// Allocate `size` bytes from the pool.  Returns `None` if out of memory.
    pub fn alloc(&mut self, size: usize, buffer_type: GpuBufferType) -> Option<GpuBufferHandle> {
        let aligned = align_up(size, 256);
        if self.next_offset + aligned > self.capacity {
            return None;
        }
        let handle = GpuBufferHandle {
            id: self.next_id,
            size,
            offset: self.next_offset,
            buffer_type,
        };
        self.next_id += 1;
        self.next_offset += aligned;
        self.allocations.insert(handle.id, handle);
        Some(handle)
    }

    /// Free a previously allocated handle.
    pub fn free(&mut self, handle: GpuBufferHandle) {
        self.allocations.remove(&handle.id);
    }

    /// Defragment: compact live allocations to the front of the pool.
    ///
    /// Allocations are re-packed in insertion order (by id).
    pub fn defragment(&mut self) {
        let mut sorted: Vec<GpuBufferHandle> = self.allocations.values().cloned().collect();
        sorted.sort_by_key(|h| h.id);

        let old_backing = self.backing.clone();
        self.backing.fill(0);
        let mut cursor = 0usize;

        for h in &mut sorted {
            let old_off = h.offset;
            let size = h.size;
            let aligned = align_up(size, 256);
            if old_off + size <= old_backing.len() {
                self.backing[cursor..cursor + size]
                    .copy_from_slice(&old_backing[old_off..old_off + size]);
            }
            h.offset = cursor;
            cursor += aligned;
        }
        self.next_offset = cursor;

        // Update allocations map with new offsets
        self.allocations.clear();
        for h in sorted {
            self.allocations.insert(h.id, h);
        }
    }

    /// Write bytes to a handle's region.
    pub fn write(&mut self, handle: &GpuBufferHandle, offset: usize, src: &[u8]) -> usize {
        let start = handle.offset + offset;
        if start >= self.backing.len() {
            return 0;
        }
        let len = src
            .len()
            .min(self.backing.len() - start)
            .min(handle.size.saturating_sub(offset));
        self.backing[start..start + len].copy_from_slice(&src[..len]);
        len
    }

    /// Read bytes from a handle's region.
    pub fn read(&self, handle: &GpuBufferHandle, offset: usize, dst: &mut [u8]) -> usize {
        let start = handle.offset + offset;
        if start >= self.backing.len() {
            return 0;
        }
        let len = dst
            .len()
            .min(self.backing.len() - start)
            .min(handle.size.saturating_sub(offset));
        dst[..len].copy_from_slice(&self.backing[start..start + len]);
        len
    }

    /// Number of live allocations.
    pub fn allocation_count(&self) -> usize {
        self.allocations.len()
    }

    /// Bytes used (sum of aligned allocations).
    pub fn used_bytes(&self) -> usize {
        self.next_offset
    }

    /// Free bytes remaining (approximate; does not account for freed holes before defrag).
    pub fn free_bytes(&self) -> usize {
        self.capacity.saturating_sub(self.next_offset)
    }
}

fn align_up(size: usize, align: usize) -> usize {
    size.div_ceil(align) * align
}

// ---------------------------------------------------------------------------
// UniformBuffer
// ---------------------------------------------------------------------------

/// A structured uniform buffer holding named `f64` fields.
#[derive(Debug, Clone)]
pub struct UniformBuffer {
    /// Name of this uniform block.
    pub name: String,
    /// Named field storage.
    fields: HashMap<String, f64>,
    /// Serialised byte cache (optional; computed on demand).
    byte_cache: Option<Vec<u8>>,
}

impl UniformBuffer {
    /// Create a new empty uniform buffer.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            fields: HashMap::new(),
            byte_cache: None,
        }
    }

    /// Set a named field value.
    pub fn set_field(&mut self, name: &str, value: f64) {
        self.fields.insert(name.to_owned(), value);
        self.byte_cache = None; // invalidate cache
    }

    /// Get a named field value, or `None` if not set.
    pub fn get_field(&self, name: &str) -> Option<f64> {
        self.fields.get(name).copied()
    }

    /// Serialise all fields to a byte vec (8 bytes per field, little-endian f64).
    ///
    /// Fields are sorted by name for determinism.
    pub fn serialize(&mut self) -> &[u8] {
        if self.byte_cache.is_none() {
            let mut names: Vec<&String> = self.fields.keys().collect();
            names.sort();
            let mut bytes = Vec::with_capacity(names.len() * 8);
            for name in names {
                let v = self.fields[name];
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            self.byte_cache = Some(bytes);
        }
        self.byte_cache.as_ref().expect("value should be Some")
    }

    /// Number of fields.
    pub fn field_count(&self) -> usize {
        self.fields.len()
    }
}

// ---------------------------------------------------------------------------
// StorageBuffer
// ---------------------------------------------------------------------------

/// A large read-write GPU storage buffer.
#[derive(Debug, Clone)]
pub struct StorageBuffer {
    /// Label for debugging.
    pub label: String,
    /// Raw byte data.
    data: Vec<u8>,
    /// Whether the buffer is mapped.
    mapped: bool,
}

impl StorageBuffer {
    /// Create a zeroed storage buffer.
    pub fn new(label: impl Into<String>, size: usize) -> Self {
        Self {
            label: label.into(),
            data: vec![0u8; size],
            mapped: false,
        }
    }

    /// Size in bytes.
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// Map for CPU access.  Returns `None` if already mapped.
    pub fn map(&mut self) -> Option<&mut [u8]> {
        if self.mapped {
            return None;
        }
        self.mapped = true;
        Some(&mut self.data)
    }

    /// Unmap the buffer.
    pub fn unmap(&mut self) {
        self.mapped = false;
    }

    /// Read `len` bytes starting at `offset`.
    pub fn read_bytes(&self, offset: usize, len: usize) -> Vec<u8> {
        if offset >= self.data.len() {
            return vec![];
        }
        let end = (offset + len).min(self.data.len());
        self.data[offset..end].to_vec()
    }

    /// Write `src` at `offset`.
    pub fn write_bytes(&mut self, offset: usize, src: &[u8]) -> usize {
        if offset >= self.data.len() {
            return 0;
        }
        let len = src.len().min(self.data.len() - offset);
        self.data[offset..offset + len].copy_from_slice(&src[..len]);
        len
    }

    /// Clear to zero.
    pub fn clear(&mut self) {
        self.data.fill(0);
    }
}

// ---------------------------------------------------------------------------
// VertexBuffer
// ---------------------------------------------------------------------------

/// Vertex layout descriptor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexLayout {
    /// Byte offset of the position field within a vertex.
    pub pos_offset: usize,
    /// Byte offset of the normal field within a vertex.
    pub normal_offset: usize,
    /// Byte offset of the UV field within a vertex.
    pub uv_offset: usize,
    /// Total vertex stride in bytes.
    pub stride: usize,
}

impl VertexLayout {
    /// Standard interleaved layout: `\[pos: f32×3, normal: f32×3, uv: f32×2\]`.
    pub fn standard() -> Self {
        Self {
            pos_offset: 0,
            normal_offset: 12,
            uv_offset: 24,
            stride: 32,
        }
    }
}

/// An interleaved vertex buffer (pos + normal + UV).
#[derive(Debug, Clone)]
pub struct VertexBuffer {
    /// Raw byte storage.
    data: Vec<u8>,
    /// Number of vertices.
    pub vertex_count: usize,
    /// Vertex layout.
    pub layout: VertexLayout,
}

impl VertexBuffer {
    /// Create an empty vertex buffer.
    pub fn new(layout: VertexLayout) -> Self {
        Self {
            data: vec![],
            vertex_count: 0,
            layout,
        }
    }

    /// Upload vertices from a slice of `f32` arrays `\[px, py, pz, nx, ny, nz, u, v\]`.
    pub fn upload_f32(&mut self, vertices: &[[f32; 8]]) {
        let stride = self.layout.stride;
        self.vertex_count = vertices.len();
        self.data = vec![0u8; self.vertex_count * stride];
        for (i, v) in vertices.iter().enumerate() {
            let base = i * stride;
            // Position (3 × f32 = 12 bytes)
            for j in 0..3 {
                let bytes = v[j].to_le_bytes();
                self.data[base + j * 4..base + j * 4 + 4].copy_from_slice(&bytes);
            }
            // Normal (3 × f32)
            let no = self.layout.normal_offset;
            for j in 0..3 {
                let bytes = v[3 + j].to_le_bytes();
                self.data[base + no + j * 4..base + no + j * 4 + 4].copy_from_slice(&bytes);
            }
            // UV (2 × f32)
            let uo = self.layout.uv_offset;
            for j in 0..2 {
                let bytes = v[6 + j].to_le_bytes();
                self.data[base + uo + j * 4..base + uo + j * 4 + 4].copy_from_slice(&bytes);
            }
        }
    }

    /// Size in bytes.
    pub fn byte_size(&self) -> usize {
        self.data.len()
    }

    /// Read position of vertex `i` as `\[f32; 3\]`.
    pub fn get_position(&self, i: usize) -> Option<[f32; 3]> {
        if i >= self.vertex_count {
            return None;
        }
        let base = i * self.layout.stride + self.layout.pos_offset;
        if base + 12 > self.data.len() {
            return None;
        }
        Some([
            f32::from_le_bytes(self.data[base..base + 4].try_into().ok()?),
            f32::from_le_bytes(self.data[base + 4..base + 8].try_into().ok()?),
            f32::from_le_bytes(self.data[base + 8..base + 12].try_into().ok()?),
        ])
    }
}

// ---------------------------------------------------------------------------
// IndexBuffer
// ---------------------------------------------------------------------------

/// Triangle or line primitive topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexTopology {
    /// Triangle list.
    Triangles,
    /// Line list.
    Lines,
}

/// An index buffer supporting u16 or u32 indices.
#[derive(Debug, Clone)]
pub struct IndexBuffer {
    /// Raw byte storage.
    data: Vec<u8>,
    /// Number of indices stored.
    pub index_count: usize,
    /// Whether indices are 32-bit (true) or 16-bit (false).
    pub wide: bool,
    /// Primitive topology.
    pub topology: IndexTopology,
}

impl IndexBuffer {
    /// Create an empty u32 triangle index buffer.
    pub fn new_u32(topology: IndexTopology) -> Self {
        Self {
            data: vec![],
            index_count: 0,
            wide: true,
            topology,
        }
    }

    /// Create an empty u16 triangle index buffer.
    pub fn new_u16(topology: IndexTopology) -> Self {
        Self {
            data: vec![],
            index_count: 0,
            wide: false,
            topology,
        }
    }

    /// Upload u32 indices.
    pub fn upload_u32(&mut self, indices: &[u32]) {
        self.wide = true;
        self.index_count = indices.len();
        self.data = indices.iter().flat_map(|&i| i.to_le_bytes()).collect();
    }

    /// Upload u16 indices.
    pub fn upload_u16(&mut self, indices: &[u16]) {
        self.wide = false;
        self.index_count = indices.len();
        self.data = indices.iter().flat_map(|&i| i.to_le_bytes()).collect();
    }

    /// Read the `i`-th index as u32.
    pub fn get_index(&self, i: usize) -> Option<u32> {
        if i >= self.index_count {
            return None;
        }
        if self.wide {
            let base = i * 4;
            let bytes: [u8; 4] = self.data[base..base + 4].try_into().ok()?;
            Some(u32::from_le_bytes(bytes))
        } else {
            let base = i * 2;
            let bytes: [u8; 2] = self.data[base..base + 2].try_into().ok()?;
            Some(u16::from_le_bytes(bytes) as u32)
        }
    }

    /// Size in bytes.
    pub fn byte_size(&self) -> usize {
        self.data.len()
    }
}

// ---------------------------------------------------------------------------
// TextureBuffer
// ---------------------------------------------------------------------------

/// Pixel format for a texture buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TextureFormat {
    /// 8-bit RGBA (4 bytes per pixel).
    Rgba8,
    /// 32-bit single-channel float (4 bytes per pixel).
    R32f,
    /// 32-bit two-channel float (8 bytes per pixel).
    Rg32f,
}

impl TextureFormat {
    /// Bytes per pixel.
    pub fn bytes_per_pixel(self) -> usize {
        match self {
            TextureFormat::Rgba8 => 4,
            TextureFormat::R32f => 4,
            TextureFormat::Rg32f => 8,
        }
    }
}

/// A 2-D texture buffer.
#[derive(Debug, Clone)]
pub struct TextureBuffer {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Pixel format.
    pub format: TextureFormat,
    /// Raw pixel data in row-major order.
    data: Vec<u8>,
}

impl TextureBuffer {
    /// Create a zero-initialised texture.
    pub fn new(width: usize, height: usize, format: TextureFormat) -> Self {
        let size = width * height * format.bytes_per_pixel();
        Self {
            width,
            height,
            format,
            data: vec![0u8; size],
        }
    }

    /// Total size in bytes.
    pub fn byte_size(&self) -> usize {
        self.data.len()
    }

    /// Write raw bytes for a row.
    pub fn write_row(&mut self, row: usize, src: &[u8]) -> usize {
        let bpp = self.format.bytes_per_pixel();
        let row_bytes = self.width * bpp;
        let offset = row * row_bytes;
        if offset + row_bytes > self.data.len() {
            return 0;
        }
        let len = src.len().min(row_bytes);
        self.data[offset..offset + len].copy_from_slice(&src[..len]);
        len
    }

    /// Read raw bytes for a row.
    pub fn read_row(&self, row: usize) -> Vec<u8> {
        let bpp = self.format.bytes_per_pixel();
        let row_bytes = self.width * bpp;
        let offset = row * row_bytes;
        if offset + row_bytes > self.data.len() {
            return vec![];
        }
        self.data[offset..offset + row_bytes].to_vec()
    }

    /// Write a single pixel at `(x, y)` from raw bytes.
    pub fn write_pixel(&mut self, x: usize, y: usize, pixel: &[u8]) -> bool {
        let bpp = self.format.bytes_per_pixel();
        let offset = (y * self.width + x) * bpp;
        if offset + bpp > self.data.len() || pixel.len() < bpp {
            return false;
        }
        self.data[offset..offset + bpp].copy_from_slice(&pixel[..bpp]);
        true
    }

    /// Read a single pixel at `(x, y)`.
    pub fn read_pixel(&self, x: usize, y: usize) -> Vec<u8> {
        let bpp = self.format.bytes_per_pixel();
        let offset = (y * self.width + x) * bpp;
        if offset + bpp > self.data.len() {
            return vec![];
        }
        self.data[offset..offset + bpp].to_vec()
    }

    /// Clear to zero.
    pub fn clear(&mut self) {
        self.data.fill(0);
    }
}

// ---------------------------------------------------------------------------
// GpuMemoryStats
// ---------------------------------------------------------------------------

/// Snapshot of GPU memory usage.
#[derive(Debug, Clone, Default)]
pub struct GpuMemoryStats {
    /// Total bytes allocated.
    pub total_allocated: usize,
    /// Total bytes freed / available.
    pub free_bytes: usize,
    /// Fragmentation ratio in [0, 1] (0 = no fragmentation).
    pub fragmentation: f64,
    /// Count of buffers by type.
    pub buffer_counts: HashMap<String, usize>,
}

impl GpuMemoryStats {
    /// Compute stats from a pool.
    pub fn from_pool(pool: &GpuBufferPool) -> Self {
        let total = pool.capacity;
        let used = pool.used_bytes();
        let free = pool.free_bytes();
        let live: usize = pool.allocations.values().map(|h| h.size).sum();
        let frag = if used > 0 {
            1.0 - (live as f64 / used as f64)
        } else {
            0.0
        };
        let mut buffer_counts = HashMap::new();
        for h in pool.allocations.values() {
            let key = format!("{:?}", h.buffer_type);
            *buffer_counts.entry(key).or_insert(0) += 1;
        }
        Self {
            total_allocated: total,
            free_bytes: free,
            fragmentation: frag.clamp(0.0, 1.0),
            buffer_counts,
        }
    }
}

// ---------------------------------------------------------------------------
// TransferQueue
// ---------------------------------------------------------------------------

/// A pending transfer operation.
#[derive(Debug, Clone)]
pub struct TransferOp {
    /// Source bytes to upload.
    pub src: Vec<u8>,
    /// Destination handle.
    pub dst_handle: GpuBufferHandle,
    /// Byte offset within the destination.
    pub dst_offset: usize,
    /// True for upload (CPU→GPU), false for download (GPU→CPU).
    pub upload: bool,
}

/// Pending download result.
#[derive(Debug, Clone)]
pub struct DownloadResult {
    /// The handle that was downloaded.
    pub handle: GpuBufferHandle,
    /// The downloaded bytes.
    pub data: Vec<u8>,
}

/// A simple upload / download transfer queue for a GPU pool.
#[derive(Debug)]
pub struct TransferQueue {
    /// Pending upload operations.
    uploads: Vec<TransferOp>,
    /// Completed download results.
    downloads: Vec<DownloadResult>,
}

impl TransferQueue {
    /// Create an empty transfer queue.
    pub fn new() -> Self {
        Self {
            uploads: vec![],
            downloads: vec![],
        }
    }

    /// Enqueue an upload of `data` to `handle` at `offset`.
    pub fn copy_to_gpu(&mut self, handle: GpuBufferHandle, offset: usize, data: Vec<u8>) {
        self.uploads.push(TransferOp {
            src: data,
            dst_handle: handle,
            dst_offset: offset,
            upload: true,
        });
    }

    /// Enqueue a download of `len` bytes from `handle` at `offset`.
    pub fn copy_from_gpu(&mut self, handle: GpuBufferHandle, _offset: usize, _len: usize) {
        // Simulate by recording a placeholder result
        self.downloads.push(DownloadResult {
            handle,
            data: vec![0u8; _len],
        });
    }

    /// Execute all pending transfers against the pool.
    ///
    /// Returns the number of upload operations executed.
    pub fn execute_transfers(&mut self, pool: &mut GpuBufferPool) -> usize {
        let n = self.uploads.len();
        for op in self.uploads.drain(..) {
            pool.write(&op.dst_handle, op.dst_offset, &op.src);
        }
        n
    }

    /// Drain completed download results.
    pub fn drain_downloads(&mut self) -> Vec<DownloadResult> {
        std::mem::take(&mut self.downloads)
    }

    /// Number of pending uploads.
    pub fn pending_uploads(&self) -> usize {
        self.uploads.len()
    }
}

impl Default for TransferQueue {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- GpuBuffer tests ---

    #[test]
    fn test_gpu_buffer_new() {
        let buf = GpuBuffer::new(GpuBufferType::Storage, 64);
        assert_eq!(buf.size(), 64);
        assert!(!buf.is_mapped());
    }

    #[test]
    fn test_gpu_buffer_write_read() {
        let mut buf = GpuBuffer::new(GpuBufferType::Storage, 16);
        let written = buf.write(4, &[1, 2, 3, 4]);
        assert_eq!(written, 4);
        let mut dst = [0u8; 4];
        buf.read(4, &mut dst);
        assert_eq!(dst, [1, 2, 3, 4]);
    }

    #[test]
    fn test_gpu_buffer_resize_grows() {
        let mut buf = GpuBuffer::new(GpuBufferType::Uniform, 8);
        buf.resize(32);
        assert_eq!(buf.size(), 32);
    }

    #[test]
    fn test_gpu_buffer_clear() {
        let mut buf = GpuBuffer::new(GpuBufferType::Vertex, 8);
        buf.write(0, &[1, 2, 3, 4, 5, 6, 7, 8]);
        buf.clear();
        let mut dst = [0xFFu8; 8];
        buf.read(0, &mut dst);
        assert_eq!(dst, [0u8; 8]);
    }

    #[test]
    fn test_gpu_buffer_map_unmap() {
        let mut buf = GpuBuffer::new(GpuBufferType::Staging, 4);
        assert!(buf.map().is_some());
        assert!(buf.is_mapped());
        assert!(buf.map().is_none()); // already mapped
        buf.unmap();
        assert!(!buf.is_mapped());
    }

    // --- GpuBufferPool tests ---

    #[test]
    fn test_pool_alloc_basic() {
        let mut pool = GpuBufferPool::new(4096);
        let h = pool.alloc(128, GpuBufferType::Storage);
        assert!(h.is_some());
        assert_eq!(pool.allocation_count(), 1);
    }

    #[test]
    fn test_pool_alloc_oom() {
        let mut pool = GpuBufferPool::new(256);
        let h = pool.alloc(512, GpuBufferType::Uniform);
        assert!(h.is_none());
    }

    #[test]
    fn test_pool_free() {
        let mut pool = GpuBufferPool::new(4096);
        let h = pool.alloc(64, GpuBufferType::Vertex).unwrap();
        pool.free(h);
        assert_eq!(pool.allocation_count(), 0);
    }

    #[test]
    fn test_pool_write_read() {
        let mut pool = GpuBufferPool::new(4096);
        let h = pool.alloc(256, GpuBufferType::Storage).unwrap();
        let src = [0xABu8; 8];
        pool.write(&h, 0, &src);
        let mut dst = [0u8; 8];
        pool.read(&h, 0, &mut dst);
        assert_eq!(dst, [0xABu8; 8]);
    }

    #[test]
    fn test_pool_defragment() {
        let mut pool = GpuBufferPool::new(4096);
        let h1 = pool.alloc(256, GpuBufferType::Storage).unwrap();
        let h2 = pool.alloc(256, GpuBufferType::Storage).unwrap();
        pool.free(h1);
        pool.defragment();
        // h2 should now be at offset 0 (or near it)
        let h2_new = pool
            .allocations
            .values()
            .find(|h| h.id == h2.id)
            .copied()
            .unwrap();
        assert!(h2_new.offset < 256 + 256);
    }

    // --- UniformBuffer tests ---

    #[test]
    fn test_uniform_set_get() {
        let mut ub = UniformBuffer::new("matrices");
        ub.set_field("time", 1.23);
        assert!((ub.get_field("time").unwrap() - 1.23).abs() < 1e-10);
    }

    #[test]
    fn test_uniform_get_missing() {
        let ub = UniformBuffer::new("test");
        assert!(ub.get_field("nonexistent").is_none());
    }

    #[test]
    fn test_uniform_serialize_length() {
        let mut ub = UniformBuffer::new("test");
        ub.set_field("a", 1.0);
        ub.set_field("b", 2.0);
        let bytes = ub.serialize();
        assert_eq!(bytes.len(), 16); // 2 × 8 bytes
    }

    #[test]
    fn test_uniform_field_count() {
        let mut ub = UniformBuffer::new("test");
        ub.set_field("x", 1.0);
        ub.set_field("y", 2.0);
        ub.set_field("z", 3.0);
        assert_eq!(ub.field_count(), 3);
    }

    // --- StorageBuffer tests ---

    #[test]
    fn test_storage_map_write() {
        let mut sb = StorageBuffer::new("particles", 64);
        {
            let slice = sb.map().unwrap();
            slice[0] = 42;
        }
        sb.unmap();
        let bytes = sb.read_bytes(0, 1);
        assert_eq!(bytes[0], 42);
    }

    #[test]
    fn test_storage_write_read_bytes() {
        let mut sb = StorageBuffer::new("buf", 32);
        sb.write_bytes(4, &[10, 20, 30]);
        let r = sb.read_bytes(4, 3);
        assert_eq!(r, vec![10, 20, 30]);
    }

    // --- VertexBuffer tests ---

    #[test]
    fn test_vertex_buffer_upload_read_pos() {
        let mut vb = VertexBuffer::new(VertexLayout::standard());
        let verts = [[0.0f32, 1.0, 2.0, 0.0, 1.0, 0.0, 0.5, 0.5]];
        vb.upload_f32(&verts);
        assert_eq!(vb.vertex_count, 1);
        let pos = vb.get_position(0).unwrap();
        assert!((pos[0]).abs() < 1e-6);
        assert!((pos[1] - 1.0).abs() < 1e-6);
        assert!((pos[2] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_vertex_buffer_out_of_bounds() {
        let vb = VertexBuffer::new(VertexLayout::standard());
        assert!(vb.get_position(0).is_none());
    }

    // --- IndexBuffer tests ---

    #[test]
    fn test_index_buffer_u32() {
        let mut ib = IndexBuffer::new_u32(IndexTopology::Triangles);
        ib.upload_u32(&[0, 1, 2, 3, 4, 5]);
        assert_eq!(ib.index_count, 6);
        assert_eq!(ib.get_index(2), Some(2));
    }

    #[test]
    fn test_index_buffer_u16() {
        let mut ib = IndexBuffer::new_u16(IndexTopology::Triangles);
        ib.upload_u16(&[0, 1, 2]);
        assert_eq!(ib.get_index(1), Some(1u32));
    }

    #[test]
    fn test_index_buffer_out_of_bounds() {
        let ib = IndexBuffer::new_u32(IndexTopology::Lines);
        assert!(ib.get_index(0).is_none());
    }

    // --- TextureBuffer tests ---

    #[test]
    fn test_texture_rgba8_size() {
        let tex = TextureBuffer::new(4, 4, TextureFormat::Rgba8);
        assert_eq!(tex.byte_size(), 4 * 4 * 4);
    }

    #[test]
    fn test_texture_write_read_pixel() {
        let mut tex = TextureBuffer::new(2, 2, TextureFormat::Rgba8);
        let pixel = [255u8, 0, 128, 255];
        assert!(tex.write_pixel(1, 0, &pixel));
        let read = tex.read_pixel(1, 0);
        assert_eq!(read, pixel);
    }

    #[test]
    fn test_texture_write_read_row() {
        let mut tex = TextureBuffer::new(4, 1, TextureFormat::Rgba8);
        let row = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
        tex.write_row(0, &row);
        let r = tex.read_row(0);
        assert_eq!(r, row);
    }

    #[test]
    fn test_texture_clear() {
        let mut tex = TextureBuffer::new(2, 2, TextureFormat::R32f);
        tex.write_pixel(0, 0, &1.0f32.to_le_bytes());
        tex.clear();
        let p = tex.read_pixel(0, 0);
        assert!(p.iter().all(|&b| b == 0));
    }

    // --- GpuMemoryStats tests ---

    #[test]
    fn test_memory_stats_from_pool() {
        let mut pool = GpuBufferPool::new(4096);
        pool.alloc(128, GpuBufferType::Vertex);
        pool.alloc(256, GpuBufferType::Storage);
        let stats = GpuMemoryStats::from_pool(&pool);
        assert_eq!(stats.total_allocated, 4096);
        assert!(stats.free_bytes < 4096);
    }

    // --- TransferQueue tests ---

    #[test]
    fn test_transfer_queue_upload() {
        let mut pool = GpuBufferPool::new(4096);
        let h = pool.alloc(16, GpuBufferType::Uniform).unwrap();
        let mut tq = TransferQueue::new();
        tq.copy_to_gpu(h, 0, vec![0xFFu8; 8]);
        assert_eq!(tq.pending_uploads(), 1);
        let executed = tq.execute_transfers(&mut pool);
        assert_eq!(executed, 1);
        assert_eq!(tq.pending_uploads(), 0);
        // Verify data landed
        let mut dst = [0u8; 8];
        pool.read(&h, 0, &mut dst);
        assert_eq!(dst, [0xFFu8; 8]);
    }

    #[test]
    fn test_transfer_queue_download() {
        let mut pool = GpuBufferPool::new(4096);
        let h = pool.alloc(64, GpuBufferType::Storage).unwrap();
        let mut tq = TransferQueue::new();
        tq.copy_from_gpu(h, 0, 32);
        let results = tq.drain_downloads();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].data.len(), 32);
    }
}
