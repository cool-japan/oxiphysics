// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-resident spatial-hash cell-list construction for the SPH pipeline.
//!
//! Every SPH step needs the particles bucketed by spatial-hash key: a
//! `sorted_ids[]` list (particle indices ordered by key) plus `cell_start[]` and
//! `cell_count[]` per bucket. This module builds that structure entirely on the
//! caller's [`WgpuBackendReal`]:
//!
//! 1. `sph_cell_list` writes the `(key, particle id)` pairs;
//! 2. a stable LSD radix sort ([`ResidentRadixPairs`]) orders the pairs by key.
//!    Each pass is a per-tile digit histogram, a multi-level GPU exclusive scan
//!    ([`ResidentScan`]) over the per-tile counts, and a stable scatter;
//! 3. `sph_cell_histogram` counts the keys per bucket and a second
//!    [`ResidentScan`] turns the counts into the bucket start offsets.
//!
//! Nothing is read back to the host between the stages, no device is created and
//! nothing is allocated per call: every buffer (including one parameter buffer
//! per dispatch shape) is created once in the constructors, so the backend's
//! buffer table stays constant while a simulation runs and the pipelines compiled
//! by the first step are reused by every later one.

use crate::compute::wgpu_backend::real::WgpuBackendReal;
use crate::compute::wgpu_backend::{WgpuBufferHandle, WgpuInitError};
use crate::kernels_wgsl::{
    RADIX_HISTOGRAM_WGSL, RADIX_SCATTER_PAIRS_WGSL, SCAN_ADD_WGSL, SCAN_WGSL,
    SPH_CELL_HISTOGRAM_WGSL, SPH_CELL_LIST_WGSL,
};

/// Width of one scan block / radix tile: the `@workgroup_size` of those kernels.
const BLOCK: usize = 256;
/// Number of digit values of one radix pass (`2^RADIX_DIGIT_BITS`).
const RADIX_BINS: usize = 256;
/// Bits consumed by one radix pass.
const RADIX_DIGIT_BITS: u32 = 8;
/// Particles handled by one workgroup of the `sph_cell_list` kernel.
const CELL_LIST_WORKGROUP: usize = 64;
/// Largest workgroup count per dispatch dimension that every wgpu backend accepts.
const MAX_WORKGROUPS: usize = 65_535;

fn ro() -> wgpu::BufferBindingType {
    wgpu::BufferBindingType::Storage { read_only: true }
}

fn rw() -> wgpu::BufferBindingType {
    wgpu::BufferBindingType::Storage { read_only: false }
}

/// Byte size of `words` 32-bit words.
fn words_to_bytes(words: usize) -> u64 {
    (words as u64) * 4
}

/// Number of `u32` words a scan input/output buffer must hold to scan `len`
/// elements: the scan kernel reads whole 256-wide blocks.
pub(super) fn scan_capacity_words(len: usize) -> usize {
    len.max(1).div_ceil(BLOCK) * BLOCK
}

/// Number of 8-bit radix passes needed to sort keys in `[0, key_limit)`.
///
/// At least one pass is always run so that the sorted data of an odd pass count
/// lands in the second buffer pair and the accessors stay uniform.
pub(super) fn radix_passes_for(key_limit: usize) -> u32 {
    let max_key = key_limit.saturating_sub(1);
    let bits = usize::BITS - max_key.leading_zeros();
    bits.div_ceil(RADIX_DIGIT_BITS).clamp(1, 4)
}

// ── Multi-level exclusive scan ───────────────────────────────────────────────

/// One level of a [`ResidentScan`].
struct ScanLevel {
    /// Data scanned by this level (level 0: the caller's input buffer; deeper
    /// levels: the block totals of the level below).
    input: WgpuBufferHandle,
    /// Exclusive scan of `input` (without the offsets of the levels above until
    /// the add pass has run).
    output: WgpuBufferHandle,
    /// Per-block totals written by the scan kernel (one word per workgroup).
    block_sums: WgpuBufferHandle,
    /// `{ n, 0, 0, 0 }` parameter words, shared by the scan and add kernels.
    params: WgpuBufferHandle,
    /// Number of 256-wide blocks (= workgroups) of this level.
    blocks: usize,
}

/// Exclusive prefix sum of a GPU-resident `u32` buffer of any length.
///
/// The kernels scan one 256-wide block per workgroup, so longer inputs are
/// handled hierarchically: the block totals of a level are themselves scanned by
/// the next level until a single block remains, and the scanned totals are then
/// added back down the hierarchy. All level buffers are allocated once by
/// [`ResidentScan::new`]; [`ResidentScan::run`] only issues dispatches.
pub(super) struct ResidentScan {
    levels: Vec<ScanLevel>,
}

impl ResidentScan {
    /// Build the level hierarchy for scanning `len` words of `input` into
    /// `output`.
    ///
    /// Both buffers must hold at least [`scan_capacity_words`]`(len)` words
    /// (`input` is read in whole blocks, so its tail must be zero or irrelevant);
    /// otherwise `InvalidHandle` is returned. `len == 0` yields a no-op scan.
    pub(super) fn new(
        backend: &mut WgpuBackendReal,
        input: WgpuBufferHandle,
        output: WgpuBufferHandle,
        len: usize,
    ) -> Result<Self, WgpuInitError> {
        if len == 0 {
            return Ok(Self { levels: Vec::new() });
        }
        let capacity = words_to_bytes(scan_capacity_words(len));
        for handle in [input, output] {
            match backend.buffer_size(handle) {
                Some(size) if size >= capacity => {}
                _ => return Err(WgpuInitError::InvalidHandle(handle.0)),
            }
        }
        if len.div_ceil(BLOCK) > MAX_WORKGROUPS {
            return Err(WgpuInitError::FeatureDisabled);
        }

        let mut levels = Vec::new();
        let mut level_len = len;
        let mut level_input = input;
        let mut level_output = output;
        loop {
            let blocks = level_len.div_ceil(BLOCK);
            // The block totals double as the next level's input, which the scan
            // kernel reads in whole blocks.
            let block_sums =
                backend.create_buffer_storage(words_to_bytes(scan_capacity_words(blocks)));
            let params = backend.create_buffer_storage(16);
            backend.queue_write_buffer_raw(
                &params,
                bytemuck::cast_slice(&[level_len as u32, 0, 0, 0]),
            );
            levels.push(ScanLevel {
                input: level_input,
                output: level_output,
                block_sums,
                params,
                blocks,
            });
            if blocks == 1 {
                break;
            }
            level_len = blocks;
            level_input = block_sums;
            level_output =
                backend.create_buffer_storage(words_to_bytes(scan_capacity_words(blocks)));
        }
        Ok(Self { levels })
    }

    /// Number of hierarchy levels (1 for up to 256 elements, 2 for up to 65536,
    /// ...).
    #[cfg(test)]
    pub(super) fn level_count(&self) -> usize {
        self.levels.len()
    }

    /// Run the scan: per-block scans bottom-up, then the block offsets added
    /// back top-down.
    pub(super) fn run(&self, backend: &WgpuBackendReal) -> Result<(), WgpuInitError> {
        for level in &self.levels {
            backend.dispatch_wgsl(
                SCAN_WGSL,
                "scan_block",
                &[
                    (level.input, ro()),
                    (level.output, rw()),
                    (level.block_sums, rw()),
                    (level.params, ro()),
                ],
                [level.blocks as u32, 1, 1],
            )?;
        }
        for pair in self.levels.windows(2).rev() {
            if let [lower, upper] = pair {
                backend.dispatch_wgsl(
                    SCAN_ADD_WGSL,
                    "add_block_offsets",
                    &[
                        (lower.output, rw()),
                        (upper.output, ro()),
                        (lower.params, ro()),
                    ],
                    [lower.blocks as u32, 1, 1],
                )?;
            }
        }
        Ok(())
    }
}

// ── Stable LSD radix sort of (key, payload) pairs ────────────────────────────

/// Stable LSD radix sort of `n` `(key, payload)` `u32` pairs held in GPU buffers.
///
/// The pairs are written into [`input_keys`](Self::input_keys) /
/// [`input_payload`](Self::input_payload) by the caller (for the SPH pipeline,
/// directly by the `sph_cell_list` kernel); [`run`](Self::run) ping-pongs them
/// between two buffer pairs, one pass per 8 key bits, and the sorted result is
/// then available through [`sorted_keys`](Self::sorted_keys) /
/// [`sorted_payload`](Self::sorted_payload). The scan of the per-tile digit
/// counts stays on the GPU ([`ResidentScan`]); there is no host readback.
pub(super) struct ResidentRadixPairs {
    num_tiles: usize,
    keys: [WgpuBufferHandle; 2],
    payload: [WgpuBufferHandle; 2],
    /// Digit-major per-tile counts: `tile_hist[digit * num_tiles + tile]`.
    tile_hist: WgpuBufferHandle,
    /// Exclusive scan of `tile_hist`: the output base of each `(digit, tile)`.
    base: WgpuBufferHandle,
    /// One `{ n, shift, num_tiles, 0 }` parameter buffer per pass.
    pass_params: Vec<WgpuBufferHandle>,
    scan: ResidentScan,
}

impl ResidentRadixPairs {
    /// Allocate everything needed to sort `n >= 1` pairs whose keys are below
    /// `key_limit` (this decides the number of passes).
    pub(super) fn new(
        backend: &mut WgpuBackendReal,
        n: usize,
        key_limit: usize,
    ) -> Result<Self, WgpuInitError> {
        if n == 0 {
            return Err(WgpuInitError::FeatureDisabled);
        }
        let num_tiles = n.div_ceil(BLOCK);
        if num_tiles > MAX_WORKGROUPS {
            return Err(WgpuInitError::FeatureDisabled);
        }
        let padded = num_tiles * BLOCK;
        let hist_len = RADIX_BINS * num_tiles;

        let keys = [
            backend.create_buffer_storage(words_to_bytes(padded)),
            backend.create_buffer_storage(words_to_bytes(padded)),
        ];
        let payload = [
            backend.create_buffer_storage(words_to_bytes(padded)),
            backend.create_buffer_storage(words_to_bytes(padded)),
        ];
        let tile_hist = backend.create_buffer_storage(words_to_bytes(hist_len));
        let base = backend.create_buffer_storage(words_to_bytes(hist_len));
        let scan = ResidentScan::new(backend, tile_hist, base, hist_len)?;

        let passes = radix_passes_for(key_limit);
        let mut pass_params = Vec::with_capacity(passes as usize);
        for pass in 0..passes {
            let params = backend.create_buffer_storage(16);
            let words = [n as u32, pass * RADIX_DIGIT_BITS, num_tiles as u32, 0];
            backend.queue_write_buffer_raw(&params, bytemuck::cast_slice(&words));
            pass_params.push(params);
        }

        Ok(Self {
            num_tiles,
            keys,
            payload,
            tile_hist,
            base,
            pass_params,
            scan,
        })
    }

    /// Buffer the caller fills with the unsorted keys.
    pub(super) fn input_keys(&self) -> WgpuBufferHandle {
        self.keys[0]
    }

    /// Buffer the caller fills with the payload that follows each key.
    pub(super) fn input_payload(&self) -> WgpuBufferHandle {
        self.payload[0]
    }

    /// Number of radix passes [`run`](Self::run) performs.
    pub(super) fn pass_count(&self) -> usize {
        self.pass_params.len()
    }

    /// Buffer holding the sorted keys after [`run`](Self::run).
    pub(super) fn sorted_keys(&self) -> WgpuBufferHandle {
        self.keys[self.pass_count() % 2]
    }

    /// Buffer holding the payload in sorted-key order after [`run`](Self::run).
    pub(super) fn sorted_payload(&self) -> WgpuBufferHandle {
        self.payload[self.pass_count() % 2]
    }

    /// Sort the pairs currently in the input buffers.
    pub(super) fn run(&self, backend: &WgpuBackendReal) -> Result<(), WgpuInitError> {
        let workgroups = [self.num_tiles as u32, 1, 1];
        for (pass, &params) in self.pass_params.iter().enumerate() {
            let src = pass % 2;
            let dst = 1 - src;
            backend.dispatch_wgsl(
                RADIX_HISTOGRAM_WGSL,
                "radix_histogram",
                &[
                    (self.keys[src], ro()),
                    (self.tile_hist, rw()),
                    (params, ro()),
                ],
                workgroups,
            )?;
            self.scan.run(backend)?;
            backend.dispatch_wgsl(
                RADIX_SCATTER_PAIRS_WGSL,
                "radix_scatter_pairs",
                &[
                    (self.keys[src], ro()),
                    (self.payload[src], ro()),
                    (self.keys[dst], rw()),
                    (self.payload[dst], rw()),
                    (self.base, ro()),
                    (params, ro()),
                ],
                workgroups,
            )?;
        }
        Ok(())
    }
}

// ── Spatial-hash cell list ───────────────────────────────────────────────────

/// Host copy of a built cell list (see [`ResidentCellList::read_back`]).
pub(super) struct CellListHost {
    /// Particle indices ordered by hash key (length = particle count).
    pub(super) sorted_ids: Vec<u32>,
    /// Exclusive-scan start offset into `sorted_ids` per bucket.
    pub(super) cell_start: Vec<u32>,
    /// Particle count per bucket.
    pub(super) cell_count: Vec<u32>,
}

/// GPU-resident spatial-hash cell list over `n` particles and `num_cells`
/// buckets: hash kernel, radix sort, bucket histogram and bucket scan, all on the
/// caller's backend with every buffer preallocated.
pub(super) struct ResidentCellList {
    n: usize,
    num_cells: usize,
    sort: ResidentRadixPairs,
    /// Per-bucket particle count (capacity padded to whole scan blocks).
    cell_count: WgpuBufferHandle,
    /// Per-bucket start offset: exclusive scan of `cell_count`.
    cell_start: WgpuBufferHandle,
    cell_scan: ResidentScan,
    /// `CellParams` words of the hash kernel (rewritten every build).
    cell_params: WgpuBufferHandle,
    /// `BinParams` words of the histogram kernel (constant).
    bin_params: WgpuBufferHandle,
    /// Host zeros used to clear `cell_count` before each histogram.
    zero_counts: Vec<u8>,
}

impl ResidentCellList {
    /// Allocate the cell list for `n >= 1` particles and `num_cells >= 1`
    /// buckets.
    pub(super) fn new(
        backend: &mut WgpuBackendReal,
        n: usize,
        num_cells: usize,
    ) -> Result<Self, WgpuInitError> {
        if n == 0 || num_cells == 0 || n.div_ceil(CELL_LIST_WORKGROUP) > MAX_WORKGROUPS {
            return Err(WgpuInitError::FeatureDisabled);
        }
        let sort = ResidentRadixPairs::new(backend, n, num_cells)?;
        let capacity = words_to_bytes(scan_capacity_words(num_cells));
        let cell_count = backend.create_buffer_storage(capacity);
        let cell_start = backend.create_buffer_storage(capacity);
        let cell_scan = ResidentScan::new(backend, cell_count, cell_start, num_cells)?;
        let cell_params = backend.create_buffer_storage(32);
        let bin_params = backend.create_buffer_storage(16);
        backend.queue_write_buffer_raw(
            &bin_params,
            bytemuck::cast_slice(&[n as u32, num_cells as u32, 0, 0]),
        );
        Ok(Self {
            n,
            num_cells,
            sort,
            cell_count,
            cell_start,
            cell_scan,
            cell_params,
            bin_params,
            zero_counts: vec![0u8; num_cells * 4],
        })
    }

    /// Number of hash buckets.
    pub(super) fn num_cells(&self) -> usize {
        self.num_cells
    }

    /// Buffer with the particle indices ordered by key (valid after a build).
    pub(super) fn sorted_ids(&self) -> WgpuBufferHandle {
        self.sort.sorted_payload()
    }

    /// Per-bucket start offsets into [`sorted_ids`](Self::sorted_ids).
    pub(super) fn cell_start(&self) -> WgpuBufferHandle {
        self.cell_start
    }

    /// Per-bucket particle counts.
    pub(super) fn cell_count(&self) -> WgpuBufferHandle {
        self.cell_count
    }

    /// Rebuild the cell list from `positions` (xyz-interleaved `f32`, `3 * n`
    /// values) on a grid of spacing `cell_size` anchored at `origin`.
    pub(super) fn build(
        &self,
        backend: &WgpuBackendReal,
        positions: WgpuBufferHandle,
        cell_size: f32,
        origin: [f32; 3],
    ) -> Result<(), WgpuInitError> {
        let words = [
            self.n as u32,
            self.num_cells as u32,
            cell_size.to_bits(),
            origin[0].to_bits(),
            origin[1].to_bits(),
            origin[2].to_bits(),
            0,
            0,
        ];
        backend.queue_write_buffer_raw(&self.cell_params, bytemuck::cast_slice(&words));
        backend.dispatch_wgsl(
            SPH_CELL_LIST_WGSL,
            "sph_cell_list",
            &[
                (positions, ro()),
                (self.sort.input_keys(), rw()),
                (self.sort.input_payload(), rw()),
                (self.cell_params, ro()),
            ],
            [self.n.div_ceil(CELL_LIST_WORKGROUP) as u32, 1, 1],
        )?;
        self.index_keys(backend)
    }

    /// Sort the `(key, id)` pairs currently in the sorter's input buffers and
    /// derive `cell_count` / `cell_start` from the sorted keys.
    ///
    /// [`build`](Self::build) calls this after the hash kernel; it is separate so
    /// that tests can feed synthetic keys through the identical sort/bucket code.
    pub(super) fn index_keys(&self, backend: &WgpuBackendReal) -> Result<(), WgpuInitError> {
        self.sort.run(backend)?;
        backend.queue_write_buffer_raw(&self.cell_count, &self.zero_counts);
        backend.dispatch_wgsl(
            SPH_CELL_HISTOGRAM_WGSL,
            "sph_cell_histogram",
            &[
                (self.sort.sorted_keys(), ro()),
                (self.cell_count, rw()),
                (self.bin_params, ro()),
            ],
            [self.n.div_ceil(BLOCK) as u32, 1, 1],
        )?;
        self.cell_scan.run(backend)
    }

    /// Download the built cell list. Returns `None` if a readback came back short.
    pub(super) fn read_back(&self, backend: &WgpuBackendReal) -> Option<CellListHost> {
        let mut sorted_ids = backend.read_buffer_u32(self.sorted_ids()).ok()?;
        let mut cell_start = backend.read_buffer_u32(self.cell_start).ok()?;
        let mut cell_count = backend.read_buffer_u32(self.cell_count).ok()?;
        if sorted_ids.len() < self.n
            || cell_start.len() < self.num_cells
            || cell_count.len() < self.num_cells
        {
            return None;
        }
        sorted_ids.truncate(self.n);
        cell_start.truncate(self.num_cells);
        cell_count.truncate(self.num_cells);
        Some(CellListHost {
            sorted_ids,
            cell_start,
            cell_count,
        })
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic LCG, one 32-bit word per call.
    fn lcg(state: &mut u32) -> u32 {
        *state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        *state
    }

    fn try_backend() -> Option<WgpuBackendReal> {
        match WgpuBackendReal::try_new() {
            Ok(backend) => Some(backend),
            Err(_) => {
                eprintln!("SKIPPED: no GPU adapter available");
                None
            }
        }
    }

    fn cpu_exclusive_scan(data: &[u32]) -> Vec<u32> {
        let mut acc = 0u32;
        data.iter()
            .map(|&v| {
                let before = acc;
                acc = acc.wrapping_add(v);
                before
            })
            .collect()
    }

    #[test]
    fn scan_capacity_rounds_up_to_whole_blocks() {
        assert_eq!(scan_capacity_words(0), 256);
        assert_eq!(scan_capacity_words(1), 256);
        assert_eq!(scan_capacity_words(256), 256);
        assert_eq!(scan_capacity_words(257), 512);
        assert_eq!(scan_capacity_words(65_536), 65_536);
    }

    #[test]
    fn radix_pass_count_follows_key_range() {
        assert_eq!(radix_passes_for(1), 1);
        assert_eq!(radix_passes_for(2), 1);
        assert_eq!(radix_passes_for(256), 1);
        assert_eq!(radix_passes_for(257), 2);
        assert_eq!(radix_passes_for(65_536), 2);
        assert_eq!(radix_passes_for(65_537), 3);
        assert_eq!(radix_passes_for(1 << 24), 3);
        assert_eq!(radix_passes_for((1 << 24) + 1), 4);
        assert_eq!(radix_passes_for(u32::MAX as usize + 1), 4);
    }

    #[test]
    fn wgsl_kernels_pass_naga_validation() {
        for (name, src) in [
            ("SPH_CELL_HISTOGRAM_WGSL", SPH_CELL_HISTOGRAM_WGSL),
            ("SPH_CELL_LIST_WGSL", SPH_CELL_LIST_WGSL),
            ("SCAN_WGSL", SCAN_WGSL),
            ("SCAN_ADD_WGSL", SCAN_ADD_WGSL),
            ("RADIX_HISTOGRAM_WGSL", RADIX_HISTOGRAM_WGSL),
            ("RADIX_SCATTER_PAIRS_WGSL", RADIX_SCATTER_PAIRS_WGSL),
        ] {
            let module = naga::front::wgsl::parse_str(src)
                .unwrap_or_else(|e| panic!("{name}: parse failed: {}", e.emit_to_string(src)));
            let mut validator = naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::all(),
            );
            validator
                .validate(&module)
                .unwrap_or_else(|e| panic!("{name}: validation failed: {e:?}"));
        }
    }

    /// The resident scan must equal the CPU scan for lengths that hit one block,
    /// a block boundary, several blocks, and the three-level hierarchy (> 65536).
    #[test]
    fn resident_scan_matches_cpu_for_all_hierarchy_depths() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        let mut rng = 0x00C0_FFEE_u32;
        for &(len, expected_levels) in &[
            (1usize, 1usize),
            (2, 1),
            (255, 1),
            (256, 1),
            (257, 2),
            (1000, 2),
            (65_536, 2),
            (65_537, 3),
            (70_000, 3),
        ] {
            let capacity = scan_capacity_words(len);
            let input = backend.create_buffer_storage(words_to_bytes(capacity));
            let output = backend.create_buffer_storage(words_to_bytes(capacity));
            let data: Vec<u32> = (0..len).map(|_| lcg(&mut rng) >> 16).collect();
            backend.queue_write_buffer_raw(&input, bytemuck::cast_slice(&data));

            let scan = ResidentScan::new(&mut backend, input, output, len).expect("scan build");
            assert_eq!(scan.level_count(), expected_levels, "levels for len {len}");
            scan.run(&backend).expect("scan run");

            let mut gpu = backend.read_buffer_u32(output).expect("scan readback");
            gpu.truncate(len);
            assert_eq!(gpu, cpu_exclusive_scan(&data), "scan mismatch at len {len}");
        }
    }

    #[test]
    fn resident_scan_rejects_undersized_buffers() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        let small = backend.create_buffer_storage(64);
        let big = backend.create_buffer_storage(words_to_bytes(512));
        assert!(matches!(
            ResidentScan::new(&mut backend, small, big, 300),
            Err(WgpuInitError::InvalidHandle(_))
        ));
        assert!(matches!(
            ResidentScan::new(&mut backend, big, small, 300),
            Err(WgpuInitError::InvalidHandle(_))
        ));
    }

    /// Stable CPU reference: payload = original index, so equal keys must keep
    /// their input order.
    fn cpu_stable_sort(keys: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let mut idx: Vec<u32> = (0..keys.len() as u32).collect();
        idx.sort_by_key(|&i| keys[i as usize]);
        let sorted_keys = idx.iter().map(|&i| keys[i as usize]).collect();
        (sorted_keys, idx)
    }

    /// Sort `n` random keys below `key_limit` on the GPU and compare keys and
    /// payload with the stable CPU sort. Covers single-tile, partial-tile,
    /// multi-tile and (n = 70000) multi-level-scan inputs and 1 to 4 passes.
    #[test]
    fn resident_radix_sort_is_stable_and_matches_cpu() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        let mut rng = 0x1234_5678_u32;
        for &(n, key_limit) in &[
            (1usize, 1usize),
            (2, 2),
            (255, 300),
            (256, 1 << 16),
            (257, 1 << 16),
            (3000, 1 << 12),
            (3000, 1 << 24),
            (5000, u32::MAX as usize + 1),
            (70_000, 1 << 16),
        ] {
            let sorter = ResidentRadixPairs::new(&mut backend, n, key_limit).expect("sorter");
            let keys: Vec<u32> = (0..n)
                .map(|_| (u64::from(lcg(&mut rng)) % key_limit as u64) as u32)
                .collect();
            let ids: Vec<u32> = (0..n as u32).collect();
            backend.queue_write_buffer_raw(&sorter.input_keys(), bytemuck::cast_slice(&keys));
            backend.queue_write_buffer_raw(&sorter.input_payload(), bytemuck::cast_slice(&ids));
            sorter.run(&backend).expect("sort run");

            let mut gpu_keys = backend
                .read_buffer_u32(sorter.sorted_keys())
                .expect("sorted key readback");
            let mut gpu_ids = backend
                .read_buffer_u32(sorter.sorted_payload())
                .expect("sorted payload readback");
            gpu_keys.truncate(n);
            gpu_ids.truncate(n);
            let (cpu_keys, cpu_ids) = cpu_stable_sort(&keys);
            assert_eq!(gpu_keys, cpu_keys, "keys differ (n={n}, limit={key_limit})");
            assert_eq!(
                gpu_ids, cpu_ids,
                "payload differs (n={n}, limit={key_limit})"
            );
        }
    }

    /// Bucket bounds produced from synthetic keys by the real index code must
    /// equal the CPU counting-sort layout, including > 256 buckets (which the
    /// 256-bin workgroup histogram cannot handle).
    #[test]
    fn index_keys_matches_cpu_bucket_layout_for_many_buckets() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        let mut rng = 0xDEAD_BEEF_u32;
        for &(n, num_cells) in &[(1usize, 1usize), (64, 64), (3000, 4096), (6000, 65_536)] {
            let cells = ResidentCellList::new(&mut backend, n, num_cells).expect("cell list");
            let keys: Vec<u32> = (0..n).map(|_| lcg(&mut rng) % num_cells as u32).collect();
            let ids: Vec<u32> = (0..n as u32).collect();
            backend.queue_write_buffer_raw(&cells.sort.input_keys(), bytemuck::cast_slice(&keys));
            backend.queue_write_buffer_raw(&cells.sort.input_payload(), bytemuck::cast_slice(&ids));
            cells.index_keys(&backend).expect("index keys");
            let host = cells.read_back(&backend).expect("read back");

            let mut expected_count = vec![0u32; num_cells];
            for &k in &keys {
                expected_count[k as usize] += 1;
            }
            let expected_start = cpu_exclusive_scan(&expected_count);
            let (_, expected_ids) = cpu_stable_sort(&keys);
            assert_eq!(
                host.cell_count, expected_count,
                "counts (n={n}, cells={num_cells})"
            );
            assert_eq!(
                host.cell_start, expected_start,
                "starts (n={n}, cells={num_cells})"
            );
            assert_eq!(
                host.sorted_ids, expected_ids,
                "ids (n={n}, cells={num_cells})"
            );
            assert_eq!(cells.num_cells(), num_cells);
        }
    }

    /// Building and running must never allocate: the live buffer count after the
    /// first and after the hundredth build are identical.
    #[test]
    fn repeated_builds_keep_the_buffer_table_constant() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        let n = 500usize;
        let cells = ResidentCellList::new(&mut backend, n, 512).expect("cell list");
        let positions = backend.create_buffer_storage((n * 3 * 4) as u64);
        let mut rng = 7u32;
        let host: Vec<f32> = (0..n * 3)
            .map(|_| (lcg(&mut rng) >> 8) as f32 / (1u32 << 24) as f32)
            .collect();
        backend.queue_write_buffer_f32(&positions, &host);

        cells
            .build(&backend, positions, 0.1, [0.0; 3])
            .expect("first build");
        let count_after_first = backend.live_buffer_count();
        let bytes_after_first = backend.live_buffer_bytes();
        for _ in 0..100 {
            cells
                .build(&backend, positions, 0.1, [0.0; 3])
                .expect("repeated build");
        }
        assert_eq!(backend.live_buffer_count(), count_after_first);
        assert_eq!(backend.live_buffer_bytes(), bytes_after_first);
    }

    #[test]
    fn constructors_reject_empty_problems() {
        let Some(mut backend) = try_backend() else {
            return;
        };
        assert!(ResidentRadixPairs::new(&mut backend, 0, 16).is_err());
        assert!(ResidentCellList::new(&mut backend, 0, 16).is_err());
        assert!(ResidentCellList::new(&mut backend, 16, 0).is_err());
    }
}
