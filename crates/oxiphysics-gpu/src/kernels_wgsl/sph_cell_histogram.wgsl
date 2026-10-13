// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH cell-list bucket histogram: cell_count[k] = number of particles whose
// spatial-hash key equals k.
//
// The hash table has up to 65536 buckets, far beyond the 256 bins the
// workgroup-local `histogram.wgsl` can hold, so every particle adds directly
// into a global atomic bin. The host zeroes `bins` before the dispatch. Keys
// that fall outside the table are clamped into the last bin (defensive only:
// the cell-list kernel always emits keys < num_cells).
//
// `params` is a read-only storage buffer to match the all-storage bind-group
// layout used by `WgpuBackendReal::dispatch_wgsl`.

struct BinParams {
    n: u32,         // number of keys
    num_bins: u32,  // number of buckets (hash-table size)
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       keys: array<u32>;
@group(0) @binding(1) var<storage, read_write> bins: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read>       params: BinParams;

@compute @workgroup_size(256)
fn sph_cell_histogram(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    // A dispatch wider than one dimension allows is folded into y (see
    // `dispatch_count_for`): linearise the invocation index.
    let i = gid.x + gid.y * nwg.x * 256u;
    if (i >= params.n) { return; }
    let k = min(keys[i], params.num_bins - 1u);
    atomicAdd(&bins[k], 1u);
}
