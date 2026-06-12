// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Radix-sort digit histogram (one workgroup per 256-element tile).
// Builds a digit-major histogram: tile_hist[d * num_tiles + tile_id] = count
// of elements in tile `tile_id` whose 8-bit digit == d. Out-of-range lanes
// (gid >= n) contribute nothing. Workgroup_size = 256 == number of digits.

struct RadixParams {
    n: u32,
    shift: u32,
    num_tiles: u32,
    pad: u32,
}

@group(0) @binding(0) var<storage, read>       keys: array<u32>;
@group(0) @binding(1) var<storage, read_write> tile_hist: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read>       params: RadixParams;

var<workgroup> local_bins: array<atomic<u32>, 256u>;

@compute @workgroup_size(256)
fn radix_histogram(
    @builtin(workgroup_id) wg_vec: vec3<u32>,
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
) {
    let lid = lid_vec.x;
    let tile_id = wg_vec.x;

    atomicStore(&local_bins[lid], 0u);
    workgroupBarrier();

    let gid = tile_id * 256u + lid;
    if (gid < params.n) {
        let digit = (keys[gid] >> params.shift) & 0xFFu;
        atomicAdd(&local_bins[digit], 1u);
    }
    workgroupBarrier();

    tile_hist[lid * params.num_tiles + tile_id] = atomicLoad(&local_bins[lid]);
}
