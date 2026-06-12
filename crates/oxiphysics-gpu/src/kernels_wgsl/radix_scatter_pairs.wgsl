// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Radix-sort stable scatter (key + u32 payload). One workgroup per 256-element
// tile. Identical stable-rank logic to radix_scatter.wgsl, but also moves the
// payload alongside the key so (key,payload) pairs stay aligned and stable.
// Out-of-range lanes (gid >= n) do not write. Workgroup_size = 256.

struct RadixParams {
    n: u32,
    shift: u32,
    num_tiles: u32,
    pad: u32,
}

@group(0) @binding(0) var<storage, read>       keys_in: array<u32>;
@group(0) @binding(1) var<storage, read>       payload_in: array<u32>;
@group(0) @binding(2) var<storage, read_write> keys_out: array<u32>;
@group(0) @binding(3) var<storage, read_write> payload_out: array<u32>;
@group(0) @binding(4) var<storage, read>       base: array<u32>;
@group(0) @binding(5) var<storage, read>       params: RadixParams;

var<workgroup> shared_digit: array<u32, 256u>;

@compute @workgroup_size(256)
fn radix_scatter_pairs(
    @builtin(workgroup_id) wg_vec: vec3<u32>,
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
) {
    let lid = lid_vec.x;
    let tile_id = wg_vec.x;
    let gid = tile_id * 256u + lid;

    var digit_self: u32;
    if (gid < params.n) {
        digit_self = (keys_in[gid] >> params.shift) & 0xFFu;
    } else {
        digit_self = 0xFFFFFFFFu;
    }
    shared_digit[lid] = digit_self;
    workgroupBarrier();

    if (gid < params.n) {
        let d = shared_digit[lid];
        var local_rank: u32 = 0u;
        var j: u32 = 0u;
        loop {
            if (j >= lid) { break; }
            if (shared_digit[j] == d) { local_rank = local_rank + 1u; }
            j = j + 1u;
        }
        let dst = base[d * params.num_tiles + tile_id] + local_rank;
        keys_out[dst] = keys_in[gid];
        payload_out[dst] = payload_in[gid];
    }
}
