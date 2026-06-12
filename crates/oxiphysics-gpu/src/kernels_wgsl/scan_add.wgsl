// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Uniform add pass: output[i] += block_offsets[workgroup_id] for the
// cross-block combination of the per-block exclusive scans.

struct Params {
    n: u32,
    num_bins: u32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read_write> output: array<u32>;
@group(0) @binding(1) var<storage, read>       block_offsets: array<u32>;
@group(0) @binding(2) var<storage, read>       params: Params;

@compute @workgroup_size(256)
fn add_block_offsets(
    @builtin(global_invocation_id) gid_vec: vec3<u32>,
    @builtin(workgroup_id) wg_vec: vec3<u32>,
) {
    let gid = gid_vec.x;
    let wg_id = wg_vec.x;
    let n = params.n;
    let offset_val = block_offsets[wg_id];
    if (gid < n) {
        output[gid] = output[gid] + offset_val;
    }
}
