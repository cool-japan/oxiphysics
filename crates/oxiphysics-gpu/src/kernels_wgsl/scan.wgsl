// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Exclusive prefix scan (Blelloch up-sweep / down-sweep) over u32.
// Fixed 256-wide intra-block scan; host handles cross-block combination.

struct Params {
    n: u32,
    num_bins: u32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       input: array<u32>;
@group(0) @binding(1) var<storage, read_write> output: array<u32>;
@group(0) @binding(2) var<storage, read_write> block_sums: array<u32>;
@group(0) @binding(3) var<storage, read>       params: Params;

var<workgroup> shared_data: array<u32, 256u>;

@compute @workgroup_size(256)
fn scan_block(
    @builtin(global_invocation_id) gid_vec: vec3<u32>,
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
    @builtin(workgroup_id) wg_vec: vec3<u32>,
) {
    let lid = lid_vec.x;
    let gid = gid_vec.x;
    let wg_id = wg_vec.x;
    let n = params.n;

    shared_data[lid] = input[gid];
    workgroupBarrier();

    var offset_val: u32 = 1u;
    var d: u32 = 256u >> 1u;
    loop {
        if (d == 0u) { break; }
        if (lid < d) {
            let ai = offset_val * (2u * lid + 1u) - 1u;
            let bi = offset_val * (2u * lid + 2u) - 1u;
            shared_data[bi] = shared_data[bi] + shared_data[ai];
        }
        offset_val = offset_val << 1u;
        workgroupBarrier();
        d = d >> 1u;
    }

    if (lid == 0u) {
        block_sums[wg_id] = shared_data[255u];
        shared_data[255u] = 0u;
    }
    workgroupBarrier();

    d = 1u;
    loop {
        if (d >= 256u) { break; }
        offset_val = offset_val >> 1u;
        if (lid < d) {
            let ai = offset_val * (2u * lid + 1u) - 1u;
            let bi = offset_val * (2u * lid + 2u) - 1u;
            let t = shared_data[ai];
            shared_data[ai] = shared_data[bi];
            shared_data[bi] = shared_data[bi] + t;
        }
        workgroupBarrier();
        d = d << 1u;
    }

    if (gid < n) {
        output[gid] = shared_data[lid];
    }
}
