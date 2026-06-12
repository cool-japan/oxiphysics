// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Tree reductions over f32: per-block sum and per-block max.
// The input buffer is host-padded to a multiple of 256 with the neutral
// element (0 for sum, -FLT_MAX for max) so all 256 lanes load in-bounds.

struct Params {
    n: u32,
    num_bins: u32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       input: array<f32>;
@group(0) @binding(1) var<storage, read_write> partials: array<f32>;
@group(0) @binding(2) var<storage, read>       params: Params;

var<workgroup> scratch: array<f32, 256u>;

@compute @workgroup_size(256)
fn reduce_sum_block(
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
    @builtin(global_invocation_id) gid_vec: vec3<u32>,
    @builtin(workgroup_id) wg_vec: vec3<u32>,
) {
    let lid = lid_vec.x;
    let gid = gid_vec.x;
    let wg_id = wg_vec.x;
    scratch[lid] = select(0.0, input[gid], gid < params.n);
    workgroupBarrier();

    var stride: u32 = 256u >> 1u;
    loop {
        if (stride == 0u) { break; }
        if (lid < stride) {
            scratch[lid] = scratch[lid] + scratch[lid + stride];
        }
        workgroupBarrier();
        stride = stride >> 1u;
    }

    if (lid == 0u) {
        partials[wg_id] = scratch[0u];
    }
}

@compute @workgroup_size(256)
fn reduce_max_block(
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
    @builtin(global_invocation_id) gid_vec: vec3<u32>,
    @builtin(workgroup_id) wg_vec: vec3<u32>,
) {
    let lid = lid_vec.x;
    let gid = gid_vec.x;
    let wg_id = wg_vec.x;
    scratch[lid] = select(-3.4028235e38, input[gid], gid < params.n);
    workgroupBarrier();

    var stride: u32 = 256u >> 1u;
    loop {
        if (stride == 0u) { break; }
        if (lid < stride) {
            scratch[lid] = max(scratch[lid], scratch[lid + stride]);
        }
        workgroupBarrier();
        stride = stride >> 1u;
    }

    if (lid == 0u) {
        partials[wg_id] = scratch[0u];
    }
}
