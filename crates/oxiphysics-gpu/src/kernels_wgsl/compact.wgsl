// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Stream-compaction scatter: writes kept values to their exclusive-scan
// positions, preserving input order (stable compaction).

struct Params {
    n: u32,
    num_bins: u32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       values: array<u32>;
@group(0) @binding(1) var<storage, read>       keep: array<u32>;
@group(0) @binding(2) var<storage, read>       positions: array<u32>;
@group(0) @binding(3) var<storage, read_write> output: array<u32>;
@group(0) @binding(4) var<storage, read>       params: Params;

@compute @workgroup_size(256)
fn scatter_kept(@builtin(global_invocation_id) gid_vec: vec3<u32>) {
    let gid = gid_vec.x;
    if (gid >= params.n) {
        return;
    }
    if (keep[gid] == 1u) {
        output[positions[gid]] = values[gid];
    }
}
