// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Histogram over u32: each value v contributes to bin min(v, num_bins-1).
// Per-workgroup atomic local bins are merged into global atomic bins.
// Constraint: num_bins <= 256 (host falls back to CPU otherwise).

struct Params {
    n: u32,
    num_bins: u32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       data: array<u32>;
@group(0) @binding(1) var<storage, read_write> global_bins: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read>       params: Params;

var<workgroup> local_bins: array<atomic<u32>, 256u>;

@compute @workgroup_size(256)
fn histogram_main(
    @builtin(global_invocation_id) gid_vec: vec3<u32>,
    @builtin(local_invocation_id) lid_vec: vec3<u32>,
) {
    let gid = gid_vec.x;
    let lid = lid_vec.x;
    let n = params.n;
    let num_bins = params.num_bins;

    var b: u32 = lid;
    loop {
        if (b >= 256u) { break; }
        atomicStore(&local_bins[b], 0u);
        b = b + 256u;
    }
    workgroupBarrier();

    if (gid < n) {
        var bin_idx = data[gid];
        let last = num_bins - 1u;
        if (bin_idx > last) {
            bin_idx = last;
        }
        atomicAdd(&local_bins[bin_idx], 1u);
    }
    workgroupBarrier();

    var m: u32 = lid;
    loop {
        if (m >= num_bins) { break; }
        let v = atomicLoad(&local_bins[m]);
        if (v > 0u) {
            atomicAdd(&global_bins[m], v);
        }
        m = m + 256u;
    }
}
