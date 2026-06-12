// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH spatial-hash cell-list key generation.
//
// For each particle, computes the integer grid cell containing it
// (floor((pos - origin) / cell_size) per axis), folds the 3-D cell
// coordinate into a 1-D spatial-hash key in [0, num_cells), and writes
// (key, particle_index) to parallel output buffers. The key buffer feeds
// `radix_sort_pairs_gpu`; the payload buffer is the particle index that
// follows the key through the sort to produce the sorted particle list.
//
// `params` is bound as a read-only storage buffer (not a uniform) to match
// the dispatch path's all-storage bind-group layout (`dispatch_wgsl`).

struct CellParams {
    n: u32,           // particle count
    num_cells: u32,   // hash table size (= nx*ny*nz, clamped > 0)
    cell_size: f32,   // grid spacing (= kernel support radius)
    origin_x: f32,    // domain-min corner
    origin_y: f32,
    origin_z: f32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read>       positions: array<f32>;   // x0,y0,z0, x1,y1,z1, ...
@group(0) @binding(1) var<storage, read_write> cell_keys: array<u32>;    // spatial hash per particle
@group(0) @binding(2) var<storage, read_write> particle_ids: array<u32>; // payload = particle index
@group(0) @binding(3) var<storage, read>       params: CellParams;

// 3-D integer cell coordinate → 1-D hash in [0, num_cells).
// Large odd primes scramble each axis; XOR mixes them; modulo folds into range.
fn cell_hash(ix: i32, iy: i32, iz: i32, num_cells: u32) -> u32 {
    let p1: u32 = 73856093u;
    let p2: u32 = 19349663u;
    let p3: u32 = 83492791u;
    let hx = bitcast<u32>(ix) * p1;
    let hy = bitcast<u32>(iy) * p2;
    let hz = bitcast<u32>(iz) * p3;
    return (hx ^ hy ^ hz) % num_cells;
}

@compute @workgroup_size(64)
fn sph_cell_list(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }

    let px = positions[i * 3u];
    let py = positions[i * 3u + 1u];
    let pz = positions[i * 3u + 2u];

    let inv_cs = 1.0 / params.cell_size;
    let ix = i32(floor((px - params.origin_x) * inv_cs));
    let iy = i32(floor((py - params.origin_y) * inv_cs));
    let iz = i32(floor((pz - params.origin_z) * inv_cs));

    cell_keys[i] = cell_hash(ix, iy, iz, params.num_cells);
    particle_ids[i] = i;
}
