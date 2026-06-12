// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH density summation over a spatial-hash cell-list.
//
// For each particle i, visits the 27 adjacent grid cells, hashes each into the
// same table the cell-list kernel used, walks the sorted particle slice in that
// bucket (cell_start[h] .. cell_start[h] + cell_count[h]), and accumulates the
// cubic-spline kernel contribution W(|r_ij|, h). Hash collisions are tolerated:
// every candidate is distance-filtered (r2 < support2) before contributing, so a
// foreign cell sharing a bucket adds nothing. To avoid double-visiting a bucket
// when two of the 27 neighbor cells collide to the same hash, a small visited
// list of already-processed hashes is kept (27 entries max).
//
// `params` and the cell-list arrays are read-only storage buffers (matching the
// all-storage dispatch bind-group layout).

struct SphGridParams {
    n: u32,           // particle count
    num_cells: u32,   // hash table size
    h: f32,           // smoothing length
    support: f32,     // kernel support radius (= 2h for cubic spline)
    mass: f32,        // particle mass
    cell_size: f32,   // grid spacing
    origin_x: f32,
    origin_y: f32,
    origin_z: f32,
    pad0: u32,
    pad1: u32,
    pad2: u32,
}

@group(0) @binding(0) var<storage, read>       positions:  array<f32>; // x,y,z interleaved
@group(0) @binding(1) var<storage, read>       sorted_ids: array<u32>; // particle list sorted by cell key
@group(0) @binding(2) var<storage, read>       cell_start: array<u32>; // exclusive-scan offset per hash bucket
@group(0) @binding(3) var<storage, read>       cell_count: array<u32>; // particle count per hash bucket
@group(0) @binding(4) var<storage, read_write> densities:  array<f32>; // output density per particle
@group(0) @binding(5) var<storage, read>       params:     SphGridParams;

fn cell_hash(ix: i32, iy: i32, iz: i32, num_cells: u32) -> u32 {
    let p1: u32 = 73856093u;
    let p2: u32 = 19349663u;
    let p3: u32 = 83492791u;
    let hx = bitcast<u32>(ix) * p1;
    let hy = bitcast<u32>(iy) * p2;
    let hz = bitcast<u32>(iz) * p3;
    return (hx ^ hy ^ hz) % num_cells;
}

// Cubic-spline kernel W3(r, h), normalized for 3-D with support 2h.
// sigma = 1/(pi h^3); piecewise on q = r/h over [0,1), [1,2).
fn w_spline3(r: f32, h: f32) -> f32 {
    let q = r / h;
    let sigma = 1.0 / (3.14159265358979 * h * h * h);
    if (q < 1.0) {
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    } else if (q < 2.0) {
        let t = 2.0 - q;
        return sigma * 0.25 * t * t * t;
    }
    return 0.0;
}

@compute @workgroup_size(64)
fn sph_density(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }

    let pix = positions[i * 3u];
    let piy = positions[i * 3u + 1u];
    let piz = positions[i * 3u + 2u];

    let inv_cs = 1.0 / params.cell_size;
    let cx = i32(floor((pix - params.origin_x) * inv_cs));
    let cy = i32(floor((piy - params.origin_y) * inv_cs));
    let cz = i32(floor((piz - params.origin_z) * inv_cs));

    let support2 = params.support * params.support;
    var density: f32 = 0.0;

    // Track up to 27 already-visited bucket hashes to avoid double counting
    // when distinct neighbor cells alias to the same hash bucket.
    var visited: array<u32, 27u>;
    var n_visited: u32 = 0u;

    for (var dz: i32 = -1; dz <= 1; dz = dz + 1) {
        for (var dy: i32 = -1; dy <= 1; dy = dy + 1) {
            for (var dx: i32 = -1; dx <= 1; dx = dx + 1) {
                let h_bucket = cell_hash(cx + dx, cy + dy, cz + dz, params.num_cells);

                // Skip if this bucket hash was already processed.
                var already: bool = false;
                for (var v: u32 = 0u; v < n_visited; v = v + 1u) {
                    if (visited[v] == h_bucket) {
                        already = true;
                        break;
                    }
                }
                if (already) { continue; }
                visited[n_visited] = h_bucket;
                n_visited = n_visited + 1u;

                let start = cell_start[h_bucket];
                let count = cell_count[h_bucket];
                for (var s: u32 = 0u; s < count; s = s + 1u) {
                    let j = sorted_ids[start + s];
                    let rx = pix - positions[j * 3u];
                    let ry = piy - positions[j * 3u + 1u];
                    let rz = piz - positions[j * 3u + 2u];
                    let r2 = rx * rx + ry * ry + rz * rz;
                    if (r2 < support2) {
                        density = density + params.mass * w_spline3(sqrt(r2), params.h);
                    }
                }
            }
        }
    }

    densities[i] = max(density, 1e-6);
}
