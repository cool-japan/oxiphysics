// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH pressure + viscosity + gravity acceleration over a spatial-hash cell-list.
//
// For each particle i: Tait EOS pressure p = k (rho/rho0 - 1), then over the
// neighbors found in the 27 adjacent hash buckets accumulates the symmetric
// Monaghan pressure acceleration  -m (p_i/rho_i^2 + p_j/rho_j^2) grad W  plus
// Monaghan artificial viscosity (active only on approaching pairs), and finally
// adds gravity (-g in Y). The result is the per-particle acceleration written
// to `accel` (xyz interleaved). This is bit-for-bit the same physical model the
// CPU reference (`sph_gpu::pressure_and_integrate`) uses, so f32-vs-f64 is the
// only divergence the parity test sees.
//
// All buffers including `params` are storage buffers (all-storage dispatch).

struct SphForceParams {
    n: u32,           // particle count
    num_cells: u32,   // hash table size
    h: f32,           // smoothing length
    support: f32,     // kernel support radius (= 2h)
    mass: f32,        // particle mass
    cell_size: f32,   // grid spacing
    rest_density: f32,// rho0
    pressure_k: f32,  // EOS stiffness k
    viscosity: f32,   // kinematic viscosity nu
    gravity: f32,     // gravitational acceleration (applied -Y)
    origin_x: f32,
    origin_y: f32,
    origin_z: f32,
    pad0: u32,
    pad1: u32,
    pad2: u32,
}

@group(0) @binding(0) var<storage, read>       positions:  array<f32>;
@group(0) @binding(1) var<storage, read>       velocities: array<f32>;
@group(0) @binding(2) var<storage, read>       densities:  array<f32>;
@group(0) @binding(3) var<storage, read>       sorted_ids: array<u32>;
@group(0) @binding(4) var<storage, read>       cell_start: array<u32>;
@group(0) @binding(5) var<storage, read>       cell_count: array<u32>;
@group(0) @binding(6) var<storage, read_write> accel:      array<f32>;
@group(0) @binding(7) var<storage, read>       params:     SphForceParams;

fn cell_hash(ix: i32, iy: i32, iz: i32, num_cells: u32) -> u32 {
    let p1: u32 = 73856093u;
    let p2: u32 = 19349663u;
    let p3: u32 = 83492791u;
    let hx = bitcast<u32>(ix) * p1;
    let hy = bitcast<u32>(iy) * p2;
    let hz = bitcast<u32>(iz) * p3;
    return (hx ^ hy ^ hz) % num_cells;
}

// |grad W| / r factor for cubic spline: returns dW/dr divided by r so the
// caller multiplies by the separation vector. Matches CPU cubic_spline_dw_dr.
// dW/dr piecewise: sigma4 = 1/(pi h^4); over q in [0,1): (-3q + 2.25 q^2),
// over [1,2): -0.75 (2-q)^2.
fn dw_dr(r: f32, h: f32) -> f32 {
    let q = r / h;
    let sigma4 = 1.0 / (3.14159265358979 * h * h * h * h);
    if (q < 1.0) {
        return sigma4 * (-3.0 * q + 2.25 * q * q);
    } else if (q < 2.0) {
        let t = 2.0 - q;
        return sigma4 * (-0.75 * t * t);
    }
    return 0.0;
}

@compute @workgroup_size(64)
fn sph_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }

    let pix = positions[i * 3u];
    let piy = positions[i * 3u + 1u];
    let piz = positions[i * 3u + 2u];
    let vix = velocities[i * 3u];
    let viy = velocities[i * 3u + 1u];
    let viz = velocities[i * 3u + 2u];

    let rho_i = densities[i];
    let p_i = params.pressure_k * (rho_i / params.rest_density - 1.0);

    let inv_cs = 1.0 / params.cell_size;
    let cx = i32(floor((pix - params.origin_x) * inv_cs));
    let cy = i32(floor((piy - params.origin_y) * inv_cs));
    let cz = i32(floor((piz - params.origin_z) * inv_cs));

    let support2 = params.support * params.support;
    let h = params.h;
    let m = params.mass;
    let nu = params.viscosity;

    var ax: f32 = 0.0;
    var ay: f32 = -params.gravity;
    var az: f32 = 0.0;

    var visited: array<u32, 27u>;
    var n_visited: u32 = 0u;

    for (var dz: i32 = -1; dz <= 1; dz = dz + 1) {
        for (var dy: i32 = -1; dy <= 1; dy = dy + 1) {
            for (var dx: i32 = -1; dx <= 1; dx = dx + 1) {
                let h_bucket = cell_hash(cx + dx, cy + dy, cz + dz, params.num_cells);

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
                    if (j == i) { continue; }

                    let rx = pix - positions[j * 3u];
                    let ry = piy - positions[j * 3u + 1u];
                    let rz = piz - positions[j * 3u + 2u];
                    let r2 = rx * rx + ry * ry + rz * rz;
                    if (r2 >= support2 || r2 <= 1e-12) { continue; }

                    let r = sqrt(r2);
                    let rho_j = densities[j];
                    let p_j = params.pressure_k * (rho_j / params.rest_density - 1.0);

                    // Symmetric Monaghan pressure acceleration.
                    let dw = dw_dr(r, h);
                    let pf = -m * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw;
                    let inv_r = 1.0 / r;
                    ax = ax + pf * rx * inv_r;
                    ay = ay + pf * ry * inv_r;
                    az = az + pf * rz * inv_r;

                    // Artificial viscosity: active only when pair approaches.
                    let vjx = velocities[j * 3u];
                    let vjy = velocities[j * 3u + 1u];
                    let vjz = velocities[j * 3u + 2u];
                    let dvx = vix - vjx;
                    let dvy = viy - vjy;
                    let dvz = viz - vjz;
                    let vdotr = dvx * rx + dvy * ry + dvz * rz;
                    if (vdotr < 0.0) {
                        let vf = nu * m / rho_j * vdotr / (r2 + 0.01 * h * h) * dw * inv_r;
                        ax = ax + vf * rx;
                        ay = ay + vf * ry;
                        az = az + vf * rz;
                    }
                }
            }
        }
    }

    accel[i * 3u] = ax;
    accel[i * 3u + 1u] = ay;
    accel[i * 3u + 2u] = az;
}
