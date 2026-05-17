// D3Q19 BGK streaming + collision
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Dispatch: (ceil(nx/8), ceil(ny/8), ceil(nz/8))
// Each workgroup processes an 8×8×8 tile of lattice sites.
//
// Buffer layout:
//   binding(0) — params (storage read, 5 x u32/f32: nx, ny, nz, omega, _pad)
//   binding(1) — f_in   (storage read,  19 * nx*ny*nz f32 values)
//   binding(2) — f_out  (storage r/w,   19 * nx*ny*nz f32 values)

struct LbmParams {
    nx:    u32,
    ny:    u32,
    nz:    u32,
    omega: f32,
    _pad:  u32,
};

@group(0) @binding(0) var<storage, read>       params: array<u32>;   // 5 words: nx, ny, nz, omega_bits, _pad
@group(0) @binding(1) var<storage, read>       f_in:   array<f32>;   // 19 * nx*ny*nz
@group(0) @binding(2) var<storage, read_write> f_out:  array<f32>;

// D3Q19 weights (compile-time constants)
const W: array<f32, 19> = array<f32, 19>(
    1.0 / 3.0,
    1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
);

// D3Q19 velocity directions (ex, ey, ez) x 19
const CX: array<i32, 19> = array<i32, 19>(0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0);
const CY: array<i32, 19> = array<i32, 19>(0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1);
const CZ: array<i32, 19> = array<i32, 19>(0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1);

// Compute flat index into the 19*N population array (q-major layout)
fn idx(x: u32, y: u32, z: u32, q: u32, nx: u32, ny: u32, nz: u32) -> u32 {
    return q * (nx * ny * nz) + z * (nx * ny) + y * nx + x;
}

@compute @workgroup_size(8, 8, 8)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    // Read params from the storage buffer (bit-cast f32 fields)
    let nx    = params[0];
    let ny    = params[1];
    let nz    = params[2];
    let omega = bitcast<f32>(params[3]);

    let x = gid.x;
    let y = gid.y;
    let z = gid.z;

    // Guard: discard threads outside the domain
    if x >= nx || y >= ny || z >= nz {
        return;
    }

    // --- Streaming (pull) -----------------------------------------------
    // Each thread collects its post-streaming populations by pulling from
    // the upstream neighbour in each direction.  Wrap-around gives periodic
    // boundary conditions on all faces.
    var f: array<f32, 19>;
    for (var q: u32 = 0u; q < 19u; q++) {
        let sx = u32((i32(x) - CX[q] + i32(nx)) % i32(nx));
        let sy = u32((i32(y) - CY[q] + i32(ny)) % i32(ny));
        let sz = u32((i32(z) - CZ[q] + i32(nz)) % i32(nz));
        f[q] = f_in[idx(sx, sy, sz, q, nx, ny, nz)];
    }

    // --- Macroscopic quantities ------------------------------------------
    var rho: f32 = 0.0;
    var ux:  f32 = 0.0;
    var uy:  f32 = 0.0;
    var uz:  f32 = 0.0;
    for (var q: u32 = 0u; q < 19u; q++) {
        rho += f[q];
        ux  += f32(CX[q]) * f[q];
        uy  += f32(CY[q]) * f[q];
        uz  += f32(CZ[q]) * f[q];
    }
    // Avoid division by zero (degenerate / empty cells get no update)
    if rho < 1e-10 {
        for (var q: u32 = 0u; q < 19u; q++) {
            f_out[idx(x, y, z, q, nx, ny, nz)] = f[q];
        }
        return;
    }
    ux /= rho;
    uy /= rho;
    uz /= rho;
    let usq = ux * ux + uy * uy + uz * uz;

    // --- BGK collision ---------------------------------------------------
    for (var q: u32 = 0u; q < 19u; q++) {
        let cu  = f32(CX[q]) * ux + f32(CY[q]) * uy + f32(CZ[q]) * uz;
        let feq = W[q] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usq);
        f_out[idx(x, y, z, q, nx, ny, nz)] = f[q] - omega * (f[q] - feq);
    }
}
