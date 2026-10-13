// D3Q19 BGK collision (Guo forcing) + push streaming with bounce-back walls
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Dispatch: (ceil(nx/4), ceil(ny/4), ceil(nz/4))
// Each workgroup processes a 4×4×4 tile of lattice sites (64 invocations,
// within the WebGPU default limit of 256 invocations per workgroup).
//
// Algorithm (identical to `LbmSimulation::step_cpu`, collide-then-push):
//   1. Solid cells copy their populations through unchanged.
//   2. Fluid cells compute ρ and the Guo-shifted velocity u = (Σ fᵢ eᵢ + F/2)/ρ,
//      relax towards f^eq(ρ, u) with BGK and add the Guo (2002) source term
//      Sᵢ = (1 − ω/2) wᵢ [3(eᵢ − u) + 9(eᵢ·u) eᵢ]·F.
//   3. Each post-collision population is pushed to the neighbour x + eᵢ
//      (periodic wrap).  When the neighbour is solid, halfway bounce-back
//      sends it back into the opposite slot of the source cell; when that
//      solid neighbour lies on the lid plane y = ny − 1 the Ladd moving-wall
//      correction −6 wᵢ ρ (eᵢ·u_lid) is applied.
//
// Every output slot is written by exactly one invocation (streaming is a
// bijection per direction on the torus, and a bounce-back slot can only be
// reached from a solid cell, which never streams), so no atomics are needed.
//
// Buffer layout:
//   binding(0) — params (storage read, 10 x u32):
//                nx, ny, nz, omega, fx, fy, fz, lid_ux, lid_uy, lid_uz
//                (the last seven are f32 bit patterns)
//   binding(1) — f_in   (storage read,  19 * nx*ny*nz f32 values, q-major)
//   binding(2) — f_out  (storage r/w,   19 * nx*ny*nz f32 values, q-major)
//   binding(3) — solid  (storage read,  nx*ny*nz u32 values, 1 = solid)

@group(0) @binding(0) var<storage, read>       params: array<u32>;
@group(0) @binding(1) var<storage, read>       f_in:   array<f32>;
@group(0) @binding(2) var<storage, read_write> f_out:  array<f32>;
@group(0) @binding(3) var<storage, read>       solid:  array<u32>;

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

// Opposite direction: e[OPP[q]] = -e[q]
const OPP: array<u32, 19> = array<u32, 19>(0u, 2u, 1u, 4u, 3u, 6u, 5u, 10u, 9u, 8u, 7u, 14u, 13u, 12u, 11u, 18u, 17u, 16u, 15u);

// Flat cell index: x + nx * (y + ny * z)
fn cell_index(x: u32, y: u32, z: u32, nx: u32, ny: u32) -> u32 {
    return x + nx * (y + ny * z);
}

// Periodic neighbour coordinate c + d in [0, n)
fn wrap(c: u32, d: i32, n: u32) -> u32 {
    let ni = i32(n);
    return u32(((i32(c) + d) % ni + ni) % ni);
}

@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let nx    = params[0];
    let ny    = params[1];
    let nz    = params[2];
    let omega = bitcast<f32>(params[3]);
    let force = vec3<f32>(
        bitcast<f32>(params[4]),
        bitcast<f32>(params[5]),
        bitcast<f32>(params[6]),
    );
    let u_lid = vec3<f32>(
        bitcast<f32>(params[7]),
        bitcast<f32>(params[8]),
        bitcast<f32>(params[9]),
    );

    let x = gid.x;
    let y = gid.y;
    let z = gid.z;

    // Guard: discard threads outside the domain
    if x >= nx || y >= ny || z >= nz {
        return;
    }

    let nc   = nx * ny * nz;
    let cell = cell_index(x, y, z, nx, ny);

    // --- Solid cells: populations are carried through unchanged ----------
    if solid[cell] != 0u {
        for (var q: u32 = 0u; q < 19u; q++) {
            f_out[q * nc + cell] = f_in[q * nc + cell];
        }
        return;
    }

    // --- Load the local (post-streaming) populations ---------------------
    var f: array<f32, 19>;
    var rho: f32 = 0.0;
    var mom = vec3<f32>(0.0, 0.0, 0.0);
    for (var q: u32 = 0u; q < 19u; q++) {
        let fq = f_in[q * nc + cell];
        f[q] = fq;
        rho += fq;
        mom += fq * vec3<f32>(f32(CX[q]), f32(CY[q]), f32(CZ[q]));
    }

    // --- BGK collision with Guo forcing ----------------------------------
    // Degenerate (empty) cells skip the collision and only stream.
    var post: array<f32, 19>;
    if rho < 1e-10 {
        for (var q: u32 = 0u; q < 19u; q++) {
            post[q] = f[q];
        }
    } else {
        let u   = (mom + 0.5 * force) / rho;
        let usq = dot(u, u);
        let guo = 1.0 - 0.5 * omega;
        for (var q: u32 = 0u; q < 19u; q++) {
            let e   = vec3<f32>(f32(CX[q]), f32(CY[q]), f32(CZ[q]));
            let eu  = dot(e, u);
            let feq = W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
            let src = guo * W[q] * (3.0 * dot(e - u, force) + 9.0 * eu * dot(e, force));
            post[q] = f[q] - omega * (f[q] - feq) + src;
        }
    }

    // --- Streaming (push) with halfway bounce-back ------------------------
    for (var q: u32 = 0u; q < 19u; q++) {
        let xd  = wrap(x, CX[q], nx);
        let yd  = wrap(y, CY[q], ny);
        let zd  = wrap(z, CZ[q], nz);
        let dst = cell_index(xd, yd, zd, nx, ny);
        if solid[dst] != 0u {
            var value = post[q];
            if yd + 1u == ny {
                // Ladd moving-wall bounce-back on the lid plane.
                let e = vec3<f32>(f32(CX[q]), f32(CY[q]), f32(CZ[q]));
                value -= 6.0 * W[q] * rho * dot(e, u_lid);
            }
            f_out[OPP[q] * nc + cell] = value;
        } else {
            f_out[q * nc + dst] = post[q];
        }
    }
}
