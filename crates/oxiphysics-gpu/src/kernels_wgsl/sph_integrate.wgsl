// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH symplectic (semi-implicit) Euler integration.
//
// Given per-particle acceleration (already including gravity from the force
// kernel), updates velocity first then position:  v += a dt ;  x += v dt.
// Position and velocity are xyz-interleaved storage buffers, updated in place.
// This is exactly the CPU reference's integration step.

struct IntegrateParams {
    n: u32,
    dt: f32,
    pad0: u32,
    pad1: u32,
}

@group(0) @binding(0) var<storage, read_write> positions:  array<f32>;
@group(0) @binding(1) var<storage, read_write> velocities: array<f32>;
@group(0) @binding(2) var<storage, read>       accel:      array<f32>;
@group(0) @binding(3) var<storage, read>       params:     IntegrateParams;

@compute @workgroup_size(64)
fn sph_integrate(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let dt = params.dt;

    let b = i * 3u;
    var vx = velocities[b]     + accel[b]      * dt;
    var vy = velocities[b + 1u] + accel[b + 1u] * dt;
    var vz = velocities[b + 2u] + accel[b + 2u] * dt;

    velocities[b]      = vx;
    velocities[b + 1u] = vy;
    velocities[b + 2u] = vz;

    positions[b]      = positions[b]      + vx * dt;
    positions[b + 1u] = positions[b + 1u] + vy * dt;
    positions[b + 2u] = positions[b + 2u] + vz * dt;
}
