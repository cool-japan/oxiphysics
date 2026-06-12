// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH axis-aligned box boundary enforcement (wall clamp + restitution bounce).
//
// For each particle, clamps the position to the domain AABB and reflects the
// outward velocity component, scaled by the restitution coefficient. The
// reflection rule matches the CPU reference's `reflect!` macro: on the min wall
// the velocity becomes |v| * e (pushed inward), on the max wall it becomes
// -|v| * e (pushed inward), per axis. Positions/velocities are xyz-interleaved
// storage buffers updated in place.

struct BoundaryParams {
    n: u32,
    restitution: f32,
    min_x: f32,
    min_y: f32,
    min_z: f32,
    max_x: f32,
    max_y: f32,
    max_z: f32,
}

@group(0) @binding(0) var<storage, read_write> positions:  array<f32>;
@group(0) @binding(1) var<storage, read_write> velocities: array<f32>;
@group(0) @binding(2) var<storage, read>       params:     BoundaryParams;

@compute @workgroup_size(64)
fn sph_boundary(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let e = params.restitution;
    let b = i * 3u;

    var px = positions[b];
    var py = positions[b + 1u];
    var pz = positions[b + 2u];
    var vx = velocities[b];
    var vy = velocities[b + 1u];
    var vz = velocities[b + 2u];

    if (px < params.min_x) { px = params.min_x; vx =  abs(vx) * e; }
    if (px > params.max_x) { px = params.max_x; vx = -abs(vx) * e; }
    if (py < params.min_y) { py = params.min_y; vy =  abs(vy) * e; }
    if (py > params.max_y) { py = params.max_y; vy = -abs(vy) * e; }
    if (pz < params.min_z) { pz = params.min_z; vz =  abs(vz) * e; }
    if (pz > params.max_z) { pz = params.max_z; vz = -abs(vz) * e; }

    positions[b]      = px;
    positions[b + 1u] = py;
    positions[b + 2u] = pz;
    velocities[b]      = vx;
    velocities[b + 1u] = vy;
    velocities[b + 2u] = vz;
}
