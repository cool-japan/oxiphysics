// BVH traversal — AABB slab-test ray traversal
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Flat BVH layout (AoS):
//   Each FlatBvhNode is 8 f32 words:
//     [0..2]  aabb.min  (3 × f32)
//     [3..5]  aabb.max  (3 × f32)
//     [6]     left_first (u32 cast to f32 via bitcast)
//     [7]     count      (u32 cast to f32 via bitcast)
//
// Each GpuRay is 8 f32 words:
//     [0..2]  origin    (3 × f32)
//     [3..5]  direction (3 × f32)
//     [6]     max_t     (f32)
//     [7]     _pad      (f32, unused)
//
// Primitive layout: prim_aabbs[prim_idx * 6 + 0..5] = [min_x, min_y, min_z, max_x, max_y, max_z]
//
// Bindings:
//   binding(0) — params        storage read  [n_nodes: u32, n_rays: u32, n_prims: u32, _pad]
//   binding(1) — bvh_nodes     storage read  [n_nodes * 8 f32]
//   binding(2) — rays          storage read  [n_rays  * 8 f32]
//   binding(3) — hit_results   storage r/w   [n_rays  i32]  (leaf object_id or -1)
//   binding(4) — prim_indices  storage read  [n_prims u32]  (leaf prim index array)
//   binding(5) — object_ids    storage read  [n_prims i32]  (object_id per primitive)
//   binding(6) — prim_aabbs    storage read  [n_prims * 6 f32] (primitive AABBs)
//
// Dispatch: (ceil(n_rays / 64), 1, 1)

@group(0) @binding(0) var<storage, read>       params:       array<u32>;
@group(0) @binding(1) var<storage, read>       bvh_nodes:    array<f32>;
@group(0) @binding(2) var<storage, read>       rays:         array<f32>;
@group(0) @binding(3) var<storage, read_write> hit_results:  array<i32>;
@group(0) @binding(4) var<storage, read>       prim_indices: array<u32>;
@group(0) @binding(5) var<storage, read>       object_ids:   array<i32>;
@group(0) @binding(6) var<storage, read>       prim_aabbs:   array<f32>;

const NODE_STRIDE: u32 = 8u;
const RAY_STRIDE:  u32 = 8u;
const PRIM_STRIDE: u32 = 6u;

// Read a node field: word offset within node (0-7)
fn node_f32(node_idx: u32, word: u32) -> f32 {
    return bvh_nodes[node_idx * NODE_STRIDE + word];
}

fn node_u32(node_idx: u32, word: u32) -> u32 {
    return bitcast<u32>(bvh_nodes[node_idx * NODE_STRIDE + word]);
}

// Ray–AABB slab test.  Returns t_near if hit (>= 0 and <= max_t), else 1e38.
fn ray_aabb_tnear(
    ox: f32, oy: f32, oz: f32,
    idx_x: f32, idx_y: f32, idx_z: f32,   // 1/dir components
    min_x: f32, min_y: f32, min_z: f32,
    max_x: f32, max_y: f32, max_z: f32,
    max_t: f32,
) -> f32 {
    let tx1 = (min_x - ox) * idx_x;
    let tx2 = (max_x - ox) * idx_x;
    var tmin = min(tx1, tx2);
    var tmax = max(tx1, tx2);

    let ty1 = (min_y - oy) * idx_y;
    let ty2 = (max_y - oy) * idx_y;
    tmin = max(tmin, min(ty1, ty2));
    tmax = min(tmax, max(ty1, ty2));

    let tz1 = (min_z - oz) * idx_z;
    let tz2 = (max_z - oz) * idx_z;
    tmin = max(tmin, min(tz1, tz2));
    tmax = min(tmax, max(tz1, tz2));

    if tmin <= tmax && tmax >= 0.0 && tmin <= max_t {
        return max(tmin, 0.0);
    }
    return 1e38;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n_nodes  = params[0];
    let n_rays   = params[1];
    let n_prims  = params[2];

    let ray_id = gid.x;
    if ray_id >= n_rays || n_nodes == 0u {
        return;
    }

    let base = ray_id * RAY_STRIDE;
    let ox    = rays[base + 0u];
    let oy    = rays[base + 1u];
    let oz    = rays[base + 2u];
    let dx    = rays[base + 3u];
    let dy    = rays[base + 4u];
    let dz    = rays[base + 5u];
    let max_t = rays[base + 6u];

    // Safe inverse direction (avoid div-by-zero with large finite sentinel)
    let eps = 1e-30;
    let idx_x = 1.0 / select(eps, dx, abs(dx) > eps);
    let idx_y = 1.0 / select(eps, dy, abs(dy) > eps);
    let idx_z = 1.0 / select(eps, dz, abs(dz) > eps);

    // Iterative BVH traversal using a stack (max depth 64)
    var stack: array<u32, 64>;
    var stack_top: i32 = 0;
    stack[0] = 0u;
    stack_top = 1;

    var best_hit: i32 = -1;
    var best_t: f32 = max_t;

    while stack_top > 0 {
        stack_top -= 1;
        let node_idx = stack[u32(stack_top)];

        if node_idx >= n_nodes {
            continue;
        }

        let min_x = node_f32(node_idx, 0u);
        let min_y = node_f32(node_idx, 1u);
        let min_z = node_f32(node_idx, 2u);
        let max_x = node_f32(node_idx, 3u);
        let max_y = node_f32(node_idx, 4u);
        let max_z = node_f32(node_idx, 5u);

        let t_node = ray_aabb_tnear(ox, oy, oz, idx_x, idx_y, idx_z,
                                    min_x, min_y, min_z,
                                    max_x, max_y, max_z,
                                    best_t);
        if t_node >= 1e37 {
            continue;
        }

        let left_first = node_u32(node_idx, 6u);
        let count      = node_u32(node_idx, 7u);

        if count > 0u {
            // Leaf node — test each primitive's AABB individually for closest hit
            let start = left_first;
            let end   = min(start + count, n_prims);
            var pi: u32 = start;
            while pi < end {
                let prim_idx = prim_indices[pi];
                let ab = prim_idx * PRIM_STRIDE;
                let p_min_x = prim_aabbs[ab + 0u];
                let p_min_y = prim_aabbs[ab + 1u];
                let p_min_z = prim_aabbs[ab + 2u];
                let p_max_x = prim_aabbs[ab + 3u];
                let p_max_y = prim_aabbs[ab + 4u];
                let p_max_z = prim_aabbs[ab + 5u];

                let t_prim = ray_aabb_tnear(ox, oy, oz, idx_x, idx_y, idx_z,
                                            p_min_x, p_min_y, p_min_z,
                                            p_max_x, p_max_y, p_max_z,
                                            best_t);
                if t_prim < best_t {
                    best_t   = t_prim;
                    best_hit = object_ids[prim_idx];
                }
                pi += 1u;
            }
        } else {
            // Internal node — push right child then left child (left processed first)
            let right_child = left_first;
            let left_child  = node_idx + 1u;
            if stack_top < 63 {
                if right_child < n_nodes {
                    stack[u32(stack_top)] = right_child;
                    stack_top += 1;
                }
                if stack_top < 63 && left_child < n_nodes {
                    stack[u32(stack_top)] = left_child;
                    stack_top += 1;
                }
            }
        }
    }

    hit_results[ray_id] = best_hit;
}
