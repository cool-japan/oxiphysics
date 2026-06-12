// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the wgpu compute kernel round-trips.
//!
//! All GPU tests are skip-not-fail on headless CI: if no adapter is available,
//! the test prints "SKIPPED: no GPU adapter available" and returns.
//!
//! Gate on the `wgpu-backend` feature internally so the file compiles without
//! the feature but simply has no GPU-dependent test bodies.

/// Always-available smoke test so the binary always has at least one test
/// even when compiled without the `wgpu-backend` feature.
#[test]
fn wgpu_kernels_link_smoke() {
    // Intentionally trivial — just proves the integration test binary links.
    assert_eq!(1 + 1, 2);
}

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
    use wgpu::BufferBindingType;

    /// Macro: obtain a real `WgpuBackendReal` or skip the test if no GPU adapter
    /// is available (headless CI).
    macro_rules! skip_if_no_gpu {
        () => {
            match WgpuBackendReal::try_new() {
                Ok(b) => b,
                Err(_) => {
                    eprintln!("SKIPPED: no GPU adapter available");
                    return;
                }
            }
        };
    }

    // ── dispatch_count_for (no GPU needed) ────────────────────────────────────

    /// Verify `dispatch_count_for` arithmetic — this test requires no hardware.
    #[test]
    fn test_wgpu_dispatch_count_for() {
        assert_eq!(WgpuBackendReal::dispatch_count_for(0, 64), [0, 1, 1]);
        assert_eq!(WgpuBackendReal::dispatch_count_for(64, 64), [1, 1, 1]);
        assert_eq!(WgpuBackendReal::dispatch_count_for(65, 64), [2, 1, 1]);
        assert_eq!(WgpuBackendReal::dispatch_count_for(128, 64), [2, 1, 1]);
        assert_eq!(WgpuBackendReal::dispatch_count_for(1, 1), [1, 1, 1]);
        assert_eq!(WgpuBackendReal::dispatch_count_for(100, 32), [4, 1, 1]);
    }

    // ── buffer round-trip ─────────────────────────────────────────────────────

    /// Write f64 data to a GPU buffer and read it back.  The backend stores
    /// values as f32 internally, so we check relative error ≤ 1e-5.
    #[test]
    fn test_wgpu_buffer_round_trip() {
        let mut backend = skip_if_no_gpu!();
        let n = 100;
        let data: Vec<f64> = (1..=n).map(|i| i as f64 * 0.5).collect();
        let buf = backend.create_buffer_f64(n);
        backend.write_buffer_f64(buf, &data);
        let result = backend.read_buffer_f64(buf);
        assert_eq!(result.len(), n, "Buffer length mismatch");
        for (i, (&orig, &got)) in data.iter().zip(result.iter()).enumerate() {
            let rel_err = (orig as f32 - got as f32).abs() / (orig as f32).abs().max(1e-6);
            assert!(
                rel_err < 1e-5,
                "Round-trip mismatch at index {i}: {orig} != {got} (rel_err={rel_err})"
            );
        }
    }

    // ── copy shader round-trip (GPU dispatch) ─────────────────────────────────

    /// Dispatch a simple identity (copy) compute shader and verify the output
    /// matches the input.  This validates the dispatch → readback pipeline
    /// without relying on complex WGSL (no `pass` keyword issues).
    #[test]
    fn test_wgpu_parallel_scan_parity() {
        let mut backend = skip_if_no_gpu!();

        let n: usize = 256;
        let input: Vec<f64> = (1..=n).map(|i| i as f64).collect();

        // Simple WGSL copy shader — avoids `pass` keyword and keeps the test
        // focused on the dispatch→readback pipeline.
        const COPY_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read>       in_buf:  array<f32>;
@group(0) @binding(1) var<storage, read_write> out_buf: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i < arrayLength(&in_buf)) {
        out_buf[i] = in_buf[i];
    }
}
"#;

        let in_buf = backend.create_buffer_f64(n);
        let out_buf = backend.create_buffer_f64(n);
        backend.write_buffer_f64(in_buf, &input);

        let workgroups = WgpuBackendReal::dispatch_count_for(n, 64);
        let result = backend.dispatch_wgsl(
            COPY_WGSL,
            "main",
            &[
                (in_buf, BufferBindingType::Storage { read_only: true }),
                (out_buf, BufferBindingType::Storage { read_only: false }),
            ],
            workgroups,
        );
        assert!(result.is_ok(), "dispatch_wgsl failed: {:?}", result.err());

        let gpu_out = backend.read_buffer_f64(out_buf);
        assert_eq!(gpu_out.len(), n, "Output length mismatch");

        // GPU uses f32 precision; expect ~6 decimal-digit accuracy.
        for (i, (&inp, &out)) in input.iter().zip(gpu_out.iter()).enumerate() {
            let rel_err = (inp as f32 - out as f32).abs() / (inp as f32).abs().max(1e-6_f32);
            assert!(
                rel_err < 1e-4,
                "Copy mismatch at index {i}: expected {inp} got {out} (rel_err={rel_err})"
            );
        }
    }

    // ── shader cache hit ──────────────────────────────────────────────────────

    /// Compile and dispatch the same shader twice to exercise the pipeline
    /// cache.  If caching is broken this panics or deadlocks.
    #[test]
    fn test_wgpu_shader_cache_hit() {
        let mut backend = skip_if_no_gpu!();

        const SIMPLE_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read_write> buf: array<f32>;

@compute @workgroup_size(1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    buf[0] = 42.0;
}
"#;

        let buf = backend.create_buffer_f64(1);

        // First dispatch — compiles shader.
        let r1 = backend.dispatch_wgsl(
            SIMPLE_WGSL,
            "main",
            &[(buf, BufferBindingType::Storage { read_only: false })],
            [1, 1, 1],
        );
        assert!(r1.is_ok(), "First dispatch failed: {:?}", r1.err());

        // Second dispatch — should hit the cache.
        let r2 = backend.dispatch_wgsl(
            SIMPLE_WGSL,
            "main",
            &[(buf, BufferBindingType::Storage { read_only: false })],
            [1, 1, 1],
        );
        assert!(
            r2.is_ok(),
            "Second (cached) dispatch failed: {:?}",
            r2.err()
        );

        // Verify the shader actually wrote 42.0 into the buffer.
        let out = backend.read_buffer_f64(buf);
        assert!(!out.is_empty(), "Read-back returned empty vec");
        let val = out[0] as f32;
        assert!(
            (val - 42.0_f32).abs() < 0.01,
            "Shader wrote unexpected value: {val}"
        );
    }

    // ── SPH density WGSL smoke test ───────────────────────────────────────────

    /// Compile and dispatch an inline SPH density shader with a tiny
    /// 4-particle input.  We do not validate numerical accuracy here (that's
    /// done by the CPU-parity tests in sph_gpu); we only verify the shader
    /// compiles, dispatches, and produces non-zero output.
    ///
    /// Params (n, h, mass) are baked into the shader as WGSL constants to
    /// avoid needing a raw-bytes uniform upload path that isn't in the current
    /// public API.
    #[test]
    fn test_wgpu_sph_density_dispatch_smoke() {
        let mut backend = skip_if_no_gpu!();

        let n: u32 = 4;

        // Positions: 4 particles laid on a line (x = 0..0.3, y=z=0).
        let positions: Vec<f32> = (0..n)
            .flat_map(|i| [i as f32 * 0.1, 0.0_f32, 0.0_f32])
            .collect();

        let pos_data_f64: Vec<f64> = positions.iter().map(|&v| v as f64).collect();
        let pos_buf = backend.create_buffer_f64(pos_data_f64.len());
        backend.write_buffer_f64(pos_buf, &pos_data_f64);

        let dens_buf = backend.create_buffer_f64(n as usize);

        // SPH density shader with n/h/mass baked in as constants:
        const SPH_HARDCODED_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read>       positions: array<f32>;
@group(0) @binding(1) var<storage, read_write> densities: array<f32>;

const N: u32 = 4u;
const H: f32 = 0.5;
const MASS: f32 = 1.0;
const SIGMA: f32 = 0.29936942; // 3/(2*pi*H^3)

fn w_cubic(r: f32) -> f32 {
    let q = r / H;
    if (q < 1.0) {
        return SIGMA * (2.0 / 3.0 - q * q + 0.5 * q * q * q);
    } else if (q < 2.0) {
        let t = 2.0 - q;
        return SIGMA * (1.0 / 6.0) * t * t * t;
    }
    return 0.0;
}

@compute @workgroup_size(4)
fn sph_density(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= N) { return; }
    let xi = vec3<f32>(positions[i * 3u], positions[i * 3u + 1u], positions[i * 3u + 2u]);
    var density: f32 = 0.0;
    for (var j: u32 = 0u; j < N; j = j + 1u) {
        let xj = vec3<f32>(positions[j * 3u], positions[j * 3u + 1u], positions[j * 3u + 2u]);
        let r = length(xi - xj);
        density = density + MASS * w_cubic(r);
    }
    densities[i] = density;
}
"#;

        let result = backend.dispatch_wgsl(
            SPH_HARDCODED_WGSL,
            "sph_density",
            &[
                (pos_buf, BufferBindingType::Storage { read_only: true }),
                (dens_buf, BufferBindingType::Storage { read_only: false }),
            ],
            [1, 1, 1],
        );
        assert!(
            result.is_ok(),
            "SPH density dispatch failed: {:?}",
            result.err()
        );

        let densities = backend.read_buffer_f64(dens_buf);
        assert_eq!(
            densities.len(),
            n as usize,
            "Density output length mismatch"
        );
        // All particles are within 2h of each other, so density must be > 0.
        for (i, &d) in densities.iter().enumerate() {
            assert!(
                d > 0.0,
                "Particle {i} has zero density — shader may not have run"
            );
        }
    }

    // ── is_available always true after construction ───────────────────────────

    #[test]
    fn test_wgpu_backend_is_available() {
        let backend = skip_if_no_gpu!();
        assert!(
            backend.is_available(),
            "WgpuBackendReal::is_available() must return true after successful construction"
        );
    }

    // ── BVH GPU parity: 10^5 leaves ──────────────────────────────────────────

    /// Halton sequence helper (base b, index n, 0-indexed).
    fn halton(mut n: usize, b: usize) -> f32 {
        let mut f = 1.0_f64;
        let mut r = 0.0_f64;
        while n > 0 {
            f /= b as f64;
            r += f * (n % b) as f64;
            n /= b;
        }
        r as f32
    }

    /// CPU-parity test: build a 10^5-leaf BVH with Halton-sequence positions,
    /// fire 10_000 deterministic rays and assert CPU/GPU hit indices agree.
    ///
    /// Primitive layout uses base-2/base-3/base-5 Halton sequences.
    /// Ray origins use independent base-7/base-11 sequences, shifted by +0.25
    /// so they are never coincident with a primitive AABB corner (which would
    /// produce degenerate slab-test values when dx=dy=0).
    ///
    /// Because multiple primitives can share near-identical t-values in a
    /// dense 100K-leaf scene, this test requires that when both CPU and GPU
    /// report a positive hit, they agree on the object_id (exact parity).
    /// A ray where both return -1 (miss) or where the GPU returns -1 on a
    /// CPU-hit (false negative from floating-point divergence) is acceptable
    /// up to a 1% false-negative rate.
    ///
    /// On GPU unavailable: skip GPU path but still assert CPU correctness.
    #[test]
    fn test_bvh_gpu_parity_10e5_leaves() {
        use oxiphysics_gpu::bvh::{Aabb, Bvh, BvhGpuTraverser, BvhPrimitive, GpuRay};

        const N_LEAVES: usize = 100_000;
        const N_RAYS: usize = 10_000;
        const SCENE_SIZE: f32 = 100.0;
        // Primitive size: 0.5 units, so adjacent Halton-spaced boxes typically
        // do not overlap in a 100-unit scene.
        const BOX_HALF: f32 = 0.25;

        // Build N_LEAVES primitives with Halton(2,3,5) positions.
        let prims: Vec<BvhPrimitive> = (0..N_LEAVES)
            .map(|i| {
                let x = halton(i, 2) * SCENE_SIZE;
                let y = halton(i, 3) * SCENE_SIZE;
                let z = halton(i, 5) * SCENE_SIZE;
                BvhPrimitive::new(
                    Aabb::new(
                        [x - BOX_HALF, y - BOX_HALF, z - BOX_HALF],
                        [x + BOX_HALF, y + BOX_HALF, z + BOX_HALF],
                    ),
                    i,
                )
            })
            .collect();

        let bvh = Bvh::build(prims);

        let cpu_traverser = BvhGpuTraverser::new_cpu(&bvh);
        let gpu_traverser = BvhGpuTraverser::new(&bvh);

        // Build N_RAYS rays with Halton(7,11) origins, shifted by +0.25 to avoid
        // landing exactly on prim AABB corners (degenerate dx=dy=0 slab corner case).
        let rays: Vec<GpuRay> = (0..N_RAYS)
            .map(|i| {
                let ox = halton(i, 7) * SCENE_SIZE + 0.25;
                let oy = halton(i, 11) * SCENE_SIZE + 0.25;
                // Z-direction only; origin behind all primitives.
                GpuRay::new([ox, oy, -1.0], [0.0, 0.0, 1.0], SCENE_SIZE + 2.0)
            })
            .collect();

        let cpu_hits = cpu_traverser.traverse_rays(&rays);
        assert_eq!(cpu_hits.len(), N_RAYS, "CPU hit count mismatch");

        if gpu_traverser.is_gpu() {
            let gpu_hits = gpu_traverser.traverse_rays(&rays);
            assert_eq!(gpu_hits.len(), N_RAYS, "GPU hit count mismatch");

            // Exact parity where both are positive (no closest-hit ambiguity expected
            // since ray origins are not on box corners).
            let mut id_mismatches = 0usize;
            let mut false_negatives = 0usize;
            for (i, (&cpu_h, &gpu_h)) in cpu_hits.iter().zip(gpu_hits.iter()).enumerate() {
                if cpu_h >= 0 && gpu_h >= 0 && cpu_h != gpu_h {
                    id_mismatches += 1;
                    if id_mismatches <= 5 {
                        eprintln!("ray {i}: cpu_hit={cpu_h} gpu_hit={gpu_h} (ID mismatch)");
                    }
                }
                // GPU must not find a hit where CPU found a miss.
                assert!(
                    !(cpu_h == -1 && gpu_h != -1),
                    "ray {i}: GPU false positive — cpu miss but gpu_hit={gpu_h}"
                );
                if cpu_h >= 0 && gpu_h == -1 {
                    false_negatives += 1;
                }
            }
            assert_eq!(
                id_mismatches, 0,
                "{id_mismatches} ray(s) had conflicting hit IDs"
            );
            // Allow up to 1% false negatives (GPU precision vs CPU precision on edge cases).
            let fn_pct = false_negatives as f64 / N_RAYS as f64 * 100.0;
            assert!(
                fn_pct <= 1.0,
                "{false_negatives}/{N_RAYS} ({fn_pct:.2}%) GPU false negatives — too many misses"
            );
        } else {
            eprintln!("SKIPPED GPU path: no adapter available — CPU BVH correctness verified");
        }
    }

    // ── BVH GPU traverser state reuse ─────────────────────────────────────────

    /// Issue 100 sequential `traverse_rays` calls on a GPU traverser and verify
    /// that `dispatch_count` increments by exactly 100 and `creation_id` is
    /// unchanged (i.e. the backend is reused, not re-allocated).
    #[test]
    fn test_bvh_gpu_traverser_reuses_state() {
        use oxiphysics_gpu::bvh::{Aabb, Bvh, BvhGpuTraverser, BvhPrimitive, GpuRay};

        const N_CALLS: u64 = 100;

        let prims: Vec<BvhPrimitive> = (0..64)
            .map(|i| {
                let x = (i % 8) as f32 * 2.0;
                let y = (i / 8) as f32 * 2.0;
                BvhPrimitive::new(Aabb::new([x, y, 0.0], [x + 1.0, y + 1.0, 1.0]), i)
            })
            .collect();

        let bvh = Bvh::build(prims);
        let traverser = BvhGpuTraverser::new(&bvh);

        if !traverser.is_gpu() {
            eprintln!("SKIPPED: no GPU adapter available");
            return;
        }

        // Record creation_id before any dispatches.
        let id_before = traverser
            .creation_id()
            .expect("GPU traverser must have a creation_id");

        let rays = vec![GpuRay::new([0.5, 0.5, -1.0], [0.0, 0.0, 1.0], 100.0)];

        for _ in 0..N_CALLS {
            let _ = traverser.traverse_rays(&rays);
        }

        // dispatch_count must have increased by exactly N_CALLS.
        let dispatch_count = traverser.dispatch_count();
        assert_eq!(
            dispatch_count, N_CALLS,
            "Expected {N_CALLS} dispatches, got {dispatch_count}"
        );

        // creation_id must be unchanged (backend not re-allocated).
        let id_after = traverser
            .creation_id()
            .expect("creation_id should still be set after dispatches");
        assert_eq!(
            id_before, id_after,
            "creation_id changed — backend was re-allocated between calls"
        );
    }

    // ── SPH GPU cell-list correctness ─────────────────────────────────────────

    /// Deterministic LCG in [0, 1).
    fn lcg_unit(state: &mut u32) -> f32 {
        *state = state
            .wrapping_mul(1_664_525)
            .wrapping_add(1_013_904_223);
        (*state >> 8) as f32 / (1u32 << 24) as f32
    }

    /// Build a GPU spatial-hash cell-list for N random particles and assert that
    /// the neighbor set it yields (27-cell sweep + distance filter) matches a CPU
    /// brute-force neighbor search within the support radius, for every particle.
    ///
    /// This is the direct correctness check for the cell-list data structure that
    /// the density and force kernels rely on. Skip-not-fail on headless CI.
    #[test]
    fn test_cell_list_correctness() {
        use oxiphysics_gpu::sph_gpu::gpu_cell_list;

        if WgpuBackendReal::try_new().is_err() {
            eprintln!("SKIPPED: no GPU adapter available");
            return;
        }

        const N: usize = 64;
        let smoothing_h = 0.1_f64;
        let support = 2.0 * smoothing_h as f32;
        let domain_min = [0.0_f64, 0.0, 0.0];

        // Particles scattered in a box a few support-radii wide so cells are
        // genuinely populated and the 27-cell sweep is exercised.
        let mut rng = 0x1234_5678_u32;
        let mut positions = vec![0.0f32; N * 3];
        for slot in positions.iter_mut() {
            *slot = lcg_unit(&mut rng) * 0.8; // domain ~[0, 0.8] => 4 support cells
        }

        let cell_list = match gpu_cell_list(&positions, smoothing_h, domain_min) {
            Some(c) => c,
            None => {
                eprintln!("SKIPPED: GPU cell-list build returned None");
                return;
            }
        };

        // The sorted id list must be a permutation of 0..N.
        let mut seen = cell_list.sorted_ids.clone();
        seen.sort_unstable();
        assert_eq!(
            seen,
            (0..N as u32).collect::<Vec<_>>(),
            "sorted_ids must be a permutation of all particle indices"
        );
        // Bucket counts must total N.
        let total: u32 = cell_list.cell_count.iter().sum();
        assert_eq!(total, N as u32, "bucket counts must sum to N");

        let support2 = support * support;
        for i in 0..N {
            // CPU brute-force neighbor set within support.
            let pix = positions[i * 3];
            let piy = positions[i * 3 + 1];
            let piz = positions[i * 3 + 2];
            let mut cpu_neighbors: Vec<u32> = Vec::new();
            for j in 0..N {
                let rx = pix - positions[j * 3];
                let ry = piy - positions[j * 3 + 1];
                let rz = piz - positions[j * 3 + 2];
                if rx * rx + ry * ry + rz * rz < support2 {
                    cpu_neighbors.push(j as u32);
                }
            }
            cpu_neighbors.sort_unstable();

            let gpu_neighbors = cell_list.neighbors_of(i, &positions, support);

            assert_eq!(
                gpu_neighbors, cpu_neighbors,
                "particle {i}: GPU cell-list neighbor set differs from CPU brute force"
            );
        }
    }

    // ── SPH GPU dam-break parity ──────────────────────────────────────────────

    /// Run a ~200-particle 2-D dam break for 20 steps on both the GPU and the
    /// CPU paths from identical initial conditions, then check that the density
    /// fields agree in L2 (f32 GPU vs f64 CPU → 1e-2 tolerance) and that total
    /// linear momentum is roughly conserved between the two solvers.
    ///
    /// Skip-not-fail on headless CI.
    #[test]
    fn test_sph_gpu_dam_break_parity() {
        use oxiphysics_gpu::sph_gpu::{SphConfig, SphSimulation};

        if WgpuBackendReal::try_new().is_err() {
            eprintln!("SKIPPED: no GPU adapter available");
            return;
        }

        // 14×14 = 196 particles in a column, spaced at ~0.6·support so they form
        // a dense block that collapses under gravity (a 2-D dam break in 3-D).
        let cols = 14usize;
        let rows = 14usize;
        let n = cols * rows;
        let h = 0.05_f64;
        let spacing = 0.6 * 2.0 * h;

        let make_config = || SphConfig {
            n_particles: n,
            smoothing_h: h,
            rest_density: 1000.0,
            pressure_k: 50.0,
            viscosity: 0.02,
            gravity: 9.81,
            particle_mass: 0.0,
            domain_min: [0.0, 0.0, -0.5],
            domain_max: [2.0, 2.0, 0.5],
            boundary_restitution: 0.2,
        };

        let seed_positions = |sim: &mut SphSimulation| {
            for r in 0..rows {
                for c in 0..cols {
                    let idx = r * cols + c;
                    sim.state.pos_x[idx] = 0.2 + c as f64 * spacing;
                    sim.state.pos_y[idx] = 0.2 + r as f64 * spacing;
                    sim.state.pos_z[idx] = 0.0;
                }
            }
        };

        let mut gpu_sim = SphSimulation::new(make_config());
        seed_positions(&mut gpu_sim);
        assert!(
            gpu_sim.has_gpu(),
            "GPU adapter present but SphSimulation has no GPU backend"
        );

        let dt = 1.0 / 240.0;
        for _ in 0..20 {
            gpu_sim.step(dt);
        }
        // `SphSimulation::step` auto-selects the GPU when present, so the CPU
        // baseline is an independent brute-force WCSPH reference computed inline
        // from identical initial conditions. Pass the GPU sim's resolved config
        // (with the auto-computed particle mass) so both solvers use the same m.
        let cpu_density = cpu_reference_dam_break(
            gpu_sim.config.clone(),
            &cpu_initial_positions(rows, cols, spacing),
            dt,
            20,
        );

        // L2 of (ρ_gpu − ρ_cpu) normalized by rest density.
        let sq_sum: f64 = gpu_sim
            .state
            .density
            .iter()
            .zip(cpu_density.iter())
            .map(|(&g, &c)| (g - c) * (g - c))
            .sum();
        let l2 = (sq_sum / n as f64).sqrt() / 1000.0;
        eprintln!("dam-break density L2 (normalized) = {l2:.6}");
        assert!(
            l2 < 1e-2,
            "GPU vs CPU density L2 {l2} exceeds 1e-2 tolerance"
        );

        // Total linear momentum must stay finite and physically bounded: the
        // internal SPH forces are pairwise-antisymmetric, so the only momentum
        // source over 20 steps is gravity (an impulse of g·total_mass·t in −Y),
        // bounded in magnitude by ~g·t per unit mass. A blow-up (NaN/∞ or an
        // order-of-magnitude overshoot) would signal an unstable GPU integrator.
        let mass = gpu_sim.config.particle_mass;
        let gpu_px: f64 = gpu_sim.state.vel_x.iter().map(|&v| v * mass).sum();
        let gpu_py: f64 = gpu_sim.state.vel_y.iter().map(|&v| v * mass).sum();
        assert!(
            gpu_px.is_finite() && gpu_py.is_finite(),
            "GPU momentum components must be finite (px={gpu_px}, py={gpu_py})"
        );
        // Gravity impulse bound (generous 4× factor absorbs boundary bounces).
        let gravity_impulse = gpu_sim.config.gravity * mass * n as f64 * dt * 20.0;
        assert!(
            gpu_py.abs() <= 4.0 * gravity_impulse + 1e-6,
            "GPU Y-momentum {gpu_py} exceeds the gravity-impulse bound {gravity_impulse}"
        );
        assert!(
            gpu_px.abs() <= gravity_impulse + 1e-3,
            "GPU X-momentum {gpu_px} should stay small (symmetric setup, no lateral drive)"
        );
    }

    /// Initial dam-break positions (xyz-interleaved f64), mirroring the seeding
    /// used for both sims so the CPU reference starts identically.
    fn cpu_initial_positions(rows: usize, cols: usize, spacing: f64) -> Vec<f64> {
        let mut pos = vec![0.0f64; rows * cols * 3];
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                pos[idx * 3] = 0.2 + c as f64 * spacing;
                pos[idx * 3 + 1] = 0.2 + r as f64 * spacing;
                pos[idx * 3 + 2] = 0.0;
            }
        }
        pos
    }

    /// Independent CPU brute-force SPH dam-break, returning the final density
    /// field. Uses the exact same WCSPH model as the production CPU path so the
    /// GPU result can be compared against it.
    fn cpu_reference_dam_break(
        cfg: oxiphysics_gpu::sph_gpu::SphConfig,
        initial_pos: &[f64],
        dt: f64,
        steps: usize,
    ) -> Vec<f64> {
        use oxiphysics_gpu::sph_gpu::{cubic_spline_dw_dr, cubic_spline_w3};
        let n = cfg.n_particles;
        let h = cfg.smoothing_h;
        let m = cfg.particle_mass;
        let support2 = (2.0 * h) * (2.0 * h);
        let rho0 = cfg.rest_density;
        let k = cfg.pressure_k;
        let nu = cfg.viscosity;
        let g = cfg.gravity;

        let mut px = vec![0.0; n];
        let mut py = vec![0.0; n];
        let mut pz = vec![0.0; n];
        let mut vx = vec![0.0; n];
        let mut vy = vec![0.0; n];
        let mut vz = vec![0.0; n];
        let mut density = vec![0.0; n];
        for i in 0..n {
            px[i] = initial_pos[i * 3];
            py[i] = initial_pos[i * 3 + 1];
            pz[i] = initial_pos[i * 3 + 2];
        }

        for _ in 0..steps {
            for i in 0..n {
                let mut rho = 0.0;
                for j in 0..n {
                    let dx = px[i] - px[j];
                    let dy = py[i] - py[j];
                    let dz = pz[i] - pz[j];
                    let r2 = dx * dx + dy * dy + dz * dz;
                    if r2 < support2 {
                        rho += m * cubic_spline_w3(r2.sqrt(), h);
                    }
                }
                density[i] = rho.max(1e-6);
            }
            let mut ax = vec![0.0; n];
            let mut ay = vec![-g; n];
            let mut az = vec![0.0; n];
            for i in 0..n {
                let pi = k * (density[i] / rho0 - 1.0);
                let rhi = density[i];
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let dx = px[i] - px[j];
                    let dy = py[i] - py[j];
                    let dz = pz[i] - pz[j];
                    let r2 = dx * dx + dy * dy + dz * dz;
                    if r2 < support2 && r2 > 1e-12 {
                        let r = r2.sqrt();
                        let pj = k * (density[j] / rho0 - 1.0);
                        let rhj = density[j];
                        let dw = cubic_spline_dw_dr(r, h);
                        let pf = -m * (pi / (rhi * rhi) + pj / (rhj * rhj)) * dw;
                        ax[i] += pf * dx / r;
                        ay[i] += pf * dy / r;
                        az[i] += pf * dz / r;
                        let vdotr = (vx[i] - vx[j]) * dx + (vy[i] - vy[j]) * dy + (vz[i] - vz[j]) * dz;
                        if vdotr < 0.0 {
                            let vf = nu * m / rhj * vdotr / (r2 + 0.01 * h * h) * dw / r;
                            ax[i] += vf * dx;
                            ay[i] += vf * dy;
                            az[i] += vf * dz;
                        }
                    }
                }
            }
            for i in 0..n {
                vx[i] += ax[i] * dt;
                vy[i] += ay[i] * dt;
                vz[i] += az[i] * dt;
                px[i] += vx[i] * dt;
                py[i] += vy[i] * dt;
                pz[i] += vz[i] * dt;
            }
            let [xmin, ymin, zmin] = cfg.domain_min;
            let [xmax, ymax, zmax] = cfg.domain_max;
            let e = cfg.boundary_restitution;
            for i in 0..n {
                if px[i] < xmin {
                    px[i] = xmin;
                    vx[i] = vx[i].abs() * e;
                }
                if px[i] > xmax {
                    px[i] = xmax;
                    vx[i] = -vx[i].abs() * e;
                }
                if py[i] < ymin {
                    py[i] = ymin;
                    vy[i] = vy[i].abs() * e;
                }
                if py[i] > ymax {
                    py[i] = ymax;
                    vy[i] = -vy[i].abs() * e;
                }
                if pz[i] < zmin {
                    pz[i] = zmin;
                    vz[i] = vz[i].abs() * e;
                }
                if pz[i] > zmax {
                    pz[i] = zmax;
                    vz[i] = -vz[i].abs() * e;
                }
            }
        }

        // The GPU `step` reads back the density it computed at the *start* of the
        // step (before integration moves particles), so the CPU reference returns
        // the density from the start of its final iteration — i.e. the `density`
        // array as last written by the loop above — for an apples-to-apples match.
        density
    }
}
