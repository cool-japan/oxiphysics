// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Folded (2-D) dispatches of the shipped 1-D WGSL kernels.
//!
//! A 1-D kernel dispatched over more than 65535 workgroups needs its dispatch
//! folded into the y dimension (`dispatch_count_for`), and the kernel must
//! linearise `gid.x + gid.y * num_workgroups.x * WG` / `wid.x + wid.y *
//! num_workgroups.x` itself.  Instead of allocating hundreds of megabytes, these
//! tests fold *small* inputs with an artificial per-dimension limit of 3
//! workgroups (`dispatch_count_with_limit`), so a handful of blocks already spans
//! a multi-row grid with extra workgroups past the end of the data, and compare
//! every kernel's result with the CPU reference.
//!
//! All GPU tests skip (pass) when no adapter is available.

/// Always-available smoke test so the binary has at least one test without the
/// `wgpu-backend` feature.
#[test]
fn folded_dispatch_binary_links() {
    assert_eq!(2 + 2, 4);
}

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::compute::timestamp::dispatch_count_with_limit;
    use oxiphysics_gpu::compute::wgpu_backend::WgpuBufferHandle;
    use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
    use oxiphysics_gpu::kernels_wgsl::{
        COMPACT_WGSL, HISTOGRAM_WGSL, RADIX_HISTOGRAM_WGSL, RADIX_SCATTER_PAIRS_WGSL,
        RADIX_SCATTER_WGSL, REDUCE_WGSL, SCAN_ADD_WGSL, SCAN_WGSL, SPH_CELL_HISTOGRAM_WGSL,
        SPH_INTEGRATE_WGSL,
    };
    use wgpu::BufferBindingType;

    /// Per-dimension workgroup limit used to force folding of small inputs.
    const TINY_LIMIT: u32 = 3;
    /// The portable WebGPU limit: small inputs are not folded at all.
    const REAL_LIMIT: u32 = 65_535;
    const BLOCK: usize = 256;

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

    fn ro() -> BufferBindingType {
        BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> BufferBindingType {
        BufferBindingType::Storage { read_only: false }
    }

    /// Deterministic pseudo-random values in `0..modulus`.
    fn lcg(n: usize, seed: u32, modulus: u32) -> Vec<u32> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (state >> 8) % modulus
            })
            .collect()
    }

    /// Dispatch shape for `groups` workgroups under `limit`; with the tiny limit it
    /// must really be a multi-row grid with at least one spare workgroup.
    fn shape(groups: usize, limit: u32) -> [u32; 3] {
        let shape = dispatch_count_with_limit(groups, 1, limit).expect("dispatch shape");
        if limit == TINY_LIMIT {
            assert!(
                shape[1] > 1,
                "{groups} workgroups must fold under a limit of {limit}"
            );
            assert!(
                u64::from(shape[0]) * u64::from(shape[1]) > groups as u64,
                "the folded grid must contain extra workgroups past the data ({shape:?})"
            );
        }
        shape
    }

    fn upload_u32(backend: &mut WgpuBackendReal, data: &[u32]) -> WgpuBufferHandle {
        let handle = backend.create_buffer_storage((data.len() * 4) as u64);
        backend.queue_write_buffer_raw(&handle, bytemuck::cast_slice(data));
        handle
    }

    fn upload_f32(backend: &mut WgpuBackendReal, data: &[f32]) -> WgpuBufferHandle {
        let handle = backend.create_buffer_storage((data.len() * 4) as u64);
        backend.queue_write_buffer_f32(&handle, data);
        handle
    }

    fn zeros_u32(backend: &mut WgpuBackendReal, words: usize) -> WgpuBufferHandle {
        upload_u32(backend, &vec![0u32; words])
    }

    fn params(backend: &mut WgpuBackendReal, words: [u32; 4]) -> WgpuBufferHandle {
        upload_u32(backend, &words)
    }

    fn cpu_exclusive_scan(data: &[u32]) -> Vec<u32> {
        let mut acc = 0u32;
        data.iter()
            .map(|&v| {
                let out = acc;
                acc = acc.wrapping_add(v);
                out
            })
            .collect()
    }

    // ── scan ──────────────────────────────────────────────────────────────────

    /// Two-level exclusive scan driven from the host, as `gpu_primitives` does.
    fn gpu_scan(backend: &mut WgpuBackendReal, data: &[u32], limit: u32) -> (Vec<u32>, Vec<u32>) {
        let n = data.len();
        let blocks = n.div_ceil(BLOCK);
        let mut padded = data.to_vec();
        padded.resize(blocks * BLOCK, 0);
        let input = upload_u32(backend, &padded);
        let output = zeros_u32(backend, padded.len());
        // Sentinel words behind the block sums detect writes past the last block.
        let mut sums_init = vec![0u32; blocks];
        sums_init.extend(std::iter::repeat_n(0xDEAD_BEEF, 16));
        let block_sums = upload_u32(backend, &sums_init);
        let p = params(backend, [n as u32, 0, 0, 0]);
        let grid = shape(blocks, limit);
        backend
            .dispatch_wgsl(
                SCAN_WGSL,
                "scan_block",
                &[(input, ro()), (output, rw()), (block_sums, rw()), (p, ro())],
                grid,
            )
            .expect("scan_block");
        let sums_back = backend.read_buffer_u32(block_sums).expect("block sums");
        let sums = sums_back[..blocks].to_vec();
        let sentinels = sums_back[blocks..].to_vec();
        let offsets = upload_u32(backend, &cpu_exclusive_scan(&sums));
        backend
            .dispatch_wgsl(
                SCAN_ADD_WGSL,
                "add_block_offsets",
                &[(output, rw()), (offsets, ro()), (p, ro())],
                grid,
            )
            .expect("add_block_offsets");
        let mut scanned = backend.read_buffer_u32(output).expect("scan output");
        scanned.truncate(n);
        for handle in [input, output, block_sums, p, offsets] {
            backend.destroy_buffer(handle);
        }
        (scanned, sentinels)
    }

    #[test]
    fn folded_scan_matches_cpu_and_leaves_block_sums_untouched() {
        let mut backend = skip_if_no_gpu!();
        for n in [2600usize, 256 * 7 + 1, 256 * 11] {
            let data = lcg(n, 0x1234, 9);
            let expected = cpu_exclusive_scan(&data);
            let (scanned, sentinels) = gpu_scan(&mut backend, &data, TINY_LIMIT);
            assert_eq!(scanned, expected, "folded scan of {n} words");
            assert!(
                sentinels.iter().all(|&w| w == 0xDEAD_BEEF),
                "extra workgroups wrote past the last block sum (n = {n})"
            );
            let (unfolded, _) = gpu_scan(&mut backend, &data, REAL_LIMIT);
            assert_eq!(unfolded, expected, "unfolded scan of {n} words");
        }
    }

    // ── reductions ────────────────────────────────────────────────────────────

    fn gpu_block_reduce(
        backend: &mut WgpuBackendReal,
        entry: &str,
        data: &[f32],
        neutral: f32,
        limit: u32,
    ) -> (Vec<f32>, Vec<f32>) {
        let n = data.len();
        let blocks = n.div_ceil(BLOCK);
        let mut padded = data.to_vec();
        padded.resize(blocks * BLOCK, neutral);
        let input = upload_f32(backend, &padded);
        let mut init = vec![0.0f32; blocks];
        init.extend(std::iter::repeat_n(-12345.0f32, 16));
        let partials = upload_f32(backend, &init);
        let p = params(backend, [n as u32, 0, 0, 0]);
        backend
            .dispatch_wgsl(
                REDUCE_WGSL,
                entry,
                &[(input, ro()), (partials, rw()), (p, ro())],
                shape(blocks, limit),
            )
            .expect("reduce block");
        let back = backend.read_buffer_f32(partials).expect("partials");
        for handle in [input, partials, p] {
            backend.destroy_buffer(handle);
        }
        (back[..blocks].to_vec(), back[blocks..].to_vec())
    }

    #[test]
    fn folded_reductions_match_cpu_and_do_not_overrun_partials() {
        let mut backend = skip_if_no_gpu!();
        let data: Vec<f32> = lcg(2600, 77, 8).into_iter().map(|v| v as f32).collect();
        let blocks = data.len().div_ceil(BLOCK);
        let expected_sum: Vec<f32> = (0..blocks)
            .map(|b| {
                data[b * BLOCK..((b + 1) * BLOCK).min(data.len())]
                    .iter()
                    .sum()
            })
            .collect();
        let expected_max: Vec<f32> = (0..blocks)
            .map(|b| {
                data[b * BLOCK..((b + 1) * BLOCK).min(data.len())]
                    .iter()
                    .copied()
                    .fold(f32::NEG_INFINITY, f32::max)
            })
            .collect();
        for limit in [TINY_LIMIT, REAL_LIMIT] {
            let (sums, sum_tail) =
                gpu_block_reduce(&mut backend, "reduce_sum_block", &data, 0.0, limit);
            assert_eq!(sums, expected_sum, "block sums under limit {limit}");
            assert!(
                sum_tail.iter().all(|&v| v == -12345.0),
                "sum overran under {limit}"
            );
            let (maxima, max_tail) = gpu_block_reduce(
                &mut backend,
                "reduce_max_block",
                &data,
                -3.4028235e38,
                limit,
            );
            assert_eq!(maxima, expected_max, "block maxima under limit {limit}");
            assert!(
                max_tail.iter().all(|&v| v == -12345.0),
                "max overran under {limit}"
            );
        }
    }

    // ── histogram and compaction ──────────────────────────────────────────────

    #[test]
    fn folded_histogram_counts_every_element_once() {
        let mut backend = skip_if_no_gpu!();
        let bins = 16usize;
        let data = lcg(2600, 5, 40);
        let mut expected = vec![0u32; bins];
        for &v in &data {
            expected[(v as usize).min(bins - 1)] += 1;
        }
        for limit in [TINY_LIMIT, REAL_LIMIT] {
            let input = upload_u32(&mut backend, &data);
            let counts = zeros_u32(&mut backend, bins);
            let p = params(&mut backend, [data.len() as u32, bins as u32, 0, 0]);
            backend
                .dispatch_wgsl(
                    HISTOGRAM_WGSL,
                    "histogram_main",
                    &[(input, ro()), (counts, rw()), (p, ro())],
                    shape(data.len().div_ceil(BLOCK), limit),
                )
                .expect("histogram_main");
            assert_eq!(
                backend.read_buffer_u32(counts).expect("bins"),
                expected,
                "histogram under limit {limit}"
            );
        }
    }

    #[test]
    fn folded_sph_cell_histogram_counts_every_key_once() {
        let mut backend = skip_if_no_gpu!();
        let bins = 300usize;
        let keys = lcg(2600, 99, 300);
        let mut expected = vec![0u32; bins];
        for &k in &keys {
            expected[k as usize] += 1;
        }
        for limit in [TINY_LIMIT, REAL_LIMIT] {
            let input = upload_u32(&mut backend, &keys);
            let counts = zeros_u32(&mut backend, bins);
            let p = params(&mut backend, [keys.len() as u32, bins as u32, 0, 0]);
            backend
                .dispatch_wgsl(
                    SPH_CELL_HISTOGRAM_WGSL,
                    "sph_cell_histogram",
                    &[(input, ro()), (counts, rw()), (p, ro())],
                    shape(keys.len().div_ceil(BLOCK), limit),
                )
                .expect("sph_cell_histogram");
            assert_eq!(
                backend.read_buffer_u32(counts).expect("bins"),
                expected,
                "cell histogram under limit {limit}"
            );
        }
    }

    #[test]
    fn folded_compaction_keeps_the_stable_order() {
        let mut backend = skip_if_no_gpu!();
        let values = lcg(2600, 31, 1000);
        let keep: Vec<u32> = lcg(2600, 32, 3)
            .into_iter()
            .map(|v| u32::from(v == 0))
            .collect();
        let positions = cpu_exclusive_scan(&keep);
        let expected: Vec<u32> = values
            .iter()
            .zip(&keep)
            .filter(|&(_, &k)| k == 1)
            .map(|(&v, _)| v)
            .collect();
        for limit in [TINY_LIMIT, REAL_LIMIT] {
            let v = upload_u32(&mut backend, &values);
            let k = upload_u32(&mut backend, &keep);
            let pos = upload_u32(&mut backend, &positions);
            let out = zeros_u32(&mut backend, expected.len());
            let p = params(&mut backend, [values.len() as u32, 0, 0, 0]);
            backend
                .dispatch_wgsl(
                    COMPACT_WGSL,
                    "scatter_kept",
                    &[(v, ro()), (k, ro()), (pos, ro()), (out, rw()), (p, ro())],
                    shape(values.len().div_ceil(BLOCK), limit),
                )
                .expect("scatter_kept");
            assert_eq!(
                backend.read_buffer_u32(out).expect("compacted"),
                expected,
                "compaction under limit {limit}"
            );
        }
    }

    // ── radix sort ────────────────────────────────────────────────────────────

    /// Host-driven LSD radix sort of `(key, payload)` pairs with the shipped
    /// kernels; `with_payload == false` runs the key-only scatter instead.
    fn gpu_radix_sort(
        backend: &mut WgpuBackendReal,
        keys: &[u32],
        payload: &[u32],
        with_payload: bool,
        limit: u32,
    ) -> (Vec<u32>, Vec<u32>) {
        let n = keys.len();
        let tiles = n.div_ceil(BLOCK);
        let padded = tiles * BLOCK;
        let hist_len = 256 * tiles;
        let mut padded_keys = keys.to_vec();
        padded_keys.resize(padded, 0);
        let mut padded_pay = payload.to_vec();
        padded_pay.resize(padded, 0);
        let mut keys_a = upload_u32(backend, &padded_keys);
        let mut keys_b = zeros_u32(backend, padded);
        let mut pay_a = upload_u32(backend, &padded_pay);
        let mut pay_b = zeros_u32(backend, padded);
        let tile_hist = zeros_u32(backend, hist_len);
        let base = zeros_u32(backend, hist_len);
        let p = zeros_u32(backend, 4);
        let grid = shape(tiles, limit);
        for pass in 0u32..4 {
            backend.queue_write_buffer_raw(
                &p,
                bytemuck::cast_slice(&[n as u32, pass * 8, tiles as u32, 0u32]),
            );
            backend
                .dispatch_wgsl(
                    RADIX_HISTOGRAM_WGSL,
                    "radix_histogram",
                    &[(keys_a, ro()), (tile_hist, rw()), (p, ro())],
                    grid,
                )
                .expect("radix_histogram");
            let hist = backend.read_buffer_u32(tile_hist).expect("tile histogram");
            backend.queue_write_buffer_raw(&base, bytemuck::cast_slice(&cpu_exclusive_scan(&hist)));
            if with_payload {
                backend
                    .dispatch_wgsl(
                        RADIX_SCATTER_PAIRS_WGSL,
                        "radix_scatter_pairs",
                        &[
                            (keys_a, ro()),
                            (pay_a, ro()),
                            (keys_b, rw()),
                            (pay_b, rw()),
                            (base, ro()),
                            (p, ro()),
                        ],
                        grid,
                    )
                    .expect("radix_scatter_pairs");
            } else {
                backend
                    .dispatch_wgsl(
                        RADIX_SCATTER_WGSL,
                        "radix_scatter",
                        &[(keys_a, ro()), (keys_b, rw()), (base, ro()), (p, ro())],
                        grid,
                    )
                    .expect("radix_scatter");
            }
            std::mem::swap(&mut keys_a, &mut keys_b);
            std::mem::swap(&mut pay_a, &mut pay_b);
        }
        let mut sorted_keys = backend.read_buffer_u32(keys_a).expect("sorted keys");
        sorted_keys.truncate(n);
        let mut sorted_pay = if with_payload {
            backend.read_buffer_u32(pay_a).expect("sorted payload")
        } else {
            Vec::new()
        };
        sorted_pay.truncate(n);
        (sorted_keys, sorted_pay)
    }

    #[test]
    fn folded_radix_sort_is_stable_and_complete() {
        let mut backend = skip_if_no_gpu!();
        for n in [2600usize, 256 * 5 + 17] {
            // Many duplicate keys, so stability is observable.
            let keys: Vec<u32> = lcg(n, 4242, 1 << 20)
                .into_iter()
                .enumerate()
                .map(|(i, k)| {
                    if i % 2 == 0 {
                        k % 5000
                    } else {
                        k.wrapping_mul(2_654_435_761)
                    }
                })
                .collect();
            let payload: Vec<u32> = (0..n as u32).collect();
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by_key(|&i| keys[i]);
            let expected_keys: Vec<u32> = order.iter().map(|&i| keys[i]).collect();
            let expected_pay: Vec<u32> = order.iter().map(|&i| payload[i]).collect();
            for limit in [TINY_LIMIT, REAL_LIMIT] {
                let (k, pay) = gpu_radix_sort(&mut backend, &keys, &payload, true, limit);
                assert_eq!(k, expected_keys, "pair sort keys, n = {n}, limit {limit}");
                assert_eq!(
                    pay, expected_pay,
                    "pair sort payload, n = {n}, limit {limit}"
                );
                let (only_keys, _) = gpu_radix_sort(&mut backend, &keys, &payload, false, limit);
                assert_eq!(only_keys, expected_keys, "key sort, n = {n}, limit {limit}");
            }
        }
    }

    // ── SPH integration ───────────────────────────────────────────────────────

    #[test]
    fn folded_sph_integrate_updates_every_particle_exactly_once() {
        let mut backend = skip_if_no_gpu!();
        let n = 700usize;
        let positions: Vec<f32> = (0..n * 3).map(|i| (i / 3) as f32).collect();
        let dt = 0.5f32;
        for limit in [TINY_LIMIT, REAL_LIMIT] {
            let pos = upload_f32(&mut backend, &positions);
            let vel = upload_f32(&mut backend, &vec![0.0f32; n * 3]);
            let acc = upload_f32(&mut backend, &vec![1.0f32; n * 3]);
            let p = params(&mut backend, [n as u32, dt.to_bits(), 0, 0]);
            backend
                .dispatch_wgsl(
                    SPH_INTEGRATE_WGSL,
                    "sph_integrate",
                    &[(pos, rw()), (vel, rw()), (acc, ro()), (p, ro())],
                    shape(n.div_ceil(64), limit),
                )
                .expect("sph_integrate");
            let new_vel = backend.read_buffer_f32(vel).expect("velocities");
            let new_pos = backend.read_buffer_f32(pos).expect("positions");
            // v = a dt = 0.5, x = x0 + v dt = x0 + 0.25 -- a second pass would double both.
            assert!(
                new_vel.iter().take(n * 3).all(|&v| v == 0.5),
                "velocity updated more or less than once under limit {limit}"
            );
            for (i, (&x, &x0)) in new_pos.iter().zip(&positions).enumerate() {
                assert_eq!(x, x0 + 0.25, "position word {i} under limit {limit}");
            }
        }
    }
}
