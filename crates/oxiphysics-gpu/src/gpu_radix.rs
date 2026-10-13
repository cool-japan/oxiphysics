// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU LSD radix sort of `u32` keys and `(key, payload)` pairs.
//!
//! The sort runs four passes of eight bits each (least-significant digit first).
//! Each pass uses a decoupled, stable per-tile scheme: a histogram kernel builds
//! a digit-major per-tile count (`tile_hist[d * num_tiles + tile]`), the host
//! runs an exclusive scan over that array (it has just been read back, so the
//! scan costs no device work) to obtain exact per-(digit, tile) output bases,
//! and a scatter kernel writes each element to its base plus its stable local
//! rank within the tile. Two buffers per stream are ping-ponged across the four
//! passes (a second buffer pair carries the payload for the pairs variant); the
//! histogram, base and parameter buffers are allocated once per sort and reused
//! by every pass.
//!
//! The sorts run on the process-wide backend shared with [`crate::gpu_primitives`]:
//! the device and the compiled pipelines are reused by every call, and the scratch
//! buffers of a call are released before it returns.
//!
//! When no GPU adapter is available (or the `wgpu-backend` feature is disabled)
//! the public entry points transparently fall back to the CPU reference sort,
//! and so do inputs of more than 16 776 960 keys (a dispatch spans at most
//! 65 535 tiles of 256 keys).  The public signatures never change between
//! feature configurations.

/// GPU LSD radix sort of u32 keys (4x8-bit passes). CPU fallback when no adapter.
pub fn radix_sort_u32_gpu(keys: &[u32]) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        if keys.len() >= 2
            && let Some(result) =
                crate::gpu_primitives::shared::run(|gpu| gpu::radix_sort_u32(gpu, keys))
        {
            return result;
        }
    }
    cpu::radix_sort_u32(keys)
}

/// GPU LSD radix sort of (key, payload) u32 pairs, stable, payload follows key.
///
/// If `keys.len() != payload.len()`, the inputs are returned unchanged as
/// `(keys.to_vec(), payload.to_vec())`.
pub fn radix_sort_pairs_gpu(keys: &[u32], payload: &[u32]) -> (Vec<u32>, Vec<u32>) {
    if keys.len() != payload.len() {
        return (keys.to_vec(), payload.to_vec());
    }
    #[cfg(feature = "wgpu-backend")]
    {
        if keys.len() >= 2
            && let Some(result) =
                crate::gpu_primitives::shared::run(|gpu| gpu::radix_sort_pairs(gpu, keys, payload))
        {
            return result;
        }
    }
    cpu::cpu_radix_pairs(keys, payload)
}

// ── CPU reference (parity oracle) ─────────────────────────────────────────────

mod cpu {
    /// Keys-only CPU fallback using the crate's reference radix sort.
    pub(super) fn radix_sort_u32(keys: &[u32]) -> Vec<u32> {
        let mut v = keys.to_vec();
        crate::parallel_sort::radix_sort_u32(&mut v);
        v
    }

    /// Stable 4x8-bit counting-sort of (key, payload) pairs (parity oracle).
    pub(super) fn cpu_radix_pairs(keys: &[u32], payload: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let n = keys.len();
        if n == 0 {
            return (Vec::new(), Vec::new());
        }
        let mut cur_keys = keys.to_vec();
        let mut cur_pay = payload.to_vec();
        let mut nxt_keys = vec![0u32; n];
        let mut nxt_pay = vec![0u32; n];
        for pass_idx in 0u32..4 {
            let shift = pass_idx * 8;
            let mut count = [0u32; 256];
            for &k in &cur_keys {
                let digit = ((k >> shift) & 0xFF) as usize;
                count[digit] += 1;
            }
            let mut offsets = [0u32; 256];
            let mut acc = 0u32;
            for (off, &c) in offsets.iter_mut().zip(count.iter()) {
                *off = acc;
                acc += c;
            }
            for (&k, &p) in cur_keys.iter().zip(cur_pay.iter()) {
                let digit = ((k >> shift) & 0xFF) as usize;
                let dst = offsets[digit] as usize;
                nxt_keys[dst] = k;
                nxt_pay[dst] = p;
                offsets[digit] += 1;
            }
            std::mem::swap(&mut cur_keys, &mut nxt_keys);
            std::mem::swap(&mut cur_pay, &mut nxt_pay);
        }
        (cur_keys, cur_pay)
    }
}

// ── GPU implementation (real wgpu) ────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
mod gpu {
    use crate::compute::wgpu_backend::WgpuInitError;
    use crate::gpu_primitives::shared::{Scratch, check_workgroups};
    use crate::kernels_wgsl::{RADIX_HISTOGRAM_WGSL, RADIX_SCATTER_PAIRS_WGSL, RADIX_SCATTER_WGSL};

    const TILE: usize = 256;

    fn ro() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: false }
    }

    fn radix_params_bytes(n: u32, shift: u32, num_tiles: u32) -> Vec<u8> {
        let p: [u32; 4] = [n, shift, num_tiles, 0];
        bytemuck::cast_slice(&p).to_vec()
    }

    /// Stable LSD radix sort of `keys` in four 8-bit passes.
    ///
    /// Besides the two ping-pong key buffers the sort owns one histogram, one base
    /// and one parameter buffer, allocated once and reused by all four passes: the
    /// histogram kernel overwrites every `(digit, tile)` entry in every pass, so the
    /// buffer needs no clearing between passes, and the per-pass scan runs on the
    /// host over the histogram that has just been read back.
    pub(super) fn radix_sort_u32(
        gpu: &mut Scratch<'_>,
        keys: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = keys.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if n == 1 {
            return Ok(keys.to_vec());
        }
        let num_tiles = n.div_ceil(TILE);
        check_workgroups(num_tiles)?;
        let padded = num_tiles * TILE;
        let hist_len = 256 * num_tiles;

        let mut padded_keys = vec![0u32; padded];
        padded_keys[..n].copy_from_slice(keys);

        let mut buf_a = gpu.alloc((padded * 4) as u64);
        let mut buf_b = gpu.alloc((padded * 4) as u64);
        let tile_hist_buf = gpu.alloc((hist_len * 4) as u64);
        let base_buf = gpu.alloc((hist_len * 4) as u64);
        let params_buf = gpu.alloc(16);
        gpu.queue_write_buffer_raw(&buf_a, bytemuck::cast_slice(&padded_keys));

        for pass_idx in 0u32..4 {
            let shift = pass_idx * 8;
            gpu.queue_write_buffer_raw(
                &params_buf,
                &radix_params_bytes(n as u32, shift, num_tiles as u32),
            );

            gpu.dispatch_wgsl(
                RADIX_HISTOGRAM_WGSL,
                "radix_histogram",
                &[(buf_a, ro()), (tile_hist_buf, rw()), (params_buf, ro())],
                [num_tiles as u32, 1, 1],
            )?;

            let tile_hist = gpu.read_u32(tile_hist_buf, hist_len)?;
            let base = crate::gpu_primitives::cpu::exclusive_scan_u32(&tile_hist);
            gpu.queue_write_buffer_raw(&base_buf, bytemuck::cast_slice(&base));

            gpu.dispatch_wgsl(
                RADIX_SCATTER_WGSL,
                "radix_scatter",
                &[
                    (buf_a, ro()),
                    (buf_b, rw()),
                    (base_buf, ro()),
                    (params_buf, ro()),
                ],
                [num_tiles as u32, 1, 1],
            )?;

            std::mem::swap(&mut buf_a, &mut buf_b);
        }

        gpu.read_u32(buf_a, n)
    }

    /// Stable LSD radix sort of `(key, payload)` pairs; see [`radix_sort_u32`] for the
    /// buffer scheme.  `payload.len()` must equal `keys.len()`.
    pub(super) fn radix_sort_pairs(
        gpu: &mut Scratch<'_>,
        keys: &[u32],
        payload: &[u32],
    ) -> Result<(Vec<u32>, Vec<u32>), WgpuInitError> {
        let n = keys.len();
        if n == 0 {
            return Ok((Vec::new(), Vec::new()));
        }
        if n == 1 {
            let k = keys.first().copied().unwrap_or_default();
            let p = payload.first().copied().unwrap_or_default();
            return Ok((vec![k], vec![p]));
        }
        let num_tiles = n.div_ceil(TILE);
        check_workgroups(num_tiles)?;
        let padded = num_tiles * TILE;
        let hist_len = 256 * num_tiles;

        let mut padded_keys = vec![0u32; padded];
        padded_keys[..n].copy_from_slice(keys);
        let mut padded_pay = vec![0u32; padded];
        padded_pay[..n].copy_from_slice(payload);

        let mut keys_a = gpu.alloc((padded * 4) as u64);
        let mut keys_b = gpu.alloc((padded * 4) as u64);
        let mut pay_a = gpu.alloc((padded * 4) as u64);
        let mut pay_b = gpu.alloc((padded * 4) as u64);
        let tile_hist_buf = gpu.alloc((hist_len * 4) as u64);
        let base_buf = gpu.alloc((hist_len * 4) as u64);
        let params_buf = gpu.alloc(16);
        gpu.queue_write_buffer_raw(&keys_a, bytemuck::cast_slice(&padded_keys));
        gpu.queue_write_buffer_raw(&pay_a, bytemuck::cast_slice(&padded_pay));

        for pass_idx in 0u32..4 {
            let shift = pass_idx * 8;
            gpu.queue_write_buffer_raw(
                &params_buf,
                &radix_params_bytes(n as u32, shift, num_tiles as u32),
            );

            gpu.dispatch_wgsl(
                RADIX_HISTOGRAM_WGSL,
                "radix_histogram",
                &[(keys_a, ro()), (tile_hist_buf, rw()), (params_buf, ro())],
                [num_tiles as u32, 1, 1],
            )?;

            let tile_hist = gpu.read_u32(tile_hist_buf, hist_len)?;
            let base = crate::gpu_primitives::cpu::exclusive_scan_u32(&tile_hist);
            gpu.queue_write_buffer_raw(&base_buf, bytemuck::cast_slice(&base));

            gpu.dispatch_wgsl(
                RADIX_SCATTER_PAIRS_WGSL,
                "radix_scatter_pairs",
                &[
                    (keys_a, ro()),
                    (pay_a, ro()),
                    (keys_b, rw()),
                    (pay_b, rw()),
                    (base_buf, ro()),
                    (params_buf, ro()),
                ],
                [num_tiles as u32, 1, 1],
            )?;

            std::mem::swap(&mut keys_a, &mut keys_b);
            std::mem::swap(&mut pay_a, &mut pay_b);
        }

        let out_keys = gpu.read_u32(keys_a, n)?;
        let out_pay = gpu.read_u32(pay_a, n)?;
        Ok((out_keys, out_pay))
    }
}

// ── Tests (CPU oracle unit tests; GPU parity exercised via the public API) ─────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_radix_pairs_stable() {
        let (k, p) = cpu::cpu_radix_pairs(&[2, 1, 2, 1, 2], &[10, 20, 30, 40, 50]);
        assert_eq!(k, vec![1, 1, 2, 2, 2]);
        assert_eq!(p, vec![20, 40, 10, 30, 50]);
    }

    #[test]
    fn cpu_radix_pairs_matches_oracle() {
        // Deterministic LCG-generated keys; payload = original index.
        let n = 1000usize;
        let mut state = 0x1234_5678u32;
        let mut keys = Vec::with_capacity(n);
        for _ in 0..n {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            keys.push(state);
        }
        let payload: Vec<u32> = (0..n as u32).collect();
        let (sorted_keys, sorted_pay) = cpu::cpu_radix_pairs(&keys, &payload);

        // Oracle: sort an index array by key (stable).
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by_key(|&i| keys[i]);
        let oracle_keys: Vec<u32> = idx.iter().map(|&i| keys[i]).collect();
        let oracle_pay: Vec<u32> = idx.iter().map(|&i| i as u32).collect();

        assert_eq!(sorted_keys, oracle_keys);
        assert_eq!(sorted_pay, oracle_pay);
    }

    #[test]
    fn empty_and_single() {
        assert_eq!(radix_sort_u32_gpu(&[]), Vec::<u32>::new());
        assert_eq!(radix_sort_u32_gpu(&[42]), vec![42]);
        assert_eq!(
            radix_sort_pairs_gpu(&[], &[]),
            (Vec::<u32>::new(), Vec::<u32>::new())
        );
        assert_eq!(radix_sort_pairs_gpu(&[7], &[9]), (vec![7], vec![9]));
    }

    #[test]
    fn unequal_pairs_unchanged() {
        let (k, p) = radix_sort_pairs_gpu(&[3, 1], &[1, 2, 3]);
        assert_eq!(k, vec![3, 1]);
        assert_eq!(p, vec![1, 2, 3]);
    }

    #[test]
    fn scatter_local_rank_stable_tiny() {
        #[cfg(feature = "wgpu-backend")]
        {
            if crate::compute::wgpu_backend::real::WgpuBackendReal::try_new().is_err() {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            }
        }
        // Public API returns correct sorted output via GPU when present, else CPU.
        let sorted = radix_sort_u32_gpu(&[5, 3, 3, 9, 1, 3, 255, 0]);
        assert_eq!(sorted, vec![0, 1, 3, 3, 3, 5, 9, 255]);

        let (k, p) = radix_sort_pairs_gpu(&[2, 1, 2, 1, 2], &[10, 20, 30, 40, 50]);
        assert_eq!(k, vec![1, 1, 2, 2, 2]);
        assert_eq!(p, vec![20, 40, 10, 30, 50]);
    }

    /// Deterministic pseudo-random keys (numerical-recipes LCG), reduced modulo
    /// `modulus` unless it is zero.
    fn lcg_keys(n: usize, modulus: u32) -> Vec<u32> {
        let mut state = 0x1234_5678u32;
        (0..n)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                if modulus == 0 { state } else { state % modulus }
            })
            .collect()
    }

    /// Stable reference sort of `(key, payload)` pairs.
    fn reference_pairs(keys: &[u32], payload: &[u32]) -> (Vec<u32>, Vec<u32>) {
        let mut order: Vec<usize> = (0..keys.len()).collect();
        order.sort_by_key(|&i| keys[i]);
        (
            order.iter().map(|&i| keys[i]).collect(),
            order.iter().map(|&i| payload[i]).collect(),
        )
    }

    #[test]
    fn public_sorts_match_the_reference_across_tile_boundaries() {
        for &n in &[0usize, 1, 2, 3, 255, 256, 257, 511, 512, 513, 1000, 4097] {
            for &modulus in &[0u32, 7, 1000] {
                let keys = lcg_keys(n, modulus);
                let mut expected = keys.clone();
                expected.sort_unstable();
                assert_eq!(
                    radix_sort_u32_gpu(&keys),
                    expected,
                    "keys n={n} mod={modulus}"
                );

                let payload: Vec<u32> = (0..n as u32).map(|i| i.wrapping_mul(3) ^ 0xA5).collect();
                let (sorted_keys, sorted_payload) = radix_sort_pairs_gpu(&keys, &payload);
                let (ref_keys, ref_payload) = reference_pairs(&keys, &payload);
                assert_eq!(sorted_keys, ref_keys, "pair keys n={n} mod={modulus}");
                assert_eq!(
                    sorted_payload, ref_payload,
                    "pair payload (stability) n={n} mod={modulus}"
                );
            }
        }
    }

    #[test]
    fn high_byte_decides_the_order_of_otherwise_close_keys() {
        // Only the high byte differs: the last radix pass must decide the order.
        let keys = [
            0x0300_0001u32,
            0x0100_0002,
            0x0200_0003,
            0x0100_0001,
            0x0300_0000,
        ];
        let mut expected = keys.to_vec();
        expected.sort_unstable();
        assert_eq!(radix_sort_u32_gpu(&keys), expected);
        let payload = [0u32, 1, 2, 3, 4];
        let (sorted_keys, sorted_payload) = radix_sort_pairs_gpu(&keys, &payload);
        assert_eq!(sorted_keys, expected);
        assert_eq!(sorted_payload, vec![3, 1, 2, 4, 0]);
    }

    #[test]
    fn concurrent_sorts_are_independent_and_correct() {
        let workers: Vec<_> = (0..4u32)
            .map(|t| {
                std::thread::spawn(move || {
                    for round in 0..3u32 {
                        let n = 700 + 111 * (t as usize) + round as usize;
                        let keys = lcg_keys(n, 50 + t * 13);
                        let mut expected = keys.clone();
                        expected.sort_unstable();
                        assert_eq!(
                            radix_sort_u32_gpu(&keys),
                            expected,
                            "thread {t} round {round}"
                        );
                        let payload: Vec<u32> = (0..n as u32).collect();
                        assert_eq!(
                            radix_sort_pairs_gpu(&keys, &payload),
                            reference_pairs(&keys, &payload),
                            "thread {t} round {round}"
                        );
                    }
                })
            })
            .collect();
        for worker in workers {
            worker.join().expect("sorting thread panicked");
        }
    }

    #[cfg(feature = "wgpu-backend")]
    mod device_reuse {
        use super::*;
        use crate::gpu_primitives::shared;
        use std::sync::PoisonError;

        /// Keys per radix tile (the workgroup size of the radix kernels).
        const TILE_KEYS: usize = 256;

        #[test]
        fn pass_buffers_are_allocated_once_not_once_per_pass() {
            let keys = lcg_keys(5000, 0);
            let payload: Vec<u32> = (0..keys.len() as u32).collect();
            let measured = shared::run(|scratch| {
                let before = scratch.live_buffer_count();
                let sorted = gpu::radix_sort_u32(scratch, &keys)?;
                let key_buffers = scratch.live_buffer_count() - before;
                let pairs = gpu::radix_sort_pairs(scratch, &keys, &payload)?;
                let pair_buffers = scratch.live_buffer_count() - before - key_buffers;
                Ok((key_buffers, pair_buffers, sorted, pairs))
            });
            let Some((key_buffers, pair_buffers, sorted, pairs)) = measured else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            // Two key buffers + histogram + base + parameters, regardless of the four passes.
            assert_eq!(key_buffers, 5);
            // Two key and two payload buffers + histogram + base + parameters.
            assert_eq!(pair_buffers, 7);
            let mut expected = keys.clone();
            expected.sort_unstable();
            assert_eq!(sorted, expected);
            assert_eq!(pairs, reference_pairs(&keys, &payload));
        }

        #[test]
        fn radix_pipelines_are_compiled_once_and_reused_by_later_sorts() {
            // One scope holds the lock of the shared backend for the whole measurement,
            // so no other test thread can compile a pipeline in between.
            let measured = shared::run(|scratch| {
                let keys = lcg_keys(4000, 0);
                let payload: Vec<u32> = (0..keys.len() as u32).collect();
                // Warm-up: one sort of each kind compiles the histogram kernel and both
                // scatter kernels.
                gpu::radix_sort_u32(scratch, &keys)?;
                gpu::radix_sort_pairs(scratch, &keys, &payload)?;
                let warm = scratch.compiled_pipeline_count();
                for round in 1..=3usize {
                    let keys = lcg_keys(1000 * round + 7, 0);
                    let payload: Vec<u32> = (0..keys.len() as u32).collect();
                    let sorted = gpu::radix_sort_u32(scratch, &keys)?;
                    assert!(sorted.windows(2).all(|w| w[0] <= w[1]), "round {round}");
                    let pairs = gpu::radix_sort_pairs(scratch, &keys, &payload)?;
                    assert_eq!(pairs, reference_pairs(&keys, &payload), "round {round}");
                }
                Ok((warm, scratch.compiled_pipeline_count()))
            });
            let Some((warm, after_reuse)) = measured else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            // radix_histogram, radix_scatter and radix_scatter_pairs.
            assert!(warm >= 3, "expected at least three pipelines, got {warm}");
            assert_eq!(
                warm, after_reuse,
                "repeated sorts must reuse the cached pipelines instead of recompiling"
            );
        }

        #[test]
        fn sorts_wider_than_one_dispatch_dimension_are_refused_before_allocating() {
            if shared::shared_backend().is_none() {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            }
            // 65 536 tiles of 256 keys: one workgroup more than a dimension allows.
            let n = TILE_KEYS * (shared::MAX_WORKGROUPS_PER_DIMENSION + 1);
            let keys: Vec<u32> = (0..n as u32)
                .map(|i| i.wrapping_mul(2_654_435_761))
                .collect();
            let payload: Vec<u32> = (0..n as u32).collect();
            let refused = shared::run(|scratch| {
                let plain = gpu::radix_sort_u32(scratch, &keys).is_err();
                let pairs = gpu::radix_sort_pairs(scratch, &keys, &payload).is_err();
                Ok((plain, pairs, scratch.live_buffer_count()))
            });
            assert_eq!(
                refused,
                Some((true, true, 0)),
                "the oversized sorts must be refused before allocating anything"
            );
        }

        #[test]
        fn repeated_sorts_reuse_the_backend_and_release_their_buffers() {
            let Some(backend) = shared::shared_backend() else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            let runs_before = shared::gpu_runs();
            let keys = lcg_keys(3000, 0);
            let payload: Vec<u32> = (0..keys.len() as u32).collect();
            for _ in 0..3 {
                let _ = radix_sort_u32_gpu(&keys);
                let _ = radix_sort_pairs_gpu(&keys, &payload);
            }
            assert!(
                shared::gpu_runs() >= runs_before + 6,
                "the sorts must have run on the GPU"
            );
            assert_eq!(
                shared::init_attempts(),
                1,
                "one backend for the whole process"
            );
            let guard = backend.lock().unwrap_or_else(PoisonError::into_inner);
            assert_eq!(guard.live_buffer_count(), 0);
            assert_eq!(guard.live_buffer_bytes(), 0);
        }
    }
}
