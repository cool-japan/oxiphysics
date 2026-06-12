// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU LSD radix sort of `u32` keys and `(key, payload)` pairs.
//!
//! The sort runs four passes of eight bits each (least-significant digit first).
//! Each pass uses a decoupled, stable per-tile scheme: a histogram kernel builds
//! a digit-major per-tile count (`tile_hist[d * num_tiles + tile]`), the host
//! runs an exclusive scan over that array via [`crate::gpu_primitives::exclusive_scan_u32`]
//! to obtain exact per-(digit, tile) output bases, and a scatter kernel writes
//! each element to its base plus its stable local rank within the tile. Two
//! buffers per stream are ping-ponged across the four passes (a second buffer
//! pair carries the payload for the pairs variant).
//!
//! When no GPU adapter is available (or the `wgpu-backend` feature is disabled)
//! the public entry points transparently fall back to the CPU reference sort.
//! The public signatures never change between feature configurations.

/// GPU LSD radix sort of u32 keys (4x8-bit passes). CPU fallback when no adapter.
pub fn radix_sort_u32_gpu(keys: &[u32]) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) = crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::radix_sort_u32(&mut backend, keys)
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
        if let Ok(mut backend) = crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::radix_sort_pairs(&mut backend, keys, payload)
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
    use crate::compute::wgpu_backend::real::WgpuBackendReal;
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

    pub(super) fn radix_sort_u32(
        backend: &mut WgpuBackendReal,
        keys: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = keys.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if n == 1 {
            return Ok(keys.first().copied().map(|k| vec![k]).unwrap_or_default());
        }
        let num_tiles = n.div_ceil(TILE);
        let padded = num_tiles * TILE;
        let hist_len = 256 * num_tiles;

        let mut padded_keys = vec![0u32; padded];
        padded_keys[..n].copy_from_slice(keys);

        let mut buf_a = backend.create_buffer_storage((padded * 4) as u64);
        let mut buf_b = backend.create_buffer_storage((padded * 4) as u64);
        backend.queue_write_buffer_raw(&buf_a, bytemuck::cast_slice(&padded_keys));

        for pass_idx in 0u32..4 {
            let shift = pass_idx * 8;

            let tile_hist_buf = backend.create_buffer_storage((hist_len * 4) as u64);
            let zero_hist = vec![0u32; hist_len];
            backend.queue_write_buffer_raw(&tile_hist_buf, bytemuck::cast_slice(&zero_hist));

            let params_buf = backend.create_buffer_storage(16);
            backend.queue_write_buffer_raw(
                &params_buf,
                &radix_params_bytes(n as u32, shift, num_tiles as u32),
            );

            backend.dispatch_wgsl(
                RADIX_HISTOGRAM_WGSL,
                "radix_histogram",
                &[(buf_a, ro()), (tile_hist_buf, rw()), (params_buf, ro())],
                [num_tiles as u32, 1, 1],
            )?;

            let mut tile_hist = backend.read_buffer_u32(tile_hist_buf);
            tile_hist.truncate(hist_len);
            let base = crate::gpu_primitives::exclusive_scan_u32(&tile_hist);

            let base_buf = backend.create_buffer_storage((hist_len * 4) as u64);
            backend.queue_write_buffer_raw(&base_buf, bytemuck::cast_slice(&base));

            backend.dispatch_wgsl(
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

        let mut out = backend.read_buffer_u32(buf_a);
        out.truncate(n);
        Ok(out)
    }

    pub(super) fn radix_sort_pairs(
        backend: &mut WgpuBackendReal,
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
        let padded = num_tiles * TILE;
        let hist_len = 256 * num_tiles;

        let mut padded_keys = vec![0u32; padded];
        padded_keys[..n].copy_from_slice(keys);
        let mut padded_pay = vec![0u32; padded];
        padded_pay[..n].copy_from_slice(payload);

        let mut keys_a = backend.create_buffer_storage((padded * 4) as u64);
        let mut keys_b = backend.create_buffer_storage((padded * 4) as u64);
        let mut pay_a = backend.create_buffer_storage((padded * 4) as u64);
        let mut pay_b = backend.create_buffer_storage((padded * 4) as u64);
        backend.queue_write_buffer_raw(&keys_a, bytemuck::cast_slice(&padded_keys));
        backend.queue_write_buffer_raw(&pay_a, bytemuck::cast_slice(&padded_pay));

        for pass_idx in 0u32..4 {
            let shift = pass_idx * 8;

            let tile_hist_buf = backend.create_buffer_storage((hist_len * 4) as u64);
            let zero_hist = vec![0u32; hist_len];
            backend.queue_write_buffer_raw(&tile_hist_buf, bytemuck::cast_slice(&zero_hist));

            let params_buf = backend.create_buffer_storage(16);
            backend.queue_write_buffer_raw(
                &params_buf,
                &radix_params_bytes(n as u32, shift, num_tiles as u32),
            );

            backend.dispatch_wgsl(
                RADIX_HISTOGRAM_WGSL,
                "radix_histogram",
                &[(keys_a, ro()), (tile_hist_buf, rw()), (params_buf, ro())],
                [num_tiles as u32, 1, 1],
            )?;

            let mut tile_hist = backend.read_buffer_u32(tile_hist_buf);
            tile_hist.truncate(hist_len);
            let base = crate::gpu_primitives::exclusive_scan_u32(&tile_hist);

            let base_buf = backend.create_buffer_storage((hist_len * 4) as u64);
            backend.queue_write_buffer_raw(&base_buf, bytemuck::cast_slice(&base));

            backend.dispatch_wgsl(
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

        let mut out_keys = backend.read_buffer_u32(keys_a);
        out_keys.truncate(n);
        let mut out_pay = backend.read_buffer_u32(pay_a);
        out_pay.truncate(n);
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
}
