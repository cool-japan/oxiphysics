// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU/CPU parity tests for the GPU LSD radix sort (`gpu_radix`).
//!
//! Skip-not-fail on headless CI: if no adapter is available the GPU-gated tests
//! print "SKIPPED: no GPU adapter available" and return.

/// Always-available smoke test so the binary links even without the feature.
#[test]
fn gpu_radix_link_smoke() {
    assert_eq!(1 + 1, 2);
}

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
    use oxiphysics_gpu::gpu_radix;

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

    /// Deterministic pseudo-random key generator (numerical-recipes LCG).
    fn lcg_keys(n: usize) -> Vec<u32> {
        let mut keys = Vec::with_capacity(n);
        let mut state: u32 = 0x1234_5678;
        for _ in 0..n {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            keys.push(state);
        }
        keys
    }

    /// CPU oracle: ascending sort (radix sort ascending == sort_unstable).
    fn cpu_sorted(keys: &[u32]) -> Vec<u32> {
        let mut v = keys.to_vec();
        v.sort_unstable();
        v
    }

    #[test]
    fn test_radix_keys_parity() {
        let _backend = skip_if_no_gpu!();
        let mut sizes = vec![255usize, 256, 257, 1024, 65536, 1_000_000];
        if std::env::var("OXIPHYSICS_GPU_BENCH").is_ok() {
            sizes.push(16_777_216);
        }
        for n in sizes {
            let keys = lcg_keys(n);
            let expected = cpu_sorted(&keys);
            let got = gpu_radix::radix_sort_u32_gpu(&keys);
            assert_eq!(got.len(), expected.len(), "length mismatch at n={n}");
            assert_eq!(got.first(), expected.first(), "first mismatch at n={n}");
            let mid = n / 2;
            assert_eq!(got.get(mid), expected.get(mid), "mid mismatch at n={n}");
            assert_eq!(got.last(), expected.last(), "last mismatch at n={n}");
            assert_eq!(got, expected, "full-vec mismatch at n={n}");
        }
    }

    #[test]
    fn test_radix_duplicate_heavy() {
        let _backend = skip_if_no_gpu!();
        let n = 65536usize;
        // Keys in a tiny range -> many duplicate digits per tile.
        let mut state: u32 = 0x0BAD_F00D;
        let mut keys = Vec::with_capacity(n);
        for _ in 0..n {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            keys.push(state % 17);
        }
        let expected = cpu_sorted(&keys);
        let got = gpu_radix::radix_sort_u32_gpu(&keys);
        assert_eq!(got.len(), expected.len(), "length mismatch (dup heavy)");
        assert_eq!(got.first(), expected.first(), "first mismatch (dup heavy)");
        assert_eq!(got.last(), expected.last(), "last mismatch (dup heavy)");
        assert_eq!(got, expected, "full-vec mismatch (dup heavy)");
    }

    #[test]
    fn test_radix_reverse_sorted() {
        let _backend = skip_if_no_gpu!();
        let n = 100_000usize;

        // Strictly descending input -> ascending 0..n.
        let descending: Vec<u32> = (0..n).map(|i| (n - 1 - i) as u32).collect();
        let ascending: Vec<u32> = (0..n as u32).collect();
        let got = gpu_radix::radix_sort_u32_gpu(&descending);
        assert_eq!(got.len(), ascending.len(), "length mismatch (reverse)");
        assert_eq!(got.first(), ascending.first(), "first mismatch (reverse)");
        assert_eq!(got.last(), ascending.last(), "last mismatch (reverse)");
        assert_eq!(got, ascending, "full-vec mismatch (reverse)");

        // Already-ascending input -> itself.
        let got_asc = gpu_radix::radix_sort_u32_gpu(&ascending);
        assert_eq!(got_asc, ascending, "ascending input should return itself");
    }

    #[test]
    fn test_radix_pairs_keeps_payload_aligned() {
        let _backend = skip_if_no_gpu!();
        let n = 50_000usize;
        let keys = lcg_keys(n);
        let payload: Vec<u32> = (0..n as u32).collect();

        let (keys_out, payload_out) = gpu_radix::radix_sort_pairs_gpu(&keys, &payload);

        // (a) Returned keys equal the sorted keys.
        let expected_keys = cpu_sorted(&keys);
        assert_eq!(keys_out, expected_keys, "pairs keys not sorted");

        // (b) Each output key is the key originally at its carried payload index,
        //     proving the payload stayed aligned to its key across all 4 passes.
        assert_eq!(keys_out.len(), n, "keys_out length mismatch");
        assert_eq!(payload_out.len(), n, "payload_out length mismatch");
        for j in 0..n {
            assert_eq!(
                keys_out[j],
                keys[payload_out[j] as usize],
                "payload misaligned at output position {j}"
            );
        }

        // (c) Payload is a permutation of 0..n.
        let mut payload_sorted = payload_out.clone();
        payload_sorted.sort_unstable();
        let identity: Vec<u32> = (0..n as u32).collect();
        assert_eq!(
            payload_sorted, identity,
            "payload is not a permutation of 0..n"
        );
    }

    #[test]
    fn test_radix_pairs_stable_duplicates() {
        let _backend = skip_if_no_gpu!();
        let keys = [2u32, 1, 2, 1, 2, 1];
        let payload = [0u32, 1, 2, 3, 4, 5];
        let (keys_out, payload_out) = gpu_radix::radix_sort_pairs_gpu(&keys, &payload);
        assert_eq!(keys_out, vec![1u32, 1, 1, 2, 2, 2], "stable keys mismatch");
        assert_eq!(
            payload_out,
            vec![1u32, 3, 5, 0, 2, 4],
            "stable payload order mismatch"
        );
    }

    #[test]
    fn test_radix_small_edge() {
        let _backend = skip_if_no_gpu!();
        assert_eq!(
            gpu_radix::radix_sort_u32_gpu(&[]),
            Vec::<u32>::new(),
            "empty slice should yield empty"
        );
        assert_eq!(
            gpu_radix::radix_sort_u32_gpu(&[7]),
            vec![7u32],
            "single element should be unchanged"
        );
        assert_eq!(
            gpu_radix::radix_sort_u32_gpu(&[2, 1]),
            vec![1u32, 2],
            "two elements should be sorted"
        );
    }
}
