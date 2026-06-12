// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU/CPU parity tests for the reduction-primitive library.
//!
//! Skip-not-fail on headless CI: if no adapter is available the GPU-gated tests
//! print "SKIPPED: no GPU adapter available" and return.

/// Always-available smoke test so the binary links even without the feature.
#[test]
fn gpu_primitives_link_smoke() {
    assert_eq!(1 + 1, 2);
}

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
    use oxiphysics_gpu::gpu_primitives;

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

    fn cpu_exclusive_scan(data: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(data.len());
        let mut acc: u32 = 0;
        for &v in data {
            out.push(acc);
            acc = acc.wrapping_add(v);
        }
        out
    }

    #[test]
    fn test_exclusive_scan_parity() {
        let _backend = skip_if_no_gpu!();
        let mut sizes = vec![1usize, 255, 256, 257, 1024, 65536, 1_000_000];
        if std::env::var("OXIPHYSICS_GPU_BENCH").is_ok() {
            sizes.push(16_777_216);
        }
        for n in sizes {
            let data: Vec<u32> = (0..n).map(|i| (i % 7) as u32).collect();
            let expected = cpu_exclusive_scan(&data);
            let got = gpu_primitives::exclusive_scan_u32(&data);
            assert_eq!(got.len(), expected.len(), "length mismatch at n={n}");
            assert_eq!(got.first(), expected.first(), "first mismatch at n={n}");
            let mid = n / 2;
            assert_eq!(got.get(mid), expected.get(mid), "mid mismatch at n={n}");
            assert_eq!(got.last(), expected.last(), "last mismatch at n={n}");
            assert_eq!(got, expected, "full-vec mismatch at n={n}");
        }
    }

    #[test]
    fn test_exclusive_scan_empty() {
        let _backend = skip_if_no_gpu!();
        assert_eq!(gpu_primitives::exclusive_scan_u32(&[]), Vec::<u32>::new());
    }

    #[test]
    fn test_reduce_sum_parity() {
        let _backend = skip_if_no_gpu!();
        for n in [1usize, 256, 257, 65536, 1_000_000] {
            let data: Vec<f32> = (0..n).map(|i| (i % 13) as f32 * 0.5).collect();
            let mut cpu_sum = 0f32;
            for &v in &data {
                cpu_sum += v;
            }
            let gpu_sum = gpu_primitives::reduce_sum_f32(&data);
            let tol = (n as f32) * f32::EPSILON * cpu_sum.abs().max(1.0) * 4.0;
            assert!(
                (gpu_sum - cpu_sum).abs() <= tol,
                "sum mismatch at n={n}: gpu={gpu_sum} cpu={cpu_sum} tol={tol}"
            );
        }
    }

    #[test]
    fn test_reduce_max_parity() {
        let _backend = skip_if_no_gpu!();
        for n in [1usize, 256, 257, 1000, 65536] {
            let mut data: Vec<f32> = (0..n).map(|i| (i % 50) as f32).collect();
            let max_idx = n / 2;
            data[max_idx] = 9999.5;
            let gpu_max = gpu_primitives::reduce_max_f32(&data);
            assert_eq!(gpu_max, 9999.5, "max mismatch at n={n}");
        }
    }

    #[test]
    fn test_reduce_empty() {
        let _backend = skip_if_no_gpu!();
        assert_eq!(gpu_primitives::reduce_sum_f32(&[]), 0.0);
        assert_eq!(gpu_primitives::reduce_max_f32(&[]), f32::NEG_INFINITY);
    }

    #[test]
    fn test_histogram_parity() {
        let _backend = skip_if_no_gpu!();
        let n = 50_000usize;
        let data: Vec<u32> = (0..n)
            .map(|i| (i as u64 * 2654435761u64 % 200) as u32)
            .collect();
        let num_bins = 200usize;
        let mut cpu_bins = vec![0u32; num_bins];
        for &v in &data {
            let idx = (v as usize).min(num_bins - 1);
            cpu_bins[idx] += 1;
        }
        let gpu_bins = gpu_primitives::histogram_u32(&data, num_bins);
        assert_eq!(gpu_bins, cpu_bins, "histogram parity mismatch");

        let small = [0u32, 1, 1, 2, 2, 2, 5];
        assert_eq!(
            gpu_primitives::histogram_u32(&small, 4),
            vec![1u32, 2, 3, 1],
            "small histogram clamp mismatch"
        );
    }

    #[test]
    fn test_compact_parity() {
        let _backend = skip_if_no_gpu!();
        let n = 10_000usize;
        let values: Vec<u32> = (0..n).map(|i| i as u32).collect();
        let keep: Vec<u32> = (0..n).map(|i| if i % 3 == 0 { 1u32 } else { 0 }).collect();
        let cpu: Vec<u32> = values
            .iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k != 0)
            .map(|(&v, _)| v)
            .collect();
        let gpu = gpu_primitives::compact_u32(&values, &keep);
        assert_eq!(gpu, cpu, "compact parity mismatch");

        let all_zero = vec![0u32; 100];
        let vals = (0..100u32).collect::<Vec<_>>();
        assert!(gpu_primitives::compact_u32(&vals, &all_zero).is_empty());

        let all_one = vec![1u32; 100];
        assert_eq!(gpu_primitives::compact_u32(&vals, &all_one), vals);
    }
}
