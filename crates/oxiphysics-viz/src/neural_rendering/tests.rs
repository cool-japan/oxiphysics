// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the neural rendering primitives.

use super::*;

// ── PositionalEncoding ────────────────────────────────────────────────

#[test]
fn test_pe_scalar_output_dim() {
    let pe = PositionalEncoding::new(4, false);
    let enc = pe.encode_scalar(0.5);
    assert_eq!(enc.len(), 8); // 2 * 4 freqs, no input
}

#[test]
fn test_pe_scalar_with_input_dim() {
    let pe = PositionalEncoding::new(4, true);
    let enc = pe.encode_scalar(0.5);
    assert_eq!(enc.len(), 9); // 1 + 2*4
}

#[test]
fn test_pe_point_output_dim() {
    let pe = PositionalEncoding::new(10, true);
    let enc = pe.encode_point([0.1, 0.2, 0.3]);
    assert_eq!(enc.len(), pe.point_output_dim());
    assert_eq!(enc.len(), 3 * (1 + 2 * 10));
}

#[test]
fn test_pe_sin_cos_pairs() {
    let pe = PositionalEncoding::new(1, false);
    let enc = pe.encode_scalar(0.0);
    // freq=1, x=0: sin(0) = 0, cos(0) = 1
    assert!(enc[0].abs() < 1e-10, "sin(0) = 0");
    assert!((enc[1] - 1.0).abs() < 1e-10, "cos(0) = 1");
}

#[test]
fn test_pe_frequency_doubling() {
    let pe = PositionalEncoding::new(3, false);
    assert_eq!(pe.frequencies.len(), 3);
    assert!((pe.frequencies[0] - 1.0).abs() < 1e-10);
    assert!((pe.frequencies[1] - 2.0).abs() < 1e-10);
    assert!((pe.frequencies[2] - 4.0).abs() < 1e-10);
}

#[test]
fn test_pe_different_inputs_produce_different_encodings() {
    let pe = PositionalEncoding::new(6, true);
    let a = pe.encode_scalar(0.1);
    let b = pe.encode_scalar(0.9);
    let diff: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum();
    assert!(
        diff > 0.01,
        "different inputs should produce different encodings"
    );
}

// ── MlpLayer ─────────────────────────────────────────────────────────

#[test]
fn test_mlp_layer_output_dim() {
    let layer = MlpLayer::new(8, 4);
    let input = vec![0.5; 8];
    let out = layer.forward_relu(&input);
    assert_eq!(out.len(), 4);
}

#[test]
fn test_mlp_relu_non_negative() {
    let layer = MlpLayer::new(4, 8);
    let input = vec![-1.0, 1.0, -0.5, 0.5];
    let out = layer.forward_relu(&input);
    for &v in &out {
        assert!(v >= 0.0, "ReLU output must be non-negative, got {}", v);
    }
}

#[test]
fn test_mlp_sigmoid_range() {
    let layer = MlpLayer::new(4, 4);
    let input = vec![1.0, -1.0, 0.0, 2.0];
    let out = layer.forward_sigmoid(&input);
    for &v in &out {
        assert!(
            (0.0..=1.0).contains(&v),
            "sigmoid output must be in [0,1], got {}",
            v
        );
    }
}

// ── NeuralRadianceField ───────────────────────────────────────────────

#[test]
fn test_nerf_query_density_nonneg() {
    let nerf = NeuralRadianceField::new(4, 2, 32, 0.1, 10.0);
    let (density, _color) = nerf.query([0.1, 0.2, 0.3], [0.0, 0.0, 1.0]);
    assert!(
        density >= 0.0,
        "density must be non-negative, got {}",
        density
    );
}

#[test]
fn test_nerf_query_color_range() {
    let nerf = NeuralRadianceField::new(4, 2, 32, 0.1, 10.0);
    let (_density, color) = nerf.query([0.5, 0.5, 0.5], [1.0, 0.0, 0.0]);
    for &c in &color {
        assert!(
            (0.0..=1.0).contains(&c),
            "color must be in [0,1], got {}",
            c
        );
    }
}

#[test]
fn test_nerf_pos_enc_dim() {
    let nerf = NeuralRadianceField::new(10, 4, 64, 0.1, 10.0);
    assert_eq!(nerf.pos_enc_dim(), 3 * (1 + 2 * 10));
}

// ── NeRfVolume ────────────────────────────────────────────────────────

#[test]
fn test_nerf_volume_stratified_sample_count() {
    let vol = NeRfVolume::new(64, 128);
    let ts = vol.stratified_sample(0.1, 10.0, 64, false);
    assert_eq!(ts.len(), 64);
}

#[test]
fn test_nerf_volume_samples_in_range() {
    let vol = NeRfVolume::new(32, 0);
    let near = 0.5;
    let far = 8.0;
    let ts = vol.stratified_sample(near, far, 32, false);
    for &t in &ts {
        assert!(t >= near && t <= far, "t={} out of [{},{}]", t, near, far);
    }
}

#[test]
fn test_nerf_volume_render_ray_alpha_range() {
    let nerf = NeuralRadianceField::new(4, 2, 16, 0.1, 5.0);
    let vol = NeRfVolume::new(16, 0);
    let result = vol.render_ray(&nerf, [0.0; 3], [0.0, 0.0, 1.0]);
    assert!(
        result.alpha >= 0.0 && result.alpha <= 1.0,
        "alpha out of range: {}",
        result.alpha
    );
}

#[test]
fn test_nerf_volume_render_ray_color_range() {
    let nerf = NeuralRadianceField::new(4, 2, 16, 0.1, 5.0);
    let vol = NeRfVolume::new(16, 0);
    let result = vol.render_ray(&nerf, [0.0; 3], [0.0, 0.0, 1.0]);
    for &c in &result.color {
        assert!((0.0..=1.0).contains(&c), "color out of range: {}", c);
    }
}

#[test]
fn test_nerf_volume_integral_empty_samples() {
    let vol = NeRfVolume::new(8, 0);
    let result = vol.volume_render_integral(&[]);
    assert_eq!(result.num_samples, 0);
    assert_eq!(result.alpha, 0.0);
}

// ── InstantNgp ────────────────────────────────────────────────────────

#[test]
fn test_instant_ngp_encode_length() {
    let ngp = InstantNgp::new(8, 16, 1.5, 512);
    let enc = ngp.encode([0.5, 0.5, 0.5]);
    assert_eq!(enc.len(), 8 * 2); // 2 features per level
}

#[test]
fn test_instant_ngp_query_density_nonneg() {
    let ngp = InstantNgp::new(4, 8, 2.0, 256);
    let (density, _) = ngp.query([0.3, 0.4, 0.5]);
    assert!(density >= 0.0, "density must be non-negative: {}", density);
}

#[test]
fn test_instant_ngp_param_count_positive() {
    let ngp = InstantNgp::new(4, 8, 2.0, 256);
    let params = ngp.param_count();
    assert!(params > 0, "param count must be positive: {}", params);
}

#[test]
fn test_hash_grid_level_hash_deterministic() {
    let level = HashGridLevel::new(16, 512);
    let h1 = level.hash_coords(3, 7, 11);
    let h2 = level.hash_coords(3, 7, 11);
    assert_eq!(h1, h2, "hash must be deterministic");
}

#[test]
fn test_hash_grid_lookup_in_range() {
    let level = HashGridLevel::new(8, 128);
    let f = level.lookup(0.5, 0.5, 0.5);
    // Table initialized to zeros so result should be zero.
    assert_eq!(f, [0.0f32; 2]);
}

// ── GaussianSplatting ─────────────────────────────────────────────────

#[test]
fn test_gaussian_covariance_symmetric() {
    let g = Gaussian3D::new(
        [0.0; 3],
        [1.0, 0.5, 0.2],
        0.8,
        [0.1, 0.2, 0.3],
        [1.0, 0.0, 0.0, 0.0],
    );
    let cov = g.covariance();
    for (i, row) in cov.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(
                (val - cov[j][i]).abs() < 1e-10,
                "covariance not symmetric at ({},{}) : {} vs {}",
                i,
                j,
                val,
                cov[j][i]
            );
        }
    }
}

#[test]
fn test_gaussian_density_max_at_center() {
    let g = Gaussian3D::new([0.0; 3], [1.0; 3], 1.0, [1.0; 3], [1.0, 0.0, 0.0, 0.0]);
    let center_density = g.eval_density([0.0; 3]);
    let offset_density = g.eval_density([2.0, 0.0, 0.0]);
    assert!(
        center_density > offset_density,
        "center density {} must exceed offset density {}",
        center_density,
        offset_density
    );
}

/// Quaternion `[w, x, y, z]` for a rotation of `angle` radians about `axis`.
fn axis_angle_quaternion(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let n = len3(axis);
    let (s, c) = (0.5 * angle).sin_cos();
    [c, s * axis[0] / n, s * axis[1] / n, s * axis[2] / n]
}

/// Inverse of a 3x3 matrix by cofactor expansion (test helper).
fn invert3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    let inv_det = 1.0 / det;
    [
        [
            c00 * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            c01 * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            c02 * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

#[test]
fn test_gaussian_eval_density_respects_rotation() {
    // A 90 degree rotation about +z carries the long local x axis onto
    // world +y, so the density must be elongated along y, not x.
    let half = std::f64::consts::FRAC_PI_4;
    let g = Gaussian3D::new(
        [0.0; 3],
        [1.0; 3],
        1.0,
        [1.0, 0.1, 0.1],
        [half.cos(), 0.0, 0.0, half.sin()],
    );
    let along_long_axis = g.eval_density([0.0, 1.0, 0.0]);
    assert!(
        (along_long_axis - (-0.5_f64).exp()).abs() < 1e-8,
        "one sigma along the rotated long axis must give exp(-0.5), got {along_long_axis}"
    );
    let across = g.eval_density([1.0, 0.0, 0.0]);
    assert!(
        across < 1e-10,
        "ten sigma along the rotated short axis must vanish, got {across}"
    );
}

#[test]
fn test_gaussian_eval_density_matches_inverse_covariance() {
    // General orientation and anisotropic scale: the density must equal
    // exp(-0.5 d^T (Sigma + eps I)^-1 d) computed from covariance().
    let g = Gaussian3D::new(
        [0.0; 3],
        [1.0; 3],
        1.0,
        [1.5, 0.4, 0.8],
        axis_angle_quaternion([1.0, 2.0, 3.0], 0.7),
    );
    let mut cov = g.covariance();
    for (i, row) in cov.iter_mut().enumerate() {
        row[i] += 1e-10;
    }
    let inv = invert3(cov);
    for offset in [
        [0.3, -0.2, 0.9],
        [1.0, 1.0, -1.0],
        [-0.7, 0.05, 0.4],
        [0.0, 0.0, 0.0],
    ] {
        let mut md2 = 0.0;
        for (i, inv_row) in inv.iter().enumerate() {
            for (j, entry) in inv_row.iter().enumerate() {
                md2 += offset[i] * entry * offset[j];
            }
        }
        let expected = (-0.5 * md2).exp();
        let actual = g.eval_density(offset);
        assert!(
            (actual - expected).abs() < 1e-9,
            "offset {offset:?}: density {actual} vs covariance-inverse {expected}"
        );
    }
}

#[test]
fn test_gaussian_rotation_matrix_is_orthonormal_and_normalises_quaternion() {
    let q = axis_angle_quaternion([0.3, -0.5, 0.8], 1.1);
    let unit = Gaussian3D::new([0.0; 3], [1.0; 3], 1.0, [1.0; 3], q);
    let scaled = Gaussian3D::new(
        [0.0; 3],
        [1.0; 3],
        1.0,
        [1.0; 3],
        [3.0 * q[0], 3.0 * q[1], 3.0 * q[2], 3.0 * q[3]],
    );
    let r = unit.rotation_matrix();
    let r_scaled = scaled.rotation_matrix();
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (r[i][j] - r_scaled[i][j]).abs() < 1e-12,
                "scaled quaternion must give the same rotation at ({i},{j})"
            );
            let dot: f64 = (0..3).map(|k| r[k][i] * r[k][j]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (dot - expected).abs() < 1e-12,
                "R^T R must be the identity at ({i},{j}), got {dot}"
            );
        }
    }
    let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    assert!((det - 1.0).abs() < 1e-12, "det(R) must be +1, got {det}");
}

#[test]
fn test_gaussian_rotation_matrix_degenerate_quaternion_is_identity() {
    for q in [
        [0.0; 4],
        [f64::NAN, 0.0, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0, 0.0],
    ] {
        let g = Gaussian3D::new([0.0; 3], [1.0; 3], 1.0, [1.0; 3], q);
        let r = g.rotation_matrix();
        for (i, row) in r.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_eq!(value, expected, "quaternion {q:?} must map to identity");
            }
        }
    }
}

#[test]
fn test_gaussian_splatting_render_ray_uses_rotation() {
    // Elongated red Gaussian whose long axis points along world +y.
    let half = std::f64::consts::FRAC_PI_4;
    let mut scene = GaussianSplatting::new();
    scene.add(Gaussian3D::new(
        [0.0; 3],
        [1.0, 0.0, 0.0],
        1.0,
        [2.0, 0.05, 0.05],
        [half.cos(), 0.0, 0.0, half.sin()],
    ));
    // The ray passes one unit above the centre, i.e. half a sigma along
    // the long axis: alpha = exp(-0.5 * 0.25).
    let colour = scene.render_ray([0.0, 1.0, -5.0], [0.0, 0.0, 1.0]);
    let alpha = (-0.125_f64).exp();
    assert!(
        (colour[0] - 1.0).abs() < 1e-8,
        "red channel stays saturated, got {}",
        colour[0]
    );
    assert!(
        (colour[1] - (1.0 - alpha)).abs() < 1e-6,
        "green channel must show the white background leaking through ({}), got {}",
        1.0 - alpha,
        colour[1]
    );
    // A ray one unit to the side lies ~20 sigma across the short axis.
    let side = scene.render_ray([1.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
    assert!(
        (side[1] - 1.0).abs() < 1e-8,
        "ray across the short axis must miss the Gaussian, got {}",
        side[1]
    );
}

#[test]
fn test_gaussian_splatting_empty_scene_white() {
    let scene = GaussianSplatting::new();
    let color = scene.render_ray([0.0; 3], [0.0, 0.0, 1.0]);
    // Empty scene = pure white background
    for &c in &color {
        assert!(
            (c - 1.0).abs() < 1e-10,
            "empty scene should be white, got {}",
            c
        );
    }
}

#[test]
fn test_gaussian_splatting_bounding_sphere_empty() {
    let scene = GaussianSplatting::new();
    let (_center, radius) = scene.bounding_sphere();
    assert_eq!(radius, 0.0);
}

#[test]
fn test_gaussian_splatting_bounding_sphere_nonempty() {
    let mut scene = GaussianSplatting::new();
    scene.add(Gaussian3D::new(
        [1.0, 0.0, 0.0],
        [1.0; 3],
        0.5,
        [0.1; 3],
        [1.0, 0.0, 0.0, 0.0],
    ));
    scene.add(Gaussian3D::new(
        [-1.0, 0.0, 0.0],
        [0.5; 3],
        0.5,
        [0.1; 3],
        [1.0, 0.0, 0.0, 0.0],
    ));
    let (center, radius) = scene.bounding_sphere();
    assert!(radius > 0.0, "bounding radius must be positive");
    let _ = center; // used
}

// ── NeuralSdf ─────────────────────────────────────────────────────────

#[test]
fn test_neural_sdf_normal_unit_length() {
    let sdf = NeuralSdf::new(4, 32);
    let n = sdf.normal([0.5, 0.5, 0.5]);
    let mag = len3(n);
    // Normal may be zero-vector if SDF gradient is zero, otherwise unit.
    if mag > 1e-8 {
        assert!(
            (mag - 1.0).abs() < 1e-4,
            "normal magnitude should be ≈1, got {:.6}",
            mag
        );
    }
}

#[test]
fn test_neural_sdf_query_fields() {
    let sdf = NeuralSdf::new(4, 32);
    let res = sdf.query([0.0, 0.0, 0.0]);
    // Just verify that is_inside matches distance sign.
    assert_eq!(res.is_inside, res.distance < 0.0);
}

#[test]
fn test_neural_sdf_sphere_trace_or_none() {
    let sdf = NeuralSdf::new(2, 16);
    // Sphere trace may or may not find a surface; just ensure no panic.
    let _result = sdf.sphere_trace([0.0; 3], [1.0, 0.0, 0.0], 32, 0.01);
}

// ── DeepRenderImage ───────────────────────────────────────────────────

#[test]
fn test_deep_render_image_init() {
    let img = DeepRenderImage::new(16, 8);
    assert_eq!(img.width, 16);
    assert_eq!(img.height, 8);
    assert_eq!(img.pixels.len(), 128);
}

#[test]
fn test_deep_render_image_set_get_pixel() {
    let mut img = DeepRenderImage::new(4, 4);
    img.set_pixel(2, 3, [0.5, 0.25, 0.75], 1.0, 0.9);
    let px = img.get_pixel(2, 3);
    assert!((px[0] - 0.5).abs() < 1e-6);
    assert!((px[1] - 0.25).abs() < 1e-6);
    assert!((px[2] - 0.75).abs() < 1e-6);
}

#[test]
fn test_deep_render_image_rgba8_length() {
    let img = DeepRenderImage::new(4, 4);
    let bytes = img.to_rgba8();
    assert_eq!(bytes.len(), 4 * 4 * 4); // 4 channels per pixel
}

#[test]
fn test_deep_render_image_mean_luminance_black() {
    let img = DeepRenderImage::new(8, 8);
    assert!(
        img.mean_luminance().abs() < 1e-6,
        "all-black image should have zero luminance"
    );
}

// ── NeuralDenoiser ────────────────────────────────────────────────────

#[test]
fn test_neural_denoiser_kernel_weights_sum_to_one() {
    let d = NeuralDenoiser::new(DenoiserKernel::Gaussian, 3);
    let sum: f64 = d.kernel_weights.iter().sum();
    assert!(
        (sum - 1.0).abs() < 1e-10,
        "kernel weights should sum to 1, got {:.6}",
        sum
    );
}

#[test]
fn test_neural_denoiser_output_same_size() {
    let d = NeuralDenoiser::new(DenoiserKernel::ATrous, 2);
    let img = DeepRenderImage::new(8, 8);
    let out = d.denoise(&img);
    assert_eq!(out.width, img.width);
    assert_eq!(out.height, img.height);
}

// ── StyleTransferViz ─────────────────────────────────────────────────

#[test]
fn test_style_transfer_gram_diagonal_nonneg() {
    let st = StyleTransferViz::new(1.0, 0.01, 100);
    let features = vec![vec![1.0, 0.0, 0.5], vec![0.0, 1.0, 0.5]];
    let gram = st.compute_gram(&features);
    assert_eq!(gram.len(), 4); // 2x2
    // Diagonal elements (self dot products) should be non-negative.
    assert!(gram[0] >= 0.0 && gram[3] >= 0.0);
}

#[test]
fn test_style_transfer_loss_zero_for_equal_grams() {
    let st = StyleTransferViz::new(1.0, 1.0, 10);
    let gram = vec![1.0, 0.0, 0.0, 1.0];
    let loss = st.style_loss(&gram, &gram);
    assert!(
        loss.abs() < 1e-10,
        "identical grams should have zero style loss"
    );
}

// ── LatentSpaceViz ────────────────────────────────────────────────────

#[test]
fn test_latent_space_viz_add_point() {
    let mut viz = LatentSpaceViz::new(30.0, 100, 5);
    viz.add_point(vec![1.0, 0.0, 0.0], 0);
    viz.add_point(vec![0.0, 1.0, 0.0], 1);
    assert_eq!(viz.points.len(), 2);
}

#[test]
fn test_latent_space_viz_embedding_distance_zero_self() {
    let viz = LatentSpaceViz::new(10.0, 50, 3);
    let emb = vec![1.0, 2.0, 3.0];
    let d = viz.embedding_distance(&emb, &emb);
    assert!(d.abs() < 1e-10, "self-distance should be zero, got {}", d);
}

#[test]
fn test_latent_space_viz_distance_matrix_diagonal_zero() {
    let mut viz = LatentSpaceViz::new(10.0, 50, 3);
    viz.add_point(vec![1.0, 0.0], 0);
    viz.add_point(vec![0.0, 1.0], 1);
    let dm = viz.distance_matrix();
    assert!(dm[0][0].abs() < 1e-10, "diagonal should be zero");
    assert!(dm[1][1].abs() < 1e-10, "diagonal should be zero");
}

#[test]
fn test_latent_space_viz_cluster_color_range() {
    let viz = LatentSpaceViz::new(10.0, 50, 8);
    for i in 0..8 {
        let c = viz.cluster_color(i);
        for &v in &c {
            assert!(
                (0.0..=1.0).contains(&v),
                "color component {} out of [0,1]",
                v
            );
        }
    }
}

#[test]
fn test_latent_space_viz_pca_projection() {
    let mut viz = LatentSpaceViz::new(10.0, 50, 3);
    viz.add_point(vec![1.0, 2.0, 3.0], 0);
    viz.add_point(vec![4.0, 5.0, 6.0], 1);
    viz.project_pca();
    // After centering, mean should be near zero.
    let mx: f64 = viz.points.iter().map(|p| p.projected[0]).sum::<f64>() / 2.0;
    let my: f64 = viz.points.iter().map(|p| p.projected[1]).sum::<f64>() / 2.0;
    assert!(mx.abs() < 1e-10, "projected x-mean should be 0, got {}", mx);
    assert!(my.abs() < 1e-10, "projected y-mean should be 0, got {}", my);
}

#[test]
fn test_hsl_to_rgb_pure_red() {
    let (r, g, b) = hsl_to_rgb(0.0, 1.0, 0.5);
    assert!((r - 1.0).abs() < 1e-6, "r should be 1 for hue=0");
    assert!(g.abs() < 1e-6, "g should be 0 for hue=0");
    assert!(b.abs() < 1e-6, "b should be 0 for hue=0");
}

#[test]
fn test_nerf_volume_stratified_monotone() {
    let vol = NeRfVolume::new(16, 0);
    let ts = vol.stratified_sample(0.1, 10.0, 16, false);
    for i in 1..ts.len() {
        assert!(
            ts[i] > ts[i - 1],
            "samples must be monotonically increasing"
        );
    }
}
