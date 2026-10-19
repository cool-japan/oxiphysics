// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the medical marching-cubes extractor
//! ([`marching_cubes`]): topology (welded, watertight), orientation (outward
//! normals, positive signed volume), accuracy and handling of anisotropic voxels.

use super::*;
use crate::mc_tables::test_support::{
    area_and_volume, edge_valence_histogram, is_closed_and_consistently_oriented,
};
use crate::mc_tables::{MC_EDGE_TABLE, case_triangles};
use std::f64::consts::PI;

/// Volume with `value(world position)` sampled at every voxel centre.
fn sampled_volume(
    dims: [usize; 3],
    voxel_size: [f64; 3],
    origin: Vec3,
    value: impl Fn(Vec3) -> f64,
) -> ScalarVolume {
    let mut volume = ScalarVolume::new(dims[0], dims[1], dims[2], voxel_size, origin);
    for iz in 0..dims[2] {
        for iy in 0..dims[1] {
            for ix in 0..dims[0] {
                let position = Vec3::new(
                    origin.x + ix as f64 * voxel_size[0],
                    origin.y + iy as f64 * voxel_size[1],
                    origin.z + iz as f64 * voxel_size[2],
                );
                volume.set(ix, iy, iz, value(position));
            }
        }
    }
    volume
}

/// Mesh as plain arrays for the shared mesh-analysis helpers.
fn as_arrays(mesh: &IsoSurface) -> Vec<[f64; 3]> {
    mesh.vertices.iter().map(|v| [v.x, v.y, v.z]).collect()
}

/// Signed volume (divergence theorem); positive when the triangles face outwards.
fn signed_volume(mesh: &IsoSurface) -> f64 {
    area_and_volume(&as_arrays(mesh), &mesh.triangles).1
}

fn assert_watertight(mesh: &IsoSurface) {
    assert_eq!(
        edge_valence_histogram(&mesh.triangles),
        vec![(2, mesh.triangles.len() * 3 / 2)],
        "every mesh edge must belong to exactly two triangles"
    );
    assert!(is_closed_and_consistently_oriented(&mesh.triangles));
}

/// Binary ball of radius `radius` mm: 1000 inside, -1000 outside.
fn binary_ball(center: Vec3, radius: f64) -> impl Fn(Vec3) -> f64 {
    move |p| {
        if p.distance_to(&center) < radius {
            1000.0
        } else {
            -1000.0
        }
    }
}

#[test]
fn binary_sphere_of_radius_ten_voxels_is_watertight_and_within_three_percent() {
    let center = Vec3::new(15.5, 15.5, 15.5);
    let volume = sampled_volume(
        [32, 32, 32],
        [1.0, 1.0, 1.0],
        Vec3::zero(),
        binary_ball(center, 10.0),
    );
    let mesh = marching_cubes(&volume, 0.0);
    assert!(
        mesh.triangle_count() > 1000,
        "{} triangles",
        mesh.triangle_count()
    );
    assert_watertight(&mesh);
    let exact = 4.0 / 3.0 * PI * 10.0_f64.powi(3);
    let volume_mm3 = signed_volume(&mesh);
    assert!(
        volume_mm3 > 0.0,
        "triangles must face outwards: {volume_mm3}"
    );
    assert!(
        (volume_mm3 / exact - 1.0).abs() < 0.03,
        "volume {volume_mm3} vs {exact}"
    );
    // `mesh_volume_mm3` reports the same magnitude.
    assert!((mesh_volume_mm3(&mesh) - volume_mm3).abs() < 1e-9 * exact);
}

#[test]
fn smooth_sphere_field_is_accurate_and_vertices_lie_on_the_sphere() {
    let center = Vec3::new(15.5, 15.5, 15.5);
    let radius = 10.0;
    let volume = sampled_volume([32, 32, 32], [1.0; 3], Vec3::zero(), |p| {
        radius - p.distance_to(&center)
    });
    let mesh = marching_cubes(&volume, 0.0);
    assert_watertight(&mesh);
    let exact = 4.0 / 3.0 * PI * radius.powi(3);
    let measured = signed_volume(&mesh);
    assert!(
        (measured / exact - 1.0).abs() < 0.01,
        "volume {measured} vs {exact}"
    );
    let exact_area = 4.0 * PI * radius * radius;
    assert!(
        (mesh.surface_area() / exact_area - 1.0).abs() < 0.02,
        "area {} vs {exact_area}",
        mesh.surface_area()
    );
    for v in &mesh.vertices {
        let r = v.distance_to(&center);
        assert!((r - radius).abs() < 0.06, "vertex radius {r}");
    }
}

#[test]
fn every_surface_edge_gets_exactly_one_welded_vertex() {
    let center = Vec3::new(11.0, 11.0, 11.0);
    let volume = sampled_volume(
        [23, 23, 23],
        [1.0; 3],
        Vec3::zero(),
        binary_ball(center, 7.3),
    );
    let mesh = marching_cubes(&volume, 0.0);
    // Count the lattice edges whose end voxels lie on different sides.
    let inside = |ix: usize, iy: usize, iz: usize| volume.data[volume.index(ix, iy, iz)] >= 0.0;
    let mut crossed = 0usize;
    for iz in 0..23 {
        for iy in 0..23 {
            for ix in 0..23 {
                crossed += usize::from(ix + 1 < 23 && inside(ix, iy, iz) != inside(ix + 1, iy, iz));
                crossed += usize::from(iy + 1 < 23 && inside(ix, iy, iz) != inside(ix, iy + 1, iz));
                crossed += usize::from(iz + 1 < 23 && inside(ix, iy, iz) != inside(ix, iy, iz + 1));
            }
        }
    }
    assert_eq!(mesh.vertices.len(), crossed);
    // No two vertices coincide.
    let mut keys: Vec<[u64; 3]> = mesh
        .vertices
        .iter()
        .map(|v| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()])
        .collect();
    keys.sort_unstable();
    keys.dedup();
    assert_eq!(keys.len(), mesh.vertices.len());
}

#[test]
fn anisotropic_voxels_and_a_shifted_origin_are_respected() {
    let spacing = [0.5, 0.8, 2.0];
    let origin = Vec3::new(-3.0, 1.0, 7.0);
    // 6 mm ball centred inside the volume, which spans [-3, 12] x [1, 33] x [7, 69].
    let center = Vec3::new(4.5, 17.0, 38.0);
    let radius = 6.0;
    let volume = sampled_volume([31, 41, 32], spacing, origin, |p| {
        radius - p.distance_to(&center)
    });
    let mesh = marching_cubes(&volume, 0.0);
    assert_watertight(&mesh);
    let exact = 4.0 / 3.0 * PI * radius.powi(3);
    let measured = signed_volume(&mesh);
    assert!(
        (measured / exact - 1.0).abs() < 0.03,
        "volume {measured} vs {exact}"
    );
    for v in &mesh.vertices {
        assert!((v.distance_to(&center) - radius).abs() < 0.3);
    }
}

#[test]
fn plane_between_voxel_layers_is_flat_with_exact_area_and_outward_normals() {
    // Voxels with z index >= 2 are inside; the surface is the plane halfway
    // between layers 1 and 2.
    let spacing = [0.7, 1.3, 2.5];
    let origin = Vec3::new(10.0, -5.0, 3.0);
    let (nx, ny, nz) = (9usize, 6usize, 5usize);
    let mut volume = ScalarVolume::new(nx, ny, nz, spacing, origin);
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                volume.set(ix, iy, iz, if iz >= 2 { 800.0 } else { -50.0 });
            }
        }
    }
    let mesh = marching_cubes(&volume, 375.0);
    let plane_z = origin.z + spacing[2] * (1.0 + (375.0 - -50.0) / (800.0 - -50.0));
    for v in &mesh.vertices {
        assert!((v.z - plane_z).abs() < 1e-12, "vertex off the plane: {v:?}");
    }
    let expected_area = (nx - 1) as f64 * spacing[0] * (ny - 1) as f64 * spacing[1];
    assert!((mesh.surface_area() - expected_area).abs() < 1e-9);
    // Outward = away from the structure, which lies at larger z: normals point to -z.
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]],
            mesh.vertices[t[1]],
            mesh.vertices[t[2]],
        );
        let normal = b.sub(&a).cross(&c.sub(&a));
        assert!(normal.z < 0.0, "triangle {t:?} faces into the structure");
    }
    // An open patch: its outline edges belong to a single triangle.
    let histogram = edge_valence_histogram(&mesh.triangles);
    assert!(histogram.iter().any(|&(valence, _)| valence == 1));
    assert!(histogram.iter().all(|&(valence, _)| valence <= 2));
}

#[test]
fn threshold_uses_greater_or_equal_for_inside_voxels() {
    // A voxel exactly at the iso value is inside, a NaN voxel outside.
    let mut volume = ScalarVolume::new(2, 2, 2, [1.0; 3], Vec3::zero());
    volume.set(0, 0, 0, 100.0);
    let mesh = marching_cubes(&volume, 100.0);
    assert_eq!(mesh.triangle_count(), 1);
    let mut nan_volume = ScalarVolume::new(2, 2, 2, [1.0; 3], Vec3::zero());
    nan_volume.set(0, 0, 0, f64::NAN);
    assert_eq!(marching_cubes(&nan_volume, 0.5).triangle_count(), 0);
    // NaN next to an inside voxel never produces a non-finite vertex.
    nan_volume.set(1, 0, 0, 10.0);
    let mesh = marching_cubes(&nan_volume, 0.5);
    assert!(mesh.triangle_count() > 0);
    assert!(
        mesh.vertices
            .iter()
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
    );
}

#[test]
fn single_inside_voxel_cuts_one_outward_triangle_at_the_edge_midpoints() {
    let mut volume = ScalarVolume::new(2, 2, 2, [2.0, 3.0, 4.0], Vec3::new(1.0, 1.0, 1.0));
    volume.set(1, 1, 1, 500.0);
    for v in volume.data.iter_mut().filter(|v| **v == 0.0) {
        *v = -500.0;
    }
    let mesh = marching_cubes(&volume, 0.0);
    assert_eq!(mesh.triangle_count(), 1);
    assert_eq!(mesh.vertex_count(), 3);
    // The inside corner is (3, 4, 5); the crossings are halfway along its edges.
    let expected = [
        Vec3::new(2.0, 4.0, 5.0),
        Vec3::new(3.0, 2.5, 5.0),
        Vec3::new(3.0, 4.0, 3.0),
    ];
    for e in expected {
        assert!(
            mesh.vertices.iter().any(|v| v.distance_to(&e) < 1e-12),
            "missing crossing {e:?}"
        );
    }
    let t = mesh.triangles[0];
    let (a, b, c) = (
        mesh.vertices[t[0]],
        mesh.vertices[t[1]],
        mesh.vertices[t[2]],
    );
    let normal = b.sub(&a).cross(&c.sub(&a));
    let centroid = a.add(&b).add(&c).scale(1.0 / 3.0);
    let to_inside_corner = Vec3::new(3.0, 4.0, 5.0).sub(&centroid);
    assert!(
        normal.dot(&to_inside_corner) < 0.0,
        "the normal must point away from the inside voxel"
    );
}

#[test]
fn every_corner_configuration_emits_the_table_triangles_facing_outwards() {
    for (case, &crossed_edges) in MC_EDGE_TABLE.iter().enumerate().take(255).skip(1) {
        let mut volume = ScalarVolume::new(2, 2, 2, [1.0; 3], Vec3::zero());
        for (corner, offset) in MC_CORNER_OFFSETS.iter().enumerate() {
            let value = if (case >> corner) & 1 == 1 {
                10.0
            } else {
                -10.0
            };
            volume.set(offset[0], offset[1], offset[2], value);
        }
        let mesh = marching_cubes(&volume, 0.0);
        let expected_triangles = case_triangles(case).count();
        assert_eq!(mesh.triangle_count(), expected_triangles, "case {case}");
        assert_eq!(
            mesh.vertex_count(),
            crossed_edges.count_ones() as usize,
            "case {case}: one vertex per crossed edge"
        );
        for v in &mesh.vertices {
            // All crossings are edge midpoints: exactly one coordinate is 0.5.
            let halves = [v.x, v.y, v.z]
                .iter()
                .filter(|c| (**c - 0.5).abs() < 1e-12)
                .count();
            assert_eq!(halves, 1, "case {case}: vertex {v:?}");
        }
        // The triangle normals oppose the trilinear gradient of the field.
        let value_of = |c: usize| if (case >> c) & 1 == 1 { 10.0 } else { -10.0 };
        let mut alignment = 0.0;
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let normal = b.sub(&a).cross(&c.sub(&a));
            let g = a.add(&b).add(&c).scale(1.0 / 3.0);
            let mut grad = [0.0; 3];
            for (corner, o) in MC_CORNER_OFFSETS.iter().enumerate() {
                let w = |axis: usize, p: f64| if o[axis] == 1 { p } else { 1.0 - p };
                let dw = |axis: usize| if o[axis] == 1 { 1.0 } else { -1.0 };
                grad[0] += value_of(corner) * dw(0) * w(1, g.y) * w(2, g.z);
                grad[1] += value_of(corner) * w(0, g.x) * dw(1) * w(2, g.z);
                grad[2] += value_of(corner) * w(0, g.x) * w(1, g.y) * dw(2);
            }
            alignment += normal.x * grad[0] + normal.y * grad[1] + normal.z * grad[2];
        }
        assert!(
            alignment < 0.0,
            "case {case}: normals face the structure ({alignment})"
        );
    }
}

#[test]
fn degenerate_volumes_give_empty_meshes() {
    let flat = ScalarVolume::new(8, 8, 1, [1.0; 3], Vec3::zero());
    assert_eq!(marching_cubes(&flat, 0.0).triangle_count(), 0);
    let empty = ScalarVolume::new(0, 0, 0, [1.0; 3], Vec3::zero());
    assert_eq!(marching_cubes(&empty, 0.0).vertex_count(), 0);
    // `data` shorter than `nx * ny * nz` must not panic.
    let mut truncated = ScalarVolume::new(4, 4, 4, [1.0; 3], Vec3::zero());
    truncated.data.truncate(10);
    let mesh = marching_cubes(&truncated, 0.0);
    assert_eq!((mesh.vertex_count(), mesh.triangle_count()), (0, 0));
}

#[test]
fn extraction_is_deterministic() {
    let center = Vec3::new(9.0, 9.0, 9.0);
    let volume = sampled_volume([19, 19, 19], [1.0; 3], Vec3::zero(), |p| {
        6.0 - p.distance_to(&center) + 0.4 * (p.x * 1.3).sin()
    });
    let first = marching_cubes(&volume, 0.0);
    let second = marching_cubes(&volume, 0.0);
    assert_eq!(first.triangles, second.triangles);
    assert_eq!(first.vertices, second.vertices);
}

#[test]
fn taubin_smoothing_keeps_the_extracted_mesh_watertight_and_nearly_volume_preserving() {
    let center = Vec3::new(15.5, 15.5, 15.5);
    let volume = sampled_volume(
        [32, 32, 32],
        [1.0; 3],
        Vec3::zero(),
        binary_ball(center, 9.0),
    );
    let mesh = marching_cubes(&volume, 0.0);
    let smoothed = taubin_smooth_medical_mesh(&mesh, 0.33, -0.34, 10);
    assert_watertight(&smoothed);
    let exact = 4.0 / 3.0 * PI * 9.0_f64.powi(3);
    assert!((signed_volume(&smoothed) / exact - 1.0).abs() < 0.04);
}
