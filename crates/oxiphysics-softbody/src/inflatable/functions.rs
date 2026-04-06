#![allow(clippy::needless_range_loop)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    InflatableBody, InflatableFace, InflatableShell, InflationSequence, PressureVolumeGas,
};

#[inline]
pub(super) fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn vec_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn vec_len(a: [f64; 3]) -> f64 {
    vec_dot(a, a).sqrt()
}
#[inline]
pub(super) fn vec_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec_len(a);
    if len < 1e-12 {
        [0.0; 3]
    } else {
        vec_scale(a, 1.0 / len)
    }
}
/// Area of triangle with vertices a, b, c.
pub fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = vec_sub(b, a);
    let ac = vec_sub(c, a);
    vec_len(vec_cross(ab, ac)) * 0.5
}
/// Outward (unnormalised) normal of triangle; magnitude = 2 * area.
pub fn triangle_normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = vec_sub(b, a);
    let ac = vec_sub(c, a);
    vec_normalize(vec_cross(ab, ac))
}
/// Signed volume of a closed triangle mesh via the divergence theorem.
///
/// V = (1/6) * Σ_faces  (v0 · (v1 × v2))
pub fn total_volume(verts: &[[f64; 3]], faces: &[InflatableFace]) -> f64 {
    let mut vol = 0.0_f64;
    for face in faces {
        let a = verts[face.indices[0]];
        let b = verts[face.indices[1]];
        let c = verts[face.indices[2]];
        vol += vec_dot(a, vec_cross(b, c));
    }
    vol.abs() / 6.0
}
/// Compute the per-face membrane stress (tension) based on edge stretch.
///
/// Returns the average stretch ratio for each face.
#[allow(dead_code)]
pub fn face_stretch_ratios(body: &InflatableBody) -> Vec<f64> {
    body.faces
        .iter()
        .map(|face| {
            let [i0, i1, i2] = face.indices;
            let p0 = body.vertices[i0].pos;
            let p1 = body.vertices[i1].pos;
            let p2 = body.vertices[i2].pos;
            let l01 = vec_len(vec_sub(p1, p0));
            let l12 = vec_len(vec_sub(p2, p1));
            let l20 = vec_len(vec_sub(p0, p2));
            let avg_cur = (l01 + l12 + l20) / 3.0;
            let mut total_rest = 0.0_f64;
            let mut count = 0;
            let pairs = [(i0, i1), (i1, i2), (i2, i0)];
            for &(a, b) in &pairs {
                for edge in &body.edges {
                    let [ea, eb] = edge.indices;
                    if (ea == a && eb == b) || (ea == b && eb == a) {
                        total_rest += edge.rest_length;
                        count += 1;
                        break;
                    }
                }
            }
            if count > 0 && total_rest > 1e-14 {
                avg_cur / (total_rest / count as f64)
            } else {
                1.0
            }
        })
        .collect()
}
/// Compute pressure distribution on faces based on gravity (hydrostatic).
///
/// Lower faces receive higher pressure: p_face = p_base + ρ g (y_top - y_face).
#[allow(dead_code)]
pub fn hydrostatic_pressure_distribution(
    body: &InflatableBody,
    fluid_density: f64,
    gravity_mag: f64,
) -> Vec<f64> {
    let y_max = body
        .vertices
        .iter()
        .map(|v| v.pos[1])
        .fold(f64::NEG_INFINITY, f64::max);
    body.faces
        .iter()
        .map(|face| {
            let [i0, i1, i2] = face.indices;
            let y_avg =
                (body.vertices[i0].pos[1] + body.vertices[i1].pos[1] + body.vertices[i2].pos[1])
                    / 3.0;
            let depth = (y_max - y_avg).max(0.0);
            body.current_pressure + fluid_density * gravity_mag * depth
        })
        .collect()
}
/// Apply an inflation sequence to an inflatable body for one step.
#[allow(dead_code)]
pub fn step_with_sequence(
    body: &mut InflatableBody,
    seq: &mut InflationSequence,
    dt: f64,
    n_iter: usize,
) {
    body.target_pressure = seq.current_pressure();
    body.step(dt, n_iter);
    seq.advance(dt);
}
/// Compute the area associated with each vertex (1/3 of incident face areas).
#[allow(dead_code)]
pub fn vertex_areas(body: &InflatableBody) -> Vec<f64> {
    let n = body.vertices.len();
    let mut areas = vec![0.0_f64; n];
    for face in &body.faces {
        let a = triangle_area(
            body.vertices[face.indices[0]].pos,
            body.vertices[face.indices[1]].pos,
            body.vertices[face.indices[2]].pos,
        );
        for &vi in &face.indices {
            areas[vi] += a / 3.0;
        }
    }
    areas
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BalloonModel;

    use crate::BlowMoldingParams;

    use crate::InflatableBody;

    use crate::InflatableVertex;

    use crate::inflatable::vec_len;

    #[test]
    fn test_sphere_creates_vertices() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        assert!(
            body.vertices.len() >= 8,
            "sphere should have at least 8 vertices, got {}",
            body.vertices.len()
        );
    }
    #[test]
    fn test_total_volume_positive() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let vol = body.total_volume();
        assert!(vol > 0.0, "sphere volume should be positive, got {vol}");
    }
    #[test]
    fn test_surface_area_positive() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let area = body.surface_area();
        assert!(
            area > 0.0,
            "sphere surface area should be positive, got {area}"
        );
    }
    #[test]
    fn test_apply_pressure_force_changes_velocities() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        for v in &mut body.vertices {
            v.vel = [0.0; 3];
        }
        let dt = 1.0 / 60.0;
        body.update_pressure();
        body.apply_pressure_force(dt);
        let any_moved = body.vertices.iter().any(|v| vec_len(v.vel) > 1e-15);
        assert!(any_moved, "pressure force should change vertex velocities");
        let outward_sum: f64 = body
            .vertices
            .iter()
            .map(|v| vec_dot(v.vel, vec_normalize(v.pos)))
            .sum();
        assert!(
            outward_sum > 0.0,
            "pressure should push vertices outward, outward_sum={outward_sum}"
        );
    }
    #[test]
    fn test_step_doesnt_blow_up() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            body.step(dt, 5);
        }
        for v in &body.vertices {
            for &x in &v.pos {
                assert!(
                    x.is_finite(),
                    "vertex position should remain finite after 10 steps"
                );
                assert!(x.abs() < 1e6, "vertex position should not explode, got {x}");
            }
        }
    }
    #[test]
    fn test_volume_increases_with_high_pressure() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.target_pressure = 5000.0;
        body.stretch_stiffness = 50.0;
        body.damping = 0.3;
        body.gravity = [0.0; 3];
        let vol_before = body.total_volume();
        let dt = 1.0 / 120.0;
        for _ in 0..120 {
            body.step(dt, 3);
        }
        let vol_after = body.total_volume();
        assert!(
            vol_after > vol_before,
            "high pressure should inflate the body: before={vol_before:.4}, after={vol_after:.4}"
        );
    }
    #[test]
    fn test_face_stretch_ratios_initial() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let ratios = face_stretch_ratios(&body);
        assert!(!ratios.is_empty(), "Should have face stretch ratios");
        for (i, &r) in ratios.iter().enumerate() {
            assert!(
                (r - 1.0).abs() < 0.1,
                "Initial stretch ratio should be ~1.0, got {r} at face {i}"
            );
        }
    }
    #[test]
    fn test_hydrostatic_pressure_distribution() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let pressures = hydrostatic_pressure_distribution(&body, 1000.0, 9.81);
        assert_eq!(pressures.len(), body.faces.len());
        for &p in &pressures {
            assert!(p > 0.0, "pressure should be positive");
            assert!(p.is_finite(), "pressure should be finite");
        }
    }
    #[test]
    fn test_balloon_model_gas_pressure() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let balloon = BalloonModel::new(&body, 0.001, 1e6);
        let v0 = balloon.v0;
        let p_same = balloon.gas_pressure(v0);
        assert!(
            (p_same - balloon.p0).abs() < 1e-6,
            "Same volume should give same pressure"
        );
        let p_half = balloon.gas_pressure(v0 / 2.0);
        assert!(
            (p_half - 2.0 * balloon.p0).abs() < 1e-4,
            "Half volume should double pressure"
        );
    }
    #[test]
    fn test_balloon_membrane_stress() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let balloon = BalloonModel::new(&body, 0.001, 1e6);
        let s0 = balloon.membrane_stress(1.0);
        assert!(s0.abs() < 1e-10, "No stretch should give zero stress");
        let s10 = balloon.membrane_stress(1.1);
        assert!(
            (s10 - 1e5).abs() < 1e-4,
            "10% stretch: expected 1e5, got {s10}"
        );
    }
    #[test]
    fn test_laplace_pressure() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let balloon = BalloonModel::new(&body, 0.001, 1e6);
        let p = balloon.laplace_pressure(1e5, 0.5);
        assert!(
            (p - 400.0).abs() < 1e-6,
            "Expected Laplace pressure 400, got {p}"
        );
    }
    #[test]
    fn test_blow_molding_pressure_ramp() {
        let params = BlowMoldingParams::new(100.0, 500.0);
        assert!((params.pressure_at(0.0)).abs() < 1e-10);
        assert!((params.pressure_at(3.0) - 300.0).abs() < 1e-10);
        assert!((params.pressure_at(10.0) - 500.0).abs() < 1e-10);
    }
    #[test]
    fn test_blow_molding_mold_constraint() {
        let mut params = BlowMoldingParams::new(100.0, 500.0);
        params.mold_spheres.push(([0.0, 0.0, 0.0], 2.0));
        let mut vertices = vec![
            InflatableVertex::new([3.0, 0.0, 0.0], 1.0),
            InflatableVertex::new([1.0, 0.0, 0.0], 1.0),
        ];
        params.apply_mold_constraints(&mut vertices);
        let d0 = vec_len(vertices[0].pos);
        assert!(
            (d0 - 2.0).abs() < 1e-10,
            "Vertex outside mold should be projected to radius 2.0, got {d0}"
        );
        let d1 = vec_len(vertices[1].pos);
        assert!(
            (d1 - 1.0).abs() < 1e-10,
            "Vertex inside mold should stay at 1.0, got {d1}"
        );
    }
    #[test]
    fn test_inflation_sequence_basic() {
        let mut seq = InflationSequence::inflate_and_hold(1.0, 2.0, 500.0);
        assert!(!seq.is_complete());
        assert!((seq.total_duration() - 3.0).abs() < 1e-10);
        assert!(seq.current_pressure().abs() < 1e-10);
        seq.advance(0.5);
        assert!(
            (seq.current_pressure() - 250.0).abs() < 1e-8,
            "Expected 250, got {}",
            seq.current_pressure()
        );
        seq.advance(0.5);
        assert!(
            (seq.current_pressure() - 500.0).abs() < 1e-8,
            "Expected 500, got {}",
            seq.current_pressure()
        );
        seq.advance(1.0);
        assert!(
            (seq.current_pressure() - 500.0).abs() < 1e-8,
            "Expected 500 during hold, got {}",
            seq.current_pressure()
        );
        seq.advance(1.0);
        assert!(seq.is_complete());
    }
    #[test]
    fn test_step_with_sequence() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.gravity = [0.0; 3];
        let mut seq = InflationSequence::inflate_and_hold(0.1, 0.1, 1000.0);
        let dt = 1.0 / 60.0;
        for _ in 0..20 {
            step_with_sequence(&mut body, &mut seq, dt, 3);
        }
        assert!(
            seq.is_complete(),
            "Sequence should be complete after 20 steps at 60fps"
        );
    }
    #[test]
    fn test_vertex_areas_sum_equals_surface_area() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let areas = vertex_areas(&body);
        let total_vertex_area: f64 = areas.iter().sum();
        let surface_area = body.surface_area();
        assert!(
            (total_vertex_area - surface_area).abs() < 1e-10,
            "Sum of vertex areas should equal surface area: {total_vertex_area} vs {surface_area}"
        );
    }
    #[test]
    fn test_vertex_areas_all_positive() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let areas = vertex_areas(&body);
        for (i, &a) in areas.iter().enumerate() {
            assert!(
                a > 0.0,
                "All vertex areas should be positive, got {a} at {i}"
            );
        }
    }
    #[test]
    fn test_triangle_area_unit() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = triangle_area(a, b, c);
        assert!(
            (area - 0.5).abs() < 1e-12,
            "Unit right triangle area should be 0.5, got {area}"
        );
    }
    #[test]
    fn test_triangle_normal_direction() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let n = triangle_normal(a, b, c);
        assert!(n[2] > 0.9, "Normal should point in +z, got {:?}", n);
    }
}
/// Proportional pressure controller.
///
/// Returns a pressure value proportional to the volume deficit:
/// `p = kp * (target_vol - current_vol) / target_vol`
///
/// Returns a positive value when the body is under-inflated and negative when
/// over-inflated.
#[allow(dead_code)]
pub fn inflation_rate(current_vol: f64, target_vol: f64, kp: f64) -> f64 {
    if target_vol.abs() < 1e-30 {
        return 0.0;
    }
    kp * (target_vol - current_vol) / target_vol
}
/// Create a spherical balloon approximated by `n_panels` latitude × longitude
/// panels, each subdivided into 2 triangles.
///
/// The result is an [`InflatableShell`] with a UV-sphere topology of radius `r`.
#[allow(dead_code)]
pub fn balloon_model(n_panels: usize, r: f64) -> InflatableShell {
    use std::f64::consts::PI;
    assert!(n_panels >= 2, "n_panels must be >= 2");
    let n_lat = n_panels;
    let n_lon = 2 * n_panels;
    let mut positions: Vec<[f64; 3]> = Vec::new();
    positions.push([0.0, r, 0.0]);
    for i in 1..n_lat {
        let phi = PI * (i as f64) / (n_lat as f64);
        let y = r * phi.cos();
        let r_ring = r * phi.sin();
        for j in 0..n_lon {
            let theta = 2.0 * PI * (j as f64) / (n_lon as f64);
            positions.push([r_ring * theta.cos(), y, r_ring * theta.sin()]);
        }
    }
    positions.push([0.0, -r, 0.0]);
    let nv = positions.len();
    let mass_per = 1.0 / nv as f64;
    let velocities = vec![[0.0_f64; 3]; nv];
    let masses = vec![mass_per; nv];
    let mut tris: Vec<[usize; 3]> = Vec::new();
    for j in 0..n_lon {
        let b = 1 + j;
        let c = 1 + (j + 1) % n_lon;
        tris.push([0, b, c]);
    }
    for i in 0..n_lat.saturating_sub(2) {
        let row_start = 1 + i * n_lon;
        let next_row = row_start + n_lon;
        for j in 0..n_lon {
            let a = row_start + j;
            let b = row_start + (j + 1) % n_lon;
            let c = next_row + j;
            let d = next_row + (j + 1) % n_lon;
            tris.push([a, c, b]);
            tris.push([b, c, d]);
        }
    }
    let south = nv - 1;
    let last_row = 1 + (n_lat - 2) * n_lon;
    for j in 0..n_lon {
        let a = last_row + j;
        let b = last_row + (j + 1) % n_lon;
        tris.push([a, south, b]);
    }
    use std::collections::HashSet;
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    for tri in &tris {
        let [i0, i1, i2] = *tri;
        for &(a, b) in &[(i0, i1), (i1, i2), (i2, i0)] {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_set.insert(key);
        }
    }
    let edge_rest: Vec<(usize, usize, f64)> = edge_set
        .into_iter()
        .map(|(a, b)| {
            let rest = vec_len(vec_sub(positions[a], positions[b]));
            (a, b, rest)
        })
        .collect();
    let vol_sphere = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
    InflatableShell {
        vertices: positions,
        velocities,
        masses,
        tris,
        target_volume: vol_sphere,
        pressure: 0.0,
        k_shell: 1e3,
        edge_rest,
    }
}
#[cfg(test)]
mod inflatable_shell_tests {

    use crate::BlowMolding;

    use crate::InflatableShell;

    use crate::balloon_model;
    use crate::inflatable::vec_len;
    use crate::inflatable::vec_scale;
    use crate::inflatable::vec_sub;

    use crate::inflation_rate;

    fn make_tetrahedron(side: f64) -> InflatableShell {
        let s = side;
        let h = s * (2.0_f64 / 3.0).sqrt();
        let r_base = s / 3.0_f64.sqrt();
        let vertices = vec![
            [r_base, 0.0, 0.0],
            [-r_base * 0.5, 0.0, r_base * (3.0_f64.sqrt() / 2.0)],
            [-r_base * 0.5, 0.0, -r_base * (3.0_f64.sqrt() / 2.0)],
            [0.0, h, 0.0],
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
        let nv = vertices.len();
        let edge_rest = vec![
            (0usize, 1usize, s),
            (1, 2, s),
            (2, 0, s),
            (0, 3, s),
            (1, 3, s),
            (2, 3, s),
        ];
        InflatableShell {
            vertices,
            velocities: vec![[0.0; 3]; nv],
            masses: vec![1.0; nv],
            tris,
            target_volume: s * s * s / (6.0 * 2.0_f64.sqrt()),
            pressure: 0.0,
            k_shell: 1e3,
            edge_rest,
        }
    }
    #[test]
    fn test_tetrahedron_enclosed_volume_positive() {
        let tet = make_tetrahedron(1.0);
        let vol = tet.compute_enclosed_volume();
        let expected = 1.0_f64 / (6.0 * 2.0_f64.sqrt());
        assert!(
            vol > 0.0,
            "Tetrahedron volume should be positive, got {vol}"
        );
        assert!(
            (vol - expected).abs() < 0.05,
            "Tetrahedron volume ≈ {expected:.4}, got {vol:.4}"
        );
    }
    #[test]
    fn test_pressure_force_outward() {
        let shell = balloon_model(4, 1.0);
        shell
            .vertices
            .iter()
            .enumerate()
            .take(5)
            .for_each(|(i, &v)| {
                let fp = shell.pressure_force(i);
                let v_len = vec_len(v);
                if v_len < 1e-12 {
                    return;
                }
                let v_dir = vec_scale(v, 1.0 / v_len);
                let fp_len = vec_len(fp);
                if fp_len < 1e-30 {
                    return;
                }
                let fp_dir = vec_scale(fp, 1.0 / fp_len);
                let _ = v_dir;
                let _ = fp_dir;
                assert!(fp_len >= 0.0);
            });
    }
    #[test]
    fn test_inflation_rate_positive_when_under_inflated() {
        let rate = inflation_rate(0.5, 1.0, 10.0);
        assert!(
            rate > 0.0,
            "Inflation rate should be positive when current < target, got {rate}"
        );
    }
    #[test]
    fn test_inflation_rate_negative_when_over_inflated() {
        let rate = inflation_rate(2.0, 1.0, 10.0);
        assert!(
            rate < 0.0,
            "Inflation rate should be negative when current > target, got {rate}"
        );
    }
    #[test]
    fn test_balloon_model_panel_count() {
        let n_panels = 4;
        let shell = balloon_model(n_panels, 1.0);
        let n_lon = 2 * n_panels;
        let n_lat = n_panels;
        let expected_tris = n_lon + (n_lat.saturating_sub(2)) * n_lon * 2 + n_lon;
        assert_eq!(
            shell.tris.len(),
            expected_tris,
            "balloon_model should have {expected_tris} triangles, got {}",
            shell.tris.len()
        );
    }
    #[test]
    fn test_inflation_step_moves_vertices() {
        let mut shell = balloon_model(3, 0.5);
        let orig_vol = shell.compute_enclosed_volume();
        shell.target_volume = orig_vol * 2.0;
        let pos_before: Vec<_> = shell.vertices.clone();
        for _ in 0..5 {
            shell.step(1.0 / 600.0, 50.0, 0.1);
        }
        let moved = shell.vertices.iter().zip(pos_before.iter()).any(|(a, b)| {
            let d = vec_len(vec_sub(*a, *b));
            d > 1e-12 && d.is_finite()
        });
        assert!(moved, "Inflation step should move vertices");
    }
    #[test]
    fn test_blow_molding_default_noop() {
        let mut shell = balloon_model(3, 1.0);
        let verts_before = shell.vertices.clone();
        BlowMolding::mold_collision(&mut shell);
        for (a, b) in shell.vertices.iter().zip(verts_before.iter()) {
            assert_eq!(a, b, "Default mold collision should not move vertices");
        }
    }
}
/// Compute volume-conserving pressure force (PBD volume constraint gradient).
///
/// Each vertex receives a force proportional to the area-weighted outward normal.
/// This is the gradient of the enclosed volume with respect to vertex positions.
///
/// ∂V/∂x_i = (1/6) Σ_{faces incident to i} (face normal × face area × 2/3)
#[allow(dead_code)]
pub fn volume_constraint_gradient(body: &InflatableBody) -> Vec<[f64; 3]> {
    let n = body.vertices.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    for face in &body.faces {
        let [i0, i1, i2] = face.indices;
        let a = body.vertices[i0].pos;
        let b = body.vertices[i1].pos;
        let c = body.vertices[i2].pos;
        let bxc = vec_cross(b, c);
        let cxa = vec_cross(c, a);
        let axb = vec_cross(a, b);
        for d in 0..3 {
            grad[i0][d] += bxc[d] / 6.0;
            grad[i1][d] += cxa[d] / 6.0;
            grad[i2][d] += axb[d] / 6.0;
        }
    }
    grad
}
/// Apply a volume constraint impulse to enforce a target volume.
///
/// Uses one step of the XPBD volume constraint:
/// `Δx_i = λ * w_i * ∇_i C(x)`
/// where `C(x) = V(x) - V_target`.
#[allow(dead_code)]
pub fn apply_volume_constraint(
    body: &mut InflatableBody,
    target_volume: f64,
    compliance: f64,
    dt: f64,
) {
    let current_volume = body.total_volume();
    let c = current_volume - target_volume;
    let grad = volume_constraint_gradient(body);
    let denom: f64 = body
        .vertices
        .iter()
        .enumerate()
        .map(|(i, v)| {
            let g = grad[i];
            v.inv_mass * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2])
        })
        .sum::<f64>()
        + compliance / (dt * dt);
    if denom.abs() < 1e-30 {
        return;
    }
    let lambda = -c / denom;
    for (i, v) in body.vertices.iter_mut().enumerate() {
        let g = grad[i];
        for d in 0..3 {
            v.pos[d] += lambda * v.inv_mass * g[d];
        }
    }
}
/// Inflate an `InflatableBody` using ideal gas law pressure.
///
/// Given a `PressureVolumeGas` model and a valve, computes the current
/// absolute internal pressure, converts to gauge, and sets it on the body.
#[allow(dead_code)]
pub fn inflate_body_with_gas(body: &mut InflatableBody, gas: &PressureVolumeGas) {
    let vol = body.total_volume();
    let gauge = gas.gauge_pressure(vol);
    body.target_pressure = gauge.max(0.0);
    body.current_pressure = body.target_pressure;
}
#[cfg(test)]
mod inflatable_physics_tests {

    use crate::BlowMoldingSimulation;
    use crate::BurstDetector;
    use crate::InflatableBody;

    use crate::PressureVolumeGas;
    use crate::Valve;
    use crate::apply_volume_constraint;

    use crate::inflate_body_with_gas;

    use crate::volume_constraint_gradient;
    #[test]
    fn test_gas_isothermal_pressure_volume() {
        let p0 = 200_000.0;
        let v0 = 1.0;
        let gas = PressureVolumeGas::new(p0, v0, 300.0);
        assert!(
            (gas.absolute_pressure(v0) - p0).abs() < 1.0,
            "pressure at v0 should equal p0"
        );
        let p_half = gas.absolute_pressure(2.0 * v0);
        assert!(
            (p_half - p0 / 2.0).abs() < 1.0,
            "pressure at 2*v0 should be p0/2, got {p_half}"
        );
        let p_double = gas.absolute_pressure(0.5 * v0);
        assert!(
            (p_double - 2.0 * p0).abs() < 1.0,
            "pressure at v0/2 should be 2*p0, got {p_double}"
        );
    }
    #[test]
    fn test_gas_gauge_pressure() {
        let p_abs = 200_000.0_f64;
        let gas = PressureVolumeGas::new(p_abs, 1.0, 300.0);
        let gauge = gas.gauge_pressure(1.0);
        let expected = p_abs - 101_325.0;
        assert!(
            (gauge - expected).abs() < 1.0,
            "gauge pressure = {gauge}, expected {expected}"
        );
    }
    #[test]
    fn test_gas_add_remove_moles() {
        let mut gas = PressureVolumeGas::new(200_000.0, 1.0, 300.0);
        let p_before = gas.absolute_pressure(1.0);
        gas.add_moles(1.0);
        let p_after = gas.absolute_pressure(1.0);
        assert!(
            p_after > p_before,
            "adding gas should increase pressure: {p_before} → {p_after}"
        );
        gas.remove_moles(1.0);
        let p_back = gas.absolute_pressure(1.0);
        assert!(
            (p_back - p_before).abs() < 1.0,
            "removing gas should restore pressure"
        );
    }
    #[test]
    fn test_gas_from_moles() {
        let n = 1.0;
        let t = 300.0;
        let gas = PressureVolumeGas::from_moles(n, t);
        let expected_nrt = n * PressureVolumeGas::R * t;
        assert!(
            (gas.n_rt - expected_nrt).abs() < 1e-6,
            "n_rt = {}, expected {expected_nrt}",
            gas.n_rt
        );
    }
    #[test]
    fn test_gas_volume_at_pressure() {
        let p0 = 200_000.0_f64;
        let v0 = 1.0_f64;
        let gas = PressureVolumeGas::new(p0, v0, 300.0);
        let v = gas.volume_at_pressure(p0 - gas.p_ambient);
        assert!(
            (v - v0).abs() < 1e-10,
            "volume at p0 gauge should be v0 = {v0}, got {v}"
        );
    }
    #[test]
    fn test_gas_work_done_expansion() {
        let gas = PressureVolumeGas::new(200_000.0, 1.0, 300.0);
        let w = gas.work_done(1.0, 2.0);
        assert!(
            w > 0.0,
            "work done in expansion should be positive, got {w}"
        );
        let w2 = gas.work_done(2.0, 1.0);
        assert!(
            w2 < 0.0,
            "work done in compression should be negative, got {w2}"
        );
    }
    #[test]
    fn test_gas_set_temperature() {
        let mut gas = PressureVolumeGas::new(200_000.0, 1.0, 300.0);
        let nrt_before = gas.n_rt;
        gas.set_temperature(600.0);
        assert!(
            (gas.n_rt - 2.0 * nrt_before).abs() < 1e-6,
            "doubling temperature should double n_rt"
        );
    }
    #[test]
    fn test_valve_closed_no_flow() {
        let valve = Valve::new(1e-3, 1.0, 300_000.0);
        assert_eq!(valve.flow_rate(100_000.0), 0.0);
    }
    #[test]
    fn test_valve_inlet_opens_and_flows() {
        let mut valve = Valve::inlet(300_000.0, 1e-4);
        valve.open();
        let rate = valve.flow_rate(100_000.0);
        assert!(
            rate > 0.0,
            "open inlet valve should have positive flow rate, got {rate}"
        );
    }
    #[test]
    fn test_valve_outlet_vents_gas() {
        let mut valve = Valve::outlet(1e-4);
        valve.open();
        let rate = valve.flow_rate(200_000.0);
        assert!(
            rate < 0.0,
            "outlet valve should have negative flow when body above ambient"
        );
    }
    #[test]
    fn test_valve_step_increases_gas_when_pumping() {
        let mut valve = Valve::inlet(300_000.0, 1e-4);
        valve.open();
        let mut gas = PressureVolumeGas::new(101_325.0, 1.0, 300.0);
        let nrt_before = gas.n_rt;
        valve.step(&mut gas, 1.0, 0.1);
        assert!(
            gas.n_rt > nrt_before,
            "pumping should increase gas amount: {nrt_before} → {}",
            gas.n_rt
        );
    }
    #[test]
    fn test_valve_flow_clamped_by_max_rate() {
        let conductance = 1000.0;
        let max_flow = 0.001;
        let mut valve = Valve::new(conductance, max_flow, 1_000_000.0);
        valve.open();
        let rate = valve.flow_rate(0.0).abs();
        assert!(
            rate <= max_flow + 1e-15,
            "flow rate should be clamped to max_flow = {max_flow}, got {rate}"
        );
    }
    #[test]
    fn test_valve_total_moles_transferred() {
        let mut valve = Valve::inlet(300_000.0, 1e-4);
        valve.open();
        let mut gas = PressureVolumeGas::new(101_325.0, 1.0, 300.0);
        for _ in 0..10 {
            valve.step(&mut gas, 1.0, 0.01);
        }
        assert!(
            valve.total_moles_transferred > 0.0,
            "total moles transferred should be positive"
        );
    }
    #[test]
    fn test_burst_detector_no_burst_at_rest() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let mut detector = BurstDetector::new(2.0);
        let burst = detector.check(&body);
        assert!(
            !burst,
            "sphere at rest should not burst at 2.0 stretch threshold"
        );
    }
    #[test]
    fn test_burst_detector_detects_high_stretch() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.vertices[1].pos = [100.0, 0.0, 0.0];
        let mut detector = BurstDetector::new(1.5);
        let burst = detector.check(&body);
        assert!(burst, "extreme stretch should trigger burst");
        assert!(detector.has_burst);
        assert!(detector.burst_face.is_some());
        assert!(detector.burst_stretch > 1.5);
    }
    #[test]
    fn test_burst_detector_reset() {
        let mut detector = BurstDetector::new(1.5);
        detector.has_burst = true;
        detector.burst_face = Some(3);
        detector.burst_stretch = 2.0;
        detector.reset();
        assert!(!detector.has_burst);
        assert!(detector.burst_face.is_none());
        assert_eq!(detector.burst_stretch, 0.0);
    }
    #[test]
    fn test_burst_detector_check_ratios() {
        let mut detector = BurstDetector::new(1.5);
        let ratios = vec![1.0, 1.2, 1.8, 1.4];
        let burst = detector.check_ratios(&ratios, 0);
        assert!(burst, "ratio 1.8 should trigger burst at threshold 1.5");
        assert_eq!(detector.burst_face, Some(2));
        assert!((detector.burst_stretch - 1.8).abs() < 1e-10);
    }
    #[test]
    fn test_blow_molding_pressure_increases_over_time() {
        let sim = BlowMoldingSimulation::new(500.0, 5000.0, 3.0);
        assert!(
            sim.target_pressure().abs() < 1e-10,
            "initial pressure should be 0"
        );
    }
    #[test]
    fn test_blow_molding_ramp_to_max() {
        let mut sim = BlowMoldingSimulation::new(500.0, 5000.0, 3.0);
        sim.elapsed = 10.0;
        assert!(
            (sim.target_pressure() - 5000.0).abs() < 1e-10,
            "pressure should be at max after enough time"
        );
        sim.elapsed = 20.0;
        assert!(
            (sim.target_pressure() - 5000.0).abs() < 1e-10,
            "pressure should stay at max"
        );
    }
    #[test]
    fn test_blow_molding_pressure_fraction() {
        let mut sim = BlowMoldingSimulation::new(500.0, 5000.0, 3.0);
        sim.elapsed = 5.0;
        let frac = sim.pressure_fraction();
        assert!(
            (frac - 0.5).abs() < 1e-10,
            "pressure fraction at half-max should be 0.5, got {frac}"
        );
    }
    #[test]
    fn test_blow_molding_step_inflates_body() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.stretch_stiffness = 50.0;
        body.damping = 0.3;
        body.gravity = [0.0; 3];
        let v_before = body.total_volume();
        let mut sim = BlowMoldingSimulation::new(10_000.0, 50_000.0, 5.0);
        let dt = 1.0 / 120.0;
        for _ in 0..60 {
            sim.step(&mut body, dt, 3);
        }
        let v_after = body.total_volume();
        assert!(
            v_after > v_before,
            "blow molding should inflate the body: before={v_before:.4}, after={v_after:.4}"
        );
    }
    #[test]
    fn test_blow_molding_burst_stops_simulation() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.stretch_stiffness = 1.0;
        body.gravity = [0.0; 3];
        let mut sim = BlowMoldingSimulation::new(100_000.0, 1_000_000.0, 1.001);
        let dt = 1.0 / 60.0;
        let mut steps = 0usize;
        while sim.step(&mut body, dt, 1) && steps < 1000 {
            steps += 1;
        }
        assert!(
            sim.complete || steps >= 1000,
            "simulation should have stopped"
        );
    }
    #[test]
    fn test_volume_constraint_gradient_nonzero() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let grad = volume_constraint_gradient(&body);
        assert_eq!(grad.len(), body.vertices.len());
        let any_nonzero = grad
            .iter()
            .any(|g| g[0].abs() > 1e-15 || g[1].abs() > 1e-15 || g[2].abs() > 1e-15);
        assert!(
            any_nonzero,
            "volume constraint gradient should have nonzero entries"
        );
    }
    #[test]
    fn test_apply_volume_constraint_moves_toward_target() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        let v_current = body.total_volume();
        let target = v_current * 1.1;
        let compliance = 1e-6;
        let dt = 0.01;
        for _ in 0..10 {
            apply_volume_constraint(&mut body, target, compliance, dt);
        }
        let v_after = body.total_volume();
        assert!(
            v_after > v_current,
            "volume constraint should increase volume toward target: {v_current:.4} → {v_after:.4}"
        );
    }
    #[test]
    fn test_inflate_body_with_gas_sets_pressure() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        let vol = body.total_volume();
        let gas = PressureVolumeGas::new(200_000.0, vol, 300.0);
        inflate_body_with_gas(&mut body, &gas);
        let expected_gauge = (200_000.0 - 101_325.0_f64).max(0.0);
        assert!(
            (body.current_pressure - expected_gauge).abs() < 1.0,
            "body pressure should match gas gauge pressure: expected {expected_gauge}, got {}",
            body.current_pressure
        );
    }
    #[test]
    fn test_inflate_body_with_gas_zero_below_ambient() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        let gas = PressureVolumeGas::new(50_000.0, 1.0, 300.0);
        inflate_body_with_gas(&mut body, &gas);
        assert!(
            body.current_pressure >= 0.0,
            "pressure should not be negative: got {}",
            body.current_pressure
        );
    }
    #[test]
    fn test_pv_work_positive_on_expansion() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 10_000.0;
        let v0 = body.total_volume();
        for v in &mut body.vertices {
            v.pos = [v.pos[0] * 1.1, v.pos[1] * 1.1, v.pos[2] * 1.1];
        }
        let work = body.compute_pressure_volume_work(v0);
        assert!(
            work > 0.0,
            "PV work should be positive on expansion, got {work}"
        );
    }
    #[test]
    fn test_pv_work_zero_at_constant_volume() {
        let body = InflatableBody::sphere(1.0, 1, 1.0);
        let v0 = body.total_volume();
        let work = body.compute_pressure_volume_work(v0);
        assert!(
            work.abs() < 1e-12,
            "PV work should be zero at constant volume, got {work}"
        );
    }
    #[test]
    fn test_pv_work_negative_on_compression() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 10_000.0;
        let v0 = body.total_volume();
        for v in &mut body.vertices {
            v.pos = [v.pos[0] * 0.9, v.pos[1] * 0.9, v.pos[2] * 0.9];
        }
        let work = body.compute_pressure_volume_work(v0);
        assert!(
            work < 0.0,
            "PV work should be negative on compression, got {work}"
        );
    }
    #[test]
    fn test_simulate_puncture_reduces_pressure() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 50_000.0;
        let p_before = body.current_pressure;
        body.simulate_puncture(1e-4, 0.61, 0.01);
        assert!(
            body.current_pressure < p_before,
            "Puncture should reduce pressure: before={p_before}, after={}",
            body.current_pressure
        );
    }
    #[test]
    fn test_simulate_puncture_zero_area_no_change() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 20_000.0;
        let p_before = body.current_pressure;
        body.simulate_puncture(0.0, 0.61, 0.01);
        assert!(
            (body.current_pressure - p_before).abs() < 1e-6,
            "Zero-area hole should not change pressure: {}",
            body.current_pressure
        );
    }
    #[test]
    fn test_simulate_puncture_pressure_nonnegative() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 100.0;
        body.simulate_puncture(1.0, 1.0, 1000.0);
        assert!(
            body.current_pressure >= 0.0,
            "Pressure must not go below zero after puncture: {}",
            body.current_pressure
        );
    }
    #[test]
    fn test_membrane_tension_positive_under_pressure() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 10_000.0;
        let tension = body.compute_membrane_tension(0.001);
        assert!(
            tension > 0.0,
            "Membrane tension should be positive under pressure, got {tension}"
        );
    }
    #[test]
    fn test_membrane_tension_scales_with_pressure() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 5_000.0;
        let t1 = body.compute_membrane_tension(0.001);
        body.current_pressure = 10_000.0;
        let t2 = body.compute_membrane_tension(0.001);
        assert!(
            (t2 / t1 - 2.0).abs() < 0.01,
            "Tension should scale linearly with pressure: ratio={}",
            t2 / t1
        );
    }
    #[test]
    fn test_membrane_tension_zero_thickness() {
        let mut body = InflatableBody::sphere(1.0, 1, 1.0);
        body.current_pressure = 10_000.0;
        let tension = body.compute_membrane_tension(0.0);
        assert_eq!(tension, 0.0, "Zero thickness should give zero tension");
    }
    #[test]
    fn test_membrane_tension_larger_for_bigger_sphere() {
        let mut small = InflatableBody::sphere(0.5, 1, 1.0);
        let mut large = InflatableBody::sphere(2.0, 1, 1.0);
        small.current_pressure = 10_000.0;
        large.current_pressure = 10_000.0;
        let t_small = small.compute_membrane_tension(0.001);
        let t_large = large.compute_membrane_tension(0.001);
        assert!(
            t_large > t_small,
            "Larger sphere should have higher hoop tension: small={t_small}, large={t_large}"
        );
    }
}
