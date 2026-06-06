//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ClippingPlane, CoupledArrowConfig, EmFieldLine, EmFieldLineConfig, FieldLineMethod, FieldNode,
    FsiSnapshot, InterfacePoint, LineSegment, MeshInterpolationMethod, Rgba, TempColorKey,
    TempColorRamp, ThermalMechanicalOverlayConfig,
};

/// Add two 3-D vectors.
#[inline]
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-D vectors.
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-D vector by a scalar.
#[inline]
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product of two 3-D vectors.
#[inline]
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean length of a 3-D vector.
#[inline]
pub(super) fn vec3_len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
/// Normalize a 3-D vector; returns `[0,0,0]` for near-zero vectors.
#[inline]
pub(super) fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}
/// Linear interpolation between two 3-D vectors.
#[inline]
pub(super) fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
/// Map a normalized `t` in `[0,1]` to a hot colormap (black → red → yellow → white).
pub(super) fn hot_colormap(t: f64) -> Rgba {
    let t = t.clamp(0.0, 1.0) as f32;
    let r = (t * 3.0).min(1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}
/// Map a normalized `t` in `[0,1]` to a cool colormap (cyan → magenta).
pub(super) fn cool_colormap(t: f64) -> Rgba {
    let t = t.clamp(0.0, 1.0) as f32;
    Rgba::new(t, 1.0 - t, 1.0, 1.0)
}
/// Normalize a scalar to `[0,1]` given min/max.
pub(super) fn normalize_scalar(val: f64, min: f64, max: f64) -> f64 {
    if (max - min).abs() < 1e-15 {
        0.5
    } else {
        ((val - min) / (max - min)).clamp(0.0, 1.0)
    }
}
/// Compute per-node colors for a coupled thermal-mechanical overlay.
///
/// For each node the thermal color (hot colormap from temperature) and
/// mechanical color (cool colormap from displacement magnitude) are blended
/// according to `config.thermal_weight`.
pub fn thermal_mechanical_overlay(
    temperature_nodes: &[FieldNode],
    displacement_nodes: &[FieldNode],
    config: &ThermalMechanicalOverlayConfig,
) -> Vec<Rgba> {
    let n = temperature_nodes.len().min(displacement_nodes.len());
    let mut colors = Vec::with_capacity(n);
    for i in 0..n {
        let t_norm = normalize_scalar(
            temperature_nodes[i].scalar,
            config.temp_min,
            config.temp_max,
        );
        let d_norm = normalize_scalar(
            vec3_len(displacement_nodes[i].vector),
            config.disp_min,
            config.disp_max,
        );
        let thermal_color = hot_colormap(t_norm);
        let mech_color = cool_colormap(d_norm);
        let w = config.thermal_weight.clamp(0.0, 1.0) as f32;
        colors.push(Rgba::lerp(mech_color, thermal_color, w));
    }
    colors
}
/// Generate arrow [`LineSegment`]s for two coupled vector fields at the same nodes.
///
/// Returns a flat list: for each node two segments are generated — one for
/// `heat_flux` (orange) and one for `displacement` (blue), filtered by
/// `config.min_magnitude`.
pub fn coupled_field_arrows(
    nodes: &[FieldNode],
    heat_flux: &[[f64; 3]],
    displacement: &[[f64; 3]],
    config: &CoupledArrowConfig,
) -> Vec<LineSegment> {
    let n = nodes.len().min(heat_flux.len()).min(displacement.len());
    let mut segments = Vec::with_capacity(n * 2);
    for i in 0..n {
        let origin = nodes[i].position;
        let hf = heat_flux[i];
        if vec3_len(hf) >= config.min_magnitude {
            let tip = vec3_add(origin, vec3_scale(hf, config.heat_flux_scale));
            segments.push(LineSegment::new(origin, tip, config.heat_flux_color));
        }
        let disp = displacement[i];
        if vec3_len(disp) >= config.min_magnitude {
            let tip = vec3_add(origin, vec3_scale(disp, config.displacement_scale));
            segments.push(LineSegment::new(origin, tip, config.displacement_color));
        }
    }
    segments
}
/// Extract interface points from a scalar volume field using marching-squares
/// in 1-D stencil (simplified): a point is on the interface if its scalar
/// value straddles `iso_value` between adjacent samples.
///
/// `scalars` and `positions` must have the same length.  Returns the
/// interface samples sorted by position along the first coordinate axis.
pub fn extract_interface_points(
    positions: &[[f64; 3]],
    scalars: &[f64],
    iso_value: f64,
) -> Vec<InterfacePoint> {
    let n = positions.len().min(scalars.len());
    let mut points = Vec::new();
    if n < 2 {
        return points;
    }
    for i in 0..n - 1 {
        let s0 = scalars[i];
        let s1 = scalars[i + 1];
        let crosses = (s0 < iso_value) != (s1 < iso_value);
        if crosses {
            let t = (iso_value - s0) / (s1 - s0 + 1e-300);
            let pos = vec3_lerp(positions[i], positions[i + 1], t);
            let dir = vec3_sub(positions[i + 1], positions[i]);
            let normal = vec3_norm([dir[1], -dir[0], dir[2]]);
            let phase = if s0 < iso_value { 0 } else { 1 };
            points.push(InterfacePoint::new(pos, normal, phase));
        }
    }
    points
}
/// Generate [`LineSegment`]s that display phase interface normals as short ticks.
pub fn interface_normal_lines(
    interface: &[InterfacePoint],
    tick_length: f64,
    color_a: Rgba,
    color_b: Rgba,
) -> Vec<LineSegment> {
    interface
        .iter()
        .map(|ip| {
            let tip = vec3_add(ip.position, vec3_scale(ip.normal, tick_length));
            let color = if ip.phase == 0 { color_a } else { color_b };
            LineSegment::new(ip.position, tip, color)
        })
        .collect()
}
/// Generate velocity arrow lines for the fluid domain.
pub fn fsi_fluid_arrows(
    snapshot: &FsiSnapshot,
    scale: f64,
    min_mag: f64,
    color: Rgba,
) -> Vec<LineSegment> {
    snapshot
        .fluid_positions
        .iter()
        .zip(snapshot.fluid_velocities.iter())
        .filter(|(_, v)| vec3_len(**v) >= min_mag)
        .map(|(p, v)| LineSegment::new(*p, vec3_add(*p, vec3_scale(*v, scale)), color))
        .collect()
}
/// Generate displacement arrow lines for the structural domain.
pub fn fsi_structural_arrows(
    snapshot: &FsiSnapshot,
    scale: f64,
    min_mag: f64,
    color: Rgba,
) -> Vec<LineSegment> {
    snapshot
        .structural_positions
        .iter()
        .zip(snapshot.structural_displacements.iter())
        .filter(|(_, d)| vec3_len(**d) >= min_mag)
        .map(|(p, d)| LineSegment::new(*p, vec3_add(*p, vec3_scale(*d, scale)), color))
        .collect()
}
/// Generate coupling lines between FSI interface nodes and their structural counterparts.
///
/// Uses the `fsi_interface_indices` in the snapshot to pair fluid boundary
/// nodes with structural nodes (by matching index into `structural_positions`).
pub fn fsi_coupling_lines(snapshot: &FsiSnapshot, color: Rgba) -> Vec<LineSegment> {
    snapshot
        .fsi_interface_indices
        .iter()
        .filter_map(|&fi| {
            let fp = snapshot.fluid_positions.get(fi)?;
            let sp = snapshot.structural_positions.get(fi)?;
            Some(LineSegment::new(*fp, *sp, color))
        })
        .collect()
}
/// Trace one electromagnetic field line starting from `seed`.
///
/// `field_fn` takes a position and returns the field vector at that point.
pub fn trace_em_field_line<F>(
    seed: [f64; 3],
    field_fn: &F,
    config: &EmFieldLineConfig,
) -> EmFieldLine
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    let mut line = EmFieldLine::new();
    let mut pos = seed;
    line.points.push(pos);
    line.magnitudes.push(vec3_len(field_fn(pos)));
    for _ in 0..config.max_steps {
        let b = field_fn(pos);
        let mag = vec3_len(b);
        if mag < config.min_magnitude {
            break;
        }
        let dir = vec3_norm(b);
        pos = match config.method {
            FieldLineMethod::Euler => vec3_add(pos, vec3_scale(dir, config.step_size)),
            FieldLineMethod::Rk4 => {
                let h = config.step_size;
                let k1 = vec3_norm(field_fn(pos));
                let k2 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k1, h * 0.5))));
                let k3 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k2, h * 0.5))));
                let k4 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k3, h))));
                let avg = [
                    (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
                    (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
                    (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
                ];
                vec3_add(pos, vec3_scale(avg, h))
            }
        };
        let b_new = field_fn(pos);
        line.points.push(pos);
        line.magnitudes.push(vec3_len(b_new));
    }
    line
}
/// Trace multiple field lines from a list of seed points.
pub fn trace_em_field_lines<F>(
    seeds: &[[f64; 3]],
    field_fn: &F,
    config: &EmFieldLineConfig,
) -> Vec<EmFieldLine>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    seeds
        .iter()
        .map(|&s| trace_em_field_line(s, field_fn, config))
        .collect()
}
/// Convert traced field lines to [`LineSegment`]s, colored by local magnitude.
///
/// `mag_min` and `mag_max` define the normalization range for coloring;
/// the cool colormap (cyan–magenta) is used.
pub fn em_field_lines_to_segments(
    lines: &[EmFieldLine],
    mag_min: f64,
    mag_max: f64,
) -> Vec<LineSegment> {
    let mut segs = Vec::new();
    for line in lines {
        let pts = &line.points;
        let mags = &line.magnitudes;
        for i in 0..pts.len().saturating_sub(1) {
            let t = normalize_scalar(mags[i], mag_min, mag_max);
            let color = cool_colormap(t);
            segs.push(LineSegment::new(pts[i], pts[i + 1], color));
        }
    }
    segs
}
/// Build a standard iron-body black-body temperature ramp
/// (black → red → orange → yellow → white).
pub fn blackbody_ramp() -> TempColorRamp {
    TempColorRamp::sorted(vec![
        TempColorKey::new(0.0, Rgba::new(0.0, 0.0, 0.0, 1.0)),
        TempColorKey::new(300.0, Rgba::new(0.1, 0.0, 0.0, 1.0)),
        TempColorKey::new(800.0, Rgba::new(0.8, 0.1, 0.0, 1.0)),
        TempColorKey::new(1200.0, Rgba::new(1.0, 0.5, 0.0, 1.0)),
        TempColorKey::new(2000.0, Rgba::new(1.0, 1.0, 0.5, 1.0)),
        TempColorKey::new(5000.0, Rgba::new(1.0, 1.0, 1.0, 1.0)),
    ])
}
/// Transfer a scalar field from a source mesh to a target mesh.
///
/// `source_positions` and `source_values` define the source field.
/// `target_positions` defines the query points.
/// Returns one interpolated scalar per target node.
pub fn interpolate_scalar_field(
    source_positions: &[[f64; 3]],
    source_values: &[f64],
    target_positions: &[[f64; 3]],
    method: MeshInterpolationMethod,
) -> Vec<f64> {
    let n_src = source_positions.len().min(source_values.len());
    target_positions
        .iter()
        .map(|&qp| match method {
            MeshInterpolationMethod::NearestNeighbor => {
                let mut best_dist = f64::INFINITY;
                let mut best_val = 0.0;
                for i in 0..n_src {
                    let d = vec3_len(vec3_sub(source_positions[i], qp));
                    if d < best_dist {
                        best_dist = d;
                        best_val = source_values[i];
                    }
                }
                best_val
            }
            MeshInterpolationMethod::InverseDistance { power, k_neighbors } => {
                let mut dist_val: Vec<(f64, f64)> = (0..n_src)
                    .map(|i| {
                        (
                            vec3_len(vec3_sub(source_positions[i], qp)),
                            source_values[i],
                        )
                    })
                    .collect();
                dist_val.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                let k = k_neighbors.min(n_src).max(1);
                let mut w_sum = 0.0;
                let mut v_sum = 0.0;
                for (d, v) in dist_val.iter().take(k) {
                    if *d < 1e-14 {
                        return *v;
                    }
                    let w = 1.0 / d.powf(power);
                    w_sum += w;
                    v_sum += w * v;
                }
                if w_sum < 1e-300 { 0.0 } else { v_sum / w_sum }
            }
            MeshInterpolationMethod::Barycentric => {
                let mut dist_val: Vec<(f64, f64)> = (0..n_src)
                    .map(|i| {
                        (
                            vec3_len(vec3_sub(source_positions[i], qp)),
                            source_values[i],
                        )
                    })
                    .collect();
                dist_val.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                let k = 3_usize.min(n_src);
                let mut w_sum = 0.0;
                let mut v_sum = 0.0;
                for (d, v) in dist_val.iter().take(k) {
                    if *d < 1e-14 {
                        return *v;
                    }
                    let w = 1.0 / (d * d);
                    w_sum += w;
                    v_sum += w * v;
                }
                if w_sum < 1e-300 { 0.0 } else { v_sum / w_sum }
            }
        })
        .collect()
}
/// Transfer a vector field from a source mesh to a target mesh using IDW-2.
pub fn interpolate_vector_field(
    source_positions: &[[f64; 3]],
    source_vectors: &[[f64; 3]],
    target_positions: &[[f64; 3]],
    k_neighbors: usize,
) -> Vec<[f64; 3]> {
    let n_src = source_positions.len().min(source_vectors.len());
    target_positions
        .iter()
        .map(|&qp| {
            let mut dist_vec: Vec<(f64, [f64; 3])> = (0..n_src)
                .map(|i| {
                    (
                        vec3_len(vec3_sub(source_positions[i], qp)),
                        source_vectors[i],
                    )
                })
                .collect();
            dist_vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let k = k_neighbors.min(n_src).max(1);
            let mut w_sum = 0.0;
            let mut v_sum = [0.0_f64; 3];
            for (d, v) in dist_vec.iter().take(k) {
                if *d < 1e-14 {
                    return *v;
                }
                let w = 1.0 / (d * d);
                w_sum += w;
                v_sum = vec3_add(v_sum, vec3_scale(*v, w));
            }
            if w_sum < 1e-300 {
                [0.0; 3]
            } else {
                vec3_scale(v_sum, 1.0 / w_sum)
            }
        })
        .collect()
}
/// Generate a color-coded scalar field visualization as colored [`LineSegment`]s
/// (vertical bars from node to colored tip), useful for 2-D scalar profiles.
pub fn scalar_profile_bars(
    positions: &[[f64; 3]],
    scalars: &[f64],
    min_val: f64,
    max_val: f64,
    bar_height_scale: f64,
) -> Vec<LineSegment> {
    positions
        .iter()
        .zip(scalars.iter())
        .map(|(&pos, &val)| {
            let t = normalize_scalar(val, min_val, max_val);
            let color = hot_colormap(t);
            let tip = [pos[0], pos[1] + t * bar_height_scale, pos[2]];
            LineSegment::new(pos, tip, color)
        })
        .collect()
}
/// Generate iso-contour lines in 2-D (XZ plane, y=0) using simple marching squares
/// on a regular grid.
///
/// `grid_values` is row-major with `nx` columns and `nz` rows.
/// `cell_size` is the spacing between grid points.
/// `origin` is the bottom-left corner of the grid.
pub fn iso_contour_lines(
    grid_values: &[f64],
    nx: usize,
    nz: usize,
    iso_value: f64,
    cell_size: f64,
    origin: [f64; 3],
    color: Rgba,
) -> Vec<LineSegment> {
    let mut segs = Vec::new();
    if nx < 2 || nz < 2 {
        return segs;
    }
    let idx = |ix: usize, iz: usize| iz * nx + ix;
    for iz in 0..nz - 1 {
        for ix in 0..nx - 1 {
            let s00 = *grid_values.get(idx(ix, iz)).unwrap_or(&0.0);
            let s10 = *grid_values.get(idx(ix + 1, iz)).unwrap_or(&0.0);
            let s01 = *grid_values.get(idx(ix, iz + 1)).unwrap_or(&0.0);
            let s11 = *grid_values.get(idx(ix + 1, iz + 1)).unwrap_or(&0.0);
            let above = [
                s00 >= iso_value,
                s10 >= iso_value,
                s11 >= iso_value,
                s01 >= iso_value,
            ];
            let code: u8 = (above[0] as u8)
                | ((above[1] as u8) << 1)
                | ((above[2] as u8) << 2)
                | ((above[3] as u8) << 3);
            if code == 0 || code == 15 {
                continue;
            }
            let x0 = origin[0] + ix as f64 * cell_size;
            let z0 = origin[2] + iz as f64 * cell_size;
            let x1 = x0 + cell_size;
            let z1 = z0 + cell_size;
            let y = origin[1];
            let interp = |a: f64, b: f64, pa: [f64; 3], pb: [f64; 3]| -> [f64; 3] {
                let denom = b - a;
                if denom.abs() < 1e-14 {
                    return pa;
                }
                let t = (iso_value - a) / denom;
                vec3_lerp(pa, pb, t)
            };
            let e_bottom = interp(s00, s10, [x0, y, z0], [x1, y, z0]);
            let e_right = interp(s10, s11, [x1, y, z0], [x1, y, z1]);
            let e_top = interp(s11, s01, [x1, y, z1], [x0, y, z1]);
            let e_left = interp(s01, s00, [x0, y, z1], [x0, y, z0]);
            let push_seg = |segs: &mut Vec<LineSegment>, a: [f64; 3], b: [f64; 3]| {
                segs.push(LineSegment::new(a, b, color));
            };
            match code {
                1 | 14 => push_seg(&mut segs, e_bottom, e_left),
                2 | 13 => push_seg(&mut segs, e_bottom, e_right),
                3 | 12 => push_seg(&mut segs, e_left, e_right),
                4 | 11 => push_seg(&mut segs, e_right, e_top),
                5 => {
                    push_seg(&mut segs, e_bottom, e_right);
                    push_seg(&mut segs, e_left, e_top);
                }
                6 | 9 => push_seg(&mut segs, e_bottom, e_top),
                7 | 8 => push_seg(&mut segs, e_left, e_top),
                10 => {
                    push_seg(&mut segs, e_bottom, e_left);
                    push_seg(&mut segs, e_right, e_top);
                }
                _ => {}
            }
        }
    }
    segs
}
/// Compute the multiphysics residual norm — the maximum absolute residual across
/// multiple field arrays.
///
/// Used to assess overall coupling convergence.
pub fn multiphysics_residual_norm(fields: &[&[f64]]) -> f64 {
    fields
        .iter()
        .flat_map(|f| f.iter())
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max)
}
/// Compute per-node coupling error between two sets of field values.
///
/// Returns the absolute pointwise difference.
pub fn coupling_error(old_values: &[f64], new_values: &[f64]) -> Vec<f64> {
    old_values
        .iter()
        .zip(new_values.iter())
        .map(|(a, b)| (b - a).abs())
        .collect()
}
/// Smooth a 1-D scalar field using a simple box filter of half-width `radius`.
pub fn smooth_scalar_field(values: &[f64], radius: usize) -> Vec<f64> {
    let n = values.len();
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(radius);
            let hi = (i + radius + 1).min(n);
            let sum: f64 = values[lo..hi].iter().sum();
            sum / (hi - lo) as f64
        })
        .collect()
}
/// Apply a clipping plane to a scalar field, returning only the visible nodes
/// and their values.
pub fn clip_scalar_field<'a>(
    nodes: &'a [[f64; 3]],
    values: &'a [f64],
    plane: &ClippingPlane,
) -> (Vec<&'a [f64; 3]>, Vec<f64>) {
    let mut out_nodes = Vec::new();
    let mut out_vals = Vec::new();
    for (p, v) in nodes.iter().zip(values.iter()) {
        if plane.visible(*p) {
            out_nodes.push(p);
            out_vals.push(*v);
        }
    }
    (out_nodes, out_vals)
}
