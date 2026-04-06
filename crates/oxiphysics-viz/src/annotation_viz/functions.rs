//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BoundingBox2D, PlacedLabel, TextAnchor};

/// Add two 3-vectors.
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract b from a.
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product.
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// L2 norm.
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Normalize a 3-vector. Returns zero vector if norm is tiny.
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}
/// Midpoint of two 3-vectors.
#[inline]
pub(super) fn midpoint3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        (a[0] + b[0]) * 0.5,
        (a[1] + b[1]) * 0.5,
        (a[2] + b[2]) * 0.5,
    ]
}
/// Distance between two 3-vectors.
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3(sub3(a, b))
}
/// Measure the Euclidean distance between two 3D points.
pub fn measure_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    dist3(a, b)
}
/// Measure the angle (in radians) at vertex B in the triangle A-B-C.
pub fn measure_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = normalize3(sub3(a, b));
    let bc = normalize3(sub3(c, b));
    let d = dot3(ba, bc).clamp(-1.0, 1.0);
    d.acos()
}
/// Measure the angle in degrees.
pub fn measure_angle_deg(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    measure_angle(a, b, c).to_degrees()
}
/// Measure the area of a triangle defined by three 3D points.
pub fn measure_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    0.5 * norm3(cross3(ab, ac))
}
/// Measure the area of a polygon defined by coplanar 3D points.
///
/// Uses the shoelace-like formula via cross products.
/// Points should be ordered (CW or CCW).
pub fn measure_polygon_area(points: &[[f64; 3]]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    let mut total = [0.0_f64; 3];
    let n = points.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let c = cross3(points[i], points[j]);
        total[0] += c[0];
        total[1] += c[1];
        total[2] += c[2];
    }
    0.5 * norm3(total)
}
/// Measure the perimeter of a polygon.
pub fn measure_polygon_perimeter(points: &[[f64; 3]]) -> f64 {
    if points.len() < 2 {
        return 0.0;
    }
    let n = points.len();
    let mut perim = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        perim += dist3(points[i], points[j]);
    }
    perim
}
/// Compute the centroid of a set of 3D points.
pub fn compute_centroid(points: &[[f64; 3]]) -> [f64; 3] {
    if points.is_empty() {
        return [0.0; 3];
    }
    let n = points.len() as f64;
    let mut sum = [0.0; 3];
    for p in points {
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
    }
    [sum[0] / n, sum[1] / n, sum[2] / n]
}
/// Snap a 3D point to the nearest grid intersection.
pub fn snap_to_grid(point: [f64; 3], spacing: f64) -> [f64; 3] {
    if spacing <= 0.0 {
        return point;
    }
    [
        (point[0] / spacing).round() * spacing,
        (point[1] / spacing).round() * spacing,
        (point[2] / spacing).round() * spacing,
    ]
}
/// Resolve label collisions using a greedy displacement algorithm.
///
/// Takes a list of ideal label positions and label sizes, and returns adjusted
/// positions that minimise overlap.
///
/// * `positions` — ideal 2D positions for each label center.
/// * `sizes` — `(width, height)` for each label.
/// * `max_displacement` — maximum distance a label can be pushed.
/// * `iterations` — number of relaxation iterations.
pub fn resolve_label_collisions(
    positions: &[[f64; 2]],
    sizes: &[(f64, f64)],
    max_displacement: f64,
    iterations: usize,
) -> Vec<PlacedLabel> {
    assert_eq!(positions.len(), sizes.len());
    let n = positions.len();
    let mut placed: Vec<PlacedLabel> = positions
        .iter()
        .zip(sizes.iter())
        .map(|(&pos, &(w, h))| PlacedLabel {
            position: pos,
            bounds: BoundingBox2D::from_center(pos, w, h),
            ideal_position: pos,
            displaced: false,
        })
        .collect();
    if n <= 1 {
        return placed;
    }
    for _iter in 0..iterations {
        let mut any_moved = false;
        for i in 0..n {
            let mut force = [0.0_f64; 2];
            for j in 0..n {
                if i == j {
                    continue;
                }
                if placed[i].bounds.overlaps(&placed[j].bounds) {
                    let dx = placed[i].position[0] - placed[j].position[0];
                    let dy = placed[i].position[1] - placed[j].position[1];
                    let d = (dx * dx + dy * dy).sqrt().max(1e-10);
                    let ov = placed[i].bounds.overlap_area(&placed[j].bounds).sqrt();
                    force[0] += ov * dx / d;
                    force[1] += ov * dy / d;
                }
            }
            let dx_ideal = placed[i].ideal_position[0] - placed[i].position[0];
            let dy_ideal = placed[i].ideal_position[1] - placed[i].position[1];
            force[0] += dx_ideal * 0.1;
            force[1] += dy_ideal * 0.1;
            let mag = (force[0] * force[0] + force[1] * force[1]).sqrt();
            if mag > 1e-10 {
                let step = mag.min(max_displacement * 0.1);
                let dx = force[0] / mag * step;
                let dy = force[1] / mag * step;
                placed[i].position[0] += dx;
                placed[i].position[1] += dy;
                let total_dx = placed[i].position[0] - placed[i].ideal_position[0];
                let total_dy = placed[i].position[1] - placed[i].ideal_position[1];
                let total_d = (total_dx * total_dx + total_dy * total_dy).sqrt();
                if total_d > max_displacement {
                    let s = max_displacement / total_d;
                    placed[i].position[0] = placed[i].ideal_position[0] + total_dx * s;
                    placed[i].position[1] = placed[i].ideal_position[1] + total_dy * s;
                }
                let (w, h) = sizes[i];
                placed[i].bounds = BoundingBox2D::from_center(placed[i].position, w, h);
                placed[i].displaced = true;
                any_moved = true;
            }
        }
        if !any_moved {
            break;
        }
    }
    placed
}
/// Count the number of overlapping label pairs.
pub fn count_overlaps(labels: &[PlacedLabel]) -> usize {
    let mut count = 0;
    for i in 0..labels.len() {
        for j in (i + 1)..labels.len() {
            if labels[i].bounds.overlaps(&labels[j].bounds) {
                count += 1;
            }
        }
    }
    count
}
/// Estimate the bounding box for a text string.
///
/// Uses a simple fixed-width character model.
///
/// * `text` — the text string.
/// * `font_size` — font size in drawing units.
/// * `char_width_ratio` — ratio of character width to font size (typically 0.6).
pub fn estimate_text_bounds(text: &str, font_size: f64, char_width_ratio: f64) -> BoundingBox2D {
    let lines: Vec<&str> = text.lines().collect();
    let max_chars = lines.iter().map(|l| l.len()).max().unwrap_or(0);
    let width = max_chars as f64 * font_size * char_width_ratio;
    let height = lines.len() as f64 * font_size * 1.2;
    BoundingBox2D::from_center([0.0, 0.0], width, height)
}
/// Position a text bounding box at a given anchor point.
pub fn position_text_at(
    text: &str,
    font_size: f64,
    position: [f64; 2],
    anchor: TextAnchor,
) -> BoundingBox2D {
    let bounds = estimate_text_bounds(text, font_size, 0.6);
    let offset = anchor.offset_to_center(bounds.width(), bounds.height());
    let center = [position[0] + offset[0], position[1] + offset[1]];
    BoundingBox2D::from_center(center, bounds.width(), bounds.height())
}
/// Compute line break positions for wrapping text within a given width.
///
/// Returns a vector of line strings.
pub fn wrap_text(text: &str, max_width: f64, font_size: f64, char_width_ratio: f64) -> Vec<String> {
    let char_width = font_size * char_width_ratio;
    let max_chars = (max_width / char_width).floor() as usize;
    if max_chars == 0 {
        return vec![text.to_string()];
    }
    let mut lines = Vec::new();
    let mut current_line = String::new();
    for word in text.split_whitespace() {
        if current_line.is_empty() {
            if word.len() > max_chars {
                let mut remaining = word;
                while remaining.len() > max_chars {
                    let (chunk, rest) = remaining.split_at(max_chars);
                    lines.push(chunk.to_string());
                    remaining = rest;
                }
                current_line = remaining.to_string();
            } else {
                current_line = word.to_string();
            }
        } else if current_line.len() + 1 + word.len() <= max_chars {
            current_line.push(' ');
            current_line.push_str(word);
        } else {
            lines.push(current_line);
            current_line = word.to_string();
        }
    }
    if !current_line.is_empty() {
        lines.push(current_line);
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation_viz::AnchorStrategy;
    use crate::annotation_viz::AngleArc;
    use crate::annotation_viz::AngleMeasurement;
    use crate::annotation_viz::AnnotationAnchor;
    use crate::annotation_viz::AnnotationLayer;
    use crate::annotation_viz::AnnotationStyle;
    use crate::annotation_viz::AreaMeasurement;
    use crate::annotation_viz::CalloutBox;
    use crate::annotation_viz::CompassRose;
    use crate::annotation_viz::CoordinateAxesWidget;
    use crate::annotation_viz::DimensionLine;
    use crate::annotation_viz::DistanceMeasurement;
    use crate::annotation_viz::GridAnnotation;
    use crate::annotation_viz::LeaderLine;
    use crate::annotation_viz::PolylineAnnotation;
    use crate::annotation_viz::ScaleBar;
    use std::f64::consts::PI;
    #[test]
    fn bbox_width_height() {
        let bb = BoundingBox2D::new([0.0, 0.0], [4.0, 3.0]);
        assert!((bb.width() - 4.0).abs() < 1e-10);
        assert!((bb.height() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn bbox_area() {
        let bb = BoundingBox2D::new([1.0, 2.0], [5.0, 6.0]);
        assert!((bb.area() - 16.0).abs() < 1e-10);
    }
    #[test]
    fn bbox_center() {
        let bb = BoundingBox2D::new([0.0, 0.0], [10.0, 6.0]);
        let c = bb.center();
        assert!((c[0] - 5.0).abs() < 1e-10);
        assert!((c[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn bbox_overlap_true() {
        let a = BoundingBox2D::new([0.0, 0.0], [3.0, 3.0]);
        let b = BoundingBox2D::new([2.0, 2.0], [5.0, 5.0]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn bbox_overlap_false() {
        let a = BoundingBox2D::new([0.0, 0.0], [1.0, 1.0]);
        let b = BoundingBox2D::new([2.0, 2.0], [3.0, 3.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn bbox_contains_point() {
        let bb = BoundingBox2D::new([0.0, 0.0], [4.0, 4.0]);
        assert!(bb.contains_point([2.0, 2.0]));
        assert!(!bb.contains_point([5.0, 5.0]));
    }
    #[test]
    fn bbox_expand() {
        let bb = BoundingBox2D::new([1.0, 1.0], [3.0, 3.0]);
        let expanded = bb.expand(0.5);
        assert!((expanded.min[0] - 0.5).abs() < 1e-10);
        assert!((expanded.max[0] - 3.5).abs() < 1e-10);
    }
    #[test]
    fn bbox_union() {
        let a = BoundingBox2D::new([0.0, 0.0], [2.0, 2.0]);
        let b = BoundingBox2D::new([1.0, 1.0], [4.0, 3.0]);
        let u = a.union(&b);
        assert!((u.min[0]).abs() < 1e-10);
        assert!((u.max[0] - 4.0).abs() < 1e-10);
        assert!((u.max[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn bbox_overlap_area() {
        let a = BoundingBox2D::new([0.0, 0.0], [3.0, 3.0]);
        let b = BoundingBox2D::new([2.0, 2.0], [5.0, 5.0]);
        assert!((a.overlap_area(&b) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn bbox_from_center() {
        let bb = BoundingBox2D::from_center([5.0, 5.0], 4.0, 2.0);
        assert!((bb.min[0] - 3.0).abs() < 1e-10);
        assert!((bb.max[0] - 7.0).abs() < 1e-10);
        assert!((bb.min[1] - 4.0).abs() < 1e-10);
        assert!((bb.max[1] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn dimension_line_length() {
        let dim = DimensionLine::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.3);
        assert!((dim.length() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn dimension_line_direction() {
        let dim = DimensionLine::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], 0.5);
        let dir = dim.direction();
        assert!((dir[0] - 1.0).abs() < 1e-10);
        assert!(dir[1].abs() < 1e-10);
    }
    #[test]
    fn dimension_line_midpoint() {
        let dim = DimensionLine::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 1.0);
        let mid = dim.midpoint();
        assert!((mid[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn dimension_line_with_label() {
        let dim = DimensionLine::with_label([0.0; 3], [1.0, 0.0, 0.0], 0.2, "100mm");
        assert_eq!(dim.display_label(), "100mm");
    }
    #[test]
    fn dimension_line_auto_label() {
        let dim = DimensionLine::new([0.0; 3], [1.0, 0.0, 0.0], 0.2);
        let label = dim.display_label();
        assert!(label.starts_with("1.0"));
    }
    #[test]
    fn leader_line_straight_length() {
        let ll = LeaderLine::new([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((ll.length() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn leader_line_with_bend() {
        let ll = LeaderLine::with_bend([0.0; 3], [2.0, 0.0, 0.0], [2.0, 3.0, 0.0]);
        assert_eq!(ll.segment_count(), 2);
        assert!((ll.length() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn leader_line_path_points() {
        let ll = LeaderLine::with_bend([0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]);
        let pts = ll.path_points();
        assert_eq!(pts.len(), 3);
    }
    #[test]
    fn callout_box_area() {
        let cb = CalloutBox::new([0.0; 3], 4.0, 2.0, "Test");
        assert!((cb.area() - 8.0).abs() < 1e-10);
    }
    #[test]
    fn callout_box_corners() {
        let cb = CalloutBox::new([5.0, 5.0, 0.0], 4.0, 2.0, "Note");
        let corners = cb.corners();
        assert!((corners[0][0] - 3.0).abs() < 1e-10);
        assert!((corners[2][0] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn callout_box_rounded_perimeter() {
        let cb = CalloutBox::rounded([0.0; 3], 4.0, 2.0, "X", 0.3);
        let rect_perim = 2.0 * (4.0 + 2.0);
        assert!(cb.perimeter() < rect_perim);
    }
    #[test]
    fn angle_arc_right_angle() {
        let arc = AngleArc::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!((arc.angle_deg() - 90.0).abs() < 1e-6);
    }
    #[test]
    fn angle_arc_zero_angle() {
        let arc = AngleArc::new([0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        assert!(arc.angle_deg().abs() < 1e-6);
    }
    #[test]
    fn angle_arc_points_count() {
        let arc = AngleArc::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        let pts = arc.arc_points();
        assert_eq!(pts.len(), arc.segments + 1);
    }
    #[test]
    fn angle_arc_length() {
        let arc = AngleArc::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        assert!((arc.arc_length() - PI).abs() < 1e-6);
    }
    #[test]
    fn coord_axes_tips() {
        let axes = CoordinateAxesWidget::new([0.0; 3], 1.0);
        let tips = axes.axis_tips();
        assert!((tips[0][0] - 1.0).abs() < 1e-10);
        assert!((tips[1][1] - 1.0).abs() < 1e-10);
        assert!((tips[2][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn coord_axes_custom_labels() {
        let axes = CoordinateAxesWidget::new([0.0; 3], 1.0).with_labels("R", "S", "T");
        assert_eq!(axes.labels[0], "R");
    }
    #[test]
    fn scale_bar_ticks() {
        let sb = ScaleBar::new([0.0; 3], 10.0, 5, "m");
        let ticks = sb.tick_positions();
        assert_eq!(ticks.len(), 6);
    }
    #[test]
    fn scale_bar_division_width() {
        let sb = ScaleBar::new([0.0; 3], 10.0, 4, "cm");
        assert!((sb.division_width() - 2.5).abs() < 1e-10);
    }
    #[test]
    fn measure_distance_3d() {
        let d = measure_distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-10);
    }
    #[test]
    fn measure_right_angle() {
        let deg = measure_angle_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((deg - 90.0).abs() < 1e-6);
    }
    #[test]
    fn measure_triangle_area_unit() {
        let area = measure_triangle_area([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((area - 0.5).abs() < 1e-10);
    }
    #[test]
    fn measure_polygon_area_square() {
        let sq = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let area = measure_polygon_area(&sq);
        assert!((area - 1.0).abs() < 1e-6);
    }
    #[test]
    fn measure_polygon_perimeter_square() {
        let sq = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let perim = measure_polygon_perimeter(&sq);
        assert!((perim - 4.0).abs() < 1e-10);
    }
    #[test]
    fn anchor_fixed_does_not_move() {
        let mut anchor = AnnotationAnchor::new([0.0; 3], [1.0, 1.0, 0.0], AnchorStrategy::Fixed);
        anchor.update([5.0, 5.0, 0.0]);
        assert!((anchor.label_position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn anchor_follow_point() {
        let mut anchor =
            AnnotationAnchor::new([0.0; 3], [1.0, 1.0, 0.0], AnchorStrategy::FollowPoint);
        anchor.update([5.0, 5.0, 0.0]);
        assert!((anchor.label_position[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn anchor_fixed_offset() {
        let mut anchor =
            AnnotationAnchor::new([0.0; 3], [2.0, 3.0, 0.0], AnchorStrategy::FixedOffset);
        anchor.update([1.0, 1.0, 0.0]);
        assert!((anchor.label_position[0] - 3.0).abs() < 1e-10);
        assert!((anchor.label_position[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn collision_avoidance_non_overlapping() {
        let positions = [[0.0, 0.0], [10.0, 10.0]];
        let sizes = [(1.0, 1.0), (1.0, 1.0)];
        let placed = resolve_label_collisions(&positions, &sizes, 5.0, 10);
        assert_eq!(placed.len(), 2);
        assert_eq!(count_overlaps(&placed), 0);
    }
    #[test]
    fn collision_avoidance_reduces_overlap() {
        let positions = [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0]];
        let sizes = [(2.0, 1.0), (2.0, 1.0), (2.0, 1.0)];
        let placed = resolve_label_collisions(&positions, &sizes, 10.0, 50);
        let initial_overlaps = {
            let initial = resolve_label_collisions(&positions, &sizes, 0.0, 0);
            count_overlaps(&initial)
        };
        let final_overlaps = count_overlaps(&placed);
        assert!(
            final_overlaps <= initial_overlaps,
            "should reduce overlaps: {} -> {}",
            initial_overlaps,
            final_overlaps
        );
    }
    #[test]
    fn annotation_layer_count() {
        let mut layer = AnnotationLayer::new();
        layer.add_dimension(DimensionLine::new([0.0; 3], [1.0, 0.0, 0.0], 0.1));
        layer.add_leader(LeaderLine::new([0.0; 3], [1.0; 3]));
        layer.add_callout(CalloutBox::new([0.0; 3], 2.0, 1.0, "X"));
        assert_eq!(layer.annotation_count(), 3);
    }
    #[test]
    fn annotation_layer_clear() {
        let mut layer = AnnotationLayer::new();
        layer.add_dimension(DimensionLine::new([0.0; 3], [1.0, 0.0, 0.0], 0.1));
        layer.clear();
        assert_eq!(layer.annotation_count(), 0);
    }
    #[test]
    fn estimate_text_bounds_non_empty() {
        let bb = estimate_text_bounds("Hello", 12.0, 0.6);
        assert!(bb.width() > 0.0);
        assert!(bb.height() > 0.0);
    }
    #[test]
    fn wrap_text_single_line() {
        let lines = wrap_text("short", 100.0, 10.0, 0.6);
        assert_eq!(lines.len(), 1);
    }
    #[test]
    fn wrap_text_multiple_lines() {
        let lines = wrap_text("this is a longer text that should wrap", 50.0, 10.0, 0.6);
        assert!(lines.len() >= 2);
    }
    #[test]
    fn annotation_style_default() {
        let style = AnnotationStyle::default();
        assert!(style.line_width > 0.0);
        assert!(style.font_size > 0.0);
    }
    #[test]
    fn annotation_style_arrow_area() {
        let style = AnnotationStyle::default();
        assert!(style.arrow_area() > 0.0);
    }
    #[test]
    fn distance_measurement_display() {
        let dm = DistanceMeasurement::new([0.0; 3], [5.0, 0.0, 0.0], 0.3, "mm");
        assert!(dm.display().contains("5.000"));
        assert!(dm.display().contains("mm"));
    }
    #[test]
    fn angle_measurement_display() {
        let am = AngleMeasurement::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!(am.display().contains("90.0"));
    }
    #[test]
    fn area_measurement_display() {
        let am = AreaMeasurement::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            "m",
        );
        let txt = am.display();
        assert!(txt.contains("m"));
    }
    #[test]
    fn polyline_length_open() {
        let pa = PolylineAnnotation::new(vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], false);
        assert!((pa.length() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn polyline_length_closed() {
        let pa = PolylineAnnotation::new(vec![[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], true);
        let expected = 2.0 + (2.0_f64).sqrt();
        assert!((pa.length() - expected).abs() < 1e-10);
    }
    #[test]
    fn grid_annotation_cells() {
        let grid = GridAnnotation::new([0.0; 3], 10.0, 10.0, 5, 5);
        assert_eq!(grid.cell_count(), 25);
    }
    #[test]
    fn grid_annotation_intersections() {
        let grid = GridAnnotation::new([0.0; 3], 10.0, 10.0, 4, 3);
        assert_eq!(grid.intersection_count(), 20);
    }
    #[test]
    fn compass_rose_direction_count() {
        let cr = CompassRose::new([0.0; 3], 1.0);
        let lines = cr.direction_lines();
        assert_eq!(lines.len(), 8);
    }
    #[test]
    fn compass_rose_north_tip() {
        let cr = CompassRose::new([0.0; 3], 1.0);
        let tip = cr.direction_point(0.0);
        assert!((tip[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn snap_to_grid_basic() {
        let snapped = snap_to_grid([0.7, 1.3, 2.6], 1.0);
        assert!((snapped[0] - 1.0).abs() < 1e-10);
        assert!((snapped[1] - 1.0).abs() < 1e-10);
        assert!((snapped[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn text_anchor_center_offset_zero() {
        let offset = TextAnchor::Center.offset_to_center(10.0, 5.0);
        assert!(offset[0].abs() < 1e-10);
        assert!(offset[1].abs() < 1e-10);
    }
    #[test]
    fn text_anchor_top_left_offset() {
        let offset = TextAnchor::TopLeft.offset_to_center(10.0, 6.0);
        assert!((offset[0] - 5.0).abs() < 1e-10);
        assert!((offset[1] - (-3.0)).abs() < 1e-10);
    }
}
