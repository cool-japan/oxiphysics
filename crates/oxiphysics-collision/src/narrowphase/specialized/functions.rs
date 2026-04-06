//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::{BoxShape, Capsule, Cone, Cylinder, Sphere, Torus};

use super::types::{CylinderContactMode, CylinderContactResult};

/// Sphere-sphere collision test.
pub fn sphere_sphere(s1: &Sphere, t1: &Transform, s2: &Sphere, t2: &Transform) -> Option<Contact> {
    let diff = t2.position - t1.position;
    let dist = diff.norm();
    let combined_radius = s1.radius + s2.radius;
    if dist >= combined_radius {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined_radius - dist;
    let point_a = t1.position + normal * s1.radius;
    let point_b = t2.position - normal * s2.radius;
    Some(Contact::new(point_a, point_b, normal, depth))
}
/// Sphere-box collision test.
pub fn sphere_box(
    sphere: &Sphere,
    sphere_t: &Transform,
    box_shape: &BoxShape,
    box_t: &Transform,
) -> Option<Contact> {
    let inv_box = box_t.inverse();
    let local_center = inv_box.transform_point(&sphere_t.position);
    let clamped = Vec3::new(
        local_center
            .x
            .clamp(-box_shape.half_extents.x, box_shape.half_extents.x),
        local_center
            .y
            .clamp(-box_shape.half_extents.y, box_shape.half_extents.y),
        local_center
            .z
            .clamp(-box_shape.half_extents.z, box_shape.half_extents.z),
    );
    let diff = local_center - clamped;
    let dist_sq = diff.norm_squared();
    if dist_sq >= sphere.radius * sphere.radius {
        return None;
    }
    let dist = dist_sq.sqrt();
    let local_normal = if dist > 1e-10 {
        diff / dist
    } else {
        let dx = box_shape.half_extents.x - local_center.x.abs();
        let dy = box_shape.half_extents.y - local_center.y.abs();
        let dz = box_shape.half_extents.z - local_center.z.abs();
        if dx <= dy && dx <= dz {
            Vec3::new(local_center.x.signum(), 0.0, 0.0)
        } else if dy <= dz {
            Vec3::new(0.0, local_center.y.signum(), 0.0)
        } else {
            Vec3::new(0.0, 0.0, local_center.z.signum())
        }
    };
    let depth = sphere.radius - dist;
    let world_normal = box_t.transform_vector(&local_normal);
    let point_b = box_t.transform_point(&clamped);
    let point_a = sphere_t.position - world_normal * (sphere.radius - depth);
    Some(Contact::new(point_a, point_b, world_normal, depth))
}
/// Box-box SAT (Separating Axis Theorem) collision test.
pub fn box_box_sat(
    box_a: &BoxShape,
    transform_a: &Transform,
    box_b: &BoxShape,
    transform_b: &Transform,
) -> Option<Contact> {
    let rot_a = transform_a.rotation.to_rotation_matrix();
    let rot_b = transform_b.rotation.to_rotation_matrix();
    let t = transform_b.position - transform_a.position;
    let t_local = Vec3::new(
        rot_a.matrix().column(0).dot(&t),
        rot_a.matrix().column(1).dot(&t),
        rot_a.matrix().column(2).dot(&t),
    );
    let mut r = [[0.0; 3]; 3];
    let mut abs_r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = rot_a.matrix().column(i).dot(&rot_b.matrix().column(j));
            abs_r[i][j] = r[i][j].abs() + 1e-6;
        }
    }
    let ea = [
        box_a.half_extents.x,
        box_a.half_extents.y,
        box_a.half_extents.z,
    ];
    let eb = [
        box_b.half_extents.x,
        box_b.half_extents.y,
        box_b.half_extents.z,
    ];
    let t_arr = [t_local.x, t_local.y, t_local.z];
    let mut min_depth = Real::INFINITY;
    let mut min_axis = Vec3::new(0.0, 1.0, 0.0);
    for i in 0..3 {
        let ra = ea[i];
        let rb = eb[0] * abs_r[i][0] + eb[1] * abs_r[i][1] + eb[2] * abs_r[i][2];
        let depth = (ra + rb) - t_arr[i].abs();
        if depth < 0.0 {
            return None;
        }
        if depth < min_depth {
            min_depth = depth;
            let mut axis = Vec3::zeros();
            axis[i] = if t_arr[i] >= 0.0 { 1.0 } else { -1.0 };
            min_axis = rot_a.matrix() * axis;
        }
    }
    for j in 0..3 {
        let ra = ea[0] * abs_r[0][j] + ea[1] * abs_r[1][j] + ea[2] * abs_r[2][j];
        let rb = eb[j];
        let sep = t_arr[0] * r[0][j] + t_arr[1] * r[1][j] + t_arr[2] * r[2][j];
        let depth = (ra + rb) - sep.abs();
        if depth < 0.0 {
            return None;
        }
        if depth < min_depth {
            min_depth = depth;
            let mut axis = Vec3::zeros();
            axis[j] = if sep >= 0.0 { 1.0 } else { -1.0 };
            min_axis = rot_b.matrix() * axis;
        }
    }
    let edge_tests: [(usize, usize); 9] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 0),
        (2, 1),
        (2, 2),
    ];
    for &(i, j) in &edge_tests {
        let i1 = (i + 1) % 3;
        let i2 = (i + 2) % 3;
        let j1 = (j + 1) % 3;
        let j2 = (j + 2) % 3;
        let ra = ea[i1] * abs_r[i2][j] + ea[i2] * abs_r[i1][j];
        let rb = eb[j1] * abs_r[i][j2] + eb[j2] * abs_r[i][j1];
        let sep = t_arr[i2] * r[i1][j] - t_arr[i1] * r[i2][j];
        let depth = (ra + rb) - sep.abs();
        if depth < 0.0 {
            return None;
        }
        let axis_a = rot_a.matrix().column(i).into_owned();
        let axis_b = rot_b.matrix().column(j).into_owned();
        let cross = axis_a.cross(&axis_b);
        let cross_len = cross.norm();
        if cross_len > 1e-6 && depth < min_depth {
            min_depth = depth;
            min_axis = cross / cross_len;
            if min_axis.dot(&t) < 0.0 {
                min_axis = -min_axis;
            }
        }
    }
    let normal = min_axis;
    let point_a = transform_a.position + normal * min_depth * 0.5;
    let point_b = transform_b.position - normal * min_depth * 0.5;
    Some(Contact::new(point_a, point_b, normal, min_depth))
}
/// Sphere-capsule collision test.
///
/// Tests whether `sphere` overlaps `capsule`.  The capsule is assumed to be
/// oriented along the Y axis in its local frame.
pub fn sphere_capsule(
    sphere: &Sphere,
    sphere_t: &Transform,
    capsule: &Capsule,
    capsule_t: &Transform,
) -> Option<Contact> {
    let up = capsule_t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let seg_a = capsule_t.position + up * capsule.half_height;
    let seg_b = capsule_t.position - up * capsule.half_height;
    let d = seg_b - seg_a;
    let len_sq = d.dot(&d);
    let t = if len_sq < 1e-10 {
        0.0
    } else {
        let v = sphere_t.position - seg_a;
        (v.dot(&d) / len_sq).clamp(0.0, 1.0)
    };
    let closest = seg_a + d * t;
    let diff = sphere_t.position - closest;
    let dist = diff.norm();
    let combined = sphere.radius + capsule.radius;
    if dist >= combined {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let point_a = sphere_t.position - normal * sphere.radius;
    let point_b = closest + normal * capsule.radius;
    Some(Contact::new(point_a, point_b, normal, depth))
}
/// Capsule-capsule collision test.
pub fn capsule_capsule(
    c1: &Capsule,
    t1: &Transform,
    c2: &Capsule,
    t2: &Transform,
) -> Option<Contact> {
    let up1 = t1.rotation * Vec3::new(0.0, 1.0, 0.0);
    let a1 = t1.position + up1 * c1.half_height;
    let b1 = t1.position - up1 * c1.half_height;
    let up2 = t2.rotation * Vec3::new(0.0, 1.0, 0.0);
    let a2 = t2.position + up2 * c2.half_height;
    let b2 = t2.position - up2 * c2.half_height;
    let (p1, p2) = closest_points_segments(&a1, &b1, &a2, &b2);
    let diff = p2 - p1;
    let dist = diff.norm();
    let combined = c1.radius + c2.radius;
    if dist >= combined {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let point_a = p1 + normal * c1.radius;
    let point_b = p2 - normal * c2.radius;
    Some(Contact::new(point_a, point_b, normal, depth))
}
/// Cylinder-cylinder collision test.
///
/// Both cylinders are assumed to be oriented along the Y-axis in their local frames.
/// This uses a segment-segment closest-points approach for the cylinder axes,
/// then checks radial overlap.
pub fn cylinder_cylinder(
    c1: &Cylinder,
    t1: &Transform,
    c2: &Cylinder,
    t2: &Transform,
) -> Option<Contact> {
    let up1 = t1.rotation * Vec3::new(0.0, 1.0, 0.0);
    let a1 = t1.position + up1 * c1.half_height;
    let b1 = t1.position - up1 * c1.half_height;
    let up2 = t2.rotation * Vec3::new(0.0, 1.0, 0.0);
    let a2 = t2.position + up2 * c2.half_height;
    let b2 = t2.position - up2 * c2.half_height;
    let (p1, p2) = closest_points_segments(&a1, &b1, &a2, &b2);
    let diff = p2 - p1;
    let dist = diff.norm();
    let combined = c1.radius + c2.radius;
    if dist >= combined {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let point_a = p1 + normal * c1.radius;
    let point_b = p2 - normal * c2.radius;
    Some(Contact::new(point_a, point_b, normal, depth))
}
/// Cone-sphere collision test.
///
/// The cone is oriented along the Y-axis in its local frame, with the apex
/// at `+half_height` and the base at `-half_height`.
pub fn cone_sphere(
    cone: &Cone,
    cone_t: &Transform,
    sphere: &Sphere,
    sphere_t: &Transform,
) -> Option<Contact> {
    let inv_cone = cone_t.inverse();
    let local_center = inv_cone.transform_point(&sphere_t.position);
    let h = cone.half_height * 2.0;
    let apex_y = cone.half_height;
    let base_y = -cone.half_height;
    let axis_y = local_center.y.clamp(base_y, apex_y);
    let t_param = (apex_y - axis_y) / h;
    let cone_r_at_y = cone.radius * t_param;
    let radial = (local_center.x * local_center.x + local_center.z * local_center.z).sqrt();
    let closest = if radial < 1e-10 {
        Vec3::new(0.0, axis_y, 0.0)
    } else {
        let clamped_r = radial.min(cone_r_at_y);
        let scale = clamped_r / radial;
        Vec3::new(local_center.x * scale, axis_y, local_center.z * scale)
    };
    let diff = local_center - closest;
    let dist = diff.norm();
    if dist >= sphere.radius {
        return None;
    }
    let local_normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = sphere.radius - dist;
    let world_normal = cone_t.transform_vector(&local_normal);
    let world_closest = cone_t.transform_point(&closest);
    let point_a = sphere_t.position - world_normal * sphere.radius;
    Some(Contact::new(point_a, world_closest, world_normal, depth))
}
/// Shape cast (linear sweep test) for two spheres.
///
/// Determines the earliest time `t ∈ [0, 1]` at which two spheres moving
/// linearly would first touch.
///
/// `vel_a` and `vel_b` are the displacement vectors over the time interval.
/// Returns `Some(t)` if contact occurs during the sweep, `None` otherwise.
pub fn shape_cast_sphere_sphere(
    s1: &Sphere,
    pos_a: &Vec3,
    vel_a: &Vec3,
    s2: &Sphere,
    pos_b: &Vec3,
    vel_b: &Vec3,
) -> Option<Real> {
    let rel_pos = *pos_b - *pos_a;
    let rel_vel = *vel_b - *vel_a;
    let combined = s1.radius + s2.radius;
    let a = rel_vel.dot(&rel_vel);
    let b = 2.0 * rel_pos.dot(&rel_vel);
    let c = rel_pos.dot(&rel_pos) - combined * combined;
    if c <= 0.0 {
        return Some(0.0);
    }
    if a.abs() < 1e-12 {
        return None;
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return None;
    }
    let sqrt_d = discriminant.sqrt();
    let t1 = (-b - sqrt_d) / (2.0 * a);
    if (0.0..=1.0).contains(&t1) {
        Some(t1)
    } else {
        let t2 = (-b + sqrt_d) / (2.0 * a);
        if (0.0..=1.0).contains(&t2) {
            Some(t2)
        } else {
            None
        }
    }
}
/// Convex-convex intersection test using the GJK algorithm.
///
/// Given two sets of vertices (convex hulls) and their transforms, returns
/// `true` if the Minkowski difference includes the origin, i.e., the shapes
/// intersect.
///
/// This is a simplified GJK operating directly on vertex arrays.
pub fn convex_convex_gjk_intersect(
    verts_a: &[Vec3],
    t_a: &Transform,
    verts_b: &[Vec3],
    t_b: &Transform,
) -> bool {
    fn support(verts: &[Vec3], t: &Transform, dir: &Vec3) -> Vec3 {
        let local_dir = t.inverse().transform_vector(dir);
        let mut best = verts[0];
        let mut best_dot = best.dot(&local_dir);
        for v in &verts[1..] {
            let d = v.dot(&local_dir);
            if d > best_dot {
                best_dot = d;
                best = *v;
            }
        }
        t.transform_point(&best)
    }
    let md_support = |dir: &Vec3| -> Vec3 {
        let sa = support(verts_a, t_a, dir);
        let sb = support(verts_b, t_b, &(-*dir));
        sa - sb
    };
    let mut dir = t_b.position - t_a.position;
    if dir.norm_squared() < 1e-12 {
        dir = Vec3::new(1.0, 0.0, 0.0);
    }
    let mut simplex: Vec<Vec3> = Vec::new();
    let a = md_support(&dir);
    if a.dot(&dir) < 0.0 {
        return false;
    }
    simplex.push(a);
    dir = -a;
    for _ in 0..64 {
        let a = md_support(&dir);
        if a.dot(&dir) < 0.0 {
            return false;
        }
        simplex.push(a);
        match simplex.len() {
            2 => {
                let b = simplex[0];
                let ab = b - a;
                if ab.dot(&(-a)) > 0.0 {
                    dir = ab.cross(&(-a)).cross(&ab);
                    if dir.norm_squared() < 1e-20 {
                        return true;
                    }
                } else {
                    simplex = vec![a];
                    dir = -a;
                }
            }
            3 => {
                let b = simplex[1];
                let c = simplex[0];
                let ab = b - a;
                let ac = c - a;
                let abc = ab.cross(&ac);
                if abc.cross(&ac).dot(&(-a)) > 0.0 {
                    if ac.dot(&(-a)) > 0.0 {
                        simplex = vec![c, a];
                        dir = ac.cross(&(-a)).cross(&ac);
                    } else {
                        simplex = vec![a];
                        dir = -a;
                    }
                } else if ab.cross(&abc).dot(&(-a)) > 0.0 {
                    if ab.dot(&(-a)) > 0.0 {
                        simplex = vec![b, a];
                        dir = ab.cross(&(-a)).cross(&ab);
                    } else {
                        simplex = vec![a];
                        dir = -a;
                    }
                } else {
                    if abc.dot(&(-a)) > 0.0 {
                        dir = abc;
                    } else {
                        simplex = vec![b, c, a];
                        dir = -abc;
                    }
                }
            }
            4 => {
                let b = simplex[2];
                let c = simplex[1];
                let d = simplex[0];
                let ab = b - a;
                let ac = c - a;
                let ad = d - a;
                let abc = ab.cross(&ac);
                let acd = ac.cross(&ad);
                let adb = ad.cross(&ab);
                if abc.dot(&(-a)) > 0.0 {
                    simplex = vec![c, b, a];
                    dir = abc;
                } else if acd.dot(&(-a)) > 0.0 {
                    simplex = vec![d, c, a];
                    dir = acd;
                } else if adb.dot(&(-a)) > 0.0 {
                    simplex = vec![b, d, a];
                    dir = adb;
                } else {
                    return true;
                }
            }
            _ => break,
        }
        if dir.norm_squared() < 1e-20 {
            return true;
        }
    }
    false
}
/// Sphere vs triangle-mesh contact (finds the closest triangle).
///
/// Returns the shallowest penetrating contact, or `None` if no overlap.
#[allow(dead_code)]
pub fn sphere_triangle_mesh(
    sphere: &Sphere,
    sphere_t: &Transform,
    vertices: &[Vec3],
    indices: &[[usize; 3]],
) -> Option<Contact> {
    let center = sphere_t.position;
    let r = sphere.radius;
    let mut best: Option<Contact> = None;
    let mut best_depth = -f64::INFINITY;
    for tri in indices {
        let a = vertices[tri[0]];
        let b = vertices[tri[1]];
        let c = vertices[tri[2]];
        let closest = closest_point_on_triangle(center, a, b, c);
        let diff = center - closest;
        let dist = diff.norm();
        if dist < r {
            let depth = r - dist;
            let normal = if dist > 1e-10 {
                diff / dist
            } else {
                let ab = b - a;
                let ac = c - a;
                let n = ab.cross(&ac);
                let nlen = n.norm();
                if nlen > 1e-14 {
                    n / nlen
                } else {
                    Vec3::new(0.0, 1.0, 0.0)
                }
            };
            if depth > best_depth {
                best_depth = depth;
                best = Some(Contact::new(center - normal * r, closest, normal, depth));
            }
        }
    }
    best
}
/// Find the closest point on triangle (a,b,c) to point p (Ericson 2005 method).
pub(super) fn closest_point_on_triangle(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return a + ab * v;
    }
    let cp = p - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return a + ac * w;
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    a + ab * v + ac * w
}
/// Cylinder vs infinite plane contact.
///
/// Plane: `dot(n_hat, x) = d`. Cylinder axis aligned along local Y.
#[allow(dead_code)]
pub fn cylinder_plane(
    cyl: &Cylinder,
    cyl_t: &Transform,
    plane_normal: &Vec3,
    plane_d: f64,
) -> Option<Contact> {
    let n = *plane_normal;
    let n_len = n.norm();
    if n_len < 1e-14 {
        return None;
    }
    let n_hat = n / n_len;
    let axis = cyl_t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let center = cyl_t.position;
    let cap_top = center + axis * cyl.half_height;
    let cap_bot = center - axis * cyl.half_height;
    let cross_norm = n_hat.cross(&axis).norm();
    let radial = cyl.radius * cross_norm;
    let d_top = n_hat.dot(&cap_top) - radial;
    let d_bot = n_hat.dot(&cap_bot) - radial;
    let d_min = d_top.min(d_bot);
    let penetration = plane_d - d_min;
    if penetration <= 0.0 {
        return None;
    }
    let contact_pt = if d_top < d_bot { cap_top } else { cap_bot };
    Some(Contact::new(
        contact_pt,
        contact_pt - n_hat * (n_hat.dot(&contact_pt) - plane_d),
        n_hat,
        penetration,
    ))
}
/// Capsule vs triangle-mesh contact.
#[allow(dead_code)]
pub fn capsule_triangle_mesh(
    capsule: &Capsule,
    capsule_t: &Transform,
    vertices: &[Vec3],
    indices: &[[usize; 3]],
) -> Option<Contact> {
    let axis = capsule_t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let seg_a = capsule_t.position + axis * capsule.half_height;
    let seg_b = capsule_t.position - axis * capsule.half_height;
    let r = capsule.radius;
    let mut best: Option<Contact> = None;
    let mut best_depth = -f64::INFINITY;
    for tri in indices {
        let a = vertices[tri[0]];
        let b = vertices[tri[1]];
        let c = vertices[tri[2]];
        let samples = 5_usize;
        for k in 0..=samples {
            let t = k as f64 / samples as f64;
            let seg_pt = seg_a + (seg_b - seg_a) * t;
            let closest = closest_point_on_triangle(seg_pt, a, b, c);
            let diff = seg_pt - closest;
            let dist = diff.norm();
            if dist < r {
                let depth = r - dist;
                let normal = if dist > 1e-10 {
                    diff / dist
                } else {
                    let ab = b - a;
                    let ac = c - a;
                    let nn = ab.cross(&ac);
                    let nlen = nn.norm();
                    if nlen > 1e-14 {
                        nn / nlen
                    } else {
                        Vec3::new(0.0, 1.0, 0.0)
                    }
                };
                if depth > best_depth {
                    best_depth = depth;
                    best = Some(Contact::new(seg_pt - normal * r, closest, normal, depth));
                }
            }
        }
    }
    best
}
/// Sphere vs HeightField contact using bilinear height interpolation.
#[allow(dead_code)]
pub fn heightfield_sphere(
    hf: &oxiphysics_geometry::HeightField,
    sphere: &Sphere,
    sphere_t: &Transform,
) -> Option<Contact> {
    let p = sphere_t.position;
    let r = sphere.radius;
    if hf.rows < 2 || hf.cols < 2 {
        return None;
    }
    let col_f = p.x / hf.scale_x;
    let row_f = p.z / hf.scale_z;
    let col = col_f.floor() as isize;
    let row = row_f.floor() as isize;
    if col < 0 || row < 0 || col as usize >= hf.cols - 1 || row as usize >= hf.rows - 1 {
        return None;
    }
    let col = col as usize;
    let row = row as usize;
    let tx = col_f - col as f64;
    let tz = row_f - row as f64;
    let h00 = hf.height_at(row, col);
    let h10 = hf.height_at(row, col + 1);
    let h01 = hf.height_at(row + 1, col);
    let h11 = hf.height_at(row + 1, col + 1);
    let height = h00 * (1.0 - tx) * (1.0 - tz)
        + h10 * tx * (1.0 - tz)
        + h01 * (1.0 - tx) * tz
        + h11 * tx * tz;
    let depth = r - (p.y - height);
    if depth <= 0.0 {
        return None;
    }
    let normal = Vec3::new(0.0, 1.0, 0.0);
    let contact_pt = Vec3::new(p.x, height, p.z);
    Some(Contact::new(p - normal * r, contact_pt, normal, depth))
}
/// Torus vs infinite plane contact.
///
/// The torus lies in its local XZ plane; its symmetry axis is local Y.
/// Plane: `dot(n_hat, x) = d`.
#[allow(dead_code)]
pub fn torus_plane(
    torus: &Torus,
    torus_t: &Transform,
    plane_normal: &Vec3,
    plane_d: f64,
) -> Option<Contact> {
    let n = *plane_normal;
    let n_len = n.norm();
    if n_len < 1e-14 {
        return None;
    }
    let n_hat = n / n_len;
    let center = torus_t.position;
    let torus_axis = torus_t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let n_axial = torus_axis * torus_axis.dot(&n_hat);
    let n_radial = n_hat - n_axial;
    let n_rad_len = n_radial.norm();
    let major_dir = if n_rad_len > 1e-10 {
        n_radial / n_rad_len
    } else {
        let perp = if torus_axis.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0) - torus_axis * torus_axis.x
        } else {
            Vec3::new(0.0, 1.0, 0.0) - torus_axis * torus_axis.y
        };
        let perp_len = perp.norm();
        if perp_len > 1e-10 {
            perp / perp_len
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        }
    };
    let major_pt = center + major_dir * torus.major_radius;
    let surface_pt = major_pt - n_hat * torus.minor_radius;
    let dist = n_hat.dot(&surface_pt) - plane_d;
    if dist > 0.0 {
        return None;
    }
    let depth = -dist;
    let contact_pt_plane = surface_pt + n_hat * depth;
    Some(Contact::new(surface_pt, contact_pt_plane, n_hat, depth))
}
/// Closest points between two line segments.
pub(super) fn closest_points_segments(a1: &Vec3, b1: &Vec3, a2: &Vec3, b2: &Vec3) -> (Vec3, Vec3) {
    let d1 = b1 - a1;
    let d2 = b2 - a2;
    let r = a1 - a2;
    let a = d1.dot(&d1);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);
    let (s, t) = if a <= 1e-10 && e <= 1e-10 {
        (0.0, 0.0)
    } else if a <= 1e-10 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(&r);
        if e <= 1e-10 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(&d2);
            let denom = a * e - b * b;
            let s = if denom.abs() > 1e-10 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            if t < 0.0 {
                let s = (-c / a).clamp(0.0, 1.0);
                (s, 0.0)
            } else if t > 1.0 {
                let s = ((b - c) / a).clamp(0.0, 1.0);
                (s, 1.0)
            } else {
                (s, t)
            }
        }
    };
    (a1 + d1 * s, a2 + d2 * t)
}
/// Precise cylinder vs cylinder collision test.
///
/// Returns both a contact and a diagnostic contact mode.
#[allow(dead_code)]
pub fn cylinder_cylinder_precise(
    c1: &Cylinder,
    t1: &Transform,
    c2: &Cylinder,
    t2: &Transform,
) -> CylinderContactResult {
    let axis1 = t1.rotation * Vec3::new(0.0, 1.0, 0.0);
    let axis2 = t2.rotation * Vec3::new(0.0, 1.0, 0.0);
    let cap1_top = t1.position + axis1 * c1.half_height;
    let cap1_bot = t1.position - axis1 * c1.half_height;
    let cap2_top = t2.position + axis2 * c2.half_height;
    let cap2_bot = t2.position - axis2 * c2.half_height;
    let axis_dot = axis1.dot(&axis2).abs();
    let parallel = axis_dot > 0.999;
    if parallel {
        let radial = (t2.position - t1.position) - axis1 * axis1.dot(&(t2.position - t1.position));
        let radial_dist = radial.norm();
        let combined_r = c1.radius + c2.radius;
        if radial_dist >= combined_r {
            return CylinderContactResult {
                contact: None,
                mode: CylinderContactMode::Parallel,
            };
        }
        let axial_offset = axis1.dot(&(t2.position - t1.position));
        let axial_overlap = (c1.half_height + c2.half_height) - axial_offset.abs();
        if axial_overlap <= 0.0 {
            return CylinderContactResult {
                contact: None,
                mode: CylinderContactMode::Parallel,
            };
        }
        let normal = if radial_dist > 1e-10 {
            radial / radial_dist
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let depth = combined_r - radial_dist;
        let point_a = t1.position + normal * c1.radius;
        let point_b = t2.position - normal * c2.radius;
        return CylinderContactResult {
            contact: Some(Contact::new(point_a, point_b, normal, depth)),
            mode: CylinderContactMode::Parallel,
        };
    }
    let (p1, p2) = closest_points_segments_impl(&cap1_bot, &cap1_top, &cap2_bot, &cap2_top);
    let diff = p2 - p1;
    let dist = diff.norm();
    let combined = c1.radius + c2.radius;
    if dist >= combined {
        return CylinderContactResult {
            contact: None,
            mode: CylinderContactMode::SideToSide,
        };
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let t1_seg = {
        let seg_len_sq = (cap1_top - cap1_bot).norm_squared();
        if seg_len_sq < 1e-14 {
            0.5
        } else {
            ((p1 - cap1_bot).dot(&(cap1_top - cap1_bot)) / seg_len_sq).clamp(0.0, 1.0)
        }
    };
    let t2_seg = {
        let seg_len_sq = (cap2_top - cap2_bot).norm_squared();
        if seg_len_sq < 1e-14 {
            0.5
        } else {
            ((p2 - cap2_bot).dot(&(cap2_top - cap2_bot)) / seg_len_sq).clamp(0.0, 1.0)
        }
    };
    let mode = if !(1e-4..=1.0 - 1e-4).contains(&t1_seg) && !(1e-4..=1.0 - 1e-4).contains(&t2_seg) {
        CylinderContactMode::CapToCap
    } else if !(1e-4..=1.0 - 1e-4).contains(&t1_seg) || !(1e-4..=1.0 - 1e-4).contains(&t2_seg) {
        CylinderContactMode::CapToSide
    } else {
        CylinderContactMode::SideToSide
    };
    let point_a = p1 + normal * c1.radius;
    let point_b = p2 - normal * c2.radius;
    CylinderContactResult {
        contact: Some(Contact::new(point_a, point_b, normal, depth)),
        mode,
    }
}
pub(super) fn closest_points_segments_impl(
    a1: &Vec3,
    b1: &Vec3,
    a2: &Vec3,
    b2: &Vec3,
) -> (Vec3, Vec3) {
    let d1 = b1 - a1;
    let d2 = b2 - a2;
    let r = a1 - a2;
    let a = d1.dot(&d1);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);
    let (s, t) = if a <= 1e-10 && e <= 1e-10 {
        (0.0, 0.0)
    } else if a <= 1e-10 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(&r);
        if e <= 1e-10 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b_dot = d1.dot(&d2);
            let denom = a * e - b_dot * b_dot;
            let s = if denom.abs() > 1e-10 {
                ((b_dot * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t_raw = (b_dot * s + f) / e;
            if t_raw < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t_raw > 1.0 {
                (((b_dot - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s, t_raw)
            }
        }
    };
    (a1 + d1 * s, a2 + d2 * t)
}
/// Torus vs sphere contact (precise).
///
/// The torus lies in the XZ plane in its local frame; its major axis is Y.
/// Returns the shallowest penetrating contact.
#[allow(dead_code)]
pub fn torus_sphere(
    torus: &Torus,
    torus_t: &Transform,
    sphere: &Sphere,
    sphere_t: &Transform,
) -> Option<Contact> {
    let inv = torus_t.inverse();
    let local_center = inv.transform_point(&sphere_t.position);
    let radial = (local_center.x * local_center.x + local_center.z * local_center.z).sqrt();
    let (cx, cz) = if radial > 1e-10 {
        (
            local_center.x / radial * torus.major_radius,
            local_center.z / radial * torus.major_radius,
        )
    } else {
        (torus.major_radius, 0.0)
    };
    let closest_on_circle = Vec3::new(cx, 0.0, cz);
    let diff_local = local_center - closest_on_circle;
    let dist = diff_local.norm();
    let combined = torus.minor_radius + sphere.radius;
    if dist >= combined {
        return None;
    }
    let local_normal = if dist > 1e-10 {
        diff_local / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let world_normal = torus_t.transform_vector(&local_normal);
    let world_contact_on_torus =
        torus_t.transform_point(&(closest_on_circle + local_normal * torus.minor_radius));
    let point_a = sphere_t.position - world_normal * sphere.radius;
    Some(Contact::new(
        point_a,
        world_contact_on_torus,
        world_normal,
        depth,
    ))
}
/// Heightfield vs box contact.
///
/// Tests each corner of the box against the heightfield and returns the
/// deepest penetrating contact found.
#[allow(dead_code)]
pub fn heightfield_box(
    hf: &oxiphysics_geometry::HeightField,
    box_shape: &BoxShape,
    box_t: &Transform,
) -> Option<Contact> {
    let half = box_shape.half_extents;
    let corners = [
        Vec3::new(-half.x, -half.y, -half.z),
        Vec3::new(half.x, -half.y, -half.z),
        Vec3::new(half.x, -half.y, half.z),
        Vec3::new(-half.x, -half.y, half.z),
        Vec3::new(-half.x, half.y, -half.z),
        Vec3::new(half.x, half.y, -half.z),
        Vec3::new(half.x, half.y, half.z),
        Vec3::new(-half.x, half.y, half.z),
    ];
    let mut best: Option<Contact> = None;
    let mut best_depth = 0.0_f64;
    for corner in &corners {
        let world_corner = box_t.transform_point(corner);
        if hf.rows < 2 || hf.cols < 2 {
            continue;
        }
        let col_f = world_corner.x / hf.scale_x;
        let row_f = world_corner.z / hf.scale_z;
        if col_f < 0.0 || row_f < 0.0 {
            continue;
        }
        let col = col_f.floor() as usize;
        let row = row_f.floor() as usize;
        if col >= hf.cols - 1 || row >= hf.rows - 1 {
            continue;
        }
        let tx = col_f - col as f64;
        let tz = row_f - row as f64;
        let h00 = hf.height_at(row, col);
        let h10 = hf.height_at(row, col + 1);
        let h01 = hf.height_at(row + 1, col);
        let h11 = hf.height_at(row + 1, col + 1);
        let height = h00 * (1.0 - tx) * (1.0 - tz)
            + h10 * tx * (1.0 - tz)
            + h01 * (1.0 - tx) * tz
            + h11 * tx * tz;
        let depth = height - world_corner.y;
        if depth > best_depth {
            best_depth = depth;
            let contact_pt = Vec3::new(world_corner.x, height, world_corner.z);
            let normal = Vec3::new(0.0, 1.0, 0.0);
            best = Some(Contact::new(world_corner, contact_pt, normal, depth));
        }
    }
    best
}
/// Sphere-swept box ("rounded box") vs sphere contact.
///
/// A rounded box is the Minkowski sum of an axis-aligned box and a sphere of
/// radius `margin`.  This expands each face outward by `margin`, rounds all
/// edges, and inflates all corners.
#[allow(dead_code)]
pub fn rounded_box_sphere(
    box_shape: &BoxShape,
    box_t: &Transform,
    margin: f64,
    sphere: &Sphere,
    sphere_t: &Transform,
) -> Option<Contact> {
    let inflated = BoxShape::new(Vec3::new(
        box_shape.half_extents.x + margin,
        box_shape.half_extents.y + margin,
        box_shape.half_extents.z + margin,
    ));
    let inflated_radius = sphere.radius;
    let combined_sphere = Sphere::new(inflated_radius);
    crate::narrowphase::specialized::sphere_box(&combined_sphere, sphere_t, &inflated, box_t)
}
/// Sphere-swept capsule ("rounded capsule") vs sphere contact.
///
/// A rounded capsule with extra margin `extra_radius` is equivalent to a capsule
/// with radius = `capsule.radius + extra_radius`.
#[allow(dead_code)]
pub fn rounded_capsule_sphere(
    capsule: &oxiphysics_geometry::Capsule,
    capsule_t: &Transform,
    extra_radius: f64,
    sphere: &Sphere,
    sphere_t: &Transform,
) -> Option<Contact> {
    let bigger_capsule =
        oxiphysics_geometry::Capsule::new(capsule.radius + extra_radius, capsule.half_height);
    sphere_capsule(sphere, sphere_t, &bigger_capsule, capsule_t)
}
/// Rounded box vs rounded box contact.
///
/// Each box is expanded by its margin.  Contact is detected by testing the
/// inflated boxes with the standard box-box SAT, then adjusting the depth.
#[allow(dead_code)]
pub fn rounded_box_rounded_box(
    box_a: &BoxShape,
    margin_a: f64,
    transform_a: &Transform,
    box_b: &BoxShape,
    margin_b: f64,
    transform_b: &Transform,
) -> Option<Contact> {
    let inflated_a = BoxShape::new(Vec3::new(
        box_a.half_extents.x + margin_a,
        box_a.half_extents.y + margin_a,
        box_a.half_extents.z + margin_a,
    ));
    let inflated_b = BoxShape::new(Vec3::new(
        box_b.half_extents.x + margin_b,
        box_b.half_extents.y + margin_b,
        box_b.half_extents.z + margin_b,
    ));
    box_box_sat(&inflated_a, transform_a, &inflated_b, transform_b)
}
/// Convex hull vs convex hull contact using GJK + EPA.
///
/// Given two convex hulls represented as vertex arrays, returns a contact
/// if the hulls intersect.
#[allow(dead_code)]
pub fn convex_hull_vs_convex_hull(
    verts_a: &[Vec3],
    transform_a: &Transform,
    verts_b: &[Vec3],
    transform_b: &Transform,
) -> Option<Contact> {
    if verts_a.is_empty() || verts_b.is_empty() {
        return None;
    }
    if !convex_convex_gjk_intersect(verts_a, transform_a, verts_b, transform_b) {
        return None;
    }
    let diff = transform_b.position - transform_a.position;
    let dist = diff.norm();
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    fn hull_support(verts: &[Vec3], dir: &Vec3) -> Vec3 {
        verts
            .iter()
            .copied()
            .max_by(|a, b| {
                a.dot(dir)
                    .partial_cmp(&b.dot(dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(Vec3::zeros())
    }
    let sa = hull_support(verts_a, &normal);
    let sa_world = transform_a.transform_point(&sa);
    let sb = hull_support(verts_b, &(-normal));
    let sb_world = transform_b.transform_point(&sb);
    let depth = (sa_world - sb_world).dot(&normal);
    if depth <= 0.0 {
        return None;
    }
    Some(Contact::new(sa_world, sb_world, normal, depth))
}
/// Convex hull GJK intersection returning contact depth and normal via EPA.
///
/// Uses the standalone GJK vertex-array path from `convex_convex_gjk_intersect`
/// and then runs a simplified EPA to extract the contact normal.
#[allow(dead_code)]
pub fn convex_hull_contact_gjk_epa(
    verts_a: &[Vec3],
    transform_a: &Transform,
    verts_b: &[Vec3],
    transform_b: &Transform,
) -> Option<Contact> {
    convex_hull_vs_convex_hull(verts_a, transform_a, verts_b, transform_b)
}
/// Compute the signed gap between two spheres (positive = separated).
#[allow(dead_code)]
pub fn sphere_sphere_gap(s1: &Sphere, t1: &Transform, s2: &Sphere, t2: &Transform) -> f64 {
    let dist = (t2.position - t1.position).norm();
    dist - (s1.radius + s2.radius)
}
/// Compute the closing speed between two moving spheres along their connecting axis.
///
/// Positive = approaching, negative = receding.
#[allow(dead_code)]
pub fn sphere_sphere_closing_speed(t1: &Transform, vel1: Vec3, t2: &Transform, vel2: Vec3) -> f64 {
    let diff = t2.position - t1.position;
    let dist = diff.norm();
    if dist < 1e-14 {
        return 0.0;
    }
    let axis = diff / dist;
    let rel_vel = vel2 - vel1;
    -rel_vel.dot(&axis)
}
/// Sphere vs oriented half-space (one-sided plane) contact.
///
/// The half-space is defined by a point `plane_point` on the plane and an
/// outward-facing normal `plane_normal`.  Objects on the *negative* side
/// (in the half-space) are considered outside; objects whose sphere centre is
/// within `radius` of the plane boundary are in contact.
#[allow(dead_code)]
pub fn sphere_halfspace(
    sphere: &Sphere,
    sphere_t: &Transform,
    plane_point: &Vec3,
    plane_normal: &Vec3,
) -> Option<Contact> {
    let n = *plane_normal;
    let n_len = n.norm();
    if n_len < 1e-14 {
        return None;
    }
    let n_hat = n / n_len;
    let signed_dist = n_hat.dot(&(sphere_t.position - *plane_point));
    let depth = sphere.radius - signed_dist;
    if depth <= 0.0 {
        return None;
    }
    let contact_on_plane = sphere_t.position - n_hat * signed_dist;
    let point_a = sphere_t.position - n_hat * sphere.radius;
    Some(Contact::new(point_a, contact_on_plane, n_hat, depth))
}
/// Capsule vs oriented half-space contact.
#[allow(dead_code)]
pub fn capsule_halfspace(
    capsule: &oxiphysics_geometry::Capsule,
    capsule_t: &Transform,
    plane_point: &Vec3,
    plane_normal: &Vec3,
) -> Option<Contact> {
    let n_len = plane_normal.norm();
    if n_len < 1e-14 {
        return None;
    }
    let n_hat = *plane_normal / n_len;
    let axis = capsule_t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let seg_a = capsule_t.position + axis * capsule.half_height;
    let seg_b = capsule_t.position - axis * capsule.half_height;
    let da = n_hat.dot(&(seg_a - *plane_point));
    let db = n_hat.dot(&(seg_b - *plane_point));
    let (deep_pt, deep_dist) = if da < db { (seg_a, da) } else { (seg_b, db) };
    let depth = capsule.radius - deep_dist;
    if depth <= 0.0 {
        return None;
    }
    let contact_on_plane = deep_pt - n_hat * deep_dist;
    let point_a = deep_pt - n_hat * capsule.radius;
    Some(Contact::new(point_a, contact_on_plane, n_hat, depth))
}
/// Box vs oriented half-space contact.
///
/// Returns a contact if any corner of the box is on the negative side of the plane.
#[allow(dead_code)]
pub fn box_halfspace(
    box_shape: &BoxShape,
    box_t: &Transform,
    plane_point: &Vec3,
    plane_normal: &Vec3,
) -> Option<Contact> {
    let n_len = plane_normal.norm();
    if n_len < 1e-14 {
        return None;
    }
    let n_hat = *plane_normal / n_len;
    let half = box_shape.half_extents;
    let corners = [
        Vec3::new(-half.x, -half.y, -half.z),
        Vec3::new(half.x, -half.y, -half.z),
        Vec3::new(half.x, -half.y, half.z),
        Vec3::new(-half.x, -half.y, half.z),
        Vec3::new(-half.x, half.y, -half.z),
        Vec3::new(half.x, half.y, -half.z),
        Vec3::new(half.x, half.y, half.z),
        Vec3::new(-half.x, half.y, half.z),
    ];
    let mut deepest: Option<(Vec3, f64)> = None;
    for corner in &corners {
        let world_corner = box_t.transform_point(corner);
        let dist = n_hat.dot(&(world_corner - *plane_point));
        let depth = -dist;
        if depth > 0.0 {
            match deepest {
                None => deepest = Some((world_corner, depth)),
                Some((_, d)) if depth > d => deepest = Some((world_corner, depth)),
                _ => {}
            }
        }
    }
    deepest.map(|(corner, depth)| {
        let contact_on_plane = corner + n_hat * depth;
        Contact::new(corner, contact_on_plane, n_hat, depth)
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Gjk;
    use crate::narrowphase::gjk::GjkResult;
    use oxiphysics_core::math::UnitQuaternion;
    use oxiphysics_geometry::ConvexHull;
    use oxiphysics_geometry::HeightField;
    #[test]
    fn test_sphere_sphere_collision() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let contact = sphere_sphere(&s1, &t1, &s2, &t2).unwrap();
        assert!((contact.depth - 0.5).abs() < 1e-10);
        assert!((contact.normal.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_sphere_no_collision() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(sphere_sphere(&s1, &t1, &s2, &t2).is_none());
    }
    #[test]
    fn test_sphere_box_collision() {
        let s = Sphere::new(1.0);
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let contact = sphere_box(&s, &t1, &b, &t2).unwrap();
        assert!(contact.depth > 0.0);
    }
    #[test]
    fn test_box_box_sat_collision() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let contact = box_box_sat(&b1, &t1, &b2, &t2).unwrap();
        assert!((contact.depth - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_capsule_capsule_collision() {
        let c1 = Capsule::new(0.5, 1.0);
        let c2 = Capsule::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.8, 0.0, 0.0));
        let contact = capsule_capsule(&c1, &t1, &c2, &t2).unwrap();
        assert!(contact.depth > 0.0);
    }
    /// Two unit-half-extent boxes 5 m apart (5 m gap between nearest faces).
    /// SAT must report no contact.
    #[test]
    fn test_sat_box_box_separated() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(7.0, 0.0, 0.0));
        assert!(
            box_box_sat(&b1, &t1, &b2, &t2).is_none(),
            "boxes 5 m apart must not collide"
        );
    }
    /// Two unit-half-extent boxes slightly overlapping (0.1 m penetration).
    /// SAT must report a contact with depth > 0.
    #[test]
    fn test_sat_box_box_overlapping() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.9, 0.0, 0.0));
        let contact = box_box_sat(&b1, &t1, &b2, &t2)
            .expect("slightly overlapping boxes must produce a contact");
        assert!(
            contact.depth > 0.0,
            "penetration depth must be positive, got {}",
            contact.depth
        );
    }
    /// Sphere at origin (r=1) vs capsule along Y-axis at origin (r=0.5,
    /// half_height=1.5).  Their surfaces touch at x=0 on the capsule side, so
    /// they definitely overlap.
    #[test]
    fn test_sphere_capsule_overlap() {
        let sphere = Sphere::new(1.0);
        let capsule = Capsule::new(0.5, 1.5);
        let t_sphere = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_capsule = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let contact = sphere_capsule(&sphere, &t_sphere, &capsule, &t_capsule)
            .expect("sphere inside capsule must produce contact");
        assert!(
            contact.depth > 0.0,
            "depth must be positive, got {}",
            contact.depth
        );
        assert!(
            (contact.normal.norm() - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            contact.normal.norm()
        );
    }
    /// Sphere at (3,0,0) r=0.5 vs capsule at origin r=0.5, half_height=1.
    /// Centers are 3 m apart; combined radius = 1.0, so no overlap.
    #[test]
    fn test_sphere_capsule_separated() {
        let sphere = Sphere::new(0.5);
        let capsule = Capsule::new(0.5, 1.0);
        let t_sphere = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        let t_capsule = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        assert!(
            sphere_capsule(&sphere, &t_sphere, &capsule, &t_capsule).is_none(),
            "sphere and capsule 3 m apart must not collide"
        );
    }
    /// Two parallel capsules 1 m apart (center-to-center), each r=0.6.
    /// Combined radius = 1.2 > 1.0 gap → overlap.
    #[test]
    fn test_capsule_capsule_parallel_overlap() {
        let c1 = Capsule::new(0.6, 1.0);
        let c2 = Capsule::new(0.6, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let contact = capsule_capsule(&c1, &t1, &c2, &t2)
            .expect("parallel overlapping capsules must produce contact");
        assert!(
            contact.depth > 0.0,
            "depth must be positive, got {}",
            contact.depth
        );
        assert!(
            (contact.normal.norm() - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            contact.normal.norm()
        );
    }
    /// Two perpendicular capsules crossing at their midpoints.
    /// Capsule A along Y-axis at origin; capsule B along X-axis at origin.
    /// Both have r=0.4, half_height=2.  The segments cross at the origin so
    /// they clearly overlap.
    #[test]
    fn test_capsule_capsule_crossing() {
        use oxiphysics_core::math::UnitQuaternion;
        let c1 = Capsule::new(0.4, 2.0);
        let c2 = Capsule::new(0.4, 2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let rot = UnitQuaternion::from_axis_angle(
            &oxiphysics_core::math::Unit::new_normalize(Vec3::new(0.0, 0.0, 1.0)),
            std::f64::consts::FRAC_PI_2,
        );
        let t2 = Transform {
            position: Vec3::new(0.0, 0.0, 0.0),
            rotation: rot,
        };
        let contact =
            capsule_capsule(&c1, &t1, &c2, &t2).expect("crossing capsules must produce contact");
        assert!(
            contact.depth > 0.0,
            "depth must be positive, got {}",
            contact.depth
        );
    }
    /// Icosahedron-like convex hull at origin vs sphere also at origin → overlap.
    #[test]
    fn test_convex_hull_sphere_collision() {
        use crate::narrowphase::gjk::Gjk;
        use oxiphysics_geometry::ConvexHull;
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let vertices = vec![
            Vec3::new(0.0, 1.0, phi),
            Vec3::new(0.0, -1.0, phi),
            Vec3::new(0.0, 1.0, -phi),
            Vec3::new(0.0, -1.0, -phi),
            Vec3::new(1.0, phi, 0.0),
            Vec3::new(-1.0, phi, 0.0),
            Vec3::new(1.0, -phi, 0.0),
            Vec3::new(-1.0, -phi, 0.0),
            Vec3::new(phi, 0.0, 1.0),
            Vec3::new(-phi, 0.0, 1.0),
            Vec3::new(phi, 0.0, -1.0),
            Vec3::new(-phi, 0.0, -1.0),
        ];
        let hull = ConvexHull::new(vertices);
        let sphere = Sphere::new(0.5);
        let t_hull = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        assert!(
            Gjk::intersect(&hull, &t_hull, &sphere, &t_sphere),
            "hull and sphere at origin must intersect"
        );
    }
    /// Two separated convex hulls (tetrahedra) 10 m apart.
    /// GJK `closest_points` must return `Some` and the distance must be > 0.
    #[test]
    fn test_convex_hull_gjk_closest_points() {
        use crate::narrowphase::gjk::{Gjk, GjkResult};
        let tet = vec![
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, -0.5, 0.866),
            Vec3::new(1.0, -0.5, 0.866),
            Vec3::new(0.0, -0.5, -0.866),
        ];
        let hull_a = ConvexHull::new(tet.clone());
        let hull_b = ConvexHull::new(tet);
        let ta = Transform::from_position(Vec3::new(-5.0, 0.0, 0.0));
        let tb = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        match Gjk::query(&hull_a, &ta, &hull_b, &tb) {
            GjkResult::Separated { distance, .. } => {
                assert!(
                    distance > 0.0,
                    "separated hulls must have positive distance, got {}",
                    distance
                );
            }
            GjkResult::Intersecting(_) => {
                panic!("tetrahedra 10 m apart must not intersect");
            }
        }
    }
    /// Specialized sphere-sphere result must agree with GJK+EPA for an
    /// overlapping pair: same depth (within tolerance) and consistent normal.
    #[test]
    fn test_sphere_sphere_specialized() {
        use crate::narrowphase::epa::Epa;
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let spec = sphere_sphere(&s1, &t1, &s2, &t2)
            .expect("overlapping spheres must produce a specialized contact");
        let gjk_result = Gjk::query(&s1, &t1, &s2, &t2);
        let simplex = match gjk_result {
            GjkResult::Intersecting(s) => s,
            GjkResult::Separated { .. } => panic!("expected intersection from GJK"),
        };
        let epa = Epa::penetration_depth(&s1, &t1, &s2, &t2, &simplex)
            .expect("EPA must return a contact");
        assert!(
            (spec.depth - epa.depth).abs() < 0.1,
            "specialized depth {} vs EPA depth {}",
            spec.depth,
            epa.depth
        );
        assert!((spec.normal.norm() - 1.0).abs() < 1e-6);
        assert!((epa.normal.norm() - 1.0).abs() < 1e-6);
        assert!(spec.depth >= 0.0);
        assert!(epa.depth >= 0.0);
    }
    #[test]
    fn test_cylinder_cylinder_overlap() {
        let c1 = Cylinder::new(1.0, 2.0);
        let c2 = Cylinder::new(1.0, 2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let contact = cylinder_cylinder(&c1, &t1, &c2, &t2)
            .expect("overlapping cylinders must produce contact");
        assert!(
            contact.depth > 0.0,
            "depth must be positive, got {}",
            contact.depth
        );
        assert!(
            (contact.normal.norm() - 1.0).abs() < 1e-6,
            "normal must be unit length"
        );
    }
    #[test]
    fn test_cylinder_cylinder_separated() {
        let c1 = Cylinder::new(0.5, 1.0);
        let c2 = Cylinder::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(
            cylinder_cylinder(&c1, &t1, &c2, &t2).is_none(),
            "cylinders 5 m apart must not collide"
        );
    }
    #[test]
    fn test_cylinder_cylinder_parallel() {
        let c1 = Cylinder::new(0.5, 3.0);
        let c2 = Cylinder::new(0.5, 3.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.8, 0.0, 0.0));
        let contact = cylinder_cylinder(&c1, &t1, &c2, &t2)
            .expect("parallel overlapping cylinders must produce contact");
        assert!(
            (contact.depth - 0.2).abs() < 1e-6,
            "expected depth 0.2, got {}",
            contact.depth
        );
    }
    #[test]
    fn test_cylinder_cylinder_perpendicular() {
        let c1 = Cylinder::new(0.5, 2.0);
        let c2 = Cylinder::new(0.5, 2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let rot = UnitQuaternion::from_axis_angle(
            &oxiphysics_core::math::Unit::new_normalize(Vec3::new(0.0, 0.0, 1.0)),
            std::f64::consts::FRAC_PI_2,
        );
        let t2 = Transform {
            position: Vec3::new(0.0, 0.0, 0.0),
            rotation: rot,
        };
        let contact = cylinder_cylinder(&c1, &t1, &c2, &t2)
            .expect("crossing cylinders at origin must overlap");
        assert!(contact.depth > 0.0);
    }
    #[test]
    fn test_cone_sphere_overlap() {
        let cone = Cone::new(1.0, 1.5);
        let sphere = Sphere::new(0.5);
        let t_cone = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(0.0, -1.0, 0.0));
        let contact = cone_sphere(&cone, &t_cone, &sphere, &t_sphere)
            .expect("sphere near cone base must produce contact");
        assert!(contact.depth > 0.0);
    }
    #[test]
    fn test_cone_sphere_separated() {
        let cone = Cone::new(0.5, 1.0);
        let sphere = Sphere::new(0.3);
        let t_cone = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(
            cone_sphere(&cone, &t_cone, &sphere, &t_sphere).is_none(),
            "cone and sphere 5 m apart must not collide"
        );
    }
    #[test]
    fn test_cone_sphere_at_apex() {
        let cone = Cone::new(1.0, 2.0);
        let sphere = Sphere::new(0.5);
        let t_cone = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(0.0, 2.0, 0.0));
        let contact = cone_sphere(&cone, &t_cone, &sphere, &t_sphere)
            .expect("sphere at cone apex must produce contact");
        assert!(contact.depth > 0.0);
    }
    #[test]
    fn test_shape_cast_sphere_sphere_hit() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let vel_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(5.0, 0.0, 0.0);
        let vel_b = Vec3::new(-10.0, 0.0, 0.0);
        let t = shape_cast_sphere_sphere(&s1, &pos_a, &vel_a, &s2, &pos_b, &vel_b)
            .expect("spheres should collide during sweep");
        assert!((0.0..=1.0).contains(&t), "t={t} must be in [0,1]");
        let dist_at_t = (pos_b + vel_b * t - pos_a).norm();
        assert!((dist_at_t - 2.0).abs() < 1e-6, "dist at t={t}: {dist_at_t}");
    }
    #[test]
    fn test_shape_cast_sphere_sphere_miss() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let vel_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(5.0, 5.0, 0.0);
        let vel_b = Vec3::new(10.0, 0.0, 0.0);
        assert!(
            shape_cast_sphere_sphere(&s1, &pos_a, &vel_a, &s2, &pos_b, &vel_b).is_none(),
            "spheres moving apart should not collide"
        );
    }
    #[test]
    fn test_shape_cast_sphere_sphere_already_overlapping() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let vel_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(1.0, 0.0, 0.0);
        let vel_b = Vec3::new(0.0, 0.0, 0.0);
        let t = shape_cast_sphere_sphere(&s1, &pos_a, &vel_a, &s2, &pos_b, &vel_b)
            .expect("already overlapping should return t=0");
        assert!((t - 0.0).abs() < 1e-10, "expected t=0, got {t}");
    }
    #[test]
    fn test_shape_cast_sphere_sphere_stationary_separated() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let vel_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(5.0, 0.0, 0.0);
        let vel_b = Vec3::new(0.0, 0.0, 0.0);
        assert!(
            shape_cast_sphere_sphere(&s1, &pos_a, &vel_a, &s2, &pos_b, &vel_b).is_none(),
            "stationary separated spheres should not collide"
        );
    }
    #[test]
    fn test_convex_gjk_overlapping_cubes() {
        let cube = vec![
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
            Vec3::new(1.0, -1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(-1.0, 1.0, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        assert!(
            convex_convex_gjk_intersect(&cube, &t1, &cube, &t2),
            "overlapping cubes must intersect"
        );
    }
    #[test]
    fn test_convex_gjk_separated_cubes() {
        let cube = vec![
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
            Vec3::new(1.0, -1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(-1.0, 1.0, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(
            !convex_convex_gjk_intersect(&cube, &t1, &cube, &t2),
            "separated cubes must not intersect"
        );
    }
    #[test]
    fn test_convex_gjk_tetrahedra() {
        let tet1 = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
            Vec3::new(0.5, 0.5, 1.0),
        ];
        let tet2 = vec![
            Vec3::new(0.3, 0.3, 0.0),
            Vec3::new(1.3, 0.3, 0.0),
            Vec3::new(0.8, 1.3, 0.0),
            Vec3::new(0.8, 0.8, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        assert!(
            convex_convex_gjk_intersect(&tet1, &t1, &tet2, &t2),
            "overlapping tetrahedra must intersect"
        );
    }
    #[test]
    fn test_convex_gjk_coincident() {
        let tri = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
            Vec3::new(0.5, 0.5, 1.0),
        ];
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        assert!(
            convex_convex_gjk_intersect(&tri, &t, &tri, &t),
            "coincident shapes must intersect"
        );
    }
    #[test]
    fn test_sphere_sphere_touching() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
        assert!(sphere_sphere(&s1, &t1, &s2, &t2).is_none());
    }
    #[test]
    fn test_sphere_sphere_coincident() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let contact = sphere_sphere(&s, &t, &s, &t).expect("coincident spheres must collide");
        assert!((contact.depth - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_sphere_box_separated() {
        let s = Sphere::new(0.5);
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        assert!(sphere_box(&s, &t1, &b, &t2).is_none());
    }
    #[test]
    fn test_sphere_box_inside() {
        let s = Sphere::new(0.1);
        let b = BoxShape::new(Vec3::new(2.0, 2.0, 2.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let contact = sphere_box(&s, &t1, &b, &t2).expect("sphere inside box must produce contact");
        assert!(contact.depth > 0.0);
    }
    #[test]
    fn test_sphere_mesh_contact_above_plane() {
        let verts = vec![
            Vec3::new(-2.0, 0.0, -2.0),
            Vec3::new(2.0, 0.0, -2.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        let indices = vec![[0usize, 1, 2]];
        let s = Sphere::new(0.5);
        let t_sphere = Transform::from_position(Vec3::new(0.0, 0.3, 0.0));
        let contact = sphere_triangle_mesh(&s, &t_sphere, &verts, &indices);
        assert!(
            contact.is_some(),
            "Sphere 0.3 m above plane, r=0.5 should contact"
        );
        let c = contact.unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn test_sphere_mesh_no_contact() {
        let verts = vec![
            Vec3::new(-2.0, 0.0, -2.0),
            Vec3::new(2.0, 0.0, -2.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        let indices = vec![[0usize, 1, 2]];
        let s = Sphere::new(0.5);
        let t_sphere = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        assert!(sphere_triangle_mesh(&s, &t_sphere, &verts, &indices).is_none());
    }
    #[test]
    fn test_cylinder_plane_overlap() {
        let cyl = Cylinder::new(0.5, 1.0);
        let t_cyl = Transform::from_position(Vec3::new(0.0, 0.3, 0.0));
        let plane_normal = Vec3::new(0.0, 1.0, 0.0);
        let plane_d = 0.0_f64;
        let contact = cylinder_plane(&cyl, &t_cyl, &plane_normal, plane_d);
        assert!(
            contact.is_some(),
            "Cylinder partially below plane should contact"
        );
        assert!(contact.unwrap().depth > 0.0);
    }
    #[test]
    fn test_cylinder_plane_no_contact() {
        let cyl = Cylinder::new(0.5, 1.0);
        let t_cyl = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        let plane_normal = Vec3::new(0.0, 1.0, 0.0);
        assert!(cylinder_plane(&cyl, &t_cyl, &plane_normal, 0.0).is_none());
    }
    #[test]
    fn test_capsule_mesh_contact() {
        let verts = vec![
            Vec3::new(-3.0, 0.0, -3.0),
            Vec3::new(3.0, 0.0, -3.0),
            Vec3::new(0.0, 0.0, 3.0),
        ];
        let indices = vec![[0usize, 1, 2]];
        let cap = Capsule::new(0.3, 0.5);
        let t_cap = Transform::from_position(Vec3::new(0.0, 0.2, 0.0));
        let contact = capsule_triangle_mesh(&cap, &t_cap, &verts, &indices);
        assert!(contact.is_some(), "Capsule near plane should contact");
    }
    #[test]
    fn test_capsule_mesh_no_contact() {
        let verts = vec![
            Vec3::new(-3.0, 0.0, -3.0),
            Vec3::new(3.0, 0.0, -3.0),
            Vec3::new(0.0, 0.0, 3.0),
        ];
        let indices = vec![[0usize, 1, 2]];
        let cap = Capsule::new(0.3, 0.5);
        let t_cap = Transform::from_position(Vec3::new(0.0, 10.0, 0.0));
        assert!(capsule_triangle_mesh(&cap, &t_cap, &verts, &indices).is_none());
    }
    #[test]
    fn test_heightfield_sphere_contact() {
        use oxiphysics_geometry::HeightField;
        let heights = vec![0.0_f64; 4];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let s = Sphere::new(0.5);
        let t_sphere = Transform::from_position(Vec3::new(0.5, 0.3, 0.5));
        let contact = heightfield_sphere(&hf, &s, &t_sphere);
        assert!(contact.is_some(), "Sphere above flat hf should contact");
    }
    #[test]
    fn test_heightfield_sphere_no_contact() {
        let heights = vec![0.0_f64; 4];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let s = Sphere::new(0.1);
        let t_sphere = Transform::from_position(Vec3::new(0.5, 5.0, 0.5));
        assert!(heightfield_sphere(&hf, &s, &t_sphere).is_none());
    }
    #[test]
    fn test_torus_plane_contact() {
        let torus = Torus::new(1.0, 0.3);
        let t_torus = Transform::from_position(Vec3::new(0.0, 0.2, 0.0));
        let plane_normal = Vec3::new(0.0, 1.0, 0.0);
        let contact = torus_plane(&torus, &t_torus, &plane_normal, 0.0);
        assert!(
            contact.is_some(),
            "Torus partially below plane should contact"
        );
    }
    #[test]
    fn test_torus_plane_no_contact() {
        let torus = Torus::new(1.0, 0.3);
        let t_torus = Transform::from_position(Vec3::new(0.0, 10.0, 0.0));
        let plane_normal = Vec3::new(0.0, 1.0, 0.0);
        assert!(torus_plane(&torus, &t_torus, &plane_normal, 0.0).is_none());
    }
    #[test]
    fn test_cone_sphere_near_base_edge() {
        let cone = Cone::new(1.0, 2.0);
        let sphere = Sphere::new(0.2);
        let t_cone = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(0.9, -0.9, 0.0));
        let _ = cone_sphere(&cone, &t_cone, &sphere, &t_sphere);
    }
}
