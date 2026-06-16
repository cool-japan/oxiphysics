//! Analytic coplanarity cubic time-of-impact (TOI) solver for continuous
//! collision detection (CCD) of deformable surfaces.
//!
//! This module implements the classic Provot (1997) / Bridson (2002) approach
//! to continuous collision detection between linearly-moving simplices. Both
//! the vertex-face and edge-edge primitive tests reduce to finding the earliest
//! time `t` in `[0, 1]` at which four moving points become *coplanar*.
//!
//! Assuming every point moves with constant velocity over the timestep,
//! `x_i(t) = p_i + t * v_i`, the coplanarity condition is the vanishing of a
//! scalar triple product of three relative position vectors. Because each
//! relative position is affine in `t`, the triple product expands into a cubic
//! polynomial `a*t^3 + b*t^2 + c*t + d = 0`. The real roots of this cubic in
//! `[0, 1]` are the candidate times of coplanarity; each is then validated
//! geometrically (barycentric containment for vertex-face, edge-parameter
//! containment for edge-edge) within a thickness tolerance.
//!
//! References:
//! - X. Provot, "Collision and self-collision handling in cloth model dedicated
//!   to design garments" (Computer Animation and Simulation, 1997).
//! - R. Bridson, R. Fedkiw, J. Anderson, "Robust treatment of collisions,
//!   contact and friction for cloth animation" (SIGGRAPH, 2002).
//!
//! The cubic itself is solved analytically using the trigonometric (Viète)
//! solution for the three-real-root case and Cardano's formula for the
//! single-real-root case, followed by a few Newton iterations to polish the
//! roots against the original (un-depressed) polynomial.
//!
//! Contacts that only graze within the collision `thickness` without ever
//! reaching exact coplanarity have no real root in the shell; for those the
//! solver falls back to a coarse prescan + bisection on the exact distance to
//! recover the earliest shell entry, so the analytic path is itself zero-miss.

use super::functions::{
    ccd_prescan_bisect, closest_point_on_triangle, closest_points_edge_edge, v3_add, v3_cross,
    v3_dot, v3_scale, v3_sub,
};
use core::f64::consts::PI;

/// Newton-polish candidate roots against the general cubic, filter to `[0, 1]`
/// (with a small tolerance), clamp, sort ascending, and deduplicate.
fn finalize_roots(mut cands: Vec<f64>, a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    for t in cands.iter_mut() {
        for _ in 0..3 {
            let f = a * (*t) * (*t) * (*t) + b * (*t) * (*t) + c * (*t) + d;
            let fp = 3.0 * a * (*t) * (*t) + 2.0 * b * (*t) + c;
            if fp.abs() > 1e-14 {
                *t -= f / fp;
            }
        }
    }

    let mut kept: Vec<f64> = cands
        .into_iter()
        .filter(|t| (-1e-8..=1.0 + 1e-8).contains(t))
        .map(|t| t.clamp(0.0, 1.0))
        .collect();

    kept.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));

    let mut out: Vec<f64> = Vec::new();
    for t in kept {
        match out.last() {
            Some(&last) if (t - last).abs() <= 1e-6 => {}
            _ => out.push(t),
        }
    }
    out
}

/// Solve `a2*t^2 + b2*t + c2 = 0` (with linear fallback), returning real roots
/// in `[0, 1]` after polishing/clamping.
fn solve_quadratic_in_01(a2: f64, b2: f64, c2: f64) -> Vec<f64> {
    if a2.abs() < 1e-12 * b2.abs().max(c2.abs()).max(1e-15) {
        if b2.abs() < 1e-15 {
            return Vec::new();
        }
        let cand = -c2 / b2;
        return finalize_roots(vec![cand], 0.0, a2, b2, c2);
    }

    let disc = b2 * b2 - 4.0 * a2 * c2;
    if disc < 0.0 {
        return Vec::new();
    }

    let sq = disc.max(0.0).sqrt();
    let q = -0.5 * (b2 + b2.signum() * sq);
    let mut cands: Vec<f64> = Vec::new();
    if q.abs() > 1e-300 {
        cands.push(q / a2);
        cands.push(c2 / q);
    } else {
        cands.push(-b2 / (2.0 * a2));
    }
    finalize_roots(cands, 0.0, a2, b2, c2)
}

/// Solve a*t^3 + b*t^2 + c*t + d = 0, returning all real roots in [0,1]
/// (after clamping/tolerance), sorted ascending and deduplicated.
pub fn solve_cubic_in_01(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    let scale = b.abs().max(c.abs()).max(d.abs()).max(1e-15);
    if a.abs() < 1e-12 * scale {
        return solve_quadratic_in_01(b, c, d);
    }

    let p = b / a;
    let q = c / a;
    let r = d / a;
    let big_p = q - p * p / 3.0;
    let big_q = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;
    let shift = -p / 3.0;
    let delta = -4.0 * big_p * big_p * big_p - 27.0 * big_q * big_q;
    let delta_eps = 1e-9 * (big_p * big_p * big_p).abs().max(big_q * big_q).max(1e-30);

    let mut ys: Vec<f64> = Vec::new();
    if delta > delta_eps {
        let m = 2.0 * (-big_p / 3.0).sqrt();
        let acos_arg = ((3.0 * big_q) / (2.0 * big_p) * (-3.0 / big_p).sqrt()).clamp(-1.0, 1.0);
        let theta = acos_arg.acos();
        for k in 0..3 {
            let yk = m * (theta / 3.0 - 2.0 * PI * (k as f64) / 3.0).cos();
            ys.push(yk);
        }
    } else if delta.abs() <= delta_eps {
        if big_p.abs() < 1e-12 {
            ys.push(0.0);
        } else {
            let y_simple = 3.0 * big_q / big_p;
            let y_double = -3.0 * big_q / (2.0 * big_p);
            ys.push(y_simple);
            ys.push(y_double);
        }
    } else {
        let half_q = big_q / 2.0;
        let inner = (big_q * big_q / 4.0 + big_p * big_p * big_p / 27.0).max(0.0);
        let cube = inner.sqrt();
        let u = (-half_q + cube).cbrt();
        let v = (-half_q - cube).cbrt();
        let y = u + v;
        ys.push(y);
    }

    let cands: Vec<f64> = ys.iter().map(|y| y + shift).collect();
    finalize_roots(cands, a, b, c, d)
}

/// Compute the cubic coefficients of the scalar triple product
/// `(r1 + t*w1) . ((r2 + t*w2) x (r3 + t*w3))` expanded in `t`.
///
/// Returns `(a, b, c, d)` for `a*t^3 + b*t^2 + c*t + d`.
fn triple_product_cubic(
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    w1: [f64; 3],
    w2: [f64; 3],
    w3: [f64; 3],
) -> (f64, f64, f64, f64) {
    let a = v3_dot(w1, v3_cross(w2, w3));
    let b = v3_dot(r1, v3_cross(w2, w3))
        + v3_dot(w1, v3_cross(r2, w3))
        + v3_dot(w1, v3_cross(w2, r3));
    let c = v3_dot(r1, v3_cross(r2, w3))
        + v3_dot(r1, v3_cross(w2, r3))
        + v3_dot(w1, v3_cross(r2, r3));
    let d = v3_dot(r1, v3_cross(r2, r3));
    (a, b, c, d)
}

/// Analytic vertex-face continuous CCD time-of-impact.
/// p4/v4 = the moving vertex; (p1,p2,p3)/(v1,v2,v3) = the moving triangle.
/// Returns earliest t in [0,1] where the vertex becomes coplanar with the
/// triangle AND projects inside it (within `thickness`), else None.
pub fn cubic_toi_vertex_face(
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    p4: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    // t=0 early-out: primitives already within thickness => immediate contact.
    let (_c0, _b0, dsq0) = closest_point_on_triangle(p4, p1, p2, p3);
    if dsq0.sqrt() < thickness {
        return Some(0.0);
    }
    let r1 = v3_sub(p1, p4);
    let r2 = v3_sub(p2, p4);
    let r3 = v3_sub(p3, p4);
    let w1 = v3_sub(v1, v4);
    let w2 = v3_sub(v2, v4);
    let w3 = v3_sub(v3, v4);
    let (a, b, c, d) = triple_product_cubic(r1, r2, r3, w1, w2, w3);
    let roots = solve_cubic_in_01(a, b, c, d);
    let bary_tol = 1e-6;
    // Primary path: earliest coplanarity crossing that projects inside the
    // triangle and is within the thickness shell (distance ~= 0 at a true root).
    for &t in &roots {
        let vp = v3_add(p4, v3_scale(v4, t));
        let at = v3_add(p1, v3_scale(v1, t));
        let bt = v3_add(p2, v3_scale(v2, t));
        let ct = v3_add(p3, v3_scale(v3, t));
        let (_closest, bary, dist_sq) = closest_point_on_triangle(vp, at, bt, ct);
        let inside = bary[0] >= -bary_tol && bary[1] >= -bary_tol && bary[2] >= -bary_tol;
        let within_shell = dist_sq.sqrt() <= thickness;
        if inside && within_shell {
            return Some(t);
        }
    }
    // Thickness-graze backstop: a vertex can dip within `thickness` of the face
    // without the relative motion ever reaching exact coplanarity (so the cubic
    // has no real root in the shell). The coplanarity cubic cannot express this
    // case analytically (the exact distance-shell boundary is degree six), so we
    // fall back to the shared analytic-free prescan+bisection on the exact
    // closest distance to recover the earliest shell entry.
    ccd_prescan_bisect(thickness, |t| {
        let vp = v3_add(p4, v3_scale(v4, t));
        let at = v3_add(p1, v3_scale(v1, t));
        let bt = v3_add(p2, v3_scale(v2, t));
        let ct = v3_add(p3, v3_scale(v3, t));
        let (_closest, _bary, dist_sq) = closest_point_on_triangle(vp, at, bt, ct);
        dist_sq.sqrt()
    })
}

/// Analytic edge-edge continuous CCD time-of-impact.
/// Edge 1 = p1->p2 (vel v1->v2), Edge 2 = p3->p4 (vel v3->v4).
pub fn cubic_toi_edge_edge(
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    p4: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    // t=0 early-out: edges already within thickness => immediate contact.
    let (_s0, _t0, _cp0, _cq0, dsq0) = closest_points_edge_edge(p1, p2, p3, p4);
    if dsq0.sqrt() < thickness {
        return Some(0.0);
    }
    let r1 = v3_sub(p2, p1);
    let w1 = v3_sub(v2, v1);
    let r2 = v3_sub(p4, p3);
    let w2 = v3_sub(v4, v3);
    let r3 = v3_sub(p3, p1);
    let w3 = v3_sub(v3, v1);
    let (a, b, c, d) = triple_product_cubic(r1, r2, r3, w1, w2, w3);
    let roots = solve_cubic_in_01(a, b, c, d);
    let param_tol = 1e-6;
    // Primary path: earliest coplanarity crossing whose closest points lie on
    // both segments and within the thickness shell.
    for &t in &roots {
        let a_pt = v3_add(p1, v3_scale(v1, t));
        let b_pt = v3_add(p2, v3_scale(v2, t));
        let c_pt = v3_add(p3, v3_scale(v3, t));
        let e_pt = v3_add(p4, v3_scale(v4, t));
        let (s, tp, _cp, _cq, dsq) = closest_points_edge_edge(a_pt, b_pt, c_pt, e_pt);
        let in_params = (-param_tol..=1.0 + param_tol).contains(&s)
            && (-param_tol..=1.0 + param_tol).contains(&tp);
        let within_shell = dsq.sqrt() <= thickness;
        if in_params && within_shell {
            return Some(t);
        }
    }
    // Thickness-graze backstop (see `cubic_toi_vertex_face`): recover shell entry
    // for edge pairs that pass within `thickness` without reaching coplanarity.
    ccd_prescan_bisect(thickness, |t| {
        let a_pt = v3_add(p1, v3_scale(v1, t));
        let b_pt = v3_add(p2, v3_scale(v2, t));
        let c_pt = v3_add(p3, v3_scale(v3, t));
        let e_pt = v3_add(p4, v3_scale(v4, t));
        let (_s, _tp, _cp, _cq, dsq) = closest_points_edge_edge(a_pt, b_pt, c_pt, e_pt);
        dsq.sqrt()
    })
}
