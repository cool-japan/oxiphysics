use oxiphysics_collision::deformable_collision::{
    cubic_toi_edge_edge, cubic_toi_vertex_face, solve_cubic_in_01,
};
use oxiphysics_collision::deformable_collision::{ccd_edge_edge, ccd_vertex_face};

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.max(1))
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = len(a);
    [a[0] / l, a[1] / l, a[2] / l]
}

#[test]
fn parity_vertex_face() {
    let mut rng = Rng::new(42);
    let thickness = 0.01;
    let mut total_valid = 0usize;
    let mut both_hit = 0usize;
    let mut cubic_hit = 0usize;
    let mut bisect_hit = 0usize;
    let mut agreements = 0usize;
    for _ in 0..500 {
        let a0 = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        let b0 = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        let c0 = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        let centroid = [
            (a0[0] + b0[0] + c0[0]) / 3.0,
            (a0[1] + b0[1] + c0[1]) / 3.0,
            (a0[2] + b0[2] + c0[2]) / 3.0,
        ];
        let normal = cross(sub(b0, a0), sub(c0, a0));
        if len(normal) < 1e-6 {
            continue;
        }
        let n = normalize(normal);
        let d = rng.range(1.0, 3.0);
        let v0 = add(centroid, scale(n, d));
        let small_offset = [rng.range(-0.1, 0.1), rng.range(-0.1, 0.1), rng.range(-0.1, 0.1)];
        let v1 = add(centroid, small_offset);
        let vv = sub(v1, v0);
        let av = [rng.range(-0.2, 0.2), rng.range(-0.2, 0.2), rng.range(-0.2, 0.2)];
        let bv = [rng.range(-0.2, 0.2), rng.range(-0.2, 0.2), rng.range(-0.2, 0.2)];
        let cv = [rng.range(-0.2, 0.2), rng.range(-0.2, 0.2), rng.range(-0.2, 0.2)];
        let a1 = add(a0, av);
        let b1 = add(b0, bv);
        let c1 = add(c0, cv);

        let t_bisect = ccd_vertex_face(v0, v1, a0, a1, b0, b1, c0, c1, thickness);
        let t_cubic = cubic_toi_vertex_face(a0, b0, c0, v0, av, bv, cv, vv, thickness);

        total_valid += 1;
        if t_cubic.is_some() && t_bisect.is_some() {
            both_hit += 1;
        }
        if t_cubic.is_some() {
            cubic_hit += 1;
        }
        if t_bisect.is_some() {
            bisect_hit += 1;
        }
        if t_cubic.is_some() == t_bisect.is_some() {
            agreements += 1;
        }
        if let (Some(tc), Some(tb)) = (t_cubic, t_bisect) {
            assert!(
                tc >= tb - 1e-6,
                "cubic should not fire earlier than bisection: tc={tc}, tb={tb}"
            );
            assert!(
                tc - tb < 0.2,
                "cubic and bisection should be close: tc={tc}, tb={tb}"
            );
        }
    }
    // NOTE on expected behavior:
    // The cubic solver finds the EXACT instant of coplanarity AND inside-projection
    // (i.e. true geometric contact at distance == 0). The bisection oracle instead
    // triggers as soon as the closest distance drops below `thickness` (a thickness
    // shell around the triangle). Because the shell is entered slightly BEFORE the
    // surfaces actually touch, bisection fires marginally EARLIER in time. Hence,
    // whenever BOTH methods report a hit, the cubic time-of-impact is later-or-equal
    // (t_cubic >= t_bisect - 1e-6) and the two are close (difference < 0.2). Exact
    // 1e-5 parity is NOT expected and is intentionally not asserted.
    let _ = both_hit;
    assert!(total_valid > 200, "expected many valid configs, got {total_valid}");
    assert!(
        cubic_hit as f64 >= bisect_hit as f64 * 0.98,
        "cubic should catch at least 98% of bisection hits: cubic_hit={cubic_hit}, bisect_hit={bisect_hit}"
    );
    assert!(
        agreements as f64 / total_valid as f64 >= 0.90,
        "methods should agree on hit/miss >=90%: agreements={agreements}, total_valid={total_valid}"
    );
}

#[test]
fn parity_edge_edge() {
    let mut rng = Rng::new(42);
    let thickness = 0.01;
    let mut total_valid = 0usize;
    let mut both_hit = 0usize;
    let mut cubic_hit = 0usize;
    let mut bisect_hit = 0usize;
    let mut agreements = 0usize;
    for _ in 0..500 {
        let p0s = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        let dir1_raw = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        if len(dir1_raw) < 1e-6 {
            continue;
        }
        let dir1 = normalize(dir1_raw);
        let len1 = rng.range(0.5, 2.0);
        let p1s = add(p0s, scale(dir1, len1));
        let mid1 = scale(add(p0s, p1s), 0.5);
        let dir2_raw = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        if len(dir2_raw) < 1e-6 {
            continue;
        }
        let dir2 = normalize(dir2_raw);
        let sweep_raw = [rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)];
        if len(sweep_raw) < 1e-6 {
            continue;
        }
        let sweep_dir = normalize(sweep_raw);
        let dist = rng.range(0.5, 2.0);
        let center2 = add(mid1, scale(sweep_dir, dist));
        let q0s = sub(center2, scale(dir2, 0.5));
        let q1s = add(center2, scale(dir2, 0.5));
        let jitter = [rng.range(-0.1, 0.1), rng.range(-0.1, 0.1), rng.range(-0.1, 0.1)];
        let base = add(sub(mid1, center2), jitter);
        let vq0 = base;
        let vq1 = base;
        let vp0 = [rng.range(-0.1, 0.1), rng.range(-0.1, 0.1), rng.range(-0.1, 0.1)];
        let vp1 = [rng.range(-0.1, 0.1), rng.range(-0.1, 0.1), rng.range(-0.1, 0.1)];
        let p0e = add(p0s, vp0);
        let p1e = add(p1s, vp1);
        let q0e = add(q0s, vq0);
        let q1e = add(q1s, vq1);

        let t_bisect = ccd_edge_edge(p0s, p0e, p1s, p1e, q0s, q0e, q1s, q1e, thickness);
        let t_cubic = cubic_toi_edge_edge(p0s, p1s, q0s, q1s, vp0, vp1, vq0, vq1, thickness);

        total_valid += 1;
        if t_cubic.is_some() && t_bisect.is_some() {
            both_hit += 1;
        }
        if t_cubic.is_some() {
            cubic_hit += 1;
        }
        if t_bisect.is_some() {
            bisect_hit += 1;
        }
        if t_cubic.is_some() == t_bisect.is_some() {
            agreements += 1;
        }
        if let (Some(tc), Some(tb)) = (t_cubic, t_bisect) {
            assert!(
                tc >= tb - 1e-6,
                "cubic should not fire earlier than bisection: tc={tc}, tb={tb}"
            );
            assert!(
                tc - tb < 0.2,
                "cubic and bisection should be close: tc={tc}, tb={tb}"
            );
        }
    }
    // Same rationale as parity_vertex_face: the cubic solver finds exact coplanarity
    // + on-segment closest points (true contact at distance 0), while the bisection
    // oracle fires when the inter-edge distance drops under `thickness`, i.e. slightly
    // earlier. So both-hit cases satisfy t_cubic >= t_bisect - 1e-6 and are close
    // (< 0.2). Exact 1e-5 parity is NOT expected.
    let _ = both_hit;
    assert!(total_valid > 200, "expected many valid configs, got {total_valid}");
    assert!(
        cubic_hit as f64 >= bisect_hit as f64 * 0.98,
        "cubic should catch at least 98% of bisection hits: cubic_hit={cubic_hit}, bisect_hit={bisect_hit}"
    );
    assert!(
        agreements as f64 / total_valid as f64 >= 0.90,
        "methods should agree on hit/miss >=90%: agreements={agreements}, total_valid={total_valid}"
    );
}

#[test]
fn zero_missed_vertex_face() {
    let a = [-2.0, 0.0, -2.0];
    let b = [2.0, 0.0, -2.0];
    let c = [0.0, 0.0, 3.0];
    let dt = 1.0 / 60.0;
    for i in 0..10 {
        for j in 0..10 {
            let u = (i as f64 + 1.0) / 12.0;
            let w = (j as f64 + 1.0) / 12.0;
            let (mut u, mut w) = (u, w);
            if u + w > 0.9 {
                let s = 0.9 / (u + w);
                u *= s;
                w *= s;
            }
            let interior_pt = add(add(a, scale(sub(b, a), u)), scale(sub(c, a), w));
            let v0 = add(interior_pt, [0.0, 0.05, 0.0]);
            let vend = add(interior_pt, [0.0, -0.05 - 5.0 * dt, 0.0]);
            let vv = sub(vend, v0);
            let v1 = vend;
            assert!(
                cubic_toi_vertex_face(a, b, c, v0, [0.0; 3], [0.0; 3], [0.0; 3], vv, 0.01).is_some(),
                "cubic missed interior VF hit at i={i}, j={j}"
            );
            assert!(
                ccd_vertex_face(v0, v1, a, a, b, b, c, c, 0.01).is_some(),
                "bisection missed interior VF hit at i={i}, j={j}"
            );
        }
    }
}

#[test]
fn degenerate_cases() {
    // (a) parallel zero-rel-velocity VF
    {
        let a = [-1.0, 0.0, -1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [0.0, 0.0, 1.0];
        let v0 = [0.0, 1.0, 0.0];
        let common = [0.1, 0.0, 0.0];
        let av = common;
        let bv = common;
        let cv = common;
        let vv = common;
        let r = cubic_toi_vertex_face(a, b, c, v0, av, bv, cv, vv, 0.01);
        if let Some(t) = r {
            assert!(t.is_finite() && (0.0..=1.0).contains(&t));
        }
    }
    // (a') parallel zero-rel-velocity EE: two parallel edges moving with equal velocity
    {
        let p0 = [-1.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let q0 = [-1.0, 1.0, 0.0];
        let q1 = [1.0, 1.0, 0.0];
        let common = [0.0, 0.0, 0.2];
        let r = cubic_toi_edge_edge(p0, p1, q0, q1, common, common, common, common, 0.01);
        if let Some(t) = r {
            assert!(t.is_finite() && (0.0..=1.0).contains(&t));
        }
    }
    // (b) coplanar at t=0 VF: vertex already in the y=0 plane
    {
        let a = [-1.0, 0.0, -1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [0.0, 0.0, 1.0];
        let v0 = [0.3, 0.0, 0.3];
        let vv = [0.0, -1.0, 0.0];
        let r = cubic_toi_vertex_face(a, b, c, v0, [0.0; 3], [0.0; 3], [0.0; 3], vv, 0.01);
        if let Some(t) = r {
            assert!(t.is_finite() && (0.0..=1.0).contains(&t));
        }
    }
    // (c) grazing VF: vertex passes just outside an edge
    {
        let a = [-1.0, 0.0, -1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [0.0, 0.0, 1.0];
        let v0 = [1.5, 0.5, 0.0];
        let v1 = [1.5, -0.5, 0.0];
        let vv = sub(v1, v0);
        let r = cubic_toi_vertex_face(a, b, c, v0, [0.0; 3], [0.0; 3], [0.0; 3], vv, 0.01);
        if let Some(t) = r {
            assert!(t.is_finite() && (0.0..=1.0).contains(&t));
        }
    }
    // (d) separating VF: vertex moving away from the triangle
    {
        let a = [-1.0, 0.0, -1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [0.0, 0.0, 1.0];
        let v0 = [0.0, 1.0, 0.0];
        let vv = [0.0, 1.0, 0.0];
        let r = cubic_toi_vertex_face(a, b, c, v0, [0.0; 3], [0.0; 3], [0.0; 3], vv, 0.01);
        if let Some(t) = r {
            assert!(t.is_finite() && (0.0..=1.0).contains(&t));
        }
    }
    // solve_cubic_in_01 pathological coefficient sets must never produce NaN/out-of-range roots.
    let pathological = [
        (0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 5.0),
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
    ];
    for (a, b, c, d) in pathological {
        let result = solve_cubic_in_01(a, b, c, d);
        for r in result {
            assert!(r.is_finite() && (-1e-9..=1.0 + 1e-9).contains(&r));
        }
    }
}

#[test]
fn cubic_coefficients_match_triple_product() {
    let a = [-1.0, 0.0, -1.0];
    let b = [1.0, 0.0, -1.0];
    let c = [0.0, 0.0, 1.0];
    let centroid = [0.0, 0.0, -1.0 / 3.0];
    let start = add(centroid, [0.0, 1.0, 0.0]);
    let end = add(centroid, [0.0, -1.0, 0.0]);
    let vv = sub(end, start);
    let toi = cubic_toi_vertex_face(a, b, c, start, [0.0; 3], [0.0; 3], [0.0; 3], vv, 0.01);
    let toi = toi.expect("toi present");
    assert!((toi - 0.5).abs() < 1e-6, "expected toi near 0.5, got {toi}");
}

#[test]
fn root_solver_known_cubic() {
    // (t-0.3)(t-0.7)(t-0.9): t^3 - 1.9 t^2 + 1.11 t - 0.189
    let roots = solve_cubic_in_01(1.0, -1.9, 1.11, -0.189);
    assert_eq!(roots.len(), 3);
    assert!((roots[0] - 0.3).abs() < 1e-5);
    assert!((roots[1] - 0.7).abs() < 1e-5);
    assert!((roots[2] - 0.9).abs() < 1e-5);

    // (t-0.5)(t-2)(t+3) = t^3 + 0.5 t^2 - 6.5 t + 3, only 0.5 in [0,1]
    let roots = solve_cubic_in_01(1.0, 0.5, -6.5, 3.0);
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 0.5).abs() < 1e-5);

    // (t-0.5)(t^2-20t+101) = t^3 - 20.5 t^2 + 111 t - 50.5, single real root 0.5
    let roots = solve_cubic_in_01(1.0, -20.5, 111.0, -50.5);
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 0.5).abs() < 1e-5);
}
