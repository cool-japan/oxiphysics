//! Triangulated test surfaces shared by the topology tests.
//!
//! Every builder returns `(vertex positions, counter-clockwise triangles)`;
//! closed convex solids are oriented outward and every surface is a valid
//! oriented 2-manifold.

use std::f64::consts::TAU;

/// Vertex positions and triangles of a test surface.
pub(super) type Fixture = (Vec<[f64; 3]>, Vec<[usize; 3]>);

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Flip triangles whose normal points towards the centroid (outward
/// orientation for a closed convex solid).
fn orient_outward(vertices: &[[f64; 3]], triangles: &mut [[usize; 3]]) {
    let n = vertices.len() as f64;
    let centroid = vertices.iter().fold([0.0; 3], |acc, v| {
        [acc[0] + v[0] / n, acc[1] + v[1] / n, acc[2] + v[2] / n]
    });
    for tri in triangles.iter_mut() {
        let [a, b, c] = [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]];
        let normal = cross(sub(b, a), sub(c, a));
        if dot(normal, sub(a, centroid)) < 0.0 {
            tri.swap(1, 2);
        }
    }
}

/// Regular tetrahedron: 4 vertices, 6 edges, 4 faces.
pub(super) fn tetrahedron() -> Fixture {
    let vertices = vec![
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let mut triangles = vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    orient_outward(&vertices, &mut triangles);
    (vertices, triangles)
}

/// Octahedron: 6 vertices (`+x -x +y -y +z -z`), 12 edges, 8 faces.
pub(super) fn octahedron() -> Fixture {
    let vertices = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let mut triangles = Vec::new();
    for x in [0, 1] {
        for y in [2, 3] {
            for z in [4, 5] {
                triangles.push([x, y, z]);
            }
        }
    }
    orient_outward(&vertices, &mut triangles);
    (vertices, triangles)
}

/// Icosahedron: 12 vertices, 30 edges, 20 faces.
pub(super) fn icosahedron() -> Fixture {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut vertices = Vec::new();
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            vertices.push([0.0, s1, s2 * phi]);
            vertices.push([s1, s2 * phi, 0.0]);
            vertices.push([s2 * phi, 0.0, s1]);
        }
    }
    let edge = 2.0;
    let adjacent = |i: usize, j: usize| {
        let d = sub(vertices[i], vertices[j]);
        (dot(d, d).sqrt() - edge).abs() < 1e-9
    };
    let mut triangles = Vec::new();
    for i in 0..vertices.len() {
        for j in (i + 1)..vertices.len() {
            for k in (j + 1)..vertices.len() {
                if adjacent(i, j) && adjacent(j, k) && adjacent(i, k) {
                    triangles.push([i, j, k]);
                }
            }
        }
    }
    orient_outward(&vertices, &mut triangles);
    (vertices, triangles)
}

/// `nu x nv` grid torus (two triangles per quad): `nu * nv` vertices,
/// `3 * nu * nv` edges, `2 * nu * nv` faces. Needs `nu, nv >= 3` to be a
/// simplicial complex; `4 x 4` is the smallest size with no short loops.
pub(super) fn grid_torus(nu: usize, nv: usize) -> Fixture {
    let (major, minor) = (3.0, 1.0);
    let mut vertices = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let theta = TAU * i as f64 / nu as f64;
        for j in 0..nv {
            let phi = TAU * j as f64 / nv as f64;
            let ring = major + minor * phi.cos();
            vertices.push([ring * theta.cos(), ring * theta.sin(), minor * phi.sin()]);
        }
    }
    let at = |i: usize, j: usize| (i % nu) * nv + (j % nv);
    let mut triangles = Vec::with_capacity(2 * nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    (vertices, triangles)
}

/// The 7-vertex Moebius-Csaszar torus: every pair of vertices is an edge
/// (7 vertices, 21 edges, 14 faces).
pub(super) fn seven_vertex_torus() -> Fixture {
    let vertices: Vec<[f64; 3]> = (0..7)
        .map(|i| {
            let angle = TAU * i as f64 / 7.0;
            [angle.cos(), angle.sin(), 0.2 * i as f64]
        })
        .collect();
    let mut triangles = Vec::new();
    for i in 0..7 {
        triangles.push([i, (i + 1) % 7, (i + 3) % 7]);
        triangles.push([i, (i + 3) % 7, (i + 2) % 7]);
    }
    (vertices, triangles)
}

/// Open tube of `rings` rings with `segments` vertices each (ring `r`, segment
/// `j` has index `r * segments + j`). Its two boundary loops are the first and
/// last ring. A sphere with two holes topologically.
pub(super) fn tube(rings: usize, segments: usize) -> Fixture {
    let mut vertices = Vec::with_capacity(rings * segments);
    for r in 0..rings {
        for j in 0..segments {
            let angle = TAU * j as f64 / segments as f64;
            vertices.push([angle.cos(), angle.sin(), r as f64]);
        }
    }
    let at = |r: usize, j: usize| r * segments + j % segments;
    let mut triangles = Vec::new();
    for r in 0..rings - 1 {
        for j in 0..segments {
            triangles.push([at(r, j), at(r, j + 1), at(r + 1, j + 1)]);
            triangles.push([at(r, j), at(r + 1, j + 1), at(r + 1, j)]);
        }
    }
    (vertices, triangles)
}

/// `rows x cols` vertex grid (`(rows - 1) * (cols - 1)` quads, two triangles
/// each). Cells alternate the diagonal direction when `alternate` is set; the
/// default diagonal runs from the lower-left to the upper-right vertex of every
/// cell, which joins two boundary vertices through an interior edge in the
/// corner cells.
pub(super) fn grid(rows: usize, cols: usize, alternate: bool) -> Fixture {
    let mut vertices = Vec::with_capacity(rows * cols);
    for r in 0..rows {
        for c in 0..cols {
            vertices.push([c as f64, r as f64, 0.0]);
        }
    }
    let at = |r: usize, c: usize| r * cols + c;
    let mut triangles = Vec::new();
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            let flip = alternate && (r + c) % 2 == 1;
            if flip {
                triangles.push([at(r, c), at(r, c + 1), at(r + 1, c)]);
                triangles.push([at(r, c + 1), at(r + 1, c + 1), at(r + 1, c)]);
            } else {
                triangles.push([at(r, c), at(r, c + 1), at(r + 1, c + 1)]);
                triangles.push([at(r, c), at(r + 1, c + 1), at(r + 1, c)]);
            }
        }
    }
    (vertices, triangles)
}

/// Disjoint union of two surfaces (vertex indices of `b` are shifted).
pub(super) fn disjoint_union(a: &Fixture, b: &Fixture) -> Fixture {
    let offset = a.0.len();
    let mut vertices = a.0.clone();
    vertices.extend(b.0.iter().map(|p| [p[0] + 20.0, p[1], p[2]]));
    let mut triangles = a.1.clone();
    triangles.extend(
        b.1.iter()
            .map(|t| [t[0] + offset, t[1] + offset, t[2] + offset]),
    );
    (vertices, triangles)
}

/// Connected sum: remove face 0 of each surface and glue the two holes along
/// their three vertices (with the orientation reversal an oriented gluing
/// needs). The vertices of `a`'s removed face are kept; `b`'s are identified
/// with them.
pub(super) fn connected_sum(a: &Fixture, b: &Fixture) -> Fixture {
    let [a0, a1, a2] = a.1[0];
    let [b0, b1, b2] = b.1[0];
    let mut vertices = a.0.clone();
    let mut map_b: Vec<Option<usize>> = vec![None; b.0.len()];
    map_b[b0] = Some(a0);
    map_b[b1] = Some(a2);
    map_b[b2] = Some(a1);
    for (v, slot) in map_b.iter_mut().enumerate() {
        if slot.is_none() {
            *slot = Some(vertices.len());
            let p = b.0[v];
            vertices.push([p[0] + 20.0, p[1], p[2]]);
        }
    }
    let mut triangles: Vec<[usize; 3]> = a.1[1..].to_vec();
    for tri in &b.1[1..] {
        let mapped: Vec<usize> = tri.iter().filter_map(|&v| map_b[v]).collect();
        triangles.push([mapped[0], mapped[1], mapped[2]]);
    }
    (vertices, triangles)
}

/// Remove the listed faces.
pub(super) fn without_faces(fixture: &Fixture, removed: &[usize]) -> Fixture {
    let triangles = fixture
        .1
        .iter()
        .enumerate()
        .filter(|(i, _)| !removed.contains(i))
        .map(|(_, t)| *t)
        .collect();
    (fixture.0.clone(), triangles)
}
