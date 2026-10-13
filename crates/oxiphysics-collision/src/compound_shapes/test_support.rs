//! Triangle meshes shared by the compound-shape tests.

/// Vertices and triangle indices of a test mesh.
pub(super) type Mesh = (Vec<[f64; 3]>, Vec<[usize; 3]>);

/// Cube `[-half, half]^3` translated by `offset`, with 12 outward-facing
/// (counter-clockwise seen from outside) triangles. Vertex `i` has the sign
/// pattern `(x > 0) + 2 (y > 0) + 4 (z > 0)`.
pub(super) fn cube_mesh(half: f64, offset: [f64; 3]) -> Mesh {
    let vertices = (0..8)
        .map(|i| {
            let sign = |bit: usize| if (i >> bit) & 1 == 1 { half } else { -half };
            [
                offset[0] + sign(0),
                offset[1] + sign(1),
                offset[2] + sign(2),
            ]
        })
        .collect();
    let indices = vec![
        [1, 3, 7],
        [1, 7, 5],
        [0, 4, 6],
        [0, 6, 2],
        [2, 6, 7],
        [2, 7, 3],
        [0, 1, 5],
        [0, 5, 4],
        [4, 5, 7],
        [4, 7, 6],
        [0, 2, 3],
        [0, 3, 1],
    ];
    (vertices, indices)
}

/// Regular icosahedron with the given circumradius and 20 outward-facing
/// triangles.
pub(super) fn icosahedron(radius: f64) -> Mesh {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let norm = (1.0 + phi * phi).sqrt();
    let s = radius / norm;
    let raw = [
        [-1.0, phi, 0.0],
        [1.0, phi, 0.0],
        [-1.0, -phi, 0.0],
        [1.0, -phi, 0.0],
        [0.0, -1.0, phi],
        [0.0, 1.0, phi],
        [0.0, -1.0, -phi],
        [0.0, 1.0, -phi],
        [phi, 0.0, -1.0],
        [phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
    ];
    let vertices = raw.iter().map(|v| [v[0] * s, v[1] * s, v[2] * s]).collect();
    let indices = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    (vertices, indices)
}

/// Icosahedron with every triangle split into four and all vertices projected
/// back onto the sphere: 42 vertices, all of them extreme points.
pub(super) fn icosphere(radius: f64) -> Mesh {
    let (mut vertices, base) = icosahedron(radius);
    let mut indices = Vec::new();
    let mut midpoints: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    let mut midpoint = |vertices: &mut Vec<[f64; 3]>, a: usize, b: usize| -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&index) = midpoints.get(&key) {
            return index;
        }
        let (pa, pb) = (vertices[a], vertices[b]);
        let m = [
            0.5 * (pa[0] + pb[0]),
            0.5 * (pa[1] + pb[1]),
            0.5 * (pa[2] + pb[2]),
        ];
        let len = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
        vertices.push([
            m[0] * radius / len,
            m[1] * radius / len,
            m[2] * radius / len,
        ]);
        let index = vertices.len() - 1;
        midpoints.insert(key, index);
        index
    };
    for tri in &base {
        let ab = midpoint(&mut vertices, tri[0], tri[1]);
        let bc = midpoint(&mut vertices, tri[1], tri[2]);
        let ca = midpoint(&mut vertices, tri[2], tri[0]);
        indices.push([tri[0], ab, ca]);
        indices.push([tri[1], bc, ab]);
        indices.push([tri[2], ca, bc]);
        indices.push([ab, bc, ca]);
    }
    (vertices, indices)
}

/// Ear-clipping triangulation of a simple counter-clockwise polygon.
fn triangulate(polygon: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| {
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    };
    let inside = |p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        cross(a, b, p) >= 0.0 && cross(b, c, p) >= 0.0 && cross(c, a, p) >= 0.0
    };
    let mut remaining: Vec<usize> = (0..polygon.len()).collect();
    let mut triangles = Vec::new();
    while remaining.len() > 3 {
        let n = remaining.len();
        let ear = (0..n).find(|&i| {
            let (ia, ib, ic) = (
                remaining[(i + n - 1) % n],
                remaining[i],
                remaining[(i + 1) % n],
            );
            let (a, b, c) = (polygon[ia], polygon[ib], polygon[ic]);
            cross(a, b, c) > 1e-12
                && !remaining
                    .iter()
                    .any(|&k| k != ia && k != ib && k != ic && inside(polygon[k], a, b, c))
        });
        let Some(i) = ear else {
            break;
        };
        triangles.push([
            remaining[(i + n - 1) % n],
            remaining[i],
            remaining[(i + 1) % n],
        ]);
        remaining.remove(i);
    }
    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }
    triangles
}

/// Prism over a simple counter-clockwise polygon in the xy-plane, from `z = 0`
/// to `z = height`, with outward-facing triangles. Vertex `i` is the polygon
/// corner `i` at the bottom, vertex `n + i` the one at the top.
pub(super) fn extrude_polygon(polygon: &[[f64; 2]], height: f64) -> Mesh {
    let n = polygon.len();
    let mut vertices: Vec<[f64; 3]> = polygon.iter().map(|p| [p[0], p[1], 0.0]).collect();
    vertices.extend(polygon.iter().map(|p| [p[0], p[1], height]));
    let mut indices = Vec::new();
    for [a, b, c] in triangulate(polygon) {
        indices.push([n + a, n + b, n + c]);
        indices.push([a, c, b]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.push([i, j, n + j]);
        indices.push([i, n + j, n + i]);
    }
    (vertices, indices)
}

/// Torus around the z axis: `major` is the radius of the centre circle of the
/// tube and `minor` the radius of the tube, with `u_steps` segments around the
/// axis and `v_steps` around the tube. It has `u_steps * v_steps` vertices and
/// twice as many outward-facing triangles; vertex `(i, j)` has the index
/// `i * v_steps + j`.
pub(super) fn torus(major: f64, minor: f64, u_steps: usize, v_steps: usize) -> Mesh {
    let mut vertices = Vec::with_capacity(u_steps * v_steps);
    for i in 0..u_steps {
        let u = std::f64::consts::TAU * i as f64 / u_steps as f64;
        for j in 0..v_steps {
            let v = std::f64::consts::TAU * j as f64 / v_steps as f64;
            let ring = major + minor * v.cos();
            vertices.push([ring * u.cos(), ring * u.sin(), minor * v.sin()]);
        }
    }
    let index = |i: usize, j: usize| (i % u_steps) * v_steps + (j % v_steps);
    let mut indices = Vec::with_capacity(2 * u_steps * v_steps);
    for i in 0..u_steps {
        for j in 0..v_steps {
            let a = index(i, j);
            let b = index(i + 1, j);
            let c = index(i + 1, j + 1);
            let d = index(i, j + 1);
            indices.push([a, b, c]);
            indices.push([a, c, d]);
        }
    }
    (vertices, indices)
}

/// The same mesh with the winding of every triangle reversed.
pub(super) fn flipped(mesh: &Mesh) -> Mesh {
    (
        mesh.0.clone(),
        mesh.1.iter().map(|t| [t[0], t[2], t[1]]).collect(),
    )
}

/// Signed volume enclosed by a mesh (positive for outward-facing triangles).
pub(super) fn signed_volume(mesh: &Mesh) -> f64 {
    mesh.1
        .iter()
        .map(|t| {
            let (a, b, c) = (mesh.0[t[0]], mesh.0[t[1]], mesh.0[t[2]]);
            let cross = [
                b[1] * c[2] - b[2] * c[1],
                b[2] * c[0] - b[0] * c[2],
                b[0] * c[1] - b[1] * c[0],
            ];
            (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]) / 6.0
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube_has_outward_triangles_and_volume() {
        let mesh = cube_mesh(1.0, [0.0; 3]);
        assert!((signed_volume(&mesh) - 8.0).abs() < 1e-12);
    }

    #[test]
    fn icosahedron_and_icosphere_are_closed_outward_meshes() {
        for mesh in [icosahedron(1.0), icosphere(1.0)] {
            assert!(signed_volume(&mesh) > 0.0);
            // Every directed edge must be matched by its reverse exactly once.
            let mut edges = std::collections::HashMap::new();
            for t in &mesh.1 {
                for k in 0..3 {
                    *edges.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
                }
            }
            for (&(a, b), &count) in &edges {
                assert_eq!(count, 1);
                assert_eq!(edges.get(&(b, a)), Some(&1));
            }
        }
        assert_eq!(icosphere(1.0).0.len(), 42);
    }

    #[test]
    fn torus_is_a_closed_outward_mesh_with_the_expected_volume() {
        let (major, minor) = (3.0, 1.0);
        let mesh = torus(major, minor, 24, 12);
        assert_eq!(mesh.0.len(), 24 * 12);
        assert_eq!(mesh.1.len(), 2 * 24 * 12);
        let mut edges = std::collections::HashMap::new();
        for t in &mesh.1 {
            for k in 0..3 {
                *edges.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
            }
        }
        for (&(a, b), &count) in &edges {
            assert_eq!(count, 1);
            assert_eq!(edges.get(&(b, a)), Some(&1));
        }
        // The smooth torus has volume 2 pi^2 R r^2; the inscribed polyhedron is
        // a little smaller.
        let smooth = 2.0 * std::f64::consts::PI.powi(2) * major * minor * minor;
        let volume = signed_volume(&mesh);
        assert!(volume > 0.9 * smooth && volume < smooth, "volume {volume}");
    }

    #[test]
    fn extruded_l_shape_has_the_expected_volume() {
        let l_shape = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let mesh = extrude_polygon(&l_shape, 1.0);
        assert_eq!(mesh.0.len(), 12);
        assert!((signed_volume(&mesh) - 3.0).abs() < 1e-12);
        assert!((signed_volume(&flipped(&mesh)) + 3.0).abs() < 1e-12);
    }
}
