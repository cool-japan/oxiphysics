// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Small reference meshes shared by the unit tests of several modules.

use std::collections::HashMap;

/// Vertices and counter-clockwise triangles of a mesh.
pub(crate) type RawMesh = (Vec<[f64; 3]>, Vec<[usize; 3]>);

/// Unit icosphere with outward, counter-clockwise triangles:
/// `20 * 4^subdivisions` triangles and `10 * 4^subdivisions + 2` vertices.
pub(crate) fn icosphere(subdivisions: usize) -> RawMesh {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let unit = |p: [f64; 3]| {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        [p[0] / l, p[1] / l, p[2] / l]
    };
    let mut vertices: Vec<[f64; 3]> = [
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ]
    .iter()
    .map(|&p| unit(p))
    .collect();
    let mut faces: Vec<[usize; 3]> = vec![
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
    for _ in 0..subdivisions {
        let mut midpoint: HashMap<(usize, usize), usize> = HashMap::new();
        let mut next = Vec::with_capacity(faces.len() * 4);
        for &[a, b, c] in &faces {
            let mut mid = |p: usize, q: usize| -> usize {
                let key = (p.min(q), p.max(q));
                *midpoint.entry(key).or_insert_with(|| {
                    let (vp, vq) = (vertices[p], vertices[q]);
                    vertices.push(unit([vp[0] + vq[0], vp[1] + vq[1], vp[2] + vq[2]]));
                    vertices.len() - 1
                })
            };
            let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
            next.push([a, ab, ca]);
            next.push([b, bc, ab]);
            next.push([c, ca, bc]);
            next.push([ab, bc, ca]);
        }
        faces = next;
    }
    (vertices, faces)
}

/// A flat `nx x ny` grid of quads with the given spacing in the plane z = 0
/// (vertex `(i, j)` has index `j * (nx + 1) + i`), two counter-clockwise
/// triangles per quad.
pub(crate) fn flat_grid(nx: usize, ny: usize, spacing: f64) -> RawMesh {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([i as f64 * spacing, j as f64 * spacing, 0.0]);
        }
    }
    let at = |i: usize, j: usize| j * (nx + 1) + i;
    let mut faces = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            faces.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            faces.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    (vertices, faces)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn icosphere_has_the_expected_size_and_orientation() {
        for s in 0..4 {
            let (v, f) = icosphere(s);
            assert_eq!(f.len(), 20 * 4usize.pow(s as u32));
            assert_eq!(v.len(), 10 * 4usize.pow(s as u32) + 2);
            let volume: f64 = f
                .iter()
                .map(|t| {
                    let (a, b, c) = (v[t[0]], v[t[1]], v[t[2]]);
                    (a[0] * (b[1] * c[2] - b[2] * c[1])
                        + a[1] * (b[2] * c[0] - b[0] * c[2])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]))
                        / 6.0
                })
                .sum();
            assert!(volume > 0.0, "triangles must face outwards: {volume}");
            for p in &v {
                let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                assert!((r - 1.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn flat_grid_counts() {
        let (v, f) = flat_grid(4, 3, 0.5);
        assert_eq!(v.len(), 5 * 4);
        assert_eq!(f.len(), 2 * 4 * 3);
        assert_eq!(v[7], [1.0, 0.5, 0.0]);
    }
}
