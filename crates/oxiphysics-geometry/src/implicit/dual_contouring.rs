// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dual contouring of a signed distance grid (Ju, Losasso, Schaefer and Warren,
//! *Dual Contouring of Hermite Data*, SIGGRAPH 2002).
//!
//! The grid samples ([`SdfGrid3D`]) are the corners of the cells.  A cell is
//! *active* when the surface crosses at least one of its twelve edges, and every
//! active cell receives one mesh vertex: the minimiser of the *quadratic error
//! function* (QEF) of the Hermite data of its edges, i.e. the point whose summed
//! squared distance to the tangent planes at the edge crossings is smallest.  Every
//! grid edge that the surface crosses then yields one quad joining the vertices of
//! the four cells around the edge.

use super::{
    IsoMesh, SdfGrid3D, surface_normal, vec3_add, vec3_dot, vec3_norm, vec3_scale, vec3_sub,
};
use crate::mc_tables::{MC_CORNER_OFFSETS, MC_EDGE_ASCENDING, case_index};
use oxiphysics_core::math::Mat3;

/// Marks a cell that has no mesh vertex.
const NO_VERTEX: usize = usize::MAX;

/// A Hermite sample: a point of the surface and the unit normal there.
type Hermite = ([f64; 3], [f64; 3]);

/// Eigenvalues of the QEF matrix below this fraction of the largest one are
/// treated as zero (singular values below a tenth of the largest, as in the
/// original paper): along such a direction the Hermite data do not constrain the
/// vertex, which then stays at the mass point of the crossings.
const QEF_TRUNCATION: f64 = 1.0e-2;

/// Dual contouring with QEF minimisation for sharp-feature preservation.
pub struct DualContouring;

impl DualContouring {
    /// Extract the isosurface `value = isovalue` of an SDF grid with dual
    /// contouring.
    ///
    /// A grid node is *inside* when its value is strictly below `isovalue`.  The
    /// mesh has one vertex per active cell, placed at the QEF minimiser of the
    /// crossings of the cell's edges (linear interpolation along the edge) and their
    /// normals (the normalised grid gradient there), regularised towards the mass
    /// point of the crossings and clamped to the cell.  Each grid edge that
    /// changes sign contributes a quad of the four surrounding cells, split along
    /// its shorter diagonal into two triangles that are wound counter-clockwise
    /// when seen from outside, so their normals point towards increasing field
    /// values.  The per-vertex normals are the normalised grid gradient at the
    /// vertex (the zero vector where the gradient is not finite).
    ///
    /// Edges on the boundary of the grid do not have four surrounding cells and
    /// produce no quad, so a surface that leaves the grid ends in an open boundary.
    /// Cells with a non-finite corner value are ignored; a grid with fewer than two
    /// nodes along an axis, a non-positive or non-finite voxel size, a non-finite
    /// origin or `isovalue`, or a data array of the wrong length gives an empty
    /// mesh.
    pub fn extract(grid: &SdfGrid3D, isovalue: f64) -> IsoMesh {
        let mut mesh = IsoMesh::new();
        if !grid.is_samplable() || !isovalue.is_finite() {
            return mesh;
        }
        let dims = [grid.nx, grid.ny, grid.nz];
        let cells = [dims[0] - 1, dims[1] - 1, dims[2] - 1];
        let cell_id = |c: [usize; 3]| c[0] + cells[0] * (c[1] + cells[1] * c[2]);
        let mut vertex_of_cell = vec![NO_VERTEX; cells[0] * cells[1] * cells[2]];

        // One vertex per active cell.
        for cz in 0..cells[2] {
            for cy in 0..cells[1] {
                for cx in 0..cells[0] {
                    let cell = [cx, cy, cz];
                    let values = cell_values(grid, cell);
                    if !is_active(&values, isovalue) {
                        continue;
                    }
                    let positions = cell_positions(grid, cell);
                    let vertex = cell_vertex(grid, &values, &positions, isovalue);
                    vertex_of_cell[cell_id(cell)] = mesh.vertices.len();
                    mesh.vertices.push(vertex);
                    mesh.normals.push(surface_normal(grid, vertex));
                }
            }
        }

        // One quad per sign-changing grid edge that has four surrounding cells.
        for axis in 0..3 {
            let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
            for nc in 1..dims[c] - 1 {
                for nb in 1..dims[b] - 1 {
                    for na in 0..dims[axis] - 1 {
                        let mut lower = [0usize; 3];
                        lower[axis] = na;
                        lower[b] = nb;
                        lower[c] = nc;
                        let mut upper = lower;
                        upper[axis] += 1;
                        let value_lower = grid.get(lower[0], lower[1], lower[2]);
                        let value_upper = grid.get(upper[0], upper[1], upper[2]);
                        let lower_inside = value_lower < isovalue;
                        if lower_inside == (value_upper < isovalue) {
                            continue;
                        }
                        // The four cells around the edge, counter-clockwise when
                        // seen from +axis (b to the right, c up): the cell offsets
                        // in (b, c) are (-1, -1), (0, -1), (0, 0), (-1, 0).
                        let mut quad = [NO_VERTEX; 4];
                        for (slot, (back_b, back_c)) in
                            [(1, 1), (0, 1), (0, 0), (1, 0)].into_iter().enumerate()
                        {
                            let mut cell = lower;
                            cell[b] -= back_b;
                            cell[c] -= back_c;
                            quad[slot] = vertex_of_cell[cell_id(cell)];
                        }
                        if quad.contains(&NO_VERTEX) {
                            continue;
                        }
                        // The field increases along +axis when the lower node is
                        // inside, and the outward normal then is +axis.
                        if !lower_inside {
                            quad.swap(1, 3);
                        }
                        push_quad(&mut mesh, quad);
                    }
                }
            }
        }
        mesh
    }
}

/// Field values at the eight corners of a cell, in the order of
/// [`MC_CORNER_OFFSETS`].
fn cell_values(grid: &SdfGrid3D, cell: [usize; 3]) -> [f64; 8] {
    std::array::from_fn(|corner| {
        let [di, dj, dk] = MC_CORNER_OFFSETS[corner];
        grid.get(cell[0] + di, cell[1] + dj, cell[2] + dk)
    })
}

/// Node positions of the eight corners of a cell, in the order of
/// [`MC_CORNER_OFFSETS`].
fn cell_positions(grid: &SdfGrid3D, cell: [usize; 3]) -> [[f64; 3]; 8] {
    std::array::from_fn(|corner| {
        let [di, dj, dk] = MC_CORNER_OFFSETS[corner];
        grid.voxel_center(cell[0] + di, cell[1] + dj, cell[2] + dk)
    })
}

/// `true` when the surface crosses at least one of the twelve edges of the cell
/// (its corners are not all on the same side) and all corner values are finite.
fn is_active(values: &[f64; 8], isovalue: f64) -> bool {
    if values.iter().any(|v| !v.is_finite()) {
        return false;
    }
    let case = case_index(values, isovalue);
    case != 0 && case != 255
}

/// Hermite data of the crossed edges of a cell: the point where the surface
/// crosses the edge (linear interpolation of the field along it) and the unit
/// surface normal there (the normalised grid gradient).
fn hermite_samples(
    grid: &SdfGrid3D,
    values: &[f64; 8],
    positions: &[[f64; 3]; 8],
    isovalue: f64,
) -> Vec<Hermite> {
    let mut samples: Vec<Hermite> = Vec::with_capacity(12);
    for &[lo, hi] in &MC_EDGE_ASCENDING {
        if (values[lo] < isovalue) == (values[hi] < isovalue) {
            continue;
        }
        let t = ((isovalue - values[lo]) / (values[hi] - values[lo])).clamp(0.0, 1.0);
        let span = vec3_sub(positions[hi], positions[lo]);
        let point = vec3_add(positions[lo], vec3_scale(span, t));
        samples.push((point, vec3_norm(grid.gradient(point))));
    }
    samples
}

/// Mesh vertex of an active cell: the QEF minimiser of the Hermite data of its
/// crossed edges, clamped to the cell.
///
/// `values` and `positions` are the field values and node positions of the eight
/// cell corners in the order of [`MC_CORNER_OFFSETS`].
fn cell_vertex(
    grid: &SdfGrid3D,
    values: &[f64; 8],
    positions: &[[f64; 3]; 8],
    isovalue: f64,
) -> [f64; 3] {
    let minimiser = qef_minimiser(&hermite_samples(grid, values, positions, isovalue));
    // Corner 0 is the lower and corner 6 the upper corner of the cell.
    let (low, high) = (positions[0], positions[6]);
    std::array::from_fn(|d| minimiser[d].max(low[d]).min(high[d]))
}

/// Mean of the sample points.
fn mass_point(samples: &[Hermite]) -> [f64; 3] {
    let sum = samples
        .iter()
        .fold([0.0; 3], |acc, (point, _)| vec3_add(acc, *point));
    vec3_scale(sum, 1.0 / samples.len().max(1) as f64)
}

/// Minimiser of `sum_i (n_i . (x - p_i))^2` over the Hermite samples `(p_i, n_i)`,
/// regularised towards the mass point of the samples.
///
/// With the mass point `c` as the origin the minimiser solves the normal equations
/// `(A^T A) d = A^T b` with the rows `n_i` of `A` and `b_i = n_i . (p_i - c)`.  The
/// symmetric matrix `A^T A` is decomposed into its eigenvectors and inverted only
/// along the eigenvectors with an eigenvalue above [`QEF_TRUNCATION`] times the
/// largest; the other directions keep the mass point.  A smooth surface (nearly
/// parallel normals) therefore only moves the vertex along the common normal, a
/// crease fixes it on the crease line, and a corner pins it down completely.
fn qef_minimiser(samples: &[Hermite]) -> [f64; 3] {
    if samples.is_empty() {
        return [0.0; 3];
    }
    let mass = mass_point(samples);
    let mut ata = [[0.0_f64; 3]; 3];
    let mut atb = [0.0_f64; 3];
    for (point, normal) in samples {
        let b = vec3_dot(*normal, vec3_sub(*point, mass));
        for (row, atb_row) in atb.iter_mut().enumerate() {
            for (entry, &n_col) in ata[row].iter_mut().zip(normal.iter()) {
                *entry += normal[row] * n_col;
            }
            *atb_row += normal[row] * b;
        }
    }
    let finite = ata
        .iter()
        .flatten()
        .chain(atb.iter())
        .all(|v| v.is_finite());
    if !finite {
        return mass;
    }
    let eigen = Mat3::from_fn(|r, c| ata[r][c]).symmetric_eigen();
    let largest = eigen.eigenvalues.iter().copied().fold(0.0_f64, f64::max);
    if largest <= 0.0 {
        return mass;
    }
    let mut offset = [0.0_f64; 3];
    for (k, &eigenvalue) in eigen.eigenvalues.iter().enumerate() {
        if eigenvalue <= QEF_TRUNCATION * largest {
            continue;
        }
        let axis = eigen.eigenvectors.column(k);
        let axis = [axis[0], axis[1], axis[2]];
        let coefficient = vec3_dot(axis, atb) / eigenvalue;
        offset = vec3_add(offset, vec3_scale(axis, coefficient));
    }
    vec3_add(mass, offset)
}

/// Append the quad `[v0, v1, v2, v3]` (counter-clockwise seen from outside) as two
/// triangles split along the shorter diagonal.
fn push_quad(mesh: &mut IsoMesh, quad: [usize; 4]) {
    let p = quad.map(|v| mesh.vertices[v]);
    let squared_length = |a: [f64; 3], b: [f64; 3]| {
        let d = vec3_sub(a, b);
        vec3_dot(d, d)
    };
    if squared_length(p[0], p[2]) <= squared_length(p[1], p[3]) {
        mesh.triangles.push([quad[0], quad[1], quad[2]]);
        mesh.triangles.push([quad[0], quad[2], quad[3]]);
    } else {
        mesh.triangles.push([quad[1], quad[2], quad[3]]);
        mesh.triangles.push([quad[1], quad[3], quad[0]]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::implicit::{ImplicitBox, ImplicitSphere, ImplicitTorus, vec3_cross};
    use crate::mc_tables::test_support::{
        area_and_volume, degenerate_triangle_count, edge_valence_histogram,
        is_closed_and_consistently_oriented,
    };
    use std::f64::consts::PI;

    /// Centre of the sphere sampled by [`sphere_grid`]: off the lattice of grid
    /// nodes, so that no node lies exactly on the surface.
    fn sphere_center(n: usize, voxel: f64) -> [f64; 3] {
        let half = 0.5 * n as f64 * voxel;
        [half + 0.013, half - 0.007, half + 0.0031]
    }

    fn sphere_grid(n: usize, voxel: f64, radius: f64) -> SdfGrid3D {
        let s = ImplicitSphere::new(sphere_center(n, voxel), radius);
        SdfGrid3D::from_sdf(n, n, n, [0.0; 3], voxel, |p| s.sdf(p))
    }

    #[test]
    fn sphere_mesh_is_watertight_outward_and_accurate() {
        let (n, voxel, radius) = (26, 0.12, 1.0);
        let grid = sphere_grid(n, voxel, radius);
        let mesh = DualContouring::extract(&grid, 0.0);
        assert!(mesh.num_triangles() > 1000, "{}", mesh.num_triangles());
        assert_eq!(mesh.normals.len(), mesh.vertices.len());
        assert_eq!(
            edge_valence_histogram(&mesh.triangles),
            vec![(2, mesh.triangles.len() * 3 / 2)],
            "every edge must be shared by exactly two triangles"
        );
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        assert_eq!(
            degenerate_triangle_count(&mesh.vertices, &mesh.triangles),
            0
        );
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        let exact_area = 4.0 * PI * radius * radius;
        let exact_volume = 4.0 / 3.0 * PI * radius.powi(3);
        assert!(
            (area / exact_area - 1.0).abs() < 0.03,
            "area {area} vs {exact_area}"
        );
        assert!(volume > 0.0, "triangles must face outwards: {volume}");
        assert!(
            (volume / exact_volume - 1.0).abs() < 0.03,
            "volume {volume} vs {exact_volume}"
        );
    }

    #[test]
    fn torus_mesh_is_a_closed_surface_of_genus_one() {
        let (n, voxel) = (56, 0.07);
        let half = 0.5 * n as f64 * voxel;
        let centre = [half + 0.013, half - 0.007, half + 0.0031];
        let torus = ImplicitTorus::new(1.0, 0.4);
        let grid =
            SdfGrid3D::from_sdf(n, n, n, [0.0; 3], voxel, |p| torus.sdf(vec3_sub(p, centre)));
        let mesh = DualContouring::extract(&grid, 0.0);
        let faces = mesh.num_triangles();
        assert!(faces > 5000, "{faces}");
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        assert_eq!(
            edge_valence_histogram(&mesh.triangles),
            vec![(2, faces * 3 / 2)],
            "every edge must be shared by exactly two triangles"
        );
        // With E = 3F/2 the Euler characteristic V - E + F is V - F/2; a torus has 0.
        let mut used = vec![false; mesh.num_vertices()];
        for &vertex in mesh.triangles.iter().flatten() {
            used[vertex] = true;
        }
        let vertices = used.iter().filter(|&&u| u).count();
        assert_eq!(2 * vertices, faces, "Euler characteristic is not zero");
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        let exact_area = 4.0 * PI * PI * 1.0 * 0.4;
        let exact_volume = 2.0 * PI * PI * 1.0 * 0.4 * 0.4;
        assert!((area / exact_area - 1.0).abs() < 0.03, "area {area}");
        assert!(
            (volume / exact_volume - 1.0).abs() < 0.03,
            "volume {volume}"
        );
    }

    #[test]
    fn the_isovalue_selects_the_level_set() {
        // The level set `sdf = 0.25` of a unit sphere is the sphere of radius 1.25.
        let (n, voxel) = (26, 0.13);
        let grid = sphere_grid(n, voxel, 1.0);
        let centre = sphere_center(n, voxel);
        let mesh = DualContouring::extract(&grid, 0.25);
        assert!(mesh.num_triangles() > 1000, "{}", mesh.num_triangles());
        for v in &mesh.vertices {
            let r = vec3_len_from(*v, centre);
            assert!((r - 1.25).abs() < 0.05 * voxel, "vertex at radius {r}");
        }
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        let exact_area = 4.0 * PI * 1.25 * 1.25;
        let exact_volume = 4.0 / 3.0 * PI * 1.25_f64.powi(3);
        assert!((area / exact_area - 1.0).abs() < 0.03, "area {area}");
        assert!(
            (volume / exact_volume - 1.0).abs() < 0.03,
            "volume {volume}"
        );
    }

    #[test]
    fn vertices_sit_on_the_sphere_not_at_the_cell_corners() {
        // The QEF moves every vertex onto the tangent planes of its cell; a
        // placement at a grid node would be off by up to a voxel.
        let (n, voxel, radius) = (26, 0.12, 1.0);
        let grid = sphere_grid(n, voxel, radius);
        let centre = sphere_center(n, voxel);
        let mesh = DualContouring::extract(&grid, 0.0);
        let mut worst = 0.0_f64;
        for v in &mesh.vertices {
            let r = vec3_len_from(*v, centre);
            worst = worst.max((r - radius).abs());
        }
        assert!(worst < 0.05 * voxel, "worst radial error {worst}");
    }

    fn vec3_len_from(a: [f64; 3], b: [f64; 3]) -> f64 {
        let d = vec3_sub(a, b);
        vec3_dot(d, d).sqrt()
    }

    #[test]
    fn normals_follow_the_gradient_and_the_winding_agrees() {
        let (n, voxel) = (22, 0.14);
        let grid = sphere_grid(n, voxel, 1.0);
        let centre = sphere_center(n, voxel);
        let mesh = DualContouring::extract(&grid, 0.0);
        for (p, nrm) in mesh.vertices.iter().zip(mesh.normals.iter()) {
            assert!((vec3_dot(*nrm, *nrm).sqrt() - 1.0).abs() < 1e-9);
            let radial = vec3_norm(vec3_sub(*p, centre));
            assert!(vec3_dot(radial, *nrm) > 0.98, "normal {nrm:?} at {p:?}");
        }
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let face = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
            let avg = vec3_add(
                vec3_add(mesh.normals[t[0]], mesh.normals[t[1]]),
                mesh.normals[t[2]],
            );
            assert!(vec3_dot(face, avg) > 0.0, "triangle {t:?} faces inwards");
        }
    }

    /// The box used by the sharp-feature tests: its faces are not aligned with the
    /// grid, so every corner and edge cuts through the cells.
    fn off_lattice_box() -> (ImplicitBox, [f64; 3], [f64; 3], SdfGrid3D) {
        let centre = [1.5237, 1.4871, 1.5413];
        let half = [0.8312, 0.6944, 0.5711];
        let b = ImplicitBox::new(centre, half);
        let n = 32;
        let grid = SdfGrid3D::from_sdf(n, n, n, [0.0; 3], 0.1, |p| b.sdf(p));
        (b, centre, half, grid)
    }

    /// Distance from every corner of the box to the nearest of `vertices`.
    fn corner_errors(vertices: &[[f64; 3]], centre: [f64; 3], half: [f64; 3]) -> Vec<f64> {
        let mut errors = Vec::new();
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                for sz in [-1.0, 1.0] {
                    let corner = [
                        centre[0] + sx * half[0],
                        centre[1] + sy * half[1],
                        centre[2] + sz * half[2],
                    ];
                    let nearest = vertices
                        .iter()
                        .map(|&p| vec3_len_from(p, corner))
                        .fold(f64::INFINITY, f64::min);
                    errors.push(nearest);
                }
            }
        }
        errors
    }

    #[test]
    fn sharp_features_are_better_preserved_than_by_the_mass_point() {
        // The vertices of the cells at the corners and edges of a box are pulled onto
        // the creases by the QEF.  The Hermite normals come from the gradient of the
        // sampled distance field, which near a crease is a blend of the adjacent face
        // normals, so the corners are recovered to about half a voxel -- but a
        // vertex at the mass point of the crossings is almost a whole voxel off.
        let (b, centre, half, grid) = off_lattice_box();
        let voxel = grid.voxel_size;
        let mut qef = Vec::new();
        let mut mass = Vec::new();
        for cz in 0..grid.nz - 1 {
            for cy in 0..grid.ny - 1 {
                for cx in 0..grid.nx - 1 {
                    let values = cell_values(&grid, [cx, cy, cz]);
                    if !is_active(&values, 0.0) {
                        continue;
                    }
                    let positions = cell_positions(&grid, [cx, cy, cz]);
                    qef.push(cell_vertex(&grid, &values, &positions, 0.0));
                    mass.push(mass_point(&hermite_samples(
                        &grid, &values, &positions, 0.0,
                    )));
                }
            }
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let qef_corners = corner_errors(&qef, centre, half);
        let mass_corners = corner_errors(&mass, centre, half);
        assert!(
            mean(&qef_corners) < 0.6 * voxel,
            "mean corner error {} voxels",
            mean(&qef_corners) / voxel
        );
        assert!(
            qef_corners.iter().all(|&e| e < 0.7 * voxel),
            "{qef_corners:?}"
        );
        assert!(
            mean(&qef_corners) < 0.75 * mean(&mass_corners),
            "QEF {} vs mass point {}",
            mean(&qef_corners),
            mean(&mass_corners)
        );
        // The same on the creases and faces: the vertices are closer to the surface.
        let worst = |vertices: &[[f64; 3]]| {
            vertices
                .iter()
                .map(|&p| b.sdf(p).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(worst(&qef) < 0.3 * voxel, "QEF {}", worst(&qef));
        assert!(worst(&qef) < worst(&mass));
    }

    #[test]
    fn box_mesh_is_closed_outward_and_has_the_right_size() {
        let (_, _, half, grid) = off_lattice_box();
        let mesh = DualContouring::extract(&grid, 0.0);
        assert!(mesh.num_triangles() > 200);
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        assert_eq!(
            degenerate_triangle_count(&mesh.vertices, &mesh.triangles),
            0
        );
        // The creases are chamfered slightly (the Hermite normals of a sampled
        // distance field are blurred next to a crease), so the box comes out a few
        // percent too small.
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        let exact_area = 8.0 * (half[0] * half[1] + half[0] * half[2] + half[1] * half[2]);
        let exact_volume = 8.0 * half[0] * half[1] * half[2];
        assert!(
            (area / exact_area - 1.0).abs() < 0.06,
            "area {area} vs {exact_area}"
        );
        assert!(
            (volume / exact_volume - 1.0).abs() < 0.04,
            "volume {volume} vs {exact_volume}"
        );
    }

    #[test]
    fn a_single_inside_node_activates_all_eight_surrounding_cells() {
        // Only the corner (1, 1, 1) of the cells differs in sign.  For the cell
        // (0, 0, 0) that is the far corner, which a test of the three edges at the
        // cell origin cannot see.
        let mut grid = SdfGrid3D::from_sdf(3, 3, 3, [0.0; 3], 1.0, |_| 1.0);
        grid.set(1, 1, 1, -1.0);
        let mesh = DualContouring::extract(&grid, 0.0);
        assert_eq!(mesh.num_vertices(), 8);
        // The six edges at the node give six quads.
        assert_eq!(mesh.num_triangles(), 12);
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        assert_eq!(
            edge_valence_histogram(&mesh.triangles),
            vec![(2, 18)],
            "a closed cube of 8 vertices: 12 cube edges + 6 diagonals"
        );
        // The crossings are the midpoints of the six edges at the node (1.5, 1.5, 1.5)
        // with axis-aligned normals, so the QEF vertices are the corners of the unit
        // cube [1, 2]^3 and the dual surface is that cube.
        for v in &mesh.vertices {
            assert!(
                v.iter()
                    .all(|&c| (c - 1.0).abs() < 1e-9 || (c - 2.0).abs() < 1e-9),
                "{v:?} is not a corner of [1, 2]^3"
            );
        }
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        assert!((area - 6.0).abs() < 1e-9, "area {area}");
        assert!((volume - 1.0).abs() < 1e-9, "volume {volume}");
    }

    #[test]
    fn plane_between_node_layers_gives_a_flat_quad_grid() {
        // x = 5.3 with unit voxels (node x are half-integers): one layer of quads.
        let (nx, ny, nz) = (9, 7, 6);
        let grid = SdfGrid3D::from_sdf(nx, ny, nz, [0.0; 3], 1.0, |p| p[0] - 5.3);
        let mesh = DualContouring::extract(&grid, 0.0);
        assert_eq!(mesh.num_vertices(), (ny - 1) * (nz - 1));
        assert_eq!(mesh.num_triangles(), 2 * (ny - 2) * (nz - 2));
        for v in &mesh.vertices {
            assert!((v[0] - 5.3).abs() < 1e-12, "{v:?}");
        }
        let (area, _) = area_and_volume(&mesh.vertices, &mesh.triangles);
        assert!(
            (area - ((ny - 2) * (nz - 2)) as f64).abs() < 1e-9,
            "area {area}"
        );
        // Quads face +x, the direction of increasing field.
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let face = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
            assert!(face[0] > 0.0, "triangle {t:?} faces {face:?}");
        }
    }

    #[test]
    fn tilted_plane_is_reproduced_by_every_vertex() {
        let normal = vec3_norm([1.0, 2.0, -2.0]);
        let offset = 0.43;
        let n = 14;
        let voxel = 0.5;
        let grid = SdfGrid3D::from_sdf(n, n, n, [0.0; 3], voxel, |p| {
            vec3_dot(normal, vec3_sub(p, [3.1, 3.2, 3.3])) - offset
        });
        let mesh = DualContouring::extract(&grid, 0.0);
        assert!(mesh.num_triangles() > 100);
        let mut checked = 0;
        for v in &mesh.vertices {
            // Away from the grid boundary, where the gradient stencil is clamped.
            if v.iter()
                .all(|&c| c > 1.5 * voxel && c < (n as f64 - 1.5) * voxel)
            {
                let d = vec3_dot(normal, vec3_sub(*v, [3.1, 3.2, 3.3])) - offset;
                assert!(d.abs() < 1e-9, "vertex {v:?} is {d} off the plane");
                checked += 1;
            }
        }
        assert!(checked > 20, "only {checked} vertices checked");
    }

    #[test]
    fn nonfinite_and_degenerate_grids_give_empty_meshes() {
        let flat = SdfGrid3D::from_sdf(1, 5, 5, [0.0; 3], 1.0, |p| p[1] - 2.0);
        assert_eq!(DualContouring::extract(&flat, 0.0).num_triangles(), 0);
        let constant = SdfGrid3D::from_sdf(5, 5, 5, [0.0; 3], 1.0, |_| 1.0);
        assert_eq!(DualContouring::extract(&constant, 0.0).num_vertices(), 0);
        let mut bad_size = SdfGrid3D::from_sdf(5, 5, 5, [0.0; 3], 1.0, |p| p[0] - 2.5);
        for size in [0.0, -1.0, f64::NAN] {
            bad_size.voxel_size = size;
            assert_eq!(DualContouring::extract(&bad_size, 0.0).num_vertices(), 0);
        }
        let grid = SdfGrid3D::from_sdf(5, 5, 5, [0.0; 3], 1.0, |p| p[0] - 2.5);
        assert_eq!(
            DualContouring::extract(&grid, f64::NAN).num_vertices(),
            0,
            "NaN isovalue"
        );
        // A NaN sample only removes the cells it touches.
        let mut holey = SdfGrid3D::from_sdf(7, 7, 7, [0.0; 3], 1.0, |p| p[0] - 3.4);
        holey.set(3, 3, 3, f64::NAN);
        let mesh = DualContouring::extract(&holey, 0.0);
        let full = DualContouring::extract(
            &SdfGrid3D::from_sdf(7, 7, 7, [0.0; 3], 1.0, |p| p[0] - 3.4),
            0.0,
        );
        assert!(mesh.num_vertices() < full.num_vertices());
        assert!(mesh.num_vertices() > 0);
        assert!(
            mesh.vertices
                .iter()
                .all(|v| v.iter().all(|c| c.is_finite()))
        );
        assert!(
            mesh.normals.iter().all(|n| n.iter().all(|c| c.is_finite())),
            "a stencil that reaches the NaN sample must not leak NaN normals"
        );
    }

    #[test]
    fn qef_minimiser_recovers_a_corner_and_ignores_parallel_planes() {
        // Three orthogonal planes through (0.3, -0.2, 0.7).
        let corner = [0.3, -0.2, 0.7];
        let samples = [
            ([0.3, 0.1, 0.5], [1.0, 0.0, 0.0]),
            ([0.0, -0.2, 0.9], [0.0, 1.0, 0.0]),
            ([0.5, 0.4, 0.7], [0.0, 0.0, 1.0]),
        ];
        let x = qef_minimiser(&samples);
        for d in 0..3 {
            assert!((x[d] - corner[d]).abs() < 1e-12, "{x:?}");
        }
        // Three parallel planes' worth of samples: only the normal direction is
        // constrained, the vertex stays at the mass point in the plane.
        let parallel = [
            ([0.0, 0.0, 0.2], [0.0, 0.0, 1.0]),
            ([1.0, 0.0, 0.2], [0.0, 0.0, 1.0]),
            ([0.0, 1.0, 0.2], [0.0, 0.0, 1.0]),
        ];
        let y = qef_minimiser(&parallel);
        assert!((y[0] - 1.0 / 3.0).abs() < 1e-12 && (y[1] - 1.0 / 3.0).abs() < 1e-12);
        assert!((y[2] - 0.2).abs() < 1e-12);
        assert_eq!(qef_minimiser(&[]), [0.0; 3]);
    }
}
