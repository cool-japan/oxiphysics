//! Simplicial homology of a [`SimplicialComplex`] over the field `Z/2Z`.
//!
//! The Betti number `beta_k` is the dimension of the `k`-th homology group,
//! `beta_k = dim C_k - rank(d_k) - rank(d_{k+1})`, where `C_k` is the vector
//! space of `k`-chains and `d_k : C_k -> C_{k-1}` is the boundary map. With
//! coefficients in `Z/2Z` the boundary of a simplex is simply the sum of its
//! facets (no orientation signs), every matrix entry is a bit, and the ranks
//! follow from Gaussian elimination in which row addition is XOR.
//!
//! Over `Z/2Z` the Betti numbers of every orientable surface agree with the
//! rational ones. Complexes with 2-torsion differ: the projective plane has
//! `(beta_0, beta_1, beta_2) = (1, 1, 1)` over `Z/2Z` but `(1, 0, 0)` over `Q`.

use std::collections::BTreeSet;

use super::types::SimplicialComplex;

/// Betti numbers `(beta_0, beta_1, beta_2)` of a simplicial complex over `Z/2Z`.
///
/// `beta_0` is the number of connected components, `beta_1` the number of
/// independent loops and `beta_2` the number of independent enclosed voids.
/// For example a filled triangle is `(1, 0, 0)`, two disjoint filled triangles
/// are `(2, 0, 0)`, the boundary of a tetrahedron (a sphere) is `(1, 0, 1)` and
/// a torus is `(1, 2, 1)`.
///
/// The complex is closed under faces before the computation, so a complex that
/// was filled in through the public `simplices` field with missing faces is
/// handled as the complex those simplices generate. Dimensions above 2 take
/// part in the elimination (a solid tetrahedron has `beta_2 = 0`) but are not
/// reported; use [`betti_numbers_by_dimension`] for the full list.
pub fn betti_numbers(sc: &SimplicialComplex) -> (usize, usize, usize) {
    let all = betti_numbers_by_dimension(sc);
    let get = |k: usize| all.get(k).copied().unwrap_or(0);
    (get(0), get(1), get(2))
}

/// Betti numbers `beta_0 ..= beta_d` over `Z/2Z`, where `d` is the highest
/// dimension that contains a simplex. An empty complex yields an empty list.
///
/// See [`betti_numbers`] for the conventions. The result satisfies the
/// Euler-Poincare identity `sum_k (-1)^k beta_k = sum_k (-1)^k |C_k|`.
pub fn betti_numbers_by_dimension(sc: &SimplicialComplex) -> Vec<usize> {
    let levels = face_closed_levels(sc);
    let top = levels.len();
    if top == 0 {
        return Vec::new();
    }
    // rank[k] is the rank of d_k : C_k -> C_{k-1}; d_0 and d_{top} (zero) are 0.
    let mut rank = vec![0usize; top + 1];
    for k in 1..top {
        rank[k] = boundary_rank(&levels[k], &levels[k - 1]);
    }
    // `d_k . d_{k+1} = 0` over Z/2Z guarantees `rank[k] + rank[k + 1] <= |C_k|`.
    (0..top)
        .map(|k| levels[k].len().saturating_sub(rank[k] + rank[k + 1]))
        .collect()
}

/// Deprecated name of [`betti_numbers`].
///
/// The old implementation guessed the numbers from the Euler characteristic
/// and was wrong even for a single filled triangle; this alias now computes
/// the exact homology.
#[deprecated(
    since = "0.1.3",
    note = "use `betti_numbers`, which computes the exact homology"
)]
pub fn betti_numbers_approx(sc: &SimplicialComplex) -> (usize, usize, usize) {
    betti_numbers(sc)
}

/// Gather the simplices of `sc` per dimension (sorted, deduplicated vertex
/// tuples), closing the collection under taking faces.
///
/// The simplices are ordered lexicographically so that indices (and therefore
/// the elimination) are deterministic regardless of `HashSet` iteration order.
/// The returned list has one entry per dimension up to the highest non-empty
/// one.
fn face_closed_levels(sc: &SimplicialComplex) -> Vec<Vec<Vec<usize>>> {
    let mut levels: Vec<BTreeSet<Vec<usize>>> = Vec::new();
    for simplex in sc.simplices.iter().flatten() {
        let mut vertices = simplex.vertices.clone();
        vertices.sort_unstable();
        vertices.dedup();
        if vertices.is_empty() {
            continue;
        }
        let dim = vertices.len() - 1;
        if levels.len() <= dim {
            levels.resize_with(dim + 1, BTreeSet::new);
        }
        levels[dim].insert(vertices);
    }
    for dim in (1..levels.len()).rev() {
        let (lower, upper) = levels.split_at_mut(dim);
        for simplex in &upper[0] {
            for skip in 0..simplex.len() {
                let face: Vec<usize> = simplex
                    .iter()
                    .enumerate()
                    .filter(|&(i, _)| i != skip)
                    .map(|(_, &v)| v)
                    .collect();
                lower[dim - 1].insert(face);
            }
        }
    }
    levels
        .into_iter()
        .map(|set| set.into_iter().collect())
        .collect()
}

/// Rank over `Z/2Z` of the boundary map from `simplices` (dimension `k`) to
/// `faces` (dimension `k - 1`).
///
/// Both slices are sorted lexicographically; the column of a simplex lists the
/// indices of its facets in ascending order.
fn boundary_rank(simplices: &[Vec<usize>], faces: &[Vec<usize>]) -> usize {
    let columns: Vec<Vec<usize>> = simplices
        .iter()
        .map(|simplex| {
            let mut column: Vec<usize> = (0..simplex.len())
                .filter_map(|skip| {
                    let facet: Vec<usize> = simplex
                        .iter()
                        .enumerate()
                        .filter(|&(i, _)| i != skip)
                        .map(|(_, &v)| v)
                        .collect();
                    faces.binary_search(&facet).ok()
                })
                .collect();
            column.sort_unstable();
            column
        })
        .collect();
    rank_gf2(columns, faces.len())
}

/// Rank over `Z/2Z` of a matrix given by its sparse columns.
///
/// Each column is a strictly ascending list of the row indices holding a 1.
/// Columns are reduced by XOR (symmetric difference) against earlier columns
/// that share the same lowest row until the lowest row is unique; every column
/// that survives non-empty contributes one to the rank.
fn rank_gf2(columns: Vec<Vec<usize>>, n_rows: usize) -> usize {
    let mut pivot_owner = vec![usize::MAX; n_rows];
    let mut reduced: Vec<Vec<usize>> = Vec::new();
    for mut column in columns {
        while let Some(&lowest) = column.last() {
            let owner = pivot_owner[lowest];
            if owner == usize::MAX {
                pivot_owner[lowest] = reduced.len();
                reduced.push(column);
                break;
            }
            column = symmetric_difference(&column, &reduced[owner]);
        }
    }
    reduced.len()
}

/// Symmetric difference of two strictly ascending index lists.
fn symmetric_difference(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::fixtures::{
        Fixture, connected_sum, disjoint_union, grid_torus, icosahedron, octahedron,
        seven_vertex_torus, tetrahedron,
    };
    use crate::topology::types::Simplex;

    fn complex_of_triangles(triangles: &[[usize; 3]]) -> SimplicialComplex {
        let mut complex = SimplicialComplex::new(2);
        for tri in triangles {
            complex.add_simplex(Simplex::new(tri.to_vec()));
        }
        complex
    }

    fn betti_of(fixture: &Fixture) -> (usize, usize, usize) {
        betti_numbers(&complex_of_triangles(&fixture.1))
    }

    /// Euler-Poincare: the alternating sum of Betti numbers equals the
    /// alternating sum of simplex counts.
    fn assert_euler_poincare(complex: &SimplicialComplex) {
        let betti = betti_numbers_by_dimension(complex);
        let alternating: i64 = betti
            .iter()
            .enumerate()
            .map(|(k, &b)| if k % 2 == 0 { b as i64 } else { -(b as i64) })
            .sum();
        assert_eq!(
            alternating,
            complex.euler_characteristic(),
            "Betti numbers {betti:?}"
        );
    }

    #[test]
    fn empty_complex_has_no_homology() {
        let complex = SimplicialComplex::new(2);
        assert_eq!(betti_numbers(&complex), (0, 0, 0));
        assert!(betti_numbers_by_dimension(&complex).is_empty());
    }

    #[test]
    fn points_and_graphs() {
        let mut points = SimplicialComplex::new(2);
        points.add_simplex(Simplex::new(vec![0]));
        assert_eq!(betti_numbers(&points), (1, 0, 0));
        points.add_simplex(Simplex::new(vec![7]));
        assert_eq!(betti_numbers(&points), (2, 0, 0));

        let mut path = SimplicialComplex::new(2);
        for edge in [[0, 1], [1, 2], [2, 3]] {
            path.add_simplex(Simplex::new(edge.to_vec()));
        }
        assert_eq!(betti_numbers(&path), (1, 0, 0));

        // A hollow triangle is a circle.
        let mut circle = SimplicialComplex::new(2);
        for edge in [[0, 1], [1, 2], [0, 2]] {
            circle.add_simplex(Simplex::new(edge.to_vec()));
        }
        assert_eq!(betti_numbers(&circle), (1, 1, 0));

        // Two circles sharing a vertex (figure eight), and two disjoint circles.
        let mut eight = circle.clone();
        for edge in [[0, 3], [3, 4], [0, 4]] {
            eight.add_simplex(Simplex::new(edge.to_vec()));
        }
        assert_eq!(betti_numbers(&eight), (1, 2, 0));
        let mut apart = circle.clone();
        for edge in [[5, 6], [6, 7], [5, 7]] {
            apart.add_simplex(Simplex::new(edge.to_vec()));
        }
        assert_eq!(betti_numbers(&apart), (2, 2, 0));
        assert_euler_poincare(&eight);
    }

    #[test]
    fn filled_triangle_is_contractible() {
        // The old Euler-characteristic heuristic reported (1, 1, 1) here.
        let complex = complex_of_triangles(&[[0, 1, 2]]);
        assert_eq!(complex.euler_characteristic(), 1);
        assert_eq!(betti_numbers(&complex), (1, 0, 0));
        assert_euler_poincare(&complex);
    }

    #[test]
    fn two_disjoint_filled_triangles() {
        // The old heuristic reported (2, 1, 1).
        let complex = complex_of_triangles(&[[0, 1, 2], [3, 4, 5]]);
        assert_eq!(betti_numbers(&complex), (2, 0, 0));
        assert_euler_poincare(&complex);
    }

    #[test]
    fn spheres() {
        for (name, fixture) in [
            ("tetrahedron", tetrahedron()),
            ("octahedron", octahedron()),
            ("icosahedron", icosahedron()),
        ] {
            assert_eq!(betti_of(&fixture), (1, 0, 1), "{name}");
            assert_euler_poincare(&complex_of_triangles(&fixture.1));
        }
    }

    #[test]
    fn solid_tetrahedron_has_no_void() {
        let mut solid = SimplicialComplex::new(3);
        solid.add_simplex(Simplex::new(vec![0, 1, 2, 3]));
        assert_eq!(betti_numbers(&solid), (1, 0, 0));
        assert_eq!(betti_numbers_by_dimension(&solid), vec![1, 0, 0, 0]);
        let hollow = complex_of_triangles(&tetrahedron().1);
        assert_eq!(betti_numbers(&hollow), (1, 0, 1));
        assert_euler_poincare(&solid);
    }

    #[test]
    fn tori() {
        assert_eq!(betti_of(&seven_vertex_torus()), (1, 2, 1));
        assert_eq!(betti_of(&grid_torus(4, 4)), (1, 2, 1));
        assert_eq!(betti_of(&grid_torus(5, 3)), (1, 2, 1));
        assert_euler_poincare(&complex_of_triangles(&grid_torus(4, 4).1));
    }

    #[test]
    fn genus_two_surfaces() {
        for (name, double) in [
            (
                "grid tori",
                connected_sum(&grid_torus(4, 4), &grid_torus(4, 4)),
            ),
            (
                "Moebius-Csaszar tori",
                connected_sum(&seven_vertex_torus(), &seven_vertex_torus()),
            ),
        ] {
            let complex = complex_of_triangles(&double.1);
            assert_eq!(complex.euler_characteristic(), -2, "{name}");
            assert_eq!(betti_numbers(&complex), (1, 4, 1), "{name}");
            assert_euler_poincare(&complex);
        }
    }

    #[test]
    fn disjoint_surfaces_add_up() {
        let both = disjoint_union(&octahedron(), &grid_torus(4, 4));
        assert_eq!(betti_of(&both), (2, 2, 2));
    }

    #[test]
    fn projective_plane_depends_on_the_coefficient_field() {
        // The 6-vertex projective plane (hemi-icosahedron): over Z/2Z both
        // beta_1 and beta_2 are 1, whereas over the rationals they vanish.
        let faces = [
            [0, 1, 2],
            [0, 1, 4],
            [0, 2, 3],
            [0, 3, 5],
            [0, 4, 5],
            [1, 2, 5],
            [1, 3, 4],
            [1, 3, 5],
            [2, 3, 4],
            [2, 4, 5],
        ];
        let complex = complex_of_triangles(&faces);
        assert_eq!(complex.euler_characteristic(), 1);
        assert_eq!(betti_numbers(&complex), (1, 1, 1));
        assert_euler_poincare(&complex);
    }

    #[test]
    fn higher_dimensional_spheres() {
        // Boundary of the 4-simplex: a 3-sphere.
        let mut sphere = SimplicialComplex::new(3);
        for skip in 0..5 {
            let tet: Vec<usize> = (0..5).filter(|&v| v != skip).collect();
            sphere.add_simplex(Simplex::new(tet));
        }
        assert_eq!(betti_numbers_by_dimension(&sphere), vec![1, 0, 0, 1]);
        assert_eq!(betti_numbers(&sphere), (1, 0, 0));
        assert_euler_poincare(&sphere);
        // Filling the 4-simplex makes it contractible.
        let mut ball = SimplicialComplex::new(4);
        ball.add_simplex(Simplex::new(vec![0, 1, 2, 3, 4]));
        assert_eq!(betti_numbers_by_dimension(&ball), vec![1, 0, 0, 0, 0]);
    }

    #[test]
    fn complexes_missing_faces_are_closed_before_counting() {
        // Fill the public field directly with a triangle and nothing else.
        let mut open = SimplicialComplex::new(2);
        open.simplices[2].insert(Simplex::new(vec![0, 1, 2]));
        assert!(!open.is_valid());
        assert_eq!(betti_numbers(&open), (1, 0, 0));
        // Unsorted vertex lists (bypassing `Simplex::new`) are normalised.
        let mut unsorted = SimplicialComplex::new(2);
        unsorted.simplices[2].insert(Simplex {
            vertices: vec![2, 0, 1],
        });
        assert_eq!(betti_numbers(&unsorted), (1, 0, 0));
    }

    #[test]
    fn result_does_not_depend_on_insertion_order() {
        let fixture = connected_sum(&grid_torus(4, 4), &grid_torus(4, 4));
        let forward = complex_of_triangles(&fixture.1);
        let reversed: Vec<[usize; 3]> = fixture.1.iter().rev().copied().collect();
        let backward = complex_of_triangles(&reversed);
        for _ in 0..4 {
            assert_eq!(betti_numbers(&forward), (1, 4, 1));
            assert_eq!(betti_numbers(&backward), (1, 4, 1));
        }
    }

    #[test]
    fn larger_complex_is_handled() {
        let complex = complex_of_triangles(&grid_torus(24, 24).1);
        assert_eq!(betti_numbers(&complex), (1, 2, 1));
    }

    #[test]
    fn gf2_rank_of_dependent_columns() {
        // Columns 0 + 1 = 2 over Z/2Z, so the rank is 2, not 3.
        let columns = vec![vec![0, 1], vec![1, 2], vec![0, 2]];
        assert_eq!(rank_gf2(columns, 3), 2);
        assert_eq!(rank_gf2(vec![vec![], vec![0]], 1), 1);
        assert_eq!(symmetric_difference(&[0, 2, 4], &[2, 3]), vec![0, 3, 4]);
    }
}
