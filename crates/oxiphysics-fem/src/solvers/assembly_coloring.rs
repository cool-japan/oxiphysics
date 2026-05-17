// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Element graph coloring for race-free parallel FEM assembly.
//!
//! Elements sharing a DOF cannot be assembled simultaneously without data races.
//! Graph coloring assigns each element a color such that no two elements of the
//! same color share any DOF.  Elements within the same color can then be
//! processed in parallel without synchronization.

#![allow(dead_code)]

use std::sync::Mutex;

use crate::parallel_solver::{AssemblyTask, CsrMatrix};

/// Result of coloring the element graph.
#[derive(Debug, Clone)]
pub struct ElementColoringResult {
    /// `colors[elem_idx]` — color (0-based) assigned to this element.
    pub colors: Vec<usize>,
    /// Total number of colors used.
    pub n_colors: usize,
    /// `buckets[color]` — list of element indices with that color.
    pub buckets: Vec<Vec<usize>>,
}

/// Maximum number of graph colors before switching to serial fallback.
const MAX_COLORS: usize = 64;

/// Assign colors to elements such that no two elements of the same color share a DOF.
///
/// Uses a greedy graph coloring algorithm:
/// 1. Build DOF → element adjacency.
/// 2. For each element, collect the set of colors used by adjacent elements (sharing a DOF).
/// 3. Assign the smallest unused color.
/// 4. If more than `MAX_COLORS` would be needed, all remaining elements go into a fallback bucket.
///
/// # Arguments
/// * `n_elements` — total number of elements.
/// * `element_dofs` — for each element, the list of global DOF indices it uses.
pub fn color_elements(n_elements: usize, element_dofs: &[Vec<usize>]) -> ElementColoringResult {
    if n_elements == 0 {
        return ElementColoringResult {
            colors: Vec::new(),
            n_colors: 0,
            buckets: Vec::new(),
        };
    }

    // Determine the maximum DOF index to size the adjacency structure
    let max_dof = element_dofs
        .iter()
        .flat_map(|dofs| dofs.iter())
        .copied()
        .max()
        .unwrap_or(0);

    // DOF → list of elements that use this DOF
    let mut dof_to_elems: Vec<Vec<usize>> = vec![Vec::new(); max_dof + 1];
    for (e, dofs) in element_dofs.iter().enumerate() {
        for &d in dofs {
            dof_to_elems[d].push(e);
        }
    }

    let mut colors = vec![MAX_COLORS; n_elements]; // MAX_COLORS = unassigned sentinel
    let mut n_colors_used = 0usize;

    for e in 0..n_elements {
        // Collect colors used by adjacent elements (those sharing a DOF with e)
        let mut neighbor_colors = std::collections::HashSet::new();
        for &d in &element_dofs[e] {
            for &adj_elem in &dof_to_elems[d] {
                if adj_elem != e && colors[adj_elem] < MAX_COLORS {
                    neighbor_colors.insert(colors[adj_elem]);
                }
            }
        }

        // Find smallest unused color
        let mut chosen_color = 0;
        while neighbor_colors.contains(&chosen_color) {
            chosen_color += 1;
        }

        if chosen_color >= MAX_COLORS {
            // Fallback: use the last bucket (MAX_COLORS - 1)
            eprintln!(
                "assembly_coloring: element {e} needs color {chosen_color} >= MAX_COLORS={MAX_COLORS}, \
                 using fallback serial bucket"
            );
            chosen_color = MAX_COLORS - 1;
        }

        colors[e] = chosen_color;
        if chosen_color + 1 > n_colors_used {
            n_colors_used = chosen_color + 1;
        }
    }

    // Build buckets
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n_colors_used];
    for (e, &c) in colors.iter().enumerate() {
        if c < n_colors_used {
            buckets[c].push(e);
        }
    }

    ElementColoringResult {
        colors,
        n_colors: n_colors_used,
        buckets,
    }
}

/// Assemble a global CSR stiffness matrix using element graph coloring.
///
/// Elements within the same color share no DOFs, so they can be processed in
/// parallel without contention.  The Mutex-per-CSR-entry approach is used for
/// correctness; since elements of the same color never share a CSR entry, no
/// Mutex is ever contended within a color.
///
/// # Arguments
/// * `ndofs` — total number of global DOFs.
/// * `tasks` — per-element stiffness contributions.
/// * `coloring` — result of `color_elements`.
pub fn assemble_colored_csr(
    ndofs: usize,
    tasks: &[AssemblyTask],
    coloring: &ElementColoringResult,
) -> CsrMatrix {
    // ── Step 1: Build sparsity pattern (sequential) ──────────────────────────
    let mut row_cols: Vec<std::collections::BTreeSet<usize>> =
        vec![std::collections::BTreeSet::new(); ndofs];

    for task in tasks {
        for &row_dof in &task.global_dofs {
            for &col_dof in &task.global_dofs {
                row_cols[row_dof].insert(col_dof);
            }
        }
    }

    let mut row_offsets = vec![0usize; ndofs + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    for (i, cols) in row_cols.iter().enumerate() {
        row_offsets[i + 1] = row_offsets[i] + cols.len();
        col_indices.extend(cols.iter().copied());
    }
    let nnz = col_indices.len();

    // ── Step 2: Per-row lookup table (col → CSR index) ───────────────────────
    let row_col_to_csr: Vec<std::collections::HashMap<usize, usize>> = row_cols
        .iter()
        .enumerate()
        .map(|(i, cols)| {
            let base = row_offsets[i];
            cols.iter()
                .enumerate()
                .map(|(j, &c)| (c, base + j))
                .collect()
        })
        .collect();

    // ── Step 3: Assemble color by color using Mutex per CSR entry ────────────
    // Since elements within the same color share no DOFs, Mutexes are never
    // contended within a single color's Rayon parallel loop.
    let values_locked: Vec<Mutex<f64>> = (0..nnz).map(|_| Mutex::new(0.0f64)).collect();

    use rayon::prelude::*;

    for bucket in &coloring.buckets {
        bucket.par_iter().for_each(|&e| {
            let task = &tasks[e];
            let ndof = task.ndof();
            for (li, &row) in task.global_dofs.iter().enumerate() {
                for (lj, &col) in task.global_dofs.iter().enumerate() {
                    let ke_val = task.ke[li * ndof + lj];
                    if let Some(&csr_idx) = row_col_to_csr[row].get(&col) {
                        // Elements in the same color never share CSR entries
                        // (guaranteed by graph coloring), so this Mutex is uncontended.
                        let mut guard = values_locked[csr_idx]
                            .lock()
                            .unwrap_or_else(|p| p.into_inner());
                        *guard += ke_val;
                    }
                }
            }
        });
    }

    let values: Vec<f64> = values_locked
        .into_iter()
        .map(|m| m.into_inner().unwrap_or_else(|p| p.into_inner()))
        .collect();

    CsrMatrix {
        nrows: ndofs,
        ncols: ndofs,
        row_offsets,
        col_indices,
        values,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel_solver::ParallelAssembler;

    fn make_test_tasks() -> (Vec<AssemblyTask>, Vec<Vec<usize>>) {
        // 4 elements in a 1D bar: DOFs [0,1], [1,2], [2,3], [3,4]
        // Element stiffness: [[1,-1],[-1,1]]
        let dofs = vec![vec![0, 1], vec![1, 2], vec![2, 3], vec![3, 4]];
        let ke = vec![1.0, -1.0, -1.0, 1.0];
        let tasks: Vec<AssemblyTask> = dofs
            .iter()
            .map(|d| AssemblyTask::new(d.clone(), ke.clone()))
            .collect();
        (tasks, dofs)
    }

    #[test]
    fn test_coloring_valid() {
        let (tasks, dofs) = make_test_tasks();
        let coloring = color_elements(tasks.len(), &dofs);

        // Verify: no two elements of the same color share a DOF
        for (color, bucket) in coloring.buckets.iter().enumerate() {
            for &e1 in bucket {
                for &e2 in bucket {
                    if e1 == e2 {
                        continue;
                    }
                    let dofs1: std::collections::HashSet<usize> =
                        dofs[e1].iter().copied().collect();
                    let dofs2: std::collections::HashSet<usize> =
                        dofs[e2].iter().copied().collect();
                    let overlap: std::collections::HashSet<usize> =
                        dofs1.intersection(&dofs2).copied().collect();
                    assert!(
                        overlap.is_empty(),
                        "Color {color}: elements {e1} and {e2} share DOFs {overlap:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_colored_assembly_matches_serial() {
        let (tasks, dofs) = make_test_tasks();
        let ndofs = 5;

        // Serial assembly via ParallelAssembler
        let asm = ParallelAssembler::new(ndofs);
        let mat_serial = asm.assemble(&tasks);

        // Colored assembly
        let coloring = color_elements(tasks.len(), &dofs);
        let mat_colored = assemble_colored_csr(ndofs, &tasks, &coloring);

        // Compare: same sparsity pattern and values
        assert_eq!(mat_serial.nrows, mat_colored.nrows);
        assert_eq!(mat_serial.ncols, mat_colored.ncols);
        assert_eq!(mat_serial.nnz(), mat_colored.nnz());

        // Compare value at each (row, col) pair
        let get_val = |mat: &CsrMatrix, row: usize, col: usize| -> f64 {
            for k in mat.row_offsets[row]..mat.row_offsets[row + 1] {
                if mat.col_indices[k] == col {
                    return mat.values[k];
                }
            }
            0.0
        };

        for i in 0..ndofs {
            for k in mat_serial.row_offsets[i]..mat_serial.row_offsets[i + 1] {
                let j = mat_serial.col_indices[k];
                let v_ser = mat_serial.values[k];
                let v_col = get_val(&mat_colored, i, j);
                assert!(
                    (v_ser - v_col).abs() < 1e-14,
                    "Mismatch at ({i},{j}): serial={v_ser}, colored={v_col}"
                );
            }
        }
    }

    #[test]
    fn test_coloring_small_mesh() {
        // 4 elements sharing nodes:
        //   elem 0: dofs [0,1]
        //   elem 1: dofs [1,2]
        //   elem 2: dofs [2,3]
        //   elem 3: dofs [0,3]  (wraps around)
        let dofs = vec![vec![0usize, 1], vec![1, 2], vec![2, 3], vec![0, 3]];
        let coloring = color_elements(4, &dofs);

        assert_eq!(coloring.colors.len(), 4);
        assert!(coloring.n_colors >= 2, "Expected at least 2 colors");

        // Elements 0 and 2 don't share DOFs → can be same color
        // Elements 1 and 3 don't share DOFs → can be same color
        // Verify no two same-color elements share a DOF
        for bucket in &coloring.buckets {
            for &e1 in bucket {
                for &e2 in bucket {
                    if e1 == e2 {
                        continue;
                    }
                    let dofs1: std::collections::HashSet<usize> =
                        dofs[e1].iter().copied().collect();
                    let dofs2: std::collections::HashSet<usize> =
                        dofs[e2].iter().copied().collect();
                    let overlap: Vec<usize> = dofs1.intersection(&dofs2).copied().collect();
                    assert!(
                        overlap.is_empty(),
                        "Same-color elements {e1} and {e2} share DOFs {overlap:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_coloring_empty() {
        let result = color_elements(0, &[]);
        assert_eq!(result.n_colors, 0);
        assert_eq!(result.buckets.len(), 0);
    }
}
