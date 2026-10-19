// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bridging-domain coupling of a meshfree (EFG) region with a finite element region.
//!
//! In the bridging-domain method the two models overlap in a strip.  Inside the
//! strip the strain energy is shared between them with a partition function
//! `α(x)` that rises from `0` (finite element model only) to `1` (meshfree model
//! only), and the two displacement fields are tied together:
//!
//! ```text
//! Π = Σ_FEM (1 − α) W_fem + Σ_MF α W_mf + ½ p Σ_pairs (u_mf − u_fem)².
//! ```
//!
//! [`setup_bridging_domain`] finds the overlap nodes, pairs every meshfree node
//! with its nearest finite element node and assigns `α` from the position
//! ([`CouplingInterface::compute_bridging_weights`]);
//! [`CouplingInterface::apply_bridging_blend`] applies the energy shares to the
//! two stiffness matrices and [`CouplingInterface::compute_penalty_coupling`] /
//! [`CouplingInterface::compute_field_penalty_coupling`] add the tie.
//!
//! The partition `α` is stored per coupling node as
//! [`CouplingNode::blend_weight`](super::types::CouplingNode::blend_weight): it is
//! always the **meshfree** share (`0` = pure FEM, `1` = pure meshfree); the finite
//! element share `1 − α` is formed where it is needed.

use super::functions::{Point2, dist2, mls_shape_2d};
use super::types::{CouplingInterface, CouplingMethod, CouplingNode, KernelType, SupportDomain};
use crate::{Error, Result};

/// Displacement components per node (plane problems).
const DOFS_PER_NODE: usize = 2;

/// Smallest overlap width, so that a zero-width overlap becomes a sharp step.
const MIN_OVERLAP_WIDTH: f64 = 1.0e-15;

/// Penalty parameter [`setup_bridging_domain`] gives the tie.
const DEFAULT_COUPLING_PENALTY: f64 = 1.0e6;

/// Fraction of a meshfree node's support radius within which a finite element
/// node counts as its coincident partner.
const PAIRING_RADIUS_FRACTION: f64 = 0.5;

/// Meshfree share `α(x) = clamp((x − x_min)/width, 0, 1)` of a point.
fn meshfree_share(x: f64, overlap_x_min: f64, overlap_width: f64) -> f64 {
    ((x - overlap_x_min) / overlap_width.max(MIN_OVERLAP_WIDTH)).clamp(0.0, 1.0)
}

impl CouplingInterface {
    /// Compute blending weights for the bridging domain method.
    ///
    /// In the overlap region, the total energy is partitioned:
    ///   E_total = w * E_meshfree + (1-w) * E_fem
    ///
    /// Every coupling node gets the meshfree share
    /// `w = clamp((x − x_min)/overlap_width, 0, 1)` from its position `x`, where
    /// `x_min` is [`CouplingInterface::overlap_x_min`] (set by
    /// [`setup_bridging_domain`]): `0` at the finite element edge of the overlap
    /// and `1` at its meshfree edge.  A non-positive width gives a sharp step at
    /// `x_min`.  The weights are consumed by
    /// [`apply_bridging_blend`](Self::apply_bridging_blend).
    pub fn compute_bridging_weights(&mut self, overlap_width: f64) {
        let x_min = self.overlap_x_min;
        for cn in &mut self.nodes {
            cn.blend_weight = meshfree_share(cn.pos[0], x_min, overlap_width);
        }
    }

    /// Per-node energy shares `(α_meshfree, α_fem)` of the two models.
    ///
    /// A node that is not part of the interface keeps the share `1` in its own
    /// model.  Every node of either model may appear in at most one coupling node.
    fn energy_shares(
        &self,
        n_meshfree_nodes: usize,
        n_fem_nodes: usize,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut meshfree = vec![1.0_f64; n_meshfree_nodes];
        let mut fem = vec![1.0_f64; n_fem_nodes];
        let mut meshfree_seen = vec![false; n_meshfree_nodes];
        let mut fem_seen = vec![false; n_fem_nodes];
        for (c, cn) in self.nodes.iter().enumerate() {
            let w = cn.blend_weight;
            if !(0.0..=1.0).contains(&w) {
                return Err(Error::InvalidInput(format!(
                    "coupling node {c}: the blend weight {w} is not in [0, 1]"
                )));
            }
            if let Some(i) = cn.meshfree_idx {
                if i >= n_meshfree_nodes {
                    return Err(Error::InvalidInput(format!(
                        "coupling node {c}: meshfree node {i} does not exist \
                         (the meshfree model has {n_meshfree_nodes} nodes)"
                    )));
                }
                if std::mem::replace(&mut meshfree_seen[i], true) {
                    return Err(Error::InvalidInput(format!(
                        "meshfree node {i} appears in more than one coupling node"
                    )));
                }
                meshfree[i] = w;
            }
            if let Some(j) = cn.fem_idx {
                if j >= n_fem_nodes {
                    return Err(Error::InvalidInput(format!(
                        "coupling node {c}: finite element node {j} does not exist \
                         (the finite element model has {n_fem_nodes} nodes)"
                    )));
                }
                if std::mem::replace(&mut fem_seen[j], true) {
                    return Err(Error::InvalidInput(format!(
                        "finite element node {j} appears in more than one coupling node"
                    )));
                }
                fem[j] = 1.0 - w;
            }
        }
        Ok((meshfree, fem))
    }

    /// Apply the bridging-domain energy partition to the two stiffness matrices.
    ///
    /// The meshfree matrix carries the share `α` of every coupling node and the
    /// finite element matrix the complementary share `1 − α`; nodes outside the
    /// interface keep the share `1` in their own model.  With the nodal shares
    /// `α_I` the entries are scaled symmetrically,
    ///
    /// `K̃_IJ = ½ (α_I + α_J) K_IJ`,
    ///
    /// the Hessian of the nodally weighted energy `½ Σ_I α_I u_I·(K u)_I`.  It
    /// reduces to `α K` where the share is uniform, keeps a symmetric matrix
    /// symmetric (scaling the rows alone would not) and, for the same `K` on
    /// coincident nodes, the two blended matrices add up to the original, so the
    /// strain energy is neither lost nor counted twice.  It is a nodal
    /// approximation of the exact scheme that weights the energy density at the
    /// quadrature points; the two agree to the variation of `α` across one node
    /// spacing.
    ///
    /// Both matrices are dense, row-major, with two DOFs per node (`ndof = 2 ×
    /// nodes`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] if a matrix does not match its DOF count, a DOF
    /// count is odd, a coupling node refers to a node that does not exist or a node
    /// appears in two coupling nodes, or a blend weight lies outside `[0, 1]`
    /// (call [`compute_bridging_weights`](Self::compute_bridging_weights) first).
    /// The matrices are left unchanged on error.
    pub fn apply_bridging_blend(
        &self,
        k_meshfree: &mut [f64],
        ndof_meshfree: usize,
        k_fem: &mut [f64],
        ndof_fem: usize,
    ) -> Result<()> {
        for (name, ndof, len) in [
            ("meshfree", ndof_meshfree, k_meshfree.len()),
            ("finite element", ndof_fem, k_fem.len()),
        ] {
            if ndof % DOFS_PER_NODE != 0 || len != ndof * ndof {
                return Err(Error::InvalidInput(format!(
                    "the {name} stiffness matrix must be {ndof} x {ndof} with an even number \
                     of DOFs, got {len} entries"
                )));
            }
        }
        let (alpha_meshfree, alpha_fem) =
            self.energy_shares(ndof_meshfree / DOFS_PER_NODE, ndof_fem / DOFS_PER_NODE)?;
        blend_symmetric(k_meshfree, ndof_meshfree, &alpha_meshfree);
        blend_symmetric(k_fem, ndof_fem, &alpha_fem);
        Ok(())
    }

    /// Penalty tie of the *field* of the meshfree model to the finite element
    /// displacement at every paired coupling node.
    ///
    /// MLS shape functions are not interpolatory, so the solution vector of the
    /// meshfree model holds parameters `û_J` and the displacement at a point is
    /// `u_h(x) = Σ_J φ_J(x) û_J`.  Tying the parameter of one node to a finite
    /// element value (as [`compute_penalty_coupling`](Self::compute_penalty_coupling)
    /// does) is therefore not the continuity condition `u_h(x_c) = u_fem`.  For
    /// every paired coupling node `c` (at position `x_c`) and component `d` this
    /// adds the penalty `½ p (Σ_J φ_J(x_c) û_{J,d} − u_{fem,d})²`, i.e. `p a aᵀ`
    /// with `a = [φ_J(x_c) at the meshfree DOFs, −1 at the FEM DOF]`.
    ///
    /// The global DOF numbering is the one of
    /// [`compute_penalty_coupling`](Self::compute_penalty_coupling): the meshfree
    /// DOFs `0..ndof_meshfree` followed by the FEM DOFs.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] if `ndof_meshfree` is not twice the node count of
    /// `domain`, a paired finite element node lies outside the FEM DOFs, or the MLS
    /// shape functions do not exist at a coupling node.
    pub fn compute_field_penalty_coupling(
        &self,
        domain: &SupportDomain,
        ktype: KernelType,
        ndof_meshfree: usize,
        ndof_fem: usize,
    ) -> Result<Vec<f64>> {
        if ndof_meshfree != DOFS_PER_NODE * domain.node_count() {
            return Err(Error::InvalidInput(format!(
                "the meshfree model has {} DOFs but {ndof_meshfree} were given",
                DOFS_PER_NODE * domain.node_count()
            )));
        }
        let ndof_total = ndof_meshfree + ndof_fem;
        let mut coupling_k = vec![0.0; ndof_total * ndof_total];
        for (c, cn) in self.nodes.iter().enumerate() {
            let (Some(_), Some(fi)) = (cn.meshfree_idx, cn.fem_idx) else {
                continue;
            };
            if DOFS_PER_NODE * (fi + 1) > ndof_fem {
                return Err(Error::InvalidInput(format!(
                    "coupling node {c}: finite element node {fi} lies outside the {ndof_fem} \
                     finite element DOFs"
                )));
            }
            let shape = mls_shape_2d(cn.pos, domain, ktype).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "coupling node {c}: the MLS shape functions are undefined at ({}, {})",
                    cn.pos[0], cn.pos[1]
                ))
            })?;
            for dof in 0..DOFS_PER_NODE {
                let mut row: Vec<(usize, f64)> = shape
                    .neighbours
                    .iter()
                    .zip(shape.phi.iter())
                    .map(|(&node, &phi)| (DOFS_PER_NODE * node + dof, phi))
                    .collect();
                row.push((ndof_meshfree + DOFS_PER_NODE * fi + dof, -1.0));
                for &(r, a_r) in &row {
                    for &(s, a_s) in &row {
                        coupling_k[r * ndof_total + s] += self.penalty * a_r * a_s;
                    }
                }
            }
        }
        Ok(coupling_k)
    }
}

/// `K_IJ ← ½ (α_I + α_J) K_IJ` for a dense row-major matrix with two DOFs per node.
fn blend_symmetric(k: &mut [f64], ndof: usize, alpha_node: &[f64]) {
    for (row, k_row) in k.chunks_exact_mut(ndof.max(1)).enumerate() {
        let alpha_row = alpha_node[row / DOFS_PER_NODE];
        for (col, entry) in k_row.iter_mut().enumerate() {
            *entry *= 0.5 * (alpha_row + alpha_node[col / DOFS_PER_NODE]);
        }
    }
}

/// Set up bridging domain coupling between a meshfree domain and FEM mesh.
///
/// Identifies the nodes of both models in the overlap strip
/// `overlap_x_min ≤ x ≤ overlap_x_max` and assigns the meshfree energy share
/// `α(x) = (x − overlap_x_min)/(overlap_x_max − overlap_x_min)`.
///
/// Every overlap meshfree node is paired with its nearest overlap FEM node when
/// that node is closer than half the (dilated) support radius of the meshfree
/// node; pairs are formed greedily by increasing distance so that each node
/// belongs to at most one pair.  A pair becomes one [`CouplingNode`] holding both
/// indices (at the position of the meshfree node), which is what the penalty tie
/// of [`CouplingInterface::compute_penalty_coupling`] acts on.  Overlap nodes
/// without a partner stay in the interface with a single index: they carry no tie
/// but still need their energy share in
/// [`CouplingInterface::apply_bridging_blend`].
///
/// The interface is returned with a penalty of `1e6`, and
/// [`CouplingInterface::overlap_x_min`] set to `overlap_x_min`.
pub fn setup_bridging_domain(
    meshfree_domain: &SupportDomain,
    fem_nodes: &[Point2],
    overlap_x_min: f64,
    overlap_x_max: f64,
) -> CouplingInterface {
    let mut interface =
        CouplingInterface::new(CouplingMethod::BridgingDomain, DEFAULT_COUPLING_PENALTY);
    interface.overlap_x_min = overlap_x_min;
    let width = (overlap_x_max - overlap_x_min).max(MIN_OVERLAP_WIDTH);
    let in_overlap = |x: f64| (overlap_x_min..=overlap_x_max).contains(&x);

    let meshfree_overlap: Vec<usize> = meshfree_domain
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, node)| in_overlap(node.pos[0]))
        .map(|(i, _)| i)
        .collect();
    let fem_overlap: Vec<usize> = fem_nodes
        .iter()
        .enumerate()
        .filter(|(_, pos)| in_overlap(pos[0]))
        .map(|(j, _)| j)
        .collect();

    // Candidate pairs (distance, meshfree node, FEM node) within the pairing radius.
    let mut candidates: Vec<(f64, usize, usize)> = Vec::new();
    for &i in &meshfree_overlap {
        let node = &meshfree_domain.nodes[i];
        let radius = PAIRING_RADIUS_FRACTION * node.support_radius * meshfree_domain.dilation;
        for &j in &fem_overlap {
            let distance = dist2(node.pos, fem_nodes[j]);
            if distance < radius {
                candidates.push((distance, i, j));
            }
        }
    }
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut partner_of_meshfree: Vec<Option<usize>> = vec![None; meshfree_domain.nodes.len()];
    let mut fem_taken = vec![false; fem_nodes.len()];
    for &(_, i, j) in &candidates {
        if partner_of_meshfree[i].is_none() && !fem_taken[j] {
            partner_of_meshfree[i] = Some(j);
            fem_taken[j] = true;
        }
    }

    for &i in &meshfree_overlap {
        let pos = meshfree_domain.nodes[i].pos;
        interface.add_node(CouplingNode {
            pos,
            meshfree_idx: Some(i),
            fem_idx: partner_of_meshfree[i],
            blend_weight: meshfree_share(pos[0], overlap_x_min, width),
        });
    }
    for &j in &fem_overlap {
        if !fem_taken[j] {
            let pos = fem_nodes[j];
            interface.add_node(CouplingNode {
                pos,
                meshfree_idx: None,
                fem_idx: Some(j),
                blend_weight: meshfree_share(pos[0], overlap_x_min, width),
            });
        }
    }
    interface
}
