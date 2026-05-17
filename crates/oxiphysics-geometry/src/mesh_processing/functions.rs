//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet, VecDeque};

use super::types::{
    AtlasPatch, BooleanResult, FeatureEdge, MeshQualityStats, ProcessMesh, Quadric,
    TriangleQuality, UvParameterization,
};

/// A 3-D vertex position stored as plain `[f64; 3]`.
pub type Vertex = [f64; 3];
/// A triangular face represented as three vertex indices.
pub type Face = [usize; 3];
/// A 2-D UV coordinate pair.
pub type Uv = [f64; 2];
#[inline]
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn vec3_len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
#[inline]
pub(super) fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / l)
    }
}
#[inline]
pub(super) fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    vec3_add(vec3_scale(a, 1.0 - t), vec3_scale(b, t))
}
/// Apply `iters` passes of uniform Laplacian smoothing with step `lambda`.
///
/// Each vertex is moved toward the average of its neighbours.
pub fn laplacian_smooth(mesh: &ProcessMesh, lambda: f64, iters: usize) -> ProcessMesh {
    let mut out = mesh.clone();
    let adj = mesh.build_adjacency();
    for _ in 0..iters {
        let old = out.verts.clone();
        for (i, nbrs) in adj.iter().enumerate() {
            if nbrs.is_empty() {
                continue;
            }
            let mut avg = [0.0_f64; 3];
            for &j in nbrs {
                avg = vec3_add(avg, old[j]);
            }
            avg = vec3_scale(avg, 1.0 / nbrs.len() as f64);
            let delta = vec3_sub(avg, old[i]);
            out.verts[i] = vec3_add(old[i], vec3_scale(delta, lambda));
        }
    }
    out
}
/// Apply Taubin λ/μ smoothing to reduce shrinkage.
///
/// Alternates a positive step (λ) with a negative step (μ < 0) so that
/// the mesh volume is approximately preserved.
pub fn taubin_smooth(mesh: &ProcessMesh, lambda: f64, mu: f64, iters: usize) -> ProcessMesh {
    let mut out = mesh.clone();
    let adj = mesh.build_adjacency();
    let one_step = |verts: &[Vertex], step: f64, a: &[Vec<usize>]| -> Vec<Vertex> {
        let mut next = verts.to_vec();
        for (i, nbrs) in a.iter().enumerate() {
            if nbrs.is_empty() {
                continue;
            }
            let mut avg = [0.0_f64; 3];
            for &j in nbrs {
                avg = vec3_add(avg, verts[j]);
            }
            avg = vec3_scale(avg, 1.0 / nbrs.len() as f64);
            let delta = vec3_sub(avg, verts[i]);
            next[i] = vec3_add(verts[i], vec3_scale(delta, step));
        }
        next
    };
    for _ in 0..iters {
        out.verts = one_step(&out.verts, lambda, &adj);
        out.verts = one_step(&out.verts, mu, &adj);
    }
    out
}
/// Simple midpoint (butterfly) subdivision: each edge is split and each
/// triangle becomes four smaller triangles.
pub fn midpoint_subdivide(mesh: &ProcessMesh) -> ProcessMesh {
    let mut verts = mesh.verts.clone();
    let mut faces = Vec::new();
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut get_mid = |a: usize, b: usize, v: &mut Vec<Vertex>| -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        *edge_mid.entry(key).or_insert_with(|| {
            let mid = vec3_scale(vec3_add(v[a], v[b]), 0.5);
            let idx = v.len();
            v.push(mid);
            idx
        })
    };
    for &[a, b, c] in &mesh.faces {
        let ab = get_mid(a, b, &mut verts);
        let bc = get_mid(b, c, &mut verts);
        let ca = get_mid(c, a, &mut verts);
        faces.push([a, ab, ca]);
        faces.push([ab, b, bc]);
        faces.push([ca, bc, c]);
        faces.push([ab, bc, ca]);
    }
    ProcessMesh::new(verts, faces)
}
/// Loop subdivision with proper weighting for interior vertices.
///
/// Implements Charles Loop's 1987 scheme for smooth subdivision surfaces.
pub fn loop_subdivide(mesh: &ProcessMesh) -> ProcessMesh {
    let mut verts = mesh.verts.clone();
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut edge_opp: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for &[a, b, c] in &mesh.faces {
        let edges = [(a, b, c), (b, c, a), (c, a, b)];
        for (i, j, k) in edges {
            let key = if i < j { (i, j) } else { (j, i) };
            edge_opp.entry(key).or_default().push(k);
        }
    }
    let mut get_mid = |a: usize, b: usize, v: &mut Vec<Vertex>| -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        *edge_mid.entry(key).or_insert_with(|| {
            let opps = edge_opp.get(&key).cloned().unwrap_or_default();
            let pos = if opps.len() == 2 {
                let sum_ab = vec3_scale(vec3_add(v[a], v[b]), 3.0 / 8.0);
                let sum_op = vec3_scale(vec3_add(v[opps[0]], v[opps[1]]), 1.0 / 8.0);
                vec3_add(sum_ab, sum_op)
            } else {
                vec3_scale(vec3_add(v[a], v[b]), 0.5)
            };
            let idx = v.len();
            v.push(pos);
            idx
        })
    };
    let adj = mesh.build_adjacency();
    let orig = mesh.verts.clone();
    let n = orig.len();
    let mut new_pos = orig.clone();
    for i in 0..n {
        let nbrs = &adj[i];
        let k = nbrs.len() as f64;
        let beta = if k <= 3.0 {
            3.0 / 16.0
        } else {
            3.0 / (8.0 * k)
        };
        let mut sum = [0.0_f64; 3];
        for &j in nbrs {
            sum = vec3_add(sum, orig[j]);
        }
        new_pos[i] = vec3_add(vec3_scale(orig[i], 1.0 - k * beta), vec3_scale(sum, beta));
    }
    let mut faces = Vec::new();
    for &[a, b, c] in &mesh.faces {
        let ab = get_mid(a, b, &mut verts);
        let bc = get_mid(b, c, &mut verts);
        let ca = get_mid(c, a, &mut verts);
        faces.push([a, ab, ca]);
        faces.push([ab, b, bc]);
        faces.push([ca, bc, c]);
        faces.push([ab, bc, ca]);
    }
    verts[..n].copy_from_slice(&new_pos[..n]);
    ProcessMesh::new(verts, faces)
}
/// Catmull-Clark subdivision adapted for triangle meshes.
///
/// For pure triangle meshes this is approximated by computing face centroids,
/// edge midpoints, and updated vertex positions following the CC weights.
pub fn catmull_clark_subdivide(mesh: &ProcessMesh) -> ProcessMesh {
    let nv = mesh.verts.len();
    let mut new_verts = mesh.verts.clone();
    let mut face_pts: Vec<Vertex> = Vec::with_capacity(mesh.faces.len());
    for &[a, b, c] in &mesh.faces {
        let fp = vec3_scale(
            vec3_add(vec3_add(mesh.verts[a], mesh.verts[b]), mesh.verts[c]),
            1.0 / 3.0,
        );
        face_pts.push(fp);
    }
    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut edge_pts: Vec<Vertex> = Vec::new();
    let mut vert_face_avg = vec![[0.0_f64; 3]; nv];
    let mut vert_face_cnt = vec![0usize; nv];
    for (fi, &[a, b, c]) in mesh.faces.iter().enumerate() {
        let fp = face_pts[fi];
        for &v in &[a, b, c] {
            vert_face_avg[v] = vec3_add(vert_face_avg[v], fp);
            vert_face_cnt[v] += 1;
        }
        for (i, j) in [(a, b), (b, c), (c, a)] {
            let key = if i < j { (i, j) } else { (j, i) };
            edge_map.entry(key).or_insert_with(|| {
                let ep = vec3_scale(vec3_add(mesh.verts[i], mesh.verts[j]), 0.5);
                let idx = edge_pts.len();
                edge_pts.push(ep);
                idx
            });
        }
    }
    let adj = mesh.build_adjacency();
    for i in 0..nv {
        let n = vert_face_cnt[i] as f64;
        if n == 0.0 {
            continue;
        }
        let f_avg = vec3_scale(vert_face_avg[i], 1.0 / n);
        let nbrs = &adj[i];
        let mut e_avg = [0.0_f64; 3];
        for &j in nbrs {
            e_avg = vec3_add(e_avg, mesh.verts[j]);
        }
        let m = nbrs.len() as f64;
        e_avg = vec3_scale(e_avg, 1.0 / m.max(1.0));
        new_verts[i] = vec3_scale(
            vec3_add(
                vec3_add(f_avg, vec3_scale(e_avg, 2.0)),
                vec3_scale(mesh.verts[i], (n - 3.0).max(0.0)),
            ),
            1.0 / n,
        );
    }
    let ep_offset = nv;
    let fp_offset = ep_offset + edge_pts.len();
    for ep in &edge_pts {
        new_verts.push(*ep);
    }
    for fp in &face_pts {
        new_verts.push(*fp);
    }
    let mut faces = Vec::new();
    for (fi, &[a, b, c]) in mesh.faces.iter().enumerate() {
        let fp_idx = fp_offset + fi;
        let ab_idx = ep_offset
            + *edge_map
                .get(&if a < b { (a, b) } else { (b, a) })
                .expect("key must exist");
        let bc_idx = ep_offset
            + *edge_map
                .get(&if b < c { (b, c) } else { (c, b) })
                .expect("key must exist");
        let ca_idx = ep_offset
            + *edge_map
                .get(&if c < a { (c, a) } else { (a, c) })
                .expect("key must exist");
        faces.push([a, ab_idx, fp_idx]);
        faces.push([ab_idx, b, fp_idx]);
        faces.push([b, bc_idx, fp_idx]);
        faces.push([bc_idx, c, fp_idx]);
        faces.push([c, ca_idx, fp_idx]);
        faces.push([ca_idx, a, fp_idx]);
    }
    ProcessMesh::new(new_verts, faces)
}
/// Build one quadric per face (the fundamental plane quadric).
pub(super) fn face_quadrics(mesh: &ProcessMesh) -> Vec<Quadric> {
    mesh.faces
        .iter()
        .map(|&[a, b, c]| {
            let va = mesh.verts[a];
            let vb = mesh.verts[b];
            let vc = mesh.verts[c];
            let ab = vec3_sub(vb, va);
            let ac = vec3_sub(vc, va);
            let n = vec3_normalize(vec3_cross(ab, ac));
            let d = -vec3_dot(n, va);
            Quadric::from_plane(n[0], n[1], n[2], d)
        })
        .collect()
}
/// Build per-vertex quadrics by summing face quadrics for all adjacent faces.
pub(super) fn vertex_quadrics(mesh: &ProcessMesh) -> Vec<Quadric> {
    let nv = mesh.verts.len();
    let fq = face_quadrics(mesh);
    let mut vq = vec![Quadric::zero(); nv];
    for (fi, &[a, b, c]) in mesh.faces.iter().enumerate() {
        for &v in &[a, b, c] {
            vq[v] = vq[v].add(&fq[fi]);
        }
    }
    vq
}
/// Decimate a mesh to at most `target_faces` triangles using quadric error
/// metrics with greedy edge collapses.
pub fn qem_decimate(mesh: &ProcessMesh, target_faces: usize) -> ProcessMesh {
    if mesh.faces.len() <= target_faces {
        return mesh.clone();
    }
    let nv = mesh.verts.len();
    let mut verts = mesh.verts.clone();
    let mut faces = mesh.faces.clone();
    let mut vq = vertex_quadrics(mesh);
    let mut removed_verts = vec![false; nv];
    let mut removed_faces = vec![false; faces.len()];
    while faces.iter().filter(|&&_f| true).count() - removed_faces.iter().filter(|&&r| r).count()
        > target_faces
    {
        let mut best_cost = f64::INFINITY;
        let mut best_edge: Option<(usize, usize)> = None;
        let mut best_pos = [0.0_f64; 3];
        let mut seen = HashSet::new();
        for &[a, b, c] in &faces {
            for (i, j) in [(a, b), (b, c), (c, a)] {
                let key = if i < j { (i, j) } else { (j, i) };
                if seen.contains(&key) {
                    continue;
                }
                seen.insert(key);
                if removed_verts[i] || removed_verts[j] {
                    continue;
                }
                let q = vq[i].add(&vq[j]);
                let mid = vec3_scale(vec3_add(verts[i], verts[j]), 0.5);
                let cost = q.evaluate(mid);
                if cost < best_cost {
                    best_cost = cost;
                    best_edge = Some((i, j));
                    best_pos = mid;
                }
            }
        }
        if let Some((vi, vj)) = best_edge {
            verts[vi] = best_pos;
            vq[vi] = vq[vi].add(&vq[vj]);
            removed_verts[vj] = true;
            for (fi, f) in faces.iter_mut().enumerate() {
                if removed_faces[fi] {
                    continue;
                }
                for idx in f.iter_mut() {
                    if *idx == vj {
                        *idx = vi;
                    }
                }
                if f[0] == f[1] || f[1] == f[2] || f[0] == f[2] {
                    removed_faces[fi] = true;
                }
            }
        } else {
            break;
        }
    }
    let remap: Vec<usize> = {
        let mut m = vec![0usize; nv];
        let mut cnt = 0;
        for i in 0..nv {
            if !removed_verts[i] {
                m[i] = cnt;
                cnt += 1;
            }
        }
        m
    };
    let new_verts: Vec<Vertex> = (0..nv)
        .filter(|&i| !removed_verts[i])
        .map(|i| verts[i])
        .collect();
    let new_faces: Vec<Face> = faces
        .iter()
        .enumerate()
        .filter(|&(fi, _)| !removed_faces[fi])
        .map(|(_, &[a, b, c])| [remap[a], remap[b], remap[c]])
        .collect();
    ProcessMesh::new(new_verts, new_faces)
}
/// Compute quality metrics for a single triangle given its three vertices.
pub fn triangle_quality(va: Vertex, vb: Vertex, vc: Vertex) -> TriangleQuality {
    let ab = vec3_sub(vb, va);
    let ac = vec3_sub(vc, va);
    let bc = vec3_sub(vc, vb);
    let la = vec3_len(bc);
    let lb = vec3_len(ac);
    let lc = vec3_len(ab);
    let s = (la + lb + lc) / 2.0;
    let area2 = vec3_len(vec3_cross(ab, ac));
    let area = area2 / 2.0;
    let inradius = if s > 1e-15 { area / s } else { 0.0 };
    let circumradius = if area > 1e-15 {
        la * lb * lc / (4.0 * area)
    } else {
        f64::INFINITY
    };
    let aspect_ratio = if inradius > 1e-15 {
        circumradius / (2.0 * inradius)
    } else {
        f64::INFINITY
    };
    let angle_a = if la > 1e-15 && lb > 1e-15 {
        let cos_a = (lb * lb + lc * lc - la * la) / (2.0 * lb * lc);
        cos_a.clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    let angle_b = if la > 1e-15 && lc > 1e-15 {
        let cos_b = (la * la + lc * lc - lb * lb) / (2.0 * la * lc);
        cos_b.clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    let angle_c = std::f64::consts::PI - angle_a - angle_b;
    let perim = 2.0 * s;
    let ideal_area = perim * perim * (3.0_f64.sqrt()) / 36.0;
    let area_quality = if ideal_area > 1e-15 {
        (area / ideal_area).min(1.0)
    } else {
        0.0
    };
    TriangleQuality {
        aspect_ratio,
        min_angle: angle_a.min(angle_b).min(angle_c),
        max_angle: angle_a.max(angle_b).max(angle_c),
        area_quality,
    }
}
/// Compute aggregate quality statistics for `mesh`.
pub fn mesh_quality_stats(mesh: &ProcessMesh) -> MeshQualityStats {
    if mesh.faces.is_empty() {
        return MeshQualityStats {
            mean_aspect_ratio: 0.0,
            max_aspect_ratio: 0.0,
            min_angle_deg: 0.0,
            mean_area_quality: 0.0,
        };
    }
    let n = mesh.faces.len() as f64;
    let mut sum_ar = 0.0;
    let mut max_ar = 0.0_f64;
    let mut min_ang = f64::INFINITY;
    let mut sum_aq = 0.0;
    for &[a, b, c] in &mesh.faces {
        let q = triangle_quality(mesh.verts[a], mesh.verts[b], mesh.verts[c]);
        sum_ar += q.aspect_ratio;
        max_ar = max_ar.max(q.aspect_ratio);
        min_ang = min_ang.min(q.min_angle);
        sum_aq += q.area_quality;
    }
    MeshQualityStats {
        mean_aspect_ratio: sum_ar / n,
        max_aspect_ratio: max_ar,
        min_angle_deg: min_ang.to_degrees(),
        mean_area_quality: sum_aq / n,
    }
}
/// Detect boundary (open) edge loops in `mesh`.
///
/// Returns a list of loops, where each loop is an ordered list of vertex
/// indices forming the boundary.
pub fn detect_boundary_loops(mesh: &ProcessMesh) -> Vec<Vec<usize>> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for &[a, b, c] in &mesh.faces {
        for (i, j) in [(a, b), (b, c), (c, a)] {
            *edge_count.entry((i, j)).or_insert(0) += 1;
        }
    }
    let mut boundary_adj: HashMap<usize, usize> = HashMap::new();
    for &(i, j) in edge_count.keys() {
        if edge_count.get(&(j, i)).copied().unwrap_or(0) == 0 {
            boundary_adj.insert(i, j);
        }
    }
    let mut visited = HashSet::new();
    let mut loops = Vec::new();
    for &start in boundary_adj.keys() {
        if visited.contains(&start) {
            continue;
        }
        let mut loop_ = Vec::new();
        let mut cur = start;
        loop {
            if visited.contains(&cur) {
                break;
            }
            visited.insert(cur);
            loop_.push(cur);
            match boundary_adj.get(&cur) {
                Some(&next) => cur = next,
                None => break,
            }
        }
        if !loop_.is_empty() {
            loops.push(loop_);
        }
    }
    loops
}
/// Fill a single boundary loop with a fan triangulation from its centroid.
///
/// Returns the new faces to append.
pub fn fill_hole_fan(verts: &mut Vec<Vertex>, loop_: &[usize]) -> Vec<Face> {
    if loop_.len() < 3 {
        return Vec::new();
    }
    let mut centroid = [0.0_f64; 3];
    for &vi in loop_ {
        centroid = vec3_add(centroid, verts[vi]);
    }
    centroid = vec3_scale(centroid, 1.0 / loop_.len() as f64);
    let center_idx = verts.len();
    verts.push(centroid);
    let mut new_faces = Vec::new();
    let n = loop_.len();
    for i in 0..n {
        new_faces.push([loop_[(i + 1) % n], loop_[i], center_idx]);
    }
    new_faces
}
/// Repair `mesh` by detecting all boundary loops and filling each with a fan.
pub fn fill_all_holes(mesh: &ProcessMesh) -> ProcessMesh {
    let mut out = mesh.clone();
    loop {
        let loops = detect_boundary_loops(&out);
        if loops.is_empty() {
            break;
        }
        for lp in &loops {
            let new_faces = fill_hole_fan(&mut out.verts, lp);
            out.faces.extend(new_faces);
        }
        if detect_boundary_loops(&out).len() >= loops.len() {
            break;
        }
    }
    out
}
/// Remove non-manifold vertices (vertices shared by more than two disconnected
/// fans) by duplicating them.
pub fn fix_non_manifold_vertices(mesh: &ProcessMesh) -> ProcessMesh {
    let nv = mesh.verts.len();
    let mut vert_faces: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (fi, &[a, b, c]) in mesh.faces.iter().enumerate() {
        vert_faces[a].push(fi);
        vert_faces[b].push(fi);
        vert_faces[c].push(fi);
    }
    let mut new_verts = mesh.verts.clone();
    let mut new_faces = mesh.faces.clone();
    for vi in 0..nv {
        let flist = &vert_faces[vi];
        if flist.len() < 2 {
            continue;
        }
        let face_set: HashSet<usize> = flist.iter().cloned().collect();
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(flist[0]);
        while let Some(fi) = queue.pop_front() {
            if visited.contains(&fi) {
                continue;
            }
            visited.insert(fi);
            let [a, b, c] = new_faces[fi];
            for &fj in &vert_faces[a] {
                if face_set.contains(&fj) && !visited.contains(&fj) {
                    queue.push_back(fj);
                }
            }
            for &fj in &vert_faces[b] {
                if face_set.contains(&fj) && !visited.contains(&fj) {
                    queue.push_back(fj);
                }
            }
            for &fj in &vert_faces[c] {
                if face_set.contains(&fj) && !visited.contains(&fj) {
                    queue.push_back(fj);
                }
            }
        }
        let extra: Vec<usize> = flist
            .iter()
            .filter(|&&f| !visited.contains(&f))
            .cloned()
            .collect();
        if !extra.is_empty() {
            let new_idx = new_verts.len();
            new_verts.push(new_verts[vi]);
            for fi in extra {
                for idx in new_faces[fi].iter_mut() {
                    if *idx == vi {
                        *idx = new_idx;
                    }
                }
            }
        }
    }
    ProcessMesh::new(new_verts, new_faces)
}
/// Detect feature edges in `mesh` whose dihedral angle exceeds
/// `threshold_deg` degrees.
pub fn detect_feature_edges(mesh: &ProcessMesh, threshold_deg: f64) -> Vec<FeatureEdge> {
    let thresh = threshold_deg.to_radians();
    let face_normals: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|&[a, b, c]| {
            let ab = vec3_sub(mesh.verts[b], mesh.verts[a]);
            let ac = vec3_sub(mesh.verts[c], mesh.verts[a]);
            vec3_normalize(vec3_cross(ab, ac))
        })
        .collect();
    let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (fi, &[a, b, c]) in mesh.faces.iter().enumerate() {
        for (i, j) in [(a, b), (b, c), (c, a)] {
            let key = if i < j { (i, j) } else { (j, i) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }
    let mut features = Vec::new();
    for ((v0, v1), flist) in &edge_faces {
        if flist.len() != 2 {
            continue;
        }
        let n0 = face_normals[flist[0]];
        let n1 = face_normals[flist[1]];
        let cos_a = vec3_dot(n0, n1).clamp(-1.0, 1.0);
        let angle = cos_a.acos();
        if angle > thresh {
            features.push(FeatureEdge {
                v0: *v0,
                v1: *v1,
                dihedral_angle: angle,
            });
        }
    }
    features
}
/// Isotropic remeshing targeting a uniform target edge length.
///
/// Performs `iters` rounds of split-long / collapse-short / flip /
/// Laplacian-smooth.
pub fn isotropic_remesh(mesh: &ProcessMesh, target_len: f64, iters: usize) -> ProcessMesh {
    let mut out = mesh.clone();
    for _ in 0..iters {
        out = split_long_edges(&out, target_len * 4.0 / 3.0);
        out = collapse_short_edges(&out, target_len * 4.0 / 5.0);
        out = laplacian_smooth(&out, 0.5, 1);
    }
    out
}
/// Split edges longer than `max_len` by inserting midpoints.
pub fn split_long_edges(mesh: &ProcessMesh, max_len: f64) -> ProcessMesh {
    let mut verts = mesh.verts.clone();
    let mut faces = Vec::new();
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut get_mid = |a: usize, b: usize, v: &mut Vec<Vertex>| -> Option<usize> {
        let d = vec3_len(vec3_sub(v[b], v[a]));
        if d > max_len {
            let key = if a < b { (a, b) } else { (b, a) };
            let idx = *edge_mid.entry(key).or_insert_with(|| {
                let mid = vec3_scale(vec3_add(v[a], v[b]), 0.5);
                let i = v.len();
                v.push(mid);
                i
            });
            Some(idx)
        } else {
            None
        }
    };
    for &[a, b, c] in &mesh.faces {
        let ab = get_mid(a, b, &mut verts);
        let bc = get_mid(b, c, &mut verts);
        let ca = get_mid(c, a, &mut verts);
        match (ab, bc, ca) {
            (None, None, None) => faces.push([a, b, c]),
            (Some(m), None, None) => {
                faces.push([a, m, c]);
                faces.push([m, b, c]);
            }
            (None, Some(m), None) => {
                faces.push([a, b, m]);
                faces.push([a, m, c]);
            }
            (None, None, Some(m)) => {
                faces.push([a, b, m]);
                faces.push([m, b, c]);
            }
            (Some(ab), Some(bc), None) => {
                faces.push([a, ab, c]);
                faces.push([ab, b, bc]);
                faces.push([ab, bc, c]);
            }
            (Some(ab), None, Some(ca)) => {
                faces.push([a, ab, ca]);
                faces.push([ab, b, c]);
                faces.push([ab, c, ca]);
            }
            (None, Some(bc), Some(ca)) => {
                faces.push([a, b, bc]);
                faces.push([a, bc, ca]);
                faces.push([ca, bc, c]);
            }
            (Some(ab), Some(bc), Some(ca)) => {
                faces.push([a, ab, ca]);
                faces.push([ab, b, bc]);
                faces.push([ca, bc, c]);
                faces.push([ab, bc, ca]);
            }
        }
    }
    ProcessMesh::new(verts, faces)
}
/// Collapse edges shorter than `min_len` by merging endpoints.
pub fn collapse_short_edges(mesh: &ProcessMesh, min_len: f64) -> ProcessMesh {
    let nv = mesh.verts.len();
    let mut verts = mesh.verts.clone();
    let mut faces = mesh.faces.clone();
    let mut remap: Vec<usize> = (0..nv).collect();
    fn root(r: &mut Vec<usize>, i: usize) -> usize {
        if r[i] == i {
            i
        } else {
            let p = root(r, r[i]);
            r[i] = p;
            p
        }
    }
    let mut collapsed = true;
    while collapsed {
        collapsed = false;
        let mut seen = HashSet::new();
        for &[a, b, c] in &faces {
            for (i, j) in [(a, b), (b, c), (c, a)] {
                let ri = root(&mut remap, i);
                let rj = root(&mut remap, j);
                if ri == rj {
                    continue;
                }
                let key = if ri < rj { (ri, rj) } else { (rj, ri) };
                if seen.contains(&key) {
                    continue;
                }
                seen.insert(key);
                let d = vec3_len(vec3_sub(verts[rj], verts[ri]));
                if d < min_len {
                    verts[ri] = vec3_scale(vec3_add(verts[ri], verts[rj]), 0.5);
                    remap[rj] = ri;
                    collapsed = true;
                }
            }
        }
        for i in 0..nv {
            let r = root(&mut remap, i);
            remap[i] = r;
        }
        for f in &mut faces {
            for idx in f.iter_mut() {
                *idx = remap[*idx];
            }
        }
        faces.retain(|&[a, b, c]| a != b && b != c && a != c);
    }
    let active: Vec<bool> = (0..nv).map(|i| remap[i] == i).collect();
    let new_idx: Vec<usize> = {
        let mut m = vec![0usize; nv];
        let mut cnt = 0;
        for i in 0..nv {
            if active[i] {
                m[i] = cnt;
                cnt += 1;
            }
        }
        m
    };
    let new_verts: Vec<Vertex> = (0..nv).filter(|&i| active[i]).map(|i| verts[i]).collect();
    let new_faces: Vec<Face> = faces
        .iter()
        .map(|&[a, b, c]| {
            let ra = root(&mut remap, a);
            let rb = root(&mut remap, b);
            let rc = root(&mut remap, c);
            [new_idx[ra], new_idx[rb], new_idx[rc]]
        })
        .collect();
    ProcessMesh::new(new_verts, new_faces)
}
/// Tutte parameterization: map boundary to a circle, solve harmonic interior.
///
/// Uses a simple iterative relaxation (not a direct solver) for portability.
pub fn tutte_parameterize(mesh: &ProcessMesh, iters: usize) -> UvParameterization {
    let nv = mesh.verts.len();
    let loops = detect_boundary_loops(mesh);
    let boundary: Vec<usize> = loops.into_iter().next().unwrap_or_default();
    let nb = boundary.len();
    let mut uv = vec![[0.0_f64; 2]; nv];
    for (k, &vi) in boundary.iter().enumerate() {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / nb as f64;
        uv[vi] = [0.5 + 0.5 * angle.cos(), 0.5 + 0.5 * angle.sin()];
    }
    let boundary_set: HashSet<usize> = boundary.iter().cloned().collect();
    let adj = mesh.build_adjacency();
    for _ in 0..iters {
        let old = uv.clone();
        for i in 0..nv {
            if boundary_set.contains(&i) {
                continue;
            }
            let nbrs = &adj[i];
            if nbrs.is_empty() {
                continue;
            }
            let mut sum = [0.0_f64; 2];
            for &j in nbrs {
                sum[0] += old[j][0];
                sum[1] += old[j][1];
            }
            let n = nbrs.len() as f64;
            uv[i] = [sum[0] / n, sum[1] / n];
        }
    }
    let mut out = mesh.clone();
    out.uvs = Some(uv.clone());
    UvParameterization { mesh: out, uvs: uv }
}
/// Least-Squares Conformal Maps (LSCM) UV parameterization.
///
/// Minimises the conformal (angle-preserving) energy over all triangles by solving
/// the Lévy 2002 sparse linear system.  Two vertices are pinned to break the
/// rotation/translation ambiguity:
/// - vertex 0 → UV (0, 0)
/// - vertex 1 → UV (1, 0)
///
/// The free-vertex system `A_free^T A_free x = A_free^T b` is solved with the
/// [`oxiphysics_fem::parallel_solver::ParallelPcgSolver`].
///
/// # Returns
/// A [`UvParameterization`] with UV coordinates clamped to \[0, 1\] per axis.
pub fn lscm_parameterize(mesh: &ProcessMesh) -> UvParameterization {
    use oxiphysics_fem::parallel_solver::{CsrMatrix, ParallelPcgSolver};
    use std::collections::BTreeMap;

    let nv = mesh.verts.len();
    let nf = mesh.faces.len();

    // Degenerate cases: fall back to Tutte
    if nv < 3 || nf == 0 {
        return tutte_parameterize(mesh, 200);
    }

    // Each triangle contributes 2 rows (Re and Im of the conformal constraint)
    // acting on 6 unknowns [u_i, v_i, u_j, v_j, u_k, v_k].
    // Total rows = 2 * nf, total columns = 2 * nv (u and v interleaved per vertex).
    //
    // The DOF ordering: column index = 2*vertex + 0 for u, 2*vertex + 1 for v.
    //
    // Two pinned DOFs (columns 0,1 for vertex 0 and columns 2,3 for vertex 1)
    // are removed from the system and their contributions moved to the RHS.
    // The pinned values are: u0=0,v0=0, u1=1,v1=0.

    // Pinned vertex indices and values
    let pin_col = [0usize, 1, 2, 3]; // u0, v0, u1, v1
    let pin_val = [0.0_f64, 0.0, 1.0, 0.0];

    // Map from original column (0..2*nv) to free column (0..2*(nv-2))
    let mut col_map = vec![usize::MAX; 2 * nv];
    let mut free_col = 0usize;
    for (c, slot) in col_map.iter_mut().enumerate() {
        if pin_col.contains(&c) {
            continue;
        }
        *slot = free_col;
        free_col += 1;
    }
    let n_free = free_col; // = 2*(nv-2)

    // Accumulate A matrix (as triplets) and rhs b
    // A is (2*nf) × n_free, but we only need A^T A and A^T b for PCG.
    // We'll build it row-by-row and form the normal equations.

    // Triplets for A_free: (row, free_col, value)
    let mut a_triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(12 * nf);
    // rhs[row] accumulates b = -A_pinned * x_pinned
    let mut b_full = vec![0.0_f64; 2 * nf];

    for (fi, &[vi, vj, vk]) in mesh.faces.iter().enumerate() {
        let pi = mesh.verts[vi];
        let pj = mesh.verts[vj];
        let pk = mesh.verts[vk];

        // Local 2-D frame for the triangle
        let d1 = vec3_sub(pj, pi); // p_j - p_i
        let d2 = vec3_sub(pk, pi); // p_k - p_i

        let x1 = vec3_len(d1); // |d1|
        if x1 < 1e-14 {
            continue; // degenerate triangle
        }
        let e1 = vec3_scale(d1, 1.0 / x1); // unit edge direction
        let n_hat = vec3_normalize(vec3_cross(d1, d2));
        let e2 = vec3_cross(n_hat, e1);

        let x2 = vec3_dot(d2, e1);
        let y2 = vec3_dot(d2, e2);
        if y2.abs() < 1e-14 {
            continue; // degenerate triangle
        }

        // Per Lévy 2002, the two LSCM constraint rows for this triangle
        // acting on [u_i, v_i, u_j, v_j, u_k, v_k] are:
        //   Row Re: [x2 - x1, -y2, -x2, y2,  x1,  0 ] / y2
        //   Row Im: [y2, x2 - x1, -y2, -x2,   0, x1 ] / y2
        // (The factor 1/y2 normalises by the triangle "height" in local coords.)

        let inv_y2 = 1.0 / y2;
        let row_re = 2 * fi;
        let row_im = 2 * fi + 1;

        // Local coefficient vectors for each DOF index in the 6-vector
        // local_dofs[0] = global col for u_i, [1] = v_i, [2] = u_j, etc.
        let local_dofs = [2 * vi, 2 * vi + 1, 2 * vj, 2 * vj + 1, 2 * vk, 2 * vk + 1];
        let re_coeffs = [
            (x2 - x1) * inv_y2,
            -1.0,
            -x2 * inv_y2,
            1.0,
            x1 * inv_y2,
            0.0,
        ];
        let im_coeffs = [
            1.0,
            (x2 - x1) * inv_y2,
            -1.0,
            -x2 * inv_y2,
            0.0,
            x1 * inv_y2,
        ];

        for k in 0..6 {
            let gc = local_dofs[k];
            if let Some(pin_pos) = pin_col.iter().position(|&p| p == gc) {
                // Pinned DOF: move to RHS
                b_full[row_re] -= re_coeffs[k] * pin_val[pin_pos];
                b_full[row_im] -= im_coeffs[k] * pin_val[pin_pos];
            } else {
                let fc = col_map[gc];
                a_triplets.push((row_re, fc, re_coeffs[k]));
                a_triplets.push((row_im, fc, im_coeffs[k]));
            }
        }
    }

    if n_free == 0 {
        return tutte_parameterize(mesh, 200);
    }

    // Build the normal equations A^T A x = A^T b
    // A^T A is n_free × n_free symmetric positive semi-definite.
    // Accumulate into a BTreeMap for sparse assembly.
    let mut ata_map: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut atb = vec![0.0_f64; n_free];

    for &(row, fc, val) in &a_triplets {
        atb[fc] += val * b_full[row];
    }

    // Group triplets by row for efficient A^T A accumulation
    // For each pair of nonzeros in the same row: (fc_a, val_a) and (fc_b, val_b)
    // contributes val_a * val_b to ATA[fc_a, fc_b].
    let mut row_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); 2 * nf];
    for &(row, fc, val) in &a_triplets {
        row_entries[row].push((fc, val));
    }
    for entries in &row_entries {
        for &(fc_a, val_a) in entries {
            for &(fc_b, val_b) in entries {
                *ata_map.entry((fc_a, fc_b)).or_insert(0.0) += val_a * val_b;
            }
        }
    }

    // Convert BTreeMap to CSR
    let mut row_offsets = vec![0usize; n_free + 1];
    for &(r, _) in ata_map.keys() {
        row_offsets[r + 1] += 1;
    }
    for i in 0..n_free {
        row_offsets[i + 1] += row_offsets[i];
    }
    let nnz = ata_map.len();
    let mut col_indices_csr = vec![0usize; nnz];
    let mut values_csr = vec![0.0_f64; nnz];
    let mut row_fill = vec![0usize; n_free];
    for (&(r, c), &v) in &ata_map {
        let pos = row_offsets[r] + row_fill[r];
        col_indices_csr[pos] = c;
        values_csr[pos] = v;
        row_fill[r] += 1;
    }

    let ata = CsrMatrix {
        nrows: n_free,
        ncols: n_free,
        row_offsets,
        col_indices: col_indices_csr,
        values: values_csr,
    };

    // Solve A^T A x = A^T b  using parallel PCG
    let solver = ParallelPcgSolver::new(2000, 1e-8);
    let mut x_free = vec![0.0_f64; n_free];
    let _stats = solver.solve(&ata, &atb, &mut x_free);

    // Reconstruct full UV array
    let mut uv = vec![[0.0_f64; 2]; nv];
    // Pinned vertices
    uv[0] = [0.0, 0.0];
    uv[1] = [1.0, 0.0];
    // Free vertices
    for (v, uv_v) in uv.iter_mut().enumerate() {
        let cu = 2 * v;
        let cv = 2 * v + 1;
        let u_val = if pin_col.contains(&cu) {
            let pos = pin_col.iter().position(|&p| p == cu).unwrap_or(0);
            pin_val[pos]
        } else {
            x_free[col_map[cu]]
        };
        let v_val = if pin_col.contains(&cv) {
            let pos = pin_col.iter().position(|&p| p == cv).unwrap_or(1);
            pin_val[pos]
        } else {
            x_free[col_map[cv]]
        };
        *uv_v = [u_val.clamp(0.0, 1.0), v_val.clamp(0.0, 1.0)];
    }

    let mut out = mesh.clone();
    out.uvs = Some(uv.clone());
    UvParameterization { mesh: out, uvs: uv }
}
/// Generate a simple texture atlas by partitioning the mesh into connected
/// components and packing their UV patches into a square atlas.
pub fn generate_texture_atlas(mesh: &ProcessMesh, atlas_size: usize) -> Vec<AtlasPatch> {
    let _ = atlas_size;
    let param = tutte_parameterize(mesh, 100);
    let uvs = param.uvs.clone();
    let bounds = uvs.iter().fold(
        [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
        |mut b, &[u, v]| {
            if u < b[0] {
                b[0] = u;
            }
            if v < b[1] {
                b[1] = v;
            }
            if u > b[2] {
                b[2] = u;
            }
            if v > b[3] {
                b[3] = v;
            }
            b
        },
    );
    vec![AtlasPatch {
        group_id: 0,
        bounds,
        uvs,
    }]
}
/// Test whether a point `p` is inside `mesh` using ray casting.
pub fn point_in_mesh(p: [f64; 3], mesh: &ProcessMesh) -> bool {
    let dir = vec3_normalize([1.0, 1e-4, 1e-4]);
    let mut crossings = 0usize;
    for &[a, b, c] in &mesh.faces {
        if let Some(_t) =
            ray_triangle_intersect(p, dir, mesh.verts[a], mesh.verts[b], mesh.verts[c])
        {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}
/// Möller-Trumbore ray-triangle intersection.
pub(super) fn ray_triangle_intersect(
    orig: [f64; 3],
    dir: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<f64> {
    let e1 = vec3_sub(v1, v0);
    let e2 = vec3_sub(v2, v0);
    let h = vec3_cross(dir, e2);
    let a = vec3_dot(e1, h);
    if a.abs() < 1e-15 {
        return None;
    }
    let f = 1.0 / a;
    let s = vec3_sub(orig, v0);
    let u = f * vec3_dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = vec3_cross(s, e1);
    let v = f * vec3_dot(dir, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * vec3_dot(e2, q);
    if t > 1e-15 { Some(t) } else { None }
}
/// Boolean union of two meshes (retains all faces from both).
///
/// This is a surface-union approximation: faces of mesh A that lie inside B
/// are removed, and vice-versa.
pub fn mesh_union(a: &ProcessMesh, b: &ProcessMesh) -> BooleanResult {
    let offset_b = a.verts.len();
    let mut verts = a.verts.clone();
    verts.extend_from_slice(&b.verts);
    let faces_a: Vec<Face> = a
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(a.verts[i], a.verts[j]), a.verts[k]),
                1.0 / 3.0,
            );
            !point_in_mesh(centroid, b)
        })
        .cloned()
        .collect();
    let faces_b: Vec<Face> = b
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(b.verts[i], b.verts[j]), b.verts[k]),
                1.0 / 3.0,
            );
            !point_in_mesh(centroid, a)
        })
        .map(|&[i, j, k]| [i + offset_b, j + offset_b, k + offset_b])
        .collect();
    let mut faces = faces_a;
    faces.extend(faces_b);
    let mut result = BooleanResult {
        mesh: ProcessMesh::new(verts, faces),
        is_exact: false,
    };
    result.is_exact = result.is_topologically_exact();
    result
}
