//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::io::{self, Write};

use super::types::{KdTree, PointCloud, Triangle};

/// Voxel map type: maps (ix, iy, iz) cell coordinates to (point indices, centroid accumulator, count).
type VoxelMap = std::collections::HashMap<(i64, i64, i64), (Vec<usize>, [f64; 3], usize)>;

/// Maximum leaf size for the KD-tree.
pub(super) const KD_LEAF_SIZE: usize = 16;
/// Squared distance between two 3D points.
pub(super) fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}
/// Euclidean distance between two 3D points.
pub(super) fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    dist_sq(a, b).sqrt()
}
/// Normalize a 3D vector. Returns `[0; 3]` if length < eps.
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
/// Cross product of two 3D vectors.
#[cfg(test)]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Downsample a point cloud using a voxel grid filter.
///
/// Each voxel cell of side `voxel_size` keeps only the centroid of points
/// falling within it.
pub fn voxel_grid_downsample(cloud: &PointCloud, voxel_size: f64) -> PointCloud {
    if cloud.is_empty() || voxel_size <= 0.0 {
        return cloud.clone();
    }
    let inv_size = 1.0 / voxel_size;
    let mut voxels: VoxelMap = VoxelMap::new();
    for (i, p) in cloud.positions.iter().enumerate() {
        let vx = (p[0] * inv_size).floor() as i64;
        let vy = (p[1] * inv_size).floor() as i64;
        let vz = (p[2] * inv_size).floor() as i64;
        let key = (vx, vy, vz);
        let entry = voxels
            .entry(key)
            .or_insert_with(|| (Vec::new(), [0.0; 3], 0));
        entry.0.push(i);
        entry.1[0] += p[0];
        entry.1[1] += p[1];
        entry.1[2] += p[2];
        entry.2 += 1;
    }
    let mut result = PointCloud::new();
    let has_normals = cloud.has_normals();
    let has_colors = cloud.has_colors();
    let has_intensities = cloud.has_intensities();
    if has_normals {
        result.normals = Some(Vec::new());
    }
    if has_colors {
        result.colors = Some(Vec::new());
    }
    if has_intensities {
        result.intensities = Some(Vec::new());
    }
    for (indices, sum, count) in voxels.values() {
        let n = *count as f64;
        let centroid = [sum[0] / n, sum[1] / n, sum[2] / n];
        result.positions.push(centroid);
        if has_normals {
            let normals = cloud.normals.as_ref().expect("value should be Some");
            let mut avg_n = [0.0; 3];
            for &i in indices {
                avg_n[0] += normals[i][0];
                avg_n[1] += normals[i][1];
                avg_n[2] += normals[i][2];
            }
            result
                .normals
                .as_mut()
                .expect("value should be Some")
                .push(normalize3(avg_n));
        }
        if has_colors {
            let colors = cloud.colors.as_ref().expect("value should be Some");
            let mut avg_c = [0.0_f64; 3];
            for &i in indices {
                avg_c[0] += colors[i][0] as f64;
                avg_c[1] += colors[i][1] as f64;
                avg_c[2] += colors[i][2] as f64;
            }
            result.colors.as_mut().expect("value should be Some").push([
                (avg_c[0] / n).round().clamp(0.0, 255.0) as u8,
                (avg_c[1] / n).round().clamp(0.0, 255.0) as u8,
                (avg_c[2] / n).round().clamp(0.0, 255.0) as u8,
            ]);
        }
        if has_intensities {
            let ints = cloud.intensities.as_ref().expect("value should be Some");
            let avg_i: f64 = indices.iter().map(|&i| ints[i]).sum::<f64>() / n;
            result
                .intensities
                .as_mut()
                .expect("value should be Some")
                .push(avg_i);
        }
    }
    result
}
/// Remove statistical outliers from a point cloud.
///
/// For each point, computes the mean distance to its `k` nearest neighbors.
/// Points whose mean distance is more than `std_dev_multiplier` standard
/// deviations above the global mean are removed.
pub fn statistical_outlier_removal(
    cloud: &PointCloud,
    k: usize,
    std_dev_multiplier: f64,
) -> PointCloud {
    if cloud.len() <= k {
        return cloud.clone();
    }
    let tree = KdTree::build(cloud);
    let mut mean_dists = Vec::with_capacity(cloud.len());
    for (i, pos) in cloud.positions.iter().enumerate() {
        let neighbors = tree.k_nearest(*pos, k + 1);
        let sum: f64 = neighbors
            .iter()
            .filter(|(idx, _)| *idx != i)
            .take(k)
            .map(|(_, d_sq)| d_sq.sqrt())
            .sum();
        mean_dists.push(sum / k as f64);
    }
    let global_mean: f64 = mean_dists.iter().sum::<f64>() / mean_dists.len() as f64;
    let variance: f64 = mean_dists
        .iter()
        .map(|d| (d - global_mean).powi(2))
        .sum::<f64>()
        / mean_dists.len() as f64;
    let std_dev = variance.sqrt();
    let threshold = global_mean + std_dev_multiplier * std_dev;
    let mut result = PointCloud::new();
    if cloud.has_normals() {
        result.normals = Some(Vec::new());
    }
    if cloud.has_colors() {
        result.colors = Some(Vec::new());
    }
    if cloud.has_intensities() {
        result.intensities = Some(Vec::new());
    }
    for (i, &md) in mean_dists.iter().enumerate() {
        if md <= threshold {
            result.positions.push(cloud.positions[i]);
            if let Some(ref normals) = cloud.normals {
                result
                    .normals
                    .as_mut()
                    .expect("value should be Some")
                    .push(normals[i]);
            }
            if let Some(ref colors) = cloud.colors {
                result
                    .colors
                    .as_mut()
                    .expect("value should be Some")
                    .push(colors[i]);
            }
            if let Some(ref ints) = cloud.intensities {
                result
                    .intensities
                    .as_mut()
                    .expect("value should be Some")
                    .push(ints[i]);
            }
        }
    }
    result
}
/// Estimate normals for a point cloud using k-nearest neighbors.
///
/// For each point, finds `k` nearest neighbors, computes the covariance
/// matrix, and takes the eigenvector with smallest eigenvalue as the normal.
/// Normals are oriented consistently by checking alignment with the centroid
/// direction.
pub fn estimate_normals(cloud: &PointCloud, k: usize) -> Vec<[f64; 3]> {
    if cloud.is_empty() {
        return Vec::new();
    }
    let tree = KdTree::build(cloud);
    let mut normals = Vec::with_capacity(cloud.len());
    let centroid = compute_centroid(cloud);
    for (i, pos) in cloud.positions.iter().enumerate() {
        let neighbors = tree.k_nearest(*pos, k + 1);
        let neighbor_positions: Vec<[f64; 3]> = neighbors
            .iter()
            .filter(|(idx, _)| *idx != i)
            .take(k)
            .map(|(idx, _)| cloud.positions[*idx])
            .collect();
        if neighbor_positions.len() < 3 {
            normals.push([0.0, 0.0, 1.0]);
            continue;
        }
        let n_neigh = neighbor_positions.len() as f64;
        let mut mean = [0.0; 3];
        for np in &neighbor_positions {
            mean[0] += np[0];
            mean[1] += np[1];
            mean[2] += np[2];
        }
        mean[0] /= n_neigh;
        mean[1] /= n_neigh;
        mean[2] /= n_neigh;
        let mut cov = [[0.0_f64; 3]; 3];
        for np in &neighbor_positions {
            let d = [np[0] - mean[0], np[1] - mean[1], np[2] - mean[2]];
            for (r, cov_row) in cov.iter_mut().enumerate() {
                for (c, cell) in cov_row.iter_mut().enumerate() {
                    *cell += d[r] * d[c];
                }
            }
        }
        for cov_row in cov.iter_mut() {
            for cell in cov_row.iter_mut() {
                *cell /= n_neigh;
            }
        }
        let normal = smallest_eigenvector_3x3(cov);
        let to_centroid = [
            centroid[0] - pos[0],
            centroid[1] - pos[1],
            centroid[2] - pos[2],
        ];
        let dot_val =
            normal[0] * to_centroid[0] + normal[1] * to_centroid[1] + normal[2] * to_centroid[2];
        if dot_val > 0.0 {
            normals.push([-normal[0], -normal[1], -normal[2]]);
        } else {
            normals.push(normal);
        }
    }
    normals
}
/// Compute the centroid of a point cloud.
pub fn compute_centroid(cloud: &PointCloud) -> [f64; 3] {
    if cloud.is_empty() {
        return [0.0; 3];
    }
    let n = cloud.len() as f64;
    let mut c = [0.0; 3];
    for p in &cloud.positions {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}
/// Approximate smallest eigenvector of a 3x3 symmetric matrix using
/// inverse power iteration.
fn smallest_eigenvector_3x3(m: [[f64; 3]; 3]) -> [f64; 3] {
    let mut v = [1.0, 0.0, 0.0];
    for _ in 0..20 {
        let mv = mat_vec_3x3(m, v);
        let len = (mv[0] * mv[0] + mv[1] * mv[1] + mv[2] * mv[2]).sqrt();
        if len < 1e-15 {
            return v;
        }
        v = [mv[0] / len, mv[1] / len, mv[2] / len];
    }
    let max_eval =
        v[0] * mat_vec_3x3(m, v)[0] + v[1] * mat_vec_3x3(m, v)[1] + v[2] * mat_vec_3x3(m, v)[2];
    let shift = max_eval + 1.0;
    let b = [
        [shift - m[0][0], -m[0][1], -m[0][2]],
        [-m[1][0], shift - m[1][1], -m[1][2]],
        [-m[2][0], -m[2][1], shift - m[2][2]],
    ];
    let mut w = [0.577, 0.577, 0.577];
    for _ in 0..30 {
        let bw = mat_vec_3x3(b, w);
        let len = (bw[0] * bw[0] + bw[1] * bw[1] + bw[2] * bw[2]).sqrt();
        if len < 1e-15 {
            return w;
        }
        w = [bw[0] / len, bw[1] / len, bw[2] / len];
    }
    w
}
/// Multiply 3x3 matrix by 3-vector.
fn mat_vec_3x3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// Merge two point clouds into one.
pub fn merge_point_clouds(a: &PointCloud, b: &PointCloud) -> PointCloud {
    let mut result = PointCloud::new();
    result.positions = a
        .positions
        .iter()
        .chain(b.positions.iter())
        .copied()
        .collect();
    if a.has_normals() || b.has_normals() {
        let mut normals = Vec::with_capacity(a.len() + b.len());
        if let Some(ref na) = a.normals {
            normals.extend_from_slice(na);
        } else {
            normals.extend(std::iter::repeat_n([0.0_f64; 3], a.len()));
        }
        if let Some(ref nb) = b.normals {
            normals.extend_from_slice(nb);
        } else {
            normals.extend(std::iter::repeat_n([0.0_f64; 3], b.len()));
        }
        result.normals = Some(normals);
    }
    if a.has_colors() || b.has_colors() {
        let mut colors = Vec::with_capacity(a.len() + b.len());
        if let Some(ref ca) = a.colors {
            colors.extend_from_slice(ca);
        } else {
            colors.extend(std::iter::repeat_n([255_u8; 3], a.len()));
        }
        if let Some(ref cb) = b.colors {
            colors.extend_from_slice(cb);
        } else {
            colors.extend(std::iter::repeat_n([255_u8; 3], b.len()));
        }
        result.colors = Some(colors);
    }
    if a.has_intensities() || b.has_intensities() {
        let mut ints = Vec::with_capacity(a.len() + b.len());
        if let Some(ref ia) = a.intensities {
            ints.extend_from_slice(ia);
        } else {
            ints.extend(std::iter::repeat_n(0.0_f64, a.len()));
        }
        if let Some(ref ib) = b.intensities {
            ints.extend_from_slice(ib);
        } else {
            ints.extend(std::iter::repeat_n(0.0_f64, b.len()));
        }
        result.intensities = Some(ints);
    }
    result
}
/// Apply a rigid transformation to a point cloud.
///
/// `rotation` is a 3x3 matrix (row-major), `translation` is `[tx, ty, tz]`.
pub fn transform_point_cloud(
    cloud: &PointCloud,
    rotation: [[f64; 3]; 3],
    translation: [f64; 3],
) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        let rotated = mat_vec_3x3(rotation, *pos);
        *pos = [
            rotated[0] + translation[0],
            rotated[1] + translation[1],
            rotated[2] + translation[2],
        ];
    }
    if let Some(ref mut normals) = result.normals {
        for n in normals.iter_mut() {
            *n = normalize3(mat_vec_3x3(rotation, *n));
        }
    }
    result
}
/// Scale a point cloud uniformly.
pub fn scale_point_cloud(cloud: &PointCloud, scale: f64) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        pos[0] *= scale;
        pos[1] *= scale;
        pos[2] *= scale;
    }
    result
}
/// Translate a point cloud.
pub fn translate_point_cloud(cloud: &PointCloud, offset: [f64; 3]) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        pos[0] += offset[0];
        pos[1] += offset[1];
        pos[2] += offset[2];
    }
    result
}
/// Center a point cloud at the origin.
pub fn center_point_cloud(cloud: &PointCloud) -> PointCloud {
    let c = compute_centroid(cloud);
    translate_point_cloud(cloud, [-c[0], -c[1], -c[2]])
}
/// Write a point cloud to PLY ASCII format.
pub fn write_ply_ascii(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    writeln!(writer, "ply")?;
    writeln!(writer, "format ascii 1.0")?;
    writeln!(writer, "element vertex {}", cloud.len())?;
    writeln!(writer, "property float x")?;
    writeln!(writer, "property float y")?;
    writeln!(writer, "property float z")?;
    if cloud.has_normals() {
        writeln!(writer, "property float nx")?;
        writeln!(writer, "property float ny")?;
        writeln!(writer, "property float nz")?;
    }
    if cloud.has_colors() {
        writeln!(writer, "property uchar red")?;
        writeln!(writer, "property uchar green")?;
        writeln!(writer, "property uchar blue")?;
    }
    if cloud.has_intensities() {
        writeln!(writer, "property float intensity")?;
    }
    writeln!(writer, "end_header")?;
    for i in 0..cloud.len() {
        let p = cloud.positions[i];
        write!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
        if let Some(ref normals) = cloud.normals {
            let n = normals[i];
            write!(writer, " {:.6} {:.6} {:.6}", n[0], n[1], n[2])?;
        }
        if let Some(ref colors) = cloud.colors {
            let c = colors[i];
            write!(writer, " {} {} {}", c[0], c[1], c[2])?;
        }
        if let Some(ref ints) = cloud.intensities {
            write!(writer, " {:.6}", ints[i])?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
/// Write a point cloud to PLY binary little-endian format (to a byte buffer).
pub fn write_ply_binary_le(cloud: &PointCloud) -> Vec<u8> {
    let mut buf = Vec::new();
    let header = format!(
        "ply\nformat binary_little_endian 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nend_header\n",
        cloud.len()
    );
    buf.extend_from_slice(header.as_bytes());
    for p in &cloud.positions {
        buf.extend_from_slice(&(p[0] as f32).to_le_bytes());
        buf.extend_from_slice(&(p[1] as f32).to_le_bytes());
        buf.extend_from_slice(&(p[2] as f32).to_le_bytes());
    }
    buf
}
/// Parse a PLY ASCII string into a point cloud.
pub fn parse_ply_ascii(data: &str) -> Option<PointCloud> {
    let mut lines = data.lines();
    let first = lines.next()?;
    if first.trim() != "ply" {
        return None;
    }
    let mut vertex_count = 0usize;
    let mut has_normals = false;
    let mut has_colors = false;
    let mut _has_intensity = false;
    let mut in_header = true;
    while in_header {
        let line = lines.next()?;
        let line = line.trim();
        if line.starts_with("element vertex") {
            vertex_count = line.split_whitespace().last()?.parse().ok()?;
        } else if line.starts_with("property") {
            if line.contains(" nx") || line.contains(" ny") || line.contains(" nz") {
                has_normals = true;
            }
            if line.contains(" red") || line.contains(" green") || line.contains(" blue") {
                has_colors = true;
            }
            if line.contains(" intensity") {
                _has_intensity = true;
            }
        } else if line == "end_header" {
            in_header = false;
        }
    }
    let mut cloud = PointCloud::new();
    cloud.positions.reserve(vertex_count);
    if has_normals {
        cloud.normals = Some(Vec::with_capacity(vertex_count));
    }
    if has_colors {
        cloud.colors = Some(Vec::with_capacity(vertex_count));
    }
    for _ in 0..vertex_count {
        let line = lines.next()?;
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }
        let x: f64 = vals[0].parse().ok()?;
        let y: f64 = vals[1].parse().ok()?;
        let z: f64 = vals[2].parse().ok()?;
        cloud.positions.push([x, y, z]);
        let mut idx = 3;
        if has_normals && vals.len() >= idx + 3 {
            let nx: f64 = vals[idx].parse().ok()?;
            let ny: f64 = vals[idx + 1].parse().ok()?;
            let nz: f64 = vals[idx + 2].parse().ok()?;
            cloud
                .normals
                .as_mut()
                .expect("value should be Some")
                .push([nx, ny, nz]);
            idx += 3;
        }
        if has_colors && vals.len() >= idx + 3 {
            let r: u8 = vals[idx].parse().ok()?;
            let g: u8 = vals[idx + 1].parse().ok()?;
            let b: u8 = vals[idx + 2].parse().ok()?;
            cloud
                .colors
                .as_mut()
                .expect("value should be Some")
                .push([r, g, b]);
        }
    }
    Some(cloud)
}
/// Write a point cloud to PCD ASCII format.
pub fn write_pcd_ascii(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    writeln!(writer, "# .PCD v0.7 - Point Cloud Data")?;
    writeln!(writer, "VERSION 0.7")?;
    writeln!(writer, "FIELDS x y z")?;
    writeln!(writer, "SIZE 4 4 4")?;
    writeln!(writer, "TYPE F F F")?;
    writeln!(writer, "COUNT 1 1 1")?;
    writeln!(writer, "WIDTH {}", cloud.len())?;
    writeln!(writer, "HEIGHT 1")?;
    writeln!(writer, "VIEWPOINT 0 0 0 1 0 0 0")?;
    writeln!(writer, "POINTS {}", cloud.len())?;
    writeln!(writer, "DATA ascii")?;
    for p in &cloud.positions {
        writeln!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
    }
    Ok(())
}
/// Parse a PCD ASCII string into a point cloud.
pub fn parse_pcd_ascii(data: &str) -> Option<PointCloud> {
    let mut lines = data.lines();
    let mut num_points = 0usize;
    let mut data_started = false;
    let mut cloud = PointCloud::new();
    for line in &mut lines {
        let line = line.trim();
        if line.starts_with("POINTS") {
            num_points = line.split_whitespace().last()?.parse().ok()?;
            cloud.positions.reserve(num_points);
        } else if line.starts_with("DATA") {
            data_started = true;
            break;
        }
    }
    if !data_started {
        return None;
    }
    for _ in 0..num_points {
        let line = lines.next()?;
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }
        let x: f64 = vals[0].parse().ok()?;
        let y: f64 = vals[1].parse().ok()?;
        let z: f64 = vals[2].parse().ok()?;
        cloud.positions.push([x, y, z]);
    }
    Some(cloud)
}
/// Write a point cloud to XYZ format (with optional normals and colors).
pub fn write_xyz(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    for i in 0..cloud.len() {
        let p = cloud.positions[i];
        write!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
        if let Some(ref normals) = cloud.normals {
            let n = normals[i];
            write!(writer, " {:.6} {:.6} {:.6}", n[0], n[1], n[2])?;
        }
        if let Some(ref colors) = cloud.colors {
            let c = colors[i];
            write!(writer, " {} {} {}", c[0], c[1], c[2])?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
/// Parse an XYZ file (space-separated x y z \[nx ny nz\] \[r g b\]).
pub fn parse_xyz(data: &str) -> PointCloud {
    let mut cloud = PointCloud::new();
    let mut first_line_fields = 0;
    for (line_idx, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }
        if line_idx == 0 || first_line_fields == 0 {
            first_line_fields = vals.len();
            if first_line_fields >= 6 {
                cloud.normals = Some(Vec::new());
            }
            if first_line_fields >= 9 {
                cloud.colors = Some(Vec::new());
            }
        }
        let x: f64 = vals[0].parse().unwrap_or(0.0);
        let y: f64 = vals[1].parse().unwrap_or(0.0);
        let z: f64 = vals[2].parse().unwrap_or(0.0);
        cloud.positions.push([x, y, z]);
        if vals.len() >= 6
            && let Some(ref mut normals) = cloud.normals
        {
            let nx: f64 = vals[3].parse().unwrap_or(0.0);
            let ny: f64 = vals[4].parse().unwrap_or(0.0);
            let nz: f64 = vals[5].parse().unwrap_or(0.0);
            normals.push([nx, ny, nz]);
        }
        if vals.len() >= 9
            && let Some(ref mut colors) = cloud.colors
        {
            let r: u8 = vals[6].parse().unwrap_or(0);
            let g: u8 = vals[7].parse().unwrap_or(0);
            let b: u8 = vals[8].parse().unwrap_or(0);
            colors.push([r, g, b]);
        }
    }
    cloud
}
/// Simplified ball-pivoting mesh generation.
///
/// Given a point cloud with normals, attempts to create a triangle mesh by
/// pivoting a ball of radius `ball_radius` over the point cloud surface.
///
/// This is a conceptual implementation suitable for small point clouds.
pub fn ball_pivoting_mesh(cloud: &PointCloud, ball_radius: f64) -> Vec<Triangle> {
    if cloud.len() < 3 {
        return Vec::new();
    }
    let tree = KdTree::build(cloud);
    let mut triangles = Vec::new();
    let mut used = vec![false; cloud.len()];
    let search_radius_sq = (2.0 * ball_radius) * (2.0 * ball_radius);
    let neighbors = tree.k_nearest(cloud.positions[0], 20.min(cloud.len()));
    let mut found_seed = false;
    let mut seed = [0usize; 3];
    'outer: for i in 0..neighbors.len() {
        for j in (i + 1)..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let a = neighbors[i].0;
                let b = neighbors[j].0;
                let c = neighbors[k].0;
                let ab = dist(cloud.positions[a], cloud.positions[b]);
                let bc = dist(cloud.positions[b], cloud.positions[c]);
                let ca = dist(cloud.positions[c], cloud.positions[a]);
                if ab < 2.0 * ball_radius && bc < 2.0 * ball_radius && ca < 2.0 * ball_radius {
                    seed = [a, b, c];
                    found_seed = true;
                    break 'outer;
                }
            }
        }
    }
    if !found_seed {
        return triangles;
    }
    triangles.push(Triangle { indices: seed });
    used[seed[0]] = true;
    used[seed[1]] = true;
    used[seed[2]] = true;
    let max_triangles = cloud.len().min(5000);
    let mut edge_queue: Vec<(usize, usize)> =
        vec![(seed[0], seed[1]), (seed[1], seed[2]), (seed[2], seed[0])];
    while let Some((a, b)) = edge_queue.pop() {
        if triangles.len() >= max_triangles {
            break;
        }
        let mid = [
            0.5 * (cloud.positions[a][0] + cloud.positions[b][0]),
            0.5 * (cloud.positions[a][1] + cloud.positions[b][1]),
            0.5 * (cloud.positions[a][2] + cloud.positions[b][2]),
        ];
        let nearby = tree.radius_search(mid, search_radius_sq);
        for (c, _) in nearby {
            if c == a || c == b || used[c] {
                continue;
            }
            let ac = dist(cloud.positions[a], cloud.positions[c]);
            let bc = dist(cloud.positions[b], cloud.positions[c]);
            if ac < 2.0 * ball_radius && bc < 2.0 * ball_radius {
                triangles.push(Triangle { indices: [a, b, c] });
                used[c] = true;
                edge_queue.push((a, c));
                edge_queue.push((c, b));
                break;
            }
        }
    }
    triangles
}
