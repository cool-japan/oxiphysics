//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use crate::xdmf;
use oxiphysics_core::math::Vec3;
use std::io::Write;

/// Compute the velocity autocorrelation function (VACF) from velocity frames.
///
/// VACF(lag) = ⟨ v(t) · v(t + lag) ⟩ / ⟨ v(0) · v(0) ⟩
#[allow(dead_code)]
pub struct VacfCalculator;
#[allow(dead_code)]
impl VacfCalculator {
    /// Compute the normalised VACF up to `max_lag` frames.
    ///
    /// Returns a vector of `(lag, vacf)` pairs.  At lag=0, VACF=1.0 by definition.
    pub fn compute_normalized(frames: &[VelocityFrame], max_lag: usize) -> Vec<(usize, f64)> {
        if frames.is_empty() {
            return Vec::new();
        }
        let effective_lag = max_lag.min(frames.len() - 1);
        let n_atoms = frames[0].n_atoms();
        let c0 = Self::vacf_at_lag(frames, 0, n_atoms);
        if c0.abs() < 1e-30 {
            return (0..=effective_lag).map(|l| (l, 0.0)).collect();
        }
        (0..=effective_lag)
            .map(|lag| {
                let c = Self::vacf_at_lag(frames, lag, n_atoms);
                (lag, c / c0)
            })
            .collect()
    }
    fn vacf_at_lag(frames: &[VelocityFrame], lag: usize, n_atoms: usize) -> f64 {
        let mut sum = 0.0_f64;
        let mut count = 0_usize;
        for t0 in 0..(frames.len().saturating_sub(lag)) {
            let t1 = t0 + lag;
            let n = n_atoms.min(frames[t0].n_atoms()).min(frames[t1].n_atoms());
            for i in 0..n {
                let v0 = frames[t0].velocities[i];
                let v1 = frames[t1].velocities[i];
                sum += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                count += 1;
            }
        }
        if count > 0 { sum / count as f64 } else { 0.0 }
    }
    /// Integrate the VACF using the trapezoidal rule to estimate the
    /// self-diffusion coefficient.  `dt` is the time-step.
    ///
    /// D = (1/3) ∫ C(t) dt  (unnormalized VACF)
    pub fn integrate_vacf(vacf: &[(usize, f64)], dt: f64) -> f64 {
        if vacf.len() < 2 || dt <= 0.0 {
            return 0.0;
        }
        let mut integral = 0.0_f64;
        for i in 0..vacf.len() - 1 {
            integral += (vacf[i].1 + vacf[i + 1].1) * 0.5 * dt;
        }
        integral / 3.0
    }
}
/// Reads particle trajectories from XYZ-format strings.
#[allow(dead_code)]
pub struct XyzReader;
#[allow(dead_code)]
impl XyzReader {
    /// Parse a single XYZ frame from the start of `data`.
    pub fn read_frame(data: &str) -> Result<TrajectoryFrame, crate::Error> {
        let mut lines = data.lines();
        let count_line = lines
            .next()
            .ok_or_else(|| crate::Error::Parse("XYZ: empty input".to_string()))?;
        let n_atoms: usize = count_line
            .trim()
            .parse()
            .map_err(|_| crate::Error::Parse(format!("XYZ: bad atom count '{}'", count_line)))?;
        let comment = lines
            .next()
            .ok_or_else(|| crate::Error::Parse("XYZ: missing comment line".to_string()))?;
        let (timestep, time) = parse_xyz_comment(comment);
        let mut positions = Vec::with_capacity(n_atoms);
        let mut atom_types = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            let line = lines.next().ok_or_else(|| {
                crate::Error::Parse(format!("XYZ: missing atom line {} of {}", i + 1, n_atoms))
            })?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(crate::Error::Parse(format!(
                    "XYZ: atom line {} has only {} fields",
                    i,
                    parts.len()
                )));
            }
            atom_types.push(parts[0].to_string());
            let x: f64 = parts[1]
                .parse()
                .map_err(|_| crate::Error::Parse(format!("XYZ: bad x coord on atom line {}", i)))?;
            let y: f64 = parts[2]
                .parse()
                .map_err(|_| crate::Error::Parse(format!("XYZ: bad y coord on atom line {}", i)))?;
            let z: f64 = parts[3]
                .parse()
                .map_err(|_| crate::Error::Parse(format!("XYZ: bad z coord on atom line {}", i)))?;
            positions.push([x, y, z]);
        }
        Ok(TrajectoryFrame {
            timestep,
            time,
            positions,
            atom_types,
        })
    }
    /// Parse all frames from a concatenated XYZ string.
    pub fn read_all_frames(data: &str) -> Result<Vec<TrajectoryFrame>, crate::Error> {
        let mut frames = Vec::new();
        let mut remaining = data;
        loop {
            remaining = remaining.trim_start_matches('\n');
            remaining = remaining.trim_start_matches('\r');
            if remaining.is_empty() {
                break;
            }
            let mut it = remaining.lines();
            let count_line = match it.next() {
                Some(l) => l,
                None => break,
            };
            let n_atoms: usize = count_line.trim().parse().map_err(|_| {
                crate::Error::Parse(format!("XYZ: bad atom count '{}'", count_line))
            })?;
            let frame_lines = 2 + n_atoms;
            let frame_text: String = remaining
                .lines()
                .take(frame_lines)
                .collect::<Vec<&str>>()
                .join("\n");
            let frame = Self::read_frame(&frame_text)?;
            frames.push(frame);
            let consumed = count_line_bytes(remaining, frame_lines);
            remaining = &remaining[consumed..];
        }
        Ok(frames)
    }
}
/// Compute statistics over trajectory frames.
#[allow(dead_code)]
pub struct TrajectoryStatistics;
#[allow(dead_code)]
impl TrajectoryStatistics {
    /// Compute the mean position of each atom across all frames.
    pub fn mean_positions(frames: &[TrajectoryFrame]) -> Vec<[f64; 3]> {
        if frames.is_empty() {
            return Vec::new();
        }
        let n_atoms = frames[0].n_atoms();
        let mut mean = vec![[0.0_f64; 3]; n_atoms];
        let mut count = vec![0_usize; n_atoms];
        for frame in frames {
            let n = n_atoms.min(frame.n_atoms());
            for i in 0..n {
                mean[i][0] += frame.positions[i][0];
                mean[i][1] += frame.positions[i][1];
                mean[i][2] += frame.positions[i][2];
                count[i] += 1;
            }
        }
        for i in 0..n_atoms {
            if count[i] > 0 {
                let c = count[i] as f64;
                mean[i][0] /= c;
                mean[i][1] /= c;
                mean[i][2] /= c;
            }
        }
        mean
    }
    /// Compute the per-atom RMSF (root mean square fluctuation).
    pub fn rmsf(frames: &[TrajectoryFrame]) -> Vec<f64> {
        if frames.is_empty() {
            return Vec::new();
        }
        let mean = Self::mean_positions(frames);
        let n_atoms = mean.len();
        let mut sum_sq = vec![0.0_f64; n_atoms];
        let mut count = vec![0_usize; n_atoms];
        for frame in frames {
            let n = n_atoms.min(frame.n_atoms());
            for i in 0..n {
                let dx = frame.positions[i][0] - mean[i][0];
                let dy = frame.positions[i][1] - mean[i][1];
                let dz = frame.positions[i][2] - mean[i][2];
                sum_sq[i] += dx * dx + dy * dy + dz * dz;
                count[i] += 1;
            }
        }
        sum_sq
            .iter()
            .zip(count.iter())
            .map(|(&sq, &c)| if c > 0 { (sq / c as f64).sqrt() } else { 0.0 })
            .collect()
    }
    /// Compute the radius of gyration for a frame.
    pub fn radius_of_gyration(frame: &TrajectoryFrame, masses: &[f64]) -> f64 {
        let n = frame.n_atoms().min(masses.len());
        if n == 0 {
            return 0.0;
        }
        let com = center_of_mass(&frame.positions[..n], &masses[..n]);
        let total_mass: f64 = masses[..n].iter().sum();
        if total_mass < 1e-30 {
            return 0.0;
        }
        let mut sum = 0.0_f64;
        for i in 0..n {
            let dx = frame.positions[i][0] - com[0];
            let dy = frame.positions[i][1] - com[1];
            let dz = frame.positions[i][2] - com[2];
            sum += masses[i] * (dx * dx + dy * dy + dz * dz);
        }
        (sum / total_mass).sqrt()
    }
    /// Compute the end-to-end distance of a chain (first to last atom).
    pub fn end_to_end_distance(frame: &TrajectoryFrame) -> f64 {
        if frame.n_atoms() < 2 {
            return 0.0;
        }
        let first = frame.positions[0];
        let last = frame.positions[frame.n_atoms() - 1];
        let dx = last[0] - first[0];
        let dy = last[1] - first[1];
        let dz = last[2] - first[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Compute the radial distribution function g(r) from one or more frames.
#[allow(dead_code)]
pub struct RdfCalculator;
#[allow(dead_code)]
impl RdfCalculator {
    /// Compute g(r) for all atom pairs within `r_max`, using `n_bins` bins.
    ///
    /// `box_size` is required for the normalisation (number density).
    /// Returns `(r_values, g_values)` both of length `n_bins`.
    pub fn compute(
        frames: &[TrajectoryFrame],
        r_max: f64,
        n_bins: usize,
        box_size: [f64; 3],
    ) -> (Vec<f64>, Vec<f64>) {
        if frames.is_empty() || n_bins == 0 || r_max <= 0.0 {
            return (Vec::new(), Vec::new());
        }
        let dr = r_max / n_bins as f64;
        let mut hist = vec![0_u64; n_bins];
        let mut n_frames = 0_u64;
        let mut n_atoms_total = 0_u64;
        for frame in frames {
            let n = frame.n_atoms();
            if n < 2 {
                continue;
            }
            n_frames += 1;
            n_atoms_total += n as u64;
            for i in 0..n {
                for j in (i + 1)..n {
                    let mut dx = frame.positions[j][0] - frame.positions[i][0];
                    let mut dy = frame.positions[j][1] - frame.positions[i][1];
                    let mut dz = frame.positions[j][2] - frame.positions[i][2];
                    for (d, l) in [
                        (&mut dx, box_size[0]),
                        (&mut dy, box_size[1]),
                        (&mut dz, box_size[2]),
                    ] {
                        if l > 1e-30 {
                            *d -= (*d / l).round() * l;
                        }
                    }
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    if r < r_max {
                        let bin = (r / dr) as usize;
                        if bin < n_bins {
                            hist[bin] += 2;
                        }
                    }
                }
            }
        }
        let vol = box_size[0] * box_size[1] * box_size[2];
        let avg_n = if n_frames > 0 {
            n_atoms_total as f64 / n_frames as f64
        } else {
            1.0
        };
        let rho = avg_n / vol;
        let r_vals: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * dr).collect();
        let g_vals: Vec<f64> = (0..n_bins)
            .map(|b| {
                let r = r_vals[b];
                let shell_vol = 4.0 * std::f64::consts::PI * r * r * dr;
                let ideal = rho * shell_vol * avg_n;
                if ideal < 1e-30 || n_frames == 0 {
                    0.0
                } else {
                    hist[b] as f64 / (n_frames as f64 * ideal)
                }
            })
            .collect();
        (r_vals, g_vals)
    }
}
/// Compute the Mean Square Displacement (MSD) as a function of lag time
/// from a trajectory.
///
/// Returns a vector of `(lag_index, msd_value)` pairs.
/// `lag_index` runs from 0 to `max_lag` (inclusive); at lag 0 the MSD is 0.
#[allow(dead_code)]
pub struct MsdCalculator;
#[allow(dead_code)]
impl MsdCalculator {
    /// Compute MSD for all atoms combined.
    ///
    /// `max_lag` is the maximum lag (number of frames) to compute.
    pub fn compute(frames: &[TrajectoryFrame], max_lag: usize) -> Vec<(usize, f64)> {
        if frames.len() < 2 || max_lag == 0 {
            return vec![(0, 0.0)];
        }
        let n_atoms = frames[0].n_atoms();
        let effective_lag = max_lag.min(frames.len() - 1);
        let mut result = Vec::with_capacity(effective_lag + 1);
        for lag in 0..=effective_lag {
            let mut sum = 0.0_f64;
            let mut count = 0_usize;
            for t0 in 0..(frames.len() - lag) {
                let t1 = t0 + lag;
                let n = n_atoms.min(frames[t0].n_atoms()).min(frames[t1].n_atoms());
                for i in 0..n {
                    let dx = frames[t1].positions[i][0] - frames[t0].positions[i][0];
                    let dy = frames[t1].positions[i][1] - frames[t0].positions[i][1];
                    let dz = frames[t1].positions[i][2] - frames[t0].positions[i][2];
                    sum += dx * dx + dy * dy + dz * dz;
                    count += 1;
                }
            }
            let msd = if count > 0 { sum / count as f64 } else { 0.0 };
            result.push((lag, msd));
        }
        result
    }
    /// Estimate the self-diffusion coefficient D from MSD data using
    /// the Einstein relation: `MSD = 6 D t` in 3-D.
    ///
    /// Performs a linear fit to `msd_data` (pairs of `(lag, msd)`),
    /// using `dt` as the frame time-step.  Returns `D` (in units of
    /// distance²/time).
    pub fn diffusion_coefficient(msd_data: &[(usize, f64)], dt: f64) -> f64 {
        if msd_data.len() < 2 || dt <= 0.0 {
            return 0.0;
        }
        let n = msd_data.len() as f64;
        let sum_x: f64 = msd_data.iter().map(|(lag, _)| *lag as f64 * dt).sum();
        let sum_y: f64 = msd_data.iter().map(|(_, msd)| *msd).sum();
        let sum_xy: f64 = msd_data
            .iter()
            .map(|(lag, msd)| *lag as f64 * dt * msd)
            .sum();
        let sum_x2: f64 = msd_data
            .iter()
            .map(|(lag, _)| (*lag as f64 * dt).powi(2))
            .sum();
        let denom = n * sum_x2 - sum_x * sum_x;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        let slope = (n * sum_xy - sum_x * sum_y) / denom;
        slope / 6.0
    }
}
/// Analyse bond lengths from trajectory data.
#[allow(dead_code)]
pub struct BondLengthAnalyser;
#[allow(dead_code)]
impl BondLengthAnalyser {
    /// Compute the distance between atoms `i` and `j` in a frame.
    pub fn bond_length(frame: &TrajectoryFrame, i: usize, j: usize) -> f64 {
        let pi = frame.positions[i];
        let pj = frame.positions[j];
        let dx = pi[0] - pj[0];
        let dy = pi[1] - pj[1];
        let dz = pi[2] - pj[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Compute bond length over all frames and return the time series.
    pub fn time_series(frames: &[TrajectoryFrame], i: usize, j: usize) -> Vec<f64> {
        frames
            .iter()
            .filter(|f| f.n_atoms() > i.max(j))
            .map(|f| Self::bond_length(f, i, j))
            .collect()
    }
    /// Compute average bond length and standard deviation over all frames.
    pub fn mean_and_std(frames: &[TrajectoryFrame], i: usize, j: usize) -> (f64, f64) {
        let series = Self::time_series(frames, i, j);
        if series.is_empty() {
            return (0.0, 0.0);
        }
        let mean = series.iter().sum::<f64>() / series.len() as f64;
        let var =
            series.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / series.len() as f64;
        (mean, var.sqrt())
    }
    /// Return the bond angle (in degrees) for atoms i–j–k in a frame.
    ///
    /// The angle is at atom `j` (the middle atom).
    pub fn bond_angle_deg(frame: &TrajectoryFrame, i: usize, j: usize, k: usize) -> f64 {
        let pi = frame.positions[i];
        let pj = frame.positions[j];
        let pk = frame.positions[k];
        let v1 = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
        let v2 = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let l1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let l2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if l1 < 1e-30 || l2 < 1e-30 {
            return 0.0;
        }
        let cos_theta = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (l1 * l2);
        cos_theta.clamp(-1.0, 1.0).acos().to_degrees()
    }
}
/// Resample a trajectory to a different time resolution.
#[allow(dead_code)]
pub struct TrajectoryResampler;
#[allow(dead_code)]
impl TrajectoryResampler {
    /// Resample frames at uniform time intervals.
    ///
    /// Linearly interpolates positions between the two nearest frames
    /// for each target time. Returns the resampled frames.
    pub fn resample_uniform(frames: &[TrajectoryFrame], target_dt: f64) -> Vec<TrajectoryFrame> {
        if frames.is_empty() || target_dt <= 0.0 {
            return Vec::new();
        }
        if frames.len() == 1 {
            return vec![frames[0].clone()];
        }
        let t_start = frames[0].time;
        let t_end = frames.last().expect("collection should not be empty").time;
        if t_end <= t_start {
            return vec![frames[0].clone()];
        }
        let n_samples = ((t_end - t_start) / target_dt).ceil() as usize + 1;
        let mut result = Vec::with_capacity(n_samples);
        for s in 0..n_samples {
            let t = t_start + s as f64 * target_dt;
            if t > t_end {
                break;
            }
            let frame = Self::interpolate_at_time(frames, t);
            result.push(frame);
        }
        result
    }
    /// Interpolate a frame at a specific time.
    fn interpolate_at_time(frames: &[TrajectoryFrame], t: f64) -> TrajectoryFrame {
        let mut lo = 0;
        let mut hi = frames.len() - 1;
        for i in 0..frames.len() - 1 {
            if frames[i].time <= t && frames[i + 1].time >= t {
                lo = i;
                hi = i + 1;
                break;
            }
        }
        if lo == hi || (frames[hi].time - frames[lo].time).abs() < 1e-30 {
            return frames[lo].clone();
        }
        let alpha = (t - frames[lo].time) / (frames[hi].time - frames[lo].time);
        let alpha = alpha.clamp(0.0, 1.0);
        let n_atoms = frames[lo].n_atoms().min(frames[hi].n_atoms());
        let mut positions = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            positions.push([
                frames[lo].positions[i][0] * (1.0 - alpha) + frames[hi].positions[i][0] * alpha,
                frames[lo].positions[i][1] * (1.0 - alpha) + frames[hi].positions[i][1] * alpha,
                frames[lo].positions[i][2] * (1.0 - alpha) + frames[hi].positions[i][2] * alpha,
            ]);
        }
        TrajectoryFrame {
            timestep: (t / 1.0) as u64,
            time: t,
            positions,
            atom_types: frames[lo].atom_types[..n_atoms].to_vec(),
        }
    }
    /// Subsample: keep every n-th frame.
    pub fn subsample(frames: &[TrajectoryFrame], every_n: usize) -> Vec<TrajectoryFrame> {
        if every_n == 0 {
            return Vec::new();
        }
        frames.iter().step_by(every_n).cloned().collect()
    }
}
/// One snapshot in a particle trajectory.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TrajectoryFrame {
    /// Integer timestep index.
    pub timestep: u64,
    /// Simulation time (e.g. in ps or fs).
    pub time: f64,
    /// Atom positions as `[x, y, z]` triples.
    pub positions: Vec<[f64; 3]>,
    /// Atom type labels (e.g. element symbols like `"C"`, `"H"`, `"O"`).
    pub atom_types: Vec<String>,
}
#[allow(dead_code)]
impl TrajectoryFrame {
    /// Create a new `TrajectoryFrame`.
    pub fn new(
        timestep: u64,
        time: f64,
        positions: Vec<[f64; 3]>,
        atom_types: Vec<String>,
    ) -> Self {
        Self {
            timestep,
            time,
            positions,
            atom_types,
        }
    }
    /// Number of atoms in this frame.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }
}
/// Writes particle trajectories in the standard XYZ file format.
///
/// Format per frame:
/// ```text
/// `N`
/// <comment line>
/// `type` `x` `y` `z`
/// ...
/// ```
#[allow(dead_code)]
pub struct XyzWriter;
#[allow(dead_code)]
impl XyzWriter {
    /// Render one frame as an XYZ-format string.
    pub fn write_frame(frame: &TrajectoryFrame) -> String {
        let mut out = String::new();
        out.push_str(&format!("{}\n", frame.n_atoms()));
        out.push_str(&format!(
            "Timestep {} time={}\n",
            frame.timestep, frame.time
        ));
        for (pos, atom_type) in frame.positions.iter().zip(frame.atom_types.iter()) {
            out.push_str(&format!("{} {} {} {}\n", atom_type, pos[0], pos[1], pos[2]));
        }
        out
    }
    /// Render multiple frames as a concatenated XYZ string.
    pub fn write_frames(frames: &[TrajectoryFrame]) -> String {
        frames.iter().map(Self::write_frame).collect()
    }
}
/// Writes `TrajectoryFrame` data in LAMMPS dump format.
///
/// Note: this is distinct from `LammpsDumpWriter` in `lammps_dump.rs`, which
/// operates on `LammpsDumpFrame` (with integer atom IDs and type indices).
/// This writer uses the string atom types from `TrajectoryFrame`.
#[allow(dead_code)]
pub struct TrajLammpsWriter;
#[allow(dead_code)]
impl TrajLammpsWriter {
    /// Render one frame in LAMMPS dump format.
    ///
    /// `box_lo` and `box_hi` are the simulation box extents as `[xlo, ylo, zlo]`
    /// and `[xhi, yhi, zhi]`.
    pub fn write_dump_frame(frame: &TrajectoryFrame, box_lo: [f64; 3], box_hi: [f64; 3]) -> String {
        let mut out = String::new();
        out.push_str("ITEM: TIMESTEP\n");
        out.push_str(&format!("{}\n", frame.timestep));
        out.push_str("ITEM: NUMBER OF ATOMS\n");
        out.push_str(&format!("{}\n", frame.n_atoms()));
        out.push_str("ITEM: BOX BOUNDS pp pp pp\n");
        for i in 0..3 {
            out.push_str(&format!("{} {}\n", box_lo[i], box_hi[i]));
        }
        out.push_str("ITEM: ATOMS id type x y z\n");
        for (i, (pos, atom_type)) in frame
            .positions
            .iter()
            .zip(frame.atom_types.iter())
            .enumerate()
        {
            out.push_str(&format!(
                "{} {} {} {} {}\n",
                i + 1,
                atom_type,
                pos[0],
                pos[1],
                pos[2]
            ));
        }
        out
    }
}
/// Convert trajectory data between different representations.
#[allow(dead_code)]
pub struct TrajectoryConverter;
#[allow(dead_code)]
impl TrajectoryConverter {
    /// Convert a `TrajectoryFrame` to a flat `[x0,y0,z0, x1,y1,z1, ...]` array.
    pub fn to_flat_xyz(frame: &TrajectoryFrame) -> Vec<f64> {
        let mut out = Vec::with_capacity(frame.n_atoms() * 3);
        for pos in &frame.positions {
            out.extend_from_slice(pos);
        }
        out
    }
    /// Build a `TrajectoryFrame` from a flat xyz array.
    pub fn from_flat_xyz(
        flat: &[f64],
        atom_types: Vec<String>,
        timestep: u64,
        time: f64,
    ) -> std::result::Result<TrajectoryFrame, crate::Error> {
        if !flat.len().is_multiple_of(3) {
            return Err(crate::Error::Parse(format!(
                "from_flat_xyz: length {} not divisible by 3",
                flat.len()
            )));
        }
        let n_atoms = flat.len() / 3;
        if atom_types.len() != n_atoms {
            return Err(crate::Error::Parse(format!(
                "from_flat_xyz: {} atom types but {} positions",
                atom_types.len(),
                n_atoms
            )));
        }
        let positions: Vec<[f64; 3]> = flat.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        Ok(TrajectoryFrame::new(timestep, time, positions, atom_types))
    }
    /// Compute the displacement between two frames (atom-by-atom).
    ///
    /// Returns displacement vectors `[[dx, dy, dz\], ...]` for each atom.
    pub fn frame_displacement(
        frame_a: &TrajectoryFrame,
        frame_b: &TrajectoryFrame,
    ) -> Vec<[f64; 3]> {
        let n = frame_a.n_atoms().min(frame_b.n_atoms());
        (0..n)
            .map(|i| {
                [
                    frame_b.positions[i][0] - frame_a.positions[i][0],
                    frame_b.positions[i][1] - frame_a.positions[i][1],
                    frame_b.positions[i][2] - frame_a.positions[i][2],
                ]
            })
            .collect()
    }
    /// Translate all atoms in a frame by vector `delta`.
    pub fn translate(frame: &TrajectoryFrame, delta: [f64; 3]) -> TrajectoryFrame {
        let positions = frame
            .positions
            .iter()
            .map(|p| [p[0] + delta[0], p[1] + delta[1], p[2] + delta[2]])
            .collect();
        TrajectoryFrame::new(
            frame.timestep,
            frame.time,
            positions,
            frame.atom_types.clone(),
        )
    }
    /// Scale all atom positions by a scalar factor.
    pub fn scale(frame: &TrajectoryFrame, factor: f64) -> TrajectoryFrame {
        let positions = frame
            .positions
            .iter()
            .map(|p| [p[0] * factor, p[1] * factor, p[2] * factor])
            .collect();
        TrajectoryFrame::new(
            frame.timestep,
            frame.time,
            positions,
            frame.atom_types.clone(),
        )
    }
}
/// Concatenate multiple trajectories into one.
#[allow(dead_code)]
pub struct TrajectoryConcatenator;
#[allow(dead_code)]
impl TrajectoryConcatenator {
    /// Concatenate trajectories, adjusting timestamps to be continuous.
    pub fn concatenate(trajectories: &[Vec<TrajectoryFrame>]) -> Vec<TrajectoryFrame> {
        let mut result = Vec::new();
        let mut time_offset = 0.0_f64;
        let mut step_offset = 0_u64;
        for traj in trajectories {
            for frame in traj {
                let mut new_frame = frame.clone();
                new_frame.time += time_offset;
                new_frame.timestep += step_offset;
                result.push(new_frame);
            }
            if let Some(last) = traj.last() {
                time_offset += last.time;
                step_offset += last.timestep + 1;
            }
        }
        result
    }
    /// Merge trajectories by interleaving based on time.
    pub fn merge_sorted(trajectories: &[Vec<TrajectoryFrame>]) -> Vec<TrajectoryFrame> {
        let mut all: Vec<&TrajectoryFrame> = trajectories.iter().flat_map(|t| t.iter()).collect();
        all.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        all.into_iter().cloned().collect()
    }
}
/// Accumulates simulation frames (time + particle positions) and writes
/// animation/trajectory files in multiple formats.
pub struct TrajectoryWriter {
    pub(super) frames: Vec<(f64, Vec<Vec3>)>,
}
impl TrajectoryWriter {
    /// Create a new empty `TrajectoryWriter`.
    pub fn new() -> Self {
        Self { frames: Vec::new() }
    }
    /// Append a frame to the trajectory.
    ///
    /// `time` is the simulation time of this frame; `positions` are the
    /// particle positions at that time.
    pub fn add_frame(&mut self, time: f64, positions: &[Vec3]) {
        self.frames.push((time, positions.to_vec()));
    }
    /// Return the number of frames currently stored.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }
    /// Write all frames as an XDMF temporal collection.
    pub fn write_xdmf<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let steps: Vec<(f64, &[Vec3])> = self
            .frames
            .iter()
            .map(|(t, positions)| (*t, positions.as_slice()))
            .collect();
        xdmf::write_xdmf_temporal(writer, &steps)
    }
    /// Write all frames in simple XYZ trajectory format.
    ///
    /// Each frame is written as:
    /// ```text
    /// `N`
    /// Frame `i` time=`t`
    /// X  x  y  z
    /// X  x  y  z
    /// ...
    /// ```
    pub fn write_xyz<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        for (i, (time, positions)) in self.frames.iter().enumerate() {
            writeln!(writer, "{}", positions.len())?;
            writeln!(writer, "Frame {} time={}", i, time)?;
            for p in positions {
                writeln!(writer, "X  {}  {}  {}", p.x, p.y, p.z)?;
            }
        }
        Ok(())
    }
}
/// High-level analysis operations on multi-frame trajectories.
#[allow(dead_code)]
pub struct Trajectory;
#[allow(dead_code)]
impl Trajectory {
    /// Compute the per-frame Root-Mean-Square Deviation (RMSD) relative to a
    /// reference frame.
    ///
    /// For each frame *i* the RMSD is:
    ///
    /// ```text
    /// RMSD_i = sqrt( (1/N) * sum_k |r_k^i - r_k^ref|^2 )
    /// ```
    ///
    /// Returns a `Vec`f64` with one entry per input frame.  If `frames` is
    /// empty or `reference` has a different atom count, the corresponding
    /// entry is `f64::NAN`.
    pub fn compute_rmsd_trajectory(
        frames: &[TrajectoryFrame],
        reference: &TrajectoryFrame,
    ) -> Vec<f64> {
        let n_ref = reference.n_atoms();
        frames
            .iter()
            .map(|frame| {
                let n = frame.n_atoms();
                if n != n_ref || n == 0 {
                    return f64::NAN;
                }
                let msd: f64 = frame
                    .positions
                    .iter()
                    .zip(reference.positions.iter())
                    .map(|(p, r)| {
                        let dx = p[0] - r[0];
                        let dy = p[1] - r[1];
                        let dz = p[2] - r[2];
                        dx * dx + dy * dy + dz * dz
                    })
                    .sum::<f64>()
                    / n as f64;
                msd.sqrt()
            })
            .collect()
    }
    /// Compute the radius of gyration for each frame.
    ///
    /// The radius of gyration is defined as:
    ///
    /// ```text
    /// Rg = sqrt( (1/N) * sum_k |r_k - r_com|^2 )
    /// ```
    ///
    /// where `r_com` is the centre of mass (equal mass assumed for all atoms).
    ///
    /// Returns a `Vec`f64` with one value per input frame.
    pub fn compute_radius_of_gyration(frames: &[TrajectoryFrame]) -> Vec<f64> {
        frames
            .iter()
            .map(|frame| {
                let n = frame.n_atoms();
                if n == 0 {
                    return 0.0;
                }
                let inv_n = 1.0 / n as f64;
                let com: [f64; 3] = {
                    let mut s = [0.0_f64; 3];
                    for p in &frame.positions {
                        s[0] += p[0];
                        s[1] += p[1];
                        s[2] += p[2];
                    }
                    [s[0] * inv_n, s[1] * inv_n, s[2] * inv_n]
                };
                let rg2 = frame
                    .positions
                    .iter()
                    .map(|p| {
                        let dx = p[0] - com[0];
                        let dy = p[1] - com[1];
                        let dz = p[2] - com[2];
                        dx * dx + dy * dy + dz * dz
                    })
                    .sum::<f64>()
                    * inv_n;
                rg2.sqrt()
            })
            .collect()
    }
    /// Translate each frame so its centre of mass coincides with that of
    /// `reference`, then apply a *simplified* least-squares rotation using
    /// the Kabsch algorithm (iterative gradient-free 3×3 SVD-free version).
    ///
    /// **Note:** This implementation uses an analytical closed-form
    /// solution for the 3-D case that is exact for rigid bodies when N ≥ 3.
    /// For very small N the result is still the minimum-RMSD superposition.
    ///
    /// Returns a new `Vec<TrajectoryFrame>` with aligned coordinates.
    /// Frames whose atom count differs from `reference` are returned unchanged.
    pub fn align_to_reference(
        frames: &[TrajectoryFrame],
        reference: &TrajectoryFrame,
    ) -> Vec<TrajectoryFrame> {
        let n_ref = reference.n_atoms();
        let ref_com = Self::centre_of_mass(&reference.positions);
        let ref_centred: Vec<[f64; 3]> = reference
            .positions
            .iter()
            .map(|p| [p[0] - ref_com[0], p[1] - ref_com[1], p[2] - ref_com[2]])
            .collect();
        frames
            .iter()
            .map(|frame| {
                if frame.n_atoms() != n_ref {
                    return frame.clone();
                }
                let com = Self::centre_of_mass(&frame.positions);
                let centred: Vec<[f64; 3]> = frame
                    .positions
                    .iter()
                    .map(|p| [p[0] - com[0], p[1] - com[1], p[2] - com[2]])
                    .collect();
                let mut h = [[0.0_f64; 3]; 3];
                for (c, r) in centred.iter().zip(ref_centred.iter()) {
                    for i in 0..3 {
                        for j in 0..3 {
                            h[i][j] += c[i] * r[j];
                        }
                    }
                }
                let rot = polar_rotation_3x3(h);
                let positions: Vec<[f64; 3]> = centred
                    .iter()
                    .map(|c| {
                        let x = rot[0][0] * c[0] + rot[0][1] * c[1] + rot[0][2] * c[2] + ref_com[0];
                        let y = rot[1][0] * c[0] + rot[1][1] * c[1] + rot[1][2] * c[2] + ref_com[1];
                        let z = rot[2][0] * c[0] + rot[2][1] * c[1] + rot[2][2] * c[2] + ref_com[2];
                        [x, y, z]
                    })
                    .collect();
                TrajectoryFrame::new(
                    frame.timestep,
                    frame.time,
                    positions,
                    frame.atom_types.clone(),
                )
            })
            .collect()
    }
    fn centre_of_mass(positions: &[[f64; 3]]) -> [f64; 3] {
        let n = positions.len();
        if n == 0 {
            return [0.0; 3];
        }
        let inv_n = 1.0 / n as f64;
        let mut s = [0.0_f64; 3];
        for p in positions {
            s[0] += p[0];
            s[1] += p[1];
            s[2] += p[2];
        }
        [s[0] * inv_n, s[1] * inv_n, s[2] * inv_n]
    }
}
/// Filter trajectory frames by various criteria.
#[allow(dead_code)]
pub struct TrajectoryFilter;
#[allow(dead_code)]
impl TrajectoryFilter {
    /// Keep only frames within a time range \[t_start, t_end\].
    pub fn time_range(
        frames: &[TrajectoryFrame],
        t_start: f64,
        t_end: f64,
    ) -> Vec<TrajectoryFrame> {
        frames
            .iter()
            .filter(|f| f.time >= t_start && f.time <= t_end)
            .cloned()
            .collect()
    }
    /// Keep only certain atom types in each frame.
    pub fn filter_atom_types(
        frames: &[TrajectoryFrame],
        keep_types: &[&str],
    ) -> Vec<TrajectoryFrame> {
        frames
            .iter()
            .map(|frame| {
                let mut new_positions = Vec::new();
                let mut new_types = Vec::new();
                for (pos, atype) in frame.positions.iter().zip(frame.atom_types.iter()) {
                    if keep_types.iter().any(|&t| t == atype) {
                        new_positions.push(*pos);
                        new_types.push(atype.clone());
                    }
                }
                TrajectoryFrame {
                    timestep: frame.timestep,
                    time: frame.time,
                    positions: new_positions,
                    atom_types: new_types,
                }
            })
            .collect()
    }
    /// Remove frames where no atoms are present.
    pub fn remove_empty(frames: &[TrajectoryFrame]) -> Vec<TrajectoryFrame> {
        frames
            .iter()
            .filter(|f| !f.positions.is_empty())
            .cloned()
            .collect()
    }
}
/// Trajectory frame with per-atom velocity data.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VelocityFrame {
    /// Timestep index.
    pub timestep: u64,
    /// Simulation time.
    pub time: f64,
    /// Per-atom velocities `[vx, vy, vz]`.
    pub velocities: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl VelocityFrame {
    /// Create a new velocity frame.
    pub fn new(timestep: u64, time: f64, velocities: Vec<[f64; 3]>) -> Self {
        Self {
            timestep,
            time,
            velocities,
        }
    }
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.velocities.len()
    }
}
/// Handle periodic boundary conditions in trajectory analysis.
#[allow(dead_code)]
pub struct PeriodicImageHandler;
#[allow(dead_code)]
impl PeriodicImageHandler {
    /// Wrap positions into the primary box \[0, box_size\].
    pub fn wrap_positions(positions: &mut [[f64; 3]], box_size: [f64; 3]) {
        for pos in positions.iter_mut() {
            for d in 0..3 {
                if box_size[d] > 1e-30 {
                    pos[d] = pos[d] - (pos[d] / box_size[d]).floor() * box_size[d];
                }
            }
        }
    }
    /// Compute the minimum image distance between two positions.
    pub fn minimum_image_distance(a: [f64; 3], b: [f64; 3], box_size: [f64; 3]) -> f64 {
        let mut r2 = 0.0_f64;
        for d in 0..3 {
            let mut dx = b[d] - a[d];
            if box_size[d] > 1e-30 {
                dx -= (dx / box_size[d]).round() * box_size[d];
            }
            r2 += dx * dx;
        }
        r2.sqrt()
    }
    /// Unwrap a trajectory to remove periodic jumps.
    ///
    /// Adjusts positions so that each atom moves continuously
    /// (no jumps across box boundaries).
    pub fn unwrap_trajectory(frames: &mut [TrajectoryFrame], box_size: [f64; 3]) {
        if frames.len() < 2 {
            return;
        }
        let n_atoms = frames[0].n_atoms();
        for i in 1..frames.len() {
            let n = n_atoms.min(frames[i].n_atoms());
            for j in 0..n {
                for d in 0..3 {
                    if box_size[d] > 1e-30 {
                        let dx = frames[i].positions[j][d] - frames[i - 1].positions[j][d];
                        let shift = (dx / box_size[d]).round() * box_size[d];
                        frames[i].positions[j][d] -= shift;
                    }
                }
            }
        }
    }
}
