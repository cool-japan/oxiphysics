//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    Atom, Bond, ColorScheme, ContourSegment, Element, MolColor, SecondaryStructure,
};

#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = len3(v);
    if n < 1e-30 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
#[inline]
pub(super) fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
#[inline]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}
pub(super) fn ss_width(ss: SecondaryStructure, helix: f64, strand: f64, loop_w: f64) -> f64 {
    match ss {
        SecondaryStructure::Helix => helix,
        SecondaryStructure::Strand => strand,
        SecondaryStructure::Loop | SecondaryStructure::Turn => loop_w,
    }
}
pub(super) fn blend_mol_color(a: MolColor, b: MolColor, t: f32) -> MolColor {
    MolColor::blend(a, b, t)
}
/// Catmull-Rom spline evaluation at parameter `t` ∈ `[0,1]` between `p1` and `p2`.
pub fn catmull_rom(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], t: f64) -> [f64; 3] {
    let t2 = t * t;
    let t3 = t2 * t;
    let f0 = -0.5 * t3 + t2 - 0.5 * t;
    let f1 = 1.5 * t3 - 2.5 * t2 + 1.0;
    let f2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
    let f3 = 0.5 * t3 - 0.5 * t2;
    [
        f0 * p0[0] + f1 * p1[0] + f2 * p2[0] + f3 * p3[0],
        f0 * p0[1] + f1 * p1[1] + f2 * p2[1] + f3 * p3[1],
        f0 * p0[2] + f1 * p1[2] + f2 * p2[2] + f3 * p3[2],
    ]
}
pub(super) fn is_favoured(phi: f64, psi: f64) -> bool {
    let helix = (-100.0..=-30.0).contains(&phi) && (-70.0..=30.0).contains(&psi);
    let strand = (-180.0..=-40.0).contains(&phi) && (psi >= 90.0 || psi <= -150.0);
    helix || strand
}
pub(super) fn is_allowed(phi: f64, psi: f64) -> bool {
    let _phi_in = (-180.0..=180.0).contains(&phi);
    let _psi_in = (-180.0..=180.0).contains(&psi);
    let steric_clash = (50.0..=180.0).contains(&phi) && (-100.0..=50.0).contains(&psi);
    !steric_clash
}
/// Automatically detect bonds between atoms based on covalent radii.
///
/// Two atoms are considered bonded when their distance is within
/// `tolerance * (r_cov_a + r_cov_b)`.
pub fn detect_bonds(atoms: &[Atom], tolerance: f64) -> Vec<Bond> {
    let mut bonds = Vec::new();
    let n = atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = len3(sub3(atoms[i].position, atoms[j].position));
            let r_sum = (atoms[i].element.covalent_radius() + atoms[j].element.covalent_radius())
                * tolerance;
            if d <= r_sum {
                bonds.push(Bond::single(atoms[i].index, atoms[j].index));
            }
        }
    }
    bonds
}
/// Apply a colour scheme to atoms, returning one colour per atom.
pub fn color_by_scheme(
    atoms: &[Atom],
    scheme: ColorScheme,
    uniform_color: MolColor,
    b_factors: Option<&[f64]>,
) -> Vec<MolColor> {
    match scheme {
        ColorScheme::Cpk => atoms.iter().map(|a| a.element.cpk_color()).collect(),
        ColorScheme::Uniform => vec![uniform_color; atoms.len()],
        ColorScheme::BFactor => {
            if let Some(bfs) = b_factors {
                let max_b = bfs.iter().cloned().fold(0.0f64, f64::max).max(1.0);
                atoms
                    .iter()
                    .enumerate()
                    .map(|(i, _)| {
                        let t = clamp(bfs.get(i).copied().unwrap_or(0.0) / max_b, 0.0, 1.0) as f32;
                        MolColor::new(t, 0.0, 1.0 - t, 1.0)
                    })
                    .collect()
            } else {
                vec![uniform_color; atoms.len()]
            }
        }
        ColorScheme::ResidueRainbow => {
            let n = atoms.len() as f32;
            atoms
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    let t = i as f32 / n;
                    rainbow_color(t)
                })
                .collect()
        }
        ColorScheme::SecondaryStructure => atoms.iter().map(|_a| MolColor::white()).collect(),
    }
}
/// Evaluate a rainbow colour map at `t ∈ [0, 1]`.
pub fn rainbow_color(t: f32) -> MolColor {
    let t = t.clamp(0.0, 1.0);
    let hue = t * 300.0;
    hue_to_rgb(hue)
}
pub(super) fn hue_to_rgb(h: f32) -> MolColor {
    let h = h % 360.0;
    let s = 1.0_f32;
    let v = 1.0_f32;
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r1, g1, b1) = if h < 60.0 {
        (c, x, 0.0)
    } else if h < 120.0 {
        (x, c, 0.0)
    } else if h < 180.0 {
        (0.0, c, x)
    } else if h < 240.0 {
        (0.0, x, c)
    } else if h < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    MolColor::new(r1 + m, g1 + m, b1 + m, 1.0)
}
pub(super) fn element_symbol(e: Element) -> &'static str {
    match e {
        Element::H => "H",
        Element::C => "C",
        Element::N => "N",
        Element::O => "O",
        Element::S => "S",
        Element::P => "P",
        Element::F => "F",
        Element::Cl => "Cl",
        Element::Br => "Br",
        Element::I => "I",
        Element::Fe => "Fe",
        Element::Ca => "Ca",
        Element::Zn => "Zn",
        Element::Mg => "Mg",
        Element::Unknown => "X",
    }
}
/// Extract contour line segments from a 2-D scalar field using the marching-squares algorithm.
///
/// - `field` — flat row-major array of size `nx × ny`.
/// - `nx`, `ny` — grid dimensions.
/// - `x_range`, `y_range` — axis extents `(min, max)`.
/// - `level` — isovalue.
pub fn marching_squares(
    field: &[f64],
    nx: usize,
    ny: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    level: f64,
) -> Vec<ContourSegment> {
    if nx < 2 || ny < 2 || field.len() < nx * ny {
        return vec![];
    }
    let dx = (x_range.1 - x_range.0) / (nx - 1) as f64;
    let dy = (y_range.1 - y_range.0) / (ny - 1) as f64;
    let mut segments = Vec::new();
    let val = |ix: usize, iy: usize| -> f64 { field[iy * nx + ix] };
    let pos = |ix: usize, iy: usize| -> [f64; 2] {
        [x_range.0 + ix as f64 * dx, y_range.0 + iy as f64 * dy]
    };
    let interp_edge = |va: f64, pa: [f64; 2], vb: f64, pb: [f64; 2]| -> [f64; 2] {
        let dv = vb - va;
        let t = if dv.abs() < 1e-30 {
            0.5
        } else {
            (level - va) / dv
        };
        let t = clamp(t, 0.0, 1.0);
        [pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])]
    };
    for iy in 0..ny - 1 {
        for ix in 0..nx - 1 {
            let v00 = val(ix, iy);
            let v10 = val(ix + 1, iy);
            let v11 = val(ix + 1, iy + 1);
            let v01 = val(ix, iy + 1);
            let p00 = pos(ix, iy);
            let p10 = pos(ix + 1, iy);
            let p11 = pos(ix + 1, iy + 1);
            let p01 = pos(ix, iy + 1);
            let mask = ((v00 >= level) as u8)
                | (((v10 >= level) as u8) << 1)
                | (((v11 >= level) as u8) << 2)
                | (((v01 >= level) as u8) << 3);
            let e_bot = interp_edge(v00, p00, v10, p10);
            let e_rig = interp_edge(v10, p10, v11, p11);
            let e_top = interp_edge(v01, p01, v11, p11);
            let e_lft = interp_edge(v00, p00, v01, p01);
            let seg = |a: [f64; 2], b: [f64; 2]| ContourSegment {
                p0: a,
                p1: b,
                level,
            };
            match mask {
                1 | 14 => segments.push(seg(e_bot, e_lft)),
                2 | 13 => segments.push(seg(e_bot, e_rig)),
                3 | 12 => segments.push(seg(e_rig, e_lft)),
                4 | 11 => segments.push(seg(e_top, e_rig)),
                5 => {
                    segments.push(seg(e_bot, e_top));
                    segments.push(seg(e_rig, e_lft));
                }
                6 | 9 => segments.push(seg(e_bot, e_top)),
                7 | 8 => segments.push(seg(e_top, e_lft)),
                10 => {
                    segments.push(seg(e_bot, e_lft));
                    segments.push(seg(e_top, e_rig));
                }
                _ => {}
            }
        }
    }
    segments
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_viz::IsoTriangle;
    use crate::molecular_viz::AtomFractional;
    use crate::molecular_viz::AtomSelector;
    use crate::molecular_viz::BallStickScene;
    use crate::molecular_viz::BondOrder;
    use crate::molecular_viz::CpkScene;
    use crate::molecular_viz::CrystalCell;
    use crate::molecular_viz::DihedralPoint;
    use crate::molecular_viz::ElectronDensityGrid;
    use crate::molecular_viz::EnergyLandscape;
    use crate::molecular_viz::LatticeParameters;
    use crate::molecular_viz::MolecularOrbital;
    use crate::molecular_viz::MolecularScene;
    use crate::molecular_viz::Molecule;
    use crate::molecular_viz::OrbitalSymmetry;
    use crate::molecular_viz::RadialDistribution;
    use crate::molecular_viz::RamachandranData;
    use crate::molecular_viz::RibbonControlPoint;
    use crate::molecular_viz::RibbonDiagram;
    use crate::molecular_viz::TrajectoryAnimation;
    use crate::molecular_viz::TrajectoryFrame;
    #[test]
    fn test_element_vdw_carbon() {
        let r = Element::C.vdw_radius();
        assert!((r - 1.70).abs() < 1e-10);
    }
    #[test]
    fn test_element_covalent_hydrogen() {
        let r = Element::H.covalent_radius();
        assert!(r < 0.5, "H covalent radius should be < 0.5 Å");
    }
    #[test]
    fn test_element_from_symbol_roundtrip() {
        assert_eq!(Element::from_symbol("C"), Element::C);
        assert_eq!(Element::from_symbol("N"), Element::N);
        assert_eq!(Element::from_symbol("Cl"), Element::Cl);
        assert_eq!(Element::from_symbol("??"), Element::Unknown);
    }
    #[test]
    fn test_element_cpk_color_oxygen_red() {
        let c = Element::O.cpk_color();
        assert!(c.r > 0.8 && c.g < 0.3, "oxygen should be reddish");
    }
    #[test]
    fn test_bond_single_constructor() {
        let b = Bond::single(0, 1);
        assert_eq!(b.atom_a, 0);
        assert_eq!(b.atom_b, 1);
        assert_eq!(b.order, BondOrder::Single);
    }
    #[test]
    fn test_ball_stick_water() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::O),
            Atom::new(1, [0.96, 0.0, 0.0], Element::H),
            Atom::new(2, [-0.25, 0.93, 0.0], Element::H),
        ];
        let bonds = vec![Bond::single(0, 1), Bond::single(0, 2)];
        let scene = BallStickScene::build(&atoms, &bonds, 0.4, 0.10);
        assert_eq!(scene.spheres.len(), 3);
        assert_eq!(scene.cylinders.len(), 4);
    }
    #[test]
    fn test_ball_stick_aabb_non_empty() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [3.0, 0.0, 0.0], Element::C),
        ];
        let scene = BallStickScene::build(&atoms, &[], 0.4, 0.1);
        let (mn, mx) = scene.aabb().unwrap();
        assert!(mn[0] < 0.0);
        assert!(mx[0] > 3.0 - 1e-10);
    }
    #[test]
    fn test_ball_stick_no_atoms_aabb_none() {
        let scene = BallStickScene::build(&[], &[], 0.4, 0.1);
        assert!(scene.aabb().is_none());
    }
    #[test]
    fn test_cpk_radius_oxygen() {
        let atoms = vec![Atom::new(0, [0.0, 0.0, 0.0], Element::O)];
        let scene = CpkScene::build(&atoms, 1.0);
        assert!((scene.spheres[0].radius - Element::O.vdw_radius()).abs() < 1e-10);
    }
    #[test]
    fn test_cpk_volume_positive() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [5.0, 0.0, 0.0], Element::O),
        ];
        let scene = CpkScene::build(&atoms, 1.0);
        assert!(scene.total_volume() > 0.0);
    }
    #[test]
    fn test_cpk_vdw_scale() {
        let atoms = vec![Atom::new(0, [0.0, 0.0, 0.0], Element::N)];
        let scene1 = CpkScene::build(&atoms, 1.0);
        let scene2 = CpkScene::build(&atoms, 2.0);
        assert!((scene2.spheres[0].radius - 2.0 * scene1.spheres[0].radius).abs() < 1e-10);
    }
    #[test]
    fn test_density_grid_index_consistency() {
        let mut g = ElectronDensityGrid::new(4, 5, 6, 0.5, [0.0; 3]);
        g.set(2, 3, 4, 1.5);
        assert!((g.get(2, 3, 4) - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_density_grid_position() {
        let g = ElectronDensityGrid::new(3, 3, 3, 1.0, [1.0, 2.0, 3.0]);
        let p = g.position(1, 1, 1);
        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!((p[1] - 3.0).abs() < 1e-10);
        assert!((p[2] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_density_grid_out_of_bounds_returns_zero() {
        let g = ElectronDensityGrid::new(2, 2, 2, 1.0, [0.0; 3]);
        assert_eq!(g.get(99, 99, 99), 0.0);
    }
    #[test]
    fn test_density_grid_max_min() {
        let mut g = ElectronDensityGrid::new(2, 2, 2, 1.0, [0.0; 3]);
        g.set(0, 0, 0, -1.0);
        g.set(1, 1, 1, 5.0);
        assert!((g.max_value() - 5.0).abs() < 1e-10);
        assert!((g.min_value() - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_density_grid_slice_z() {
        let mut g = ElectronDensityGrid::new(2, 2, 2, 1.0, [0.0; 3]);
        g.set(0, 0, 1, 3.0);
        let s = g.slice_z(1);
        assert_eq!(s.len(), 4);
        assert!(s.contains(&3.0));
    }
    #[test]
    fn test_density_grid_interpolate_corner() {
        let mut g = ElectronDensityGrid::new(2, 2, 2, 1.0, [0.0; 3]);
        g.set(0, 0, 0, 8.0);
        let v = g.interpolate([0.0, 0.0, 0.0]);
        assert!((v - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_catmull_rom_t0_equals_p1() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [2.0, 0.0, 0.0];
        let p3 = [3.0, 0.0, 0.0];
        let pt = catmull_rom(p0, p1, p2, p3, 0.0);
        assert!((pt[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_catmull_rom_t1_equals_p2() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [2.0, 0.0, 0.0];
        let p3 = [3.0, 0.0, 0.0];
        let pt = catmull_rom(p0, p1, p2, p3, 1.0);
        assert!((pt[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_ribbon_empty_input() {
        let r = RibbonDiagram::build(&[], 4, 1.5, 2.5, 0.5);
        assert!(r.segments.is_empty());
    }
    #[test]
    fn test_ribbon_two_residues_produces_segments() {
        let pts = vec![
            RibbonControlPoint {
                ca_position: [0.0, 0.0, 0.0],
                ribbon_normal: [0.0, 1.0, 0.0],
                ss: SecondaryStructure::Helix,
                residue_seq: 1,
            },
            RibbonControlPoint {
                ca_position: [3.8, 0.0, 0.0],
                ribbon_normal: [0.0, 1.0, 0.0],
                ss: SecondaryStructure::Helix,
                residue_seq: 2,
            },
        ];
        let r = RibbonDiagram::build(&pts, 4, 1.5, 2.5, 0.5);
        assert!(!r.segments.is_empty());
    }
    #[test]
    fn test_mo_density_at_zero_occupation() {
        let mo = MolecularOrbital::new(
            "HOMO",
            -5.0,
            0.0,
            OrbitalSymmetry::Sigma,
            2,
            2,
            2,
            1.0,
            [0.0; 3],
        );
        assert_eq!(mo.density_at(0, 0, 0), 0.0);
    }
    #[test]
    fn test_mo_phase_seeds_empty_at_high_threshold() {
        let mo = MolecularOrbital::new(
            "LUMO",
            -3.0,
            2.0,
            OrbitalSymmetry::Pi,
            2,
            2,
            2,
            1.0,
            [0.0; 3],
        );
        let (pos, neg) = mo.phase_seeds(100.0);
        assert!(pos.is_empty());
        assert!(neg.is_empty());
    }
    #[test]
    fn test_trajectory_advance_wraps() {
        let mut traj = TrajectoryAnimation::new(3, 25.0);
        for i in 0..3 {
            traj.push_frame(TrajectoryFrame {
                time_ps: i as f64,
                positions: vec![[0.0; 3]; 3],
                velocities: None,
                potential_energy: None,
                kinetic_energy: None,
                temperature: None,
            });
        }
        traj.advance();
        traj.advance();
        let wrapped = traj.advance();
        assert!(wrapped);
        assert_eq!(traj.current_frame, 0);
    }
    #[test]
    fn test_trajectory_rmsd_same_is_zero() {
        let pos = vec![[1.0, 2.0, 3.0]; 5];
        let frame = TrajectoryFrame {
            time_ps: 0.0,
            positions: pos.clone(),
            velocities: None,
            potential_energy: Some(-100.0),
            kinetic_energy: Some(50.0),
            temperature: Some(300.0),
        };
        assert!((frame.rmsd(&pos)).abs() < 1e-10);
    }
    #[test]
    fn test_trajectory_total_energy() {
        let frame = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![],
            velocities: None,
            potential_energy: Some(-200.0),
            kinetic_energy: Some(75.0),
            temperature: None,
        };
        assert!((frame.total_energy().unwrap() - (-125.0)).abs() < 1e-10);
    }
    #[test]
    fn test_rdf_empty_gives_zero_g() {
        let rdf = RadialDistribution::compute(&[], &[], &[], 5.0, 50, 0.0335, "O-H");
        assert!(rdf.g_r.iter().all(|&g| g == 0.0));
    }
    #[test]
    fn test_rdf_bin_count_matches_n_bins() {
        let rdf = RadialDistribution::compute(&[], &[], &[], 10.0, 100, 0.033, "X-X");
        assert_eq!(rdf.r.len(), 100);
        assert_eq!(rdf.g_r.len(), 100);
    }
    #[test]
    fn test_rdf_first_peak_r_some() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        let rdf = RadialDistribution::compute(&positions, &[0], &[1], 5.0, 20, 0.033, "C-C");
        assert!(rdf.first_peak_r().is_some());
    }
    #[test]
    fn test_energy_landscape_minimum_zero() {
        let cv1: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let cv2: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let pes = EnergyLandscape::from_samples(
            &cv1,
            &cv2,
            (0.0, 9.0),
            (0.0, 9.0),
            5,
            5,
            300.0,
            "x",
            "y",
        );
        if let Some((emin, _, _)) = pes.minimum() {
            assert!(emin.abs() < 1e-6, "min energy should be 0 after shifting");
        }
    }
    #[test]
    fn test_energy_landscape_contour_levels_count() {
        let pes = EnergyLandscape::from_samples(
            &[0.0],
            &[0.0],
            (-1.0, 1.0),
            (-1.0, 1.0),
            2,
            2,
            300.0,
            "a",
            "b",
        );
        let levels = pes.contour_levels(5, 20.0);
        assert_eq!(levels.len(), 5);
    }
    #[test]
    fn test_ramachandran_favoured_fraction_helical() {
        let pts = vec![
            DihedralPoint {
                phi: -60.0,
                psi: -45.0,
                residue_seq: 1,
                ss: SecondaryStructure::Helix,
            },
            DihedralPoint {
                phi: -65.0,
                psi: -40.0,
                residue_seq: 2,
                ss: SecondaryStructure::Helix,
            },
        ];
        let r = RamachandranData::new(pts, 36);
        assert!(
            r.favoured_fraction() > 0.9,
            "helical residues should be in favoured region"
        );
    }
    #[test]
    fn test_ramachandran_by_ss() {
        let pts = vec![
            DihedralPoint {
                phi: -60.0,
                psi: -45.0,
                residue_seq: 1,
                ss: SecondaryStructure::Helix,
            },
            DihedralPoint {
                phi: -120.0,
                psi: 130.0,
                residue_seq: 2,
                ss: SecondaryStructure::Strand,
            },
        ];
        let r = RamachandranData::new(pts, 36);
        assert_eq!(r.by_ss(SecondaryStructure::Helix).len(), 1);
        assert_eq!(r.by_ss(SecondaryStructure::Strand).len(), 1);
        assert_eq!(r.by_ss(SecondaryStructure::Loop).len(), 0);
    }
    #[test]
    fn test_ramachandran_density_map_size() {
        let r = RamachandranData::new(vec![], 36);
        assert_eq!(r.density_map.len(), 36 * 36);
    }
    #[test]
    fn test_cubic_cell_volume() {
        let a = 5.0;
        let params = LatticeParameters::cubic(a);
        let v = params.volume();
        assert!(
            (v - a * a * a).abs() < 1e-8,
            "cubic volume should be a³ = {}",
            a * a * a
        );
    }
    #[test]
    fn test_ortho_cell_wireframe_12_edges() {
        let params = LatticeParameters::orthorhombic(4.0, 5.0, 6.0);
        let cell = CrystalCell::new(params, [0.0; 3], vec![]);
        assert_eq!(cell.wireframe_edges().len(), 12);
    }
    #[test]
    fn test_crystal_supercell_atom_count() {
        let params = LatticeParameters::cubic(4.0);
        let atoms = vec![AtomFractional {
            fa: 0.0,
            fb: 0.0,
            fc: 0.0,
            element: Element::C,
            occupancy: 1.0,
        }];
        let cell = CrystalCell::new(params, [0.0; 3], atoms);
        let sup = cell.supercell_atoms(2, 2, 2);
        assert_eq!(sup.len(), 8);
    }
    #[test]
    fn test_lattice_vectors_cubic_orthogonal() {
        let params = LatticeParameters::cubic(3.0);
        let (va, vb, vc) = params.lattice_vectors();
        assert!(dot3(va, vb).abs() < 1e-8);
        assert!(dot3(va, vc).abs() < 1e-8);
        assert!(dot3(vb, vc).abs() < 1e-8);
    }
    #[test]
    fn test_detect_bonds_close_atoms() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [1.5, 0.0, 0.0], Element::C),
        ];
        let bonds = detect_bonds(&atoms, 1.15);
        assert_eq!(bonds.len(), 1);
    }
    #[test]
    fn test_detect_bonds_far_atoms() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [10.0, 0.0, 0.0], Element::C),
        ];
        let bonds = detect_bonds(&atoms, 1.15);
        assert!(bonds.is_empty());
    }
    #[test]
    fn test_selector_select_all() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [1.0, 0.0, 0.0], Element::N),
        ];
        let sel = AtomSelector::new();
        let indices = sel.select(&atoms);
        assert_eq!(indices.len(), 2);
    }
    #[test]
    fn test_selector_by_element() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [1.0, 0.0, 0.0], Element::N),
            Atom::new(2, [2.0, 0.0, 0.0], Element::C),
        ];
        let sel = AtomSelector::new().by_elements(&[Element::C]);
        let indices = sel.select(&atoms);
        assert_eq!(indices.len(), 2);
        assert!(indices.contains(&0) && indices.contains(&2));
    }
    #[test]
    fn test_selector_within_distance() {
        let atoms = vec![
            Atom::new(0, [0.0, 0.0, 0.0], Element::C),
            Atom::new(1, [1.0, 0.0, 0.0], Element::C),
            Atom::new(2, [10.0, 0.0, 0.0], Element::C),
        ];
        let sel = AtomSelector::new().within([0.0, 0.0, 0.0], 2.0);
        let indices = sel.select(&atoms);
        assert_eq!(indices.len(), 2);
        assert!(!indices.contains(&2));
    }
    #[test]
    fn test_molecule_add_atom_increments_count() {
        let mut mol = Molecule::new("test");
        mol.add_atom([0.0; 3], Element::C);
        mol.add_atom([1.0, 0.0, 0.0], Element::O);
        assert_eq!(mol.atoms.len(), 2);
    }
    #[test]
    fn test_molecule_formula_water() {
        let mut mol = Molecule::new("water");
        mol.add_atom([0.0, 0.0, 0.0], Element::O);
        mol.add_atom([0.96, 0.0, 0.0], Element::H);
        mol.add_atom([-0.25, 0.93, 0.0], Element::H);
        let formula = mol.formula();
        assert!(
            formula.contains("H2"),
            "formula should contain H2: got {}",
            formula
        );
        assert!(
            formula.contains("O1"),
            "formula should contain O1: got {}",
            formula
        );
    }
    #[test]
    fn test_color_by_scheme_cpk_count() {
        let atoms = vec![
            Atom::new(0, [0.0; 3], Element::C),
            Atom::new(1, [1.0, 0.0, 0.0], Element::O),
        ];
        let colors = color_by_scheme(&atoms, ColorScheme::Cpk, MolColor::white(), None);
        assert_eq!(colors.len(), 2);
    }
    #[test]
    fn test_color_by_scheme_uniform() {
        let atoms = vec![
            Atom::new(0, [0.0; 3], Element::C),
            Atom::new(1, [1.0, 0.0, 0.0], Element::N),
        ];
        let u = MolColor::new(0.5, 0.5, 0.5, 1.0);
        let colors = color_by_scheme(&atoms, ColorScheme::Uniform, u, None);
        for c in &colors {
            assert!((c.r - 0.5).abs() < 1e-6);
        }
    }
    #[test]
    fn test_marching_squares_uniform_above_no_contour() {
        let field = vec![2.0f64; 4];
        let segs = marching_squares(&field, 2, 2, (0.0, 1.0), (0.0, 1.0), 1.0);
        assert!(
            segs.is_empty(),
            "uniform-above field should have no contour at level 1.0"
        );
    }
    #[test]
    fn test_marching_squares_step_function_produces_segment() {
        let field = vec![0.0, 2.0, 0.0, 2.0];
        let segs = marching_squares(&field, 2, 2, (0.0, 1.0), (0.0, 1.0), 1.0);
        assert!(
            !segs.is_empty(),
            "step function should produce contour segment"
        );
    }
    #[test]
    fn test_mol_color_blend_midpoint() {
        let a = MolColor::new(0.0, 0.0, 0.0, 1.0);
        let b = MolColor::new(1.0, 1.0, 1.0, 1.0);
        let mid = MolColor::blend(a, b, 0.5);
        assert!((mid.r - 0.5).abs() < 1e-6);
        assert!((mid.g - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_molecular_scene_empty_primitive_count() {
        let s = MolecularScene::empty();
        assert_eq!(s.primitive_count(), 0);
    }
    #[test]
    fn test_molecular_scene_with_cpk() {
        let atoms = vec![Atom::new(0, [0.0; 3], Element::C)];
        let mut s = MolecularScene::empty();
        s.cpk = Some(CpkScene::build(&atoms, 1.0));
        assert_eq!(s.primitive_count(), 1);
    }
    #[test]
    fn test_iso_triangle_face_normal_unit_length() {
        let tri = IsoTriangle {
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normals: [[0.0, 0.0, 1.0]; 3],
            level: 0.5,
        };
        let n = tri.face_normal();
        let l = len3(n);
        assert!(
            (l - 1.0).abs() < 1e-10,
            "face normal should be unit length, got {}",
            l
        );
    }
    #[test]
    fn test_iso_triangle_centroid() {
        let tri = IsoTriangle {
            vertices: [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            normals: [[0.0, 0.0, 1.0]; 3],
            level: 1.0,
        };
        let c = tri.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_normalize3_unit_vector_unchanged() {
        let v = [1.0, 0.0, 0.0];
        let n = normalize3(v);
        assert!((len3(n) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_normalize3_zero_vector() {
        let v = [0.0, 0.0, 0.0];
        let n = normalize3(v);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_cross3_orthogonal() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = cross3(x, y);
        assert!((z[0] - 0.0).abs() < 1e-10);
        assert!((z[1] - 0.0).abs() < 1e-10);
        assert!((z[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rainbow_color_red_at_zero() {
        let c = rainbow_color(0.0);
        assert!(c.r > 0.8, "rainbow at 0 should be reddish");
    }
}
