// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the AMBER NetCDF trajectory convention writer and reader.

use crate::netcdf::{
    AMBER_VELOCITY_SCALE_FACTOR, CdfVersion, NcType, NcValues, NetcdfTrajectoryBuilder,
    TrajectoryFrame, create_amber_trajectory, read_amber_trajectory, read_cdf, read_cdf_dataset,
    write_cdf,
};

fn sample_builder(nm: bool) -> NetcdfTrajectoryBuilder {
    let mut b = NetcdfTrajectoryBuilder::new().with_title("water box");
    if nm {
        b = b.in_nanometres();
    }
    for i in 0..3 {
        let s = i as f64;
        b.add_frame(TrajectoryFrame {
            time_ps: 0.5 * s,
            positions: vec![[1.0 + s, 2.0, 3.0], [4.0, 5.0 + s, 6.25]],
            velocities: Some(vec![[0.1 * s, -0.2, 0.3], [0.0, 0.5, -1.0]]),
            box_lengths: Some([20.0, 21.0, 22.5]),
        });
    }
    b
}

fn amber_bytes(b: &NetcdfTrajectoryBuilder) -> Vec<u8> {
    let mut buf = Vec::new();
    b.write_amber_netcdf(&mut buf).unwrap();
    buf
}

#[test]
fn builder_writes_a_spec_compliant_64bit_offset_file() {
    let bytes = amber_bytes(&sample_builder(false));
    assert_eq!(&bytes[..4], b"CDF\x02");
    let ds = read_cdf_dataset(&bytes).unwrap();
    assert_eq!(ds.text_attribute("Conventions").as_deref(), Some("AMBER"));
    assert_eq!(
        ds.text_attribute("ConventionVersion").as_deref(),
        Some("1.0")
    );
    assert_eq!(ds.text_attribute("program").as_deref(), Some("OxiPhysics"));
    assert_eq!(
        ds.text_attribute("programVersion").as_deref(),
        Some(env!("CARGO_PKG_VERSION"))
    );
    assert_eq!(ds.text_attribute("title").as_deref(), Some("water box"));
    let frame = ds.dimension("frame").unwrap();
    assert!(frame.unlimited);
    assert_eq!(frame.len, 3);
    assert_eq!(ds.dimension("spatial").map(|d| d.len), Some(3));
    assert_eq!(ds.dimension("label").map(|d| d.len), Some(5));
    assert_eq!(
        ds.variable("spatial").unwrap().values,
        NcValues::text("xyz")
    );
    assert_eq!(
        ds.variable("cell_spatial").unwrap().values,
        NcValues::text("abc")
    );
    assert_eq!(
        ds.variable("cell_angular").unwrap().values,
        NcValues::text("alphabeta gamma")
    );
    for (name, ty, units) in [
        ("coordinates", NcType::Float, "angstrom"),
        ("time", NcType::Float, "picosecond"),
        ("cell_lengths", NcType::Float, "angstrom"),
        ("cell_angles", NcType::Float, "degree"),
        ("velocities", NcType::Float, "angstrom/picosecond"),
    ] {
        let var = ds.variable(name).unwrap();
        assert_eq!(var.nc_type(), ty, "{name}");
        assert_eq!(
            var.text_attribute("units").as_deref(),
            Some(units),
            "{name}"
        );
        assert!(ds.is_record_variable(var), "{name}");
    }
    let cell = ds.variable("cell_lengths").unwrap();
    let dim_names: Vec<&str> = cell
        .dim_ids
        .iter()
        .map(|&id| ds.dimensions[id].name.as_str())
        .collect();
    assert_eq!(dim_names, ["frame", "cell_spatial"]);
    assert_eq!(
        ds.variable("velocities")
            .unwrap()
            .attribute("scale_factor")
            .map(|a| a.values.clone()),
        Some(NcValues::Double(vec![AMBER_VELOCITY_SCALE_FACTOR]))
    );
}

#[test]
fn builder_round_trips_through_the_reader() {
    let builder = sample_builder(false);
    let traj = read_amber_trajectory(&amber_bytes(&builder)).unwrap();
    assert_eq!(traj.n_atoms, 2);
    assert_eq!(traj.n_frames(), 3);
    assert_eq!(traj.title.as_deref(), Some("water box"));
    assert_eq!(traj.time_ps, Some(vec![0.0, 0.5, 1.0]));
    assert_eq!(traj.coordinates[2][0], [3.0, 2.0, 3.0]);
    assert_eq!(traj.coordinates[1][1], [4.0, 6.0, 6.25]);
    assert_eq!(traj.cell_lengths.as_ref().unwrap()[0], [20.0, 21.0, 22.5]);
    assert_eq!(traj.cell_angles.as_ref().unwrap()[2], [90.0, 90.0, 90.0]);
    // Velocities: nm/ps in the builder, Å/ps after scale_factor in the file.
    let vels = traj.velocities.as_ref().unwrap();
    let expected = [[0.0, 5.0, -10.0], [2.0, -2.0, 3.0]];
    for (got, want) in [vels[0][1], vels[2][0]].iter().zip(expected) {
        for k in 0..3 {
            assert!((got[k] - want[k]).abs() < 1e-5, "{got:?} vs {want:?}");
        }
    }
    let frames = traj.to_frames();
    assert_eq!(frames.len(), 3);
    assert_eq!(frames[1].box_lengths, Some([20.0, 21.0, 22.5]));
    let v = frames[2].velocities.as_ref().unwrap()[0];
    assert!((v[0] - 0.2).abs() < 1e-6 && (v[1] + 0.2).abs() < 1e-6);
    let rebuilt = traj.to_builder();
    assert_eq!(rebuilt.frame_count(), 3);
    assert_eq!(rebuilt.frames[0].positions, builder.frames[0].positions);
}

#[test]
fn nanometre_builder_is_converted_to_angstrom() {
    let traj = read_amber_trajectory(&amber_bytes(&sample_builder(true))).unwrap();
    assert_eq!(traj.coordinates[0][0], [10.0, 20.0, 30.0]);
    assert_eq!(
        traj.cell_lengths.as_ref().unwrap()[0],
        [200.0, 210.0, 225.0]
    );
}

#[test]
fn optional_quantities_are_omitted_when_some_frames_lack_them() {
    let mut b = sample_builder(false);
    b.frames[1].velocities = None;
    b.frames[2].box_lengths = None;
    let ds = b.to_amber_dataset().unwrap();
    assert!(ds.variable("velocities").is_none());
    assert!(ds.variable("cell_lengths").is_none());
    assert!(ds.dimension("cell_spatial").is_none());
    let traj = read_amber_trajectory(&write_cdf_bytes(&ds)).unwrap();
    assert!(traj.velocities.is_none());
    assert!(traj.cell_lengths.is_none());
}

fn write_cdf_bytes(ds: &crate::netcdf::CdfDataset) -> Vec<u8> {
    crate::netcdf::write_cdf_dataset(ds, CdfVersion::Offset64).unwrap()
}

#[test]
fn ragged_frames_are_rejected() {
    let mut b = sample_builder(false);
    b.frames[1].positions.pop();
    assert!(b.to_amber_dataset().is_err());
    let mut buf = Vec::new();
    assert!(b.write_amber_netcdf(&mut buf).is_err());
    assert!(buf.is_empty());
}

#[test]
fn create_amber_trajectory_follows_the_convention() {
    let file = create_amber_trajectory(
        2,
        2,
        (0..12).map(f64::from).collect(),
        vec![0.0, 2.0],
        Some(vec![30.0; 6]),
        None,
    );
    assert_eq!(
        file.get_variable("coordinates").unwrap().nc_type,
        NcType::Float
    );
    assert_eq!(
        file.get_variable("spatial").unwrap().as_text().as_deref(),
        Some("xyz")
    );
    assert_eq!(
        file.get_variable("cell_lengths").unwrap().dims,
        vec!["frame".to_string(), "cell_spatial".to_string()]
    );
    // cell_angles accompanies cell_lengths; unknown angles are zero.
    assert_eq!(file.get_variable("cell_angles").unwrap().data, vec![0.0; 6]);
    let bytes = write_cdf(&file, CdfVersion::Offset64).unwrap();
    assert_eq!(&bytes[..4], b"CDF\x02");
    let traj = read_amber_trajectory(&bytes).unwrap();
    assert_eq!(traj.coordinates[1][1], [9.0, 10.0, 11.0]);
    assert_eq!(traj.time_ps, Some(vec![0.0, 2.0]));
    assert_eq!(traj.cell_lengths.as_ref().unwrap()[1], [30.0, 30.0, 30.0]);
    let back = read_cdf(&bytes).unwrap();
    assert_eq!(
        back.get_attribute("programVersion"),
        Some(env!("CARGO_PKG_VERSION"))
    );
    // Mismatched data length is reported when encoding.
    let bad = create_amber_trajectory(2, 2, vec![0.0; 5], vec![0.0, 2.0], None, None);
    assert!(write_cdf(&bad, CdfVersion::Offset64).is_err());
}

#[test]
fn reader_rejects_non_amber_files() {
    let mut file = create_amber_trajectory(1, 1, vec![0.0; 3], vec![0.0], None, None);
    file.attributes.retain(|(k, _)| k != "Conventions");
    file.add_attribute("Conventions", "CF-1.8");
    let bytes = write_cdf(&file, CdfVersion::Offset64).unwrap();
    assert!(read_amber_trajectory(&bytes).is_err());
    assert!(read_amber_trajectory(b"CDF\x02").is_err());
    assert!(read_amber_trajectory(b"not netcdf").is_err());
}

#[test]
fn triclinic_cells_are_not_reported_as_orthorhombic_boxes() {
    let file = create_amber_trajectory(
        1,
        1,
        vec![0.0; 3],
        vec![0.0],
        Some(vec![10.0, 11.0, 12.0]),
        Some(vec![60.0, 90.0, 90.0]),
    );
    let traj = read_amber_trajectory(&write_cdf(&file, CdfVersion::Classic).unwrap()).unwrap();
    assert_eq!(traj.cell_angles.as_ref().unwrap()[0], [60.0, 90.0, 90.0]);
    assert_eq!(traj.to_frames()[0].box_lengths, None);
}
