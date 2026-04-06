// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{DataType, Dataset, ShdfGroup, XdmfParams};

pub(super) const MAGIC: &[u8; 4] = b"SHDF";
pub(super) const VERSION: u32 = 1;
/// Encode a UTF-8 string as a length-prefixed byte sequence.
///
/// Format: `[length: u32 LE][UTF-8 bytes]`
#[allow(dead_code)]
pub fn encode_string(s: &str) -> Vec<u8> {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(4 + bytes.len());
    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
    out.extend_from_slice(bytes);
    out
}
/// Decode a length-prefixed UTF-8 string from `data` starting at `*offset`.
///
/// Advances `*offset` past the consumed bytes.
#[allow(dead_code)]
pub fn decode_string(data: &[u8], offset: &mut usize) -> Result<String, String> {
    let len = read_u32(data, offset)? as usize;
    if *offset + len > data.len() {
        return Err(format!(
            "string data out of bounds: need {} bytes at offset {}",
            len, offset
        ));
    }
    let s = std::str::from_utf8(&data[*offset..*offset + len])
        .map_err(|e| format!("invalid UTF-8 in string: {e}"))?
        .to_string();
    *offset += len;
    Ok(s)
}
pub(super) fn read_u8(data: &[u8], pos: &mut usize) -> Result<u8, String> {
    if *pos >= data.len() {
        return Err(format!("unexpected EOF reading u8 at offset {pos}"));
    }
    let v = data[*pos];
    *pos += 1;
    Ok(v)
}
pub(super) fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, String> {
    require_bytes(data, *pos, 4)?;
    let v = u32::from_le_bytes(
        data[*pos..*pos + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 4;
    Ok(v)
}
pub(super) fn read_u64(data: &[u8], pos: &mut usize) -> Result<u64, String> {
    require_bytes(data, *pos, 8)?;
    let v = u64::from_le_bytes(
        data[*pos..*pos + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 8;
    Ok(v)
}
pub(super) fn read_f32(data: &[u8], pos: &mut usize) -> Result<f32, String> {
    require_bytes(data, *pos, 4)?;
    let v = f32::from_le_bytes(
        data[*pos..*pos + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 4;
    Ok(v)
}
pub(super) fn read_f64(data: &[u8], pos: &mut usize) -> Result<f64, String> {
    require_bytes(data, *pos, 8)?;
    let v = f64::from_le_bytes(
        data[*pos..*pos + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 8;
    Ok(v)
}
pub(super) fn read_i32(data: &[u8], pos: &mut usize) -> Result<i32, String> {
    require_bytes(data, *pos, 4)?;
    let v = i32::from_le_bytes(
        data[*pos..*pos + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 4;
    Ok(v)
}
pub(super) fn read_i64(data: &[u8], pos: &mut usize) -> Result<i64, String> {
    require_bytes(data, *pos, 8)?;
    let v = i64::from_le_bytes(
        data[*pos..*pos + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 8;
    Ok(v)
}
pub(super) fn require_bytes(data: &[u8], pos: usize, n: usize) -> Result<(), String> {
    if pos + n > data.len() {
        Err(format!(
            "unexpected EOF: need {n} bytes at offset {pos}, have {}",
            data.len()
        ))
    } else {
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::hdf5_simple::AttributeHelper;

    use crate::hdf5_simple::ChunkingConfig;

    use crate::hdf5_simple::CompressionAlgorithm;

    use crate::hdf5_simple::CompressionSettings;

    use crate::hdf5_simple::ShdfFile;
    use crate::hdf5_simple::ShdfSchema;

    use crate::hdf5_simple::types::*;
    #[test]
    fn test_roundtrip_f64() {
        let mut file = ShdfFile::new();
        let original_data: Vec<f64> = vec![1.0, 2.5, -3.15625, 0.0, 1e100];
        file.add_dataset_f64("temperatures", vec![5], original_data.clone());
        let bytes = file.to_bytes();
        let recovered = ShdfFile::from_bytes(&bytes).expect("from_bytes failed");
        let got = recovered
            .get_f64("temperatures")
            .expect("dataset not found");
        assert_eq!(got.len(), original_data.len());
        for (a, b) in original_data.iter().zip(got.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "f64 value mismatch: {a} vs {b}");
        }
    }
    #[test]
    fn test_roundtrip_i32() {
        let mut file = ShdfFile::new();
        let original_data: Vec<i32> = vec![0, 1, -1, i32::MAX, i32::MIN];
        file.add_dataset_i32("indices", vec![5], original_data.clone());
        let bytes = file.to_bytes();
        let recovered = ShdfFile::from_bytes(&bytes).expect("from_bytes failed");
        let got = recovered.get_i32("indices").expect("dataset not found");
        assert_eq!(got, original_data.as_slice());
    }
    #[test]
    fn test_multiple_datasets() {
        let mut file = ShdfFile::new();
        file.add_dataset_f64("velocities", vec![3, 3], vec![0.1; 9]);
        file.add_dataset_i32("labels", vec![3], vec![10, 20, 30]);
        file.add_dataset_f64("pressure", vec![1], vec![101325.0]);
        let bytes = file.to_bytes();
        let recovered = ShdfFile::from_bytes(&bytes).expect("from_bytes failed");
        let vel = recovered.get_f64("velocities").expect("velocities missing");
        assert_eq!(vel.len(), 9);
        let lbl = recovered.get_i32("labels").expect("labels missing");
        assert_eq!(lbl, &[10, 20, 30]);
        let pres = recovered.get_f64("pressure").expect("pressure missing");
        assert!((pres[0] - 101325.0).abs() < 1e-6);
        assert_eq!(recovered.datasets[0].shape, vec![3, 3]);
        assert_eq!(recovered.datasets[1].shape, vec![3]);
        assert_eq!(recovered.datasets[2].shape, vec![1]);
    }
    #[test]
    fn test_empty_file_roundtrip() {
        let file = ShdfFile::new();
        let bytes = file.to_bytes();
        let recovered = ShdfFile::from_bytes(&bytes).expect("from_bytes failed");
        assert!(recovered.datasets.is_empty());
        assert!(recovered.global_attributes.is_empty());
    }
    #[test]
    fn test_attribute_roundtrip() {
        let mut file = ShdfFile::new();
        file.add_global_attr("author", "Team KitaSan");
        file.add_global_attr("date", "2026-03-14");
        file.add_dataset_f64("energy", vec![2], vec![1.0, 2.0]);
        file.datasets[0]
            .attributes
            .push(("units".to_string(), "Joules".to_string()));
        let bytes = file.to_bytes();
        let recovered = ShdfFile::from_bytes(&bytes).expect("from_bytes failed");
        assert_eq!(recovered.global_attributes.len(), 2);
        assert_eq!(
            recovered.global_attributes[0],
            ("author".to_string(), "Team KitaSan".to_string())
        );
        assert_eq!(
            recovered.global_attributes[1],
            ("date".to_string(), "2026-03-14".to_string())
        );
        assert_eq!(recovered.datasets[0].attributes.len(), 1);
        assert_eq!(recovered.datasets[0].attributes[0].0, "units");
        assert_eq!(recovered.datasets[0].attributes[0].1, "Joules");
    }
    #[test]
    fn test_write_to_text() {
        let mut file = ShdfFile::new();
        file.add_global_attr("title", "Demo");
        file.add_dataset_f64("pos", vec![2, 3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let text = file.write_to_text();
        assert!(text.contains("SHDF"));
        assert!(text.contains("title = Demo"));
        assert!(text.contains("[pos]"));
        assert!(text.contains("2×3"));
        assert!(text.contains("Float64"));
    }
    #[test]
    fn test_encode_decode_string() {
        let original = "hello, 世界!";
        let encoded = encode_string(original);
        let mut offset = 0usize;
        let decoded = decode_string(&encoded, &mut offset).expect("decode failed");
        assert_eq!(decoded, original);
        assert_eq!(offset, encoded.len());
    }
    #[test]
    fn test_group_empty() {
        let group = ShdfGroup::new("root");
        assert_eq!(group.name, "root");
        assert!(group.datasets.is_empty());
        assert!(group.children.is_empty());
        assert_eq!(group.total_datasets(), 0);
    }
    #[test]
    fn test_group_add_dataset() {
        let mut group = ShdfGroup::new("root");
        group.add_dataset_f64("temperatures", vec![10], vec![0.0; 10]);
        assert_eq!(group.datasets.len(), 1);
        assert_eq!(group.total_datasets(), 1);
        assert!(group.get_dataset("temperatures").is_some());
    }
    #[test]
    fn test_group_nested() {
        let mut root = ShdfGroup::new("root");
        let mut child = ShdfGroup::new("particles");
        child.add_dataset_f64("positions", vec![100, 3], vec![0.0; 300]);
        child.add_dataset_i32("types", vec![100], vec![0; 100]);
        root.add_child(child);
        assert_eq!(root.total_datasets(), 2);
        let particles = root.get_child("particles").unwrap();
        assert_eq!(particles.datasets.len(), 2);
    }
    #[test]
    fn test_group_attributes() {
        let mut group = ShdfGroup::new("simulation");
        group.add_attribute("timestep", "0.001");
        group.add_attribute("units", "SI");
        assert_eq!(group.attributes.len(), 2);
    }
    #[test]
    fn test_group_summary() {
        let mut root = ShdfGroup::new("root");
        root.add_dataset_f64("energy", vec![1], vec![42.0]);
        let summary = root.summary(0);
        assert!(summary.contains("root"));
        assert!(summary.contains("energy"));
    }
    #[test]
    fn test_group_get_missing_dataset() {
        let group = ShdfGroup::new("root");
        assert!(group.get_dataset("nonexistent").is_none());
    }
    #[test]
    fn test_group_get_missing_child() {
        let group = ShdfGroup::new("root");
        assert!(group.get_child("nonexistent").is_none());
    }
    #[test]
    fn test_chunking_n_chunks_1d() {
        let config = ChunkingConfig::new(vec![10]);
        let n = config.n_chunks(&[100]);
        assert_eq!(n, 10);
    }
    #[test]
    fn test_chunking_n_chunks_2d() {
        let config = ChunkingConfig::new(vec![10, 5]);
        let n = config.n_chunks(&[20, 15]);
        assert_eq!(n, 2 * 3);
    }
    #[test]
    fn test_chunking_n_chunks_remainder() {
        let config = ChunkingConfig::new(vec![10]);
        let n = config.n_chunks(&[15]);
        assert_eq!(n, 2);
    }
    #[test]
    fn test_chunking_chunk_index() {
        let config = ChunkingConfig::new(vec![10]);
        assert_eq!(config.chunk_index(&[5], &[100]), 0);
        assert_eq!(config.chunk_index(&[15], &[100]), 1);
    }
    #[test]
    fn test_chunking_default() {
        let config = ChunkingConfig::default_for_shape(&[1000, 3]);
        assert_eq!(config.chunk_dims, vec![64, 3]);
    }
    #[test]
    fn test_chunking_mismatched_dims() {
        let config = ChunkingConfig::new(vec![10, 10]);
        assert_eq!(config.n_chunks(&[100]), 0);
    }
    #[test]
    fn test_delta_encode_decode_f64() {
        let data = vec![1.0, 3.0, 6.0, 10.0, 15.0];
        let encoded = CompressionSettings::delta_encode_f64(&data);
        let decoded = CompressionSettings::delta_decode_f64(&encoded);
        for (a, b) in data.iter().zip(decoded.iter()) {
            assert!((a - b).abs() < 1e-14, "Mismatch: {a} vs {b}");
        }
    }
    #[test]
    fn test_delta_encode_decode_i32() {
        let data = vec![10, 20, 30, 25, 35];
        let encoded = CompressionSettings::delta_encode_i32(&data);
        let decoded = CompressionSettings::delta_decode_i32(&encoded);
        assert_eq!(data, decoded);
    }
    #[test]
    fn test_delta_encode_empty() {
        assert!(CompressionSettings::delta_encode_f64(&[]).is_empty());
        assert!(CompressionSettings::delta_decode_f64(&[]).is_empty());
        assert!(CompressionSettings::delta_encode_i32(&[]).is_empty());
        assert!(CompressionSettings::delta_decode_i32(&[]).is_empty());
    }
    #[test]
    fn test_delta_encode_single() {
        let data = vec![42.0];
        let encoded = CompressionSettings::delta_encode_f64(&data);
        let decoded = CompressionSettings::delta_decode_f64(&encoded);
        assert_eq!(decoded, data);
    }
    #[test]
    fn test_compression_none() {
        let settings = CompressionSettings::none();
        assert_eq!(settings.algorithm, CompressionAlgorithm::None);
        assert_eq!(settings.level, 0);
    }
    #[test]
    fn test_compression_delta() {
        let settings = CompressionSettings::delta();
        assert_eq!(settings.algorithm, CompressionAlgorithm::Delta);
    }
    #[test]
    fn test_attribute_value_string_roundtrip() {
        let val = AttributeValue::String("hello".to_string());
        let s = AttributeHelper::to_string(&val);
        let recovered = AttributeHelper::from_string(&s);
        assert_eq!(recovered, val);
    }
    #[test]
    fn test_attribute_value_float_roundtrip() {
        let val = AttributeValue::Float64(3.125);
        let s = AttributeHelper::to_string(&val);
        let recovered = AttributeHelper::from_string(&s);
        if let AttributeValue::Float64(f) = recovered {
            assert!((f - 3.125).abs() < 1e-10);
        } else {
            panic!("Expected Float64");
        }
    }
    #[test]
    fn test_attribute_value_int_roundtrip() {
        let val = AttributeValue::Int32(42);
        let s = AttributeHelper::to_string(&val);
        let recovered = AttributeHelper::from_string(&s);
        assert_eq!(recovered, val);
    }
    #[test]
    fn test_attribute_value_bool_roundtrip() {
        let val = AttributeValue::Bool(true);
        let s = AttributeHelper::to_string(&val);
        let recovered = AttributeHelper::from_string(&s);
        assert_eq!(recovered, val);
        let val2 = AttributeValue::Bool(false);
        let s2 = AttributeHelper::to_string(&val2);
        let recovered2 = AttributeHelper::from_string(&s2);
        assert_eq!(recovered2, val2);
    }
    #[test]
    fn test_attribute_value_unknown_prefix() {
        let recovered = AttributeHelper::from_string("unknown_value");
        assert_eq!(
            recovered,
            AttributeValue::String("unknown_value".to_string())
        );
    }
    #[test]
    fn test_schema_valid_file() {
        let mut schema = ShdfSchema::new();
        schema.expect_dataset("positions", DataType::Float64);
        schema.require_attribute("author");
        let mut file = ShdfFile::new();
        file.add_dataset_f64("positions", vec![10, 3], vec![0.0; 30]);
        file.add_global_attr("author", "Test");
        let errors = schema.validate(&file);
        assert!(errors.is_empty(), "Expected no errors, got: {:?}", errors);
    }
    #[test]
    fn test_schema_missing_dataset() {
        let mut schema = ShdfSchema::new();
        schema.expect_dataset("positions", DataType::Float64);
        let file = ShdfFile::new();
        let errors = schema.validate(&file);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("Missing dataset: positions"));
    }
    #[test]
    fn test_schema_wrong_dtype() {
        let mut schema = ShdfSchema::new();
        schema.expect_dataset("indices", DataType::Float64);
        let mut file = ShdfFile::new();
        file.add_dataset_i32("indices", vec![10], vec![0; 10]);
        let errors = schema.validate(&file);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("expected"));
    }
    #[test]
    fn test_schema_missing_attribute() {
        let mut schema = ShdfSchema::new();
        schema.require_attribute("version");
        let file = ShdfFile::new();
        let errors = schema.validate(&file);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("Missing global attribute: version"));
    }
    #[test]
    fn test_schema_empty_valid() {
        let schema = ShdfSchema::new();
        let file = ShdfFile::new();
        let errors = schema.validate(&file);
        assert!(errors.is_empty());
    }
    #[test]
    fn test_datatype_equality() {
        assert_eq!(DataType::Float64, DataType::Float64);
        assert_ne!(DataType::Float64, DataType::Float32);
        assert_ne!(DataType::Int32, DataType::Int64);
    }
    #[test]
    fn test_from_bytes_bad_magic() {
        let data = b"BAAD\x01\x00\x00\x00";
        let result = ShdfFile::from_bytes(data);
        assert!(result.is_err());
    }
    #[test]
    fn test_from_bytes_too_short() {
        let result = ShdfFile::from_bytes(&[0u8; 3]);
        assert!(result.is_err());
    }
    #[test]
    fn test_from_bytes_bad_version() {
        let mut data = Vec::new();
        data.extend_from_slice(b"SHDF");
        data.extend_from_slice(&99u32.to_le_bytes());
        let result = ShdfFile::from_bytes(&data);
        assert!(result.is_err());
    }
}
/// Write an XDMF XML file that points to data in an HDF5/SHDF file.
#[allow(dead_code)]
pub fn write_xdmf(path: &str, params: &XdmfParams) -> std::io::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, r#"<?xml version="1.0" ?>"#)?;
    writeln!(w, r#"<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>"#)?;
    writeln!(w, r#"<Xdmf Version="2.0">"#)?;
    writeln!(w, r#"  <Domain>"#)?;
    writeln!(w, r#"    <Grid Name="Mesh" GridType="Uniform">"#)?;
    writeln!(
        w,
        r#"      <Topology TopologyType="{}" NumberOfElements="{}">"#,
        params.topology.as_str(),
        params.n_elements
    )?;
    writeln!(
        w,
        r#"        <DataItem Dimensions="{} {}" NumberType="Int" Format="HDF">"#,
        params.n_elements, params.nodes_per_element
    )?;
    writeln!(
        w,
        r#"          {}:{}"#,
        params.hdf5_path, params.connectivity_dataset
    )?;
    writeln!(w, r#"        </DataItem>"#)?;
    writeln!(w, r#"      </Topology>"#)?;
    writeln!(w, r#"      <Geometry GeometryType="XYZ">"#)?;
    writeln!(
        w,
        r#"        <DataItem Dimensions="{} 3" NumberType="Float" Precision="8" Format="HDF">"#,
        params.n_nodes
    )?;
    writeln!(
        w,
        r#"          {}:{}"#,
        params.hdf5_path, params.coords_dataset
    )?;
    writeln!(w, r#"        </DataItem>"#)?;
    writeln!(w, r#"      </Geometry>"#)?;
    for (attr_name, ds_path) in &params.attributes {
        writeln!(
            w,
            r#"      <Attribute Name="{attr_name}" AttributeType="Scalar" Center="Node">"#
        )?;
        writeln!(
            w,
            r#"        <DataItem Dimensions="{}" NumberType="Float" Precision="8" Format="HDF">"#,
            params.n_nodes
        )?;
        writeln!(w, r#"          {}:{}"#, params.hdf5_path, ds_path)?;
        writeln!(w, r#"        </DataItem>"#)?;
        writeln!(w, r#"      </Attribute>"#)?;
    }
    writeln!(w, r#"    </Grid>"#)?;
    writeln!(w, r#"  </Domain>"#)?;
    writeln!(w, r#"</Xdmf>"#)?;
    w.flush()?;
    Ok(())
}
#[cfg(test)]
mod tests_hdf5_extended {
    use super::*;

    use crate::hdf5_simple::ChunkedDataset;

    use crate::hdf5_simple::CompoundDataset;
    use crate::hdf5_simple::CompoundField;

    use crate::hdf5_simple::CompressionLevel;

    use crate::hdf5_simple::DeflateMetadata;

    use crate::hdf5_simple::XdmfTopologyType;
    use crate::hdf5_simple::types::*;
    #[test]
    fn test_chunked_dataset_n_elements() {
        let ds = ChunkedDataset::new("data", vec![100, 50], vec![10, 10]);
        assert_eq!(ds.n_elements(), 5000);
    }
    #[test]
    fn test_chunked_dataset_n_chunks_1d() {
        let ds = ChunkedDataset::new("v", vec![100], vec![32]);
        let nchunks = ds.n_chunks_per_dim();
        assert_eq!(nchunks[0], 4);
    }
    #[test]
    fn test_chunked_dataset_total_chunks() {
        let ds = ChunkedDataset::new("v", vec![100, 50], vec![10, 10]);
        assert_eq!(ds.total_chunks(), 50);
    }
    #[test]
    fn test_chunked_dataset_write_read_chunk() {
        let mut ds = ChunkedDataset::new("v", vec![10], vec![5]);
        ds.write_chunk_1d(0, &[1.0, 2.0, 3.0, 4.0, 5.0]);
        ds.write_chunk_1d(1, &[6.0, 7.0, 8.0, 9.0, 10.0]);
        assert!((ds.data[0] - 1.0).abs() < 1e-12);
        assert!((ds.data[9] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_chunked_dataset_to_bytes_nonempty() {
        let mut ds = ChunkedDataset::new("pressure", vec![4], vec![2]);
        ds.data = vec![1.0, 2.0, 3.0, 4.0];
        let bytes = ds.to_bytes();
        assert!(!bytes.is_empty());
    }
    #[test]
    fn test_chunked_dataset_attrs() {
        let mut ds = ChunkedDataset::new("temp", vec![5], vec![5]);
        ds.add_attr("units", "K");
        ds.add_attr("long_name", "Temperature");
        assert_eq!(ds.attrs.len(), 2);
    }
    #[test]
    fn test_shdf_group_new() {
        let g = ShdfGroup::new("root");
        assert_eq!(g.name, "root");
        assert!(g.datasets.is_empty());
        assert!(g.children.is_empty());
    }
    #[test]
    fn test_shdf_group_add_dataset() {
        let mut g = ShdfGroup::new("data");
        g.add_dataset_f64("velocity", vec![3], vec![1.0, 2.0, 3.0]);
        assert_eq!(g.datasets.len(), 1);
        assert!(g.get_dataset("velocity").is_some());
    }
    #[test]
    fn test_shdf_group_hierarchy() {
        let mut root = ShdfGroup::new("root");
        let child1 = ShdfGroup::new("geometry");
        let child2 = ShdfGroup::new("fields");
        root.add_child(child1);
        root.add_child(child2);
        assert_eq!(root.children.len(), 2);
        assert!(root.get_child("geometry").is_some());
        assert!(root.get_child("missing").is_none());
    }
    #[test]
    fn test_shdf_group_total_datasets_recursive() {
        let mut root = ShdfGroup::new("root");
        let mut child = ShdfGroup::new("child");
        child.add_dataset_f64("x", vec![1], vec![0.0]);
        root.add_child(child);
        assert_eq!(root.total_datasets(), 1);
    }
    #[test]
    fn test_shdf_group_name() {
        let g = ShdfGroup::new("geometry");
        assert_eq!(g.name, "geometry");
    }
    #[test]
    fn test_write_xdmf_creates_file() {
        let params = XdmfParams {
            hdf5_path: "sim.shdf".to_string(),
            coords_dataset: "/coordinates".to_string(),
            connectivity_dataset: "/connectivity".to_string(),
            n_nodes: 100,
            n_elements: 50,
            nodes_per_element: 4,
            topology: XdmfTopologyType::Tetrahedron,
            attributes: vec![("pressure".to_string(), "/fields/pressure".to_string())],
        };
        let path = "/tmp/test_oxiphysics_xdmf.xdmf";
        write_xdmf(path, &params).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("Xdmf"));
        assert!(content.contains("Tetrahedron"));
        assert!(content.contains("pressure"));
        assert!(content.contains("sim.shdf"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_xdmf_topology_type_str() {
        assert_eq!(XdmfTopologyType::Triangle.as_str(), "Triangle");
        assert_eq!(XdmfTopologyType::Hexahedron.as_str(), "Hexahedron");
    }
    #[test]
    fn test_compression_level_values() {
        assert_eq!(CompressionLevel::None.level(), 0);
        assert_eq!(CompressionLevel::Fast.level(), 1);
        assert_eq!(CompressionLevel::Balanced.level(), 5);
        assert_eq!(CompressionLevel::Maximum.level(), 9);
    }
    #[test]
    fn test_compression_level_is_compressed() {
        assert!(!CompressionLevel::None.is_compressed());
        assert!(CompressionLevel::Fast.is_compressed());
        assert!(CompressionLevel::Maximum.is_compressed());
    }
    #[test]
    fn test_deflate_metadata_uncompressed() {
        let meta = DeflateMetadata::uncompressed(1000);
        assert_eq!(meta.compression_ratio(), 1.0);
        assert_eq!(meta.space_savings(), 0.0);
    }
    #[test]
    fn test_deflate_metadata_compression_ratio() {
        let meta = DeflateMetadata {
            level: CompressionLevel::Balanced,
            shuffle: true,
            chunk_shape: vec![100],
            compressed_size: 400,
            uncompressed_size: 1000,
        };
        assert!(
            (meta.compression_ratio() - 2.5).abs() < 1e-10,
            "ratio={}",
            meta.compression_ratio()
        );
        assert!(
            (meta.space_savings() - 0.6).abs() < 1e-10,
            "savings={}",
            meta.space_savings()
        );
    }
    #[test]
    fn test_deflate_metadata_zero_compressed() {
        let meta = DeflateMetadata {
            level: CompressionLevel::Balanced,
            shuffle: false,
            chunk_shape: vec![],
            compressed_size: 0,
            uncompressed_size: 1000,
        };
        assert_eq!(meta.compression_ratio(), 1.0);
    }
    #[test]
    fn test_compound_dataset_new() {
        let ds = CompoundDataset::new("particles", 10);
        assert_eq!(ds.n_records, 10);
        assert_eq!(ds.n_fields(), 0);
    }
    #[test]
    fn test_compound_dataset_add_field() {
        let mut ds = CompoundDataset::new("atoms", 3);
        ds.add_field(CompoundField::new(
            "x",
            DataType::Float64,
            vec![1.0, 2.0, 3.0],
        ));
        ds.add_field(CompoundField::new(
            "y",
            DataType::Float64,
            vec![4.0, 5.0, 6.0],
        ));
        assert_eq!(ds.n_fields(), 2);
        let xs = ds.get_field("x").unwrap();
        assert_eq!(xs.len(), 3);
        assert!((xs[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_compound_dataset_get_missing_field() {
        let ds = CompoundDataset::new("particles", 5);
        assert!(ds.get_field("missing").is_none());
    }
    #[test]
    fn test_compound_dataset_to_csv_bytes() {
        let mut ds = CompoundDataset::new("data", 2);
        ds.add_field(CompoundField::new("a", DataType::Float64, vec![1.0, 2.0]));
        ds.add_field(CompoundField::new("b", DataType::Float32, vec![3.0, 4.0]));
        let csv = ds.to_csv_bytes();
        let s = String::from_utf8(csv).unwrap();
        assert!(s.contains("a,b"));
        assert!(s.contains("1.000000"));
    }
    #[test]
    fn test_compound_dataset_attrs() {
        let mut ds = CompoundDataset::new("traj", 0);
        ds.add_attr("source", "MD simulation");
        ds.add_attr("units", "nm");
        assert_eq!(ds.attrs.len(), 2);
    }
    #[test]
    #[should_panic]
    fn test_compound_dataset_wrong_field_length() {
        let mut ds = CompoundDataset::new("atoms", 5);
        ds.add_field(CompoundField::new(
            "x",
            DataType::Float64,
            vec![1.0, 2.0, 3.0],
        ));
    }
}
/// Generate a CDL (Common Data form Language) text representation of an
/// [`ShdfGroup`] hierarchy for debugging purposes.
pub fn cdl_dump(root: &ShdfGroup) -> String {
    let mut out = String::new();
    cdl_dump_group(root, 0, &mut out);
    out
}
pub(super) fn cdl_dump_group(group: &ShdfGroup, depth: usize, out: &mut String) {
    let indent = "  ".repeat(depth);
    out.push_str(&format!("{}group: {} {{\n", indent, group.name));
    for ds in &group.datasets {
        cdl_dump_dataset(ds, depth + 1, out);
    }
    for child in &group.children {
        cdl_dump_group(child, depth + 1, out);
    }
    out.push_str(&format!("{}}} // group: {}\n", indent, group.name));
}
pub(super) fn cdl_dump_dataset(ds: &Dataset, depth: usize, out: &mut String) {
    let indent = "  ".repeat(depth);
    let type_str = match ds.dtype {
        DataType::Float64 => "double",
        DataType::Float32 => "float",
        DataType::Int32 => "int",
        DataType::Int64 => "int64",
    };
    let shape_str: Vec<String> = ds.shape.iter().map(|d| d.to_string()).collect();
    out.push_str(&format!(
        "{}{}  {}({}) ;\n",
        indent,
        type_str,
        ds.name,
        shape_str.join(", ")
    ));
    for (k, v) in &ds.attributes {
        out.push_str(&format!("{}    {}:{} = \"{}\" ;\n", indent, ds.name, k, v));
    }
    let preview_len = ds.data_f64.len().min(8);
    if preview_len > 0 {
        let vals: Vec<String> = ds.data_f64[..preview_len]
            .iter()
            .map(|v| format!("{:.6}", v))
            .collect();
        let ellipsis = if ds.data_f64.len() > 8 { ", ..." } else { "" };
        out.push_str(&format!(
            "{}    // data = {}{} ;\n",
            indent,
            vals.join(", "),
            ellipsis
        ));
    }
}
#[cfg(test)]
mod tests_hdf5_additions {
    use super::*;

    use crate::hdf5_simple::DatasetStats;

    use crate::hdf5_simple::GroupNavigator;

    use crate::hdf5_simple::TimeSeriesAppender;
    use crate::hdf5_simple::VirtualLink;

    use crate::hdf5_simple::types::*;
    #[test]
    fn test_navigator_get_dataset_simple() {
        let mut root = ShdfGroup::new("root");
        root.add_dataset_f64("pos", vec![3], vec![1.0, 2.0, 3.0]);
        let nav = GroupNavigator::new(root);
        assert!(nav.get_dataset("/root/pos").is_some());
        assert!(nav.get_dataset("/root/missing").is_none());
    }
    #[test]
    fn test_navigator_nested_path() {
        let mut root = ShdfGroup::new("root");
        let mut sim = ShdfGroup::new("simulation");
        let mut atoms = ShdfGroup::new("atoms");
        atoms.add_dataset_f64("positions", vec![10, 3], vec![0.0; 30]);
        sim.add_child(atoms);
        root.add_child(sim);
        let nav = GroupNavigator::new(root);
        assert!(
            nav.get_dataset("/root/simulation/atoms/positions")
                .is_some()
        );
        assert_eq!(nav.total_datasets(), 1);
    }
    #[test]
    fn test_navigator_all_paths() {
        let mut root = ShdfGroup::new("root");
        root.add_dataset_f64("a", vec![1], vec![1.0]);
        root.add_dataset_f64("b", vec![1], vec![2.0]);
        let nav = GroupNavigator::new(root);
        let paths = nav.all_paths();
        assert_eq!(paths.len(), 2);
        assert!(paths.iter().any(|p| p.ends_with("/a")));
        assert!(paths.iter().any(|p| p.ends_with("/b")));
    }
    #[test]
    fn test_navigator_empty_group() {
        let root = ShdfGroup::new("empty");
        let nav = GroupNavigator::new(root);
        assert_eq!(nav.all_paths().len(), 0);
        assert_eq!(nav.total_datasets(), 0);
    }
    #[test]
    fn test_navigator_leading_slash_optional() {
        let mut root = ShdfGroup::new("root");
        root.add_dataset_f64("v", vec![2], vec![1.0, 2.0]);
        let nav = GroupNavigator::new(root);
        assert!(nav.get_dataset("/root/v").is_some());
    }
    #[test]
    fn test_timeseries_appender_basic() {
        let mut app = TimeSeriesAppender::new("temperature", 3);
        app.append(&[300.0, 301.0, 302.0]);
        app.append(&[303.0, 304.0, 305.0]);
        assert_eq!(app.n_frames, 2);
        assert_eq!(app.total_samples(), 6);
    }
    #[test]
    fn test_timeseries_get_frame() {
        let mut app = TimeSeriesAppender::new("vel", 3);
        app.append(&[1.0, 2.0, 3.0]);
        app.append(&[4.0, 5.0, 6.0]);
        let frame = app.get_frame(1).unwrap();
        assert_eq!(frame, &[4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_timeseries_get_frame_out_of_bounds() {
        let app = TimeSeriesAppender::new("x", 2);
        assert!(app.get_frame(0).is_none());
    }
    #[test]
    fn test_timeseries_to_dataset_shape() {
        let mut app = TimeSeriesAppender::new("pressure", 1);
        app.append(&[1.0]);
        app.append(&[2.0]);
        app.append(&[3.0]);
        let ds = app.to_dataset();
        assert_eq!(ds.shape, vec![3, 1]);
        assert!(ds.attributes.iter().any(|(k, _)| k == "n_frames"));
    }
    #[test]
    #[should_panic]
    fn test_timeseries_wrong_frame_width() {
        let mut app = TimeSeriesAppender::new("x", 3);
        app.append(&[1.0, 2.0]);
    }
    #[test]
    fn test_cdl_dump_contains_group_name() {
        let mut root = ShdfGroup::new("simulation");
        root.add_dataset_f64("energy", vec![5], vec![1.0; 5]);
        let cdl = cdl_dump(&root);
        assert!(cdl.contains("simulation"), "CDL should contain group name");
        assert!(cdl.contains("energy"), "CDL should contain dataset name");
        assert!(cdl.contains("double"), "CDL should contain type");
    }
    #[test]
    fn test_cdl_dump_nested() {
        let mut root = ShdfGroup::new("root");
        let mut child = ShdfGroup::new("atoms");
        child.add_dataset_f64("x", vec![4], vec![1.0, 2.0, 3.0, 4.0]);
        root.add_child(child);
        let cdl = cdl_dump(&root);
        assert!(cdl.contains("atoms"));
        assert!(cdl.contains("x"));
    }
    #[test]
    fn test_cdl_dump_with_attrs() {
        let mut root = ShdfGroup::new("meta");
        let mut ds = Dataset {
            name: "temp".to_string(),
            dtype: DataType::Float32,
            shape: vec![2],
            data_f64: vec![300.0, 310.0],
            data_i32: Vec::new(),
            attributes: vec![("units".to_string(), "K".to_string())],
        };
        ds.attributes
            .push(("source".to_string(), "sensor".to_string()));
        root.datasets.push(ds);
        let cdl = cdl_dump(&root);
        assert!(cdl.contains("units"));
        assert!(cdl.contains("source"));
    }
    #[test]
    fn test_cdl_dump_int32_type_label() {
        let mut root = ShdfGroup::new("ids");
        root.datasets.push(Dataset {
            name: "atom_id".to_string(),
            dtype: DataType::Int32,
            shape: vec![3],
            data_f64: vec![1.0, 2.0, 3.0],
            data_i32: Vec::new(),
            attributes: vec![],
        });
        let cdl = cdl_dump(&root);
        assert!(cdl.contains("int"), "CDL should label Int32 as 'int'");
    }
    #[test]
    fn test_virtual_link_no_slice() {
        let link = VirtualLink::new("/virtual/pos", "sim.shdf", "/atoms/positions");
        let cdl = link.to_cdl();
        assert!(cdl.contains("/virtual/pos"));
        assert!(cdl.contains("sim.shdf"));
        assert!(cdl.contains("(:)"), "no-slice should show '(:)'");
    }
    #[test]
    fn test_virtual_link_with_slice() {
        let link = VirtualLink::new("/v/vel", "traj.shdf", "/vel")
            .with_slice(0, 100, 1)
            .with_slice(0, 3, 1);
        let cdl = link.to_cdl();
        assert!(
            cdl.contains("0:100:1"),
            "should contain first dimension slice"
        );
        assert!(
            cdl.contains("0:3:1"),
            "should contain second dimension slice"
        );
    }
    #[test]
    fn test_virtual_link_clone() {
        let link = VirtualLink::new("/a", "f.shdf", "/b");
        let cloned = link.clone();
        assert_eq!(cloned.virtual_path, link.virtual_path);
    }
    #[test]
    fn test_dataset_stats_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = DatasetStats::from_slice(&data).unwrap();
        assert!((stats.min - 1.0).abs() < 1e-12);
        assert!((stats.max - 5.0).abs() < 1e-12);
        assert!((stats.mean - 3.0).abs() < 1e-12);
        assert!((stats.range() - 4.0).abs() < 1e-12);
        assert_eq!(stats.count, 5);
    }
    #[test]
    fn test_dataset_stats_single_element() {
        let stats = DatasetStats::from_slice(&[42.0]).unwrap();
        assert!((stats.min - 42.0).abs() < 1e-12);
        assert!((stats.max - 42.0).abs() < 1e-12);
        assert!((stats.mean - 42.0).abs() < 1e-12);
        assert!((stats.variance).abs() < 1e-12);
        assert!((stats.std_dev()).abs() < 1e-12);
    }
    #[test]
    fn test_dataset_stats_empty() {
        assert!(DatasetStats::from_slice(&[]).is_none());
    }
    #[test]
    fn test_dataset_stats_std_dev_constant() {
        let data = vec![5.0; 10];
        let stats = DatasetStats::from_slice(&data).unwrap();
        assert!(
            stats.std_dev().abs() < 1e-12,
            "constant data has zero std dev"
        );
    }
    #[test]
    fn test_dataset_stats_variance_known() {
        let data = vec![0.0, 2.0, 4.0];
        let stats = DatasetStats::from_slice(&data).unwrap();
        let expected_var = 8.0 / 3.0;
        assert!(
            (stats.variance - expected_var).abs() < 1e-10,
            "variance={}",
            stats.variance
        );
    }
}
