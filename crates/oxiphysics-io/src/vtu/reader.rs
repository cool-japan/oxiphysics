// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reader for VTK XML unstructured grid (`.vtu`) files.
//!
//! The XML is tokenised with the pure-Rust `oxixml-xml` pull parser and every
//! `DataArray` is interpreted through its attributes (`type`,
//! `NumberOfComponents`, `format`, `Name`), never through the textual layout of
//! the file: files with all content on one line, `Float32` or `Int64` arrays,
//! several `Piece`s and inline base64 `binary` arrays are all read correctly.
//!
//! Supported array encodings:
//!
//! * `format="ascii"` — whitespace separated numbers;
//! * `format="binary"` — inline base64 with an unsigned 32- or 64-bit length
//!   header (`header_type` attribute of `VTKFile`), little or big endian
//!   (`byte_order`), written without compression.  Both the VTK layout (header
//!   and payload base64-encoded separately) and a single contiguous base64
//!   stream are accepted.
//!
//! Encodings that are *not* supported — `format="appended"` data and
//! compressed binary data (`compressor` attribute) — are reported as
//! [`Error::Parse`] instead of silently producing an empty mesh.

use crate::mesh_io::base64_decode_buffer;
use crate::{Error, Result};
use oxixml_xml::XmlVersion;
use oxixml_xml::events::{BytesStart, Event};
use oxixml_xml::reader::Reader;
use std::path::Path;

/// A named point-data or cell-data array read from a VTU file.
#[derive(Debug, Clone, PartialEq)]
pub struct VtuDataArray {
    /// Array name (`Name` attribute).
    pub name: String,
    /// Number of components per tuple (`NumberOfComponents`, default 1).
    pub num_components: usize,
    /// All values, tuple after tuple (`num_tuples() * num_components` values).
    pub values: Vec<f64>,
}

impl VtuDataArray {
    /// Number of tuples (points or cells) the array holds.
    pub fn num_tuples(&self) -> usize {
        self.values
            .len()
            .checked_div(self.num_components)
            .unwrap_or(0)
    }
}

/// Reader for VTK XML unstructured grids: points, cell connectivity, cell
/// types and the point/cell data arrays.
///
/// When the file holds several `Piece`s they are concatenated: connectivity
/// indices and offsets of later pieces are shifted so that they address the
/// combined point list, and data arrays are appended (every piece must provide
/// the same arrays).
#[derive(Debug, Clone)]
pub struct VtuReader {
    pub(super) points: Vec<[f64; 3]>,
    pub(super) num_cells: usize,
    pub(super) cell_types: Vec<u8>,
    pub(super) connectivity: Vec<usize>,
    pub(super) offsets: Vec<usize>,
    pub(super) point_scalar_names: Vec<String>,
    pub(super) cell_scalar_names: Vec<String>,
    pub(super) point_arrays: Vec<VtuDataArray>,
    pub(super) cell_arrays: Vec<VtuDataArray>,
}

impl VtuReader {
    /// Parse a VTU XML string.
    ///
    /// Reads the point coordinates, the cell connectivity (`connectivity`,
    /// `offsets`, `types`) and all named `PointData`/`CellData` arrays, and
    /// validates them against the `NumberOfPoints` / `NumberOfCells` counts of
    /// each `Piece`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the XML is malformed, the document is not an
    /// `UnstructuredGrid` VTK file, an array uses an unsupported encoding
    /// (`appended`, compressed binary), numbers cannot be decoded, array
    /// lengths contradict the piece counts, or a cell refers to a point that
    /// does not exist.
    pub fn from_xml(data: &str) -> Result<Self> {
        let mut reader = Reader::from_str(data);
        let mut parser = VtuParser::new();
        loop {
            let event = reader.read_event().map_err(|err| {
                Error::Parse(format!(
                    "VTU XML error at byte {}: {err}",
                    reader.error_position()
                ))
            })?;
            match event {
                Event::Eof => break,
                Event::Start(element) => parser.start_element(&element)?,
                Event::Empty(element) => {
                    parser.start_element(&element)?;
                    parser.end_element()?;
                }
                Event::End(_) => parser.end_element()?,
                Event::Text(text) if parser.collecting_text() => {
                    let text = text
                        .decode()
                        .map_err(|err| Error::Parse(format!("invalid text encoding: {err}")))?;
                    parser.push_text(&text);
                }
                Event::CData(text) if parser.collecting_text() => {
                    let text = text
                        .decode()
                        .map_err(|err| Error::Parse(format!("invalid text encoding: {err}")))?;
                    parser.push_text(&text);
                }
                Event::GeneralRef(reference) if parser.collecting_text() => {
                    match reference.resolve_char_ref() {
                        Ok(Some(c)) => parser.push_text(c.encode_utf8(&mut [0u8; 4])),
                        Ok(None) => {
                            return Err(Error::Parse(
                                "entity reference inside numeric DataArray content".to_string(),
                            ));
                        }
                        Err(err) => {
                            return Err(Error::Parse(format!(
                                "invalid character reference in DataArray content: {err}"
                            )));
                        }
                    }
                }
                _ => {}
            }
        }
        parser.finish()
    }

    /// Read and parse the VTU file at `path`.
    ///
    /// Files that are not valid UTF-8 text — notably ParaView's default
    /// `appended` raw binary layout — are rejected with [`Error::Parse`].
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let bytes = std::fs::read(path.as_ref())?;
        match String::from_utf8(bytes) {
            Ok(text) => Self::from_xml(&text),
            Err(err) => {
                let raw = err.into_bytes();
                let has_appended = raw
                    .windows(b"<AppendedData".len())
                    .any(|window| window == b"<AppendedData");
                if has_appended {
                    Err(Error::Parse(
                        "VTU file stores its arrays as raw appended binary data, which is not \
                         supported; write it with ascii or inline binary data arrays"
                            .to_string(),
                    ))
                } else {
                    Err(Error::Parse("VTU file is not valid UTF-8 text".to_string()))
                }
            }
        }
    }

    /// Return parsed point coordinates.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }
    /// Return the number of points.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }
    /// Return the number of cells reported.
    pub fn num_cells(&self) -> usize {
        self.num_cells
    }
    /// Return the cell types parsed from the file.
    pub fn cell_types(&self) -> &[u8] {
        &self.cell_types
    }
    /// Return the flat connectivity array: the point indices of all cells,
    /// concatenated in cell order.
    pub fn connectivity(&self) -> &[usize] {
        &self.connectivity
    }
    /// Return the offsets array: `offsets()[c]` is the end of cell `c` inside
    /// [`connectivity`](Self::connectivity).
    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }
    /// Return the point indices of cell `cell`, or `None` if it does not exist.
    pub fn cell_nodes(&self, cell: usize) -> Option<&[usize]> {
        let end = *self.offsets.get(cell)?;
        let start = if cell == 0 {
            0
        } else {
            *self.offsets.get(cell - 1)?
        };
        self.connectivity.get(start..end)
    }
    /// Return point data array names found in the file.
    pub fn point_data_names(&self) -> &[String] {
        &self.point_scalar_names
    }
    /// Return cell data array names found in the file.
    pub fn cell_data_names(&self) -> &[String] {
        &self.cell_scalar_names
    }
    /// Return all point data arrays.
    pub fn point_data_arrays(&self) -> &[VtuDataArray] {
        &self.point_arrays
    }
    /// Return all cell data arrays.
    pub fn cell_data_arrays(&self) -> &[VtuDataArray] {
        &self.cell_arrays
    }
    /// Return the point data array called `name`.
    pub fn point_data_array(&self, name: &str) -> Option<&VtuDataArray> {
        self.point_arrays.iter().find(|array| array.name == name)
    }
    /// Return the cell data array called `name`.
    pub fn cell_data_array(&self, name: &str) -> Option<&VtuDataArray> {
        self.cell_arrays.iter().find(|array| array.name == name)
    }

    /// Read a cell data array from a VTU XML string by name.
    ///
    /// Returns the values of the requested array as a `Vec<f64>` (all
    /// components of a tuple interleaved), or an error message if the document
    /// cannot be parsed (see [`VtuReader::from_xml`]) or has no array with that
    /// name in a `CellData` section.
    pub fn read_cell_data(data: &str, array_name: &str) -> std::result::Result<Vec<f64>, String> {
        let reader = Self::from_xml(data).map_err(|err| err.to_string())?;
        reader
            .cell_data_array(array_name)
            .map(|array| array.values.clone())
            .ok_or_else(|| format!("cell data array '{}' not found", array_name))
    }

    /// Return the bounding box of the parsed points.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.points.is_empty() {
            return None;
        }
        let mut min = self.points[0];
        let mut max = self.points[0];
        for p in &self.points[1..] {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        Some((min, max))
    }
}

// ---------------------------------------------------------------------------
// Parsing machinery
// ---------------------------------------------------------------------------

/// Scalar element type of a `DataArray` (`type` attribute).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScalarType {
    Int8,
    UInt8,
    Int16,
    UInt16,
    Int32,
    UInt32,
    Int64,
    UInt64,
    Float32,
    Float64,
}

impl ScalarType {
    fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "Int8" => Self::Int8,
            "UInt8" => Self::UInt8,
            "Int16" => Self::Int16,
            "UInt16" => Self::UInt16,
            "Int32" => Self::Int32,
            "UInt32" => Self::UInt32,
            "Int64" => Self::Int64,
            "UInt64" => Self::UInt64,
            "Float32" => Self::Float32,
            "Float64" => Self::Float64,
            _ => return None,
        })
    }

    fn byte_size(self) -> usize {
        match self {
            Self::Int8 | Self::UInt8 => 1,
            Self::Int16 | Self::UInt16 => 2,
            Self::Int32 | Self::UInt32 | Self::Float32 => 4,
            Self::Int64 | Self::UInt64 | Self::Float64 => 8,
        }
    }

    /// Decode one little- or big-endian element from `chunk` (exactly
    /// `byte_size()` bytes) into an `f64`.
    fn decode(self, chunk: &[u8], big_endian: bool) -> Option<f64> {
        macro_rules! read {
            ($ty:ty, $n:literal) => {{
                let bytes: [u8; $n] = *chunk.first_chunk::<$n>()?;
                if big_endian {
                    <$ty>::from_be_bytes(bytes) as f64
                } else {
                    <$ty>::from_le_bytes(bytes) as f64
                }
            }};
        }
        Some(match self {
            Self::Int8 => read!(i8, 1),
            Self::UInt8 => read!(u8, 1),
            Self::Int16 => read!(i16, 2),
            Self::UInt16 => read!(u16, 2),
            Self::Int32 => read!(i32, 4),
            Self::UInt32 => read!(u32, 4),
            Self::Int64 => read!(i64, 8),
            Self::UInt64 => read!(u64, 8),
            Self::Float32 => read!(f32, 4),
            Self::Float64 => read!(f64, 8),
        })
    }
}

/// How the values of a `DataArray` are stored (`format` attribute).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArrayFormat {
    Ascii,
    Binary,
    Appended,
}

/// Width of the unsigned length header in front of binary data
/// (`header_type` attribute of `VTKFile`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HeaderType {
    UInt32,
    UInt64,
}

impl HeaderType {
    fn byte_size(self) -> usize {
        match self {
            Self::UInt32 => 4,
            Self::UInt64 => 8,
        }
    }
}

/// File-wide settings from the `VTKFile` element.
#[derive(Debug, Clone)]
struct FileSettings {
    big_endian: bool,
    header_type: HeaderType,
    compressor: Option<String>,
}

/// What a `DataArray` is used for, derived from its parent element and name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArraySlot {
    Points,
    Connectivity,
    Offsets,
    Types,
    PointData,
    CellData,
    /// Not needed (e.g. `FieldData`, `faces`); its content is skipped.
    Ignored,
}

/// Parsed attributes of a `DataArray` whose content is being collected.
#[derive(Debug, Clone)]
struct PendingArray {
    slot: ArraySlot,
    name: Option<String>,
    scalar: ScalarType,
    components: usize,
    format: ArrayFormat,
    tuples: Option<usize>,
}

/// Element kinds tracked on the element stack.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Context {
    VtkFile,
    Piece,
    Points,
    Cells,
    PointData,
    CellData,
    DataArray,
    Other,
}

/// Arrays and counts of the `Piece` currently being read.
#[derive(Debug, Default)]
struct PieceBuilder {
    declared_points: usize,
    declared_cells: usize,
    points: Vec<f64>,
    connectivity: Vec<f64>,
    offsets: Vec<f64>,
    types: Vec<f64>,
    point_arrays: Vec<VtuDataArray>,
    cell_arrays: Vec<VtuDataArray>,
}

/// Event-driven state machine behind [`VtuReader::from_xml`].
struct VtuParser {
    settings: FileSettings,
    seen_vtk_file: bool,
    stack: Vec<Context>,
    piece: Option<PieceBuilder>,
    pending: Option<PendingArray>,
    text: String,
    result: VtuReader,
}

impl VtuParser {
    fn new() -> Self {
        Self {
            settings: FileSettings {
                big_endian: false,
                header_type: HeaderType::UInt32,
                compressor: None,
            },
            seen_vtk_file: false,
            stack: Vec::new(),
            piece: None,
            pending: None,
            text: String::new(),
            result: VtuReader {
                points: Vec::new(),
                num_cells: 0,
                cell_types: Vec::new(),
                connectivity: Vec::new(),
                offsets: Vec::new(),
                point_scalar_names: Vec::new(),
                cell_scalar_names: Vec::new(),
                point_arrays: Vec::new(),
                cell_arrays: Vec::new(),
            },
        }
    }

    /// Whether character data is currently part of a `DataArray` that has to
    /// be decoded.
    fn collecting_text(&self) -> bool {
        matches!(&self.pending, Some(array) if array.slot != ArraySlot::Ignored)
    }

    fn push_text(&mut self, text: &str) {
        self.text.push_str(text);
    }

    fn start_element(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let context = match element.local_name().as_ref() {
            b"VTKFile" => {
                self.read_vtk_file(element)?;
                Context::VtkFile
            }
            b"Piece" => {
                self.begin_piece(element)?;
                Context::Piece
            }
            b"Points" => Context::Points,
            b"Cells" => Context::Cells,
            b"PointData" => Context::PointData,
            b"CellData" => Context::CellData,
            b"DataArray" => {
                self.begin_array(element)?;
                Context::DataArray
            }
            _ => Context::Other,
        };
        self.stack.push(context);
        Ok(())
    }

    fn end_element(&mut self) -> Result<()> {
        match self.stack.pop() {
            Some(Context::DataArray) => self.finish_array(),
            Some(Context::Piece) => self.finish_piece(),
            _ => Ok(()),
        }
    }

    /// Parse the attributes of the `VTKFile` root element.
    fn read_vtk_file(&mut self, element: &BytesStart<'_>) -> Result<()> {
        if self.seen_vtk_file {
            return Err(Error::Parse("more than one <VTKFile> element".to_string()));
        }
        self.seen_vtk_file = true;
        let attrs = Attributes::collect(element)?;
        match attrs.get("type") {
            Some("UnstructuredGrid") => {}
            Some(other) => {
                return Err(Error::Parse(format!(
                    "unsupported VTK XML file type '{other}': only UnstructuredGrid (.vtu) is supported"
                )));
            }
            None => {
                return Err(Error::Parse(
                    "<VTKFile> element has no 'type' attribute".to_string(),
                ));
            }
        }
        self.settings.big_endian = match attrs.get("byte_order") {
            None | Some("LittleEndian") => false,
            Some("BigEndian") => true,
            Some(other) => {
                return Err(Error::Parse(format!("unknown byte_order '{other}'")));
            }
        };
        self.settings.header_type = match attrs.get("header_type") {
            None | Some("UInt32") => HeaderType::UInt32,
            Some("UInt64") => HeaderType::UInt64,
            Some(other) => {
                return Err(Error::Parse(format!("unsupported header_type '{other}'")));
            }
        };
        self.settings.compressor = attrs
            .get("compressor")
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        Ok(())
    }

    fn begin_piece(&mut self, element: &BytesStart<'_>) -> Result<()> {
        if self.piece.is_some() {
            return Err(Error::Parse("nested <Piece> elements".to_string()));
        }
        let attrs = Attributes::collect(element)?;
        let declared_points = attrs.require_count("NumberOfPoints", "Piece")?;
        let declared_cells = attrs.require_count("NumberOfCells", "Piece")?;
        self.piece = Some(PieceBuilder {
            declared_points,
            declared_cells,
            ..PieceBuilder::default()
        });
        Ok(())
    }

    fn begin_array(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let attrs = Attributes::collect(element)?;
        let name = attrs.get("Name").map(str::to_string);
        let slot = match self.stack.last() {
            Some(Context::Points) => ArraySlot::Points,
            Some(Context::Cells) => match name.as_deref() {
                Some("connectivity") => ArraySlot::Connectivity,
                Some("offsets") => ArraySlot::Offsets,
                Some("types") => ArraySlot::Types,
                _ => ArraySlot::Ignored,
            },
            Some(Context::PointData) => ArraySlot::PointData,
            Some(Context::CellData) => ArraySlot::CellData,
            _ => ArraySlot::Ignored,
        };
        // Unnamed point/cell data cannot be looked up; skip it.
        let slot = if matches!(slot, ArraySlot::PointData | ArraySlot::CellData) && name.is_none() {
            ArraySlot::Ignored
        } else {
            slot
        };
        self.text.clear();

        if slot == ArraySlot::Ignored {
            self.pending = Some(PendingArray {
                slot,
                name,
                scalar: ScalarType::Float64,
                components: 1,
                format: ArrayFormat::Ascii,
                tuples: None,
            });
            return Ok(());
        }

        let label = name.clone().unwrap_or_else(|| format!("{slot:?}"));
        let type_name = attrs
            .get("type")
            .ok_or_else(|| Error::Parse(format!("DataArray '{label}' has no 'type' attribute")))?;
        let scalar = ScalarType::parse(type_name).ok_or_else(|| {
            Error::Parse(format!(
                "DataArray '{label}' has unsupported type '{type_name}'"
            ))
        })?;
        let default_components = if slot == ArraySlot::Points { 3 } else { 1 };
        let components = match attrs.get("NumberOfComponents") {
            None => default_components,
            Some(text) => text.trim().parse::<usize>().map_err(|_| {
                Error::Parse(format!(
                    "DataArray '{label}' has invalid NumberOfComponents '{text}'"
                ))
            })?,
        };
        if components == 0 {
            return Err(Error::Parse(format!(
                "DataArray '{label}' declares zero components"
            )));
        }
        if slot == ArraySlot::Points && components != 3 {
            return Err(Error::Parse(format!(
                "Points array must have 3 components, found {components}"
            )));
        }
        let format = match attrs.get("format") {
            None | Some("ascii") => ArrayFormat::Ascii,
            Some("binary") => ArrayFormat::Binary,
            Some("appended") => ArrayFormat::Appended,
            Some(other) => {
                return Err(Error::Parse(format!(
                    "DataArray '{label}' has unknown format '{other}'"
                )));
            }
        };
        if format == ArrayFormat::Appended {
            return Err(Error::Parse(format!(
                "DataArray '{label}' uses format=\"appended\", which is not supported; \
                 write the file with ascii or inline binary data arrays"
            )));
        }
        let tuples = match attrs.get("NumberOfTuples") {
            None => None,
            Some(text) => Some(text.trim().parse::<usize>().map_err(|_| {
                Error::Parse(format!(
                    "DataArray '{label}' has invalid NumberOfTuples '{text}'"
                ))
            })?),
        };
        self.pending = Some(PendingArray {
            slot,
            name,
            scalar,
            components,
            format,
            tuples,
        });
        Ok(())
    }

    fn finish_array(&mut self) -> Result<()> {
        let Some(array) = self.pending.take() else {
            return Ok(());
        };
        if array.slot == ArraySlot::Ignored {
            self.text.clear();
            return Ok(());
        }
        let label = array
            .name
            .clone()
            .unwrap_or_else(|| format!("{:?}", array.slot));
        let values = match array.format {
            ArrayFormat::Ascii => decode_ascii(&self.text),
            ArrayFormat::Binary => {
                if let Some(compressor) = &self.settings.compressor {
                    Err(Error::Parse(format!(
                        "binary data compressed with '{compressor}' is not supported"
                    )))
                } else {
                    decode_binary(&self.text, array.scalar, &self.settings)
                }
            }
            ArrayFormat::Appended => Err(Error::Parse(
                "appended data arrays are not supported".to_string(),
            )),
        }
        .map_err(|err| match err {
            Error::Parse(message) => Error::Parse(format!("DataArray '{label}': {message}")),
            other => other,
        })?;
        self.text.clear();

        if let Some(tuples) = array.tuples
            && tuples.checked_mul(array.components) != Some(values.len())
        {
            return Err(Error::Parse(format!(
                "DataArray '{label}' declares {tuples} tuples of {} components but holds {} values",
                array.components,
                values.len()
            )));
        }
        let Some(piece) = self.piece.as_mut() else {
            return Err(Error::Parse(format!(
                "DataArray '{label}' appears outside of a <Piece>"
            )));
        };
        match array.slot {
            ArraySlot::Points => piece.points = values,
            ArraySlot::Connectivity => piece.connectivity = values,
            ArraySlot::Offsets => piece.offsets = values,
            ArraySlot::Types => piece.types = values,
            ArraySlot::PointData | ArraySlot::CellData => {
                let data = VtuDataArray {
                    name: array.name.unwrap_or_default(),
                    num_components: array.components,
                    values,
                };
                if array.slot == ArraySlot::PointData {
                    piece.point_arrays.push(data);
                } else {
                    piece.cell_arrays.push(data);
                }
            }
            ArraySlot::Ignored => {}
        }
        Ok(())
    }

    /// Validate the finished piece and append it to the result.
    fn finish_piece(&mut self) -> Result<()> {
        let Some(piece) = self.piece.take() else {
            return Ok(());
        };
        let n_points = piece.declared_points;
        let n_cells = piece.declared_cells;

        if n_points.checked_mul(3) != Some(piece.points.len()) {
            return Err(Error::Parse(format!(
                "Piece declares {n_points} points but the Points array holds {} coordinates \
                 (3 per point expected)",
                piece.points.len()
            )));
        }

        let types = piece
            .types
            .iter()
            .map(|&value| to_cell_type(value))
            .collect::<Result<Vec<u8>>>()?;
        if types.len() != n_cells {
            return Err(Error::Parse(format!(
                "Piece declares {n_cells} cells but the 'types' array has {} entries",
                types.len()
            )));
        }

        let mut offsets = piece
            .offsets
            .iter()
            .map(|&value| to_index(value, "offsets"))
            .collect::<Result<Vec<usize>>>()?;
        // Some writers store the n_cells + 1 cell boundaries starting at 0.
        if n_cells.checked_add(1) == Some(offsets.len()) && offsets.first() == Some(&0) {
            offsets.remove(0);
        }
        if offsets.len() != n_cells {
            return Err(Error::Parse(format!(
                "Piece declares {n_cells} cells but the 'offsets' array has {} entries",
                offsets.len()
            )));
        }
        let connectivity = piece
            .connectivity
            .iter()
            .map(|&value| to_index(value, "connectivity"))
            .collect::<Result<Vec<usize>>>()?;
        let mut previous = 0usize;
        for (cell, &end) in offsets.iter().enumerate() {
            if end < previous {
                return Err(Error::Parse(format!(
                    "'offsets' decreases at cell {cell} ({previous} -> {end})"
                )));
            }
            previous = end;
        }
        if previous != connectivity.len() {
            return Err(Error::Parse(format!(
                "last offset is {previous} but the connectivity array has {} entries",
                connectivity.len()
            )));
        }
        if let Some(&bad) = connectivity.iter().find(|&&index| index >= n_points) {
            return Err(Error::Parse(format!(
                "cell connectivity refers to point {bad} but the piece has {n_points} points"
            )));
        }

        let points_before = self.result.points.len();
        let connectivity_before = self.result.connectivity.len();
        let cells_before = self.result.cell_types.len();

        self.result.points.extend(
            piece
                .points
                .chunks_exact(3)
                .map(|xyz| [xyz[0], xyz[1], xyz[2]]),
        );
        self.result
            .connectivity
            .extend(connectivity.iter().map(|&index| index + points_before));
        self.result
            .offsets
            .extend(offsets.iter().map(|&end| end + connectivity_before));
        self.result.cell_types.extend_from_slice(&types);
        self.result.num_cells += n_cells;

        merge_arrays(
            &mut self.result.point_arrays,
            piece.point_arrays,
            points_before,
            n_points,
            "point data",
        )?;
        merge_arrays(
            &mut self.result.cell_arrays,
            piece.cell_arrays,
            cells_before,
            n_cells,
            "cell data",
        )?;
        self.result.point_scalar_names = self
            .result
            .point_arrays
            .iter()
            .map(|array| array.name.clone())
            .collect();
        self.result.cell_scalar_names = self
            .result
            .cell_arrays
            .iter()
            .map(|array| array.name.clone())
            .collect();
        Ok(())
    }

    fn finish(self) -> Result<VtuReader> {
        if !self.seen_vtk_file {
            return Err(Error::Parse(
                "not a VTK XML file: no <VTKFile> element found".to_string(),
            ));
        }
        if !self.stack.is_empty() {
            return Err(Error::Parse(
                "unexpected end of document: unclosed XML elements".to_string(),
            ));
        }
        Ok(self.result)
    }
}

/// Append the data arrays of a new piece to the arrays collected so far.
///
/// `tuples_before` is the number of points (or cells) of the preceding pieces,
/// `tuples` the number in the new piece.  Every array must cover all pieces.
fn merge_arrays(
    target: &mut Vec<VtuDataArray>,
    incoming: Vec<VtuDataArray>,
    tuples_before: usize,
    tuples: usize,
    what: &str,
) -> Result<()> {
    for array in incoming {
        match target.iter_mut().find(|known| known.name == array.name) {
            Some(known) => {
                if known.num_components != array.num_components {
                    return Err(Error::Parse(format!(
                        "{what} array '{}' has {} components in one piece and {} in another",
                        array.name, known.num_components, array.num_components
                    )));
                }
                known.values.extend(array.values);
            }
            None => {
                if tuples_before > 0 {
                    return Err(Error::Parse(format!(
                        "{what} array '{}' is missing from the first piece(s)",
                        array.name
                    )));
                }
                target.push(array);
            }
        }
    }
    for known in target.iter() {
        let expected = tuples_before
            .checked_add(tuples)
            .and_then(|total| total.checked_mul(known.num_components));
        if expected != Some(known.values.len()) {
            return Err(Error::Parse(format!(
                "{what} array '{}' holds {} values, which does not match {} tuples of {} \
                 components (missing from a piece, or wrong tuple count)",
                known.name,
                known.values.len(),
                tuples_before.saturating_add(tuples),
                known.num_components
            )));
        }
    }
    Ok(())
}

/// Attributes of one XML start tag as owned `(name, value)` pairs.
struct Attributes {
    pairs: Vec<(String, String)>,
}

impl Attributes {
    fn collect(element: &BytesStart<'_>) -> Result<Self> {
        let mut pairs = Vec::new();
        for attribute in element.attributes() {
            let attribute =
                attribute.map_err(|err| Error::Parse(format!("malformed XML attribute: {err}")))?;
            let value = attribute
                .normalized_value(XmlVersion::Implicit1_0)
                .map_err(|err| Error::Parse(format!("invalid attribute value: {err}")))?;
            pairs.push((
                String::from_utf8_lossy(attribute.key.as_ref()).into_owned(),
                value.into_owned(),
            ));
        }
        Ok(Self { pairs })
    }

    fn get(&self, key: &str) -> Option<&str> {
        self.pairs
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// A required non-negative integer attribute.
    fn require_count(&self, key: &str, element: &str) -> Result<usize> {
        let text = self
            .get(key)
            .ok_or_else(|| Error::Parse(format!("<{element}> has no '{key}' attribute")))?;
        text.trim().parse::<usize>().map_err(|_| {
            Error::Parse(format!(
                "<{element}> attribute {key}=\"{text}\" is not a non-negative integer"
            ))
        })
    }
}

/// Convert a decoded value to a non-negative integer index.
fn to_index(value: f64, what: &str) -> Result<usize> {
    if value.is_finite() && value >= 0.0 && value.fract() == 0.0 && value < 9.007_199_254_740_992e15
    {
        Ok(value as usize)
    } else {
        Err(Error::Parse(format!(
            "'{what}' contains the value {value}, which is not a valid index"
        )))
    }
}

/// Convert a decoded value to a VTK cell type id.
fn to_cell_type(value: f64) -> Result<u8> {
    if value.is_finite() && (0.0..=255.0).contains(&value) && value.fract() == 0.0 {
        Ok(value as u8)
    } else {
        Err(Error::Parse(format!(
            "'types' contains the value {value}, which is not a VTK cell type"
        )))
    }
}

/// Decode a whitespace separated list of numbers.
fn decode_ascii(text: &str) -> Result<Vec<f64>> {
    text.split_ascii_whitespace()
        .map(|token| {
            token.parse::<f64>().map_err(|_| {
                Error::Parse(format!("'{token}' is not a number in an ascii DataArray"))
            })
        })
        .collect()
}

/// Decode an inline base64 `binary` array into values.
fn decode_binary(text: &str, scalar: ScalarType, settings: &FileSettings) -> Result<Vec<f64>> {
    let bytes = binary_payload(text, settings)?;
    let size = scalar.byte_size();
    if !bytes.len().is_multiple_of(size) {
        return Err(Error::Parse(format!(
            "binary payload of {} bytes is not a multiple of the element size {size}",
            bytes.len()
        )));
    }
    bytes
        .chunks_exact(size)
        .map(|chunk| {
            scalar
                .decode(chunk, settings.big_endian)
                .ok_or_else(|| Error::Parse("truncated binary element".to_string()))
        })
        .collect()
}

/// Extract the payload bytes of an uncompressed inline binary array and check
/// them against the length header.
///
/// VTK base64-encodes the header and the payload separately (so the header
/// carries its own padding); a single contiguous base64 stream is accepted as
/// well.  Either way the decoded payload must be exactly as long as the header
/// states.
fn binary_payload(text: &str, settings: &FileSettings) -> Result<Vec<u8>> {
    let compact: String = text.chars().filter(|c| !c.is_whitespace()).collect();
    if compact.is_empty() {
        return Ok(Vec::new());
    }
    if !compact
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'+' | b'/' | b'='))
    {
        return Err(Error::Parse(
            "invalid base64 data: unexpected character".to_string(),
        ));
    }
    let header_len = settings.header_type.byte_size();
    let header_chars = header_len.div_ceil(3) * 4;

    // Layout 1 (VTK): base64(header) immediately followed by base64(payload).
    if compact.len() >= header_chars
        && let Some(header_bytes) = base64_decode_buffer(&compact[..header_chars])
        && header_bytes.len() == header_len
        && let Some(declared) = read_header(&header_bytes, settings)
        && let Some(payload) = base64_decode_buffer(&compact[header_chars..])
        && payload.len() == declared
    {
        return Ok(payload);
    }

    // Layout 2: header and payload in one base64 stream.
    if let Some(all) = base64_decode_buffer(&compact)
        && all.len() >= header_len
        && let Some(declared) = read_header(&all[..header_len], settings)
        && all.len() - header_len == declared
    {
        return Ok(all[header_len..].to_vec());
    }
    Err(Error::Parse(
        "binary data does not match its length header (truncated, or compressed data \
         without a 'compressor' attribute)"
            .to_string(),
    ))
}

/// Read the unsigned length header.
fn read_header(bytes: &[u8], settings: &FileSettings) -> Option<usize> {
    let value = match settings.header_type {
        HeaderType::UInt32 => {
            let raw = *bytes.first_chunk::<4>()?;
            u64::from(if settings.big_endian {
                u32::from_be_bytes(raw)
            } else {
                u32::from_le_bytes(raw)
            })
        }
        HeaderType::UInt64 => {
            let raw = *bytes.first_chunk::<8>()?;
            if settings.big_endian {
                u64::from_be_bytes(raw)
            } else {
                u64::from_le_bytes(raw)
            }
        }
    };
    usize::try_from(value).ok()
}

#[cfg(test)]
mod tests;
