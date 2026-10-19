// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The finite-element, Fluent and OpenFOAM readers must meet damaged input
//! with an error, or with whatever they can still read, and never with a
//! panic or a hang.
//!
//! Every reader is fed thousands of mutations (deleted, duplicated and
//! overwritten byte ranges, inserted delimiters, truncations) of a valid
//! file.  The mutations come from a fixed seed, so a failure reproduces.

use std::panic::{AssertUnwindSafe, catch_unwind};

use oxiphysics_io::finite_element_io::{AbaqusInpIo, GmshMesh, MapdlDeck, NastranBulkData};
use oxiphysics_io::fluent_format::FluentReader;
use oxiphysics_io::foam_io::{parse_foam_scalar_field_internal, parse_foam_vector_field};

/// How many mutations of each seed file each reader gets, unless the
/// environment variable [`ROUNDS_VARIABLE`] asks for another number (for long
/// soak runs).
const DEFAULT_MUTATIONS_PER_SEED: usize = 4000;

/// Environment variable that overrides [`DEFAULT_MUTATIONS_PER_SEED`].
const ROUNDS_VARIABLE: &str = "OXIPHYSICS_MESH_FUZZ_ROUNDS";

fn mutations_per_seed() -> usize {
    std::env::var(ROUNDS_VARIABLE)
        .ok()
        .and_then(|text| text.trim().parse().ok())
        .unwrap_or(DEFAULT_MUTATIONS_PER_SEED)
}

/// Bytes that change the structure of the files under test.
const DELIMITERS: &[u8] = b"()\"{};,$*+-.eE0123456789 \n\t\0\xff/abcxyz";

/// xorshift64*: small, deterministic and good enough for choosing edits.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, bound: usize) -> usize {
        let bound = u64::try_from(bound.max(1)).unwrap_or(1);
        usize::try_from(self.next() % bound).unwrap_or(0)
    }
}

/// `seed` with one to four random edits applied.
fn mutate(seed: &[u8], rng: &mut Rng) -> Vec<u8> {
    let mut data = seed.to_vec();
    for _ in 0..=rng.below(4) {
        let len = data.len();
        match rng.below(6) {
            0 => data.truncate(rng.below(len + 1)),
            1 => {
                let start = rng.below(len + 1);
                let end = (start + rng.below(24)).min(len);
                data.drain(start..end);
            }
            2 => {
                let start = rng.below(len + 1);
                let end = (start + rng.below(24)).min(len);
                let copy = data[start..end].to_vec();
                let at = rng.below(len + 1);
                data.splice(at..at, copy);
            }
            3 => {
                if len > 0 {
                    let at = rng.below(len);
                    data[at] = DELIMITERS[rng.below(DELIMITERS.len())];
                }
            }
            4 => {
                let at = rng.below(len + 1);
                data.insert(at, DELIMITERS[rng.below(DELIMITERS.len())]);
            }
            _ => {
                let start = rng.below(len + 1);
                let end = (start + rng.below(64)).min(len);
                let fill = DELIMITERS[rng.below(DELIMITERS.len())];
                data[start..end].fill(fill);
            }
        }
    }
    data
}

/// Feed `reader` the seed itself and many mutations of it.
fn hammer(name: &str, seed: &str, reader: impl Fn(&[u8])) {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ seed.len() as u64);
    reader(seed.as_bytes());
    for round in 0..mutations_per_seed() {
        let input = mutate(seed.as_bytes(), &mut rng);
        if catch_unwind(AssertUnwindSafe(|| reader(&input))).is_err() {
            panic!(
                "{name} panicked on mutation {round}: {:?}",
                String::from_utf8_lossy(&input)
            );
        }
    }
}

fn as_text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

const CDB: &str = "\
/PREP7
ET,       1,185
ET,       2,187
MP,EX,1,2.0E11
NBLOCK,6,SOLID,         4,         4
(3i9,6e21.13e3)
        1        0        0 0.0000000000000E+000 0.0000000000000E+000 0.0000000000000E+000
        2        0        0 1.0000000000000E+000-2.0000000000000E+000 3.0000000000000E+000
        3        0        0-4.0000000000000E+001 5.0000000000000E+000 6.0000000000000E+002
        4        0        0 0.0000000000000E+000 0.0000000000000E+000 1.0000000000000E+000
-1
EBLOCK,19,SOLID,         2,         2
(19i9)
        1        1        1        1        0        0        0        0        4        0        1        1        2        3        4
        1        2        1        1        0        0        0        0       10        0        2        1        2        3        4        5        6        7        8
        9       10
-1
D,2,UY,0.5
F,4,FZ,-100.0
";

const GMSH_V4: &str = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
3 1 \"cube\"
$EndPhysicalNames
$Entities
0 0 0 1
1 0 0 0 1 1 1 0 0
$EndEntities
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 2 1 2
3 1 4 2
1 1 2 3 4
2 4 3 2 1
$EndElements
";

const GMSH_V2: &str = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
2
1 4 2 11 7 1 2 3 4
2 2 0 1 2 3
$EndElements
";

const NASTRAN: &str = "\
$ Nastran bulk data
SOL 101
CEND
SPC = 1
BEGIN BULK
GRID    1               0.0     0.0     0.0
GRID*   2                               1.0             0.0
*       0.0
GRID,3,0,0.0,1.0,0.0
CTETRA,1,1,1,2,3,4
CHEXA   1       1       1       2       3       4       5       6       +C1
+C1     7       8
CHEXA,2,3,1,2,3,4,5,6,+
+,7,8,9,10,11,12,13,14,+
+,15,16,17,18,19,20
MAT1,1,200E9,77E9,0.3,7800
SPC,1,1,123,0.0
FORCE,1,4,0,100.0,0.0,0.0,1.0
ENDDATA
";

const ABAQUS: &str = "\
*HEADING
cantilever
*NODE, NSET=ALL
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0
*ELEMENT, TYPE=C3D8R, ELSET=BEAM
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D20
2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
16, 17, 18, 19, 20
*NSET, NSET=FIXED, GENERATE
1, 10, 3
*NSET, NSET=TIP
2, 3
*ELSET, ELSET=ALL, GENERATE
1, 2
*MATERIAL, NAME=Steel
*ELASTIC
200.0E9, 0.3
*BOUNDARY
FIXED, ENCASTRE
2, 1, 3, 0.5
*STEP, NAME=Load
*STATIC
*CLOAD
TIP, 3, -250.0
*END STEP
";

const FLUENT_2D: &str = "\
(0 \"2-D mesh (two triangles)\")
(2 2)
(10 (0 1 4 0 2))
(10 (1 1 4 1 2)(
0.0 0.0
1.0 0.0
1.0 1.0
0.0 1.0
))
(12 (0 1 2 0))
(12 (3 1 2 1 1))
(13 (0 1 5 0))
(13 (4 1 1 2 2)(
1 3 1 2
))
(13 (5 2 5 3 0)(
2 1 2 1 0
2 2 3 1 0
2 3 4 2 0
2 4 1 2 0
))
(12 (7 3 5 1 0)(2 5 6))
(45 (3 fluid fluid-3)())
(39 (5 wall wall-5 1)((a . 1)))
";

const FOAM_VELOCITY: &str = "\
/* banner (with parentheses) */
FoamFile
{
    format      ascii;
    class       volVectorField;
    object      U;
}
dimensions      [0 1 -1 0 0 0 0];
internalField   nonuniform List<vector>
3
(
(0.5 0 0)
(1 -0.25 0)
(1.5e-1 0 2)
)
;
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    walls { type fixedValue; value nonuniform List<vector> 2((0 0 0) (7 8 9)); }
}
";

const FOAM_PRESSURE: &str = "\
FoamFile { class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField   nonuniform List<scalar> 3(1 2.5 -3);
boundaryField
{
    inlet { type fixedValue; value nonuniform List<scalar> 2(5 6); }
    outlet { type fixedValue; value $internalField; }
}
";

#[test]
fn mapdl_reader_survives_damaged_decks() {
    hammer("MapdlDeck::parse", CDB, |bytes| {
        let text = as_text(bytes);
        let _ = MapdlDeck::parse(&text);
        let _ = MapdlDeck::try_parse(&text);
    });
}

#[test]
fn gmsh_reader_survives_damaged_files() {
    for (name, seed) in [("Gmsh 4.1", GMSH_V4), ("Gmsh 2.2", GMSH_V2)] {
        hammer(name, seed, |bytes| {
            let _ = GmshMesh::parse(&as_text(bytes));
        });
    }
}

#[test]
fn nastran_reader_survives_damaged_bulk_data() {
    hammer("NastranBulkData::parse", NASTRAN, |bytes| {
        let data = NastranBulkData::parse(&as_text(bytes));
        // Whatever was read must be writable again.
        let _ = data.write();
    });
}

#[test]
fn abaqus_reader_survives_damaged_input_files() {
    hammer("AbaqusInpIo", ABAQUS, |bytes| {
        let text = as_text(bytes);
        let mesh = AbaqusInpIo::parse(&text);
        let _ = AbaqusInpIo::write(&mesh);
        let _ = AbaqusInpIo::try_parse(&text);
    });
}

#[test]
fn fluent_reader_survives_damaged_meshes() {
    hammer("FluentReader::parse", FLUENT_2D, |bytes| {
        let _ = FluentReader::parse(bytes);
    });
}

#[test]
fn openfoam_field_readers_survive_damaged_files() {
    hammer("parse_foam_vector_field", FOAM_VELOCITY, |bytes| {
        let _ = parse_foam_vector_field(&as_text(bytes));
        let _ = parse_foam_scalar_field_internal(&as_text(bytes));
    });
    hammer("parse_foam_scalar_field_internal", FOAM_PRESSURE, |bytes| {
        let _ = parse_foam_scalar_field_internal(&as_text(bytes));
        let _ = parse_foam_vector_field(&as_text(bytes));
    });
}

#[test]
fn the_undamaged_seed_files_are_valid() {
    let deck = MapdlDeck::try_parse(CDB).expect("the CDB seed is valid");
    assert_eq!(deck.nodes.len(), 4);
    assert_eq!(deck.elements.len(), 2);
    assert_eq!(deck.nodes[1].coords, [1.0, -2.0, 3.0]);
    assert!(GmshMesh::parse(GMSH_V4).is_ok());
    assert!(GmshMesh::parse(GMSH_V2).is_ok());
    let nastran = NastranBulkData::parse(NASTRAN);
    assert_eq!(nastran.elements.len(), 3);
    assert_eq!(nastran.grids.len(), 3);
    let abaqus = AbaqusInpIo::try_parse(ABAQUS).expect("the Abaqus seed is valid");
    assert_eq!(abaqus.elements.len(), 2);
    let fluent = FluentReader::parse(FLUENT_2D.as_bytes()).expect("the Fluent seed is valid");
    assert_eq!(fluent.node_count(), 4);
    assert_eq!(fluent.cell_count(), 5);
    assert_eq!(
        parse_foam_vector_field(FOAM_VELOCITY)
            .expect("the U seed is valid")
            .len(),
        3
    );
    assert_eq!(
        parse_foam_scalar_field_internal(FOAM_PRESSURE).expect("the p seed is valid"),
        vec![1.0, 2.5, -3.0]
    );
}
