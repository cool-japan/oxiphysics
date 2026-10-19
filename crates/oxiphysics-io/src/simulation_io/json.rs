// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! JSON encoding helpers for the scene and checkpoint documents.
//!
//! Floating-point numbers are written with the shortest decimal form that
//! parses back to the identical `f64`, and the crate builds `serde_json` with
//! its `float_roundtrip` feature, so state survives a write/read cycle
//! bit-exactly (including `-0.0` and subnormals).
//!
//! JSON has no spelling for NaN and the infinities, so they are written as the
//! strings `"NaN"`, `"Infinity"` and `"-Infinity"` and accepted in the same
//! form on read; this is the convention of the restart-file JSON format. Without
//! it a checkpoint of a diverged simulation would be written as `null` and could
//! never be read back.

use std::fmt;

use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::{SerializeMap, SerializeSeq, Serializer};
use serde::{Deserialize, Serialize};

/// Encode `value` as compact JSON.
pub(super) fn encode<T: Serialize>(value: &T) -> std::io::Result<String> {
    serde_json::to_string(value).map_err(std::io::Error::from)
}

fn serialize_f64<S: Serializer>(value: f64, serializer: S) -> Result<S::Ok, S::Error> {
    if value.is_finite() {
        serializer.serialize_f64(value)
    } else if value.is_nan() {
        serializer.serialize_str("NaN")
    } else if value > 0.0 {
        serializer.serialize_str("Infinity")
    } else {
        serializer.serialize_str("-Infinity")
    }
}

struct F64Visitor;

impl Visitor<'_> for F64Visitor {
    type Value = f64;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("a number, \"NaN\", \"Infinity\" or \"-Infinity\"")
    }

    fn visit_f64<E: de::Error>(self, value: f64) -> Result<f64, E> {
        Ok(value)
    }

    fn visit_i64<E: de::Error>(self, value: i64) -> Result<f64, E> {
        Ok(value as f64)
    }

    fn visit_u64<E: de::Error>(self, value: u64) -> Result<f64, E> {
        Ok(value as f64)
    }

    fn visit_str<E: de::Error>(self, value: &str) -> Result<f64, E> {
        match value {
            "NaN" => Ok(f64::NAN),
            "Infinity" | "inf" => Ok(f64::INFINITY),
            "-Infinity" | "-inf" => Ok(f64::NEG_INFINITY),
            other => Err(E::invalid_value(de::Unexpected::Str(other), &self)),
        }
    }
}

/// An `f64` that is encoded with the non-finite convention described in the
/// module documentation.
pub(super) struct SpecialF64(pub(super) f64);

impl Serialize for SpecialF64 {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serialize_f64(self.0, serializer)
    }
}

impl<'de> Deserialize<'de> for SpecialF64 {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(F64Visitor).map(SpecialF64)
    }
}

/// `#[serde(with = ...)]` adapter for a single `f64`.
pub(super) mod f64_value {
    use super::*;

    pub fn serialize<S: Serializer>(value: &f64, serializer: S) -> Result<S::Ok, S::Error> {
        serialize_f64(*value, serializer)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<f64, D::Error> {
        deserializer.deserialize_any(F64Visitor)
    }
}

/// `#[serde(with = ...)]` adapter for a fixed-size `f64` array.
///
/// Deserialization requires exactly `N` elements.
pub(super) mod f64_array {
    use super::*;

    pub fn serialize<S: Serializer, const N: usize>(
        values: &[f64; N],
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(N))?;
        for &value in values {
            seq.serialize_element(&SpecialF64(value))?;
        }
        seq.end()
    }

    struct ArrayVisitor<const N: usize>;

    impl<'de, const N: usize> Visitor<'de> for ArrayVisitor<N> {
        type Value = [f64; N];

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(formatter, "an array of {N} numbers")
        }

        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<[f64; N], A::Error> {
            let mut out = [0.0f64; N];
            for (index, slot) in out.iter_mut().enumerate() {
                let element = seq.next_element::<SpecialF64>()?;
                match element {
                    Some(SpecialF64(value)) => *slot = value,
                    None => return Err(<A::Error as de::Error>::invalid_length(index, &self)),
                }
            }
            if seq.next_element::<de::IgnoredAny>()?.is_some() {
                return Err(<A::Error as de::Error>::invalid_length(N + 1, &self));
            }
            Ok(out)
        }
    }

    pub fn deserialize<'de, D: Deserializer<'de>, const N: usize>(
        deserializer: D,
    ) -> Result<[f64; N], D::Error> {
        deserializer.deserialize_seq(ArrayVisitor::<N>)
    }
}

/// `#[serde(with = ...)]` adapter for an ordered list of string pairs, stored
/// as a JSON object.
///
/// Going through a map type would sort the keys and drop duplicates; this
/// keeps the entries in file order.
pub(super) mod ordered_string_map {
    use super::*;

    pub fn serialize<S: Serializer>(
        entries: &[(String, String)],
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        let mut map = serializer.serialize_map(Some(entries.len()))?;
        for (key, value) in entries {
            map.serialize_entry(key, value)?;
        }
        map.end()
    }

    struct MapVisitor;

    impl<'de> Visitor<'de> for MapVisitor {
        type Value = Vec<(String, String)>;

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("an object whose values are strings")
        }

        fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
            let mut entries = Vec::with_capacity(map.size_hint().unwrap_or(0).min(1024));
            while let Some((key, value)) = map.next_entry::<String, String>()? {
                entries.push((key, value));
            }
            Ok(entries)
        }
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Vec<(String, String)>, D::Error> {
        deserializer.deserialize_map(MapVisitor)
    }
}

/// Borrowed metadata list that serializes as a JSON object (key order kept).
pub(super) struct MetadataRef<'a>(pub(super) &'a [(String, String)]);

impl Serialize for MetadataRef<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        ordered_string_map::serialize(self.0, serializer)
    }
}
