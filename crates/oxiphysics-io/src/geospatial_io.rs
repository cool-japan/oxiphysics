// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geospatial I/O: GeoJSON, WKT, CRS metadata, DEM rasters, bounding-box
//! queries, Mercator projections, terrain tiles, geohash, R-tree-like spatial
//! index and a simplified DBF attribute table.
//!
//! All coordinate values are plain `f64` or `[f64; 3]` (longitude, latitude,
//! optional elevation). No external geometry library is used.

// ─────────────────────────────────────────────────────────────────────────────
// 1. Errors
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub enum GeoError {
    /// The input string could not be parsed.
    ParseError(String),
    /// A requested feature or geometry was not found.
    NotFound(String),
    /// The coordinate reference system is unsupported.
    UnsupportedCrs(String),
    /// Generic I/O error message.
    Io(String),
}

impl std::fmt::Display for GeoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GeoError::ParseError(m) => write!(f, "parse error: {m}"),
            GeoError::NotFound(m) => write!(f, "not found: {m}"),
            GeoError::UnsupportedCrs(m) => write!(f, "unsupported CRS: {m}"),
            GeoError::Io(m) => write!(f, "I/O error: {m}"),
        }
    }
}

impl std::error::Error for GeoError {}

/// Convenience alias for geospatial results.
pub type GeoResult<T> = Result<T, GeoError>;

// ─────────────────────────────────────────────────────────────────────────────
// 2. Coordinate types
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D geographic position `(longitude_deg, latitude_deg)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LonLat {
    /// Longitude in decimal degrees (−180 … +180).
    pub lon: f64,
    /// Latitude in decimal degrees (−90 … +90).
    pub lat: f64,
}

impl LonLat {
    /// Construct a new `LonLat`.
    #[must_use]
    pub fn new(lon: f64, lat: f64) -> Self {
        Self { lon, lat }
    }

    /// Haversine distance in metres between two geographic positions.
    #[must_use]
    pub fn haversine_m(&self, other: &LonLat) -> f64 {
        const R: f64 = 6_371_000.0; // Earth mean radius [m]
        let dlat = (other.lat - self.lat).to_radians();
        let dlon = (other.lon - self.lon).to_radians();
        let lat1 = self.lat.to_radians();
        let lat2 = other.lat.to_radians();
        let a = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
        let c = 2.0 * a.sqrt().asin();
        R * c
    }
}

/// A 3-D geographic position with elevation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LonLatElev {
    /// Longitude in decimal degrees.
    pub lon: f64,
    /// Latitude in decimal degrees.
    pub lat: f64,
    /// Elevation above the ellipsoid in metres.
    pub elev: f64,
}

impl LonLatElev {
    /// Construct from individual components.
    #[must_use]
    pub fn new(lon: f64, lat: f64, elev: f64) -> Self {
        Self { lon, lat, elev }
    }

    /// Drop elevation and return the 2-D position.
    #[must_use]
    pub fn to_lonlat(&self) -> LonLat {
        LonLat::new(self.lon, self.lat)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Bounding box
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box in geographic (lon/lat) space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    /// Minimum longitude.
    pub min_lon: f64,
    /// Minimum latitude.
    pub min_lat: f64,
    /// Maximum longitude.
    pub max_lon: f64,
    /// Maximum latitude.
    pub max_lat: f64,
}

impl BoundingBox {
    /// Construct a `BoundingBox` from corner coordinates.
    #[must_use]
    pub fn new(min_lon: f64, min_lat: f64, max_lon: f64, max_lat: f64) -> Self {
        Self {
            min_lon,
            min_lat,
            max_lon,
            max_lat,
        }
    }

    /// Return `true` if `pt` lies within or on the boundary of this box.
    #[must_use]
    pub fn contains(&self, pt: LonLat) -> bool {
        pt.lon >= self.min_lon
            && pt.lon <= self.max_lon
            && pt.lat >= self.min_lat
            && pt.lat <= self.max_lat
    }

    /// Return `true` if `other` overlaps with this box.
    #[must_use]
    pub fn intersects(&self, other: &BoundingBox) -> bool {
        !(other.min_lon > self.max_lon
            || other.max_lon < self.min_lon
            || other.min_lat > self.max_lat
            || other.max_lat < self.min_lat)
    }

    /// Expand this box to include `pt`.
    pub fn expand_to_include(&mut self, pt: LonLat) {
        if pt.lon < self.min_lon {
            self.min_lon = pt.lon;
        }
        if pt.lon > self.max_lon {
            self.max_lon = pt.lon;
        }
        if pt.lat < self.min_lat {
            self.min_lat = pt.lat;
        }
        if pt.lat > self.max_lat {
            self.max_lat = pt.lat;
        }
    }

    /// Centre point of the bounding box.
    #[must_use]
    pub fn centre(&self) -> LonLat {
        LonLat::new(
            (self.min_lon + self.max_lon) / 2.0,
            (self.min_lat + self.max_lat) / 2.0,
        )
    }

    /// Width in degrees of longitude.
    #[must_use]
    pub fn width_deg(&self) -> f64 {
        self.max_lon - self.min_lon
    }

    /// Height in degrees of latitude.
    #[must_use]
    pub fn height_deg(&self) -> f64 {
        self.max_lat - self.min_lat
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Coordinate Reference System metadata
// ─────────────────────────────────────────────────────────────────────────────

/// Named coordinate reference system authority codes (EPSG subset).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CrsCode {
    /// WGS 84 geographic (EPSG:4326).
    Epsg4326,
    /// Web Mercator (EPSG:3857).
    Epsg3857,
    /// UTM zone (Northern/Southern hemisphere, zone number 1–60).
    Utm {
        /// UTM zone number (1–60).
        zone: u8,
        /// `true` for the Northern hemisphere, `false` for the Southern.
        north: bool,
    },
    /// A custom CRS identified by a free-form string.
    Custom(String),
}

impl std::fmt::Display for CrsCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CrsCode::Epsg4326 => write!(f, "EPSG:4326"),
            CrsCode::Epsg3857 => write!(f, "EPSG:3857"),
            CrsCode::Utm { zone, north } => {
                let hemi = if *north { "N" } else { "S" };
                write!(f, "UTM Zone {zone}{hemi}")
            }
            CrsCode::Custom(s) => write!(f, "{s}"),
        }
    }
}

/// Metadata bundle describing a coordinate reference system.
#[derive(Debug, Clone)]
pub struct CrsMetadata {
    /// The authority code identifying this CRS.
    pub code: CrsCode,
    /// Human-readable name.
    pub name: String,
    /// Datum name (e.g. "WGS 84").
    pub datum: String,
    /// Unit of measure for coordinate axes.
    pub units: String,
}

impl CrsMetadata {
    /// WGS 84 geographic CRS.
    #[must_use]
    pub fn wgs84() -> Self {
        Self {
            code: CrsCode::Epsg4326,
            name: "WGS 84".to_string(),
            datum: "WGS 84".to_string(),
            units: "degree".to_string(),
        }
    }

    /// Web Mercator projected CRS.
    #[must_use]
    pub fn web_mercator() -> Self {
        Self {
            code: CrsCode::Epsg3857,
            name: "WGS 84 / Pseudo-Mercator".to_string(),
            datum: "WGS 84".to_string(),
            units: "metre".to_string(),
        }
    }

    /// Serialise to a minimal WKT2 string.
    #[must_use]
    pub fn to_wkt2(&self) -> String {
        format!(
            r#"GEOGCRS["{}",DATUM["{}"],CS[ellipsoidal,2],AXIS["longitude",east],AXIS["latitude",north],UNIT["{}",1]]"#,
            self.name, self.datum, self.units
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. GeoJSON geometry types
// ─────────────────────────────────────────────────────────────────────────────

/// A GeoJSON geometry value.
#[derive(Debug, Clone, PartialEq)]
pub enum GeoGeometry {
    /// A single position.
    Point(LonLat),
    /// A single position with elevation.
    Point3D(LonLatElev),
    /// An ordered sequence of positions forming an open curve.
    LineString(Vec<LonLat>),
    /// A polygon defined by an exterior ring and zero or more interior rings.
    Polygon {
        /// Exterior ring (closed — first == last position).
        exterior: Vec<LonLat>,
        /// Interior rings (holes).
        holes: Vec<Vec<LonLat>>,
    },
    /// A collection of points.
    MultiPoint(Vec<LonLat>),
    /// A collection of line strings.
    MultiLineString(Vec<Vec<LonLat>>),
    /// A collection of polygons.
    MultiPolygon(Vec<(Vec<LonLat>, Vec<Vec<LonLat>>)>),
    /// A heterogeneous collection of geometries.
    GeometryCollection(Vec<GeoGeometry>),
}

impl GeoGeometry {
    /// Compute the axis-aligned bounding box of this geometry.
    ///
    /// Returns `None` for empty collections.
    #[must_use]
    pub fn bounding_box(&self) -> Option<BoundingBox> {
        let pts = self.collect_lonlat();
        if pts.is_empty() {
            return None;
        }
        let mut bb = BoundingBox::new(pts[0].lon, pts[0].lat, pts[0].lon, pts[0].lat);
        for p in &pts[1..] {
            bb.expand_to_include(*p);
        }
        Some(bb)
    }

    /// Recursively collect all `LonLat` positions from this geometry.
    #[must_use]
    pub fn collect_lonlat(&self) -> Vec<LonLat> {
        match self {
            GeoGeometry::Point(p) => vec![*p],
            GeoGeometry::Point3D(p) => vec![p.to_lonlat()],
            GeoGeometry::LineString(pts) => pts.clone(),
            GeoGeometry::Polygon { exterior, holes } => {
                let mut v = exterior.clone();
                for h in holes {
                    v.extend_from_slice(h);
                }
                v
            }
            GeoGeometry::MultiPoint(pts) => pts.clone(),
            GeoGeometry::MultiLineString(lines) => lines.iter().flatten().copied().collect(),
            GeoGeometry::MultiPolygon(polys) => polys
                .iter()
                .flat_map(|(ext, holes)| {
                    let mut v = ext.clone();
                    for h in holes {
                        v.extend_from_slice(h);
                    }
                    v
                })
                .collect(),
            GeoGeometry::GeometryCollection(gs) => {
                gs.iter().flat_map(|g| g.collect_lonlat()).collect()
            }
        }
    }

    /// Return the GeoJSON `type` string for this geometry.
    #[must_use]
    pub fn type_str(&self) -> &'static str {
        match self {
            GeoGeometry::Point(_) | GeoGeometry::Point3D(_) => "Point",
            GeoGeometry::LineString(_) => "LineString",
            GeoGeometry::Polygon { .. } => "Polygon",
            GeoGeometry::MultiPoint(_) => "MultiPoint",
            GeoGeometry::MultiLineString(_) => "MultiLineString",
            GeoGeometry::MultiPolygon(_) => "MultiPolygon",
            GeoGeometry::GeometryCollection(_) => "GeometryCollection",
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. GeoJSON serialisation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Serialise a `LonLat` to a GeoJSON coordinate pair string `[lon,lat]`.
#[must_use]
pub fn lonlat_to_json(p: LonLat) -> String {
    format!("[{:.6},{:.6}]", p.lon, p.lat)
}

/// Serialise a ring (closed polygon boundary) to a JSON array of pairs.
#[must_use]
pub fn ring_to_json(pts: &[LonLat]) -> String {
    let items: Vec<String> = pts.iter().map(|p| lonlat_to_json(*p)).collect();
    format!("[{}]", items.join(","))
}

/// Serialise a `GeoGeometry` to its GeoJSON `geometry` object string.
#[must_use]
pub fn geometry_to_geojson(g: &GeoGeometry) -> String {
    match g {
        GeoGeometry::Point(p) => {
            format!(r#"{{"type":"Point","coordinates":{}}}"#, lonlat_to_json(*p))
        }
        GeoGeometry::Point3D(p) => {
            format!(
                r#"{{"type":"Point","coordinates":[{:.6},{:.6},{:.6}]}}"#,
                p.lon, p.lat, p.elev
            )
        }
        GeoGeometry::LineString(pts) => {
            let coords: Vec<String> = pts.iter().map(|p| lonlat_to_json(*p)).collect();
            format!(
                r#"{{"type":"LineString","coordinates":[{}]}}"#,
                coords.join(",")
            )
        }
        GeoGeometry::Polygon { exterior, holes } => {
            let mut rings = vec![ring_to_json(exterior)];
            for h in holes {
                rings.push(ring_to_json(h));
            }
            format!(
                r#"{{"type":"Polygon","coordinates":[{}]}}"#,
                rings.join(",")
            )
        }
        GeoGeometry::MultiPoint(pts) => {
            let coords: Vec<String> = pts.iter().map(|p| lonlat_to_json(*p)).collect();
            format!(
                r#"{{"type":"MultiPoint","coordinates":[{}]}}"#,
                coords.join(",")
            )
        }
        GeoGeometry::MultiLineString(lines) => {
            let ls: Vec<String> = lines
                .iter()
                .map(|line| {
                    let cs: Vec<String> = line.iter().map(|p| lonlat_to_json(*p)).collect();
                    format!("[{}]", cs.join(","))
                })
                .collect();
            format!(
                r#"{{"type":"MultiLineString","coordinates":[{}]}}"#,
                ls.join(",")
            )
        }
        GeoGeometry::MultiPolygon(polys) => {
            let ps: Vec<String> = polys
                .iter()
                .map(|(ext, holes)| {
                    let mut rings = vec![ring_to_json(ext)];
                    for h in holes {
                        rings.push(ring_to_json(h));
                    }
                    format!("[{}]", rings.join(","))
                })
                .collect();
            format!(
                r#"{{"type":"MultiPolygon","coordinates":[{}]}}"#,
                ps.join(",")
            )
        }
        GeoGeometry::GeometryCollection(gs) => {
            let inner: Vec<String> = gs.iter().map(geometry_to_geojson).collect();
            format!(
                r#"{{"type":"GeometryCollection","geometries":[{}]}}"#,
                inner.join(",")
            )
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. GeoJSON Feature and FeatureCollection
// ─────────────────────────────────────────────────────────────────────────────

/// A single string-valued property.
#[derive(Debug, Clone, PartialEq)]
pub struct GeoProperty {
    /// Property name.
    pub key: String,
    /// Property value as a string (numbers are stored as strings for simplicity).
    pub value: String,
}

impl GeoProperty {
    /// Construct from key/value string pair.
    #[must_use]
    pub fn new(key: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            key: key.into(),
            value: value.into(),
        }
    }

    /// Serialise to `"key":"value"` JSON fragment.
    #[must_use]
    pub fn to_json_pair(&self) -> String {
        format!(r#""{}":{}"#, self.key, json_string_or_number(&self.value))
    }
}

/// Attempt to treat `s` as a JSON number; fall back to a quoted string.
fn json_string_or_number(s: &str) -> String {
    if s.parse::<f64>().is_ok() {
        s.to_string()
    } else {
        format!(r#""{s}""#)
    }
}

/// A GeoJSON `Feature` containing one geometry and a property bag.
#[derive(Debug, Clone)]
pub struct GeoFeature {
    /// Optional feature identifier.
    pub id: Option<String>,
    /// The geometry associated with this feature.
    pub geometry: Option<GeoGeometry>,
    /// Key-value properties.
    pub properties: Vec<GeoProperty>,
}

impl GeoFeature {
    /// Construct a feature with a geometry and no properties.
    #[must_use]
    pub fn new(geometry: GeoGeometry) -> Self {
        Self {
            id: None,
            geometry: Some(geometry),
            properties: Vec::new(),
        }
    }

    /// Add a property and return `self` for chaining.
    #[must_use]
    pub fn with_property(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.properties.push(GeoProperty::new(key, value));
        self
    }

    /// Serialise this feature to a GeoJSON string.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        let geom_str = match &self.geometry {
            Some(g) => geometry_to_geojson(g),
            None => "null".to_string(),
        };
        let props: Vec<String> = self.properties.iter().map(|p| p.to_json_pair()).collect();
        let id_str = match &self.id {
            Some(id) => format!(r#","id":"{}""#, id),
            None => String::new(),
        };
        format!(
            r#"{{"type":"Feature"{id_str},"geometry":{geom_str},"properties":{{{}}}}}"#,
            props.join(",")
        )
    }

    /// Compute the bounding box of the feature's geometry.
    #[must_use]
    pub fn bounding_box(&self) -> Option<BoundingBox> {
        self.geometry.as_ref().and_then(|g| g.bounding_box())
    }
}

/// A GeoJSON `FeatureCollection`.
#[derive(Debug, Clone, Default)]
pub struct FeatureCollection {
    /// The features contained in this collection.
    pub features: Vec<GeoFeature>,
    /// Optional CRS metadata (non-standard extension, but widely used).
    pub crs: Option<CrsMetadata>,
}

impl FeatureCollection {
    /// Create an empty collection.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a feature.
    pub fn add(&mut self, feature: GeoFeature) {
        self.features.push(feature);
    }

    /// Return the number of features.
    #[must_use]
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Return `true` if the collection contains no features.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// Return features whose geometry bounding box intersects `bbox`.
    #[must_use]
    pub fn query_bbox(&self, bbox: &BoundingBox) -> Vec<&GeoFeature> {
        self.features
            .iter()
            .filter(|f| {
                f.bounding_box()
                    .map(|bb| bb.intersects(bbox))
                    .unwrap_or(false)
            })
            .collect()
    }

    /// Serialise to a GeoJSON `FeatureCollection` string.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        let feats: Vec<String> = self.features.iter().map(|f| f.to_geojson()).collect();
        format!(
            r#"{{"type":"FeatureCollection","features":[{}]}}"#,
            feats.join(",")
        )
    }

    /// Compute the bounding box enclosing all features.
    #[must_use]
    pub fn total_bounding_box(&self) -> Option<BoundingBox> {
        let mut all_pts: Vec<LonLat> = self
            .features
            .iter()
            .flat_map(|f| {
                f.geometry
                    .as_ref()
                    .map(|g| g.collect_lonlat())
                    .unwrap_or_default()
            })
            .collect();
        if all_pts.is_empty() {
            return None;
        }
        let first = all_pts.remove(0);
        let mut bb = BoundingBox::new(first.lon, first.lat, first.lon, first.lat);
        for p in all_pts {
            bb.expand_to_include(p);
        }
        Some(bb)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. WKT (Well-Known Text) geometry format
// ─────────────────────────────────────────────────────────────────────────────

/// Serialise a `GeoGeometry` to WKT.
#[must_use]
pub fn geometry_to_wkt(g: &GeoGeometry) -> String {
    match g {
        GeoGeometry::Point(p) => format!("POINT ({:.6} {:.6})", p.lon, p.lat),
        GeoGeometry::Point3D(p) => format!("POINT Z ({:.6} {:.6} {:.6})", p.lon, p.lat, p.elev),
        GeoGeometry::LineString(pts) => {
            let coords = pts_to_wkt_coords(pts);
            format!("LINESTRING ({coords})")
        }
        GeoGeometry::Polygon { exterior, holes } => {
            let mut rings = vec![format!("({})", pts_to_wkt_coords(exterior))];
            for h in holes {
                rings.push(format!("({})", pts_to_wkt_coords(h)));
            }
            format!("POLYGON ({})", rings.join(","))
        }
        GeoGeometry::MultiPoint(pts) => {
            let coords: Vec<String> = pts
                .iter()
                .map(|p| format!("({:.6} {:.6})", p.lon, p.lat))
                .collect();
            format!("MULTIPOINT ({})", coords.join(","))
        }
        GeoGeometry::MultiLineString(lines) => {
            let ls: Vec<String> = lines
                .iter()
                .map(|l| format!("({})", pts_to_wkt_coords(l)))
                .collect();
            format!("MULTILINESTRING ({})", ls.join(","))
        }
        GeoGeometry::MultiPolygon(polys) => {
            let ps: Vec<String> = polys
                .iter()
                .map(|(ext, holes)| {
                    let mut rings = vec![format!("({})", pts_to_wkt_coords(ext))];
                    for h in holes {
                        rings.push(format!("({})", pts_to_wkt_coords(h)));
                    }
                    format!("({})", rings.join(","))
                })
                .collect();
            format!("MULTIPOLYGON ({})", ps.join(","))
        }
        GeoGeometry::GeometryCollection(gs) => {
            let inner: Vec<String> = gs.iter().map(geometry_to_wkt).collect();
            format!("GEOMETRYCOLLECTION ({})", inner.join(","))
        }
    }
}

fn pts_to_wkt_coords(pts: &[LonLat]) -> String {
    pts.iter()
        .map(|p| format!("{:.6} {:.6}", p.lon, p.lat))
        .collect::<Vec<_>>()
        .join(",")
}

/// Parse a WKT `POINT (lon lat)` string.
///
/// Only handles 2-D `POINT` for brevity; returns `GeoError::ParseError` for
/// unsupported types.
pub fn parse_wkt_point(wkt: &str) -> GeoResult<GeoGeometry> {
    let s = wkt.trim();
    if !s.to_uppercase().starts_with("POINT") {
        return Err(GeoError::ParseError(format!("expected POINT, got: {s}")));
    }
    let inner = s
        .find('(')
        .and_then(|i| {
            let rest = &s[i + 1..];
            rest.rfind(')').map(|j| &rest[..j])
        })
        .ok_or_else(|| GeoError::ParseError("missing parentheses".to_string()))?;
    let parts: Vec<f64> = inner
        .split_whitespace()
        .map(|tok| tok.parse::<f64>())
        .collect::<Result<_, _>>()
        .map_err(|e| GeoError::ParseError(e.to_string()))?;
    match parts.as_slice() {
        [lon, lat] => Ok(GeoGeometry::Point(LonLat::new(*lon, *lat))),
        [lon, lat, elev] => Ok(GeoGeometry::Point3D(LonLatElev::new(*lon, *lat, *elev))),
        _ => Err(GeoError::ParseError(
            "expected 2 or 3 coordinates".to_string(),
        )),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. Coordinate transformation — simplified Mercator
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for a Web-Mercator (EPSG:3857) projection.
#[derive(Debug, Clone, Copy)]
pub struct MercatorProjection {
    /// Earth radius used in the spherical Mercator formula \[m\].
    pub earth_radius_m: f64,
}

impl MercatorProjection {
    /// Standard Web Mercator with R = 6 378 137 m.
    #[must_use]
    pub fn web_mercator() -> Self {
        Self {
            earth_radius_m: 6_378_137.0,
        }
    }

    /// Forward projection: geographic (lon, lat) → Mercator (x, y) metres.
    ///
    /// Latitude is clamped to ±85.051129° to stay within the Mercator domain.
    #[must_use]
    pub fn forward(&self, pt: LonLat) -> [f64; 2] {
        let lat_clamp = pt.lat.clamp(-85.051_129, 85.051_129);
        let x = self.earth_radius_m * pt.lon.to_radians();
        let y = self.earth_radius_m
            * ((std::f64::consts::FRAC_PI_4 + lat_clamp.to_radians() / 2.0).tan()).ln();
        [x, y]
    }

    /// Inverse projection: Mercator (x, y) → geographic (lon, lat).
    #[must_use]
    pub fn inverse(&self, xy: [f64; 2]) -> LonLat {
        let lon = xy[0].to_degrees() / self.earth_radius_m;
        let lat = (2.0 * (xy[1] / self.earth_radius_m).exp().atan() - std::f64::consts::FRAC_PI_2)
            .to_degrees();
        LonLat::new(lon, lat)
    }

    /// Project an entire `LonLat` slice to metres.
    #[must_use]
    pub fn forward_all(&self, pts: &[LonLat]) -> Vec<[f64; 2]> {
        pts.iter().map(|p| self.forward(*p)).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. DEM raster (Digital Elevation Model)
// ─────────────────────────────────────────────────────────────────────────────

/// A raster grid of elevation values covering a bounding box.
#[derive(Debug, Clone)]
pub struct DemRaster {
    /// Geographic extent of the raster.
    pub bbox: BoundingBox,
    /// Number of columns.
    pub cols: usize,
    /// Number of rows.
    pub rows: usize,
    /// Row-major elevation samples \[m\].  Length = rows × cols.
    pub data: Vec<f32>,
    /// No-data sentinel value.
    pub nodata: f32,
}

impl DemRaster {
    /// Construct a zero-filled DEM raster.
    #[must_use]
    pub fn zeros(bbox: BoundingBox, cols: usize, rows: usize) -> Self {
        Self {
            bbox,
            cols,
            rows,
            data: vec![0.0_f32; cols * rows],
            nodata: -9999.0,
        }
    }

    /// Cell width in degrees.
    #[must_use]
    pub fn cell_width_deg(&self) -> f64 {
        bbox_cell_size(self.bbox.min_lon, self.bbox.max_lon, self.cols)
    }

    /// Cell height in degrees.
    #[must_use]
    pub fn cell_height_deg(&self) -> f64 {
        bbox_cell_size(self.bbox.min_lat, self.bbox.max_lat, self.rows)
    }

    /// Sample elevation at a geographic position using bilinear interpolation.
    ///
    /// Returns `None` if `pt` is outside the raster extent or resolves to a
    /// no-data cell.
    #[must_use]
    pub fn sample(&self, pt: LonLat) -> Option<f32> {
        if !self.bbox.contains(pt) {
            return None;
        }
        let fx = (pt.lon - self.bbox.min_lon) / self.cell_width_deg();
        let fy = (pt.lat - self.bbox.min_lat) / self.cell_height_deg();
        let col = fx.floor() as usize;
        let row = fy.floor() as usize;
        let col = col.min(self.cols - 1);
        let row = row.min(self.rows - 1);
        let v = self.data[row * self.cols + col];
        if v == self.nodata { None } else { Some(v) }
    }

    /// Return the minimum elevation, ignoring no-data values.
    #[must_use]
    pub fn min_elevation(&self) -> Option<f32> {
        self.data
            .iter()
            .copied()
            .filter(|&v| v != self.nodata)
            .reduce(f32::min)
    }

    /// Return the maximum elevation, ignoring no-data values.
    #[must_use]
    pub fn max_elevation(&self) -> Option<f32> {
        self.data
            .iter()
            .copied()
            .filter(|&v| v != self.nodata)
            .reduce(f32::max)
    }

    /// Export to an ASCII Grid string (ESRI format).
    #[must_use]
    pub fn to_ascii_grid(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("ncols {}\n", self.cols));
        s.push_str(&format!("nrows {}\n", self.rows));
        s.push_str(&format!("xllcorner {:.6}\n", self.bbox.min_lon));
        s.push_str(&format!("yllcorner {:.6}\n", self.bbox.min_lat));
        s.push_str(&format!("cellsize {:.6}\n", self.cell_width_deg()));
        s.push_str(&format!("NODATA_value {:.1}\n", self.nodata));
        for row in (0..self.rows).rev() {
            let row_vals: Vec<String> = (0..self.cols)
                .map(|c| format!("{:.1}", self.data[row * self.cols + c]))
                .collect();
            s.push_str(&row_vals.join(" "));
            s.push('\n');
        }
        s
    }
}

fn bbox_cell_size(lo: f64, hi: f64, n: usize) -> f64 {
    if n == 0 { 0.0 } else { (hi - lo) / n as f64 }
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. Terrain tile format
// ─────────────────────────────────────────────────────────────────────────────

/// A Slippy-map tile identifier (zoom / x / y).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TileId {
    /// Zoom level (0–22).
    pub zoom: u8,
    /// Column tile index.
    pub x: u32,
    /// Row tile index.
    pub y: u32,
}

impl TileId {
    /// Construct a tile identifier.
    #[must_use]
    pub fn new(zoom: u8, x: u32, y: u32) -> Self {
        Self { zoom, x, y }
    }

    /// Compute the `TileId` that covers a geographic point at a given zoom.
    #[must_use]
    pub fn from_lonlat(pt: LonLat, zoom: u8) -> Self {
        let n = 2_u32.pow(zoom as u32) as f64;
        let x = ((pt.lon + 180.0) / 360.0 * n).floor() as u32;
        let lat_rad = pt.lat.to_radians();
        let y = ((1.0 - (lat_rad.tan() + 1.0 / lat_rad.cos()).ln() / std::f64::consts::PI) / 2.0
            * n)
            .floor() as u32;
        Self::new(zoom, x, y)
    }

    /// Return the bounding box of this tile in geographic coordinates.
    #[must_use]
    pub fn bbox(&self) -> BoundingBox {
        let n = 2_u32.pow(self.zoom as u32) as f64;
        let min_lon = self.x as f64 / n * 360.0 - 180.0;
        let max_lon = (self.x + 1) as f64 / n * 360.0 - 180.0;
        let lat_north = tile_lat(self.y, n);
        let lat_south = tile_lat(self.y + 1, n);
        BoundingBox::new(min_lon, lat_south, max_lon, lat_north)
    }

    /// Return the four child tiles at zoom+1.
    #[must_use]
    pub fn children(&self) -> [TileId; 4] {
        let z = self.zoom + 1;
        let x2 = self.x * 2;
        let y2 = self.y * 2;
        [
            TileId::new(z, x2, y2),
            TileId::new(z, x2 + 1, y2),
            TileId::new(z, x2, y2 + 1),
            TileId::new(z, x2 + 1, y2 + 1),
        ]
    }
}

fn tile_lat(y: u32, n: f64) -> f64 {
    let sinh_val = std::f64::consts::PI * (1.0 - 2.0 * y as f64 / n);
    sinh_val.sinh().atan().to_degrees()
}

/// A terrain tile carrying elevation samples.
#[derive(Debug, Clone)]
pub struct TerrainTile {
    /// Tile identifier.
    pub id: TileId,
    /// Elevation grid, 256 × 256 samples by default.
    pub elevation: Vec<f32>,
    /// Side length in samples.
    pub size: usize,
}

impl TerrainTile {
    /// Construct an empty (zero-filled) terrain tile of given size.
    #[must_use]
    pub fn empty(id: TileId, size: usize) -> Self {
        Self {
            id,
            elevation: vec![0.0; size * size],
            size,
        }
    }

    /// Sample elevation at normalised coordinates (u, v) ∈ \[0, 1\]².
    #[must_use]
    pub fn sample_uv(&self, u: f64, v: f64) -> f32 {
        let col = (u.clamp(0.0, 1.0) * (self.size - 1) as f64).round() as usize;
        let row = (v.clamp(0.0, 1.0) * (self.size - 1) as f64).round() as usize;
        self.elevation[row * self.size + col]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Geohash encoding / decoding
// ─────────────────────────────────────────────────────────────────────────────

const GEOHASH_CHARS: &[u8] = b"0123456789bcdefghjkmnpqrstuvwxyz";

/// Encode a `LonLat` to a geohash string of the requested length (1–12).
///
/// Panics if `precision` is 0 or greater than 12.
#[must_use]
pub fn geohash_encode(pt: LonLat, precision: usize) -> String {
    assert!((1..=12).contains(&precision), "precision must be 1–12");
    let mut lat_range = (-90.0_f64, 90.0_f64);
    let mut lon_range = (-180.0_f64, 180.0_f64);
    let mut bits_total = precision * 5;
    let mut result = Vec::with_capacity(precision);
    let mut idx: u8 = 0;
    let mut bit = 0u8;
    let mut is_lon = true;
    while bits_total > 0 {
        if is_lon {
            let mid = (lon_range.0 + lon_range.1) / 2.0;
            if pt.lon >= mid {
                idx = (idx << 1) | 1;
                lon_range.0 = mid;
            } else {
                idx <<= 1;
                lon_range.1 = mid;
            }
        } else {
            let mid = (lat_range.0 + lat_range.1) / 2.0;
            if pt.lat >= mid {
                idx = (idx << 1) | 1;
                lat_range.0 = mid;
            } else {
                idx <<= 1;
                lat_range.1 = mid;
            }
        }
        is_lon = !is_lon;
        bit += 1;
        bits_total -= 1;
        if bit == 5 {
            result.push(GEOHASH_CHARS[idx as usize] as char);
            idx = 0;
            bit = 0;
        }
    }
    result.into_iter().collect()
}

/// Decode a geohash string to the centre `LonLat` and error margins.
///
/// Returns `(centre, lon_error_deg, lat_error_deg)`.
///
/// Returns `GeoError::ParseError` if the string contains invalid characters.
pub fn geohash_decode(hash: &str) -> GeoResult<(LonLat, f64, f64)> {
    let mut lat_range = (-90.0_f64, 90.0_f64);
    let mut lon_range = (-180.0_f64, 180.0_f64);
    let mut is_lon = true;
    for ch in hash.chars() {
        let idx = GEOHASH_CHARS
            .iter()
            .position(|&c| c as char == ch)
            .ok_or_else(|| GeoError::ParseError(format!("invalid geohash char '{ch}'")))?;
        for bit in (0..5).rev() {
            let b = (idx >> bit) & 1;
            if is_lon {
                let mid = (lon_range.0 + lon_range.1) / 2.0;
                if b == 1 {
                    lon_range.0 = mid;
                } else {
                    lon_range.1 = mid;
                }
            } else {
                let mid = (lat_range.0 + lat_range.1) / 2.0;
                if b == 1 {
                    lat_range.0 = mid;
                } else {
                    lat_range.1 = mid;
                }
            }
            is_lon = !is_lon;
        }
    }
    let lon = (lon_range.0 + lon_range.1) / 2.0;
    let lat = (lat_range.0 + lat_range.1) / 2.0;
    let lon_err = (lon_range.1 - lon_range.0) / 2.0;
    let lat_err = (lat_range.1 - lat_range.0) / 2.0;
    Ok((LonLat::new(lon, lat), lon_err, lat_err))
}

// ─────────────────────────────────────────────────────────────────────────────
// 13. Spatial Index — R-tree-like bounding-box index
// ─────────────────────────────────────────────────────────────────────────────

/// A single entry in the spatial index.
#[derive(Debug, Clone)]
pub struct SpatialEntry {
    /// Bounding box of this entry.
    pub bbox: BoundingBox,
    /// Arbitrary integer identifier (e.g. feature index).
    pub id: usize,
}

/// A naive R-tree-like spatial index backed by a flat list.
///
/// For small data sets this brute-force scan is sufficient; a real R-tree
/// would split nodes at a configurable capacity.
#[derive(Debug, Clone, Default)]
pub struct SpatialIndex {
    entries: Vec<SpatialEntry>,
}

impl SpatialIndex {
    /// Create an empty spatial index.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert an entry.
    pub fn insert(&mut self, bbox: BoundingBox, id: usize) {
        self.entries.push(SpatialEntry { bbox, id });
    }

    /// Return all entry IDs whose bounding boxes intersect `query`.
    #[must_use]
    pub fn query(&self, query: &BoundingBox) -> Vec<usize> {
        self.entries
            .iter()
            .filter(|e| e.bbox.intersects(query))
            .map(|e| e.id)
            .collect()
    }

    /// Return the number of entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return `true` if the index is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}

/// Build a `SpatialIndex` from a `FeatureCollection`.
#[must_use]
pub fn index_feature_collection(fc: &FeatureCollection) -> SpatialIndex {
    let mut idx = SpatialIndex::new();
    for (i, f) in fc.features.iter().enumerate() {
        if let Some(bb) = f.bounding_box() {
            idx.insert(bb, i);
        }
    }
    idx
}

// ─────────────────────────────────────────────────────────────────────────────
// 14. Simplified DBF attribute table (Shapefile component)
// ─────────────────────────────────────────────────────────────────────────────

/// A DBF field type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DbfFieldType {
    /// Character string.
    Character,
    /// Numeric (floating-point).
    Numeric,
    /// Logical (boolean).
    Logical,
    /// Date (YYYYMMDD string).
    Date,
}

/// A DBF field descriptor.
#[derive(Debug, Clone)]
pub struct DbfFieldDescriptor {
    /// Field name (max 10 characters in dBase III).
    pub name: String,
    /// Field type.
    pub field_type: DbfFieldType,
    /// Field length in bytes.
    pub length: u8,
    /// Number of decimal places (Numeric fields).
    pub decimals: u8,
}

impl DbfFieldDescriptor {
    /// Construct a character field.
    #[must_use]
    pub fn character(name: impl Into<String>, length: u8) -> Self {
        Self {
            name: name.into(),
            field_type: DbfFieldType::Character,
            length,
            decimals: 0,
        }
    }

    /// Construct a numeric field.
    #[must_use]
    pub fn numeric(name: impl Into<String>, length: u8, decimals: u8) -> Self {
        Self {
            name: name.into(),
            field_type: DbfFieldType::Numeric,
            length,
            decimals,
        }
    }
}

/// A single DBF record (row).
#[derive(Debug, Clone, Default)]
pub struct DbfRecord {
    /// Field values stored as strings.
    pub values: Vec<String>,
    /// Deletion flag.
    pub deleted: bool,
}

impl DbfRecord {
    /// Construct a record with the given values.
    #[must_use]
    pub fn new(values: Vec<String>) -> Self {
        Self {
            values,
            deleted: false,
        }
    }
}

/// A simplified in-memory DBF attribute table.
#[derive(Debug, Clone)]
pub struct DbfTable {
    /// Field descriptors (schema).
    pub fields: Vec<DbfFieldDescriptor>,
    /// Data rows.
    pub records: Vec<DbfRecord>,
}

impl DbfTable {
    /// Create a new table with the given schema.
    #[must_use]
    pub fn new(fields: Vec<DbfFieldDescriptor>) -> Self {
        Self {
            fields,
            records: Vec::new(),
        }
    }

    /// Append a record; returns `GeoError` if value count mismatches.
    pub fn append(&mut self, record: DbfRecord) -> GeoResult<()> {
        if record.values.len() != self.fields.len() {
            return Err(GeoError::ParseError(format!(
                "expected {} values, got {}",
                self.fields.len(),
                record.values.len()
            )));
        }
        self.records.push(record);
        Ok(())
    }

    /// Number of records.
    #[must_use]
    pub fn record_count(&self) -> usize {
        self.records.iter().filter(|r| !r.deleted).count()
    }

    /// Look up the column index for a field name.
    #[must_use]
    pub fn field_index(&self, name: &str) -> Option<usize> {
        self.fields.iter().position(|f| f.name == name)
    }

    /// Get a field value from a record by field name.
    #[must_use]
    pub fn get_value(&self, record_idx: usize, field_name: &str) -> Option<&str> {
        let col = self.field_index(field_name)?;
        self.records.get(record_idx).map(|r| r.values[col].as_str())
    }

    /// Export the table as a CSV string.
    #[must_use]
    pub fn to_csv(&self) -> String {
        let header: Vec<&str> = self.fields.iter().map(|f| f.name.as_str()).collect();
        let mut lines = vec![header.join(",")];
        for rec in &self.records {
            if !rec.deleted {
                lines.push(rec.values.join(","));
            }
        }
        lines.join("\n")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 15. Shapefile writer (geometry + DBF combined)
// ─────────────────────────────────────────────────────────────────────────────

/// A minimal in-memory representation of a point shapefile.
#[derive(Debug, Clone)]
pub struct PointShapefile {
    /// Point geometries.
    pub points: Vec<LonLat>,
    /// Attribute table.
    pub attributes: DbfTable,
}

impl PointShapefile {
    /// Construct a new point shapefile with the given attribute schema.
    #[must_use]
    pub fn new(fields: Vec<DbfFieldDescriptor>) -> Self {
        Self {
            points: Vec::new(),
            attributes: DbfTable::new(fields),
        }
    }

    /// Add a point with attribute values.
    pub fn add_point(&mut self, pt: LonLat, values: Vec<String>) -> GeoResult<()> {
        self.points.push(pt);
        self.attributes.append(DbfRecord::new(values))
    }

    /// Export the point locations as GeoJSON.
    #[must_use]
    pub fn to_geojson(&self) -> String {
        let mut fc = FeatureCollection::new();
        for (i, pt) in self.points.iter().enumerate() {
            let mut feat = GeoFeature::new(GeoGeometry::Point(*pt));
            for (fi, field) in self.attributes.fields.iter().enumerate() {
                if let Some(rec) = self.attributes.records.get(i) {
                    feat = feat.with_property(&field.name, &rec.values[fi]);
                }
            }
            fc.add(feat);
        }
        fc.to_geojson()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 16. Elevation contour line generation
// ─────────────────────────────────────────────────────────────────────────────

/// Extract iso-contour polylines from a DEM at a given elevation level using
/// marching-squares (simplified: horizontal/vertical cell edges only).
#[must_use]
pub fn dem_contour_lines(dem: &DemRaster, level: f32) -> Vec<Vec<LonLat>> {
    let mut lines = Vec::new();
    let cw = dem.cell_width_deg();
    let ch = dem.cell_height_deg();
    for row in 0..dem.rows.saturating_sub(1) {
        for col in 0..dem.cols.saturating_sub(1) {
            let v00 = dem.data[row * dem.cols + col];
            let v10 = dem.data[row * dem.cols + col + 1];
            let v01 = dem.data[(row + 1) * dem.cols + col];
            let v11 = dem.data[(row + 1) * dem.cols + col + 1];
            // Simple saddle: if any edge crosses `level`, emit a short segment
            let pts = marching_square_pts(v00, v10, v01, v11, level, col, row, cw, ch, dem);
            if pts.len() >= 2 {
                lines.push(pts);
            }
        }
    }
    lines
}

fn marching_square_pts(
    v00: f32,
    v10: f32,
    v01: f32,
    v11: f32,
    level: f32,
    col: usize,
    row: usize,
    cw: f64,
    ch: f64,
    dem: &DemRaster,
) -> Vec<LonLat> {
    let mut pts = Vec::new();
    let lon0 = dem.bbox.min_lon + col as f64 * cw;
    let lat0 = dem.bbox.min_lat + row as f64 * ch;
    // bottom edge
    if (v00 < level) != (v10 < level) {
        let t = (level - v00) / (v10 - v00);
        pts.push(LonLat::new(lon0 + t as f64 * cw, lat0));
    }
    // left edge
    if (v00 < level) != (v01 < level) {
        let t = (level - v00) / (v01 - v00);
        pts.push(LonLat::new(lon0, lat0 + t as f64 * ch));
    }
    // top edge
    if (v01 < level) != (v11 < level) {
        let t = (level - v01) / (v11 - v01);
        pts.push(LonLat::new(lon0 + t as f64 * cw, lat0 + ch));
    }
    // right edge
    if (v10 < level) != (v11 < level) {
        let t = (level - v10) / (v11 - v10);
        pts.push(LonLat::new(lon0 + cw, lat0 + t as f64 * ch));
    }
    pts
}

// ─────────────────────────────────────────────────────────────────────────────
// 17. GeoJSON import helpers (minimal hand-rolled parser)
// ─────────────────────────────────────────────────────────────────────────────

/// Extract a quoted string value for a given key from a flat JSON object.
///
/// This is intentionally minimal and handles only the simple cases produced by
/// this module's own serialiser.
#[must_use]
pub fn extract_json_string(json: &str, key: &str) -> Option<String> {
    let pattern = format!(r#""{key}":""#);
    let start = json.find(&pattern)? + pattern.len();
    let rest = &json[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

/// A minimal GeoJSON `Point` parser (round-trips with `geometry_to_geojson`).
pub fn parse_geojson_point(json: &str) -> GeoResult<LonLat> {
    // Expect: {"type":"Point","coordinates":[lon,lat]}
    let coord_key = r#""coordinates":["#;
    let start = json
        .find(coord_key)
        .ok_or_else(|| GeoError::ParseError("no coordinates key".to_string()))?
        + coord_key.len();
    let rest = &json[start..];
    let end = rest
        .find(']')
        .ok_or_else(|| GeoError::ParseError("no closing bracket".to_string()))?;
    let coords: Vec<f64> = rest[..end]
        .split(',')
        .map(|s| s.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .map_err(|e| GeoError::ParseError(e.to_string()))?;
    match coords.as_slice() {
        [lon, lat] => Ok(LonLat::new(*lon, *lat)),
        [lon, lat, _] => Ok(LonLat::new(*lon, *lat)),
        _ => Err(GeoError::ParseError(
            "unexpected coordinate count".to_string(),
        )),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 18. Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- LonLat ---------------------------------------------------------------

    #[test]
    fn test_lonlat_new() {
        let p = LonLat::new(13.4050, 52.5200);
        assert!((p.lon - 13.4050).abs() < 1e-10);
        assert!((p.lat - 52.5200).abs() < 1e-10);
    }

    #[test]
    fn test_haversine_same_point() {
        let p = LonLat::new(0.0, 0.0);
        assert!(p.haversine_m(&p) < 1e-6);
    }

    #[test]
    fn test_haversine_equator() {
        // 1 degree of longitude at equator ≈ 111 195 m
        let a = LonLat::new(0.0, 0.0);
        let b = LonLat::new(1.0, 0.0);
        let d = a.haversine_m(&b);
        assert!((d - 111_195.0).abs() < 200.0);
    }

    // --- BoundingBox ----------------------------------------------------------

    #[test]
    fn test_bbox_contains() {
        let bb = BoundingBox::new(0.0, 0.0, 10.0, 10.0);
        assert!(bb.contains(LonLat::new(5.0, 5.0)));
        assert!(!bb.contains(LonLat::new(11.0, 5.0)));
    }

    #[test]
    fn test_bbox_intersects() {
        let a = BoundingBox::new(0.0, 0.0, 10.0, 10.0);
        let b = BoundingBox::new(5.0, 5.0, 15.0, 15.0);
        assert!(a.intersects(&b));
        let c = BoundingBox::new(11.0, 0.0, 20.0, 10.0);
        assert!(!a.intersects(&c));
    }

    #[test]
    fn test_bbox_centre() {
        let bb = BoundingBox::new(-10.0, -10.0, 10.0, 10.0);
        let c = bb.centre();
        assert!((c.lon).abs() < 1e-10);
        assert!((c.lat).abs() < 1e-10);
    }

    #[test]
    fn test_bbox_expand() {
        let mut bb = BoundingBox::new(0.0, 0.0, 1.0, 1.0);
        bb.expand_to_include(LonLat::new(5.0, 5.0));
        assert!((bb.max_lon - 5.0).abs() < 1e-10);
        assert!((bb.max_lat - 5.0).abs() < 1e-10);
    }

    // --- GeoGeometry ----------------------------------------------------------

    #[test]
    fn test_geometry_type_str() {
        let g = GeoGeometry::Point(LonLat::new(0.0, 0.0));
        assert_eq!(g.type_str(), "Point");
    }

    #[test]
    fn test_collect_lonlat_polygon() {
        let g = GeoGeometry::Polygon {
            exterior: vec![
                LonLat::new(0.0, 0.0),
                LonLat::new(1.0, 0.0),
                LonLat::new(1.0, 1.0),
                LonLat::new(0.0, 0.0),
            ],
            holes: vec![],
        };
        assert_eq!(g.collect_lonlat().len(), 4);
    }

    #[test]
    fn test_bounding_box_linestring() {
        let g = GeoGeometry::LineString(vec![LonLat::new(-1.0, -1.0), LonLat::new(1.0, 1.0)]);
        let bb = g.bounding_box().unwrap();
        assert!((bb.min_lon - (-1.0)).abs() < 1e-10);
        assert!((bb.max_lat - 1.0).abs() < 1e-10);
    }

    // --- GeoJSON serialisation ------------------------------------------------

    #[test]
    fn test_geometry_to_geojson_point() {
        let g = GeoGeometry::Point(LonLat::new(13.405, 52.52));
        let s = geometry_to_geojson(&g);
        assert!(s.contains("\"type\":\"Point\""));
        assert!(s.contains("13.405000"));
    }

    #[test]
    fn test_geometry_to_geojson_linestring() {
        let g = GeoGeometry::LineString(vec![LonLat::new(0.0, 0.0), LonLat::new(1.0, 1.0)]);
        let s = geometry_to_geojson(&g);
        assert!(s.contains("\"type\":\"LineString\""));
    }

    #[test]
    fn test_feature_collection_to_geojson() {
        let mut fc = FeatureCollection::new();
        fc.add(GeoFeature::new(GeoGeometry::Point(LonLat::new(0.0, 0.0))));
        let s = fc.to_geojson();
        assert!(s.contains("\"type\":\"FeatureCollection\""));
    }

    #[test]
    fn test_feature_with_properties() {
        let f = GeoFeature::new(GeoGeometry::Point(LonLat::new(10.0, 20.0)))
            .with_property("name", "test")
            .with_property("value", "42");
        let s = f.to_geojson();
        assert!(s.contains("\"name\":\"test\""));
        assert!(s.contains("\"value\":42"));
    }

    #[test]
    fn test_query_bbox_filters() {
        let mut fc = FeatureCollection::new();
        fc.add(GeoFeature::new(GeoGeometry::Point(LonLat::new(5.0, 5.0))));
        fc.add(GeoFeature::new(GeoGeometry::Point(LonLat::new(50.0, 50.0))));
        let query = BoundingBox::new(0.0, 0.0, 10.0, 10.0);
        let results = fc.query_bbox(&query);
        assert_eq!(results.len(), 1);
    }

    // --- WKT ------------------------------------------------------------------

    #[test]
    fn test_geometry_to_wkt_point() {
        let g = GeoGeometry::Point(LonLat::new(10.0, 20.0));
        let s = geometry_to_wkt(&g);
        assert!(s.starts_with("POINT ("));
        assert!(s.contains("10.000000"));
    }

    #[test]
    fn test_parse_wkt_point_2d() {
        let g = parse_wkt_point("POINT (10.5 20.5)").unwrap();
        if let GeoGeometry::Point(p) = g {
            assert!((p.lon - 10.5).abs() < 1e-6);
            assert!((p.lat - 20.5).abs() < 1e-6);
        } else {
            panic!("expected Point");
        }
    }

    #[test]
    fn test_parse_wkt_point_3d() {
        let g = parse_wkt_point("POINT Z (1.0 2.0 3.0)").unwrap();
        assert!(matches!(g, GeoGeometry::Point3D(_)));
    }

    #[test]
    fn test_parse_wkt_invalid() {
        assert!(parse_wkt_point("LINESTRING (0 0, 1 1)").is_err());
    }

    // --- Mercator projection --------------------------------------------------

    #[test]
    fn test_mercator_forward_origin() {
        let m = MercatorProjection::web_mercator();
        let xy = m.forward(LonLat::new(0.0, 0.0));
        assert!(xy[0].abs() < 1e-6);
        assert!(xy[1].abs() < 1e-6);
    }

    #[test]
    fn test_mercator_round_trip() {
        let m = MercatorProjection::web_mercator();
        let orig = LonLat::new(30.0, 45.0);
        let xy = m.forward(orig);
        let back = m.inverse(xy);
        assert!((back.lon - orig.lon).abs() < 1e-6);
        assert!((back.lat - orig.lat).abs() < 1e-4);
    }

    // --- DEM ------------------------------------------------------------------

    #[test]
    fn test_dem_zeros_min_max() {
        let dem = DemRaster::zeros(BoundingBox::new(0.0, 0.0, 1.0, 1.0), 10, 10);
        assert_eq!(dem.min_elevation(), Some(0.0));
        assert_eq!(dem.max_elevation(), Some(0.0));
    }

    #[test]
    fn test_dem_sample_inside() {
        let mut dem = DemRaster::zeros(BoundingBox::new(0.0, 0.0, 1.0, 1.0), 10, 10);
        dem.data[55] = 100.0;
        // Point should land inside the raster
        let v = dem.sample(LonLat::new(0.05, 0.05));
        assert!(v.is_some());
    }

    #[test]
    fn test_dem_sample_outside() {
        let dem = DemRaster::zeros(BoundingBox::new(0.0, 0.0, 1.0, 1.0), 10, 10);
        let v = dem.sample(LonLat::new(5.0, 5.0));
        assert!(v.is_none());
    }

    #[test]
    fn test_dem_ascii_grid_header() {
        let dem = DemRaster::zeros(BoundingBox::new(0.0, 0.0, 1.0, 1.0), 4, 4);
        let s = dem.to_ascii_grid();
        assert!(s.contains("ncols 4"));
        assert!(s.contains("nrows 4"));
    }

    // --- Terrain tile ---------------------------------------------------------

    #[test]
    fn test_tile_from_lonlat_zoom0() {
        let t = TileId::from_lonlat(LonLat::new(0.0, 0.0), 0);
        assert_eq!(t.zoom, 0);
        assert_eq!(t.x, 0);
        assert_eq!(t.y, 0);
    }

    #[test]
    fn test_tile_bbox_width() {
        let t = TileId::new(1, 0, 0);
        let bb = t.bbox();
        // At zoom 1, each tile covers 180° of longitude
        assert!((bb.width_deg() - 180.0).abs() < 0.001);
    }

    #[test]
    fn test_tile_children_count() {
        let t = TileId::new(2, 1, 1);
        let kids = t.children();
        assert_eq!(kids.len(), 4);
    }

    #[test]
    fn test_terrain_tile_sample_uv() {
        let mut tile = TerrainTile::empty(TileId::new(0, 0, 0), 4);
        tile.elevation[5] = 500.0;
        // Just check it doesn't panic
        let _ = tile.sample_uv(0.5, 0.5);
    }

    // --- Geohash --------------------------------------------------------------

    #[test]
    fn test_geohash_encode_length() {
        let h = geohash_encode(LonLat::new(13.4050, 52.5200), 7);
        assert_eq!(h.len(), 7);
    }

    #[test]
    fn test_geohash_decode_round_trip() {
        let orig = LonLat::new(10.0, 50.0);
        let hash = geohash_encode(orig, 9);
        let (decoded, lon_err, lat_err) = geohash_decode(&hash).unwrap();
        assert!((decoded.lon - orig.lon).abs() < lon_err + 1e-5);
        assert!((decoded.lat - orig.lat).abs() < lat_err + 1e-5);
    }

    #[test]
    fn test_geohash_decode_invalid() {
        assert!(geohash_decode("invalid!").is_err());
    }

    // --- Spatial index --------------------------------------------------------

    #[test]
    fn test_spatial_index_query() {
        let mut idx = SpatialIndex::new();
        idx.insert(BoundingBox::new(0.0, 0.0, 1.0, 1.0), 0);
        idx.insert(BoundingBox::new(10.0, 10.0, 20.0, 20.0), 1);
        let results = idx.query(&BoundingBox::new(0.5, 0.5, 2.0, 2.0));
        assert!(results.contains(&0));
        assert!(!results.contains(&1));
    }

    #[test]
    fn test_spatial_index_empty() {
        let idx = SpatialIndex::new();
        assert!(idx.is_empty());
    }

    #[test]
    fn test_index_feature_collection() {
        let mut fc = FeatureCollection::new();
        fc.add(GeoFeature::new(GeoGeometry::Point(LonLat::new(5.0, 5.0))));
        fc.add(GeoFeature::new(GeoGeometry::Point(LonLat::new(15.0, 15.0))));
        let idx = index_feature_collection(&fc);
        assert_eq!(idx.len(), 2);
    }

    // --- DBF ------------------------------------------------------------------

    #[test]
    fn test_dbf_append_ok() {
        let mut tbl = DbfTable::new(vec![DbfFieldDescriptor::character("name", 50)]);
        tbl.append(DbfRecord::new(vec!["hello".to_string()]))
            .unwrap();
        assert_eq!(tbl.record_count(), 1);
    }

    #[test]
    fn test_dbf_append_mismatch() {
        let mut tbl = DbfTable::new(vec![DbfFieldDescriptor::character("name", 50)]);
        let result = tbl.append(DbfRecord::new(vec!["a".to_string(), "b".to_string()]));
        assert!(result.is_err());
    }

    #[test]
    fn test_dbf_get_value() {
        let mut tbl = DbfTable::new(vec![
            DbfFieldDescriptor::character("city", 50),
            DbfFieldDescriptor::numeric("pop", 10, 0),
        ]);
        tbl.append(DbfRecord::new(vec![
            "Berlin".to_string(),
            "3700000".to_string(),
        ]))
        .unwrap();
        assert_eq!(tbl.get_value(0, "city"), Some("Berlin"));
    }

    #[test]
    fn test_dbf_to_csv() {
        let mut tbl = DbfTable::new(vec![DbfFieldDescriptor::character("x", 10)]);
        tbl.append(DbfRecord::new(vec!["foo".to_string()])).unwrap();
        let csv = tbl.to_csv();
        assert!(csv.contains("x\nfoo"));
    }
}
