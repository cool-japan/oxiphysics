//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{GevParameters, MetarReport, MetarWind};

/// Parse the wind group from a METAR token (e.g. `"25015G25KT"` or `"VRB05KT"`).
pub fn parse_metar_wind(token: &str) -> Option<MetarWind> {
    if token.len() < 7 {
        return None;
    }
    let kt_pos = token.find("KT")?;
    let body = &token[..kt_pos];
    let variable = body.starts_with("VRB");
    let (dir, rest) = if variable {
        (0u16, &body[3..])
    } else if body.len() >= 5 {
        let dir: u16 = body[..3].parse().ok()?;
        (dir, &body[3..])
    } else {
        return None;
    };
    let (speed, gust) = if let Some(g_pos) = rest.find('G') {
        let spd: u16 = rest[..g_pos].parse().ok()?;
        let gst: u16 = rest[g_pos + 1..].parse().ok()?;
        (spd, gst)
    } else {
        let spd: u16 = rest.parse().ok()?;
        (spd, 0u16)
    };
    Some(MetarWind {
        direction: dir,
        speed,
        gust,
        variable,
    })
}
/// Parse a basic METAR string into a [`MetarReport`].
///
/// This handles the most common METAR format:
/// `METAR KJFK 271856Z 25015G25KT 10SM ... T/DP Axxxx`
pub fn parse_metar(raw: &str) -> Option<MetarReport> {
    let tokens: Vec<&str> = raw.split_whitespace().collect();
    if tokens.len() < 4 {
        return None;
    }
    let start = if tokens[0] == "METAR" || tokens[0] == "SPECI" {
        1
    } else {
        0
    };
    let station_id = tokens[start].to_string();
    let time_tok = tokens[start + 1];
    if time_tok.len() < 7 {
        return None;
    }
    let day: u8 = time_tok[..2].parse().ok()?;
    let hour: u8 = time_tok[2..4].parse().ok()?;
    let minute: u8 = time_tok[4..6].parse().ok()?;
    let wind = parse_metar_wind(tokens[start + 2])?;
    let mut visibility = 9999.0;
    let mut temperature = f64::NAN;
    let mut dewpoint = f64::NAN;
    let mut altimeter = 1013.25;
    for &tok in &tokens[start + 3..] {
        if let Some(vis_str) = tok.strip_suffix("SM")
            && let Ok(v) = vis_str.parse::<f64>()
        {
            visibility = v * 1609.34;
        }
        if tok.contains('/') && !tok.starts_with('A') && !tok.starts_with('Q') {
            let parts: Vec<&str> = tok.split('/').collect();
            if parts.len() == 2
                && let (Some(t), Some(d)) = (parse_metar_temp(parts[0]), parse_metar_temp(parts[1]))
            {
                temperature = t;
                dewpoint = d;
            }
        }
        if tok.len() == 5
            && let Some(rest) = tok.strip_prefix('A')
            && let Ok(a) = rest.parse::<f64>()
        {
            altimeter = a / 100.0 * 33.8639;
        }
        if tok.len() >= 4
            && let Some(rest) = tok.strip_prefix('Q')
            && let Ok(q) = rest.parse::<f64>()
        {
            altimeter = q;
        }
    }
    Some(MetarReport {
        station_id,
        day,
        hour,
        minute,
        wind,
        visibility,
        temperature,
        dewpoint,
        altimeter,
        raw: raw.to_string(),
    })
}
/// Parse a METAR temperature token like `"15"` or `"M02"`.
pub(super) fn parse_metar_temp(s: &str) -> Option<f64> {
    if let Some(rest) = s.strip_prefix('M') {
        rest.parse::<f64>().ok().map(|v| -v)
    } else {
        s.parse::<f64>().ok()
    }
}
/// Approximate gamma function using Stirling's approximation.
pub(super) fn gamma_approx(x: f64) -> f64 {
    if x < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma_approx(1.0 - x))
    } else {
        let t = x - 1.0 + 7.5;
        let coeffs = [
            0.999_999_999_999_809_9,
            676.520_368_121_885_1,
            -1_259.139_216_722_403,
            771.323_428_777_653_1,
            -176.615_029_162_140_6,
            12.507_343_278_686_905,
            -0.138_571_095_265_720_12,
            9.984_369_578_019_572e-6,
            1.505_632_735_149_311_6e-7,
        ];
        let mut sum = coeffs[0];
        for (i, &c) in coeffs[1..].iter().enumerate() {
            sum += c / (x + i as f64);
        }
        (2.0 * std::f64::consts::PI).sqrt() * t.powf(x - 0.5) * (-t).exp() * sum
    }
}
/// Fit GEV distribution using L-moments method.
pub(super) fn fit_gev_lmoments(data: &[f64]) -> GevParameters {
    let n = data.len();
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let nf = n as f64;
    let mut b0 = 0.0_f64;
    let mut b1 = 0.0_f64;
    let mut b2 = 0.0_f64;
    for (i, &x) in sorted.iter().enumerate() {
        let j = i as f64;
        b0 += x;
        b1 += j / (nf - 1.0) * x;
        if n > 2 {
            b2 += j * (j - 1.0) / ((nf - 1.0) * (nf - 2.0)) * x;
        }
    }
    b0 /= nf;
    b1 /= nf;
    b2 /= nf;
    let l1 = b0;
    let l2 = 2.0 * b1 - b0;
    let l3 = 6.0 * b2 - 6.0 * b1 + b0;
    if l2.abs() < 1e-15 {
        return GevParameters::new(l1, 1.0, 0.0);
    }
    let t3 = l3 / l2;
    let c = 2.0 / (3.0 + t3) - std::f64::consts::LN_2 / (2.0_f64.ln() + 3.0_f64.ln());
    let shape = 7.859 * c + 2.9554 * c * c;
    let shape = shape.clamp(-0.5, 0.5);
    let scale = if shape.abs() < 1e-10 {
        l2 / std::f64::consts::LN_2
    } else {
        l2 * shape / (gamma_approx(1.0 + shape) * (1.0 - 2.0_f64.powf(-shape)))
    };
    let location = if shape.abs() < 1e-10 {
        l1 - scale * 0.5772
    } else {
        l1 - scale * (gamma_approx(1.0 + shape) - 1.0) / shape
    };
    GevParameters::new(location, scale.abs(), shape)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::weather_data_io::CellMethod;
    use crate::weather_data_io::CfStandardName;
    use crate::weather_data_io::CfTimeAxis;
    use crate::weather_data_io::CfVariable;
    use crate::weather_data_io::ClimateRecord;
    use crate::weather_data_io::ClimateTimeSeries;
    use crate::weather_data_io::DataRepresentation;
    use crate::weather_data_io::GribMessage;
    use crate::weather_data_io::GribReader;
    use crate::weather_data_io::GribSectionType;
    use crate::weather_data_io::GridDefinition;
    use crate::weather_data_io::NetCdfCfReader;
    use crate::weather_data_io::PrecipitationAnalysis;
    use crate::weather_data_io::SynopObservation;
    use crate::weather_data_io::TimeResolution;
    use crate::weather_data_io::WeatherStation;
    use crate::weather_data_io::WeatherWriter;
    use crate::weather_data_io::WindRose;
    #[test]
    fn test_grib_message_grid_point_count() {
        let msg = GribMessage::new(
            0,
            0,
            0,
            GridDefinition::LatLon {
                ni: 360,
                nj: 181,
                lat_first: -90.0,
                lat_last: 90.0,
                lon_first: 0.0,
                lon_last: 359.0,
                di: 1.0,
                dj: 1.0,
            },
            DataRepresentation::SimplePacking,
        );
        assert_eq!(msg.grid_point_count(), 360 * 181);
    }
    #[test]
    fn test_grib_message_valid_data_no_bitmap() {
        let mut msg = GribMessage::new(
            0,
            0,
            0,
            GridDefinition::LatLon {
                ni: 2,
                nj: 2,
                lat_first: 0.0,
                lat_last: 1.0,
                lon_first: 0.0,
                lon_last: 1.0,
                di: 1.0,
                dj: 1.0,
            },
            DataRepresentation::SimplePacking,
        );
        msg.data = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(msg.valid_data().len(), 4);
    }
    #[test]
    fn test_grib_message_valid_data_with_bitmap() {
        let mut msg = GribMessage::new(
            0,
            0,
            0,
            GridDefinition::LatLon {
                ni: 2,
                nj: 2,
                lat_first: 0.0,
                lat_last: 1.0,
                lon_first: 0.0,
                lon_last: 1.0,
                di: 1.0,
                dj: 1.0,
            },
            DataRepresentation::SimplePacking,
        );
        msg.data = vec![1.0, 2.0, 3.0, 4.0];
        msg.bitmap = Some(vec![true, false, true, false]);
        let valid = msg.valid_data();
        assert_eq!(valid, vec![1.0, 3.0]);
    }
    #[test]
    fn test_grib_parameter_description() {
        let msg = GribMessage::new(
            0,
            0,
            0,
            GridDefinition::LatLon {
                ni: 1,
                nj: 1,
                lat_first: 0.0,
                lat_last: 0.0,
                lon_first: 0.0,
                lon_last: 0.0,
                di: 1.0,
                dj: 1.0,
            },
            DataRepresentation::SimplePacking,
        );
        assert_eq!(msg.parameter_description(), "Temperature");
    }
    #[test]
    fn test_grib_section_parsing_minimal() {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"GRIB");
        buf.extend_from_slice(&[0, 0]);
        buf.push(0);
        buf.push(2);
        buf.extend_from_slice(&20u64.to_be_bytes());
        buf.extend_from_slice(b"7777");
        let mut reader = GribReader::new(&buf);
        let msg = reader.next_message();
        assert!(msg.is_some());
        let msg = msg.unwrap();
        assert_eq!(msg.discipline, 0);
        let sec_types: Vec<_> = reader.sections().iter().map(|s| s.section_type).collect();
        assert!(sec_types.contains(&GribSectionType::Indicator));
        assert!(sec_types.contains(&GribSectionType::End));
    }
    #[test]
    fn test_grib_lambert_grid_point_count() {
        let msg = GribMessage::new(
            0,
            2,
            2,
            GridDefinition::LambertConformal {
                nx: 100,
                ny: 100,
                lat_ref: 45.0,
                lon_ref: -100.0,
                latin1: 30.0,
                latin2: 60.0,
                dx: 3000.0,
                dy: 3000.0,
            },
            DataRepresentation::ComplexPacking,
        );
        assert_eq!(msg.grid_point_count(), 10000);
    }
    #[test]
    fn test_metar_wind_simple() {
        let wind = parse_metar_wind("25015KT").unwrap();
        assert_eq!(wind.direction, 250);
        assert_eq!(wind.speed, 15);
        assert_eq!(wind.gust, 0);
        assert!(!wind.variable);
    }
    #[test]
    fn test_metar_wind_with_gust() {
        let wind = parse_metar_wind("18020G35KT").unwrap();
        assert_eq!(wind.direction, 180);
        assert_eq!(wind.speed, 20);
        assert_eq!(wind.gust, 35);
    }
    #[test]
    fn test_metar_wind_variable() {
        let wind = parse_metar_wind("VRB05KT").unwrap();
        assert!(wind.variable);
        assert_eq!(wind.speed, 5);
    }
    #[test]
    fn test_metar_parse_full() {
        let raw = "METAR KJFK 271856Z 25015G25KT 10SM FEW250 15/M02 A3012";
        let report = parse_metar(raw).unwrap();
        assert_eq!(report.station_id, "KJFK");
        assert_eq!(report.day, 27);
        assert_eq!(report.hour, 18);
        assert_eq!(report.minute, 56);
        assert_eq!(report.wind.direction, 250);
        assert_eq!(report.wind.speed, 15);
        assert!((report.temperature - 15.0).abs() < 0.1);
        assert!((report.dewpoint - (-2.0)).abs() < 0.1);
    }
    #[test]
    fn test_metar_invalid_short() {
        assert!(parse_metar("METAR").is_none());
    }
    #[test]
    fn test_cf_time_axis_base_date() {
        let axis = CfTimeAxis::new("days since 1850-01-01", "standard", vec![0.0, 365.0]);
        let (y, m, d) = axis.base_date().unwrap();
        assert_eq!(y, 1850);
        assert_eq!(m, 1);
        assert_eq!(d, 1);
    }
    #[test]
    fn test_cf_time_axis_seconds_per_unit() {
        let axis = CfTimeAxis::new("hours since 2000-01-01", "standard", vec![0.0]);
        assert!((axis.seconds_per_unit() - 3600.0).abs() < 1e-10);
    }
    #[test]
    fn test_cf_time_axis_to_seconds() {
        let axis = CfTimeAxis::new("days since 2000-01-01", "standard", vec![1.0, 2.5]);
        let secs = axis.to_seconds();
        assert!((secs[0] - 86400.0).abs() < 1e-6);
        assert!((secs[1] - 216000.0).abs() < 1e-6);
    }
    #[test]
    fn test_cf_standard_name_roundtrip() {
        let name = CfStandardName::from_keyword("air_temperature");
        assert_eq!(name.as_str(), "air_temperature");
        assert_eq!(name, CfStandardName::AirTemperature);
    }
    #[test]
    fn test_netcdf_reader_find_variable() {
        let mut reader = NetCdfCfReader::new();
        reader.add_variable(CfVariable {
            name: "tas".to_string(),
            standard_name: CfStandardName::AirTemperature,
            units: "K".to_string(),
            cell_method: CellMethod::Mean,
            dimensions: vec!["time".to_string(), "lat".to_string(), "lon".to_string()],
            fill_value: -9999.0,
            data: vec![280.0, 285.0],
        });
        let found = reader.find_by_standard_name(&CfStandardName::AirTemperature);
        assert!(found.is_some());
        assert_eq!(found.unwrap().name, "tas");
    }
    #[test]
    fn test_climate_mean_temperature() {
        let mut ts = ClimateTimeSeries::new("test", TimeResolution::Monthly);
        for i in 1..=12 {
            ts.add_record(ClimateRecord {
                year: 2020,
                month: i,
                day: 0,
                temperature: i as f64,
                precipitation: 50.0,
            });
        }
        let mean = ts.mean_temperature().unwrap();
        assert!((mean - 6.5).abs() < 1e-10);
    }
    #[test]
    fn test_climate_temperature_anomalies() {
        let mut ts = ClimateTimeSeries::new("test", TimeResolution::Annual);
        ts.set_baseline(2000, 2002);
        for y in 2000..=2005 {
            ts.add_record(ClimateRecord {
                year: y,
                month: 0,
                day: 0,
                temperature: 10.0 + (y - 2000) as f64,
                precipitation: 0.0,
            });
        }
        let anomalies = ts.temperature_anomalies();
        assert!((anomalies[0] - (-1.0)).abs() < 1e-10);
        assert!((anomalies[3] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_climate_temperature_trend() {
        let mut ts = ClimateTimeSeries::new("test", TimeResolution::Annual);
        for i in 0..10 {
            ts.add_record(ClimateRecord {
                year: 2000 + i,
                month: 0,
                day: 0,
                temperature: 10.0 + 0.5 * i as f64,
                precipitation: 0.0,
            });
        }
        let (slope, _intercept) = ts.temperature_trend().unwrap();
        assert!((slope - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_climate_monthly_means() {
        let mut ts = ClimateTimeSeries::new("test", TimeResolution::Monthly);
        for _yr in 0..2 {
            for m in 1..=12 {
                ts.add_record(ClimateRecord {
                    year: 2020,
                    month: m,
                    day: 0,
                    temperature: m as f64 * 2.0,
                    precipitation: 0.0,
                });
            }
        }
        let means = ts.monthly_means();
        assert!((means[0] - 2.0).abs() < 1e-10);
        assert!((means[11] - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_wind_rose_frequency_sum() {
        let mut wr = WindRose::new(16, vec![0.0, 5.0, 10.0, 15.0], 0.5);
        wr.add_observation(0.0, 3.0);
        wr.add_observation(90.0, 7.0);
        wr.add_observation(180.0, 12.0);
        wr.add_observation(270.0, 0.1);
        assert_eq!(wr.frequency_sum(), 3);
        assert_eq!(wr.calm_count, 1);
        assert_eq!(wr.total_count, 4);
    }
    #[test]
    fn test_wind_rose_direction_frequency() {
        let mut wr = WindRose::new(4, vec![0.0, 10.0, 20.0], 0.5);
        for _ in 0..10 {
            wr.add_observation(0.0, 5.0);
        }
        let freq = wr.direction_frequency(0);
        assert!((freq - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_wind_rose_weibull_fit() {
        let wr = WindRose::new(16, vec![0.0, 5.0, 10.0], 0.5);
        let obs: Vec<(f64, f64)> = (0..100).map(|i| (0.0, 2.0 + (i as f64) * 0.1)).collect();
        let (k, c) = wr.weibull_fit(&obs);
        assert!(k > 0.0);
        assert!(c > 0.0);
    }
    #[test]
    fn test_wind_rose_mean_direction() {
        let mut wr = WindRose::new(4, vec![0.0, 20.0], 0.0);
        for _ in 0..10 {
            wr.add_observation(90.0, 5.0);
        }
        let mean_dir = wr.mean_direction();
        assert!((mean_dir - 90.0).abs() < 1.0);
    }
    #[test]
    fn test_gev_quantile_gumbel() {
        let params = GevParameters::new(100.0, 10.0, 0.0);
        let q = params.return_level(100.0);
        assert!(q > 100.0);
    }
    #[test]
    fn test_gev_pdf_positive() {
        let params = GevParameters::new(50.0, 10.0, 0.1);
        let p = params.pdf(55.0);
        assert!(p > 0.0);
    }
    #[test]
    fn test_idf_curve_monotonicity() {
        let mut analysis = PrecipitationAnalysis::new();
        for &dur in &[15u32, 30, 60, 120, 360] {
            for i in 0..20 {
                let depth = (dur as f64).sqrt() * (10.0 + i as f64 * 0.5);
                analysis.add_annual_max(dur, depth);
            }
        }
        analysis.fit_all();
        analysis.compute_idf(&[10.0, 50.0, 100.0]);
        for &rp in &[10.0, 50.0, 100.0] {
            assert!(analysis.check_idf_monotonicity(rp));
        }
    }
    #[test]
    fn test_gev_fit_basic() {
        let data = vec![10.0, 12.0, 15.0, 11.0, 13.0, 14.0, 16.0, 10.5, 12.5, 13.5];
        let params = fit_gev_lmoments(&data);
        assert!(params.scale > 0.0);
        let mean: f64 = data.iter().sum::<f64>() / data.len() as f64;
        assert!((params.location - mean).abs() < 10.0);
    }
    #[test]
    fn test_weather_station_mean_temperature() {
        let mut station = WeatherStation::new("TEST01", "Test Station", 35.0, 139.0, 40.0);
        for t in [280.0, 285.0, 290.0] {
            station.add_observation(SynopObservation {
                time: 0,
                temperature: t,
                dewpoint: 270.0,
                pressure: 101325.0,
                wind_direction: 180.0,
                wind_speed: 5.0,
                weather_code: 0,
                cloud_cover: 4,
                visibility: 10000.0,
                precipitation: 0.0,
            });
        }
        let mean = station.mean_temperature().unwrap();
        assert!((mean - 285.0).abs() < 1e-10);
    }
    #[test]
    fn test_weather_station_max_wind() {
        let mut station = WeatherStation::new("TEST02", "Windy", 50.0, 0.0, 10.0);
        for ws in [5.0, 15.0, 10.0] {
            station.add_observation(SynopObservation {
                time: 0,
                temperature: 280.0,
                dewpoint: 270.0,
                pressure: 101325.0,
                wind_direction: 270.0,
                wind_speed: ws,
                weather_code: 0,
                cloud_cover: 2,
                visibility: 10000.0,
                precipitation: 0.0,
            });
        }
        assert!((station.max_wind_speed().unwrap() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_writer_station_csv_has_header() {
        let station = WeatherStation::new("S01", "Station One", 10.0, 20.0, 100.0);
        let writer = WeatherWriter::new();
        let csv = writer.write_station_csv(&station);
        assert!(csv.contains("Station One"));
        assert!(csv.contains("temperature_K"));
    }
    #[test]
    fn test_writer_format_nan() {
        let writer = WeatherWriter::new();
        assert_eq!(writer.format_value(f64::NAN), "NaN");
        assert_eq!(writer.format_value(2.54159), "2.54");
    }
    #[test]
    fn test_writer_gridded_output() {
        let mut msg = GribMessage::new(
            0,
            0,
            0,
            GridDefinition::LatLon {
                ni: 2,
                nj: 2,
                lat_first: 0.0,
                lat_last: 1.0,
                lon_first: 0.0,
                lon_last: 1.0,
                di: 1.0,
                dj: 1.0,
            },
            DataRepresentation::SimplePacking,
        );
        msg.data = vec![10.0, 20.0, 30.0, 40.0];
        let writer = WeatherWriter::new();
        let out = writer.write_gridded(&msg);
        assert!(out.contains("Temperature"));
        assert!(out.contains("10.00"));
    }
    #[test]
    fn test_writer_wind_rose_csv() {
        let mut wr = WindRose::new(4, vec![0.0, 5.0, 10.0], 0.5);
        wr.add_observation(0.0, 3.0);
        let writer = WeatherWriter::new();
        let csv = writer.write_wind_rose_csv(&wr);
        assert!(csv.contains("Wind Rose"));
        assert!(csv.contains("direction_deg"));
    }
}
