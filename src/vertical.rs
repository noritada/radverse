use std::{ops::Range, vec::IntoIter};

use itertools::Itertools;

use crate::{
    calc_distance_and_direction, AxisTransformation, LatLonInDegrees, LatLonInRadians,
    RadarCenteredPoint, RadarObsCell, RadarObsCellVertical, RadarSite, Xyz, HALF_PI, TWO_PI,
};

#[derive(Debug, PartialEq)]
pub struct VerticalCrossSection {
    pub path_points: Vec<VerticalCrossSectionHorizontalPoint>,
    v_axis: VerticalCrossSectionVerticalAxis,
    site: RadarSite,
    pub shape: (usize, usize),
    pub max_distance_meter: f64,
    pub max_alt_km: u8,
}

impl VerticalCrossSection {
    pub fn new(
        path: &[LatLonInRadians],
        max_alt_km: u8,
        (width, height): (usize, usize),
        site: &RadarSite,
    ) -> Option<Self> {
        let earth_radius = crate::earth::calc_earth_radius(site.lat_deg.to_degrees());
        let h_axis = VerticalCrossSectionHorizontalAxis::from(&path, width, site)?;
        let v_axis = VerticalCrossSectionVerticalAxis::from(max_alt_km, height);
        let h_cells = h_axis
            .cells
            .iter()
            .map(|point| {
                let mut new_point = point.clone();
                new_point.site_distance = point.site_distance * earth_radius;
                new_point
            })
            .collect::<Vec<_>>();

        let [s_start, s_end] = h_axis.phi_bounds;
        let max_distance_meter = (s_start - s_end).abs() * earth_radius;
        Some(Self {
            path_points: h_cells,
            v_axis,
            site: site.clone(),
            shape: (width, height),
            max_distance_meter,
            max_alt_km,
        })
    }

    pub fn cells(&self) -> Vec<RadarObsCell> {
        let VerticalCrossSectionVerticalAxis(ref v_cells) = self.v_axis;
        v_cells
            .iter()
            .cartesian_product(self.path_points.iter())
            .map(|(&alt_meter, point)| {
                let cell = RadarCenteredPoint {
                    alt_meter,
                    dist_meter: point.site_distance,
                };
                let cell = RadarObsCellVertical::from((&cell, &self.site));
                let az_deg = point.site_direction.to_degrees();
                RadarObsCell::new(cell.r_meter, az_deg, cell.el_deg)
            })
            .collect()
    }
}

#[derive(Debug, PartialEq)]
pub struct VerticalCrossSectionHorizontalAxis {
    cells: Vec<VerticalCrossSectionHorizontalPoint>,
    phi_bounds: [f64; 2],
}

impl VerticalCrossSectionHorizontalAxis {
    // `n` is a number of pixels between `loc1` and `loc2`.
    // `cells` will have `n + 1` length.
    pub fn from(path: &[LatLonInRadians], n: usize, site: &RadarSite) -> Option<Self> {
        let path = path
            .into_iter()
            .map(|waypoint| Xyz::from(waypoint))
            .tuple_windows::<(_, _)>();
        let mut phi_total_start = 0_f64;
        let segments = path
            .map(|(start, end)| {
                let segment = PathSegment::from(&start, &end);
                let phi_total_end = phi_total_start + segment.len();
                phi_total_start = phi_total_end;
                (phi_total_end, segment)
            })
            .collect::<Vec<_>>();

        let points = PathPoints::new(&segments, n, site)?;
        let length = points.length;

        Some(Self {
            cells: points.collect(),
            phi_bounds: [0., length],
        })
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct VerticalCrossSectionHorizontalPoint {
    phi: f64,
    pub site_distance: f64,
    pub site_direction: f64,
}

#[derive(Debug, PartialEq)]
pub struct PathInDegrees(pub Vec<LatLonInDegrees>);

impl std::str::FromStr for PathInDegrees {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lines = s
            .lines()
            .map(|line| {
                let (lat, lon) = line.trim_end().split_once(",")?;
                Some(LatLonInDegrees(
                    lat.trim().parse::<f64>().ok()?,
                    lon.trim().parse::<f64>().ok()?,
                ))
            })
            .collect::<Option<Vec<_>>>();
        lines.map(|vec| Self(vec)).ok_or_else(|| {
            "invalid format; path should be specified as lines of comma-separated lat/lon values"
        })
    }
}

#[derive(Debug, Clone)]
struct PathSegment {
    tx: AxisTransformation,
    phi_start: f64,
    phi_end: f64,
}

impl PathSegment {
    fn from(start: &Xyz, end: &Xyz) -> Self {
        let normal_vec = start.unit_normal_vector(&end);
        let tx = AxisTransformation::new(normal_vec);
        let project = |xyz: &Xyz| {
            let xyz = xyz.transform(&tx);
            LatLonInRadians::from(&xyz)
        };

        let LatLonInRadians(_theta, phi_start) = project(&start);
        let LatLonInRadians(_theta, phi_end) = project(&end);
        Self {
            tx,
            phi_start,
            phi_end,
        }
    }

    fn xyz(&self, d_phi: f64) -> Xyz {
        let phi = self.phi_start + d_phi;
        Xyz::from(&LatLonInRadians(0.0, phi)).transform_inverse(&self.tx)
    }

    fn len(&self) -> f64 {
        self.phi_end - self.phi_start
    }
}

struct PathPoints<'s> {
    segments: std::slice::Iter<'s, (f64, PathSegment)>,
    current: Option<&'s (f64, PathSegment)>,
    length: f64,
    points: IntoIter<f64>,
    site: LatLonInRadians,
}

impl<'s> PathPoints<'s> {
    fn new(segments: &'s [(f64, PathSegment)], n: usize, site: &RadarSite) -> Option<Self> {
        let (length, _) = segments.last()?;
        if *length == 0. {
            return None;
        }

        let mut segments = segments.iter();
        let current = segments.next();

        let phi_inc = length / n as f64;
        let points = (0..=n)
            .map(|k| phi_inc * k as f64)
            .collect::<Vec<_>>()
            .into_iter();

        let site = LatLonInRadians::from(&LatLonInDegrees(site.lat_deg, site.lon_deg));

        Some(Self {
            segments,
            current,
            length: *length,
            points,
            site,
        })
    }
}

impl<'s> Iterator for PathPoints<'s> {
    type Item = VerticalCrossSectionHorizontalPoint;

    fn next(&mut self) -> Option<Self::Item> {
        let (mut end, mut segment) = self.current.unwrap().clone(); // safely unwrapped
        let phi = self.points.next()?;
        while phi > end + f64::EPSILON {
            self.current = self.segments.next();
            (end, segment) = self.current.unwrap().clone(); // safely unwrapped
        }

        let d_phi = (segment.phi_end - segment.phi_start) - (end - phi);
        let xyz = segment.xyz(d_phi);
        let (site_distance, site_direction) =
            calc_distance_and_direction(&self.site, &LatLonInRadians::from(&xyz));
        // convert "anti-clockwise from east" to "clockwise from north"
        let site_direction = HALF_PI - site_direction;
        // normalize into `[0, TWO_PI)`
        let site_direction = if site_direction < 0. {
            site_direction + TWO_PI
        } else if site_direction < TWO_PI {
            site_direction
        } else {
            site_direction - TWO_PI
        };
        Some(VerticalCrossSectionHorizontalPoint {
            phi,
            site_distance,
            site_direction,
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.points.size_hint()
    }
}

#[derive(Debug, PartialEq)]
pub struct VerticalCrossSectionVerticalAxis(Vec<f64>);

impl VerticalCrossSectionVerticalAxis {
    // `n` is a number of pixels between the lower and upper altitude.
    // `cells` will have `n + 1` length.
    pub fn from(max_alt_km: u8, n: usize) -> Self {
        let alt_inc_m = max_alt_km as f64 * 1000_f64 / n as f64;
        let cells = (0..=n).rev().map(|k| alt_inc_m * k as f64).collect();
        Self(cells)
    }
}

#[derive(Debug, PartialEq)]
pub struct ElevationRanges(Vec<(Range<f64>, Option<Elevation>)>);

impl ElevationRanges {
    pub fn from(elevations_sorted: &[Elevation]) -> Self {
        let mut iter = elevations_sorted
            .into_iter()
            .map(|el| (el, el.angle - el.width, el.angle + el.width))
            .peekable();
        let mut vec = Vec::new();
        let mut boundary = None;

        while let Some((el, start, end)) = iter.next() {
            let start = boundary.unwrap_or(start);
            let end = if let Some((_, next_start, _)) = iter.peek() {
                if end < *next_start {
                    boundary = None;
                    end
                } else {
                    let new_end = (end + next_start) / 2.;
                    boundary = Some(new_end);
                    new_end
                }
            } else {
                end
            };
            vec.push((start..end, Some(el.clone())));

            if boundary.is_none() {
                if let Some((_, next_start, _)) = iter.peek() {
                    vec.push((end..*next_start, None));
                }
            }
        }

        Self(vec)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Elevation {
    angle: f64,
    width: f64,
}

impl Elevation {
    pub fn new(angle: f64, width: f64) -> Self {
        Self { angle, width }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! test_parsing_path_specifiers {
        ($((
            $name:ident,
            $spec:expr,
            $expected:expr,
        ),)*) => ($(
            #[test]
            fn $name() {
                let spec = $spec;
                let actual = spec.parse::<PathInDegrees>();
                let expected = $expected;
                assert_eq!(actual, expected);
            }
        )*);
    }

    test_parsing_path_specifiers! {
        (
            parsing_path_specifier_succeeding_for_empty_string,
            "",
            Ok(PathInDegrees(Vec::new())),
        ),
        (
            parsing_path_specifier_succeeding,
            "36.0,140.0\n37.0,141.0",
            Ok(PathInDegrees(vec![
                LatLonInDegrees(36.0, 140.0),
                LatLonInDegrees(37.0, 141.0),
            ])),
        ),
        (
            parsing_path_specifier_succeeding_with_trailing_newline,
            "36.0,140.0\n37.0,141.0\n",
            Ok(PathInDegrees(vec![
                LatLonInDegrees(36.0, 140.0),
                LatLonInDegrees(37.0, 141.0),
            ])),
        ),
        (
            parsing_path_specifier_succeeding_with_whitespaces,
            " 36.0 , 140.0 \n 37.0 , 141.0 \n",
            Ok(PathInDegrees(vec![
                LatLonInDegrees(36.0, 140.0),
                LatLonInDegrees(37.0, 141.0),
            ])),
        ),
        // should be fixed?
        (
            parsing_path_specifier_failing_for_wrong_lat_lon_values,
            "136.0,140.0\n137.0,141.0\n",
            Ok(PathInDegrees(vec![
                LatLonInDegrees(136.0, 140.0),
                LatLonInDegrees(137.0, 141.0),
            ])),
        ),
        (
            parsing_path_specifier_failing_for_nonnumerical_values,
            "foo\nbar",
            Err(
                "invalid format; path should be specified as lines of comma-separated lat/lon values"
            ),
        ),
    }

    #[test]
    fn regression_test_for_crashes_caused_by_rounding_errors() {
        let path_spec =
            "-30.61788354116184,-40.608893002640784\n-30.87973113791563,-40.8162709218129";
        let PathInDegrees(path) = path_spec.parse::<PathInDegrees>().unwrap();
        let path: Vec<LatLonInRadians> = path.iter().map(|loc| loc.into()).collect();
        let width = 689;
        let site = RadarSite {
            lat_deg: 34.17444568059641,
            lon_deg: 139.17029552581388,
            alt_meter: 87.9,
        };
        let h_axis = VerticalCrossSectionHorizontalAxis::from(&path, width, &site);
        assert!(h_axis.is_some())
    }

    macro_rules! test_elevation_ranges {
        ($((
            $name:ident,
            $elevations:expr,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let elevations = $elevations;
                let actual = ElevationRanges::from(&elevations);
                let expected = ElevationRanges($expected);
                assert_eq!(actual, expected);
            }
        )*);
    }

    test_elevation_ranges! {
        (
            elevation_ranges_for_nonoverlapping_elevations,
            vec![Elevation::new(2., 1.), Elevation::new(5., 1.)],
            vec![
                (1.0..3.0, Some(Elevation::new(2., 1.))),
                (3.0..4.0, None),
                (4.0..6.0, Some(Elevation::new(5., 1.))),
            ]
        ),
        (
            elevation_ranges_for_overlapping_elevations,
            vec![Elevation::new(2., 2.), Elevation::new(5., 2.)],
            vec![
                (0.0..3.5, Some(Elevation::new(2., 2.))),
                (3.5..7.0, Some(Elevation::new(5., 2.))),
            ]
        ),
    }
}
