use std::{f64::consts::PI, vec::IntoIter};

use itertools::Itertools;

use crate::{RadarCenteredPoint, RadarObsCell, RadarObsCellVertical, RadarSite};

pub(crate) const HALF_PI: f64 = PI / 2.0;
const TWO_PI: f64 = PI + PI;

#[derive(Debug, PartialEq, Clone)]
pub struct LatLonInDegrees(pub f64, pub f64);

impl From<&LatLonInRadians> for LatLonInDegrees {
    fn from(value: &LatLonInRadians) -> Self {
        let LatLonInRadians(lat, lon) = value;
        Self(lat.to_degrees(), lon.to_degrees())
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Xyz(f64, f64, f64);

impl Xyz {
    fn unit_normal_vector(&self, other: &Self) -> Self {
        let vec = self.cross_product(other);
        let len = vec.norm();
        Self(vec.0 / len, vec.1 / len, vec.2 / len)
    }

    fn cross_product(&self, other: &Self) -> Self {
        Self(
            self.1 * other.2 - self.2 * other.1,
            self.2 * other.0 - self.0 * other.2,
            self.0 * other.1 - self.1 * other.0,
        )
    }

    fn dot_product(&self, other: &Self) -> f64 {
        self.0 * other.0 + self.1 * other.1 + self.2 * other.2
    }

    fn norm(&self) -> f64 {
        (self.0 * self.0 + self.1 * self.1 + self.2 * self.2).sqrt()
    }

    pub(crate) fn rotate_around_x_axis(&self, sin_theta: f64, cos_theta: f64) -> Self {
        let Self(x, y, z) = &self;
        Self(
            *x,
            y * cos_theta - z * sin_theta,
            y * sin_theta + z * cos_theta,
        )
    }

    fn rotate_around_y_axis(&self, sin_theta: f64, cos_theta: f64) -> Self {
        let Self(x, y, z) = &self;
        Self(
            z * sin_theta + x * cos_theta,
            *y,
            z * cos_theta - x * sin_theta,
        )
    }

    fn rotate_around_z_axis(&self, sin_theta: f64, cos_theta: f64) -> Self {
        let Self(x, y, z) = &self;
        Self(
            x * cos_theta - y * sin_theta,
            x * sin_theta + y * cos_theta,
            *z,
        )
    }

    fn transform(&self, ax: &AxisTransformation) -> Self {
        ax.transform(&self)
    }

    fn transform_inverse(&self, ax: &AxisTransformation) -> Self {
        ax.transform_inverse(&self)
    }
}

impl From<&LatLonInRadians> for Xyz {
    fn from(value: &LatLonInRadians) -> Self {
        let LatLonInRadians(lat, lon) = value;
        Self(lon.cos() * lat.cos(), lon.sin() * lat.cos(), lat.sin())
    }
}

impl From<&LatLonInDegrees> for Xyz {
    fn from(value: &LatLonInDegrees) -> Self {
        let value = LatLonInRadians::from(value);
        Self::from(&value)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct LatLonInRadians(pub f64, pub f64);

impl From<&LatLonInDegrees> for LatLonInRadians {
    fn from(value: &LatLonInDegrees) -> Self {
        let LatLonInDegrees(lat, lon) = value;
        Self(lat.to_radians(), lon.to_radians())
    }
}

impl From<&Xyz> for LatLonInRadians {
    fn from(value: &Xyz) -> Self {
        let Xyz(x, y, z) = value;
        // assuming that sqrt(x * x + y * y + z * z) == 1
        let theta = z.asin();
        let phi = y.atan2(*x);
        Self(theta, phi)
    }
}

// Transformation of the polar axis from z-axis to the axis defined.
#[derive(Debug, Clone)]
struct AxisTransformation {
    theta: f64,
    phi: f64,
}

impl AxisTransformation {
    fn new(axis_vec: Xyz) -> Self {
        let LatLonInRadians(theta, phi) = LatLonInRadians::from(&axis_vec);
        let theta = HALF_PI - theta; // from equator -> from north pole
        let (theta, phi) = (-theta, -phi); // for inverse transformation
        Self { theta, phi }
    }

    fn transform(&self, xyz: &Xyz) -> Xyz {
        xyz.rotate_around_z_axis(self.phi.sin(), self.phi.cos())
            .rotate_around_y_axis(self.theta.sin(), self.theta.cos())
    }

    fn transform_inverse(&self, xyz: &Xyz) -> Xyz {
        xyz.rotate_around_y_axis(-self.theta.sin(), self.theta.cos())
            .rotate_around_z_axis(-self.phi.sin(), self.phi.cos())
    }
}

pub(crate) struct LocalXyzTransformation([Xyz; 3]);

impl LocalXyzTransformation {
    pub(crate) fn from(center: &LatLonInRadians) -> Self {
        let LatLonInRadians(_lat, lon) = center;
        let x = Xyz::from(&LatLonInRadians(0., lon + HALF_PI));
        let z: Xyz = center.into();
        let y = z.cross_product(&x);
        Self([x, y, z])
    }

    pub(crate) fn transform(&self, xyz: &Xyz) -> Xyz {
        let Self([ax, ay, az]) = self;
        Xyz(
            ax.dot_product(&xyz),
            ay.dot_product(&xyz),
            az.dot_product(&xyz),
        )
    }

    pub(crate) fn transform_inverse(&self, xyz: &Xyz) -> Xyz {
        let Self([ax, ay, az]) = self;
        let inv_c = 1.
            / (ax.0 * ay.1 * az.2 + ax.1 * ay.2 * az.0 + ax.2 * ay.0 * az.1
                - ax.2 * ay.1 * az.0
                - ax.1 * ay.0 * az.2
                - ax.0 * ay.2 * az.1);
        let tmp_x = ay.cross_product(&az);
        let tmp_y = az.cross_product(&ax);
        let tmp_z = ax.cross_product(&ay);
        let inv_x = Xyz(tmp_x.0, tmp_y.0, tmp_z.0);
        let inv_y = Xyz(tmp_x.1, tmp_y.1, tmp_z.1);
        let inv_z = Xyz(tmp_x.2, tmp_y.2, tmp_z.2);
        Xyz(
            inv_x.dot_product(&xyz) * inv_c,
            inv_y.dot_product(&xyz) * inv_c,
            inv_z.dot_product(&xyz) * inv_c,
        )
    }
}

#[derive(Debug, PartialEq)]
pub struct VerticalCrossSection {
    pub path_points: Vec<VerticalCrossSectionHorizontalPoint>,
    pub cells: Vec<RadarObsCell>,
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
        let VerticalCrossSectionVerticalAxis(v_cells) = v_axis;
        let cells = v_cells
            .iter()
            .cartesian_product(h_cells.iter())
            .map(|(&alt_meter, point)| {
                let cell = RadarCenteredPoint {
                    alt_meter,
                    dist_meter: point.site_distance,
                };
                let cell = RadarObsCellVertical::from((&cell, site));
                let az_deg = point.site_direction.to_degrees();
                RadarObsCell::new(cell.r_meter, az_deg, cell.el_deg)
            })
            .collect();

        let [s_start, s_end] = h_axis.phi_bounds;
        let max_distance_meter = (s_start - s_end).abs() * earth_radius;
        Some(Self {
            path_points: h_cells,
            cells,
            shape: (width, height),
            max_distance_meter,
            max_alt_km,
        })
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

/// Calculates distance and direction from `loc1` to `loc2` in radians.
/// Direction is anti-clockwise from x axis.
pub fn calc_distance_and_direction(loc1: &LatLonInRadians, loc2: &LatLonInRadians) -> (f64, f64) {
    let xyz1 = Xyz::from(loc1);
    let xyz2 = Xyz::from(loc2);
    let normal_vec = xyz1.unit_normal_vector(&xyz2);

    let tx = AxisTransformation::new(normal_vec);
    let z_axis = Xyz(0., 0., 1.);
    let LatLonInRadians(z_axis_coordp_lat, z_axis_coordp_lon) =
        LatLonInRadians::from(&z_axis.transform(&tx));
    let LatLonInRadians(_zero, loc1_coordp_lon) = LatLonInRadians::from(&xyz1.transform(&tx));
    let LatLonInRadians(_zero, loc2_coordp_lon) = LatLonInRadians::from(&xyz2.transform(&tx));

    let distance = loc2_coordp_lon - loc1_coordp_lon;
    // normalize into `[-PI, PI)`
    let distance = if distance < -PI {
        distance + TWO_PI
    } else if distance < PI {
        distance
    } else {
        distance - TWO_PI
    };

    // Angle of a vector loc1->loc2 from parallel is equal to angle of z_axis_coordp
    // from current z axis. So, calculate how z_axis_coordp is projected.
    let z_axis_coordp_lon_rel = z_axis_coordp_lon - loc1_coordp_lon;
    let xx = z_axis_coordp_lat.cos() * z_axis_coordp_lon_rel.sin();
    let yy = z_axis_coordp_lat.sin();

    (distance.abs(), xx.atan2(yy))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::*;

    macro_rules! test_xyz_to_lat_lon_rad {
        ($((
            $name:ident,
            $x:expr,
            $y:expr,
            $z:expr,
            $lat:expr,
            $lon:expr,
        ),)*) => ($(
            #[test]
            fn $name() {
                let xyz = Xyz($x, $y, $z);
                let actual = LatLonInRadians::from(&xyz);
                let expected = LatLonInRadians($lat, $lon);
                assert_eq!(actual, expected);
            }
        )*);
    }

    test_xyz_to_lat_lon_rad! {
        (xyz_to_lat_lon_rad_x, 1.0, 0.0, 0.0, 0.0, 0.0,),
        (xyz_to_lat_lon_rad_y, 0.0, 1.0, 0.0, 0.0, HALF_PI,),
        (xyz_to_lat_lon_rad_z, 0.0, 0.0, 1.0, HALF_PI, 0.0,),
    }

    macro_rules! test_lat_lon_deg_to_xyz {
        ($((
            $name:ident,
            $lat:expr,
            $lon:expr,
            $x:expr,
            $y:expr,
            $z:expr,
        ),)*) => ($(
            #[test]
            fn $name() {
                let ll = LatLonInDegrees($lat, $lon);
                let actual = Xyz::from(&ll);
                let expected = Xyz($x, $y, $z);
                assert_almost_eq!(actual.0, expected.0, 1.0e-16);
                assert_almost_eq!(actual.1, expected.1, 1.0e-16);
                assert_almost_eq!(actual.2, expected.2, 1.0e-16);
            }
        )*);
    }

    test_lat_lon_deg_to_xyz! {
        (lat_lon_deg_to_lat_lon_rad_x, 0.0, 0.0, 1.0, 0.0, 0.0,),
        (lat_lon_deg_to_lat_lon_rad_y, 0.0, 90.0, 0.0, 1.0, 0.0,),
        (lat_lon_deg_to_lat_lon_rad_z, 90.0, 0.0, 0.0, 0.0, 1.0,),
    }

    #[test]
    fn unit_normal_vector() {
        let actual = Xyz(1., 1., 0.).unit_normal_vector(&Xyz(0., 1., 1.));
        let element = (1_f64 / 3.).sqrt();
        let expected = Xyz(element, -element, element);
        assert_almost_eq!(actual.0, expected.0, 1.0e-15);
        assert_almost_eq!(actual.1, expected.1, 1.0e-15);
        assert_almost_eq!(actual.2, expected.2, 1.0e-15);
    }

    #[test]
    fn cross_product() {
        let actual = Xyz(1., 1., 0.).cross_product(&Xyz(0., 1., 1.));
        let expected = Xyz(1., -1., 1.);
        assert_eq!(actual, expected)
    }

    #[test]
    fn dot_product() {
        let actual = Xyz(1., 1., 0.).dot_product(&Xyz(0., 1., 1.));
        let expected = 1.0;
        assert_eq!(actual, expected)
    }

    macro_rules! test_rotation {
        ($((
            $name:ident,
            $origin:expr,
            $op:ident,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let theta = HALF_PI;
                let actual = $origin.$op(theta.sin(), theta.cos());
                let expected = $expected;
                assert_almost_eq!(actual.0, expected.0, 1.0e-15);
                assert_almost_eq!(actual.1, expected.1, 1.0e-15);
                assert_almost_eq!(actual.2, expected.2, 1.0e-15);
            }
        )*);
    }

    test_rotation! {
        (rotation_of_x_around_x_axis, Xyz(1., 0., 0.), rotate_around_x_axis, Xyz(1., 0., 0.)),
        (rotation_of_y_around_x_axis, Xyz(0., 1., 0.), rotate_around_x_axis, Xyz(0., 0., 1.)),
        (rotation_of_z_around_x_axis, Xyz(0., 0., 1.), rotate_around_x_axis, Xyz(0., -1., 0.)),
        (rotation_of_x_around_y_axis, Xyz(1., 0., 0.), rotate_around_y_axis, Xyz(0., 0., -1.)),
        (rotation_of_y_around_y_axis, Xyz(0., 1., 0.), rotate_around_y_axis, Xyz(0., 1., 0.)),
        (rotation_of_z_around_y_axis, Xyz(0., 0., 1.), rotate_around_y_axis, Xyz(1., 0., 0.)),
        (rotation_of_x_around_z_axis, Xyz(1., 0., 0.), rotate_around_z_axis, Xyz(0., 1., 0.)),
        (rotation_of_y_around_z_axis, Xyz(0., 1., 0.), rotate_around_z_axis, Xyz(-1., 0., 0.)),
        (rotation_of_z_around_z_axis, Xyz(0., 0., 1.), rotate_around_z_axis, Xyz(0., 0., 1.)),
    }

    macro_rules! test_axis_transformation {
        ($((
            $name:ident,
            $axis:expr,
            $input:expr,
            $op:ident,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let tx = AxisTransformation::new($axis);
                let input = $input;
                let actual = tx.$op(&input);
                let expected = $expected;
                assert_almost_eq!(actual.0, expected.0, 1.0e-15);
                assert_almost_eq!(actual.1, expected.1, 1.0e-15);
                assert_almost_eq!(actual.2, expected.2, 1.0e-15);
            }
        )*);
    }

    test_axis_transformation! {
        (
            axis_transformation_of_x_using_x_axis,
            Xyz(1., 0., 0.), Xyz(1., 0., 0.), transform, Xyz(0., 0., 1.)
        ),
        (
            axis_transformation_of_y_using_x_axis,
            Xyz(1., 0., 0.), Xyz(0., 1., 0.), transform, Xyz(0., 1., 0.)
        ),
        (
            axis_transformation_of_z_using_x_axis,
            Xyz(1., 0., 0.), Xyz(0., 0., 1.), transform, Xyz(-1., 0., 0.)
        ),
        (
            axis_transformation_of_x_using_y_axis,
            Xyz(0., 1., 0.), Xyz(1., 0., 0.), transform, Xyz(0., -1., 0.)
        ),
        (
            axis_transformation_of_y_using_y_axis,
            Xyz(0., 1., 0.), Xyz(0., 1., 0.), transform, Xyz(0., 0., 1.)
        ),
        (
            axis_transformation_of_z_using_y_axis,
            Xyz(0., 1., 0.), Xyz(0., 0., 1.), transform, Xyz(-1., 0., 0.)
        ),
        (
            axis_transformation_of_x_using_z_axis,
            Xyz(0., 0., 1.), Xyz(1., 0., 0.), transform, Xyz(1., 0., 0.)
        ),
        (
            axis_transformation_of_y_using_z_axis,
            Xyz(0., 0., 1.), Xyz(0., 1., 0.), transform, Xyz(0., 1., 0.)
        ),
        (
            axis_transformation_of_z_using_z_axis,
            Xyz(0., 0., 1.), Xyz(0., 0., 1.), transform, Xyz(0., 0., 1.)
        ),
        (
            inverse_axis_transformation_of_x_using_x_axis,
            Xyz(1., 0., 0.), Xyz(1., 0., 0.), transform_inverse, Xyz(0., 0., -1.)
        ),
        (
            inverse_axis_transformation_of_y_using_x_axis,
            Xyz(1., 0., 0.), Xyz(0., 1., 0.), transform_inverse, Xyz(0., 1., 0.)
        ),
        (
            inverse_axis_transformation_of_z_using_x_axis,
            Xyz(1., 0., 0.), Xyz(0., 0., 1.), transform_inverse, Xyz(1., 0., 0.)
        ),
        (
            inverse_axis_transformation_of_x_using_y_axis,
            Xyz(0., 1., 0.), Xyz(1., 0., 0.), transform_inverse, Xyz(0., 0., -1.)
        ),
        (
            inverse_axis_transformation_of_y_using_y_axis,
            Xyz(0., 1., 0.), Xyz(0., 1., 0.), transform_inverse, Xyz(-1., 0., 0.)
        ),
        (
            inverse_axis_transformation_of_z_using_y_axis,
            Xyz(0., 1., 0.), Xyz(0., 0., 1.), transform_inverse, Xyz(0., 1., 0.)
        ),
        (
            inverse_axis_transformation_of_x_using_z_axis,
            Xyz(0., 0., 1.), Xyz(1., 0., 0.), transform_inverse, Xyz(1., 0., 0.)
        ),
        (
            inverse_axis_transformation_of_y_using_z_axis,
            Xyz(0., 0., 1.), Xyz(0., 1., 0.), transform_inverse, Xyz(0., 1., 0.)
        ),
        (
            inverse_axis_transformation_of_z_using_z_axis,
            Xyz(0., 0., 1.), Xyz(0., 0., 1.), transform_inverse, Xyz(0., 0., 1.)
        ),
    }

    macro_rules! test_local_xyz_transformation {
        ($((
            $name:ident,
            $center:expr,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let center = $center;
                let tx = LocalXyzTransformation::from(&LatLonInRadians::from(&center));
                let LocalXyzTransformation([actual_ax, actual_ay, actual_az]) = tx;
                let [expected_ax, expected_ay, expected_az] = $expected;
                assert_almost_eq!(actual_ax.0, expected_ax.0, 1.0e-15);
                assert_almost_eq!(actual_ax.1, expected_ax.1, 1.0e-15);
                assert_almost_eq!(actual_ax.2, expected_ax.2, 1.0e-15);
                assert_almost_eq!(actual_ay.0, expected_ay.0, 1.0e-15);
                assert_almost_eq!(actual_ay.1, expected_ay.1, 1.0e-15);
                assert_almost_eq!(actual_ay.2, expected_ay.2, 1.0e-15);
                assert_almost_eq!(actual_az.0, expected_az.0, 1.0e-15);
                assert_almost_eq!(actual_az.1, expected_az.1, 1.0e-15);
                assert_almost_eq!(actual_az.2, expected_az.2, 1.0e-15);
            }
        )*);
    }

    test_local_xyz_transformation! {
        (
            local_xyz_transformation_with_x_axis_centered,
            LatLonInDegrees(0., 0.), [Xyz(0., 1., 0.), Xyz(0., 0., 1.), Xyz(1., 0., 0.)]
        ),
        (
            local_xyz_transformation_with_y_axis_centered,
            LatLonInDegrees(0., 90.), [Xyz(-1., 0., 0.), Xyz(0., 0., 1.), Xyz(0., 1., 0.)]
        ),
        (
            local_xyz_transformation_with_z_axis_centered,
            LatLonInDegrees(90., 0.), [Xyz(0., 1., 0.), Xyz(-1., 0., 0.), Xyz(0., 0., 1.)]
        ),
    }

    macro_rules! test_local_xyz_transformation_forward_operations {
        ($((
            $name:ident,
            $center:expr,
            $input:expr,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let center = $center;
                let tx = LocalXyzTransformation::from(&LatLonInRadians::from(&center));
                let input = $input;
                let actual = tx.transform(&input);
                let expected = $expected;
                assert_almost_eq!(actual.0, expected.0, 1.0e-15);
                assert_almost_eq!(actual.1, expected.1, 1.0e-15);
                assert_almost_eq!(actual.2, expected.2, 1.0e-15);
            }
        )*);
    }

    test_local_xyz_transformation_forward_operations! {
        (
            local_xyz_forward_transformation_of_x_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(1., 0., 0.), Xyz(0., 0., 1.)
        ),
        (
            local_xyz_forward_transformation_of_y_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(0., 1., 0.), Xyz(1., 0., 0.)
        ),
        (
            local_xyz_forward_transformation_of_z_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(0., 0., 1.), Xyz(0., 1., 0.)
        ),
        (
            local_xyz_forward_transformation_of_x_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(1., 0., 0.), Xyz(-1., 0., 0.)
        ),
        (
            local_xyz_forward_transformation_of_y_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(0., 1., 0.), Xyz(0., 0., 1.)
        ),
        (
            local_xyz_forward_transformation_of_z_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(0., 0., 1.), Xyz(0., 1., 0.)
        ),
        (
            local_xyz_forward_transformation_of_x_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(1., 0., 0.), Xyz(0., -1., 0.)
        ),
        (
            local_xyz_forward_transformation_of_y_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(0., 1., 0.), Xyz(1., 0., 0.)
        ),
        (
            local_xyz_forward_transformation_of_z_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(0., 0., 1.), Xyz(0., 0., 1.)
        ),
    }

    macro_rules! test_local_xyz_transformation_inverse_operations {
        ($((
            $name:ident,
            $center:expr,
            $input:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let center = $center;
                let tx = LocalXyzTransformation::from(&LatLonInRadians::from(&center));
                let input = $input;
                let actual = tx.transform_inverse(&tx.transform(&input));
                let expected = input;
                assert_almost_eq!(actual.0, expected.0, 1.0e-15);
                assert_almost_eq!(actual.1, expected.1, 1.0e-15);
                assert_almost_eq!(actual.2, expected.2, 1.0e-15);
            }
        )*);
    }

    test_local_xyz_transformation_inverse_operations! {
        (
            local_xyz_inverse_transformation_of_x_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(1., 0., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_y_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(0., 1., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_z_axis_to_center_on_x_axis,
            LatLonInDegrees(0., 0.), Xyz(0., 0., 1.)
        ),
        (
            local_xyz_inverse_transformation_of_x_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(1., 0., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_y_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(0., 1., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_z_axis_to_center_on_y_axis,
            LatLonInDegrees(0., 90.), Xyz(0., 0., 1.)
        ),
        (
            local_xyz_inverse_transformation_of_x_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(1., 0., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_y_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(0., 1., 0.)
        ),
        (
            local_xyz_inverse_transformation_of_z_axis_to_center_on_z_axis,
            LatLonInDegrees(90., 0.), Xyz(0., 0., 1.)
        ),
    }

    macro_rules! test_distance_and_direction_calculation {
        ($((
            $name:ident,
            $loc1:expr,
            $loc2:expr,
            $expected_distance:expr,
            $expected_direction:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let loc1 = $loc1;
                let loc2 = $loc2;
                let (distance, direction) = calc_distance_and_direction(
                    &LatLonInRadians::from(&loc1),
                    &LatLonInRadians::from(&loc2),
                );

                let lat_center = (loc1.0 + loc2.0) / 2.0;
                let earth_radius = crate::earth::calc_earth_radius(lat_center);
                let actual_distance = distance * earth_radius;
                let actual_direction = direction.to_degrees();
                assert_almost_eq!(actual_distance, $expected_distance, 1.0e+3);
                assert_almost_eq!(actual_direction, $expected_direction, 1.0e+1);
            }
        )*);
    }

    test_distance_and_direction_calculation! {
        (
            distance_and_direction_calculation_for_the_same_loc,
            LatLonInDegrees(0., 0.), LatLonInDegrees(0., 0.), 0., 0.
        ),
        (
            distance_and_direction_calculation_along_equator,
            LatLonInDegrees(0., 0.), LatLonInDegrees(0., 1.), 111_000., 0.
        ),
        (
            distance_and_direction_calculation_along_prime_meridian,
            LatLonInDegrees(0., 0.), LatLonInDegrees(1., 0.), 111_000., 90.
        ),
        (
            distance_and_direction_calculation_along_parallel,
            LatLonInDegrees(36., 140.), LatLonInDegrees(36., 141.),
            111_000_f64 * 36_f64.to_radians().cos(), 0.
        ),
        (
            distance_and_direction_calculation_along_meridian,
            LatLonInDegrees(36., 140.), LatLonInDegrees(37., 140.), 111_000., 90.
        ),
    }

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
}
