const HALF_PI: f64 = std::f64::consts::PI / 2.0;

#[derive(Debug, PartialEq)]
pub struct LatLonInDegrees(pub f64, pub f64);

#[derive(Debug, PartialEq)]
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

    fn rotate_around_x_axis(&self, sin_theta: f64, cos_theta: f64) -> Self {
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

#[derive(Debug, PartialEq)]
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

fn compute_distance_from_great_circle_route<I>(
    latlons: I,
    loc1: &LatLonInRadians,
    loc2: &LatLonInRadians,
) -> std::iter::Map<I, impl FnMut(LatLonInRadians) -> f64>
where
    I: Iterator<Item = LatLonInRadians>,
{
    let lat_center = (loc1.0 + loc2.0) / 2.0;
    let earth_radius = crate::earth::calc_earth_radius(lat_center);

    let xyz1 = Xyz::from(loc1);
    let xyz2 = Xyz::from(loc2);
    let normal_vec = xyz1.unit_normal_vector(&xyz2);

    latlons.map(move |ll_rad| {
        let xyz = Xyz::from(&ll_rad);
        let dot_product = normal_vec.dot_product(&xyz);
        // arccos returns the angle between the normal vector and vector to the point
        let distance_rad = dot_product.acos() - HALF_PI;
        distance_rad * earth_radius
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_almost_eq {
        ($a1:expr, $a2:expr, $d:expr) => {
            if $a1 - $a2 > $d || $a2 - $a1 > $d {
                panic!();
            }
        };
    }

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

    #[test]
    fn computation_of_distance_from_great_circle_route() {
        let latlons = [
            (35.0_f64, 140.0_f64),
            (36.0, 140.0),
            (35.0, 141.0),
            (36.0, 141.0),
        ];
        let loc1 = LatLonInDegrees(35.0, 140.0);
        let loc2 = LatLonInDegrees(36.0, 141.0);
        let actual_distances = compute_distance_from_great_circle_route(
            latlons
                .into_iter()
                .map(|(lat, lon)| LatLonInRadians::from(&LatLonInDegrees(lat, lon))),
            &LatLonInRadians::from(&loc1),
            &LatLonInRadians::from(&loc2),
        )
        .map(|d| d.abs())
        .collect::<Vec<_>>();
        let expected_distances = vec![0., 70_000., 70_000., 0.];
        assert_almost_eq!(actual_distances[0], expected_distances[0], 1.0e+3);
        assert_almost_eq!(actual_distances[1], expected_distances[1], 1.0e+4);
        assert_almost_eq!(actual_distances[2], expected_distances[2], 1.0e+4);
        assert_almost_eq!(actual_distances[3], expected_distances[3], 1.0e+3);
    }
}
