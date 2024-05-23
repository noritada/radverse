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

fn compute_distance_from_great_circle_route<I1, I2>(
    lats_rad: I1,
    lons_rad: I2,
    loc1: &LatLonInRadians,
    loc2: &LatLonInRadians,
) -> std::iter::Map<std::iter::Zip<I1, I2>, impl FnMut((f64, f64)) -> f64>
where
    I1: Iterator<Item = f64>,
    I2: Iterator<Item = f64>,
{
    let lat_center = (loc1.0 + loc2.0) / 2.0;
    let earth_radius = crate::earth::calc_earth_radius(lat_center);

    let xyz1 = Xyz::from(loc1);
    let xyz2 = Xyz::from(loc2);
    let normal_vec = xyz1.unit_normal_vector(&xyz2);

    lats_rad.zip(lons_rad).map(move |(lat, lon)| {
        let xyz = Xyz::from(&LatLonInRadians(lat, lon));
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

    #[test]
    fn computation_of_distance_from_great_circle_route() {
        let lats = [35.0_f64, 36.0, 35.0, 36.0];
        let lons = [140.0_f64, 140.0, 141.0, 141.0];
        let loc1 = LatLonInDegrees(35.0, 140.0);
        let loc2 = LatLonInDegrees(36.0, 141.0);
        let actual_distances = compute_distance_from_great_circle_route(
            lats.into_iter().map(|f| f.to_radians()),
            lons.into_iter().map(|f| f.to_radians()),
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
