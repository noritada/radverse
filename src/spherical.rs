#[derive(Debug, PartialEq)]
pub struct LatLonInDegrees(pub f64, pub f64);

#[derive(Debug, PartialEq)]
pub struct Xyz(f64, f64, f64);

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

    const HALF_PI: f64 = std::f64::consts::PI / 2.0;

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
}
