use std::f64::consts::PI;

pub(crate) const HALF_PI: f64 = PI / 2.0;
pub(crate) const TWO_PI: f64 = PI + PI;

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
    pub(crate) fn unit_normal_vector(&self, other: &Self) -> Self {
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

    pub(crate) fn transform(&self, ax: &AxisTransformation) -> Self {
        ax.transform(&self)
    }

    pub(crate) fn transform_inverse(&self, ax: &AxisTransformation) -> Self {
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
pub(crate) struct AxisTransformation {
    theta: f64,
    phi: f64,
}

impl AxisTransformation {
    pub(crate) fn new(axis_vec: Xyz) -> Self {
        let LatLonInRadians(theta, phi) = LatLonInRadians::from(&axis_vec);
        let theta = HALF_PI - theta; // from equator -> from north pole
        let (theta, phi) = (-theta, -phi); // for inverse transformation
        Self { theta, phi }
    }

    pub(crate) fn transform(&self, xyz: &Xyz) -> Xyz {
        xyz.rotate_around_z_axis(self.phi.sin(), self.phi.cos())
            .rotate_around_y_axis(self.theta.sin(), self.theta.cos())
    }

    pub(crate) fn transform_inverse(&self, xyz: &Xyz) -> Xyz {
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
}
