#[derive(Debug)]
pub struct RadarSite {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub alt_meter: f64,
}

#[derive(Debug, Clone)]
pub struct RadarObsCell {
    pub r_meter: f64,
    pub az_deg: f64,
    pub el_deg: f64,
}

impl RadarObsCell {
    pub fn new(r_meter: f64, az_deg: f64, el_deg: f64) -> Self {
        Self {
            r_meter,
            az_deg,
            el_deg,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct RadarObsCellVertical {
    pub r_meter: f64,
    pub el_deg: f64,
}

impl From<(&RadarCenteredPoint, &RadarSite)> for RadarObsCellVertical {
    fn from(value: (&RadarCenteredPoint, &RadarSite)) -> Self {
        let (point, site) = value;
        calc_altitude_and_distance_on_sphere_inverse(point, site)
    }
}

#[derive(Debug)]
pub struct RadarCenteredPoint {
    pub alt_meter: f64,
    pub dist_meter: f64,
}

impl From<(&RadarObsCellVertical, &RadarSite)> for RadarCenteredPoint {
    fn from(value: (&RadarObsCellVertical, &RadarSite)) -> Self {
        let (cell, site) = value;
        calc_altitude_and_distance_on_sphere(
            cell.r_meter,
            cell.el_deg,
            site.lat_deg,
            site.alt_meter,
        )
    }
}

// References:
//
// - https://docs.wradlib.org/en/stable/generated/wradlib.georef.misc.site_distance.html
// - https://docs.wradlib.org/en/stable/generated/wradlib.georef.misc.bin_altitude.html
fn calc_altitude_and_distance_on_sphere(
    r_meter: f64,
    el_deg: f64,
    lat_deg: f64,
    alt_meter: f64,
) -> RadarCenteredPoint {
    let el = el_deg.to_radians();
    let r_earth = crate::earth::calc_earth_radius(lat_deg.to_radians());
    let r_eff = r_earth * 4_f64 / 3_f64;
    let sr = r_eff + alt_meter;
    let z = (r_meter * r_meter + sr * sr + r_meter * sr * 2_f64 * el.sin()).sqrt() - r_eff;
    let s = r_eff * ((r_meter * el.cos()) / (r_eff + z)).asin();
    RadarCenteredPoint {
        alt_meter: z,
        dist_meter: s,
    }
}

// Inverse function:
// ```
// X = Z * sin(s / r_eff)
// Y = sqrt(Z * Z - X * X) - sr
// Z = r_eff + z
// X = r * cos(el)
// Y = r * sin(el)
// ```
fn calc_altitude_and_distance_on_sphere_inverse(
    point: &RadarCenteredPoint,
    site: &RadarSite,
) -> RadarObsCellVertical {
    let r_earth = crate::earth::calc_earth_radius(site.lat_deg.to_radians());
    let r_eff = r_earth * 4_f64 / 3_f64;
    let sr = r_eff + site.alt_meter;
    let z = r_eff + point.alt_meter;
    let x = z * (point.dist_meter / r_eff).sin();
    let y = (z * z - x * x).sqrt() - sr;
    let r_meter = (x * x + y * y).sqrt();
    let el_deg = y.atan2(x).to_degrees();
    RadarObsCellVertical { r_meter, el_deg }
}

// Formula:
// ```
// Z = r_eff + z
// theta = s / r_eff
// X = Z * sin(theta)
// Y = Z * cos(theta) - sr
// X = r * cos(el)
// Y = r * sin(el)
// r = sr * sin(theta) / cos(el + theta)
// Z = sr * cos(el) / cos(el + theta)
// ```
pub fn calc_r_and_z_from_s_and_el(s_meter: f64, el_deg: f64, site: &RadarSite) -> (f64, f64) {
    let el = el_deg.to_radians();
    let r_earth = crate::earth::calc_earth_radius(site.lat_deg.to_radians());
    let r_eff = r_earth * 4_f64 / 3_f64;
    let sr = r_eff + site.alt_meter;
    let theta = s_meter / r_eff;
    let c = sr / (el + theta).cos();
    let r = c * theta.sin();
    let z = c * el.cos() - r_eff;
    (r, z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::*;

    #[test]
    fn altitude_and_distance_calculation() {
        let site = RadarSite {
            lat_deg: 36.,
            lon_deg: 140.,
            alt_meter: 30.,
        };
        let point_original = RadarObsCellVertical {
            r_meter: 50_000.,
            el_deg: 0.,
        };
        let radar_centered = RadarCenteredPoint::from((&point_original, &site));
        let point_calculated = RadarObsCellVertical::from((&radar_centered, &site));
        assert_eq!(point_calculated, point_original);
    }

    #[test]
    fn r_and_z_calculation() {
        let site = RadarSite {
            lat_deg: 36.,
            lon_deg: 140.,
            alt_meter: 30.,
        };
        let (s_meter, el_deg) = (50_000., 10.);
        let (r_meter, z_meter) = calc_r_and_z_from_s_and_el(s_meter, el_deg, &site);

        let point = RadarObsCellVertical { r_meter, el_deg };
        let radar_centered = RadarCenteredPoint::from((&point, &site));
        assert_almost_eq!(radar_centered.dist_meter, s_meter, 1.0e-10);

        let radar_centered = RadarCenteredPoint {
            alt_meter: z_meter,
            dist_meter: s_meter,
        };
        let point_calculated = RadarObsCellVertical::from((&radar_centered, &site));
        assert_almost_eq!(point_calculated.el_deg, el_deg, 1.0e-10);
    }
}
