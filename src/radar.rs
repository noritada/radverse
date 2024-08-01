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

#[derive(Debug)]
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
// X = (r_eff + z) * sin(s / r_eff)
// Y = sqrt(z * z - X * X) - sr
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
    let x = (r_eff + point.alt_meter) * (point.dist_meter / r_eff).sin();
    let y = (point.alt_meter * point.alt_meter - x * x).sqrt() - sr;
    let r_meter = (x * x + y * y).sqrt();
    let el_deg = y.atan2(x).to_degrees();
    RadarObsCellVertical { r_meter, el_deg }
}
