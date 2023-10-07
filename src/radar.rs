#[derive(Debug)]
pub(crate) struct ObservationPoint {
    pub(crate) alt_meter: f64,
    pub(crate) dist_meter: f64,
}

// References:
//
// - https://docs.wradlib.org/en/stable/generated/wradlib.georef.misc.site_distance.html
// - https://docs.wradlib.org/en/stable/generated/wradlib.georef.misc.bin_altitude.html
pub(crate) fn calc_altitude_and_distance_on_sphere(
    r_meter: f64,
    el_deg: f64,
    lat_deg: f64,
    alt_meter: f64,
) -> ObservationPoint {
    let el = el_deg.to_radians();
    let r_earth = calc_earth_radius(lat_deg);
    let r_eff = r_earth * 4_f64 / 3_f64;
    let sr = r_eff + alt_meter;
    let z = (r_meter * r_meter + sr * sr + r_meter * sr * 2_f64 * el.sin()).sqrt() - r_eff;
    let s = r_eff * ((r_meter * el.cos()) / (r_eff + z)).asin();
    ObservationPoint {
        alt_meter: z,
        dist_meter: s,
    }
}

pub(crate) const WGS84_RADIUS_EARTH_MAJOR: f64 = 6_378_137_f64;

pub(crate) const WGS84_INV_FLATTENING: f64 = 298.257_223_563;

pub(crate) fn calc_earth_radius(lat_deg: f64) -> f64 {
    let wgs84_radius_earth_minor =
        WGS84_RADIUS_EARTH_MAJOR * (1_f64 - 1_f64 / WGS84_INV_FLATTENING);
    let lat = lat_deg.to_radians();
    let cos_lat_2 = lat.cos().powf(2.into());
    let sin_lat_2 = lat.sin().powf(2.into());
    ((WGS84_RADIUS_EARTH_MAJOR.powf(4.into()) * cos_lat_2
        + wgs84_radius_earth_minor.powf(4.into()) * sin_lat_2)
        / (WGS84_RADIUS_EARTH_MAJOR.powf(2.into()) * cos_lat_2
            + wgs84_radius_earth_minor.powf(2.into()) * sin_lat_2))
        .sqrt()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    #[test]
    fn earth_radius_major() {
        assert_eq!(calc_earth_radius(0_f64), WGS84_RADIUS_EARTH_MAJOR);
    }

    #[test]
    fn earth_radius_minor() {
        assert_eq!(
            calc_earth_radius(90_f64),
            WGS84_RADIUS_EARTH_MAJOR * (1_f64 - 1_f64 / WGS84_INV_FLATTENING)
        );
    }
}
