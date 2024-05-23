pub(crate) const WGS84_RADIUS_EARTH_MAJOR: f64 = 6_378_137_f64;

pub(crate) const WGS84_INV_FLATTENING: f64 = 298.257_223_563;

pub(crate) fn calc_earth_radius(lat_rad: f64) -> f64 {
    let wgs84_radius_earth_minor =
        WGS84_RADIUS_EARTH_MAJOR * (1_f64 - 1_f64 / WGS84_INV_FLATTENING);
    let cos_lat = lat_rad.cos();
    let cos_lat_2 = cos_lat * cos_lat;
    let sin_lat = lat_rad.sin();
    let sin_lat_2 = sin_lat * sin_lat;
    let wgs84_radius_earth_minor_2 = wgs84_radius_earth_minor * wgs84_radius_earth_minor;
    let wgs84_radius_earth_minor_4 = wgs84_radius_earth_minor_2 * wgs84_radius_earth_minor_2;
    let wgs84_radius_earth_major_2 = WGS84_RADIUS_EARTH_MAJOR * WGS84_RADIUS_EARTH_MAJOR;
    let wgs84_radius_earth_major_4 = wgs84_radius_earth_major_2 * wgs84_radius_earth_major_2;
    ((wgs84_radius_earth_major_4 * cos_lat_2 + wgs84_radius_earth_minor_4 * sin_lat_2)
        / (wgs84_radius_earth_major_2 * cos_lat_2 + wgs84_radius_earth_minor_2 * sin_lat_2))
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
            calc_earth_radius(std::f64::consts::PI / 2.0),
            WGS84_RADIUS_EARTH_MAJOR * (1_f64 - 1_f64 / WGS84_INV_FLATTENING)
        );
    }
}
