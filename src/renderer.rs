use itertools::Itertools;

use crate::{
    color_map::{ListedColorMap, RgbColor, RgbColorMap},
    ElevationRanges, Ppi, RadarCenteredPoint, RadarObsCell, RadarObsCellVertical,
    VerticalCrossSection,
};

pub trait Render {
    type Color;

    fn render(&self) -> impl Iterator<Item = Self::Color>;
}

pub struct FixedElevationScanVerticalRenderer<'a, 'v> {
    obs: &'a [Ppi<'v>],
    vcs: &'a VerticalCrossSection,
    color_map: &'a ListedColorMap,
    default_color: RgbColor,
}

impl<'a, 'v> FixedElevationScanVerticalRenderer<'a, 'v> {
    pub fn new(
        obs: &'a [Ppi<'v>],
        vcs: &'a VerticalCrossSection,
        color_map: &'a ListedColorMap,
        default_color: RgbColor,
    ) -> Self {
        Self {
            obs,
            vcs,
            color_map,
            default_color,
        }
    }
}

impl<'a, 'v> Render for FixedElevationScanVerticalRenderer<'a, 'v> {
    type Color = RgbColor;

    fn render(&self) -> impl Iterator<Item = Self::Color> {
        let el_sorted = self
            .obs
            .iter()
            .map(|ppi| (ppi.el.clone(), ppi))
            .collect::<Vec<_>>();
        let el_ranges = ElevationRanges::from(&el_sorted);
        let h_items = self
            .vcs
            .horizontal_iter()
            .map(|point| {
                (
                    point,
                    el_ranges.z_ranges(point.site_distance, &self.vcs.site),
                )
            })
            .collect::<Vec<_>>();

        let h_iter = h_items.iter();
        let v_iter = self.vcs.vertical_iter();
        let pixels = v_iter
            .cartesian_product(h_iter)
            .map(|(alt_meter, (point, z_ranges))| {
                let value = z_ranges.find(alt_meter).and_then(|(_z, ppi)| {
                    ppi.and_then(|ppi| {
                        let cell = RadarCenteredPoint {
                            alt_meter: *alt_meter,
                            dist_meter: point.site_distance,
                        };
                        let cell = RadarObsCellVertical::from((&cell, &self.vcs.site));
                        let az_deg = point.site_direction.to_degrees();
                        let cell = RadarObsCell::new(cell.r_meter, az_deg, cell.el_deg);

                        ppi.value_at(&cell)
                    })
                });
                value
                    .and_then(|value| self.color_map.get_rgb(value).cloned())
                    .unwrap_or(self.default_color.clone())
            })
            .collect::<Vec<_>>();
        pixels.into_iter()
    }
}
