use crate::{
    color_map::{ListedColorMap, RgbColor, RgbColorMap},
    Ppi, VerticalCrossSection,
};

pub trait Render {
    type Color;

    fn render(&self) -> impl Iterator<Item = Self::Color>;
}

pub struct FixedElevationScanVerticalRenderer<'a, 'v> {
    ppi: &'a Ppi<'v>,
    vcs: &'a VerticalCrossSection,
    color_map: &'a ListedColorMap,
    default_color: RgbColor,
}

impl<'a, 'v> FixedElevationScanVerticalRenderer<'a, 'v> {
    pub fn new(
        ppi: &'a Ppi<'v>,
        vcs: &'a VerticalCrossSection,
        color_map: &'a ListedColorMap,
        default_color: RgbColor,
    ) -> Self {
        Self {
            ppi,
            vcs,
            color_map,
            default_color,
        }
    }
}

impl<'a, 'v> Render for FixedElevationScanVerticalRenderer<'a, 'v> {
    type Color = RgbColor;

    fn render(&self) -> impl Iterator<Item = Self::Color> {
        let cells = self.vcs.cells();
        cells.into_iter().map(|cell| {
            self.ppi
                .value_at(&cell)
                .and_then(|value| self.color_map.get_rgb(value).cloned())
                .unwrap_or(self.default_color.clone())
        })
    }
}
