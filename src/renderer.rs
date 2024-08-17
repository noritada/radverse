use crate::{
    color_map::{ColorMap, ListedColorMap, RgbColor},
    RadarObsCell, VerticalCrossSection,
};

pub trait Render {
    type Color;

    fn render(&self) -> impl Iterator<Item = Self::Color>;
}

pub struct FixedElevationScanVerticalRenderer<'a> {
    // N * M length
    values: &'a [f64],
    // N length
    r_start_meter: &'a [f64],
    // N length
    r_end_meter: &'a [f64],
    // M length
    az_deg: &'a [f64],
    // Fixed
    el_deg: f64,
    half_beam_width_deg: f64,
    vcs: &'a VerticalCrossSection,
    color_map: &'a ListedColorMap,
    default_color: RgbColor,
}

impl<'a> Render for FixedElevationScanVerticalRenderer<'a> {
    type Color = RgbColor;

    fn render(&self) -> impl Iterator<Item = Self::Color> {
        self.vcs.cells.iter().map(|cell| {
            get_point_value(
                cell,
                self.values,
                self.r_start_meter,
                self.r_end_meter,
                self.az_deg,
                self.el_deg,
                self.half_beam_width_deg,
            )
            .and_then(|value| self.color_map.get_color(value).cloned())
            .unwrap_or(self.default_color.clone())
        })
    }
}

fn get_point_value(
    cell: &RadarObsCell,
    values: &[f64],
    r_start_meter: &[f64],
    r_end_meter: &[f64],
    az_deg: &[f64],
    el_deg: f64,
    half_beam_width_deg: f64,
) -> Option<f64> {
    if (cell.el_deg - el_deg).abs() > half_beam_width_deg {
        return None;
    }

    let az_index = az_deg
        .iter()
        .position(|az| (cell.az_deg - az).abs() <= half_beam_width_deg)?;
    let r_index = r_start_meter
        .iter()
        .zip(r_end_meter.iter())
        .position(|(&start, &end)| cell.r_meter >= start && cell.r_meter <= end)?;
    let index = az_index * r_start_meter.len() + r_index;
    Some(values[index])
}
