use crate::{
    color_map::{ColorMap, ListedColorMap, RgbColor},
    RadarObsCell, RangeGateSpecInMeter, VerticalCrossSection,
};

pub trait Render {
    type Color;

    fn render(&self) -> impl Iterator<Item = Self::Color>;
}

pub struct FixedElevationScanVerticalRenderer<'a> {
    // N * M length
    values: &'a [f64],
    // n = N
    range_gate_spec: &'a RangeGateSpecInMeter,
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
                self.range_gate_spec,
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
    range_gate_spec: &RangeGateSpecInMeter,
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
    let r_index = range_gate_spec.find_index(cell.r_meter)?;
    let index = az_index * range_gate_spec.len() + r_index;
    Some(values[index])
}
