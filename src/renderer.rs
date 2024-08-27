use crate::{
    color_map::{ColorMap, ListedColorMap, RgbColor},
    PpiAngleSpecInDegrees, RadarObsCell, RangeGateSpecInMeter, VerticalCrossSection,
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
    angle_spec: &'a PpiAngleSpecInDegrees,
    // M length
    az_deg: &'a [f64],
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
                self.angle_spec,
                self.az_deg,
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
    angle_spec: &PpiAngleSpecInDegrees,
    az_deg: &[f64],
) -> Option<f64> {
    if (cell.el_deg - angle_spec.el).abs() > angle_spec.half_el_beam_width {
        return None;
    }

    let az_index = az_deg
        .iter()
        .position(|az| (cell.az_deg - az).abs() <= angle_spec.half_az_beam_width)?;
    let r_index = range_gate_spec.find_index(cell.r_meter)?;
    let index = az_index * range_gate_spec.len() + r_index;
    Some(values[index])
}
