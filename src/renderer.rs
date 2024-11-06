use crate::{
    color_map::{ListedColorMap, RgbColor, RgbColorMap},
    Azimuth, PpiElevationSpecInDegrees, RadarObsCell, RangeGateSpecInMeter, VerticalCrossSection,
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
    el_spec: &'a PpiElevationSpecInDegrees,
    // M length
    az: Azimuth,
    vcs: &'a VerticalCrossSection,
    color_map: &'a ListedColorMap,
    default_color: RgbColor,
}

impl<'a> Render for FixedElevationScanVerticalRenderer<'a> {
    type Color = RgbColor;

    fn render(&self) -> impl Iterator<Item = Self::Color> {
        let cells = self.vcs.cells();
        cells.into_iter().map(|cell| {
            get_point_value(
                &cell,
                self.values,
                self.range_gate_spec,
                self.el_spec,
                &self.az,
            )
            .and_then(|value| self.color_map.get_rgb(value).cloned())
            .unwrap_or(self.default_color.clone())
        })
    }
}

fn get_point_value(
    cell: &RadarObsCell,
    values: &[f64],
    range_gate_spec: &RangeGateSpecInMeter,
    el_spec: &PpiElevationSpecInDegrees,
    az: &Azimuth,
) -> Option<f64> {
    if (cell.el_deg - el_spec.angle).abs() > el_spec.half_beam_width {
        return None;
    }

    let az_index = az.position(cell.az_deg)?;
    let r_index = range_gate_spec.find_index(cell.r_meter)?;
    let index = az_index * range_gate_spec.len() + r_index;
    Some(values[index])
}
