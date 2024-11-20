use crate::RadarObsCell;

#[derive(Debug, PartialEq, Clone)]
pub struct Ppi {
    values: Vec<f64>,
    range: RangeGateSpecInMeter,
    az: Azimuth,
    pub(crate) el: PpiElevationSpecInDegrees,
}

impl Ppi {
    pub fn new(
        values: Vec<f64>,
        range: RangeGateSpecInMeter,
        az: Azimuth,
        el: PpiElevationSpecInDegrees,
    ) -> Self {
        Self {
            values,
            range,
            az,
            el,
        }
    }

    pub fn value_at(&self, cell: &RadarObsCell) -> Option<f64> {
        if (cell.el_deg - self.el.angle).abs() > self.el.half_beam_width {
            return None;
        }

        let az_index = self.az.position(cell.az_deg)?;
        let r_index = self.range.find_index(cell.r_meter)?;
        let index = az_index * self.range.len() + r_index;
        Some(self.values[index])
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct RangeGateSpecInMeter {
    start: f64,
    gate_spacing: f64,
    n_gates: usize,
}

impl RangeGateSpecInMeter {
    pub fn new(start: f64, gate_spacing: f64, n_gates: usize) -> Self {
        Self {
            start,
            gate_spacing,
            n_gates,
        }
    }

    pub fn find_index(&self, r_meter: f64) -> Option<usize> {
        if r_meter < self.start {
            return None;
        }
        let index = ((r_meter - self.start) / self.gate_spacing) as usize;
        if index > self.n_gates {
            None
        } else {
            Some(index)
        }
    }

    pub fn len(&self) -> usize {
        self.n_gates
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum Azimuth {
    Simple(usize),
    Degrees(Vec<f64>, f64),
}

impl Azimuth {
    pub fn position(&self, target_deg: f64) -> Option<usize> {
        if target_deg < 0. || target_deg >= 360. {
            return None;
        }
        let index = match self {
            Self::Simple(n_rays) => (target_deg / 360. * *n_rays as f64) as usize,
            Self::Degrees(degrees, half_beam_width) => degrees
                .iter()
                .position(|az| (target_deg - az).abs() <= *half_beam_width)?,
        };
        Some(index)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct PpiElevationSpecInDegrees {
    pub angle: f64,
    pub half_beam_width: f64,
}

impl PpiElevationSpecInDegrees {
    pub fn new(angle: f64, half_beam_width: f64) -> Self {
        Self {
            angle,
            half_beam_width,
        }
    }
}
