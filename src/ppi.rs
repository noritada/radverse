#[derive(Debug, PartialEq)]
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
