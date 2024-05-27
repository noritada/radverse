pub trait ColorMap {
    fn get_color(&self, value: f64) -> Option<&RgbColor>;
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RgbColor {
    pub red: u8,
    pub green: u8,
    pub blue: u8,
}

impl RgbColor {
    pub fn new(red: u8, green: u8, blue: u8) -> Self {
        Self { red, green, blue }
    }
}

pub(crate) struct ListedColorMap {
    pub(crate) name: String,
    // Colors are assumed to be sorted.
    colors: Vec<(f64, RgbColor)>,
}

impl ListedColorMap {
    pub(crate) fn new(name: String, colors: Vec<(f64, RgbColor)>) -> Self {
        Self { name, colors }
    }
}

impl ColorMap for ListedColorMap {
    fn get_color(&self, value: f64) -> Option<&RgbColor> {
        self.colors.iter().rev().find_map(|(threshold, color)| {
            if value > *threshold {
                Some(color)
            } else {
                None
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn listed_color_map() {
        let color_zero = RgbColor::new(0, 0, 0);
        let color_one = RgbColor::new(1, 1, 1);
        let color_ten = RgbColor::new(10, 10, 10);
        let color_hundred = RgbColor::new(100, 100, 100);
        let colors = vec![
            (0.0, color_zero.clone()),
            (1.0, color_one.clone()),
            (10.0, color_ten.clone()),
            (100.0, color_hundred.clone()),
        ];
        let color_map = ListedColorMap::new("foo".to_owned(), colors);
        assert_eq!(color_map.get_color(-1.0), None);
        assert_eq!(color_map.get_color(-0.0), None);
        assert_eq!(color_map.get_color(0.0), None);
        assert_eq!(color_map.get_color(0.1), Some(&color_zero));
        assert_eq!(color_map.get_color(1.0), Some(&color_zero));
        assert_eq!(color_map.get_color(1.5), Some(&color_one));
        assert_eq!(color_map.get_color(2.0), Some(&color_one));
        assert_eq!(color_map.get_color(20.0), Some(&color_ten));
        assert_eq!(color_map.get_color(200.0), Some(&color_hundred));
    }
}
