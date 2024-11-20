use std::{cell::LazyCell, str::FromStr};

use itertools::Itertools;

pub trait RgbColorMap {
    fn get_rgb(&self, value: f64) -> Option<&RgbColor>;
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

impl From<&[u8; 3]> for RgbColor {
    fn from(value: &[u8; 3]) -> Self {
        let [red, green, blue] = *value;
        Self { red, green, blue }
    }
}

impl From<&RgbColor> for [u8; 3] {
    fn from(value: &RgbColor) -> Self {
        [value.red, value.green, value.blue]
    }
}

impl FromStr for RgbColor {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        const MESSAGE: &str =
            r##"color must be in specified in hexadecimal numbers like "#ffffff""##;

        let mut chars = s.chars();
        if s.len() != 7 || chars.next().unwrap() != '#' {
            return Err(MESSAGE);
        }
        let chars: Option<Vec<_>> = chars.map(|c| c.to_digit(16)).collect();
        let mut chars = chars
            .ok_or_else(|| MESSAGE)?
            .into_iter()
            .map(|i| i as u8)
            .tuples()
            .map(|(a, b)| (a << 4) | b);
        Ok(Self::new(
            chars.next().unwrap(),
            chars.next().unwrap(),
            chars.next().unwrap(),
        ))
    }
}

// Thresholds are assumed to be sorted.
#[derive(Debug, PartialEq, Clone)]
pub struct ListedColorMap(pub Vec<(f64, RgbColor)>);

impl FromStr for ListedColorMap {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lines: Result<Vec<_>, _> = s
            .lines()
            .map(|line| {
                let (key, value) = line
                    .split_once(":")
                    .ok_or_else(|| "key/value must be separated by a colon")?;
                let key = key
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| "parsing a float value failed")?;
                let value = value.trim().parse::<RgbColor>()?;
                Ok((key, value))
            })
            .collect();
        Ok(Self(lines?))
    }
}

impl RgbColorMap for ListedColorMap {
    fn get_rgb(&self, value: f64) -> Option<&RgbColor> {
        if value.is_nan() {
            return None;
        }
        let Self(inner) = self;
        inner.iter().rev().find_map(|(threshold, color)| {
            if value > *threshold {
                Some(color)
            } else {
                None
            }
        })
    }
}

pub const JMA_RADAR_COLOR_MAP: LazyCell<ListedColorMap> = LazyCell::new(|| {
    ListedColorMap(vec![
        (-9999., RgbColor::from(&[0xff, 0xff, 0xff])),
        (0.01, RgbColor::from(&[0x99, 0xff, 0xff])),
        (1., RgbColor::from(&[0x66, 0xcc, 0xff])),
        (5., RgbColor::from(&[0x21, 0x98, 0xff])),
        (10., RgbColor::from(&[0x00, 0x38, 0xff])),
        (20., RgbColor::from(&[0xfa, 0xf5, 0x00])),
        (30., RgbColor::from(&[0xff, 0x99, 0x00])),
        (50., RgbColor::from(&[0xe7, 0x28, 0x00])),
        (80., RgbColor::from(&[0x9a, 0x00, 0x79])),
        (200., RgbColor::from(&[0x00, 0x00, 0x00])),
    ])
});

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! test_parsing_rgb_color {
        ($((
            $name:ident,
            $input:expr,
            $expected:expr
        ),)*) => ($(
            #[test]
            fn $name() {
                let s = $input;
                let actual = s.parse::<RgbColor>().ok();
                let expected = $expected;
                assert_eq!(actual, expected)
            }
        )*);
    }

    test_parsing_rgb_color! {
        (
            parsing_rgb_color_succeeds_with_numerical_chars,
            "#012345", Some(RgbColor::new(0x01, 0x23, 0x45))
        ),
        (
            parsing_rgb_color_succeeds_with_a_to_f_chars,
            "#abcdef", Some(RgbColor::new(0xab, 0xcd, 0xef))
        ),
        (parsing_rgb_color_fails, "foobar!", None),
        (parsing_rgb_color_fails_for_longer_text, "#ffffffff", None),
    }

    #[test]
    fn parsing_color_map() {
        let s = "0.0:#000000
            1.0 : #111111
            2.0 : #222222";
        let actual = s.parse::<ListedColorMap>();
        let expected = Ok(ListedColorMap(vec![
            (0.0, RgbColor::new(0x00, 0x00, 0x00)),
            (1.0, RgbColor::new(0x11, 0x11, 0x11)),
            (2.0, RgbColor::new(0x22, 0x22, 0x22)),
        ]));
        assert_eq!(actual, expected)
    }

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
        let color_map = ListedColorMap(colors);
        assert_eq!(color_map.get_rgb(-1.0), None);
        assert_eq!(color_map.get_rgb(-0.0), None);
        assert_eq!(color_map.get_rgb(0.0), None);
        assert_eq!(color_map.get_rgb(0.1), Some(&color_zero));
        assert_eq!(color_map.get_rgb(1.0), Some(&color_zero));
        assert_eq!(color_map.get_rgb(1.5), Some(&color_one));
        assert_eq!(color_map.get_rgb(2.0), Some(&color_one));
        assert_eq!(color_map.get_rgb(20.0), Some(&color_ten));
        assert_eq!(color_map.get_rgb(200.0), Some(&color_hundred));
    }
}
