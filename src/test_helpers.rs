macro_rules! assert_almost_eq {
    ($a1:expr, $a2:expr, $d:expr) => {
        if $a1 - $a2 > $d {
            panic!(
                "assertion a1 - a2 <= delta failed\n a1 - a2: {} - {}\n   delta: {}",
                $a1, $a2, $d
            );
        } else if $a2 - $a1 > $d {
            panic!(
                "assertion a2 - a1 <= delta failed\n a2 - a1: {} - {}\n   delta: {}",
                $a2, $a1, $d
            );
        }
    };
}
pub(crate) use assert_almost_eq;
