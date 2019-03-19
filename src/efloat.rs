use crate::math;
use crate::math::consts;
use crate::math::Sqrt;
use crate::types::Float;
use std::ops;
use std::cmp;

#[derive(Clone, Copy, Debug, Default)]
pub struct EFloat {
    pub value: Float,
    pub low: Float,
    pub high: Float,
}

impl EFloat {
    pub fn new(value: Float, err: Float) -> EFloat {
        let ef = EFloat {
            value,
            low: if err != 0. {
                math::next_float_down(value - err)
            } else {
                0.
            },
            high: if err != 0. {
                math::next_float_up(value + err)
            } else {
                0.
            },
        };
        ef.check();
        ef
    }

    pub fn check(&self) {
        if self.low == consts::INFINITY
            || self.low == consts::NAN
            || self.high == consts::INFINITY
            || self.high == consts::NAN
        {
            panic!("EFloat lower bound is {:?}", self.low);
        }

        if self.low > self.high {
            panic!("EFloat low is greater than high: {:?}", self);
        }
    }

    pub fn quadratic(a: Self, b: Self, c: Self) -> Option<(Self, Self)> {
        let discriminant = b * b - (a * c * 4.);
        if discriminant.value < 0. {
            return None;
        }
        let root_discriminant = discriminant.sqrt();
        let float_root_discriminant =
            EFloat::new(root_discriminant.value, (root_discriminant * consts::EPSILON).value);

        let q = if b.value < 0. {
            (b - float_root_discriminant) * -0.5
        } else {
            (b + float_root_discriminant) * -0.5
        };

        let t0 = q / a;
        let t1 = c / q;

        if t0 <= t1 {
            Some((t0, t1))
        } else {
            Some((t1, t0))
        }
    }

    pub fn sqrt(&self) -> Self {
        let out = Self {
            value: self.value.sqrt(),
            low: math::next_float_down(self.low.sqrt()),
            high: math::next_float_up(self.high.sqrt()),
        };
        out.check();
        out
    }
}

impl cmp::PartialEq for EFloat {
    fn eq(&self, other: &Self) -> bool {
        self.value.eq(&other.value)
    }
}

impl cmp::PartialOrd for EFloat {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

impl ops::Add for EFloat {
    type Output = EFloat;
    fn add(self, rhs: Self) -> Self::Output {
        let ef = EFloat {
            value: self.value + rhs.value,
            low: math::next_float_down(self.low + rhs.low),
            high: math::next_float_up(self.high + rhs.high),
        };
        ef.check();
        ef
    }
}

impl ops::Add<Float> for EFloat {
    type Output = EFloat;
    fn add(self, rhs: Float) -> Self::Output {
        self + EFloat::new(rhs, 0.)
    }
}

impl ops::Div for EFloat {
    type Output = EFloat;
    fn div(self, rhs: Self) -> Self {
        let mut out = EFloat::default();
        out.value = self.value / rhs.value;
        if rhs.low < 0. && rhs.high > 0. {
            out.low = -consts::INFINITY;
            out.high = consts::INFINITY;
        } else {
            let div = [
                self.low / rhs.low,
                self.high / rhs.low,
                self.low / rhs.high,
                self.high / rhs.high,
            ];
            out.low = math::next_float_down(math::min(
                math::min(div[0], div[1]),
                math::min(div[2], div[3]),
            ));
            out.high = math::next_float_up(math::max(
                math::max(div[0], div[1]),
                math::max(div[2], div[3]),
            ));
        }
        out.check();
        out
    }
}

impl ops::Div<Float> for EFloat {
    type Output = EFloat;
    fn div(self, rhs: Float) -> Self::Output {
        self / EFloat::new(rhs, 0.)
    }
}

impl ops::Mul for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: Self) -> Self {
        let mut out = EFloat::default();
        out.value = self.value * rhs.value;
        let div = [
            self.low * rhs.low,
            self.high * rhs.low,
            self.low * rhs.high,
            self.high * rhs.high,
        ];
        out.low = math::next_float_down(math::min(
            math::min(div[0], div[1]),
            math::min(div[2], div[3]),
        ));
        out.high = math::next_float_up(math::max(
            math::max(div[0], div[1]),
            math::max(div[2], div[3]),
        ));
        out.check();
        out
    }
}

impl ops::Mul<Float> for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: Float) -> Self::Output {
        self * EFloat::new(rhs, 0.)
    }
}

impl ops::Sub for EFloat {
    type Output = EFloat;
    fn sub(self, rhs: Self) -> Self::Output {
        let ef = EFloat {
            value: self.value - rhs.value,
            low: math::next_float_down(self.low - rhs.low),
            high: math::next_float_up(self.high - rhs.high),
        };
        ef.check();
        ef
    }
}

impl ops::Sub<Float> for EFloat {
    type Output = EFloat;
    fn sub(self, rhs: Float) -> Self::Output {
        self - EFloat::new(rhs, 0.)
    }
}
