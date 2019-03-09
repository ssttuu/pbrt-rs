use crate::types::Float;
use crate::types::Number;
use std::cmp;
use std::ops;

const N_SPECTRAL_SAMPLES: usize = 60;

pub enum SpectrumType {
    Sampled,
    RGB,
}

pub static DEFAULT_SPECTRUM_TYPE: SpectrumType = SpectrumType::Sampled;

#[derive(Debug, Clone)]
pub struct Spectrum {
    data: Vec<Float>,
}

impl Spectrum {
    fn new_with(v: Float, size: usize) -> Self {
        Self {
            data: vec![v; size],
        }
    }

    pub fn new(v: Float) -> Self {
        match DEFAULT_SPECTRUM_TYPE {
            SpectrumType::Sampled => Self {
                data: vec![v; N_SPECTRAL_SAMPLES],
            },
            SpectrumType::RGB => Self { data: vec![v; 3] },
        }
    }

    pub fn set_all(&mut self, v: Float) {
        for c in &mut self.data {
            *c = v
        }
    }

    pub fn max(&self) -> Float {
        self.data.iter().cloned().fold(0. / 0.0, Float::max)
    }

    pub fn clamp(&mut self, low: Float, high: Float) {
        for i in 0..self.data.len() {
            self[i] = num::clamp(self[i], low, high);
        }
    }

    pub fn is_black(&self) -> bool {
        for i in 0..self.data.len() {
            if self[i] != 0.0 {
                return false;
            }
        }
        true
    }

    pub fn has_nans(&self) -> bool {
        for i in 0..self.data.len() {
            if self[i].is_nan() {
                return true;
            }
        }
        false
    }
}

impl cmp::PartialEq for Spectrum {
    fn eq(&self, rhs: &Self) -> bool {
        self.data == rhs.data
    }
}

impl ops::Index<usize> for Spectrum {
    type Output = Float;

    fn index(&self, index: usize) -> &Float {
        &self.data[index]
    }
}

impl ops::IndexMut<usize> for Spectrum {
    fn index_mut(&mut self, index: usize) -> &mut Float {
        &mut self.data[index]
    }
}

impl ops::Add for Spectrum {
    type Output = Spectrum;

    fn add(self, rhs: Self) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] += rhs[i];
        }
        out
    }
}

impl ops::Add<Float> for Spectrum {
    type Output = Spectrum;

    fn add(self, rhs: Float) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] += rhs;
        }
        out
    }
}

impl ops::AddAssign for Spectrum {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.data.len() {
            self[i] += rhs[i];
        }
    }
}

impl ops::Div for Spectrum {
    type Output = Spectrum;

    fn div(self, rhs: Self) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] /= rhs[i];
        }
        out
    }
}

impl ops::Div<Float> for Spectrum {
    type Output = Spectrum;

    fn div(self, rhs: Float) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] /= rhs;
        }
        out
    }
}

impl ops::DivAssign for Spectrum {
    fn div_assign(&mut self, rhs: Self) {
        for i in 0..self.data.len() {
            self[i] /= rhs[i];
        }
    }
}

impl ops::Mul for Spectrum {
    type Output = Spectrum;

    fn mul(self, rhs: Self) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] *= rhs[i];
        }
        out
    }
}

impl ops::Mul<Float> for Spectrum {
    type Output = Spectrum;

    fn mul(self, rhs: Float) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] *= rhs;
        }
        out
    }
}

impl ops::MulAssign for Spectrum {
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..self.data.len() {
            self[i] *= rhs[i];
        }
    }
}

impl ops::Sub for Spectrum {
    type Output = Spectrum;

    fn sub(self, rhs: Self) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] -= rhs[i];
        }
        out
    }
}

impl ops::Sub<Float> for Spectrum {
    type Output = Spectrum;

    fn sub(self, rhs: Float) -> Self {
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out[i] -= rhs;
        }
        out
    }
}

impl ops::SubAssign for Spectrum {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..self.data.len() {
            self[i] -= rhs[i];
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_deref() {
        let mut r = Spectrum::new(1.0);
        r.set_all(2.0);

        assert_eq!(r, Spectrum::new(2.0));
    }
}
