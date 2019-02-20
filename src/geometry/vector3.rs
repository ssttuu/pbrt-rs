use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use num::Num;
use num::ToPrimitive;

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vector3<T> where T: Num + Copy + ToPrimitive + PartialOrd {
    fn new(x: T, y: T, z: T) -> Self {
        Vector3 {
            x,
            y,
            z,
        }
    }

    fn cross(&self, other: Self) -> Self {
        Self {
            x: (self.y - other.z) - (self.z * other.y),
            y: (self.z * other.x) - (self.x * other.z),
            z: (self.x * other.y) - (self.y * other.x),
        }
    }
}

impl<T: Add<Output=T> + Copy> Add<Vector3<T>> for Vector3<T> {
    type Output = Vector3<T>;

    fn add(self, other: Vector3<T>) -> Vector3<T> {
        Vector3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Add<Output=T> + Copy> Add<T> for Vector3<T> {
    type Output = Vector3<T>;

    fn add(self, other: T) -> Vector3<T> {
        Vector3 {
            x: self.x + other,
            y: self.y + other,
            z: self.z + other,
        }
    }
}

impl<T: AddAssign + Copy> AddAssign<Self> for Vector3<T> {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<T: AddAssign + Copy> AddAssign<T> for Vector3<T> {
    fn add_assign(&mut self, other: T) {
        self.x += other;
        self.y += other;
        self.z += other;
    }
}

impl<T: Div<Output=T> + Copy> Div<T> for Vector3<T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<T: DivAssign + Copy> DivAssign<T> for Vector3<T> {
    fn div_assign(&mut self, other: T) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

impl<T> Index<i32> for Vector3<T> {
    type Output = T;

    fn index(&self, index: i32) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            1 => &self.z,
            _ => panic!("index out of range"),
        }
    }
}

impl<T> IndexMut<i32> for Vector3<T> {
    fn index_mut(&mut self, index: i32) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index out of range"),
        }
    }
}

impl<T: Mul<Output=T> + Copy> Mul<Self> for Vector3<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl<T: MulAssign + Copy> MulAssign<Self> for Vector3<T> {
    fn mul_assign(&mut self, other: Self) {
        self.x *= other.x;
        self.y *= other.y;
        self.z *= other.z;
    }
}


impl<T: MulAssign + Copy> MulAssign<T> for Vector3<T> {
    fn mul_assign(&mut self, other: T) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl<T: Sub<Output=T> + Copy> Sub<Self> for Vector3<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T: Sub<Output=T> + Copy> Sub<T> for Vector3<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self {
        Self {
            x: self.x - rhs,
            y: self.y - rhs,
            z: self.z - rhs,
        }
    }
}

impl<T: SubAssign + Copy> SubAssign<T> for Vector3<T> {
    fn sub_assign(&mut self, rhs: T) {
        self.x -= rhs;
        self.y -= rhs;
        self.z -= rhs;
    }
}
