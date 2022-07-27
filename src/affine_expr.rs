use crate::var::{Environment, ModelVariable, Variable, VariableCollection};
use num::ToPrimitive;
use std::cell::{Ref, RefCell, RefMut};
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::{Add, AddAssign, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::rc::Rc;
use uuid::Uuid;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AffineExpression {
    pub(crate) coeffs: HashMap<Variable, f64>,
    pub(crate) constant: f64,
}

impl fmt::Display for AffineExpression {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let len = self.coeffs.len();
        let mult = "\u{00D7}";
        for (i, (var, coeff)) in self.coeffs.iter().enumerate() {
            if (i + 1 < len) | (self.constant != 0.0_f64) {
                write!(f, "{}{}{} + ", coeff, mult, var)?;
            } else {
                write!(f, "{}{}{}", coeff, mult, var)?;
            };
        }
        if (self.constant != 0.0_f64) | (len == 0) {
            write!(f, "{}", self.constant)?;
        }
        Ok(())
    }
}

impl From<&Variable> for AffineExpression {
    fn from(var: &Variable) -> Self {
        let mut coeffs = HashMap::new();
        coeffs.insert(var.clone(), 1.0_f64);
        let constant = 0.0_f64;
        let env = var.env().clone();

        Self {
            coeffs,
            constant,
        }
    }
}

impl <T: ToPrimitive> From<T> for AffineExpression {
    
    fn from(num: T) -> Self {
        let n = num.to_f64().unwrap();

        let mut expr = Self::default();
        expr.constant = n;
        expr
    }
}

impl AffineExpression {
    pub fn new(
        coeffs: HashMap<Variable, f64>,
        constant: f64,
    ) -> Self {
        Self {
            coeffs,
            constant,
        }
    }

    pub fn new_empty(env: Environment) -> Self {
        Self {
            coeffs: HashMap::new(),
            constant: 0.0_f64,
        }
    }

    pub(crate) fn clear(&mut self) {
        self.coeffs.clear();
        self.constant = 0.0_f64;
    }
    pub fn variables(&self) -> Vec<Variable> {
        self.coeffs.iter().map(|(k, _v)| k.clone()).collect()
    }

    pub fn constant(&self) -> f64 {
        self.constant
    }

    pub fn constant_mut(&mut self) -> &mut f64 {
        &mut self.constant
    }

    pub fn env(&self) -> Environment {
        self.variables()[0].env()
    }

    pub fn contains_var(&self, var: &Variable) -> bool {
        self.coeffs.contains_key(var)
    }

    pub fn replace_var(&mut self, var: Variable, expr: AffineExpression) -> bool {
        //check if variable in self
        //  replace with expr if exists

        if let Occupied(coeff) = self.coeffs.entry(var) {
            //remove variable from expression
            let (var, c) = coeff.remove_entry();
            //add expresiion
            let mut expr = expr;
            expr *= c;
            *self += expr;
            return true;
        }
        false
    }

    pub fn eval(&self, values: &HashMap<Variable, f64>) -> f64 {
        let mut val = 0.0;
        
        for (var, coeff) in self.coeffs.iter() {
            val += coeff * values[var];
        }
        val
    }
}

//AF + AF -> AF
impl Add for AffineExpression {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert!(self.env() == rhs.env());
        let mut exp = self;
        rhs.coeffs.iter().for_each(|(key, coeff)| {
            let c = exp.coeffs.entry(key.clone()).or_insert(0.0_f64);
            *c += coeff;
        });
        exp.coeffs.retain(|_, c| *c != 0.0_f64); 
        exp.constant += rhs.constant;
        exp
    }
}

//AF - AF -> AF
impl Sub for AffineExpression {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        assert!(self.env() == rhs.env());
        let mut exp = self;
        rhs.coeffs.iter().for_each(|(key, coeff)| {
            let c = exp.coeffs.entry(key.clone()).or_insert(0.0_f64);
            *c -= coeff;
        });
        exp.coeffs.retain(|_, c| *c != 0.0_f64); 
        exp.constant -= rhs.constant;
        exp
    }
}

//AF + V -> AF
impl Add<&Variable> for AffineExpression {
    type Output = Self;

    fn add(mut self, rhs: &Variable) -> Self::Output {
        //check var from same collection
        assert!(self.env() == rhs.env());

        let c = self.coeffs.entry(rhs.clone()).or_insert(0.0_f64);
        *c += 1.0_f64;

        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self
    }
}

//AF - V -> AF
impl Sub<&Variable> for AffineExpression {
    type Output = Self;

    fn sub(mut self, rhs: &Variable) -> Self::Output {
        //check var from same collection
        assert!(self.env() == rhs.env());

        let c = self.coeffs.entry(rhs.clone()).or_insert(0.0_f64);
        *c -= 1.0_f64;
        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self
    }
}

//V + AF -> AF
impl Add<AffineExpression> for &Variable {
    type Output = AffineExpression;

    fn add(self, rhs: AffineExpression) -> Self::Output {
        //check var from same collection
        assert!(self.env() == rhs.env());

        let mut _rhs = rhs;
        let c = _rhs.coeffs.entry(self.clone()).or_insert(0.0_f64);
        *c += 1.0_f64;
        _rhs.coeffs.retain(|_, c| *c != 0.0_f64); 
        _rhs
    }
}

//V - AF -> AF
impl Sub<AffineExpression> for &Variable {
    type Output = AffineExpression;

    fn sub(self, rhs: AffineExpression) -> Self::Output {
        //check var from same collection
        assert!(self.env() == rhs.env());

        let mut _rhs = -rhs;
        let c = _rhs.coeffs.entry(self.clone()).or_insert(0.0_f64);
        *c += 1.0_f64;

        _rhs.coeffs.retain(|_, c| *c != 0.0_f64); 
        _rhs
    }
}

//V + V -> AF
impl Add for &Variable {
    type Output = AffineExpression;

    fn add(self, rhs: Self) -> Self::Output {
        //check vars from same collection
        assert!(self.env() == rhs.env());

        let mut coeffs = HashMap::new();

        for var in [self, rhs] {
            let c = coeffs.entry(var.clone()).or_insert(0.0_f64);
            *c += 1.0_f64;

        }

        coeffs.retain(|_, c| *c != 0.0_f64); 
        Self::Output {
            coeffs,
            constant: 0.0_f64,
        }
    }
}

//V - V -> AF
impl Sub for &Variable {
    type Output = AffineExpression;

    fn sub(self, rhs: Self) -> Self::Output {
        //check vars from same collection
        assert!(self.env == rhs.env);

        let mut coeffs = HashMap::new();

        //add lhs
        coeffs.insert(self.clone(), 1.0_f64);

        //subtract rhs
        let c = coeffs.entry(rhs.clone()).or_insert(0.0_f64);
        *c -= 1.0_f64;

        coeffs.retain(|_, c| *c != 0.0_f64); 
        Self::Output {
            coeffs,
            constant: 0.0_f64,
        }
    }
}

//C + AF -> AF
impl<T: ToPrimitive> Add<T> for AffineExpression {
    type Output = AffineExpression;

    fn add(self, rhs: T) -> Self::Output {
        let mut lhs = self;
        lhs.constant += rhs.to_f64().unwrap();
        lhs
    }
}
//AF + C -> AF
macro_rules! var_left_scalar_add_af_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Add<AffineExpression> for $T {
            type Output = AffineExpression;

            fn add(self, rhs: AffineExpression) -> Self::Output {
                let mut rhs = rhs;
                rhs.constant += self.to_f64().unwrap();
                rhs
            }
        }
    )*}
);

var_left_scalar_add_af_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

//AF - C -> AF
impl<T: ToPrimitive> Sub<T> for AffineExpression {
    type Output = AffineExpression;

    fn sub(self, rhs: T) -> Self::Output {
        let mut lhs = self;

        lhs.constant -= rhs.to_f64().unwrap();
        lhs
    }
}

//C - AF -> AF
macro_rules! var_left_scalar_sub_af_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Sub<AffineExpression> for $T {
            type Output = AffineExpression;

            fn sub(self, rhs: AffineExpression) -> Self::Output {
                //negate expr
                let mut rhs = -rhs;
                rhs.constant += self.to_f64().unwrap();
                rhs
            }
        }
    )*}
);

var_left_scalar_sub_af_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

//AF * C -> AF
impl<T: ToPrimitive> Mul<T> for AffineExpression {
    type Output = AffineExpression;

    fn mul(self, rhs: T) -> Self::Output {
        let mut lhs = self;
        lhs.coeffs.iter_mut().for_each(|(_,v)| *v *= rhs.to_f64().unwrap());
        lhs.constant *= rhs.to_f64().unwrap();
        lhs
    }
}
//C * AF -> AF
macro_rules! var_left_scalar_mul_af_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Mul<AffineExpression> for $T {
            type Output = AffineExpression;

            fn mul(self, rhs: AffineExpression) -> Self::Output {
                let mut rhs = rhs;
                rhs.coeffs.iter_mut().for_each(|(_, v)| *v *= self.to_f64().unwrap());
                rhs.constant *= self.to_f64().unwrap();
                rhs
            }
        }
    )*}
);

var_left_scalar_mul_af_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

//C + V -> AF
impl<T: ToPrimitive> Add<T> for &Variable {
    type Output = AffineExpression;

    fn add(self, rhs: T) -> Self::Output {
        let mut lhs = AffineExpression::from(self);
        lhs.constant += rhs.to_f64().unwrap();
        lhs
    }
}
//V + C -> AF
macro_rules! var_left_scalar_add_variable_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Add<&Variable> for $T {
            type Output = AffineExpression;

            fn add(self, rhs: &Variable) -> Self::Output {
                let mut rhs = AffineExpression::from(rhs);
                rhs.constant += self.to_f64().unwrap();
                rhs
            }
        }
    )*}
);

var_left_scalar_add_variable_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

//V * C -> AF
impl<T: ToPrimitive> Mul<T> for &Variable {
    type Output = AffineExpression;

    fn mul(self, rhs: T) -> Self::Output {
        let mut lhs = AffineExpression::from(self);
        lhs* rhs
    }
}
//C * V -> AF
macro_rules! var_left_scalar_mul_v_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Mul<&Variable> for $T {
            type Output = AffineExpression;

            fn mul(self, rhs: &Variable) -> Self::Output {
                let mut rhs = AffineExpression::from(rhs);
                rhs * self.to_f64().unwrap()
            }
        }
    )*}
);

var_left_scalar_mul_v_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);
//V - C -> AF
impl<T: ToPrimitive> Sub<T> for &Variable {
    type Output = AffineExpression;

    fn sub(self, rhs: T) -> Self::Output {
        let mut lhs = AffineExpression::from(self);
        lhs.constant -= rhs.to_f64().unwrap();
        lhs
    }
}
//C - V -> AF
macro_rules! var_left_scalar_sub_variable_impl(
    ($($T: ty), *$(, )*) => {$(
        impl Sub<&Variable> for $T {
            type Output = AffineExpression;

            fn sub(self, rhs: &Variable) -> Self::Output {
                let mut rhs = -AffineExpression::from(rhs);
                rhs.constant += self.to_f64().unwrap();
                rhs
            }
        }
    )*}
);

var_left_scalar_sub_variable_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

//AF += AF
impl AddAssign for AffineExpression {
    fn add_assign(&mut self, rhs: Self) {
        rhs.coeffs.into_iter().for_each(|(rhs_k, rhs_v)| {
            self.coeffs
                .entry(rhs_k)
                .and_modify(|lhs_v| *lhs_v += rhs_v)
                .or_insert(rhs_v);
        });
        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self.constant += rhs.constant;
    }
}
//AF -= AF
impl SubAssign for AffineExpression {
    fn sub_assign(&mut self, rhs: Self) {
        rhs.coeffs.into_iter().for_each(|(rhs_k, rhs_v)| {
            self.coeffs
                .entry(rhs_k)
                .and_modify(|lhs_v| *lhs_v -= rhs_v)
                .or_insert(-rhs_v);
        });

        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self.constant -= rhs.constant;
    }
}
//AF += V
impl AddAssign<&Variable> for AffineExpression {
    fn add_assign(&mut self, rhs: &Variable) {
        let mut _rhs = AffineExpression::from(rhs);
        _rhs.coeffs.into_iter().for_each(|(rhs_k, rhs_v)| {
            self.coeffs
                .entry(rhs_k)
                .and_modify(|lhs_v| *lhs_v += rhs_v)
                .or_insert(rhs_v);
        });

        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self.constant += _rhs.constant;
    }
}
//AF -= V
impl SubAssign<&Variable> for AffineExpression {
    fn sub_assign(&mut self, rhs: &Variable) {
        let mut _rhs = AffineExpression::from(rhs);
        _rhs.coeffs.into_iter().for_each(|(rhs_k, rhs_v)| {
            self.coeffs
                .entry(rhs_k)
                .and_modify(|lhs_v| *lhs_v -= rhs_v)
                .or_insert(-rhs_v);
        });

        self.coeffs.retain(|_, c| *c != 0.0_f64); 
        self.constant -= _rhs.constant;
    }
}
//AF += C
impl<T: ToPrimitive> AddAssign<T> for AffineExpression {
    fn add_assign(&mut self, rhs: T) {
        self.constant += rhs.to_f64().unwrap();
    }
}

//AF -= C
impl<T: ToPrimitive> SubAssign<T> for AffineExpression {
    fn sub_assign(&mut self, rhs: T) {
        self.constant -= rhs.to_f64().unwrap();
    }
}

//AF *= C
impl<T: ToPrimitive> MulAssign<T> for AffineExpression {
    fn mul_assign(&mut self, rhs: T) {
        self.coeffs
            .iter_mut()
            .for_each(|(key, val)| *val *= rhs.to_f64().unwrap());
        self.constant *= rhs.to_f64().unwrap();
    }
}

//AF /= C
impl<T: ToPrimitive> DivAssign<T> for AffineExpression {
    fn div_assign(&mut self, rhs: T) {
        self.coeffs
            .iter_mut()
            .for_each(|(key, val)| *val /= rhs.to_f64().unwrap());
        self.constant /= rhs.to_f64().unwrap();
    }
}

//-AF
impl Neg for AffineExpression {
    type Output = AffineExpression;

    fn neg(self) -> Self::Output {
        let mut af = self.clone();
        af *= -1;
        af
    }
}

//-V
impl Neg for Variable {
    type Output = AffineExpression;

    fn neg(self) -> Self::Output {
        let mut af = AffineExpression::from(&self);
        -af
    }
}
impl Neg for &Variable {
    type Output = AffineExpression;

    fn neg(self) -> Self::Output {
        let mut af = AffineExpression::from(self);
        -af
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var::{VariableDefinition, VarType};

    #[test]
    fn af_add_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        let mut af2 = AffineExpression::from(&vb);
        
        *af1.constant_mut() = 1.0_f64;
        *af2.constant_mut() = 2.0_f64;
        //(a + 1) + (b+2)
        let af3 = af1 + af2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 3.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af3 == af3_comp);

    }
    #[test]
    fn af_sub_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        let mut af2 = AffineExpression::from(&vb);
        
        *af1.constant_mut() = 1.0_f64;
        *af2.constant_mut() = 2.0_f64;
        //(a + 1) - (b+2)
        let af3 = af1 - af2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = -1.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af3 == af3_comp);

    }

    #[test]
    fn af_add_v() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)+b
        let af3 = af1 + &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 1.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af3 == af3_comp);
    }

    #[test]
    fn af_sub_v() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)-b
        let af3 = af1 - &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = 1.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af3 == af3_comp);
    }

    #[test]
    fn v_add_af() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)+b
        let af3 = &vb + af1;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 1.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af3 == af3_comp);
    }

    #[test]
    fn v_sub_af() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)-b
        let af3 = &vb - af1;

        let coeffs = HashMap::from([(va.clone(), -1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = -1.0_f64;
        let af3_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af3 == af3_comp);
    }

    #[test]
    fn v_add_v() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

       
        //a+b
        let af1 = &va + &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 0.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af1 == af1_comp);
    }

    #[test]
    fn v_sub_v() {
        
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

       
        //(a+1)-b
        let af1 = &va - &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = 0.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af1 == af1_comp);
    }

    #[test]
    fn af_add_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)+2
        let af2 = af1 + 2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = 3.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn c_add_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)+2
        let af2 = 2 + af1;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = 3.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn af_sub_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)-2
        let af2 = af1 - 2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = -1.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)-2 = a-1
        println!("{}", af2);
        println!("{}", af2_comp);
        assert!(af2 == af2_comp);
    }

    #[test]
    fn c_sub_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //2-(a+1)
        let af2 = 2 - af1;

        let coeffs = HashMap::from([(va.clone(), -1.0_f64)]);
        let constant = 1.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn v_add_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

       
        //a+1
        let af1 = &va + 1;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //a+1= a+1
        assert!(af1 == af1_comp);
    }

    #[test]
    fn c_add_v() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
       
        //1+a
        let af1 = 1 + &va;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //1+a= a+1
        assert!(af1 == af1_comp);
    }

    #[test]
    fn v_sub_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        //a-1
        let af1 = &va - 1;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = -1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //a-1 = a-1
        assert!(af1 == af1_comp);
    }

    #[test]
    fn c_sub_v() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
       
        //1-a
        let af1 = 1 - &va;

        let coeffs = HashMap::from([(va.clone(), -1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //1-a = 1-a
        assert!(af1 == af1_comp);
    }

    #[test]
    fn af_addassign_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        let mut af2 = AffineExpression::from(&vb);
        
        *af1.constant_mut() = 1.0_f64;
        *af2.constant_mut() = 2.0_f64;
        //(a + 1) + (b+2)
        af1 += af2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 3.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_subassign_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        let mut af2 = AffineExpression::from(&vb);
        
        *af1.constant_mut() = 1.0_f64;
        *af2.constant_mut() = 2.0_f64;
        //(a + 1) - (b+2)
        af1 -= af2;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = -1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_addassign_v() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        
        *af1.constant_mut() = 1.0_f64;
        //(a + 1) + b
        af1 += &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), 1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_subassign_v() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        
        *af1.constant_mut() = 1.0_f64;
        //(a + 1) - b
        af1 -= &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+(b+2) = a+b+3
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_addassign_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        
        *af1.constant_mut() = 1.0_f64;
        //(a + 1) + 1
        af1 += 1;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64)]);
        let constant = 2.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+1 = a+2
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_subassign_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());
        let vdb = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);
        let vb = Variable::new(&mut env, vdb);

        let mut af1 = AffineExpression::from(&va);
        
        *af1.constant_mut() = 1.0_f64;
        //(a + 1) - b
        af1 -= &vb;

        let coeffs = HashMap::from([(va.clone(), 1.0_f64), (vb.clone(), -1.0_f64)]);
        let constant = 1.0_f64;
        let af1_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)-b = a - b + 1
        assert!(af1 == af1_comp);

    }

    #[test]
    fn af_mul_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //(a+1)-2
        let af2 = af1 * 2;

        let coeffs = HashMap::from([(va.clone(), 2.0_f64)]);
        let constant = 2.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)-2 = a-1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn c_mul_af() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        let mut af1 = AffineExpression::from(&va);
        *af1.constant_mut() = 1.0_f64;
       
        //2-(a+1)
        let af2 = 2 * af1;

        let coeffs = HashMap::from([(va.clone(), 2.0_f64)]);
        let constant = 2.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn v_mul_c() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

       
        //(a+1)-2
        let af2 = &va * 2;

        let coeffs = HashMap::from([(va.clone(), 2.0_f64)]);
        let constant = 0.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)-2 = a-1
        assert!(af2 == af2_comp);
    }

    #[test]
    fn c_mul_v() {
        let vda = VariableDefinition::new(VarType::Float).with_lb(0).with_name("a".to_string());

        let mut env = Environment::new();

        let va = Variable::new(&mut env, vda);

        //2-(a+1)
        let af2 = 2 * &va;

        let coeffs = HashMap::from([(va.clone(), 2.0_f64)]);
        let constant = 0.0_f64;
        let af2_comp = AffineExpression::new(coeffs, constant, env.clone());
        //(a+1)+b= a+b+1
        assert!(af2 == af2_comp);
    }
}
