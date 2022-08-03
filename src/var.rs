use num::ToPrimitive;
use std::{
    cell::{Ref, RefCell, RefMut},
    hash::{Hash, Hasher},
    ops::{
        Add, AddAssign, Deref, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
    },
    rc::Rc,
};

use std::fmt;

use uuid::Uuid; //used for variables ID

use crate::affine_expr::AffineExpression;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VarType {
    Float,
    Int,
    Bool,
}

impl Default for VarType {
    fn default() -> Self {
        VarType::Float
    }
}

#[derive(Clone, Debug, Default)]
pub struct VariableDefinition {
    ty: VarType,
    lb: Option<f64>,
    ub: Option<f64>,
    name: String,
}

impl VariableDefinition {
    pub fn new(ty: VarType) -> Self {
        Self {
            ty,
            lb: None,
            ub: None,
            name: String::from(""),
        }
    }

    pub fn with_lb<T: ToPrimitive>(mut self, lb: T) -> Self {
        let _lb: Option<f64> = lb.to_f64();
        self.lb = _lb;
        assert!(self.valid_bounds());
        self
    }

    pub fn with_ub<T: ToPrimitive>(mut self, ub: T) -> Self {
        let _ub: Option<f64> = ub.to_f64();
        assert!(self.valid_bounds());
        self.ub = _ub;
        self
    }

    pub fn  with_name<T:ToString>(mut self, name: T) -> Self {
        let name = name.to_string();
        self.name = name;
        self
    }

    fn valid_bounds(&self) -> bool {
        //ensure if upper bound exists, lb <= ub
        if let Some(lb) = self.lb {
            if let Some(ub) = self.ub {
                return lb <= ub;
            }
        }
        true
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RawModelVariable {
    pub ty: VarType,
    pub lb: Option<f64>,
    pub ub: Option<f64>,
    pub name: String,
    id: Uuid,
}

impl Eq for RawModelVariable {}

impl Hash for RawModelVariable {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl RawModelVariable {
    pub(crate) fn new(var_def: VariableDefinition) -> Self {
        Self {
            ty: var_def.ty,
            lb: var_def.lb,
            ub: var_def.ub,
            name: var_def.name,
            id: Uuid::new_v4(),
        }
    }

    pub fn ty(&self) -> VarType {
        self.ty
    }
    pub fn lb(&self) -> Option<f64> {
        self.lb
    }
    pub fn ub(&self) -> Option<f64> {
        self.ub
    }
    pub fn name(&self) -> &str {
        self.name.as_str()
    }
    pub fn id(&self) -> Uuid {
        self.id
    }
    pub fn ty_mut(&mut self) -> &mut VarType {
        &mut self.ty
    }
    pub fn lb_mut(&mut self) -> &mut Option<f64> {
        &mut self.lb
    }
    pub fn ub_mut(&mut self) -> &mut Option<f64> {
        &mut self.ub
    }
    pub fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelVariable {
    var: Rc<RefCell<RawModelVariable>>,
}

impl Hash for ModelVariable {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id().hash(state);
    }
}

impl ModelVariable {
    pub(crate) fn new(var_def: VariableDefinition) -> Self {
        Self {
            var: Rc::new(RefCell::new(RawModelVariable::new(var_def))),
        }
    }

    pub(crate) fn ty(&self) -> VarType {
        self.var.borrow().ty()
    }

    pub(crate) fn lb(&self) -> Option<f64> {
        self.var.borrow().lb()
    }

    pub(crate) fn ub(&self) -> Option<f64> {
        self.var.borrow().ub()
    }

    pub(crate) fn name(&self) -> Ref<str> {
        //self.var.borrow().name()
        Ref::map(self.var.borrow(), |borrow| borrow.name())
    }

    pub(crate) fn id(&self) -> Uuid {
        self.var.borrow().id()
    }

    pub fn ty_mut(&mut self) -> RefMut<VarType> {
        RefMut::map(self.var.borrow_mut(), |borrow| borrow.ty_mut())
    }
    pub fn lb_mut(&mut self) -> RefMut<Option<f64>> {
        RefMut::map(self.var.borrow_mut(), |borrow| borrow.lb_mut())
    }
    pub fn ub_mut(&mut self) -> RefMut<Option<f64>> {
        RefMut::map(self.var.borrow_mut(), |borrow| borrow.ub_mut())
    }
    pub fn name_mut(&mut self) -> RefMut<String> {
        RefMut::map(self.var.borrow_mut(), |borrow| borrow.name_mut())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct RawVariableCollection {
    vars: Vec<ModelVariable>,
}

impl RawVariableCollection {
    pub(crate) fn new() -> Self {
        Self { vars: Vec::new() }
    }
    pub(crate) fn add_var(&mut self, var: ModelVariable) {
        self.vars.push(var);
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub(crate) struct VariableCollection {
    raw_collection: Rc<RefCell<RawVariableCollection>>,
}

impl VariableCollection {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn add_var(&mut self, var: ModelVariable) {
        self.raw_collection.borrow_mut().add_var(var)
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Environment {
    collection: VariableCollection,
    env_id: Uuid,
}

impl Environment {
    pub fn new() -> Self {
        Self {
            collection: VariableCollection::default(),
            env_id: Uuid::new_v4(),
        }
    }
    pub fn add_var(&mut self, var_def: VariableDefinition) -> Variable {
        let var = ModelVariable::new(var_def);
        self.collection.add_var(var.clone());
        Variable {
            mv: var,
            env: self.clone(),
        }
    }

    pub(crate) fn collection(&self) -> &VariableCollection {
        &self.collection
    }
}
//should this be an affine expression? -> easier to deal with unbounded variables
#[derive(Clone, Debug, PartialEq)]
pub struct Variable {
    pub(crate) mv: ModelVariable,
    pub(crate) env: Environment,
}

impl Eq for Variable {}

impl Hash for Variable {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.mv.id().hash(state);
    }
}

impl fmt::Display for Variable {
    
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

pub struct VariableTransformationInfo {
    pub expr: AffineExpression,
    pub added_vars: Vec<Variable>,
}

impl Variable {
    pub fn new(env: &mut Environment, variable_definition: VariableDefinition) -> Self {
        env.add_var(variable_definition)
    }

    pub fn ty(&self) -> VarType {
        self.mv.ty()
    }

    pub fn lb(&self) -> Option<f64> {
        self.mv.lb()
    }

    pub fn ub(&self) -> Option<f64> {
        self.mv.ub()
    }

    pub fn id(&self) -> Uuid {
        self.mv.id()
    }

    pub(crate) fn env(&self) -> Environment {
        self.env.clone()
    }

    pub fn name(&self) -> Ref<str> {
        self.mv.name()
    }

    pub fn ty_mut(&mut self) -> RefMut<VarType> {
        self.mv.ty_mut()
    }
    pub fn lb_mut(&mut self) -> RefMut<Option<f64>> {
        self.mv.lb_mut()
    }
    pub fn ub_mut(&mut self) -> RefMut<Option<f64>> {
        self.mv.ub_mut()
    }
    pub fn name_mut(&mut self) -> RefMut<String> {
        self.mv.name_mut()
    }

    pub(crate) fn is_standard(&self) -> bool {
        if let Some(lb) = self.lb() {
            if lb == 0.0_f64 {
                return true;
            }
        }
        false
    }

    pub fn as_standard_form(&self, bound_var: Option<&Variable>) -> Option<VariableTransformationInfo> {
        match self.lb() {
            //lower bound at zero -> do nothing
            Some(lb) if lb == 0.0_f64 => None,

            //bounded case -> add constant
            Some(lb) => {
                //add var
                let mut env = self.env.clone();
                let vd = VariableDefinition::new(VarType::Float).with_lb(0).with_name(format!("{}_o", self.name().to_string()));
                let var = Variable::new(&mut env, vd);
                let t_info = VariableTransformationInfo {
                    expr: &var + lb,
                    added_vars: vec![var.clone()],
                };
                Some(t_info)
            }
            //unbounded case -> sum of bounded variables
            None => {
                let mut env = self.env.clone();
                let vd1 = VariableDefinition::new(VarType::Float).with_lb(0).with_name(format!("{}_p", self.name().to_string()));
                let var1 = Variable::new(&mut env, vd1);
                let t_info = if let Some(bvar) = bound_var {
                    VariableTransformationInfo {
                        expr: &var1 + bvar,
                        added_vars: vec![var1.clone()],
                    }
                } else {
                    let vd2 = VariableDefinition::new(VarType::Float).with_lb(0).with_name(format!("{}_n", self.name().to_string()));
                    let var2 = Variable::new(&mut env, vd2);
                    VariableTransformationInfo {
                        expr: &var1 + &var2,
                        added_vars: vec![var1.clone(), var2.clone()],
                    }
                };
                Some(t_info)
            }
        }
    }
}
