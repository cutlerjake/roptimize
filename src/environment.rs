use crate::var::{Var, VarType};

use std::cell::RefCell;
use std::rc::Rc;

#[derive(Default, Debug)]
pub struct Env {
    vars: Vec<Rc<RefCell<Var>>>,
}

impl Env {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_var(&mut self, ty: VarType) -> Rc<RefCell<Var>> {
        let v = Rc::new(RefCell::new(Var::new(ty, None, None)));
        self.vars.push(v.clone());
        v
    }
}
