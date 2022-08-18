use std::collections::HashSet;
use std::fmt;

use crate::affine_expr::AffineExpression;
use crate::var::{VarType, Variable, VariableDefinition};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comp {
    Le, // <=
    Ge, // >=
    Eq, // ==
}

impl fmt::Display for Comp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Comp::Le => write!(f, "\u{2264}"),
            Comp::Eq => write!(f, "="),
            Comp::Ge => write!(f, "\u{2265}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Constraint {
    pub lhs: AffineExpression,
    pub comp: Comp,
    pub rhs: AffineExpression,
}

impl fmt::Display for Constraint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} {} {}", self.lhs, self.comp, self.rhs)
    }
}

pub struct ConstraintTransformationInfo {
    pub cons: Constraint,
    pub slack_vars: Vec<Variable>,
    pub artificial_vars: Vec<Variable>,
}

impl Constraint {
    pub fn new<T: Into<AffineExpression>, U: Into<AffineExpression>>(
        lhs: T,
        comp: Comp,
        rhs: U,
    ) -> Self {
        Self {
            lhs: lhs.into(),
            comp,
            rhs: rhs.into(),
        }
    }

    pub fn lhs(&self) -> &AffineExpression {
        &self.lhs
    }

    pub fn rhs(&self) -> &AffineExpression {
        &self.rhs
    }

    pub fn comp(&self) -> Comp {
        self.comp
    }

    pub fn variables(&self) -> Vec<Variable> {
        let lhs_vars = self.lhs.variables();
        let rhs_vars = self.rhs.variables();

        let vars: HashSet<Variable> = lhs_vars.into_iter().chain(rhs_vars.into_iter()).collect();
        vars.into_iter().collect()
    }

    pub fn as_standard_form(&self, suffix: Option<String>) -> ConstraintTransformationInfo {
        let mut rhs = self.rhs().clone();
        let mut lhs = self.lhs().clone();
        let mut comp = self.comp();
        let mut env = self.variables()[0].env();

        let suffix = match suffix {
            Some(s) => s,
            None => "".to_string(),
        };
        //move all variables to lhs
        lhs -= rhs.clone(); // todo: update when subassin by ref implemented
        rhs.clear();

        //move constant to right side
        rhs -= lhs.constant();
        *lhs.constant_mut() = 0.0_f64;

        //ensure positive rhs
        if rhs.constant() < 0.0_f64 {
            lhs *= -1;
            rhs *= -1;
            match comp {
                Comp::Le => comp = Comp::Ge,
                Comp::Eq => {}
                Comp::Ge => comp = Comp::Le,
            }
        }

        //convert to equality, add slack/artificial variables as required
        let mut slack_vars = Vec::new();
        let mut artificial_vars = Vec::new();
        match comp {
            Comp::Le => {
                //add slack variable to lhs
                let vd = VariableDefinition::new(VarType::Float)
                    .with_lb(0.0)
                    .with_name("S".to_string() + suffix.as_str());
                let svar = Variable::new(&mut env, vd);
                slack_vars.push(svar.clone());
                lhs += &svar;
            }
            Comp::Eq => {
                //add artificial variable to rhs
                let vd = VariableDefinition::new(VarType::Float)
                    .with_lb(0.0)
                    .with_name("A".to_string() + suffix.as_str());
                let avar = Variable::new(&mut env, vd);
                artificial_vars.push(avar.clone());
                lhs += &avar;
            }
            Comp::Ge => {
                //add artificial and slack variable to rhs
                let vd = VariableDefinition::new(VarType::Float)
                    .with_lb(0.0)
                    .with_name("A".to_string() + suffix.as_str());
                let avar = Variable::new(&mut env, vd);
                artificial_vars.push(avar.clone());
                //add slack varibale (surplus)
                let vd = VariableDefinition::new(VarType::Float)
                    .with_lb(0.0)
                    .with_name("S".to_string() + suffix.as_str());
                let svar = Variable::new(&mut env, vd);
                slack_vars.push(svar.clone());
                lhs += &avar;
                lhs -= &svar;
            }
        }

        let cons = Constraint::new(lhs, Comp::Eq, rhs);

        ConstraintTransformationInfo {
            cons,
            slack_vars,
            artificial_vars,
        }
    }

    pub fn contains_var(&self, var: &Variable) -> bool {
        self.lhs.contains_var(var) || self.rhs.contains_var(var)
    }
    pub fn replace_var(&mut self, var: &Variable, expr: &AffineExpression) {
        self.lhs.replace_var(var.clone(), expr.clone());
        self.rhs.replace_var(var.clone(), expr.clone());
    }
}
