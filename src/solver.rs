use crate::affine_expr::AffineExpression;
use crate::model::Model;

use crate::var::Variable;

use std::collections::HashMap;

pub trait SolveAlgorithm {
    fn solve(&self, mdl: &Model) -> Solution;
}

pub trait SolveSolution {
    fn var_value(&self, variable: &Variable) -> f64;

    fn expr_value(&self, expr: &AffineExpression) -> f64;
}

pub struct Solution {
    obj_fn_val: f64,
    var_map: HashMap<Variable, AffineExpression>,
    var_values: HashMap<Variable, f64>,
}

impl Solution {
    pub fn new(
        obj_fn_val: f64,
        var_map: HashMap<Variable, AffineExpression>,
        var_values: HashMap<Variable, f64>,
    ) -> Self {
        Self {
            obj_fn_val,
            var_map,
            var_values,
        }
    }

    pub fn var_value(&self, var: &Variable) -> f64 {
        if let Some(expr) = self.var_map.get(var) {
            expr.eval(&self.var_values)
        } else {
            self.var_values[var]
        }
    }
}

//impl Algorithm for Simplex {}

pub struct Solver<T: SolveAlgorithm> {
    mdl: Model,
    solve_algo: T,
}

impl<T: SolveAlgorithm> Solver<T> {
    fn new(mdl: Model, solve_algo: T) -> Self {
        Self {
            mdl,
            solve_algo,
        }
    }

    fn solve(&self) {
        self.solve_algo.solve(&self.mdl);
    }
}
