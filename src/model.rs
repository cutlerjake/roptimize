use colored::*;
use ndarray::{s, Array2};
use tabular::{Row, Table};
use uuid::Uuid; //used for unique variable ID

use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::affine_expr::AffineExpression;
use crate::constraint::{Comp, Constraint};
use crate::tableau::Tableau;
use crate::var::{Environment, VarType, Variable, VariableDefinition, VariableTransformationInfo};

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum OptDir {
    Max,
    Min,
}

impl fmt::Display for OptDir {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            OptDir::Max => write!(f, "Max"),
            OptDir::Min => write!(f, "Min"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    pub(crate) obj_fn: AffineExpression,
    pub(crate) opt_dir: OptDir,
    pub(crate) constraints: Vec<Constraint>,
    pub(crate) env: Environment,
}

pub struct ModelTransformationInfo {
    pub mdl: Model,
    pub var_map: HashMap<Variable, VariableTransformationInfo>,
    pub artificial_vars: Vec<Variable>,
    pub flipped_obj_fn: bool,
}

impl Model {
    pub fn new(env: Environment) -> Self {
        Self {
            obj_fn: AffineExpression::default(),
            constraints: Vec::new(),
            opt_dir: OptDir::Min,
            env,
        }
    }

    //set objective function and optimization direction
    pub fn set_obj_fn(&mut self, opt_dir: OptDir, obj_fn: AffineExpression) {
        self.obj_fn = obj_fn;
        self.opt_dir = opt_dir;
    }

    //set constraints
    pub fn set_constraints(&mut self, constraints: Vec<Constraint>) {
        self.constraints = constraints;
    }

    //add a constrain to model
    pub fn add_constraint(&mut self, constraint: Constraint) {
        self.constraints.push(constraint);
    }

    //retrieve variables from the model
    pub fn variables(&self) -> Vec<Variable> {
        let mut vars: HashSet<Variable> =
            self.obj_fn.variables().iter().map(|v| v.clone()).collect();
        self.constraints.iter().for_each(|c| {
            c.variables().into_iter().for_each(|v| {
                vars.insert(v);
            })
        });
        vars.into_iter().collect()
    }

    pub fn variable_index_map(&self) -> HashMap<Variable, usize> {
        let vars = self.variables();
        vars.into_iter()
            .enumerate()
            .map(|(i, var)| (var, i))
            .collect::<HashMap<Variable, usize>>()
    }

    pub fn as_standard_form(&self, bound_var: bool) -> ModelTransformationInfo {
        //mut convert all variables to non-negative
        //all constraints must be in standard form (positive constant rhs, initial basic feasible
        //solution possible)
        //step 1 -> convert all variables to standard form
        //   how to track mapping?
        //   how to handle unbounded variables -> pass var and add constraints?
        let mut mdl = self.clone();
        let bvar = if bound_var {
            let vd = VariableDefinition::new(VarType::Float)
                .with_lb(0)
                .with_name("BV".to_string());
            let var = Variable::new(&mut mdl.env, vd);
            Some(var)
        } else {
            None
        };
        let mut var_map: HashMap<Variable, VariableTransformationInfo> = HashMap::new();
        for ref var in mdl.variables() {
            //standardize and add to var map
            if let Some(info) = var.as_standard_form(bvar.as_ref()) {
                // mdl.var_map.insert(var.clone(), info);
                var_map.insert(var.clone(), info);
                // if let Some(ref bvar) = bvar {
                //     // build constraint
                //     let lhs = AffineExpression::from(var);
                //     // let rhs = mdl.var_map[var].expr.clone();
                //     let rhs = var_map[var].expr.clone();

                //     //Current approach: N constraint -> consider converting constraint to sum
                //     //require bvar greater than all other variables
                //     let cons = Constraint::new(lhs, Comp::Ge, rhs);

                //     //add constraint to model
                //     mdl.add_constraint(cons);
                // }
            } if let Some(ub) = var.ub() {
                let cons = Constraint::new(var, Comp::Le, ub);
                mdl.add_constraint(cons);
            }
        }
        //step 2 update all constraints to replace non-standard variables
        for constraint in &mut mdl.constraints {
            for (var, info) in &var_map {
                constraint.replace_var(var, &info.expr);
            }
        }

        //step 3 update obj-fn to replace non-standard variables
        for (var, info) in &var_map {
            mdl.obj_fn.replace_var(var.clone(), info.expr.clone());
        }

        //step 4 convert all constraints to standard form
        let mut artificial_vars = Vec::new();
        let mut slack_vars = Vec::new();
        for (i, constraint) in mdl.constraints.iter_mut().enumerate() {
            let info = constraint.as_standard_form(Some("_".to_string() + i.to_string().as_str()));
            *constraint = info.cons;
            artificial_vars.extend(info.artificial_vars);
            slack_vars.extend(info.slack_vars);
        }

        //convert to minimization problem
        let mut flipped_obj_fn = false;
        match mdl.opt_dir {
            OptDir::Max => {
                //set z = -obj_expr
                //min z
                mdl.obj_fn *= -1;
                mdl.opt_dir = OptDir::Min;
                flipped_obj_fn = true;
            }
            OptDir::Min => {
                //Do nothing
            }
        }
        ModelTransformationInfo {
            mdl,
            var_map,
            artificial_vars,
            flipped_obj_fn,
        }
    }

    pub fn as_tableau(&self) -> Tableau {
        //
        //map all variables to unique index
        let var_ind_map = self.variable_index_map();
        //
        //constraint rows followed by obj_fn row
        let mut tbl = Array2::<f64>::zeros((self.constraints.len() + 1, var_ind_map.len() + 2));

        //populate constraint rows
        for (i, con) in self.constraints.iter().enumerate() {
            //ensure constraints in standard form
            assert!(
                con.rhs().coeffs.len() == 0,
                "Variable(s) on rhs of constraint"
            );
            assert!(con.comp() == Comp::Eq, "Non equality constraint");
            assert!(
                con.lhs().constant == 0.0_f64,
                "Constant on lhs of constraint"
            );

            //pupulate lhs
            for (var, coeff) in &con.lhs().coeffs {
                tbl[[i, var_ind_map[var]]] = *coeff;
            }

            //populate rhs
            tbl[[i, var_ind_map.len() + 1]] = con.rhs().constant();
        }

        //populate obj_fn
        let mut obj_fn = vec![0.0; var_ind_map.len() + 2];
        //obj_fn[var_ind_map.len()] = 1.0_f64;
        tbl[[self.constraints.len(), var_ind_map.len()]] = 1.0_f64;
        for (var, coeff) in self.obj_fn.coeffs.iter() {
            //obj_fn[var_ind_map[var]] = *coeff;
            tbl[[self.constraints.len(), var_ind_map[var]]] = *coeff;
        }

        //find basic vars
        let mut basic_vars = vec![0; self.constraints.len()];
        for col_idx in 0..var_ind_map.len() {
            let col = tbl.slice(s![..-1, col_idx]);

            let mut basic = true;
            let mut found_one = false;
            let mut ind = 0;
            for (i, &val) in col.iter().enumerate() {
                if (val == 1.0_f64) {
                    if found_one {
                        basic = false;
                        break;
                    }
                    found_one = true;
                    ind = i;
                } else if (val != 0.0_f64) {
                    basic = false;
                    break;
                }
            }

            if basic {
                basic_vars[ind] = col_idx;
                //basic_vars.push(col_idx);
            }
        }

        //create tableau
        let mut tableau = Tableau::new(tbl, var_ind_map, basic_vars);

        tableau
    }

    pub fn tbl_variable_columns_string(&self, n_vars: usize) -> String {
        "{:>}".repeat(n_vars)
    }

    pub fn table_column_string(&self, n_vars: usize) -> String {
        let mut columns = "{:<}{:^}".to_string(); //OptDir, colon
        columns += &self.tbl_variable_columns_string(2 * n_vars + 2); //lhs var, constant
        columns += "{:^}"; //Comp
        columns += &self.tbl_variable_columns_string(2 * n_vars + 2); //rhs var, constant
        columns
    }

}

impl fmt::Display for Model {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //build column spec:
        //  OptDir: ObjFn
        //  Subject To:
        //      C1
        //      C2
        //      ...
        //OptDir, colon, LHSVars, Comp, RHSVars
        let plus = "+".to_string();
        let minus = "-".to_string();
        let vars = self.variables();

        let columns = self.table_column_string(vars.len());
        let mut table = Table::new(columns.as_str());

        let mut row = Row::new();
        // row.add_cell(self.opt_dir);
        // row.add_cell(":");

        // for _ in 0..vars.len() + 2 {
        //     row.add_cell("");
        // }

        let var_ind_map = self.variable_index_map();
        let mut row_vec = vec!["".to_string(); 4 * vars.len() + 7];
        row_vec[0] = format!("{}", self.opt_dir);
        row_vec[1] = ":".to_string();
        for (var, coeff) in &self.obj_fn.coeffs {
            let ind = var_ind_map[var];
            row_vec[2 * ind + 2] = if *coeff >= 0.0 {
                plus.clone()
            } else {
                minus.clone()
            };
            row_vec[2 * ind + 3] = format!("{}*{}", coeff.abs(), var.name());
        }
        if (self.obj_fn.constant != 0.0_f64) | (self.obj_fn.coeffs.len() == 0) {
            row_vec[2 + 2 * vars.len()] = if self.obj_fn.constant >= 0.0 {
                plus.clone()
            } else {
                minus.clone()
            };
            row_vec[3 + 2 * vars.len()] = format!("{}", self.obj_fn.constant.abs());
        }

        //format objfn expression to remove uneccesary first "+" if present
        for s in row_vec[2..2 * vars.len() + 4].iter_mut() {
            if (*s == plus) {
                *s = "".to_string();
                break;
            }
            if (s != "") {
                break;
            }
        }

        row_vec.iter().for_each(|cell| {
            row.add_cell(cell);
        });

        table.add_row(row);

        let mut row_vec2 = vec!["".to_string(); 4 * vars.len() + 7];
        row_vec2[0] = "Subject to".to_string();
        row_vec2[1] = ":".to_string();
        table.add_row(Row::from_cells(row_vec2));

        for constraint in &self.constraints {
            let mut cons_vec = vec!["".to_string(); 4 * vars.len() + 7];
            // lhs
            for (var, coeff) in &constraint.lhs().coeffs {
                let ind = var_ind_map[var];
                cons_vec[2 * ind + 2] = if *coeff >= 0.0 {
                    plus.clone()
                } else {
                    minus.clone()
                };
                cons_vec[2 * ind + 3] = format!("{}*{}", coeff.abs(), var.name());
            }
            if (constraint.lhs().constant != 0.0_f64) | (constraint.lhs().coeffs.len() == 0) {
                cons_vec[2 + 2 * vars.len()] = if constraint.lhs().constant >= 0.0 {
                    plus.clone()
                } else {
                    minus.clone()
                };
                cons_vec[3 + 2 * vars.len()] = format!("{}", constraint.lhs().constant.abs());
            }

            //format lhs expression to remove uneccesary first "+" if present
            for s in cons_vec[2..2 * vars.len() + 4].iter_mut() {
                if (*s == plus) {
                    *s = "".to_string();
                    break;
                }
                if (s != "") {
                    break;
                }
            }

            // comparison
            cons_vec[2 * vars.len() + 4] = format!(" {} ", constraint.comp());

            // rhs
            for (var, coeff) in &constraint.rhs().coeffs {
                let ind = var_ind_map[var];
                cons_vec[2 * ind + 5 + 2 * vars.len()] = if *coeff >= 0.0 {
                    plus.clone()
                } else {
                    minus.clone()
                };
                cons_vec[2 * ind + 6 + 2 * vars.len()] = format!("{}*{}", coeff.abs(), var.name());
            }
            if (constraint.rhs().constant != 0.0_f64) | (constraint.rhs().coeffs.len() == 0) {
                cons_vec[5 + 4 * vars.len()] = if constraint.rhs().constant >= 0.0 {
                    plus.clone()
                } else {
                    minus.clone()
                };
                cons_vec[6 + 4 * vars.len()] = format!("{}", constraint.rhs().constant.abs());
            }

            //format rhs expression to remove uneccesary first "+" if present
            for s in cons_vec[2 * vars.len() + 5..].iter_mut() {
                if (*s == plus) {
                    *s = "".to_string();
                    break;
                }
                if (s != "") {
                    break;
                }
            }

            table.add_row(Row::from_cells(cons_vec));
        }

        write!(f, "{}", table)
    }
}
