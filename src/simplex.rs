use ndarray::{s, Array2};

use crate::affine_expr::AffineExpression;
use crate::constraint::{Comp, Constraint};
use crate::model::Model;
use crate::solver::Solution;
use crate::tableau::{Tableau, TableauIx};
use crate::var::Variable;

use std::collections::HashMap;

pub struct Simplex {}

impl Simplex {
    pub fn new() -> Self {
        Self {}
    }

    // pub fn tbl(&self) -> &Tableau {
    //     &self.tbl
    // }

    // pub fn basic_vars(&self) -> &Vec<usize> {
    //     &self.basic_vars
    // }

    #[inline(always)]
    fn pivot_ind(&self, tbl: &Tableau) -> Option<TableauIx> {
        let j = tbl
            .tbl()
            .slice(s![-1, ..])
            .iter()
            .enumerate()
            .filter(|(_i, v)| **v < 0_f64)
            .min_by(|(_i1, v1), (_i2, v2)| v1.partial_cmp(v2).expect("Nan encountered"))
            .map(|(i, _v)| i)?;

        let i = tbl
            .tbl()
            .slice(s![.., j])
            .iter()
            .enumerate()
            .zip(tbl.tbl().slice(s![.., -1]))
            .filter(|((_i, a), b)| **a > 0_f64 && **b > 0_f64)
            .min_by(|((_i1, a1), b1), ((_i2, a2), b2)| {
                (*b1 / *a1)
                    .partial_cmp(&(*b2 / *a2))
                    .expect("Nan encountered")
            })
            .map(|((i, _a), _b)| i)?;

        Some(TableauIx::new(i, j))
    }

    pub fn solve(&mut self, mdl: &Model) -> Solution {
        println!("Original model:\n{}", &mdl);

        let std_mdl_info = mdl.as_standard_form(true);

        println!("Standardized model:\n{}", &std_mdl_info.mdl);

        let (mut tableau, mut basic_vars) = std_mdl_info.mdl.as_tableau();
        let mut pvt_cnt = 0;

        println!("Original tableau.");
        println!("{}", tableau.tbl());

        while let Some(ix) = self.pivot_ind(&tableau) {
            tableau.pivot(&ix);
            basic_vars[ix.i()] = ix.j();

            pvt_cnt += 1;
            println!("Pivot {}", pvt_cnt);
            println!("{}", tableau.tbl());
        }
        println!(
            "Objective function value: {}",
            tableau.tbl().slice(s![-1, -1])
        );

        let tableau_shape = tableau.tbl().shape();

        //get obj fn value
        let mut obj_fn_val = tableau.tbl()[[tableau_shape[0] - 1, tableau_shape[1] - 1]];
        if std_mdl_info.flipped_obj_fn {
            obj_fn_val *= -1.0_f64;
        }

        //create variable map
        let var_map = std_mdl_info
            .var_map
            .into_iter()
            .map(|(var, vt_info)| (var, vt_info.expr))
            .collect::<HashMap<Variable, AffineExpression>>();

        //get variable_values
        let var_index_map = std_mdl_info.mdl.variable_index_map();
        let index_var_map = var_index_map
            .iter()
            .map(|(var, usize)| (*usize, var.clone()))
            .collect::<HashMap<usize, Variable>>();
        let mut var_values = var_index_map
            .iter()
            .map(|(var, usize)| (var.clone(), 0.0_f64))
            .collect::<HashMap<Variable, f64>>();
        for (row, var_ind) in basic_vars.iter().enumerate() {
            let var_val = tableau.tbl()[[row, tableau_shape[1] - 1]];
            let entry = var_values
                .entry(index_var_map[var_ind].clone())
                .or_insert(0.0_f64);
            *entry - var_val;
        }

        Solution::new(obj_fn_val, var_map, var_values)
    }
}
