use ndarray::{s, Array, Array2};

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
        println!("{}", tbl.tbl);
        let j = tbl
            .tbl()
            .slice(s![-1, ..-1])
            .iter()
            .enumerate()
            .filter(|(_i, v)| **v < 0_f64)
            .min_by(|(_i1, v1), (_i2, v2)| v1.partial_cmp(v2).expect("Nan encountered"))
            .map(|(i, _v)| i)?;

        println!("J: {}", j);
        let i = tbl
            .tbl()
            .slice(s![..-1, j])
            .iter()
            .enumerate()
            .zip(tbl.tbl().slice(s![.., -1]))
            .filter(|((_i, a), b)| (a.signum() == b.signum()) && (**a != 0.0_f64)) //**a > 0_f64 && **b > 0_f64)
            .min_by(|((_i1, a1), b1), ((_i2, a2), b2)| {
                (*b1 / *a1)
                    .partial_cmp(&(*b2 / *a2))
                    .expect("Nan encountered")
            })
            .map(|((i, _a), _b)| i)?;

        println!("i: {}", i);
        Some(TableauIx::new(i, j))
    }

    fn _solve(&mut self, tbl: &mut Tableau ) {
        // println!("Hello from _solve");
        let mut pvt_cnt = 0;
        while let Some(ix) = self.pivot_ind(&tbl) {
            tbl.pivot(&ix);
            // println!("pivot: {}", pvt_cnt);
            pvt_cnt += 1;
        }
    }

    pub fn solve(&mut self, mdl: &Model) -> Solution {

        let std_mdl_info = mdl.as_standard_form(true);

        let mut tableau = std_mdl_info.mdl.as_tableau();

        //solve first stage -> find inital basic feasible solution
        let mut second_stage_obj = tableau.tbl().slice(s![-1, ..]).to_owned();
        let mut first_stage_obj_vec = vec![0.0_f64; tableau.tbl().shape()[1]];

        std_mdl_info.artificial_vars.iter().for_each(|var| {
            first_stage_obj_vec[tableau.vars[var]] = 1.0_f64;
        });

        //set first stage obj
        let mut first_stage_obj =
            Array::from_shape_vec((first_stage_obj_vec.len()), first_stage_obj_vec).unwrap();

        //perform elementary row operations to set coeffs of artificial vars to 0
        std_mdl_info.artificial_vars.iter().for_each(|var| {
            //find col associated with artificial var
            let ix_j = tableau.vars[var];

            //find first row with non-zero entry in column
            let ix_i = tableau
                .tbl()
                .slice(s![..-1, ix_j])
                .iter()
                .position(|&v| v == 1.0_f64)
                .unwrap();
            if ix_i == tableau.tbl().shape()[0] {
                panic!("Encountered artificial var not included in constraints");
            }

            //get obj and constraint row
            let cons = tableau.tbl().slice(s![ix_i, ..]);

            first_stage_obj -= &cons;
        });

        //set first stage objectve
        tableau.tbl.slice_mut(s![-1, ..]).assign(&first_stage_obj);

        self._solve(&mut tableau);

        //make sure all artificial variables 0
        if tableau.tbl()[[tableau.tbl().shape()[0] - 1, tableau.tbl().shape()[1] - 1]] != 0.0_f64 {
            panic!("Infeasible")
        };

        //restore proper form through gaussian elimination
        tableau.basic_vars.iter().for_each(|&col| {
            //col associated with basic var -> only one non-zero value
            let ix_i = tableau
                .tbl()
                .slice(s![..-1, col])
                .iter()
                .position(|&v| v == 1.0_f64)
                .unwrap();

            //set coeff of obj fn to zero
            let obj_coeff = tableau.tbl()[[tableau.tbl().shape()[0] - 1, col]];
            let con_coeff = tableau.tbl()[[ix_i, col]];
            let ratio = obj_coeff / con_coeff;

            let temp = &tableau.tbl().slice(s![ix_i, ..]) * ratio;

            second_stage_obj -= &temp;
        });

        //set second stage obj
        tableau.tbl.slice_mut(s![-1, ..]).assign(&second_stage_obj);

        //remove artificial variables
        let mut col_mask = vec![true; tableau.tbl.shape()[1]];
        std_mdl_info.artificial_vars.iter().for_each(|var| {
            col_mask[tableau.vars[var]] = false;
        });

        tableau.filter_cols(col_mask);

        //solve second stage problem
        self._solve(&mut tableau);

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
            .map(|(var, ind)| (*ind, var.clone()))
            .collect::<HashMap<usize, Variable>>();
        let mut var_values = var_index_map
            .iter()
            .map(|(var, ind)| (var.clone(), 0.0_f64))
            .collect::<HashMap<Variable, f64>>();
        for (row, var_ind) in tableau.basic_vars.iter().enumerate() {
            let var_val = tableau.tbl()[[row, tableau_shape[1] - 1]];
            let entry = var_values
                .entry(index_var_map[var_ind].clone())
                .or_insert(0.0_f64);
            *entry = var_val;
        }

        Solution::new(obj_fn_val, var_map, var_values)
    }
}
