use ndarray::{s, Array};
use sprs::CsVec;

use crate::affine_expr::AffineExpression;

use crate::model::Model;
use crate::solver::Solution;
use crate::tableau::{SparseTableau, Tableau, TableauIx, TblPrintInfo};
use crate::var::Variable;

use std::collections::HashMap;

pub struct Simplex {}

impl Simplex {
    pub fn new() -> Self {
        Self {}
    }

    #[inline(always)]
    fn pivot_ind(&self, tbl: &SparseTableau) -> Option<TableauIx> {
        let j = tbl
            .nb_r_costs()
            .filter(|(_i, v)| *v < 0_f64)
            .min_by(|(_i1, v1), (_i2, v2)| v1.partial_cmp(v2).expect("Nan encountered"))
            .map(|(i, _v)| i)?;

        let a_j = tbl.col(j);

        let i = tbl
            .rhs()
            .iter()
            .filter(|(i, &b)| b.is_sign_positive())
            .map(|(i, &b)| {
                let a = a_j.get(i).unwrap_or(&0.0);
                (i, b / a)
            })
            .filter(|(i, ratio)| ratio.is_sign_positive() && ratio.is_finite())
            .min_by(|(i1, ratio1), (i2, ratio2)| {
                ratio1.partial_cmp(&ratio2).expect("Encountered Nan")
            })
            .map(|(i, _)| i)?;

        Some(TableauIx::new(i, j))
    }

    fn _solve_with_print(&mut self, tbl: &mut SparseTableau) {
        let mut pvt_cnt = 0;

        while let Some(ix) = self.pivot_ind(tbl) {
            let table = Tableau::from(&*tbl).as_table(TblPrintInfo::Pivot(ix));

            println!("pivot: {}", pvt_cnt);
            println!("{}\n", table);
            tbl.pivot(ix);
            pvt_cnt += 1;

            if pvt_cnt > 5 {
                panic!()
            }
        }
    }

    fn _solve(&mut self, tbl: &mut SparseTableau) {
        // println!("Hello from _solve");
        while let Some(ix) = self.pivot_ind(&tbl) {
            tbl.pivot(ix);
        }
    }

    pub fn solve(&mut self, mdl: &Model, print: bool) -> Solution {
        let std_mdl_info = mdl.as_standard_form(true);

        //let mut tableau = std_mdl_info.mdl.as_tableau();
        let mut tableau = SparseTableau::from(&std_mdl_info.mdl);

        //solve first stage -> find inital basic feasible solution
        let mut second_stage_obj = tableau.obj_fn.clone();

        let (inds, dstorage): (Vec<usize>, Vec<f64>) = std_mdl_info
            .artificial_vars
            .iter()
            .map(|var| (tableau.var_map[var], 1.0))
            .unzip();

        let first_stage_obj =
            CsVec::new_from_unsorted(tableau.constraints.cols(), inds, dstorage).unwrap();

        tableau.obj_fn = first_stage_obj;

        //perform elementary row operations to set coeffs of artificial vars to 0
        std_mdl_info.artificial_vars.iter().for_each(|var| {
            //find col associated with artificial var
            let ix_j = tableau.var_map[var];

            //find first row with non-zero entry in column
            let ix_i = tableau
                .col(ix_j)
                .iter()
                .find(|(ind, &val)| val == 1.0)
                .map(|(i, _)| i)
                .unwrap();

            tableau.pivot(TableauIx::new(ix_i, ix_j));
        });

        match print {
            true => self._solve_with_print(&mut tableau),
            false => self._solve(&mut tableau),
        }

        //make sure all artificial variables 0
        if tableau.z() != 0.0 {
            let table = Tableau::from(&tableau).as_table(TblPrintInfo::Default);
            println!("{}", table);
            panic!("Infeasible")
        }

        if print {
            println!(
                "Final first stage tableau:\n{}\n",
                Tableau::from(&tableau).as_table(TblPrintInfo::Default)
            );
        }

        //set second stage obj
        tableau.obj_fn = second_stage_obj;

        //restore proper form through gaussian elimination
        let bvars = tableau.basic_vars.clone();

        bvars.iter().enumerate().for_each(|(ix_i, &col)| {
            tableau.pivot(TableauIx::new(ix_i, col));
        });

        if print {
            println!(
                "Second stage after fixing:\n{}\n",
                Tableau::from(&tableau).as_table(TblPrintInfo::Default)
            );
        }

        //remove artificial variables
        tableau.remove_vars(std_mdl_info.artificial_vars);

        if print {
            println!(
                "Second stage after removal:\n{}\n",
                Tableau::from(&tableau).as_table(TblPrintInfo::Default)
            );
        }

        //solve second stage problem
        match print {
            true => self._solve_with_print(&mut tableau),
            false => self._solve(&mut tableau),
        }

        if print {
            println!(
                "Final tableau:\n{}",
                Tableau::from(&tableau).as_table(TblPrintInfo::Default)
            );
        }

        let tableau_shape = tableau.constraints.shape();

        //get obj fn value
        let mut obj_fn_val = tableau.z();
        //let mut obj_fn_val = tableau.tbl()[[tableau_shape[0] - 1, tableau_shape[1] - 1]];
        let mut obj_constant = std_mdl_info.mdl.obj_fn.constant();
        if std_mdl_info.flipped_obj_fn {
            //obj_fn_val *= -1.0_f64;
            obj_constant *= -1.0_f64;
        }

        obj_fn_val += obj_constant;

        println!("Obj fn value: {}", obj_fn_val);

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
            .map(|(var, _ind)| (var.clone(), 0.0_f64))
            .collect::<HashMap<Variable, f64>>();
        // for (row, var_ind) in tableau.basic_vars.iter().enumerate() {
        //     let var_val = tableau.tbl()[[row, tableau_shape[1] - 1]];
        //     let entry = var_values
        //         .entry(index_var_map[var_ind].clone())
        //         .or_insert(0.0_f64);
        //     *entry = var_val;
        // }

        Solution::new(obj_fn_val, var_map, var_values)
    }
}
