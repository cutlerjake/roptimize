use ndarray::{s, Array};

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

        let rhs = tbl.rhs();
        let col_j = tbl.col(j);

        let a_j = tbl.col(j);

        let i = tbl
            .rhs()
            .iter()
            .filter(|(i, &b)| b > 0.0)
            .map(|(i, &b)| {
                let a = a_j.get(i).unwrap_or(&0.0);
                (i, b / a)
            })
            .filter(|(i, ratio)| *ratio > 0.0)
            .min_by(|(i1, ratio1), (i2, ratio2)| {
                ratio1.partial_cmp(&ratio2).expect("Encountered Nan")
            })
            // .min_by(|(i1, b1), (i2, b2)| {
            //     let a1 = a_j.get(*i1).unwrap_or(&0.0);
            //     let a2 = a_j.get(*i2).unwrap_or(&0.0);
            //     (*b1 / *a1)
            //         .partial_cmp(&(*b2 / *a2))
            //         .expect("Nan encountered")
            // })
            .map(|(i, _)| i)?;
        // let i = tbl
        //     .rhs()
        //     .iter()
        //     .zip(tbl.col(j))
        //     .enumerate()
        //     .filter(|(_, (&b, _))| b > 0.0)
        //     .min_by(|(_i1, (a1, b1)), (_i2, (a2, b2))| {
        //         (*b1 / *a1)
        //             .partial_cmp(&(*b2 / *a2))
        //             .expect("Nan encountered")
        //     })
        //     .map(|(i, (a, b))| i)?;

        println!("{}, {}", i, j);
        println!("{}", Tableau::from(tbl));

        // let j = tbl
        //     .tbl()
        //     .slice(s![-1, ..-1])
        //     .iter()
        //     .enumerate()
        //     .filter(|(_i, v)| **v < 0_f64)
        //     .min_by(|(_i1, v1), (_i2, v2)| v1.partial_cmp(v2).expect("Nan encountered"))
        //     .map(|(i, _v)| i)?;

        // let i = tbl
        //     .tbl()
        //     .slice(s![..-1, j])
        //     .iter()
        //     .enumerate()
        //     .zip(tbl.tbl().slice(s![.., -1]))
        //     .filter(|((_i, a), b)| (a.signum() == b.signum()) && (**a != 0.0_f64)) //**a > 0_f64 && **b > 0_f64)
        //     .min_by(|((_i1, a1), b1), ((_i2, a2), b2)| {
        //         (*b1 / *a1)
        //             .partial_cmp(&(*b2 / *a2))
        //             .expect("Nan encountered")
        //     })
        //     .map(|((i, _a), _b)| i)?;

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

        let mut tableau = std_mdl_info.mdl.as_tableau();

        //solve first stage -> find inital basic feasible solution
        let mut second_stage_obj = tableau.tbl().slice(s![-1, ..]).to_owned();
        let mut first_stage_obj_vec = vec![0.0_f64; tableau.tbl().shape()[1]];
        first_stage_obj_vec[tableau.tbl.shape()[1] - 2] = 1.0_f64;

        std_mdl_info.artificial_vars.iter().for_each(|var| {
            first_stage_obj_vec[tableau.vars[var]] = 1.0_f64;
        });

        //set first stage obj
        let mut first_stage_obj =
            Array::from_shape_vec(first_stage_obj_vec.len(), first_stage_obj_vec).unwrap();

        tableau.tbl.slice_mut(s![-1, ..]).assign(&first_stage_obj);

        //perform elementary row operations to set coeffs of artificial vars to 0
        std_mdl_info.artificial_vars.iter().for_each(|var| {
            //find col associated with artificial var
            let ix_j = tableau.vars[var];

            //find first row with non-zero entry in column
            let ix_i = tableau
                .tbl()
                .slice(s![..-1, ix_j])
                .iter()
                .position(|&v| v != 0.0_f64)
                .unwrap();
            if ix_i == tableau.tbl().shape()[0] {
                panic!("Encountered artificial var not included in constraints");
            }

            tableau.pivot(&TableauIx::new(ix_i, ix_j));

            //get obj and constraint row
            // let mut cons = tableau.tbl().slice(s![ix_i, ..]).to_owned();
            // let ratio = first_stage_obj[[ix_j]] / cons[[ix_j]];
            // first_stage_obj -= &(ratio * &cons);
        });

        //set first stage objectve
        //tableau.tbl.slice_mut(s![-1, ..]).assign(&first_stage_obj);
        //println!("DENSE:\n{}", tableau);

        let mut sp = SparseTableau::from(tableau.clone());
        //println!("SPARSE:\n{}", Tableau::from(&sp));

        match print {
            true => {
                println!("{}", tableau);
                self._solve_with_print(&mut sp)
            }
            false => self._solve(&mut sp),
        }
        tableau = Tableau::from(&sp);

        //make sure all artificial variables 0
        if tableau.tbl()[[tableau.tbl().shape()[0] - 1, tableau.tbl().shape()[1] - 1]] != 0.0_f64 {
            println!("{}", tableau);
            panic!("Infeasible")
        };

        if print {
            println!("Final first stage tableau:\n{}\n", tableau);
        }

        tableau.tbl.slice_mut(s![-1, ..]).assign(&second_stage_obj);

        //restore proper form through gaussian elimination
        let bvars = tableau.basic_vars.clone();

        bvars.iter().enumerate().for_each(|(ix_i, &col)| {
            tableau.pivot(&TableauIx::new(ix_i, col));
            //set coeff of obj fn to zero
            // let obj_coeff = second_stage_obj[[col]];
            // let con_coeff = tableau.tbl()[[ix_i, col]];
            // let ratio = obj_coeff / con_coeff;
            // //let ratio = con_coeff/obj_coeff;

            // let temp = &tableau.tbl().slice(s![ix_i, ..]) * ratio;

            // second_stage_obj -= &temp;
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

        let mut sp = SparseTableau::from(tableau.clone());

        match print {
            true => self._solve_with_print(&mut sp),
            false => self._solve(&mut sp),
        }

        tableau = Tableau::from(&sp);

        // match print {
        //     true => {
        //         println!("Basic feasible solution found. Starting solve.");
        //         self._solve_with_print(&mut SparseTableau::from(tableau));
        //     }
        //     false => self._solve(&mut SparseTableau::from(tableau)),
        // }

        if print {
            println!("Final tableau:\n{}", tableau);
        }

        let tableau_shape = tableau.tbl().shape();

        //get obj fn value
        let mut obj_fn_val = tableau.tbl()[[tableau_shape[0] - 1, tableau_shape[1] - 1]];
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
