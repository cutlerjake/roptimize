use crate::model::Model;
use crate::var::{self, Variable};

use ndarray::linalg::Dot;
use ndarray::{s, Array, Array1, Array2, ArrayView, Axis, Slice, Zip};

use sprs::prod::{csc_mulacc_dense_colmaj, csc_mulacc_dense_rowmaj, mul_acc_mat_vec_csc};
use sprs::smmp::mul_csr_csr;
use sprs::{CompressedStorage, CsMat, CsVec};

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::{cmp::Ordering, iter::Enumerate};

use rustc_hash::FxHashSet;

use tabled::{
    builder::Builder,
    object::{Cell, Columns, Object, Rows},
    style::Border,
    Modify, Style, Table,
};

#[derive(Copy, Clone, Debug, Hash)]
pub enum TblPrintInfo {
    Default,
    Pivot(TableauIx),
}

impl TblPrintInfo {
    pub fn get_padding(&self) -> [usize; 4] {
        match self {
            TblPrintInfo::Default => [0, 0, 0, 0],
            TblPrintInfo::Pivot(_) => [1, 0, 1, 1],
        }
    }
}

#[derive(Copy, Clone, Debug, Hash)]
pub struct TableauIx {
    i: usize,
    j: usize,
}

impl TableauIx {
    pub fn new(i: usize, j: usize) -> Self {
        Self { i, j }
    }

    pub fn i(&self) -> usize {
        self.i
    }

    pub fn j(&self) -> usize {
        self.j
    }

    pub fn i_mut(&mut self) -> &mut usize {
        &mut self.i
    }

    fn j_mut(&mut self) -> &mut usize {
        &mut self.j
    }
}

#[derive(Copy, Clone, Debug)]
struct CooElem {
    pub row: usize,
    pub col: usize,
    pub data: f64,
}

#[derive(Debug, Clone)]
pub struct SparseTableau {
    pub constraints: CsMat<f64>,
    pub obj_fn: CsVec<f64>,
    pub basis: CsMat<f64>,
    pub basic_vars: Vec<usize>,
    pub non_basic_vars: FxHashSet<usize>,
    pub var_map: HashMap<Variable, usize>,
}

impl From<&Model> for SparseTableau {
    fn from(mdl: &Model) -> Self {
        Self::new(mdl)
    }
}

impl From<Tableau> for SparseTableau {
    fn from(tbl: Tableau) -> Self {
        // pub struct Tableau {
        //     pub(crate) tbl: Array2<f64>,
        //     pub(crate) vars: HashMap<Variable, usize>,
        //     pub(crate) basic_vars: Vec<usize>,
        // }

        //extract data from tableau
        let mut constraints = CsMat::csc_from_dense(tbl.tbl.slice(s![..-1, ..-2]), 1e-6);
        constraints = constraints.append_outer(tbl.tbl.slice(s![..-1, -1]).to_vec().as_slice());
        let var_map = tbl.vars;
        let basic_vars = tbl.basic_vars;

        //create obj_fn
        let (inds, dstorage): (Vec<usize>, Vec<f64>) = tbl
            .tbl
            .slice(s![-1, ..])
            .iter()
            .enumerate()
            .filter(|(i, val)| val.abs() > 1e-6 && *i != tbl.tbl.shape()[1] - 2)
            .map(|(i, val)| {
                if i < tbl.tbl.shape()[1] - 2 {
                    (i, val)
                } else {
                    (i - 1, val)
                }
            })
            .unzip();

        let obj_fn = CsVec::new(constraints.cols(), inds, dstorage);

        //construct basis
        let mut basis_dense =
            Array2::<f64>::from_elem((constraints.rows(), constraints.rows()), 0.0);

        basic_vars.iter().enumerate().for_each(|(i, &bvar_ind)| {
            basis_dense
                .slice_mut(s![.., i])
                .assign(&tbl.tbl.slice(s![..-1, bvar_ind]));
        });
        let basis = CsMat::csc_from_dense(basis_dense.view(), 1e-6);

        //get non basic inds
        let non_basic_vars = (0..constraints.cols() - 1)
            .into_iter()
            .filter(|i| !basic_vars.contains(i))
            .collect::<FxHashSet<usize>>();

        Self {
            constraints,
            obj_fn,
            basis,
            basic_vars,
            non_basic_vars,
            var_map,
        }
    }
}

impl SparseTableau {
    pub fn constraints(&self) -> &CsMat<f64> {
        &self.constraints
    }
    pub fn obj_fn(&self) -> &CsVec<f64> {
        &self.obj_fn
    }
    pub fn basis(&self) -> &CsMat<f64> {
        &self.basis
    }
    pub fn basic_vars(&self) -> &Vec<usize> {
        &self.basic_vars
    }
    pub fn non_basic_vars(&self) -> &FxHashSet<usize> {
        &self.non_basic_vars
    }
    pub fn var_map(&self) -> &HashMap<Variable, usize> {
        &self.var_map
    }

    pub fn new(mdl: &Model) -> Self {
        let var_index_map = mdl.variable_index_map();
        let nrows = mdl.constraints.len();
        let ncols = var_index_map.len() + 1;

        //create constraint matrix in COO format
        let mut coo = Vec::new();

        for (row, constraint) in mdl.constraints.iter().enumerate() {
            //populate lhs of constraint
            for (var, coeff) in &constraint.lhs.coeffs {
                if *coeff != 0.0 {
                    coo.push(CooElem {
                        row,
                        col: var_index_map[var],
                        data: *coeff,
                    });
                }
            }
            //populate rhs of constrain
            if constraint.rhs.constant != 0.0 {
                coo.push(CooElem {
                    row,
                    col: var_index_map.len(),
                    data: constraint.rhs.constant,
                })
            }
        }

        //sort by column, then row
        coo.sort_unstable_by_key(|elem| (elem.col, elem.row));

        //create column indices (cumulative nnz by beginning of each col)
        let (indptr, _) =
            coo.iter()
                .enumerate()
                .fold((vec![0; 2], 0), |(mut indptr, mut cur_col), (i, elem)| {
                    let mut back = indptr.len() - 1;

                    while cur_col < elem.col {
                        indptr.push(indptr[back]);
                        back += 1;
                        cur_col += 1;
                    }
                    indptr[back] += 1;
                    (indptr, cur_col)
                });

        //create row inds and data storage
        let (indices, dstorage): (Vec<usize>, Vec<f64>) =
            coo.into_iter().map(|elem| (elem.row, elem.data)).unzip();

        //create constraint matrix in CSC format
        let constraints =
            CsMat::new_from_unsorted_csc((nrows, ncols), indptr, indices, dstorage).unwrap();

        //create obj vec
        let (inds, data): (Vec<usize>, Vec<f64>) = mdl
            .obj_fn
            .coeffs
            .iter()
            .map(|(var, coeff)| (var_index_map[var], coeff))
            .unzip();

        let obj_fn = CsVec::new_from_unsorted(ncols, inds, data).unwrap();

        //Create basis
        let basis = CsMat::eye_csc(nrows);
        let mut basic_vars = vec![None; nrows];
        let mut non_basic_vars = FxHashSet::default();
        for (i, col) in constraints
            .slice_outer(..constraints.cols() - 1)
            .outer_iterator()
            .enumerate()
        {
            if col.data().len() == 1 && col.data()[0] == 1.0 && basic_vars[col.indices()[0]] == None
            {
                basic_vars[col.indices()[0]] = Some(i);
            } else {
                non_basic_vars.insert(i);
            }
        }
        let basic_vars = basic_vars.iter().map(|i| i.unwrap()).collect();

        //Construct
        Self {
            constraints,
            basis,
            obj_fn,
            basic_vars,
            non_basic_vars,
            var_map: var_index_map,
        }
    }

    pub fn pivot(&mut self, ind: TableauIx) {
        //update basic and non-basic vars
        self.non_basic_vars.insert(self.basic_vars[ind.i()]);
        self.basic_vars[ind.i()] = ind.j();
        self.non_basic_vars.remove(&ind.j());

        let size = self.basic_vars.len();

        //let denom = self.constraints.get(ind.i(), ind.j()).unwrap_or(&0.0);
        let mut e = CsMat::empty(CompressedStorage::CSC, size);
        let c = self.col(ind.j());
        let denom = c.get(ind.i()).unwrap_or(&0.0);
        for j in 0..size {
            if j == ind.i() {
                for i in 0..size {
                    if i == ind.i() {
                        e.insert(i, j, 1.0 / denom);
                    } else {
                        //let numerator = self.constraints.get(i, ind.j()).unwrap_or(&0.0);
                        let numerator = c.get(i).unwrap_or(&0.0);
                        e.insert(i, j, -numerator / denom);
                    }
                }
            } else {
                e.insert(j, j, 1.0);
            }
        }
        self.basis = &e * &self.basis;
        //self.basis = mul_csr_csr(e.to_csr().view(), self.basis.to_csr().view());
        //self.basis = self.basis.to_csc();
    }

    pub fn r_costs(&self) -> CsVec<f64> {
        let (inds, dstorage) = self.nb_r_costs().unzip();

        CsVec::new_from_unsorted(self.constraints.cols(), inds, dstorage).unwrap()
    }

    pub fn nb_r_costs(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        //want to create an iterable over non basic var coeffs as efficiently as possible\
        let (inds, dstorage): (Vec<usize>, Vec<f64>) = self
            .basic_vars
            .iter()
            .enumerate()
            .filter(|(i, col)| self.obj_fn.get(**col) != None)
            .map(|(i, col)| (i, *self.obj_fn.get(*col).unwrap()))
            .unzip();

        let cb = CsVec::new(self.basis.rows(), inds, dstorage);

        //let mut cb_dot_ab = Array1::from_elem(self.basic_vars.len(), 0.0); //.dot(&self.basis.to_dense().view());
        let cb_dot_ab = &cb * &self.basis;

        self.non_basic_vars.iter().map(move |&i| {
            let nb = self.constraints.outer_view(i).unwrap();
            let cn_i = self.obj_fn.get(i).unwrap_or(&0.0);
            let delta_i = &cb_dot_ab.dot(&nb);
            //let delta_i = nb.dot(&cb_dot_ab);
            //(i, cn_i - delta_i)
            (i, cn_i - delta_i)
        })
    }

    pub fn z(&self) -> f64 {
        let (inds, dstorage): (Vec<usize>, Vec<f64>) = self
            .basic_vars
            .iter()
            .enumerate()
            .filter(|(i, col)| self.obj_fn.get(**col) != None)
            .map(|(i, col)| (i, *self.obj_fn.get(*col).unwrap()))
            .unzip();

        let cb = CsVec::new_from_unsorted(self.constraints.rows(), inds, dstorage).unwrap();

        let cb_dot_ab = &cb * &self.basis; //= self.basis.dot(&cb);

        let b = self
            .constraints
            .outer_view(self.constraints.cols() - 1)
            .unwrap();

        //let b = self.rhs();

        self.obj_fn.get(self.obj_fn.dim() - 1).unwrap_or(&0.0) - cb_dot_ab.dot(&b)
        //cb_dot_ab.dot(&b)
    }

    #[inline]
    pub fn rhs(&self) -> CsVec<f64> {
        self.col(self.constraints.cols() - 1)
    }

    #[inline]
    pub fn col(&self, col: usize) -> CsVec<f64> {
        &self.basis * &self.constraints.outer_view(col).unwrap()

    }

    pub fn remove_vars(&mut self, vars: Vec<Variable>) {
        //removed vars cannot be basic
        assert!(vars
            .iter()
            .map(|v| {
                let ind = self.var_map[v];
                self.basic_vars.contains(&ind)
            })
            .all(|v| !v));

        //build remove var index vec
        let rm_var_inds = vars
            .iter()
            .map(|var| self.var_map[var])
            .collect::<Vec<usize>>();

        //remove vars from map
        vars.iter().for_each(|var| {
            self.var_map.remove(var);
        });

        //remove vars from non-basic set
        rm_var_inds.iter().for_each(|ind| {
            self.non_basic_vars.remove(&ind);
        });

        //reindex var_map
        let mut reindexed_map = HashMap::new();
        self.var_map.iter_mut().for_each(|(var, ind)| {
            let offset = rm_var_inds.iter().filter(|&i| i < ind).count();

            if offset > 0 {
                *ind -= offset;
                reindexed_map.insert(*ind + offset, ind);
            }
        });

        //reindex basic vars
        self.basic_vars.iter_mut().for_each(|ind| {
            if let Some(new_ind) = reindexed_map.get(ind) {
                *ind = **new_ind;
            }
        });

        //reindex non-basic vars
        let non_basic_vars = self
            .non_basic_vars
            .iter()
            .map(|ind| {
                if let Some(new_ind) = reindexed_map.get(ind) {
                    **new_ind
                } else {
                    *ind
                }
            })
            .collect::<FxHashSet<usize>>();
        self.non_basic_vars = non_basic_vars;

        //build new constraint matrix
        let mut constraints = CsMat::empty(CompressedStorage::CSC, self.constraints.rows());

        for (i, col) in self.constraints.outer_iterator().enumerate() {
            if rm_var_inds.contains(&i) {
                continue;
            }

            constraints = constraints.append_outer_csvec(col);
        }
        self.constraints = constraints;

        //build new obj fn
        let (inds, dstorage): (Vec<usize>, Vec<f64>) = self
            .obj_fn
            .iter()
            .filter(|(i, coeff)| !rm_var_inds.contains(i))
            .map(|(i, coeff)| {
                if let Some(ind) = reindexed_map.get(&i) {
                    (**ind, *coeff)
                } else {
                    (i, *coeff)
                }
            })
            .unzip();

        self.obj_fn = CsVec::new_from_unsorted(self.constraints.cols(), inds, dstorage).unwrap();
    }
}

#[derive(Debug, Clone)]
pub struct Tableau {
    pub(crate) tbl: Array2<f64>,
    pub(crate) vars: HashMap<Variable, usize>,
    pub(crate) basic_vars: Vec<usize>,
}

impl From<&SparseTableau> for Tableau {
    fn from(sp: &SparseTableau) -> Self {
        let rows = sp.constraints().rows() + 1;
        let cols = sp.constraints().cols() + 1;
        let mut tbl = Array2::from_elem((rows, cols), 0.0);
        //populate A
        for col in (0..cols - 2).into_iter() {
            tbl.slice_mut(s![..-1, col]).assign(&sp.col(col).to_dense());
        }
        //populate b
        tbl.slice_mut(s![..-1, cols - 1])
            .assign(&sp.col(cols - 2).to_dense());
        //populate obj fn
        tbl.slice_mut(s![-1, ..-1]).assign(&sp.r_costs().to_dense());
        //populate z
        tbl[[rows - 1, cols - 2]] = 1.0;
        //populate obj_fn val
        tbl[[rows - 1, cols - 1]] = sp.z();
        Self {
            tbl,
            vars: sp.var_map().clone(),
            basic_vars: sp.basic_vars().clone(),
        }
    }
}

impl Tableau {
    //constructor
    pub fn new(tbl: Array2<f64>, vars: HashMap<Variable, usize>, basic_vars: Vec<usize>) -> Self {
        Self {
            tbl,
            vars,
            basic_vars,
        }
    }

    //TODO!
    pub fn set_obj_fn(&mut self, _row: Vec<f64>) {}
    pub fn add_constraint(&mut self, _row: Vec<f64>) {}
    pub fn remove_constraint(&mut self, _row: usize) {}

    pub fn append_row(&mut self, row: Vec<f64>) {
        let arr = ArrayView::from(row.as_slice())
            .into_shape((1, row.len()))
            .unwrap();
        println!("Pre append:\n{}", self.tbl);
        self.tbl.append(Axis(0), arr).unwrap();
        println!("Post append:\n{}", self.tbl);
    }

    pub fn tbl(&self) -> &Array2<f64> {
        &self.tbl
    }

    pub fn filter_rows(&mut self, mask: Vec<bool>) {
        assert!(self.tbl.shape()[0] == mask.len(), "Incorrect mask length");
        let elems = self
            .tbl
            .axis_iter(Axis(0))
            .zip(mask.iter())
            .filter(|(_row, keep)| **keep)
            .flat_map(|(row, _keep)| row.to_vec());

        let new_n_rows = mask.len() - mask.iter().filter(|m| !**m).count();

        self.tbl = Array::from_iter(elems)
            .into_shape((new_n_rows, self.tbl.shape()[1]))
            .unwrap();
    }

    pub fn filter_cols(&mut self, mask: Vec<bool>) {
        assert!(self.tbl.shape()[1] == mask.len(), "Incorrect mask length");

        //create index map
        let _ind_var_map = self
            .vars
            .iter()
            .map(|(var, ind)| (*ind, var.clone()))
            .collect::<HashMap<usize, Variable>>();

        //update basic var indexes
        let _basic_offset = vec![0; self.basic_vars.len()];
        let _var_offset = HashMap::<Variable, usize>::new();
        let mask_indices = mask
            .iter()
            .enumerate()
            .filter(|(_i, &flag)| !flag)
            .map(|(i, &_flag)| i)
            .collect::<Vec<usize>>();
        self.basic_vars.iter_mut().for_each(|bvar_ind| {
            *bvar_ind -= mask_indices
                .iter()
                .filter(|mask_ind| *mask_ind < bvar_ind)
                .count();
        });

        self.vars.iter_mut().for_each(|(_var, ind)| {
            *ind -= mask_indices
                .iter()
                .filter(|mask_ind| *mask_ind < ind)
                .count();
        });

        //can do better than double transpose
        self.tbl.swap_axes(0, 1);
        self.filter_rows(mask);
        self.tbl.swap_axes(0, 1);

       
    }

    pub fn rref<T: Into<Slice>, U: Into<Slice>>(&mut self, rows: T, cols: U) {
        //convert matrix into reduced row echelon form
        //can this be faster?
        let tol = 1.0e-6;

        //for now create a copy
        let mut view = self.tbl.slice_mut(s![rows.into(), cols.into()]);

        let m = view.shape()[0]; //rows
        let n = view.shape()[1]; //cols

        let mut pivots = Vec::with_capacity(std::cmp::min(m, n));

        let mut i = 0;
        let mut j = 0;

        while (i < m) && (j < n) {
            //find pivot row
            if let Some(pivot) = view
                .slice(s![i.., j])
                .iter()
                .enumerate()
                .filter(|(_, v)| v.abs() >= tol)
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap_or(Ordering::Equal))
                .map(|(index, _)| index)
            {
                let k = i + pivot;
                //swap row i and k
                if k != i {
                    let (row_k, row_i) = view.multi_slice_mut((s![k, ..], s![i, ..]));
                    Zip::from(row_k).and(row_i).for_each(::std::mem::swap);
                    self.basic_vars.swap(i, k);
                }

                //normalize pivot row
                let div = view[[i, j]];
                let mut row_r = view.slice_mut(s![i, j..]);
                row_r /= div;

                //eliminate column
                let col = view.slice(s![.., j]).to_owned().into_shape((m, 1)).unwrap();
                let row = view
                    .slice(s![i, j..])
                    .to_owned()
                    .into_shape((1, n - j))
                    .unwrap();

                let dot = col.dot(&row);

                let row_i = view.slice(s![i, j..]).to_owned();
                let mut tmp = view.slice_mut(s![.., j..]);
                tmp -= &dot;

                let tmp_row_i = view.slice_mut(s![i, j..]);
                row_i.assign_to(tmp_row_i);

                pivots.push(j);
                i += 1;
                j += 1;
            } else {
                let mut slice = view.slice_mut(s![i..m, j]);
                slice.fill(0.0_f64);
                j += 1;
            }
        }
        self.basic_vars = pivots;
    }

    //pivot
    #[inline(always)]
    pub fn pivot(&mut self, pivot_ind: &TableauIx) {
        //assert row and col in valid range
        assert!(pivot_ind.i() < self.tbl.shape()[0] - 1);
        assert!(pivot_ind.j() < self.tbl.shape()[1]);

        //set coefficients in pivot row
        let div = self.tbl[[pivot_ind.i(), pivot_ind.j()]];
        for j in 0..self.tbl.shape()[1] {
            self.tbl[[pivot_ind.i(), j]] /= div;
        }

        //pivot
        for i in 0..self.tbl.shape()[0] {
            //skip pivot row
            if i == pivot_ind.i() {
                continue;
            }
            let ratio = self.tbl[[i, pivot_ind.j()]];
            for j in 0..self.tbl.shape()[1] {
                //calc new coefficients
                self.tbl[[i, j]] -= self.tbl[[pivot_ind.i(), j]] * ratio;
            }
        }

        //update basic vars
        self.basic_vars[pivot_ind.i()] = pivot_ind.j();
    }

    pub fn tbl_variable_columns_string(&self) -> String {
        "{:>}".repeat(self.tbl.shape()[1])
    }

    pub fn row_strings<T: Into<Slice>, U: Into<Slice>>(
        &self,
        rows: T,
        cols: U,
    ) -> Vec<Vec<String>> {
        let tbl_slice = self.tbl.slice(s![rows.into(), cols.into()]);
        let mut row_strings =
            vec![vec!["".to_string(); tbl_slice.shape()[1]]; tbl_slice.shape()[0]];

        for (ci, col) in tbl_slice.columns().into_iter().enumerate() {
            let (max_whole_len, max_decimal_len) = col.iter().fold((0, 0), |len, num| {
                let s_num = num.to_string();
                let mut s_num_it = s_num.split('.');
                (
                    std::cmp::max(len.0, s_num_it.next().unwrap_or("").len()),
                    std::cmp::max(len.1, s_num_it.next().unwrap_or("").len()),
                )
            });

            col.iter().enumerate().for_each(|(ri, coeff)| {
                // if ci != 0 && coeff.signum() == 1.0 {
                //     row_strings[ri][ci] += &" + ";
                // } else if coeff.signum() != 1.0 {
                //     row_strings[ri][ci] += &" - ";
                // }

                let sign = if coeff.signum() != 1.0 {
                    "-".to_string()
                } else {
                    "".to_string()
                };

                let s_num = coeff.abs().to_string();
                let mut s_num_it = s_num.split('.');

                let s_whole_num = match s_num_it.next() {
                    Some(whole_num) => {
                        format!("{:>width$}", sign + whole_num, width = max_whole_len)
                    }
                    _ => {
                        format!("{:>width$}", sign, width = max_whole_len)
                    }
                };
                let mut s_decimal_num = match s_num_it.next() {
                    Some(decimal_num) => {
                        let s = ".".to_string();
                        s + format!("{:<width$}", decimal_num, width = max_decimal_len).as_str()
                    }
                    _ => {
                        format!("{:<width$}", "", width = max_decimal_len + 1)
                    }
                };

                //trim s_decimal_num
                if s_decimal_num.len() > 2 {
                    s_decimal_num.replace_range(2.., "");
                }
                row_strings[ri][ci] += s_whole_num.as_str();
                row_strings[ri][ci] += s_decimal_num.as_str();
            });
        }
        row_strings
    }

    pub fn as_tabular(&self) -> tabular::Table {
        let mut table = tabular::Table::new(self.tbl_variable_columns_string().as_str());

        self.row_strings(.., ..).iter().for_each(|row_vec| {
            let row = tabular::Row::from_cells(row_vec);
            table.add_row(row);
        });

        table
    }

    pub fn print_pivot(&self, _pivot_ind: TableauIx) {}

    pub fn as_table_builder(&self, info: TblPrintInfo) -> Builder {
        let mut builder = Builder::default();
        let _offset = match info {
            TblPrintInfo::Default => 0,
            TblPrintInfo::Pivot(_) => 1,
        };

        let [top_padding, bottom_padding, left_padding, right_padding] = info.get_padding();

        //build header row
        let mut header =
            vec!["".to_string(); self.tbl.shape()[1] + 1 + left_padding + right_padding];
        header[left_padding] = "Basic Vars".to_string();
        self.vars
            .iter()
            .for_each(|(var, ind)| header[ind + 1 + left_padding] = var.name().to_string());
        header[self.tbl.shape()[1] - 1 + left_padding] = "Z".to_string();
        header[self.tbl.shape()[1] + left_padding] = "rhs".to_string();

        //add header and required empty rows
        match info {
            TblPrintInfo::Default => {
                builder.add_record(&header);
            }
            TblPrintInfo::Pivot(_) => {
                for _ in 0..top_padding {
                    builder.add_record(&vec!["".to_string(); header.len()]);
                }
                builder.add_record(&header);
            }
        }
        //data rows
        self.row_strings(.., ..)
            .into_iter()
            .enumerate()
            .for_each(|(i, row_coeffs)| {
                let mut row = match info {
                    TblPrintInfo::Default => Vec::new(),
                    TblPrintInfo::Pivot(_) => vec!["".to_string(); 1],
                };
                //name of basic var for row
                if i < self.tbl.shape()[0] - 1 {
                    row.push(header[self.basic_vars[i] + 1 + left_padding].clone());
                } else {
                    row.push("Z".to_string());
                };

                row.extend(row_coeffs);
                builder.add_record(row);
            });

        match info {
            TblPrintInfo::Default => {}
            TblPrintInfo::Pivot(_) => {
                for _ in 0..bottom_padding {
                    builder.add_record(&vec!["".to_string(); header.len()]);
                }
            }
        }

        builder
    }

    pub fn as_table(&self, info: TblPrintInfo) -> Table {
        let [top_padding, bottom_padding, left_padding, right_padding] = info.get_padding();

        //build table and get dims
        let mut table = self.as_table_builder(info).build().with(Style::empty());
        let (nrows, ncols) = table.shape();

        //format table borders
        table = table
            .with(
                Modify::new(
                    Rows::single(top_padding)
                        .not(Columns::new(0..left_padding))
                        .not(Columns::new(ncols - right_padding..ncols)),
                )
                .with(Border::default().bottom('─')),
            )
            .with(
                Modify::new(
                    Rows::single(nrows - 1 - bottom_padding)
                        .not(Columns::new(0..left_padding))
                        .not(Columns::new(ncols - right_padding..ncols)),
                )
                .with(Border::default().top('─')),
            )
            .with(
                Modify::new(
                    Columns::single(left_padding)
                        .not(Rows::new(0..top_padding))
                        .not(Rows::new(nrows - bottom_padding..nrows)),
                )
                .with(Border::default().right('│')),
            )
            .with(
                Modify::new(
                    Columns::single(ncols - 1 - right_padding)
                        .not(Rows::new(0..top_padding))
                        .not(Rows::new(nrows - bottom_padding..nrows)),
                )
                .with(Border::default().left('│')),
            )
            .with(
                Modify::new(Cell(top_padding, left_padding))
                    .with(Border::default().bottom_right_corner('┼')),
            )
            .with(
                Modify::new(Cell(nrows - 1 - bottom_padding, left_padding))
                    .with(Border::default().top_right_corner('┼')),
            )
            .with(
                Modify::new(Cell(nrows - 1 - bottom_padding, ncols - 1 - right_padding))
                    .with(Border::default().top_left_corner('┼')),
            )
            .with(
                Modify::new(Cell(top_padding, ncols - 1 - right_padding))
                    .with(Border::default().bottom_left_corner('┼')),
            );

        //annotate table if info is pivot variant
        match info {
            TblPrintInfo::Default => table,
            TblPrintInfo::Pivot(TableauIx { i, j }) => {
                table = table
                    .with(
                        Modify::new(Cell(i + 1 + top_padding, 0))
                            .with(|_s: &str| "Leaving".to_string()),
                    )
                    .with(
                        Modify::new(Cell(0, j + 1 + left_padding))
                            .with(|_s: &str| "Entering".to_string()),
                    );

                for (row, val) in self.tbl.slice(s![..-1, j]).iter().enumerate() {
                    let denominator = self.tbl[[i, j]];
                    let sign = if val.signum() == denominator.signum() {
                        "-"
                    } else {
                        "+"
                    };
                    let annotation = format!(
                        "R{} = R{} {} {:.2}/{:.2}*R{}",
                        row,
                        row,
                        sign,
                        val.abs(), //numerator.abs(),
                        denominator.abs(),
                        i
                    );
                    table = table.with(
                        Modify::new(Cell(row + 1 + top_padding, ncols - 1))
                            .with(|_s: &str| annotation.to_string()),
                    );
                }

                table
            }
        }
    }
}

impl fmt::Display for Tableau {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let table = self.as_table(TblPrintInfo::Default);
        writeln!(f, "{}", table)

        //let table = self.as_tabular();

        // write!(f, "{}", table)
    }
}
