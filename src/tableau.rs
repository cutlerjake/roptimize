use crate::var::Variable;

use ndarray::{s, stack, Array, Array1, Array2, ArrayView, Axis, Slice, Zip};

use std::cmp::Ordering;
use std::collections::HashMap;
use std::ops::RangeBounds;

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

#[derive(Debug, Clone)]
pub struct Tableau {
    pub(crate) tbl: Array2<f64>,
    pub(crate) vars: HashMap<Variable, usize>,
    pub(crate) basic_vars: Vec<usize>,
}

impl Tableau {
    //constructor
    pub fn new(tbl: Array2<f64>, vars: HashMap<Variable, usize>, basic_vars: Vec<usize>) -> Self {
        return Self {
            tbl,
            vars,
            basic_vars,
        };
    }

    //TODO!
    pub fn set_obj_fn(&mut self, row: Vec<f64>) {}
    pub fn add_constraint(&mut self, row: Vec<f64>) {}
    pub fn remove_constraint(&mut self, row: usize) {}

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

        //update basic var indexes
        mask.iter().enumerate().for_each(|(var_ind, var)| {
            self.basic_vars.iter_mut().for_each(|bvar_ind| {
                if *bvar_ind > var_ind {
                    *bvar_ind -= 1;
                }
            });
        });
        //can do better than double transpose
        self.tbl.swap_axes(0, 1);
        self.filter_rows(mask);
        self.tbl.swap_axes(0, 1);

        // let elems = self
        // .tbl
        // .axis_iter(Axis(1))
        // .zip(mask.iter())
        // .filter(|(_col, keep)| **keep)
        // .flat_map(|(col, _keep)| col.to_vec());

        // let new_n_cols = mask.len() - mask.iter().filter(|m| !**m).count();

        // self.tbl = Array::from_iter(elems)
        //     .into_shape((self.tbl.shape()[0], new_n_cols))
        //     .unwrap();
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
}
