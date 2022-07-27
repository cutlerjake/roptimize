use ndarray::{s, ArrayView, Array2, Axis, Zip, stack, Array1};

use std::cmp::Ordering;

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

#[derive(Debug)]
pub struct Tableau {
    pub(crate) tbl: Array2<f64>,
}

impl Tableau {
    //constructor
    pub fn new(tbl: Array2<f64>) -> Self {
        return Self { tbl };
    }

    pub fn append_row(&mut self, row: Vec<f64>) {
        let arr = ArrayView::from(row.as_slice()).into_shape((1, row.len())).unwrap();
        println!("Pre append:\n{}", self.tbl);
        self.tbl.append(Axis(0), arr).unwrap();
        println!("Post append:\n{}", self.tbl);
        
    }

    pub fn tbl(&self) -> &Array2<f64> {
        &self.tbl
    }

    pub fn rref(&mut self) -> Vec<usize> {
        //convert matrix into reduced row echelon form
        //can this be faster?
        let tol = 1.0e-6;
        let m = self.tbl.shape()[0]; //rows
        let n = self.tbl.shape()[1]; //cols

        let mut pivots = Vec::with_capacity(std::cmp::min(m, n));

        let mut i = 0;
        let mut j = 0;

        while (i < m) && (j < n) {
            //find pivot row
            if let Some(pivot) = self
                .tbl
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
                    let (row_k, row_i) = self.tbl.multi_slice_mut((s![k, ..], s![i, ..]));
                    Zip::from(row_k).and(row_i).for_each(::std::mem::swap);
                }

                //normalize pivot row
                let div = self.tbl[[i, j]];
                let mut row_r = self.tbl.slice_mut(s![i, j..]);
                row_r /= div;
                //eliminate column
                let col = self
                    .tbl
                    .slice(s![.., j])
                    .to_owned()
                    .into_shape((m, 1))
                    .unwrap();
                let row = self
                    .tbl
                    .slice(s![i, j..])
                    .to_owned()
                    .into_shape((1, n - j))
                    .unwrap();

                let mut dot = col.dot(&row);

                let row_i = self.tbl.slice(s![i, j..]).to_owned();
                let mut tmp = self.tbl.slice_mut(s![.., j..]);
                tmp -= &dot;

                let tmp_row_i = self.tbl.slice_mut(s![i, j..]);
                row_i.assign_to(tmp_row_i);

                pivots.push(j);
                i += 1;
                j += 1;
                
            } else {
                let mut slice = self.tbl.slice_mut(s![i..m, j]);
                slice.fill(0.0_f64);
                j += 1;
            }
        }
        pivots
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
    }
}
