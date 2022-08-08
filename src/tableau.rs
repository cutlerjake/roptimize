use crate::var::Variable;

use ndarray::{s, stack, Array, Array1, Array2, ArrayView, Axis, Slice, Zip};

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;
use std::ops::RangeBounds;

use tabled::{
    builder::Builder,
    object::{Cell, Columns, Object, Rows, Segment},
    style::Border,
    Highlight, Modify, Style, Table,
};
use tabular;

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

        //create index map
        let ind_var_map = self
            .vars
            .iter()
            .map(|(var, ind)| (*ind, var.clone()))
            .collect::<HashMap<usize, Variable>>();

        //update basic var indexes
        let mut basic_offset = vec![0; self.basic_vars.len()];
        let mut var_offset = HashMap::<Variable, usize>::new();
        let mask_indices = mask
            .iter()
            .enumerate()
            .filter(|(i, &flag)| !flag)
            .map(|(i, &flag)| i)
            .collect::<Vec<usize>>();
        self.basic_vars.iter_mut().for_each(|bvar_ind| {
            *bvar_ind -= mask_indices
                .iter()
                .filter(|mask_ind| *mask_ind < bvar_ind)
                .count();
        });

        self.vars.iter_mut().for_each(|(var, ind)| {
            *ind -= mask_indices
                .iter()
                .filter(|mask_ind| *mask_ind < ind)
                .count();
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
                let mut s_num_it = s_num.split(".");
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
                let mut s_num_it = s_num.split(".");

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

    pub fn print_pivot(&self, pivot_ind: TableauIx) {}

    pub fn as_table_builder(&self, info: TblPrintInfo) -> Builder {
        let mut builder = Builder::default();
        let offset = match info {
            TblPrintInfo::Default => 0,
            TblPrintInfo::Pivot(_) => 1,
        };

        let [top_padding, bottom_padding, left_padding, right_padding] = info.get_padding();

        //build header row
        let mut header =
            vec!["".to_string(); self.tbl.shape()[1] + 1 + left_padding + right_padding];
        header[0 + left_padding] = "Basic Vars".to_string();
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
                    Rows::single(0 + top_padding)
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
                    Columns::single(0 + left_padding)
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
                Modify::new(Cell(0 + top_padding, 0 + left_padding))
                    .with(Border::default().bottom_right_corner('┼')),
            )
            .with(
                Modify::new(Cell(nrows - 1 - bottom_padding, 0 + left_padding))
                    .with(Border::default().top_right_corner('┼')),
            )
            .with(
                Modify::new(Cell(nrows - 1 - bottom_padding, ncols - 1 - right_padding))
                    .with(Border::default().top_left_corner('┼')),
            )
            .with(
                Modify::new(Cell(0 + top_padding, ncols - 1 - right_padding))
                    .with(Border::default().bottom_left_corner('┼')),
            );

        //annotate table if info is pivot variant
        match info {
            TblPrintInfo::Default => table,
            TblPrintInfo::Pivot(TableauIx { i, j }) => {
                table = table
                    .with(
                        Modify::new(Cell(i + 1 + top_padding, 0))
                            .with(|s: &str| "Leaving".to_string()),
                    )
                    .with(
                        Modify::new(Cell(0, j + 1 + left_padding))
                            .with(|s: &str| "Entering".to_string()),
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
                        val.abs(),//numerator.abs(),
                        denominator.abs(),
                        i
                    );
                    table = table.with(
                        Modify::new(Cell(row + 1 + top_padding, ncols - 1))
                            .with(|s: &str| annotation.to_string()),
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
        write!(f, "{}\n", table)

        //let table = self.as_tabular();

        // write!(f, "{}", table)
    }
}
