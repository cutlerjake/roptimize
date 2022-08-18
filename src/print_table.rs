use ndarray::{
    s, Array2, Dim, DimAdd, SliceArg, SliceInfo, SliceInfoElem, SliceNextDim,
};


use unicode_width::UnicodeWidthStr;

use std::fmt;

pub struct BorderStyle {
    top: char,
    right: char,
    bottom: char,
    left: char,

    down_and_left: char,
    down_and_right: char,

    up_and_right: char,
    up_and_left: char,

    vertical_and_right: char,
    vertical_and_left: char,

    up_and_horizontal: char,
    down_and_horizontal: char,

    vertical_and_horizontal: char,
}

impl Default for BorderStyle {
    fn default() -> Self {
        Self {
            top: '─',
            right: '│',
            bottom: '─',
            left: '│',

            down_and_left: '┐',
            down_and_right: '┌',

            up_and_right: '┘',
            up_and_left: '└',

            vertical_and_right: '├',
            vertical_and_left: '┤',

            up_and_horizontal: '┴',
            down_and_horizontal: '┬',

            vertical_and_horizontal: '┼',
        }
    }
}

pub struct PrintTableCell {
    data: String,
    width: usize,
    border: BorderStyle
}

impl PrintTableCell {
    pub fn new(data: String) -> Self {
        let width = UnicodeWidthStr::width(data.as_str());
        Self { data, width, border: BorderStyle::default() }
    }

    pub fn data(&self) -> &str {
        &self.data
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn set_data(self, data: String) -> Self {
        Self::new(data)
    }
}

#[derive(Copy, Clone, Debug, Hash)]
pub struct PrintTableIx {
    pub i: usize,
    pub j: usize,
}

#[derive(Copy, Clone, Debug, Hash)]
pub struct CellBorderIxs {
    //horizontal borders
    pub top: PrintTableIx,
    pub bottom: PrintTableIx,

    //vertical borders
    pub left: PrintTableIx,
    pub right: PrintTableIx,

    //corners
    pub top_left: PrintTableIx,
    pub top_right: PrintTableIx,
    pub bottom_right: PrintTableIx,
    pub bottom_left: PrintTableIx,
}

pub struct PrintTable {
    cells: Array2<PrintTableCell>,
    borders: Array2<char>,
}

impl PrintTable {
    pub fn new(cells: Array2<PrintTableCell>) -> Self {
        let borders = Array2::from_elem((cells.shape()[0], cells.shape()[1]), ' ');
        Self { cells, borders }
    }

    pub fn col_widths(&self) -> Vec<usize> {
        self.cells
            .columns()
            .into_iter()
            .map(|col| col.iter().fold(0, |l, cell| std::cmp::max(l, cell.width())))
            .collect::<Vec<usize>>()
    }

    pub fn row_heights(&self) -> Vec<usize> {
        self.cells
            .rows()
            .into_iter()
            .map(|row| {
                row.iter()
                    .fold(0, |l, cell| std::cmp::max(l, cell.data.split('\n').count()))
            })
            .collect::<Vec<usize>>()
    }

    pub fn cell_border_inds(&self, ind: PrintTableIx) -> CellBorderIxs {
        assert!(ind.i < self.cells.shape()[0] && ind.j < self.cells.shape()[1]);

        CellBorderIxs {
            top: PrintTableIx {
                i: ind.i - 1,
                j: ind.j,
            },
            bottom: PrintTableIx {
                i: ind.i + 1,
                j: ind.j,
            },

            left: PrintTableIx {
                i: ind.i,
                j: ind.j - 1,
            },
            right: PrintTableIx {
                i: ind.i,
                j: ind.j + 1,
            },

            top_left: PrintTableIx {
                i: ind.i - 1,
                j: ind.j - 1,
            },
            top_right: PrintTableIx {
                i: ind.i - 1,
                j: ind.j + 1,
            },
            bottom_right: PrintTableIx {
                i: ind.i + 1,
                j: ind.j + 1,
            },
            bottom_left: PrintTableIx {
                i: ind.i + 1,
                j: ind.j - 1,
            },
        }
    }

    pub fn fill_border<T, U>(&mut self, b: char, rows: T, cols: U)
    where
        SliceInfoElem: From<U>,
        SliceInfoElem: From<T>,
        T: SliceNextDim,
        U: SliceNextDim,
        <T as SliceNextDim>::OutDim: DimAdd<<U as SliceNextDim>::OutDim>,
        <T as SliceNextDim>::InDim: DimAdd<<U as SliceNextDim>::InDim>,
        SliceInfo<
            [SliceInfoElem; 2],
            <<T as SliceNextDim>::InDim as DimAdd<<U as SliceNextDim>::InDim>>::Output,
            <<T as SliceNextDim>::OutDim as DimAdd<<U as SliceNextDim>::OutDim>>::Output,
        >: SliceArg<Dim<[usize; 2]>>,
    {
        self.borders.slice_mut(s![rows, cols]).fill(b);
    }

    pub fn testy_boy(&mut self) {
        self.fill_border('-', .., 1);
    }
}

impl fmt::Display for PrintTable {
    fn fmt(&self, _f: &mut fmt::Formatter) -> fmt::Result {
        //print table row by row
        let _col_widths = self.col_widths();
        let _row_heights = self.row_heights();

        for(i, _row) in self.borders.rows().into_iter().enumerate() {
            if i%2 != 0 {
                //print cells
            }
            //then print bottom border

        }

        for (_i, _row) in self.cells.rows().into_iter().enumerate() {}
        Ok(())
    }
}
