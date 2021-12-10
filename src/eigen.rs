use nalgebra::DMatrix;

pub fn add_sub_matrix(
    mat: &mut DMatrix<f64>,
    (row, col): (usize, usize),
    submat: &DMatrix<f64>,
) -> Result<(), String> {
    //check if col , row + size(submat) does not overpass size(mat)
    if (submat.nrows() + row) > mat.nrows() || (submat.ncols() + col) > mat.ncols() {
        return Err("Size submatrix and position are incompatible with matrix".to_string());
    }

    for (i, col_submat) in submat.column_iter().enumerate() {
        let pos_col_mat = i + col;
        for (j, value) in col_submat.iter().enumerate() {
            let pos_row_mat = j + row;
            mat[(pos_row_mat, pos_col_mat)] += *value;
        }
    }
    Ok(())
}
