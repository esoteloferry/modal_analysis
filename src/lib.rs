extern crate nalgebra as na;
use std::ops::{DivAssign, Mul};

use na::DMatrix;

#[derive(Debug)]
pub struct Mode {
    pub frequency: Vec<f64>,
    pub eigenvalues: DMatrix<f64>,
    pub eigenvectors_normalized: DMatrix<f64>,
    pub eigenvectors: DMatrix<f64>,
}

pub fn eigen(mass_matrix: &DMatrix<f64>, stiffness_matrix: DMatrix<f64>) -> Result<Mode, String> {
    //Checks
    if !mass_matrix.is_square() {
        return Err(String::from("Mass matrix is not square"));
    }
    if !stiffness_matrix.is_square() {
        return Err(String::from("Stiffness matrix is not square"));
    }
    let shape_mass = mass_matrix.shape();
    let shape_stiff = stiffness_matrix.shape();
    if shape_mass != shape_stiff {
        return Err(String::from(format!(
            "Matrix and stiffness matrix have different dimensions. Mass={}, Stiffness={}",
            shape_mass.0, shape_stiff.0
        )));
    }

    let dim = shape_mass.0;

    println!("M_matrix : {}", mass_matrix);
    // let inv_m_matrix = m_matrix.clone().try_inverse();
    let inv_m_matrix = match mass_matrix.clone().try_inverse() {
        Some(i) => i,
        None => return Err(String::from("Matrix could not be inversed")),
    };
    let a_matrix = inv_m_matrix * stiffness_matrix;

    // let eigen = a_matrix.eigenvalues().expect("Error getting eigenvalues");
    let eigen = match a_matrix.eigenvalues() {
        Some(i) => i,
        None => return Err(String::from("Error getting eigenvalues")),
    };
    println!("A_matrix :{}", a_matrix);

    let mut eigen_vect: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]);

    for (i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec = vec![0.0; dim];
        // Determinant is zero, which implies that there is only one independent equation for the two unknowns.
        // Hence, we can use either equation to solve for the ratios of the unknowns, therefore we
        // set b_vec[i]=1.0
        b_vec[i] = 1.0;
        let b = DMatrix::from_vec(dim, 1, b_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        // let mut x = decomp.solve(&b).expect("Linear resolution failed.");
        let mut x = match decomp.solve(&b) {
            Some(i) => i,
            None => return Err(String::from("Linear resolution failed")),
        };
        x.div_assign(x[i]);
        // add_sub_matrix(&mut eigen_vect, (0, i), &x).expect("Error adding sub_matrix");
        match add_sub_matrix(&mut eigen_vect, (0, i), &x) {
            Ok(_) => (),
            Err(e) => return Err(e),
        };
    }
    println!("---------------------------------------------------------");

    let omega = eigen.clone().map(|e| e.sqrt()).as_slice().to_vec();
    // let om2 = eigen.clone().map(|e| e.sqrt());
    // println!("Omega : {:?}", omega);
    // println!("EigenValues : {}", eigen);
    // println!("EigenVectors : {}", eigen_vect);

    println!("---------------------------------------------------------");

    // Get modal mass and modal stiffness matrix
    let eigen_vect_trans = eigen_vect.transpose();

    let modal_mass = eigen_vect_trans
        .clone()
        .mul(mass_matrix.clone())
        .mul(eigen_vect.clone());

    // println!("m_se: {}", m_se);
    // println!("{}", eigen_vect);
    // println!("{}", eigen_vect_trans);
    // println!("modal mass: {}", modal_mass);

    let mut eigen_vect_norm = eigen_vect.clone();
    // println!("{}", eigen_vect_norm[(0, 0)]);
    for (i, mut col) in eigen_vect_norm.column_iter_mut().enumerate() {
        col.div_assign(modal_mass[(i, i)].powf(0.5));
    }

    // println!("eigenvectors norm : {}", eigen_vect_norm);
    let modal_mass_norm = eigen_vect_norm
        .transpose()
        .clone()
        .mul(mass_matrix.clone())
        .mul(eigen_vect_norm.clone());
    println!("modal mass norm : {}", modal_mass_norm);
    //
    // let diag_c = eigen.clone().div_assign(2.0);
    // println!("Diag_c : {}", diag_c);
    // let mut c_natural = DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]);
    // for (i, _) in c_natural.diagonal().iter_mut().enumerate() {
    //     c_natural[(i, i)] = 2.0 * 0.0025 * omega[(i)];
    // }
    // println!("{}", c_natural);
    //
    // let inv_eigen_vect_trans = match eigen_vect_trans.clone().try_inverse() {
    //     Some(i) => i,
    //     None => panic!("Matrix could not be inversed"),
    // };
    // let inv_eigen_vect = match eigen_vect.clone().try_inverse() {
    //     Some(i) => i,
    //     None => panic!("Matrix could not be inversed"),
    // };
    //
    // let c_damping = inv_eigen_vect_trans.mul(c_natural).mul(inv_eigen_vect);
    //
    // println!("C damping : {}", c_damping);
    // let mut omega_vec: Vec<f64> = Vec::with_capacity(dim);
    // for om in omega.iter() {
    //     omega_vec.push(om.to_owned())
    // }
    // let mut eigen_vect_norm_vec: Vec<Vec<f64>> = Vec::new();
    // for col in eigen_vect_norm.column_iter() {
    //     eigen_vect_norm_vec.push(col.as_slice().to_owned())
    // }
    // let mut modes: Vec<Mode> = Vec::with_capacity(dim);
    // TODO: what effect has a change in column position of eigenvector matrix? to x variables?
    // for (omega, shape) in omega_vec.iter().zip(eigen_vect_norm_vec) {
    //     modes.push(Mode {
    //         frequency: *omega,
    //         mode: shape,
    //     });
    // }
    // modes.sort_by(|a, b| a.frequency.partial_cmp(&b.frequency).unwrap());
    let mode = Mode {
        frequency: omega,
        eigenvalues: DMatrix::from_diagonal(&eigen),
        eigenvectors_normalized: eigen_vect_norm,
        eigenvectors: eigen_vect,
    };
    return Ok(mode);
}

fn add_sub_matrix(
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

#[cfg(test)]
mod tests {

    use na::DVector;

    use super::*;
    #[test]
    fn it_works() {
        let mass_mat_vec: Vec<f64> = vec![1.0, 0.0, 0.0, 2.0];
        let stiff_mat_vec: Vec<f64> = vec![2.0, -1.0, -1.0, 2.0];

        let dim: usize = 2;
        let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
        let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);

        let modes = eigen(&m_matrix, k_matrix);
        assert!(modes.is_ok());
        let modes = modes.unwrap();
        //Print
        println!("Frequencies : {:?}", modes.frequency);
        println!("Eigenvalues : {}", modes.eigenvalues);
        println!("Eigenvectors : {}", modes.eigenvectors);
        println!("Eigenvectors normalized: {}", modes.eigenvectors_normalized);

        let (w1, w2) = (1.53819, 0.79623);
        let sol_freq_vec = DVector::from_vec(vec![w1, w2]);
        let sol_eigenval = DMatrix::from_vec(dim, dim, vec![w1 * w1, 0.0, 0.0, w2 * w2]);
        let sol_eigenvec_norm =
            DMatrix::from_vec(dim, dim, vec![0.88807, -0.32506, 0.4597, 0.62796]);
        // println!("{:?}", modes.eigenvectors.clone().mul(mass_mat_vec));
        // assert_eq!(2 + 2, 4);
        let freq_vec_na = DVector::from_vec(modes.frequency);
        let eps = 0.0001;
        assert!(freq_vec_na.relative_eq(&sol_freq_vec, eps, eps));
        assert!(modes.eigenvalues.relative_eq(&sol_eigenval, eps, eps));
        assert!(modes
            .eigenvectors_normalized
            .relative_eq(&sol_eigenvec_norm, eps, eps));
    }
}
