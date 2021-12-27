extern crate nalgebra as na;
use std::ops::{DivAssign, Mul};

use na::DMatrix;

pub fn eigen(dim_: usize, m_mat_vec: Vec<f64>, k_mat_vec: Vec<f64>) -> Vec<f64> {
    //Input:
    let dim = dim_;
    let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, m_mat_vec);

    let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, k_mat_vec);

    // ----------------------------------------------

    println!("M_matrix : {}", m_matrix);
    let inv_m_matrix = match m_matrix.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Matrix could not be inversed"),
    };
    let a_matrix = inv_m_matrix * k_matrix;

    let eigen = a_matrix.eigenvalues().expect("Error getting eigenvalues");
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
        let mut x = decomp.solve(&b).expect("Linear resolution failed.");
        x.div_assign(x[i]);
        add_sub_matrix(&mut eigen_vect, (0, i), &x).expect("Error adding sub_matrix");
    }
    println!("---------------------------------------------------------");

    let omega = eigen.clone().map(|e| e.sqrt());
    println!("Omega : {}", omega);
    println!("EigenValues : {}", eigen);
    println!("EigenVectors : {}", eigen_vect);

    println!("---------------------------------------------------------");

    // Get modal mass and modal stiffness matrix
    let eigen_vect_trans = eigen_vect.transpose();

    let modal_mass = eigen_vect_trans
        .clone()
        .mul(m_matrix.clone())
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

    println!("eigenvectors norm : {}", eigen_vect_norm);
    let modal_mass_norm = eigen_vect_norm
        .transpose()
        .clone()
        .mul(m_matrix.clone())
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
    let mut omega_vec = vec![];
    for om in omega.iter() {
        omega_vec.push(om.to_owned())
    }
    omega_vec
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
