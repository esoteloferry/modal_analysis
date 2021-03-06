mod eigen;
mod tower;

extern crate nalgebra as na;
use std::ops::{AddAssign, DivAssign, Mul, MulAssign};

use na::DMatrix;

fn main() {
    let dim = 2;
    let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![1.0, 0.0, 0.0, 2.0]);

    let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![2.0, -1.0, -1.0, 2.0]);

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
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(2 * dim, 2 * dim).mul(*lambda);
        let a_mat_for_eigen_vec = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec = vec![0.0; dim * 2];
        // Determinant is zero, which implies that there is only one independent equation for the two unknowns.
        // Hence, we can use either equation to solve for the ratios of the unknowns, therefore we
        // set b_vec[i]=1.0
        b_vec[i] = 1.0;
        let b = DMatrix::from_vec(dim * 2, 1, b_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        let mut x = decomp.solve(&b).expect("Linear resolution failed.");
        x.div_assign(x[i]);
        eigen::add_sub_matrix(&mut eigen_vect, (0, i), &x).expect("Error adding sub_matrix");
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
    // println!("eigenvectors norm : {}", eigen_vect_norm);
    // let modal_mass_norm = eigen_vect_norm
    //     .transpose()
    //     .clone()
    //     .mul(m_se.clone())
    //     .mul(eigen_vect_norm.clone());
    // println!("modal mass norm : {}", modal_mass_norm);
    //
    // let diag_c = eigen.clone().div_assign(2.0);
    // println!("Diag_c : {}", diag_c);
    let mut c_natural = DMatrix::from_vec(dim * 2, dim * 2, vec![0.0; dim * dim * 4]);
    for (i, _) in c_natural.diagonal().iter_mut().enumerate() {
        c_natural[(i, i)] = 2.0 * 0.0025 * omega[(i)];
    }
    println!("{}", c_natural);

    let inv_eigen_vect_trans = match eigen_vect_trans.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Matrix could not be inversed"),
    };
    let inv_eigen_vect = match eigen_vect.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Matrix could not be inversed"),
    };

    let c_damping = inv_eigen_vect_trans.mul(c_natural).mul(inv_eigen_vect);

    println!("C damping : {}", c_damping);
}
