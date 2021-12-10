mod eigen;
mod tower;

extern crate nalgebra as na;
use std::ops::{DivAssign, Mul};

use na::DMatrix;

fn main() {
    // let tower = tower::load_tower().expect("error");
    // println!("Tower : {:?}", tower)
    let dim = 2;
    let m_se: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![1.0, 0.0, 0.0, 2.0]);
    let m_tow: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![1.0, 0.0, 0.0, 2.0]);
    let mut m_matrix: DMatrix<f64> = DMatrix::from_vec(dim * 2, dim * 2, vec![0.0; dim * dim * 4]);

    eigen::add_sub_matrix(&mut m_matrix, (0, 0), &m_se).expect("Error adding sub_matrix");
    eigen::add_sub_matrix(&mut m_matrix, (0, dim), &m_tow).expect("Error adding sub_matrix");

    let k_se: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![2.0, -1.0, -1.0, 2.0]);
    let k_tow: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![5.0, 6.0, 7.0, 8.0]);

    let mut k_matrix: DMatrix<f64> = DMatrix::from_vec(dim * 2, dim * 2, vec![0.0; dim * dim * 4]);

    eigen::add_sub_matrix(&mut k_matrix, (0, 0), &k_se).expect("Error adding sub_matrix");
    eigen::add_sub_matrix(&mut k_matrix, (0, dim), &k_tow).expect("Error adding sub_matrix");

    let inv_m_matrix = match m_se.clone().try_inverse() {
        Some(i) => i,
        None => panic!("Matrix could not be inversed"),
    };
    let a_matrix = inv_m_matrix * k_se;

    let eigen = a_matrix.eigenvalues().expect("Error getting eigenvalues");
    println!("A_matrix :{}", a_matrix);

    // let norm_a = a_matrix.clone().normalize();
    // println!("Norm_A_matrix :{}", norm_a);
    // a_matrix.mul_assign(DMatrix::from_vec(dim, dim, vec![1.0, 1.0, 2.0, 2.0]));
    // println!("A_matrix :{}", a_matrix);

    // let aaa = a_matrix.clone().symmetric_eigen();
    // println!("Eigen values :{}", aaa.eigenvalues);
    // println!("Eigen vectors :{}", aaa.eigenvectors);

    let mut eigen_vect: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]);
    // DMatrix::from_vec(dim * 2, dim * 2, vec![0.0; dim * dim * 4]);

    for (i, lambda) in eigen.iter().enumerate() {
        // println!("---------------------------------------------------------");
        // println!("Lambda : {}", lambda);
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec = a_matrix.clone() - &eigen_identity_matrix;
        // println!("ident : {}", eigen_identity_matrix);
        // println!("a_mat_for_eigen_vec : {}", a_mat_for_eigen_vec);
        // let x = a_mat_for_eigen_vec.symmetric_eigen();
        let mut b_vec = vec![0.0; dim];
        // Determinant is zero, which implies that there is only one independent equation for the two unknowns.
        // Hence, we can use either equation to solve for the ratios of the unknowns, therefore we
        // set b_vec[i]=1.0
        b_vec[i] = 1.0;
        let b = DMatrix::from_vec(dim, 1, b_vec);
        // println!("b : {}", b);
        let decomp = a_mat_for_eigen_vec.lu();
        let mut x = decomp.solve(&b).expect("Linear resolution failed.");
        x.div_assign(x[i]);
        // println!("Eigenvector : {}", x);
        eigen::add_sub_matrix(&mut eigen_vect, (0, i), &x).expect("Error adding sub_matrix");
    }
    println!("---------------------------------------------------------");

    let omega = eigen.map(|e| e.sqrt());
    println!("Omega : {}", omega);
    println!("EigenValues : {}", eigen);
    println!("EigenVectors : {}", eigen_vect);

    println!("---------------------------------------------------------");

    // Get modal mass and modal stiffness matrix
    let eigen_vect_trans = eigen_vect.transpose();

    let modal_mass = eigen_vect_trans
        .clone()
        .mul(m_se.clone())
        .mul(eigen_vect.clone());

    println!("m_se: {}", m_se);
    println!("{}", eigen_vect);
    println!("{}", eigen_vect_trans);
    println!("modal mass: {}", modal_mass);

    let mut eigen_vect_norm = eigen_vect.clone();
    // println!("{}", eigen_vect_norm[(0, 0)]);
    for (i, mut col) in eigen_vect_norm.column_iter_mut().enumerate() {
        col.div_assign(modal_mass[(i, i)].powf(0.5));
    }
    println!("eigenvectors norm : {}", eigen_vect_norm);
    let modal_mass_norm = eigen_vect_norm
        .transpose()
        .clone()
        .mul(m_se.clone())
        .mul(eigen_vect_norm.clone());
    println!("modal mass norm : {}", modal_mass_norm);
}
