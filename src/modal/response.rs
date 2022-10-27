extern crate nalgebra as na;
use na::{DMatrix, DVector};
use std::ops::{AddAssign, Mul};

pub fn get_free_vibration(
    eigenvec_norm: &DMatrix<f64>,
    mass_mat: &DMatrix<f64>,
    init_cond: &DVector<f64>,
) -> DMatrix<f64> {
    // q init = EigenVec_trans * Mass_matrix * InitCond_Pos
    let q_0 = eigenvec_norm.transpose().mul(mass_mat).mul(init_cond);
    // So converting from modal space to physics space
    //  [x] = [[eigen_vec_norm]] * [q]
    //  Where [q] = [[q_0]] * [[cos(omega_i * t)]], each of matrix in right are diagonal
    //  By association:
    //  [x] = ( [[eigen_vec_norm]] * [[q_0]] ) * [[cos(omega_i * t)]]
    //          \--------------------------/
    //               this will not change
    //  [x] = [[constant]] * [[cos(omega_i * t)]]
    let init_mod_space_mat = DMatrix::from_diagonal(&q_0);
    let eigen_vec_mul_q_0 = eigenvec_norm.mul(init_mod_space_mat);
    eigen_vec_mul_q_0
}

pub fn solve_with_form_solution<F: Fn(f64, f64) -> f64>(
    time: f64,
    const_mat: &DMatrix<f64>,
    omegas: &DVector<f64>,
    form_solution: F,
) -> DVector<f64> {
    let mut sum: DVector<f64> = DVector::from_vec(vec![0.0; const_mat.nrows()]);
    for (i, col) in const_mat.column_iter().enumerate() {
        let cos_at_time = form_solution(time, omegas[i]);
        sum.add_assign(col.mul(cos_at_time));
    }
    sum
}

// Solution form : x= A cos (omega * t)
pub fn cos_sol_form(time: f64, omega: f64) -> f64 {
    (omega * time).cos()
}
