extern crate nalgebra as na;
use std::ops::{AddAssign, DivAssign, Mul};

use na::{DMatrix, DVector};

#[derive(Debug)]
pub struct Mode {
    pub frequency: DVector<f64>,               //rad/s
    pub eigenvalues: DMatrix<f64>,             // Also known as Spectral Matrix
    pub eigenvectors_normalized: DMatrix<f64>, //Also known as Modal Matrix
    pub eigenvectors: DMatrix<f64>,
}
//In general, it is assumed that the diagonal entries of the spectral
//matrix and the columns in the modal matrix are ordered according to the
//increasing magnitudes of the eigenvalues BUT in this case is not returned in that way

// Properties:

// Orthonality property
// For an eigenvalue (lambda_i) and eigenvector (x_i), the eigenvalue problem equation [K]*x = lambda*[M]*x can be
// rewriten as [K]*x_a = lambda_a*[M]*x_a , another pair can be [K]*x_b = lambda_b*[M]*x_b
// The property makes:
//  x_b_trans * [M] * x_a = 0
//  x_b_trans * [K] * x_a = 0
pub fn check_orthogonality_property(
    mass: DMatrix<f64>,
    _stiff: DMatrix<f64>,
    eigenvectors: DMatrix<f64>,
) -> bool {
    let dim = eigenvectors.ncols();
    let vec_dim_num: Vec<usize> = (0..dim).collect();

    for (i, col) in eigenvectors.column_iter().enumerate() {
        let mut other_cols = vec_dim_num.clone();
        other_cols.retain(|&x| !x.eq(&i));
        for other_col_idx in other_cols {
            let mul_vecb_mass_veca = col
                .transpose()
                .mul(&mass)
                .mul(eigenvectors.column(other_col_idx));
            // println!("Check {}", mul_vecb_mass_veca);
            if mul_vecb_mass_veca.sum() != 0.0 {
                return false;
            }
        }
    }
    true
}
// Expansion theorem or Principle of modal superposition

// Lib functions
fn print_type_of<T>(_: T) {
    println!("{}", std::any::type_name::<T>())
}
pub fn eigen(mass_matrix: &DMatrix<f64>, stiffness_matrix: &DMatrix<f64>) -> Result<Mode, String> {
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

    let inv_m_matrix = match mass_matrix.clone().try_inverse() {
        Some(i) => i,
        None => return Err(String::from("Matrix could not be inversed")),
    };
    let a_matrix = inv_m_matrix * stiffness_matrix;

    let eigen = match a_matrix.eigenvalues() {
        Some(i) => i,
        None => return Err(String::from("Error getting eigenvalues")),
    };
    println!("A_matrix :{}", a_matrix);

    let mut eigen_vect: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]);

    // println!("HEllloo {}", eigen);
    for (i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec: DMatrix<f64> = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec: Vec<f64> = vec![0.0; dim];
        // Determinant is zero, which implies that there is only one independent equation for the two unknowns.
        // Hence, we can use either equation to solve for the ratios of the unknowns, therefore we
        // set b_vec[i]=1.0
        b_vec[i] = 1.0;
        // Other way can be to set always the first value to 1.0
        // b_vec[0] = 1.0;
        // println!("HEllloo {}", lambda);
        let b: DMatrix<f64> = DMatrix::from_vec(dim, 1, b_vec);
        // let b: DVector<f64> = DVector::from_vec(b_vec);
        println!("helllooa aftehet {}", a_mat_for_eigen_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        let mut x = match decomp.solve(&b) {
            Some(i) => i,
            None => return Err(String::from("Linear resolution failed")),
        };
        x.div_assign(x[i]);
        // This way corresponds to setting first value to 1.0 in b_vec
        // x.div_assign(x[0]);
        match add_sub_matrix(&mut eigen_vect, (0, i), &x) {
            Ok(_) => (),
            Err(e) => return Err(e),
        };
    }

    let omega = eigen.clone().map(|e| e.sqrt());

    // Get modal mass and modal stiffness matrix
    let modal_mass = eigen_vect.transpose().mul(mass_matrix).mul(&eigen_vect);

    let mut eigen_vect_norm = eigen_vect.clone();
    for (i, mut col) in eigen_vect_norm.column_iter_mut().enumerate() {
        col.div_assign(modal_mass[(i, i)].powf(0.5));
    }

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

//Functions to get vibration motion
//
// x = x_free  +  x_forced
//
//
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

// pub fn get_free_vibration(
//     eigenvec_norm: &DMatrix<f64>,
//     omegas: &Vec<f64>,
//     mass_mat: &DMatrix<f64>,
//     init_cond: &DVector<f64>,
// ) -> Result<impl Fn(f64) -> DVector<f64>, String> {
//     // q init = EigenVec_trans * Mass_matrix * InitCond_Pos
//     let q_0 = eigenvec_norm.transpose().mul(mass_mat).mul(init_cond);
//     // So converting from modal space to physics space
//     //  [x] = [[eigen_vec_norm]] * [q]
//     //  Where [q] = [[q_0]] * [[cos(omega_i * t)]], each of matrix in right are diagonal
//     //  By association:
//     //  [x] = ( [[eigen_vec_norm]] * [[q_0]] ) * [[cos(omega_i * t)]]
//     //          \--------------------------/
//     //               this will not change
//     //  [x] = [[constant]] * [[cos(omega_i * t)]]
//     let init_mod_space_mat = DMatrix::from_diagonal(&q_0);
//     let eigen_vec_mul_q_0 = eigenvec_norm.mul(init_mod_space_mat);
//     Ok(solve_with_form_solution(
//         eigen_vec_mul_q_0,
//         omegas.to_owned(),
//         cos_sol_form,
//     ))
// }
// Solution form : x= A cos (omega * t)
pub fn cos_sol_form(time: f64, omega: f64) -> f64 {
    (omega * time).cos()
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
// fn solve_with_form_solution<F: Fn(f64, f64) -> f64>(
//     const_mat: DMatrix<f64>,
//     omegas: Vec<f64>,
//     form_solution: F,
// ) -> impl Fn(f64) -> DVector<f64> {
//     move |time| {
//         let mut sum: DVector<f64> = DVector::from_vec(vec![0.0; const_mat.nrows()]);
//         for (i, col) in const_mat.column_iter().enumerate() {
//             let cos_at_time = form_solution(time, omegas[i]);
//             sum.add_assign(col.mul(cos_at_time));
//         }
//         sum
//     }
// }

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
pub mod tests {

    //This allows me to test private functions
    use super::*;
    use std::ops::Div;

    use na::{DVector, Matrix2x4, Vector2};

    pub const EPS: f64 = 0.0001;
    pub struct Modal {
        pub mass_matrix: DMatrix<f64>,
        pub stiff_matrix: DMatrix<f64>,
        frequencies: DVector<f64>,
        eigenvalues: DMatrix<f64>,
        eigenvectors_normalized: DMatrix<f64>,
        pub free_transient: Transient,
    }

    pub struct Transient {
        pub init_cond: DVector<f64>,
        pub response: Vec<Response>,
    }
    pub struct Response {
        pub time: f64,
        pub values: DVector<f64>,
    }

    pub struct Setup {
        pub _2d: Modal,
    }

    impl Setup {
        pub fn new() -> Self {
            let dim = 2;
            let mass_mat_vec: Vec<f64> = vec![1.0, 0.0, 0.0, 2.0]; //kg
            let stiff_mat_vec: Vec<f64> = vec![2.0, -1.0, -1.0, 2.0]; //N/m

            let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
            let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);
            //Solution
            let (w1_2d, w2_2d) = (1.53819, 0.79623); //rad/s
            let sol_freq_vec_2d = DVector::from_vec(vec![w1_2d, w2_2d]);
            let sol_eigenval_2d =
                DMatrix::from_vec(dim, dim, vec![w1_2d * w1_2d, 0.0, 0.0, w2_2d * w2_2d]);
            let sol_eigenvec_norm_2d =
                DMatrix::from_vec(dim, dim, vec![0.88807, -0.32506, 0.4597, 0.62796]);
            let response_t5s = Response {
                time: 5.0,
                values: DVector::from_vec(vec![-0.0131029, -0.2396485]).div(100.0),
            };
            let response_t10s = Response {
                time: 10.0,
                values: DVector::from_vec(vec![-0.7700099, 0.24229929]).div(100.0),
            };
            let hand_calc_sol = Transient {
                init_cond: DVector::from_vec(vec![1.0, 0.0]).div(100.0),
                response: vec![response_t5s, response_t10s],
            }; // X1 and X2 [m] (real)
            Self {
                _2d: Modal {
                    mass_matrix: m_matrix,
                    stiff_matrix: k_matrix,
                    frequencies: sol_freq_vec_2d,
                    eigenvalues: sol_eigenval_2d,
                    eigenvectors_normalized: sol_eigenvec_norm_2d,
                    free_transient: hand_calc_sol,
                },
            }
        }
    }

    #[test]
    fn eigen_analysis() {
        let setup = Setup::new();
        let m_matrix = setup._2d.mass_matrix;
        let k_matrix = setup._2d.stiff_matrix;
        let modes = eigen(&m_matrix, &k_matrix);
        assert!(modes.is_ok());
        let modes = modes.unwrap();
        //Print
        println!("Frequencies : {:?}", modes.frequency);
        println!("Eigenvalues : {}", modes.eigenvalues);
        println!("Eigenvectors : {}", modes.eigenvectors);
        println!("Eigenvectors normalized: {}", modes.eigenvectors_normalized);

        // println!("{:?}", modes.eigenvectors.clone().mul(mass_mat_vec));
        // assert_eq!(2 + 2, 4);
        let freq_vec_na = modes.frequency;
        assert!(freq_vec_na.relative_eq(&setup._2d.frequencies, EPS, EPS));
        assert!(modes
            .eigenvalues
            .relative_eq(&setup._2d.eigenvalues, EPS, EPS));
        assert!(modes.eigenvectors_normalized.relative_eq(
            &setup._2d.eigenvectors_normalized,
            EPS,
            EPS
        ));

        println!(
            "Xt * M * X : {}",
            modes
                .eigenvectors_normalized
                .transpose()
                .mul(&m_matrix)
                .mul(&modes.eigenvectors_normalized)
        );
        let generalized_mass = modes
            .eigenvectors
            .transpose()
            .mul(&m_matrix)
            .mul(&modes.eigenvectors);
        println!("Xt * M * X : {}", generalized_mass);
        let check_orthog_prop =
            check_orthogonality_property(m_matrix, k_matrix, modes.eigenvectors_normalized);
        assert!(check_orthog_prop);
    }
    #[test]
    fn test_case_runs_ok() {
        let mass_mat_vec: Vec<f64> = vec![2.0, 0.0, 0.0, 1.0];
        let stiff_mat_vec: Vec<f64> = vec![4000.0, -3000.0, -3000.0, 5000.0];

        let dim: usize = 2;
        let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
        let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);

        let modes = eigen(&m_matrix, &k_matrix);
        assert!(modes.is_ok());
        // let modes = modes.unwrap();
        //Print
        // println!("Frequencies : {:?}", modes.frequency);
        // println!("Eigenvalues : {}", modes.eigenvalues);
        // println!("Eigenvectors : {}", modes.eigenvectors);
        // println!("Eigenvectors normalized: {}", modes.eigenvectors_normalized);
    }
    #[test]
    fn rigid_body() {
        let mass_mat_vec: Vec<f64> = vec![150.0, 0.0, 0.0, 100.0];
        let stiff_mat_vec: Vec<f64> = vec![15000.0, -15000.0, -15000.0, 15000.0];

        let dim: usize = 2;
        let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
        let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);

        let modes = eigen(&m_matrix, &k_matrix);
        assert!(modes.is_err());
        // TODO: check eigenvectors and eigenvalues of elastic modes
        println!("{}", modes.expect_err(""))
    }

    #[test]
    fn test_get_free_vibration() {
        let setup = Setup::new();
        let init_cond_modal = get_free_vibration(
            &setup._2d.eigenvectors_normalized,
            &setup._2d.mass_matrix,
            &setup._2d.free_transient.init_cond,
        );

        for trans in &setup._2d.free_transient.response {
            let sol = solve_with_form_solution(
                trans.time,
                &init_cond_modal,
                &setup._2d.frequencies,
                cos_sol_form,
            );
            assert!(trans.values.relative_eq(&sol, EPS, EPS));
        }
        // let sol_t5 =
        //     solve_with_form_solution(5.0, &init_cond_modal, &setup._2d.frequencies, cos_sol_form);
        // let sol_t10 =
        //     solve_with_form_solution(10.0, &init_cond_modal, &setup._2d.frequencies, cos_sol_form);
        // // assert!(fnvib.is_ok());
        // // let fnvib = fnvib.unwrap();
        //
        // println!("t5: {}", sol_t5);
        // println!("t10: {}", sol_t10);
        // assert!(sol_t5.relative_eq(&setup._2d.free_transient.response_t5s, eps, eps));
        // assert!(sol_t10.relative_eq(&setup._2d.free_transient.response_t10s, eps, eps));

        // println!("{}", fnvib(0.0));
        // println!("{}", fnvib(2.0));
        // println!("{}", fnvib(4.0));,
    }
    // #[test]
    // fn vec_mul_mat() {
    //     let vec: DVector<f64> = DVector::from_vec(vec![1.0, 1.0]);
    //     let matrix: DMatrix<f64> = DMatrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, 2.0]);
    //     //     println!("{}", vec.transpose());
    //     //     println!("{}", matrix);
    //     //     println!("{}", vec.transpose().mul(matrix));
    // }
    // #[test]
    // fn vec_mul_mat_compile() {
    //     let vec: Vector2<f64> = Vector2::from_vec(vec![1.0, 1.0]);
    //     let matrix: Matrix2x4<f64> =
    //         Matrix2x4::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    //     // println!("{}", vec.mul(matrix))
    //     // println!("{}", vec.transpose());
    //     // println!("{}", matrix);
    //     // println!("{}", vec.transpose().mul(matrix));
    // }
    //
    // #[test]
    // fn vec_add() {
    //     let mut vec_empty: Vector2<f64> = Vector2::from_vec(vec![0.0; 2]);
    //     let vec_to_add: Vector2<f64> = Vector2::from_vec(vec![1.0, 1.0]);
    //     vec_empty.add_assign(vec_to_add);
    //     println!("{}", vec_empty);
    // }
}
