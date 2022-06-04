extern crate nalgebra as na;
use std::ops::Mul;

use na::{DMatrix, DVector, Dynamic, SymmetricEigen};

#[derive(Debug)]
pub struct Modes {
    pub frequencies: DVector<f64>, //rad/s
    pub eigenvalues: DVector<f64>,
    pub eigenvectors_normalized: DMatrix<f64>,
    pub physic2modal: DMatrix<f64>,
    pub modal2physic: DMatrix<f64>,
}

impl Modes {
    fn report(&self) {
        println!("Frequencies : {}", self.frequencies);
        println!("Eigenvalues : {}", self.eigenvalues);
        println!(
            "Eigenvectors (normalized): {}",
            self.eigenvectors_normalized
        );
        println!("Physic2Modal: {}", self.physic2modal);
        println!("Modal2Physic (mode shapes): {}", self.modal2physic);
    }
}
pub fn eigen2(
    mass_matrix: &DMatrix<f64>,
    stiffness_matrix: &DMatrix<f64>,
) -> Result<Modes, String> {
    // let a =
    // mass_matrix.is
    // 1. Calculate M^-1/2 and M^1/2(assuming mass_matrix is diagonal)
    let diag_mass_sqroot = mass_matrix.diagonal().map(|e| e.powf(0.5));
    let mat_diag_mass_sqroot = DMatrix::from_diagonal(&diag_mass_sqroot);
    let inv_diag_mass_sqroot = mass_matrix.diagonal().map(|e| e.powf(-0.5));
    let mat_inv_diag_mass_sqroot = DMatrix::from_diagonal(&inv_diag_mass_sqroot);

    // 2. Calculate ~K = M^-1/2 * K * M^-1/2  ( mass normalized stiffness matrix )
    let kk = mat_inv_diag_mass_sqroot
        .clone()
        .mul(stiffness_matrix.mul(mat_inv_diag_mass_sqroot.clone()));

    // 3. Calculate the symmetric eigenvalue problem of ~K to get omega_i^2 and v_i
    let eigen = match SymmetricEigen::try_new(kk, 0.001, 100) {
        Some(i) => i,
        None => {
            return Err(String::from(
                "Symmetric eigenvalue problem could not be solved",
            ))
        }
    };
    // println!("Bef eigenvalues {}", eigen.eigenvalues);
    // println!("Bef eigenvectors {}", eigen.eigenvectors);
    // println!(
    //     "Bef Pt * P = I {}",
    //     eigen
    //         .eigenvectors
    //         .clone()
    //         .transpose()
    //         .mul(eigen.eigenvectors.clone())
    // );

    // 3.a.  Eigen formatting :
    // - Sort eigenvectors and eigenvalues by eigenvalues ascending
    // - If eigenvector has all values negative then change it to positive
    let (eigenvalues, mut eigenvectors) = sort_eigenvectors(eigen);
    prefer_all_positive_column(&mut eigenvectors);
    // println!(
    //     "After Pt * P = I {}",
    //     eigenvectors.clone().transpose().mul(eigenvectors.clone())
    // );

    // 4. Normalize v_i and form the matrix P
    // eigen.eigenvectors is already normalized in matrix form

    // 5. Calculate S = M^-0.5 * P
    let ss = mat_inv_diag_mass_sqroot.clone().mul(eigenvectors.clone());
    let ss_inv = eigenvectors.clone().transpose().mul(mat_diag_mass_sqroot);

    // println!("Verify S * S^-1 = I : {}", ss.clone().mul(ss_inv.clone()));
    // 6. Calculate the modal initial conditions: r(0)= S^-1 * x_0 , r_dot(0) = S^-1 * x_dot_0

    // let modal_init_cond_pos = ss_inv.clone().mul()

    let mut modes = Modes {
        frequencies: eigenvalues.map(|e| e.sqrt()),
        eigenvalues: eigenvalues,
        eigenvectors_normalized: eigenvectors,
        physic2modal: ss_inv,
        modal2physic: ss,
    };
    return Ok(modes);
}

fn prefer_all_positive_column(eigenvectors: &mut DMatrix<f64>) {
    for mut col in eigenvectors.column_iter_mut() {
        let any_positive = col.iter().any(|x| x > &0.0);
        if (!any_positive) {
            col.iter_mut().for_each(|x| *x *= -1.0);
        }
    }
}
fn sort_eigenvectors(eigen: SymmetricEigen<f64, Dynamic>) -> (DVector<f64>, DMatrix<f64>) {
    let mut eig_values = eigen.eigenvalues.clone();
    let mut eig_vectors = eigen.eigenvectors.clone();
    // Sort
    struct Mode {
        eigenvalue: f64,
        eigenvector: Vec<f64>,
    }

    let num_eigenvalues = eig_values.len();
    let mut vec_modes = vec![];
    for i in 0..num_eigenvalues {
        let newMode = Mode {
            eigenvalue: eig_values[i],
            eigenvector: eig_vectors.column(i).as_slice().to_vec(),
        };
        vec_modes.push(newMode);
    }
    vec_modes.sort_by(|a, b| a.eigenvalue.total_cmp(&b.eigenvalue));
    for i in 0..num_eigenvalues {
        eig_values[i] = vec_modes[i].eigenvalue;
        eig_vectors.set_column(i, &DVector::from_vec(vec_modes[i].eigenvector.clone()));
    }

    (eig_values, eig_vectors)
}

#[cfg(test)]
pub mod modal_tests {

    //This allows me to test private functions
    use super::*;
    use std::ops::Div;

    use na::DVector;

    // pub const EPS: f64 = 0.0001;
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
        // pub _2d_with_rigid_mode: Modal,
    }

    impl Setup {
        pub fn new() -> Self {
            let dim = 2;
            let mass_mat_vec: Vec<f64> = vec![9.0, 0.0, 0.0, 1.0]; //kg
            let stiff_mat_vec: Vec<f64> = vec![27.0, -3.0, -3.0, 3.0]; //N/m

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
            let _2d: Modal = Modal {
                mass_matrix: m_matrix,
                stiff_matrix: k_matrix,
                frequencies: sol_freq_vec_2d,
                eigenvalues: sol_eigenval_2d,
                eigenvectors_normalized: sol_eigenvec_norm_2d,
                free_transient: hand_calc_sol,
            };

            Self { _2d }
        }
    }

    #[test]
    fn eigen2_test() {
        let setup = Setup::new();
        let sol = eigen2(&setup._2d.mass_matrix, &setup._2d.stiff_matrix);
        assert!(sol.is_ok());
        let sol_data = sol.unwrap();
        sol_data.report();
    }
}
