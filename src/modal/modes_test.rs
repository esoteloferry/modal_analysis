#[cfg(test)]
pub mod test {
    extern crate nalgebra as na;
    use std::ops::MulAssign;

    use super::super::modes::*;
    use na::{DMatrix, DVector};

    pub const EPS: f64 = 0.0001;

    pub struct Modal {
        pub mass_matrix: DMatrix<f64>,
        pub stiff_matrix: DMatrix<f64>,
        pub frequencies: DVector<f64>,
        pub eigenvalues: DVector<f64>,
        pub eigenvectors_normalized: DMatrix<f64>,
    }

    pub struct Setup {
        pub _2d: Modal,
        pub _2d_with_rigid_mode: Modal,
    }

    impl Setup {
        pub fn new() -> Self {
            let dim = 2;

            // ***** Rigid modes *****
            // Example 4.6.1 - Engineering vibration 2nd edition - Daniel J Inman
            let mass_mat_vec: Vec<f64> = vec![9.0, 0.0, 0.0, 1.0]; //kg
            let stiff_mat_vec: Vec<f64> = vec![27.0, -3.0, -3.0, 3.0]; //N/m

            let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
            let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);
            //Solution
            let lambda1_2d: f64 = 2.0;
            let lambda2_2d: f64 = 4.0;
            let sol_freq_vec_2d = DVector::from_vec(vec![lambda1_2d.sqrt(), lambda2_2d.sqrt()]);
            let sol_eigenval_2d = DVector::from_vec(vec![lambda1_2d, lambda2_2d]);
            let mut sol_eigenvec_norm_2d = DMatrix::from_vec(dim, dim, vec![1.0, 1.0, 1.0, -1.0]); // X1 and X2 [m] (real)
            sol_eigenvec_norm_2d.mul_assign(0.7071);
            let _2d: Modal = Modal {
                mass_matrix: m_matrix,
                stiff_matrix: k_matrix,
                frequencies: sol_freq_vec_2d,
                eigenvalues: sol_eigenval_2d,
                eigenvectors_normalized: sol_eigenvec_norm_2d,
            };
            // ***** Rigid modes *****
            // Example 4.4.4 - Engineering vibration 2nd edition - Daniel J Inman
            let mass_mat_vec: Vec<f64> = vec![1.0, 0.0, 0.0, 4.0]; //kg
            let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
            let mut k_matrix: DMatrix<f64> =
                DMatrix::from_vec(dim, dim, vec![1.0, -1.0, -1.0, 1.0]);
            k_matrix.mul_assign(400.0);
            // Solution
            let lambda1_2d: f64 = 0.0;
            let lambda2_2d: f64 = 5.0;
            let sol_freq_vec_2d = DVector::from_vec(vec![lambda1_2d.sqrt(), lambda2_2d.sqrt()]);
            let sol_eigenval_2d = DVector::from_vec(vec![lambda1_2d, lambda2_2d]);
            let mut sol_eigenvec_norm_2d = DMatrix::from_vec(dim, dim, vec![1.0, 2.0, -2.0, 1.0]);
            sol_eigenvec_norm_2d.mul_assign(0.4472);
            let _2d_with_rigid_mode: Modal = Modal {
                mass_matrix: m_matrix,
                stiff_matrix: k_matrix,
                frequencies: sol_freq_vec_2d,
                eigenvalues: sol_eigenval_2d,
                eigenvectors_normalized: sol_eigenvec_norm_2d,
            };

            //
            Self {
                _2d,
                _2d_with_rigid_mode,
            }
        }
    }
    #[test]
    fn eigen_2d() {
        let setup = Setup::new();
        let sol = eigen(&setup._2d.mass_matrix, &setup._2d.stiff_matrix);
        assert!(sol.is_ok());
        let sol_data = sol.unwrap();

        assert!(sol_data
            .frequencies
            .relative_eq(&setup._2d.frequencies, EPS, EPS));
        assert!(sol_data
            .eigenvalues
            .relative_eq(&setup._2d.eigenvalues, EPS, EPS));
        assert!(sol_data.eigenvectors_normalized.relative_eq(
            &setup._2d.eigenvectors_normalized,
            EPS,
            EPS
        ));
    }

    #[test]
    fn basic() {
        let setup = Setup::new();
        let sol = eigen(&setup._2d.mass_matrix, &setup._2d.stiff_matrix);
        assert!(sol.is_ok());
        let sol_data = sol.unwrap();
        sol_data.report();
    }
    // #[test]
    // fn matrix_order() {
    //     let mut mat = DMatrix::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0]);
    //     mat.mul_assign(2.0);
    //     println!("Hello {}", mat);
    // }
}
