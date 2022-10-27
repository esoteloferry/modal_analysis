#[cfg(test)]
pub mod test {
    extern crate nalgebra as na;
    use super::super::modes::*;
    use na::{DMatrix, DVector};

    pub const EPS: f64 = 0.0001;

    pub struct Modal {
        pub mass_matrix: DMatrix<f64>,
        pub stiff_matrix: DMatrix<f64>,
        pub frequencies: DVector<f64>,
        pub eigenvalues: DMatrix<f64>,
        pub eigenvectors_normalized: DMatrix<f64>,
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
                DMatrix::from_vec(dim, dim, vec![0.88807, -0.32506, 0.4597, 0.62796]); // X1 and X2 [m] (real)
            let _2d: Modal = Modal {
                mass_matrix: m_matrix,
                stiff_matrix: k_matrix,
                frequencies: sol_freq_vec_2d,
                eigenvalues: sol_eigenval_2d,
                eigenvectors_normalized: sol_eigenvec_norm_2d,
            };

            Self { _2d }
        }
    }

    #[test]
    fn basic() {
        let setup = Setup::new();
        let sol = eigen(&setup._2d.mass_matrix, &setup._2d.stiff_matrix);
        assert!(sol.is_ok());
        let sol_data = sol.unwrap();
        sol_data.report();
    }
}
