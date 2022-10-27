extern crate nalgebra as na;
use na::{DMatrix, DVector, Dynamic, SymmetricEigen};
use std::ops::Mul;

#[derive(Debug)]
pub struct Modes {
    pub frequencies: DVector<f64>, //rad/s
    pub eigenvalues: DVector<f64>,
    pub eigenvectors_normalized: DMatrix<f64>,
    pub physic2modal: DMatrix<f64>,
    pub modal2physic: DMatrix<f64>,
}

impl Modes {
    pub fn report(&self) {
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
pub fn eigen(mass_matrix: &DMatrix<f64>, stiffness_matrix: &DMatrix<f64>) -> Result<Modes, String> {
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

    let modes = Modes {
        frequencies: eigenvalues.map(|e| e.sqrt()),
        eigenvalues,
        eigenvectors_normalized: eigenvectors,
        physic2modal: ss_inv,
        modal2physic: ss,
    };
    return Ok(modes);
}

fn prefer_all_positive_column(eigenvectors: &mut DMatrix<f64>) {
    for mut col in eigenvectors.column_iter_mut() {
        let any_positive = col.iter().any(|x| x > &0.0);
        if !any_positive {
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
        let new_mode = Mode {
            eigenvalue: eig_values[i],
            eigenvector: eig_vectors.column(i).as_slice().to_vec(),
        };
        vec_modes.push(new_mode);
    }
    vec_modes.sort_by(|a, b| a.eigenvalue.total_cmp(&b.eigenvalue));
    for i in 0..num_eigenvalues {
        eig_values[i] = vec_modes[i].eigenvalue;
        eig_vectors.set_column(i, &DVector::from_vec(vec_modes[i].eigenvector.clone()));
    }

    (eig_values, eig_vectors)
}
