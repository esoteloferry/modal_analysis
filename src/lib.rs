extern crate nalgebra as na;
use std::ops::{DivAssign, Mul};

use na::DMatrix;

#[derive(Debug)]
pub struct Mode {
    pub frequency: Vec<f64>,
    pub eigenvalues: DMatrix<f64>, // Also known as Spectral Matrix
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
    stiff: DMatrix<f64>,
    eigenvectors: DMatrix<f64>,
) -> bool {
    let dim = eigenvectors.ncols();
    let vec_dim_num: Vec<usize> = (0..dim).collect();

    println!("Inputs: {}{}{}", mass, stiff, eigenvectors);
    for (i, col) in eigenvectors.column_iter().enumerate() {
        // col.mul_to(&mass, out)
        // let other_cols: Vec<usize> = vec_dim_num.iter().clone().filter(|&x| i.eq(x)).collect();
        let mut other_cols = vec_dim_num.clone();
        other_cols.retain(|&x| !x.eq(&i));
        for other_col_idx in other_cols {
            println!("{:?} {:?}", i, other_col_idx);
            println!(
                "{}{}{}",
                col.transpose(),
                mass,
                eigenvectors.column(other_col_idx)
            );
            println!(
                "Check {}",
                col.transpose()
                    .mul(&mass)
                    .mul(eigenvectors.column(other_col_idx))
            )
        }
        // let other_cols = vec_dim_num.clone().retain(|&x| x != i);
        // println!("HEHETEETETEE{} {}", col.transpose().mul(&mass), i);
        // println!("ASEHNNTE{}, {:?}, {:?}", i, vec_dim_num, other_cols);
    }
    true
}
// Expansion theorem or Principle of modal superposition

// Lib functions

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

    // println!("M_matrix : {}", mass_matrix);
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
    // println!("A_matrix :{}", a_matrix);

    let mut eigen_vect: DMatrix<f64> = DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]);

    for (i, lambda) in eigen.iter().enumerate() {
        let eigen_identity_matrix: DMatrix<f64> = DMatrix::identity(dim, dim).mul(*lambda);
        let a_mat_for_eigen_vec = a_matrix.clone() - &eigen_identity_matrix;
        let mut b_vec = vec![0.0; dim];
        // Determinant is zero, which implies that there is only one independent equation for the two unknowns.
        // Hence, we can use either equation to solve for the ratios of the unknowns, therefore we
        // set b_vec[i]=1.0
        b_vec[i] = 1.0;
        // Other way can be to set always the first value to 1.0
        // b_vec[0] = 1.0;
        let b = DMatrix::from_vec(dim, 1, b_vec);
        let decomp = a_mat_for_eigen_vec.lu();
        // let mut x = decomp.solve(&b).expect("Linear resolution failed.");
        let mut x = match decomp.solve(&b) {
            Some(i) => i,
            None => return Err(String::from("Linear resolution failed")),
        };
        x.div_assign(x[i]);
        // This way corresponds to setting first value to 1.0 in b_vec
        // x.div_assign(x[0]);
        // add_sub_matrix(&mut eigen_vect, (0, i), &x).expect("Error adding sub_matrix");
        match add_sub_matrix(&mut eigen_vect, (0, i), &x) {
            Ok(_) => (),
            Err(e) => return Err(e),
        };
    }
    // println!("---------------------------------------------------------");

    let omega = eigen.clone().map(|e| e.sqrt()).as_slice().to_vec();
    // let om2 = eigen.clone().map(|e| e.sqrt());
    // println!("Omega : {:?}", omega);
    // println!("EigenValues : {}", eigen);
    // println!("EigenVectors : {}", eigen_vect);

    // println!("---------------------------------------------------------");

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
    // let modal_mass_norm = eigen_vect_norm
    //     .transpose()
    //     .clone()
    //     .mul(mass_matrix.clone())
    //     .mul(eigen_vect_norm.clone());
    // println!("modal mass norm : {}", modal_mass_norm);
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

    use na::{DVector, Matrix2x4, Vector2};

    use super::*;
    #[test]
    fn it_works() {
        let mass_mat_vec: Vec<f64> = vec![1.0, 0.0, 0.0, 2.0];
        let stiff_mat_vec: Vec<f64> = vec![2.0, -1.0, -1.0, 2.0];

        let dim: usize = 2;
        let m_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, mass_mat_vec);
        let k_matrix: DMatrix<f64> = DMatrix::from_vec(dim, dim, stiff_mat_vec);

        let modes = eigen(&m_matrix, &k_matrix);
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
        // check_orthogonality_property(m_matrix, k_matrix, modes.eigenvectors_normalized);
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
}
