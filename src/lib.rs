extern crate nalgebra as na;

mod modal;
use na::{DMatrix, DVector};

#[derive(Debug)]
pub struct Structure {
    pub mass: DMatrix<f64>,
    pub stiffness: DMatrix<f64>,
    pub modes: modal::Mode,
}

impl Structure {
    fn new(mass: DMatrix<f64>, stiffness: DMatrix<f64>) -> Result<Structure, String> {
        let modes = match modal::eigen(&mass, &stiffness) {
            Ok(val) => val,
            Err(err) => return Err(format!("Error getting modes : {}", err)),
        };

        Ok(Structure {
            mass,
            stiffness,
            modes,
        })
    }

    //TODO:_
    fn get_total_mass() {}

    fn get_frequencies() {}
}

#[derive(Debug)]
pub struct StructureSim {
    pub structure: Structure,
    pub init_cond_modal: DMatrix<f64>,
}

impl StructureSim {
    fn new(structure: Structure, init_position: &DVector<f64>) -> StructureSim {
        let init_modal = modal::get_free_vibration(
            &structure.modes.eigenvectors_normalized,
            &structure.mass,
            init_position,
        );
        StructureSim {
            structure: structure,
            init_cond_modal: init_modal,
        }
    }
    fn step(&self, time: f64) -> DVector<f64> {
        modal::solve_with_form_solution(
            time,
            &self.init_cond_modal,
            &self.structure.modes.frequency,
            modal::cos_sol_form,
        )
    }
}
// fn set_initial_conditions(
//     structure: &Structure,
//     init_position: &DVector<f64>,
// ) -> impl Fn(f64) -> DVector<f64> {
//     let free_vib = modal::get_free_vibration(
//         &structure.modes.eigenvectors_normalized,
//         &structure.modes.frequency,
//         &structure.mass,
//         &init_position,
//     )
//     .unwrap();
//     free_vib
// }
#[derive(Debug)]
pub struct ConfigSim {
    pub timestep: f64,
}

#[derive(Debug)]
pub struct Simulation {
    pub time: f64,
    pub config: ConfigSim,
    pub position: Vec<DVector<f64>>, //this will be solution variables
    pub structure: StructureSim,
}

impl Simulation {
    fn new(structure: Structure, init_position: DVector<f64>) -> Simulation {
        let struc_sim = StructureSim::new(structure, &init_position);
        let mut position: Vec<DVector<f64>> = Vec::with_capacity(100);
        position.push(init_position);

        Simulation {
            time: 0.0,
            config: ConfigSim { timestep: 0.05 },
            position: position,
            structure: struc_sim,
        }
    }

    fn step(&mut self) {
        self.time += self.config.timestep;
        self.position.push(self.structure.step(self.time));
    }
}

#[cfg(test)]
mod struct_sim_tests {
    use super::*;

    use modal::tests::Setup;

    #[test]
    fn simulation_step() {
        let setup = Setup::new();
        let struct_ = Structure::new(setup._2d.mass_matrix, setup._2d.stiff_matrix);
        assert!(struct_.is_ok());
        let struct_ = struct_.unwrap();
        let init_position = DVector::from_vec(vec![1.0, 0.0]);
        let mut sim = Simulation::new(struct_, init_position);

        sim.step();
        sim.step();

        println!("sim position : {:?}", sim.position);
    }
}
