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

// fn get_num_timesteps(total_time: f64, timestep: f64) -> usize {
//     total_time % timestep
// }
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
    pub total_time: f64,
}
const TIME_STEP: f64 = 0.1;
const TOTAL_TIME: f64 = 150.0;

#[derive(Debug)]
pub struct Simulation {
    pub time: Vec<f64>,
    pub config: ConfigSim,
    pub position: Vec<DVector<f64>>, //this will be solution variables
    pub structure: StructureSim,
}

impl Simulation {
    pub fn new(structure: Structure, init_position: DVector<f64>) -> Simulation {
        let struc_sim = StructureSim::new(structure, &init_position);
        let time_vect = get_time_array(TOTAL_TIME, TIME_STEP);
        let mut position: Vec<DVector<f64>> = vec![];
        position.push(init_position);

        Simulation {
            time: time_vect,
            config: ConfigSim {
                timestep: TIME_STEP,
                total_time: TOTAL_TIME,
            },
            position: position,
            structure: struc_sim,
        }
    }

    pub fn run(&mut self) {
        for time in self.time.iter() {
            let new_position = self.structure.step(time.clone());
            self.position.push(new_position);
        }
    }
}
fn get_time_array(total_time: f64, time_step: f64) -> Vec<f64> {
    let num_steps = (total_time / time_step).round() as usize;
    let mut time = vec![0.0; num_steps + 1];
    let mut current_time = 0.0;
    for i in 0..(num_steps + 1) {
        time[i] = if current_time > total_time {
            total_time
        } else {
            current_time
        };
        current_time += time_step;
    }
    return time;
}

#[cfg(test)]
mod tests {
    use crate::get_time_array;

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

        sim.run();

        println!("time : {:?}", sim.time);
        // for pos in sim.position {
        //     println!("sim position : {}", pos);
        // }
    }

    #[test]
    fn test_get_time_array() {
        let time_vec1 = get_time_array(10.0, 1.0);
        assert!(time_vec1.len().eq(&11));
        assert!(time_vec1[time_vec1.len() - 1].eq(&10.0));

        let time_vec2 = get_time_array(10.0, 20.0);
        assert!(time_vec2.len().eq(&2));
        assert!(time_vec2[time_vec2.len() - 1].eq(&10.0));

        let time_vec3 = get_time_array(10.0, 3.3);
        assert!(time_vec3.len().eq(&4));
        assert!(time_vec3[time_vec3.len() - 1].lt(&10.0));
    }
}
