extern crate nalgebra as na;

pub mod modal;
use na::{DMatrix, DVector};

#[derive(Debug)]
pub struct Structure {
    pub gen_coord: Vec<SensorName>,
    pub mass: DMatrix<f64>,
    pub stiffness: DMatrix<f64>,
    pub modes: modal::Mode,
}

impl Structure {
    pub fn new(mass: DMatrix<f64>, stiffness: DMatrix<f64>) -> Result<Self, String> {
        let modes = match modal::eigen(&mass, &stiffness) {
            Ok(val) => val,
            Err(err) => return Err(format!("Error getting modes : {}", err)),
        };
        let mut sensor_names: Vec<SensorName> = vec![];
        for (i, _) in modes.frequency.iter().enumerate() {
            sensor_names.push(SensorName {
                name: format!("X{}", i.to_string()),
                units: "m".to_string(),
            })
        }

        Ok(Structure {
            gen_coord: sensor_names,
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
    pub sensors: Vec<Sensor>,
    pub init_cond_modal: DMatrix<f64>,
    time: f64,
    config: ConfigSim,
}

impl StructureSim {
    pub fn new(structure: Structure, config: ConfigSim) -> Self {
        let dim = structure.mass.nrows();
        let mut sensors = vec![];
        for coord in &structure.gen_coord {
            sensors.push(Sensor {
                info: coord.clone(),
                timeseries: vec![],
            })
        }
        // let sensors = structure
        //     .gen_coord
        //     .iter()
        //     .map(|n| Sensor {
        //         info: n.clone(),
        //         timeseries: vec![],
        //     })
        //     .collect();
        StructureSim {
            structure,
            sensors,
            init_cond_modal: DMatrix::from_vec(dim, dim, vec![0.0; dim * dim]),
            time: 0.0,
            config,
        }
    }
    pub fn set_initial_conditions(&mut self, init_position: &DVector<f64>) {
        let init_modal = modal::get_free_vibration(
            &self.structure.modes.eigenvectors_normalized,
            &self.structure.mass,
            init_position,
        );
        self.init_cond_modal = init_modal;
        //Reset timeseseries
        for sensor in self.sensors.iter_mut() {
            sensor.timeseries = vec![]
        }
        // for (sensor, ini) in self.sensors.iter_mut().zip(init_position) {
        //     sensor.timeseries = vec![*ini]
        // }
    }
    //Step will start from time =0.0 s
    fn step(&mut self) {
        let latest_pos_all = modal::solve_with_form_solution(
            self.time,
            &self.init_cond_modal,
            &self.structure.modes.frequency,
            modal::cos_sol_form,
        );
        for (i, latest_pos) in latest_pos_all.iter().enumerate() {
            self.sensors[i].timeseries.push(*latest_pos)
        }

        self.time += self.config.timestep;
    }
}

#[derive(Debug, Clone)]
pub struct ConfigSim {
    pub timestep: f64,
    pub total_time: f64,
}
const TIME_STEP: f64 = 0.1;
const TOTAL_TIME: f64 = 150.0;

#[derive(Debug, Clone)]
pub struct SensorName {
    pub name: String,
    pub units: String,
}

#[derive(Debug)]
pub struct Sensor {
    pub info: SensorName,
    pub timeseries: Vec<f64>,
}

#[derive(Debug)]
pub struct Simulation {
    pub time: Vec<f64>,
    pub config: ConfigSim,
    pub structure: StructureSim,
}

impl Simulation {
    pub fn new(structure: StructureSim, config: ConfigSim) -> Self {
        // let struc_sim = StructureSim::new(structure, &init_position);
        let time_vect = get_time_array(config.total_time, config.timestep);

        Simulation {
            time: time_vect,
            config,
            structure,
        }
    }

    pub fn run(&mut self) {
        for _ in 0..(self.time.len()) {
            self.structure.step();
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

fn find_time_pos_in_vec(time_vec: &Vec<f64>, time_to_find: f64) -> i32 {
    let (mut t_found, mut pos_time) = (false, 0);
    for (i, t) in time_vec.iter().enumerate() {
        if (t - time_to_find).abs() < 0.0001 {
            t_found = true;
            pos_time = i;
        }
    }

    if !t_found {
        return -1;
    }
    return pos_time as i32;
}
#[cfg(test)]
mod tests {
    use core::panic;

    use crate::get_time_array;

    use super::*;

    use modal::tests::Setup;

    #[test]
    fn simulation_step() {
        let setup = Setup::new();
        let struct_ = Structure::new(setup._2d.mass_matrix, setup._2d.stiff_matrix);
        assert!(struct_.is_ok());
        let struct_ = struct_.unwrap();
        let config = ConfigSim {
            timestep: 1.0,
            total_time: 10.0,
        };
        let mut struct_sim = StructureSim::new(struct_, config.clone());
        struct_sim.set_initial_conditions(&setup._2d.free_transient.init_cond);
        let mut sim = Simulation::new(struct_sim, config.clone());

        sim.run();

        // println!("time : {:?}", sim.time.len());

        // println!("position {:?}", sim.structure.sensors[0]);
        for trans in &setup._2d.free_transient.response {
            let pos_time = find_time_pos_in_vec(&sim.time, trans.time);
            if pos_time == -1 {
                panic!("Time from test not found in sim.time");
            }
            for (i, x_sens) in sim.structure.sensors.iter().enumerate() {
                let pos_from_sensor = x_sens.timeseries[pos_time as usize];
                assert!((pos_from_sensor - trans.values[i]).abs() < modal::tests::EPS)
            }
        }
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
