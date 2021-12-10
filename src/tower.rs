use std::{fs::File, io, path::Path};

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Tower {
    density: f32,
    emodulus: f32,
    sections: Vec<Section>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Section {
    h: f32,
    t: f32,
    d: f32,
    mlump: f32,
}

pub fn load_tower() -> Result<Tower, io::Error> {
    let json_file_path = Path::new("./src/tower.gf.json");
    // let file = File::open(json_file_path).expect("error opening file");
    let file = File::open(json_file_path)?;
    // let tower: Tower = serde_json::from_reader(file).expect("error while reading or parsing");
    let tower: Tower = serde_json::from_reader(file)?;
    // println!("tower = {:?}", tower);
    Ok(tower)
}
