extern crate serde_yaml;
extern crate yaml_rust;
use std::process;
use crate::num_core::{
    smp_vector::SmpVector,
    n_body::NBody,
    io_config::{self, *},
};

mod num_core;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // read file
    let args: Vec<String> = std::env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        eprintln!("Problems with parsing command line arguments: {}", err);
        process::exit(1);
    });
    // read yaml
    let yaml_file = io_config::read_setting_yaml(config).unwrap();
    let n_body: NBody = serde_yaml::from_str(&yaml_file).unwrap();
    n_body.check().unwrap();

    println!("{:?}", n_body);
    
    Ok(())
}
