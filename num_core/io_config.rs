//
// Around IO and YAML parsing
//
use std::error::Error;
use std::fs::File;
use std::io::Read;

// The following code is based on Chapter 12 of The Book.
// https://doc.rust-jp.rs/book-ja/ch12-03-improving-error-handling-and-modularity.html
pub struct Config {
    pub build_command: String,
    pub file_name: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 2 {
            return Err("Not enough arguments");
        }

        let build_command = args[0].clone();
        let file_name = args[1].clone();

        Ok(Config {
            build_command,
            file_name,
        })
    }
}

// Extract the contents of the setting file from config.
pub fn read_setting_yaml(config: Config) -> Result<String, Box<dyn Error>> {
    let mut f = File::open(config.file_name)?;

    let mut contents = String::new();
    f.read_to_string(&mut contents)?;

    Ok(contents)
}
