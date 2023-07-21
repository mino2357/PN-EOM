use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(rename_all = "PascalCase")]
pub struct NBody {
    pub setting_name: String,
    pub number_of_bodies: usize,
    pub mass: Vec<f64>,
    pub position: Vec<Vec<f64>>,
}

impl NBody {
    pub fn check(&self) -> Result<&NBody, &'static str> {
        if self.mass.len() != self.number_of_bodies {
            return Err("mass.len() != number_of_bodies");
        }
        if self.position.len() != self.number_of_bodies {
            return Err("position.len() != number_of_bodies");
        }

        Ok(self)
    }
}