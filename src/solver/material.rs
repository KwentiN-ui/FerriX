//! This module contains data structures for material definitions. A Material is defined using the `*MATERIAL` card
//! globally. To use a material, it has to be referenced in a section definition.
#![allow(unused)]

use crate::solver::inp::InpFile;

#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub density: Option<f64>,
    pub elastic: Option<(f64, f64)>,
}

impl Material {
    pub fn from_input(input: &InpFile) -> Vec<Self> {
        let mut materials = Vec::new();

        for (nr, line) in input.0.lines().enumerate() {
            if line.starts_with("*MATERIAL") {
                materials.push(Self::from_definition(input, nr));
            }
        }

        materials
    }
    pub fn from_definition(input: &InpFile, line_nr: usize) -> Self {
        let name = input
            .0
            .lines()
            .nth(line_nr)
            .unwrap()
            .split('=')
            .next_back()
            .unwrap()
            .to_string();

        let mut material = Self {
            name,
            density: None,
            elastic: None,
        };

        let mut lines = input.0.lines();
        while let Some(line) = lines.next() {
            if line.starts_with("*DENSITY") {
                if let Some(def) = lines.next() {
                    material.density = def.parse().ok();
                }
            } else if line.starts_with("*ELASTIC") {
                if let Some(def) = lines.next() {
                    let args: Vec<&str> = def.split(',').map(str::trim).collect();
                    material.elastic = Some((args[0].parse().unwrap(), args[1].parse().unwrap()));
                }
            }
        }

        material
    }
}
