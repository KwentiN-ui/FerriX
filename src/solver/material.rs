//! This module contains data structures for material definitions. A Material is defined using the `*MATERIAL` card
//! globally. To use a material, it has to be referenced in a section definition.
#![allow(unused)]

use ndarray::Array2;

use crate::solver::inp::InpFile;

#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub density: Option<f64>,
    /// Elastic Modulus, Poisson
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

    /// constructs linear elastic materialmatrix D (6x6 for 3D)
    /// Voigt-Notation: xx, yy, zz, xy, yz, zx
    pub fn build_elastic_d_matrix(&self) -> Array2<f64> {
        let (e, nu) = self.elastic.expect("No *ELASTIC Card definition found!");

        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c1 = 1.0 - nu;
        let c2 = nu;
        let c3 = (1.0 - 2.0 * nu) / 2.0;

        let data = vec![
            c1, c2, c2, 0., 0., 0., c2, c1, c2, 0., 0., 0., c2, c2, c1, 0., 0., 0., 0., 0., 0., c3,
            0., 0., 0., 0., 0., 0., c3, 0., 0., 0., 0., 0., 0., c3,
        ];

        Array2::from_shape_vec((6, 6), data).expect("Matrix shape error") * factor
    }
}
