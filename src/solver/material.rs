//! Material property definitions.
//!
//! This module contains structures for defining physical properties like density
//! and elasticity, and provides methods for generating material law matrices.

use crate::solver::error::{FerrixError, Result};
use nalgebra::DMatrix;

/// A Temperature Dependent Look-Up Table (LUT) for scalar values.
///
/// Stores data as a list of (temperature, value) pairs, sorted by temperature.
/// Provides linear interpolation between data points.
#[derive(Debug, Clone, Default)]
pub struct TemperatureDependentLUT {
    pub data: Vec<(f64, f64)>,
}

impl TemperatureDependentLUT {
    /// Creates a new LUT from a list of (temperature, value) pairs.
    ///
    /// # Panics
    /// Panics if any temperature is NaN.
    #[must_use]
    pub fn new(mut data: Vec<(f64, f64)>) -> Self {
        data.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("NaN temperature in LUT"));
        Self { data }
    }

    /// Interpolates the value at a given temperature.
    ///
    /// Clamps to the first or last value if the temperature is out of range.
    ///
    /// # Panics
    /// Panics if the LUT is empty.
    #[must_use]
    pub fn interpolate(&self, temp: f64) -> f64 {
        assert!(
            !self.data.is_empty(),
            "Attempted to interpolate an empty LUT"
        );

        if self.data.len() == 1 || temp <= self.data[0].0 {
            return self.data[0].1;
        }

        if temp >= self.data.last().unwrap().0 {
            return self.data.last().unwrap().1;
        }

        // Binary search to find the interval [low, high]
        let mut low = 0;
        let mut high = self.data.len() - 1;

        while high - low > 1 {
            let mid = low + (high - low) / 2;
            if self.data[mid].0 <= temp {
                low = mid;
            } else {
                high = mid;
            }
        }

        let (t0, v0) = self.data[low];
        let (t1, v1) = self.data[high];

        let factor = (temp - t0) / (t1 - t0);
        v0 + (v1 - v0) * factor
    }
}

/// Represents the solution-dependent state variables at a single material point.
#[derive(Debug, Clone, Default)]
pub struct MaterialPointState {
    pub variables: Vec<f64>,
}

impl std::ops::Index<usize> for MaterialPointState {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.variables[index]
    }
}

impl std::ops::IndexMut<usize> for MaterialPointState {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.variables[index]
    }
}

/// Represents the solution-dependent state variables for all integration points of an element.
#[derive(Debug, Clone, Default)]
pub struct ElementMaterialState {
    pub ip_states: Vec<MaterialPointState>,
}

impl ElementMaterialState {
    /// Creates a new element state with specified number of integration points and variables.
    #[must_use]
    pub fn new(num_ips: usize, num_vars: usize) -> Self {
        Self {
            ip_states: vec![
                MaterialPointState {
                    variables: vec![0.0; num_vars]
                };
                num_ips
            ],
        }
    }
}

impl std::ops::Index<usize> for ElementMaterialState {
    type Output = MaterialPointState;
    fn index(&self, index: usize) -> &Self::Output {
        &self.ip_states[index]
    }
}

impl std::ops::IndexMut<usize> for ElementMaterialState {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.ip_states[index]
    }
}

/// Defines the physical and mechanical properties of a material.
///
/// Properties return `Option<f64>` to allow for sparse material definitions.
/// New properties can be added with a default `None` return to maintain backward compatibility.
pub trait Material: std::fmt::Debug + Send + Sync {
    /// Returns the unique name of the material.
    fn name(&self) -> &str;

    /// Returns the mass density of the material at a given temperature.
    fn density(&self, temp: f64) -> Option<f64> {
        let _ = temp;
        None
    }

    /// Returns Young's Modulus (E) at a given temperature.
    fn youngs_modulus(&self, temp: f64) -> Option<f64> {
        let _ = temp;
        None
    }

    /// Returns Poisson's ratio (nu) at a given temperature.
    fn poisson_ratio(&self, temp: f64) -> Option<f64> {
        let _ = temp;
        None
    }

    /// Returns the thermal expansion coefficient (alpha) at a given temperature.
    fn thermal_expansion(&self, temp: f64) -> Option<f64> {
        let _ = temp;
        None
    }

    /// Returns the reference temperature for thermal expansion.
    fn reference_temperature(&self) -> f64 {
        0.0
    }

    /// Returns the number of solution-dependent state variables (SDVs) for this material.
    fn num_state_variables(&self) -> usize {
        0
    }

    /// Updates the state-dependent variables and returns the tangent stiffness and stress.
    ///
    /// This is a simplified version of UMAT.
    /// # Arguments
    /// * `temp` - Temperature at the integration point.
    /// * `strain` - Total strain at the end of the increment.
    /// * `state_old` - State variables at the start of the increment.
    /// * `dtime` - Time increment.
    ///
    /// # Returns
    /// A tuple of (Tangent Matrix, Stress Vector, Updated State Variables).
    ///
    /// # Errors
    /// Returns `FerrixError` if the update fails.
    fn update_state(
        &self,
        temp: f64,
        strain: &nalgebra::DVector<f64>,
        state_old: &MaterialPointState,
        _dtime: f64,
    ) -> Result<(
        nalgebra::DMatrix<f64>,
        nalgebra::DVector<f64>,
        MaterialPointState,
    )> {
        // Default implementation for linear material
        let d_mat = self.build_elastic_d_matrix(temp)?;
        let stress = &d_mat * strain;
        Ok((d_mat, stress, state_old.clone()))
    }

    /// Builds the elastic constitutive matrix (D-matrix) for the material (6x6 for 3D).
    /// Uses Voigt notation: [xx, yy, zz, xy, yz, zx].
    ///
    /// # Errors
    /// Returns `FerrixError::InvalidModelState` if Young's Modulus or Poisson's Ratio is missing.
    fn build_elastic_d_matrix(&self, temp: f64) -> Result<DMatrix<f64>> {
        let e = self.youngs_modulus(temp).ok_or_else(|| {
            FerrixError::InvalidModelState(format!(
                "No Young's Modulus definition found for material '{}'!",
                self.name()
            ))
        })?;
        let nu = self.poisson_ratio(temp).ok_or_else(|| {
            FerrixError::InvalidModelState(format!(
                "No Poisson's Ratio definition found for material '{}'!",
                self.name()
            ))
        })?;

        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c1 = 1.0 - nu;
        let c2 = nu;
        let c3 = (1.0 - 2.0 * nu) / 2.0;

        Ok(DMatrix::from_row_slice(
            6,
            6,
            &[
                c1, c2, c2, 0.0, 0.0, 0.0, c2, c1, c2, 0.0, 0.0, 0.0, c2, c2, c1, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, c3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, c3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                c3,
            ],
        ) * factor)
    }
}

/// Base implementation of the `Material` trait using LUTs for property storage.
#[derive(Debug, Clone)]
pub struct BaseMaterial {
    pub name: String,
    pub density: Option<TemperatureDependentLUT>,
    pub youngs_modulus: Option<TemperatureDependentLUT>,
    pub poisson_ratio: Option<TemperatureDependentLUT>,
    pub thermal_expansion: Option<TemperatureDependentLUT>,
    pub reference_temperature: f64,
    pub num_depvars: usize,
}

impl Material for BaseMaterial {
    fn name(&self) -> &str {
        &self.name
    }

    fn num_state_variables(&self) -> usize {
        self.num_depvars
    }

    fn density(&self, temp: f64) -> Option<f64> {
        self.density.as_ref().map(|lut| lut.interpolate(temp))
    }

    fn youngs_modulus(&self, temp: f64) -> Option<f64> {
        self.youngs_modulus
            .as_ref()
            .map(|lut| lut.interpolate(temp))
    }

    fn poisson_ratio(&self, temp: f64) -> Option<f64> {
        self.poisson_ratio.as_ref().map(|lut| lut.interpolate(temp))
    }

    fn thermal_expansion(&self, temp: f64) -> Option<f64> {
        self.thermal_expansion
            .as_ref()
            .map(|lut| lut.interpolate(temp))
    }

    fn reference_temperature(&self) -> f64 {
        self.reference_temperature
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_lut_interpolation() {
        let lut = TemperatureDependentLUT::new(vec![(0.0, 100.0), (100.0, 200.0), (200.0, 400.0)]);

        assert_eq!(lut.interpolate(-10.0), 100.0);
        assert_eq!(lut.interpolate(0.0), 100.0);
        assert_eq!(lut.interpolate(50.0), 150.0);
        assert_eq!(lut.interpolate(100.0), 200.0);
        assert_eq!(lut.interpolate(150.0), 300.0);
        assert_eq!(lut.interpolate(200.0), 400.0);
        assert_eq!(lut.interpolate(250.0), 400.0);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_base_material_elastic() {
        let e_lut = TemperatureDependentLUT::new(vec![(0.0, 210_000.0), (500.0, 150_000.0)]);
        let nu_lut = TemperatureDependentLUT::new(vec![(0.0, 0.3), (500.0, 0.35)]);

        let material = BaseMaterial {
            name: "Steel".to_string(),
            density: None,
            youngs_modulus: Some(e_lut),
            poisson_ratio: Some(nu_lut),
            thermal_expansion: None,
            reference_temperature: 0.0,
            num_depvars: 0,
        };

        assert_eq!(material.youngs_modulus(0.0), Some(210_000.0));
        assert_eq!(material.poisson_ratio(0.0), Some(0.3));

        assert_eq!(material.youngs_modulus(250.0), Some(180_000.0));
        let nu_interp = material.poisson_ratio(250.0).unwrap();
        assert!((nu_interp - 0.325).abs() < 1e-12);
    }
}
