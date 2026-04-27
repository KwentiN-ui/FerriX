//! Time-dependent amplitude definitions.
//!
//! This module handles "Amplitudes," which are used to scale loads or boundary
//! conditions over time, either linearly (ramp) or via a user-defined time series.

use crate::solver::time::SolverTime;

/// Defines how a value (like a load) scales over simulation time.
#[derive(Debug, Clone)]
pub struct Amplitude {
    /// If true, the time-values are relative to total simulation time;
    /// otherwise, they are relative to the start of the current step.
    pub total_time: bool,
    /// X-axis shift (time offset) for the amplitude data.
    pub shift_x: f64,
    /// Y-axis shift (value offset) for the amplitude data.
    pub shift_y: f64,
    /// Optional piecewise linear time-value pairs.
    pub data: Option<TimeSeries>,
}

impl Default for Amplitude {
    fn default() -> Self {
        Self {
            total_time: false,
            shift_x: 0.,
            shift_y: 0.,
            data: None,
        }
    }
}

impl Amplitude {
    /// Calculates the scaling factor for a given simulation time.
    ///
    /// If no data series is provided, it defaults to a linear ramp from 0 to 1
    /// during the step where it was first defined, and remains 1 thereafter.
    #[must_use]
    pub fn y(&self, time: &SolverTime, origin_step: usize, current_step: usize) -> f64 {
        match &self.data {
            Some(series) => {
                if self.total_time {
                    interpolate(time.global_time(), series, self.shift_x, self.shift_y)
                } else {
                    interpolate(time.local_time(), series, self.shift_x, self.shift_y)
                }
            }
            None => {
                if current_step > origin_step {
                    1.0
                } else {
                    // apply a ramp local to the step
                    time.local_time() / time.local_max_time()
                }
            }
        }
    }
}

/// A collection of time-value pairs used for interpolation.
///
/// The first vector contains the time points, and the second contains the corresponding values.
#[derive(Debug, Clone)]
pub struct TimeSeries(pub Vec<f64>, pub Vec<f64>);

fn interpolate(t: f64, data: &TimeSeries, shift_x: f64, shift_y: f64) -> f64 {
    let times = &data.0;
    let values = &data.1;

    if times.is_empty() {
        return 0.0;
    }

    // shift time-value
    let t_target = t - shift_x;

    // use binary search for interval
    let idx = times.partition_point(|&x| x < t_target);

    // clamp at borders
    if idx == 0 {
        return values[0] + shift_y;
    }
    if idx >= times.len() {
        return values[values.len() - 1] + shift_y;
    }

    let t0 = times[idx - 1];
    let t1 = times[idx];
    let v0 = values[idx - 1];
    let v1 = values[idx];

    // interpolate
    let slope = (v1 - v0) / (t1 - t0);
    let interpolated_v = v0 + slope * (t_target - t0);

    interpolated_v + shift_y
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_data() -> TimeSeries {
        TimeSeries(vec![0.0, 10.0, 20.0], vec![0.0, 100.0, 200.0])
    }

    #[test]
    fn test_exact_match() {
        let data = setup_data();
        let result = interpolate(10.0, &data, 0.0, 0.0);
        assert!((result - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_linear_interpolation() {
        let data = setup_data();
        // Midpoint between 0.0 and 10.0 -> 50.0
        let result = interpolate(5.0, &data, 0.0, 0.0);
        assert!((result - 50.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_shifts() {
        let data = setup_data();
        // t=15 with shift_x=5 corresponds to t_target=10 -> value=100
        // adding shift_y=20 results in 120
        let result = interpolate(15.0, &data, 5.0, 20.0);
        assert!((result - 120.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_clamping() {
        let data = setup_data();
        // Below range
        let result_low = interpolate(-10.0, &data, 0.0, 0.0);
        assert!((result_low - 0.0).abs() < f64::EPSILON);
        // Above range
        let result_high = interpolate(30.0, &data, 0.0, 0.0);
        assert!((result_high - 200.0).abs() < f64::EPSILON);
    }
}
