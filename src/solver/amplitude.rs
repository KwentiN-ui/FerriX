use crate::solver::time::SolverTime;

#[derive(Debug, Clone)]
pub struct Amplitude {
    pub name: String,
    /// Defines, if the time-values are defined in reference to the total simulation time (over multiple steps)
    pub total_time: bool,
    pub shift_x: f64,
    pub shift_y: f64,
    pub data: Option<TimeSeries>,
}

impl Default for Amplitude {
    fn default() -> Self {
        Self {
            name: "default".into(),
            total_time: false,
            shift_x: 0.,
            shift_y: 0.,
            data: None,
        }
    }
}

impl Amplitude {
    pub fn y(&self, time: &SolverTime) -> f64 {
        match &self.data {
            Some(series) => {
                if self.total_time {
                    interpolate(time.global_time(), series, self.shift_x, self.shift_y)
                } else {
                    interpolate(time.local_time(), series, self.shift_x, self.shift_y)
                }
            }
            None => {
                // apply a ramp local to the step
                time.local_time() / time.local_max_time()
            }
        }
    }
}

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
        assert_eq!(interpolate(10.0, &data, 0.0, 0.0), 100.0);
    }

    #[test]
    fn test_linear_interpolation() {
        let data = setup_data();
        // Midpoint between 0.0 and 10.0 -> 50.0
        assert_eq!(interpolate(5.0, &data, 0.0, 0.0), 50.0);
    }

    #[test]
    fn test_shifts() {
        let data = setup_data();
        // t=15 with shift_x=5 corresponds to t_target=10 -> value=100
        // adding shift_y=20 results in 120
        assert_eq!(interpolate(15.0, &data, 5.0, 20.0), 120.0);
    }

    #[test]
    fn test_clamping() {
        let data = setup_data();
        // Below range
        assert_eq!(interpolate(-10.0, &data, 0.0, 0.0), 0.0);
        // Above range
        assert_eq!(interpolate(30.0, &data, 0.0, 0.0), 200.0);
    }
}
