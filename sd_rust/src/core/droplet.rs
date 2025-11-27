//! Superdroplet data structure and methods.
//!
//! A superdroplet represents a collection of identical cloud droplets,
//! characterized by their multiplicity and physical properties.

use crate::core::constants::{FOUR_THIRD, PI, RHO_WATER};

/// Represents a superdroplet in the simulation.
///
/// A superdroplet is a computational particle that represents many
/// identical physical droplets. The `multi` field indicates how many
/// real droplets this superdroplet represents.
#[derive(Clone)]
pub struct Droplet {
    /// Multiplicity - number of real droplets represented
    pub multi: usize,

    /// Radius cubed (m³) - stored for efficiency
    pub rcubed: f64,

    /// Solute mass (kg)
    pub solute: f64,

    /// Density (kg/m³) - private, typically constant
    _density: f64,

    /// Unique identifier - private, for debugging
    _id: usize,

    /// Cached radius (m)
    radius: f64,

    /// Cached volume (m³)
    volume: f64,

    /// Cached mass (kg)
    mass: f64,

    /// Cached terminal velocity (m/s)
    terminal_velocity: f64,
}

impl Droplet {
    /// Creates a new droplet with the given multiplicity and radius.
    ///
    /// # Arguments
    /// * `multi` - Number of real droplets this superdroplet represents
    /// * `radius` - Radius of the droplet in meters
    ///
    /// # Example
    /// ```
    /// use sd_rust::core::Droplet;
    /// let droplet = Droplet::new(1000, 1e-6);  // 1000 droplets of 1 micron radius
    /// ```
    pub fn new(multi: usize, radius: f64) -> Droplet {
        let rcubed = radius.powi(3);
        let volume = rcubed * PI * FOUR_THIRD;
        let mass = volume * RHO_WATER;
        let terminal_velocity = Self::compute_terminal_velocity(radius, mass);

        Droplet {
            multi,
            rcubed,
            solute: 0.,
            _density: RHO_WATER,
            _id: 0,
            radius,
            volume,
            mass,
            terminal_velocity,
        }
    }

    /// Updates the radius cubed and recalculates cached properties.
    ///
    /// This method should be called whenever the droplet size changes
    /// (e.g., after coalescence).
    ///
    /// # Arguments
    /// * `new_rcubed` - New value for radius cubed (m³)
    pub fn update_rcubed(&mut self, new_rcubed: f64) {
        self.rcubed = new_rcubed;
        self.radius = new_rcubed.cbrt();
        self.volume = new_rcubed * PI * FOUR_THIRD;
        self.mass = self.volume * RHO_WATER;
        self.terminal_velocity = Self::compute_terminal_velocity(self.radius, self.mass);
    }

    /// Returns the volume of the droplet (m³).
    #[inline(always)]
    pub fn volume(&self) -> f64 {
        self.volume
    }

    /// Returns the radius of the droplet (m).
    #[inline(always)]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the mass of the droplet (kg).
    #[inline(always)]
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Returns the terminal velocity of the droplet (m/s).
    #[inline(always)]
    pub fn terminal_velocity(&self) -> f64 {
        self.terminal_velocity
    }

    /// Computes terminal velocity using Beard (1976) parameterization.
    ///
    /// # Arguments
    /// * `radius` - Droplet radius in meters
    /// * `mass` - Droplet mass in kilograms
    ///
    /// # Returns
    /// Terminal velocity in m/s
    fn compute_terminal_velocity(radius: f64, mass: f64) -> f64 {
        let d = 2.0 * radius * 1e6; // diameter, m -> μm
        let x = mass * 1e3; // mass, kg -> g

        let (alpha, x_to_beta) = match d {
            d if d <= 134.43 => {
                // x^(2/3) = cbrt(x)^2
                let cbrt_x = x.cbrt();
                (4.5795e5, cbrt_x * cbrt_x)
            }
            d if d <= 1511.64 => {
                // x^(1/3) = cbrt(x)
                (4962.0, x.cbrt())
            }
            d if d <= 3477.84 => {
                // x^(1/6) = sqrt(cbrt(x))
                (1732.0, x.cbrt().sqrt())
            }
            _ => (917.0, 1.0),
        };

        1e-2 * alpha * x_to_beta // cm/s -> m/s
    }
}

/// Calculates the total water mass in a collection of droplets.
///
/// # Arguments
/// * `droplets` - Slice of droplets to sum
///
/// # Returns
/// Total water mass in kilograms
#[inline(always)]
pub fn total_water(droplets: &[Droplet]) -> f64 {
    droplets
        .iter()
        .map(|droplet| droplet.mass() * (droplet.multi as f64))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_droplet_creation() {
        let d = Droplet::new(1000, 1e-6);
        assert_eq!(d.multi, 1000);
        assert!(d.radius() > 0.0);
        assert!(d.volume() > 0.0);
        assert!(d.mass() > 0.0);
    }

    #[test]
    fn test_droplet_volume() {
        let d = Droplet::new(1, 1e-6);
        let expected = 4.0 / 3.0 * PI * 1e-18;
        assert!((d.volume() - expected).abs() < 1e-25);
    }

    #[test]
    fn test_update_rcubed() {
        let mut d = Droplet::new(1, 1e-6);
        let old_mass = d.mass();

        d.update_rcubed(8e-18); // Double the radius
        assert!(d.mass() > old_mass);
    }

    #[test]
    fn test_total_water() {
        let droplets = vec![Droplet::new(1000, 1e-6), Droplet::new(2000, 2e-6)];
        let total = total_water(&droplets);
        assert!(total > 0.0);
    }
}
