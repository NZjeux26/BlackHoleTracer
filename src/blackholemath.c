#include "blackholemath.h"

BlackHoleParams init_BH_params(double mass, double spin, double observer_distance){
    BlackHoleParams params = {0}; // Initialize all fields to zero
    
    // Core properties
    params.mass = mass;
    params.spin = fmax(-0.97, fmin(0.97, spin)); // Clamp to physical range I clamp to 0.97 since higher gives artifacting
    params.spin_dimensional = params.spin * mass;   // a = (a/M) * M
    
    // Basic Schwarzschild quantities
    params.schwarzschild_radius = 2.0 * mass;
    params.c_squared = 1.0; // Geometric units where c = 1
    
    // Observer parameters
    params.observer_distance = observer_distance;
    
    // Accretion disk parameters (adjusted for Kerr geometry)(Needs fixed up for Hamiltonian)
    params.disk.inner_radius = params.isco_radius_prograde; // Start at ISCO
    params.disk.outer_radius = 15.0 * params.schwarzschild_radius;
    params.disk.thickness = 0.2 * params.schwarzschild_radius;
    params.disk.opacity = 0.2;
    params.disk.temperature_factor = 1.0;
    
    // Integration parameters (may need adjustment for Kerr)
    params.eps = 0.01;
    params.dtau = 0.05; 
    params.max_steps = 2500;
    
    return params;
}


