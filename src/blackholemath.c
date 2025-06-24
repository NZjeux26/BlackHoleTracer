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
    
    // Initialize default disk
    init_default_disk(&params.disk);

    // Calculate Kerr-specific quantities
    calculate_kerr_horizons(&params);
    calculate_critical_orbits(&params);
    calculate_recommended_disk_params(&params);

    // Integration parameters (may need adjustment for Kerr)
    params.eps = 0.01;
    params.dtau = 0.05; 
    params.max_steps = 2000;
    
    return params;
}

void init_default_disk(AccretionDisk* disk) {
     // Set reasonable default values
    disk->opacity = 0.2;
    disk->brightness = 0.01;
    disk->turbulence_strength = 0.9;
    disk->spiral_arms = 4;
    disk->spiral_tightness = 0.5;
    disk->rotation_speed = 4.0;
    disk->turbulence_speed = 0.4;
   // disk->volume_samples = 64; // Good for offline rendering
    disk->self_shadowing = 0.5;
    disk->magnetic_field_strength = 2.0;
    disk->viscosity_alpha = 0.8; 
}

void set_disk_turbulence_params(AccretionDisk* disk, double strength, double speed) {
    disk->turbulence_strength = fmax(0.0, fmin(1.0, strength));
    disk->turbulence_speed = fmax(0.0, speed);
}

void set_disk_spiral_params(AccretionDisk* disk, int arms, double tightness) {
    disk->spiral_arms = fmax(1, fmin(8, arms));
    disk->spiral_tightness = fmax(0.1, fmin(2.0, tightness));
}

void set_disk_rendering_params(AccretionDisk* disk, int volume_samples, double brightness) {
    disk->volume_samples = fmax(16, fmin(256, volume_samples));
    disk->brightness = fmax(0.1, brightness);
}

void calculate_recommended_disk_params(BlackHoleParams* params) {
    double safety_margin = 1.15;
    double spin_correction = 1.0 + 0.2 * fabs(params->spin);
    params->disk.inner_radius = safety_margin * spin_correction * params->isco_radius_prograde;
   
    // Set outer radius based on observer distance
    params->disk.outer_radius = 25;
    
    // Scale thickness with inner radius
    params->disk.thickness = 0.1 + 0.4 * pow(params->isco_radius_prograde / 6.0, 0.5); // Thicker for slower-spinning BHs
    
    // Adjust temperature
    params->disk.temperature_factor = 8000;
    
    // More turbulence for higher spin
    params->disk.turbulence_strength = 0.2 + 0.3 * fabs(params->spin);
    
    // Spiral tightness increases with spin
    params->disk.spiral_tightness = 0.3 + 0.4 * fabs(params->spin);
    
}
// Calculate Kerr horizons and ergosphere
void calculate_kerr_horizons(BlackHoleParams* params) {
    double M = params->mass;
    double a = params->spin_dimensional;
    double a_squared = a * a;
    double M_squared = M * M;
    
    // Event horizons
    double delta_discriminant = M_squared - a_squared;
    if (delta_discriminant >= 0) {
        params->kerr_radius_outer = M + sqrt(delta_discriminant);
        params->kerr_radius_inner = M - sqrt(delta_discriminant);
    } else {
        // Naked singularity case (a > M) - shouldn't happen in practice
        params->kerr_radius_outer = M;
        params->kerr_radius_inner = 0.0;
        fprintf(stderr, "Warning: Naked singularity detected (a > M)\n");
    }
    
    // Ergosphere boundaries (static limit surface)
    // At equator (θ = π/2): r = M + sqrt(M² - a²cos²θ) = M + sqrt(M² - 0) = M + M = 2M
    // At poles (θ = 0, π): r = M + sqrt(M² - a²) = r+
    params->ergosphere_radius_eq = 2.0 * M; // = 2M at equator //GPT said it would be better as params->ergosphere_radius_eq = M + sqrt(M*M - a*a);
    params->ergosphere_radius_pole = params->kerr_radius_outer;
    params->static_limit_eq = params->ergosphere_radius_eq;
}

// Calculate critical orbital radii for Kerr geometry
void calculate_critical_orbits(BlackHoleParams* params) {
    double M = params->mass;
    double a_norm = params->spin;  // Use dimensionless spin consistently
    
    // ISCO calculation with correct dimensionless spin
    double Z1 = 1.0 + pow(1.0 - a_norm*a_norm, 1.0/3.0) * (pow(1.0 + a_norm, 1.0/3.0) + pow(1.0 - a_norm, 1.0/3.0));
    double Z2 = sqrt(3.0*a_norm*a_norm + Z1*Z1);

    // Prograde ISCO (co-rotating)
    params->isco_radius_prograde = M * (3.0 + Z2 - sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2)));
    
    // Retrograde ISCO (counter-rotating)  
    params->isco_radius_retrograde = M * (3.0 + Z2 + sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2)));
    
    // Photon sphere radii (unstable circular photon orbits)
    // Simplified approximations - exact calculation requires solving quartic
    double rho_plus = 2.0*M*(1.0 + cos(2.0/3.0 * acos(-a_norm)));
    double rho_minus = 2.0*M*(1.0 + cos(2.0/3.0 * acos(a_norm)));
    
    params->photon_sphere_prograde = rho_plus;
    params->photon_sphere_retrograde = rho_minus;
    
    // Ensure physical consistency
    if (params->photon_sphere_prograde < params->kerr_radius_outer) {
        params->photon_sphere_prograde = 1.5 * params->kerr_radius_outer;
    }
    if (params->photon_sphere_retrograde < params->kerr_radius_outer) {
        params->photon_sphere_retrograde = 1.5 * params->kerr_radius_outer;
    }
}

