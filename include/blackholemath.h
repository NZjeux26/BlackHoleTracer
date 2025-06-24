#ifndef BLACKHOLEMATH_H
#define BLACKHOLEMATH_H

#include "3DMathtypes.h"
#include <sys/types.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct {
    //Basic Params
    double inner_radius;       // Inner stable orbit radius
    double outer_radius;       // Outer extent of disk
    double thickness;          // Vertical thickness of disk
    double opacity;            // Transparency factor
    double temperature_factor; // Controls disk temperature/brightness

    //Volumetric params
    double brightness;          // Overall disk brightness multiplier
    double turbulence_strength; // How much turbulence affects density (0-1)
    int spiral_arms;            // Number of spiral arms (typically 2-4)
    double spiral_tightness;    // How tightly wound the spirals are

    // Animation parameters
    double rotation_speed;      // Disk rotation speed (for animation)
    double turbulence_speed;    // Speed of turbulence evolution

    // Rendering parameters
    int volume_samples;         // Number of samples for volumetric rendering
    double self_shadowing;      // Strength of self-shadowing effects

    // Physical model parameters
    double magnetic_field_strength; // For future magnetic field visualization
    double viscosity_alpha;         // Alpha parameter for viscous heating
} AccretionDisk;

// Black hole and simulation parameters
typedef struct {
    // Core black hole properties
    double mass;                    // Black hole mass (M)
    double spin;                    // Dimensionless spin parameter (a/M, where -1 ≤ a ≤ 1)
    double spin_dimensional;        // Dimensional spin parameter (a = J/Mc)
    
    // Derived Kerr quantities
    double schwarzschild_radius;    // rs = 2M
    double kerr_radius_outer;       // r+ = M + sqrt(M² - a²) (outer event horizon)
    double kerr_radius_inner;       // r- = M - sqrt(M² - a²) (inner event horizon)
    double ergosphere_radius_eq;    // Ergosphere radius at equator
    double ergosphere_radius_pole;  // Ergosphere radius at poles
    double static_limit_eq;         // Static limit surface at equator
    
    // Critical orbits (depend on spin)
    double isco_radius_prograde;     // Innermost stable circular orbit (prograde)
    double isco_radius_retrograde;   // Innermost stable circular orbit (retrograde)
    double photon_sphere_prograde;   // Photon sphere radius (prograde)
    double photon_sphere_retrograde; // Photon sphere radius (retrograde)
    
    // Observer parameters
    double observer_distance;       // Distance from black hole center
    
    // Physical constants (geometric units)
    double c_squared;               // Speed of light squared

    // Accretion disk properties
    AccretionDisk disk;
    
    // Integration parameters
    double dtau;                    // Integration step size
    double eps;                     // Numerical precision parameter
    int max_steps;                  // Maximum ray tracing steps
} BlackHoleParams;

BlackHoleParams init_BH_params(double mass, double spin, double observer_distance);
void calculate_kerr_horizons(BlackHoleParams* params);
void calculate_critical_orbits(BlackHoleParams* params);

void init_default_disk(AccretionDisk* disk);
void set_disk_turbulence_params(AccretionDisk* disk, double strength, double speed);
void set_disk_spiral_params(AccretionDisk* disk, int arms, double tightness);
void set_disk_rendering_params(AccretionDisk* disk, int volume_samples, double brightness);

// Helper function to calculate recommended parameters
void calculate_recommended_disk_params(BlackHoleParams* params);

#endif // BLACKHOLEMATH_H