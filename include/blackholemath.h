#ifndef BLACKHOLEMATH_H
#define BLACKHOLEMATH_H

#include "3DMathtypes.h"
#include <sys/types.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct {
    double inner_radius; // Inner stable orbit radius
    double outer_radius; // Outer extent of disk
    double thickness; // Vertical thickness of disk
    double opacity; // Transparency factor
    double temperature_factor; // Controls disk temperature/brightness
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
    double isco_radius_prograde;    // Innermost stable circular orbit (prograde)
    double isco_radius_retrograde;  // Innermost stable circular orbit (retrograde)
    double photon_sphere_prograde;  // Photon sphere radius (prograde)
    double photon_sphere_retrograde; // Photon sphere radius (retrograde)
    
    // Observer parameters
    double observer_distance;       // Distance from black hole center
    double observer_inclination;    // Viewing angle (0 = pole-on, π/2 = edge-on)
    
    // Physical constants (geometric units)
    double c_squared;               // Speed of light squared

    // Accretion disk properties
    AccretionDisk disk;
    
    // Integration parameters
    double dtau;                      // Integration step size
    double eps;
    int max_steps;                  // Maximum ray tracing steps
} BlackHoleParams;

BlackHoleParams init_BH_params(double mass, double spin, double observer_distance);
void calculate_kerr_horizons(BlackHoleParams* params);
void calculate_critical_orbits(BlackHoleParams* params);
#endif // BLACKHOLEMATH_H