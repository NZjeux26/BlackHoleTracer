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

// Boyer-Lindquist coordinate structure
typedef struct {
    double r;      // radial coordinate
    double theta;  // polar angle
    double phi;    // azimuthal angle
} BoyerLindquistCoords;

// Calculate metric components in Boyer-Lindquist coordinates
typedef struct {
    double g_tt, g_tr, g_rr, g_theta_theta, g_phi_phi, g_t_phi;
    double rho2, delta;
} KerrMetric;

typedef struct {
    double x;      // Cartesian x-coordinate
    double y;      // Cartesian y-coordinate
    double z;      // Cartesian z-coordinate
} CartesianCoords;

// 4-velocity structures
typedef struct {
    double dt;     // Time component
    double dr;     // Radial component
    double dtheta; // Polar component
    double dphi;   // Azimuthal component
} FourVelocityBL;

typedef struct {
    double dt;     // Time component
    double dx;     // x-component
    double dy;     // y-component
    double dz;     // z-component
} FourVelocityCartesian;

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
    double min_dt;                      // Integration step size
    double max_dt;
    int max_steps;                  // Maximum ray tracing steps
} BlackHoleParams;

BlackHoleParams init_BH_params(double mass, double spin, double observer_distance);
void calculate_kerr_horizons(BlackHoleParams* params);
void calculate_critical_orbits(BlackHoleParams* params);
void print_kerr_properties(const BlackHoleParams* params);
// Coordinate conversion functions
BoyerLindquistCoords cartesian_to_boyer_lindquist(double x, double y, double z, double a);
CartesianCoords boyer_lindquist_to_cartesian(double r, double theta, double phi, double a);
// Kerr metric helper functions
double calculate_rho_squared(double r, double theta, double a);
double calculate_delta(double r, double M, double a);
double calculate_sigma(double r, double theta, double M, double a);
// Spacetime region checking functions
bool is_inside_ergosphere(double r, double theta, double M, double a);
bool is_inside_event_horizon(double r, double M, double a);
// 4-velocity conversion functions
FourVelocityBL cartesian_4vel_to_boyer_lindquist(
    double x, double y, double z,
    double vx, double vy, double vz,
    double a);

FourVelocityCartesian boyer_lindquist_4vel_to_cartesian(
    double r, double theta, double phi,
    double vt, double vr, double vtheta, double vphi,
    double a);
#endif // BLACKHOLEMATH_H