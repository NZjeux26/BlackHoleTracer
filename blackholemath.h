#ifndef BLACKHOLEMATH_H
#define BLACKHOLEMATH_H

#include "types.h"
#include <math.h>

// Black hole and simulation parameters
typedef struct {
    double mass;             // Black hole mass
    double schwarzschild_radius; // 2*G*M/c^2
    double accretion_disk_inner_radius;
    double accretion_disk_outer_radius;
    double accretion_disk_thickness;
    double observer_distance;
} BlackHoleParams;

// Structure to represent the state of a photon in Schwarzschild spacetime
// Using spherical coordinates (t, r, θ, φ) and their derivatives with respect to affine parameter λ
typedef struct {
    double r;       // radial coordinate
    double theta;   // polar angle
    double phi;     // azimuthal angle
    double pr;      // radial momentum component
    double ptheta;  // θ momentum component
    double pphi;    // φ momentum component (angular momentum, constant of motion)
    double E;       // Energy (constant of motion)
} PhotonState;

void cartesian_to_schwarzschild(Vec3 pos, double* r, double* theta, double* phi);
Vec3 schwarzschild_to_cartesian(double r, double theta, double phi);
PhotonState init_photon_state(Vec3 pos, Vec3 dir, double schwarzschild_radius);
void geodesic_derivatives(PhotonState state, double rs, PhotonState* derivatives);
void rk4_step(PhotonState* state, double rs, double stepsize);
double calculate_adaptive_step(PhotonState state, double rs);
bool trace_schwarzschild_geodesic(Vec3* pos, Vec3* dir, BlackHoleParams params, double stepsize);
#endif // BLACKHOLEMATH_H