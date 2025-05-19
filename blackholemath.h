#ifndef BLACKHOLEMATH_H
#define BLACKHOLEMATH_H

#include "3DMathtypes.h"
#include <math.h>

typedef struct {
    double inner_radius; // Inner stable orbit radius
    double outer_radius; // Outer extent of disk
    double thickness; // Vertical thickness of disk
    double opacity; // Transparency factor
    double temperature_factor; // Controls disk temperature/brightness
} AccretionDisk;

// Black hole and simulation parameters
typedef struct {
    double mass;             // Black hole mass
    double schwarzschild_radius; // 2*G*M/c^2
    AccretionDisk disk; // Accretion disk properties
    double c_squard; // Speed of light squared (for calculations)
    double dt;  // Integration time step
    int max_steps;  // Maximum integration steps
    double observer_distance; // Distance of camera from black hole
} BlackHoleParams;
#endif // BLACKHOLEMATH_H