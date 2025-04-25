#ifndef BLACKHOLEMATH_H
#define BLACKHOLEMATH_H

// Include necessary standard or project-specific headers here
#include <cmath>
#include <iostream>

// Black hole and simulation parameters
typedef struct {
    double mass;             // Black hole mass
    double schwarzschild_radius; // 2*G*M/c^2
    double accretion_disk_inner_radius;
    double accretion_disk_outer_radius;
    double observer_distance;
} BlackHoleParams;

#endif // BLACKHOLEMATH_H