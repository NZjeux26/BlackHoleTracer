#include "blackholemath.h"

BlackHoleParams init_BH_params(double mass, double observer_distance){
    BlackHoleParams params;

    //core properties
    params.mass = mass;
    params.schwarzschild_radius = 2.0 * mass; //In geometric units where G=c=1
    params.c_squard = 1.0; // Speed of light (in geometric units)
    params.observer_distance = observer_distance;

    //Accretion Disk
    params.disk.inner_radius = 3.0 * params.schwarzschild_radius;
    params.disk.outer_radius = 15 * params.schwarzschild_radius;
    params.disk.thickness = 0.2 * params.schwarzschild_radius;
    params.disk.opacity = 0.2;
    params.disk.temperature_factor = 1.0;

    //intergration params
    params.dt = 0.05 * params.schwarzschild_radius;
    params.max_steps = 2000;

    return params;
}