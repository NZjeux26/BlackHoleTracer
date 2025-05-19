#ifndef RAYTRACER_H
#define RAYTRACER_H

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include "3DMathtypes.h"
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL.h>

typedef struct{
    Vec3 position; //current pos
    Vec3 direction; //
    double redshift;
    double intensity;
} RayState;

void raytrace_blackhole(BlackHoleParams params, SDL_Surface* surface);
void apply_strippling_effect(SDL_Surface* surface);
bool trace_rayStep(RayState* ray, BlackHoleParams params);
bool accretion_disk_intersection(Vec3 pos, Vec3 dir, BlackHoleParams params, 
                                        double* intersection_distance, Vec3* intersection_point, Vec3* normal);
Uint32 cal_accretion_disk_colour(Vec3 position, Vec3 view_direction, BlackHoleParams params);
Uint32 trace_black_hole_ray(Vec3 ray_origin, Vec3 ray_dir, BlackHoleParams Params);
BlackHoleParams init_BH_params(double mass, double observer_distance);
double calculate_doppler(Vec3 disk_velocity, Vec3 view_direction);
Vec3 calculate_orbital_velocity(Vec3 position, double schwarzschild_radius);
#endif // RAYTRACER_H
