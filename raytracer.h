#ifndef RAYTRACER_H
#define RAYTRACER_H

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include "types.h"
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL.h>

void raytrace_blackhole(BlackHoleParams params, SDL_Surface* surface);
void apply_strippling_effect(SDL_Surface* surface);
bool trace_rayStep(Vec3* pos, Vec3* dir, BlackHoleParams params, double step_size);
bool check_accretion_disk_intersection(Vec3 pos, Vec3 dir, BlackHoleParams params, 
                                        double* intersection_distance, Vec3* intersection_point);
Uint32 cal_accretion_disk_colour(Vec3 intersection_point, BlackHoleParams params);
Uint32 trace_black_hole_ray(Vec3 ray_origin, Vec3 ray_dir, BlackHoleParams Params);

#endif // RAYTRACER_H
