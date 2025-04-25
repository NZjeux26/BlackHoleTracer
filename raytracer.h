#ifndef RAYTRACER_H
#define RAYTRACER_H

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <blackholemath.h>
#include <SDL2/SDL_surface.h>

typedef struct {
    double x, y, z;
} Vec3;

Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_scale(Vec3 v, double s);
double vec3_length(Vec3 v);
Vec3 vec3_normalise(Vec3 v);    

void raytrace_blackhole(BlackHoleParams params, SDL_Surface* surface);
#endif // RAYTRACER_H
