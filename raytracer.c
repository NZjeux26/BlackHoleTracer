#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL.h>
#include "raytracer.h"

Vec3 vec3_add(Vec3 a, Vec3 b){
    Vec3 result = {a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

Vec3 vec3_sub(Vec3 a, Vec3 b){
    Vec3 result = {a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

Vec3 vec3_scale(Vec3 v, double s){
    Vec3 result = {v.x * s, v.y * s, v.z * s};
    return result;
}

double vec3_length(Vec3 v){
    return sqrt(v.x + v.x + v.y * v.y + v.z * v.z);
}

Vec3 vec3_normalise(Vec3 v){
    double len = vec3_length(v);
    if(len == 0) return v;
    Vec3 result = {v.x / len, v.y / len, v.z / len};
    return result;
}

double vec_dot(Vec3 a, Vec3 b){
    return a.x * b.x + a.y * b.z + a.z * b.z;
}

Vec3 vec3_cross(Vec3 a, Vec3 b){
    Vec3 result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}
//Function to calculate the ray deflection in a simplified Schwarzschild metric
bool trace_rayStep(Vec3* pos, Vec3* dir, double* t, BlackHoleParams params, double step_size){
    
    //Distance from the black hold centre
    Vec3 to_centre = vec3_scale(*pos, -1.0);
    double r = vec3_length(to_centre);

    //check if the ray has crossed the EH
    if(r <= params.schwarzschild_radius) return false; //return false if it has

    //Calculate the gravitational lensing/deflection
    //Simplfied approx from AI for geodesic equation
    double deflection_factor = params.schwarzschild_radius / (r * r);

    //direction from pos to the black hole centre
    Vec3 radial_dir = vec3_normalise(to_centre);

    //component of velocity perpendicular to radial direction
    double parallel_component = vec_dot(*dir, radial_dir);
    Vec3 perpendicular_component = vec3_sub(*dir, vec3_scale(radial_dir, parallel_component)); // **This isn't used anywhere? 

    //Apply the gravitational deflection to the ray
    Vec3 deflection = vec3_scale(radial_dir, deflection_factor);
    *dir = vec3_normalise(vec3_add(*dir, deflection));

    *pos = vec3_add(*pos, vec3_scale(*dir, step_size));

    return true;
}

bool check_accretion_disk_intersection(Vec3 pos, Vec3 dir, BlackHoleParams params, 
                                        double* intersection_distance, Vec3* intersection_point){
        //Accretion disk is in XY plane
        if(fabs(dir.z) < 1e-10) return false;

        //distance to plane intersection
        double t = -pos.z / dir.z;
        if(t < 0) return false;

        //Intersection point
        Vec3 intersect = vec3_add(pos, vec3_scale(dir, t));
        double distance_from_centre = sqrt(intersect.x * intersect.x + intersect.y * intersect.y);

        //check if the point is within the accetion disk bounds
        if(distance_from_centre >= params.accretion_disk_inner_radius &&
            distance_from_centre <= params.accretion_disk_outer_radius &&
            fabs(intersect.z <= params.accretion_disk_thickness)){
                *intersection_distance = t;
                *intersection_point = intersect;
                return true;
        }
        return false;    
    }