#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL.h>
#include "raytracer.h"
#include "camera.h"

//Function to calculate the ray deflection in a simplified Schwarzschild metric
bool trace_rayStep(Vec3* pos, Vec3* dir, BlackHoleParams params, double step_size){
    
    //Distance from the black hold centre
    Vec3 to_centre = vec3_scale(*pos, -1.0);
    double r = vec3_length(to_centre);

    //check if the ray has crossed the EH
    if(r <= params.schwarzschild_radius) return false; //return false if it has

    //Calculate the gravitational lensing/deflection
    //Simplfied approx from AI for geodesic equation
    double deflection_factor = 2.5 * params.schwarzschild_radius / (r * r);

    //direction from pos to the black hole centre
    Vec3 radial_dir = vec3_normalise(to_centre);

    //component of velocity perpendicular to radial direction
    double parallel_component = vec_dot(*dir, radial_dir);
    Vec3 parallel_vector = vec3_scale(radial_dir, parallel_component);
    Vec3 perpendicular_vector = vec3_sub(*dir, parallel_vector);

    //Apply the gravitational deflection to the ray
    Vec3 deflection_perpendicular = vec3_scale(perpendicular_vector, 1.0 - deflection_factor);
    
    *dir = vec3_normalise(vec3_add(parallel_vector, deflection_perpendicular));

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
            fabs(intersect.z) <= params.accretion_disk_thickness){
                *intersection_distance = t;
                *intersection_point = intersect;
                return true;
        }
        return false;    
    }
//calculate the colour based on the temp/brightness of the disk
Uint32 cal_accretion_disk_colour(Vec3 intersection_point, BlackHoleParams params){
    //distance from the centre
    double r = sqrt(intersection_point.x * intersection_point.x + 
                    intersection_point.y * intersection_point.y);

    //Cal the angualar velocity of the disk(simple model)
    //double omega = sqrt(params.mass / (r * r * r)); //this is never used(in this simplfied model)

    //Calculate the reativistic effects(simplified)
    //Doppler effect based on whether material is moving away to towards the observer
    double angle = atan2(intersection_point.y, intersection_point.x);

    double doppler_factor;

    if(sin(angle) < 0){
        doppler_factor = 1.5;
    }else{
        doppler_factor = 0.7;
    }

    doppler_factor += (1.0 + 0.2 * cos(angle));

    double base_intensity = 3.0 * params.schwarzschild_radius / r;
    double intensity = base_intensity * doppler_factor;
    
    //Clamp the intensity to a max value
    if(intensity > 1.0) intensity = 1.0;
    
    //Convert to a colour value
    Uint32 colour = SDL_MapRGB(SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888), 
        (Uint8)(intensity * 200), 
        (Uint8)(intensity * 130), 
        (Uint8)(intensity * 240));
    return colour;
}

Uint32 trace_black_hole_ray(Vec3 ray_origin, Vec3 ray_dir, BlackHoleParams Params){
    Vec3 pos = ray_origin;
    Vec3 dir = ray_dir;
    
    //BG colour
    Uint32 colour = 0x001020; //default blue to help with debuging

    //RT parameters
    double base_step_size = Params.schwarzschild_radius * 0.01; //step size for ray tracing
    int max_steps = 1000; //max number of steps
    double max_distance = Params.observer_distance * 10;
    double distance_travelled = 0;

    for (int step = 0; step < max_steps; step++){
        double intersection_distance = 0.0;
        Vec3 intersection_point = {0.0, 0.0, 0.0};
        
        if(check_accretion_disk_intersection(pos, dir, Params, &intersection_distance, &intersection_point)){
            //Calculate the colour based on the intersection point
            colour = cal_accretion_disk_colour(intersection_point, Params);
            break;
        }
        // Adaptive step size: shrinks near BH to capture lensing, grows far away for speed
        double r = vec3_length(pos);
        double epsilon = 1e-6;
        double step_size = base_step_size * (1.0 + (Params.schwarzschild_radius / (r + epsilon)) * (Params.schwarzschild_radius / (r + epsilon)));
        
        //Move ray forward and apply gravitional lensing
        if(!trace_rayStep(&pos, &dir, Params, step_size)){
            //Ray has crossed the event horizon
            colour = 0x000000; //black hole colour
            break;
        }
        //Check if the ray has moved too far
        distance_travelled += step_size;
        if (distance_travelled > max_distance) {
            // Simple starfield effect
            double d = fabs(dir.y) * 0.7 + fabs(dir.x * dir.z) * 0.3;
            unsigned char intensity = (unsigned char)(d * 64);
            colour = (intensity << 16) | (intensity << 8) | intensity;
            break;
        }
        //Check if the ray has moved too far
        // if(vec3_length(vec3_sub(pos, ray_origin)) > max_distance){
        //     //Ray has moved too far
        //     colour = 0x000000; //black hole colour
        //     break;
        // }
    }
    return colour;
}

//Ray tracing function for the black hole
//This function will be called from the main function
void raytrace_blackhole(BlackHoleParams params, SDL_Surface* surface){
    int width = surface->w;
    int height = surface->h;
    //setup camera
    double aspect_ratio = (double)width / (double)height;
    double fov = 50.0 * M_PI / 180.0; //Field of view in radians

    //fixed position camera
    Vec3 cam_pos = {0.0,0.0, -params.observer_distance};
    Vec3 cam_target = {0.0,0.0,0.0};
    Vec3 cam_up = {0.0, 1.0, 0.0};

    Camera cam = make_camera(cam_pos, cam_target, cam_up, fov, aspect_ratio);

    // Print camera details
    // printf("Camera Position: (%f, %f, %f)\n", cam.position.x, cam.position.y, cam.position.z);
    // printf("Camera Forward: (%f, %f, %f)\n", cam.forward.x, cam.forward.y, cam.forward.z);
    // printf("Camera Up: (%f, %f, %f)\n", cam.up.x, cam.up.y, cam.up.z);
    // printf("Camera Right: (%f, %f, %f)\n", cam.right.x, cam.right.y, cam.right.z);

    double scale = tan(cam.fov / 2.0); // precompute once

    SDL_LockSurface(surface);
    Uint32* pixels = (Uint32*)surface->pixels;

    memset(pixels, 0, width * height * sizeof(Uint32)); //clear the surface
    
    //Loop through each pixel in the surface    
    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            //Calculate the ray direction
            double screen_x = (2.0 * (x + 0.5) / width - 1.0) * cam.aspect * scale;
            double screen_y = (1.0 - 2.0 * (y + 0.5) / height) * scale;

            //Ray origin (comes from the back towards the observer)
            Vec3 ray_origin = cam.position; //ray origin

            //Ray direction
            Vec3 ray_dir = vec3_normalise(vec3_add3(
                vec3_scale(cam.right, screen_x),
                vec3_scale(cam.up, screen_y),
                cam.forward
            ));
            
            //Trace the ray
            Uint32 colour = trace_black_hole_ray(ray_origin, ray_dir, params);
            
            //Set the pixel colour
            pixels[y * width + x] = colour;
        }
    }
    SDL_UnlockSurface(surface);
}