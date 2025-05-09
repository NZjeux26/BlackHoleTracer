#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL.h>
#include "raytracer.h"
#include "camera.h"

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
    params.disk.thickness = 0.1 * params.schwarzschild_radius;
    params.disk.opacity = 0.8;
    params.disk.temperature_factor = 1.0;

    //intergration params
    params.dt = 0.1 * params.schwarzschild_radius;
    params.max_steps = 2000;

    return params;
}

double calculate_doppler(Vec3 disk_velocity, Vec3 view_direction){
     // Compute velocity component along viewing direction (positive is receding)
    double vel_component = vec_dot(disk_velocity, view_direction);

    //TODO Apply relativisitic doppler formula: sqrt((1-v/c)/(1+v/c))
    // Simplified for low velocities
    return 1.0 / (1.0 + vel_component);
}

Vec3 calculate_orbital_velocity(Vec3 pos, double schwarzschild_radius){
    //Distance from black hole centre
    double r = vec3_length(pos);

    //Orbital velocity (Keplerian Approx)
    double v = sqrt(schwarzschild_radius / (2.0 / r));

    //Get unit vector perpendicular to radius in the xy-plane
    Vec3 radial = {pos.x, pos.y, 0.0};
    radial = vec3_normalise(radial);

    //Orbital velocity is perpendicular to radial direction
    Vec3 velocity = {-radial.y, radial.x, 0.0};

    return vec3_scale(velocity, v);
}

//calculate the colour based on the temp/brightness of the disk
Uint32 cal_accretion_disk_colour(Vec3 position, Vec3 view_direction, BlackHoleParams params){
    //distance from the centre in XY plane
    double r_xy = sqrt(position.x * position.x + position.y * position.y);

    // Calculate temperature/brightness based on distance
    // T ∝ r^(-3/4) for standard accretion disk
    double temp_factor = pow(params.disk.inner_radius / r_xy, 0.75);
    temp_factor *= params.disk.temperature_factor;

    // Calculate orbital velocity and apply Doppler shift
    Vec3 orbital_velocity = calculate_orbital_velocity(position, params.schwarzschild_radius);
    double doppler = calculate_doppler(orbital_velocity, view_direction);

    // Calculate gravitational redshift factor
    // g = sqrt(1 - 2M/r) in Schwarzschild metric
    double grav_shift = sqrt(1.0 - params.schwarzschild_radius / r_xy);
    double total_shift = doppler * grav_shift;

    // Apply redshift to temperature (for visualization)
    double apparent_temp = temp_factor * total_shift;
    
     // Convert temperature to RGB (blackbody approximation)
    // Hotter = more blue/white, cooler = more red/yellow
    double intensity = fmin(apparent_temp * 3.0, 1.0);
    
    //base colours
    double rvalue = fmin(intensity * 1.2, 1.0);
    double gvalue = fmin(intensity * 0.9, 1.0);
    double bvalue = fmin(intensity * 0.6, 1.0);
    
    // Blueshift effect (hotter regions are more blue/white)
    if(total_shift > 1.0){
        double blue_factor = fmin((total_shift - 1.0) * 5.0, 1.0);
        bvalue = fmin(bvalue + blue_factor * 0.4, 1.0);
        gvalue = fmin(gvalue + blue_factor * 0.2, 1.0);
    }
    //Convert to a colour value
    Uint32 colour = SDL_MapRGB(SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888), 
        (Uint8)(intensity * 200), 
        (Uint8)(intensity * 130), 
        (Uint8)(intensity * 240));
    return colour;
}

bool accretion_disk_intersection(Vec3 pos, Vec3 dir, BlackHoleParams params, 
                                double* intersection_distance, Vec3* intersection_point, Vec3* normal){
    // For simplicity, we'll check intersection with a cylinder and then verify z bounds

    //Project Ray onto XY Plane for disk intersection

    Vec3 pos_xy = {pos.x, pos.y, 0.0};
    Vec3 dir_xy = {dir.x, dir.y, 0.0};
    dir_xy = vec3_normalise(dir_xy);

    //handle caes where the ray is parallel to the Z axis
    if (vec3_length(dir_xy) < 1e-10){
        //Ray is vertical, check if we're within the disk radius
        double r = vec3_length(pos_xy);
        if (r >= params.disk.inner_radius && r <= params.disk.outer_radius) {
            //Check if we're above/below the disk and moving toward it
            if((pos.z > 0 && dir.z < 0) || (pos.z < 0 && dir.z > 0)){
                double t = -pos.z / dir.z;
                *intersection_distance = t;
                *intersection_point = vec3_add(pos, vec3_scale(dir, t));
                *normal = (Vec3){0,0,(pos.z > 0) ? 1.0 : - 1.0}; //normal is in the Z direction
                return true; //intersection found
            }
        }
        return false;  
    }
    //General Case: Ray not vertical
    //Solve uadratic equation for intersection with cylinder

    //Check intersection with plane
    if(fabs(dir.z) > 1e-10){
        double t_plane = -pos.z / dir.z;
        if(t_plane > 0){
            Vec3 plane_hit = vec3_add(pos, vec3_scale(dir, t_plane));
            double r = sqrt(plane_hit.x * plane_hit.x + plane_hit.y * plane_hit.y);     
            if(r >= params.disk.inner_radius && r <= params.disk.outer_radius){
                //Check if the ray is moving towards the disk
                if((pos.z > 0 && dir.z < 0) || (pos.z < 0 && dir.z > 0)){
                    *intersection_distance = t_plane;
                    *intersection_point = plane_hit;
                    *normal = (Vec3){0,0,(pos.z > 0) ? 1.0 : -1.0}; //normal is in the Z direction
                    return true; //intersection found
                }
            }
        } 
    }
}

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
    double theta = 60 * M_PI / 180; //Camera Angle
    double r = params.observer_distance;
    
    //Movable camera
    Vec3 cam_pos = {
        0.0,
        r * sin(theta),   // height
        -r * cos(theta)   // distance back
    };
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