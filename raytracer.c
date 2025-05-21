#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "blackholemath.h"
#include <SDL_surface.h>
#include <SDL.h>
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
    params.disk.inner_radius = 3.83 * params.schwarzschild_radius;
    params.disk.outer_radius = 15 * params.schwarzschild_radius;
    params.disk.thickness = 0.2 * params.schwarzschild_radius;
    params.disk.opacity = 0.7;
    params.disk.temperature_factor = 1.0;

    //intergration params
    params.dt = 0.05 * params.schwarzschild_radius;
    params.max_steps = 1000;

    return params;
}

double calculate_doppler(Vec3 disk_velocity, Vec3 view_direction){
     // Compute velocity component along viewing direction (positive is receding)
    double v = vec3_dot(disk_velocity, view_direction);

    //TODO Apply relativisitic doppler formula: sqrt((1-v/c)/(1+v/c))
    if (v > 0.99) v = 0.99;  // prevent division by zero or near-zero
    if (v < -0.99) v = -0.99;
    return sqrt((1 - v) / (1 + v));
}

Vec3 calculate_orbital_velocity(Vec3 pos, double schwarzschild_radius){
    //Distance from black hole centre
    double r = vec3_length(pos);

    //Orbital velocity (Keplerian Approx)
    double v = sqrt(schwarzschild_radius / (2.0 * r));

    //Get unit vector perpendicular to radius in the xy-plane
    Vec3 radial = {pos.x, pos.y, 0.0};
    radial = vec3_normalise(radial);

    //Orbital velocity is perpendicular to radial direction
    Vec3 velocity = {-radial.y, radial.x, 0.0};

    return vec3_scale(velocity, v);
}

//calculate the colour based on the temp/brightness of the disk
Uint32 cal_accretion_disk_colour(Vec3 position, Vec3 view_direction, BlackHoleParams params){
    static SDL_PixelFormat* format = NULL;
    if(!format){
        format = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
    }

    //distance from the centre in XY plane
    double r_xy = sqrt(position.x * position.x + position.y * position.y);
    if (r_xy < 1e-6) r_xy = 1e-6;
    
    // Calculate temperature/brightness based on distance
    // T ∝ r^(-3/4) for standard accretion disk
    double temp_factor = pow(params.disk.inner_radius / r_xy, 0.75);
    temp_factor *= params.disk.temperature_factor;

    // // Add temperature variation based on angle (creates disk features)
    // double angle = atan2(position.y, position.x);
    // double angle_var = 0.1 * sin(6.0 * angle) + 0.05 * sin(12.0 * angle + 0.7);
    // temp_factor *= (1.0 + angle_var);

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
    double intensity = apparent_temp * 2.0; // Scale for visibility
    intensity = intensity / (1.0 + intensity); // Normalize to [0, 1]
    
    double rvalue, gvalue, bvalue;
    if(apparent_temp < 0.6){
         // Cooler regions (reddish)
        rvalue = fmin(intensity * 2.0, (1.0 + intensity)); //Normalised red
        gvalue = fmin(intensity * 0.7, (1.0 + intensity * 0.65)); //reduced green
        bvalue = fmin(intensity * 0.4, (1.0 + intensity * 1.1)); //minmal blue
    }
    else{
        //Hotter Regions
        double norm_factor = 1.0 / (1.0 + 1.2 * intensity);
        rvalue = intensity * norm_factor;
        gvalue = intensity * norm_factor;
        bvalue = intensity * 1.3 * norm_factor;
    }
    
    // Blueshift effect (hotter regions are more blue/white)
    if(total_shift > 1.0){
        double blue_factor = fmin((total_shift - 1.0) * 5.0, 1.0);
        bvalue = fmin(bvalue + blue_factor * 0.5, 1.0);
        gvalue = fmin(gvalue + blue_factor * 0.3, 1.0);
    }

    if(total_shift < 1.0){
        double red_factor = fmin((1.0 - total_shift) * 5.0, 1.0);
        rvalue = fmin(rvalue + red_factor * 0.4, 1.0);
        gvalue = fmin(gvalue - red_factor * 0.2, fmax(0.0, gvalue));
        bvalue = fmin(bvalue - red_factor * 0.3, fmax(0.0, bvalue));
    }
    //Convert to a colour value
    Uint32 colour = SDL_MapRGB(format, 
         (Uint8)(rvalue * 255.0), //ravlue times white
         (Uint8)(gvalue * 255.0), 
         (Uint8)(bvalue * 255.0));
    return colour;
}

bool accretion_disk_intersection(Vec3 pos, Vec3 dir, BlackHoleParams params, 
                                double* intersection_distance, Vec3* intersection_point, Vec3* normal){
    // For simplicity, we'll check intersection with a cylinder and then verify z bounds
    *intersection_distance = -1.0;
    dir = vec3_normalise(dir);

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

    // 2. Check intersection with the inner/outer cylindrical boundaries
    // For both inner and outer radii
    double radii[2] = {params.disk.inner_radius, params.disk.outer_radius};

    for (int i = 0; i < 2; i++){
        double radius = radii[i];

        //Quadratric equation coefficents for cylinder intersection
        double a = dir_xy.x * dir_xy.x + dir_xy.y * dir_xy.y;
        double b = 2.0 * (pos_xy.x * dir_xy.x + pos_xy.y * dir_xy.y);
        double c = pos_xy.x * pos_xy.x + pos_xy.y * pos_xy.y - radius * radius;

        double discriminant = b * b - 4 * a * c;

        if (discriminant >= 0){
            // We have intersection with the infinite cylinder
            double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
            
            // Use the closest positive intersection
            double t = (t1 > 0) ? t1 : ((t2 > 0)) ? t2 : -1.0;

            if (t > 0){
                Vec3 cylinder_hit = vec3_add(pos, vec3_scale(dir, t));

                // Check if the intersection is within the z-bounds of the disk
                if (fabs(cylinder_hit.z) <= params.disk.thickness) {
                    // If this is closer than any previously found intersection
                    if (*intersection_distance < 0 || t < *intersection_distance) {
                        *intersection_distance = t;
                        *intersection_point = cylinder_hit;
                        
                        // Normal points outward for outer radius, inward for inner radius
                        Vec3 radial = {cylinder_hit.x, cylinder_hit.y, 0.0};
                        radial = vec3_normalise(radial);
                        if (i == 0) { // Inner radius - normal points inward
                            radial = vec3_scale(radial, -1.0);
                        }
                        *normal = radial;
                        return true;
                    }//if4   
                }//if3 
            }//if2
        }//if1
    }//for

    //3 checks of intersection with the top and bottom faces of the disk
    
    double z_faces[2] = {params.disk.thickness, -params.disk.thickness};

    for(int i = 0; i < 2; i++){
        double z_face = z_faces[i];

        //check if the ray is moving toward face
        if((z_face > pos.z && dir.z > 0) || (z_face < pos.z && dir.z < 0)){
            double t = (z_face - pos.z) / dir.z;

            if(t > 0){
                Vec3 face_hit = vec3_add(pos, vec3_scale(dir, t));
                double r = sqrt(face_hit.x * face_hit.x + face_hit.y * face_hit.y);

                if(r >= params.disk.inner_radius && r <=params.disk.outer_radius){
                    if(*intersection_distance < 0 || t < *intersection_distance){
                        *intersection_distance = t;
                        *intersection_point = face_hit;
                        *normal = (Vec3){0,0,(i == 0) ? 1.0 : -1.0};
                        return true;
                    }
                }
            }
        }
    }
    return *intersection_distance > 0;
}//function

//Function to calculate the ray deflection in a simplified Schwarzschild metric
bool trace_rayStep(RayState* ray, BlackHoleParams params){
    
    //Distance from the black hold centre
    double r = vec3_length(ray->position);
    if(r < 1e-6) r = 1e-6;

    //check if the ray has crossed the EH
    if(r <= params.schwarzschild_radius * 1.01) {
        ray->intensity = 0.0; //set intensity to zero
        return false; //return false if it has
    }

    //direction from pos to the black hole centre
    Vec3 radial_dir = vec3_scale(ray->position, -1.0 / r);//check this

    // Calculate the gravitational deflection based on Schwarzschild metric
    // Using a more accurate deflection formula based on GR
    double impact_parameter = vec3_length(vec3_cross(ray->position, ray->direction));
    double rs = params.schwarzschild_radius;

    double critial_impact = 2.6 * rs; //critical impact parameter for photon sphere

    double impact_proxmity = fabs(impact_parameter - critial_impact);
    
    //calculate the gravitational deflection strength
    double deflection_factor = 1.5 * rs / (r * r);

    if(impact_proxmity < 0.5){
        // For rays near critical impact, limit maximum deflection more strictly
        // This prevents the chaotic behavior causing white artifacts
        deflection_factor = fmin(deflection_factor, 0.2 / (impact_proxmity + 0.1));
    }

    // Stronger deflection near the photon sphere (approx. 1.5 * schwarzschild_radius)
    if(r < 3.0 * rs){
        //enhanced deflection near photon sphere
        double photon_sphere_factor = 1.0 + 2.0 * rs / (r - 1.5 * rs);
        deflection_factor *= fmin(photon_sphere_factor, 5.0);//limit the max deflection
    }

    //Apply gravitational deflection - pull direction toward black hole
    Vec3 deflection = vec3_scale(radial_dir, deflection_factor);//this was -deflection
    ray->direction = vec3_normalise(vec3_add(ray->direction, deflection));

    //alculate redshift
    // g = sqrt(1 - 2M/r) for static observer in Schwarzschild metric
    double redshift_factor = sqrt(fmax(0.1, 1.0 - params.schwarzschild_radius / r));
    ray->redshift *= redshift_factor;
    
    //Add stronger intensity damping
    //Dot product between direction and radial vector will be small when moving away from the black hole
    //double radial_component = fabs(vec3_dot(ray->direction, radial_dir));
    double radial_component = fabs(vec3_dot(ray->direction, radial_dir));
    double tangential_damping = 1.0 - 0.1 * (1.0 - radial_component);

    //ray->intensity *= (1.0 - params.dt / (r * r + 1.0)); //attenuate intensity based on distance
    ray->intensity *= fmax(0.97, tangential_damping) * (1.0 - params.dt / (r * r + 5.0));

    double adaptive_dt = params.dt * fmin(1.0, r / (5.0 * params.schwarzschild_radius));
    ray->position = vec3_add(ray->position, vec3_scale(ray->direction, adaptive_dt)); //move the ray forward
    
    return true;
}

Uint32 trace_black_hole_ray(Vec3 ray_origin, Vec3 ray_dir, BlackHoleParams Params){
    static SDL_PixelFormat* format = NULL;
    if (!format) {
        format = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
    }

    //BG colour
   // Uint32 colour = 0x001020;//set to light blue for debbugging
    Uint8 bg_r = 0, bg_g = 2, bg_b = 8;  // Very dark blue background
    Uint32 colour = SDL_MapRGB(format, bg_r, bg_g, bg_b);
    
    RayState ray;
    ray.position = ray_origin;
    ray.direction = ray_dir;
    ray.redshift = 1.0; //initial redshift
    ray.intensity = 1.0; //initial intensity
    
    // Calculate impact parameter for this ray
    Vec3 h = vec3_cross(ray_origin,ray_dir); //Angular Momentum
    double impact_params = vec3_length(h);
    double rs = Params.schwarzschild_radius;

    double photon_ring_factor = 0;//= exp(-pow(fabs(impact_params - 2.6 * rs) / (0.05 * rs), 2));
    double einstein_ring_factor = 0; //= exp(-pow(fabs(impact_params - 2.8 * rs) / (0.15 * rs), 2));

    double accumulated_opacity = 0.0;
    double accum_r = bg_r / 255.0;
    double accum_g = bg_g / 255.0;
    double accum_b = bg_b / 255.0;

    bool hit_disk = false;
    int num_disk_hits = 0;
    //Ray tracing loop
    for (int step = 0; step < Params.max_steps && accumulated_opacity < 0.99; step++){
        double distance = -1.0;
        Vec3 hit_point = {0.0, 0.0, 0.0};
        Vec3 normal = {0.0, 0.0, 0.0};
        
        if(accretion_disk_intersection(ray.position, ray.direction, Params, &distance, &hit_point, &normal)){
            //We hit the disk
            hit_disk = true;
            num_disk_hits++;

            //Move to the hit point
            ray.position = hit_point;

            Uint32 disk_colour = cal_accretion_disk_colour(hit_point, ray.direction, Params);

            //Apply acumlated redshift to the colour
            Uint8 r = ((disk_colour >> 16) & 0xFF);
            Uint8 g = ((disk_colour >> 8) & 0xFF);
            Uint8 b = (disk_colour & 0xFF);

            // Apply redshift with smoother transition
            double redshift_factor = ray.redshift;
           // Apply intensity with smoother function to prevent sudden changes
            double intensity_factor = ray.intensity;

           // Convert to float for better precision during calculations
            double rf = (double)r / 255.0 * redshift_factor * intensity_factor;
            double gf = (double)g / 255.0 * redshift_factor * intensity_factor;
            double bf = (double)b / 255.0 * redshift_factor * intensity_factor;

            double base_opacity = Params.disk.opacity;
            if(num_disk_hits > 1){
                // Apply a smoother transition for the opacity
                base_opacity *= (0.7 / num_disk_hits);
            }

            //calculate this segments contribution based on opacity
            double segment_opacity = base_opacity * (1.0 - accumulated_opacity);
            accumulated_opacity = 1.0 - (1.0 - accumulated_opacity) * (1.0 - segment_opacity);

            //blend accumulated colour with the disk colour
            accum_r = accum_r * (1.0 - segment_opacity) + rf * segment_opacity;
            accum_g = accum_g * (1.0 - segment_opacity) + gf * segment_opacity;
            accum_b = accum_b * (1.0 - segment_opacity) + bf * segment_opacity;
            
            // Calculate reflection angle - simplified for now
            double reflection_factor = 0.3; //reflection factor
            
            // Adjust ray for partial reflection/transmission
            Vec3 reflection = vec3_reflect(ray.direction, normal); //niether of these are used

            //Blend between perfect reflection and transmission
            ray.direction = vec3_normalise(
                vec3_add(
                    vec3_scale(reflection, reflection_factor),
                    vec3_scale(ray.direction, 1.0 - reflection_factor)
                )
            );

            ray.intensity *= (1.0 - Params.disk.opacity); //attenuate intensity based on reflection

            //ray.position = vec3_add(ray.position, vec3_scale(ray.direction, distance)); //move the ray forward
            ray.position = vec3_add(ray.position, vec3_scale(normal, 1e-3));//move ray slightly forward
        }   
        if(!trace_rayStep(&ray, Params)){
            //Ray has crossed the event horizon
            break;
        }

        if(vec3_length(ray.position) > Params.observer_distance * 10){
            //claude has adding starfield here, not adding for now.

            //Ray has moved outside the observer distance
            break;
        }
    }

    // Smoother version of your Einstein ring calculation
    if (fabs(impact_params - 2.8 * rs) < 0.3 * rs) {
        // Smooth falloff based on distance from ideal radius
        double dist = fabs(impact_params - 2.8 * rs) / (0.3 * rs);
        einstein_ring_factor = fmax(0.0, 1.0 - dist * dist); // Quadratic falloff
    }
    
    // Smoother version of your photon ring calculation
    if (fabs(impact_params - 2.6 * rs) < 0.12 * rs) {
        // Smooth falloff based on distance from ideal radius
        double dist = fabs(impact_params - 2.6 * rs) / (0.12 * rs);
        photon_ring_factor = fmax(0.0, 1.0 - dist * dist); // Quadratic falloff
    }

    // Apply ring effects with smooth transitions
    if(!hit_disk){

        Uint8 r = (Uint8)(fmin(accum_r * 255.0, 255.0));
        Uint8 g = (Uint8)(fmin(accum_g * 255.0, 255.0));
        Uint8 b = (Uint8)(fmin(accum_b * 255.0, 255.0));

        //Einstein ring effect - Blueish glow
        if(einstein_ring_factor > 0.01){
            // Light blue glow for Einstein ring
            r = (Uint8)fmin(r + 60 * einstein_ring_factor, 255);
            g = (Uint8)fmin(g + 110 * einstein_ring_factor, 255);
            b = (Uint8)fmin(b + 170 * einstein_ring_factor, 255);
        }
    
        // Apply photon ring effect
        if (photon_ring_factor > 0.01) {
            // Bright yellow ring for photon ring
            r = (Uint8)fmin(r + 180 * photon_ring_factor, 255);
            g = (Uint8)fmin(g + 170 * photon_ring_factor, 255);
            b = (Uint8)fmin(b + 120 * photon_ring_factor, 255);
        }
        // Update the accumulated color values
        accum_r = (double)r / 255.0;
        accum_g = (double)g / 255.0;
        accum_b = (double)b / 255.0;
    }
    // Final color mapping with normalization
    Uint8 final_r = (Uint8)(fmin(accum_r * 255.0, 255.0));
    Uint8 final_g = (Uint8)(fmin(accum_g * 255.0, 255.0));
    Uint8 final_b = (Uint8)(fmin(accum_b * 255.0, 255.0));
    
    colour = SDL_MapRGB(format, final_r, final_g, final_b);
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
    double theta = 30 * M_PI / 180.0; //Camera Angle
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