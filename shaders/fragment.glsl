#version 410 core

out vec4 FragColour;
in vec2 TexCoord;

// Black hole parameters
uniform float u_mass;
uniform float u_schwarzschild_radius;
uniform float u_observer_distance;
uniform float u_dt;
uniform int u_max_steps;

// Camera parameters
uniform vec3 u_cam_pos;
uniform vec3 u_cam_forward;
uniform vec3 u_cam_up;
uniform vec3 u_cam_right;
uniform float u_fov;
uniform float u_aspect;

// Screen parameters
uniform vec2 u_resolution;

// Skybox texture
uniform samplerCube u_skybox;

// Ray state structure (simulated with separate variables)
struct RayState {
    vec3 position;
    vec3 direction;
    float redshift;
    float intensity; //this is redudant atm with no disk
};

// Vector math helpers
float vec3_length_safe(vec3 v) {
    return max(length(v), 1e-10);
}

vec3 vec3_normalise_safe(vec3 v) {
    float len = length(v);
    return len > 1e-10 ? v / len : vec3(0.0, 0.0, 1.0);
}

// Sample the spherical skybox using ray direction
vec3 sample_skybox(vec3 direction) {
    // Debug: Return pure red to verify function is being called
    //return vec3(1.0, 0.0, 0.0);//this is actually happening, or something is overiding it.
    return texture(u_skybox, normalize(direction)).rgb;
}

// Schwarzschild acceleration (geodesic equation)
vec3 schwarzschild_acceleration(vec3 pos, vec3 vel) {
      float r = vec3_length_safe(pos);
    float rs = u_schwarzschild_radius;
    
    if (r <= rs * 1.01) {
        return vec3(0.0); // Near singularity
    }
    
    vec3 r_hat = pos / r;
    
    // Schwarzschild metric acceleration
    float factor = -0.5 * rs / (r * r);
    
    // Radial component
    float vr = dot(vel, r_hat);
    vec3 radial_accel = factor * (1.0 / (1.0 - rs/r)) * r_hat;
    
    // Tangential component  
    vec3 tangential_vel = vel - vr * r_hat;
    vec3 tangential_accel = factor * (2.0 / r) * tangential_vel;
    
    return radial_accel + tangential_accel;
}

// Trace a single ray step through curved spacetime
bool trace_rayStep(inout RayState ray) {
    float r = vec3_length_safe(ray.position);
    float rs = u_schwarzschild_radius;
    
    // Check termination conditions first
    if (r <= 1.5 * rs * 1.01) {
        return false;
    }
    
    // Convert to spherical-like coordinates for geodesic calculation
    vec3 pos = ray.position;
    vec3 vel = ray.direction;
    
    // RK4 step
    vec3 k1_pos = vel;
    vec3 k1_vel = schwarzschild_acceleration(pos, vel);
    
    vec3 k2_pos = vel + 0.5 * u_dt * k1_vel;
    vec3 k2_vel = schwarzschild_acceleration(pos + 0.5 * u_dt * k1_pos, vel + 0.5 * u_dt * k1_vel);
    
    vec3 k3_pos = vel + 0.5 * u_dt * k2_vel;
    vec3 k3_vel = schwarzschild_acceleration(pos + 0.5 * u_dt * k2_pos, vel + 0.5 * u_dt * k2_vel);
    
    vec3 k4_pos = vel + u_dt * k3_vel;
    vec3 k4_vel = schwarzschild_acceleration(pos + u_dt * k3_pos, vel + u_dt * k3_vel);
    
    // Update position and velocity
    ray.position += (u_dt / 6.0) * (k1_pos + 2.0 * k2_pos + 2.0 * k3_pos + k4_pos);
    ray.direction += (u_dt / 6.0) * (k1_vel + 2.0 * k2_vel + 2.0 * k3_vel + k4_vel);
    
    // Normalise direction to maintain unit vector
    ray.direction = vec3_normalise_safe(ray.direction);
    
    // Update redshift and intensity (simplified)
    float new_r = vec3_length_safe(ray.position);
    float redshift_factor = sqrt(max(0.1, 1.0 - rs / new_r));
    ray.redshift *= redshift_factor;
    ray.intensity *= 0.995; // Small decay
    
    return true;
}
// Main raytracing function
vec3 trace_black_hole_ray(vec3 ray_origin, vec3 ray_dir) {
    // Debug: Return pure green to verify function is being called
    //return vec3(0.0, 1.0, 0.0);
    
    RayState ray;
    ray.position = ray_origin;
    ray.direction = ray_dir;
    ray.redshift = 1.0;
    ray.intensity = 1.0;
    
    // Store initial ray direction for skybox sampling if ray escapes
    vec3 final_direction = ray_dir;
    float accumulated_opacity = 0.8; //this is a left over from the accrention disk, I don't think it's required anymore in calculating the EH colour
    
    // Calculate impact parameter for ring effects
    vec3 h = cross(ray_origin, ray_dir);
    float impact_params = length(h);
    float rs = u_schwarzschild_radius;

    float photon_ring_factor = 0.0;
    float einstein_ring_factor = 0.0;
    
    // Ray tracing loop
    for (int step = 0; step < u_max_steps; step++) {

        if (!trace_rayStep(ray)) {
            // Ray crossed event horizon shadow
            return vec3(0.0, 0.0, 0.0);
        }
        
        final_direction = ray.direction;
         
        if (length(ray.position) > u_observer_distance * 10) {
            // Ray moved too far away
            break;
        }
    }

    // Sample skybox with the final (potentially bent) ray direction
    vec3 skybox_colour = sample_skybox(final_direction);

    // Calculate ring effects with smooth transitions and redshift
    float doppler_shift = ray.redshift;

    // // Einstein ring effect - smooth falloff
    // if (abs(impact_params - 2.8 * rs) < 0.3 * rs) {
    //     float dist = abs(impact_params - 2.8 * rs) / (0.3 * rs);
    //     einstein_ring_factor = max(0.0, 1.0 - dist * dist);
    // }
    
    // Photon ring effect - smooth falloff  
    if (abs(impact_params - 2.6 * rs) < 0.35 * rs) { //the 0.2 is the width of the ring, with the side getting pregressivly redder further out
        float dist = abs(impact_params - 2.6 * rs) / (0.35 * rs);
        
        // Combine multiple smooth functions for very gradual falloff **IDK if this really is right but it's smoother thans for sure.
        float smooth_factor1 = 1.0 - smoothstep(0.0, 1.0, dist);           // Smooth hermite
        float smooth_factor2 = exp(-pow(dist, 1.5) * 2.0);                  // Gentler exponential
        float smooth_factor3 = 1.0 / (1.0 + pow(dist * 3.0, 4.0));         // Rational function
        
        // Blend the smooth functions for ultra-smooth transition
        photon_ring_factor = mix(smooth_factor1, smooth_factor2, 0.4);
        photon_ring_factor = mix(photon_ring_factor, smooth_factor3, 0.3);
        
        // Additional outer glow for even smoother edges
        float glow_dist = abs(impact_params - 2.6 * rs) / (0.2 * rs);
        if (glow_dist < 1.0) {
            float glow_factor = pow(1.0 - glow_dist, 3.0) * 0.3;
            photon_ring_factor = max(photon_ring_factor, glow_factor);
        }
    }
    
    // Apply Einstein ring effect - blueish glow 
    //**These two are adding ontop of eachother meaning that the colours are wrong.
    // if (einstein_ring_factor > 0.01) {
    //     skybox_colour += vec3(0.235, 0.431, 0.667) * einstein_ring_factor; // Light blue
    // }
    
    // // Apply photon ring effect 
    if (photon_ring_factor > 0.01) {
        /*The previous vversion without the velocity calcculations looked almost the same*/
        
        //Start with the concentrated skybox light at this direction
        vec3 concentrated_light = sample_skybox(final_direction);
        
        // Calculate orbital Doppler effect for photon sphere
        // Photons at r = 1.5 * rs orbit with velocity v = c/sqrt(3) ≈ 0.577c
        float photon_sphere_radius = 1.5 * rs;
        vec3 to_observer = normalize(u_cam_pos - ray.position);
        
        // Calculate orbital velocity direction (tangent to sphere)
        vec3 radial_dir = normalize(ray.position);
        vec3 orbital_velocity_dir = normalize(cross(cross(radial_dir, to_observer), radial_dir));
        
        // Orbital speed at photon sphere (relativistic)
        float orbital_speed = 0.577; // v/c = 1/sqrt(3)
        
        // Calculate angle between orbital motion and line of sight
        float cos_angle = dot(orbital_velocity_dir, to_observer);
        
        // Relativistic Doppler formula: f_obs/f_source = sqrt((1-β)/(1+β)) * (1 + β*cos(θ))/(1 - β*cos(θ))
        // Simplified for moderate velocities: factor ≈ 1 + β*cos(θ)
        float orbital_doppler = 1.0 + orbital_speed * cos_angle;
        
        // Combine gravitational redshift with orbital Doppler
        float total_shift = max(0.3, doppler_shift) * orbital_doppler;//doppler_shift * orbital_doppler;
        
        // Apply gravitational lensing amplification
        float lensing_amplification = 3.0 + 5.0 / max(0.05, abs(impact_params - 2.6 * rs));
        concentrated_light *= lensing_amplification;
        
        // Apply spectral shifting based on combined effects
        vec3 shifted_colour = concentrated_light;
        
        if (total_shift > 1.2) { //this is neevr triggered
            // Blueshift - approaching side of orbit
            float blue_factor = min(3.0, total_shift);
            shifted_colour.b *= (1.0 + (blue_factor - 1.0) * 2.0);  
            shifted_colour.g *= (1.0 + (blue_factor - 1.0) * 0.5);  
            shifted_colour.r *= (1.0 - (blue_factor - 1.0) * 0.8);  
        } else {
            // Redshift - receding side of orbit
            float red_factor = max(0.1, total_shift);
            shifted_colour.r *= (1.0 + (1.0 - red_factor) * 3.0);   
            shifted_colour.g *= red_factor;                          
            shifted_colour.b *= (red_factor * red_factor);           
        }
     
        // Apply energy conservation 
        float intensity_factor = pow(max(0.13, total_shift), 1.8); //the top end affects intensity, higher = less intense (bright)
        shifted_colour *= intensity_factor;
        
        // Apply relativistic beaming - approaching side appears brighter
        float beaming_factor = pow(orbital_doppler, 3.0); // Intensity ∝ δ^3 for synchrotron this is the beam around the rim.
        shifted_colour *= beaming_factor;
        
        // Make the effect visible
        skybox_colour += shifted_colour * photon_ring_factor * 0.83;
    }
    
    /* ****DEBUG RINGS**** */
    // float eps = 0.05 * u_schwarzschild_radius;
    // // Event horizon ring (b ≤ r_s = 2M)
    // if (impact_params < 2.0 * u_mass + eps) {
    //     return vec3(1.0, 0.0, 0.0); // Red
    // }

    // //Photon sphere radius (r = 3M, not b)
    // if (abs(impact_params - 3.0 * u_mass) < eps) {
    //     return vec3(0.0, 0.0, 1.0); // Blue
    // }

    // //Shadow edge: b_crit = 3√3 M ≈ 5.196M
    // if (abs(impact_params - 3.0 * sqrt(3.0) * u_mass) < eps) {
    //     return vec3(0.0, 1.0, 0.0); // Green
    // }

    // //Einstein ring region (approx): b ≈ 6.0–7.0 M
    // if (abs(impact_params - 6.0 * u_mass) < eps) {
    //     return vec3(1.0, 1.0, 0.0); // Yellow inside edge of ring
    // }
    // if (abs(impact_params - 7.0 * u_mass) < eps) {
    //     return vec3(1.0, 0.5, 0.0); // Orange outside egde
    // }

    // //Outer escape region (reference): b = 8.0 M
    // if (abs(impact_params - 8.0 * u_mass) < eps) {
    //     return vec3(1.0, 1.0, 1.0); // White
    // }
    return clamp(skybox_colour, 0.0, 1.0);
}//function

void main() {
    // Debug: Return pure blue to verify main() is running
    // FragColour = vec4(0.0, 0.0, 1.0, 1.0);
    // return;

    // Calculate ray direction
    vec2 screen_pos = TexCoord * 2.0 - 1.0;
    screen_pos.y = -screen_pos.y;
    
    float scale = tan(u_fov / 2.0);
    float screen_x = screen_pos.x * u_aspect * scale;
    float screen_y = screen_pos.y * scale;
    
    vec3 ray_dir = normalize(
        u_cam_right * screen_x + 
        u_cam_up * screen_y + 
        u_cam_forward
    );
    
    vec3 colour = trace_black_hole_ray(u_cam_pos, ray_dir);
    
    FragColour = vec4(colour, 1.0);
}