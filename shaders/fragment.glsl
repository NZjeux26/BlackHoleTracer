#version 410 core

out vec4 FragColour;
in vec2 TexCoord;

// Black hole parameters
uniform float u_mass;
uniform float u_schwarzschild_radius;
uniform float u_observer_distance;
uniform float u_dt;
uniform int u_max_steps;
uniform float u_spin;  // Kerr spin parameter (a)

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

// Ray state structure
struct RayState {
    vec3 position;
    vec3 direction;
    float redshift;
    float intensity;
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
    return texture(u_skybox, normalize(direction)).rgb;
}

// Kerr metric critical radii calculations
float kerr_outer_horizon(float M, float a) {
    return M + sqrt(M*M - a*a);
}

float kerr_inner_horizon(float M, float a) {
    return M - sqrt(M*M - a*a);
}

float kerr_ergosphere_radius(float M, float a, float theta) {
    float cos_theta = cos(theta);
    return M + sqrt(M*M - a*a * cos_theta * cos_theta);
}

// Kerr photon sphere radii (depends on inclination and direction)
float kerr_photon_sphere_prograde(float M, float a) {
    return 2.0 * M * (1.0 + cos(2.0/3.0 * acos(-a/M)));
}

float kerr_photon_sphere_retrograde(float M, float a) {
    return 2.0 * M * (1.0 + cos(2.0/3.0 * acos(a/M)));
}

// Calculate impact parameters for Kerr geodesics
vec2 calculate_kerr_impact_parameters(vec3 pos, vec3 dir, float a) {
    // Convert to Boyer-Lindquist-like coordinates
    float r = length(pos);
    float cos_theta = pos.z / r;
    float sin_theta = sqrt(max(0.0, 1.0 - cos_theta * cos_theta));
    
    // Angular momentum and Carter constant approximations
    vec3 angular_momentum = cross(pos, dir);
    float L_z = angular_momentum.z;  // z-component of angular momentum
    
    // Carter constant (simplified - exact calculation is complex)
    float carter_Q = dot(angular_momentum, angular_momentum) - L_z * L_z;
    
    return vec2(abs(L_z), sqrt(max(0.0, carter_Q)));
}

// Kerr geodesic acceleration (simplified version)
vec3 kerr_acceleration(vec3 pos, vec3 vel, float M, float a) {
    float r = vec3_length_safe(pos);
    float rs = 2.0 * M;
    
    if (r <= kerr_outer_horizon(M, a) * 1.01) {
        return vec3(0.0); // Near horizon
    }
    
    // Convert to spherical coordinates
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    
    // Kerr metric components (simplified)
    float rho2 = r*r + a*a * cos_theta*cos_theta;
    float delta = r*r - rs*r + a*a;
    
    // Effective potential modifications for spin
    vec3 r_hat = pos / r;
    float vr = dot(vel, r_hat);
    
    // Radial acceleration with Kerr corrections
    float radial_factor = -M / (rho2 * rho2) * (rs * rho2 - 2.0 * M * r*r);
    vec3 radial_accel = radial_factor * r_hat;
    
    // Frame dragging effect
    vec3 spin_axis = vec3(0.0, 0.0, 1.0);
    vec3 frame_drag = 2.0 * M * a * r / (rho2 * rho2) * cross(spin_axis, vel);
    
    // Tangential component with spin coupling
    vec3 tangential_vel = vel - vr * r_hat;
    float tangential_factor = -M * rs / (rho2 * rho2 * r);
    vec3 tangential_accel = tangential_factor * tangential_vel;
    
    return radial_accel + tangential_accel + frame_drag;
}

// Trace a single ray step through Kerr spacetime
bool trace_kerr_rayStep(inout RayState ray, float M, float a) {
    float r = vec3_length_safe(ray.position);
    float r_outer = kerr_outer_horizon(M, a);
    
    // Check termination conditions
    if (r <= r_outer * 1.05) {
        return false;
    }
    
    vec3 pos = ray.position;
    vec3 vel = ray.direction;
    
    // RK4 integration with Kerr geodesics
    vec3 k1_pos = vel;
    vec3 k1_vel = kerr_acceleration(pos, vel, M, a);
    
    vec3 k2_pos = vel + 0.5 * u_dt * k1_vel;
    vec3 k2_vel = kerr_acceleration(pos + 0.5 * u_dt * k1_pos, vel + 0.5 * u_dt * k1_vel, M, a);
    
    vec3 k3_pos = vel + 0.5 * u_dt * k2_vel;
    vec3 k3_vel = kerr_acceleration(pos + 0.5 * u_dt * k2_pos, vel + 0.5 * u_dt * k2_vel, M, a);
    
    vec3 k4_pos = vel + u_dt * k3_vel;
    vec3 k4_vel = kerr_acceleration(pos + u_dt * k3_pos, vel + u_dt * k3_vel, M, a);
    
    // Update position and velocity
    ray.position += (u_dt / 6.0) * (k1_pos + 2.0 * k2_pos + 2.0 * k3_pos + k4_pos);
    ray.direction += (u_dt / 6.0) * (k1_vel + 2.0 * k2_vel + 2.0 * k3_vel + k4_vel);
    
    // Maintain unit vector
    ray.direction = vec3_normalise_safe(ray.direction);
    
    // Kerr redshift calculation (simplified)
    float new_r = vec3_length_safe(ray.position);
    float rho2 = new_r*new_r + a*a * pow(ray.position.z / new_r, 2.0);
    float delta = new_r*new_r - 2.0*M*new_r + a*a;
    
    float redshift_factor = sqrt(max(0.1, delta / rho2));
    ray.redshift *= redshift_factor;
    ray.intensity *= 0.995;
    
    return true;
}

// Main Kerr raytracing function
vec3 trace_kerr_ray(vec3 ray_origin, vec3 ray_dir, float M, float a) {
    RayState ray;
    ray.position = ray_origin;
    ray.direction = ray_dir;
    ray.redshift = 1.0;
    ray.intensity = 1.0;
    
    vec3 final_direction = ray_dir;
    
    // Calculate Kerr impact parameters
    vec2 impact_params = calculate_kerr_impact_parameters(ray_origin, ray_dir, a);
    float L_z = impact_params.x;  // Angular momentum
    float carter_Q = impact_params.y; // Carter constant
    
    // Determine if orbit is prograde or retrograde
    bool is_prograde = dot(cross(ray_origin, ray_dir), vec3(0.0, 0.0, 1.0)) > 0.0;
    
    // Ray tracing loop
    for (int step = 0; step < u_max_steps; step++) {
        if (!trace_kerr_rayStep(ray, M, a)) {
            // Ray crossed event horizon
            return vec3(0.0, 0.0, 0.0);
        }
        
        final_direction = ray.direction;
        
        if (length(ray.position) > u_observer_distance * 10.0) {
            break;
        }
    }
    
    // Sample skybox with final ray direction
    vec3 skybox_colour = sample_skybox(final_direction);
    
    // Calculate Kerr-specific ring effects
    float doppler_shift = ray.redshift;
    
    // Photon sphere effects (asymmetric for Kerr)
    float r_ph_pro = kerr_photon_sphere_prograde(M, a);
    float r_ph_retro = kerr_photon_sphere_retrograde(M, a);
    
    float photon_ring_factor = 0.0;
    float ring_width = 0.4 * M;
    
    if (is_prograde) {
        // Prograde photon ring (closer to black hole, more intense)
        if (abs(L_z - r_ph_pro) < ring_width) {
            float dist = abs(L_z - r_ph_pro) / ring_width;
            photon_ring_factor = exp(-pow(dist * 2.0, 2.0)) * 1.5; // More intense
        }
    } else {
        // Retrograde photon ring (farther out, less intense)
        if (abs(L_z - r_ph_retro) < ring_width) {
            float dist = abs(L_z - r_ph_retro) / ring_width;
            photon_ring_factor = exp(-pow(dist * 2.0, 2.0)) * 0.8; // Less intense
        }
    }
    
    // Apply photon ring effects with Kerr-specific features
    if (photon_ring_factor > 0.01) {
        vec3 concentrated_light = sample_skybox(final_direction);
        
        // Frame dragging Doppler effect
        vec3 to_observer = normalize(u_cam_pos - ray.position);
        vec3 spin_axis = vec3(0.0, 1.0, 0.0);
        
        // Calculate frame dragging velocity
        float r = length(ray.position);
        float omega_frame = 2.0 * M * a * r / (pow(r*r + a*a, 2.0));
        vec3 drag_velocity = omega_frame * cross(spin_axis, normalize(ray.position));
        
        // Combined Doppler effect (orbital + frame dragging)
        float orbital_speed = is_prograde ? 0.6 : 0.4; // Different speeds for pro/retrograde
        vec3 total_velocity = orbital_speed * normalize(cross(spin_axis, to_observer)) + drag_velocity;
        float total_doppler = 1.0 + dot(total_velocity, to_observer);
        
        // Gravitational lensing (stronger for prograde orbits)
        float lensing_factor = is_prograde ? 4.0 : 2.5;
        float lensing_amplification = lensing_factor + 3.0 / max(0.05, abs(L_z - (is_prograde ? r_ph_pro : r_ph_retro)));
        concentrated_light *= lensing_amplification;
        
        // Spectral shifting with asymmetric effects
        vec3 shifted_colour = concentrated_light;
        float combined_shift = doppler_shift * total_doppler;
        
        if (is_prograde) {
            // Prograde: stronger redshift on inner edge, blueshift on outer
            if (L_z < r_ph_pro) {
                // Inner edge - strong redshift
                float red_factor = max(0.2, combined_shift);
                shifted_colour.r *= (1.0 + (1.0 - red_factor) * 4.0);
                shifted_colour.g *= red_factor * 0.7;
                shifted_colour.b *= red_factor * red_factor * 0.5;
            } else {
                // Outer edge - moderate blueshift
                float blue_factor = min(2.0, combined_shift);
                shifted_colour.b *= (1.0 + (blue_factor - 1.0) * 1.5);
                shifted_colour.g *= (1.0 + (blue_factor - 1.0) * 0.3);
            }
        } else {
            // Retrograde: more uniform redshift
            float red_factor = max(0.3, combined_shift);
            shifted_colour.r *= (1.0 + (1.0 - red_factor) * 2.5);
            shifted_colour.g *= red_factor * 0.8;
            shifted_colour.b *= red_factor * red_factor * 0.7;
        }
        
        // Intensity with relativistic beaming
        float beaming_power = is_prograde ? 3.5 : 2.5;
        float intensity_factor = pow(max(0.1, combined_shift), 1.5);
        float beaming_factor = pow(max(0.5, total_doppler), beaming_power);
        
        shifted_colour *= intensity_factor * beaming_factor;
        
        // Apply the ring effect
        float ring_intensity = is_prograde ? 1.2 : 0.9;
        skybox_colour += shifted_colour * photon_ring_factor * ring_intensity;
    }
    
    // Ergosphere effects (subtle glow for rays passing through)
    float r = length(ray.position);
    float theta = acos(clamp(ray.position.z / r, -1.0, 1.0));
    float r_ergo = kerr_ergosphere_radius(M, a, theta);
    
    if (r < r_ergo && r > kerr_outer_horizon(M, a)) {
        // Add subtle frame-dragging glow
        float ergo_factor = (r_ergo - r) / (r_ergo - kerr_outer_horizon(M, a));
        vec3 ergo_color = vec3(0.3, 0.1, 0.4) * ergo_factor * 0.2;
        skybox_colour += ergo_color;
    }
    
    /* DEBUG RINGS FOR KERR CRITICAL SURFACES */
    // float eps = 0.05 * M;
    
    // // Event horizon
    // if (abs(L_z - kerr_outer_horizon(M, a)) < eps) {
    //     return vec3(1.0, 0.0, 0.0); // Red - outer horizon
    // }
    
    // // Prograde photon sphere
    // if (abs(L_z - r_ph_pro) < eps) {
    //     return vec3(0.0, 1.0, 0.0); // Green - prograde photon sphere
    // }
    
    // // Retrograde photon sphere  
    // if (abs(L_z - r_ph_retro) < eps) {
    //     return vec3(0.0, 0.0, 1.0); // Blue - retrograde photon sphere
    // }
    
    // // Ergosphere boundary
    // if (abs(r - r_ergo) < eps) {
    //     return vec3(1.0, 1.0, 0.0); // Yellow - ergosphere
    // }
    
    return clamp(skybox_colour, 0.0, 1.0);
}

void main() {
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
    
    // Use Kerr raytracing with mass and spin
    float M = u_mass;
    float a = u_spin;
    
    vec3 colour = trace_kerr_ray(u_cam_pos, ray_dir, M, a);
    
    FragColour = vec4(colour, 1.0);
}