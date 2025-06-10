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
// Hamilton-Jacobi separation for null geodesics
vec2 calculate_kerr_impact_parameters(vec3 pos, vec3 dir, float a) {
    float r = length(pos);
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    // Normalised position and direction
    vec3 x_hat = normalize(pos);
    vec3 v_hat = normalize(dir);

    // z-component of angular momentum (around spin axis)
    float Lz = r * (pos.x * dir.y - pos.y * dir.x) / r;

    // Estimate p_theta component
    float p_theta = r * (cos_theta * dir.z - sin_theta * (dir.x * pos.x + dir.y * pos.y) / r);

    // Approximate Carter constant Q
    float Q = p_theta * p_theta + pow(cos_theta, 2.0) * (a * a + Lz * Lz / pow(sin_theta, 2.0));

    return vec2(Lz, sqrt(max(Q, 0.0)));  // Return Lz and sqrt(Q)
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

// Helper function to calculate critical impact parameter for Kerr photon sphere
float kerr_critical_impact_parameter(float M, float a, float theta, bool is_prograde) {
    float r_ph = is_prograde ? kerr_photon_sphere_prograde(M, a) : kerr_photon_sphere_retrograde(M, a);
    
    // Simplified approximation - for exact calculation you'd need to solve
    // the full geodesic equations, but this gives good visual results
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    
    // Base impact parameter (would be exact for Schwarzschild)
    float b_base = r_ph * sqrt(27.0) / sqrt(2.0);  // ~3.464 * r_ph for circular photon orbits
    
    // Kerr corrections for spin effects
    float spin_factor = 1.0 + (a / M) * cos_theta * (is_prograde ? -0.3 : 0.3);
    
    return b_base * spin_factor;
}

// Main Kerr raytracing function
vec3 trace_kerr_ray(vec3 ray_origin, vec3 ray_dir, float M, float a) {
    RayState ray;
    ray.position = ray_origin;
    ray.direction = ray_dir;
    ray.redshift = 1.0;
    ray.intensity = 1.0;
    
    vec3 final_direction = ray_dir;
    
    // Calculate impact parameter (perpendicular distance from ray to black hole)
    vec3 ray_to_bh = -ray_origin;  // Vector from ray origin to black hole at origin
    vec3 ray_dir_norm = normalize(ray_dir);
    float impact_param = length(cross(ray_to_bh, ray_dir_norm));
    
    // Calculate angle with respect to spin axis for Kerr corrections
    float theta = acos(clamp(ray_origin.z / length(ray_origin), -1.0, 1.0));
    
    // Determine if orbit is prograde or retrograde based on angular momentum
    vec3 angular_momentum = cross(ray_origin, ray_dir);
    bool is_prograde = angular_momentum.z > 0.0;
    
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
    
    // Calculate critical impact parameters for photon spheres
    float b_crit_pro = kerr_critical_impact_parameter(M, a, theta, true);
    float b_crit_retro = kerr_critical_impact_parameter(M, a, theta, false);
    
    // Choose the appropriate critical impact parameter
    float b_crit = is_prograde ? b_crit_pro : b_crit_retro;
    
    // Photon ring effects based on impact parameter
    float photon_ring_factor = 0.0;
    float ring_width = 0.35 * M * sqrt(27.0) / sqrt(2.0); // Scale ring width with critical impact parameter
    
    // Photon sphere radii - prograde and retrograde
    float r_ph_pro = kerr_photon_sphere_prograde(M, a);
    float r_ph_retro = kerr_photon_sphere_retrograde(M, a);

    float photon_sphere_radius = is_prograde ? r_ph_pro : r_ph_retro;
    if (abs(impact_param - b_crit) < ring_width) {
        float dist = abs(impact_param - b_crit) / ring_width;
        
        // Blend multiple falloffs for smooth ring
        float smooth_factor1 = 1.0 - smoothstep(0.0, 1.0, dist);
        float smooth_factor2 = exp(-pow(dist, 1.5) * 2.0);
        float smooth_factor3 = 1.0 / (1.0 + pow(dist * 3.0, 4.0));
        
        photon_ring_factor = mix(smooth_factor1, smooth_factor2, 0.4);
        photon_ring_factor = mix(photon_ring_factor, smooth_factor3, 0.3);
        
        // Outer glow based on impact parameter
        float glow_dist = abs(impact_param - b_crit) / (0.2 * M);
        if (glow_dist < 1.0) {
            float glow_factor = pow(1.0 - glow_dist, 3.0) * 0.3;
            photon_ring_factor = max(photon_ring_factor, glow_factor);
        }
    }
    
    // // Apply photon ring visual effects
    // if (photon_ring_factor > 0.01) {
    //     vec3 concentrated_light = sample_skybox(final_direction);
    //     vec3 to_observer = normalize(u_cam_pos - ray.position);
    //     vec3 radial_dir = normalize(ray.position);
        
    //     // Orbital velocity direction (perpendicular to radial and toward observer)
    //     vec3 orbital_velocity_dir = normalize(cross(cross(radial_dir, to_observer), radial_dir));
        
    //     // Calculate orbital frequency at photon sphere
    //     float r_ph = is_prograde ? kerr_photon_sphere_prograde(M, a) : kerr_photon_sphere_retrograde(M, a);
    //     float omega = is_prograde
    //         ? 1.0 / (a + pow(r_ph, 1.5) / M)
    //         : 1.0 / (-a + pow(r_ph, 1.5) / M);
        
    //     float orbital_speed = clamp(r_ph * omega, 0.0, 0.9999);  // Avoid v > c
    //     float cos_angle = dot(orbital_velocity_dir, to_observer);
    //     float orbital_doppler = 1.0 + orbital_speed * cos_angle;
        
    //     // Combine with gravitational redshift
    //     float total_shift = max(0.3, ray.redshift) * orbital_doppler;
        
    //     // Lensing amplification (depends on how close to critical impact parameter)
    //     float lensing_amplification = (is_prograde ? 3.5 : 2.2) + 5.0 / max(0.05, abs(impact_param - b_crit));
    //     concentrated_light *= lensing_amplification;
        
    //     vec3 shifted_colour = concentrated_light;
        
    //     // Spectral shifting based on total Doppler effect
    //     if (total_shift > 1.2) {
    //         // Blueshift regime
    //         float blue_factor = min(3.0, total_shift);
    //         shifted_colour.b *= (1.0 + (blue_factor - 1.0) * 2.0);
    //         shifted_colour.g *= (1.0 + (blue_factor - 1.0) * 0.5);
    //         shifted_colour.r *= (1.0 - (blue_factor - 1.0) * 0.8);
    //     } else {
    //         // Redshift regime
    //         float red_factor = max(0.1, total_shift);
    //         shifted_colour.r *= (1.0 + (1.0 - red_factor) * 3.0);
    //         shifted_colour.g *= red_factor;
    //         shifted_colour.b *= red_factor * red_factor;
    //     }
        
    //     // Energy scaling with relativistic beaming
    //     float intensity_factor = pow(max(0.13, total_shift), 1.8);
    //     float beaming_factor = pow(max(0.5, orbital_doppler), is_prograde ? 3.0 : 2.5);
        
    //     shifted_colour *= intensity_factor * beaming_factor;
        
    //     // Add to skybox with asymmetric intensity
    //     skybox_colour += shifted_colour * photon_ring_factor * (is_prograde ? 1.1 : 0.9);
    // }
    
    // Ergosphere effects (subtle glow for rays passing through)
    float r = length(ray.position);
    float current_theta = acos(clamp(ray.position.z / r, -1.0, 1.0));
    float r_ergo = kerr_ergosphere_radius(M, a, current_theta);
    
    if (r < r_ergo && r > kerr_outer_horizon(M, a)) {
        // Add subtle frame-dragging glow
        float ergo_factor = (r_ergo - r) / (r_ergo - kerr_outer_horizon(M, a));
        vec3 ergo_colour = vec3(0.3, 0.1, 0.4) * ergo_factor * 0.2;
        skybox_colour += ergo_colour;
    }
    
    // /* DEBUG RINGS FOR KERR CRITICAL SURFACES */
     //float eps = 0.05 * M;
    
    // // Event horizon rings - outer and inner horizons for Kerr
    // float r_plus = kerr_outer_horizon(M, a);   // r+ = M + sqrt(M^2 - a^2)
    // float r_minus = kerr_inner_horizon(M, a);  // r- = M - sqrt(M^2 - a^2)
    
    // if (impact_param < r_plus + eps) {
    //     return vec3(1.0, 0.0, 0.0); // Red - outer event horizon
    // }
    // if (abs(impact_param - r_minus) < eps && r_minus > 0.0) {
    //     return vec3(0.8, 0.0, 0.0); // Dark red - inner event horizon (Cauchy horizon)
    // }
    
    // // Ergosphere boundary
    // if (abs(impact_param - r_ergo) < eps) {
    //     return vec3(1.0, 0.0, 1.0); // Magenta - ergosphere boundary
    // }
    
    // if (abs(impact_param - r_ph_pro) < eps) {
    //     return vec3(0.0, 0.0, 1.0); // Blue - prograde photon sphere
    // }
    // if (abs(impact_param - r_ph_retro) < eps) {
    //     return vec3(0.0, 0.5, 1.0); // Light blue - retrograde photon sphere
    // }
    
    // // Critical impact parameters - shadow edge
    // if (abs(impact_param - b_crit_pro) < eps) {
    //     return vec3(0.0, 1.0, 0.0); // Green - prograde critical impact parameter (shadow edge)
    // }
    // if (abs(impact_param - b_crit_retro) < eps) {
    //     return vec3(0.5, 1.0, 0.0); // Light green - retrograde critical impact parameter
    // }
    
    // // Einstein ring region approximations for Kerr (these are more complex than Schwarzschild)
    // // For Kerr, Einstein ring depends on source position, but we'll use rough estimates
    // float einstein_ring_inner = 6.0 * M;  // Approximate inner Einstein ring
    // float einstein_ring_outer = 8.0 * M;  // Approximate outer Einstein ring
    
    if (abs(impact_param - einstein_ring_inner) < eps) {
        return vec3(1.0, 1.0, 0.0); // Yellow - inner Einstein ring region
    }
    if (abs(impact_param - einstein_ring_outer) < eps) {
        return vec3(1.0, 0.5, 0.0); // Orange - outer Einstein ring region
    }
    
    // // Stable circular orbit radius (ISCO - Innermost Stable Circular Orbit)
    // // float r_isco = kerr_isco_radius(M, a, is_prograde);
    // // if (abs(impact_param - r_isco) < eps) {
    // //     return vec3(0.0, 1.0, 1.0); // Cyan - ISCO
    // // }
    
    // // Marginally bound orbit
    // // float r_mb = kerr_marginally_bound_radius(M, a, is_prograde);
    // // if (abs(impact_param - r_mb) < eps) {
    // //     return vec3(1.0, 1.0, 1.0); // White - marginally bound orbit
    // // }
    
    // // Additional reference rings
    // if (abs(impact_param - 10.0 * M) < eps) {
    //     return vec3(0.5, 0.5, 0.5); // Gray - 10M reference ring
    // }
    // if (abs(impact_param - 15.0 * M) < eps) {
    //     return vec3(0.3, 0.3, 0.3); // Dark gray - 15M reference ring
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