#version 410 core

out vec4 FragColor;
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
uniform sampler2D u_skybox;

// Ray state structure (simulated with separate variables)
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

vec3 vec3_normalize_safe(vec3 v) {
    float len = length(v);
    return len > 1e-10 ? v / len : vec3(0.0, 0.0, 1.0);
}

// Sample the spherical skybox using ray direction
vec3 sample_skybox(vec3 direction) {
    // Normalize direction vector\n"
    vec3 dir = normalize(direction);
   
    // Convert to spherical coordinates\n"
    // theta: azimuthal angle (0 to 2π) - horizontal rotation\n"
    // phi: polar angle (0 to π) - vertical angle from north pole\n"
  
    float theta = atan(dir.z, dir.x) + 3.14159265; // [0, 2π]\n"
    float phi = acos(clamp(dir.y, -1.0, 1.0));      // [0, π]\n"

    // Convert to UV coordinates for equirectangular mapping\n"
    // U maps theta (horizontal): 0 to 1 across image width\n"
    // V maps phi (vertical): 0 to 1 from top to bottom\n"
    vec2 uv = vec2(
        theta / (2.0 * 3.14159265),  // u = θ/(2π)\n"
        phi / 3.14159265             // v = φ/π\n"
    );
   
    // Sample the skybox texture\n"
   return texture(u_skybox, uv).rgb;
}

// Trace a single ray step through curved spacetime
bool trace_rayStep(inout RayState ray) {
    float r = vec3_length_safe(ray.position);
    
    // Check if ray crossed event horizon
    if (r <= u_schwarzschild_radius * 1.01) {
        ray.intensity = 0.0;
        return false;
    }
    
    // Direction from position to black hole center
    vec3 radial_dir = -ray.position / r;
    
    // Calculate impact parameter
    vec3 h = cross(ray.position, ray.direction);
    float impact_parameter = length(h);
    float rs = u_schwarzschild_radius;
    
    float critical_impact = 2.6 * rs;
    float impact_proximity = abs(impact_parameter - critical_impact);
    
    // Calculate gravitational deflection strength
    float deflection_factor = 1.5 * rs / (r * r);
    
    if (impact_proximity < 0.5) {
        deflection_factor = min(deflection_factor, 0.2 / (impact_proximity + 0.1));
    }
    
    // Enhanced deflection near photon sphere
    if (r < 3.0 * rs) {
        float photon_sphere_factor = 1.0 + 2.0 * rs / (r - 1.5 * rs);
        deflection_factor *= min(photon_sphere_factor, 5.0);
    }
    
    // Apply gravitational deflection
    vec3 deflection = radial_dir * deflection_factor;
    ray.direction = vec3_normalize_safe(ray.direction + deflection);
    
    // Calculate redshift
    float redshift_factor = sqrt(max(0.1, 1.0 - u_schwarzschild_radius / r));
    ray.redshift *= redshift_factor;
    
    // Intensity damping
    float radial_component = abs(dot(ray.direction, radial_dir));
    float tangential_damping = 1.0 - 0.1 * (1.0 - radial_component);
    ray.intensity *= max(0.97, tangential_damping) * (1.0 - u_dt / (r * r + 5.0));
    
    // Move ray forward with adaptive step size
    float adaptive_dt = u_dt * min(1.0, r / (5.0 * u_schwarzschild_radius));
    ray.position += ray.direction * adaptive_dt;
    
    return true;
}

// Main raytracing function
vec3 trace_black_hole_ray(vec3 ray_origin, vec3 ray_dir) {
    vec3 bg_color = vec3(0.0, 2.0/255.0, 8.0/255.0); // Dark blue background
    
    RayState ray;
    ray.position = ray_origin;
    ray.direction = ray_dir;
    ray.redshift = 1.0;
    ray.intensity = 1.0;
    
    // Calculate impact parameter for ring effects
    vec3 h = cross(ray_origin, ray_dir);
    float impact_params = length(h);
    float rs = u_schwarzschild_radius;
    
    float photon_ring_factor = 0.0;
    float einstein_ring_factor = 0.0;
    
    float accumulated_opacity = 0.0;
    vec3 accum_color = bg_color;
    
    // Ray tracing loop
    for (int step = 0; step < u_max_steps && accumulated_opacity < 0.99; step++) {
        if (!trace_rayStep(ray)) {
            // Ray crossed event horizon
            break;
        }
        
        if (length(ray.position) > u_observer_distance * 10.0) {
            // Ray moved too far away
            break;
        }
    }
    
    // Calculate ring effects with smooth transitions
    
    // Einstein ring effect - smooth falloff
    if (abs(impact_params - 2.8 * rs) < 0.3 * rs) {
        float dist = abs(impact_params - 2.8 * rs) / (0.3 * rs);
        einstein_ring_factor = max(0.0, 1.0 - dist * dist);
    }
    
    // Photon ring effect - smooth falloff  
    if (abs(impact_params - 2.6 * rs) < 0.12 * rs) {
        float dist = abs(impact_params - 2.6 * rs) / (0.12 * rs);
        photon_ring_factor = max(0.0, 1.0 - dist * dist);
    }
    
    vec3 final_color = accum_color;
    
    // Apply Einstein ring effect - blueish glow
    if (einstein_ring_factor > 0.01) {
        final_color.r = min(final_color.r + 60.0/255.0 * einstein_ring_factor, 1.0);
        final_color.g = min(final_color.g + 110.0/255.0 * einstein_ring_factor, 1.0);
        final_color.b = min(final_color.b + 170.0/255.0 * einstein_ring_factor, 1.0);
    }
    
    // Apply photon ring effect - bright yellow ring
    if (photon_ring_factor > 0.01) {
        final_color.r = min(final_color.r + 180.0/255.0 * photon_ring_factor, 1.0);
        final_color.g = min(final_color.g + 170.0/255.0 * photon_ring_factor, 1.0);
        final_color.b = min(final_color.b + 120.0/255.0 * photon_ring_factor, 1.0);
    }
    
    return final_color;
}

void main() {
    // Calculate screen coordinates (from texture coordinates)
    vec2 screen_pos = TexCoord * 2.0 - 1.0; // Convert from [0,1] to [-1,1]
    screen_pos.y = -screen_pos.y; // Flip Y coordinate
    
    // Calculate ray direction using camera parameters
    float scale = tan(u_fov / 2.0);
    float screen_x = screen_pos.x * u_aspect * scale;
    float screen_y = screen_pos.y * scale;
    
    // Construct ray direction in world space
    vec3 ray_dir = normalize(
        u_cam_right * screen_x + 
        u_cam_up * screen_y + 
        u_cam_forward
    );
    
    // Trace the ray
    vec3 color = trace_black_hole_ray(u_cam_pos, ray_dir);
    
    FragColor = vec4(color, 1.0);
}