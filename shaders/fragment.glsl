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

vec3 vec3_normalize_safe(vec3 v) {
    float len = length(v);
    return len > 1e-10 ? v / len : vec3(0.0, 0.0, 1.0);
}

// Sample the spherical skybox using ray direction
vec3 sample_skybox(vec3 direction) {
    // Debug: Return pure red to verify function is being called
    //return vec3(1.0, 0.0, 0.0);//this is actually happening, or something is overiding it.
    return texture(u_skybox, normalize(direction)).rgb;
}

// Trace a single ray step through curved spacetime
bool trace_rayStep(inout RayState ray) {
    float r = vec3_length_safe(ray.position);
    
    // Check if ray crossed event horizon
    if (r <= u_schwarzschild_radius * 1.01) {
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
        deflection_factor *= min(photon_sphere_factor, 5.0); //not as numerical robustness as it could be
    }
    
    // Apply gravitational deflection
    vec3 deflection = radial_dir * deflection_factor;
    ray.direction = vec3_normalize_safe(ray.direction + deflection);
    
    // Calculate redshift
    float redshift_factor = sqrt(max(0.1, 1.0 - u_schwarzschild_radius / r));
    ray.redshift *= redshift_factor; //nothing is actually done with this red shifting.
    
    // Intensity damping
    float radial_component = abs(dot(ray.direction, radial_dir));
    float tangential_damping = 1.0 - 0.1 * (1.0 - radial_component);
    ray.intensity *= max(0.97, tangential_damping) * (1.0 - u_dt / (r * r + 5.0));
    
    // Move ray forward with adaptive step size
    float adaptive_dt = u_dt * min(1.0, r / (5.0 * u_schwarzschild_radius));
    ray.position += ray.direction * adaptive_dt; //this is basicly Eular intergration and should be upgraded to RK4 at some point
    
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
    
    // Ray tracing loop **Keeping in the accumulated_opacity for now, May be needed later and doesn't affect anything so far
    for (int step = 0; step < u_max_steps; step++) {
        
        if (!trace_rayStep(ray)) {
            // Ray crossed event horizon
            return vec3(0.0, 1.0, 0.0);//debug colour for now

        }//this isn't being correctly triggered or the skybox is overriding it as a layer 
        
        final_direction = ray.direction;
        
        if (length(ray.position) > u_observer_distance * 10.0) {
            // Ray moved too far away
            break;
        }
    }
    
    // Sample skybox with the final (potentially bent) ray direction
    vec3 skybox_colour = sample_skybox(final_direction);
    
    // Apply redshift and intensity effects to skybox
    //skybox_colour *= ray.redshift * ray.intensity;//this line is the issue

    // Calculate ring effects with smooth transitions
    
    // Einstein ring effect - smooth falloff
    if (abs(impact_params - 2.8 * rs) < 0.3 * rs) {
        float dist = abs(impact_params - 2.8 * rs) / (0.3 * rs);
        einstein_ring_factor = max(0.0, 1.0 - dist * dist);
    }
    
    // Photon ring effect - smooth falloff  
    if (abs(impact_params - 2.6 * rs) < 0.12 * rs) {
        float dist = abs(impact_params - 2.6 * rs) / (0.12 * rs);
        //photon_ring_factor = max(0.0, 1.0 - dist * dist);
        photon_ring_factor = exp(-pow(dist, 2.0) * 5.0); //Gaussian instead
    }
    
    // Apply Einstein ring effect - blueish glow **These two are adding ontop of eachother meaning that the colours are wrong.
    if (einstein_ring_factor > 0.01) {
        skybox_colour += vec3(0.235, 0.431, 0.667) * einstein_ring_factor; // Light blue
    }
    
    // Apply photon ring effect - bright yellow ring
    if (photon_ring_factor > 0.01) {
       skybox_colour += vec3(0.706, 0.667, 0.471) * photon_ring_factor; // Bright yellow
    }
    
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