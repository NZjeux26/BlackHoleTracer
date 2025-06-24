#version 410 core

out vec4 FragColour;
in vec2 TexCoord;
#define PI 3.1415926538

// Black hole parameters
uniform float u_mass;
uniform float u_schwarzschild_radius; //not needed anymore
uniform float u_observer_distance;
uniform float u_dtau;
uniform float u_eps;
uniform int u_max_steps;
uniform float u_spin;  // Kerr spin parameter (a)

// Accretion disk parameters
uniform float u_disk_inner_radius;
uniform float u_disk_outer_radius;
uniform float u_disk_opacity;
uniform float u_disk_temperature_factor;
uniform float u_disk_thickness;

// Enhanced disk parameters
float u_time = 1.3;                     // For animation (not used, set to 0.0)
float u_doppler_factor  = 5.0;           // Doppler factor for relativistic effects
uniform float u_disk_turbulence;         // Turbulence strength
uniform float u_disk_spiral_arms;        // Number of spiral arms
uniform float u_disk_spiral_tightness;   // How tight the spirals are
uniform float u_disk_brightness;         // Overall disk brightness
uniform int u_disk_volume_samples;       // Number of volume samples

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
uniform sampler2D u_disk_texture;  // Accretion disk texture

//Cconstants
float M = u_mass; // Mass of the black hole
float a = u_spin * M; //Spin parameter spin amount * Mass of the black hole

//////////////////////////////////////////////////////////////
// Kerr Metric in Kerr-Schild Coordinates
//////////////////////////////////////////////////////////////

///////////////
// Utility Functions
/////////////////

/*Computes the inverse of a 4x4 matrix*/
mat4 diag(vec4 v) {
    return mat4(v.x,0,0,0, 0,v.y,0,0, 0,0,v.z,0, 0,0,0,v.w);
}

///Normalizes a 4D vector using the metric tensor g, ensuring it has unit length
vec4 unit(vec4 v, mat4 g) {
    float norm2 = dot(g * v, v);
    return (norm2 != 0.0) ? v / sqrt(abs(norm2)) : v;
}

// Enhanced noise functions for procedural disk structure
float hash(float n) {
    return fract(sin(n) * 43758.5453123);
}

float noise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    
    float n = i.x + i.y * 57.0 + 113.0 * i.z;
    return mix(
        mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
            mix(hash(n + 57.0), hash(n + 58.0), f.x), f.y),
        mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
            mix(hash(n + 170.0), hash(n + 171.0), f.x), f.y), f.z);
}

/////////////////////
// Core Physics Functions
/////////////////////

/*Converts 4D Kerr-Schild coordinates to the radial coordinate r used in the Kerr metric*/
float rFromCoords(vec4 pos) {
    vec3 p = pos.yzw;
    float rho2 = dot(p,p) - a*a;
    float r2 = 0.5 * (rho2 + sqrt(rho2*rho2 + 4.0*a*a*p.z*p.z));
    return sqrt(r2);
}
/*Computes the Kerr spacetime metric tensor at a given 4D position,
describing the curvature of spacetime around the rotating black hole*/
mat4 metric(vec4 pos) {
    float r = rFromCoords(pos);
    vec4 k = vec4(-1.0, (r*pos.y - a*pos.z)/(r*r + a*a),
                        (r*pos.z + a*pos.y)/(r*r + a*a),
                         pos.a / r);
    float f = 2.0 * M * r / (r*r + a*a * pos.a*pos.a / (r*r));
    return f * mat4(k.x*k, k.y*k, k.z*k, k.w*k) + diag(vec4(-1,1,1,1));
}
/*Calculates the Hamiltonian (total energy) of a photon at position x with momentum p, 
used for geodesic integration*/
float hamiltonian(vec4 x, vec4 p) {
    return 0.5 * dot(inverse(metric(x)) * p, p);
}

// Enhanced 4-velocity calculation for orbiting disk material
vec4 getDiskFourVelocity(vec4 pos) {
    float r = rFromCoords(pos);
    float r2 = r * r;
    float a2 = a * a;  // Uses your existing a = u_spin * M
    float delta = r2 - 2.0 * M * r + a2;
    float sigma = r2 + a2 * pos.w * pos.w / r2;
    float A = (r2 + a2) * (r2 + a2) - a2 * delta * pos.w * pos.w / r2;
    
    // Frame dragging frequency (depends on your u_spin parameter)
    float omega_drag = 2.0 * M * a * r / A;
    
    // Keplerian frequency
    float omega_k = sqrt(M) / pow(r, 1.5);
    
    // Final orbital frequency (frame dragging dominates as u_spin increases)
    float orbital_freq = omega_k + omega_drag;
    
    // Four-velocity in Kerr-Schild coordinates
    vec4 u_disk = vec4(1.0, 0.0, 0.0, orbital_freq);
    
    // Normalize the four-velocity
    mat4 g = metric(pos);
    return unit(u_disk, g);
}

// Calculate Doppler shift factor
float calculateDopplerShift(vec4 pos, vec4 photon_momentum, vec4 observer_pos) {
    mat4 g = metric(pos);
    
    // Get disk material four-velocity
    vec4 u_disk = getDiskFourVelocity(pos);
    
    // Observer four-velocity (assumed stationary at camera position)
    vec4 u_obs = vec4(1.0, 0.0, 0.0, 0.0);
    u_obs = unit(u_obs, g);
    
    // Calculate relative velocity between disk material and observer
    // Project disk velocity onto the line of sight
    vec3 disk_pos_3d = pos.yzw;
    vec3 obs_pos_3d = observer_pos.yzw;
    vec3 line_of_sight = normalize(obs_pos_3d - disk_pos_3d);
    
    // Disk orbital velocity (tangential to radius)
    float r = rFromCoords(pos);
    vec3 disk_center = disk_pos_3d;
    disk_center.z = 0.0; // Project to disk plane
    vec3 radial_dir = normalize(disk_center);
    vec3 orbital_dir = vec3(-radial_dir.y, radial_dir.x, 0.0); // Perpendicular to radial
    
    // Calculate orbital speed (from four-velocity)
    float orbital_speed = length(u_disk.yzw);
    vec3 velocity_3d = orbital_dir * orbital_speed;
    
    // Doppler factor based on line-of-sight velocity
    float v_los = dot(velocity_3d, line_of_sight);
    
    // Relativistic Doppler formula: f_obs/f_emit = sqrt((1-β)/(1+β)) where β = v/c
    // For small velocities: f_obs/f_emit ≈ 1 + v/c
    float beta = v_los; // v/c (in natural units)
    float doppler_factor;
    
    if (abs(beta) < 0.1) {
        // Non-relativistic approximation
        doppler_factor = 1.0 + beta;
    } else {
        // Full relativistic formula
        doppler_factor = sqrt((1.0 - beta) / (1.0 + beta));
    }
    
    // Apply user-controlled Doppler strength
    return mix(1.0, doppler_factor, u_doppler_factor);
}

// Fractal Brownian Motion for multi-scale turbulence
float fbm(vec3 p, int octaves) {
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;
    
    for (int i = 0; i < octaves; i++) {
        value += amplitude * noise(p * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }
    return value;
}

// Calculate spiral density pattern
float spiralDensity(vec2 pos, float r) {
    float angle = atan(pos.y, pos.x);
    float spiral_angle = u_disk_spiral_arms * (angle - u_disk_spiral_tightness * log(r / u_disk_inner_radius));
    return 0.5 + 0.5 * cos(spiral_angle);
}
/////////////////
// Geodesic Integration Functions
/////////////////

/*Computes the numerical gradient of the Hamiltonian 
using finite differences for the equations of motion*/
vec4 hamiltonianGradient(vec4 x, vec4 p) {
    float eps = u_eps;
    return (vec4(
        hamiltonian(x + vec4(eps,0,0,0), p),
        hamiltonian(x + vec4(0,eps,0,0), p),
        hamiltonian(x + vec4(0,0,eps,0), p),
        hamiltonian(x + vec4(0,0,0,eps), p)
    ) - hamiltonian(x, p)) / eps;
}

/*Performs one step of Hamiltonian integration to advance the photon's 
position and momentum along its geodesic path*/
void transportStep(inout vec4 x, inout vec4 p) {
    float dtau = u_dtau;
    p -= dtau * hamiltonianGradient(x, p);
    x += dtau * inverse(metric(x)) * p;
}

/*Determines when to stop raytracing (either the photon hits the event horizon or escapes to infinity)*/
bool stopCondition(vec4 pos) {
    float r = rFromCoords(pos);
    float horizon = M + sqrt(M*M - a*a);
    return r < horizon || r > u_observer_distance * 10.0;
}

/////////////// 
// Rendering Functions & Coordinate Systems
/////////////////

/*Constructs a tetrad (orthonormal basis) for the observer's frame in Kerr spacetime,
given the observer's position, time direction, aim direction, and vertical direction.
This tetrad is used to convert between 4D spacetime coordinates and 3D spatial coordinates.*/
mat4 tetrad(vec4 x, vec4 time, vec4 aim, vec4 vert) {
    mat4 g = metric(x);
    vec4 E0 = unit(time, g);
    vec4 E1 = unit(aim + dot(g*aim, E0) * E0, g);
    vec4 E3 = unit(vert - dot(g*vert,E1)*E1 + dot(g*vert,E0)*E0, g);
    vec4 E2 = unit(inverse(g) * vec4(
        dot(E0.yzw, cross(E1.yzw, E3.yzw)),
        -dot(E0.zwx, cross(E1.zwx, E3.zwx)),
         dot(E0.wxy, cross(E1.wxy, E3.wxy)),
        -dot(E0.xyz, cross(E1.xyz, E3.xyz))), g);
    mat4 basis;
    basis[0] = E0;
    basis[1] = E1;
    basis[2] = E2;
    basis[3] = E3;
    return basis;
}

// Sample the spherical skybox using ray direction
vec3 sample_skybox(vec3 direction) {
    // Debug: Return pure red to verify function is being called
    //return vec3(1.0, 0.0, 0.0);//this is actually happening, or something is overiding it.
    return texture(u_skybox, normalize(direction)).rgb;
}

vec3 toSRGB(vec3 linear) {
    return pow(linear, vec3(1.0 / 2.2));
}

vec3 toneMap(vec3 colour, float exposure) {
    return vec3(1.0) - exp(-colour * exposure);
}

// Enhanced disk density function with turbulence and spiral structure
float getDiskDensity(vec4 pos) {
    float r = rFromCoords(pos);
    float r_norm = r / u_disk_inner_radius;
    
    // Base density calculation (your existing code)
    float core_density = 1.0 / pow(r_norm, 0.8);
    float mid_density = 0.4 / pow(r_norm, 0.5);
    float outer_density = 0.15 / pow(r_norm, 0.2);
    
    float base_density = core_density + mid_density + outer_density;
    
    // Smooth outer edge fade
    float edge_fade = 1.0 - smoothstep(0.6, 1.0, (r - u_disk_inner_radius) / (u_disk_outer_radius - u_disk_inner_radius));
    base_density *= edge_fade;
    
    // Vertical density profile
    float z_scale = u_disk_thickness * (0.3 + 0.7 * sqrt(r_norm));
    float vertical_density = exp(-0.5 * pow(pos.w / z_scale, 2.0));
    
    // Frame dragging enhancement (depends on u_spin)
    float r2 = r * r;
    float a2 = a * a;  // Uses your existing a = u_spin * M
    
    // Frame dragging factor (stronger with higher u_spin)
    float frame_drag_factor = 1.0 + 0.5 * u_spin * M / r;
    frame_drag_factor = clamp(frame_drag_factor, 0.8, 2.5);
    
    // Spiral structure modified by frame dragging
    vec2 disk_pos = pos.yz;
    
    // Frame dragging affects spiral tightness (more spin = tighter spirals)
    float dragged_spiral_tightness = u_disk_spiral_tightness * (1.0 + u_spin);
    float angle = atan(disk_pos.y, disk_pos.x);
    float enhanced_spiral_angle = u_disk_spiral_arms * (angle - dragged_spiral_tightness * log(r / u_disk_inner_radius));
    float enhanced_spiral = 0.5 + 0.5 * cos(enhanced_spiral_angle);
    
    // Multi-scale turbulence
    vec3 turb_pos = vec3(pos.yz * (0.1 + 0.05 * r_norm), u_time * 0.02);
    float turbulence = fbm(turb_pos, 4);
    
    // Combine all effects
    float structure_modulation = (0.7 + 0.3 * enhanced_spiral) * (0.6 + 0.4 * turbulence * u_disk_turbulence);
    
    // Apply frame dragging to final density
    float final_density = base_density * vertical_density * structure_modulation * frame_drag_factor * u_disk_opacity;
    
    return final_density;
}

// Enhanced temperature calculation with turbulence
float getDiskTemperature(vec4 pos, float density, vec4 observer_pos) {
    float r = rFromCoords(pos);
    float r_norm = r / u_disk_inner_radius;
    
    // Base temperature profile
    float base_temp = u_disk_temperature_factor * pow(r_norm, -0.8);
    
    // Calculate Doppler shift
    float doppler_shift = calculateDopplerShift(pos, vec4(0.0), observer_pos);
    
    // Apply Doppler shift to temperature (frequency shift affects blackbody temperature)
    float doppler_temp = base_temp * doppler_shift;
    
    // Temperature fluctuations from turbulence
    vec3 temp_pos = vec3(pos.yz * 0.05, u_time * 0.01);
    float temp_turbulence = fbm(temp_pos, 3);
    
    // Hot spots correlate with density
    float temp_enhancement = 1.0 + 0.5 * density * temp_turbulence;
    
    // Final temperature with Doppler effect
    float final_temp = doppler_temp * temp_enhancement;
    
    return clamp(final_temp, 200.0, 1e7);
}

// Accurate blackbody colour from temperature in Kelvin (1000K–40000K)
// Returns linear RGB (not gamma corrected)
vec3 blackbodycolour(float temperature) {
    float t = clamp(temperature, 1000.0, 40000.0);
    float x, y;

    // Planckian locus approximation (CIE 1960 UCS to CIE xy)
    if (t <= 4000.0) {
        x = -0.2661239e9 / (t*t*t) - 0.2343580e6 / (t*t) + 0.8776956e3 / t + 0.179910;
    } else {
        x = -3.0258469e9 / (t*t*t) + 2.1070379e6 / (t*t) + 0.2226347e3 / t + 0.240390;
    }

    y = -1.1063814 * x * x * x - 1.34811020 * x * x + 2.18555832 * x - 0.20219683;

    // Convert xyY to XYZ (Y = 1.0)
    float Y = 1.0;
    float X = Y * x / y;
    float Z = Y * (1.0 - x - y) / y;

    // Convert XYZ to linear sRGB
    vec3 rgb;
    rgb.r =  3.2406 * X - 1.5372 * Y - 0.4986 * Z;
    rgb.g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
    rgb.b =  0.0557 * X - 0.2040 * Y + 1.0570 * Z;

    // Clamp to valid range (avoid negative RGBs)
    return clamp(rgb, 0.0, 1.0);
}

// Volumetric ray marching through the disk
vec3 volumetricDiskRender(vec4 start_pos, vec4 ray_dir, vec4 momentum, float max_distance, vec4 observer_pos) {
   // Debug: Return fixed colour to test if function is being called
    //return vec3(0.5, 0.2, 0.8); // Uncomment this line for basic test
    
    vec3 accumulated_colour = vec3(0.0);
    float accumulated_opacity = 0.0;
    
    int samples = max(24, min(u_disk_volume_samples, 128));
    float step_size = max_distance / float(samples);
    
    vec4 current_pos = start_pos;

    for (int i = 0; i < samples; i++) {
        if (accumulated_opacity > 0.98) break;

        float r = rFromCoords(current_pos);
        if (r >= u_disk_inner_radius && r <= u_disk_outer_radius &&
            abs(current_pos.w) < u_disk_thickness * 2.0) {

            // Enhanced density with frame dragging
            float density = getDiskDensity(current_pos);
            
            if (density < 0.0005) {
                current_pos += ray_dir * step_size;
                continue;
            }

            // Enhanced temperature with Doppler shifting
            float temp = getDiskTemperature(current_pos, density, observer_pos);

            // Emission color (linear RGB)
            vec3 emission = blackbodycolour(temp);

            // Get Doppler factor for beaming
            float doppler_factor = calculateDopplerShift(current_pos, vec4(0.0), observer_pos);
            
            // Apply Doppler beaming to intensity (affects brightness)
            float beaming_enhancement = pow(abs(doppler_factor), 2.0); // Relativistic beaming ∝ δ²
            beaming_enhancement = clamp(beaming_enhancement, 0.1, 8.0);
            
            // Intensity scaling with temperature
            float temp_ratio = temp / 5000.0;
            float brightness;

            if (temp_ratio > 2.0) {
                brightness = 4.0 + 2.0 * log(temp_ratio / 2.0);
            } else {
                brightness = pow(temp_ratio, 1.8);
            }

            brightness = clamp(brightness, 0.1, 12.0);
            emission *= brightness * beaming_enhancement;

            // Opacity calculation
            float opacity_per_step = density * step_size * 0.8;
            opacity_per_step = min(opacity_per_step, 0.3);

            float extinction = exp(-accumulated_opacity * 1.5);

            accumulated_colour += emission * opacity_per_step * extinction * u_disk_brightness;
            accumulated_opacity += opacity_per_step * (1.0 - accumulated_opacity * 0.8);
        }

        current_pos += ray_dir * step_size;
    }

    return accumulated_colour;
}

/*Main raytracing function that follows a light ray through curved spacetime 
and returns the colour from the skybox or black for captured rays*/
vec3 trace_kerr_ray(vec3 dir, vec4 camPos, mat4 axes) {
     vec4 pos = camPos;
    vec4 dir4D = -axes[0] + vec4(0.0, dir.x, dir.y, dir.z);
    vec4 p = metric(pos) * dir4D;

    vec3 disk_contribution = vec3(0.0);
    vec4 final_pos;

    bool in_disk_region = false;
    vec4 disk_entry_pos;
    vec4 disk_entry_momentum;

    for (int i = 0; i < u_max_steps; i++) {
        vec4 last_pos = pos;
        transportStep(pos, p);

        float r = rFromCoords(pos);
        bool currently_in_disk = (r >= u_disk_inner_radius &&
                                  r <= u_disk_outer_radius &&
                                  abs(pos.w) < u_disk_thickness);

        if (currently_in_disk && !in_disk_region) {
            in_disk_region = true;
            disk_entry_pos = pos;
            disk_entry_momentum = p;
        } else if (!currently_in_disk && in_disk_region) {
            float disk_distance = length(pos.yzw - disk_entry_pos.yzw);
            if (disk_distance > 0.001 && disk_distance < 1000.0) {
                // Use enhanced volumetric rendering with Doppler effects
                disk_contribution += volumetricDiskRender(
                    disk_entry_pos,
                    normalize(pos - disk_entry_pos),
                    disk_entry_momentum,
                    disk_distance,
                    camPos  // Pass observer position for Doppler calculation
                );
            }
            in_disk_region = false;
        }

        if (stopCondition(pos)) {
            final_pos = pos;
            break;
        }

        final_pos = pos;
    }

    // Handle final disk contribution if still inside
    if (in_disk_region) {
        float disk_distance = length(final_pos.yzw - disk_entry_pos.yzw);
        if (disk_distance > 0.001 && disk_distance < 1000.0) {
            disk_contribution += volumetricDiskRender(
                disk_entry_pos,
                normalize(final_pos - disk_entry_pos),
                disk_entry_momentum,
                disk_distance,
                camPos  // Pass observer position for Doppler calculation
            );
        }
    }

    // Background determination
    float final_r = rFromCoords(final_pos);
    bool captured = final_r < (M + sqrt(M * M - a * a));
    vec3 background_colour = vec3(0.0);

    if (!captured) {
        vec4 out_dir = inverse(metric(final_pos)) * p;
        vec3 cube_dir = normalize(vec3(-out_dir.y, out_dir.w, -out_dir.z));
        background_colour = sample_skybox(cube_dir);
    }

    return background_colour + disk_contribution;
}

void main() {
    // Calculate ray direction
    vec2 screen_pos = TexCoord * 2.0 - 1.0;
    screen_pos.y = -screen_pos.y;
    
    float scale = tan(u_fov * 0.5);
    float screen_x = screen_pos.x * u_aspect * scale;
    float screen_y = screen_pos.y * scale;
    
    vec3 ray_dir = normalize(
        u_cam_right * screen_x + 
        u_cam_up * screen_y + 
        u_cam_forward
    );
    
    // Observer 4-position and initial time direction
    // Note: u_cam_pos is in the format (t, x, y, z) where t is time
    // and (x, y, z) are the spatial coordinates.
    vec4 camPos = vec4(0.0, u_cam_pos.x, u_cam_pos.y, u_cam_pos.z);
    // Map 3D camera vectors to 4D space
    vec4 aim = vec4(0.0, u_cam_forward.x, u_cam_forward.y, u_cam_forward.z);
    vec4 vert = vec4(0.0, u_cam_up.x, u_cam_up.y, u_cam_up.z);
    vec4 timeDir = vec4(1.0, 0.0, 0.0, 0.0);
    
    mat4 camFrame = tetrad(camPos, timeDir, aim, vert);

    //This area is controlling the colour and the background saturation, In reality the brightness of the disk would overwhelm 
    //the background I'm not happy with the saturated HDR look.
    vec3 colour = trace_kerr_ray(ray_dir, camPos, camFrame);
    vec3 mapped = toneMap(colour, 1.0 / 5.0); // <-- Try tuning this number (smaller = darker)
    vec3 final = toSRGB(mapped);
    
    // Colour debugging
    // float temp = mix(1000.0, 40000.0, TexCoord.x); // left to right: 1000K to 40000K
    // vec3 final = blackbodycolour(temp);
    // Temporary test - remove after debugging
    FragColour = vec4(final, 1.0);
}