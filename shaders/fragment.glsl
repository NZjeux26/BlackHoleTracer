#version 410 core

out vec4 FragColour;
in vec2 TexCoord;
#define PI 3.1415926538

// Black hole parameters
uniform float u_mass;
uniform float u_spin;  // Kerr spin parameter (a)
uniform float u_observer_distance;
uniform float u_dtau;
uniform float u_eps;
uniform int u_max_steps;

// SPH particle data
uniform sampler2D u_particle_positions;    // RGBA32F texture: xyz = position, w = mass
uniform sampler2D u_particle_velocities;   // RGBA32F texture: xyz = velocity, w = density
uniform sampler2D u_particle_properties;   // RGBA32F texture: x = temperature, y = pressure, z = smoothing_length, w = flags
uniform sampler2D u_particle_thermal;      // RGBA32F texture: x = thermal_energy, y = radiative_cooling, z = heating_rate, w = unused
uniform int u_particle_count;              // Total number of active particles
uniform int u_particle_texture_size;       // Size of particle texture (sqrt of max particles)

//SPH hash table data
uniform sampler2D u_hash_indices;
uniform sampler2D u_hash_counts;
uniform float u_grid_cell_size;
uniform int u_hash_table_size;

// Accretion disk parameters
uniform float u_disk_inner_radius;
uniform float u_disk_outer_radius;
uniform float u_disk_opacity;
uniform float u_disk_thickness;

// Enhanced disk parameters
uniform float u_doppler_factor;           // Doppler factor for relativistic effects
uniform float u_disk_brightness;         // Overall disk brightness

// Camera parameters
uniform vec3 u_cam_pos;
uniform vec3 u_cam_forward;
uniform vec3 u_cam_up;
uniform vec3 u_cam_right;
uniform float u_fov;
uniform float u_aspect;

// Skybox texture
uniform samplerCube u_skybox;

//Constants
float M = u_mass; // Mass of the black hole
float a = u_spin * M; //Spin parameter spin amount * Mass of the black hole

//////////////////////////////////////////////////////////////
// Kerr Metric in Kerr-Schild Coordinates
//////////////////////////////////////////////////////////////

/*Computes the inverse of a 4x4 matrix*/
mat4 diag(vec4 v) {
    return mat4(v.x,0,0,0, 0,v.y,0,0, 0,0,v.z,0, 0,0,0,v.w);
}

///Normalizes a 4D vector using the metric tensor g, ensuring it has unit length
vec4 unit(vec4 v, mat4 g) {
    float norm2 = dot(g * v, v);
    return (norm2 != 0.0) ? v / sqrt(abs(norm2)) : v;
}

// Hash function matching your CPU implementation
uint hashPosition(vec3 position, float cell_size) {
    ivec3 cell = ivec3(floor(position / cell_size));
    
    // Use the same hash function as your CPU code
    uint hash = uint(cell.x * 73856093) ^ uint(cell.y * 19349663) ^ uint(cell.z * 83492791);
    return hash % uint(u_hash_table_size);
}

// Convert hash index to texture coordinates
vec2 hashIndexToTexCoord(uint hash_index) {
    int hash_tex_size = int(sqrt(float(u_hash_table_size)) + 0.5);
    int x = int(hash_index) % hash_tex_size;
    int y = int(hash_index) / hash_tex_size;
    return vec2(float(x) + 0.5, float(y) + 0.5) / float(hash_tex_size);
}

// Get particles from hash cell
void getParticlesFromHashCell(vec3 position, inout int particle_indices[9], inout int particle_count) {
    uint hash_index = hashPosition(position, u_grid_cell_size);
    vec2 tex_coord = hashIndexToTexCoord(hash_index);
    
    // Get particle count for this hash cell
    vec4 count_data = texture(u_hash_counts, tex_coord);
    int cell_count = int(count_data.r);
    
    // Get particle indices (up to 3 per texel for now)
    vec4 indices_data = texture(u_hash_indices, tex_coord);
    
    for (int i = 0; i < 3 && i < cell_count && particle_count < 9; i++) {
        particle_indices[particle_count] = int(indices_data[i]);
        particle_count++;
    }
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
    //return vec3(1.0, 0.0, 0.0);
    return texture(u_skybox, normalize(direction)).rgb;
}

vec3 toSRGB(vec3 linear) {
    return pow(linear, vec3(1.0 / 2.2));
}

vec3 toneMap(vec3 colour) {
     // ACES approximation by Krzysztof Narkowicz
    return clamp((colour * (2.51 * colour + 0.03)) / (colour * (2.43 * colour + 0.59) + 0.14), 0.0, 1.0);
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

//////////////////////////////////////////////////////////////
// SPH PARTICLE FUNCTIONS
//////////////////////////////////////////////////////////////

// Sample particle data from textures
vec4 getParticlePosition(int index) {
    if (index < 0 || index >= u_particle_count) return vec4(0.0);
    ivec2 coord = ivec2(index % u_particle_texture_size, index / u_particle_texture_size);
    return texelFetch(u_particle_positions, coord, 0);
}

vec4 getParticleVelocity(int index) {
    if (index < 0 || index >= u_particle_count) return vec4(0.0);
    ivec2 coord = ivec2(index % u_particle_texture_size, index / u_particle_texture_size);
    return texelFetch(u_particle_velocities, coord, 0);
}

vec4 getParticleProperties(int index) {
    if (index < 0 || index >= u_particle_count) return vec4(0.0);
    ivec2 coord = ivec2(index % u_particle_texture_size, index / u_particle_texture_size);
    return texelFetch(u_particle_properties, coord, 0);
}

vec4 getParticleThermal(int index) {
    if (index < 0 || index >= u_particle_count) return vec4(0.0);
    ivec2 coord = ivec2(index % u_particle_texture_size, index / u_particle_texture_size);
    return texelFetch(u_particle_thermal, coord, 0);
}

// SPH Wendland C2 kernel (matches your header)
float wendlandC2Kernel(float r, float h) {
    if (r >= h) return 0.0;
    float q = r / h;
    float factor = 1.0 - q;
    
     return factor * factor * factor * factor * (1.0 + 4.0 * q);
}

// Calculate Doppler shift for particle-based rendering
float calculateParticleDopplerShift(vec3 particle_pos, vec3 particle_velocity, vec4 observer_pos) {
    vec3 observer_3d = observer_pos.yzw;
    vec3 line_of_sight = normalize(observer_3d - particle_pos);
    float v_los = dot(particle_velocity, line_of_sight);
    
    // Relativistic Doppler formula
    float beta = v_los; // v/c in natural units
    float doppler_factor;
    
    if (abs(beta) < 0.1) {
        doppler_factor = 1.0 + beta;
    } else {
        doppler_factor = sqrt((1.0 - beta) / (1.0 + beta));
    }
    
    return mix(1.0, doppler_factor, u_doppler_factor);
}

// Calculate particle contribution at a given point
vec3 calculateParticleContribution(vec3 point, int particle_index, vec4 observer_pos) {
     // Debug: Return pure red to verify function is being called
    //return vec3(0.0, 0.0, 1.0);
    
    vec4 pos_mass = getParticlePosition(particle_index);
    vec3 particle_pos = pos_mass.xyz;
    float mass = pos_mass.w;
    
    vec4 vel_density = getParticleVelocity(particle_index);
    vec3 velocity = vel_density.xyz;
    float density = vel_density.w;
    
    vec4 properties = getParticleProperties(particle_index);
    float temperature = properties.x;
    float pressure = properties.y;
    float smoothing_length = properties.z;
    uint flags = uint(properties.w);
    
    vec4 thermal = getParticleThermal(particle_index);
    float thermal_energy = thermal.x;
    
    // Check if point is within particle's influence
    float dist = length(point - particle_pos);
    if (dist >= smoothing_length) {
        return vec3(0.0);
    }
    
    // Calculate kernel weight
    float kernel_weight = wendlandC2Kernel(dist, smoothing_length);
    if (kernel_weight <= 1e-6) {
        return vec3(0.0);
    }
    
    // Check particle flags for rendering
    bool is_active = (flags & 1u) != 0u;
    bool is_hot = (flags & 16u) != 0u;
    bool is_accreting = (flags & 32u) != 0u;
    
    if (!is_active) {
        return vec3(0.0);
    }
    
    // Enhanced temperature based on particle state
    float effective_temp = temperature;
    if (is_hot) {
        effective_temp *= 1.5; // Hot particles are hotter
    }
    if (is_accreting) {
        effective_temp *= 1.2; // Accreting particles have additional heating
    }
    
    // Add thermal energy contribution
    effective_temp += thermal_energy * 100.0; // Scale thermal energy to temperature
    
    // Calculate Doppler shift
    float doppler_factor = calculateParticleDopplerShift(particle_pos, velocity, observer_pos);
    effective_temp *= doppler_factor;
    
    // Blackbody emission
    vec3 emission = blackbodycolour(effective_temp);
    
    // Intensity based on density and kernel weight
    float intensity = density * kernel_weight * mass;
    
    // Apply Doppler beaming
    float beaming_factor = pow(abs(doppler_factor), 2.0);
    beaming_factor = clamp(beaming_factor, 0.1, 8.0);
    
    intensity *= beaming_factor;
    
    // Temperature-based brightness scaling
    float temp_ratio = effective_temp / 5000.0;
    float brightness;
    if (temp_ratio > 2.0) {
        brightness = 4.0 + 2.0 * log(temp_ratio / 2.0);
    } else {
        brightness = pow(temp_ratio, 1.8);
    }
    brightness = clamp(brightness, 0.1, 12.0);
    
    return emission * intensity * brightness * u_disk_brightness * u_disk_opacity;
}

// Optimised particle-based volumetric rendering
vec3 particleVolumetricRender(vec4 start_pos, vec4 ray_dir, float max_distance, vec4 observer_pos, mat4 frame) {
      // Debug: Return pure red to verify function is being called
    //return vec3(0.0, 1.0, 0.0);
    // // this commented out section is a red bebug that allows you to see all particles without any filtering based on temp. 
    // int samples = 8;
    // float step_size = max_distance / float(samples);
    // vec4 current_pos = start_pos;

    // for (int i = 0; i < samples; ++i) {
    //     vec3 sample_point = current_pos.yzw;
        
    //     // Get particles from current hash cell and neighbors
    //     int particle_indices[9]; // Max particles to check
    //     int particle_count = 0;
        
    //     // Check current cell
    //     getParticlesFromHashCell(sample_point, particle_indices, particle_count);
        
    //     // Check neighboring cells (simplified - just check a few key neighbors)
    //     vec3 cell_offset = vec3(u_grid_cell_size * 0.5);
    //     getParticlesFromHashCell(sample_point + vec3(cell_offset.x, 0, 0), particle_indices, particle_count);
    //     getParticlesFromHashCell(sample_point + vec3(0, cell_offset.y, 0), particle_indices, particle_count);
    //     getParticlesFromHashCell(sample_point + vec3(0, 0, cell_offset.z), particle_indices, particle_count);
        
    //     // DEBUG: If any particles found, just return red immediately
    //     if (particle_count > 0) {
    //             for (int p = 0; p < particle_count; p++) {
    //             vec4 pos_mass = getParticlePosition(particle_indices[p]);
    //             vec4 properties = getParticleProperties(particle_indices[p]);
    //             float dist = length(sample_point - pos_mass.xyz);
    //             float smoothing_length = properties.z;
                
    //             if (dist < smoothing_length) {
    //                 return vec3(1.0, 0.0, 0.0); // Red if within influence
    //             }
    //         }
    //         return vec3(0.0, 1.0, 0.0); // Green if particles found but outside influence
    //     }
        
    //     current_pos += ray_dir * step_size;
    // }

    //  return vec3(0.0, 0.0, 0.0); // Black if no particles found
    vec3 accumulated_colour = vec3(0.0);
    float accumulated_opacity = 0.0;
    
    int samples = 24;
    float step_size = max_distance / float(samples);
    vec4 current_pos = start_pos;
    
    for (int i = 0; i < samples; ++i) {
        if (accumulated_opacity > 0.9999) break; // Early termination
        
        vec3 sample_point = current_pos.yzw;
        vec3 step_contribution = vec3(0.0);
        float step_opacity = 0.0;
        
        // Get particles from current hash cell and neighbors
        int particle_indices[9]; // Max particles to check chnaging this crashes the programs
        int particle_count = 0;
        
        // Check current cell
        getParticlesFromHashCell(sample_point, particle_indices, particle_count);

        // Check all 26 neighboring cells (3x3x3 - 1)
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) continue; // Skip center (already checked)
                    if (particle_count >= 9) break; // Array limit
                    
                    vec3 neighbor_point = sample_point + vec3(dx, dy, dz) * u_grid_cell_size;
                    getParticlesFromHashCell(neighbor_point, particle_indices, particle_count);
                }
            }
        }
        
        // Process only the nearby particles
        for (int p = 0; p < particle_count; p++) {
            int particle_idx = particle_indices[p];
            
            vec4 pos_mass = getParticlePosition(particle_idx);
            vec3 particle_pos = pos_mass.xyz;
            
            float dist = length(sample_point - particle_pos);
            
            // Check if particle is within influence radius
            vec4 properties = getParticleProperties(particle_idx);
            float smoothing_length = properties.z;

            if (dist < smoothing_length) {
                vec3 particle_contrib = calculateParticleContribution(sample_point, particle_idx, observer_pos);
                step_contribution += particle_contrib;
                
                // FIXED: Boost the opacity contribution significantly
                vec4 vel_density = getParticleVelocity(particle_idx);
                float density = vel_density.w;
                float kernel_weight = wendlandC2Kernel(dist, smoothing_length);
                
                // Much larger multiplier instead of 0.1
                float particle_opacity = density * kernel_weight * step_size * u_disk_opacity; // Increased from 0.1 to 5.0
                step_opacity += particle_opacity;
            }
        }
        
        // Limit opacity per step
        step_opacity = min(step_opacity, 0.9);
        
        // Apply extinction
        float extinction = exp(-accumulated_opacity * 1.5);
        
        accumulated_colour += step_contribution * extinction;
        accumulated_opacity += step_opacity * (1.0 - accumulated_opacity * 0.25);
        
        current_pos += ray_dir * step_size;
    }
    
    return accumulated_colour;
}

// Main ray tracing function adapted for particles
vec3 trace_kerr_ray_particles(vec3 dir, vec4 camPos, mat4 axes) {
    // Debug: Return fixed colour to test if function is being called
    //return vec3(0.5, 0.2, 0.8); // Uncomment this line for basic test
    
    vec4 pos = camPos;
    vec4 dir4D = -axes[0] + vec4(0.0, dir.x, dir.y, dir.z);
   
    vec4 p = metric(pos) * dir4D;

    vec3 disk_contribution = vec3(0.0);
    vec4 final_pos;

    bool in_disk_region = false;
    vec4 disk_entry_pos;

    for (int i = 0; i < u_max_steps; i++) {
        transportStep(pos, p);
        //remember that in Kerr-Schild coordinates, the time coordinate is the first component so it's x,y,z,w 
        float r = rFromCoords(pos);
        bool currently_in_disk = (r >= u_disk_inner_radius &&
                                  r <= u_disk_outer_radius &&
                                  abs(pos.w) < u_disk_thickness);

        if (currently_in_disk && !in_disk_region) {
            in_disk_region = true;
            disk_entry_pos = pos;
        } else if (!currently_in_disk && in_disk_region) {
            float disk_distance = length(pos.yzw - disk_entry_pos.yzw);
            if (disk_distance > 0.001 && disk_distance < 1000.0) {
                // Use particle-based volumetric rendering
                disk_contribution += particleVolumetricRender(
                    disk_entry_pos,
                    normalize(pos - disk_entry_pos),
                    disk_distance,
                    camPos,
                    axes  // Observer position THIS IS PROB FUCKED AGAIN CAM IS 90 OFF
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
            disk_contribution += particleVolumetricRender(
                disk_entry_pos,
                normalize(final_pos - disk_entry_pos),
                disk_distance,
                camPos,
                axes
            );
        }
    }

    // Background determination
    float final_r = rFromCoords(final_pos);
    bool captured = final_r < (M + sqrt(M * M - a * a));
    vec3 background_colour = vec3(0.0);

    if (!captured) {
        vec4 out_dir = inverse(metric(final_pos)) * p;
        vec3 cube_dir = normalize(vec3(out_dir.y, out_dir.z, out_dir.w));
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
    vec4 camPos = vec4(0.0, u_cam_pos.x, u_cam_pos.y, u_cam_pos.z);
    
    // Map 3D camera vectors to 4D space
    vec4 aim = vec4(0.0, u_cam_forward.x, u_cam_forward.y, u_cam_forward.z);
    vec4 vert = vec4(0.0, u_cam_up.x, u_cam_up.y, u_cam_up.z);
    vec4 timeDir = vec4(1.0, 0.0, 0.0, 0.0);
    
    mat4 camFrame = tetrad(camPos, timeDir, aim, vert);

    vec3 colour = trace_kerr_ray_particles(ray_dir, camPos, camFrame);

    colour *= 1.0 / 5.0;
    vec3 mapped = toneMap(colour);
    vec3 final = toSRGB(mapped);
    FragColour = vec4(final, 1.0);
}