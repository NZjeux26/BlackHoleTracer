#version 410

layout(local_size_x = 64) in;

// SPH particle structure
struct SPHParticle {
    vec4 position;       // (t, x, y, z) in Kerr-Schild coordinates
    vec4 velocity;       // 4-velocity
    float mass;
    float density;
    float pressure;
    float temperature;
    float smoothing_h;   // Smoothing length
    uint hash_key;       // For spatial hashing
    uint neighbors[32];  // Neighbor indices (limited for GPU memory)
    uint neighbor_count;
    float padding[3];    // Alignment padding
};

// Compute shader storage buffers
layout(std430, binding = 0) restrict buffer ParticleBuffer {
    SPHParticle particles[];
};

layout(std430, binding = 1) restrict buffer HashBuffer {
    uint hash_table[];
};

layout(std430, binding = 2) restrict buffer CellBuffer {
    uint cell_start[];
    uint cell_count[];
};

// Uniforms
uniform float u_dt;
uniform float u_mass;           // Black hole mass
uniform float u_spin;           // Black hole spin
uniform int u_particle_count;
uniform float u_smoothing_h;
uniform float u_gas_constant;
uniform float u_viscosity;
uniform vec3 u_hash_grid_size;
uniform float u_cell_size;

// SPH kernel functions
float poly6_kernel(float r, float h) {
    if (r >= h) return 0.0;
    float h2 = h * h;
    float h9 = h2 * h2 * h2 * h2 * h;
    float q = (h2 - r * r);
    return (315.0 / (64.0 * 3.14159265359 * h9)) * q * q * q;
}

vec3 spiky_gradient(vec3 r, float r_mag, float h) {
    if (r_mag >= h || r_mag < 1e-6) return vec3(0.0);
    float h6 = h * h * h * h * h * h;
    float coeff = -45.0 / (3.14159265359 * h6);
    float term = (h - r_mag) * (h - r_mag);
    return coeff * term * (r / r_mag);
}

float viscosity_laplacian(float r, float h) {
    if (r >= h) return 0.0;
    float h6 = h * h * h * h * h * h;
    return (45.0 / (3.14159265359 * h6)) * (h - r);
}

// Kerr metric computation (same as fragment shader)
float rFromCoords(vec4 pos) {
    float a = u_spin * u_mass;
    vec3 p = pos.yzw;
    float rho2 = dot(p, p) - a * a;
    float r2 = 0.5 * (rho2 + sqrt(rho2 * rho2 + 4.0 * a * a * p.z * p.z));
    return sqrt(r2);
}

mat4 metric(vec4 pos) {
    float M = u_mass;
    float a = u_spin * M;
    float r = rFromCoords(pos);
    vec4 k = vec4(-1.0, (r * pos.y - a * pos.z) / (r * r + a * a),
                        (r * pos.z + a * pos.y) / (r * r + a * a),
                         pos.w / r);
    float f = 2.0 * M * r / (r * r + a * a * pos.w * pos.w / (r * r));
    mat4 diag_metric = mat4(
        -1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1
    );
    return f * mat4(k.x * k, k.y * k, k.z * k, k.w * k) + diag_metric;
}

// Spatial hashing
uint hash_position(vec3 pos) {
    ivec3 grid_pos = ivec3(floor(pos / u_cell_size));
    const uint p1 = 73856093u;
    const uint p2 = 19349663u;
    const uint p3 = 83492791u;
    return (uint(grid_pos.x) * p1) ^ (uint(grid_pos.y) * p2) ^ (uint(grid_pos.z) * p3);
}

void find_neighbors(uint particle_idx) {
    SPHParticle p = particles[particle_idx];
    vec3 pos = p.position.yzw;
    
    particles[particle_idx].neighbor_count = 0;
    uint count = 0;
    
    // Search in 3x3x3 grid around particle
    for (int dx = -1; dx <= 1 && count < 32; dx++) {
        for (int dy = -1; dy <= 1 && count < 32; dy++) {
            for (int dz = -1; dz <= 1 && count < 32; dz++) {
                vec3 search_pos = pos + vec3(dx, dy, dz) * u_cell_size;
                uint hash = hash_position(search_pos) % hash_table.length();
                
                uint cell_idx = hash_table[hash];
                if (cell_idx == 0xFFFFFFFF) continue;
                
                for (uint i = 0; i < cell_count[cell_idx] && count < 32; i++) {
                    uint neighbor_idx = cell_start[cell_idx] + i;
                    if (neighbor_idx >= u_particle_count || neighbor_idx == particle_idx) continue;
                    
                    vec3 neighbor_pos = particles[neighbor_idx].position.yzw;
                    float dist = length(pos - neighbor_pos);
                    
                    if (dist < p.smoothing_h) {
                        particles[particle_idx].neighbors[count] = neighbor_idx;
                        count++;
                    }
                }
            }
        }
    }
    
    particles[particle_idx].neighbor_count = count;
}

// SPH density calculation
void compute_density(uint particle_idx) {
    SPHParticle p = particles[particle_idx];
    float density = 0.0;
    
    // Self contribution
    density += p.mass * poly6_kernel(0.0, p.smoothing_h);
    
    // Neighbor contributions
    for (uint i = 0; i < p.neighbor_count; i++) {
        uint neighbor_idx = p.neighbors[i];
        SPHParticle neighbor = particles[neighbor_idx];
        
        vec3 r = p.position.yzw - neighbor.position.yzw;
        float r_mag = length(r);
        
        density += neighbor.mass * poly6_kernel(r_mag, p.smoothing_h);
    }
    
    particles[particle_idx].density = max(density, 1e-6);
    
    // Calculate pressure using ideal gas law
    // P = ρ * k * T, where k is gas constant and T is temperature
    particles[particle_idx].pressure = particles[particle_idx].density * u_gas_constant * p.temperature;
}

// SPH force calculation
vec4 compute_sph_forces(uint particle_idx) {
    SPHParticle p = particles[particle_idx];
    vec4 force = vec4(0.0);
    
    for (uint i = 0; i < p.neighbor_count; i++) {
        uint neighbor_idx = p.neighbors[i];
        SPHParticle neighbor = particles[neighbor_idx];
        
        vec3 r = p.position.yzw - neighbor.position.yzw;
        float r_mag = length(r);
        
        if (r_mag < 1e-6) continue;
        
        // Pressure force
        vec3 pressure_gradient = spiky_gradient(r, r_mag, p.smoothing_h);
        float pressure_term = (p.pressure / (p.density * p.density)) + 
                             (neighbor.pressure / (neighbor.density * neighbor.density));
        vec3 pressure_force = -neighbor.mass * pressure_term * pressure_gradient;
        
        // Viscosity force
        vec3 vel_diff = p.velocity.yzw - neighbor.velocity.yzw;
        float viscosity_term = viscosity_laplacian(r_mag, p.smoothing_h);
        vec3 viscosity_force = u_viscosity * neighbor.mass * vel_diff * viscosity_term / neighbor.density;
        
        force.yzw += pressure_force + viscosity_force;
    }
    
    return force;
}

// Kerr gravitational force (geodesic acceleration)
vec4 compute_kerr_gravity(uint particle_idx) {
    SPHParticle p = particles[particle_idx];
    vec4 pos = p.position;
    vec4 vel = p.velocity;
    
    // Compute Christoffel symbols numerically for geodesic equation
    // This is a simplified version - full implementation would need all connection coefficients
    mat4 g = metric(pos);
    float eps = 1e-6;
    
    vec4 gravity = vec4(0.0);
    
    // Simplified radial gravitational acceleration
    float r = rFromCoords(pos);
    float M = u_mass;
    float a = u_spin * M;
    
    // Basic Newtonian-like acceleration modified by metric
    vec3 radial_dir = normalize(pos.yzw);
    float r2 = r * r;
    float a2 = a * a;
    
    // Effective potential approach
    float L = length(cross(pos.yzw, vel.yzw)); // Angular momentum approximation
    float effective_mass_factor = 1.0 - 2.0 * M / r + a2 / r2;
    
    vec3 radial_accel = -M / (r2 * effective_mass_factor) * radial_dir;
    
    // Frame dragging effect
    vec3 frame_drag = 2.0 * M * a / (r * r2) * vec3(-pos.z, 0.0, pos.y);
    
    gravity.yzw = radial_accel + frame_drag;
    
    return gravity;
}

// Main compute shader stages
void main() {
    uint particle_idx = gl_GlobalInvocationID.x;
    if (particle_idx >= u_particle_count) return;
    
    // Multi-pass approach using gl_LocalInvocationIndex to determine stage
    uint stage = gl_LocalInvocationIndex % 4;
    
    barrier();
    
    switch (stage) {
        case 0: // Spatial hashing
            {
                vec3 pos = particles[particle_idx].position.yzw;
                uint hash = hash_position(pos) % hash_table.length();
                particles[particle_idx].hash_key = hash;
            }
            break;
            
        case 1: // Find neighbors
            find_neighbors(particle_idx);
            break;
            
        case 2: // Compute density and pressure
            compute_density(particle_idx);
            break;
            
        case 3: // Integrate forces
            {
                vec4 sph_force = compute_sph_forces(particle_idx);
                vec4 gravity_force = compute_kerr_gravity(particle_idx);
                vec4 total_force = sph_force + gravity_force;
                
                // Leapfrog integration
                particles[particle_idx].velocity += total_force * u_dt / particles[particle_idx].mass;
                particles[particle_idx].position += particles[particle_idx].velocity * u_dt;
                
                // Update temperature based on compression/expansion
                float compression_rate = -particles[particle_idx].density / max(particles[particle_idx].density, 1e-6);
                particles[particle_idx].temperature += compression_rate * u_dt * 1000.0;
                particles[particle_idx].temperature = clamp(particles[particle_idx].temperature, 1000.0, 100000.0);
            }
            break;
    }
}