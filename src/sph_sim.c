#include "sph_sim.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

//========================================
// System Management Functions
//========================================

SPHSystem* sph_create_system(int max_particles, BlackHoleParams* black_hole) {
    // Validate input parameters
    if (max_particles <= 0 || max_particles > MAX_PARTICLES) {
        fprintf(stderr, "Error: Invalid max_particles (%d). Must be between 1 and %d\n", 
                max_particles, MAX_PARTICLES);
        return NULL;
    }

    // Allocate system structure
    SPHSystem* system = (SPHSystem*)malloc(sizeof(SPHSystem));
    if (!system) {
        fprintf(stderr, "Error: Failed to allocate memory for SPHSystem\n");
        return NULL;
    }

    // initialise basic parameters
    system->max_particles = max_particles;
    system->particle_count = 0;
    system->hash_table_size = HASH_TABLE_SIZE;
    system->grid_cell_size = GRID_CELL_SIZE;

    // Allocate particle array
    system->particles = (SPHParticle*)calloc(max_particles, sizeof(SPHParticle));
    if (!system->particles) {
        fprintf(stderr, "Error: Failed to allocate memory for particles\n");
        free(system);
        return NULL;
    }

    // Allocate spatial hash table
    system->hash_table = (HashCell*)calloc(system->hash_table_size, sizeof(HashCell));
    if (!system->hash_table) {
        fprintf(stderr, "Error: Failed to allocate memory for hash table\n");
        free(system->particles);
        free(system);
        return NULL;
    }

    // initialise SPH parameters with defaults from header
    system->kernel_radius = SPH_KERNEL_RADIUS;
    system->kernel_radius_sq = SPH_KERNEL_RADIUS_SQ;
    system->rest_density = REST_DENSITY;
    system->gas_constant = GAS_CONSTANT;
    system->viscosity_coeff = VISCOSITY_COEFF;
    system->surface_tension = SURFACE_TENSION;
    system->thermal_conductivity = THERMAL_CONDUCTIVITY;

    // initialise time stepping parameters
    system->dt = 0.001;  // 1ms default time step
    system->current_time = 0.0;
    system->max_iterations = 1000;
    system->adaptive_time_step = true;

    // Set black hole reference
    system->black_hole = black_hole;

    // initialise world bounds (will be updated as particles are added)
    system->world_min = (Vec3){-10.0, -10.0, -10.0};
    system->world_max = (Vec3){10.0, 10.0, 10.0};

    // initialise OpenGL render buffer (will be created later if needed)
    system->render_vbo = 0;

    // initialise performance counters
    system->last_update_time = 0.0;
    system->neighbor_searches = 0;
    system->density_calculations = 0;

    printf("SPH System created successfully:\n");
    printf("  Max particles: %d\n", max_particles);
    printf("  Hash table size: %d\n", system->hash_table_size);
    printf("  Kernel radius: %.3f\n", system->kernel_radius);
    printf("  Rest density: %.3f\n", system->rest_density);

    return system;
}

void sph_destroy_system(SPHSystem* system) {
    if (!system) {
        return;
    }

    printf("Destroying SPH system...\n");

    // Clean up OpenGL resources if they exist
    if (system->render_vbo != 0) {
        glDeleteBuffers(1, &system->render_vbo);
        system->render_vbo = 0;
    }

    // Free allocated memory
    if (system->particles) {
        free(system->particles);
        system->particles = NULL;
    }

    if (system->hash_table) {
        free(system->hash_table);
        system->hash_table = NULL;
    }

    // Note: We don't free black_hole as it's owned by the caller
    system->black_hole = NULL;

    // Free the system structure itself
    free(system);

    printf("SPH system destroyed successfully\n");
}

void sph_reset_system(SPHSystem* system) {
    if (!system) {
        fprintf(stderr, "Error: Cannot reset NULL system\n");
        return;
    }

    printf("Resetting SPH system...\n");

    // Reset particle count and clear all particle data
    system->particle_count = 0;
    memset(system->particles, 0, system->max_particles * sizeof(SPHParticle));

    // Clear spatial hash table
    sph_clear_hash_table(system);

    // Reset time stepping
    system->current_time = 0.0;
    system->dt = 0.001;  // Reset to default time step

    // Reset world bounds
    system->world_min = (Vec3){-10.0, -10.0, -10.0};
    system->world_max = (Vec3){10.0, 10.0, 10.0};

    // Reset performance counters
    system->last_update_time = 0.0;
    system->neighbor_searches = 0;
    system->density_calculations = 0;

    printf("SPH system reset complete\n");
}

void sph_update_system(SPHSystem* system, double dt) {
    if (!system) {
        fprintf(stderr, "Error: Cannot update NULL system\n");
        return;
    }

    if (system->particle_count == 0) {
        // Nothing to update
        return;
    }

    if (dt <= 0.0) {
        fprintf(stderr, "Warning: Invalid time step %.6f, skipping update\n", dt);
        return;
    }

    // Store the time step
    double original_dt = dt; //this isn't used
    
    // Use adaptive time stepping if enabled
    if (system->adaptive_time_step) {
        double adaptive_dt = sph_calculate_adaptive_time_step(system);
        if (adaptive_dt > 0.0 && adaptive_dt < dt) {
            dt = adaptive_dt;
        }
    }

    // Store time step in system
    system->dt = dt;

    // Clear and rebuild spatial hash table for neighbor finding
    sph_clear_hash_table(system);
    sph_build_hash_table(system);

    // Find neighbors for all particles
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            sph_find_neighbors(system, i);
        }
    }
    system->neighbor_searches += system->particle_count;

    // Update smoothing lengths based on local particle density
    sph_update_smoothing_lengths(system);

    // Calculate fluid properties
    sph_calculate_density(system);
    system->density_calculations += system->particle_count;
    
    sph_calculate_pressure(system);
    sph_calculate_temperature(system);

    // Calculate forces
    sph_calculate_forces(system);
    sph_calculate_viscosity_forces(system);
    sph_calculate_surface_tension_forces(system);

    // Apply thermal effects
    sph_calculate_thermal_diffusion(system);
    sph_apply_radiative_cooling(system, dt);
    sph_apply_viscous_heating(system, dt);

    // Integrate equations of motion
    // Use RK4 for better accuracy with black hole dynamics
    sph_integrate_rk4(system, dt);

    // Apply boundary conditions
    sph_apply_boundary_conditions(system);

    // Check for particles crossing event horizon
    if (system->black_hole) {
        for (int i = system->particle_count - 1; i >= 0; i--) {
            if (system->particles[i].flags & PARTICLE_ACTIVE) {
                if (sph_check_event_horizon_crossing(system->particles[i].position, 
                                                   system->black_hole)) {
                    // Mark particle as accreted
                    system->particles[i].flags |= PARTICLE_ACCRETING;
                    system->particles[i].flags &= ~PARTICLE_ACTIVE;
                }
            }
        }
    }

    // Update simulation time
    system->current_time += dt;
    system->last_update_time = dt;

    // Validate system state in debug builds
    #ifdef DEBUG
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            if (!sph_validate_particle_state(system, i)) {
                fprintf(stderr, "Warning: Particle %d failed validation at time %.6f\n", 
                        i, system->current_time);
            }
        }
    }
    
    // Check conservation laws periodically
    static int conservation_check_counter = 0;
    if (++conservation_check_counter >= 100) {
        sph_check_conservation_laws(system);
        sph_detect_numerical_instabilities(system);
        conservation_check_counter = 0;
    }
    #endif

    // Update render buffer if VBO exists
    if (system->render_vbo != 0) {
        // This would typically update the GPU buffer
        // Implementation depends on rendering system
    }
}

//========================================
// Particle Management Functions
//========================================

/**
 * Adds a new particle to the SPH system
 * @param system Pointer to the SPH system
 * @param position Initial position of the particle
 * @param velocity Initial velocity of the particle
 * @param mass Mass of the particle
 * @return Index of the newly added particle, or -1 if failed
 */
int sph_add_particle(SPHSystem* system, Vec3 position, Vec3 velocity, double mass) {
    // Validate input parameters
    if (!system || !system->particles) {
        fprintf(stderr, "Error: Invalid SPH system pointer\n");
        return -1;
    }
    
    if (system->particle_count >= system->max_particles) {
        fprintf(stderr, "Error: Maximum particle count reached (%d)\n", system->max_particles);
        return -1;
    }
    
    if (mass <= 0.0) {
        fprintf(stderr, "Warning: Particle mass should be positive, using default mass\n");
        mass = 1.0;
    }
    
    // Get the index for the new particle
    int index = system->particle_count;
    SPHParticle* particle = &system->particles[index];
    
    // initialise kinematic properties
    particle->position = position;
    particle->velocity = velocity;
    particle->acceleration = (Vec3){0.0, 0.0, 0.0};
    particle->prev_position = position;
    particle->prev_velocity = velocity;
    
    // initialise physical properties
    particle->mass = mass;
    particle->density = system->rest_density;  // Start with rest density
    particle->pressure = 0.0;
    particle->temperature = 1000.0;  // Default temperature in Kelvin
    particle->specific_heat = 1000.0;  // J/(kg·K) - typical for gas
    particle->viscosity = system->viscosity_coeff;
    
    // initialise SPH-specific properties
    particle->smoothing_length = system->kernel_radius;
    particle->color_gradient = (Vec3){0.0, 0.0, 0.0};
    particle->color_laplacian = 0.0;
    particle->surface_normal = (Vec3){0.0, 0.0, 0.0};
    particle->curvature = 0.0;
    
    // initialise thermal properties
    particle->thermal_energy = particle->specific_heat * particle->temperature;
    particle->thermal_diffusion = 0.0;
    particle->radiative_cooling = 0.0;
    particle->heating_rate = 0.0;
    
    // initialise simulation properties
    particle->flags = PARTICLE_ACTIVE;
    particle->id = (uint32_t)index;  // Simple ID assignment
    particle->neighbor_count = 0;
    
    // Clear neighbor list
    memset(particle->neighbors, 0, MAX_NEIGHBORS * sizeof(uint32_t));
    
    // initialise Kerr black hole coordinates (Boyer-Lindquist)
    if (system->black_hole) {
        double r = vec3_length(position);
        particle->r_coordinate = r;
        particle->theta_coordinate = acos(position.z / r);
        particle->phi_coordinate = atan2(position.y, position.x);
        
        // Calculate orbital properties
        particle->orbital_velocity = sph_calculate_orbital_velocity_kerr(
            particle->r_coordinate, particle->theta_coordinate, system->black_hole);
        particle->angular_momentum = particle->mass * particle->orbital_velocity * 
                                   particle->r_coordinate * sin(particle->theta_coordinate);
        particle->binding_energy = sph_calculate_binding_energy(position, velocity, system->black_hole);
    } else {
        particle->orbital_velocity = 0.0;
        particle->angular_momentum = 0.0;
        particle->binding_energy = 0.0;
        particle->r_coordinate = vec3_length(position);
        particle->theta_coordinate = 0.0;
        particle->phi_coordinate = 0.0;
    }
    
    // Increment particle count
    system->particle_count++;
    
    return index;
}

/**
 * Removes a particle from the SPH system
 * @param system Pointer to the SPH system
 * @param index Index of the particle to remove
 */
void sph_remove_particle(SPHSystem* system, int index) {
    // Validate input parameters
    if (!system || !system->particles) {
        fprintf(stderr, "Error: Invalid SPH system pointer\n");
        return;
    }
    
    if (index < 0 || index >= system->particle_count) {
        fprintf(stderr, "Error: Invalid particle index %d (valid range: 0-%d)\n", 
                index, system->particle_count - 1);
        return;
    }
    
    // Mark particle as inactive instead of actually removing it
    // This prevents issues with neighbor lists and indexing
    system->particles[index].flags &= ~PARTICLE_ACTIVE;
    system->particles[index].mass = 0.0;
    system->particles[index].density = 0.0;
    
    // If removing the last particle, we can actually decrease the count
    if (index == system->particle_count - 1) {
        system->particle_count--;
        
        // Check if there are more inactive particles at the end
        while (system->particle_count > 0 && 
               !(system->particles[system->particle_count - 1].flags & PARTICLE_ACTIVE)) {
            system->particle_count--;
        }
    }
    
    // Alternative approach: Swap with last particle and decrease count
    // This is more memory efficient but changes particle indices
    /*
    if (index != system->particle_count - 1) {
        // Copy last particle to the removed particle's position
        system->particles[index] = system->particles[system->particle_count - 1];
        // Update the moved particle's ID to reflect new index
        system->particles[index].id = (uint32_t)index;
    }
    system->particle_count--;
    */
}

/**
 * Sets physical properties of a specific particle
 * @param system Pointer to the SPH system
 * @param index Index of the particle to modify
 * @param density Density value to set
 * @param temperature Temperature value to set
 * @param flags Particle flags to set
 */
void sph_set_particle_properties(SPHSystem* system, int index, double density,
                                 double temperature, uint32_t flags) {
    // Validate input parameters
    if (!system || !system->particles) {
        fprintf(stderr, "Error: Invalid SPH system pointer\n");
        return;
    }
    
    if (index < 0 || index >= system->particle_count) {
        fprintf(stderr, "Error: Invalid particle index %d (valid range: 0-%d)\n", 
                index, system->particle_count - 1);
        return;
    }
    
    SPHParticle* particle = &system->particles[index];
    
    // Validate and set density
    if (density > 0.0) {
        if (density < MIN_DENSITY) {
            fprintf(stderr, "Warning: Density %f below minimum threshold %f\n", 
                    density, MIN_DENSITY);
            particle->density = MIN_DENSITY;
        } else if (density > MAX_DENSITY) {
            fprintf(stderr, "Warning: Density %f above maximum threshold %f\n", 
                    density, MAX_DENSITY);
            particle->density = MAX_DENSITY;
        } else {
            particle->density = density;
        }
        
        // Recalculate pressure based on new density
        particle->pressure = system->gas_constant * (particle->density - system->rest_density);
    } else {
        fprintf(stderr, "Warning: Invalid density %f, keeping current value\n", density);
    }
    
    // Validate and set temperature
    if (temperature > 0.0) {
        particle->temperature = temperature;
        // Update thermal energy based on new temperature
        particle->thermal_energy = particle->specific_heat * particle->temperature;
    } else {
        fprintf(stderr, "Warning: Invalid temperature %f, keeping current value\n", temperature);
    }
    
    // Set particle flags
    particle->flags = flags;
    
    // Update thermal properties based on flags
    if (flags & PARTICLE_HOT) {
        particle->heating_rate = 1000.0;  // Increased heating for hot particles
    }
    
    if (flags & PARTICLE_BOUNDARY) {
        // Boundary particles have fixed properties
        particle->velocity = (Vec3){0.0, 0.0, 0.0};
        particle->acceleration = (Vec3){0.0, 0.0, 0.0};
    }
    
    if (flags & PARTICLE_ACCRETING) {
        // Accreting particles may have enhanced thermal properties
        particle->radiative_cooling = 10.0;  // Enhanced cooling
    }
    
    // Update smoothing length based on density (adaptive smoothing)
    if (particle->density > 0.0) {
        // Smoothing length inversely proportional to density^(1/3) in 3D
        particle->smoothing_length = system->kernel_radius * 
                                   pow(system->rest_density / particle->density, 1.0/3.0);
        
        // Clamp smoothing length to reasonable bounds
        double min_h = system->kernel_radius * 0.5;
        double max_h = system->kernel_radius * 2.0;
        
        if (particle->smoothing_length < min_h) {
            particle->smoothing_length = min_h;
        } else if (particle->smoothing_length > max_h) {
            particle->smoothing_length = max_h;
        }
    }
}

/**
 * Helper function to validate particle properties after modification
 * @param system Pointer to the SPH system
 * @param index Index of the particle to validate
 * @return true if particle properties are valid, false otherwise
 */
bool sph_validate_particle_properties(SPHSystem* system, int index) {
    if (!system || !system->particles || index < 0 || index >= system->particle_count) {
        return false;
    }
    
    SPHParticle* particle = &system->particles[index];
    
    // Check for NaN or infinite values
    if (!isfinite(particle->position.x) || !isfinite(particle->position.y) || 
        !isfinite(particle->position.z)) {
        fprintf(stderr, "Error: Particle %d has invalid position\n", index);
        return false;
    }
    
    if (!isfinite(particle->velocity.x) || !isfinite(particle->velocity.y) || 
        !isfinite(particle->velocity.z)) {
        fprintf(stderr, "Error: Particle %d has invalid velocity\n", index);
        return false;
    }
    
    if (!isfinite(particle->density) || particle->density <= 0.0) {
        fprintf(stderr, "Error: Particle %d has invalid density: %f\n", index, particle->density);
        return false;
    }
    
    if (!isfinite(particle->temperature) || particle->temperature <= 0.0) {
        fprintf(stderr, "Error: Particle %d has invalid temperature: %f\n", index, particle->temperature);
        return false;
    }
    
    if (!isfinite(particle->mass) || particle->mass <= 0.0) {
        fprintf(stderr, "Error: Particle %d has invalid mass: %f\n", index, particle->mass);
        return false;
    }
    
    return true;
}

//========================================
// Initialisation Functions
//========================================

// initialise accretion disk with particles distributed in a ring
void sph_initialise_accretion_disk(SPHSystem* system, double inner_radius, 
                                  double outer_radius, int num_particles) {
    if (!system || num_particles <= 0 || inner_radius >= outer_radius) {
        fprintf(stderr, "Invalid parameters for accretion disk initialization\n");
        return;
    }
    
    // Limit particles to system capacity
    if (num_particles > system->max_particles) {
        num_particles = system->max_particles;
        printf("Warning: Clamping particle count to %d\n", num_particles);
    }
    
    // Clear existing particles
    system->particle_count = 0;
    
    // Disk parameters
    double disk_height = (outer_radius - inner_radius) * 0.1; // Thin disk
    double mass_per_particle = 1.0; // Base mass per particle
    
    // Generate particles in cylindrical coordinates
    for (int i = 0; i < num_particles; i++) {
        // Radial distribution - more particles at inner edge (steeper density gradient)
        double u = (double)i / (double)(num_particles - 1);
        double r = inner_radius + (outer_radius - inner_radius) * sqrt(u);
        
        // Angular distribution - uniform
        double phi = 2.0 * M_PI * ((double)rand() / RAND_MAX);
        
        // Vertical distribution - Gaussian around midplane
        double z_scale = disk_height * 0.5;
        double z = z_scale * (2.0 * ((double)rand() / RAND_MAX) - 1.0) * 
                   exp(-0.5 * pow(2.0 * ((double)rand() / RAND_MAX) - 1.0, 2));
        
        // Convert to Cartesian coordinates
        Vec3 position = {
            .x = r * cos(phi),
            .y = r * sin(phi), 
            .z = z
        };
        
        // Initial velocity - will be set by Keplerian initialization
        Vec3 velocity = {0.0, 0.0, 0.0};
        
        // Add some radial velocity perturbation for turbulence
        double radial_perturbation = 0.01 * ((double)rand() / RAND_MAX - 0.5);
        velocity.x += radial_perturbation * cos(phi);
        velocity.y += radial_perturbation * sin(phi);
        
        // Add particle to system
        int particle_id = sph_add_particle(system, position, velocity, mass_per_particle);
        
        if (particle_id >= 0) {
            // Set additional properties
            SPHParticle* particle = &system->particles[particle_id];
            
            // Store cylindrical coordinates for reference
            particle->r_coordinate = r;
            particle->theta_coordinate = M_PI/2.0; // Disk plane
            particle->phi_coordinate = phi;
            
            // initialise thermal properties
            particle->temperature = 1000.0; // Initial temperature in K
            particle->specific_heat = 1.0;
            particle->thermal_energy = particle->specific_heat * particle->temperature;
            
            // Set particle flags
            particle->flags = PARTICLE_ACTIVE;
            if (r < inner_radius * 1.1) {
                particle->flags |= PARTICLE_HOT; // Inner particles are hotter
                particle->temperature *= 2.0;
            }
            
            // initialise smoothing length based on local density estimate
            particle->smoothing_length = SPH_KERNEL_RADIUS;
        }
    }
    
    printf("initialised accretion disk with %d particles (r: %.2f - %.2f)\n", 
           system->particle_count, inner_radius, outer_radius);
}

// Set Keplerian orbital velocities for all particles
void sph_initialise_keplerian_velocities(SPHSystem* system) {
    if (!system || !system->black_hole) {
        fprintf(stderr, "Invalid system or missing black hole for Keplerian initialization\n");
        return;
    }
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        // Skip inactive particles
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        // Calculate cylindrical radius
        double r_cyl = sqrt(particle->position.x * particle->position.x + 
                           particle->position.y * particle->position.y);
        
        if (r_cyl < 1e-10) continue; // Skip particles at origin
        
        // Calculate Keplerian orbital velocity for Kerr black hole
        double v_orbital = sph_calculate_orbital_velocity_kerr(r_cyl, M_PI/2.0, system->black_hole);
        
        // Unit vector in azimuthal direction (perpendicular to radial)
        Vec3 radial_unit = {
            .x = particle->position.x / r_cyl,
            .y = particle->position.y / r_cyl,
            .z = 0.0
        };
        
        Vec3 azimuthal_unit = {
            .x = -radial_unit.y,
            .y = radial_unit.x,
            .z = 0.0
        };
        
        // Set orbital velocity
        particle->velocity = vec3_scale(azimuthal_unit, v_orbital);
        
        // Add small random perturbations for turbulence
        double perturbation_magnitude = 0.05 * v_orbital;
        Vec3 perturbation = {
            .x = perturbation_magnitude * (2.0 * ((double)rand() / RAND_MAX) - 1.0),
            .y = perturbation_magnitude * (2.0 * ((double)rand() / RAND_MAX) - 1.0),
            .z = perturbation_magnitude * (2.0 * ((double)rand() / RAND_MAX) - 1.0) * 0.1
        };
        
        particle->velocity = vec3_add(particle->velocity, perturbation);
        
        // Store orbital properties
        particle->orbital_velocity = v_orbital;
        particle->angular_momentum = r_cyl * v_orbital * particle->mass;
        particle->binding_energy = sph_calculate_binding_energy(particle->position, 
                                                               particle->velocity, 
                                                               system->black_hole);
        
        // Copy current state to previous state for integration
        particle->prev_position = particle->position;
        particle->prev_velocity = particle->velocity;
    }
    
    printf("initialised Keplerian velocities for %d particles\n", system->particle_count);
}

// initialise thermal equilibrium based on accretion physics
void sph_initialise_thermal_equilibrium(SPHSystem* system) {
    if (!system || !system->black_hole) {
        fprintf(stderr, "Invalid system for thermal equilibrium initialization\n");
        return;
    }
    
    // Physical constants (in simulation units)
    const double stefan_boltzmann = 5.67e-8; // Stefan-Boltzmann constant
    const double opacity = 0.1; // Simplified opacity ***This should come from the BlackholeParams
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        // Skip inactive particles
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        // Calculate distance from black hole
        double r = vec3_length(particle->position);
        if (r < 1e-10) continue;
        
        // Temperature profile for thin accretion disk 
        // T ~ (GM*mdot/(8*pi*sigma*r^3))^(1/4)
        double r_cyl = sqrt(particle->position.x * particle->position.x + 
                           particle->position.y * particle->position.y);
        
        // Simplified temperature profile (decreases with radius) ***This or maybe the blackbody code I have
        double base_temp = 1000.0; // K
        double temp_scale = pow(r_cyl / 10.0, -0.75); // r^(-3/4) scaling
        particle->temperature = base_temp * temp_scale;
        
        // Clamp temperature to reasonable range
        if (particle->temperature < 100.0) particle->temperature = 100.0;
        if (particle->temperature > 10000.0) particle->temperature = 10000.0;
        
        // Calculate thermal energy
        particle->thermal_energy = particle->specific_heat * particle->temperature;
        
        // Estimate viscous heating rate based on shear
        double shear_rate = particle->orbital_velocity / r_cyl;
        particle->heating_rate = VISCOSITY_COEFF * particle->density * 
                                shear_rate * shear_rate;
        
        // Radiative cooling rate (simplified blackbody) ***This or maybe the blackbody code I have
        particle->radiative_cooling = 4.0 * stefan_boltzmann * 
                                     pow(particle->temperature, 4) * opacity;
        
        // Thermal diffusion coefficient
        particle->thermal_diffusion = THERMAL_CONDUCTIVITY / 
                                     (particle->density * particle->specific_heat);
        
        // Set thermal flags based on temperature
        if (particle->temperature > 5000.0) {
            particle->flags |= PARTICLE_HOT;
        }
        
        // Inner particles are actively accreting
        if (r_cyl < 5.0) {
            particle->flags |= PARTICLE_ACCRETING;
        }
    }
    
    // Calculate average thermal energy for diagnostics
    double total_thermal = 0.0;
    int active_count = 0;
    
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            total_thermal += system->particles[i].thermal_energy;
            active_count++;
        }
    }
    
    double avg_thermal = (active_count > 0) ? total_thermal / active_count : 0.0;
    
    printf("initialised thermal equilibrium for %d particles\n", active_count);
    printf("Average thermal energy: %.3f, Average temperature: %.1f K\n", 
           avg_thermal, avg_thermal / 1.0); // Assuming specific_heat = 1.0
}

//========================================
// Spatial Hashing Functions
//========================================

// Helper function to get grid coordinates
static inline void get_grid_coords(Vec3 position, double cell_size, int* x, int* y, int* z) {
    *x = (int)floor(position.x / cell_size);
    *y = (int)floor(position.y / cell_size);
    *z = (int)floor(position.z / cell_size);
}

// Improved hash function with better distribution
uint32_t sph_hash_position(Vec3 position, double cell_size) {
    int x, y, z;
    get_grid_coords(position, cell_size, &x, &y, &z);
    
    // Better hash function with prime coefficients
    uint32_t hash = 0;
    hash = hash * 73856093U + (uint32_t)x;
    hash = hash * 19349663U + (uint32_t)y;
    hash = hash * 83492791U + (uint32_t)z;
    
    return hash;
}

void sph_clear_hash_table(SPHSystem* system) {
    if (!system || !system->hash_table) return;
    
    // Clear all hash cells efficiently
    for (int i = 0; i < system->hash_table_size; i++) {
        system->hash_table[i].count = 0;
    }
}

void sph_build_hash_table(SPHSystem* system) {
    if (!system || !system->particles || !system->hash_table) return;
    
    // Clear table first
    sph_clear_hash_table(system);
    
    // Update world bounds for better spatial distribution
    if (system->particle_count > 0) {
        system->world_min = system->particles[0].position;
        system->world_max = system->particles[0].position;
        
        for (int i = 1; i < system->particle_count; i++) {
            if (system->particles[i].flags & PARTICLE_ACTIVE) {
                Vec3 pos = system->particles[i].position;
                if (pos.x < system->world_min.x) system->world_min.x = pos.x;
                if (pos.y < system->world_min.y) system->world_min.y = pos.y;
                if (pos.z < system->world_min.z) system->world_min.z = pos.z;
                if (pos.x > system->world_max.x) system->world_max.x = pos.x;
                if (pos.y > system->world_max.y) system->world_max.y = pos.y;
                if (pos.z > system->world_max.z) system->world_max.z = pos.z;
            }
        }
    }
    
    // Insert each active particle into hash table
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            uint32_t hash = sph_hash_position(system->particles[i].position, system->grid_cell_size);
            hash = hash % system->hash_table_size;
            
            HashCell* cell = &system->hash_table[hash];
            if (cell->count < 32) {  // Max particles per cell
                cell->particle_indices[cell->count] = i;
                cell->count++;
            } else {
                // Warning: cell overflow - might need larger hash table or smaller cells
                static int overflow_count = 0;
                if (overflow_count < 10) {  // Limit warning spam
                    printf("Warning: Hash cell overflow at hash %u\n", hash);
                    overflow_count++;
                }
            }
        }
    }
}

void sph_find_neighbors(SPHSystem* system, int particle_index) {
    if (!system || particle_index < 0 || particle_index >= system->particle_count) return;
    
    SPHParticle* particle = &system->particles[particle_index];
    particle->neighbor_count = 0;
    
    Vec3 pos = particle->position;
    double h = system->kernel_radius;
    double cell_size = system->grid_cell_size;
    
    // Calculate how many cells we need to search in each direction
    int search_radius = (int)ceil(h / cell_size);
    
    // Get the grid coordinates of the current particle
    int base_x, base_y, base_z;
    get_grid_coords(pos, cell_size, &base_x, &base_y, &base_z);
    
    // Search neighboring cells
    for (int dx = -search_radius; dx <= search_radius; dx++) {
        for (int dy = -search_radius; dy <= search_radius; dy++) {
            for (int dz = -search_radius; dz <= search_radius; dz++) {
                // Calculate the position of this cell
                Vec3 cell_pos = {
                    (base_x + dx) * cell_size,
                    (base_y + dy) * cell_size,
                    (base_z + dz) * cell_size
                };
                
                uint32_t hash = sph_hash_position(cell_pos, cell_size);
                hash = hash % system->hash_table_size;
                
                HashCell* cell = &system->hash_table[hash];
                
                // Check all particles in this cell
                for (int i = 0; i < cell->count; i++) {
                    int neighbor_idx = cell->particle_indices[i];
                    if (neighbor_idx == particle_index) continue;
                    
                    // Skip inactive particles
                    if (!(system->particles[neighbor_idx].flags & PARTICLE_ACTIVE)) continue;
                    
                    Vec3 neighbor_pos = system->particles[neighbor_idx].position;
                    Vec3 diff = vec3_sub(pos, neighbor_pos);
                    double dist_sq = vec3_dot(diff, diff);
                    
                    // Check if within kernel radius
                    if (dist_sq < h * h && particle->neighbor_count < MAX_NEIGHBORS) {
                        particle->neighbors[particle->neighbor_count] = neighbor_idx;
                        particle->neighbor_count++;
                    }
                }
                
                // Early exit if we've found enough neighbors
                if (particle->neighbor_count >= MAX_NEIGHBORS) {
                    return;
                }
            }
        }
    }
}

//========================================
// Kernel Functions
//========================================

// Wendland C2 kernel function
// This is a compact support kernel with C2 continuity, excellent for SPH
double sph_wendland_c2_kernel(double r, double h) {
    if (h <= 0.0) return 0.0;
    
    double q = r / h;
    
    // Wendland C2 kernel has compact support: W = 0 for q >= 2
    if (q >= 2.0) return 0.0;
    
    // Normalisation constant for 3D
    double alpha_3d = 21.0 / (2.0 * M_PI * h * h * h);
    
    // Wendland C2 kernel: W(q) = α * (1 - q/2)^4 * (2q + 1) for q < 2
    double factor = 1.0 - 0.5 * q;
    double factor4 = factor * factor * factor * factor;
    
    return alpha_3d * factor4 * (2.0 * q + 1.0);
}

// Gradient of Wendland C2 kernel (returns magnitude of gradient)
double sph_wendland_c2_gradient(double r, double h) {
    if (h <= 0.0 || r <= 0.0) return 0.0;
    
    double q = r / h;
    
    // Compact support: gradient = 0 for q >= 2
    if (q >= 2.0) return 0.0;
    
    // Normalisation constant for 3D
    double alpha_3d = 21.0 / (2.0 * M_PI * h * h * h);
    
    // Derivative of Wendland C2 kernel
    // dW/dr = (α/h) * d/dq[(1-q/2)^4 * (2q+1)]
    // = (α/h) * [(1-q/2)^3 * (-2) * (2q+1) + (1-q/2)^4 * 2]
    double factor = 1.0 - 0.5 * q;
    double factor3 = factor * factor * factor;
    double factor4 = factor3 * factor;
    
    double dw_dq = -2.0 * factor3 * (2.0 * q + 1.0) + 2.0 * factor4;
    
    return (alpha_3d / h) * dw_dq;
}

// Laplacian of Wendland C2 kernel
double sph_wendland_c2_laplacian(double r, double h) {
    if (h <= 0.0) return 0.0;
    
    double q = r / h;
    
    // Compact support
    if (q >= 2.0) return 0.0;
    
    // For very small r, use analytical limit
    if (r < 1e-10) {
        double alpha_3d = 21.0 / (2.0 * M_PI * h * h * h);
        return (alpha_3d / (h * h)) * (-10.0); // Limit as r->0
    }
    
    // Normalisation constant for 3D
    double alpha_3d = 21.0 / (2.0 * M_PI * h * h * h);
    
    // Laplacian in 3D: ∇²W = d²W/dr² + (2/r) * dW/dr
    double factor = 1.0 - 0.5 * q;
    double factor2 = factor * factor;
    double factor3 = factor2 * factor;
    
    // First derivative terms
    double dw_dq = -2.0 * factor3 * (2.0 * q + 1.0) + 2.0 * factor2 * factor;
    
    // Second derivative
    double d2w_dq2 = 6.0 * factor2 * (2.0 * q + 1.0) - 8.0 * factor3;
    
    double d2w_dr2 = (alpha_3d / (h * h)) * d2w_dq2;
    double dw_dr = (alpha_3d / h) * dw_dq;
    
    return d2w_dr2 + (2.0 / r) * dw_dr;
}

// Calculate kernel gradient vector between two particles
Vec3 sph_calculate_kernel_gradient(Vec3 ri, Vec3 rj, double h, int kernelType) {
    Vec3 r_ij = vec3_sub(ri, rj);
    double r = vec3_length(r_ij);
    
    // Handle zero distance case
    if (r < 1e-10) {
        Vec3 zero = {0.0, 0.0, 0.0};
        return zero;
    }
    
    double grad_magnitude = 0.0;
    
    switch (kernelType) {
        case 0: // Wendland C2 (default)
        default:
            grad_magnitude = sph_wendland_c2_gradient(r, h);
            break;
        
        case 1: // Cubic spline (alternative implementation)
        {
            double q = r / h;
            if (q >= 2.0) {
                grad_magnitude = 0.0;
            } else {
                double alpha_3d = 8.0 / (M_PI * h * h * h);
                if (q <= 1.0) {
                    grad_magnitude = (alpha_3d / h) * (-3.0 + 2.25 * q);
                } else {
                    double factor = 2.0 - q;
                    grad_magnitude = (alpha_3d / h) * (-0.75 * factor * factor / q);
                }
            }
            break;
        }
    }
    
    // Convert to vector: ∇W = (dW/dr) * (r_ij / |r_ij|)
    Vec3 unit_vector = vec3_scale(r_ij, 1.0 / r);
    return vec3_scale(unit_vector, grad_magnitude);
}

// Calculate kernel Laplacian between two particles
double sph_calculate_kernel_laplacian(Vec3 ri, Vec3 rj, double h, int kernelType) {
    Vec3 r_ij = vec3_sub(ri, rj);
    double r = vec3_length(r_ij);
    
    switch (kernelType) {
        case 0: // Wendland C2 (default)
        default:
            return sph_wendland_c2_laplacian(r, h);
        
        case 1: // Cubic spline Laplacian
        {
            if (r >= 2.0 * h) return 0.0;
            
            double q = r / h;
            double alpha_3d = 8.0 / (M_PI * h * h * h);
            
            if (r < 1e-10) {
                // Analytical limit as r->0
                return (alpha_3d / (h * h)) * (-6.0);
            }
            
            if (q <= 1.0) {
                return (alpha_3d / (h * h)) * (-6.0 + 9.0 * q);
            } else {
                double factor = 2.0 - q;
                return (alpha_3d / (h * h)) * (-3.0 * factor / (2.0 * q));
            }
        }
    }
}

//========================================
// Physical calculations
//========================================

// Calculate density for all particles using SPH kernel summation
void sph_calculate_density(SPHSystem* system) {
    if (!system) return;
    
    // Reset performance counter
    system->density_calculations = 0;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        // Skip inactive particles
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        double density_sum = 0.0;
        double h = particle_i->smoothing_length;
        
        // Self-contribution
        density_sum += particle_i->mass * sph_wendland_c2_kernel(0.0, h);
        
        // Neighbor contributions
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            // Distance between particles
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < h) {
                density_sum += particle_j->mass * sph_wendland_c2_kernel(r, h);
            }
        }
        
        // Set density with bounds checking
        particle_i->density = density_sum;
        if (particle_i->density < MIN_DENSITY) particle_i->density = MIN_DENSITY;
        if (particle_i->density > MAX_DENSITY) particle_i->density = MAX_DENSITY;
        
        system->density_calculations++;
    }
}

// Calculate pressure from density using equation of state
void sph_calculate_pressure(SPHSystem* system) {
    if (!system) return;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        // Skip inactive particles
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        // Tait equation of state: P = k * ((ρ/ρ₀)^γ - 1)
        // For weakly compressible fluid, use gas constant directly
        double density_ratio = particle->density / system->rest_density;
        
        if (density_ratio > 1.0) {
            // Standard gas law for compressible regions
            particle->pressure = system->gas_constant * (density_ratio - 1.0);
        } else {
            // Small negative pressure allowed for slight under-density
            particle->pressure = system->gas_constant * (density_ratio - 1.0) * 0.1;
        }
        
        // Add thermal pressure component
        if (particle->temperature > 0.0) {
            // P_thermal = ρ * R * T (simplified ideal gas)
            double thermal_pressure = particle->density * 0.1 * particle->temperature / 1000.0;
            particle->pressure += thermal_pressure;
        }
        
        // Ensure non-negative pressure for stability
        if (particle->pressure < 0.0) particle->pressure = 0.0;
    }
}

// Calculate pressure and viscosity forces on all particles
void sph_calculate_forces(SPHSystem* system) {
    if (!system) return;
    
    // Reset accelerations
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            system->particles[i].acceleration = (Vec3){0.0, 0.0, 0.0};
        }
    }
    
    // Calculate pairwise forces
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        double h_i = particle_i->smoothing_length;
        
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i >= j) continue; // Avoid double counting
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue; // Avoid singularity
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                Vec3 dr_unit = vec3_scale(dr, 1.0 / r);
                
                // Pressure force (symmetric)
                double pressure_term = (particle_i->pressure / (particle_i->density * particle_i->density)) +
                                     (particle_j->pressure / (particle_j->density * particle_j->density));
                
                double kernel_grad = sph_wendland_c2_gradient(r, h_avg);
                Vec3 pressure_force = vec3_scale(dr_unit, -particle_j->mass * pressure_term * kernel_grad);
                
                // Apply Newton's 3rd law
                particle_i->acceleration = vec3_add(particle_i->acceleration, 
                                                  vec3_scale(pressure_force, 1.0 / particle_i->mass));
                particle_j->acceleration = vec3_sub(particle_j->acceleration, 
                                                  vec3_scale(pressure_force, 1.0 / particle_j->mass));
            }
        }
        
        // Add black hole gravitational acceleration
        if (system->black_hole) {
            Vec3 bh_accel = sph_calculate_kerr_acceleration(particle_i->position, 
                                                           particle_i->velocity, 
                                                           system->black_hole);
            particle_i->acceleration = vec3_add(particle_i->acceleration, bh_accel);
        }
    }
}

// Calculate viscosity forces separately for clarity
void sph_calculate_viscosity_forces(SPHSystem* system) {
    if (!system || system->viscosity_coeff <= 0.0) return;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        Vec3 viscosity_accel = {0.0, 0.0, 0.0};
        double h_i = particle_i->smoothing_length;
        
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue;
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                Vec3 dv = vec3_sub(particle_j->velocity, particle_i->velocity);
                
                // Artificial viscosity (Monaghan form)
                double v_dot_r = vec3_dot(dv, dr);
                
                if (v_dot_r < 0.0) { // Only apply if particles approaching
                    double rho_avg = (particle_i->density + particle_j->density) * 0.5;
                    double c_s = sqrt(system->gas_constant); // Sound speed approximation
                    
                    double alpha = system->viscosity_coeff;
                    double beta = 2.0 * alpha;
                    
                    double mu = h_avg * v_dot_r / (r * r + 0.01 * h_avg * h_avg);
                    double pi_visc = (-alpha * c_s * mu + beta * mu * mu) / rho_avg;
                    
                    double kernel_grad = sph_wendland_c2_gradient(r, h_avg);
                    Vec3 viscous_force = vec3_scale(dr, -particle_j->mass * pi_visc * kernel_grad / r);
                    
                    viscosity_accel = vec3_add(viscosity_accel, 
                                             vec3_scale(viscous_force, 1.0 / particle_i->mass));
                }
            }
        }
        
        particle_i->acceleration = vec3_add(particle_i->acceleration, viscosity_accel);
    }
}

// Calculate surface tension forces
void sph_calculate_surface_tension_forces(SPHSystem* system) {
    if (!system || system->surface_tension <= 0.0) return;
    
    // First pass: calculate color field gradient and laplacian
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        particle_i->color_gradient = (Vec3){0.0, 0.0, 0.0};
        particle_i->color_laplacian = 0.0;
        
        double h_i = particle_i->smoothing_length;
        
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue;
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                double mass_over_density = particle_j->mass / particle_j->density;
                
                // Color field gradient
                Vec3 grad_kernel = sph_calculate_kernel_gradient(particle_i->position, 
                                                               particle_j->position, h_avg, 0);
                Vec3 contrib = vec3_scale(grad_kernel, mass_over_density);
                particle_i->color_gradient = vec3_add(particle_i->color_gradient, contrib);
                
                // Color field laplacian
                double lapl_kernel = sph_calculate_kernel_laplacian(particle_i->position, 
                                                                  particle_j->position, h_avg, 0);
                particle_i->color_laplacian += mass_over_density * lapl_kernel;
            }
        }
        
        // Calculate surface normal and curvature
        double grad_magnitude = vec3_length(particle_i->color_gradient);
        if (grad_magnitude > 1e-10) {
            particle_i->surface_normal = vec3_scale(particle_i->color_gradient, 1.0 / grad_magnitude);
            particle_i->curvature = -particle_i->color_laplacian / grad_magnitude;
            
            // Mark as surface particle if gradient is significant
            if (grad_magnitude > 0.1) {
                particle_i->flags |= PARTICLE_SURFACE;
            } else {
                particle_i->flags &= ~PARTICLE_SURFACE;
                particle_i->flags |= PARTICLE_INTERIOR;
            }
        } else {
            particle_i->surface_normal = (Vec3){0.0, 0.0, 0.0};
            particle_i->curvature = 0.0;
            particle_i->flags |= PARTICLE_INTERIOR;
        }
    }
    
    // Second pass: apply surface tension forces
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE) || 
            !(particle_i->flags & PARTICLE_SURFACE)) continue;
        
        // Surface tension force = σ * κ * n
        Vec3 surface_force = vec3_scale(particle_i->surface_normal, 
                                       system->surface_tension * particle_i->curvature);
        
        Vec3 surface_accel = vec3_scale(surface_force, 1.0 / particle_i->mass);
        particle_i->acceleration = vec3_add(particle_i->acceleration, surface_accel);
    }
}

// Calculate thermal diffusion
void sph_calculate_thermal_diffusion(SPHSystem* system) {
    if (!system || system->thermal_conductivity <= 0.0) return;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        double thermal_diffusion_rate = 0.0;
        double h_i = particle_i->smoothing_length;
        
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue;
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                // Temperature difference
                double dT = particle_j->temperature - particle_i->temperature;
                
                // Thermal diffusion coefficient (average)
                double k_avg = (particle_i->thermal_diffusion + particle_j->thermal_diffusion) * 0.5;
                
                // Laplacian of temperature approximation
                double kernel_lapl = sph_wendland_c2_laplacian(r, h_avg);
                
                thermal_diffusion_rate += k_avg * particle_j->mass * dT * kernel_lapl / particle_j->density;
            }
        }
        
        // Update thermal diffusion rate
        particle_i->thermal_diffusion = thermal_diffusion_rate;
    }
}
//========================================
// Black Hole Physics
//========================================

// Calculate gravitational acceleration from Kerr black hole
Vec3 sph_calculate_kerr_acceleration(Vec3 position, Vec3 velocity, BlackHoleParams* bh) {
    if (!bh) return (Vec3){0.0, 0.0, 0.0};
    
    double x = position.x, y = position.y, z = position.z;
    double vx = velocity.x, vy = velocity.y, vz = velocity.z;
    
    // Convert to spherical coordinates
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1e-10) return (Vec3){0.0, 0.0, 0.0};
    
    double theta = acos(z / r);
    double phi = atan2(y, x);
    
    // Kerr metric parameters
    double M = bh->mass;
    double a = bh->spin;  // Dimensionless spin parameter
    double a2 = a * a;
    
    // Boyer-Lindquist coordinates
    double r2 = r * r;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double cos2_theta = cos_theta * cos_theta;
    double sin2_theta = sin_theta * sin_theta;
    
    // Kerr metric functions
    double rho2 = r2 + a2 * cos2_theta;
    double delta = r2 - 2.0 * M * r + a2;
    double sigma = (r2 + a2) * (r2 + a2) - a2 * delta * sin2_theta;
    
    // Avoid singularities
    if (rho2 < 1e-10 || delta < 1e-10) {
        // Fall back to Schwarzschild near singularities
        double r3 = r2 * r;
        double accel_mag = M / r3;
        return vec3_scale(position, -accel_mag / r);
    }
    
    // Convert velocity to Boyer-Lindquist coordinates
    double vr = (x * vx + y * vy) / r;
    double vtheta = (z * (x * vx + y * vy) - r2 * vz) / (r * sqrt(r2 - z*z));
    double vphi = (x * vy - y * vx) / (x*x + y*y);
    
    if (sin_theta < 1e-10) vtheta = 0.0; // Handle pole
    if (x*x + y*y < 1e-10) vphi = 0.0;   // Handle z-axis
    
    // Kerr geodesic acceleration components
    double ar, atheta, aphi;
    
    // Radial acceleration
    ar = -M * (r2 - a2 * cos2_theta) / (rho2 * rho2 * rho2) * (r2 + a2) +
         r * (r2 + a2 * cos2_theta) / (rho2 * rho2) +
         2.0 * M * r * a2 * sin2_theta / (rho2 * rho2 * rho2) * (r2 + a2) +
         (r - M) * delta / (rho2 * rho2) * vr * vr +
         a2 * sin_theta * cos_theta / (rho2 * rho2) * vtheta * vtheta +
         2.0 * M * a * r * sin2_theta / (rho2 * rho2) * vphi * vphi;
    
    // Theta acceleration  
    atheta = a2 * sin_theta * cos_theta / (rho2 * rho2) * (1.0 - 2.0 * M * r / rho2) -
             2.0 * r / rho2 * vr * vtheta +
             2.0 * M * a * r * sin_theta * cos_theta / (rho2 * rho2) * vphi * vphi;
    
    // Phi acceleration
    aphi = -2.0 * M * a * r / (rho2 * rho2 * sigma) * (r2 + a2) +
           2.0 * a * (M * r - rho2) / (rho2 * sigma) * vr * vphi +
           2.0 * a * r * cos_theta / (rho2 * sin_theta) * vtheta * vphi;
    
    // Convert back to Cartesian coordinates
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    
    // Transformation matrix elements
    double ax = ar * sin_theta * cos_phi + atheta * cos_theta * cos_phi - aphi * sin_phi;
    double ay = ar * sin_theta * sin_phi + atheta * cos_theta * sin_phi + aphi * cos_phi;
    double az = ar * cos_theta - atheta * sin_theta;
    
    return (Vec3){ax, ay, az};
}

// Calculate Keplerian orbital velocity for Kerr metric
double sph_calculate_orbital_velocity_kerr(double r, double theta, BlackHoleParams* bh) {
    if (!bh || r <= 0.0) return 0.0;
    
    double M = bh->mass;
    double a = bh->spin;
    double a2 = a * a;
    
    // For circular orbits in Kerr metric (Boyer-Lindquist coordinates)
    // This is the ISCO-stable circular orbit velocity
    
    // Innermost stable circular orbit (ISCO) radius
    double Z1 = 1.0 + pow(1.0 - a2, 1.0/3.0) * (pow(1.0 + a, 1.0/3.0) + pow(1.0 - a, 1.0/3.0));
    double Z2 = sqrt(3.0 * a2 + Z1 * Z1);
    double r_isco = 3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
    
    // Don't allow orbits inside ISCO
    if (r < r_isco) r = r_isco;
    
    // Keplerian velocity with relativistic corrections
    double r2 = r * r;
    double r3 = r2 * r;
    
    // Leading term (Keplerian)
    double v_kepler = sqrt(M / r);
    
    // First-order relativistic correction
    double correction = 1.0 - 1.5 * M / r + a * sqrt(M) / (r * sqrt(r));
    
    // Ensure positive correction
    if (correction <= 0.0) correction = 0.1;
    
    return v_kepler * sqrt(correction);
}

// Calculate specific binding energy
double sph_calculate_binding_energy(Vec3 position, Vec3 velocity, BlackHoleParams* bh) {
    if (!bh) return 0.0;
    
    double r = vec3_length(position);
    if (r < 1e-10) return 0.0;
    
    double M = bh->mass;
    double a = bh->spin;
    
    // Convert to spherical coordinates
    double theta = acos(position.z / r);
    double phi = atan2(position.y, position.x);
    
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double cos2_theta = cos_theta * cos_theta;
    double sin2_theta = sin_theta * sin_theta;
    
    // Kerr metric functions
    double r2 = r * r;
    double a2 = a * a;
    double rho2 = r2 + a2 * cos2_theta;
    double delta = r2 - 2.0 * M * r + a2;
    double sigma = (r2 + a2) * (r2 + a2) - a2 * delta * sin2_theta;
    
    if (rho2 < 1e-10 || sigma < 1e-10) {
        // Fallback to Newtonian
        double v2 = vec3_dot(velocity, velocity);
        return 0.5 * v2 - M / r;
    }
    
    // Convert velocity to Boyer-Lindquist coordinates
    double vr = (position.x * velocity.x + position.y * velocity.y) / r;
    double vtheta = (position.z * (position.x * velocity.x + position.y * velocity.y) - 
                    r2 * velocity.z) / (r * sqrt(r2 - position.z * position.z));
    double vphi = (position.x * velocity.y - position.y * velocity.x) / (position.x * position.x + position.y * position.y);
    
    // Handle singularities
    if (sin_theta < 1e-10) vtheta = 0.0;
    if (position.x * position.x + position.y * position.y < 1e-10) vphi = 0.0;
    
    // Kerr metric components
    double g_tt = -(1.0 - 2.0 * M * r / rho2);
    double g_rr = rho2 / delta;
    double g_theta_theta = rho2;
    double g_phi_phi = sin2_theta * sigma / rho2;
    double g_t_phi = -2.0 * M * a * r * sin2_theta / rho2;
    
    // Specific energy (conserved quantity)
    double ut = 1.0; // Normalized time component
    double energy = -g_tt * ut * ut - 2.0 * g_t_phi * ut * vphi - g_phi_phi * vphi * vphi;
    energy += g_rr * vr * vr + g_theta_theta * vtheta * vtheta;
    
    return -sqrt(fabs(energy)); // Binding energy is negative
}

// Check if particle has crossed the event horizon
bool sph_check_event_horizon_crossing(Vec3 position, BlackHoleParams* bh) {
    if (!bh) return false;
    
    double r = vec3_length(position);
    double M = bh->mass;
    double a = bh->spin;
    
    // Event horizon radius for Kerr black hole
    double r_plus = M + sqrt(M * M - a * a);
    
    return (r <= r_plus);
}

//========================================
// RK Intergration
//========================================

// RK4 integration for the SPH system
void sph_integrate_rk4(SPHSystem* system, double dt) {
    if (!system || dt <= 0.0) return;
    
    // Temporary storage for RK4 stages
    Vec3* k1_pos = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k1_vel = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k2_pos = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k2_vel = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k3_pos = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k3_vel = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k4_pos = malloc(system->particle_count * sizeof(Vec3));
    Vec3* k4_vel = malloc(system->particle_count * sizeof(Vec3));
    
    Vec3* original_pos = malloc(system->particle_count * sizeof(Vec3));
    Vec3* original_vel = malloc(system->particle_count * sizeof(Vec3));
    
    if (!k1_pos || !k1_vel || !k2_pos || !k2_vel || !k3_pos || !k3_vel || 
        !k4_pos || !k4_vel || !original_pos || !original_vel) {
        fprintf(stderr, "Memory allocation failed in RK4 integration\n");
        goto cleanup;
    }
    
    // Store original state
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            original_pos[i] = system->particles[i].position;
            original_vel[i] = system->particles[i].velocity;
        }
    }
    
    // Stage 1: k1 = f(t, y)
    sph_build_hash_table(system);
    sph_calculate_density(system);
    sph_calculate_pressure(system);
    sph_calculate_forces(system);
    sph_calculate_viscosity_forces(system);
    
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            k1_pos[i] = system->particles[i].velocity;
            k1_vel[i] = system->particles[i].acceleration;
        }
    }
    
    // Stage 2: k2 = f(t + dt/2, y + k1*dt/2)
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            system->particles[i].position = vec3_add(original_pos[i], 
                                                    vec3_scale(k1_pos[i], dt * 0.5));
            system->particles[i].velocity = vec3_add(original_vel[i], 
                                                    vec3_scale(k1_vel[i], dt * 0.5));
        }
    }
    
    sph_build_hash_table(system);
    sph_calculate_density(system);
    sph_calculate_pressure(system);
    sph_calculate_forces(system);
    sph_calculate_viscosity_forces(system);
    
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            k2_pos[i] = system->particles[i].velocity;
            k2_vel[i] = system->particles[i].acceleration;
        }
    }
    
    // Stage 3: k3 = f(t + dt/2, y + k2*dt/2)
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            system->particles[i].position = vec3_add(original_pos[i], 
                                                    vec3_scale(k2_pos[i], dt * 0.5));
            system->particles[i].velocity = vec3_add(original_vel[i], 
                                                    vec3_scale(k2_vel[i], dt * 0.5));
        }
    }
    
    sph_build_hash_table(system);
    sph_calculate_density(system);
    sph_calculate_pressure(system);
    sph_calculate_forces(system);
    sph_calculate_viscosity_forces(system);
    
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            k3_pos[i] = system->particles[i].velocity;
            k3_vel[i] = system->particles[i].acceleration;
        }
    }
    
    // Stage 4: k4 = f(t + dt, y + k3*dt)
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            system->particles[i].position = vec3_add(original_pos[i], 
                                                    vec3_scale(k3_pos[i], dt));
            system->particles[i].velocity = vec3_add(original_vel[i], 
                                                    vec3_scale(k3_vel[i], dt));
        }
    }
    
    sph_build_hash_table(system);
    sph_calculate_density(system);
    sph_calculate_pressure(system);
    sph_calculate_forces(system);
    sph_calculate_viscosity_forces(system);
    
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            k4_pos[i] = system->particles[i].velocity;
            k4_vel[i] = system->particles[i].acceleration;
        }
    }
    
    // Final update: y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    for (int i = 0; i < system->particle_count; i++) {
        if (system->particles[i].flags & PARTICLE_ACTIVE) {
            Vec3 pos_increment = vec3_add3(k1_pos[i], 
                                          vec3_scale(k2_pos[i], 2.0),
                                          vec3_scale(k3_pos[i], 2.0));
            pos_increment = vec3_add(pos_increment, k4_pos[i]);
            pos_increment = vec3_scale(pos_increment, dt / 6.0);
            
            Vec3 vel_increment = vec3_add3(k1_vel[i], 
                                          vec3_scale(k2_vel[i], 2.0),
                                          vec3_scale(k3_vel[i], 2.0));
            vel_increment = vec3_add(vel_increment, k4_vel[i]);
            vel_increment = vec3_scale(vel_increment, dt / 6.0);
            
            // Update particle state
            system->particles[i].prev_position = system->particles[i].position;
            system->particles[i].prev_velocity = system->particles[i].velocity;
            
            system->particles[i].position = vec3_add(original_pos[i], pos_increment);
            system->particles[i].velocity = vec3_add(original_vel[i], vel_increment);
            
            // Apply damping
            system->particles[i].velocity = vec3_scale(system->particles[i].velocity, 
                                                      DAMPING_FACTOR);
            
            // Check for event horizon crossing
            if (sph_check_event_horizon_crossing(system->particles[i].position, 
                                                system->black_hole)) {
                system->particles[i].flags &= ~PARTICLE_ACTIVE; // Remove particle
                printf("Particle %d crossed event horizon\n", i);
            }
        }
    }
    
    // Update simulation time
    system->current_time += dt;
    
cleanup:
    free(k1_pos);
    free(k1_vel);
    free(k2_pos);
    free(k2_vel);
    free(k3_pos);
    free(k3_vel);
    free(k4_pos);
    free(k4_vel);
    free(original_pos);
    free(original_vel);
}

//========================================
// Adaptive Parameters
//========================================

// Update smoothing lengths based on local density
void sph_update_smoothing_lengths(SPHSystem* system) {
    if (!system) return;
    
    const double TARGET_NEIGHBORS = 32.0;  // Target number of neighbors **Why 32 when MAX_Neighbours is 64?
    const double MIN_H = 0.1;              // Minimum smoothing length
    const double MAX_H = 2.0;              // Maximum smoothing length
    const double H_TOLERANCE = 0.05;       // Convergence tolerance
    const int MAX_ITERATIONS = 10;         // Maximum iterations per particle
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        double h_old = particle->smoothing_length;
        double h_new = h_old;
        
        // Iterative smoothing length adjustment
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            // Count actual neighbors within current smoothing length
            int actual_neighbors = 0;
            double volume_sum = 0.0;
            
            for (int j = 0; j < system->particle_count; j++) {
                if (i == j) continue;
                
                SPHParticle* neighbor = &system->particles[j];
                if (!(neighbor->flags & PARTICLE_ACTIVE)) continue;
                
                Vec3 dr = vec3_sub(particle->position, neighbor->position);
                double r = vec3_length(dr);
                
                if (r < h_new) {
                    actual_neighbors++;
                    // Accumulate kernel volume for density estimation
                    volume_sum += sph_wendland_c2_kernel(r, h_new);
                }
            }
            
            // Add self-contribution
            volume_sum += sph_wendland_c2_kernel(0.0, h_new);
            
            // Calculate density-based smoothing length
            // h ∝ (m/ρ)^(1/d) where d is dimension (3D)
            double density_estimate = particle->mass * volume_sum;
            if (density_estimate > MIN_DENSITY) {
                double h_density = pow(particle->mass / density_estimate, 1.0/3.0);
                h_density *= 1.2; // Scale factor for SPH kernel support
                
                // Blend density-based and neighbor-count-based estimates
                double neighbor_ratio = (double)actual_neighbors / TARGET_NEIGHBORS;
                double h_neighbors = h_new * pow(neighbor_ratio, -1.0/3.0);
                
                // Weighted average
                h_new = 0.7 * h_density + 0.3 * h_neighbors;
            } else {
                // Fallback: adjust based on neighbor count only
                double neighbor_ratio = (double)actual_neighbors / TARGET_NEIGHBORS;
                if (neighbor_ratio > 0.1) {
                    h_new = h_new * pow(neighbor_ratio, -1.0/3.0);
                } else {
                    h_new *= 1.1; // Increase if too few neighbors
                }
            }
            
            // Apply bounds
            if (h_new < MIN_H) h_new = MIN_H;
            if (h_new > MAX_H) h_new = MAX_H;
            
            // Check convergence
            double h_change = fabs(h_new - h_old) / h_old;
            if (h_change < H_TOLERANCE) break;
            
            h_old = h_new;
        }
        
        particle->smoothing_length = h_new;
        
        // Update neighbor count for performance tracking
        particle->neighbor_count = 0;
        for (int j = 0; j < system->particle_count; j++) {
            if (i == j) continue;
            
            SPHParticle* neighbor = &system->particles[j];
            if (!(neighbor->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle->position, neighbor->position);
            double r = vec3_length(dr);
            
            if (r < h_new && particle->neighbor_count < MAX_NEIGHBORS) {
                particle->neighbors[particle->neighbor_count] = j;
                particle->neighbor_count++;
            }
        }
    }
}

// Calculate adaptive time step based on CFL condition and force constraints
double sph_calculate_adaptive_time_step(SPHSystem* system) {
    if (!system) return 0.001; // Default fallback
    
    double dt_min = 1e6; // Start with large value
    
    // CFL condition parameters
    const double CFL_FACTOR = 0.3;      // CFL safety factor
    const double FORCE_FACTOR = 0.25;   // Force-based time step factor
    const double VISCOUS_FACTOR = 0.125; // Viscous time step factor
    const double MIN_DT = 1e-6;         // Minimum allowed time step
    const double MAX_DT = 0.01;         // Maximum allowed time step
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        double h = particle->smoothing_length;
        double rho = particle->density;
        
        // 1. CFL condition: dt < CFL * h / (c_s + |v|)
        double sound_speed = sqrt(system->gas_constant * particle->pressure / rho);
        double velocity_magnitude = vec3_length(particle->velocity);
        double signal_speed = sound_speed + velocity_magnitude;
        
        if (signal_speed > 1e-10) {
            double dt_cfl = CFL_FACTOR * h / signal_speed;
            if (dt_cfl < dt_min) dt_min = dt_cfl;
        }
        
        // 2. Force-based constraint: dt < sqrt(h / |a|)
        double acceleration_magnitude = vec3_length(particle->acceleration);
        if (acceleration_magnitude > 1e-10) {
            double dt_force = FORCE_FACTOR * sqrt(h / acceleration_magnitude);
            if (dt_force < dt_min) dt_min = dt_force;
        }
        
        // 3. Viscous constraint: dt < 0.125 * h^2 * ρ / μ
        if (system->viscosity_coeff > 1e-10) {
            double dt_viscous = VISCOUS_FACTOR * h * h * rho / system->viscosity_coeff;
            if (dt_viscous < dt_min) dt_min = dt_viscous;
        }
        
        // 4. Thermal diffusion constraint
        if (particle->thermal_diffusion > 1e-10) {
            double dt_thermal = 0.5 * h * h / particle->thermal_diffusion;
            if (dt_thermal < dt_min) dt_min = dt_thermal;
        }
        
        // 5. Black hole constraint (prevent particles from moving too far per step)
        if (system->black_hole) {
            double r = vec3_length(particle->position);
            if (r > 1e-10) {
                // Limit to moving at most 5% of current radius per time step
                double dt_bh = 0.05 * r / (velocity_magnitude + 1e-10);
                if (dt_bh < dt_min) dt_min = dt_bh;
                
                // Additional constraint near event horizon
                double r_s = 2.0 * system->black_hole->mass; // Schwarzschild radius approximation
                if (r < 5.0 * r_s) {
                    double proximity_factor = r / (5.0 * r_s);
                    dt_bh *= proximity_factor;
                    if (dt_bh < dt_min) dt_min = dt_bh;
                }
            }
        }
    }
    
    // Apply global bounds
    if (dt_min < MIN_DT) dt_min = MIN_DT;
    if (dt_min > MAX_DT) dt_min = MAX_DT;
    
    return dt_min;
}

// Apply boundary conditions to particles
void sph_apply_boundary_conditions(SPHSystem* system) {
    if (!system) return;
    
    // Boundary parameters
    const double RESTITUTION_COEFF = 0.1;    // Energy loss on collision
    const double FRICTION_COEFF = 0.9;       // Tangential friction
    const double BOUNDARY_DAMPING = 0.95;    // General damping near boundaries
    
    // Define simulation boundaries (adjust as needed)
    const double X_MIN = -50.0, X_MAX = 50.0;
    const double Y_MIN = -50.0, Y_MAX = 50.0;
    const double Z_MIN = -10.0, Z_MAX = 10.0;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        bool boundary_collision = false;
        Vec3 new_position = particle->position;
        Vec3 new_velocity = particle->velocity;
        
        // X boundaries
        if (particle->position.x < X_MIN) {
            new_position.x = X_MIN;
            if (particle->velocity.x < 0.0) {
                new_velocity.x = -particle->velocity.x * RESTITUTION_COEFF;
                new_velocity.y *= FRICTION_COEFF;
                new_velocity.z *= FRICTION_COEFF;
                boundary_collision = true;
            }
        } else if (particle->position.x > X_MAX) {
            new_position.x = X_MAX;
            if (particle->velocity.x > 0.0) {
                new_velocity.x = -particle->velocity.x * RESTITUTION_COEFF;
                new_velocity.y *= FRICTION_COEFF;
                new_velocity.z *= FRICTION_COEFF;
                boundary_collision = true;
            }
        }
        
        // Y boundaries
        if (particle->position.y < Y_MIN) {
            new_position.y = Y_MIN;
            if (particle->velocity.y < 0.0) {
                new_velocity.y = -particle->velocity.y * RESTITUTION_COEFF;
                new_velocity.x *= FRICTION_COEFF;
                new_velocity.z *= FRICTION_COEFF;
                boundary_collision = true;
            }
        } else if (particle->position.y > Y_MAX) {
            new_position.y = Y_MAX;
            if (particle->velocity.y > 0.0) {
                new_velocity.y = -particle->velocity.y * RESTITUTION_COEFF;
                new_velocity.x *= FRICTION_COEFF;
                new_velocity.z *= FRICTION_COEFF;
                boundary_collision = true;
            }
        }
        
        // Z boundaries (vertical confinement for disk)
        if (particle->position.z < Z_MIN) {
            new_position.z = Z_MIN;
            if (particle->velocity.z < 0.0) {
                new_velocity.z = -particle->velocity.z * RESTITUTION_COEFF;
                new_velocity.x *= FRICTION_COEFF;
                new_velocity.y *= FRICTION_COEFF;
                boundary_collision = true;
            }
        } else if (particle->position.z > Z_MAX) {
            new_position.z = Z_MAX;
            if (particle->velocity.z > 0.0) {
                new_velocity.z = -particle->velocity.z * RESTITUTION_COEFF;
                new_velocity.x *= FRICTION_COEFF;
                new_velocity.y *= FRICTION_COEFF;
                boundary_collision = true;
            }
        }
        
        // Apply boundary conditions
        if (boundary_collision) {
            particle->position = new_position;
            particle->velocity = vec3_scale(new_velocity, BOUNDARY_DAMPING);
            particle->flags |= PARTICLE_BOUNDARY;
            
            // Cool down particles that hit boundaries (energy dissipation)
            particle->temperature *= 0.9;
            particle->thermal_energy = particle->specific_heat * particle->temperature;
        } else {
            particle->flags &= ~PARTICLE_BOUNDARY;
        }
        
        // Black hole event horizon boundary
        if (system->black_hole) {
            double r = vec3_length(particle->position);
            double r_s = 2.0 * system->black_hole->mass; // Schwarzschild radius
            
            // Check for event horizon crossing
            if (sph_check_event_horizon_crossing(particle->position, system->black_hole)) {
                // Mark particle as accreted (remove from simulation)
                particle->flags &= ~PARTICLE_ACTIVE;
                particle->flags |= PARTICLE_ACCRETING;
                
                // Optional: add particle mass to black hole
                // system->black_hole->mass += particle->mass;
                
                continue;
            }
            
            // Apply strong damping very close to black hole
            if (r < 3.0 * r_s) {
                double proximity = r / (3.0 * r_s);
                double damping = 0.9 + 0.1 * proximity; // Stronger damping closer to BH
                particle->velocity = vec3_scale(particle->velocity, damping);
            }
        }
        
        // Periodic boundaries (optional - comment out if using reflecting boundaries)
        /*
        // X periodic
        if (particle->position.x < X_MIN) {
            particle->position.x += (X_MAX - X_MIN);
        } else if (particle->position.x > X_MAX) {
            particle->position.x -= (X_MAX - X_MIN);
        }
        
        // Y periodic
        if (particle->position.y < Y_MIN) {
            particle->position.y += (Y_MAX - Y_MIN);
        } else if (particle->position.y > Y_MAX) {
            particle->position.y -= (Y_MAX - Y_MIN);
        }
        */
        
        // Density-based boundary particles (for complex geometries)
        if (particle->flags & PARTICLE_BOUNDARY) {
            // Boundary particles have higher density and zero velocity
            particle->density = system->rest_density * 2.0;
            particle->velocity = (Vec3){0.0, 0.0, 0.0};
            particle->acceleration = (Vec3){0.0, 0.0, 0.0};
        }
    }
}

//========================================
// Thermal Dynamics
//========================================

// Calculate temperature from thermal energy and apply thermal diffusion
void sph_calculate_temperature(SPHSystem* system) {
    if (!system) return;
    
    const double MIN_TEMPERATURE = 10.0;    // Minimum temperature (K)
    const double MAX_TEMPERATURE = 50000.0; // Maximum temperature (K)
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        // Calculate temperature from thermal energy
        // E_thermal = c_p * T (for ideal gas)
        if (particle->specific_heat > 1e-10) {
            particle->temperature = particle->thermal_energy / particle->specific_heat;
        } else {
            // Fallback if specific heat is not set
            particle->temperature = particle->thermal_energy;
        }
        
        // Apply bounds
        if (particle->temperature < MIN_TEMPERATURE) {
            particle->temperature = MIN_TEMPERATURE;
            particle->thermal_energy = particle->specific_heat * particle->temperature;
        }
        if (particle->temperature > MAX_TEMPERATURE) {
            particle->temperature = MAX_TEMPERATURE;
            particle->thermal_energy = particle->specific_heat * particle->temperature;
        }
        
        // Apply thermal diffusion (from previous calculation)
        if (fabs(particle->thermal_diffusion) > 1e-10) {
            double dT_dt = particle->thermal_diffusion / particle->specific_heat;
            particle->temperature += dT_dt * system->dt;
            
            // Update thermal energy consistently
            particle->thermal_energy = particle->specific_heat * particle->temperature;
        }
        
        // Calculate thermal pressure contribution
        // P_thermal = ρ * R_specific * T
        double thermal_pressure = particle->density * 0.287 * particle->temperature / 1000.0; // Simplified gas constant
        
        // Update total pressure (already calculated in sph_calculate_pressure, but add thermal component)
        particle->pressure += thermal_pressure;
        
        // Update particle thermal flags
        if (particle->temperature > 5000.0) {
            particle->flags |= PARTICLE_HOT;
        } else {
            particle->flags &= ~PARTICLE_HOT;
        }
        
        // Calculate thermal velocity for sound speed
        double thermal_velocity = sqrt(particle->temperature / 1000.0); // Simplified
        
        // Update smoothing length based on thermal state (hotter = larger h)
        if (particle->flags & PARTICLE_HOT) {
            particle->smoothing_length *= 1.05; // Slightly larger smoothing for hot particles
        }
    }
}

// Apply radiative cooling to particles
void sph_apply_radiative_cooling(SPHSystem* system, double dt) {
    if (!system || dt <= 0.0) return;
    
    // Physical constants (in simulation units)
    const double STEFAN_BOLTZMANN = 5.67e-8;  // Stefan-Boltzmann constant
    const double OPACITY_BASE = 0.1;          // Base opacity
    const double COOLING_EFFICIENCY = 1.0;    // Cooling efficiency factor
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle = &system->particles[i];
        
        if (!(particle->flags & PARTICLE_ACTIVE)) continue;
        
        double temperature = particle->temperature;
        double density = particle->density;
        
        // Skip cooling for very cold particles
        if (temperature < 100.0) continue;
        
        // Calculate opacity (temperature and density dependent)
        double opacity = OPACITY_BASE * pow(density / system->rest_density, 0.5) * 
                        pow(temperature / 1000.0, -0.5);
        
        // Optically thin cooling: L = 4π * σ * T^4 * V * opacity
        // Optically thick cooling: L = σ * T^4 * A (surface area limited)
        
        // Calculate optical depth (simplified)
        double optical_depth = opacity * density * particle->smoothing_length;
        
        double cooling_rate;
        
        if (optical_depth < 1.0) {
            // Optically thin regime - free-free radiation
            cooling_rate = COOLING_EFFICIENCY * STEFAN_BOLTZMANN * pow(temperature, 4) * 
                          density * density * opacity;
        } else {
            // Optically thick regime - blackbody radiation from surface
            double surface_area = 4.0 * M_PI * particle->smoothing_length * particle->smoothing_length;
            cooling_rate = COOLING_EFFICIENCY * STEFAN_BOLTZMANN * pow(temperature, 4) * 
                          surface_area / particle->mass;
        }
        
        // Additional cooling mechanisms for accretion disks
        
        // Bremsstrahlung cooling (free-free radiation)
        if (temperature > 1000.0) {
            double bremsstrahlung_rate = 1.4e-27 * density * density * sqrt(temperature) / 1000.0;
            cooling_rate += bremsstrahlung_rate;
        }
        
        // Line cooling (simplified metal line cooling)
        if (temperature > 10000.0 && temperature < 100000.0) {
            double line_cooling_rate = 1e-22 * density * density * 
                                      exp(-11600.0 / temperature); // Approximate line cooling
            cooling_rate += line_cooling_rate;
        }
        
        // Synchrotron cooling (for very hot particles in magnetic field)
        if (particle->flags & PARTICLE_HOT && temperature > 10000.0) {
            // Simplified synchrotron cooling (assumes magnetic field)
            double magnetic_field_strength = 1.0; // Simplified B field
            double synchrotron_rate = 6.3e-18 * magnetic_field_strength * magnetic_field_strength * 
                                     temperature * temperature / (1000.0 * 1000.0);
            cooling_rate += synchrotron_rate;
        }
        
        // Apply cooling
        double energy_loss = cooling_rate * dt;
        
        // Limit energy loss to prevent negative temperatures
        double max_energy_loss = 0.9 * particle->thermal_energy;
        if (energy_loss > max_energy_loss) {
            energy_loss = max_energy_loss;
        }
        
        particle->thermal_energy -= energy_loss;
        
        // Update temperature
        if (particle->specific_heat > 1e-10) {
            particle->temperature = particle->thermal_energy / particle->specific_heat;
        } else {
            particle->temperature = particle->thermal_energy;
        }
        
        // Store cooling rate for diagnostics
        particle->radiative_cooling = cooling_rate;
        
        // Ensure minimum temperature
        if (particle->temperature < 10.0) {
            particle->temperature = 10.0;
            particle->thermal_energy = particle->specific_heat * particle->temperature;
        }
    }
}

// Apply viscous heating to particles
void sph_apply_viscous_heating(SPHSystem* system, double dt) {
    if (!system || dt <= 0.0) return;
    
    for (int i = 0; i < system->particle_count; i++) {
        SPHParticle* particle_i = &system->particles[i];
        
        if (!(particle_i->flags & PARTICLE_ACTIVE)) continue;
        
        double total_heating = 0.0;
        double h_i = particle_i->smoothing_length;
        
        // Calculate viscous heating from velocity shear
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue;
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                Vec3 dv = vec3_sub(particle_j->velocity, particle_i->velocity);
                Vec3 dr_unit = vec3_scale(dr, 1.0 / r);
                
                // Velocity shear components
                double v_radial = vec3_dot(dv, dr_unit);
                Vec3 v_tangential = vec3_sub(dv, vec3_scale(dr_unit, v_radial));
                double shear_rate = vec3_length(v_tangential) / r;
                
                // Viscous heating rate: Q = μ * (∇v)^2
                double viscosity = system->viscosity_coeff * particle_i->density;
                double heating_contribution = viscosity * shear_rate * shear_rate;
                
                // Weight by kernel
                double kernel_weight = sph_wendland_c2_kernel(r, h_avg);
                total_heating += heating_contribution * kernel_weight * particle_j->mass / particle_j->density;
            }
        }
        
        // Additional sources of viscous heating
        
        // Bulk viscosity heating (from compression/expansion)
        Vec3 velocity_gradient = {0.0, 0.0, 0.0};
        double divergence = 0.0;
        
        for (int n = 0; n < particle_i->neighbor_count; n++) {
            int j = particle_i->neighbors[n];
            if (j >= system->particle_count || i == j) continue;
            
            SPHParticle* particle_j = &system->particles[j];
            if (!(particle_j->flags & PARTICLE_ACTIVE)) continue;
            
            Vec3 dr = vec3_sub(particle_i->position, particle_j->position);
            double r = vec3_length(dr);
            
            if (r < 1e-10) continue;
            
            double h_avg = (h_i + particle_j->smoothing_length) * 0.5;
            
            if (r < h_avg) {
                Vec3 dv = vec3_sub(particle_j->velocity, particle_i->velocity);
                Vec3 grad_kernel = sph_calculate_kernel_gradient(particle_i->position, 
                                                               particle_j->position, h_avg, 0);
                
                // Velocity divergence
                divergence += vec3_dot(dv, grad_kernel) * particle_j->mass / particle_j->density;
            }
        }
        
        // Bulk viscosity heating
        double bulk_viscosity = system->viscosity_coeff * 0.5; // Typically smaller than shear viscosity
        double bulk_heating = bulk_viscosity * particle_i->density * divergence * divergence;
        total_heating += bulk_heating;
        
        // Tidal heating (for particles in strong gravitational field)
        if (system->black_hole) {
            double r = vec3_length(particle_i->position);
            if (r > 1e-10) {
                // Simplified tidal heating rate
                double tidal_field = system->black_hole->mass / (r * r * r);
                double tidal_heating = 0.01 * particle_i->density * tidal_field * tidal_field * 
                                      particle_i->smoothing_length * particle_i->smoothing_length;
                total_heating += tidal_heating;
            }
        }
        
        // Magnetic field heating (simplified MHD heating)
        if (particle_i->flags & PARTICLE_ACCRETING) {
            // Assume magnetic reconnection and field dissipation in accreting regions
            double magnetic_heating = 0.1 * particle_i->density * 
                                    vec3_length(particle_i->velocity) * vec3_length(particle_i->velocity);
            total_heating += magnetic_heating;
        }
        
        // Apply heating
        double energy_gain = total_heating * dt;
        
        particle_i->thermal_energy += energy_gain;
        
        // Update temperature
        if (particle_i->specific_heat > 1e-10) {
            particle_i->temperature = particle_i->thermal_energy / particle_i->specific_heat;
        } else {
            particle_i->temperature = particle_i->thermal_energy;
        }
        
        // Store heating rate for diagnostics
        particle_i->heating_rate = total_heating;
        
        // Limit maximum temperature to prevent numerical issues
        if (particle_i->temperature > 50000.0) {
            particle_i->temperature = 50000.0;
            particle_i->thermal_energy = particle_i->specific_heat * particle_i->temperature;
        }
        
        // Update thermal flags
        if (particle_i->temperature > 5000.0) {
            particle_i->flags |= PARTICLE_HOT;
        }
    }
}

//========================================
// Stub implementations for missing SPH functions
// These provide basic functionality to make the code compile and link
//========================================

