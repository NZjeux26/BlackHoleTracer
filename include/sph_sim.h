#ifndef SPH_SIM_H
#define SPH_SIM_H

#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <GL/glew.h>
#include <SDL_opengl.h>
#include "3DMathtypes.h"
#include "blackholemath.h"

/*Look at breaking this up into multiple smaller headers eg SPH_physics, SPH_Core, SPH_Blackhole, SPH_Utils etc*/

// SPH Configuration constants
#define MAX_PARTICLES            65536       // Maximum number of particles **This will need to be dropped in size heavily for testing
#define MAX_NEIGHBORS            64          // Maximum neighbors per particle
#define HASH_TABLE_SIZE          262144      // Spatial hash table size (power of 2)
#define GRID_CELL_SIZE           0.5         // Size of spatial grid cells
#define SPH_KERNEL_RADIUS        1.0
#define SPH_KERNEL_RADIUS_SQ     (SPH_KERNEL_RADIUS * SPH_KERNEL_RADIUS)
#define REST_DENSITY             1.0         // Rest density of fluid
#define GAS_CONSTANT             20.0        // Gas constant for pressure calculation
#define VISCOSITY_COEFF          0.1         // Viscosity coefficient
#define DAMPING_FACTOR           0.99        // Velocity damping factor
#define MIN_DENSITY              0.01        // Minimum density threshold
#define MAX_DENSITY              10.0        // Maximum density threshold
#define SURFACE_TENSION          0.0728      // Surface tension coefficient
#define THERMAL_CONDUCTIVITY     0.1         // Thermal conductivity coefficient

// Particle state flags
typedef enum {
    PARTICLE_ACTIVE     = 1 << 0,
    PARTICLE_BOUNDARY   = 1 << 1,
    PARTICLE_SURFACE    = 1 << 2,
    PARTICLE_INTERIOR   = 1 << 3,
    PARTICLE_HOT        = 1 << 4,
    PARTICLE_ACCRETING  = 1 << 5
} ParticleFlags;

// CPU-side SPH particle structure
typedef struct {
    // Kinematics
    Vec3 position;
    Vec3 velocity;
    Vec3 acceleration;
    Vec3 prev_position;
    Vec3 prev_velocity;

    // Physical properties
    double mass;
    double density;
    double pressure;
    double temperature;
    double specific_heat;
    double viscosity;

    // SPH-specific
    double smoothing_length;
    Vec3 color_gradient;
    double color_laplacian;
    Vec3 surface_normal;
    double curvature;

    // Thermal
    double thermal_energy;
    double thermal_diffusion;
    double radiative_cooling;
    double heating_rate;

    // Simulation properties
    uint32_t flags;
    uint32_t id;
    int neighbor_count;
    uint32_t neighbors[MAX_NEIGHBORS];

    // Kerr black hole
    double orbital_velocity;
    double angular_momentum;
    double binding_energy;
    double r_coordinate;
    double theta_coordinate;
    double phi_coordinate;
} SPHParticle;

// Renderer-side packed particle for VBO upload
typedef struct {
    float x, y, z;      // World position
    float size;         // Point size or quad scale
    float r, g, b, a;   // Color
} RenderParticle;

// Spatial hash cell
typedef struct {
    uint32_t particle_indices[32];
    int count;
} HashCell;

// SPH system
typedef struct {
    // Simulation data
    SPHParticle* particles;
    int particle_count;
    int max_particles;

    // Spatial hashing
    HashCell* hash_table;
    int hash_table_size;
    double grid_cell_size;
    Vec3 world_min, world_max;

    // SPH params
    double kernel_radius;
    double kernel_radius_sq;
    double rest_density;
    double gas_constant;
    double viscosity_coeff;
    double surface_tension;
    double thermal_conductivity;

    // Time stepping
    double dt;
    double current_time;
    int max_iterations;
    bool adaptive_time_step;

    // Black hole reference
    BlackHoleParams* black_hole;

    // OpenGL render buffer
    GLuint render_vbo;

    // Performance counters
    double last_update_time;
    int neighbor_searches;
    int density_calculations;
} SPHSystem;

//========================================
// Function declarations
//========================================

// System management
SPHSystem* sph_create_system(int max_particles, BlackHoleParams* black_hole);
void sph_destroy_system(SPHSystem* system);
void sph_reset_system(SPHSystem* system);
void sph_update_system(SPHSystem* system, double dt);

// Particle management
int sph_add_particle(SPHSystem* system, Vec3 position, Vec3 velocity, double mass);
void sph_remove_particle(SPHSystem* system, int index);
void sph_set_particle_properties(SPHSystem* system, int index, double density,
                                 double temperature, uint32_t flags);

// Initialisation
void sph_initialise_accretion_disk(SPHSystem* system, double inner_radius,
                                  double outer_radius, int num_particles);
void sph_initialise_keplerian_velocities(SPHSystem* system);
void sph_initialise_thermal_equilibrium(SPHSystem* system);

// Spatial hashing
void sph_clear_hash_table(SPHSystem* system);
void sph_build_hash_table(SPHSystem* system);
uint32_t sph_hash_position(Vec3 position, double cell_size);
void sph_find_neighbors(SPHSystem* system, int particle_index);

// Kernel functions
double sph_wendland_c2_kernel(double r, double h);
double sph_wendland_c2_gradient(double r, double h);
double sph_wendland_c2_laplacian(double r, double h);
Vec3 sph_calculate_kernel_gradient(Vec3 ri, Vec3 rj, double h, int kernelType);
double sph_calculate_kernel_laplacian(Vec3 ri, Vec3 rj, double h, int kernelType);

// Physical calculations
void sph_calculate_density(SPHSystem* system);
void sph_calculate_pressure(SPHSystem* system);
void sph_calculate_forces(SPHSystem* system);
void sph_calculate_viscosity_forces(SPHSystem* system);
void sph_calculate_surface_tension_forces(SPHSystem* system);
void sph_calculate_thermal_diffusion(SPHSystem* system);

// Black hole physics
Vec3 sph_calculate_kerr_acceleration(Vec3 position, Vec3 velocity, BlackHoleParams* bh);
double sph_calculate_orbital_velocity_kerr(double r, double theta, BlackHoleParams* bh);
double sph_calculate_binding_energy(Vec3 position, Vec3 velocity, BlackHoleParams* bh);
bool sph_check_event_horizon_crossing(Vec3 position, BlackHoleParams* bh);

// Integration
void sph_integrate_rk4(SPHSystem* system, double dt);

// Adaptive parameters
void sph_update_smoothing_lengths(SPHSystem* system);
double sph_calculate_adaptive_time_step(SPHSystem* system);
void sph_apply_boundary_conditions(SPHSystem* system);

// Thermal dynamics
void sph_calculate_temperature(SPHSystem* system);
void sph_apply_radiative_cooling(SPHSystem* system, double dt);
void sph_apply_viscous_heating(SPHSystem* system, double dt);

// Rendering
void sph_initialize_renderer(SPHSystem* system);
void sph_update_render_buffer(SPHSystem* system, RenderParticle* buffer);
void sph_render_particles(SPHSystem* system);

// Coordinate transforms
Vec3 sph_cartesian_to_boyer_lindquist(Vec3 cartesian);
Vec3 sph_boyer_lindquist_to_cartesian(Vec3 bl_coords);
Vec3 sph_velocity_cartesian_to_bl(Vec3 vel_cart, Vec3 position);
Vec3 sph_velocity_bl_to_cartesian(Vec3 vel_bl, Vec3 position);

// Utility
void sph_print_system_stats(SPHSystem* system);
void sph_export_particles_to_file(SPHSystem* system, const char* filename);
void sph_import_particles_from_file(SPHSystem* system, const char* filename);
double sph_get_kinetic_energy(SPHSystem* system);
double sph_get_potential_energy(SPHSystem* system);
double sph_get_thermal_energy(SPHSystem* system);

// Debug
bool sph_validate_particle_state(SPHSystem* system, int particle_index);
void sph_check_conservation_laws(SPHSystem* system);
void sph_detect_numerical_instabilities(SPHSystem* system);

#endif // SPH_H
