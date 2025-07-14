#ifndef SHADERUTILS_H
#define SHADERUTILS_H

#include <SDL_image.h>
#include <GL/glew.h>
#include <SDL_opengl.h>
#include "blackholemath.h"
#include "sph_sim.h"

typedef struct {
    GLuint fbo;
    GLuint texture;
    GLuint rbo;
    int width, height;
} SSAAFramebuffer;

// Structure to hold particle GPU textures
typedef struct {
    GLuint positions_texture;     // RGBA32F: xyz = position, w = mass
    GLuint velocities_texture;    // RGBA32F: xyz = velocity, w = density  
    GLuint properties_texture;    // RGBA32F: x = temperature, y = pressure, z = smoothing_length, w = flags
    GLuint thermal_texture;       // RGBA32F: x = thermal_energy, y = radiative_cooling, z = heating_rate, w = unused
    int texture_size;             // Size of square texture (must be >= sqrt(max_particles))
} SPHGPUData;

char* load_shader_source(const char* filepath);
GLuint compile_shader(const char* source, GLenum shader_type);
GLuint create_shader_program(const char* vertex_path, const char* fragment_path);
GLuint create_simple_disk_texture();
void setup_fullscreen_quad(GLuint* VAO, GLuint* VBO);
void set_shader_uniforms(GLuint program, BlackHoleParams params, int width, int height);
void save_framebuffer_to_png(int width, int height, const char* filename);
SSAAFramebuffer create_ssaa_fbo(int base_width, int base_height, int scale);
// Create GPU textures for SPH particle data
SPHGPUData* create_sph_gpu_data(int max_particles);
// Upload SPH particle data to GPU textures
void upload_sph_particles_to_gpu(SPHSystem* sph_system, SPHGPUData* gpu_data);
// Bind SPH textures to shader
void bind_sph_textures_to_shader(GLuint shader_program, SPHGPUData* gpu_data, int particle_count);
void cleanup_sph_gpu_data(SPHGPUData* gpu_data);
void upload_hash_table_to_gpu(SPHSystem* sph_system, SPHGPUData* gpu_data);
#endif 
