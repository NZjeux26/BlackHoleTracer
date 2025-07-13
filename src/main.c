#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <SDL.h>
#include <SDL_image.h>
#include <sys/stat.h> 
#include <sys/types.h>
#include <GL/glew.h>
#include <SDL_opengl.h>
#include <omp.h>
#include "shaderutils.h"
#include "skybox.h"
#include "blackholemath.h"
#include "sph_sim.h"

int main() {
    int width = 1200; // Set the width of the window
    int height = 900; // Set the height of the window

    #ifdef _OPENMP
    printf("OpenMP IS available - compiled with OpenMP support\n");
    printf("Max threads: %d\n", omp_get_max_threads());
    #else
    printf("OpenMP NOT available - not compiled with OpenMP support\n");
    #endif

    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL initialization failed: %s\n", SDL_GetError());
        return 1;
    }

    //explictly set the openGL version to 4.10
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    
    //4x MSAA
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);

    SDL_Window* window = SDL_CreateWindow("Black Hole Raytracer - OpenGL", 
                                  SDL_WINDOWPOS_UNDEFINED, 
                                  SDL_WINDOWPOS_UNDEFINED, 
                                  width, height, 
                                  SDL_WINDOW_OPENGL);
    if (!window) {
        fprintf(stderr, "Could not create window: %s\n", SDL_GetError());
        return 1;
    }

    // Create OpenGL context
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (!gl_context) {
        fprintf(stderr, "Could not create OpenGL context: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    glEnable(GL_MULTISAMPLE);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "GLEW initialization failed\n");
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    printf("OpenGL Version: %s\n", glGetString(GL_VERSION));
    printf("GLSL Version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Replace the skybox loading code with:
    const char* cubemap_faces[6] = {
        "textures/px.jpg",  // positive x
        "textures/nx.jpg",   // negative x
        "textures/py.jpg",    // positive y
        "textures/ny.jpg", // negative y
        "textures/pz.jpg",  // positive z
        "textures/nz.jpg"    // negative z
    };
    
    GLuint skybox_texture = create_cubemap_texture(cubemap_faces);
    if (skybox_texture == 0) {
        fprintf(stderr, "Failed to create cubemap texture\n");
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    printf("Cubemap texture created successfully\n");
    
    // Load and compile shaders
    GLuint shader_program = create_shader_program("shaders/vertex.glsl", "shaders/fragment.glsl");
    if (!shader_program) {
        fprintf(stderr, "Failed to create shader program\n");
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Setup geometry
    GLuint VAO, VBO;
    setup_fullscreen_quad(&VAO, &VBO);

    //setup the BH parameters Mass(geometrix units), Spin(% speed of C), Distance (gemoetric units)
    // Note: In geometric units, mass is in terms of Schwarzschild radius (M = 1)
    // Spin is dimensionless (a/M), where -1 ≤ a/M ≤ 1, and distance is in terms of Schwarzschild radius.
    BlackHoleParams params = init_BH_params(1.0, 0.2, 30.0); // Mass and distance from black hole

    // Initialize SPH system for accretion disk
    SPHSystem* sph_system = sph_create_system(16348, &params);  // 8192 particles
    if (!sph_system) {
        fprintf(stderr, "Failed to create SPH system\n");
        // ... existing cleanup code ...
        return 1;
    }

    // Set up accretion disk particles
    printf("Initialising accretion disk particles...\n");
    sph_initialise_accretion_disk(sph_system, 6.0, 20.0, 16348);  // Inner radius: 6, Outer: 20, 4096 particles
    sph_initialise_keplerian_velocities(sph_system);
    sph_initialise_thermal_equilibrium(sph_system);

    printf("SPH system initialised with %d particles\n", sph_system->particle_count);

    for (int i = 0; i < 5 && i < sph_system->particle_count; i++) {
        printf("Particle %d: pos(%.2f,%.2f,%.2f) vel(%.2f,%.2f,%.2f)\n", 
            i, sph_system->particles[i].position.x, 
            sph_system->particles[i].position.y, 
            sph_system->particles[i].position.z,
            sph_system->particles[i].velocity.x,
            sph_system->particles[i].velocity.y,
            sph_system->particles[i].velocity.z);
    }

    SPHGPUData* gpu_data = create_sph_gpu_data(sph_system->max_particles);
    if (!gpu_data) {
        fprintf(stderr, "Failed to create SPH GPU data\n");
        // ... cleanup and return
    }

    // Upload particle data to GPU
    //upload_sph_particles_to_gpu(sph_system, gpu_data);
    
    // Initialize SDL_image for PNG saving
    if (IMG_Init(IMG_INIT_PNG) == 0) {
        fprintf(stderr, "IMG_Init failed: %s\n", IMG_GetError());
    }

    printf("Starting raytracing.....\n");
    clock_t start_time, end_time;
    double cpu_time_used;

    start_time = clock();

    // Set viewport
    glViewport(0, 0, width, height);
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        fprintf(stderr, "OpenGL Error after glViewport: 0x%X\n", err);
    }

    // Main loop
    SDL_Event event;
    bool running = true;
    bool save_image = true;
    bool should_render = true;
    bool particle_render = true;
    int particle_steps = 1000;
    float dt = 0.001f;
    
    // Main loop
    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
            if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_ESCAPE) {
                    running = false;
                }
            }
        }
        
        if(particle_render){
            printf("Particle Simulation Running!\n");
            
            #pragma omp parallel for schedule(static)
            for(int i = 0; i < particle_steps; i++){
                sph_update_system(sph_system,dt);
            }
            
            particle_render = false;
            printf("Particle Simulation has ran for %d steps\n", particle_steps);
            
            // Upload particle data to GPU
            upload_sph_particles_to_gpu(sph_system, gpu_data);
            printf("Uploaded %d particles to GPU\n", sph_system->particle_count);
            printf("GPU texture size: %d x %d\n", gpu_data->texture_size, gpu_data->texture_size);

            printf("New Particle data:\n");
            //debug
            for (int i = 0; i < 5 && i < sph_system->particle_count; i++) {
                printf("Particle %d: pos(%.2f,%.2f,%.2f) vel(%.2f,%.2f,%.2f)\n", 
                    i, sph_system->particles[i].position.x, 
                    sph_system->particles[i].position.y, 
                    sph_system->particles[i].position.z,
                    sph_system->particles[i].velocity.x,
                    sph_system->particles[i].velocity.y,
                    sph_system->particles[i].velocity.z);
            }
        }

        if(should_render) {
            // Clear screen
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            
            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after glClearColor: 0x%X\n", err);
            }

            // Use shader and set uniforms
            glUseProgram(shader_program);
            
            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after glUseProgram: 0x%X\n", err);
            }

            // Set shader uniforms
            set_shader_uniforms(shader_program, params, width, height);

            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after set_shader_uniforms: 0x%X\n", err);
            }

            // Bind skybox texture
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_texture);
            glUniform1i(glGetUniformLocation(shader_program, "u_skybox"), 0);

            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after glActiveTexture(GL_TEXTURE0): 0x%X\n", err);
            }

            glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_texture);
            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after binding skybox texture: 0x%X\n", err);
            }

            glUniform1i(glGetUniformLocation(shader_program, "u_skybox"), 0);
            err = glGetError();
            if (err != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error after setting u_skybox uniform: 0x%X\n", err);
            }

            // Bind disk texture
            bind_sph_textures_to_shader(shader_program, gpu_data, sph_system->particle_count);
            
            printf("Bound SPH textures to shader\n");
            GLint active_texture;
            glGetIntegerv(GL_ACTIVE_TEXTURE, &active_texture);
            printf("Active texture unit: GL_TEXTURE%d\n", active_texture - GL_TEXTURE0);
            
            // Render fullscreen quad
            glBindVertexArray(VAO);
            glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

            GLenum err;
            while ((err = glGetError()) != GL_NO_ERROR) {
                fprintf(stderr, "OpenGL Error: 0x%X\n", err);
            }

            // Save image if requested
            if (save_image) {
                save_framebuffer_to_png(width, height, "Images/blackhole_gpu.png");
                save_image = false;
                should_render = false; // Stop rendering after saving the image
    
                end_time = clock();
                cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
                printf("GPU raytracing completed in: %f seconds\n", cpu_time_used);
            }
            
            // Swap buffers
            SDL_GL_SwapWindow(window);
        }
        else SDL_Delay(60); // Wait if not rendering
         
    }
    // Cleanup
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shader_program);
    // Cleanup textures
    glDeleteTextures(1, &skybox_texture);
    //glDeleteTextures(1, &disk_texture);
     // Cleanup SPH system
    cleanup_sph_gpu_data(gpu_data);
    sph_destroy_system(sph_system);

    IMG_Quit();
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}