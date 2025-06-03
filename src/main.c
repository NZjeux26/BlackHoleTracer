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
#include "shaderutils.h"
#include "skybox.h"
#include "blackholemath.h"

int main() {
    int width = 1920; // Set the width of the window
    int height = 1200; // Set the height of the window

    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL initialization failed: %s\n", SDL_GetError());
        return 1;
    }

    //explictly set the openGL version to 4.10
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

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

    //setup the BH parameters
    BlackHoleParams params = init_BH_params(1.0, 30.0); // Mass and distance from black hole

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
    
    // Main loop
    SDL_Event event;
    bool running = true;
    bool save_image = true;
    bool shouldd_render = true;
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
        if(shouldd_render) {
            // Clear screen
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            
            // Use shader and set uniforms
            glUseProgram(shader_program);
            // Set shader uniforms
            set_shader_uniforms(shader_program, params, width, height);
            // Bind skybox texture
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_CUBE_MAP, skybox_texture);
            glUniform1i(glGetUniformLocation(shader_program, "u_skybox"), 0);

            // Render fullscreen quad
            glBindVertexArray(VAO);
            glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
            
            // Save image if requested
            if (save_image) {
                save_framebuffer_to_png(width, height, "Images/blackhole_gpu.png");
                save_image = false;
                shouldd_render = false; // Stop rendering after saving the image

                end_time = clock();
                cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
                printf("GPU raytracing completed in: %f seconds\n", cpu_time_used);
            }
            
            // Swap buffers
            SDL_GL_SwapWindow(window);
        }
        else SDL_Delay(100); // Wait if not rendering
         
    }
    // Cleanup
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shader_program);
    
    IMG_Quit();
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}