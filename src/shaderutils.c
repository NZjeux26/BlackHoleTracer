#include <SDL_image.h>
#include <GL/glew.h>
#include <SDL_opengl.h>
#include "blackholemath.h"
#include "shaderutils.h"

/*Loads the contents of a shader file into a dynamically allocated C string.
  Handles file I/O and ensures null-terminated source for GLSL compatibility.*/
char* load_shader_source(const char* filepath) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Could not open shader file: %s\n", filepath);
        return NULL;
    }
    
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    char* source = malloc(length + 1);
    if (!source) {
        fclose(file);
        return NULL;
    }
    
    fread(source, 1, length, file);
    source[length] = '\0';
    fclose(file);
    return source;
}

/*Compiles a shader of the given type from a source string.
Logs any compilation errors and returns 0 on failure.*/
GLuint compile_shader(const char* source, GLenum shader_type) {
    GLuint shader = glCreateShader(shader_type);
    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);
    
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char info_log[512];
        glGetShaderInfoLog(shader, 512, NULL, info_log);
        fprintf(stderr, "Shader compilation failed:\n%s\n", info_log);
        glDeleteShader(shader);
        return 0;
    }
    
    return shader;
}

/*Loads, compiles, and links a vertex and fragment shader into a GLSL program.
Cleans up shader resources and returns the final program object or 0 on error.*/
GLuint create_shader_program(const char* vertex_path, const char* fragment_path) {
    char* vertex_source = load_shader_source(vertex_path);
    char* fragment_source = load_shader_source(fragment_path);
    
    if (!vertex_source || !fragment_source) {
        free(vertex_source);
        free(fragment_source);
        return 0;
    }
    
    GLuint vertex_shader = compile_shader(vertex_source, GL_VERTEX_SHADER);
    GLuint fragment_shader = compile_shader(fragment_source, GL_FRAGMENT_SHADER);
    
    free(vertex_source);
    free(fragment_source);
    
    if (!vertex_shader || !fragment_shader) {
        if (vertex_shader) glDeleteShader(vertex_shader);
        if (fragment_shader) glDeleteShader(fragment_shader);
        return 0;
    }
    
    GLuint program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    
    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char info_log[512];
        glGetProgramInfoLog(program, 512, NULL, info_log);
        fprintf(stderr, "Shader program linking failed:\n%s\n", info_log);
        glDeleteProgram(program);
        program = 0;
    }
    
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    
    return program;
}

/*Initializes VAO, VBO, and EBO for a fullscreen quad with texture coordinates.
Binds vertex attributes for position and texture input to the shader.*/
void setup_fullscreen_quad(GLuint* VAO, GLuint* VBO) {
    // Fullscreen quad vertices (position + texture coordinates)
    float vertices[] = {
        // Positions   // Texture Coords
        -1.0f,  1.0f,  0.0f, 1.0f,  // Top Left
        -1.0f, -1.0f,  0.0f, 0.0f,  // Bottom Left
         1.0f, -1.0f,  1.0f, 0.0f,  // Bottom Right
         1.0f,  1.0f,  1.0f, 1.0f   // Top Right
    };
    
    unsigned int indices[] = {
        0, 1, 2,
        0, 2, 3
    };
    
    GLuint EBO;
    glGenVertexArrays(1, VAO);
    glGenBuffers(1, VBO);
    glGenBuffers(1, &EBO);
    
    glBindVertexArray(*VAO);
    
    glBindBuffer(GL_ARRAY_BUFFER, *VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Texture coordinate attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(0);
}

/*Sets uniforms for black hole rendering including camera position, orientation, and parameters.
Calculates view vectors and passes simulation settings to the shader.*/
void set_shader_uniforms(GLuint program, BlackHoleParams params, int width, int height) {
    glUseProgram(program);
    
    // Black hole parameters
    glUniform1f(glGetUniformLocation(program, "u_mass"), (float)params.mass);
    glUniform1f(glGetUniformLocation(program, "u_schwarzschild_radius"), (float)params.schwarzschild_radius);
    glUniform1f(glGetUniformLocation(program, "u_observer_distance"), (float)params.observer_distance);
    glUniform1f(glGetUniformLocation(program, "u_dt"), (float)params.dt);
    glUniform1i(glGetUniformLocation(program, "u_max_steps"), params.max_steps);
    
    // Screen parameters
    glUniform2f(glGetUniformLocation(program, "u_resolution"), (float)width, (float)height);
    
    // Camera setup
    double aspect_ratio = (double)width / (double)height;
    double fov = 50.0 * M_PI / 180.0; //FOV
    double theta = 0.0 * M_PI / 180.0;//camera angle
    double r = params.observer_distance;
    
    // Camera position
    float cam_pos[3] = {
        0.0f,
        (float)(r * sin(theta)),
        (float)(-r * cos(theta))
    };
    
    // Camera vectors (simplified - forward points toward origin)
    float cam_target[3] = {0.0f, 0.0f, 0.0f};
    float cam_up[3] = {0.5f, 1.0f, 0.0f};
    
    // Calculate forward vector (from camera to target)
    float forward[3] = {
        cam_target[0] - cam_pos[0],
        cam_target[1] - cam_pos[1],
        cam_target[2] - cam_pos[2]
    };
    
    // Normalize forward
    float forward_len = sqrt(forward[0]*forward[0] + forward[1]*forward[1] + forward[2]*forward[2]);
    forward[0] /= forward_len;
    forward[1] /= forward_len;
    forward[2] /= forward_len;
    
    // Calculate right vector (forward × up)
    float right[3] = {
        forward[1] * cam_up[2] - forward[2] * cam_up[1],
        forward[2] * cam_up[0] - forward[0] * cam_up[2],
        forward[0] * cam_up[1] - forward[1] * cam_up[0]
    };
    
    // Normalize right
    float right_len = sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);
    right[0] /= right_len;
    right[1] /= right_len;
    right[2] /= right_len;
    
    // Recalculate up (right × forward)
    float up[3] = {
        right[1] * forward[2] - right[2] * forward[1],
        right[2] * forward[0] - right[0] * forward[2],
        right[0] * forward[1] - right[1] * forward[0]
    };
    
    // Set camera uniforms
    glUniform3f(glGetUniformLocation(program, "u_cam_pos"), cam_pos[0], cam_pos[1], cam_pos[2]);
    glUniform3f(glGetUniformLocation(program, "u_cam_forward"), forward[0], forward[1], forward[2]);
    glUniform3f(glGetUniformLocation(program, "u_cam_up"), up[0], up[1], up[2]);
    glUniform3f(glGetUniformLocation(program, "u_cam_right"), right[0], right[1], right[2]);
    glUniform1f(glGetUniformLocation(program, "u_fov"), (float)fov);
    glUniform1f(glGetUniformLocation(program, "u_aspect"), (float)aspect_ratio);
}

/*Reads the framebuffer into memory, flips it vertically, and saves it as a PNG using SDL.
Handles memory allocation, surface manipulation, and error logging.*/
void save_framebuffer_to_png(int width, int height, const char* filename) {
    // Allocate buffer for pixel data
    unsigned char* pixels = malloc(width * height * 4);
    if (!pixels) {
        fprintf(stderr, "Failed to allocate memory for screenshot\n");
        return;
    }
    
    // Read pixels from framebuffer
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    
    // Create SDL surface from pixel data
    SDL_Surface* surface = SDL_CreateRGBSurfaceFrom(
        pixels, width, height, 32, width * 4,
        0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000
    );
    
    if (surface) {
        // Flip vertically (OpenGL has origin at bottom-left, images at top-left)
        SDL_Surface* flipped = SDL_CreateRGBSurface(0, width, height, 32,
            0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
        
        if (flipped) {
            SDL_LockSurface(surface);
            SDL_LockSurface(flipped);
            
            Uint32* src_pixels = (Uint32*)surface->pixels;
            Uint32* dst_pixels = (Uint32*)flipped->pixels;
            
            for (int y = 0; y < height; y++) {
                memcpy(&dst_pixels[y * width], 
                       &src_pixels[(height - 1 - y) * width], 
                       width * sizeof(Uint32));
            }
            
            SDL_UnlockSurface(flipped);
            SDL_UnlockSurface(surface);
            
            // Save flipped surface
            if (IMG_SavePNG(flipped, filename) != 0) {
                fprintf(stderr, "IMG_SavePNG failed: %s\n", IMG_GetError());
            } else {
                printf("Image saved as %s\n", filename);
            }
            
            SDL_FreeSurface(flipped);
        }
        SDL_FreeSurface(surface);
    }
    
    free(pixels);
}