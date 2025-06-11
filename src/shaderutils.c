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
    
    // Kerr black hole parameters
    glUniform1f(glGetUniformLocation(program, "u_mass"), (float)params.mass);
    glUniform1f(glGetUniformLocation(program, "u_spin"), (float)params.spin);
    glUniform1f(glGetUniformLocation(program, "u_spin_dimensional"), (float)params.spin_dimensional);
    
    // Event horizons
    glUniform1f(glGetUniformLocation(program, "u_schwarzschild_radius"), (float)params.schwarzschild_radius);
    glUniform1f(glGetUniformLocation(program, "u_kerr_radius_outer"), (float)params.kerr_radius_outer);
    glUniform1f(glGetUniformLocation(program, "u_kerr_radius_inner"), (float)params.kerr_radius_inner);
    
    // Ergosphere
    glUniform1f(glGetUniformLocation(program, "u_ergosphere_eq"), (float)params.ergosphere_radius_eq);
    glUniform1f(glGetUniformLocation(program, "u_ergosphere_pole"), (float)params.ergosphere_radius_pole);
    
    // Critical orbits
    glUniform1f(glGetUniformLocation(program, "u_isco_prograde"), (float)params.isco_radius_prograde);
    glUniform1f(glGetUniformLocation(program, "u_isco_retrograde"), (float)params.isco_radius_retrograde);
    glUniform1f(glGetUniformLocation(program, "u_photon_sphere_prograde"), (float)params.photon_sphere_prograde);
    glUniform1f(glGetUniformLocation(program, "u_photon_sphere_retrograde"), (float)params.photon_sphere_retrograde);
    
    // Observer parameters
    glUniform1f(glGetUniformLocation(program, "u_observer_distance"), (float)params.observer_distance);
    glUniform1f(glGetUniformLocation(program, "u_observer_inclination"), (float)params.observer_inclination);
    
    // Integration parameters
    glUniform1f(glGetUniformLocation(program, "u_eps"), (float)params.eps);
    glUniform1f(glGetUniformLocation(program, "u_dtau"), (float)params.dtau);
    glUniform1i(glGetUniformLocation(program, "u_max_steps"), params.max_steps);
    
    // Screen parameters
    glUniform2f(glGetUniformLocation(program, "u_resolution"), (float)width, (float)height);
    
    // Camera setup
    double aspect_ratio = (double)width / (double)height;
    double fov = 60.0 * M_PI / 180.0;
    double theta = 0.0 * M_PI / 180.0;  // Angle from +x axis
    double r = params.observer_distance;
    
   // Camera position in spherical coordinates -> Cartesian
     float cam_pos[3] = {
        0.0f,
        (float)(r * sin(theta)),
        (float)(-r * cos(theta))
    };
    
    // Camera vectors (simplified - forward points toward origin)
    float cam_target[3] = {0.0f, 0.0f, 0.0f};
    float cam_up[3] = {0.0f, 1.0f, 0.0f}; //camera from direction x y z
    
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

    // Calculate up (right × forward) 
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

    printf("Camera pos: (%.2f, %.2f, %.2f)\n", cam_pos[0], cam_pos[1], cam_pos[2]);
    printf("Forward: (%.2f, %.2f, %.2f)\n", forward[0], forward[1], forward[2]);
    printf("Right: (%.2f, %.2f, %.2f)\n", right[0], right[1], right[2]);
    printf("Up: (%.2f, %.2f, %.2f)\n", up[0], up[1], up[2]);
    printf("theta: %.2f degrees\n", theta * 180.0 / M_PI);
    printf("Camera pos: (%.2f, %.2f, %.2f)\n", cam_pos[0], cam_pos[1], cam_pos[2]);
    printf("Camera distance from origin: %.2f\n", sqrt(cam_pos[0]*cam_pos[0] + cam_pos[1]*cam_pos[1] + cam_pos[2]*cam_pos[2]));
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

SSAAFramebuffer create_ssaa_fbo(int base_width, int base_height, int scale) {
    SSAAFramebuffer ssaa = {0};
    ssaa.width = base_width * scale;
    ssaa.height = base_height * scale;

    glGenFramebuffers(1, &ssaa.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, ssaa.fbo);

    glGenTextures(1, &ssaa.texture);
    glBindTexture(GL_TEXTURE_2D, ssaa.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, ssaa.width, ssaa.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaa.texture, 0);

    glGenRenderbuffers(1, &ssaa.rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, ssaa.rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, ssaa.width, ssaa.height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, ssaa.rbo);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        fprintf(stderr, "SSAA FBO is not complete!\n");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0); // Unbind
    return ssaa;
}