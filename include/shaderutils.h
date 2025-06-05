#ifndef SHADERUTILS_H
#define SHADERUTILS_H

#include <SDL_image.h>
#include <GL/glew.h>
#include <SDL_opengl.h>
#include "blackholemath.h"

typedef struct {
    GLuint fbo;
    GLuint texture;
    GLuint rbo;
    int width, height;
} SSAAFramebuffer;

char* load_shader_source(const char* filepath);
GLuint compile_shader(const char* source, GLenum shader_type);
GLuint create_shader_program(const char* vertex_path, const char* fragment_path);
void setup_fullscreen_quad(GLuint* VAO, GLuint* VBO);
void set_shader_uniforms(GLuint program, BlackHoleParams params, int width, int height);
void save_framebuffer_to_png(int width, int height, const char* filename);
SSAAFramebuffer create_ssaa_fbo(int base_width, int base_height, int scale);

#endif 