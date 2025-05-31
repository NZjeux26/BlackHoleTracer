#ifndef SKYBOX_H
#define SKYBOX_H

#include <SDL.h>
#include <SDL_image.h>
#include <GL/glew.h>

GLuint create_cubemap_texture(const char* faces[6]);

#endif

