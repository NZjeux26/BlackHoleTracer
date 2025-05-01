#ifndef SKYBOX_H
#define SKYBOX_H

#include <SDL2/SDL.h>
#include "raytracer.h"
// Struct to hold the skybox texture
typedef struct {
    SDL_Surface* surface;
    int width;
    int height;
} Skybox;

Skybox load_skybox(const char* filename);
Uint32 sample_skybox(Skybox skybox, int x, int y);
void direction_to_texcoord(Vec3 dir, int* tex_x, int* tex_y, int tex_width, int tex_height);

#endif 