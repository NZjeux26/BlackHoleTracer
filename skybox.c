#include <SDL2/SDL.h>
#include "raytracer.h"
#include "skybox.h"

Skybox load_skybox(const char* filename){
    Skybox skybox;

    skybox.surface = SDL_LoadBMP(filename);
    if(!skybox.surface){
        fprintf(stderr, "Failed to load skybox BMP image %s! SDL Error: %s\n", 
            filename, SDL_GetError());
        exit(1);
    }

    skybox.width = skybox.surface->w;
    skybox.height = skybox.surface->h;
}