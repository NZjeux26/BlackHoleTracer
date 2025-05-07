#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <SDL2/SDL.h>
#include "raytracer.h"

int main() {
    int width = 1200;
    int height = 800;

    SDL_Window* window = SDL_CreateWindow("Black Hole Raytracer", 
                                  SDL_WINDOWPOS_UNDEFINED, 
                                  SDL_WINDOWPOS_UNDEFINED, 
                                  width, height, 0);
    if (!window) {
        fprintf(stderr, "Could not create window: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        fprintf(stderr, "Could not create renderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        return 1;
    }

    SDL_Surface* surface = SDL_CreateRGBSurface(0, width, height, 32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000); //Co-Pilot wanted 0

    BlackHoleParams params = {
        .mass = 1.0,
        .schwarzschild_radius = 1.0,
        .accretion_disk_inner_radius = 3.0,
        .accretion_disk_outer_radius = 15.0,
        .accretion_disk_thickness = 0.5,
        .observer_distance = 100.0,
    };

    printf("Starting raytracing.....\n");
    clock_t start_time, end_time;
    double cpu_time_used;

    start_time = clock();

    raytrace_blackhole(params, surface);

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Raytracing finished.\n");
    printf("Time taken for raytrace_blackhole: %f seconds\n", cpu_time_used);

    // Create a texture from the surface
    SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);
    if (!texture) {
        fprintf(stderr, "Could not create texture: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        return 1;
    }
    // Main loop
    SDL_Event event;
    bool running = true;
   
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

        // Clear the screen
        SDL_RenderClear(renderer);

        // Render the texture
        SDL_RenderCopy(renderer, texture, NULL, NULL);

        // Present the back buffer
        SDL_RenderPresent(renderer);
    }
    // Clean up
    SDL_DestroyTexture(texture);
    SDL_FreeSurface(surface);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}