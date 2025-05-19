#include "3DMathtypes.h"
#include "camera.h"
#include <math.h>
#include <stdio.h>

Camera make_camera(Vec3 position, Vec3 target, Vec3 world_up, double fov, double aspect) {
    Camera cam;
    cam.position = position;

    // Forward points from position to target
    cam.forward = vec3_normalise(vec3_sub(target, position));

    // If forward is nearly parallel to world_up, pick an alternative up vector
    if (fabs(vec_dot(cam.forward, world_up)) > 0.999) {
        world_up = (Vec3){1.0, 0.0, 0.0}; // Switch to X-axis as new up
    }

    // Build camera basis
    cam.right = vec3_normalise(vec3_cross(cam.forward, world_up));
    cam.up = vec3_normalise(vec3_cross(cam.right, cam.forward)); // Now orthogonal

    cam.fov = fov;
    cam.aspect = aspect;
    return cam;
}