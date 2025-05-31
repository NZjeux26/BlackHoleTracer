#ifndef CAMERA_H
#define CAMERA_H

#include "3DMathtypes.h"

typedef struct {
    Vec3 position;
    Vec3 forward;
    Vec3 up;
    Vec3 right;
    double fov;     // Field of view (in radians)
    double aspect;  // Aspect ratio = width / height
} Camera;

Camera make_camera(Vec3 position, Vec3 target, Vec3 world_up, double fov, double aspect);

#endif 