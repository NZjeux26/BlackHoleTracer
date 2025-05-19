#ifndef TYPES_H
#define TYPES_H

typedef struct {
    double x, y, z;
} Vec3;

Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_add3(Vec3 a, Vec3 b, Vec3 c);
Vec3 vec3_sub(Vec3 a, Vec3 b);
Vec3 vec3_scale(Vec3 v, double s);
double vec3_length(Vec3 v);
Vec3 vec3_normalise(Vec3 v);
double vec_dot(Vec3 a, Vec3 b);
Vec3 vec3_cross(Vec3 a, Vec3 b);   
Vec3 vec3_reflect(Vec3 i, Vec3 normal);

#endif // TYPES_H