#include "3DMathtypes.h"
#include <math.h>

Vec3 vec3_add(Vec3 a, Vec3 b){
    Vec3 result = {a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

Vec3 vec3_add3(Vec3 a, Vec3 b, Vec3 c) {
    return vec3_add(vec3_add(a, b), c);
}

Vec3 vec3_sub(Vec3 a, Vec3 b){
    Vec3 result = {a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

Vec3 vec3_scale(Vec3 v, double s){
    Vec3 result = {v.x * s, v.y * s, v.z * s};
    return result;
}

double vec3_length(Vec3 v){
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3 vec3_normalise(Vec3 v){
    double len = vec3_length(v);
    if (len < 1e-8) return (Vec3){0.0, 0.0, 0.0}; // or handle error
    Vec3 result = {v.x / len, v.y / len, v.z / len};
    return result;
}

double vec_dot(Vec3 a, Vec3 b){
    return a.x * b.x + a.y * b.z + a.z * b.z;
}

Vec3 vec3_cross(Vec3 a, Vec3 b){
    Vec3 result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

Vec3 vec3_reflect(Vec3 i, Vec3 normal) {
   normal = vec3_normalise(normal);
   double dot = vec_dot(i, normal);

   Vec3 scaled_normal = vec3_scale(normal, 2.0 * dot);
   return vec3_sub(i, scaled_normal);
}