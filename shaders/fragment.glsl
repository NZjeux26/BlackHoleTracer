#version 410 core

out vec4 FragColour;
in vec2 TexCoord;
#define PI 3.1415926538

// Black hole parameters
uniform float u_mass;
uniform float u_schwarzschild_radius; //not needed anymore
uniform float u_observer_distance;
uniform float u_dtau;
uniform float u_eps;
uniform int u_max_steps;
uniform float u_spin;  // Kerr spin parameter (a)

// Accretion disk parameters
uniform float u_disk_inner_radius;
uniform float u_disk_outer_radius;
uniform float u_disk_opacity;
uniform float u_disk_temperature_factor;
uniform float u_disk_thickness;

// Camera parameters
uniform vec3 u_cam_pos;
uniform vec3 u_cam_forward;
uniform vec3 u_cam_up;
uniform vec3 u_cam_right;
uniform float u_fov;
uniform float u_aspect;

// Screen parameters
uniform vec2 u_resolution;

// Skybox texture
uniform samplerCube u_skybox;
uniform sampler2D u_disk_texture;  // Accretion disk texture

//Cconstants
float M = u_mass; // Mass of the black hole
float a = u_spin * M; //Spin parameter spin amount * Mass of the black hole

//////////////////////////////////////////////////////////////
// Kerr Metric in Kerr-Schild Coordinates
//////////////////////////////////////////////////////////////

///////////////
// Utility Functions
/////////////////

/*Computes the inverse of a 4x4 matrix*/
mat4 diag(vec4 v) {
    return mat4(v.x,0,0,0, 0,v.y,0,0, 0,0,v.z,0, 0,0,0,v.w);
}
///Normalizes a 4D vector using the metric tensor g, ensuring it has unit length
vec4 unit(vec4 v, mat4 g) {
    float norm2 = dot(g * v, v);
    return (norm2 != 0.0) ? v / sqrt(abs(norm2)) : v;
}

/////////////////////
// Core Physics Functions
/////////////////////

/*Converts 4D Kerr-Schild coordinates to the radial coordinate r used in the Kerr metric*/
float rFromCoords(vec4 pos) {
    vec3 p = pos.yzw;
    float rho2 = dot(p,p) - a*a;
    float r2 = 0.5 * (rho2 + sqrt(rho2*rho2 + 4.0*a*a*p.z*p.z));
    return sqrt(r2);
}
/*Computes the Kerr spacetime metric tensor at a given 4D position,
describing the curvature of spacetime around the rotating black hole*/
mat4 metric(vec4 pos) {
    float r = rFromCoords(pos);
    vec4 k = vec4(-1.0, (r*pos.y - a*pos.z)/(r*r + a*a),
                        (r*pos.z + a*pos.y)/(r*r + a*a),
                         pos.a / r);
    float f = 2.0 * M * r / (r*r + a*a * pos.a*pos.a / (r*r));
    return f * mat4(k.x*k, k.y*k, k.z*k, k.w*k) + diag(vec4(-1,1,1,1));
}
/*Calculates the Hamiltonian (total energy) of a photon at position x with momentum p, 
used for geodesic integration*/
float hamiltonian(vec4 x, vec4 p) {
    return 0.5 * dot(inverse(metric(x)) * p, p);
}

/////////////////
// Geodesic Integration Functions
/////////////////

/*Computes the numerical gradient of the Hamiltonian 
using finite differences for the equations of motion*/
vec4 hamiltonianGradient(vec4 x, vec4 p) {
    float eps = u_eps;
    return (vec4(
        hamiltonian(x + vec4(eps,0,0,0), p),
        hamiltonian(x + vec4(0,eps,0,0), p),
        hamiltonian(x + vec4(0,0,eps,0), p),
        hamiltonian(x + vec4(0,0,0,eps), p)
    ) - hamiltonian(x, p)) / eps;
}

/*Performs one step of Hamiltonian integration to advance the photon's 
position and momentum along its geodesic path*/
void transportStep(inout vec4 x, inout vec4 p) {
    float dtau = u_dtau;
    p -= dtau * hamiltonianGradient(x, p);
    x += dtau * inverse(metric(x)) * p;
}

/*Determines when to stop raytracing (either the photon hits the event horizon or escapes to infinity)*/
bool stopCondition(vec4 pos) {
    float r = rFromCoords(pos);
    float horizon = M + sqrt(M*M - a*a);
    return r < horizon || r > u_observer_distance * 10.0;
}

/////////////// 
// Rendering Functions & Coordinate Systems
/////////////////

/*Constructs a tetrad (orthonormal basis) for the observer's frame in Kerr spacetime,
given the observer's position, time direction, aim direction, and vertical direction.
This tetrad is used to convert between 4D spacetime coordinates and 3D spatial coordinates.*/
mat4 tetrad(vec4 x, vec4 time, vec4 aim, vec4 vert) {
    mat4 g = metric(x);
    vec4 E0 = unit(time, g);
    vec4 E1 = unit(aim + dot(g*aim, E0) * E0, g);
    vec4 E3 = unit(vert - dot(g*vert,E1)*E1 + dot(g*vert,E0)*E0, g);
    vec4 E2 = unit(inverse(g) * vec4(
        dot(E0.yzw, cross(E1.yzw, E3.yzw)),
        -dot(E0.zwx, cross(E1.zwx, E3.zwx)),
         dot(E0.wxy, cross(E1.wxy, E3.wxy)),
        -dot(E0.xyz, cross(E1.xyz, E3.xyz))), g);
    mat4 basis;
    basis[0] = E0;
    basis[1] = E1;
    basis[2] = E2;
    basis[3] = E3;
    return basis;
}

// Sample the spherical skybox using ray direction
vec3 sample_skybox(vec3 direction) {
    // Debug: Return pure red to verify function is being called
    //return vec3(1.0, 0.0, 0.0);//this is actually happening, or something is overiding it.
    return texture(u_skybox, normalize(direction)).rgb;
}

// Calculate disk colour with temperature and Doppler effects
// Calculate disk color with temperature and Doppler effects
vec3 calculateDiskColour(vec2 diskUV, float blueshift, float r) {
    // Base disk texture
    vec3 baseColour = texture(u_disk_texture, diskUV).rgb;

    // Temperature-based emission (simplified blackbody)
    float temp_factor = u_disk_temperature_factor / pow(r / u_disk_inner_radius, 0.75);
    vec3 thermal_colour = vec3(1.0, 0.8, 0.6) * temp_factor;

    // Apply Doppler shift (relativistic beaming)
    vec3 disk_colour = mix(baseColour, thermal_colour, 0.7);
    disk_colour *= pow(abs(blueshift), 3.0); // Relativistic beaming

    return disk_colour;
}

// Check if a position is within the accretion disk
// Uses your existing AccretionDisk parameters
bool isInDisk(vec4 pos) {
    float r = rFromCoords(pos);
    return abs(pos.w) < u_disk_thickness * 0.5 && 
           r >= u_disk_inner_radius && 
           r <= u_disk_outer_radius;
}

/*Main raytracing function that follows a light ray through curved spacetime 
and returns the colour from the skybox or black for captured rays*/
vec3 trace_kerr_ray(vec3 dir, vec4 camPos, mat4 axes) {
    vec4 pos = camPos;
    vec4 dir4D = -axes[0] + vec4(0.0, dir.x, dir.y, dir.z); //changing to this fixed the rotation issue with the camera
    vec4 p = metric(pos) * dir4D;

    bool captured = false;
    bool hitDisk = false;
    vec4 intersectPos;
    vec2 diskUV;
    float blueshift;
    vec4 final_pos;

    for (int i = 0; i < u_max_steps; i++) {
        vec4 last_pos = pos;
        transportStep(pos, p);

        if(pos.w * last_pos.w < 0.0){
            // Interpolate to find exact intersection point
            intersectPos = (pos * abs(last_pos.w) + last_pos * abs(pos.w)) / (abs(last_pos.w) + abs(pos.w));
            float r = rFromCoords(intersectPos);

            // Check if intersection is within disk bounds
            if (r > u_disk_inner_radius && r < u_disk_outer_radius) {
                hitDisk = true;

                // Calculate UV coordinates for disk texture
                diskUV = (intersectPos.yz / u_disk_outer_radius + 1.0) * 0.5;

                // Calculate disk velocity for Doppler shift
                // Simplified Keplerian velocity in Kerr spacetime
                float sqrt_r = sqrt(r);
                vec4 diskVel = vec4(r + a/sqrt_r, 
                                   vec3(-intersectPos.z, intersectPos.y, 0.0) * sign(a) / sqrt_r) / 
                               sqrt(r*r - 3.0*r + 2.0*a*sqrt_r);

                // Calculate Doppler factor
                blueshift = 1.0 / dot(p, diskVel);
                break;
            }
        }
        if (stopCondition(pos)) {
            float r = rFromCoords(pos);
            captured = r < M + sqrt(M*M - a*a);
            break;
        }
        final_pos = pos;
    }

    if(hitDisk){
        float r = rFromCoords(intersectPos);
        return calculateDiskColour(diskUV, blueshift, r);
    }

    else{
        // Sample skybox direction
        vec4 out_dir = inverse(metric(final_pos)) * p;
        vec3 cube_dir = normalize(vec3(-out_dir.y, out_dir.w, -out_dir.z));//this controls the background skybox

        // In trace_kerr_ray, before the return:
        return captured ? vec3(0.0) : sample_skybox(cube_dir);
    }
}

void main() {
    // Calculate ray direction
    vec2 screen_pos = TexCoord * 2.0 - 1.0;
    screen_pos.y = -screen_pos.y;
    
    float scale = tan(u_fov * 0.5);
    float screen_x = screen_pos.x * u_aspect * scale;
    float screen_y = screen_pos.y * scale;
    
    vec3 ray_dir = normalize(
        u_cam_right * screen_x + 
        u_cam_up * screen_y + 
        u_cam_forward
    );
    
    // Observer 4-position and initial time direction
    // Note: u_cam_pos is in the format (t, x, y, z) where t is time
    // and (x, y, z) are the spatial coordinates.
    vec4 camPos = vec4(0.0, u_cam_pos.x, u_cam_pos.y, u_cam_pos.z);
    // Map 3D camera vectors to 4D space
    vec4 aim = vec4(0.0, u_cam_forward.x, u_cam_forward.y, u_cam_forward.z);
    vec4 vert = vec4(0.0, u_cam_up.x, u_cam_up.y, u_cam_up.z);
    vec4 timeDir = vec4(1.0, 0.0, 0.0, 0.0);
    
    mat4 camFrame = tetrad(camPos, timeDir, aim, vert);

    vec3 colour = trace_kerr_ray(ray_dir, camPos, camFrame);
    FragColour = vec4(clamp(colour, 0.0, 1.0), 1.0);
}