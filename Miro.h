#ifndef __MIRO_H__
#define __MIRO_H__

// #ifndef max
// #define max(a,b) ((a>b)?a:b)
// #endif

const float MIRO_TMAX = 1e12f;
const float epsilon   = 1e-4f;
const float PI = 3.1415926535897932384626433832795028841972f;
const float DegToRad = PI/180.0f;
const float RadToDeg = 180.0f/PI;
const float TRACE_DEPTH = 5;
const float TRACE_SAMPLES = 20;
const float PHOTON_MAX_DIST = 1e3;
const float PHOTON_SAMPLES = 100.f;
const float DOF_APERTURE = 0.15f;
const float DOF_FOCUS_PLANE = 20.23f;
//const float DOF_FOCUS_PLANE = 25.23f;


#include <stdlib.h>
#include "OpenGL.h"
#include <stdio.h>
#include <iostream>
#include <limits>

typedef struct tex_coord2d_s
{
	tex_coord2d_s() : u(0), v(0) {}
	tex_coord2d_s(float _u, float _v) : u(_u), v(_v) {}
	float u, v;
} tex_coord2d_t;

typedef struct tex_coord3d_s
{
	tex_coord3d_s() : u(0), v(0), w(0) {}
	tex_coord3d_s(float _u, float _v, float _w) : u(_u), v(_v), w(_w) {}
	float u, v, w;
} tex_coord3d_t;

class Ray;
class HitInfo;

class Object;
class Sphere;
class Triangle;
class TriangleMesh;
class Instance;

class PointLight;

class Camera;
class Image;
class Scene;
class Material;

extern void ParseFile(FILE* fp);
extern void initOpenGL();
extern Camera* g_camera;
extern Scene* g_scene;
extern Image* g_image;

const float infinity = std::numeric_limits<float>::infinity();

#endif
