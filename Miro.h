#ifndef __MIRO_H__
#define __MIRO_H__

// #ifndef max
// #define max(a,b) ((a>b)?a:b)
// #endif

const float MIRO_TMAX = 1e12f;
const float epsilon   = 1e-5;
const float PI = 3.1415926535897932384626433832795028841972f;
const float DegToRad = PI/180.0f;
const float RadToDeg = 180.0f/PI;
const float TRACE_DEPTH = 5;

#include <stdlib.h>
#include "OpenGL.h"
#include <stdio.h>
#include <iostream>

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

#endif
