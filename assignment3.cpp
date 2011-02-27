#include "assignment3.h"
#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Sphere.h"
#include "PointLight.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"

// local helper function declarations
namespace
{
void addMeshTrianglesToScene(TriangleMesh * mesh, Material * material);
}


void
makeTestPetalScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(2, 3, -2));
    g_camera->setLookAt(Vector3(0, 0, 6));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(10, 10, 10));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(700);
    g_scene->addLight(light);

	 // create and place a point light source
    PointLight * amlight = new PointLight;
    amlight->setPosition(Vector3(-10, 10, 10));
    amlight->setColor(Vector3(1, 1, 1));
    amlight->setWattage(300);
    g_scene->addLight(amlight);

	Material* material = new Lambert(Vector3(1.0f, 0.3f, 0.3f));
    TriangleMesh * teapot = new TriangleMesh;
    teapot->load("models/petal.obj");
	addMeshTrianglesToScene(teapot, material);
    
    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(-10, -1, -10));
    floor->setV2(Vector3(  0, -1,  10));
    floor->setV3(Vector3( 10, -1, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
	Material* floorMaterial = new Phong(Vector3(0.5f));
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(floorMaterial); 
    g_scene->addObject(t);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

// local helper function definitions
namespace
{

void
addMeshTrianglesToScene(TriangleMesh * mesh, Material * material)
{
    // create all the triangles in the bunny mesh and add to the scene
    for (int i = 0; i < mesh->numTris(); ++i)
    {
        Triangle* t = new Triangle;
        t->setIndex(i);
        t->setMesh(mesh);
        t->setMaterial(material); 
        g_scene->addObject(t);
    }
}
}