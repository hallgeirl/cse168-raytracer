#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"

#include "PointLight.h"
#include "Sphere.h"
#include "Plane.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"
#include "Phong.h"
#include "MiroWindow.h"

using namespace std;

void
makeSpiralScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    //g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setEye(Vector3(0, 0, -5));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, -15));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

    // create a spiral of spheres

    const int maxI = 150;
    const float a = 0.15f;
    for (int i = 1; i < maxI; i+=1)
    {
        float t = i/float(maxI);
        float theta = 4*PI*t;
        float r = a*theta;
        float x = r*cos(theta);
        float y = r*sin(theta);
        float z = 2*(2*PI*a - r);
        Material* mat = new Phong(Vector3(1.0f, t, i%2));
        Sphere * sphere = new Sphere;
        sphere->setCenter(Vector3(x,y,z));
        sphere->setRadius(r/10);
        sphere->setMaterial(mat);
        g_scene->addObject(sphere);
    }

    Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -2, 0));
    plane->setMaterial(new Lambert(Vector3(1.0, 0, 0)));
    g_scene->addObject(plane);

    TriangleMesh *mesh = new TriangleMesh();
    mesh->createSingleTriangle();

    mesh->setV1(Vector3(0,0,0));
    mesh->setV2(Vector3(0,3,0));
    mesh->setV3(Vector3(5,5,0));
    mesh->setN1(Vector3(0,0,-1));
    mesh->setN2(Vector3(0.1,0.1,-1).normalize());
    mesh->setN3(Vector3(-0.1,-0.2,-1).normalize());


    Triangle * triangle = new Triangle();
    triangle->setMesh(mesh);
    triangle->setIndex(0);
    triangle->setMaterial(new Lambert(Vector3(0,1,0)));
    g_scene->addObject(triangle);

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
makeSphereScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(-5, 1, 3));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 10));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

    Material* mat = new Phong(Vector3(1.0f, 0.5f, 0.25f), Vector3(0.1, 0.1, 0.1), Vector3(1, 1, 1), 10);
	mat->SetRefraction(1.0f, 1.5);
    Sphere * sphere = new Sphere;
    sphere->setCenter(Vector3(0,0,0));
    sphere->setRadius(1.5);
    sphere->setMaterial(mat);
    g_scene->addObject(sphere);

    Material* mat2 = new Phong(Vector3(0.25f, 0.5f, 0.75f), Vector3(0.1, 0.1, 0.1), Vector3(1, 1, 1), 20);
	mat2->SetReflection(0.25f);
	Sphere * sphere2 = new Sphere;
    sphere2->setCenter(Vector3(5,0,-1));
    sphere2->setRadius(2);
    sphere2->setMaterial(mat2);
    g_scene->addObject(sphere2);

	Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -3, 0));
    plane->setMaterial(new Lambert(Vector3(0.8, 0.8, 0.8), Vector3(0.1, 0.1, 0.1), 0.1f));
	mat->SetReflection(0.25f);
    g_scene->addObject(plane);

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void addModel(char* filename, Material *mat, Scene* scene)
{
	TriangleMesh * mesh = new TriangleMesh();
	mesh->load(filename);

    for (int i = 0; i < mesh->numTris(); i++)
    {
		Triangle *tri = new Triangle();
		tri->setMesh(mesh);
		tri->setIndex(i);
		tri->setMaterial(mat);
		g_scene->addObject(tri);
    }
}

void makeModelsScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(-5, 1, 3));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 10));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

    Material* bunnyMat = new Phong(Vector3(0.25f, 0.5f, 0.75f), Vector3(0.1, 0.1, 0.1), Vector3(1, 1, 1), 20);
	bunnyMat->SetReflection(0.25f);
	//bunnyMat->SetRefraction(1.0f, 1.5);
	addModel("models/bunny.obj", bunnyMat, g_scene);

/*
    Material* teapotMat = new Phong(Vector3(0.25f, 0.5f, 0.75f), Vector3(0.1, 0.1, 0.1), Vector3(1, 1, 1), 20);
	teapotMat->SetReflection(0.25f);
	addModel("models/teapot.obj", teapotMat, g_scene);
*/

	Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -3, 0));
    plane->setMaterial(new Lambert(Vector3(0.8, 0.8, 0.8), Vector3(0.1, 0.1, 0.1), 0.1f));
    g_scene->addObject(plane);

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

int
main(int argc, char*argv[])
{
    // create a scene
    //makeSpiralScene();
	//makeSphereScene();
	makeModelsScene();
    MiroWindow miro(&argc, argv);
    #ifndef NO_GFX
    miro.mainLoop();
    #else
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();
    #endif

    return 0; 
}

