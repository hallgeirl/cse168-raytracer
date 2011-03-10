#include "assignment1.h"
#include <math.h>
#include "Miro.h"
#include "includes.h"

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



//Make a scene with two spheres and a teapot. 
//Teapot is fully reflective, one sphere is texture mapped with the stone texture and one is part refractive, part reflective.
void makeScene1()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(256, 256);

    // set up the camera
    float viewAngleXZ = -PI;
    float pitch = -0.1;
    Vector3 dir(std::sin(viewAngleXZ),sin(pitch),std::cos(viewAngleXZ));
    Vector3 eye(0, 3, 2);
    
    g_scene->setBgColor(Vector3(.0f));
    g_camera->setEye(eye);
    g_camera->setLookAt(eye+dir);
    //g_camera->setEye(Vector3(-5, 20, 7));
    //g_camera->setLookAt(Vector3(-4, 1, 7));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(60);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-2, 3, -6));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    light = new PointLight();
    light->setPosition(Vector3(2, 4.5, -6.5));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 20, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 5, -7));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    Sphere  *sphere = new Sphere;
    sphere->setCenter(Vector3(-2,2.5,-9));
    sphere->setRadius(1.5);
    sphere->setMaterial(new Phong(Vector3(0.0f, 1.0f, 0.0f)));
    //sphere->setMaterial(new TexturedPhong(new StoneTexture(20), Vector3(0.1,0.1,0.1), Vector3(0,0,0), Vector3(0), 3));
    g_scene->addObject(sphere);

    sphere = new Sphere;
    sphere->setCenter(Vector3(2,2.5,-9));
    sphere->setRadius(1.5);
    sphere->setMaterial(new Phong(Vector3(1.0f, 0.0f, 0.f), Vector3(0.0), Vector3(0.0), 3, 1.5));
    g_scene->addObject(sphere);
    
    sphere = new Sphere;
    sphere->setCenter(Vector3(0,4.5,-10));
    sphere->setRadius(1.5);
    sphere->setMaterial(new Phong(Vector3(0,0,1.0), Vector3(0), Vector3(0.0, 0, 0.0), 3, 1.5));
    g_scene->addObject(sphere);

    Material* teapotMat = new Phong(Vector3(1, 1, 1), Vector3(0), Vector3(0), 3, 1.5);
    addModel("models/teapot.obj", teapotMat, g_scene, Vector3(0,0,-5));

	//g_scene->setEnvironment(autumnHDR);

	/*Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -0.5, 0));

    //plane->setMaterial(new TexturedPhong(new StoneTexture(3), Vector3(0.1,0.1,0.1), Vector3(0, 0, 0)));
    plane->
    plane->setMaterial(new Phong(Vector3(1, 0, 0)));
    g_scene->addObject(plane);*/
    addModel("models/square.obj", new Phong(Vector3(1,0,0)), g_scene, Vector3(0, 0, -8), 0, Vector3(6,6,6));

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

//Make a scene with different refractive spheres.
void makeScene2()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(2048, 2048);

    // set up the camera
    float viewAngleXZ = -PI;
    float pitch = 0;
    Vector3 dir(std::sin(viewAngleXZ),sin(pitch),std::cos(viewAngleXZ));
    Vector3 eye(0, 4, 2);
    
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(eye);
    g_camera->setLookAt(eye+dir);
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(60);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-2, 3, -6));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    light = new PointLight();
    light->setPosition(Vector3(2, 4.5, -4));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 20, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 5, -4));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    for (int y = 0; y < 3; y++)
    {
        for (int x = 0; x < 3; x++)
        {
            Sphere *sphere = new Sphere;
            sphere->setCenter(Vector3(3*(x-1), 3*y+1.5,-9));
            sphere->setRadius(1.5);
            sphere->setMaterial(new Phong(Vector3(), Vector3(0), Vector3(1), 10, 1.0+((float)y*3.0+(float)x*2.0)/20));
            g_scene->addObject(sphere);
        }
    }

    g_scene->setEnvironment(autumnHDR);

    Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -0.5, 0));
    plane->setMaterial(new TexturedPhong(new StoneTexture(3), Vector3(0, 0, 0)));
    g_scene->addObject(plane);

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}
 
//Make a scene with different refractive spheres.
void makeBUNNIZ()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    // set up the camera
    float viewAngleXZ = -PI-0.1;
    float pitch = -0.3;
    Vector3 dir(std::sin(viewAngleXZ),sin(pitch),std::cos(viewAngleXZ));
    Vector3 eye(0, 4, 3);
    
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(eye);
    g_camera->setLookAt(eye+dir);
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(60);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-2, 3, -6));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    light = new PointLight();
    light->setPosition(Vector3(2, 4.5, -4));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 20, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);
    
    light = new PointLight();
    light->setPosition(Vector3(0, 5, -4));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(30);
    g_scene->addLight(light);

    for (int x = 0; x < 3; x++)
    {
        Material* bunnyMat;
        float z = 0;
        switch (x)
        {
            case 0:
                bunnyMat = new Phong(Vector3(0.25, 0.5, 0.75), Vector3(0), Vector3(1), 3, 1.5);
                z += 2;
                break;
            case 1:
                bunnyMat = new Phong(Vector3(0.25, 0.5, 0.75), Vector3(0.7), Vector3(0.3), 3, 1.5);
                z -= 1;
                break;
            case 2:
                bunnyMat = new Phong(Vector3(0.25, 0.5, 0.75), Vector3(0), Vector3(0), 3, 1.5);
                z += 4;
                break;
        }
	    addModel("models/bunny.obj", bunnyMat, g_scene, Vector3(2*(x-1),0,-7+z), x*(PI/3));
    }

    g_scene->setEnvironment(autumnHDR);

    Plane * plane = new Plane();
    plane->setNormal(Vector3(0, 1, 0));
    plane->setOrigin(Vector3(0, -0.5, 0));
    plane->setMaterial(new TexturedPhong(new StoneTexture(3), Vector3(0, 0, 0)));
    g_scene->addObject(plane);

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
A1makeBunnyScene()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 3, 5));
    g_camera->setLookAt(Vector3(-.5, 1, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    g_scene->setEnvironment(autumnHDR);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(3500);
    g_scene->addLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    TriangleMesh * bunny = new TriangleMesh;
    bunny->load("models/bunny.obj");
    
    // create all the triangles in the bunny mesh and add to the scene
    for (int i = 0; i < bunny->numTris(); ++i)
    {
        Triangle* t = new Triangle;
        t->setIndex(i);
        t->setMesh(bunny);
        t->setMaterial(mat); 
        g_scene->addObject(t);
    }
    
    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  10));
    floor->setV2(Vector3( 10, 0, -10));
    floor->setV3(Vector3(-10, 0, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    //g_scene->addObject(t);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
A1makeSphereScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(2048, 2048);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 1, 5));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(500);
    g_scene->addLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    Sphere *sphere = new Sphere;
    sphere->setRadius(1.5);
    sphere->setMaterial(mat);

    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, -1.5,  10));
    floor->setV2(Vector3( 10, -1.5, -10));
    floor->setV3(Vector3(-10, -1.5, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    g_scene->addObject(t);
    g_scene->addObject(sphere);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
A1makeTeapotScene()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(256, 256);
    g_scene->setEnvironment(autumnHDR);
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 3, 5));
    g_camera->setLookAt(Vector3(0, 1, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(500);
    g_scene->addLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    TriangleMesh * teapot = new TriangleMesh;
    teapot->load("models/teapot.obj");
    
    // create all the triangles in the triangle mesh and add to the scene
    for (int i = 0; i < teapot->numTris(); ++i)
    {
        Triangle* t = new Triangle;
        t->setIndex(i);
        t->setMesh(teapot);
        t->setMaterial(mat); 
        g_scene->addObject(t);
    }
    
    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  10));
    floor->setV2(Vector3( 10, 0, -10));
    floor->setV3(Vector3(-10, 0, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    //g_scene->addObject(t);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}


