#include <math.h>
#include <string>
#ifdef OPENMP
#include <omp.h>
#endif
#include "assignment1.h"
#include "assignment2.h"
#include "assignment3.h"
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
#include "Texture.h"
#include <FreeImage.h>
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

void addModel(const char* filename, Material *mat, Scene* scene, Vector3 position, float rotY=0, Vector3 scale=Vector3(1,1,1))
{
	TriangleMesh * mesh = new TriangleMesh();
    Matrix4x4 m_rot(Vector4(cos(rotY),0,-sin(rotY),0), Vector4(0,1,0,0), Vector4(sin(rotY),0,cos(rotY),0), Vector4(0, 0, 0, 1));
    Matrix4x4 m_trans(Vector4(1,0,0,0), Vector4(0,1,0,0), Vector4(0,0,1,0), Vector4(position.x, position.y, position.z, 1));
    Matrix4x4 m_scale(Vector4(scale.x,0,0,0), Vector4(0,scale.y,0,0), Vector4(0,0,scale.z,0), Vector4(0, 0, 0, 1));
	mesh->load(filename, m_trans*m_rot*m_scale);
    for (int i = 0; i < mesh->numTris(); i++)
    {
		Triangle *tri = new Triangle();
		tri->setMesh(mesh);
		tri->setIndex(i);
		tri->setMaterial(mat);
		g_scene->addObject(tri);
    }
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
makeTestSphereScene()
{

	Texture *autumnHDR = new CloudTexture(3.0f,  //Scale (higher is smaller)
	                                        0.1f,  //Cloud size
	                                        0.2f, //Cloud density
	                                        50.0f,  //Sharpness
	                                        0.4f,  //Ambient (higher number means lower contrast)
	                                        0.35f,  //Shadow threshold
	                                        0.5f,  //Shadow magnitude
	                                        0.3f   //Shadow sharpness
	                                        ); 
    Material *m;
    cout << "Test scene" << endl;
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    //g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setEye(Vector3(9, 1, 0));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(90);

    

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(0, 5, -5));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

    light = new PointLight;
    light->setPosition(Vector3(0, 5, -25));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1500);
    g_scene->addLight(light);

	// focused sphere

    Sphere  *sphere = new Sphere;
    sphere->setCenter(Vector3(0,0.5,0));
    sphere->setRadius(3);
    sphere->setMaterial(m = new Phong(Vector3(0.0f, 1.0f, 0.0f), Vector3(0.1f), Vector3(0.0f), 10, 1.5));
    
    //g_scene->addObject(sphere);

	//close sphere

	Sphere  *sphereClose = new Sphere;
    sphereClose->setCenter(Vector3(2.5,1,-6));
    sphereClose->setRadius(1.5);
    sphereClose->setMaterial(m = new Phong(Vector3(1.0f, 0.0f, 0.0f), Vector3(0.1f), Vector3(0.0f), 10, 1.5));

   // g_scene->addObject(sphereClose);
	
	//far sphere

	Sphere  *sphereFar = new Sphere;
    sphereFar->setCenter(Vector3(-5,1,5));
    sphereFar->setRadius(5);
    sphereFar->setMaterial(m = new Phong(Vector3(0.0f, 0.0f, 1.0f), Vector3(0.1f), Vector3(0.0f), 10, 1.5));

 //   g_scene->addObject(sphereFar);


	sphere = new Sphere;
    sphere->setCenter(Vector3(0,0.5,0));
    sphere->setRadius(3);
    sphere->setMaterial(m = new Phong(Vector3(0.0f, 1.0f, 0.0f), Vector3(1.f), Vector3(0.0f), 10, 1.5));
    
    g_scene->addObject(sphere);

    Plane *plane = new Plane;
    plane->setMaterial(m = new TexturedPhong(new CheckerBoardTexture(Vector3(1), Vector3(0), 1), Vector3(0)));
    plane->setOrigin(Vector3(0,-1,0));
    g_scene->addObject(plane);

    g_scene->setEnvironment(autumnHDR);
    g_scene->preCalc();
}

void
makeTestScene()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    cout << "Test scene" << endl;
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(4096, 4096);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    //g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setEye(Vector3(0, 0, -8));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, -15));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

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
    addModel("models/testobj.obj", new Lambert(Vector3(0,1,0)), g_scene, Vector3(0), 0);

    g_scene->preCalc();
}

int
main(int argc, char*argv[])
{
//    srand(time(0));
    //Initialize FreeImage
    FreeImage_Initialise();
#ifdef __SSE4_1__
    cout << "Using SSE" << endl;
#endif
#ifdef OPENMP
    cout << "Using OpenMP with up to " << omp_get_max_threads() << " threads." << endl;
#endif
    //makeTestScene();
    //makeSpiralScene();

    //Assignment 1 scenes
    //makeScene1();
    //makeScene2();
    //makeBUNNIZ();
    //A1makeTeapotScene();
    //A1makeSphereScene();
    //A1makeBunnyScene();

    //Assignment 2 scenes
    //makeBunny1Scene();
    //makeBunny20Scene();
    //makeSponzaScene();
    //makeCornellScene();
    //makeTeapotScene();

    //Assignment 3 scenes
    //makeTestSphereScene(); 
    makeTestPetalScene();
	//makeTestSphereTextureScene();


    MiroWindow miro(&argc, argv);
#ifndef NO_GFX
    miro.mainLoop();
#else
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();
#endif


	FreeImage_DeInitialise();

    return 0;
}

