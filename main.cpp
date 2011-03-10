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

#include "DirectionalAreaLight.h"
using namespace std;

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
	//LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/forrest_salzburg02_big.hdr"));
    cout << "Test scene" << endl;
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
    g_scene->setEnvironment(autumnHDR);
    //g_scene->setEnvironment(new CheckerBoardTexture(Vector3(1), Vector3(0), 100));

    g_image->resize(2048, 2048);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    //g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setEye(Vector3(0, 6, -10));
    g_camera->setLookAt(Vector3(0, 3, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(60);


    // create and place a point light source
   /* PointLight * light = new PointLight;
    light->setPosition(Vector3(0, 10, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(3000);
    g_scene->addLight(light);*/

    // create and place a point light source
    DirectionalAreaLight * light = new DirectionalAreaLight(2);
    light->setPosition(Vector3(0, 10, 0));
    Vector3 lightDir = -light->position();
    lightDir.normalize();
    light->setNormal(lightDir);
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1);
    g_scene->addLight(light);

    Sphere *sp = new Sphere;
    sp->setCenter(Vector3(0,5,0));
    sp->setRadius(1.5);
    sp->setMaterial(new Phong(Vector3(0), Vector3(0), Vector3(1), infinity, 1.5));
    g_scene->addObject(sp);
    
    
//    addModel("models/teapot.obj", new Phong(Vector3(1)), g_scene, Vector3(0.f));
//    addModel("models/teapot.obj", new Phong(Vector3(1)), g_scene, Vector3(0.f));

    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  30));
    floor->setV2(Vector3( 30, 0, -30));
    floor->setV3(Vector3(-30, 0, -30));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(new TexturedPhong(new CheckerBoardTexture(Vector3(1), Vector3(0), 10))); 
    g_scene->addObject(t);

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
//    makeTestScene();
    
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
    //makeTestTextureScene();

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

