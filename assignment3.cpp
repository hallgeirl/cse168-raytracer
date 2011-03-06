#include "assignment3.h"
#include <math.h>
#include <string>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Sphere.h"
#include "PointLight.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"

using namespace std;

void addFlowerModel(const char* filename, Material *mat, Scene* scene, Vector3 position, float rotY, Vector3 scale)
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

void
makeTestPetalScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(2, 4, 0));
    g_camera->setLookAt(Vector3(0, 0, 6));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(10, 10, 10));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(400);
    g_scene->addLight(light);

	 // create and place a point light source
    PointLight * amlight = new PointLight;
    amlight->setPosition(Vector3(-10, 10, 10));
    amlight->setColor(Vector3(1, 1, 1));
    amlight->setWattage(300);
    g_scene->addLight(amlight);

	//Material* material = new Lambert(Vector3(1.0f, 0.3f, 0.3f));
	Material* material = new TexturedPhong(new PetalTexture(Vector3(0.f), 10), Vector3(0), Vector3(0), 5);
	addFlowerModel("models/petal.obj", material, g_scene, Vector3(0.f));

	Material* water = new Phong(Vector3(1.f), Vector3(0), Vector3(1.0f), 5, 1.33);
	addFlowerModel("models/WaterDrops.obj", water, g_scene, Vector3(0.f));
    
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

void 
makeTestSphereTextureScene()
{
	g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-10, 4, 0));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(10, 10, 10));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(5000);
    g_scene->addLight(light);

	 // create and place a point light source
    PointLight * amlight = new PointLight;
    amlight->setPosition(Vector3(-10, 10, 10));
    amlight->setColor(Vector3(1, 1, 1));
    amlight->setWattage(5000);
    g_scene->addLight(amlight);

	LoadedTexture *earth = new LoadedTexture(std::string("gfx/earth.jpg"));
//	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"), FIF_HDR);
	Material* material = new TexturedPhong(earth, Vector3(0), Vector3(0), 5);
	addFlowerModel("models/TexturedSphere.obj", material, g_scene, Vector3(0.f));

	// create the floor triangle
    /*TriangleMesh * floor = new TriangleMesh;
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
    g_scene->addObject(t);*/
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}


