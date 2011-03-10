#include "assignment3.h"
#include <math.h>
#include <string>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Sphere.h"
#include "DirectionalAreaLight.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"
#include "Plane.h"

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
	Texture *cloudTexture = new CloudTexture(3.0f,  //Scale (higher is smaller)
	                                        0.1f,  //Cloud size
	                                        0.2f, //Cloud density
	                                        50.0f,  //Sharpness
	                                        0.4f,  //Ambient (higher number means lower contrast)
	                                        0.35f,  //Shadow threshold
	                                        0.5f,  //Shadow magnitude
	                                        0.3f   //Shadow sharpness
	                                        ); 

    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
    g_scene->setEnvironment(new LoadedTexture("gfx/forrest_salzburg02_big.hdr"));
    g_scene->setEnvironmentRotation(PI/3+0.05, PI/8);
    g_scene->setBgColor(Vector3(1));

    float aspect = 1.5;
    int res = 4096;

    g_image->resize(res, int((float)res/aspect));
    
    Vector3 lightPos = Vector3(50,10000,40);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));

//Don't delete    
g_camera->setEye(Vector3(2, 4.4, 16.8));
//Don't delete    
g_camera->setLookAt(Vector3(3, 0.0, 4));

    //g_camera->setLookAt(Vector3(-3, 0.0, 4));

    /*g_camera->setEye(Vector3(2, 4.4, 10));
    g_camera->setLookAt(Vector3(2, 0.0, 4));*/
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(30);

    // create and place a point light source
    DirectionalAreaLight * light = new DirectionalAreaLight(7);
    light->setPosition(lightPos);
    Vector3 lightNormal = -light->position();
    lightNormal.normalize();
    light->setNormal(lightNormal);
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(4);
    g_scene->addLight(light);

	 // create and place a point light source
    PointLight * amlight = new PointLight;
    amlight->setPosition(Vector3(-20, 20, 20));
    amlight->setColor(Vector3(1, 1, 1));
    amlight->setWattage(30000);
    //g_scene->addLight(amlight);

	Material* material = new TexturedPhong(new PetalTexture(Vector3(0.f), 7), Vector3(0.), Vector3(0), infinity, 1.5);
	addFlowerModel("models/Petals2.obj", material, g_scene, Vector3(0.f), 0.0);

	//material = new Phong(Vector3(0,1,0));
	material = new TexturedPhong(new StemTexture(30));
	addFlowerModel("models/Stem.obj", material, g_scene, Vector3(0.f), 0.0);

	material = new Phong(Vector3(0,1,0));
	//material = new TexturedPhong(new LeafTexture(Vector3(0,0,0), Vector3(1,0,0)));
	//material = new TexturedPhong(new LoadedTexture("LeafTexture.jpg"));
	addFlowerModel("models/Leaf.obj", material, g_scene, Vector3(0.f, 0.5f, 0.f), 0.0);

	material = new TexturedPhong(new FlowerCenterTexture(Vector3(-0.1,-0.35,0), 1.1));
	addFlowerModel("models/FlowerCenter.obj", material, g_scene, Vector3(0.f));//Vector3(-0.05f, 0.25, 0.32f), 0);

	Material* water = new Phong(Vector3(1.f), Vector3(0), Vector3(1.0f), infinity, 1.33);
	//Material* water = new Phong(Vector3(1.f), Vector3(0), Vector3(1.0f), infinity, 1);
	addFlowerModel("models/WaterDropsMany.obj", water, g_scene, Vector3(0.f, 0.f, 0));
	//addFlowerModel("models/WaterDrops.obj", water, g_scene, Vector3(0.f, 0.f, 0));
	//addFlowerModel("models/WaterDropsHiRes.obj", water, g_scene, Vector3(0.f));
    
    /*Sphere *sp = new Sphere();
    sp->setCenter(Vector3(-0.1,-0.35,0));
    sp->setRadius(0.5);
    sp->setMaterial(new Phong(Vector3(1)));
    g_scene->addObject(sp);*/
    
   
    Plane *p = new Plane;
    p->setMaterial(new Phong(Vector3(1)));
    p->setNormal(Vector3(0,-1,0));
    p->setOrigin(Vector3(0, 51, 0));
  //  g_scene->addObject(p);
    
    p = new Plane;
    p->setMaterial(new Phong(Vector3(0.5)));
    p->setNormal(Vector3(1,0,0));
    p->setOrigin(Vector3(3, 0, 0));
    //g_scene->addObject(p);
    
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


void 
makeTestTextureScene()
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

	Texture *tex = new StemTexture();
	//Texture *tex = new LeafTexture(Vector3(0,0,0), Vector3(1,0,0));
	Material* material = new TexturedPhong(tex);
	Plane *p = new Plane;
	p->setMaterial(material);
	g_scene->addObject(p);

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
