#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "Phong.h"
#include <Perlin.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <FreeImage.h>
#include <Worley.h>

//Generate turbulent perlin noise at (x,y,z).
//The returned number is between -1 and 1.
//frequencyIncrease indicates how quickly the frequency increases over multiple turbulence iterations.
//amplitudeFalloff indicates how quickly the amplitude of the higher frequency noise values decrease.
//iterations determines the number of iterations to use for turbulence.
inline float generateNoise(float x, float y, float z, float initialFrequency, float frequencyIncrease, float amplitudeFalloff, int iterations)
{
    float amplitude = 1;
    float frequency = initialFrequency;
    float value = 0;
    float max_val = 0;
    
    for (int i = 0; i < iterations; i++)
    {
        value += amplitude * PerlinNoise::noise(x*frequency, y*frequency, z*frequency);
        max_val += amplitude;
        
        frequency *= frequencyIncrease;
        amplitude *= amplitudeFalloff;
    }
  
    return value/max_val;
}

//Used for the grid in cellular textures
typedef struct gridcell_s
{
    std::vector<tex_coord2d_t> points;
    void addPoint(tex_coord2d_t point) { points.push_back(point); }
} gridcell_t;

class Grid
{
public:
    Grid(int gridWidth, int gridHeight);
    void addPoint(tex_coord2d_t point);

    const std::vector<tex_coord2d_t> &getPoints(int cellI, int cellJ) { return m_grid[cellI][cellJ].points; }
    int getWidth() { return m_gridWidth; }
    int getHeight() { return m_gridHeight; }
private:
    gridcell_t ** m_grid;
    int m_gridWidth, m_gridHeight;
};

class Texture
{
public:
    virtual float bumpHeight2D(const tex_coord2d_t & coords) const { return 0; }
    virtual float bumpHeight3D(const tex_coord3d_t & coords) const { return 0; }
	virtual LookupCoordinates GetLookupCoordinates() const = 0;
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) { return Vector3(0,0,0); } // Look up the color value for a specified position
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) { return Vector3(0,0,0); } // For 3D textures
};

class Texture2D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UV; }
};

class Texture3D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UVW; }
};

class CellularTexture2D : public Texture2D
{
public:
    CellularTexture2D(int points, int gridWidth, int gridHeight);
    float * getClosestDistances(const tex_coord2d_t &point, int n); //Get the n closest distances from a given point
    
    //Populate the grid with grid points. Override to control distribution of points.
    virtual void populateGrid(int points); 
    //Typically, this function combines the values gotten from getClosestDistances() and maps it to a color value.
    virtual Vector3 lookup2D(const tex_coord2d_t & coords); 
    
protected:
    std::vector<tex_coord2d_t> m_points;
    
    Grid m_grid;
};

class StoneTexture : public Texture2D
{
protected:
    float m_scale;
public:
    StoneTexture(float scale=1) { m_scale = scale; }
    //StoneTexture(int points, int gridWidth, int gridHeight) : CellularTexture2D(points, gridWidth, gridHeight) {}
    virtual float bumpHeight2D(const tex_coord2d_t & coords) const;
    virtual Vector3 lookup2D(const tex_coord2d_t & coords);
};


//Texture loaded from a file
class LoadedTexture : public Texture2D
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();

    virtual Vector3 lookup2D(const tex_coord2d_t & coords);
    
protected:
    Vector3 getPixel(int x, int y);     //Get the tonemapped pixel
    FIRGBF getFloatPixel(int x, int y); //Get a pixel from the HDR
    float tonemapValue(float val);

	FIBITMAP* m_bitmap; 
	FIBITMAP* m_toneMapped;
	
	float* m_toneMappedValues;
	
	float m_maxIntensity;
};

//Shading model that also does textures
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
			            const Vector3 & ka = Vector3(0),
			            const Vector3 & ks = Vector3(0),
			            const Vector3 & kt = Vector3(0),
			            const float shinyness = 1.0f,
			            const float refractIndex = 1);
	virtual LookupCoordinates GetLookupCoordinates() const { return m_texture->GetLookupCoordinates(); }

    virtual Vector3 diffuse2D(const tex_coord2d_t & texture_coords) const;
    virtual Vector3 diffuse3D(const tex_coord3d_t & texture_coords) const;
    virtual float   bumpHeight2D(const tex_coord2d_t & texture_coords) const;
    virtual float   bumpHeight3D(const tex_coord3d_t & texture_coords) const;

protected:
    Texture * m_texture;
};

#endif
