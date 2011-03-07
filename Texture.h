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
#include <cmath>
#include "Utility.h"

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

    const std::vector<tex_coord2d_t> &getPoints(int cellI, int cellJ) const { return m_grid[cellI][cellJ].points; }
    int getWidth() const { return m_gridWidth; }
    int getHeight() const { return m_gridHeight; }
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
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const { return Vector3(0,0,0); } // Look up the color value for a specified position
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) const { return Vector3(0,0,0); } // For 3D textures
    virtual Vector3 lowresLookup2D(const tex_coord2d_t & coords) const { return lookup2D(coords); }
    virtual Vector3 lowresLookup3D(const tex_coord2d_t & coords) const { return lookup2D(coords); }
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
    float * getClosestDistances(const tex_coord2d_t &point, int n) const; //Get the n closest distances from a given point
    
    //Populate the grid with grid points. Override to control distribution of points.
    virtual void populateGrid(int points); 
    //Typically, this function combines the values gotten from getClosestDistances() and maps it to a color value.
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const; 
    
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
    
    virtual float bumpHeight2D(const tex_coord2d_t & coords) const;
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const;
};

class CheckerBoardTexture : public Texture2D
{
protected:
    float m_scale;
    Vector3 m_color1, m_color2;
public:
    CheckerBoardTexture(Vector3 color1 = Vector3(1), Vector3 color2 = Vector3(0), float scale=1) 
    { 
        m_color1 = color1; 
        m_color2 = color2; 
        m_scale = scale; 
    }

    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const
    {
        float u = std::abs(m_scale*(coords.u));
        float v = std::abs(m_scale*(coords.v));
        if (coords.u < 0) u += m_scale;
        if (coords.v < 0) v += m_scale;
        return ((int)((int)u+(int)v) % 2 == 0 ? m_color1 : m_color2); 
    }
};

class CloudTexture : public Texture2D
{
public:
    CloudTexture(float scale = 3.0f, float cloudSize = 0.3f, float density = -0.1f, float sharpness = 20.0f, float ambient = 0.3f, 
                 float shadowThreshold = 0.3, float shadowMagnitude = 0.5, float shadowSharpness = 0.3)
    { 
        m_scale = scale;
        
        m_cloudSize = cloudSize;
        m_density = density;
        m_sharpness = sharpness;
        m_ambient = ambient;
        m_shadowThreshold = shadowThreshold;
        m_shadowMagnitude = shadowMagnitude;
        m_shadowSharpness = shadowSharpness;
    }

    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const
    {
        float u = m_scale * coords.u, 
              v = m_scale * coords.v;
        
        float val = generateNoise(u, v, 0, 1.0f/m_cloudSize, 2, 0.5, 15);
        
        float cloudResult = std::min(1.0f, m_ambient+sigmoid(m_sharpness*(val + m_density)));
        float shadowResult = m_shadowMagnitude*sigmoid(m_shadowSharpness*m_sharpness*(val - m_shadowThreshold));
        
        return Vector3(cloudResult, cloudResult, 1) - 
               Vector3(shadowResult);
    }
    
private:
    float m_scale;
    float m_cloudSize, m_density, m_sharpness, m_ambient, m_shadowThreshold, m_shadowMagnitude, m_shadowSharpness;
};

class PetalTexture : public Texture3D
{
protected:
    float m_scale;
	float m_radius;
	Vector3 m_pivot;
public:
    PetalTexture(const Vector3 & Pivot, float Radius=1, float scale=1) { m_scale = scale; m_pivot = Pivot; m_radius = Radius; }
    virtual float bumpHeight3D(const tex_coord3d_t & coords) const;
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) const;
};

class FlowerCenterTexture : public Texture3D
{
protected:
    float m_scale;
	float m_radius;
	Vector3 m_pivot;
public:
    FlowerCenterTexture(const Vector3 & Pivot, float Radius=1, float scale=1) { m_scale = scale; m_pivot = Pivot; m_radius = Radius; }
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) const
    {
        float x = coords.u, y = coords.v, z = coords.w;
        Vector3 point(x,y,z);
        float dist = (point-m_pivot).length();
        float fraction = std::max(std::min((float)pow(dist / m_radius, 30.0f), 1.0f), 0.0f);
        
        float maxRed = 0.92f, maxGreen = 0.71f,
              minRed = 0.31f, minGreen = 0.18f;

        float red = std::min((1.0f-fraction)*minRed+fraction*maxRed, 1.0f);
        float green = std::min((1.0f-fraction)*minGreen+fraction*maxGreen, 1.0f);
        
       
        return Vector3(red, green, 0.1f);
    };
};

//Texture loaded from a file
class LoadedTexture : public Texture2D
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();

    Vector3 lookup(const tex_coord2d_t & coords, bool lowres) const;
    Vector3 lookup2D(const tex_coord2d_t & coords) const { return lookup(coords, false); }
    Vector3 lowresLookup2D(const tex_coord2d_t & coords) const { return lookup(coords, true); }
    
protected:
    static Vector3 getPixel(FIBITMAP* bm, int x, int y);     //Get the tonemapped pixel
    static void setPixel(FIBITMAP* bm, Vector3& value, int x, int y);
    float tonemapValue(float val) const;

	FIBITMAP* m_bitmap; 
    FIBITMAP* m_lowres;
    static const int LOWRES_WIDTH = 24;	
	float m_maxIntensity;
};

//Shading model that also does textures
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
            			const Vector3 & specularColor = Vector3(0),    //Reflectivity for each color. 1 is fully reflective, 0 is fully non-reflective.
			            const Vector3 & transparentColor = Vector3(0), //Transparency for each color. 1 is fully transparent (refracting according to refractIndex), 0 is fully opaque.
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
