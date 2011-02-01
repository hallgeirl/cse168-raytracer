#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "Phong.h"
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <FreeImage.h>
#include <noise/noisegen.h>

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

class StoneTexture : public CellularTexture2D
{
    
};

//Generates random 3D noise. For testing.
class TestTexture3D : public Texture3D
{
public:
	virtual Vector3 lookup3D(const tex_coord3d_t & coords)
	{
		float x = fabs(coords.u) * 1.0, y = fabs(coords.v) * 1.0, z = fabs(coords.w) * 1.0;
		//std::cout << x << " " << y << " " << z << std::endl;
		//std::cout << Vector3((float)(x%256) / 255.0, (float)(y%256) / 255.0, (float)(z%256) / 255.0) << std::endl;
		//return Vector3(0,0,0);
		float out = 0;
		float ampl = 1;
		float freq = 1;
		for (int i = 0; i < 6; i++)
		{
			out += ampl * (noise::GradientCoherentNoise3D (x*freq, y*freq, z*freq)/2.0 + 0.5);
			ampl *= 0.7;
			freq *= 2;
		}

		return Vector3(out);
	}
};

//Texture loaded from a file
class LoadedTexture : public Texture2D
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();

    virtual Vector3 lookup2D(const tex_coord2d_t & coords);
protected:
	FIBITMAP* m_bitmap;
};

//Shading model that also does textures
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
			            const Vector3 & ka = Vector3(0),
			            const Vector3 & ks = Vector3(1),
			            const float shinyness = 1.0f,
			            const float reflect = 0,
			            const float refract = 0,
			            const float refractIndex = 1);
	virtual LookupCoordinates GetLookupCoordinates() const { return m_texture->GetLookupCoordinates(); }

    virtual Vector3 diffuse2D(const tex_coord2d_t & texture_coords) const;
    virtual Vector3 diffuse3D(const tex_coord3d_t & texture_coords) const;

protected:
    Texture * m_texture;
};

#endif
