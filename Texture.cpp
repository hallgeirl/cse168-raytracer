#include "Texture.h"
#include <FreeImage.h>
#include <algorithm>
#include <memory.h>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>

using namespace std;


/*******************************
class LoadedTexture
Texture loaded from a image file.
********************************/
LoadedTexture::LoadedTexture(std::string filename)
{
	FIBITMAP* tmp = FreeImage_Load(FIF_HDR, filename.c_str());

	m_bitmap = FreeImage_TmoDrago03(tmp); //Do tonemapping to get pixel values in a sensible range
}
LoadedTexture::~LoadedTexture()
{
	free(m_bitmap->data);
	free(m_bitmap);
}

//Aux. function to retrieve pixel values from a freeimage bitmap and convert them to vectors.
Vector3 getPixel(FIBITMAP* img, int x, int y)
{
	RGBQUAD color;
	FreeImage_GetPixelColor(img, x, y, &color);
	return Vector3((float)color.rgbRed/255.0f,(float)color.rgbGreen/255.0f,(float)color.rgbBlue/255.0f);
}

Vector3 LoadedTexture::lookup2D(const tex_coord2d_t & texture_coords)
{
	float u = texture_coords.u, v = texture_coords.v;

	//Image dimensions
	int w = FreeImage_GetWidth(m_bitmap), h = FreeImage_GetHeight(m_bitmap);

	//Do bilinear filtering
	float px_real = (float)w*u, py_real = (float)h*v; //Point in image we really want

	int x1 = (int)px_real, x2 = x1+1;
	x1 %= w; x2 %= w;
	float x1_error = px_real - (float)x1;

	int y1 = (int)py_real, y2 = y1+1;
	y1 %= h; y2 %= h;
	float y1_error = py_real - (float)y1;

	//Get pixel values and average them
	Vector3 f = (getPixel(m_bitmap, x1, y1) * (1-x1_error) + getPixel(m_bitmap, x2, y1) * x1_error) * (1 - y1_error) + (getPixel(m_bitmap, x1, y2) * (1-x1_error) + getPixel(m_bitmap, x2, y2) * x1_error) * y1_error;

	return f;
}

/***********************************
class Grid
Grid for use with cellular textures.
************************************/
Grid::Grid(int gridWidth, int gridHeight) : m_gridWidth(gridWidth), m_gridHeight(gridHeight) 
{
    m_grid = new gridcell_t*[m_gridHeight];
    for (int i = 0; i < gridHeight; i++)
    {
        m_grid[i] = new gridcell_t[m_gridWidth];
    }
}

void Grid::addPoint(tex_coord2d_t point)
{
    if (point.u < 0 || point.u > 1.0f || point.v < 0 || point.v > 1.0)
    {
        ostringstream ss;
        ss << "Error adding point to grid. " << point.u << "," << point.v << "is outside the grid.";
        throw runtime_error(ss.str());
    }
    

    int i = (int)(point.v*(float)m_gridHeight), j = (int)(point.u*(float)m_gridWidth);
    m_grid[i][j].addPoint(point);
}


/************************************************************************************
class CellularTexture2D
Base class for all 2D cellular textures. Has methods for finding closest points, etc.
*************************************************************************************/
CellularTexture2D::CellularTexture2D(int points, int gridWidth, int gridHeight) : m_grid(gridWidth, gridHeight)
{
    populateGrid(points);
}

void CellularTexture2D::populateGrid(int points)
{
    //Generate some random points and put them in the grid
    for (int k = 0; k < points; k++)
    {
        tex_coord2d_t point((float)rand() / (float)RAND_MAX, (float)rand() / (float)RAND_MAX);
        m_grid.addPoint(point);
    }
}



Vector3 CellularTexture2D::lookup2D(const tex_coord2d_t & coords)
{
    float *f = getClosestDistances(coords, 4);
    
    float comb = exp(-(f[1]-f[0]+f[2])*100)*2;
    //float out = exp(-(comb)*100);
    //float out = comb > 0.8 ? 0.8 : 0.2;
    float out = comb;
    
    //add some noise for good measure
    //out += (rand()/(float)RAND_MAX)*0.4-0.2; 
   
    return Vector3(out);    
}

//Find the n closest distances to the given point
float* CellularTexture2D::getClosestDistances(const tex_coord2d_t & coords, int n)
{
    //Make sure coords is between 0 and 1
    float u = coords.u - int(coords.u), 
          v = coords.v - int(coords.v);
    
    float cellWidth  = 1.0f/(float)(m_grid.getWidth()),
          cellHeight = 1.0f/(float)(m_grid.getHeight()); 
    
    int gridWidth = m_grid.getWidth(),
        gridHeight = m_grid.getHeight();
        
    //State of cells (visited, added). Removes the need for a linear search.
    //Bit 1 indicates visited or not, bit 2 indicates if it's added to the search queue.
    map<pair<int,int>, int> state;
    
    //The best n points.
    float *f = new float[n];
    
    for (int i = 0; i < n; i++)
        f[i] = 2;//Initialize the distances to a number higher than the maximum distance (which is sqrt(2))
    
    vector<pair<int,int> > searchQueue;
   
    //Start in current cell.
    pair<int,int> currentPos((int)(v*(float)gridHeight), (int)(u*(float)gridWidth));
    searchQueue.push_back(currentPos);
    state[currentPos] = 0;
    int it = 0;
    
    while (searchQueue.size() > 0)
    {
        it++;
        int _i = searchQueue.back().first, _j = searchQueue.back().second; // Current cell position
        state[searchQueue.back()] |= 1; //Mark as visited 
        searchQueue.pop_back();
        const std::vector<tex_coord2d_t> &points = m_grid.getPoints(_i, _j);
        
        //Check the points in this grid cell
        for (int i = 0; i < points.size(); i++)
        {
            tex_coord2d_t diff(fabs(u - points[i].u), fabs(v - points[i].v));
            if (diff.u > 0.5) diff.u = 1-diff.u;
            if (diff.v > 0.5) diff.v = 1-diff.v;
            float dist = sqrt(diff.u*diff.u + diff.v*diff.v);
            
            //Check if we might get a better value for one of our distances. If not, just continue.
            if (f[n-1] > dist)
            {
                for (int j = 0; j < n; j++)
                {
                    if (dist < f[j])
                    {
                        for (int k = n-1; k > j; k--)
                        {
                            f[k] = f[k-1];
                        }
                        f[j] = dist;
                        break;
                    }
                }
            }
        }

        //Add any neighboring cells that might have a better match
        for (int i = _i-1; i <= _i+1; i++)
        {
            //Check the vertical distance to the grid's borders. If it's more than the distance to the best hit so far,
            //we know that we won't get a better hit in that cell and can safely ignore it.
            float y1 = (float)i*cellHeight, y2 = (float)(i+1)*cellHeight;
            //If it's impossible to get a better 2nd. best distance, we can't get a best distance either so just continue.
            //The check on i is done because of possible rounding errors when the point is near the center of the grid square.
            if (fmin(fabs(v-y1), fabs(v-y2)) > f[n-1] && i != _i) continue;

            int i_proper = i + (i<0?gridHeight:0); // Wrapped position
            i_proper %= gridHeight;
            
            //printf("y1 %f\t y2 %f dist1 %f dist2 %f pointx %f pointy %f\n", y1, y2, fabs(v-y1), fabs(v-y2), );
            for (int j = _j-1; j <= _j+1; j++)
            {
                int j_proper = j + (j<0?gridWidth:0);
                j_proper %= gridWidth;
                
                if (i_proper == _i && j_proper == _j) continue; //Ignore the center cell.
                float x1 = (float)j*cellWidth, x2 = (float)(j+1)*cellWidth;
                
                //Again, we can ignore a cell whose horizontal distance to the point is more than the current best 2nd hit. 
                if (fmin(fabs(u-x1), fabs(u-x2)) > f[n-1] && j != _j) continue;
                
                pair<int, int> neighbor(i_proper, j_proper);
                
                //If we haven't searched cell (i,j) yet, add it to the queue to be checked
                //if (find(searchQueue.begin(), searchQueue.end(), pair<int,int>(i_proper,j_proper)) == searchQueue.end() && find(visited.begin(), visited.end(), pair<int,int>(i_proper,j_proper)) == visited.end())
                if (state.find(neighbor) == state.end() || ((state[neighbor] & 1 == 0) && (state[neighbor] & 2 == 0)))
                {
                    searchQueue.push_back(neighbor);
                    state[neighbor] = 2;
                }
            }
        }
    }

    return f;
}

/**************************************************************
class TexturedPhong
Material that looks up diffuse colors from a texture (2D or 3D)
***************************************************************/
TexturedPhong::TexturedPhong(Texture * texture, const Vector3 & ka, const Vector3 & ks, const float shinyness, const float reflect, const float refract, const float refractIndex)
 : Phong(Vector3(0,0,0), ka, ks, shinyness, reflect, refract, refractIndex), m_texture(texture)
{

}

Vector3 TexturedPhong::diffuse2D(const tex_coord2d_t & texture_coords) const
{
    return m_texture->lookup2D(texture_coords);
}

Vector3 TexturedPhong::diffuse3D(const tex_coord3d_t & texture_coords) const
{
    return m_texture->lookup3D(texture_coords);
}
