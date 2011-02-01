#include "Texture.h"
#include <FreeImage.h>
#include <algorithm>
#include <memory.h>
#include <cstdlib>
#include <cstdio>
#include <map>

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

/*****************
class StoneTexture
Stone texture.
******************/
bool pointComparer (const tex_coord2d_t &p1, const tex_coord2d_t &p2) { return (p1.u < p2.u); }

StoneTexture::StoneTexture(int points)
{
    int _k;
    
    //Generate some points and put them in the grid
    for (int k = 0; k < points; k++)
    {
        tex_coord2d_t point((float)rand() / (float)RAND_MAX, (float)rand() / (float)RAND_MAX);
        int i = (int)(point.v*(float)GRID_HEIGHT), j = (int)(point.u*(float)GRID_WIDTH);
        m_grid[i][j].addPoint(point);
    }
}

Vector3 StoneTexture::lookup2D(const tex_coord2d_t & coords)
{
    //Make sure coords is between 0 and 1
    float u = coords.u - int(coords.u), 
          v = coords.v - int(coords.v);
    
    float cellWidth  = 1.0f/(float)GRID_WIDTH,
          cellHeight = 1.0f/(float)GRID_HEIGHT; 
          
    //State of cells (visited, added). Removes the need for a linear search.
    //Bit 1 indicates visited or not, bit 2 indicates if it's added to the search queue.
    map<pair<int,int>, int> state;
    
    //The two best points.
    float f[2] = {2, 2}; 
    
    vector<pair<int,int> > searchQueue;
   
    //Start in current cell.
    pair<int,int> currentPos((int)(v*(float)GRID_HEIGHT), (int)(u*(float)GRID_WIDTH));
    searchQueue.push_back(currentPos);
    state[currentPos] = 0;
    int it = 0;
    
    while (searchQueue.size() > 0)
    {
        it++;
        int _i = searchQueue.back().first, _j = searchQueue.back().second; // Current cell position
        state[searchQueue.back()] |= 1; //Mark as visited 
        searchQueue.pop_back();
        
        //Check the points in this grid cell
        for (int i = 0; i < m_grid[_i][_j].points.size(); i++)
        {
            tex_coord2d_t diff(fabs(u - m_grid[_i][_j].points[i].u), fabs(v - m_grid[_i][_j].points[i].v));
            if (diff.u > 0.5) diff.u = 1-diff.u;
            if (diff.v > 0.5) diff.v = 1-diff.v;
            float dist = sqrt(diff.u*diff.u + diff.v*diff.v);
            
            if (dist < f[1])
            {
                if (dist < f[0])
                {
                    //Found new best hit.
                    f[1] = f[0]; f[0] = dist;
                }
                else
                {
                    f[1] = dist;
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
            if (fmin(fabs(v-y1), fabs(v-y2)) > f[1] && i != _i) continue;

            int i_proper = i + (i<0?GRID_HEIGHT:0); // Wrapped position
            i_proper %= GRID_HEIGHT;
            
            //printf("y1 %f\t y2 %f dist1 %f dist2 %f pointx %f pointy %f\n", y1, y2, fabs(v-y1), fabs(v-y2), );
            for (int j = _j-1; j <= _j+1; j++)
            {
                int j_proper = j + (j<0?GRID_WIDTH:0);
                j_proper %= GRID_WIDTH;
                
                if (i_proper == _i && j_proper == _j) continue; //Ignore the center cell.
                float x1 = (float)j*cellWidth, x2 = (float)(j+1)*cellWidth;
                
                //Again, we can ignore a cell whose horizontal distance to the point is more than the current best 2nd hit. 
                if (fmin(fabs(u-x1), fabs(u-x2)) > f[1] && j != _j) continue;
                
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
    //printf("%d\n", it);
//    float out = 2.0 * f[0] / (f[1]+f[0]);
    
    float out = exp(-(f[1]-f[0])*100)/2.72;
    
    //if (f[0] > f[1]) printf("FOOOO");
    return Vector3(out);    
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
