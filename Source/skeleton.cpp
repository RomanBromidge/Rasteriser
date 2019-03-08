#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

float focalLength = SCREEN_HEIGHT;
vec4 cameraPos( 0, 0, -3.001,1 );
mat4 cameraRotMatrix;

vector<Triangle> triangles;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
void TransformationMatrix(mat4 Tranformation, vec3 angles, vec3 translation);
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );
void Rotate(mat3 rotation);
mat3 RotMatrixX(float angle);
mat3 RotMatrixY(float angle);
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels);


int main( int argc, char* argv[] ){

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles);
  
  vector<ivec2> vertexPixels(3);
  vertexPixels[0] = ivec2(10, 5);
  vertexPixels[1] = ivec2(5, 10);
  vertexPixels[2] = ivec2(15, 15);
  vector<ivec2> leftPixels;
  vector<ivec2> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  for (int row = 0; row<leftPixels.size(); ++row)
  {
	  cout << "Start: ("
		  << leftPixels[row].x << ","
		  << leftPixels[row].y << "). "
		  << "End: ("
		  << rightPixels[row].x << ","
		  << rightPixels[row].y << "). " << endl;
  }

  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for( uint32_t i=0; i<triangles.size(); ++i ) {
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    DrawPolygonEdges(vertices, screen);
	//DrawPolygon(screen, vertices);
  }
}

/*Place updates of parameters here*/
bool Update(){

	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;

	float angle = 0.01f;
	float step = 0.2f;

	SDL_Event e;
	while(SDL_PollEvent(&e)){
		if (e.type == SDL_QUIT){
			return false;
		}
		else
			if (e.type == SDL_KEYDOWN){
				int key_code = e.key.keysym.sym;
				switch(key_code){

					//Camera Movement
					case SDLK_UP:
						/* Move camera forward */
						cameraPos.z += step;
						break;
					case SDLK_DOWN:
						/* Move camera backwards */
						cameraPos.z -= step;
						break;
					case SDLK_LEFT:
						/* Move camera left */
						cameraPos.x -= step;
						break;
					case SDLK_RIGHT:
						/* Move camera right */
						cameraPos.x += step;
						break;

					//Camera Rotation
					case SDLK_k:
						/*Rotate camera downwards on X axis */
						Rotate(RotMatrixX(-angle));
						break;
					case SDLK_i:
						/*Rotate camera upwards on X axis */
						Rotate(RotMatrixX(angle));
						break;
					case SDLK_l:
						/*Rotate camera left on Y axis */
						Rotate(RotMatrixY(angle));
						break;
					case SDLK_j:
						/*Rotate camera right on Y axis */
						Rotate(RotMatrixY(-angle));
						break;

					case SDLK_ESCAPE:
						/* Move camera quit */
						return false;
				}
			}
    }
  return true;
}

//Interpolate a 2d vector between two 2d points
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ) {
  int N = result.size();
  vec2 step = vec2(b-a) / float(fmax(N-1,1));
  vec2 current( a );
  for( int i=0; i<N; ++i ) {
    result[i] = current;
    current += step;
  }
}

//Draw a line between 2 points on the screen
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color ) {

  //If a and b are the start and end of the line segment then compute the number of pixels to draw
  //Find the absolute difference
  ivec2 delta = glm::abs( a - b );
  int pixels = glm::max( delta.x, delta.y ) + 1;

  //Get the pixel positions of the line by calling the Interpolation function
  //Create a vector of vectors
  vector<ivec2> line( pixels );
  Interpolate( a, b, line );

  for( uint32_t i=0; i<line.size(); ++i ) {
    PutPixelSDL( screen, line[i].x, line[i].y, color );
  }
}

//Draw the edges of each polygon in the scene
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen ) {

  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<ivec2> projectedVertices( V );
  for( int i=0; i<V; ++i ) {
    VertexShader( vertices[i], projectedVertices[i] );
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for( int i=0; i<V; ++i ) {
    int j = (i+1)%V; // The next vertex
    vec3 color( 1, 1, 1 );
    DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
  }
}

//Determine the position of the projected vertices
void VertexShader( const vec4& v, ivec2& p ){
  vec4 n = v - cameraPos;
  p.x = (int)((focalLength * (n.x/n.z)) + (SCREEN_WIDTH * 0.5));
  p.y = (int)((focalLength * (n.y/n.z)) + (SCREEN_HEIGHT * 0.5));
}

//Transform the geometry depending on some translation and rotation vector
void TransformationMatrix(mat4 Transformation, vec3 angles, vec3 translation){
  float a = angles[0];
  float b = angles[1];
  float c = angles[2];

  vec4 row1(cosf(b) * cosf(c), cosf(b)* sinf(c), -sinf(b), 0);
  vec4 row2((sinf(a) * sinf(b)* cosf(c)) - (cosf(a) * sinf(c)), (sinf(a) * sinf(b) * sinf(c)) + (cosf(a) * cosf(c)), sinf(a) * cosf(b), 0);
  vec4 row3((cosf(a) * sinf(b)* cosf(c)) + (sinf(a) * sinf(c)), (cosf(a) * sinf(b) * sinf(c)) - (sinf(a) * cosf(c)), cosf(a) * cosf(b), 0);
  vec4 row4(0,0,0,1);

  mat4 Rotation(row1, row2, row3, row4);

  vec4 zeros(0,0,0,0);
  mat4 Positive(zeros, zeros, zeros, row4);
  Positive[0][3] = translation[0];
  Positive[1][3] = translation[1];
  Positive[2][3] = translation[2];
  mat4 Negative(zeros, zeros, zeros, row4);
  Negative[0][3] = -translation[0];
  Negative[1][3] = -translation[1];
  Negative[2][3] = -translation[2];

  Transformation = Positive * Rotation * Negative;
}

mat3 RotMatrixX(float angle) {
	vec3 x0(1, 0, 0);
	vec3 x1(0, cosf(angle), sinf(angle));
	vec3 x2(0, -sinf(angle), cosf(angle));

	return mat3(x0, x1, x2);
}

mat3 RotMatrixY(float angle) {
	vec3 y1(cosf(angle), 0, -sinf(angle));
	vec3 y2(0, 1, 0);
	vec3 y3(sinf(angle), 0, cosf(angle));

	return mat3(y1, y2, y3);
}

void Rotate(mat3 rotation) {
	vec3 loc = cameraPos * rotation;
	cameraPos = vec4(loc, 1);

	/*mat3 cameraRotExtract(cameraRotMatrix[0][0], cameraRotMatrix[0][1], cameraRotMatrix[0][2],
		cameraRotMatrix[1][0], cameraRotMatrix[1][1], cameraRotMatrix[1][2],
		cameraRotMatrix[2][0], cameraRotMatrix[2][1], cameraRotMatrix[2][2]);

	cameraRotExtract = rotation * cameraRotExtract;

	cameraRotMatrix[0][0] = cameraRotExtract[0][0];
	cameraRotMatrix[0][1] = cameraRotExtract[0][1];
	cameraRotMatrix[0][2] = cameraRotExtract[0][2];
	cameraRotMatrix[1][0] = cameraRotExtract[1][0];
	cameraRotMatrix[1][1] = cameraRotExtract[1][1];
	cameraRotMatrix[1][2] = cameraRotExtract[1][2];
	cameraRotMatrix[2][0] = cameraRotExtract[2][0];
	cameraRotMatrix[2][1] = cameraRotExtract[2][1];
	cameraRotMatrix[2][2] = cameraRotExtract[2][2];*/
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels){

	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.

	int maxYValue = vertexPixels[0].y;
	if (maxYValue < vertexPixels[1].y) {
		maxYValue = vertexPixels[1].y;
	}
	if (maxYValue < vertexPixels[2].y) {
		maxYValue = vertexPixels[2].y;
	}

	int minYValue = vertexPixels[0].y;
	if (minYValue > vertexPixels[1].y) {
		minYValue = vertexPixels[1].y;
	}
	if (minYValue > vertexPixels[2].y) {
		minYValue = vertexPixels[2].y;
	}

	int ROWS = maxYValue - minYValue + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.

	leftPixels.resize(ROWS);
	rightPixels.resize(ROWS);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.

	for (int i = 0; i<ROWS; ++i){
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.

	for (int i = 0; i<vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex
		ivec2 delta = glm::abs(vertexPixels[i] - vertexPixels[j]);

		int pixels = glm::max(delta.x, delta.y) + 1;
		vector<ivec2> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (ivec2 pixel : result) {
			int loc = pixel.y - minYValue;
			if (pixel.x < leftPixels[loc].x) {
				leftPixels[loc].x = pixel.x;
				leftPixels[loc].y = pixel.y;
			}
			if (pixel.x > rightPixels[loc].x) {
				rightPixels[loc].x = pixel.x;
				rightPixels[loc].y = pixel.y;
			}
		}
	}
}

