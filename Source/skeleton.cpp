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

//Camera variables
float focalLength = SCREEN_HEIGHT;
vec4 cameraPos( 0, 0, -3.001,1 );
mat4 cameraRotMatrix;

//Light variables
vec4 lightPos(0,-0.5,-0.7, 1);
vec3 lightPower = 1.1f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

vector<Triangle> triangles;
vec3 currentColor;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader( screen* screen, const Pixel& p );
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void Rotate(mat3 rotation);
mat3 RotMatrixX(float angle);
mat3 RotMatrixY(float angle);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices, screen* screen);
vec4 NormaliseVec(vec4 vec);
double Find3Distance(vec3 vectorOne, vec3 vectorTwo);

float **malloc2dArray(float dim1, float dim2)
{
	float i, j;

	float **array = (float **)malloc(dim1 * sizeof(float *));

	for (int i = 0; i < dim1; i++) {

		array[i] = (float *)malloc(dim2 * sizeof(float));
		for (int j = 0; j < dim2; j++) {
			array[i][j] = 0;
		}

	}
	return array;

}

//Storing the inverse depth 1/z for each pixel
float **depthBuffer;

int main( int argc, char* argv[] ){

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles);

  depthBuffer = malloc2dArray(SCREEN_WIDTH, SCREEN_HEIGHT);
  /*for (int i=0; i < SCREEN_HEIGHT; i++) {
	  for (int j = 0; j < SCREEN_WIDTH; j++) {
		  depthBuffer[i][j] = 0;
	  }
  }*/

 /* vector<ivec2> vertexPixels(3);
  vertexPixels[0] = ivec2(10, 5);
  vertexPixels[1] = ivec2(5, 10);
  vertexPixels[2] = ivec2(15, 15);
  vector<ivec2> leftPixels;
  vector<ivec2> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  for (int row = 0; row<leftPixels.size(); ++row) {
	  cout << "Start: ("
		  << leftPixels[row].x << ","
		  << leftPixels[row].y << "). "
		  << "End: ("
		  << rightPixels[row].x << ","
		  << rightPixels[row].y << "). " << endl;
  }*/

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

  /* Clear depth buffer */
  //memset(depthBuffer, 0, sizeof(uint32_t)*SCREEN_WIDTH*SCREEN_HEIGHT);
  for (int i = 0; i < SCREEN_WIDTH; i++) {
	  for (int j = 0; j < SCREEN_HEIGHT; j++) {
		  depthBuffer[i][j] = 0;
	  }
  }

  for( uint32_t i=0; i<triangles.size(); ++i ) {
	vector<Vertex> vertices(3);
	vec3 norm3 = triangles[i].normal;
	vec4 norm = vec4(norm3, 1);
	vec3 refl = triangles[i].color;
	vertices[0].position = triangles[i].v0;
	vertices[0].normal = norm;
	vertices[0].reflectance = refl;

	vertices[1].position = triangles[i].v1;
	vertices[1].normal = norm;
	vertices[1].reflectance = refl;

	vertices[2].position = triangles[i].v2;
	vertices[2].normal = norm;
	vertices[2].reflectance = refl;

	currentColor = triangles[i].color;

	//DrawPolygonEdges(vertices, screen);
	DrawPolygon(vertices, screen);
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
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
	int N = result.size();

	Pixel step;
	float stepX = (b.x - a.x) / (float)(fmax(N - 1, 1));
	float stepY = (b.y - a.y) / (float)(fmax(N - 1, 1));

	float d = sqrt((a.x-b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));

	float m = sqrt((stepX * stepX) + (stepY * stepY));
	float lamdaStep = m / d;

	vec3 aLum = a.illumination;
	vec3 bLum = b.illumination;
	vec3 illumination;

	float stepXill = (aLum.x - bLum.x)/ (float)(fmax(N - 1, 1));
	float stepYill = (aLum.y - bLum.y) / (float)(fmax(N - 1, 1));
	float stepZill = (aLum.z - bLum.z) / (float)(fmax(N - 1, 1));


	for (int i = 0; i < N; ++i) {
		result[i].x = round(a.x + i * stepX);
		result[i].y = round(a.y + i * stepY);
		float lamda = i * lamdaStep;
		result[i].zinv = (a.zinv*(1-lamda) + b.zinv*(lamda));

		illumination.x = (a.x + i * stepXill);
		illumination.y = (a.y + i * stepYill);
		illumination.z = (a.y + i * stepZill);
		result[i].illumination = illumination;
	}
}

//Draw a line between 2 points on the screen
//void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color ) {
//
//  //If a and b are the start and end of the line segment then compute the number of pixels to draw
//  //Find the absolute difference
//  ivec2 delta = glm::abs( a - b );
//  int pixels = glm::max( delta.x, delta.y ) + 1;
//
//  //Get the pixel positions of the line by calling the Interpolation function
//  //Create a vector of vectors
//  vector<ivec2> line( pixels );
//  Interpolate( a, b, line );
//
//  for( uint32_t i=0; i<line.size(); ++i ) {
//    PutPixelSDL( screen, line[i].x, line[i].y, color );
//  }
//}

//Draw the edges of each polygon in the scene
//void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen ) {
//
//  int V = vertices.size();
//  // Transform each vertex from 3D world position to 2D image position:
//  vector<ivec2> projectedVertices( V );
//  for( int i=0; i<V; ++i ) {
//    VertexShader( vertices[i], projectedVertices[i] );
//  }
//  // Loop over all vertices and draw the edge from it to the next vertex:
//  for( int i=0; i<V; ++i ) {
//    int j = (i+1)%V; // The next vertex
//    vec3 color( 1, 1, 1 );
//    DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
//  }
//}

void VertexShader(const Vertex& v, Pixel& p) {
	vec4 pos = v.position;
	vec4 n = pos - cameraPos;

	//Calculate z inverse before anything else
	p.zinv = 1 / n.z;
	//Compute projected position
	p.x = (int)((focalLength * (n.x / n.z)) + (SCREEN_WIDTH * 0.5));
	p.y = (int)((focalLength * (n.y / n.z)) + (SCREEN_HEIGHT * 0.5));

	//Compute illumination of vertex
	vec4 norm = v.normal;
	vec4 r = v.position - lightPos;
	norm = NormaliseVec(n);
	r = NormaliseVec(r);
	double d = Find3Distance(lightPos, vec3(pos));
	double dotProduct = (norm.x * r.x) + (norm.y * r.y) + (norm.z * r.z);
	vec3 D (0,0,0);
	if (dotProduct > 0) {
		D.x = dotProduct * lightPower.x / (4 * M_PI * d * d);
		D.y = dotProduct * lightPower.y / (4 * M_PI * d * d);
		D.z = dotProduct * lightPower.z / (4 * M_PI * d * d);
		cout << "positive" << 1 << endl;
	}
	else {
		cout << "negative" << endl;
	}
	p.illumination = (D + indirectLightPowerPerArea);

}

void PixelShader( screen* screen, const Pixel& p ) {
	int x = p.x;
	int y = p.y;
	if( p.zinv > depthBuffer[x][y] ) {
		depthBuffer[x][y] = p.zinv;
		PutPixelSDL( screen, x, y, p.illumination );
	}
}

//Transform the geometry depending on some translation and rotation vector
//void TransformationMatrix(mat4 Transformation, vec3 angles, vec3 translation){
//  float a = angles[0];
//  float b = angles[1];
//  float c = angles[2];
//
//  vec4 row1(cosf(b) * cosf(c), cosf(b)* sinf(c), -sinf(b), 0);
//  vec4 row2((sinf(a) * sinf(b)* cosf(c)) - (cosf(a) * sinf(c)), (sinf(a) * sinf(b) * sinf(c)) + (cosf(a) * cosf(c)), sinf(a) * cosf(b), 0);
//  vec4 row3((cosf(a) * sinf(b)* cosf(c)) + (sinf(a) * sinf(c)), (cosf(a) * sinf(b) * sinf(c)) - (sinf(a) * cosf(c)), cosf(a) * cosf(b), 0);
//  vec4 row4(0,0,0,1);
//
//  mat4 Rotation(row1, row2, row3, row4);
//
//  vec4 zeros(0,0,0,0);
//  mat4 Positive(zeros, zeros, zeros, row4);
//  Positive[0][3] = translation[0];
//  Positive[1][3] = translation[1];
//  Positive[2][3] = translation[2];
//  mat4 Negative(zeros, zeros, zeros, row4);
//  Negative[0][3] = -translation[0];
//  Negative[1][3] = -translation[1];
//  Negative[2][3] = -translation[2];
//
//  Transformation = Positive * Rotation * Negative;
//}

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
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
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
	for (int i = 0; i<ROWS; ++i) {
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i<vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex

		ivec2 delta;
		delta.x = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
		delta.y = glm::abs(vertexPixels[i].y - vertexPixels[j].y);

		int pixels = glm::max(delta.x, delta.y) + 1;

		vector<Pixel> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (Pixel pixel : result) {
			int loc = pixel.y - minYValue;
			if (pixel.x < leftPixels[loc].x) {
				leftPixels[loc].x = pixel.x;
				leftPixels[loc].y = pixel.y;
				leftPixels[loc].zinv = pixel.zinv;
				leftPixels[loc].illumination = pixel.illumination;
			}
			if (pixel.x > rightPixels[loc].x) {
				rightPixels[loc].x = pixel.x;
				rightPixels[loc].y = pixel.y;
				rightPixels[loc].zinv = pixel.zinv;
				rightPixels[loc].illumination = pixel.illumination;
			}
		}
	}
}

void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
	int length = leftPixels.size();
	//Call PutPixelSDL for each pixel between the start and end for each row
	for (int i = 0; i < length; ++i) {

		ivec2 delta;
		delta.x = glm::abs(leftPixels[i].x - rightPixels[i].x);
		delta.y = glm::abs(leftPixels[i].y - rightPixels[i].y);
		int pixels = glm::max(delta.x, delta.y) +1 ;

		//int pixels = glm::abs(leftPixels[i].x - rightPixels[i].x);

		vector<Pixel> results(pixels);

		Interpolate(leftPixels[i], rightPixels[i], results);

		for (Pixel pixel : results) {
			//Check if pixel.x and pixel.y are in the screen
			if (pixel.x < SCREEN_WIDTH && pixel.x >= 0 && pixel.y < SCREEN_HEIGHT && pixel.y >= 0) {
				//Draw pixel based on zinv
				PixelShader(screen, pixel);
			}
		}
	}
}

void DrawPolygon(const vector<Vertex>& vertices, screen* screen) {
	int V = vertices.size();

	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i) {
		VertexShader(vertices[i], vertexPixels[i]);
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(screen, leftPixels, rightPixels);
}

vec4 NormaliseVec(vec4 vec) {
	vec4 normalisedPoint;
	double length = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);

	normalisedPoint.x = vec.x / length;
	normalisedPoint.y = vec.y / length;
	normalisedPoint.z = vec.z / length;
	normalisedPoint.w = vec.w;

	return normalisedPoint;
}

double Find3Distance(vec3 vectorOne, vec3 vectorTwo) {

	double sqrDifferenceX = (vectorOne.x - vectorTwo.x)*(vectorOne.x - vectorTwo.x);
	double sqrDifferenceY = (vectorOne.y - vectorTwo.y)*(vectorOne.y - vectorTwo.y);
	double sqrDifferenceZ = (vectorOne.z - vectorTwo.z)*(vectorOne.z - vectorTwo.z);
	double root = sqrt(sqrDifferenceX + sqrDifferenceY + sqrDifferenceZ);

	return root;
}
