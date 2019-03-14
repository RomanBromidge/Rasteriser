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
mat4 transformMatrix;
float angleStep = 0.01f;
float transStep = 0.2f;

vec4 lightPos(0, -0.5, -0.7, 1);
vec3 lightPower = 14.1f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);

vector<Triangle> triangles;

vec3 currentColor;
vec4 currentNormal;
vec3 currentReflectance;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader(const vec4& v, Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void TransformationMatrix(mat4 r, vec3 t);
mat4 RotMatrixX(float angle);
mat4 RotMatrixY(float angle);
mat4 RotMatrixZ(float angle);
mat4 Ident();
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void PixelShader(screen* screen, const Pixel& p);
vec3 calculateIllumination(const Pixel& p);
void DrawPolygon(const vector<vec4>& vertices, screen* screen);
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
	transformMatrix = Ident();
	TransformationMatrix(Ident(), vec3(-cameraPos));

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
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

	currentColor = triangles[i].color;

	currentNormal = triangles[i].normal;
	currentReflectance = currentColor;

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
						TransformationMatrix(Ident(), vec3(0, 0, -transStep));
						break;
					case SDLK_DOWN:
						/* Move camera backwards */
						TransformationMatrix(Ident(), vec3(0, 0, transStep));
						break;
					case SDLK_LEFT:
						/* Move camera left */
						TransformationMatrix(Ident(), vec3(transStep, 0, 0));
						break;
					case SDLK_RIGHT:
						/* Move camera right */
						TransformationMatrix(Ident(), vec3(-transStep, 0, 0));
						break;

					//Camera Rotation
					case SDLK_k:
						/*Rotate camera downwards on X axis */
						TransformationMatrix(RotMatrixX(-angleStep), vec3(0, 0, 0));
						break;
					case SDLK_i:
						/*Rotate camera upwards on X axis */
						TransformationMatrix(RotMatrixX(angleStep), vec3(0, 0, 0));
						break;
					case SDLK_l:
						/*Rotate camera left on Y axis */
						TransformationMatrix(RotMatrixY(angleStep), vec3(0, 0, 0));
						break;
					case SDLK_j:
						/*Rotate camera right on Y axis */
						TransformationMatrix(RotMatrixY(-angleStep), vec3(0, 0, 0));
						break;

						//Light Rotation
					case SDLK_d:
						lightPos.x += transStep;
						break;
					case SDLK_a:
						lightPos.x -= transStep;
						break;
					case SDLK_s:
						lightPos.y += transStep;
						break;
					case SDLK_w:
						lightPos.y -= transStep;
						break;

					case SDLK_ESCAPE:
						/* Move camera quit */
						return false;
				}
			}
    }
  return true;
}

//Return 4x4 identity matrix
mat4 Ident() {
	return mat4(1.0f);
}
mat4 RotMatrixX(float angle) {
	vec4 c0(1, 0, 0, 0);
	vec4 c1(0, cosf(angle), sinf(angle), 0);
	vec4 c2(0, -sinf(angle), cosf(angle), 0);
	vec4 c3(0, 0, 0, 1);
	return mat4(c0, c1, c2, c3);
}
mat4 RotMatrixY(float angle) {
	vec4 c1(cosf(angle), 0, -sinf(angle), 0);
	vec4 c2(0, 1, 0, 0);
	vec4 c3(sinf(angle), 0, cosf(angle), 0);
	vec4 c4(0, 0, 0, 1);
	return mat4(c1, c2, c3, c4);
}
mat4 RotMatrixZ(float angle) {
	vec4 c1(cosf(angle), sinf(angle), 0, 0);
	vec4 c2(-sinf(angle), cosf(angle), 0, 0);
	vec4 c3(0, 0, 1, 0);
	vec4 c4(0, 0, 0, 1);
	return mat4(c1, c2, c3, c4);
}

//Transform geometry by some amount
void TransformationMatrix(mat4 r, vec3 t) {
	vec4 cT = transformMatrix[3];
	mat4 newTrans = r * transformMatrix;
	vec4 c3(cT.x + t.x, cT.y + t.y, cT.z + t.z, 1);
	newTrans[3] = c3;
	transformMatrix = newTrans ;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
	int N = result.size();

	float stepX = (b.x - a.x) / (float)(fmax(N - 1, 1));
	float stepY = (b.y - a.y) / (float)(fmax(N - 1, 1));

	float d = sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));

	float m = sqrt((stepX * stepX) + (stepY * stepY));
	float lamdaStep = m / d;

	vec3 step3d = (b.pos3d*b.zinv - a.pos3d*a.zinv) / (float)(fmax(N - 1, 1));

	for (int i = 0; i < N; ++i) {
		result[i].x = round(a.x + i * stepX);
		result[i].y = round(a.y + i * stepY);

		float lamda = i * lamdaStep;
		result[i].zinv = (a.zinv*(1-lamda) + b.zinv*(lamda));

		result[i].pos3d = vec4(((vec3)a.pos3d*a.zinv + (float)i * step3d),1) / result[i].zinv;
	}
}

void VertexShader(const vec4& v, Pixel& p) {
	//vec4 n = v - cameraPos;
	vec4 n = transformMatrix * v;

	//Calculate z inverse before anything else:
	p.zinv = 1 / n.z;

	p.x = (int)((focalLength * (n.x / n.z)) + (SCREEN_WIDTH * 0.5));
	p.y = (int)((focalLength * (n.y / n.z)) + (SCREEN_HEIGHT * 0.5));

	vec4 pixelPos(v.x, v.y, v.z, 1);
	p.pos3d = pixelPos;
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
				leftPixels[loc].pos3d = pixel.pos3d;
			}
			if (pixel.x > rightPixels[loc].x) {
				rightPixels[loc].x = pixel.x;
				rightPixels[loc].y = pixel.y;
				rightPixels[loc].zinv = pixel.zinv;
				rightPixels[loc].pos3d = pixel.pos3d;
			}
		}
	}
}

void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
	int length = leftPixels.size();
	//Call PutPixelSDL for each pixel between the start and end for each row
	for (int i = 0; i < length; ++i) {

		int pixels = glm::abs(leftPixels[i].x - rightPixels[i].x) + 1;

		vector<Pixel> results(pixels);

		Interpolate(leftPixels[i], rightPixels[i], results);

		for (Pixel pixel : results) {
			//Check if pixel.x and pixel.y are in the screen
			if (pixel.x < SCREEN_WIDTH && pixel.x >= 0 && pixel.y < SCREEN_HEIGHT && pixel.y >= 0) {
				PixelShader(screen, pixel);
			}
		}
	}
}

void PixelShader(screen* screen, const Pixel& p) {
	int x = p.x;
	int y = p.y;

	vec3 illumination = calculateIllumination(p);

	//Before drawing a new pixel check if it is in front of the one that is currently drawn there by comparing the value of the depthBuffer at that pixel
	if (p.zinv > depthBuffer[p.x][p.y]) {
		depthBuffer[p.x][p.y] = p.zinv;
		PutPixelSDL(screen, x, y, currentColor * illumination);
	}
}

vec3 calculateIllumination(const Pixel& p) {
	//Compute illumination of vertex
	vec4 norm = currentNormal;
	vec4 r = lightPos - p.pos3d;
	norm = NormaliseVec(norm);
	r = NormaliseVec(r);
	double d = Find3Distance(lightPos, vec3(p.pos3d));
	double dotProduct = (norm.x * r.x) + (norm.y * r.y) + (norm.z * r.z);
	vec3 D(0, 0, 0);
	if (dotProduct > 0) {
		D.x = dotProduct * lightPower.x / (4 * M_PI * d * d);
		D.y = dotProduct * lightPower.y / (4 * M_PI * d * d);
		D.z = dotProduct * lightPower.z / (4 * M_PI * d * d);
	}
	return  (D + indirectLightPowerPerArea);
}

void DrawPolygon(const vector<vec4>& vertices, screen* screen) {
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
