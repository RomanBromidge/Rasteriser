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


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles);
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

    // for(int v=0; v<3; ++v){
    //   ivec2 projPos;
    //   VertexShader( vertices[v], projPos );
    //   vec3 color(1,1,1);
    //   PutPixelSDL( screen, projPos.x, projPos.y, color );
    // }

    DrawPolygonEdges(vertices, screen);
  }
}

/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code)
	      {
	      case SDLK_UP:
		/* Move camera forward */
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
		break;
	      case SDLK_LEFT:
		/* Move camera left */
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
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
  vec2 step = vec2(b-a) / float(max(N-1,1));
  vec2 current( a );
  for( int i=0; i<N; ++i ) {
    result[i] = current;
    current += step;
  }
}

//Draw a line between 2 points on the screen
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color ) {
  //Find the absolute difference
  ivec2 delta = glm::abs( a - b );
  int pixels = glm::max( delta.x, delta.y ) + 1;
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
  for( int i=0; i<V; ++i )
  {
    VertexShader( vertices[i], projectedVertices[i] );
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for( int i=0; i<V; ++i )
  {
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
