#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::ivec2;
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

    for(int v=0; v<3; ++v){
      ivec2 projPos;
      VertexShader( vertices[v], projPos );
      vec3 color(1,1,1);
      PutPixelSDL( screen, projPos.x, projPos.y, color );
    }
  }

  // vec3 colour(1.0,0.0,0.0);
  // for(int i=0; i<1000; i++){
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  // }
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

void VertexShader( const vec4& v, ivec2& p ){
  vec4 n = v - cameraPos;
  p.x = (int)((focalLength * (n.x/n.z)) + (SCREEN_WIDTH * 0.5));
  p.y = (int)((focalLength * (n.y/n.z)) + (SCREEN_HEIGHT * 0.5));
}

//Rotatin about the origin
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
