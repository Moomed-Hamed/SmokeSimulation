#pragma comment(lib, "opengl32")
#pragma comment(lib, "external/GLEW/glew32s")
#pragma comment(lib, "external/GLFW/glfw3")
#pragma comment(lib, "glu32.lib")

#define GLEW_STATIC
#include <external/GLEW\glew.h> // OpenGL functions
#include <external/GLFW\glfw3.h>// window & input

#define GLM_ENABLE_EXPERIMENTAL
#include <external/GLM/glm.hpp> //for math

using glm::vec3;
using glm::ivec3;
using glm::uvec3;

#include <iostream>
#include <vector>
#include <algorithm>

#define Alloc(type, count) (type *)calloc(count, sizeof(type))

typedef unsigned char byte;
typedef unsigned int  uint;

// randomness & noise

#define BIT_NOISE_1 0xB5297A4D;
#define BIT_NOISE_2 0x68E31DA4;
#define BIT_NOISE_3 0x1B56C4E9;

uint random_uint()
{
	uint seed = __rdtsc();
	seed *= BIT_NOISE_1;
	seed *= seed; // helps avoid linearity
	seed ^= (seed >> 8);
	seed += BIT_NOISE_2;
	seed ^= (seed >> 8);
	seed *= BIT_NOISE_3;
	seed ^= (seed >> 8);
	return seed;
}

float random_normalized_float() // random float between 0 and 1
{
	uint seed = random_uint();
	return (float)seed / (float)UINT_MAX; // is there a better way to do this?
}

struct Window
{
	char window_title[64];
	GLFWwindow* handle;
	int width, height;
};

void init(Window* c)
{
	c->width  = 1500;
	c->height = 1000;

	if (!glfwInit()) { std::cerr << "glfwInit() failed!\n"; exit(-1); }

	c->handle = glfwCreateWindow(c->width, c->height, "Smoke", NULL, NULL);
	if (!c->handle) { std::cerr << "Create Window failed\n"; exit(-1); }

	glfwMakeContextCurrent(c->handle);
	//glfwSetWindowPos(c->handle, 0, 0);

	glClearColor(0, 0, 0, 1);

	if (glewInit() != GLEW_OK) { std::cerr << "glewInit() failed!\n"; exit(-1); }
}