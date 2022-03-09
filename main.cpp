#include "fluid.h"

struct Timer
{
	double start_time;
	double stop_time;

	void start()
	{
		start_time = glfwGetTime();
	}
	double restart()
	{
		stop_time = glfwGetTime();
		double elapsed = stop_time - start_time;
		start_time = stop_time;
		return elapsed;
	}
};

int main(int argc, char **argv)
{
	Window* window  = Alloc(Window, 1);
	Renderer* renderer = Alloc(Renderer, 1);
	Fluid* fluid = Alloc(Fluid , 1);

	init(window);
	init(renderer, fluid->_density);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	float x_angle =  10;
	float y_angle = -20;

	glTranslatef(0, 0, -3.8);
	glRotatef(x_angle, 1, 0, 0);
	glRotatef(y_angle, 0, 1, 0);
	
	float dt = 1.f / 20;
	vec3 light_pos = vec3{ -1.2f, .2f, 1.2f };

	Timer timer = {};
	timer.start();

	while (!glfwWindowShouldClose(window->handle))
	{
		set_light_position(renderer, light_pos);

		snprintf(window->window_title, 64, "Smoke | FPS: %4f", 1.f / timer.restart());
		glfwSetWindowTitle(window->handle, window->window_title);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		gluPerspective(50.f, (GLfloat)window->width / window->height, 0.1f, 100.f);
		glViewport(0, 0, window->width, window->height);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		GenerateSmoke(fluid, vec3(0, 1, 0));
		VelocityStep (fluid, dt);
		DensityStep  (fluid, dt);

		FillTexture(renderer->texture_data, renderer->volume_data, -light_pos, renderer->_rayTemplate);

		GLdouble mvMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

		glDisable(GL_DEPTH_TEST);

		draw_cube(renderer->_cubeVertices);
		draw_smoke(mvMatrix, renderer->_cubeVertices, renderer->_cubeEdges, false, vec3(1));

		glEnable(GL_DEPTH_TEST);

		glfwSwapBuffers(window->handle);
		glfwPollEvents();
	}

	return 0;
}