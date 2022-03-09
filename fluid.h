#include "renderer.h"

#define IDX(x,y,z) (((x)*(RES)*(RES))+((y)*(RES))+(z)) // index into 3D array

void Advect(float* quantityTmp, float* quantity, float* velX, float* velY, float* velZ, float dt)
{
	float dt0 = dt * N;

	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		float xx = i - dt0 * velX[IDX(i, j, k)];
		float yy = j - dt0 * velY[IDX(i, j, k)];
		float zz = k - dt0 * velZ[IDX(i, j, k)];

		if (xx < 0.5) xx = 0.5f; if (xx > N + 0.5) xx = N + 0.5f; int i0 = (int)xx; int i1 = i0 + 1;
		if (yy < 0.5) yy = 0.5f; if (yy > N + 0.5) yy = N + 0.5f; int j0 = (int)yy; int j1 = j0 + 1;
		if (zz < 0.5) zz = 0.5f; if (zz > N + 0.5) zz = N + 0.5f; int k0 = (int)zz; int k1 = k0 + 1;

		float sx1 = xx - i0; float sx0 = 1 - sx1;
		float sy1 = yy - j0; float sy0 = 1 - sy1;
		float sz1 = zz - k0; float sz0 = 1 - sz1;

		float v0 = sx0 * (sy0 * quantityTmp[IDX(i0, j0, k0)] + sy1 * quantityTmp[IDX(i0, j1, k0)]) + sx1 * (sy0 * quantityTmp[IDX(i1, j0, k0)] + sy1 * quantityTmp[IDX(i1, j1, k0)]);
		float v1 = sx0 * (sy0 * quantityTmp[IDX(i0, j0, k1)] + sy1 * quantityTmp[IDX(i0, j1, k1)]) + sx1 * (sy0 * quantityTmp[IDX(i1, j0, k1)] + sy1 * quantityTmp[IDX(i1, j1, k1)]);
		
		quantity[IDX(i, j, k)] = sz0 * v0 + sz1 * v1;
	} } }
}
void Diffuse(float* velTmp, float* vel, float diff, float dt)
{
	float a = dt * diff * N * N * N;

	for (int l = 0; l < 20; l++)
	{
		for (int i = 1; i <= (N); i++) {
		for (int j = 1; j <= (N); j++) {
		for (int k = 1; k <= (N); k++)
		{
			vel[IDX(i, j, k)] = (velTmp[IDX(i, j, k)] + a * (
				vel[IDX(i - 1, j    , k    )] + vel[IDX(i + 1, j    , k    )] +
				vel[IDX(i    , j - 1, k    )] + vel[IDX(i    , j + 1, k    )] +
				vel[IDX(i    , j    , k - 1)] + vel[IDX(i    , j    , k + 1)])) / (1 + 6 * a);
		} } } // EnforceBoundary(b, vel);
	}
}

struct Fluid
{
	struct info {
		float x[SIZE] = {}, y[SIZE] = {}, z[SIZE] = {};
	} vel, prev_vel;

	vec3 velocities[SIZE], prev_velocities[SIZE];

	float _density[SIZE] = {}, _densityTmp[SIZE] = {};
};

void GenerateSmoke(Fluid* f, vec3 vel)
{
	const int centerY = RES / 8;
	const int centerZ = RES / 8;
	float dens = random_normalized_float();

	for (int i = 1; i <= N; i++) {
	for (int j = 1; j <= N; j++)
	{
		if (hypot(i - centerY, j - centerZ) < RES / 10)
		{
			f->_density[IDX(5, i, j)] = dens;
			f->vel.x[IDX(5, i, j)] = vel.x;
			f->vel.y[IDX(5, i, j)] = vel.y;
			f->vel.z[IDX(5, i, j)] = vel.z;
		}
	} }
}

void Project(Fluid* f)
{
	float* p   = f->prev_vel.x;
	float* div = f->prev_vel.y;	// temporary buffers, use old velocity buffers
	int l;
	float h;
	h = 1.0f / N;

	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		div[IDX(i, j, k)] = -h * (
			f->vel.x[IDX(i + 1, j, k)] - f->vel.x[IDX(i - 1, j, k)] +
			f->vel.y[IDX(i, j + 1, k)] - f->vel.y[IDX(i, j - 1, k)] +
			f->vel.z[IDX(i, j, k + 1)] - f->vel.z[IDX(i, j, k - 1)]) / 3;
		p[IDX(i, j, k)] = 0;
	} } }

	for (l = 0; l < 20; l++)
	{
		for (int i = 1; i <= (N); i++) {
		for (int j = 1; j <= (N); j++) {
		for (int k = 1; k <= (N); k++)
		{
			p[IDX(i, j, k)] = (div[IDX(i, j, k)] +
				p[IDX(i - 1, j, k)] + p[IDX(i + 1, j, k)] +
				p[IDX(i, j - 1, k)] + p[IDX(i, j + 1, k)] +
				p[IDX(i, j, k - 1)] + p[IDX(i, j, k + 1)]) / 6;
		} } }
	}

	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		f->vel.x[IDX(i, j, k)] -= (p[IDX(i + 1, j, k)] - p[IDX(i - 1, j, k)]) / 3 / h;
		f->vel.y[IDX(i, j, k)] -= (p[IDX(i, j + 1, k)] - p[IDX(i, j - 1, k)]) / 3 / h;
		f->vel.z[IDX(i, j, k)] -= (p[IDX(i, j, k + 1)] - p[IDX(i, j, k - 1)]) / 3 / h;
	} } }
}

void VorticityConfinement(Fluid* f, float dt0)
{
	//temp buffers
	float* curlX = f->prev_vel.x, * curlY = f->prev_vel.y, * curlZ = f->prev_vel.z, * curl = f->_densityTmp;

	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		int ijk = IDX(i, j, k);
		// curlx = dw/dy - dv/dz
		curlX[ijk] = (f->vel.z[IDX(i, j + 1, k)] - f->vel.z[IDX(i, j - 1, k)]) * 0.5f - (f->vel.y[IDX(i, j, k + 1)] - f->vel.y[IDX(i, j, k - 1)]) * 0.5f;

		// curly = du/dz - dw/dx
		curlY[ijk] = (f->vel.x[IDX(i, j, k + 1)] - f->vel.x[IDX(i, j, k - 1)]) * 0.5f - (f->vel.z[IDX(i + 1, j, k)] - f->vel.z[IDX(i - 1, j, k)]) * 0.5f;

		// curlz = dv/dx - du/dy
		curlZ[ijk] = (f->vel.y[IDX(i + 1, j, k)] - f->vel.y[IDX(i - 1, j, k)]) * 0.5f - (f->vel.x[IDX(i, j + 1, k)] - f->vel.x[IDX(i, j - 1, k)]) * 0.5f;

		// curl = |curl|
		curl[ijk] = sqrtf(curlX[ijk] * curlX[ijk] + curlY[ijk] * curlY[ijk] + curlZ[ijk] * curlZ[ijk]);
	} } }

	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		int ijk = IDX(i, j, k);

		float nX = (curl[IDX(i + 1, j, k)] - curl[IDX(i - 1, j, k)]) * 0.5f;
		float nY = (curl[IDX(i, j + 1, k)] - curl[IDX(i, j - 1, k)]) * 0.5f;
		float nZ = (curl[IDX(i, j, k + 1)] - curl[IDX(i, j, k - 1)]) * 0.5f;

		float len1 = 1.0f / (sqrtf(nX * nX + nY * nY + nZ * nZ) + 0.0000001f);

		nX *= len1;
		nY *= len1;
		nZ *= len1;

		f->vel.x[ijk] += (nY * curlZ[ijk] - nZ * curlY[ijk]) * dt0;
		f->vel.y[ijk] += (nZ * curlX[ijk] - nX * curlZ[ijk]) * dt0;
		f->vel.z[ijk] += (nX * curlY[ijk] - nY * curlX[ijk]) * dt0;
	} } }
}

void VelocityStep(Fluid* f, float dt)
{
	float viscosity = 0.f;
	float buoyancy  = 4.0f;

	// buoyancy
	for (int i = 0; i < SIZE; i++) f->vel.y[i] += f->_density[i] * buoyancy * dt; // FIXME

	float vc_eps = 5.0f;
	VorticityConfinement(f, vc_eps * dt);

	// diffuse

	Fluid::info* temp = Alloc(Fluid::info, 1);
	
	*temp = f->vel;
	f->vel = f->prev_vel;
	f->prev_vel = *temp;

	Diffuse(f->prev_vel.x, f->vel.x, viscosity, dt);
	Diffuse(f->prev_vel.y, f->vel.y, viscosity, dt);
	Diffuse(f->prev_vel.z, f->vel.z, viscosity, dt);
	Project(f);

	// advect

	*temp = f->vel;
	f->vel = f->prev_vel;
	f->prev_vel = *temp;

	Advect(f->prev_vel.x, f->vel.x, f->prev_vel.x, f->prev_vel.y, f->prev_vel.z, dt);
	Advect(f->prev_vel.y, f->vel.y, f->prev_vel.x, f->prev_vel.y, f->prev_vel.z, dt);
	Advect(f->prev_vel.z, f->vel.z, f->prev_vel.x, f->prev_vel.y, f->prev_vel.z, dt);
	Project(f);

	free(temp);
}
void DensityStep(Fluid* f, float dt)
{
	float diffusion = 0.00001f;

	// diffuse
	std::swap(f->_densityTmp, f->_density);
	Diffuse(f->_densityTmp, f->_density, diffusion, dt);

	// advect
	std::swap(f->_densityTmp, f->_density);
	Advect(f->_densityTmp, f->_density, f->vel.x, f->vel.y, f->vel.z, dt);

	//decrease density
	for (int i = 1; i <= (N); i++) {
	for (int j = 1; j <= (N); j++) {
	for (int k = 1; k <= (N); k++)
	{
		f->_density[IDX(k, j, i)] -= 0.001;
		if (f->_density[IDX(k, j, i)] < 0)
			f->_density[IDX(k, j, i)] = 0;
	} } }
}