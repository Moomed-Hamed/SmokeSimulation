#include "window.h"

#define RES 36      // box resolution
#define N ((RES)-2) // valid simulation area
#define SIZE ((RES)*(RES)*(RES))

#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)
#define SLICE_NUM 64.0f

struct Convexcomp
{
	const vec3& p0, &up;
	Convexcomp(const vec3& p0, const vec3& up) : p0(p0), up(up) {}
	bool operator()(const vec3& a, const vec3& b) const
	{
		vec3 va = a - p0, vb = b - p0;
		return dot(up, cross(va, vb)) >= 0;
	}
};

struct Renderer
{
	byte* texture_data;
	GLuint texture;

	// lighting infomation
	int _rayTemplate[4096][3];
	float* volume_data;

	GLfloat _cubeVertices[8][3];
	GLfloat _cubeEdges[12][2][3];
};

void init(Renderer* r, float* volumeData)
{
	r->texture_data = Alloc(byte, SIZE * 4);
	r->volume_data = volumeData;

	// cube vertices
	GLfloat cv[][3] = {
		{1.0f, 1.0f,  1.0f}, {-1.0f, 1.0f,  1.0f}, {-1.0f, -1.0f,  1.0f}, {1.0f, -1.0f,  1.0f},
		{1.0f, 1.0f, -1.0f}, {-1.0f, 1.0f, -1.0f}, {-1.0f, -1.0f, -1.0f}, {1.0f, -1.0f, -1.0f}
	};

	// cube edges have the form edges[n][0][xyz] + t*edges[n][1][xyz]
	GLfloat ce[12][2][3] = {
		{{ 1.0f,  1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f,  1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{ 1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},

		{{ 1.0f, -1.0f,  1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f,  1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},
		{{ 1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},

		{{-1.0f,  1.0f,  1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f,  1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f,  1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}}
	};

	// Init GL
	memcpy(r->_cubeVertices, cv, sizeof(r->_cubeVertices));
	memcpy(r->_cubeEdges   , ce, sizeof(r->_cubeEdges   ));

	glEnable(GL_TEXTURE_3D);
	glDisable(GL_DEPTH_TEST);
	glCullFace(GL_FRONT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glGenTextures(1, &r->texture);
	glBindTexture(GL_TEXTURE_3D, r->texture);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
}

void light_ray(int x, int y, int z, int n, float decay, const float* dens, byte* intensity, int _rayTemplate[4096][3])
{
	const int AMBIENT = 100;

	int xx = x, yy = y, zz = z;
	int offset;

	int l = 200;
	float d;

	int i = 0;
	do {
		offset = ((xx * n) + yy) * n + zz; // FIXME

		if (intensity[offset] > 0)
			intensity[offset] = (byte)((intensity[offset] + l) * 0.5f);
		else
			intensity[offset] = (byte)l;

		d = dens[offset] * 255.0f;

		if (l > AMBIENT)
		{
			l -= d * decay;
			if (l < AMBIENT)
				l = AMBIENT;
		}

		i++;
		xx = x + _rayTemplate[i][0];
		yy = y + _rayTemplate[i][1];
		zz = z + _rayTemplate[i][2];
	} while ((xx >= 0) && (xx < n) && (yy >= 0) && (yy < n) && (zz >= 0) && (zz < n));
}
void cast_light(int edge_len, const float* dens, byte* intensity, vec3 light_dir, int _rayTemplate[4096][3])
{
	const float DECAY = 0.06f;

	int sx = (light_dir[0] > 0) ? 0 : edge_len - 1;
	int sy = (light_dir[1] > 0) ? 0 : edge_len - 1;
	int sz = (light_dir[2] > 0) ? 0 : edge_len - 1;

	float decay = 1.f / (edge_len * DECAY);

	for (int i = 0; i < edge_len; i++) {
	for (int j = 0; j < edge_len; j++)
	{
		if (abs(light_dir.x) > 0) light_ray(sx, i , j , edge_len, decay, dens, intensity, _rayTemplate);
		if (abs(light_dir.y) > 0) light_ray(i , sy, j , edge_len, decay, dens, intensity, _rayTemplate);
		if (abs(light_dir.z) > 0) light_ray(i , j , sz, edge_len, decay, dens, intensity, _rayTemplate);
	} }
}
void generate_ray_template(int edgeLen, int _rayTemplate[4096][3], vec3 light_dir)
{
	_rayTemplate[0][0] = _rayTemplate[0][2] = _rayTemplate[0][2] = 0;

	vec3 f = {};
	vec3 t = {};
	vec3 l = light_dir + +0.000001f;
	ivec3 u = {};
	uvec3 inc = { (l.x > 0) ? 1 : -1 , (l.y > 0) ? 1 : -1 , (l.z > 0) ? 1 : -1 };
	int i = 1;
	int len = 0;
	int maxlen = 3 * edgeLen * edgeLen;

	while (len <= maxlen)
	{
		// f.x + t * l.x = (x + 1) -> t = (x + 1 - f.x) / l.x
		t.x = (u.x + inc.x - f.x) / l.x;
		t.y = (u.y + inc.y - f.y) / l.y;
		t.z = (u.z + inc.z - f.z) / l.z;

		if ((t.x <= t.y) && (t.x <= t.z))
		{
			_rayTemplate[i][0] = _rayTemplate[i - 1][0] + inc.x;
			u.x = inc.x; // is this an error?
			f.x = u.x;

			if (ALMOST_EQUAL(t.y, t.x)) {
				_rayTemplate[i][1] = _rayTemplate[i - 1][1] + inc.y;
				u.y += inc.y;
				f.y = u.y;
			}
			else {
				_rayTemplate[i][1] = _rayTemplate[i - 1][1];
				f.y += t.x * l.y;
			}

			if (ALMOST_EQUAL(t.z, t.x)) {
				_rayTemplate[i][2] = _rayTemplate[i - 1][2] + inc.z;
				u.z += inc.z;
				f.z = u.z;
			}
			else {
				_rayTemplate[i][2] = _rayTemplate[i - 1][2];
				f.z += t.x * l.z;
			}
		}
		else if ((t.y < t.x) && (t.y <= t.z)) {
			_rayTemplate[i][0] = _rayTemplate[i - 1][0];
			f.x += t.y * l.x;

			_rayTemplate[i][1] = _rayTemplate[i - 1][1] + inc.y;
			u.y += inc.y;
			f.y = u.y;

			if (ALMOST_EQUAL(t.z, t.y)) {
				_rayTemplate[i][2] = _rayTemplate[i - 1][2] + inc.z;
				u.z += inc.z;
				f.z = u.z;
			}
			else {
				_rayTemplate[i][2] = _rayTemplate[i - 1][2];
				f.z += t.y * l.z;
			}
		}
		else {
			assert((t.z < t.x) && (t.z < t.y));
			_rayTemplate[i][0] = _rayTemplate[i - 1][0];
			f.x += t.z * l.x;
			_rayTemplate[i][1] = _rayTemplate[i - 1][1];
			f.y += t.z * l.y;
			_rayTemplate[i][2] = _rayTemplate[i - 1][2] + inc.z;
			u.z += inc.z;
			f.z = u.z;
		}

		len = _rayTemplate[i][0] * _rayTemplate[i][0]
			 + _rayTemplate[i][1] * _rayTemplate[i][1]
			 + _rayTemplate[i][2] * _rayTemplate[i][2];
		i++;
	}
}
void set_light_position(Renderer* r, vec3 light_pos)
{
	vec3 light_dir = -light_pos;
	generate_ray_template(RES, r->_rayTemplate, light_dir);
}

// intersect a plane with the cube; helper function for drawing slices
// plane equation is Ax + By + Cz + D = 0
typedef std::vector<vec3> Point;
Point IntersectEdges(float A, float B, float C, float D, GLfloat cube_edges[12][2][3])
{
	Point ret;
	GLfloat(*edges)[2][3] = cube_edges;

	for (int i = 0; i < 12; i++)
	{
		float t = -(A * edges[i][0][0] + B * edges[i][0][1] + C * edges[i][0][2] + D) / 
			        (A * edges[i][1][0] + B * edges[i][1][1] + C * edges[i][1][2]);

		vec3 p = {};

		if (t > 0 && t < 2)
		{
			p[0] = edges[i][0][0] + edges[i][1][0] * t;
			p[1] = edges[i][0][1] + edges[i][1][1] * t;
			p[2] = edges[i][0][2] + edges[i][1][2] * t;
			ret.push_back(p);
		}
	}

	return ret;
}

// draw the slices, mvMatrix must be the MODELVIEW_MATRIX
void draw_smoke(GLdouble mvMatrix[16], GLfloat _cubeVertices[8][3], GLfloat _cubeEdges[12][2][3], bool draw_slices, vec3 color = {1,1,1})
{
	vec3 viewdir = normalize(vec3(-mvMatrix[2], -mvMatrix[6], -mvMatrix[10]));	//FIXME

	// find cube vertex that is closest to the viewer
	GLfloat(*cv)[3] = _cubeVertices;

	int i;
	for (i = 0; i < 8; i++)
	{
		float x = cv[i][0] + viewdir.x;
		float y = cv[i][1] + viewdir.y;
		float z = cv[i][2] + viewdir.z;

		if ((x >= -1.0f) && (x <= 1.0f) && (y >= -1.0f) && (y <= 1.0f) && (z >= -1.0f) && (z <= 1.0f))
		{
			break;
		}
	} assert(i != 8);

	// slices are defined by the plane equation A*x + B*y + C*z + D = 0
	// (a,b,c), the plane normal, are given by viewdir
	// d is the parameter along the view direction
	// the first d is given by inserting previously found vertex into the plane equation
	vec3 c = vec3(cv[i][0], cv[i][1], cv[i][2]);
	float d0 = dot(c, viewdir) * -1.f;
	float dStep = 2 * d0 / SLICE_NUM;

	for (float d = -d0; d < d0; d += dStep)
	{
		// IntersectEdges returns the intersection points of all cube edges with the given plane that lie within the cube
		Point pt = IntersectEdges(viewdir.x, viewdir.y, viewdir.z, d, _cubeEdges);

		if (pt.size() < 3) continue;

		// sort points to get a convex polygon
		std::sort(pt.begin() + 1, pt.end(), Convexcomp(pt[0], viewdir));

		glEnable(GL_TEXTURE_3D);
		glBegin(GL_POLYGON);
		for (i = 0; i < pt.size(); i++)
		{
			vec3 pos = vec3(pt[i][0], pt[i][1], pt[i][2]);
			vec3 coord = (pos + vec3(1)) / vec3(2);

			glColor3f(1, (float)i / pt.size(), (float)i / pt.size() * 2);
			//glColor3f(color.x, color.y, color.z);
			glTexCoord3f(coord.x, coord.y, coord.z); // FIXME
			glVertex3f(pos.x, pos.y, pos.z);
		}
		glEnd();

		if (draw_slices)
		{
			glDisable(GL_TEXTURE_3D);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glBegin(GL_POLYGON);
			for (i = 0; i < pt.size(); i++)
			{
				glColor3f(0, 1, 1);
				glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
			}
			glEnd();
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
	}

	glDisable(GL_TEXTURE_3D);
}

// generate texture from smoke density
void FillTexture(byte* texture_data, float* volume_data, vec3 light_dir, int _rayTemplate[4096][3])
{
	float* density = volume_data;
	byte* intensity = Alloc(byte, SIZE);

	cast_light(RES, density, intensity, light_dir, _rayTemplate);

	// FIXME: It is important to beware that texture coordinate is in reverse order of the simulation coordinate
	for (int i = 0; i < RES; i++) {
	for (int j = 0; j < RES; j++) {
	for (int k = 0; k < RES; k++)
	{
		int texIndex  = i * RES * RES + j * RES + k; // reverse order
		int densIndex = k * RES * RES + j * RES + i;

		byte c = intensity[densIndex];
		float d = density[densIndex];

		texture_data[texIndex * 4 + 0] = c;
		texture_data[texIndex * 4 + 1] = c;
		texture_data[texIndex * 4 + 2] = c;
		texture_data[texIndex * 4 + 3] = (d < 0.1f) ? (d * 2550.f) : 255;
	} } }

	free(intensity);

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, RES, RES, RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, texture_data);
}

// draw outline of the cube containing the smoke
void draw_cube(GLfloat _cubeVertices[8][3])
{
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_FRAGMENT_PROGRAM_ARB);

	glEnable(GL_CULL_FACE);
	vec3 color = 0.1f * vec3(.9, .8, .4);
	glColor4f(color.x, color.y, color.z, .5f);
	GLfloat(*cv)[3] = _cubeVertices;

	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]); glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[4]); glVertex3fv(cv[5]); glVertex3fv(cv[1]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[3]); glVertex3fv(cv[7]); glVertex3fv(cv[4]);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3fv(cv[7]); glVertex3fv(cv[6]); glVertex3fv(cv[5]); glVertex3fv(cv[4]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]); glVertex3fv(cv[7]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]); glVertex3fv(cv[6]); glVertex3fv(cv[2]);
	glEnd();

	glDisable(GL_CULL_FACE);

	glBegin(GL_LINES);
	glColor4f(.6f, .6f, .6f, .6f);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]);
	glVertex3fv(cv[1]); glVertex3fv(cv[2]);
	glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glVertex3fv(cv[3]); glVertex3fv(cv[0]);

	glVertex3fv(cv[4]); glVertex3fv(cv[5]);
	glVertex3fv(cv[5]); glVertex3fv(cv[6]);
	glVertex3fv(cv[6]); glVertex3fv(cv[7]);
	glVertex3fv(cv[7]); glVertex3fv(cv[4]);

	glVertex3fv(cv[0]); glVertex3fv(cv[4]);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]);
	glVertex3fv(cv[3]); glVertex3fv(cv[7]);
	glEnd();
}
