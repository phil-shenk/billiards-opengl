#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "initShader.h"
// phil's linear albegra library, github.com/phil-shenk/linalglib
#include "../linalglib/linalglib.h"

#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))

// 2D vector struct to represent UV texture coordinates
typedef struct
{
    GLfloat x;
    GLfloat y;
} vec2;

// number of total vertices (billiard balls, light sphere, floor)
int num_vertices;

// number of billiard balls to generate
const int numBalls = 16;

// number of vertices per billiard ball
int numBallVerts;

// number of vertices for the floor
int numFloorVerts;

// how many sections on the edge of the floor (number of floor tiles = floorDivisions^2)
int floorDivisions;

// boolean keydown u,r,d,l,+,-
int keysdown[6] = {0,0,0,0,0,0};

// length (and width) of the square ground object
const float groundSize = 10.0;

// vectors to store camera orientation information
vec4 eyePoint;
vec4 relEyePoint;
vec4 atPoint;
vec4 upVector;
int following = -1;
int orbitFollow = 0; // enable/disable orbiting mode for following a ball
vec4 followMotion;


// LIGHTING INFORMATION
vec4 light_position = {0.0, 1.0, 0.0, 1.0};

float attenuation_constant  = 0.1;
float attenuation_linear    = 0.1;
float attenuation_quadratic = 0.5;

// shininess of the ball material
GLfloat shininess = 555.2;

// stores the state of the program, used for different phases of animation
int mode = 0;

// use same ctm location, but update which CTM is stored in there 16 times per display() call
GLuint ctm_location;
GLuint use_texture_location;
GLuint light_position_location;
GLuint is_shadow_location;

// need a different ctm for each ball
mat4 *ctms;
vec4 *ballpositions;

GLuint model_view_location;
mat4 model_view_matrix = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1}
};

GLuint projection_location;
mat4 projection_matrix = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1}
};

int window_size = 512;

/*
 * get the 'polar indices' of theta and phi (thetai, phii, respectively) from the given index
 * values are returned by populating the supplied array
 */
void getPolarIndicesFromIndex(int i, int nNS, int nWE, int polarIndices[2])
{
	int thetai;
	int phii;

	// num of verts per vertical sliver (time zone)
	int nSV = (nNS - 1)*6;
	// num of this sliver
	int sn = i / nSV;

	// index within this sliver
	int si = i % nSV;

	// check handle case of north and south pole
	if (si == 0){
		// set theta
		polarIndices[0] = 0;
		// set phi
		polarIndices[1] = 0;
		return;
	}
	if (si == nSV - 3){
		// set theta
		polarIndices[0] = nNS; // had -1
		// set phi
		polarIndices[1] = 0;
		return;
	}
	
	// calculate phi and theta (see onenote)
	if (si % 2 == 0){
		phii = sn + 1;
		int k = (si/2)-1;
		thetai = (k/3)+1;
	}else{
		phii = sn; 
		int k = (si-1)/2;
		thetai = (k/3) + 1;
		if (k%3 == 1)
			thetai++;
	}

	polarIndices[0] = thetai;
	polarIndices[1] = phii;

	return;
}
	
/* 
 * calculate polar coordinate representation vector of a given index of a sphere
 * sphere has nNS 'rows' (stacking from north to south) and nWE 'columns' (slivers)
 * stacking from west to east
 * see onenote for details
 */
vec4 getPolarFromIndex(int i, float ballRadius, int nNS, int nWE)
{
	float dtheta = M_PI / (float)nNS;
	float dphi = 2*M_PI / (float)nWE;
	float R = ballRadius;

	int polarIndices[2] = {0,0};
	getPolarIndicesFromIndex(i, nNS, nWE, polarIndices);

	float theta = dtheta * polarIndices[0];
	float phi   = dphi * polarIndices[1];

	vec4 result = {R, theta, phi, 1};
	return result;
}

/*
 * return cartesian representation of polar vector (or point, if w=1)
 */
vec4 getCartesianFromPolar(vec4 polarvector)
{
	// extract polar values from recycled 'cartesianly-named' vec4 struct..
	float r     = polarvector.x; 
	float theta = polarvector.y;
	float phi   = polarvector.z;
	float w     = polarvector.w;
	float x =    r*sin(theta)*cos(phi);
	float y =    r*cos(theta);
	float z = -1*r*sin(theta)*sin(phi);
	vec4 result = {x,y,z,w};
	return result;
}
	

/*
 * get the ith vertex of a ball with nNS altitudinal segments and nWE azimuthal segments
 */
vec4 getBallVertex(int i, float ballRadius, int nNS, int nWE)
{

	vec4 result = getCartesianFromPolar( getPolarFromIndex(i, ballRadius, nNS, nWE) );
	return result;
}

/*
 * get the ith normal of a ball with nNS altitudinal segments and nWE azimuthal segments
 */
vec4 getBallNormal(int i, float ballRadius, int nNS, int nWE)
{
	// the normal of a sphere at a point is just the normalized vector from the center of the sphere to the point
	vec4 point = getCartesianFromPolar( getPolarFromIndex(i, ballRadius, nNS, nWE) );
	// make the point a vector, it's already relative to the origin so we needn't worry about that
	point.w = 0;
	vec4 normal = normalize_v4( point );
	return normal;
}

/*
 * get the ith texture coordinate of ball 'ballnumber' with nNS altitudinal segments and nWE azimuthal segments
 */
vec2 getBallTextureCoordinate(int i, int ballnumber, int nNS, int nWE)
{
	//get polar indices from index
	int polarIndices[2] = {0,0};
	getPolarIndicesFromIndex(i, nNS, nWE, polarIndices); 

	int thetai = polarIndices[0];
	int phii   = polarIndices[1];

	// mapping of ball textures from the source texture image
	vec2 topLeftCoordinates[16] = {
		{0.76, 0.7625},
		{0.01, 0.01},
		{0.26, 0.01},
		{0.51, 0.01},
		{0.76, 0.01},
		{0.01, 0.26},
		{0.26, 0.26},
		{0.51, 0.26},
		{0.76, 0.26},
		{0.01, 0.511},
		{0.26, 0.511},
		{0.51, 0.511},
		{0.76, 0.511},
		{0.01, 0.7625},
		{0.26, 0.7625},
		{0.51, 0.7625}
	};
	vec2 bottomRightCoordinates[16] = {
		{0.99, 0.99},
		{0.24, 0.24},
		{0.49, 0.24},
		{0.74, 0.24},
		{0.99, 0.24},
		{0.24, 0.49},
		{0.49, 0.49},
		{0.74, 0.49},
		{0.99, 0.49},
		{0.24, 0.741},
		{0.49, 0.741},
		{0.74, 0.741},
		{0.99, 0.741},
		{0.24, 0.99},
		{0.49, 0.99},
		{0.74, 0.99}
	};

	vec2 TLC = topLeftCoordinates[ballnumber];
	vec2 BRC = bottomRightCoordinates[ballnumber];

	float dx = (BRC.x - TLC.x)*2 / (float)nWE;
	float dy = (BRC.y - TLC.y)   / (float)nNS;

	int nSV = (nNS - 1)*6;
	int si = i % nSV;

	int ti = phii % (nWE/2);
	// correct even seam cases
	if(ti == 0 && (si%2 == 0)){
		ti=phii;
		if(phii == nWE){
			ti = phii/2;
		}
	}

	int tj = thetai;

	vec2 result = {ti*dx+TLC.x, tj*dy+TLC.y};
	return result;
}

/*
 * populate the vertices that represent the floor
 */
void buildFloor(int vertoffset, vec4 vertices[], vec4 colors[], vec4 normals[], vec2 tex_coords[]){
	for(int i=0; i<floorDivisions; i++){
		for(int j=0; j<floorDivisions; j++){
			//populate vertices
			float dx = groundSize / ((float)floorDivisions);
			float dz = dx;
			float floorOffset = -groundSize/2.0;
			int tileOffset = 6*(i*floorDivisions + j);
			vec4 A = {j*dx+floorOffset, 0, i*dz+floorOffset, 1};
		       	vec4 B = {(j+1)*dx+floorOffset, 0, i*dz+floorOffset, 1};
			vec4 C = {j*dx+floorOffset, 0, (i+1)*dz+floorOffset, 1};
			vec4 D = {(j+1)*dx+floorOffset, 0, (i+1)*dx+floorOffset, 1};
			vertices[tileOffset + 0] = A;
			vertices[tileOffset + 1] = C;
			vertices[tileOffset + 2] = B;
			vertices[tileOffset + 3] = B;
			vertices[tileOffset + 4] = C;
			vertices[tileOffset + 5] = D;

			//populate colors and normals and bogus tex_coords
			vec4 n = {0.0, 1.0, 0.0, 0.0};
			for(int v=0; v<6; v++){
				vec4 green = {0.0, 1.0, 0.0, 1.0};
				colors[tileOffset + v] = green;
				normals[tileOffset + v] = n;
				vec2 t = {0.0,0.0};
				tex_coords[tileOffset + v] = t;// not used, since the floor is a solid color
			}
		}
	}
}

/*
 * build ball and place values into vertices, colors, normals, and tex_coords to transfer to the graphics card 
 * (nNS = num segments north to south, i.e.numAltitudinalSegments, nWE = num segments west to east, i.e. numAzimuthalSegments)
 */
void buildBall(int vertoffset, int numBallVerts, int nNS, int nWE, int ballnumber, int flipNormals, float ballRadius, vec4 vertices[], vec4 colors[], vec4 normals[], vec2 tex_coords[]){
	// i = relative index of vertex in supplied vertices[], colors[], normals[], and tex_coords[] arrays
	for(int i=0; i<numBallVerts; i++){
		// get ball vertex value
		vertices[i+vertoffset] = getBallVertex(i, ballRadius, nNS, nWE);

		// and populate the color values into the 'colors' array (not used, since textures are used)
		vec4 c = {1.0, 1.0, 1.0, 1.0};
		colors[i+vertoffset] = c;

		// populate the normals
		normals[i+vertoffset] = getBallNormal(i, ballRadius, nNS, nWE);
		if(flipNormals == 1){
			normals[i+vertoffset] = scale_v4(-1.0, getBallNormal(i, ballRadius, nNS, nWE));
		}

		// and populate the tex_coords, attempting to get fancy with my approach
		// format vec2 = {0,0.5}
		tex_coords[i+vertoffset] = getBallTextureCoordinate(i, ballnumber, nNS, nWE);
	}
}

/*
 * initialize relevant transformation matrices, vertex counts, textures, and graphics memory sectors
 */
void init()
{
	// INITIAL MODEL VIEW
	vec4 e = {0.0, 2.0, 5.0, 1.0};
	vec4 a = {0, 0, 0, 1};
	vec4 u = {0, 1, 0, 0};
	eyePoint = e;
	relEyePoint = e;
	atPoint  = a;
	upVector = u;
	model_view_matrix = look_at(eyePoint, atPoint, upVector);

	vec4 fm = {0.0, 0.0, 1.0, 0.0};
	followMotion = fm;

	float w = -5.0;
	float h = 5.0;
	float tnear = 0.5;
	float tfar = 20.0;

	// get projection matrix from frustum description
 	projection_matrix = frustum(-w, w, -h, h, tfar, tnear);

	// floorDivisions**2 = number of tiles on the floor
	floorDivisions = 100;
	numFloorVerts = floorDivisions*floorDivisions*2*3;

	// generate 16 billiard balls

	// parameters PER SPHERE generated
	int nWE = 24;	// E/W longitude resolution for spheres
	int nNS = 12;	// N/S latitude  resolution for spheres
	numBallVerts = (nWE)*((nNS-1)*6);

	//calculate how many vertices will be necessary
	//             floor vertices  light verts    all billiard ball verts
	num_vertices = numFloorVerts + numBallVerts + numBallVerts*numBalls;

	vec4 vertices[num_vertices];
	vec4 colors[num_vertices];
	vec4 normals[num_vertices];
	vec2 tex_coords[num_vertices];

	// build the floor
	buildFloor(0, vertices, colors, normals, tex_coords);
	
	//build  the light (with flipped normals)
	buildBall(numFloorVerts, numBallVerts, nNS, nWE, 0, 1, 0.05, vertices, colors, normals, tex_coords);

	// allocate the ctms[] array
	ctms = malloc((numBalls + 2) * sizeof(mat4));
	ballpositions = malloc((numBalls + 1) * sizeof(vec4)); // include light

	// put vertices, colors, normals, and tex_coords for each ball into the large arrays
	for(int b=0; b<numBalls; b++){
		printf("b is %d\n", b);
		// make sure vertex offset is calculated based on number of balls already constructed. Also, don't flip vertices
		buildBall(numFloorVerts + numBallVerts + b*numBallVerts, numBallVerts, nNS, nWE, (b+1)%numBalls, 0, 0.1, vertices, colors, normals, tex_coords);
		
		// initialize CTM for each ball as a simple translation
		float x = ((float)((b % 4))-1.5) * 0.2;
		float z = ((float)((b / 4))-1.5)* 0.2;

		mat4 ballctm = translation(x, 0.1, z);
		vec4 ballpos = {x, 0.1, z, 1.0};
		ballpositions[b] = ballpos;
		ctms[b] = arb_rotation_com_world(ballctm, -M_PI/2.0, 1.0, 0.0, 0.0);
		//ctms[b] = translation(x, 0.1, z); // identity is unnecessary here..
	}
	// also want a ctm for the floor
	ctms[numBalls] = translation(0.0, 0.0, 0.0);
	
	// also want a ctm for the light
	ctms[numBalls+1] = translation(0.0, 1.0, 0.0);

	// LOAD TEXTURE
	int width  = 512;
	int height = 512;
	GLubyte my_texels[width][height][3];

	FILE *fp = fopen("pb_512_512.raw", "r");

	fread(my_texels, width * height * 3, 1, fp);

	fclose(fp);

	GLuint program = initShader("vshader.glsl", "fshader.glsl");
	glUseProgram(program);

	GLuint mytex[1];
	glGenTextures(1, mytex);
	glBindTexture(GL_TEXTURE_2D, mytex[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, my_texels);
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

	int param;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &param);

	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors) + sizeof(normals) + sizeof(tex_coords), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices), sizeof(colors), colors);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors), sizeof(normals), normals);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors) + sizeof(normals), sizeof(tex_coords), tex_coords);

	GLuint vPosition = glGetAttribLocation(program, "vPosition");
	glEnableVertexAttribArray(vPosition);
	glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));

	GLuint vColor = glGetAttribLocation(program, "vColor");
	glEnableVertexAttribArray(vColor);
	glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0 + sizeof(vertices));

	GLuint vNormal = glGetAttribLocation(program, "vNormal");
	glEnableVertexAttribArray(vNormal);
	glVertexAttribPointer(vNormal, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0 + sizeof(vertices) + sizeof(colors));

	GLuint vTexCoord = glGetAttribLocation(program, "vTexCoord");
	glEnableVertexAttribArray(vTexCoord);
	glVertexAttribPointer(vTexCoord, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0 + (sizeof(vertices) + sizeof(colors) + sizeof(normals)));

	GLuint texture_location = glGetUniformLocation(program, "texture");
	glUniform1i(texture_location, 0);
	printf("texture_location: %i\n", texture_location);

	use_texture_location = glGetUniformLocation(program, "use_texture");
	glUniform1i(use_texture_location, 0);
	printf("use_texture_location: %i\n", use_texture_location);

	GLuint shininess_location;
	shininess_location = glGetUniformLocation(program, "shininess");
	glUniform1f(shininess_location, shininess);
	printf("shininess_location: %i\n", shininess_location);

	GLuint att_c_location;
	att_c_location = glGetUniformLocation(program, "attenuation_constant");
	printf("att_c_loc: %i\n", att_c_location);
	glUniform1f(att_c_location, attenuation_constant);

	GLuint att_l_location;
	att_l_location = glGetUniformLocation(program, "attenuation_linear");
	printf("att_l_loc: %i\n", att_l_location);
	glUniform1f(att_l_location, attenuation_linear);

	GLuint att_q_location;
	att_q_location = glGetUniformLocation(program, "attenuation_quadratic");
	printf("att_q_loc: %i\n", att_q_location);
	glUniform1f(att_q_location, attenuation_quadratic);

	// make room for light position
	light_position_location = glGetUniformLocation(program, "light_position");
	printf("light_pos_loc %d\n", light_position_location);
	//glUniform4fv(light_position_location, 1, (GLfloat *) &light_position);	
	glUniform4f(light_position_location, light_position.x, light_position.y, light_position.z, light_position.w);

	is_shadow_location = glGetUniformLocation(program, "is_shadow");
	glUniform1i(is_shadow_location, 0);
	printf("is_shadow_loc %d\n", is_shadow_location);

	model_view_location = glGetUniformLocation(program, "model_view_matrix");
	projection_location = glGetUniformLocation(program, "projection_matrix");
	printf("model_view_location: %i\n", model_view_location);
	printf("projection_location: %i\n", projection_location);

	ctm_location = glGetUniformLocation(program, "ctm");

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.7, 0.7, 1.0, 1.0);
	glDepthRange(1,0);
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPolygonMode(GL_BACK, GL_LINE);
	glPolygonMode(GL_FRONT, GL_FILL);

	glUniformMatrix4fv(model_view_location, 1, GL_FALSE, (GLfloat *) &model_view_matrix);
	glUniformMatrix4fv(projection_location, 1, GL_FALSE, (GLfloat *) &projection_matrix);

	// pass the location of the light every frame
	glUniform4f(light_position_location, light_position.x, light_position.y, light_position.z, light_position.w);

	// don't use textures for floor and light
	glUniform1i(use_texture_location, 0);

	glUniform1i(is_shadow_location, 0);

	// draw the floor
	glUniformMatrix4fv(ctm_location, 1, GL_FALSE, (GLfloat *) &ctms[numBalls]);
	glDrawArrays(GL_TRIANGLES, 0, numFloorVerts);

	// draw the light
	glUniformMatrix4fv(ctm_location, 1, GL_FALSE, (GLfloat *) &ctms[numBalls+1]);
	glDrawArrays(GL_TRIANGLES, numFloorVerts, numBallVerts);

	// load up a fresh CTM for each ball, and draw it (with textures)
	glUniform1i(use_texture_location, 1);
	for(int i=0; i<numBalls; i++){
		glUniform1i(is_shadow_location, 0);
		glUniformMatrix4fv(ctm_location, 1, GL_FALSE, (GLfloat *) &ctms[i]);
		//           mode,     first,        count
		//                                   extra numBallVerts is for the light, which is a Similar sphere
		glDrawArrays(GL_TRIANGLES, numFloorVerts + numBallVerts + i*numBallVerts, numBallVerts);

		// put into shadow mode
		glUniform1i(is_shadow_location, 1);
		glDrawArrays(GL_TRIANGLES, numFloorVerts + numBallVerts + i*numBallVerts, numBallVerts);
	}

	glutSwapBuffers();
}

/*
 * follow ball 'ballnumber' with camera hovering behind it
 */
void followBall(int ballnumber){
	vec4 rep = {0.0, 0.3, 1.5, 1.0};
	relEyePoint = rep;

	printf("follow %d\n",ballnumber);
	if(ballnumber == -1){
		following = -1;
		vec4 ap = {0.0, 0.0, 0.0, 1.0};
		atPoint = ap;
		vec4 ep = {0.0, 2.0, 5.0, 1.0};
		relEyePoint = ep;
		eyePoint = ep;
		model_view_matrix = look_at(eyePoint, atPoint, upVector);
		return;
	}
	following = ballnumber;
}


void keyboardDown(unsigned char key, int mousex, int mousey){
	if(key == 'q')
	{
		glutLeaveMainLoop();
	}

	if(key == 'n')
	{
		mode++;
	}
	if(key == '+' || key == '=')
		keysdown[4] = 1;
	if(key == '-')
		keysdown[5] = 1;

	float dl = 0.05;
	if(key == 'i') // -z
		light_position = multiply_m4v4(translation( 0.0,  0.0,  -dl), light_position);
	if(key == 'k') // +z
		light_position = multiply_m4v4(translation( 0.0,  0.0,   dl), light_position);
	if(key == 'l') // +x
		light_position = multiply_m4v4(translation(  dl,  0.0,  0.0), light_position);
	if(key == 'j') // -x
		light_position = multiply_m4v4(translation( -dl,  0.0,  0.0), light_position);
	if(key == 'o') // +y
		light_position = multiply_m4v4(translation( 0.0,   dl,  0.0), light_position);
	if(key == 'u') // -y
		light_position = multiply_m4v4(translation( 0.0,  -dl,  0.0), light_position);


	if(key>=49 && key <= 57)
		followBall(key-48);
	if(key == 48)
		followBall(10);
	if(key>=97 && key<=102)
		followBall(key-86);
	if(key == ' ')
		followBall(-1);

	if(key == 't'){
		if(orbitFollow == 0){
			orbitFollow = 1;
		}else if(orbitFollow == 1){
			orbitFollow = 0;
		}

		printf("toggle orbits %d\n", orbitFollow);
	}
}

void keyboardUp(unsigned char key, int mousex, int mousey){
	if(key == '+' || key == '=')
		keysdown[4] = 0;
	if(key == '-')
		keysdown[5] = 0;

}

void keyboardSpecialDown(int key, int mousex, int mousey)
{
	// boolean keydown u,r,d,l,+,-
	if(key == GLUT_KEY_UP)
		keysdown[0] = 1;
	if(key == GLUT_KEY_RIGHT){
		keysdown[1] = 1;
	}
	if(key == GLUT_KEY_DOWN)
		keysdown[2] = 1;
	if(key == GLUT_KEY_LEFT)
		keysdown[3] = 1;

}
void keyboardSpecialUp(int key, int mousex, int mousey){
	// boolean keydown u,r,d,l,+,-
	if(key == GLUT_KEY_UP)
		keysdown[0] = 0;
	if(key == GLUT_KEY_RIGHT)
		keysdown[1] = 0;
	if(key == GLUT_KEY_DOWN)
		keysdown[2] = 0;
	if(key == GLUT_KEY_LEFT)
		keysdown[3] = 0;
}

void mouse(int button, int state, int x, int y)
{
		if(button == 3)
		{
			//ctm = multiply_m4(ctm, scaling(1.02,1.02,1.02));
			glutPostRedisplay();
		}
		if(button == 4)
		{
			//ctm = multiply_m4(ctm, scaling(1/1.02,1/1.02,1/1.02));
			glutPostRedisplay();
		}

}

void reshape(int width, int height)
{
    glViewport(0, 0, window_size, window_size);
}

/*
* return the ctm for a time 0<t<1 rolling linearly from pos vi to pos vf
*/
mat4 linear_roll_ctm_interpolation(mat4 ctm, vec4 pi, vec4 pf, float radius, float t, float dt){

	// new ctm = old ctm translated by (pf-pi)*t and rotated by ___

	vec4 rollVec = subtract_v4(pf, pi);
	vec4 drollVec = scale_v4(dt, rollVec);
	vec4 rollHat = normalize_v4(rollVec);
	vec4 yHat = {0.0, 1.0, 0.0, 0.0};
	vec4 aboutHat = cross_v4(yHat, rollHat);

	float d = magnitude_v4(rollVec) * dt;
	float theta = d / radius;	

	mat4 rotated_ctm = arb_rotation_com_world(ctm, theta, aboutHat.x, aboutHat.y, aboutHat.z);
	mat4 translated_rotated_ctm = multiply_m4(translation(drollVec.x, drollVec.y, drollVec.z), rotated_ctm);
	return translated_rotated_ctm;
}

/*
 * return new ctm after a step dt of radial rolling
 */
mat4 radial_roll_step(mat4 ctm, float ballradius, float dt){
	// get position of ball (using ctm)
	vec4 bpos = {0,0,0,1};
	bpos = multiply_m4v4(ctm, bpos);

	float Rbig = sqrt(bpos.x*bpos.x + bpos.z*bpos.z);

	dt *= sin(1.0*M_PI*Rbig/(0.2*16));
	float dTheta = dt / Rbig;

	vec4 newPos = multiply_m4v4(y_rotation(dTheta), bpos);
	vec4 tr = subtract_v4(newPos, bpos);
	
	vec4 aboutHat = {bpos.x, bpos.y, bpos.z, 0};
	aboutHat = normalize_v4(aboutHat);
	vec4 yHat = {0.0, 1.0, 0.0, 0.0};	
	float cosGamma = dot_v4(aboutHat, yHat);
	float rho = ballradius*sin(acos(cosGamma));
	
	float dPhi = magnitude_v4(tr) / rho;

	mat4 rotated_ctm = arb_rotation_com_world(ctm, dPhi, -aboutHat.x, -aboutHat.y, -aboutHat.z);
	mat4 translated_rotated_ctm = multiply_m4(translation(tr.x, tr.y, tr.z), rotated_ctm);
	return translated_rotated_ctm;
}


int count = 0;
int ticks1 = 500;
int ticks2 = 10000;
vec4 prevBallpoint;
void idle(void){
	// follow ball (if applicable)
	if(following != -1){
		// get point of ball
		vec4 ballpoint = {0.0, 0.0, 0.0, 1.0};
		ballpoint = multiply_m4v4(ctms[following-1], ballpoint);
		//update followMotion
		vec4 dmotion = subtract_v4(ballpoint, prevBallpoint);
		if(magnitude_v4(dmotion) > 0.00000001){
			followMotion = normalize_v4(subtract_v4(ballpoint, prevBallpoint));
			prevBallpoint = ballpoint;
			// look at ball	
			atPoint = ballpoint;

			// if allowing orbit around the followed ball
			if(orbitFollow == 1){
				ballpoint.w = 0.0;
				eyePoint = add_v4(relEyePoint, ballpoint);

			}else{ // if following directly behind the followed ball
				eyePoint = add_v4(ballpoint, scale_v4(-1.5, followMotion));
				eyePoint.y = eyePoint.y + 0.3;
			}
			model_view_matrix = look_at(eyePoint, atPoint, upVector);
		}
		
	}
	//update light ctm
	ctms[numBalls+1] = translation(light_position.x, light_position.y, light_position.z);

	if(keysdown[0] == 1){
		vec4 lookVector = normalize_v4(subtract_v4(atPoint, eyePoint));
		//if(dot_v4(lookVector, upVector) > -0.98){
			vec4 aboutVector = normalize_v4(cross_v4(lookVector, upVector));
			relEyePoint = multiply_m4v4(arb_rotation_origin(-0.01, aboutVector.x, aboutVector.y, aboutVector.z), relEyePoint);
			if(following == -1)
				eyePoint = relEyePoint;
			model_view_matrix = look_at(eyePoint, atPoint, upVector);
		//}
	}
	if(keysdown[1] == 1){
		relEyePoint = multiply_m4v4(y_rotation(0.01), relEyePoint);
		if(following == -1)
			eyePoint = relEyePoint;
		model_view_matrix = look_at(eyePoint, atPoint, upVector);
	}
	if(keysdown[2] == 1){
		vec4 lookVector = normalize_v4(subtract_v4(atPoint, eyePoint));
		//if(dot_v4(lookVector, upVector) < 0.0){
			vec4 aboutVector = normalize_v4(cross_v4(lookVector, upVector));
			relEyePoint = multiply_m4v4(arb_rotation_origin( 0.01, aboutVector.x, aboutVector.y, aboutVector.z), relEyePoint);
			if(following == -1)
				eyePoint = relEyePoint;
			model_view_matrix = look_at(eyePoint, atPoint, upVector);
		//}
	}
	if(keysdown[3] == 1){
		relEyePoint = multiply_m4v4(y_rotation(-0.01), relEyePoint);
		if(following == -1)
			eyePoint = relEyePoint;
		model_view_matrix = look_at(eyePoint, atPoint, upVector);
	}
	if(keysdown[4] == 1){
		relEyePoint = scale_v4(1.0/1.01, relEyePoint);
		if(following == -1)
			eyePoint = relEyePoint;
		model_view_matrix = look_at(eyePoint, atPoint, upVector);
	}
	if(keysdown[5] == 1){
		relEyePoint = scale_v4(1.01, relEyePoint);
		if(following == -1)
			eyePoint = relEyePoint;
		model_view_matrix = look_at(eyePoint, atPoint, upVector);
	}

	if(mode == 0){
	}
	else if(mode == 1){
		float t = (float)count/(float)ticks1;
		float dt = 1.0/(float)ticks1;

		for(int i=0; i<numBalls; i++){
			vec4 targetpos = {0.0, 0.1, -0.2*(numBalls-i-1), 1.0};
			ctms[i] = linear_roll_ctm_interpolation(ctms[i], ballpositions[i], targetpos, 0.1, t, dt);
		}
		
		if(count>ticks1){
			mode++;
			count = 0;
		}else{
			count++;
		}
	}
	else if(mode == 3){
		float dt = 1.0/(float)500.0;

		for(int i=0; i<numBalls; i++){
			ctms[i] = radial_roll_step(ctms[i], 0.1, dt);
		}
		
		if(count>ticks2){
			mode++;
			count = 0;
		}else{
			count++;
		}
	}

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_size, window_size);
    glutCreateWindow("billiard");
    glewInit();
    init();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboardDown);
    glutKeyboardUpFunc(keyboardUp);
    glutSpecialFunc(keyboardSpecialDown);
    glutSpecialUpFunc(keyboardSpecialUp);
    glutMainLoop();

    return 0;
}
