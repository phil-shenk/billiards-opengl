#version 120

attribute vec4 vPosition;
attribute vec4 vColor;
attribute vec4 vNormal;
attribute vec2 vTexCoord;

varying vec2 texCoord;
varying vec4 color;
varying vec4 vp;
varying vec4 normal;

uniform mat4 model_view_matrix;
uniform mat4 projection_matrix;
uniform mat4 ctm;

uniform vec4 light_position;

uniform int is_shadow;

varying vec3 N;
varying vec3 L;
varying vec3 L_temp;
varying vec3 V;
varying float distance;


void main()
{
	if(is_shadow == 0){
		gl_Position = projection_matrix * model_view_matrix * ctm * vPosition;

		// N = normal vector at p
		N = (model_view_matrix * ctm * vNormal).xyz;

		// L = direction of a line from point to light source
		L = (model_view_matrix * (light_position - (ctm * vPosition))).xyz; // tan has this as L_temp

		// V is direction from point to the viewer (COP)
		V = -1 * (model_view_matrix * ctm * vPosition).xyz;

		color = vColor;
		texCoord = vTexCoord;

		distance = length(L);
	}else{
		
		vp = ctm * vPosition;
		vec4 lp = light_position;
		float x = lp.x - lp.y*(lp.x - vp.x)/(lp.y - vp.y);
		float z = lp.z - lp.y*(lp.z - vp.z)/(lp.y - vp.y);
		
		gl_Position = projection_matrix * model_view_matrix * vec4(x, 0.001, z, 1.0);
		color = vec4(0, 0, 0, 1);
	}
		
}
