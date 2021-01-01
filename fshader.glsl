#version 120

varying vec3 N, L, V;
varying vec4 color;
varying float distance;
varying vec2 texCoord;

uniform sampler2D texture;
uniform float shininess;
uniform float attenuation_constant;
uniform float attenuation_linear;
uniform float attenuation_quadratic;

uniform int use_texture;
uniform int is_shadow;

vec4 ambient, diffuse, specular;

void main()
{
	vec4 the_color = color;
	if(is_shadow == 0){
		if(use_texture == 1){
			the_color = texture2D(texture, texCoord);
		}
		vec3 NN = normalize(N);
		vec3 VV = normalize(V);
		vec3 LL = normalize(L);

		ambient = the_color * 0.2;
		vec3 H = normalize(LL + VV);
		diffuse = max(dot(LL, NN), 0.0) * 0.8 * the_color;
		specular = pow(max(dot(NN, H), 0.0), shininess) * vec4(1.0, 1.0, 1.0, 1.0);

		float attenuation = 1/(attenuation_constant + (attenuation_linear * distance) + (attenuation_quadratic * distance * distance));

		gl_FragColor = ambient + attenuation * (diffuse + specular);

		// USEFUL CALLS FOR DEBUGGING:

		// brightness as inverse of distance from point to light
		//gl_FragColor = vec4(1.0/length(L), 1.0/length(L), 1.0/length(L), 1.0);
		
		// brightness as inverse of distance from point to viewer
		//gl_FragColor = vec4(1.0/length(V), 1.0/length(V), 1.0/length(V), 1.0);

		// color as normal (just to check normals)
		//gl_FragColor = vec4(N, 1.0);
	}else{
		gl_FragColor = the_color;
	}
}
