#version 330
in vec4 oColor;
in vec3 oNormal;
in vec3 ePos;

uniform mat4 camMat;
uniform mat4 perMat;

out vec4 fragColor;

// fixed point light properties
vec3 light_position_world  = vec3 (0.0, 0.0, 2.0);
vec3 Ls = vec3 (1.0, 1.0, 1.0); // white specular colour
vec3 Ld = vec3 (0.7, 0.7, 0.7); // dull white diffuse light colour
vec3 La = vec3 (0.2, 0.2, 0.2); // grey ambient colour
  
// surface reflectance
vec3 Ks = vec3 (1.0, 1.0, 1.0); // fully reflect specular light
vec3 Kd = vec3 (1.0, 1.0, 1.0); // orange diffuse surface reflectance
vec3 Ka = vec3 (1.0, 1.0, 1.0); // fully reflect ambient light
float specular_exponent = 100.0; // specular 'power'

void main()
{
	// ambient intensity
	vec3 Ia = La * Ka;

	// diffuse intensity
	vec3 light_position_eye = vec3 (camMat * vec4 (light_position_world, 1.0));
	vec3 distance_to_light_eye = light_position_eye - ePos;
	vec3 direction_to_light_eye = normalize (distance_to_light_eye);
	float dot_prod = dot (direction_to_light_eye, oNormal);
	dot_prod = max (dot_prod, 0.0);
	vec3 Id = Ld * Kd * dot_prod; // final diffuse intensity
	
	fragColor = oColor * vec4(Id + Ia, 1.0);
}