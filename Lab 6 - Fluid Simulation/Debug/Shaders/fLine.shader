#version 330
in vec3 oColor;
out vec4 fragColor;

void main()
{
	fragColor = vec4(oColor, 1.0);
}