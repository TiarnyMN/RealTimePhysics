#version 330
layout (location = 0) in vec3 Position;
layout (location = 3) in vec3 Color;

uniform mat4 gWVP;
out vec3 oColor;

void main()
{
	gl_Position = gWVP * vec4(Position, 1.0);
	oColor = Color;
}