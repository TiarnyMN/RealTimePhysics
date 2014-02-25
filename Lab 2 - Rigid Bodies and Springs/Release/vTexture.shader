#version 330
layout (location = 0) in vec3 Position;
layout (location = 1) in vec2 TexUV;
uniform mat4 gWVP;
out vec2 oTexUV;

void main()
{
	gl_Position = gWVP * vec4(Position, 1.0);
	oTexUV = TexUV;
}