#version 330
in vec3 vPosition;
layout(location = 3) in vec4 vColor;

uniform mat4 camMat;
uniform mat4 perMat;

out vec4 fragColor;

void main()
{
	fragColor = vColor;
	gl_Position = perMat * camMat * vec4(vPosition, 1.0);
}