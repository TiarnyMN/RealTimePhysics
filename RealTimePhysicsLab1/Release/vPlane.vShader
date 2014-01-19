#version 330
layout (location = 0) in vec3 vPosition;
layout (location = 2) in vec3 vNormal;
layout(location = 3) in vec4 vColor;

uniform mat4 camMat;
uniform mat4 perMat;

out vec4 oColor;
out vec3 oNormal;
out vec3 ePos;

void main()
{
	oColor = vColor;
	ePos = vec3(camMat * vec4(vPosition, 1.0));
	gl_Position = perMat * vec4(ePos, 1.0);

	oNormal = (camMat * vec4(vNormal, 0.0)).xyz;
}