#version 330
layout (location = 0) in vec3 vPosition;
layout (location = 2) in vec3 vNormal;
layout (location = 8) in vec2 vTexCoord;

uniform mat4 mView, mProjection, mTransform;

out vec3 oNormal;
out vec3 ePos;

void main()
{
	ePos = vec3(mView * mTransform * vec4(vPosition, 1.0));
	gl_Position = mProjection * vec4(ePos, 1.0);

	oNormal = (mView * vec4(vNormal, 0.0)).xyz;
}