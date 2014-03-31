#version 330
layout (location = 0) in vec3 vPosition;

uniform mat4 mTransform, mView, mProjection;

void main()
{	
	gl_Position = mProjection * mView * mTransform * vec4(vPosition, 1.0);
}