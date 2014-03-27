#version 330
layout (location = 0) in vec3 vPosition;
layout (location = 2) in vec3 vNormal;
layout (location = 8) in vec2 vTexCoord;

uniform mat4 mTransform, mView, mProjection;

//out vec2 oTexUV;
//out vec3 oVertexPos, oNormDir;

void main()
{	
	vec3 oVertexPos = vec3(mView * mTransform * vec4(vPosition, 1.0));
	//oNormDir = normalize(mat3(mView * mTransform) * vNormal);

	//oTexUV = vTexCoord;		
	gl_Position = mProjection * vec4(oVertexPos, 1.0);
}