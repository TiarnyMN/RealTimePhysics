#version 330
in vec2 oTexUV;
in vec3 oVertexPos, oNormDir;

uniform sampler2D texSampler;

out vec4 fColor; 

void main()
{
	vec3 lightDir = vec3(1.0, 1.0, 0.0);
	float normDotLight = clamp(dot(oNormDir, lightDir), 0, 1);

	vec2 texUV = oTexUV;

	vec3 ambiColor = vec3(0.2, 0.2, 0.2);
	vec3 diffuseColor = vec3(0.0, 0.0, 1.0);

	if(normDotLight > 0.0)
		diffuseColor = diffuseColor * normDotLight;  
	
	fColor = vec4(ambiColor + diffuseColor, 1.0);
}