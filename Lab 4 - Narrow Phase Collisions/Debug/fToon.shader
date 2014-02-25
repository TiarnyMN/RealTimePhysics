#version 330
in vec2 oTexUV;
uniform sampler2D texSampler;
out vec4 fragColor;

void main()
{
	fragColor = texture2D(texSampler, oTexUV.xy);
}