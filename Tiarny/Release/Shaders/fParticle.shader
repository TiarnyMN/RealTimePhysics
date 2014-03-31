#version 330
in vec4 fragColor;
out vec4 finalColor; 

void main()
{	
	if(fragColor.a <= 0.1f)
		discard;

	finalColor = fragColor;
}