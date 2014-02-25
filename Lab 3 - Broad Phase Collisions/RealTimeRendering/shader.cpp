#include "shader.h"

Shader::Shader(void)
{
}

Shader::Shader(const string &shaderVertexPath, const string& shaderFragmentPath)
{
	CreateShaderFromFile(shaderVertexPath, shaderFragmentPath);
}

GLuint Shader::GetShaderID()
{
	return shaderID;
}

bool Shader::CreateShaderFromFile(const string &shaderVertexPath, const string& shaderFragmentPath)
{
	shaderID = glCreateProgram();
    if (shaderID == 0) 
	{
        fprintf(stderr, "ERROR - Unable to create shader program");
        exit(1);
    }

	if(!LoadShaderFile(shaderVertexPath, vertexShaderText))
	{
		fprintf(stderr, "ERROR: Shader file cannot be loaded - ", shaderVertexPath);
        exit(0);
	}

	if(!LoadShaderFile(shaderFragmentPath, fragmentShaderText))
	{
		fprintf(stderr, "ERROR: Shader file cannot be loaded - ", shaderFragmentPath);
        exit(0);
	}

	// Create two shader objects, one for the vertex, and one for the fragment shader
	AddShader(shaderID, vertexShaderText.c_str(), GL_VERTEX_SHADER);
	AddShader(shaderID, fragmentShaderText.c_str(), GL_FRAGMENT_SHADER);

	GLint Success = 0;
    GLchar ErrorLog[1024] = { 0 };
	
    glLinkProgram(shaderID);													// After compiling all shader objects and attaching them to the program, we can finally link it

    glGetProgramiv(shaderID, GL_LINK_STATUS, &Success);							// Check for program related errors using glGetProgramiv
	
	if (Success == 0) 
	{
		glGetProgramInfoLog(shaderID, sizeof(ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "ERROR: Unable to link shader program - ", ErrorLog);
        exit(1);
	}


    glValidateProgram(shaderID);												// Program has been successfully linked but needs to be validated to check whether the program can execute given the current pipeline state
    glGetProgramiv(shaderID, GL_VALIDATE_STATUS, &Success);						// Check for program related errors using glGetProgramiv
    
	if (!Success) 
	{
        glGetProgramInfoLog(shaderID, sizeof(ErrorLog), NULL, ErrorLog);
        fprintf(stderr, "ERROR: Shader program is invalid - ", ErrorLog);
        exit(1);
    }
	
    glUseProgram(shaderID);														// Finally, use the linked shader program
	
	return true;
}

bool Shader::LoadShaderFile(const string &shaderPath, string &shaderTextLocation)
{
	//Load a shader, read the input, and save it to the given location.
	std::ifstream file (shaderPath.c_str());

	if(!file)
		return false;
	
	std :: stringstream stream;
	stream << file.rdbuf();

	file.close();

	shaderTextLocation = stream.str();
	stream.clear();

	return true;
}

bool Shader::AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType)
{
    GLuint ShaderObj = glCreateShader(ShaderType);		//Create a new shader object.

    if (ShaderObj == 0)
	{
        fprintf(stderr, "ERROR: Unable to create shader type - ", ShaderType);
        exit(0);
    }
	
	glShaderSource(ShaderObj, 1, (const GLchar**)&pShaderText, NULL);		// Bind the source code to the shader, this happens before compilation.
    glCompileShader(ShaderObj);												// Compile the shader and check for errors.
   
	GLint success;
    glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);					// Check for shader related errors using glGetShaderiv.

    if (!success) 
	{
        GLchar InfoLog[1024];
        glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
        fprintf(stderr, "ERROR: Unable to compile shader type - ", ShaderType, InfoLog);
        exit(1);
    }
	
    glAttachShader(ShaderProgram, ShaderObj);								// Attach the compiled shader object to the program object.

	return true;
}