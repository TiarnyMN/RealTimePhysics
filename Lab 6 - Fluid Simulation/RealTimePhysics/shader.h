#pragma once
#include "shared.h"

#include <iostream>
#include <sstream>
#include <filesystem>

using namespace std;

class Shader
{
	public:
		Shader(void);
		Shader(const string &shaderPath, const string& shaderFragmentPath);
		GLuint GetShaderID();		//Gets the ID for the given shader.
 
	private:
		bool CreateShaderFromFile(const string &shaderPath, const string& shaderFragmentPath);	//Generates, and compiles, a vertex and fragment shader located at the given paths.
		bool LoadShaderFile(const string &shaderPath, string &shaderTextLocation);				//Loads a shader.
		bool AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType);		//Adds a shader.

		GLuint shaderID;				//Stores current shader ID.
		string vertexShaderText;		//Stores current shader vertex text, for processing.
		string fragmentShaderText;		//Stores the current shader fragment text, for processing.
};