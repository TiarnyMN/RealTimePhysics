#pragma once
#include "shared.h"

class SimpleCamera
{
	public:
		SimpleCamera(void);
		SimpleCamera(const glm::vec3 &pos, const glm::vec3 &tar, const glm::vec3 &up);

		void Update(float elapsedTimeStep);
		bool HandleKeyboard(int curKey);

		glm::mat4 GetViewMatrix(void);	//Generates a view matrix using current parameters.
		glm::mat4 GetProjMatrix(void);	//Gets the current projection matrix.

		glm::vec3 GetPos(void);

	private:
		glm::vec3 position, target, upVec;				//Camera variables to render the scene correctly
		GLfloat distance, verticalOffset;
		float aspectRatio;

		glm::mat4 perspectiveMatrix;
		glm::mat4 lookAtMatrix;
	
		GLfloat camFOV, camNear, camFar, camSpin;	//Variables used for the matrices. CamSpin controls the rate at which the camera spins around the object.
		GLboolean spinEnabled;						//Spinning can be enabled / disabled using SPACE.

		void ToggleSpin();
		void Initialise(void);
};

