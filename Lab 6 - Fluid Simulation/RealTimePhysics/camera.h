#pragma once
#include "shared.h"

class Camera
{
	public:
		Camera(void);				//Constructs the camera with default values.

		void Move(glm::vec3 direction);		//Moves the character in the given direction.
		glm::mat4 GetViewMatrix();		//Generates a view matrix using current parameters.
		glm::mat4 GetProjMatrix();	//Generates a projection matrix for the given resolution.

		void Update(GLfloat gameTime);	//Updates the camera based on elapsed frame time.
		void ToggleSpin();				//Toggles spin on / off.
		void SetTarget(glm::vec3 newTarget);
		void AimAtTarget(glm::vec3 newTarget);

		void RotateX(float deltaTimeStep);
		void RotateY(float deltaTimeStep);
		void MoveForwards(float deltaTimeStep);
		void MoveBackwards(float deltaTimeStep);
	
		void AlignWithTarget(glm::vec3 newTarget);
		glm::vec3 position, target, upVec, forwardVec, rotVec;				//Camera variables to render the scene correctly
		float distance, verticalOffset, speedConstant;

	private:
		GLfloat camFOV, camNear, camFar, camSpin;	//Variables used for the matrices. CamSpin controls the rate at which the camera spins around the object.
		GLboolean spinEnabled;						//Spinning can be enabled / disabled using SPACE.
};