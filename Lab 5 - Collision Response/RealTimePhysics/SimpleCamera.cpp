#include "SimpleCamera.h"

SimpleCamera::SimpleCamera(void)
{
	position = glm::vec3(0.0f, 0.0f, -10.0f);
	target = glm::vec3(0.0f, 0.0f, 0.0f);
	upVec = glm::vec3(0.0f, 1.0f, 0.0f);

	Initialise();
}

SimpleCamera::SimpleCamera(const glm::vec3 &pos, const glm::vec3 &tar, const glm::vec3 &up)
{
    position = pos;
	target = tar;
	upVec = up;

	target = glm::normalize(target);
	upVec = glm::normalize(upVec);

	Initialise();
}

void SimpleCamera::SetNewTarget(glm::vec3 newTarget)
{
	target = newTarget;
}

glm::vec3 SimpleCamera::GetPos(void)
{
	return position;
}

void SimpleCamera::Initialise(void)
{
	spinEnabled = false;

	camFOV = 45.0f;
	camNear = 1.0f;
	camFar = 1000.0f;
	camSpin = 0.0f;

	distance = 10.0f;
	verticalOffset = 5.0f;

	aspectRatio = (float)glutGet(GLUT_SCREEN_WIDTH) / (float)glutGet(GLUT_SCREEN_HEIGHT);

	perspectiveMatrix = glm::perspective(camFOV, aspectRatio, camNear, camFar);
	lookAtMatrix = glm::lookAt(position, target, upVec);
}

glm::mat4 SimpleCamera::GetViewMatrix(void)
{
	return glm::lookAt(position, target, upVec);
}

glm::mat4 SimpleCamera::GetProjMatrix(void)
{
	return perspectiveMatrix;
}

void SimpleCamera::ToggleSpin()
{
	spinEnabled = !spinEnabled;
}

bool SimpleCamera::HandleKeyboard(int curKey)
{
	switch(curKey)
	{
		case ' ':
			ToggleSpin();
			return true;
	}

	return false;
}

void SimpleCamera::Update(float elapsedTimeStep)
{
	if(spinEnabled)
		camSpin += elapsedTimeStep/2.0f;

	GLfloat camX = target.x - (sin(camSpin) * distance);
	GLfloat camZ = target.z - (cos(camSpin) * distance);

	position.x = camX;
	position.y = target.y + verticalOffset;
	position.z = camZ;

	lookAtMatrix = glm::lookAt(position, target, upVec);
}