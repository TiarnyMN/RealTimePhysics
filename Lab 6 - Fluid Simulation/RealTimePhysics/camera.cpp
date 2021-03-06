#include "camera.h"

Camera::Camera(void)
{
	position = glm::vec3(0.0f, 0.0f, -10.0f);
	target = glm::vec3(0.0f, 0.0f, 0.0f);
	forwardVec = glm::normalize(target - position);

	rotVec = glm::vec3(1.57079633f, 0.0f, 0.0f);
	upVec = glm::vec3(0.0f, 1.0f, 0.0f);

	camFOV = 90.0f;
	camNear = 0.1f;
	camFar = 100000.0f;
	camSpin = 0.0f;
	distance = 20.0f;
	verticalOffset = 5.0f;

	speedConstant = 2.0f;

	position.y = target.y + verticalOffset;

	spinEnabled = false;
}

void Camera::RotateX(float deltaTimeStep)
{
	rotVec.x += speedConstant * deltaTimeStep;
	target.x = position.x + cos(rotVec.x) * 10;
	target.z = position.z - sin(rotVec.x) * 10;

	forwardVec = glm::normalize(target - position);
}

void Camera::RotateY(float deltaTimeStep)
{
	rotVec.y += speedConstant * deltaTimeStep;
	target.z = position.z - cos(rotVec.y) * 10;
	target.y = position.y + sin(rotVec.y) * 10;

	forwardVec = glm::normalize(target - position);
}

void Camera::MoveForwards(float deltaTimeStep)
{
	position += forwardVec * speedConstant * 100.0f * deltaTimeStep;
}

void Camera::MoveBackwards(float deltaTimeStep)
{
	position -= forwardVec * speedConstant * 100.0f * deltaTimeStep;
}

void Camera::SetTarget(glm::vec3 newTarget)
{
	target = newTarget;
	camSpin = 0.0f;
	verticalOffset = 10.0f;
	distance = 2.0f;
}

glm::mat4 Camera::GetViewMatrix()
{
	return glm::lookAt(position, target, upVec);
}

glm::mat4 Camera::GetProjMatrix()
{
	return glm::perspective(45.0f, ((float)glutGet(GLUT_WINDOW_WIDTH) / (float)glutGet(GLUT_WINDOW_HEIGHT)), 0.001f, 1000.0f);
}

void Camera::Update(GLfloat gameTime, glm::vec3 newTarget)
{
	if(spinEnabled)
		camSpin += gameTime/2.0f;

	target = newTarget;

	GLfloat camX = target.x + (sin(camSpin) * distance);
	GLfloat camZ = target.z + (cos(camSpin) * distance);
	
	position.x = camX;
	//position.y = target.y + verticalOffset;
	position.z = camZ;
}

void Camera::ToggleSpin()
{
	spinEnabled = !spinEnabled;
}