#include "plane.h"

void Plane::GenerateVertices(glm::vec3 pos, glm::vec3 norm)
{
	position = pos;
	normal = glm::normalize(norm);

	glm::vec3 firstOrtho = glm::vec3(0, 0, 1);
	if(glm::dot(firstOrtho, normal))
		firstOrtho = glm::vec3(1.0f, 0.0f, 0.0f);

	firstOrtho = glm::normalize(glm::cross(normal, firstOrtho));
	glm::vec3 secondOrtho = glm::normalize(glm::cross(normal, firstOrtho));

	vOne = position + firstOrtho * 100.0f + secondOrtho * 100.0f;
	vTwo = position + firstOrtho * 100.0f - secondOrtho * 100.0f;
	vThree = position - firstOrtho * 100.0f - secondOrtho * 100.0f;
	vFour = position - firstOrtho * 100.0f + secondOrtho * 100.0f;
}

void Plane::Render(void)
{
	//glNormal3d(0.0f, 1.0f, 0.0f);
	glNormal3d(normal.x, normal.y, normal.z);
	glVertex4fv(glm::value_ptr(glm::vec4(vOne, 1)));
	glVertex4fv(glm::value_ptr(glm::vec4(vTwo, 1)));
	glVertex4fv(glm::value_ptr(glm::vec4(vThree, 1)));
	glVertex4fv(glm::value_ptr(glm::vec4(vFour, 1)));
}