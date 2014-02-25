#pragma once
#include "shared.h"
#include "BasicCube.h"

/*Implementation of the Rigid Body, the Physics Body State, and the Physics Derivatives comes from a mixture of lecture notes,
 as well as this article from Glenn Fiedler - http://gafferongames.com/game-physics/physics-in-3d/*/

struct BodyDerivative
{
	glm::vec3 velocity;	// Derivative of position.
	glm::vec3 force;	// Derivative of momentum.

	glm::quat spin;		// Derivative of orientation.
	glm::vec3 torque;	// Derivative of angular momentum.
};

struct BodyState
{
	// Primary Values
	glm::vec3 position;
	glm::vec3 momentum;

	glm::quat orientation;
	glm::vec3 angularMomentum;

	// Secondary Values
	glm::vec3 velocity;

	glm::quat spin;
	glm::vec3 angularVelocity;

	glm::mat4 worldMat;				//Converts body coordinates to world coordinates
	glm::mat4 inverseWorldMat;		//Converts world coordinates to body coordinates

	// Fixed Values
	float mass;
	float inverseMass;

	glm::mat3 inertialTensor;
	glm::mat3 inverseInertialTensor;

	BodyDerivative evaluate(BodyState initial, float t, float dt, const BodyDerivative &derivative, const glm::vec3 &springPosition, float kVal, float bVal);
	void updateForces(const BodyState &state, float t, glm::vec3 &force, glm::vec3 &torque, const glm::vec3 &springPosition, float kVal, float bVal);
	void integrate(BodyState &state, float t, float dt, const glm::vec3 &springPosition, float kVal, float bVal);
	void recalculate(void);
};

class RigidBody
{
public:
	RigidBody(void);
	
	void Update(float t, float dt, glm::vec3 sPos, float kVal, float bVal);
	void Catchup(const float alpha);
	void Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix);

	BodyState curState;
	BasicCube* cube;

private:
	BodyState prevState, catchupState;
	glm::vec3 springPosition;
};

