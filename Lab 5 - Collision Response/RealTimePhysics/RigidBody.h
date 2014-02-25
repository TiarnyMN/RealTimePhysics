#pragma once
#include "shared.h"

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
	glm::vec3 position;			//The position of the body (its origin) in the world.
	glm::vec3 momentum;			//The momentum of the body.

	glm::quat orientation;		//The orientation of the rigid body, as a normalised quaternion.
	glm::vec3 angularMomentum;	//The angular momentum of the body.

	// Secondary Values - Calculated from the primary values.
	glm::vec3 velocity;			//The velocity of the body, calculated from its momentum.

	glm::quat spin;				//The speed at which the body's orientation changes, as a quaternion.
	glm::vec3 angularVelocity;	//The angular velocity of the body, calculated from its angular momentum.

	glm::mat4 worldMat;				//Converts body coordinates to world coordinates

	// Fixed Values
	float mass;						//The mass of the rigid body.
	float inverseMass;				// 1 / mass

	glm::mat3 inertialTensor;			//The interia tensor of the cube, representing how it reacts to forces.
	glm::mat3 inverseInertialTensor;	//The inverse tensor, which converts the angular momentum to an angular velocity.

	void RecalculateValues(void);
};

class RigidBody
{
public:
	RigidBody(void);
	RigidBody(float mass, glm::vec3 pos, glm::vec3 mom, glm::quat orient, glm::vec3 angMom, glm::mat3 inTensor);

	void Update(float t, float dt);	
	BodyState* GetCurrentState(void);
	BodyState GetInterpolatedState(const float alpha);

	bool HasGravity(void) const;
	bool HasDamping(void) const;

	void SetGravity(bool enabled);
	void SetDamping(bool enabled);

private:
	BodyState curState, prevState;

	//Our static classes, for all Rigid Bodies
	static BodyState CreateInitialState(const float mass, const glm::vec3 position, const glm::vec3 momentum, const glm::quat orientation, const glm::vec3 angularMomentum, const glm::mat3 inertialTensor);
	static BodyDerivative Evaluate(const RigidBody* rigidBody, BodyState state, float t, float dt, const BodyDerivative &derivative);
	static void UpdateForces(const RigidBody* rigidBody, const BodyState &state, glm::vec3 &force, glm::vec3 &torque);
	static void Integrate(const RigidBody* rigidBody, BodyState &state, float t, float dt);

	bool gravity;
	bool damping;

	bool fixedObject;
};

