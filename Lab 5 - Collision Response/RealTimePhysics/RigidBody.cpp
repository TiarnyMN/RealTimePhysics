#include "RigidBody.h"

//Recalculate all secondary values (velocity, angularVelocity, normalised orientation, spin etc.) from primary values.
void BodyState::RecalculateValues(void)
{
	velocity = momentum * inverseMass;

	glm::mat3 rotationMat = glm::mat3_cast(orientation);
	glm::mat3 inverseTensor = glm::transpose(rotationMat) * inverseInertialTensor * rotationMat;

	angularVelocity = angularMomentum * inverseTensor;

	orientation = glm::normalize(orientation);
	spin = 0.5f * glm::quat(0.0f, angularVelocity.x, angularVelocity.y, angularVelocity.z) * orientation;

	glm::mat4 translationMat = glm::translate(position);	
	worldMat = translationMat * glm::mat4_cast(orientation);
}

//An empty constructor creates default values which assume it is a RB for a simple cube.
RigidBody::RigidBody(void)
{
	float height, depth, width, size, mass;
	height = depth = width = size = 2.0f;
	mass = 1.0f;

	glm::mat3 approxCubeTensor = glm::mat3(
		(1.0f /12.0f) * mass * (height * height + depth * depth), 0.0f, 0.0f,
		0.0f, (1.0f / 12.0f) * mass * (width * width + depth * depth), 0.0f,
		0.0f, 0.0f, (1.0f / 12.0f) * mass * (width * width + height * height));

	curState = CreateInitialState(mass, glm::vec3(), glm::vec3(), glm::quat(), glm::vec3(), approxCubeTensor);
}

//Takes in values from its parent object, whatever that might be.
RigidBody::RigidBody(const float mass, const glm::vec3 position, const glm::vec3 momentum, const glm::quat orientation, const glm::vec3 angularMomentum, const glm::mat3 inertialTensor)
{
	if(mass == 0.0f)
		fixedObject = true;

	curState = CreateInitialState(mass, position, momentum, orientation, angularMomentum, inertialTensor);
}

BodyState RigidBody::CreateInitialState(const float mass, const glm::vec3 position, const glm::vec3 momentum, const glm::quat orientation, const glm::vec3 angularMomentum, const glm::mat3 inertialTensor)
{
	BodyState initialState = BodyState();

	if(mass == 0.0f)
	{
		initialState.mass = 0.0f;
		initialState.inverseMass = 0.0f;

		initialState.inertialTensor = glm::mat3(0.0f);
		initialState.inverseInertialTensor = glm::mat3(0.0f);
	}
	else
	{
		initialState.mass = mass;
		initialState.inverseMass = 1.0f / initialState.mass;

		initialState.inertialTensor = inertialTensor;
		initialState.inverseInertialTensor = glm::inverse(initialState.inertialTensor);
	}

	initialState.position = position;
	initialState.momentum = momentum;

	initialState.orientation = orientation;
	initialState.angularMomentum = angularMomentum;

	return initialState;
}

BodyState* RigidBody::GetCurrentState(void)
{
	return &curState;
}

bool RigidBody::HasGravity(void) const
{
	return gravity;
}

void RigidBody::SetGravity(bool enabled)
{
	gravity = enabled;
}

bool RigidBody::HasDamping(void) const
{
	return damping;
}

void RigidBody::SetDamping(bool enabled)
{
	damping = enabled;
}

void RigidBody::Update(float t, float dt)
{
	prevState = curState;	//Move our curState to the previous state, for interpolation.
	Integrate(this, curState, t, dt);	//Integrate using RK4 for our new state.
}

//Uses RK4 to integrate the physics state forward by "dt" seconds.
void RigidBody::Integrate(const RigidBody* rigidBody, BodyState &state, float t, float dt)
{
	BodyDerivative a = Evaluate(rigidBody, state, t, 0.0f, BodyDerivative());
	BodyDerivative b = Evaluate(rigidBody, state, t, dt * 0.5f, a);
	BodyDerivative c = Evaluate(rigidBody, state, t, dt * 0.5f, b);
	BodyDerivative d = Evaluate(rigidBody, state, t, dt, c);

	state.position += (1.0f / 6.0f) * dt * (a.velocity + 2.0f * (b.velocity + c.velocity) + d.velocity);
	state.momentum += (1.0f / 6.0f) * dt * (a.force + 2.0f * (b.force + c.force) + d.force);
	state.orientation = state.orientation + ((1.0f / 6.0f) * dt * (a.spin + 2.0f * (b.spin + c.spin) + d.spin));
	state.angularMomentum += (1.0f / 6.0f) * dt * (a.torque + 2.0f * (b.torque + c.torque) + d.torque);

	state.RecalculateValues();
}

// Evaluates all the derivative values for the current physics state at a given time t;
BodyDerivative RigidBody::Evaluate(const RigidBody* rigidBody, BodyState state, float t, float dt, const BodyDerivative &derivative)
{
	//We use a copy of state, not a reference, as we only want to update the "true" state at the end of the RK4 integrator.
	state.position += derivative.velocity * dt;
	state.momentum += derivative.force * dt;

	state.orientation = state.orientation + (derivative.spin * dt);
	state.angularMomentum += derivative.torque * dt;

	state.RecalculateValues();

	BodyDerivative updatedOutput;
	updatedOutput.velocity = state.velocity;
	updatedOutput.spin = state.spin;

	UpdateForces(rigidBody, state, updatedOutput.force, updatedOutput.torque);
	return updatedOutput;
}

void RigidBody::UpdateForces(const RigidBody* rigidBody, const BodyState &state, glm::vec3 &force, glm::vec3 &torque)
{
	//First start the initial force and torque at zero.
	force = glm::vec3();
	torque = glm::vec3();

	if(rigidBody->HasGravity())
		force.y -= 9.81f * state.mass;

	const float linearDamping = 0.25f;
	const float angularDamping = 0.25f;

	if(rigidBody->HasDamping())
	{
		force -= linearDamping * state.velocity;
		torque -= angularDamping * state.angularVelocity;
	}
}

//Interpolates between the previous state and the next state for the current frame, to create a smooth intermediary step.
BodyState RigidBody::GetInterpolatedState(const float alpha)
{
	// This gives us a state that approximates the physics at the exact current frame time,
	// while still allowing the physics system to update at a fixed rate. 
	// This is used for rendering, then discarded next frame by the next physics update step.
	BodyState smoothedState = curState;

	if(alpha <= 0.0f)
		return prevState;

	smoothedState.position = curState.position * alpha + prevState.position * (1.0f - alpha);
	smoothedState.momentum = curState.momentum * alpha + prevState.momentum * (1.0f - alpha);
	smoothedState.orientation = glm::slerp(curState.orientation, prevState.orientation, alpha);
	smoothedState.angularMomentum = curState.angularMomentum * alpha + prevState.angularMomentum * (1.0f - alpha);

	smoothedState.RecalculateValues();

	return smoothedState;
}

