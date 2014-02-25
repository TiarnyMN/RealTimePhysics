#include "RigidBody.h"

//Recalculate our secondary values, from our new primary values.
void BodyState::recalculate(void)
{
	velocity = momentum * inverseMass;

	glm::mat3 rotationMat = glm::mat3_cast(orientation);
	glm::mat3 inverseTensor = glm::transpose(rotationMat) * inverseInertialTensor * rotationMat;

	angularVelocity = angularMomentum * inverseTensor;

	orientation = glm::normalize(orientation);
	spin = 0.5f * glm::quat(0.0f, angularVelocity.x, angularVelocity.y, angularVelocity.z) * orientation;

	glm::mat4 translationMat = glm::translate(position);
	
	worldMat = translationMat * glm::mat4_cast(orientation);
	inverseWorldMat = glm::inverse(worldMat);
}

//Evaluate the new state using our derivatives.
BodyDerivative BodyState::evaluate(BodyState initial, float t, float dt, const BodyDerivative &derivative, const glm::vec3 &springPosition, float kVal, float bVal)
{
	initial.position += derivative.velocity * dt;
	initial.momentum += derivative.force * dt;

	initial.orientation = initial.orientation + (derivative.spin * dt);
	initial.angularMomentum += derivative.torque * dt;

	initial.recalculate();

	BodyDerivative updatedOut;
	updatedOut.velocity = initial.velocity;
	updatedOut.spin = initial.spin;

	updateForces(initial, t + dt, updatedOut.force, updatedOut.torque, springPosition, kVal, bVal);

	return updatedOut;
}

void BodyState::updateForces(const BodyState &state, float t, glm::vec3 &force, glm::vec3 &torque, const glm::vec3 &springPosition, float kVal, float bVal)
{
	//Update force and torque.
	if(springPosition.x <= -100.0f)	//For demo purposes, this is how we detect if a rigid body is attached to a spring or not.
		return;

	force = glm::vec3(0.0f, -9.81f, 0.0f);

	glm::vec3 x = glm::vec3();	// Vector displacement of the end of the spring from its equilibrium position.
	x = state.position - springPosition;
	glm::vec3 v;				// Relative velocity between the two points connected by the spring.

	glm::vec3 p = glm::vec3(1.0f, 1.0f, -1.0f);				// Point on the object where the spring is attached.
	p = worldMat * p;

	glm::vec3 com = glm::vec3(0.0f, 0.0f, 0.0f);				// Centre of Mass
	com = worldMat * com;

	p = p - com;

	v = state.velocity + glm::cross(state.angularVelocity, p);	//Calculates the point velocity.

	glm::vec3 springLinearForce = -kVal * (x + p) - bVal * v;
	glm::vec3 springTorqueForce = glm::cross(p, springLinearForce);

	force += springLinearForce;
	torque += springTorqueForce;
}

//Integrate, using RK4.
void BodyState::integrate(BodyState &state, float t, float dt, const glm::vec3 &springPosition, float kVal, float bVal)
{
	BodyDerivative a = evaluate(state, t, 0.0f, BodyDerivative(), springPosition, kVal, bVal);
	BodyDerivative b = evaluate(state, t, dt * 0.5f, a, springPosition, kVal, bVal);
	BodyDerivative c = evaluate(state, t, dt * 0.5f, b, springPosition, kVal, bVal);
	BodyDerivative d = evaluate(state, t, dt, c, springPosition, kVal, bVal);

	state.position += (1.0f / 6.0f) * dt * (a.velocity + 2.0f * (b.velocity + c.velocity) + d.velocity);
	state.momentum += (1.0f / 6.0f) * dt * (a.force + 2.0f * (b.force + c.force) + d.force);
	state.orientation = state.orientation + ((1.0f / 6.0f) * dt * (a.spin + 2.0f * (b.spin + c.spin) + d.spin));
	state.angularMomentum += (1.0f / 6.0f) * dt * (a.torque + 2.0f * (b.torque + c.torque) + d.torque);

	state.recalculate();
}

RigidBody::RigidBody(void)
{
	cube = new BasicCube();

	float height, depth, width, size;
	height = depth = width = size = 2.0f;

	curState = BodyState();
	curState.mass = 1.0f;
	curState.inverseMass = 1.0f / curState.mass;

	curState.position = glm::vec3(0.0f, 0.0f, 0.0f);
	curState.momentum = glm::vec3(0.0f, 0.0f, 0.0f);

	curState.orientation = glm::quat();
	curState.angularMomentum = glm::vec3(0.0f, 0.0f, 0.0f);

	//This is an intertial tensor for a cube, with given values.
	curState.inertialTensor =  glm::mat3(
		(1.0f /12.0f) * curState.mass * (height * height + depth * depth), 0.0f, 0.0f,
		0.0f, (1.0f / 12.0f) * curState.mass * (width * width + depth * depth), 0.0f,
		0.0f, 0.0f, (1.0f / 12.0f) * curState.mass * (width * width + height * height));

	curState.inverseInertialTensor = glm::inverse(curState.inertialTensor);
}

void RigidBody::Update(float t, float dt, glm::vec3 sPos, float kVal, float bVal)
{
	springPosition = sPos;

	prevState = curState;	//Move our curState to the previous state, for interpolation.
	curState.integrate(curState, t, dt, sPos, kVal, bVal);	//Integrate using RK4 for our new state.

	cube->TransformVertices(curState.worldMat);	//Get our vertices, transformed and rotated, so their position is in world space.
}

void RigidBody::Catchup(const float alpha)
{
	//This gives us a state that approximates the physics at the exact current frame time,
	// while still allowing the physics system to update at a fixed rate. 
	// This is used for rendering, then discard next frame by the next physics update step.
	catchupState = curState;

	if(alpha <= 0.0f)	//We are rendering  in sync with an update, so no need to interpolate.
		return;

	catchupState.position = curState.position * alpha + prevState.position * (1.0f - alpha);
	catchupState.momentum = curState.momentum * alpha + prevState.momentum * (1.0f - alpha);
	catchupState.orientation = glm::slerp(curState.orientation, prevState.orientation, alpha);
	catchupState.angularMomentum = curState.angularMomentum * alpha + prevState.angularMomentum * (1.0f - alpha);

	catchupState.recalculate();
}

void RigidBody::Render(glm::mat4 &viewMatrix, glm::mat4 &projMatrix)
{
	glm::mat4 transformationMat = projMatrix * viewMatrix * catchupState.worldMat;	//Use our interpolated world matrix, to keep rendering smooth.
	cube->Render(transformationMat);	//Render our cube.
}