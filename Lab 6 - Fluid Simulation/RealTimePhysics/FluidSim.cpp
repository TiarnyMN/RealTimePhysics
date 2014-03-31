#include "FluidSim.h"

//Generates a random float between minVal and maxVal
float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
}

bool FluidSim::HandleRegularInput(unsigned char curKey)
{
	switch(curKey)
	{
		case 'w':
			containerVelocity += glm::vec3(0.0f, 0.0f, 0.05f);
			containerVelocity.z = glm::min(0.25f, containerVelocity.z);
			return true;

		case 's':
			containerVelocity -= glm::vec3(0.0f, 0.0f, 0.05f);
			containerVelocity.z = glm::max(-0.25f, containerVelocity.z);
			return true;

		case 'a':
			containerVelocity += glm::vec3(0.05f, 0.0f, 0.0f);
			containerVelocity.x = glm::min(0.25f, containerVelocity.x);
			return true;

		case 'd':
			containerVelocity -= glm::vec3(0.05f, 0.0f, 0.0f);
			containerVelocity.x = glm::max(-0.25f, containerVelocity.x);
			return true;

		case 'i':
			container.CalculateCapsulePosition(glm::vec3(0.1f, 0.0f, 0.0f));
			innerContainer.CalculateCapsulePosition(glm::vec3(0.1f, 0.0f, 0.0f));
			return true;

		case 'k':
			container.CalculateCapsulePosition(glm::vec3(-0.1f, 0.0f, 0.0f));
			innerContainer.CalculateCapsulePosition(glm::vec3(-0.1f, 0.0f, 0.0f));
			return true;

		case ' ':
			if(!crashCapsule && !brokenCapsule)
			{
				crashCapsule = true;
				brokenCapsule = false;
				//container.centre.y = container.length * 2.0f;
				//containerVelocity.y = 0.5f;
				//innerContainer.centre.y = container.length * 2.0f;
			}
			return true;

		case '0':
			showSurfaceParticles = false;
			showNeighbours = false;
			hideSurfaceParticles = false;	
			return true;

		case '1':
			showSurfaceParticles = true;
			showNeighbours = false;
			hideSurfaceParticles = false;	
			return true;

		case '2':
			hideSurfaceParticles = true;
			showSurfaceParticles = false;
			showNeighbours = false;
			return true;

		case '3':
			showNeighbours = true;
			showSurfaceParticles = false;
			hideSurfaceParticles = false;			
			return true;		

		case '4':
			innerRadius = glm::max(innerRadius - 0.1f, 0.0f);
			innerContainer.radius = innerRadius;
			innerContainer.CalculateCapsulePosition(glm::vec3());
			return true;

		case '5':
			innerRadius = glm::min(innerRadius + 0.1f, HorizontalBounds);
			innerContainer.radius = innerRadius;
			innerContainer.CalculateCapsulePosition(glm::vec3());
			return true;

		case '6':
			innerLength = glm::max(innerLength - 0.1f, 0.0f);
			innerContainer.length = innerLength;
			innerContainer.CalculateCapsulePosition(glm::vec3());
			return true;

		case '7':
			innerLength = glm::min(innerLength + 0.1f, HorizontalBounds);
			innerContainer.length = innerLength;
			innerContainer.CalculateCapsulePosition(glm::vec3());
			return true;

		case '-':
			if(showNeighbours)
				particleHighlightID = glm::max(particleHighlightID - 1, 0);
			return true;

		case '=':
			if(showNeighbours)
				particleHighlightID = glm::min(particleHighlightID + 1, (int)(particles.size() - 1));
			return true;
	}

	return false;
}

void FluidSim::InitialiseParticles(void)
{
	KernelPow9 = glm::pow(SmoothingKernelSize, 9.0f);
	KernelPow6 = glm::pow(SmoothingKernelSize, 6.0f);
	KernelPow2 = glm::pow(SmoothingKernelSize, 2.0f);

	float pi = glm::pi<float>();

	DefaultKernelConstant = 315.0f / (64.0f * pi * KernelPow9);
	NormalKernelConstant = -945.0f / (32.0f * pi * KernelPow9);
	SurfaceKernelConstant = NormalKernelConstant;
	PressureKernelConstant = -45.0f  / (pi * KernelPow6);
	ViscosityKernelConstant = 45.0f / (pi * KernelPow6);

	particles.reserve(ParticleCount);

	for(int x = 0; x < ParticlesPerDimension; x++)
		for(int y = 0; y < ParticlesPerDimension; y++)
			for(int z = 0; z < ParticlesPerDimension; z++)
			{
				Particle* particle = new Particle(ParticleMass);
				particle->SetPosition(GetRandomValue(-HorizontalBounds, HorizontalBounds), GetRandomValue(0.0f, HorizontalBounds), GetRandomValue(-HorizontalBounds, HorizontalBounds));
				//particle->SetPosition(0.0f, 0.0f, 0.0f);
				particle->SetVelocity(0.0f, 0.0f, 0.0f);
				particles.push_back(particle);
			}

	container = Capsule();
	container.radius = HorizontalBounds;
	container.length = HorizontalBounds;
	container.centre = glm::vec3();

	innerContainer = Capsule();

	innerRadius = HorizontalBounds / 2.0f;
	innerLength = HorizontalBounds / 2.0f;

	innerContainer.radius = innerRadius;
	innerContainer.length = innerLength;
	innerContainer.centre = glm::vec3();

	container.CalculateCapsulePosition(glm::vec3());
	innerContainer.CalculateCapsulePosition(glm::vec3());

	crashCapsule = brokenCapsule = showNeighbours = showSurfaceParticles = hideSurfaceParticles = false;
	particleHighlightID = 0;
	wallEnabled = true;
}

//Why is this broken...!?
inline bool FluidSim::CalculateInternalCapsuleCollision(Capsule* container, Particle* particle)
{
	float t = -glm::dot((container->startPoint - particle->position), (container->endPoint - container->startPoint) / glm::distance2(container->endPoint, container->startPoint));
	t = glm::min(1.0f, glm::max(0.0f, t));

	container->q = container->startPoint + (container->endPoint - container->startPoint) * t;
	
	glm::vec3 distanceFromContainerToParticle = particle->position - container->q;
	float distanceMagnitude = glm::length(distanceFromContainerToParticle);

	if(distanceMagnitude >= container->radius)
	{
		//Particle is either touching, or outside, the container.
		glm::vec3 contactPoint = container->q + container->radius * (distanceFromContainerToParticle / distanceMagnitude);	

		float penetrationDepth = glm::length(CalculateCapsuleFunc(container, particle));
		glm::vec3 unitSurfaceNormal = glm::sign(CalculateCapsuleFunc(container, particle)) * ((container->q - particle->position) / glm::length(container->q - particle->position));

		//particle->position = contactPoint;	//We set the particle to be touching the inner surface, if it was outside.

		glm::vec3 relPos = particle->position - contactPoint;
		float dotProd = glm::dot(relPos, unitSurfaceNormal);		//(x-p).n
		float dotProd2 = glm::dot(unitSurfaceNormal, particle->velocity);		//(n.v)

		float epsilon = 0.1f;	//A collision smoothing / buffer value, how close a particle can actually get to a plane before we consider it a collision.

		if(dotProd < epsilon && dotProd2 < epsilon)
		{
			glm::vec3 vN = glm::dot(unitSurfaceNormal, particle->velocity) * (container->q - particle->position);
			glm::vec3 vT = particle->velocity - vN;
			vT *= FrictionStrength;	//"Frictional" force

			particle->velocity = vT - vN * DampeningStrength;
			particle->oldVelocity = particle->velocity;
		}

		particle->position = contactPoint;	//We set the particle to be touching the inner surface, if it was outside.
		return true;
	}

	return false;
}

inline bool FluidSim::CalculateExternalCapsuleCollision(Capsule* container, Particle* particle)
{
	float t = -glm::dot((container->startPoint - particle->position), (container->endPoint - container->startPoint) / glm::distance2(container->endPoint, container->startPoint));
	t = glm::min(1.0f, glm::max(0.0f, t));

	container->q = container->startPoint + (container->endPoint - container->startPoint) * t;
	
	glm::vec3 distanceFromContainerToParticle = particle->position - container->q;
	float distanceMagnitude = glm::length(distanceFromContainerToParticle);

	if(distanceMagnitude <= container->radius)
	{
		//Particle is either touching, or outside, the container.
		glm::vec3 contactPoint = container->q + container->radius * (distanceFromContainerToParticle / distanceMagnitude);	

		float penetrationDepth = glm::length(CalculateCapsuleFunc(container, particle));
		glm::vec3 unitSurfaceNormal = glm::sign(CalculateCapsuleFunc(container, particle)) * ((particle->position - container->q) / glm::length(particle->position - container->q));

		glm::vec3 relPos = particle->position - contactPoint;
		float dotProd = glm::dot(relPos, unitSurfaceNormal);		//(x-p).n
		float dotProd2 = glm::dot(unitSurfaceNormal, particle->velocity);		//(n.v)

		float epsilon = 0.1f;	//A collision smoothing / buffer value, how close a particle can actually get to a plane before we consider it a collision.

		if(dotProd < epsilon && dotProd2 < epsilon)
		{
			glm::vec3 vN = glm::dot(unitSurfaceNormal, particle->velocity) * (particle->position - container->q);
			glm::vec3 vT = particle->velocity - vN;
			vT *= FrictionStrength;	//"Frictional" force

			particle->velocity = vT - vN * DampeningStrength;
			particle->oldVelocity = particle->velocity;
		}

		particle->position = contactPoint;	//We set the particle to be touching the inner surface, if it was outside.
		return true;
	}

	return false;
}

inline float FluidSim::CalculateCapsulePoint(Capsule* container, Particle* particle)
{
	float t = -glm::dot((container->startPoint - particle->position), (container->endPoint - container->startPoint) / glm::distance2(container->endPoint, container->startPoint));
	return glm::min(1.0f, glm::max(0.0f, t));
}

inline float FluidSim::CalculateCapsuleFunc(Capsule* container, Particle* particle)
{
	return glm::length(container->q - particle->position) - container->radius;
}

inline float FluidSim::CalculateKernal(Particle* curParticle, Particle* compareParticle)
{
	return DefaultKernelConstant * glm::pow(KernelPow2 - glm::distance2(curParticle->position, compareParticle->position), 3.0f);
}

inline glm::vec3 FluidSim::CalculateNormalKernal(glm::vec3 &r, float &rLength)
{
	if(rLength == 0.0f)
		return glm::vec3();

	return NormalKernelConstant * (r / rLength) * glm::pow(KernelPow2 - glm::pow(rLength, 2.0f), 2.0f);
}

inline float FluidSim::CalculateSurfaceKernal(glm::vec3 &r, float &rLength)
{
	if(rLength == 0.0f)
		return 0.0f;

	float rLengthSquared =  glm::pow(rLength, 2.0f);
	return SurfaceKernelConstant * (KernelPow2 - rLengthSquared) * (3.0f * KernelPow2 - 7.0f * rLengthSquared);
}

inline glm::vec3 FluidSim::CalculatePressureKernal(glm::vec3 &r, float &rLength)
{
	if(glm::length(r) == 0.0f)
		return glm::vec3();

	return PressureKernelConstant * ((r / rLength) * glm::pow((SmoothingKernelSize - rLength), 2.0f));
}

inline float FluidSim::CalculateViscosityKernal(glm::vec3 &r, float &rLength)
{
	return	ViscosityKernelConstant * (SmoothingKernelSize - rLength);
}

void FluidSim::Update(void)
{
	hashMap.clear();

	if(!crashCapsule && ! brokenCapsule)
	{
		containerVelocity *= 0.95f;
		container.centre += containerVelocity;
		container.CalculateCapsulePosition(glm::vec3());
		innerContainer.centre += containerVelocity;
		innerContainer.CalculateCapsulePosition(glm::vec3());
	}

	if(crashCapsule)
	{
		containerVelocity *= 0.95f;
		containerVelocity -= glm::vec3(0.0f, 0.05f, 0.0f);
		container.centre += containerVelocity;
		container.CalculateCapsulePosition(glm::vec3());
		innerContainer.centre += containerVelocity;
		innerContainer.CalculateCapsulePosition(glm::vec3());

		if((container.centre - container.length).y <= -20.0f)
		{
			crashCapsule = false;
			brokenCapsule = true;
		}
	}
  
	for(int i = 0; i < ParticleCount; i++)
	{
		hashMap.insert(std::pair<glm::vec3, Particle*>(glm::floor(particles[i]->position / SmoothingKernelSize), particles[i]));
		Particle* curParticle = particles[i];
		curParticle->neighbours.clear();

		for(int x = -1; x < 2; x++)
			for(int y = -1; y < 2; y++)
				for(int z = -1; z < 2; z++)
				{
					auto iterator = hashMap.equal_range(glm::vec3(x, y, z) + glm::floor(particles[i]->position / SmoothingKernelSize));
					for(auto it = iterator.first; it != iterator.second; it++)
					{
						Particle* compParticle = (*it).second;

						if(glm::distance(curParticle->position, compParticle->position) <= SmoothingKernelSize)
						{
							curParticle->neighbours.push_back((*it).second);

							if(compParticle != curParticle)
								compParticle->neighbours.push_back(curParticle);
						}
					}
				}
	}

	for(int i = 0; i < ParticleCount; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->density = 0.0f;
		curParticle->pressure = 0.0f;

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];
			curParticle->density += compareParticle->mass * CalculateKernal(curParticle, compareParticle);
		}

		curParticle->pressure = StiffnessConstant * (curParticle->density);	//Compute pressure p(i).
	}

	//cout << endl;

	//5.6.3 (ii) Compute the pressure force density
	for(int i = 0; i < ParticleCount; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->forcePressure = glm::vec3();
		curParticle->forceViscosity = glm::vec3();
		curParticle->forceSurface = glm::vec3();
		curParticle->color = Color;

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];

			if(curParticle == compareParticle)
				continue;

			glm::vec3 r = curParticle->position - compareParticle->position;
			float rLength = glm::length(r);

			float firstPart = (curParticle->mass / curParticle->density) * ((curParticle->pressure + compareParticle->pressure) * 0.5f);
			curParticle->forcePressure -= firstPart * CalculatePressureKernal(r, rLength);
			
			glm::vec3 firstHalf = (compareParticle->velocity - curParticle->velocity) * (compareParticle->mass / compareParticle->density);
			curParticle->forceViscosity += firstHalf * CalculateViscosityKernal(r, rLength);
		}

		curParticle->forceViscosity *= ViscosityConstant;
		curParticle->forceInternal = curParticle->forcePressure + curParticle->forceViscosity;

		glm::vec3 surfaceNormal = glm::vec3();	//Inward Surface Normal

		for(int j = 0; j < curParticle->neighbours.size(); j++)
		{
			Particle* compareParticle = curParticle->neighbours[j];

			glm::vec3 r = curParticle->position - compareParticle->position;
			float rLength = glm::length(r);
			surfaceNormal += (compareParticle->mass / compareParticle->density) * CalculateNormalKernal(r, rLength);
		}

		if(glm::length(surfaceNormal) > 0.6f)
		{
			//Colour the particle to indicate it is on the surface.
			if(showSurfaceParticles || hideSurfaceParticles)
				curParticle->color = glm::vec3(1.0f, 1.0f, 1.0f);

			float laplacianC = 0.0f;
			for(int j = 0; j < curParticle->neighbours.size(); j++)
			{
				Particle* compareParticle = curParticle->neighbours[j];
				glm::vec3 r = curParticle->position - compareParticle->position;
				float rLength =  glm::length(r);
				laplacianC += CalculateSurfaceKernal(r, rLength);
			}

			curParticle->forceSurface = -Sigma * laplacianC * glm::normalize(surfaceNormal);
		}
		
	}

	//5.6.4 (i) + 5.6.5 (ii) Add gravitational force and update particle acceleration. 
	for(int i = 0; i < ParticleCount; i++)
	{
		Particle* curParticle = particles[i];
		curParticle->forceExternal = glm::vec3(0.0f, -9.81f, 0.0f) * curParticle->density;

		curParticle->acceleration = (curParticle->forceInternal + curParticle->forceExternal + curParticle->forceSurface) / curParticle->density;
	}

	//5.6.5 (iiu) Use leapfrog integration to update particle position.
	for(int i = 0; i < ParticleCount; i++)
	{
		Particle* curParticle = particles[i];
		glm::vec3 newVel = curParticle->oldVelocity + FluidStep * curParticle->acceleration;
		//glm::vec3 halfVel = (curParticle->oldVelocity + newVel) / 2.0f;

		glm::vec3 oldPos = curParticle->position;
		curParticle->position = curParticle->position + newVel * FluidStep;

		curParticle->velocity = newVel + (FluidStep / 2.0f) * curParticle->acceleration;
		curParticle->oldVelocity = newVel;
	}

	

	for(int i = 0; i < ParticleCount; i++)
	{
		Particle* curParticle = particles[i];
		if(!brokenCapsule)
		{
			CalculateInternalCapsuleCollision(&container, curParticle);
			CalculateExternalCapsuleCollision(&innerContainer, curParticle);
		}
		else
		{
			if(curParticle->position.y < -20.0f && curParticle->velocity.y < 0.0f)
			{
				curParticle->position.y = -20.0f;
				curParticle->velocity.y = -curParticle->velocity.y * DampeningStrength;
				curParticle->oldVelocity.y = -curParticle->oldVelocity.y * DampeningStrength;
			}
		}
	}
}

vector<Particle*> FluidSim::GetParticles(void)
{
	return particles;
}

void FluidSim::PrintDebug(std::stringstream& ss)
{
	ss << "Particle Count: " << ParticleCount << endl;
	ss << "Particle Mass: " << ParticleMass << endl << endl;
	ss << "Stiffness: " << StiffnessConstant << endl;
	ss << "Viscosity: " << ViscosityConstant << endl;
	ss << "Particles Per Kernel: " << AverageParticlesPerKernel << endl;
	ss << "Smoothing Kernel Size: " << SmoothingKernelSize << endl;

	ss << "Dampening: " << DampeningStrength << endl;
	ss << "Friction: " << FrictionStrength << endl;

	if(!brokenCapsule)
	{
		ss << endl << "INNER CAPSULE" << endl;
		ss << "Length: " << innerContainer.length * 2.0f << endl;
		ss << "Radius: " << innerContainer.radius << endl;
	}
}

void FluidSim::Render(GLuint capsuleShaderID, GLuint particleShaderID, glm::mat4 &viewMat, glm::mat4 &projMat)
{
	glm::mat4 finalTransform;
	int transformID, colorID;

	GLUquadricObj* refSphere = gluNewQuadric();
	gluQuadricTexture(refSphere, TRUE);
	gluQuadricNormals(refSphere, TRUE);

	if(!brokenCapsule)
	{
		GLUquadricObj* refCylinder = gluNewQuadric();
		gluQuadricTexture(refSphere, TRUE);
		gluQuadricNormals(refSphere, TRUE);

		glUseProgram(capsuleShaderID);
		glUniformMatrix4fv(glGetUniformLocation(capsuleShaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
		glUniformMatrix4fv(glGetUniformLocation(capsuleShaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));

		transformID = glGetUniformLocation(capsuleShaderID, "mTransform");
		colorID = glGetUniformLocation(capsuleShaderID, "vColor");

		float initialRotationOffset = glm::pi<float>() / 2.0f;

		glDepthFunc(GL_LEQUAL);
		glDepthMask(GL_FALSE);
		glUniform3f(colorID, 0.3f, 0.3f, 0.3f);	
		finalTransform = glm::translate(container.centre) * glm::eulerAngleYXZ(initialRotationOffset + container.rotation.y, initialRotationOffset + container.rotation.x, container.rotation.z) * glm::translate(0.0f, 0.0f, container.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluSphere(refSphere, container.radius * 1.2f, 16, 16);

		finalTransform = glm::translate(container.centre) * glm::eulerAngleYXZ(initialRotationOffset + container.rotation.y, initialRotationOffset + container.rotation.x, container.rotation.z) * glm::translate(0.0f, 0.0f, -container.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluSphere(refSphere, container.radius * 1.2f, 16, 16);

		finalTransform = glm::translate(container.centre) * glm::eulerAngleYXZ(initialRotationOffset + container.rotation.y, initialRotationOffset + container.rotation.x, container.rotation.z) * glm::translate(0.0f, 0.0f, -container.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluCylinder(refCylinder, container.radius * 1.2f, container.radius * 1.2f, container.length * 2.0, 16, 16);
		glDepthMask(GL_TRUE);

		glUseProgram(particleShaderID);
		glUniformMatrix4fv(glGetUniformLocation(particleShaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
		glUniformMatrix4fv(glGetUniformLocation(particleShaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));

		transformID = glGetUniformLocation(particleShaderID, "mTransform");
		colorID = glGetUniformLocation(particleShaderID, "vColor");

		glUniform3f(colorID, 0.0f, 1.0f, 0.0f);
		finalTransform = glm::translate(innerContainer.centre) * glm::eulerAngleYXZ(initialRotationOffset + innerContainer.rotation.y, initialRotationOffset + innerContainer.rotation.x, innerContainer.rotation.z) * glm::translate(0.0f, 0.0f, innerContainer.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluSphere(refSphere, innerContainer.radius, 16, 16);

		finalTransform = glm::translate(innerContainer.centre) * glm::eulerAngleYXZ(initialRotationOffset + innerContainer.rotation.y, initialRotationOffset + innerContainer.rotation.x, innerContainer.rotation.z) * glm::translate(0.0f, 0.0f, -innerContainer.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluSphere(refSphere, innerContainer.radius, 16, 16);

		finalTransform = glm::translate(innerContainer.centre) * glm::eulerAngleYXZ(initialRotationOffset + innerContainer.rotation.y, initialRotationOffset + innerContainer.rotation.x, innerContainer.rotation.z) * glm::translate(0.0f, 0.0f, -innerContainer.length);
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		gluCylinder(refCylinder, innerContainer.radius, innerContainer.radius, innerContainer.length * 2.0, 16, 16);
	}
	else
	{
		glUseProgram(particleShaderID);
		glUniformMatrix4fv(glGetUniformLocation(particleShaderID, "mView"), 1, GL_FALSE, glm::value_ptr(viewMat));
		glUniformMatrix4fv(glGetUniformLocation(particleShaderID, "mProjection"), 1, GL_FALSE, glm::value_ptr(projMat));

		transformID = glGetUniformLocation(particleShaderID, "mTransform");
		colorID = glGetUniformLocation(particleShaderID, "vColor");
	}

	if(showNeighbours)
	{
		Particle* curParticle = particles[particleHighlightID];
		finalTransform = glm::translate(curParticle->GetPosition());
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		glUniform3f(colorID, curParticle->color.x, curParticle->color.y, curParticle->color.z);

		float particleRenderSize = glm::pow((3.0f * curParticle->mass) / (4.0f * glm::pi<float>() * curParticle->density), 1.0f / 3.0f);	
		gluSphere(refSphere, particleRenderSize, 4, 4);

		for(int i = 0; i < curParticle->neighbours.size(); i++)
		{
			Particle* neighbourParticle = curParticle->neighbours[i];

			if(curParticle == neighbourParticle)
				continue;

			finalTransform = glm::translate(neighbourParticle->GetPosition());
			glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
			glUniform3f(colorID, 1.0f, 0.0f, 0.0f);

			float particleRenderSize = glm::pow((3.0f * neighbourParticle->mass) / (4.0f * glm::pi<float>() * neighbourParticle->density), 1.0f / 3.0f);	
			gluSphere(refSphere, particleRenderSize, 4, 4);
		}
	}
	else
	{  
		for(int i = 0; i < particles.size(); i++)
		{
			finalTransform = glm::translate(particles[i]->GetPosition());
			   
			if(hideSurfaceParticles && particles[i]->color == glm::vec3(1.0f, 1.0f, 1.0f))
				continue;

			glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
			glUniform3f(colorID, particles[i]->color.x, particles[i]->color.y, particles[i]->color.z);

			float particleRenderSize = glm::pow((3.0f * particles[i]->mass) / (4.0f * glm::pi<float>() * particles[i]->density), 1.0f / 3.0f);	
			gluSphere(refSphere, particleRenderSize, 4, 4);
		}
	}
}