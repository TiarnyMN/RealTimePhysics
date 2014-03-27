#include "FluidSim.h"
#include <gtx\fast_exponential.hpp>

//Generates a random float between minVal and maxVal
float GetRandomValue(float minVal, float maxVal)
{
	float rndVal = ((float)rand()) / (float)RAND_MAX;
	return (rndVal * (maxVal - minVal)) + minVal;
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
				particle->SetVelocity(0.0f, 0.0f, 0.0f);
				particles.push_back(particle);
			}

	wallEnabled = true;
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
		curParticle->color = glm::vec3(0.0f, 0.0f, 1.0f);

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
			surfaceNormal += (compareParticle->inverseMass * compareParticle->density) * CalculateNormalKernal(r, rLength);
		}

		if(glm::length(surfaceNormal) > 0.5f)
		{
			//Colour the particle to indicate it is on the surface.
			//curParticle->color = glm::vec3(1.0f, 1.0f, 1.0f);

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
		curParticle->forceExternal = glm::vec3(0.0f, -9.81f, 0.0f) * curParticle->density;	//Fix this later.

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
		if(curParticle->position.y < -HorizontalBounds && curParticle->velocity.y < 0.0f)
		{
			curParticle->position.y = -HorizontalBounds;
			curParticle->velocity.y = -curParticle->velocity.y * DampeningStrength;
			curParticle->oldVelocity.y = -curParticle->oldVelocity.y * DampeningStrength;
		}

		if(wallEnabled && curParticle->position.x < -HorizontalBounds)
		{
			curParticle->position.x = -HorizontalBounds;
			curParticle->velocity.x = -curParticle->velocity.x * DampeningStrength;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DampeningStrength;
		}
		else if(curParticle->position.x < -HorizontalBounds * 3.0f)
		{
			curParticle->position.x = -HorizontalBounds * 3.0f;
			curParticle->velocity.x = -curParticle->velocity.x * DampeningStrength;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DampeningStrength;
		}
		else if(curParticle->position.x > HorizontalBounds)
		{
			curParticle->position.x = HorizontalBounds;
			curParticle->velocity.x = -curParticle->velocity.x * DampeningStrength;
			curParticle->oldVelocity.x = -curParticle->oldVelocity.x * DampeningStrength;
		}

		if(curParticle->position.z < -HorizontalBounds)
		{
			curParticle->position.z = -HorizontalBounds;
			curParticle->velocity.z = -curParticle->velocity.z * DampeningStrength;
			curParticle->oldVelocity.z = -curParticle->oldVelocity.z * DampeningStrength;
		}
		else if(curParticle->position.z > HorizontalBounds)
		{
			curParticle->position.z = HorizontalBounds; 
			curParticle->velocity.z = -curParticle->velocity.z * DampeningStrength;
			curParticle->oldVelocity.z = -curParticle->oldVelocity.z * DampeningStrength;
		}
	}
}

vector<Particle*> FluidSim::GetParticles(void)
{
	return particles;
}

void FluidSim::Render(GLuint shaderID)
{
	GLUquadricObj* refSphere = gluNewQuadric();
	gluQuadricTexture(refSphere, TRUE);
	gluQuadricNormals(refSphere, TRUE);

	int transformID = glGetUniformLocation(shaderID, "mTransform");
	int colorID = glGetUniformLocation(shaderID, "vColor");

	for(int i = 0; i < particles.size(); i++)
	{
		glm::mat4 finalTransform = glm::translate(particles[i]->GetPosition());
		glUniformMatrix4fv(transformID, 1, GL_FALSE, &finalTransform[0][0]);
		glUniform3f(colorID, particles[i]->color.x, particles[i]->color.y, particles[i]->color.z);

		float particleRenderSize = glm::pow((3.0f * particles[i]->mass) / (4.0f * glm::pi<float>() * particles[i]->density), 1.0f / 3.0f);	
		//gluSphere(refSphere, particleRenderSize, 5, 5);

		gluSphere(refSphere, particleRenderSize, 2, 2);

	}
}