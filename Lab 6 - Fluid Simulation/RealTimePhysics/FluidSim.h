#pragma once
#include "shared.h";
#include <map>
#include <unordered_map>
#include "Hasher.h"

class FluidSim
{
public:
	virtual void Update(void);
	virtual void Render(GLuint shaderID);

	vector<Particle*> GetParticles(void);

private:
	float CalculateKernal(Particle* curParticle, Particle* compareParticle);
	glm::vec3 CalculateNormalKernal(glm::vec3 &r, float &rLength);
	float CalculateSurfaceKernal(glm::vec3 &r, float &rLength);
	glm::vec3 CalculatePressureKernal(glm::vec3 &r, float &rLength);
	float CalculateViscosityKernal(glm::vec3 &r, float &rLength);

	//std::multimap<glm::vec3, Particle*, Hasher> hashMap;
	std::unordered_multimap<glm::vec3, Particle*, Hasher, Hasher> hashMap;

	vector<Particle*> particles;

protected:
	void InitialiseParticles(void);

	float RestDensity;
	float ParticlesPerDimension;
	float ParticleCount;
	float SurfaceThreshold;

	float StiffnessConstant;
	float ViscosityConstant;
	float AverageParticlesPerKernel;
	float FluidStep;

	float VerticalBounds;
	float HorizontalBounds;

	float BoxVolume;

	float DampeningStrength;

	float ParticleMass;

	float SmoothingKernelSize;
	float ParticleSize;

	float Sigma;

	bool wallEnabled;

	float KernelPow9, KernelPow2, KernelPow6;
	float DefaultKernelConstant, NormalKernelConstant, SurfaceKernelConstant, PressureKernelConstant, ViscosityKernelConstant;
};

//	#define PARTICLES_PER_DIMENSION 10.0f
//#define PARTICLE_COUNT (PARTICLES_PER_DIMENSION * PARTICLES_PER_DIMENSION * PARTICLES_PER_DIMENSION)
//
//#define STIFFNESS_CONSTANT 100.0f
//#define VISCOSITY_CONSTANT 50.0f
//#define AVERAGE_KERNEL_PARTICLES 20.0f
//#define SURFACE_THRESHOLD (glm::pow(REST_DENSITY / (AVERAGE_KERNEL_PARTICLES * 2.0f), 0.5f))
//
//#define SIGMA 0.01f
//
//#define FLUID_TIME_STEP 0.005f
//
//#define VERTICAL_BOUNDS PARTICLES_PER_DIMENSION * 2.0f
//#define HORIZONTAL_BOUNDS 2.0f
//#define BOX_VOLUME (HORIZONTAL_BOUNDS * HORIZONTAL_BOUNDS * HORIZONTAL_BOUNDS)
//#define DAMPENING_STRENGTH 0.2f
//#define PARTICLE_MASS ((1.0f * REST_DENSITY * BOX_VOLUME) / PARTICLE_COUNT)
//
//#define SMOOTHING_KERNEL_SIZE (glm::pow((3.0f * BOX_VOLUME * AVERAGE_KERNEL_PARTICLES) / (4.0f * glm::pi<float>() * PARTICLE_COUNT), 1.0f/3.0f))
//#define PARTICLE_SIZE (5.0f * glm::pow((3.0f * PARTICLE_MASS) / (4.0f * glm::pi<float>() * REST_DENSITY), 1.0f / 3.0f))

