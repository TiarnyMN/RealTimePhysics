#include "WaterSim.h"

WaterSim::WaterSim(void)
{
	RestDensity = 1000.0f;
	ParticlesPerDimension = 10.0f;
	ParticleCount = ParticlesPerDimension * ParticlesPerDimension * ParticlesPerDimension;

	StiffnessConstant = 100.0f;
	ViscosityConstant = 50.0f;

	AverageParticlesPerKernel = 20.0f;
	SurfaceThreshold = glm::pow(RestDensity / (AverageParticlesPerKernel * 2.0f), 0.5f);

	Sigma = 0.01f;
	FluidStep = 0.0075f;

	VerticalBounds = ParticlesPerDimension * 2.0f;
	HorizontalBounds = 2.0f;
	BoxVolume = HorizontalBounds * HorizontalBounds * HorizontalBounds;
	
	DampeningStrength = 0.1f;
	FrictionStrength = 0.95f;
	ParticleMass = (1.0f * RestDensity * BoxVolume) / ParticleCount;

	SmoothingKernelSize = glm::pow((3.0f * BoxVolume * AverageParticlesPerKernel) / (4.0f * glm::pi<float>() * ParticleCount), 1.0f/3.0f);
	ParticleSize = 5.0f * glm::pow((3.0f * ParticleMass) / (4.0f * glm::pi<float>() * RestDensity), 1.0f / 3.0f);

	Color = glm::vec3(0.0f, 0.0f, 1.0f);
	InitialiseParticles();
}


WaterSim::~WaterSim(void)
{
}
