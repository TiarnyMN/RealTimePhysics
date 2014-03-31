#include "MucusSim.h"

MucusSim::MucusSim(void)
{
	RestDensity = 2000.0f;
	ParticlesPerDimension = 10.0f;
	ParticleCount = ParticlesPerDimension * ParticlesPerDimension * ParticlesPerDimension;

	StiffnessConstant = 50.0f;
	ViscosityConstant = 700.0f;

	AverageParticlesPerKernel = 40.0f;
	SurfaceThreshold = glm::pow(RestDensity / (AverageParticlesPerKernel * 2.0f), 0.5f);

	Sigma = 0.1f;
	FluidStep = 0.005f;

	VerticalBounds = ParticlesPerDimension * 2.0f;
	HorizontalBounds = 2.0f;
	BoxVolume = HorizontalBounds * HorizontalBounds * HorizontalBounds;
	
	DampeningStrength = 0.1f;
	FrictionStrength = 0.4f;
	ParticleMass = (1.0f * RestDensity * BoxVolume) / ParticleCount;

	SmoothingKernelSize = glm::pow((3.0f * BoxVolume * AverageParticlesPerKernel) / (4.0f * glm::pi<float>() * ParticleCount), 1.0f/3.0f);
	ParticleSize = 5.0f * glm::pow((3.0f * ParticleMass) / (4.0f * glm::pi<float>() * RestDensity), 1.0f / 3.0f);

	Color = glm::vec3(0.0f, 1.0f, 0.0f);

	InitialiseParticles();
}


MucusSim::~MucusSim(void)
{
}
