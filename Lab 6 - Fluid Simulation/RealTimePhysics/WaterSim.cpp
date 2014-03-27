#include "WaterSim.h"

//float ParticlePerDimension;
//	float ParticleCount;
//	float SurfaceThreshold;
//
//	float StiffnessConstant;
//	float ViscosityConstant;
//	float AverageParticlesPerKernel;
//	float FluidStep;
//
//	float VerticalBounds;
//	float HorizontalBounds;
//
//	float BoxVolume;
//
//	float DampeningStrength;
//
//	float ParticleMass;
//
//	float SmoothingKernelSize;
//	float ParticleSize;
//
//	float Sigma;
//
//	bool wallEnabled;

//#define REST_DENSITY 1000.0f
//
//#define PARTICLES_PER_DIMENSION 10.0f
//#define PARTICLE_COUNT (PARTICLES_PER_DIMENSION * PARTICLES_PER_DIMENSION * PARTICLES_PER_DIMENSION)
//
//#define STIFFNESS_CONSTANT 100.0f
//
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

WaterSim::WaterSim(void)
{
	RestDensity = 1000.0f;
	ParticlesPerDimension = 25.0f;
	ParticleCount = ParticlesPerDimension * ParticlesPerDimension * ParticlesPerDimension;

	StiffnessConstant = 100.0f;
	ViscosityConstant = 50.0f;

	AverageParticlesPerKernel = 20.0f;
	SurfaceThreshold = glm::pow(RestDensity / (AverageParticlesPerKernel * 2.0f), 0.5f);

	Sigma = 0.01f;
	FluidStep = 0.005f;

	VerticalBounds = ParticlesPerDimension * 2.0f;
	HorizontalBounds = 2.0f;
	BoxVolume = HorizontalBounds * HorizontalBounds * HorizontalBounds;
	
	DampeningStrength = 0.2f;
	ParticleMass = (1.0f * RestDensity * BoxVolume) / ParticleCount;

	SmoothingKernelSize = glm::pow((3.0f * BoxVolume * AverageParticlesPerKernel) / (4.0f * glm::pi<float>() * ParticleCount), 1.0f/3.0f);
	ParticleSize = 5.0f * glm::pow((3.0f * ParticleMass) / (4.0f * glm::pi<float>() * RestDensity), 1.0f / 3.0f);

	InitialiseParticles();
}


WaterSim::~WaterSim(void)
{
}
