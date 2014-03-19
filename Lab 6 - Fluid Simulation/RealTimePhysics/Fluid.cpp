#include "fluid.h"

using namespace Fluid;

Solver::Solver(const Vector3D size, VecInt3D cells)
{
	mViscosity = 0.5f;
	mDiffusion = 0.5f;
	mDissipation = 0.0f;

	mSize = size;
	mCells = cells;

	if ((mCells(0) == 0) && (mCells(1) == 0))
	{
		mCells[0] =	static_cast < int >(static_cast < float >(mCells(2)) * size(0) / size(2));
		mCells[1] =	static_cast < int >(static_cast < float >(mCells(2)) * size(1) / size(2));
	}

	if ((mCells(0) == 0) && (mCells(2) == 0))
	{
		mCells[0] =	static_cast < int >(static_cast < float >(mCells(1)) * size(0) / size(1));
		mCells[2] =	static_cast < int >(static_cast < float >(mCells(1)) * size(2) / size(1));
	}

	if ((mCells(1) == 0) && (mCells(2) == 0))
	{
		mCells[1] =	static_cast < int >(static_cast < float >(mCells(0)) * size(1) / size(0));
		mCells[2] =	static_cast < int >(static_cast < float >(mCells(0)) * size(2) / size(0));
	}

	assert((mCells(0) > 0) && (mCells(1) > 0) && (mCells(2) > 0));

	mDensity[0] = new Grid < float >(mCells);
	mDensity[1] = new Grid < float >(mCells);
	mVelocity[0] = new Grid < Vector3D > (mCells);
	mVelocity[1] = new Grid < Vector3D > (mCells);

	mSource = 0;
	mDest = 1;
}

Solver::~Solver()
{
	delete mDensity[0];
	delete mDensity[1];
	delete mVelocity[0];
	delete mVelocity[1];
}

void Solver::SetViscosity(float viscosity)
{
	float h = GridSpace(1.0f);
	mViscosity = viscosity / (h * h);
}

float Solver::GetViscosity()
{
	return mViscosity;
}

void Solver::SetDiffusion(float diffusion)
{
	float h = GridSpace(1.0f);
	mDiffusion = diffusion / (h * h);
}

float Solver::GetDiffusion()
{
	return mDiffusion;
}


void Solver::SetDissipation(float dissipation)
{
	mDissipation = dissipation;
}

float Solver::GetDissipation()
{
	return mDissipation;
}

const Grid < float >*Solver::GetDensity() const
{
	return mDensity[0];
}

const Grid < Vector3D > *Solver::GetVelocity() const
{
	return mVelocity[mSource];
}

void Solver::SetBoundedGrid(const bool bounded)
{
	mDensity[0]->SetBounded(bounded);
	mDensity[1]->SetBounded(bounded);
	mVelocity[0]->SetBounded(bounded);
	mVelocity[1]->SetBounded(bounded);
}

void Solver::SetObject(const Object * object, const bool update_solid)
{
	mDensity[0]->SetObject(object, update_solid);
	mDensity[1]->SetObject(object, update_solid);
	mVelocity[0]->SetObject(object, update_solid);
	mVelocity[1]->SetObject(object, update_solid);
}

Vector3D Solver::GridSpace(const Vector3D v) const
{
	Vector3D nv;

	nv[0] = static_cast < float >(mCells(0)) * v(0) / mSize(0);
	nv[1] = static_cast < float >(mCells(1)) * v(1) / mSize(1);
	nv[2] = static_cast < float >(mCells(2)) * v(2) / mSize(2);

	return nv;
}

float Solver::GridSpace(const float f) const
{
	return static_cast < float >(mCells(0)) * f / mSize(0);
}

VecInt3D Solver::GridSpaceInt(const Vector3D v) const
{
	VecInt3D nv;

	nv[0] = static_cast < int >(static_cast < float >(mCells(0)) * v(0) / mSize(0));
	nv[1] =	static_cast < int >(static_cast < float >(mCells(1)) * v(1) / mSize(1));
	nv[2] =	static_cast < int >(static_cast < float >(mCells(2)) * v(2) / mSize(2));

	return nv;
}

int Solver::GridSpaceInt(const float f) const
{
	return static_cast < int >(static_cast < float >(mCells(0)) * f / mSize(0));
}

void Solver::DoStep(float time)
{
	// density
	mDensity[0]->AddSource(time);
	mDensity[0]->Diffuse(mDensity[1], mDiffusion, time);
	mDensity[1]->Move(mDensity[0], mVelocity[mSource], time);
	*mDensity[0] /= 1 + mDissipation * time;

	// velocity
	mVelocity[mSource]->AddForce(time);
	mVelocity[mSource]->Diffuse(mVelocity[mDest], mViscosity, time);
	mVelocity[mDest]->Move(mVelocity[mSource], mVelocity[mDest], time);
	mVelocity[mSource]->ConserveMass(mVelocity[mDest]);

	if (mSource == 0)
	{
		mSource = 1;
		mDest = 0;
	} else
	{
		mSource = 0;
		mDest = 1;
	}
}