#include "object.h"

using namespace Fluid;

int Object::mMaxID = 0;	// changed by tobias (23.01.03)

Object::Object() :
	mForce(0.0f),
	mSource(0.0f),
	mSolid(false),
	mChanged(false),
	mpSolver(0)
{
	mMaxID++;							// changed by tobias (23.01.03)
	mID = mMaxID;		// changed by tobias (23.01.03)
}

void Object::SetForce(const Vector3D force)
{
	mForce = force;
	if( mpSolver != 0 )
	{
		//float h = mpSolver->GridSpace(1.0f);
		mForceGrid = mForce * 1.0f;// * (1 / h);
	}
}

void Object::SetSource(const float source)
{
	mSource = source;
	mSourceGrid = mSource;
}

void Object::Init(const Solver & solver)
{
	mpSolver = &solver;
	//float h = mpSolver->GridSpace(1.0f);
	mForceGrid = mForce * 1.0f;// * (1 / h);
	mSourceGrid = mSource;
}

void Object::SetSolid(const bool solid)
{
	mSolid = solid;
}

bool Object::InsideSolid(const VecInt3D pos) const
{
	return (mSolid && Inside(pos));
}

const Vector3D Object::GetForce(const VecInt3D pos) const
{
	if (Inside(pos))
		return mForceGrid;
	else
		return 0;
}

float Object::GetSource(const VecInt3D pos) const
{
	if (Inside(pos))
		return mSourceGrid;
	else
		return 0;
}

void Object::SetChanged(bool changed)
{
	mChanged = changed;
}

bool Object::GetChanged() const
{
	bool changed = mChanged;
	return changed;
}

Sphere::Sphere(const Vector3D center, const float radius):
mCenter(center)
{
	mRadius = radius;
}

bool Sphere::Inside(const VecInt3D pos) const
{
	return ((pos - mCenterGrid).Norm2_2() <= mRadius2Grid);
}

void Sphere::Init(const Solver & solver)
{
	Object::Init(solver);

	mCenterGrid = solver.GridSpaceInt(mCenter);
	mRadius2Grid =
			static_cast <
			int >(solver.GridSpace(mRadius) * solver.GridSpace(mRadius));
}

void Sphere::Translate(float x, float y, float z) // changed by tobias (23.01.03)
{
	mCenter[0] += x; mCenter[1] +=y; mCenter[2] += z;
}

Box::Box(const Vector3D lt, const Vector3D gt):mCornerLT(lt), mCornerGT(gt)
{
	assert(mCornerLT <= mCornerGT);
}

bool Box::Inside(const VecInt3D pos) const
{
	return ((mCornerLTGrid <= pos) && (pos <= mCornerGTGrid));
}

void Box::Init(const Solver & solver)
{
	Object::Init(solver);

	mCornerLTGrid = solver.GridSpaceInt(mCornerLT);
	mCornerGTGrid = solver.GridSpaceInt(mCornerGT);
}

void Box::Translate(float x, float y, float z) // changed by tobias (23.01.03)
{
	mCornerLT[0] += x; mCornerLT[1] +=y; mCornerLT[2] += z;
	mCornerGT[0] += x; mCornerGT[1] +=y; mCornerGT[2] += z;
}

void Group::Init(const Solver & solver)
{
	Object::Init(solver);

	for (ObjList::iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		(*i)->Init(solver);
	}
}

void Group::Add(Object * object)
{
	mElements.push_back(object);
	SetChanged(true);
}

bool Group::Inside(const VecInt3D pos) const
{
	for (ObjList::const_iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		if ((*i)->Inside(pos))
			return true;
	}

	return false;
}

bool Group::InsideSolid(const VecInt3D pos) const
{
	for (ObjList::const_iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		if ((*i)->InsideSolid(pos))
			return true;
	}

	return false;
}

const Vector3D Group::GetForce(const VecInt3D pos) const
{
	Vector3D force = Object::GetForce(pos);

	for (ObjList::const_iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		force += (*i)->GetForce(pos);
	}

	return force;
}

float Group::GetSource(const VecInt3D pos) const
{
	float source = Object::GetSource(pos);

	for (ObjList::const_iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		source += (*i)->GetSource(pos);
	}

	return source;
}

bool Group::GetChanged() const
{
	bool changed = Object::GetChanged();

	for (ObjList::const_iterator i = mElements.begin(); i != mElements.end(); i++)
	{
		changed = changed || (*i)->GetChanged();
	}

	return changed;
}