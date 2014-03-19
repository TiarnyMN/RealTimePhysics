#include "vecint3d.h"

#include <cmath>

using namespace Fluid;

VecInt3D::VecInt3D()
{
	mValue[0] = 0;
	mValue[1] = 0;
	mValue[2] = 0;
}

VecInt3D::VecInt3D(int x, int y, int z)
{
	mValue[0] = x;
	mValue[1] = y;
	mValue[2] = z;
}

VecInt3D::VecInt3D(const VecInt3D & v)
{
	for (int i = 0; i < 3; i++)
		mValue[i] = v.mValue[i];
}

int &VecInt3D::operator[] (int i)
{
	assert((0 <= i) && (i < 3));

	return mValue[i];
}

const int &VecInt3D::operator() (int i)
const
{
	assert((0 <= i) && (i < 3));

	return mValue[i];
}

bool VecInt3D::operator==(const VecInt3D & v) const
{
	return (mValue[0] == v(0)) && (mValue[1] == v(1)) && (mValue[2] == v(2));
}

bool VecInt3D::operator<(const VecInt3D & v) const
{
	return (mValue[0] < v(0)) && (mValue[1] < v(1)) && (mValue[2] < v(2));
}

bool VecInt3D::operator<=(const VecInt3D & v) const
{
	return (mValue[0] <= v(0)) && (mValue[1] <= v(1)) && (mValue[2] <= v(2));
}

bool VecInt3D::operator>(const VecInt3D & v) const
{
	return (mValue[0] > v(0)) && (mValue[1] > v(1)) && (mValue[2] > v(2));
}

bool VecInt3D::operator>=(const VecInt3D & v) const
{
	return (mValue[0] >= v(0)) && (mValue[1] >= v(1)) && (mValue[2] >= v(2));
}

VecInt3D & VecInt3D::operator =(const VecInt3D & v)
{
	mValue[0] = v(0);
	mValue[1] = v(1);
	mValue[2] = v(2);

	return *this;
}

VecInt3D & VecInt3D::operator +=(const VecInt3D & v)
{
	mValue[0] += v(0);
	mValue[1] += v(1);
	mValue[2] += v(2);

	return *this;
}

VecInt3D & VecInt3D::operator -=(const VecInt3D & v)
{
	mValue[0] -= v(0);
	mValue[1] -= v(1);
	mValue[2] -= v(2);

	return *this;
}

VecInt3D VecInt3D::operator +(const VecInt3D & v) const
{
	VecInt3D nv = *this;
	nv += v;
	return nv;
}

VecInt3D VecInt3D::operator -(const VecInt3D & v) const
{
	VecInt3D nv = *this;
	nv -= v;
	return nv;
}

int VecInt3D::operator *(const VecInt3D & v) const
{
	return (*this) (0) * v(0) + (*this) (1) * v(1) + (*this) (2) * v(2);
}

int VecInt3D::Norm2_2() const
{
	return mValue[0] * mValue[0] + mValue[1] * mValue[1] + mValue[2] * mValue[2];
}

float VecInt3D::Norm2() const
{
	return sqrtf(static_cast < float >(Norm2_2()));
}