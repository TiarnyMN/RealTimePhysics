#include "vector3d.h"
#include "vecint3d.h"

#include <cmath>

using namespace Fluid;

Vector3D::Vector3D()
{
	mValue[0] = 0.0f;
	mValue[1] = 0.0f;
	mValue[2] = 0.0f;
}

Vector3D::Vector3D(const float f)
{
	mValue[0] = f;
	mValue[1] = f;
	mValue[2] = f;
}

Vector3D::Vector3D(const float x, const float y, const float z)
{
	mValue[0] = x;
	mValue[1] = y;
	mValue[2] = z;
}

Vector3D::Vector3D(const int x, const int y, const int z)
{
	mValue[0] = static_cast < float >(x);
	mValue[1] = static_cast < float >(y);
	mValue[2] = static_cast < float >(z);
}

Vector3D::Vector3D(const Vector3D & v)
{
	for (int i = 0; i < 3; i++)
		mValue[i] = v(i);
}

float &Vector3D::operator[] (const int i)
{
	assert((0 <= i) && (i < 3));

	return mValue[i];
}

const float &Vector3D::operator() (const int i)
const
{
	assert((0 <= i) && (i < 3));

	return mValue[i];
}

bool Vector3D::operator==(const Vector3D & v) const
{
	return (mValue[0] == v(0)) && (mValue[1] == v(1)) && (mValue[2] == v(2));
}

bool Vector3D::operator<(const Vector3D & v) const
{
	return (mValue[0] < v(0)) && (mValue[1] < v(1)) && (mValue[2] < v(2));
}

bool Vector3D::operator<=(const Vector3D & v) const
{
	return (mValue[0] <= v(0)) && (mValue[1] <= v(1)) && (mValue[2] <= v(2));
}

bool Vector3D::operator>(const Vector3D & v) const
{
	return (mValue[0] > v(0)) && (mValue[1] > v(1)) && (mValue[2] > v(2));
}

bool Vector3D::operator>=(const Vector3D & v) const
{
	return (mValue[0] >= v(0)) && (mValue[1] >= v(1)) && (mValue[2] >= v(2));
}

Vector3D & Vector3D::operator =(const Vector3D & v)
{
	mValue[0] = v(0);
	mValue[1] = v(1);
	mValue[2] = v(2);

	return *this;
}

Vector3D & Vector3D::operator =(const float f)
{
	mValue[0] = f;
	mValue[1] = f;
	mValue[2] = f;

	return *this;
}

Vector3D & Vector3D::operator +=(const Vector3D & v)
{
	mValue[0] += v(0);
	mValue[1] += v(1);
	mValue[2] += v(2);

	return *this;
}

Vector3D & Vector3D::operator -=(const Vector3D & v)
{
	mValue[0] -= v(0);
	mValue[1] -= v(1);
	mValue[2] -= v(2);

	return *this;
}

Vector3D & Vector3D::operator *=(const float f)
{
	mValue[0] *= f;
	mValue[1] *= f;
	mValue[2] *= f;

	return *this;
}

Vector3D & Vector3D::operator /=(const float f)
{
	mValue[0] /= f;
	mValue[1] /= f;
	mValue[2] /= f;

	return *this;
}

Vector3D Vector3D::operator *(const float f) const
{
	Vector3D nv = *this;
	nv *= f;
	return nv;
}

Vector3D Vector3D::operator +(const Vector3D & v) const
{
	Vector3D nv = *this;
	nv += v;
	return nv;
}

Vector3D Vector3D::operator -(const Vector3D & v) const
{
	Vector3D nv = *this;
	nv -= v;
	return nv;
}

float Vector3D::operator *(const Vector3D & v) const
{
	return (*this) (0) * v(0) + (*this) (1) * v(1) + (*this) (2) * v(2);
}

float Vector3D::Norm2_2() const
{
	return mValue[0] * mValue[0] + mValue[1] * mValue[1] + mValue[2] * mValue[2];
}

float Vector3D::Norm2() const
{
	return sqrtf(Norm2_2());
}

VecInt3D Vector3D::GetVecInt3D() const
{
	return VecInt3D( static_cast<int>(mValue[0]+0.5f),static_cast<int>(mValue[1]+0.5f),static_cast<int>(mValue[2]+0.5f));
}
