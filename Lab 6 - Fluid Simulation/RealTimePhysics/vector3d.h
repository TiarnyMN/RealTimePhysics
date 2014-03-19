#ifndef __VECTOR3D_H__
#define __VECTOR3D_H__

#include <assert.h>

// forward declaration
namespace Fluid
{
	class VecInt3D;
	class Vector3D;
}

namespace Fluid
{

/** A general purpose 3D vector.
 * This vector class is used for all kinds of 3D mathematics. You can perform
 * basic operations like multiplication or dot products on it.
 *
 * \ingroup core
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2003/01/10 21:27:31 $
 *
 * \author malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	class Vector3D
	{
	public:
	/** Create a zero vector.
	 */
		Vector3D();

	/** Create a vector initialized by a scalar.
	 * All coordinates are set to the scalar value
	 * 
	 * \param f scalar initializer
	 */
		Vector3D(const float f);

	/** Create a vector using coordinates.
	 * 
	 * \param x x coordinate
	 * \param y y coordinate
	 * \param z z coordinate
	 */
		 Vector3D(const float x, const float y, const float z);

	/** Create a vector using coordinates.
	 * 
	 * \param x x coordinate
	 * \param y y coordinate
	 * \param z z coordinate
	 */
		 Vector3D(const int x, const int y, const int z);

	/** Create a copy of another vector.
	 * 
	 * \param v source vector
	 */
		 Vector3D(const Vector3D & v);

	/** Retrieve value.
	 * Use this operator for write access. See the () operator for
	 * read only retrieval.
	 * 
	 * \param i element index
	 * \return reference to the value
	 */
		float &operator[] (const int i);

	/** Retrieve value.
	 * Use this operator for read only access. See the [] operator for
	 * write access.
	 * 
	 * \param i  element index
	 * \return const reference to the value
	 */
		const float &operator() (const int i) const;

	/** Check wether vector equals to another vector.
	 * 
	 * \param v second vector
	 * \return true if equal, false otherwise
	 */
		bool operator==(const Vector3D & v) const;

	/** Check wether all three coordinates are lower than those of another vector
	 * 
	 * \param v second vector
	 * \return true if lower, false otherwise
	 */
		bool operator<(const Vector3D & v) const;

	/** Check wether all three coordinates are lower or equal than those of another vector
	 * 
	 * \param v second vector
	 * \return true if lower or equal, false otherwise
	 */
		bool operator<=(const Vector3D & v) const;

	/** Check wether all three coordinates are greater than those of another vector
	 * 
	 * \param v second vector
	 * \return true if greater, false otherwise
	 */
		bool operator>(const Vector3D & v) const;

	/** Check wether all three coordinates are greater or equal than those of another vector
	 * 
	 * \param v second vector
	 * \return true if greater or equal, false otherwise
	 */
		bool operator>=(const Vector3D & v) const;

	/** Copy from vector.
	 * Copies the values from another vector.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator=(const Vector3D & v);

	/** Copy from float
	 * Initializes all elements with a float value
	 * 
	 * \param f initial value
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator=(const float f);

	/** Add vector.
	 * Add vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator+=(const Vector3D & v);

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator-=(const Vector3D & v);

	/** Multiply with float.
	 * Multiplies all values by a scalar.
	 * 
	 * \param v source scalar
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator*=(const float f);

	/** Divide by float.
	 * Divides all values by a scalar.
	 * 
	 * \param v source scalar
	 * \return reference to this vector (use to chain operations)
	 */
		 Vector3D & operator/=(const float f);

	/** Multiply with float.
	 * Multiplies all values by a scalar.
	 * 
	 * \param v source scalar
	 * \return result vector
	 */
		Vector3D operator*(const float f) const;

	/** Dot product.
	 * Calculates the dot product with another vector.
	 * 
	 * \param v source vector
	 * \return dot product
	 */
		float operator*(const Vector3D & v) const;

	/** Add vector.
	 * Add vectors value by value.
	 * 
	 * \param v source vector
	 * \return result vector
	 */
		Vector3D operator+(const Vector3D & v) const;

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return result vector
	 */
		Vector3D operator-(const Vector3D & v) const;

	/** Euclidean norm squared.
	 * Calculates the squared euclidean norm. This operation is less expensive
	 * than the usual euclidean norm Norm2().
	 * 
	 * \return euclidean norm squared
	 */
		float Norm2_2() const;

	/** Euclidean norm.
	 * Calculates the euclidean norm. This is the square root of the dot
	 * product of this vector with itself. It is the length of the vector
	 * in n dimensional space. You can use this as error measurement if you
	 * calculate the length of the distance between two vectors.
	 * 
	 * \return euclidean norm
	 */
		float Norm2() const;

		VecInt3D GetVecInt3D() const;

	private:
		float mValue[3];						///< coordinates
	};

};															// namespace Fluid

#endif													// __VECTOR3D_H__