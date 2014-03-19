#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <assert.h>
#include <cmath>

// forward declaration
namespace Fluid
{
	template < class T > class Vector;
}

namespace Fluid
{

/** A general purpose mathematical vector.
 * This vector provides all functionality needed to be used as container for
 * the grid cells and as mathematical vector for the linear equations solver
 * in the matrix class. It has many basic operations such as euclidean norm
 * and algebraic functions as multiplication, scaling etc.
 *
 * \ingroup core
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2003/01/10 21:27:30 $
 *
 * \author malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	template < class T > class Vector
	{
	public:
	/** Create a new vector.
	 * 
	 * \param size length (number of grid cells)
	 */
		Vector(int size)
		{
			assert(size > 0);

			mSize = size;
			mValue = new T[mSize];
		}

	/** Create a copy.
	 * This is called when creating a new vector as copy of an existing vector.
	 * 
	 * \param v source vector
	 */
		Vector(const Vector & v)
		{
			mSize = v.GetSize();
			mValue = new T[mSize];
			for (int i = 0; i < mSize; i++)
				mValue[i] = v(i);
		}

	/** Destroy vector.
	 * Deletes the allocated storage.
	 * 
	 */
		~Vector()
		{
			delete[]mValue;
		}

	/** Check wether vector equals to another vector.
	 * 
	 * \param v second vector
	 * \return true if equal, false otherwise
	 */
		bool operator==(const Vector & v)
		{
			bool equal = true;

			equal = equal && (mSize == v.GetSize());

			for (int i = 0; i < mSize; i++)
				equal = equal && (mValue[i] == v(i));

			return equal;
		}

	/** Copy from vector.
	 * Copies the values from another vector.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		Vector & operator=(const Vector & v)
		{
			assert(mSize == v.GetSize());

			for (int i = 0; i < mSize; i++)
				mValue[i] = v(i);

			return *this;
		}

	/** Copy from float
	 * Initializes all elements with a float value
	 * 
	 * \param f initial value
	 * \return reference to this vector (use to chain operations)
	 */
		Vector & operator=(const float f)
		{
			assert(mSize == v.GetSize());

			for (int i = 0; i < mSize; i++)
				mValue[i] = f;

			return *this;
		}

	/** Retrieve value.
	 * Use this operator for write access. See the () operator for
	 * read only retrieval.
	 * 
	 * \param i element index
	 * \return reference to the value
	 */
		T & operator[](int i)
		{
			assert((0 <= i) && (i < mSize));

			return mValue[i];
		}

	/** Retrieve value.
	 * Use this operator for read only access. See the [] operator for
	 * write access.
	 * 
	 * \param i  element index
	 * \return const reference to the value
	 */
		const T & operator() (int i) const
		{
			assert((0 <= i) && (i < mSize));

			return mValue[i];
		}

	/** Scale using another vector.
	 * Multiplies all values by the corresponding values of the source vector.
	 * Note that this method does not modify this vector, it returns the result
	 * in a new vector.
	 * 
	 * \param v source vector
	 * \return result of the scaling operation
	 */
		Vector < T > Scale(const Vector < float > & v) const
		{
			Vector < T > nv(*this);

			for (int i = 0; i < mSize; i++)
				nv[i] *= v(i);

			return nv;
		}

	/** Add vector.
	 * Add vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		Vector < T > &operator +=(const Vector < T > & v)
		{
			for (int i = 0; i < mSize; i++)
				(*this)[i] += v(i);

			return (*this);
		}

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		Vector < T > &operator -=(const Vector < T > & v)
		{
			for (int i = 0; i < mSize; i++)
				(*this)[i] -= v(i);

			return (*this);
		}

	/** Multiply with float.
	 * Multiplies all values by a scalar.
	 * 
	 * \param v source scalar
	 * \return reference to this vector (use to chain operations)
	 */
		Vector < T > &operator *=(const float f)
		{
			for (int i = 0; i < mSize; i++)
				(*this)[i] *= f;

			return (*this);
		}

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return result vector
	 */
		Vector < T > operator -(const Vector < T > & v) const
		{
			Vector < T > nv(*this);
			return nv -= v;
		}

	/** Multiply with float.
	 * Multiplies all values by a scalar.
	 * 
	 * \param v source scalar
	 * \return result vector
	 */
		Vector < T > operator *(const float f) const
		{
			Vector < T > nv(*this);
			return nv *= f;
		}

	/** Dot product.
	 * Calculates the dot product with another vector.
	 * 
	 * \param v source vector
	 * \return dot product
	 */
		float operator *(const Vector < T > & v) const
		{
			float sum = 0.0f;

			for (int i = 0; i < mSize; i++)
				 sum += (*this) (i) * v(i);

			 return sum;
		}

	/** Euclidean norm squared.
	 * Calculates the squared euclidean norm. This operation is less expensive
	 * than the usual euclidean norm Norm2().
	 * 
	 * \return euclidean norm squared
	 */
		float Norm2_2() const
		{
			float sum = 0.0f;

			for (int i = 0; i < mSize; i++)
				 sum += (*this) (i) * (*this) (i);

			 return sum;
		}

	/** Euclidean norm.
	 * Calculates the euclidean norm. This is the square root of the dot
	 * product of this vector with itself. It is the length of the vector
	 * in n dimensional space. You can use this as error measurement if you
	 * calculate the length of the distance between two vectors.
	 * 
	 * \return euclidean norm
	 */
		float Norm2() const
		{
			return sqrtf(Norm2_2());
		}

	/** Gets the size of the vector.
	 * 
	 * \return size
	 */
		int GetSize() const
		{
			return mSize;
		}

	private:
		int mSize;									///< number of elements
		T *mValue;									///< element array
	};

};															// namespace Fluid

#endif													// __VECTOR_H__