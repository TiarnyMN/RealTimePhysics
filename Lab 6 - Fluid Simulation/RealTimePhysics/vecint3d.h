#ifndef __VECINT3D_H__
#define __VECINT3D_H__

#include <assert.h>

// forward declaration
namespace Fluid
{
	class VecInt3D;
}

namespace Fluid
{

/** Vertex position.
 * This class is used to store vertex positions. Since the grid consists
 * of discreet cells, the vertices can be stored using integer coordinates.
 *
 * \ingroup core
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2003/01/10 21:27:30 $
 *
 * \author Malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	class VecInt3D
	{
	public:
		VecInt3D();

	/** Creates a vertex position using seperated coordinates.
	 * 
	 * \param x x coordinate 
	 * \param y y coordinate
	 * \param z z coordinate
	 */
		VecInt3D(int x, int y, int z);

	/** Creates a copy of another vertex position.
	 * 
	 * \param v source vector
	 */
		 VecInt3D(const VecInt3D & v);

	/** Retrieve coordinate.
	 * 
	 * \param i coordinate index
	 * \return coordinate
	 */
		const int &operator() (int i) const;

	/** Retrieve coordinate for write access.
	 * 
	 * \param i coordinate index
	 * \return reference to coordinate
	 */
		int &operator[] (int i);

	/** Compare with another vector.
	 * 
	 * \param v second vector
	 * \return true if equal, false otherwise
	 */
		bool operator==(const VecInt3D & v) const;

	/** Check wether all three coordinates are lower than those of another vector
	 * 
	 * \param v second vector
	 * \return true if lower, false otherwise
	 */
		bool operator<(const VecInt3D & v) const;

	/** Check wether all three coordinates are lower or equal than those of another vector
	 * 
	 * \param v second vector
	 * \return true if lower or equal, false otherwise
	 */
		bool operator<=(const VecInt3D & v) const;

	/** Check wether all three coordinates are greater than those of another vector
	 * 
	 * \param v second vector
	 * \return true if greater, false otherwise
	 */
		bool operator>(const VecInt3D & v) const;

	/** Check wether all three coordinates are greater or equal than those of another vector
	 * 
	 * \param v second vector
	 * \return true if greater or equal, false otherwise
	 */
		bool operator>=(const VecInt3D & v) const;

	/** Copy from vector.
	 * Copies the values from another vector.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 VecInt3D & operator=(const VecInt3D & v);

	/** Add vector.
	 * Add vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 VecInt3D & operator+=(const VecInt3D & v);

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return reference to this vector (use to chain operations)
	 */
		 VecInt3D & operator-=(const VecInt3D & v);

	/** Dot product.
	 * Calculates the dot product with another vector.
	 * 
	 * \param v source vector
	 * \return dot product
	 */
		int operator*(const VecInt3D & v) const;

	/** Add vector.
	 * Add vectors value by value.
	 * 
	 * \param v source vector
	 * \return result vector
	 */
		VecInt3D operator+(const VecInt3D & v) const;

	/** Subtract vector.
	 * Subtract vectors value by value.
	 * 
	 * \param v source vector
	 * \return result vector
	 */
		VecInt3D operator-(const VecInt3D & v) const;

	/** Euclidean norm squared.
	 * Calculates the squared euclidean norm. This operation is less expensive
	 * than the usual euclidean norm Norm2().
	 * 
	 * \return euclidean norm squared
	 */
		int Norm2_2() const;

	/** Euclidean norm.
	 * Calculates the euclidean norm. This is the square root of the dot
	 * product of this vector with itself. It is the length of the vector
	 * in n dimensional space. You can use this as error measurement if you
	 * calculate the length of the distance between two vectors.
	 * 
	 * \return euclidean norm
	 */
		float Norm2() const;

	private:
		int mValue[3];							///< coordinates
	};

};															// namespace Fluid

#endif													// __VECINT3D_H__