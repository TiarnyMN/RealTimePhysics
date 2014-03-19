#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <assert.h>

// forward declaration
namespace Fluid
{
	template < class T > class Matrix;
}

#include "config.h"

#include "vector.h"
#include "vecint3d.h"
#include "object.h"

namespace Fluid
{

/** A special sparse matrix.
 * This matrix is specialized to the needs in a fluid solver. It is a sparse
 * symmetric positive definite matrix that contains only entries on the main
 * diagonal and at most six other positions per row that represent the
 * adjacent grid cells in 3D. Because of this special structure, there's no
 * memory required to store the matrix itself, it can be calculated on the fly.
 * It's only purpose is to solve linear equations when calculating implicit
 * steps in the grid.
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
	template < class T > class Matrix
	{
	/** A row of the matrix.
	 * The matrix is sparse and contains 7 entries per row at most: One for 
	 * the center cell (the main diagonal) and one for each adjacent cell.
	 * This structure can be used to limit the rows to 6 columns.
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
		class Row
		{
		public:
			int mRowCount;
			int mColumn[6];
			int mCenterCount;
		};

	public:
	/** Creates and initializes the matrix to the given size.
	 * The matrix size is calculated automatically from the grid size.
	 * The grid size is also needed to calculate the proper positions
	 * of the adjacent fields.
	 * 
	 * \param size Specifies the size of the grid in three dimensions.
	 */
		 Matrix(const VecInt3D & gridsize):mGridSize(gridsize)
		{
			mSize = mGridSize(0) * mGridSize(1) * mGridSize(2);

			mRow = new Row[mSize];
		}

		~Matrix()
		{
			delete[]mRow;
		}

	private:

		float GetValueCenter(int x) const
		{
			return mCenter + mCenterAdj * mRow[x].mCenterCount;
		}

		float GetValue(int x, int y) const
		{
			if (x == y)
				GetValueCenter(x);

			for (int i = 1; i < mRow[x].mRowCount; i++)
			{
				if (mRow[x].mColumn[i] == y)
					return mAdjacent;
			}

			return 0.0f;
		}

	/** Multiply with a vector.
	 * Perform matrix-vector multiplication and return the result
	 *
	 * \param v source vector
	 * \return result vector
	 */
		Vector < T > operator*(const Vector < T > &v) const
		{
			Vector < T > nv(mSize);

			for (int x = 0; x < mSize; x++)
			{
				nv[x] = 0.0f;

				for (int y = 0; y < mRow[x].mRowCount; y++)
				{
					nv[x] += v(mRow[x].mColumn[y]) * mAdjacent;
				}
				nv[x] += v(x) * (mCenter + mCenterAdj * mRow[x].mCenterCount );
			}

			return nv;
		}

	public:
	/** Initialize matrix values.
	 * Sets the values for center and adjacent cells.
	 * 
	 * \param center base value for the center cell
	 * \param centeradj additional center value, gets multiplied by the number of adjacent cells
	 * \param adjacent value for the neighbours
	 */
		void InitValues(const float center, const float centeradj,
										const float adjacent)
		{
			mCenter = center;
			mCenterAdj = centeradj;
			mAdjacent = adjacent;
		}

	/** Initialize matrix structure.
	 * Initialize the matrix structure according to the object hierarchy.
	 * 
	 * \param object object hierarchy, used to exclude solid grid cells
	 */
		void InitNeighbour(const Object * object, const bool bounded)
		{
			int w = mGridSize(0);
			int wh = w * mGridSize(1);

			int i = 0;
			for (int z = 0; z < mGridSize(2); z++)
			{
				for (int y = 0; y < mGridSize(1); y++)
				{
					for (int x = 0; x < mGridSize(0); x++)
					{
						int col = 0;
						int cc_add = 0;

						if ((object == NULL) || !(object->InsideSolid(VecInt3D(x, y, z))))
						{
							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x, y, z - 1))))
							{
								if (z >= 1)
								{
									mRow[i].mColumn[col] = i - wh;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}
							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x, y, z + 1))))
							{
								if (z < mGridSize(2) - 1)
								{
									mRow[i].mColumn[col] = i + wh;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}

							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x, y - 1, z))))
							{
								if (y >= 1)
								{
									mRow[i].mColumn[col] = i - w;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}
							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x, y + 1, z))))
							{
								if (y < mGridSize(1) - 1)
								{
									mRow[i].mColumn[col] = i + w;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}

							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x - 1, y, z))))
							{
								if (x >= 1)
								{
									mRow[i].mColumn[col] = i - 1;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}
							if ((object == NULL)
									|| !(object->InsideSolid(VecInt3D(x + 1, y, z))))
							{
								if (x < mGridSize(0) - 1)
								{
									mRow[i].mColumn[col] = i + 1;
									col++;
								}
								else if (!bounded)
									cc_add++;
							}
						}

						mRow[i].mRowCount = col;
						mRow[i].mCenterCount = mRow[i].mRowCount + cc_add;

						if( mRow[i].mCenterCount == 0 ) // conjugate gradient hack: if 0, the matrix won't be spd
							mRow[i].mCenterCount = 1;

						i++;
					}
				}
			}
			assert(i == mSize);
		}

	/** Solve linear equations.
	 * Solves the linear equations given by the matrix and a right side vector.
	 * 
	 * \param target receives the result vector
	 * \param v source vector (right side)
	 * \return reference to the target vector (use to chain operations)
	 */
		Vector < T > &Solve(Vector < T > &target, const Vector < T > &v)
		{
			// initial x0
			for (int i = 0; i < mSize; i++)
				target[i] = 0.0f;

			if (v.Norm2_2() < gc_epsilon)	// nothing to do
				return target;

#ifdef CONJUGATE_GRADIENT
			// conjugate gradient

			Vector < T > r(mSize);
			Vector < T > z(mSize);
			Vector < T > p(mSize);
			Vector < T > q(mSize);

			Vector < float >m(mSize);
			for (int i = 0; i < mSize; i++)
				m[i] = 1.0f / GetValueCenter(i);

			float rho;								// rho^(i-1)
			float rhol = 0.0f;				// rho^(i-2)

			// r^0 = b - A * x^0
			r = v;										// x^0 = 0

			for (int i = 0; i < gc_max_iteration; i++)
			{
				// solve M * z^(i-1) = r^(i-1)
				z = r.Scale(m);

				// rho_(i-1) = r^(i-1)^T * z^(i-1)
				rho = r * z;

				if (i == 0)
				{
					// p^(i) = z^(0)
					p = z;
				} else
				{
					// beta_(i-1) = rho_(i-1) / rho_(i-2)
					float beta = rho / rhol;

					// p^(i) = z^(i-1) + beta_(i-1) * p^(i-1)
					//p = z + p * beta;
					p *= beta;
					p += z;
				}

				// q^(i) = A * p^(i)
				q = (*this) * p;

				// alpha_i = rho_(i-1) / p^(i)^T * q^(i)
				float alpha = rho / (p * q);

				// x^(i) = x^(i-1) + alpha_i * p^(i)
				target += p * alpha;

				// r^(i) = r^(i-1) - alpha_i * q^(i)
				r -= q * alpha;

				// convergence check
				if (r.Norm2_2() < gc_epsilon)
					break;

				rhol = rho;
			}
#else
			// gauss-seidel
			for (int s = 0; s < gc_max_iteration; s++)
			{
				for (int i = 0; i < mSize; i++)
				{
					target[i] = v(i);

					for (int j = 0; j < mSize; j++)
					{
						if (j != i)
						{
							target[i] -= target(j) * GetValue(i, j);
						}
					}

					target[i] /= GetValueCenter(i);
				}

				// check convergence
				Vector < T > r(mSize);
				r = (*this) * target;
				r -= v;
				if (r.Norm2_2() < gc_epsilon)
					break;
			}
#endif													// CONJUGATE GRADIENT
			return target;
		}

	private:
		int mSize;									///< height/width of the matrix (2D)
		VecInt3D mGridSize;					///< size of the grid the matrix is applied to

		float mCenter;							///< base value for the center cell
		float mCenterAdj;						///< additional center value, gets multiplied by the number of adjacent cells
		float mAdjacent;						///< value for the neighbours

		Row *mRow;
	};

};															// namespace Fluid

#endif													// __MATRIX_H____