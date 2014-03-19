#ifndef __GRID_H__
#define __GRID_H__

#include <assert.h>
#include <iostream>

// forward declaration
namespace Fluid
{
	template < class T > class Grid;
}

#include "config.h"
#include "vector.h"
#include "vector3d.h"
#include "vecint3d.h"
#include "matrix.h"
#include "object.h"

namespace Fluid
{

/** The simulation grid.
 * 
 * The grid contains the vertices the simulation consists of. You can apply
 * the basic simulation methods to it, for example Diffuse(), and it will
 * update the vertices accordingly. It is implemented as a template due to
 * it's generic nature: You can use the same basic operations for both
 * scalar and vector fields. The template argument is the type of the
 * vertices. You might want to use float for scalar grids (particle density 
 * grids) and Vector3D for vector grids (velocity grids).
 *
 * \ingroup core
 *
 * \version $Revision: 1.3 $
 *
 * \date $Date: 2003/01/13 17:36:32 $
 *
 * \author malte
 *
 * \todo
 * Add vortice reinforcement
 * Exclude solid cells from calculations.
 *
 * \bug 
 *
 */
	template < class T > class Grid
	{
	public:
	/** Creates and initializes the grid to the given size.
	 * The grid size cannot be changed later.
	 * 
	 * \param size Specifies the size of the grid in three dimensions.
	 */
	Grid(VecInt3D size):
		mSize(size),
				mSize1D(mSize(0) * mSize(1) * mSize(2)),
				mValue(mSize1D), mMatrixDiffusion(mSize), mMatrixConserveMass(mSize)
		{
			Reset();
			mBounded = true;
			mpObject = NULL;
		}

		// ==================== access ====================
	public:
	/** Determines wether a vertex belongs to the fluid simulation.
	 * At this moment all vertices inside the grid are valid. This may change
	 * to allow obstacles inside the grid. Vertices inside obstacles will
	 * be invalid.
	 * 
	 * \param x x coordinate
	 * \param y y coordinate
	 * \param z z coordinate
	 * \return true if vertex is valid, false otherwise
	 */
		 bool ValidPos(const int x, const int y, const int z) const
		{
			return ((x >= 0) && (x < mSize(0))) && ((y >= 0) && (y < mSize(1)))
					&& ((z >= 0) && (z < mSize(2)));
		}

	/** Retrieves a vertex value.
	 * 
	 * \param x x coordinate
	 * \param y y coordinate
	 * \param z z coordinate
	 * \return the value of the vertex or zero if the coordinates specify an invalid vertex
	 */
		T operator() (const int x, const int y, const int z) const
		{
			if (ValidPos(x, y, z))
				return mValue((z * mSize(1) + y) * mSize(0) + x);
			else
			{
				T t = 0.0f;
				 return t;
			}
		}

		T & Ref(const int x, const int y, const int z)
		{
			assert(ValidPos(x, y, z));
			return mValue[(z * mSize(1) + y) * mSize(0) + x];
		}

	/** Retrieves a vertex value
	 * This function provides the same functionality as the () operator
	 * 
	 * \param pos vertex position
	 * \return the value of the vertex or zero if the coordinates specify an invalid vertex
	 */
		T GetValue(const VecInt3D pos) const
		{
			return (*this) (pos(0), pos(1), pos(2));
		}

	/** Retrieve vertex value in linear fashion.
	 * Use this if you want to access all vertices in a for-loop or similar.
	 * DO NOT assume a specified grid layout, use the explicit 3D functions
	 * if you have to know the grid position of the vertex.
	 * 
	 * \param i linear position
	 * \return vertex value
	 */
		const T & operator() (int i) const
		{
			return mValue(i);
		}

	/** Retrieve vertex value in linear fashion.
	 * Use this if you want to access all vertices in a for-loop or similar.
	 * DO NOT assume a specified grid layout, use the explicit 3D functions
	 * if you have to know the grid position of the vertex.
	 * 
	 * \param i linear position
	 * \return vertex value
	 */
		const T & GetValue(const int i) const
		{
			return mValue(i);
		}

	/** Set a vertex to a given value.
	 * Use this function to set the value of a specific vertex.
	 * 
	 * \param pos 
	 * \param value 
	 */
		void SetValue(const VecInt3D pos, const T value)
		{
			mValue[(pos(2) * mSize(1) + pos(1)) * mSize(0) + pos(0)] = value;
		}

	/** Retrieves the vertices as one vector.
	 * The main purpose of this method is to allow easy matrix operations
	 * using the grid vertices.
	 * 
	 * \return A vector with all vertices
	 */
		Vector < T > &GetVector()
		{
			return mValue;
		}

	/** Resets the grid to zero.
	 * Use this to reset the simulation. All vertices are set back to zero.
	 * 
	 */
		void Reset()
		{
			for (int i = 0; i < mSize1D; i++)
				mValue[i] *= 0;
		}
	private:
	/** Returns the size of the grid.
	 * 
	 * \return The size that was passed to the constructor
	 */
		const VecInt3D & GetSize() const
		{
			return mSize;
		}

	/** Perform a linear interpolation between to values.
	 * 
	 * \param t0 first value
	 * \param t1 second value
	 * \param pos position between t0 (= 0.0f) and t1 (= 1.0f)
	 * \return interpolated calue
	 */
		T Interpolate(const T & t0, const T & t1, const float pos) const
		{
			return t0 * (1 - pos) + t1 * pos;
		}

	public:
	/** Retrieves an interpolated value.
	 * Use this function to determine values on positions that do not lie on 
	 * the grid. The surrounding vertices are used to interpolate the wanted
	 * value.
	 * 
	 * \param v sample position
	 * \return interpolated value
	 */
		 T GetValue(const Vector3D & v) const
		{
			// interpolate
			T i[4];

			// x
			 i[0] = Interpolate((*this) ((int) (v(0)), (int) (v(1)), (int) (v(2))),
													(*this) ((int) (v(0) + 1), (int) (v(1)),
																	 (int) (v(2))), v(0) - floorf(v(0)));

			 i[1] =
					Interpolate((*this) ((int) (v(0)), (int) (v(1) + 1), (int) (v(2))),
											(*this) ((int) (v(0) + 1), (int) (v(1) + 1),
															 (int) (v(2))), v(0) - floorf(v(0)));

			 i[2] =
					Interpolate((*this) ((int) (v(0)), (int) (v(1)), (int) (v(2) + 1)),
											(*this) ((int) (v(0) + 1), (int) (v(1)),
															 (int) (v(2) + 1)), v(0) - floorf(v(0)));

			 i[3] =
					Interpolate((*this)
											((int) (v(0)), (int) (v(1) + 1), (int) (v(2) + 1)),
											(*this) ((int) (v(0) + 1), (int) (v(1) + 1),
															 (int) (v(2) + 1)), v(0) - floorf(v(0)));

			// y
			 i[0] = Interpolate(i[0], i[1], v(1) - floorf(v(1)));
			 i[1] = Interpolate(i[2], i[3], v(1) - floorf(v(1)));

			// z
			 i[0] = Interpolate(i[0], i[1], v(2) - floorf(v(2)));

			 return i[0];
		}

	/** Retrieves a vertex value.
	 * 
	 * \param x x coordinate
	 * \param y y coordinate
	 * \param z z coordinate
	 * \return the value of the vertex or zero if the coordinates specify an invalid vertex
	 */
		T operator() (const float x, const float y, const float z) const
		{
			return GetValue( Vector3D(x,y,z) );
		}

	/** Determines wether coordinates are outside grid or inside a solid object.
		* \param pos position
		* \return true if valid
		*/
	bool IsValidCell( int x, int y, int z )
	{
		return( ( !mBounded || ( ( x >= 0 ) && ( x < mSize(0) ) && ( y >= 0 ) && ( y < mSize(1) ) && ( z >= 0 ) && ( z < mSize(2) ) ) ) && !mpObject->InsideSolid( VecInt3D(x,y,z) ) );
	}
	bool IsValidCell( VecInt3D pos )
	{
		return IsValidCell( pos[0],pos[1],pos[2] );
	}


	private:
	/** Retrieve vertex value in linear fashion.
	 * Use this if you want to access all vertices in a for-loop or similar
	 * for write access. If you only want to read, you should go for
	 * GetValue().
	 * DO NOT assume a specified grid layout, use the explicit 3D functions
	 * if you have to know the grid position of the vertex.
	 *
	 * \param i linear position
	 * \return vertex value
	 */
		 T & operator[] (const int i)
		{
			return mValue[i];
		}

		// ==================== calculate ====================
	public:
	/** Calculate Diffusion.
	 * Calculates the diffusion of the grid at a given rate for the specified
	 * time step and stores the result in a target grid.
	 *
	 * \param target receives the result
	 * \param rate diffusion rate, a non-negative value
	 * \param time delta time
	 * \return a pointer to the target grid (use to chain operations)
	 */
		Grid < T > *Diffuse(Grid < T > *target, const float rate, const float time)
		{
			mMatrixDiffusion.InitValues(1.0f, rate * time, -rate * time);

			// implicit step
			mMatrixDiffusion.Solve(target->GetVector(), mValue);

			return target;
		}
	
		/** Trace Particle.
		 * Used to approximate the position a particle at position pos would 
		 * have after the specified time.
		 * 
		 * \param pos position, is also used to return the new position
		 * \param time time step, can be negative
		 * \return same as pos
		 */
		Vector3D* TraceParticle( Vector3D* pos, const float time ) const
		{
#if defined( EULER_FIRST_ORDER )
			// euler first order
			Vector3D np = *pos + GetValue(*pos) * time;
#elif defined( RUNGEKUTTA_SECOND_ORDER )
			// runge-kutta second order
			Vector3D k1 = GetValue(*pos) * time;
			Vector3D k2 = GetValue(*pos - k1) * time;
			Vector3D np = *pos + (k1 + k2) * 0.5f;
#else														// RUNGEKUTTA_FOURTH_ORDER
			// runge-kutta fourth order
			Vector3D k1 = GetValue(*pos) * time;
			Vector3D k2 = GetValue(*pos - k1 * 0.5f) * time;
			Vector3D k3 = GetValue(*pos - k2 * 0.5f) * time;
			Vector3D k4 = GetValue(*pos - k3) * time;
			Vector3D np = *pos + (k1 + (k2 + k3) * 2.0f + k4) * (1.0f / 6.0f);
#endif
			if( !mpObject->InsideSolid( np.GetVecInt3D() ) )
				*pos = np;

			return pos;
		}

	/** Calculate Movement.
	 * Calculates the movement of the grid according to the passed velocity grid and stores
	 * the result in a target grid.
	 * 
	 * \param target receives the result
	 * \param velocity velocity grid
	 * \param time delta time
	 * \return a pointer to the target grid (use to chain operations)
	 *
	 * \todo improve particle tracer (use runge-kutta instead of explicit euler)
	 *
	 */
		Grid < T > *Move(Grid < T > *target, const Grid < Vector3D > *velocity,
										 const float time) const
		{
			for (int z = 0; z < mSize(2); z++)
			{
				for (int y = 0; y < mSize(1); y++)
				{
					for (int x = 0; x < mSize(0); x++)
					{
						Vector3D pos(static_cast < float >(x), static_cast < float >(y),
												 static_cast < float >(z));

						if (!mpObject->InsideSolid(VecInt3D(x, y, z)))
						{
							// -- trace particle

							// step backwards in time
							velocity->TraceParticle( &pos, -time );

							 if (mBounded)
							{
								if (pos(0) < 0.0f)
									pos[0] = 0.0f;
								if (pos(0) >= mSize(0))
									pos[0] = static_cast < float >(mSize(0) - 1);

								if (pos(1) < 0.0f)
									 pos[1] = 0.0f;
								if (pos(1) >= mSize(1))
									 pos[1] = static_cast < float >(mSize(1) - 1);

								if (pos(2) < 0.0f)
									 pos[2] = 0.0f;
								if (pos(2) >= mSize(2))
									 pos[2] = static_cast < float >(mSize(2) - 1);
							}
							// -- interpolate
							target->SetValue(VecInt3D(x, y, z), GetValue(pos));
						} else
						{
							target->SetValue(VecInt3D(x, y, z), 0.0f);
						}
					}
				}
			}

			return target;
		}

	/** Calculate mass conservation.
	 * Calculates the mass conservation of the grid. This means that the sum
	 * of incoming and outgoing flows in each grid cell has to be zero.
	 * This is a necessary step in velocity fields since you cannot store
	 * velocity in grid cells. Mass conversation provides nice vortices 
	 * instead.
	 *
	 * \todo bear in mind the mass conservation at object boundaries. At the
	 * moment, it appears as if some of the velocity is lost at boundaries.
	 * 
	 * \param target receives the result
	 * \return a pointer to the target grid (use to chain operations)
	 */
		Grid < T > *ConserveMass(Grid < T > *target)
		{
			Vector < float >b(mSize1D);
			int i = 0;
			for (int z = 0; z < mSize(2); z++)
			{
				for (int y = 0; y < mSize(1); y++)
				{
					for (int x = 0; x < mSize(0); x++)
					{
						if( IsValidCell(x,y,z) )
						{
							b[i] = 0.0f;

							if( IsValidCell(x-1,y,z) ) b[i] += ( (*this) (x - 1, y, z) (0) + (*this) (x, y, z) (0) ) * 0.5f;
							if( IsValidCell(x+1,y,z) ) b[i] -= ( (*this) (x + 1, y, z) (0) + (*this) (x, y, z) (0) ) * 0.5f;
							if( IsValidCell(x,y-1,z) ) b[i] += ( (*this) (x, y - 1, z) (1) + (*this) (x, y, z) (1) ) * 0.5f;
							if( IsValidCell(x,y+1,z) ) b[i] -= ( (*this) (x, y + 1, z) (1) + (*this) (x, y, z) (1) ) * 0.5f;
							if( IsValidCell(x,y,z-1) ) b[i] += ( (*this) (x, y, z - 1) (2) + (*this) (x, y, z) (2) ) * 0.5f;
							if( IsValidCell(x,y,z+1) ) b[i] -= ( (*this) (x, y, z + 1) (2) + (*this) (x, y, z) (2) ) * 0.5f;
						}
						else
							b[i] = 0.0f;

						i++;
					}
				}
			}


			Grid < float >gradient(mSize);
			mMatrixConserveMass.Solve(gradient.GetVector(), b);

			i = 0;
			for (int z = 0; z < mSize(2); z++)
			{
				for (int y = 0; y < mSize(1); y++)
				{
					for (int x = 0; x < mSize(0); x++)
					{
						if( IsValidCell(x,y,z ) )
						{
							if( !IsValidCell(x-1,y,z) )
								(*target)[i][0] = (*this) (i) (0) - ( gradient(x, y, z) - gradient(x + 1, y, z) ) * 0.5f;
							else if( !IsValidCell(x+1,y,z ) )
								(*target)[i][0] = (*this) (i) (0) - ( gradient(x - 1, y, z) - gradient(x, y, z) ) * 0.5f;
							else
								(*target)[i][0] = (*this) (i) (0) - (gradient(x - 1, y, z) - gradient(x + 1, y, z)) * 0.5f;

							if( !IsValidCell(x,y-1,z) )
								(*target)[i][1] = (*this) (i) (1) - ( gradient(x, y, z) - gradient(x, y + 1, z) ) * 0.5f;
							else if( !IsValidCell(x,y+1,z ) )
								(*target)[i][1] = (*this) (i) (1) - ( gradient(x, y - 1, z) - gradient(x, y, z) ) * 0.5f;
							else
								(*target)[i][1] = (*this) (i) (1) - (gradient(x, y - 1, z) - gradient(x, y + 1, z)) * 0.5f;

							if( !IsValidCell(x,y,z-1) )
								(*target)[i][2] = (*this) (i) (2) - ( gradient(x, y, z) - gradient(x, y, z + 1) ) * 0.5f;
							else if( !IsValidCell(x,y,z+1 ) )
								(*target)[i][2] = (*this) (i) (2) - ( gradient(x, y, z - 1) - gradient(x, y, z) ) * 0.5f;
							else
								(*target)[i][2] = (*this) (i) (2) - (gradient(x, y, z - 1) - gradient(x, y, z + 1)) * 0.5f;
						}
						else
						{
							(*target)[i][0] = 0.0f;
							(*target)[i][1] = 0.0f;
							(*target)[i][2] = 0.0f;
						}

						i++;
					}
				}
			}

			return target;
		}

	/** Add a grid.
	 * Adds the values of another grid of the same size to the values of this grid.
	 * Use this function to add user defined forces or density sources.
	 *
	 * \param g source grid
	 * \return reference to this grid (use to chain operations)
	 */
		Grid < T > &operator+=(const Grid < T > &g)
		{
			assert(mSize == g.GetSize());

			for (int i = 0; i < mSize1D; i++)
				mValue[i] += g(i);

			return *this;
		}


	/** Scale.
	 * Divides all values by a scalar. Use this function for dissipation.
	 *
	 * \param f scaling value
	 * \return reference to this grid (use to chain operations)
	 */
		Grid < T > &operator/=(const float f)
		{
			for (int i = 0; i < mSize1D; i++)
				mValue[i] /= f;

			return *this;
		}

	/** Copy
	 * Create a copy of another grid.
	 *
	 * \param g source grid
	 * \return reference to this grid (use to chain operations)
	 */
		Grid < T > &operator=(const Grid < T > &g)
		{
			assert(mSize == g.GetSize());

			for (int i = 0; i < mSize1D; i++)
				mValue[i] = g(i);

			return *this;
		}

	/** Scale.
	 * Multiplies all values with a scalar. Use this function for time dependent calculations.
	 *
	 * \param f scaling value
	 * \return reference to this grid (use to chain operations)
	 */
		Grid < T > operator*(const float f) const
		{
			Grid < T > g = *this;

			for (int i = 0; i < mSize1D; i++)
				g[i] *= f;

			return g;
		}

	/** Add object forces to the grid.
	 *
	 * \param time time factor, velocity = force * time
	 */
		void AddForce(const float time)
		{
			if (mpObject != NULL)
				for (int z = 0; z < mSize(2); z++)
					for (int y = 0; y < mSize(1); y++)
						for (int x = 0; x < mSize(0); x++)
						{
							Ref(x, y, z) += mpObject->GetForce(VecInt3D(x, y, z)) * time;
						}
		}

	/** Add object sources to the grid.
	 *
	 * \param time  time factor, density = source rate * time
	 */
		void AddSource(const float time)
		{
			if (mpObject != NULL)
				for (int z = 0; z < mSize(2); z++)
					for (int y = 0; y < mSize(1); y++)
						for (int x = 0; x < mSize(0); x++)
						{
							Ref(x, y, z) += mpObject->GetSource(VecInt3D(x, y, z)) * time;
						}
		}

private:
	/** Update changed in boundary conditions.
	 */
		void UpdateBoundary()
		{
			//mMatrixConserveMass.InitValues(-6, 0, 1);
			mMatrixConserveMass.InitValues(0, -1, 1);

			mMatrixConserveMass.InitNeighbour(mpObject, mBounded);
			mMatrixDiffusion.InitNeighbour(mpObject, mBounded);
		}



public:

	/** Use solid walls around the simulation grid.
	 *
	 * \param bounded true if bounded.
	 */
		void SetBounded(const bool bounded)
		{
			mBounded = bounded;
			UpdateBoundary();
		}

	/** Set the object hierarchy used to interact with the simulation.
	 *
	 * \param object pointer to the object hierarchy
	 */
		void SetObject(const Object * object, const bool update_solid)
		{
			mpObject = object;

			if (update_solid)
				UpdateBoundary();
		}


	private:
		VecInt3D mSize;							///< grid size
		int mSize1D;								///< number of vertices (linear size)
		Vector < T > mValue;				///< vertex values
		bool mBounded;							///< bounded grid, uses solid walls around the grid
		const Object *mpObject;			///< object hierarchy, used to exclude solid grid cells and add sources/forces

		Matrix < T > mMatrixDiffusion;	///< cached matrix, used in Diffuse()
		Matrix < float >mMatrixConserveMass;	///< cached matrix, used in ConserveMass()
	};

};															// namespace Fluid

#endif													// __GRID_H__