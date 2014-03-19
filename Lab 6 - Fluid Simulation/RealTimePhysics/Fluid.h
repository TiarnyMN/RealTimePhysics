#ifndef __SOLVER_H__
#define __SOLVER_H__

//http://fluids.cvs.sourceforge.net/viewvc/fluids/

// forward declaration
namespace Fluid
{
	class Solver;
}

#include "config.h"
#include "vector.h"
#include "vector3d.h"
#include "vecint3d.h"
#include "grid.h"
#include "object.h"

namespace Fluid
{

/** The fluid simulation solver. 
 *
 * This is the fluid simulation class. You can specify all simulation 
 * parameters here and do discrete time steps. The system is then solved
 * using the algorithm described in Jos Stam's paper "Stable Fluids".
 * Note that the simulation is influenced by numerical errors. Flows and
 * concentrations will show diffusion and dissipation even when you set
 * the according parameters to zero. You can improve the details by using
 * shorter time steps or higher grid resolutions.
 *
 * \ingroup core
 *
 * \version $Revision: 1.2 $
 *
 * \date $Date: 2003/02/05 15:51:27 $
 *
 * \author malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	class Solver
	{
	public:
	/** Creates and initializes the simulation to the given size.
	 * The simulation size cannot be changed later.
	 * 
	 * \param size Specifies the size of the box in world space. This
	 * is also the space the objects are defined in.
	 *
	 * \param cells Specifies the number of simulation grid cells in three 
	 * dimensions. Do not use more than 10 grid cells per MHz if you want 
	 * to work in realtime (rule of thumb). Set only one coordinate if you
	 * want the others to be calculated automatically. This way, you can
	 * be sure to get cubic voxels which are required for correct simulations.
	 */
		Solver(const Vector3D size, const VecInt3D cells);
		~Solver();

	/** Sets the viscosity of the fluid.
	 * The viscosity describes the diffusion of the velocity inside the fluid. 
	 * A low viscosity means that the fluid has a low internal friction, so
	 * flows are more straight forward and to not spread.
	 * 
	 * \param viscosity positive float value
	 */
		void SetViscosity(float viscosity);
		float GetViscosity();

	/** Set the diffusion of the fluid.
	 * The diffusion describes the rate of distribution of the particles inside
	 * the fluid. A low diffusion means that the particles will stay in their
	 * places, keeping concentration peaks. A high diffusion causes the
	 * particles to distribute equally among the fluid.
	 *
	 * \param diffusion non-negative scalar
	 */
		void SetDiffusion(float diffusion);
		float GetDiffusion();

	/** Sets the dissipation of the fluid
	 * The dissipation describes the dampening factor of the particle
	 * concentration. A high dissipation rate means that the particles
	 * disappear soon whereas a low rate causes them to stay quite long.
	 * This parameter is used to keep control over the system. It's not
	 * required by natural phenomena, but a simulation without particle
	 * dissipation is quite hard to maintain from the artistic side.
	 *
	 * \param dissipation non-negative scalar
	 */
		void SetDissipation(float dissipation);
		float GetDissipation();

	/** Does one step ahead in time.
	 * This function calculates the changes of the simulation after one
	 * discrete time step of the specified size. Use it in your main
	 * loop when doing realtime graphics:
	 * while( true )
	 * {
	 *   // process user input
	 *   solver.DoStep( delta_time );
	 *   // update visualization
     * }
	 *
	 * \param time The time difference to the last time DoStep was called in seconds
	 */
		void DoStep(float time);

	/** Gets the density grid.
	 * This function returns a const pointer to the density grid. Use it to
	 * visualize the simulation if you want to show the particles flowing
	 * through the fluid (added by SetSource()).
	 * DO NOT assume that this pointer does not change after a call to DoStep().
	 *
	 * \return The density grid
	 */
		const Grid < float >*GetDensity() const;
	/** Gets the velocity grid
	 * This function returns a const pointer to the velocity grid. Use it to
	 * visualize the simulation if you want to show streamlines of the flow
	 * itself. If you prefer showing the particles instead, use GetDensity().
	 * DO NOT assume that this pointer does not change after a call to DoStep().
	 * 
	 * \return  The velocity grid
	 */
		const Grid < Vector3D > *GetVelocity() const;

	/** Use solid walls around the simulation grid.
	 * 
	 * \param bounded true if bounded
	 */
		void SetBoundedGrid(const bool bounded);

	/** Set the object hierarchy used to interact with the simulation.
	 * If you don't change the solid structure, set update_solid to false.
	 * This way the solver can omit the update of the solid matrix setup
	 * which gives a noticable performance boost.
	 * 
	 * \param object pointer to the object hierarchy
	 * \param update_solid true if the solid objects changed
	 */
		void SetObject(const Object * object, const bool update_solid);

	/** Convert world space vector to grid space.
	 * 
	 * \param v world space vector
	 * \return grid space vector
	 */
		Vector3D GridSpace(const Vector3D v) const;

	/** Convert world space vector to grid space.
	 * 
	 * \param v world space vector
	 * \return grid space vector
	 */
		VecInt3D GridSpaceInt(const Vector3D v) const;

	/** Convert world space scalar to grid space.
	 * 
	 * \param v world space scalar
	 * \return grid space scalar
	 */
		float GridSpace(const float f) const;

	/** Convert world space scalar to grid space.
	 * 
	 * \param v world space scalar
	 * \return grid space scalar
	 */
		int GridSpaceInt(const float f) const;

	private:
		float mViscosity;						///< See SetViscosity()
		float mDiffusion;						///< See SetDiffusion()
		float mDissipation;					///< See SetDissipation()

		int mSource;								///< Used to switch between the two buffers while calculating a new time step
		int mDest;									///< Used to switch between the two buffers while calculating a new time step

		 Grid < float >*mDensity[2];	///< The density grid now and one time step before
		 Grid < Vector3D > *mVelocity[2];	///< The velocity grid now and one time step before

		Vector3D mSize;							///< simulation size in world space
		VecInt3D mCells;						///< number of grid cells
	};

}

#endif													// __SOLVER_H__


//#pragma once
//#include <math.h>
//#include <stdlib.h>
//
////http://www.mat.ucsb.edu/~wakefield/594cm/assignment.htm
//
///* Look at the below source, and the elements from the initial tutorial.
// * Then look at the code linked to by John, find and extract the common elements. */
//
//struct Fluid
//{
//	unsigned int bits, dim, dim2, dim3, dimwrap;		//Discrete field dimensions.
//	unsigned int scale;	// division of region dimension to fluid dimension (pow of 2)
//
//	double *pfield, *pfield0;		//Pressure fields.
//	double *dfield, *dfield0;		//Density fields.
//	double *uxfield, *uxfield0;		//Velocity fields.
//	double *uyfield, *uyfield0;
//	double *uzfield, *uzfield0;
//
//	double *noisefield;			//A source of turbulence.
//	double noise;
//
//	double viscosity, diffusion, decay;		//Global fluid parameters for the cube.
//
//	Fluid(unsigned int dim, double viscosity = 0.00001, double diffusion = 0.001, double decay = 0.99)
//		: dim (dim), viscosity(viscosity), diffusion(diffusion), decay(decay)
//	{
//		dim2 = dim * dim;			//dim squared.
//		dim3 = dim * dim * dim;		//dim cube.
//		dimwrap = dim - 1;
//
//		scale = 1;
//
//		//dim = DIMENSION
//
//		noisefield = new double[dim3];
//		noise = 1;
//
//		//Initialising all the arrays.
//		dfield = new double[dim3];
//
//		pfield = new double[dim3];
//		pfield0 = new double[dim3];
//
//		uxfield = new double[dim3];
//		uyfield = new double[dim3];
//		uzfield = new double[dim3];
//		dfield0 = new double[dim3];
//		uxfield0 = new double[dim3];
//		uyfield0 = new double[dim3];
//		uzfield0 = new double[dim3];
//
//		reset();
//
//		////Initialising the density field to 0 "empty" initially,
//		////and adding random velocities to the velocity fields.
//		//for (int i=0; i<dim3; i++) 
//		//{
//		//	dfield0[i] = dfield[i] = 0;
//		//	uxfield0[i] = uxfield[i] = ((double) rand() / (RAND_MAX));
//		//	uyfield0[i] = uyfield[i] = ((double) rand() / (RAND_MAX));
//		//	uzfield0[i] = uzfield[i] = ((double) rand() / (RAND_MAX));
//		//}
//	}
//
//	~Fluid()
//	{
//		delete[] dfield;
//		delete[] uxfield;
//		delete[] uyfield;
//		delete[] uzfield;
//		delete[] dfield0;
//		delete[] uxfield0;
//		delete[] uyfield0;
//		delete[] uzfield0;
//	}
//
//	//inline void addDensity(glm::vec3 &v, double d = 1)
//	//{
//	//	dfield[index(v.x, v.y, v.z)] = d;
//	//}
//
//	//inline void addVelocity(glm::vec3 &v, glm::vec3 &d)
//	//{
//	//	uxfield[index(v.x, v.y, v.z)] += d.x;
//	//	uyfield[index(v.x, v.y, v.z)] += d.y;
//	//	uzfield[index(v.x, v.y, v.z)] += d.z;
//	//}
//
//	//inline void advet(double *f, double *f0, double *ux, double *uy, double *uz, double rate = 1)
//	//{
//	//	for(int x = 0; x < dim; x++)
//	//	{
//	//		for(int y = 0; y < dim; y++)
//	//		{
//	//			for(int z = 0; z < dim; z++)
//	//			{
//	//				int i = index(x, y, z);
//
//	//				//Get co-ordinates back in time from current cell:
//	//				double x0 = x - rate * ux[i];
//	//				double y0 = y - rate * uy[i];
//	//				double z0 = z - rate * uz[i];
//
//	//				//Get interpolated value at x0, y0, z0:
//	//				f[i] = interp(f0, x0, y0, z0);
//	//			}
//	//		}
//	//	}
//	//}
//
//	//inline void project(int N, double *vx, double *vy, double *vz, double *g, double *g0)
//	//{
//	//	double h = 1.0 / N;
//
//	//	// g and g0 represent the velocity gradient
//	//	BEGIN_PER_CELL
//	//		// previous instantaneous magnitude of velocity gradient 
//	//		//		= (sum of velocity gradients per axis)/2N:
//	//		g0[i] = -0.5 * h * (
//	//								vx[INDEX(x+1, y, z)]-vx[INDEX(x-1, y, z)]+	// velocity gradients
//	//								vy[INDEX(x, y+1, z)]-vy[INDEX(x, y-1, z)]+	// velocity gradients
//	//								vz[INDEX(x, y, z+1)]-vz[INDEX(x, y, z-1)]	// velocity gradients
//	//							);
//	//		// zero out the present velocity gradient:
//	//		g[i] = 0;
//	//	END_PER_CELL
//
//	//	// reuse the Gauss-Seidel relaxation solver to safely diffuse the velocity gradients from g0 to g:
//	//	stable_solve(N, g, g0, 1, 6);
//	//
//	//	// now subtract this gradient from our current velocity field:
//	//	BEGIN_PER_CELL
//	//		vx[i] -= 0.5 * N * (g[INDEX(x+1, y, z)]-g[INDEX(x-1, y, z)]); // gradient calculated by neighbors
//	//		vy[i] -= 0.5 * N * (g[INDEX(x, y+1, z)]-g[INDEX(x, y-1, z)]);
//	//		vz[i] -= 0.5 * N * (g[INDEX(x, y, z+1)]-g[INDEX(x, y, z-1)]);
//	//	END_PER_CELL
//	//}
//
//	inline unsigned int index(int x, int y, int z) 
//	{
//		int xVal = ((x/scale) % dimwrap) * dim2;
//		int yVal = ((y/scale) % dimwrap) * dim;
//		int zVal = ((z/scale) % dimwrap);
//
//		return xVal + yVal + zVal;
//
//		//return ((x/scale) & dimwrap) * dim2 
//		//	 + ((y/scale) & dimwrap) * dim 
//		//	 + ((z/scale) & dimwrap);
//	}
//
//	inline unsigned int index(double x, double y, double z) 
//	{
//		return index((int)floor(x), (int)floor(y), (int)floor(z));
//	}
//
//	inline double fRand(double fMin, double fMax)
//	{
//		double f = (double)rand() / RAND_MAX;
//		return fMin + f * (fMax - fMin);
//	}
//
//	inline void reset()
//	{
//		for(int i = 0; i < dim3; i++)
//		{
//			dfield0[i] = dfield[i] = 0;
//			pfield0[i] = pfield[i] = 0;
//			
//			noisefield[i] = 0;
//
//			//uxfield0[i] = uxfield[i] = 0.0;
//			//uyfield0[i] = uyfield[i] = 0.0;
//			//uzfield0[i] = uzfield[i] = 0.0;
//			
//			uxfield0[i] = uxfield[i] = fRand(-1.0, 1.0);
//			uyfield0[i] = uyfield[i] = fRand(-1.0, 1.0);
//			uzfield0[i] = uzfield[i] = fRand(-1.0, 1.0);
//
//			uxfield0[i] = uxfield[i] = 1.0;
//			uyfield0[i] = uyfield[i] = 0.0;
//			uzfield0[i] = uzfield[i] = 0.0;
//		}
//	}
//
//	inline void getFlow(glm::vec3 &pos, glm::vec3 &flow)
//	{
//		int i = index(pos.x, pos.y, pos.z);
//		flow.x = uxfield[i];
//		flow.y = uyfield[i];
//		flow.z = uzfield[i];
//	}
//
//	inline void getFlowInterp(glm::vec3 &pos, glm::vec3 &flow)
//	{
//		double x = pos.x/(double)scale; 
//		double y = pos.y/(double)scale; 
//		double z = pos.z/(double)scale;
//		
//		// get the integer components, and normalized 0..1 interp factors, of x,y,z:
//		int x0 = (int)x; double xbf = x-(double)x0; double xaf = 1.-xbf; 
//		int y0 = (int)y; double ybf = y-(double)y0; double yaf = 1.-ybf;
//		int z0 = (int)z; double zbf = z-(double)z0; double zaf = 1.-zbf;
//		
//		// find the 8 cube corners by getting the corners' array indices:
//		int xa = (x0 % (dimwrap))*dim2;	int xb = ((x0+1) % (dimwrap))*dim2;
//		int ya = (y0 % (dimwrap))*dim;	int yb = ((y0+1) % (dimwrap))*dim;
//		int za = (z0 % (dimwrap));		int zb = ((z0+1) % (dimwrap));
//		
//		int iaaa = xa + ya + za; double faaa = xaf * yaf * zaf;
//		int ibaa = xb + ya + za; double fbaa = xbf * yaf * zaf;
//		int iaba = xa + yb + za; double faba = xaf * ybf * zaf;
//		int iaab = xa + ya + zb; double faab = xaf * yaf * zbf;
//		int ibab = xb + ya + zb; double fbab = xbf * yaf * zbf;
//		int iabb = xa + yb + zb; double fabb = xaf * ybf * zbf;
//		int ibba = xb + yb + za; double fbba = xbf * ybf * zaf;
//		int ibbb = xb + yb + zb; double fbbb = xbf * ybf * zbf;
//		
//		// do the interpolation:
//		flow.x =(uxfield[iaaa] * faaa) +
//				(uxfield[ibaa] * fbaa) + 
//				(uxfield[iaba] * faba) + 
//				(uxfield[iaab] * faab) +
//				(uxfield[ibab] * fbab) + 
//				(uxfield[iabb] * fabb) + 
//				(uxfield[ibba] * fbba) + 
//				(uxfield[ibbb] * fbbb);
//				
//		flow.y =(uyfield[iaaa] * faaa) +
//				(uyfield[ibaa] * fbaa) + 
//				(uyfield[iaba] * faba) + 
//				(uyfield[iaab] * faab) +
//				(uyfield[ibab] * fbab) + 
//				(uyfield[iabb] * fabb) + 
//				(uyfield[ibba] * fbba) + 
//				(uyfield[ibbb] * fbbb);
//				
//		flow.z =(uzfield[iaaa] * faaa) +
//				(uzfield[ibaa] * fbaa) + 
//				(uzfield[iaba] * faba) + 
//				(uzfield[iaab] * faab) +
//				(uzfield[ibab] * fbab) + 
//				(uzfield[iabb] * fabb) + 
//				(uzfield[ibba] * fbba) + 
//				(uzfield[ibbb] * fbbb);
//	}
//
//	inline void mixFlow(glm::vec3 &pos, glm::vec3 vel, double mix)
//	{
//		int i = index(pos.x, pos.y, pos.z);
//		uxfield[i] += mix * (2.*vel.x - uxfield0[i]);
//		uyfield[i] += mix * (2.*vel.y - uyfield0[i]);
//		uzfield[i] += mix * (2.*vel.z - uzfield0[i]);
//	}
//
//	inline void addNoise(glm::vec3 &pos, double amt)
//	{
//		int i = index(pos.x, pos.y, pos.z);
//		pfield[i] += amt;
//	}
//
//	inline void addNoise(float x, float y, float z, double amt)
//	{
//		int i = index(x, y, z);
//		pfield[i] += amt;
//	}
//
//	inline double getNoise(glm::vec3 &pos)
//	{
//		int i = index(pos.x, pos.y, pos.z);
//		return pfield[i];
//	}
//
//	inline double getNoise(int x, int y, int z)
//	{
//		int i = index(x, y, z);
//		//cout << i << " - " << pfield[i] << endl;
//		return pfield[i];
//	}
//
//	inline void step();
//};
//
//namespace stable {
//
//	#define INDEX(x, y, z) ((z % (N-1)) + (y % (N-1))*N + (x % (N-1))*N*N)
//	#define BEGIN_PER_CELL	for (int x=0 ; x<N ; x++ ) { \
//							for (int y=0 ; y<N ; y++ ) { \
//							for (int z=0 ; z<N ; z++ ) { \
//								int i = INDEX(x, y, z);
//	#define END_PER_CELL	}}}
//	#define SWAP_PTR(x0,x) {double * tmp=x0;x0=x;x=tmp;}
//
//	/*
//		Stable diffusion solver
//			for dimension N, diffuses field p (previous p0)
//		
//		This method avoids oscillation and divergence in diffusion, 
//			granting a *stable* result instead:
//			an iterative technique to solve a sparse linear system
//			(the Gauss-Seidel relaxation method)	
//	*/
//
//	inline void set_bnd(int b, double *x, int N)
//	{
//		for(int j = 1; j < N - 1; j++) {
//			for(int i = 1; i < N - 1; i++) {
//				x[INDEX(i, j, 0  )] = b == 3 ? -x[INDEX(i, j, 1  )] : x[INDEX(i, j, 1  )];
//				x[INDEX(i, j, N-1)] = b == 3 ? -x[INDEX(i, j, N-2)] : x[INDEX(i, j, N-2)];
//			}
//		}
//		for(int k = 1; k < N - 1; k++) {
//			for(int i = 1; i < N - 1; i++) {
//				x[INDEX(i, 0  , k)] = b == 2 ? -x[INDEX(i, 1  , k)] : x[INDEX(i, 1  , k)];
//				x[INDEX(i, N-1, k)] = b == 2 ? -x[INDEX(i, N-2, k)] : x[INDEX(i, N-2, k)];
//			}
//		}
//		for(int k = 1; k < N - 1; k++) {
//			for(int j = 1; j < N - 1; j++) {
//				x[INDEX(0  , j, k)] = b == 1 ? -x[INDEX(1  , j, k)] : x[INDEX(1  , j, k)];
//				x[INDEX(N-1, j, k)] = b == 1 ? -x[INDEX(N-2, j, k)] : x[INDEX(N-2, j, k)];
//			}
//		}
//    
//		x[INDEX(0, 0, 0)]       = 0.33 * (x[INDEX(1, 0, 0)]
//										+ x[INDEX(0, 1, 0)]
//										+ x[INDEX(0, 0, 1)]);
//		x[INDEX(0, N-1, 0)]     = 0.33 * (x[INDEX(1, N-1, 0)]
//										+ x[INDEX(0, N-2, 0)]
//										+ x[INDEX(0, N-1, 1)]);
//		x[INDEX(0, 0, N-1)]     = 0.33 * (x[INDEX(1, 0, N-1)]
//										+ x[INDEX(0, 1, N-1)]
//										+ x[INDEX(0, 0, N)]);
//		x[INDEX(0, N-1, N-1)]   = 0.33 * (x[INDEX(1, N-1, N-1)]
//										+ x[INDEX(0, N-2, N-1)]
//										+ x[INDEX(0, N-1, N-2)]);
//		x[INDEX(N-1, 0, 0)]     = 0.33 * (x[INDEX(N-2, 0, 0)]
//										+ x[INDEX(N-1, 1, 0)]
//										+ x[INDEX(N-1, 0, 1)]);
//		x[INDEX(N-1, N-1, 0)]   = 0.33 * (x[INDEX(N-2, N-1, 0)]
//										+ x[INDEX(N-1, N-2, 0)]
//										+ x[INDEX(N-1, N-1, 1)]);
//		x[INDEX(N-1, 0, N-1)]   = 0.33 * (x[INDEX(N-2, 0, N-1)]
//										+ x[INDEX(N-1, 1, N-1)]
//										+ x[INDEX(N-1, 0, N-2)]);
//		x[INDEX(N-1, N-1, N-1)] = 0.33 * (x[INDEX(N-2, N-1, N-1)]
//										+ x[INDEX(N-1, N-2, N-1)]
//										+ x[INDEX(N-1, N-1, N-2)]);
//	}
//
//	inline void stable_solve(int b, int N, double * p, double * p0, double diffusion, double divisor ) {
//		static int iterations = 16; //20; 
//		double div = 1.0/divisor;
//
//		for(int n = 0; n < iterations; n++)
//		{
//			for(int x = 1; x < N-1; x++)
//			{
//				for(int y = 1; y < N-1; y++)
//				{
//					for(int z = 1; z < N-1; z++)
//					{
//						int i = INDEX(x, y, z);
//						p[i] = (
//								p0[i] + 
//								diffusion*(
//									p[INDEX(x-1, y, z)] + p[INDEX(x+1, y, z)] +
//									p[INDEX(x, y-1, z)] + p[INDEX(x, y+1, z)] +
//									p[INDEX(x, y, z-1)] + p[INDEX(x, y, z+1)]
//								)
//							) * div;
//					}
//				}
//			}
//		}
//
//		set_bnd(b, p, N);
//
//
//
//	//	#define BEGIN_PER_CELL	for (int x=0 ; x<N ; x++ ) { \
//	//						for (int y=0 ; y<N ; y++ ) { \
//	//						for (int z=0 ; z<N ; z++ ) { \
//	//							int i = INDEX(x, y, z);
//	//#define END_PER_CELL	}}}
//	//	
//	//	for (int n=0 ; n<iterations ; n++) {	
//	//		BEGIN_PER_CELL	
//	//			p[i] = (
//	//						p0[i] + 
//	//						diffusion*(
//	//							p[INDEX(x-1, y, z)] + p[INDEX(x+1, y, z)] +
//	//							p[INDEX(x, y-1, z)] + p[INDEX(x, y+1, z)] +
//	//							p[INDEX(x, y, z-1)] + p[INDEX(x, y, z+1)]
//	//						)
//	//					) * div;
//	//		END_PER_CELL
//	//	}
//	}
//
//	inline void diffuse(int b, int N, double * p, double * p0, double diffusion)
//	{
//		double a = diffusion; // * (N * N * N);
//		stable_solve(b, N, p, p0, a, (1+6*a) );
//	}
//
//	// (linear) interpolate a value at the given coords:
//	inline double interp(int N, double * p0, double x, double y, double z) 
//	{
//		int dimwrap = N-1; // assumes N is a power of 2!!
//		
//		// get the integer components, and normalized 0..1 interp factors, of x,y,z:
//		int x0 = (int)x; double xbf = x-(double)x0; double xaf = 1-xbf; 
//		int y0 = (int)y; double ybf = y-(double)y0; double yaf = 1-ybf;
//		int z0 = (int)z; double zbf = z-(double)z0; double zaf = 1-zbf;
//		
//		// find the 8 cube corners by getting the corners' array indices:
//		int xa = (x0 % (dimwrap))*N*N;	int xb = ((x0+1) % (dimwrap))*N*N;
//		int ya = (y0 % (dimwrap))*N;	int yb = ((y0+1) % (dimwrap))*N;
//		int za = (z0 % (dimwrap));		int zb = ((z0+1) % (dimwrap));
//		
//		// do the interpolation:
//		return  (p0[xa + ya + za] * xaf * yaf * zaf) +
//				(p0[xb + ya + za] * xbf * yaf * zaf) + 
//				(p0[xa + yb + za] * xaf * ybf * zaf) + 
//				(p0[xa + ya + zb] * xaf * yaf * zbf) +
//				(p0[xb + ya + zb] * xbf * yaf * zbf) + 
//				(p0[xa + yb + zb] * xaf * ybf * zbf) + 
//				(p0[xb + yb + za] * xbf * ybf * zaf) + 
//				(p0[xb + yb + zb] * xbf * ybf * zbf);
//	}
//
//	inline void advect(int b, int N, double *d, double *d0,  double *velocX, double *velocY, double *velocZ, double dt)
//	{
//		double i0, i1, j0, j1, k0, k1;
//    
//		double dtx = dt * (N - 2);
//		double dty = dt * (N - 2);
//		double dtz = dt * (N - 2);
//    
//		double s0, s1, t0, t1, u0, u1;
//		double tmp1, tmp2, tmp3, x, y, z;
//    
//		double Nfloat = N;
//		double ifloat, jfloat, kfloat;
//		int i, j, k;
//    
//		for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
//			for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
//				for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
//					tmp1 = dtx * velocX[INDEX(i, j, k)];
//					tmp2 = dty * velocY[INDEX(i, j, k)];
//					tmp3 = dtz * velocZ[INDEX(i, j, k)];
//					x    = ifloat - tmp1; 
//					y    = jfloat - tmp2;
//					z    = kfloat - tmp3;
//                
//					if(x < 0.5f) x = 0.5f; 
//					if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
//					i0 = floorf(x); 
//					i1 = i0 + 1.0f;
//					if(y < 0.5f) y = 0.5f; 
//					if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
//					j0 = floorf(y);
//					j1 = j0 + 1.0f; 
//					if(z < 0.5f) z = 0.5f;
//					if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
//					k0 = floorf(z);
//					k1 = k0 + 1.0f;
//                
//					s1 = x - i0; 
//					s0 = 1.0f - s1; 
//					t1 = y - j0; 
//					t0 = 1.0f - t1;
//					u1 = z - k0;
//					u0 = 1.0f - u1;
//                
//					int i0i = i0;
//					int i1i = i1;
//					int j0i = j0;
//					int j1i = j1;
//					int k0i = k0;
//					int k1i = k1;
//                
//					d[INDEX(i, j, k)] = 
//                
//						s0 * ( t0 * (u0 * d0[INDEX(i0i, j0i, k0i)]
//									+u1 * d0[INDEX(i0i, j0i, k1i)])
//							+( t1 * (u0 * d0[INDEX(i0i, j1i, k0i)]
//									+u1 * d0[INDEX(i0i, j1i, k1i)])))
//					   +s1 * ( t0 * (u0 * d0[INDEX(i1i, j0i, k0i)]
//									+u1 * d0[INDEX(i1i, j0i, k1i)])
//							+( t1 * (u0 * d0[INDEX(i1i, j1i, k0i)]
//									+u1 * d0[INDEX(i1i, j1i, k1i)])));
//				}
//			}
//		}
//		set_bnd(b, d, N);
//	}
//
//	/*
//		Translate densities etc. (p, p0 = previous) through the vector field (vx, vy, vz) over dim N
//			Does a linear back-trace of center point through vector field, 
//			using linear interpolation to gather p0 at this source and place in new p
//	*/
//	//inline void advect(int b, int N, double * p, double * p0, double * vx, double * vy, double * vz)
//	//{
//	//	static double rate = 1.0;
//	//	for(int x = 1; x < N-1; x++)
//	//	{
//	//		for(int y = 1; y < N-1; y++)
//	//		{
//	//			for(int z = 1; z < N-1; z++)
//	//			{
//	//				int i = INDEX(x, y, z);
//	//				//back trace:
//	//				double x0 = x - rate * vx[i]; 
//	//				double y0 = y - rate * vy[i]; 
//	//				double z0 = z - rate * vz[i]; 
//	//				// trilinearinterp p0 at (x0, y0, z0):
//	//				p[i] = interp(N, p0, x0, y0, z0);
//	//			}
//	//		}
//	//	}
//
//	//	set_bnd(b, p, N);
//
//	//	//BEGIN_PER_CELL
//	//	//	// back trace:
//	//	//	double x0 = x - rate * vx[i]; 
//	//	//	double y0 = y - rate * vy[i]; 
//	//	//	double z0 = z - rate * vz[i]; 
//	//	//	// trilinearinterp p0 at (x0, y0, z0):
//	//	//	p[i] = interp(N, p0, x0, y0, z0);
//	//	//END_PER_CELL
//	//}
//
//	
//
//	/*
//		Clever part of Jos Stam's work. 
//			A velocity field can become divergent (have regions that are purely emanating or aggregating)
//				violating the definition of an incompressible fluid
//			But, since a velocity field can be seen as an incompressible velocity field + a gradient field,
//				we can subtract a gradient field from our bad velocity field to get an incompressible one
//			To calculate this gradient field and then subtract it, we use this function:
//	*/
//	inline void project(int N, double * vx, double * vy, double * vz, double * g, double * g0 )
//	{
//		double h = 1.0/N;
//
//		for(int x = 1; x < N-1; x++)
//		{
//			for(int y = 1; y < N-1; y++)
//			{
//				for(int z = 1; z < N-1; z++)
//				{
//					int i = INDEX(x, y, z);
//					g0[i] = -0.5 * h * (
//									vx[INDEX(x+1, y, z)]-vx[INDEX(x-1, y, z)]+	// velocity gradients
//									vy[INDEX(x, y+1, z)]-vy[INDEX(x, y-1, z)]+	// velocity gradients
//									vz[INDEX(x, y, z+1)]-vz[INDEX(x, y, z-1)]	// velocity gradients
//								);
//					// zero out the present velocity gradient:
//					g[i] = 0;
//				}
//			}
//		}
//		
//		//// g and g0 represent the velocity gradient
//		//BEGIN_PER_CELL
//		//	// previous instantaneous magnitude of velocity gradient 
//		//	//		= (sum of velocity gradients per axis)/2N:
//		//	g0[i] = -0.5 * h * (
//		//							vx[INDEX(x+1, y, z)]-vx[INDEX(x-1, y, z)]+	// velocity gradients
//		//							vy[INDEX(x, y+1, z)]-vy[INDEX(x, y-1, z)]+	// velocity gradients
//		//							vz[INDEX(x, y, z+1)]-vz[INDEX(x, y, z-1)]	// velocity gradients
//		//						);
//		//	// zero out the present velocity gradient:
//		//	g[i] = 0;
//		//END_PER_CELL
//		
//		// reuse the Gauss-Seidel relaxation solver to safely diffuse the velocity gradients from g0 to g:
//		set_bnd(0, g, N); 
//		set_bnd(0, g0, N);
//		stable_solve(0, N, g, g0, 1, 6);		
//
//		for(int x = 1; x < N-1; x++)
//		{
//			for(int y = 1; y < N-1; y++)
//			{
//				for(int z = 1; z < N-1; z++)
//				{
//					int i = INDEX(x, y, z);
//
//					vx[i] -= 0.5 * N * (g[INDEX(x+1, y, z)]-g[INDEX(x-1, y, z)]); // gradient calculated by neighbors
//					vy[i] -= 0.5 * N * (g[INDEX(x, y+1, z)]-g[INDEX(x, y-1, z)]);
//					vz[i] -= 0.5 * N * (g[INDEX(x, y, z+1)]-g[INDEX(x, y, z-1)]);
//				}
//			}
//		}
//
//		//// now subtract this gradient from our current velocity field:
//		//BEGIN_PER_CELL
//		//	vx[i] -= 0.5 * N * (g[INDEX(x+1, y, z)]-g[INDEX(x-1, y, z)]); // gradient calculated by neighbors
//		//	vy[i] -= 0.5 * N * (g[INDEX(x, y+1, z)]-g[INDEX(x, y-1, z)]);
//		//	vz[i] -= 0.5 * N * (g[INDEX(x, y, z+1)]-g[INDEX(x, y, z-1)]);
//		//END_PER_CELL
//	}
//
//	void dens_step (int N, double * p, double * p0, double * vx, double * vy, double * vz, double diffusion )
//	{
//		SWAP_PTR( p0, p ); 
//		diffuse(0, N, p, p0, diffusion);
//		SWAP_PTR( p0, p ); 
//		//advect(0, N, p, p0, vx, vy, vz, 1.0);
//		advect(0, N, p, p0, vx, vy, vz, 1.0);
//	}
//
//	//    diffuse(1, Vx0, Vx, visc, dt, 4, N);	//Diffuse in the x direction.
////    diffuse(2, Vy0, Vy, visc, dt, 4, N);	//Diffuse in the y direction.
////    diffuse(3, Vz0, Vz, visc, dt, 4, N);	//Diffuse in the z direction.
////    
////    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
////    
////    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
////    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
////    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
////    
////    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
//
//	void vel_step (int N, double * vx, double * vy, double * vz, double * vx0, double * vy0, double * vz0, double viscosity)
//	{
//		// diffuse the velocity field (per axis):
//		SWAP_PTR ( vx0, vx ); 
//		diffuse (1, N, vx, vx0, viscosity );
//		SWAP_PTR ( vy0, vy ); 
//		diffuse (2, N, vy, vy0, viscosity );
//		SWAP_PTR ( vz0, vz ); 
//		diffuse (3, N, vz, vz0, viscosity );
//		// stabilize it: (vx0, vy0 are whatever, being used as temporaries to store gradient field)
//		project ( N, vx, vy, vz, vx0, vy0 );
//		
//		// advect the velocity field (per axis):
//		SWAP_PTR (vx0, vx); 
//		SWAP_PTR (vy0, vy);
//		SWAP_PTR (vz0, vz);
//		advect (1, N, vx, vx0, vx0, vy0, vz0, 1.0 ); 
//		advect (2, N, vy, vy0, vx0, vy0, vz0, 1.0 );
//		advect (3, N, vz, vz0, vx0, vy0, vz0, 1.0 );
////		// stabilize it: (vx0, vy0 are whatever, being used as temporaries to store gradient field)
//		project ( N, vx, vy, vz, vx0, vy0 );
//	}
//
//	#undef BEGIN_PER_CELL
//	#undef END_PER_CELL
//	#undef INDEX
//	#undef SWAP_PTR
//
//} // namespace
//
//inline void Fluid :: step() 
//{		
//	// diffuse velocity:
//	stable::vel_step(dim, 
//		uxfield0, uyfield0, uzfield0, 
//		uxfield, uyfield, uzfield, 		
//		viscosity);
//
//		// diffuse pressure:
//	stable::dens_step(dim, 
//		pfield, pfield0, 
//		uxfield, uyfield, uzfield, 
//		diffusion );
//		
//	// distort velocity:
//	for (int i=0; i<dim3; i++) 
//	{
//		//Change happened here!
//		/*uxfield[i] += fRand(-1.0, 1.0) * noise;
//		uyfield[i] += fRand(-1.0, 1.0) * noise;
//		uzfield[i] += fRand(-1.0, 1.0) * noise;*/
//
//		//uxfield[i] += fRand(-1.0, 1.0) * noise;
//		uyfield[i] += 0.1;
//		//uzfield[i] += fRand(-1.0, 1.0) * noise;
//	}
//	//	
//	// decay velocity:
//	for (int i=0; i<dim3; i++) 
//	{
//		uxfield[i] *= decay;
//		uyfield[i] *= decay;
//		uzfield[i] *= decay;
//	}
//
//	//
//	// decay noise:
//	//for (int i=0; i<dim3; i++)
//	//{
//	//	pfield[i] *= 0.9;
//	//}
// }
//
//
////#pragma region Original
////#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)
////
////struct FluidCube
////{
////    int size;		//The size of the cube, e.g. size * size * size
////    float dt;		//The delta time to advance the state by during each update.
////    float diff;		//The diffusion amount (?) :: How fast stuff spreads out in the fluid.
////    float visc;		//The viscosity of the cube (?) :: How thick the fluid is.
////    
////    float* s;			//Scratch space array for density. (Keeps the old / previous values while the new ones are computed.)
////    float* density;		//Density array.
////    
////    float* Vx;			//Velocity arrays for X, Y and Z.
////    float* Vy;
////    float* Vz;
////
////    float* Vx0;			//Scratch space arrays for X, Y and Z.
////    float* Vy0;
////    float* Vz0;
////};
////
////
////FluidCube* FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
////{
////	//Allocates memory for the cube, and creates arrays of size * size * size for each value.
////    FluidCube* cube = (FluidCube*)malloc(sizeof(*cube));
////    int N = size;
////    
////    cube->size = size;
////    cube->dt = dt;
////    cube->diff = diffusion;
////    cube->visc = viscosity;
////    
////    cube->s = (float*)calloc(N * N * N, sizeof(float));
////    cube->density = (float*)calloc(N * N * N, sizeof(float));
////    
////    cube->Vx = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vy = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vz = (float*)calloc(N * N * N, sizeof(float));
////    
////    cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vz0 = (float*)calloc(N * N * N, sizeof(float));
////    
////    return cube;
////}
////
////void FluidCubeFree(FluidCube *cube)
////{
////    free(cube->s);
////    free(cube->density);
////    
////    free(cube->Vx);
////    free(cube->Vy);
////    free(cube->Vz);
////    
////    free(cube->Vx0);
////    free(cube->Vy0);
////    free(cube->Vz0);
////    
////    free(cube);
////}
////
////
////static void set_bnd(int b, float *x, int N)
////{
////    for(int j = 1; j < N - 1; j++) {
////        for(int i = 1; i < N - 1; i++) {
////            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
////            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
////        }
////    }
////    for(int k = 1; k < N - 1; k++) {
////        for(int i = 1; i < N - 1; i++) {
////            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
////            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
////        }
////    }
////    for(int k = 1; k < N - 1; k++) {
////        for(int j = 1; j < N - 1; j++) {
////            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
////            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
////        }
////    }
////    
////    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
////                                  + x[IX(0, 1, 0)]
////                                  + x[IX(0, 0, 1)]);
////    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
////                                  + x[IX(0, N-2, 0)]
////                                  + x[IX(0, N-1, 1)]);
////    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
////                                  + x[IX(0, 1, N-1)]
////                                  + x[IX(0, 0, N)]);
////    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
////                                  + x[IX(0, N-2, N-1)]
////                                  + x[IX(0, N-1, N-2)]);
////    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
////                                  + x[IX(N-1, 1, 0)]
////                                  + x[IX(N-1, 0, 1)]);
////    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
////                                  + x[IX(N-1, N-2, 0)]
////                                  + x[IX(N-1, N-1, 1)]);
////    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
////                                  + x[IX(N-1, 1, N-1)]
////                                  + x[IX(N-1, 0, N-2)]);
////    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
////                                  + x[IX(N-1, N-2, N-1)]
////                                  + x[IX(N-1, N-1, N-2)]);
////}
////
////static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
////{
////    float cRecip = 1.0 / c;
////    for (int k = 0; k < iter; k++) {
////        for (int m = 1; m < N - 1; m++) {
////            for (int j = 1; j < N - 1; j++) {
////                for (int i = 1; i < N - 1; i++) {
////                    x[IX(i, j, m)] =
////                        (x0[IX(i, j, m)]
////                            + a*(    x[IX(i+1, j  , m  )]
////                                    +x[IX(i-1, j  , m  )]
////                                    +x[IX(i  , j+1, m  )]
////                                    +x[IX(i  , j-1, m  )]
////                                    +x[IX(i  , j  , m+1)]
////                                    +x[IX(i  , j  , m-1)]
////                           )) * cRecip;
////                }
////            }
////        }
////        set_bnd(b, x, N);
////    }
////}
////
////static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N)
////{
////    float a = dt * diff * (N - 2) * (N - 2);
////    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
////}
////
////static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt, int N)
////{
////    float i0, i1, j0, j1, k0, k1;
////    
////    float dtx = dt * (N - 2);
////    float dty = dt * (N - 2);
////    float dtz = dt * (N - 2);
////    
////    float s0, s1, t0, t1, u0, u1;
////    float tmp1, tmp2, tmp3, x, y, z;
////    
////    float Nfloat = N;
////    float ifloat, jfloat, kfloat;
////    int i, j, k;
////    
////    for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
////        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
////            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
////                tmp1 = dtx * velocX[IX(i, j, k)];
////                tmp2 = dty * velocY[IX(i, j, k)];
////                tmp3 = dtz * velocZ[IX(i, j, k)];
////                x    = ifloat - tmp1; 
////                y    = jfloat - tmp2;
////                z    = kfloat - tmp3;
////                
////                if(x < 0.5f) x = 0.5f; 
////                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
////                i0 = floorf(x); 
////                i1 = i0 + 1.0f;
////                if(y < 0.5f) y = 0.5f; 
////                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
////                j0 = floorf(y);
////                j1 = j0 + 1.0f; 
////                if(z < 0.5f) z = 0.5f;
////                if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
////                k0 = floorf(z);
////                k1 = k0 + 1.0f;
////                
////                s1 = x - i0; 
////                s0 = 1.0f - s1; 
////                t1 = y - j0; 
////                t0 = 1.0f - t1;
////                u1 = z - k0;
////                u0 = 1.0f - u1;
////                
////                int i0i = i0;
////                int i1i = i1;
////                int j0i = j0;
////                int j1i = j1;
////                int k0i = k0;
////                int k1i = k1;
////                
////                d[IX(i, j, k)] = 
////                
////                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i)]
////                                +u1 * d0[IX(i0i, j0i, k1i)])
////                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i)]
////                                +u1 * d0[IX(i0i, j1i, k1i)])))
////                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i)]
////                                +u1 * d0[IX(i1i, j0i, k1i)])
////                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i)]
////                                +u1 * d0[IX(i1i, j1i, k1i)])));
////            }
////        }
////    }
////    set_bnd(b, d, N);
////}
////
////static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N)
////{
////    for (int k = 1; k < N - 1; k++) {
////        for (int j = 1; j < N - 1; j++) {
////            for (int i = 1; i < N - 1; i++) {
////                div[IX(i, j, k)] = -0.5f*(
////                         velocX[IX(i+1, j  , k  )]
////                        -velocX[IX(i-1, j  , k  )]
////                        +velocY[IX(i  , j+1, k  )]
////                        -velocY[IX(i  , j-1, k  )]
////                        +velocZ[IX(i  , j  , k+1)]
////                        -velocZ[IX(i  , j  , k-1)]
////                    )/N;		//A change based on the comments, from /N to *N.
////                p[IX(i, j, k)] = 0;
////            }
////        }
////    }
////    set_bnd(0, div, N); 
////    set_bnd(0, p, N);
////    lin_solve(0, p, div, 1, 6, iter, N);
////    
////    for (int k = 1; k < N - 1; k++) {
////        for (int j = 1; j < N - 1; j++) {
////            for (int i = 1; i < N - 1; i++) {
////                velocX[IX(i, j, k)] -= 0.5f * (  p[IX(i+1, j, k)]
////                                                -p[IX(i-1, j, k)]) * N;
////                velocY[IX(i, j, k)] -= 0.5f * (  p[IX(i, j+1, k)]
////                                                -p[IX(i, j-1, k)]) * N;
////                velocZ[IX(i, j, k)] -= 0.5f * (  p[IX(i, j, k+1)]
////                                                -p[IX(i, j, k-1)]) * N;
////            }
////        }
////    }
////    set_bnd(1, velocX, N);
////    set_bnd(2, velocY, N);
////    set_bnd(3, velocZ, N);
////}
////
////void FluidCubeStep(FluidCube *cube)
////{
////	int N          = cube->size;
////
////	free(cube->Vx0);
////    free(cube->Vy0);
////    free(cube->Vz0);
////
////	cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
////    cube->Vz0 = (float*)calloc(N * N * N, sizeof(float));
////
////    float visc     = cube->visc;
////    float diff     = cube->diff;
////    float dt       = cube->dt;
////    float *Vx      = cube->Vx;
////    float *Vy      = cube->Vy;
////    float *Vz      = cube->Vz;
////    float *Vx0     = cube->Vx0;
////    float *Vy0     = cube->Vy0;
////    float *Vz0     = cube->Vz0;
////    float *s       = cube->s;
////    float *density = cube->density;
////    
////    diffuse(1, Vx0, Vx, visc, dt, 4, N);	//Diffuse in the x direction.
////    diffuse(2, Vy0, Vy, visc, dt, 4, N);	//Diffuse in the y direction.
////    diffuse(3, Vz0, Vz, visc, dt, 4, N);	//Diffuse in the z direction.
////    
////    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
////    
////    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
////    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
////    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
////    
////    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
////    
////	diffuse(0, s, density, diff, dt, 4, N);
////    advect(0, density, s, Vx, Vy, Vz, dt, N);
////}
////
//////Adds "amount" to that cell in the cube, making that cell more dense.
////void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount)
////{
////    int N = cube->size;
////    cube->density[IX(x, y, z)] += amount;
////}
////
//////Gets the density at that cell.
////float FluidCubeGetDensity(FluidCube *cube, int x, int y, int z)
////{
////	int N = cube->size;
////	return cube->density[IX(x, y, z)];
////}
////
//////Adds some velocity, in all three directions, in the given cell in the cube.
////void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ)
////{
////    int N = cube->size;
////    int index = IX(x, y, z);
////    
////    cube->Vx[index] += amountX;
////    cube->Vy[index] += amountY;
////    cube->Vz[index] += amountZ;
////}
////#pragma endregion