#pragma once
#include <math.h>
#include <stdlib.h>

//http://www.mat.ucsb.edu/~wakefield/594cm/assignment.htm

struct Fluid
{
	unsigned int bits, dim, dim2, dim3, dimwrap;		//Discrete field dimensions.
	unsigned int scale;	// division of region dimension to fluid dimension (pow of 2)

	double *pfield, *pfield0;		//Pressure fields.
	double *dfield, *dfield0;		//Density fields.
	double *uxfield, *uxfield0;		//Velocity fields.
	double *uyfield, *uyfield0;
	double *uzfield, *uzfield0;

	double viscosity, diffusion, decay;		//Global fluid parameters for the cube.

	Fluid(unsigned int bits, double viscosity = 0.00001, double diffusion = 0.001, double decay = 0.99)
		: bits(bits), dim (1 << bits), viscosity(viscosity), diffusion(diffusion), decay(decay)
	{
		dim2 = dim * dim;			//dim squared.
		dim3 = dim * dim * dim;		//dim cube.
		dimwrap = dim - 1;

		//dim = DIMENSION

		//Initialising all the arrays.
		dfield = new double[dim3];
		uxfield = new double[dim3];
		uyfield = new double[dim3];
		uzfield = new double[dim3];
		dfield0 = new double[dim3];
		uxfield0 = new double[dim3];
		uyfield0 = new double[dim3];
		uzfield0 = new double[dim3];

		//Initialising the density field to 0 "empty" initially,
		//and adding random velocities to the velocity fields.
		for (int i=0; i<dim3; i++) 
		{
			dfield0[i] = dfield[i] = 0;
			uxfield0[i] = uxfield[i] = ((double) rand() / (RAND_MAX));
			uyfield0[i] = uyfield[i] = ((double) rand() / (RAND_MAX));
			uzfield0[i] = uzfield[i] = ((double) rand() / (RAND_MAX));
		}
	}

	~Fluid()
	{
		delete[] dfield;
		delete[] uxfield;
		delete[] uyfield;
		delete[] uzfield;
		delete[] dfield0;
		delete[] uxfield0;
		delete[] uyfield0;
		delete[] uzfield0;
	}

	void step();

	inline void addDensity(glm::vec3 &v, double d = 1)
	{
		dfield[index(v.x, v.y, v.z)] = d;
	}

	inline void addVelocity(glm::vec3 &v, glm::vec3 &d)
	{
		uxfield[index(v.x, v.y, v.z)] += d.x;
		uyfield[index(v.x, v.y, v.z)] += d.y;
		uzfield[index(v.x, v.y, v.z)] += d.z;
	}

	inline void advet(double *f, double *f0, double *ux, double *uy, double *uz, double rate = 1)
	{
		for(int x = 0; x < dim; x++)
		{
			for(int y = 0; y < dim; y++)
			{
				for(int z = 0; z < dim; z++)
				{
					int i = index(x, y, z);

					//Get co-ordinates back in time from current cell:
					double x0 = x - rate * ux[i];
					double y0 = y - rate * uy[i];
					double z0 = z - rate * uz[i];

					//Get interpolated value at x0, y0, z0:
					f[i] = interp(f0, x0, y0, z0);
				}
			}
		}
	}

	inline void project(int N, double *vx, double *vy, double *vz, double *g, double *g0)
	{
		double h = 1.0 / N;

		// g and g0 represent the velocity gradient
		BEGIN_PER_CELL
			// previous instantaneous magnitude of velocity gradient 
			//		= (sum of velocity gradients per axis)/2N:
			g0[i] = -0.5 * h * (
									vx[INDEX(x+1, y, z)]-vx[INDEX(x-1, y, z)]+	// velocity gradients
									vy[INDEX(x, y+1, z)]-vy[INDEX(x, y-1, z)]+	// velocity gradients
									vz[INDEX(x, y, z+1)]-vz[INDEX(x, y, z-1)]	// velocity gradients
								);
			// zero out the present velocity gradient:
			g[i] = 0;
		END_PER_CELL

		// reuse the Gauss-Seidel relaxation solver to safely diffuse the velocity gradients from g0 to g:
		stable_solve(N, g, g0, 1, 6);
	
		// now subtract this gradient from our current velocity field:
		BEGIN_PER_CELL
			vx[i] -= 0.5 * N * (g[INDEX(x+1, y, z)]-g[INDEX(x-1, y, z)]); // gradient calculated by neighbors
			vy[i] -= 0.5 * N * (g[INDEX(x, y+1, z)]-g[INDEX(x, y-1, z)]);
			vz[i] -= 0.5 * N * (g[INDEX(x, y, z+1)]-g[INDEX(x, y, z-1)]);
		END_PER_CELL
	}

	inline unsigned int index(int x, int y, int z) 
	{
		return ((x/scale) & dimwrap) * dim2 
			 + ((y/scale) & dimwrap) * dim 
			 + ((z/scale) & dimwrap);
	}
	inline unsigned int index(double x, double y, double z) 
	{
		return index((int)floor(x), (int)floor(y), (int)floor(z));
	}
}


#pragma region Original
#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)

struct FluidCube
{
    int size;		//The size of the cube, e.g. size * size * size
    float dt;		//The delta time to advance the state by during each update.
    float diff;		//The diffusion amount (?) :: How fast stuff spreads out in the fluid.
    float visc;		//The viscosity of the cube (?) :: How thick the fluid is.
    
    float* s;			//Scratch space array for density. (Keeps the old / previous values while the new ones are computed.)
    float* density;		//Density array.
    
    float* Vx;			//Velocity arrays for X, Y and Z.
    float* Vy;
    float* Vz;

    float* Vx0;			//Scratch space arrays for X, Y and Z.
    float* Vy0;
    float* Vz0;
};


FluidCube* FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
{
	//Allocates memory for the cube, and creates arrays of size * size * size for each value.
    FluidCube* cube = (FluidCube*)malloc(sizeof(*cube));
    int N = size;
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = (float*)calloc(N * N * N, sizeof(float));
    cube->density = (float*)calloc(N * N * N, sizeof(float));
    
    cube->Vx = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy = (float*)calloc(N * N * N, sizeof(float));
    cube->Vz = (float*)calloc(N * N * N, sizeof(float));
    
    cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vz0 = (float*)calloc(N * N * N, sizeof(float));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    
    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vz);
    
    free(cube->Vx0);
    free(cube->Vy0);
    free(cube->Vz0);
    
    free(cube);
}


static void set_bnd(int b, float *x, int N)
{
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0  )] = b == 3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
            x[IX(i, j, N-1)] = b == 3 ? -x[IX(i, j, N-2)] : x[IX(i, j, N-2)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k)] = b == 2 ? -x[IX(i, 1  , k)] : x[IX(i, 1  , k)];
            x[IX(i, N-1, k)] = b == 2 ? -x[IX(i, N-2, k)] : x[IX(i, N-2, k)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k)] = b == 1 ? -x[IX(1  , j, k)] : x[IX(1  , j, k)];
            x[IX(N-1, j, k)] = b == 1 ? -x[IX(N-2, j, k)] : x[IX(N-2, j, k)];
        }
    }
    
    x[IX(0, 0, 0)]       = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
    x[IX(0, N-1, 0)]     = 0.33f * (x[IX(1, N-1, 0)]
                                  + x[IX(0, N-2, 0)]
                                  + x[IX(0, N-1, 1)]);
    x[IX(0, 0, N-1)]     = 0.33f * (x[IX(1, 0, N-1)]
                                  + x[IX(0, 1, N-1)]
                                  + x[IX(0, 0, N)]);
    x[IX(0, N-1, N-1)]   = 0.33f * (x[IX(1, N-1, N-1)]
                                  + x[IX(0, N-2, N-1)]
                                  + x[IX(0, N-1, N-2)]);
    x[IX(N-1, 0, 0)]     = 0.33f * (x[IX(N-2, 0, 0)]
                                  + x[IX(N-1, 1, 0)]
                                  + x[IX(N-1, 0, 1)]);
    x[IX(N-1, N-1, 0)]   = 0.33f * (x[IX(N-2, N-1, 0)]
                                  + x[IX(N-1, N-2, 0)]
                                  + x[IX(N-1, N-1, 1)]);
    x[IX(N-1, 0, N-1)]   = 0.33f * (x[IX(N-2, 0, N-1)]
                                  + x[IX(N-1, 1, N-1)]
                                  + x[IX(N-1, 0, N-2)]);
    x[IX(N-1, N-1, N-1)] = 0.33f * (x[IX(N-2, N-1, N-1)]
                                  + x[IX(N-1, N-2, N-1)]
                                  + x[IX(N-1, N-1, N-2)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m)] =
                        (x0[IX(i, j, m)]
                            + a*(    x[IX(i+1, j  , m  )]
                                    +x[IX(i-1, j  , m  )]
                                    +x[IX(i  , j+1, m  )]
                                    +x[IX(i  , j-1, m  )]
                                    +x[IX(i  , j  , m+1)]
                                    +x[IX(i  , j  , m-1)]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j, k;
    
    for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k)];
                tmp2 = dty * velocY[IX(i, j, k)];
                tmp3 = dtz * velocZ[IX(i, j, k)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                z    = kfloat - tmp3;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f; 
                if(z < 0.5f) z = 0.5f;
                if(z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                k0 = floorf(z);
                k1 = k0 + 1.0f;
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                u1 = z - k0;
                u0 = 1.0f - u1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                int k0i = k0;
                int k1i = k1;
                
                d[IX(i, j, k)] = 
                
                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                                +u1 * d0[IX(i0i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                                +u1 * d0[IX(i0i, j1i, k1i)])))
                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                                +u1 * d0[IX(i1i, j0i, k1i)])
                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                                +u1 * d0[IX(i1i, j1i, k1i)])));
            }
        }
    }
    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N)
{
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k)] = -0.5f*(
                         velocX[IX(i+1, j  , k  )]
                        -velocX[IX(i-1, j  , k  )]
                        +velocY[IX(i  , j+1, k  )]
                        -velocY[IX(i  , j-1, k  )]
                        +velocZ[IX(i  , j  , k+1)]
                        -velocZ[IX(i  , j  , k-1)]
                    )/N;		//A change based on the comments, from /N to *N.
                p[IX(i, j, k)] = 0;
            }
        }
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k)] -= 0.5f * (  p[IX(i+1, j, k)]
                                                -p[IX(i-1, j, k)]) * N;
                velocY[IX(i, j, k)] -= 0.5f * (  p[IX(i, j+1, k)]
                                                -p[IX(i, j-1, k)]) * N;
                velocZ[IX(i, j, k)] -= 0.5f * (  p[IX(i, j, k+1)]
                                                -p[IX(i, j, k-1)]) * N;
            }
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void FluidCubeStep(FluidCube *cube)
{
	int N          = cube->size;

	free(cube->Vx0);
    free(cube->Vy0);
    free(cube->Vz0);

	cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vz0 = (float*)calloc(N * N * N, sizeof(float));

    float visc     = cube->visc;
    float diff     = cube->diff;
    float dt       = cube->dt;
    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);	//Diffuse in the x direction.
    diffuse(2, Vy0, Vy, visc, dt, 4, N);	//Diffuse in the y direction.
    diffuse(3, Vz0, Vz, visc, dt, 4, N);	//Diffuse in the z direction.
    
    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    
    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
    
	diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, Vz, dt, N);
}

//Adds "amount" to that cell in the cube, making that cell more dense.
void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount)
{
    int N = cube->size;
    cube->density[IX(x, y, z)] += amount;
}

//Gets the density at that cell.
float FluidCubeGetDensity(FluidCube *cube, int x, int y, int z)
{
	int N = cube->size;
	return cube->density[IX(x, y, z)];
}

//Adds some velocity, in all three directions, in the given cell in the cube.
void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{
    int N = cube->size;
    int index = IX(x, y, z);
    
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
    cube->Vz[index] += amountZ;
}
#pragma endregion