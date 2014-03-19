#ifndef __CONFIG_H__
#define __CONFIG_H__

// ----- linear equation solver

/** Global constant: epsilon.
 * This constant is used to determine convergence in the linear equation
 * solver. Lower values result in increased accuracy whereas higher values
 * reduce the calculation time.
 */
const float gc_epsilon = 0.00001f;

/** Global constant: maximum number of iterations.
 * This constant limits the number of iterations in the linear equation
 * solver. The solver usually terminates when the error falls below
 * gc_epsilon, but in some cases (mainly while debugging) the solver
 * converges very slow or not at all. In this case, gc_max_iteration
 * limits the number of steps, keeping the simulation running although
 * the results may be incorrect.
 */
const int gc_max_iteration = 50;

// ----- algorithms

/** Use conjugate gradient to solve linear equations.
 * Conjugate gradient is a fast algorithm to solve linear equations. If you
 * deactivate it, the slower Gauss-Seidel relaxation is used. This should
 * be required for debug purposes only.
 */
#define CONJUGATE_GRADIENT

/** Numerical method used to trace particles.
 * The particle tracer approximates a differential equation to determine the
 * origin of each particle according to the velocity field. You can choose
 * between three methods:
 * The explicit Euler first order is quite inaccurate but very fast. The 
 * explicit Runge-Kutta second order gives good looking results but takes
 * a bit more cpu power. The explicit Runge-Kutta fourth order (also known
 * as classic Runge-Kutta) provides the best results, but takes quite a lot
 * of cpu power, it's not recommended in realtime-environments.
 */
//#define EULER_FIRST_ORDER
#define RUNGEKUTTA_SECOND_ORDER
//#define RUNGEKUTTA_FOURTH_ORDER

// ----- features

#endif // __CONFIG_H__