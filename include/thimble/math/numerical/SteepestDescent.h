/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
 *
 *  Copyright 2014 Benjamin Tams
 *
 *  THIMBLE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  THIMBLE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with THIMBLE. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file SteepestDescent.h
 *
 * @brief
 *            Provides a mechanism for minimizing a real-valued
 *            <i>n</i>-dimensional differentiable functional via a steepest
 *            descent method
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_STEEPESTDESCENT_H_
#define THIMBLE_STEEPESTDESCENT_H_

#include <climits>

#include <thimble/dllcompat.h>
#include <thimble/math/numerical/RealFunctional.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Enumeration of codes indicating whether a minimization
	 *            attempt was successful or which kind of error occurred
	 *            during the attempt.
	 */
	typedef enum {

		/**
		 * @brief
		 *            Indicates that a minimization attempt successfully
		 *            found a reasonable approximation to a (local)
		 *            minimum.
		 */
		SD_CONVERGED = 0 ,

		/**
		 *  @brief
		 *            Indicates that the gradient is too small to further
		 *            minimize which, however, may still indicate that
		 *            the result is close to a (local) minimum.
		 */
		SD_GRAD_TOO_SMALL = 1 ,

		/**
		 * @brief
		 *            Indicates the number of maximal minimization steps
		 *            has been exceeded without the result being reasonable
		 *            close to a (local) minimum.
		 */
		SD_DIVERGED = 2 ,

		/**
		 * @brief
		 *            Indicates that an exception occurred during a
		 *            minimization attempt.
		 *
		 * @details
		 *            An exception might be thrown by calling
		 *            the \link RealFunctional::eval()\endlink
		 *            or \link RealFunctional::derive()\endlink
		 *            of which implementation is due to the programmer
		 *            of subclasses of the
		 *            abstract \link RealFunctional\endlink class.
		 */
		SD_EXCEPTION = 3

	} SD_STATE_T;

	/**
	 * @brief
	 *            Objects of this class represent minimizers for
	 *            real-valued <i>n</i>-dimensional differentiable functionals.
	 *
	 * @see RealFunctional
	 */
	class THIMBLE_DLL SteepestDescent {

	public:

		/**
		 * @brief
		 *            Creates a minimizer for real-valued <i>n</i>-dimensional
		 *            differentiable functionals.
		 *
		 * @param maxIts
		 *            The number of maximal iteration steps before a
		 *            minimization attempt via the \link minimize()\endlink
		 *            function will be interrupted and cause it to
		 *            return \link SD_DIVERGED\endlink.
		 *
		 * @param tol
		 *            The tolerance controlling stopping criteria for
		 *            this minimization object.
		 *
		 * @param sigma
		 *            A parameter for choosing the step size within each
		 *            iteration via the \link stepSize()\endlink function.
		 *
		 * @see RealFunctional
		 *
		 * @warning
		 *            If <code>maxIts</code> is smaller than or equals 0 or
		 *            if <code>tol</code> is smaller than or equals 0.0, then
		 *            an error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		SteepestDescent
			( int maxIts = INT_MAX-1 ,
			  double tol = 1e-6 ,
			  double sigma = 0.25 );

		/**
		 * @brief
		 *            Virtual destructor.
		 *
		 * @details
		 *            A virtual destructor is provided since a SteepestDescent
		 *            object has virtual protected functions.
		 *
		 * @see finished()
		 * @see stepSize()
		 */
		virtual ~SteepestDescent();

		/**
		 * @brief
		 *            Attempts to minimizes a real-valued <i>n</i>-dimensional
		 *            differentiable function via this minimizer.
		 *
		 * @details
		 *            Assume that this minimizer has been created via the
		 *            constructor call \link SteepestDescent(int,double,double)
		 *            SteepestDescent(maxIts,tol,sigma)\endlink and that we
		 *            want to minimize the <i>n</i>-dimensional real functional
		 *            \f$f\f$ (represented by the argument <code>f</code>) for
		 *            which the function \link RealFunctional::eval()
		 *            f.eval()\endlink and the
		 *            method \link RealFunctional::derive() f.derive()\endlink
		 *            has been reasonably implemented,
		 *            i.e., \link RealFunctional::derive() f.derive()\endlink
		 *            implements the gradient
		 *            function \f$\nabla~f\f$ for \f$f\f$ in a proper manner.
		 *
		 *            Given an initial vector
		 *            \f[
		 *             x = ( x[0], ... ,x[n-1] )
		 *            \f]
		 *            this function attempts to generate a sequence
		 *            \f[
		 *             x_i = ( x_i[0] , ... , x[n-1] )
		 *            \f]
		 *            where \f$x_0=x\f$ such that
		 *            \f[
		 *             f(x_i)\geq f(x_{i+1})
		 *            \f]
		 *            until the series fulfills certain convergence properties
		 *            or a stopping criteria causes to interrupt the
		 *            minimization attempt. The latest \f$x_i\f$ will be
		 *            stored in <code>xmin</code> independently from the
		 *            return value:
		 *            <table border="0">
		 *             <tr>
		 *              <td>\link SD_CONVERGED\endlink:</td>
		 *              <td>
		 *               if \f$f(x_i)-f(x_{i+1})\leq tol\cdot(1+|f(x_i)|)\f$
		 *               indicating that a reasonable approximation to a
		 *               local minimum of \f$f\f$ is stored in
		 *               <code>xmin</code>.
		 *              </td>
		 *             </tr>
		 *             <tr>
		 *              <td>\link SD_DIVERGED\endlink:</td>
		 *              <td>
		 *               if \f$i>maxIts\f$
		 *              </td>
		 *             </tr>
		 *             <tr>
		 *              <td>\link SD_GRAD_TOO_SMALL\endlink: </td>
		 *              <td>
		 *               If the gradient of \f$f\f$ at \f$x_i\f$ is smaller
		 *               than or equals \f$tol\f$, i.e.,
		 *               if \f$\nabla f(x_i)\leq tol\f$.
		 *              </td>
		 *             </tr>
		 *             <tr>
		 *              <td>\link SD_EXCEPTION\endlink: </td>
		 *              <td>
		 *               If an exception is thrown during the minimization
		 *               attempt, e.g., caused by an internal call
		 *               of \link RealFunctional::eval()
		 *               f.eval()\endlink
		 *               or \link RealFunctional::derive()
		 *               f.derive()\endlink
		 *              </td>
		 *             </tr>
		 *            </table>
		 *            If we want to modify the stopping criteria
		 *            \f$f(x_i)-f(x_{i+1})\leq tol\cdot(1+|f(x_i)|)\f$
		 *            returning \link SD_CONVERGED\endlink, we may
		 *            override the virtual \link finished()\endlink function.
		 *
		 * @param xmin
		 *            An array capable to hold
		 *            <i>n</i>=\link RealFunctional::getPreimageDimension() f.getPreimageDimension()\endlink
		 *            double values in which the approximation of the minimum
		 *            will be stored.
		 *
		 * @param f
		 *            A real function providing a proper implementation of
		 *            its virtual \link RealFunctional::eval() eval()\endlink
		 *            function
		 *            and \link RealFunction::derive() derive()\endlink
		 *            method.
		 *
		 * @param x
		 *            An array holding
		 *            <i>n</i>=\link RealFunctional::getPreimageDimension() f.getPreimageDimension()\endlink
		 *            double values representing the initial vector from which
		 *            as series of descending vectors is generated.
		 *
		 * @return
		 *            A state indicating how to interpret the result stored
		 *            in <code>xmin</code>
		 */
		SD_STATE_T minimize
			( double *xmin , const RealFunctional & f , const double *x ) const;

	protected:

		/**
		 * @brief
		 *            Tests whether the solution <code>x1</code> is
		 *            reasonably close to a local minimum of the
		 *            real functional <code>f</code>.
		 *
		 * @details
		 *            This function is called by the \link minimize()\endlink
		 *            function to test whether an element \f$x_{i+1}\f$ of the
		 *            descending sequence \f$x_0,...,x_i,x_{i+1}\f$ where
		 *            \f$f(x_0)\geq ...\geq f(x_i)\geq f(x_{i+1})\f$ should be
		 *            considered a reasonable approximation to a local minimum
		 *            and return \link SD_CONVERGED\endlink correspondingly.
		 *
		 *            In this implementation, the function tests
		 *            whether \f$f(x_0)-f(x_1)\leq tol\cdot (1+|f(x_1)|)\f$
		 *            and returns <code>true</code> if the inequality
		 *            is fulfilled and <code>false</code> otherwise;
		 *            thereby, the function only accounts the parameters
		 *            <code>f</code>, <code>x0</code> and <code>x1</code>.
		 *
		 *            This function can be overriden to modify the stopping
		 *            criteria in the \link minimize()\endlink function
		 *            causing it to return \link SD_CONVERGED\endlink
		 *            which indicates a successful minimization. Thereby
		 *            we can also account for the gradient \f$\nabla f(x_1)\f$
		 *            and its norm \f$\|\nabla f(x_1)\|_2\f$ represented
		 *            by <code>grad1</code> and <code>gradNorm1</code>,
		 *            respectively.
		 *
		 * @param f
		 *            The real functional to be minimized.
		 *
		 * @param x0
		 *            One element in the sequence of a minimization
		 *            attempt.
		 *
		 * @param x1
		 *            The element after <code>x0</code> in the sequence
		 *            of a minimization attempt such
		 *            that \f$f(x_1)\leq f(x_0)\f$ should be fulfilled.
		 *
		 * @param grad1
		 *            The gradient of \f$f\f$ at \f$x_1\f$, i.e.,
		 *            \f$\nabla f(x_1)\f$ as it is the result
		 *            of \link RealFunctional::derive()
		 *            f.derive(grad1,x1)\endlink.
		 *
		 * @param gradNorm1
		 *            The norm \f$\|\nabla f(x_1)\|_2\f$.
		 *
		 * @return
		 *            <code>true</code> if the stopping criteria is
		 *            fulfilled and, otherwise, <code>false</code>;
		 *            also see the details for this function.
		 */
		virtual bool finished
			( const RealFunctional & f ,
			  double *x0 ,
			  double *x1 ,
			  double *grad1 ,
			  double gradNorm1 ) const;

		/**
		 * @brief
		 *            Chooses the step size controlling how far the
		 *            next element in the descending minimization sequence
		 *            on base of its previous element is to be chosen along
		 *            the direction of the to-be-minimized functionals
		 *            steepest descent.
		 *
		 * @param f
		 *            The real functional that is attempted to be minimized.
		 *
		 * @param x
		 *            An element in the minimization sequence from which
		 *            the next element is going to be selected.
		 *
		 * @param negGrad
		 *            The negative gradient of \f$f\f$ at \f$x\f$, i.e.,
		 *            \f$-\nabla f(x)\f$.
		 *
		 * @param tmp1
		 *            An array that is capable to hold <i>n</i> double values
		 *            and that is used as a temporary array by the
		 *            implementation of this function. The array is not
		 *            allocated within this function for internal performance
		 *            reasons but may be ignored if we override this function
		 *            to modify the step size selection.
		 *
		 * @param tmp2
		 *            see <code>tmp1</code>.
		 *
		 * @return
		 *            The step size for the next element in the minimization
		 *            sequence on base of \f$x\f$ and its negative gradient
		 *            \f$-\nabla f(x)\f$.
		 */
		virtual double stepSize
			( const RealFunctional & f ,
			  double *x ,
			  double *negGrad ,
			  double *tmp1 ,
			  double *tmp2 ) const;

	private:

		/**
		 * @brief
		 *            The number of maximal iteration steps before a
		 *            minimization attempt via the \link minimize()\endlink
		 *            function will be interrupted and cause it to
		 *            return \link SD_DIVERGED\endlink.
		 *
		 * @see \link SteepestDescent() SteepestDescent(double,int,double)\endlink
		 */
		int maxIts;

		/**
		 * @brief
		 *            The tolerance controlling stopping criteria for
		 *            this minimization object.
		 *
		 * @see \link SteepestDescent() SteepestDescent(double,int,double)\endlink
		 */
		double tol;

		/**
		 * @brief
		 *            A parameter for choosing the step size within each
		 *            iteration via the \link stepSize()\endlink function.
		 *
		 * @see stepSize()
		 * @see \link SteepestDescent() SteepestDescent(double,int,double)\endlink
		 */
		double sigma;
	};
}

#endif /* THIMBLE_STEEPESTDESCENT_H_ */
