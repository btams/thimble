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
 * @file SteepestDescent.cpp
 *
 * @brief
 *            Implementation of a steepest descent method for minimizing
 *            a real-valued <i>n</i>-dimensional differentiable functional
 *            as provided by the 'SteepestDescent.h' header file.
 *
 * @author Benjamin Tams
 */

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/math/numerical/RealFunctional.h>
#include <thimble/math/numerical/SteepestDescent.h>

using namespace std;

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	SteepestDescent::SteepestDescent( int maxIts , double tol , double sigma ) {

		if ( maxIts <= 0 ) {
			cerr << "SteepestDescent: number of iterations must be greater "
				 << "than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( tol <= 0.0 ) {
			cerr << "SteepestDescent: tolerance must be greater than 0.0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->maxIts = maxIts;
		this->tol    = tol;
		this->sigma  = sigma;
	}

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
	SteepestDescent::~SteepestDescent() {
	}

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
	SD_STATE_T SteepestDescent::minimize
		( double *xmin, const RealFunctional & f , const double *x ) const {

		SD_STATE_T state = SD_DIVERGED;
		double *x1 = xmin;

		int n = f.getPreimageDimension();

		// Allocate memory
		double *x0 , *grad , *tmp1 , *tmp2;
		x0   = (double*)malloc( n * sizeof(double) );
		grad = (double*)malloc( n * sizeof(double) );
		tmp1 = (double*)malloc( n * sizeof(double) );
		tmp2 = (double*)malloc( n * sizeof(double) );
		if ( x0 == NULL || grad == NULL || tmp1 == NULL || tmp2 == NULL ) {
			cerr << "GradientDescent::minimize: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		memcpy(x0,x,n*sizeof(double));
		if ( x1 != x ) {
			memcpy(x1,x,n*sizeof(double));
		}

		try {
			// Compute steepest descent
			f.derive(grad,x0);
			for ( int i = 0 ; i < n ; i++ ) {
				grad[i] = -grad[i];
			}

			for ( int it = 0 ; it < this->maxIts ; it++ ) {

				// Step in the direction of the steepest descent to
				// furthermore make the function smaller.
				double t = stepSize(f,x0,grad,tmp1,tmp2);
				for ( int i = 0 ; i < n ; i++ ) {
					x1[i] = x0[i]+t*grad[i];
				}

				// Compute next steepest descent
				f.derive(grad,x1);
				for ( int i = 0 ; i < n ; i++ ) {
					grad[i] = -grad[i];
				}

				// Compute the norm of the gradient ...
				double gradNorm = 0.0;
				for ( int i = 0 ; i < n ; i++ ) {
					gradNorm += grad[i]*grad[i];
				}
				gradNorm = sqrt(gradNorm);

				// ... and check whether it is too small.
				if ( gradNorm <= DBL_EPSILON ) {
					memcpy(x1,x0,n*sizeof(double));
					state = SD_GRAD_TOO_SMALL;
					break;
				}

				if ( finished(f,x0,x1,grad,gradNorm) ) {
					memcpy(x1,x0,n*sizeof(double));
					state = SD_CONVERGED;
					break;
				}

				memcpy(x0,x1,n*sizeof(double));

			}
		} catch(...) {
			state = SD_EXCEPTION;
		}

		free(x0);
		free(grad);
		free(tmp1);
		free(tmp2);

		return state;
	}

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
	bool SteepestDescent::finished
		( const RealFunctional & f ,
		  double *x0 ,
		  double *x1 ,
		  double *grad1 ,
		  double gradNorm1 ) const {

		double y0 , y1;

		y0 = f.eval(x0);
		y1 = f.eval(x1);

		return y0-y1 <= tol*(1.0+abs(y1));
	}

	/**
	 * @brief
	 *           Helper class for choosing step sizes.
	 *
	 * @details
	 *           Given the negative gradient \f$d=-\nabla f(x)\f$ of a real
	 *           functional \f$f\f$ at \f$x\f$
	 *           instances of this class represent the real-valued
	 *           one-dimensional function
	 *           \f[
	 *            \Phi(t)=f(x+t\cdot d).
	 *           \f]
	 *
	 * @see SteepestDescent::stepSize()
	 */
	class Phi {

	public:

		/**
		 * @brief
		 *           Pointer to the to-be-minimized real functional
		 *           of dimension <i>n</i>.
		 */
		const RealFunctional *fPtr;

		/**
		 * @brief
		 *           <i>n</i>-dimensional vector from which the next
		 *           element furthermore decreasing the value of \f$f\f$
		 *           is to be selected.
		 */
		double *x;

		/**
		 * @brief
		 *           Negative gradient of the to-be-minimized real
		 *           functional at \f$x\f$, i.e., \f$d=-\nabla f(x)\f$.
		 */
		double *negGrad;

		/**
		 * @brief
		 *           The same parameter <code>%sigma</code> as the
		 *           private member of a \link SteepestDescent\endlink
		 *           object.
		 */
		double sigma;

		/**
		 * @brief
		 *           The value \f$\Phi(0)=f(x)\f$.
		 */
		double phi0;

		/**
		 * @brief
		 *           The derivation \f$\frac{\partial\Phi(0)}{\partial t}\f$.
		 */
		double dphi0;

		/**
		 * @brief
		 *           Temporary array capabale of holding <i>n</i>
		 *           double values.
		 */
		double *tmp1;

		/**
		 * @brief
		 *           Temporary array capabale of holding <i>n</i>
		 *           double values.
		 */
		double *tmp2;

		/**
		 * @brief
		 *           Creates a helper object for choosing step sizes
		 *           via the \link SteepestDescent::stepSize()\endlink
		 *           function.
		 *
		 * @param f
		 *           The function that is to be minimized.
		 *
		 * @param x
		 *           An array holding <i>n</i> double values.
		 *
		 * @param negGrad
		 *           The negative gradient of \f$f\f$ at \f$x\f$, i.e.,
		 *           \f$-\nabla f(x)\f$.
		 *
		 * @param tmp1
		 *            An array that is capable to hold <i>n</i> double values
		 *            and that is used as a temporary array by the
		 *            implementation of this function. The array is not
		 *            allocated within this function for internal performance
		 *            reasons but may be ignored if we override this function
		 *            to modify the step size selection.
		 *
		 * @param sigma
		 *            The same parameter <code>%sigma</code> as the
		 *            private member of the \link SteepestDescent\endlink
		 *            object that is attempting to minimize \f$f\f$.
		 *
		 * @param tmp2
		 *            see <code>tmp1</code>.
		 */
		Phi( const RealFunctional & f ,
			 double *x , double *negGrad , double sigma ,
			 double *tmp1 , double *tmp2 );


		/**
		 * @brief
		 *            Evaluates this helper function \f$\Phi(\cdot)\f$
		 *            at \f$t\f$.
		 *
		 * @param t
		 *            The real position at which this helper function
		 *            is evaluated.
		 *
		 * @return
		 *            The evaluation of this helper function \f$\Phi(\cdot)\f$
		 *            at \f$t\f$, i.e., \f$\Phi(t)\f$.
		 */
		double eval( double t );

		/**
		 * @brief
		 *            Determine the deriviation of this helper function
		 *            \f$\Phi(\cdot)\f$ at \f$t\f$, i.e.,
		 *            \f$\frac{\partial\Phi(t)}{\partial t}\f$.
		 *
		 * @param t
		 *            The position at which the derivation of this helper
		 *            function is evaluation.
		 *
		 * @return
		 *            \f$\frac{\partial\Phi(t)}{\partial t}\f$
		 */
		double deriv( double t );

		/**
		 * @brief
		 *            Returns \f$\Phi(0)+\sigma\cdot t\cdot \Phi'(0)\f$
		 *            where \f$\Phi'\f$ denotes the derivation of \f$\Phi\f$.
		 *
		 * @param t
		 *            see above.
		 *
		 * @return
		 *            see above.
		 */
		inline double eval_o( double t ) {
			return phi0+sigma*t*dphi0;
		}

		/**
		 * @brief
		 *            Returns \f$\Phi(0)+(1-\sigma)\cdot t\cdot \Phi'(0)\f$
		 *            where \f$\Phi'\f$ denotes the derivation of \f$\Phi\f$.
		 *
		 * @param t
		 *            see above.
		 *
		 * @return
		 *            see above.
		 */
		inline double eval_u( double t ) {
			return phi0+(1.0-sigma)*t*dphi0;
		}
	};

	/**
	 * @brief
	 *           Creates a helper object for choosing step sizes
	 *           via the \link SteepestDescent::stepSize()\endlink
	 *           function.
	 *
	 * @param f
	 *           The function that is to be minimized.
	 *
	 * @param x
	 *           An array holding <i>n</i> double values.
	 *
	 * @param negGrad
	 *           The negative gradient of \f$f\f$ at \f$x\f$, i.e.,
	 *           \f$-\nabla f(x)\f$.
	 *
	 * @param tmp1
	 *            An array that is capable to hold <i>n</i> double values
	 *            and that is used as a temporary array by the
	 *            implementation of this function. The array is not
	 *            allocated within this function for internal performance
	 *            reasons but may be ignored if we override this function
	 *            to modify the step size selection.
	 *
	 * @param sigma
	 *            The same parameter <code>%sigma</code> as the
	 *            private member of the \link SteepestDescent\endlink
	 *            object that is attempting to minimize \f$f\f$.
	 *
	 * @param tmp2
	 *            see <code>tmp1</code>.
	 */
	Phi::Phi
	( const RealFunctional & f ,
	  double *x , double *negGrad , double sigma ,
	  double *tmp1 , double *tmp2 ) {

		this->tmp1 = tmp1;
		this->tmp2 = tmp2;

		this->fPtr = &f;
		this->x = x;
		this->negGrad = negGrad;

		this->sigma = sigma;

		this->phi0 = eval(0);
		this->dphi0 = deriv(0);

	}

	/**
	 * @brief
	 *            Evaluates this helper function \f$\Phi(\cdot)\f$
	 *            at \f$t\f$.
	 *
	 * @param t
	 *            The real position at which this helper function
	 *            is evaluated.
	 *
	 * @return
	 *            The evaluation of this helper function \f$\Phi(\cdot)\f$
	 *            at \f$t\f$, i.e., \f$\Phi(t)\f$.
	 */
	double Phi::eval( double t ) {

		int n = fPtr->getPreimageDimension();

		for ( int i = 0 ; i < n; i++ ) {
			this->tmp1[i] = this->x[i]+t*this->negGrad[i];
		}

		return this->fPtr->eval(this->tmp1);
	}

	/**
	 * @brief
	 *            Determine the deriviation of this helper function
	 *            \f$\Phi(\cdot)\f$ at \f$t\f$, i.e.,
	 *            \f$\frac{\partial\Phi(t)}{\partial t}\f$.
	 *
	 * @param t
	 *            The position at which the derivation of this helper
	 *            function is evaluation.
	 *
	 * @return
	 *            \f$\frac{\partial\Phi(t)}{\partial t}\f$
	 */
	double Phi::deriv( double t ) {

		// As Phi(t)=f(x+t*d) we can apply the chain rule
		// for deriving real functionals.

		int n = fPtr->getPreimageDimension();

		for ( int i = 0 ; i < n ; i++ ) {
			this->tmp1[i] = this->x[i]+t*this->negGrad[i];
		}

		fPtr->derive(this->tmp2,this->tmp1);

		double diff = 0.0;
		for ( int i = 0 ; i < n ; i++ ) {
			diff += this->tmp2[i]*this->negGrad[i];
		}

		return diff;
	}

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
	double SteepestDescent::stepSize
		( const RealFunctional & f ,
	      double *x ,
	      double *negGrad ,
	      double *tmp1 ,
	      double *tmp2 ) const {

		double tu = 0.0 , to = 1.0;

		Phi phi(f,x,negGrad,this->sigma,tmp1,tmp2);

		while ( phi.eval(to) < phi.eval_u(to) ) {
			tu = to;
			to += to;
		}

		if ( phi.eval(to) < phi.eval_o(to) ) {
			return to;
		}

		double phit , t , old_t = -1.0;

		do {
			t = 0.5*(tu+to);
			if ( t == old_t ) {
				break;
			}
			old_t = t;
			phit = phi.eval(t);
			if ( phit < phi.eval_u(t) ) {
				tu = t;
			}
			if ( phit > phi.eval_o(t) ) {
				to = t;
			}
		} while ( !(phi.eval_u(t) <= phit && phit <= phi.eval_o(t)) );

		return t;
	}
}




