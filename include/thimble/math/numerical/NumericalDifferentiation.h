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
 * @file NumericalDifferentiation.h
 *
 * @brief
 *            Provides approximation functions for evaluating a real
 *            functional's gradients.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_NUMERICALDIFFERENTIATION_H_
#define THIMBLE_NUMERICALDIFFERENTIATION_H_

#include <thimble/dllcompat.h>
#include <thimble/math/numerical/RealFunctional.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Class with static functions estimating point-wise gradients
	 *            of real functionals.
	 */
	class THIMBLE_DLL NumericalDifferentiation {

	private:

		/**
		 * @brief
		 *            Private standard constructor.
		 *
		 * @details
		 *            The constructor is provided as private to prevent users
		 *            from creating objects of this class.
		 */
		inline NumericalDifferentiation() {
		}

	public:

		/**
		 * @brief
		 *            Given a
		 *            real functional \f$f:R^n\rightarrow R\f$, this method
		 *            approximates its gradient.
		 *
		 * @details
		 *            This method computes for each \f$i=0,...,n-1\f$
		 *            \f[
		 *             grad[i]=\frac{f(x[0],...,x[i-1]~,~x[i]-h~,~x[i+1],...,x[n-1])-f(x[0],...,x[i-1]~,~x[i]+h~,~x[i+1],...,x[n-1])}{2\cdot h}
		 *            \f]
		 *            such that, for \f$0<h\rightarrow 0\f$, the
		 *            vector \f$(grad[0],...,grad[n-1])\f$ approaches the
		 *            gradient of \f$f\f$ at \f$(x[0],...,x[n-1])\f$
		 *            if \f$f\f$ at \f$(x[0],...,x[n-1])\f$ is differentiable.
		 *
		 *            This method changes the content of the array <code>x</code>
		 *            but on output the content of <code>x</code> will be
		 *            set to the content of <code>x</code> as it has been
		 *            input to this method.
		 *
		 * @attention
		 *            This method calls
		 *            the \link RealFunctional::eval() f.eval()\endlink method
		 *            and thus the behavior of this method depends on the
		 *            implementation of
		 *            the \link RealFunctional::eval() f.eval()\endlink as
		 *            well.
		 *
		 * @param grad
		 *            An array capabale of storing \f$n\f$ double values in
		 *            which the gradient of \f$f\f$ at \f$(x[0],...,x[n-1])\f$
		 *            will be stored.
		 *
		 * @param f
		 *            A reference to an object representing a real functional.
		 *
		 * @param x
		 *            An array containing \f$n\f$ double values being a valid
		 *            argument
		 *            to \link RealFunctional::eval() f.eval()\endlink.
		 *
		 * @param h
		 *            A double value greater than 0.0.
		 *
		 * @warning
		 *            If <code>h</code> is smaller than or equals 0.0, then
		 *            an error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		static void jacobi
			( double *grad ,
			  const RealFunctional & f ,
			  double *x ,
			  double h );

		/**
		 * @brief
		 *            Given a
		 *            real functional \f$f:R^n\rightarrow R\f$, this method
		 *            approximates its gradient.
		 *
		 * @details
		 *            This method computes for each \f$i=0,...,n-1\f$
		 *            \f[
		 *             grad[i]=\frac{f(x[0],...,x[i-1]~,~x[i]-h~,~x[i+1],...,x[n-1])-f(x[0],...,x[i-1]~,~x[i]+h~,~x[i+1],...,x[n-1])}{2\cdot h}
		 *            \f]
		 *            such that, for \f$0<h\rightarrow 0\f$, the
		 *            vector \f$(grad[0],...,grad[n-1])\f$ approaches the
		 *            gradient of \f$f\f$ at \f$(x[0],...,x[n-1])\f$
		 *            if \f$f\f$ at \f$(x[0],...,x[n-1])\f$ is differentiable.
		 *
		 * @attention
		 *            This method calls
		 *            the \link RealFunctional::eval() f.eval()\endlink method
		 *            and thus the behavior of this method depends on the
		 *            implementation of
		 *            the \link RealFunctional::eval() f.eval()\endlink as
		 *            well.
		 *
		 * @param grad
		 *            An array capabale of storing \f$n\f$ double values in
		 *            which the gradient of \f$f\f$ at \f$(x[0],...,x[n-1])\f$
		 *            will be stored.
		 *
		 * @param f
		 *            A reference to an object representing a real functional.
		 *
		 * @param x
		 *            An array containing \f$n\f$ double values being a valid
		 *            argument
		 *            to \link RealFunctional::eval() f.eval()\endlink.
		 *
		 * @param h
		 *            A double value greater than 0.0.
		 *
		 * @warning
		 *            If <code>h</code> is smaller than or equals 0.0, then
		 *            an error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error message
		 *            will be printed to <code>stderr</code> and then the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		static void jacobi
			( double *grad ,
			  const RealFunctional & f ,
			  const double *x ,
			  double h );
	};
}


#endif /* THIMBLE_NUMERICALDIFFERENTIATION_H_ */
