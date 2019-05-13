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
 * @file NumericalDifferentiation.cpp
 *
 * @brief
 *            Implements approximation functions for evaluating a real
 *            functional's gradients.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/math/numerical/RealFunctional.h>
#include <thimble/math/numerical/NumericalDifferentiation.h>

using namespace std;

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

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
	void NumericalDifferentiation::jacobi
		( double *grad , const RealFunctional & f , double *x , double h ) {

		if ( h <= 0.0 ) {
			cerr << "NumericalDifferentiation::jacobi: Bad arguments; "
				 << "'h' must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		int n = f.getPreimageDimension();

		for ( int i = 0 ; i < n ; i++ ) {

			double xi = x[i];

			double y0;
			x[i] = xi - h;
			y0 = f.eval(x);

			double y1;
			x[i] = xi + h;
			y1 = f.eval(x);

			x[i] = xi;

			grad[i] = (y1-y0) / (h+h);
		}
	}

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
	void NumericalDifferentiation::jacobi
		( double *grad ,
		  const RealFunctional & f ,
		  const double *x ,
		  double h ) {

		int n = f.getPreimageDimension();

		double *tmp = (double*)malloc( n * sizeof(double) );
		if ( tmp == NULL ) {
			cerr << "NumericalDifferentiation::jacobi: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}
		memcpy(tmp,x,n*sizeof(double));

		jacobi(grad,f,tmp,h);

		free(tmp);
	}

}



