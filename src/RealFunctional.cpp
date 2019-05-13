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
 * @file RealFunctional.cpp
 *
 * @brief
 *            Implementation file of an abstract purely virtual class
 *            representing real-valued <i>n</i>-dimensional functions, i.e.,
 *            <em>real functionals</em>.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <iostream>

#include <thimble/math/numerical/RealFunctional.h>

using namespace std;

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Virtual destructor.
	 *
	 * @details
	 *            Essentially, this destructor does nothing and
	 *            is just provided.
	 */
	RealFunctional::~RealFunctional() {
	}

	/**
	 * @brief
	 *            Estimates the gradient of this functional at
	 *            the specified <i>n</i>-dimensional vector.
	 *
	 * @details
	 *            Unless overridden, this methods prints an
	 *            message to <code>stderr</code> indicating
	 *            that the method is not implemented and then
	 *            exits with status 'EXIT_FAILURE'. The implementation
	 *            of this method is due to a subclass which may make
	 *            use of
	 *            the \link NumericalDifferentiation::jacobi()\endlink
	 *            method.
	 *
	 * @param grad
	 *            An array capabale of storing <i>n</i> double values in
	 *            which the gradient of this funcational at
	 *            <i>(x[0],...,x[n-1])</i> will be stored where
	 *            <i>n</i>=\link getPreimageDimension()\endlink
	 *
	 * @param x
	 *            An array containing <i>n</i> double values being a valid
	 *            argument to the the \link eval()\endlink function.
	 *
	 * @see NumericalDifferentiation::jacobi()
	 */
	void RealFunctional::derive( double *grad , double *x ) const {
		cerr << "RealFunction::deriv: Not implemented." << endl;
		exit(EXIT_FAILURE);
	}
}




