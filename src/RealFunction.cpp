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
 * @file RealFunction.cpp
 *
 * @brief
 *            Implements an abstract class for representing real-valued
 *            one-dimensional functionals.
 *
 * @author Benjamin Tams
 */

#include <thimble/math/numerical/RealFunctional.h>
#include <thimble/math/numerical/RealFunction.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Returns the preimage dimension of this real functional
	 *            which is 1 for real functions.
	 *
	 * @return
	 *            The preimage dimension of this real functional
	 *            which is 1.
	 */
	int RealFunction::getPreimageDimension() const {
		return 1;
	}

	/**
	 * @brief
	 *            Evaluates this real function at the
	 *            specified one-dimensional array.
	 *
	 * @details
	 *            The implementation of this function is wrapped around
	 *            the purely virtual \link val()\endlink function of this
	 *            object.
	 *
	 * @param x
	 *            An array containing (at least) one valid double
	 *            value that can be passed as argument to
	 *            the \link val()\endlink function.
	 *
	 * @return
	 *            <code>val(x[0])</code>.
	 */
	double RealFunction::eval( double *x ) const {
		return val(x[0]);
	}
}
