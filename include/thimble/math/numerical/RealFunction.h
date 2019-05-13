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
 * @file RealFunction.h
 *
 * @brief
 *            Provides an abstract class for representing real-valued
 *            one-dimensional functions.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_REALFUNCTION_H_
#define THIMBLE_REALFUNCTION_H_

#include <thimble/dllcompat.h>
#include <thimble/math/numerical/RealFunctional.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Objects of this class are deemed to represent real-valued
	 *            one-dimensional functionals which map single real values
	 *            to single real values.
	 *
	 * @details
	 *            A real function is a special case of a real functional which
	 *            is modeled in this programming paradigm by deriving it
	 *            from the corresponding class from
	 *            the \link RealFunctional\endlink class and implementing some
	 *            of its purely virtual functions analogous.
	 */
	class THIMBLE_DLL RealFunction : public RealFunctional {

	public:

		/**
		 * @brief
		 *            Purely virtual function that is deemed to define the
		 *            real function.
		 *
		 * @details
		 *            If this real function maps the real value <i>x</i> to
		 *            the real value <i>y</i>, i.e., \f$x\mapsto y\f$, then
		 *            this function should return <i>y</i> as <i>x</i> has
		 *            been passed to an implementation of this function.
		 *
		 * @param x
		 *            The real value at which this real function is evaluated.
		 *
		 * @return
		 *            The evaluation of this real function at <i>x</i>.
		 */
		virtual double val( double x ) const = 0;

		/**
		 * @brief
		 *            Returns the preimage dimension of this real functional
		 *            which is 1 for real functions.
		 *
		 * @return
		 *            The preimage dimension of this real functional
		 *            which is 1.
		 */
		virtual int getPreimageDimension() const;

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
		virtual double eval( double *x ) const;
	};

}

#endif /* THIMBLE_REALFUNCTION_H_ */
