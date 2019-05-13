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
 * @file RealFunctional.h
 *
 * @brief
 *            Provides an abstract purely virtual class representing
 *            real-valued <i>n</i>-dimensional functions, i.e.,
 *            <em>real functionals</em>.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_REALFUNCTIONAL_H_
#define THIMBLE_REALFUNCTIONAL_H_

#include <thimble/dllcompat.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Abstract class of which subclass should implement
	 *           its purely virtual member functions properly to
	 *           let its objects represent real functionals, i.e.,
	 *           real-valued <i>n</i>-dimensional maps.
	 */
	class THIMBLE_DLL RealFunctional {

	public:

		/**
		 * @brief
		 *            Virtual destructor.
		 *
		 * @details
		 *            Essentially, this destructor does nothing and
		 *            is just provided.
		 */
		virtual ~RealFunctional();

		/**
		 * @brief
		 *            Returns the preimage dimension of this real functional.
		 *
		 * @details
		 *            A real function maps real vectors of fixed
		 *            length <i>n</i> to a single real value. This
		 *            function is deemed to return the value of
		 *            <i>n</i>.
		 *
		 * @return
		 *            The preimage dimension of this real functional.
		 */
		virtual int getPreimageDimension() const = 0;

		/**
		 * @brief
		 *            Evaluates this real functional at the specified
		 *            vector containing <i>n</i> double values.
		 *
		 * @param x
		 *            A vector of <i>n</i> valid double values at which
		 *            this real functional is to be evaluated where
		 *            <i>n</i>=\link getPreimageDimension()\endlink.
		 *
		 * @return
		 *            The evaluation of this functional at
		 *            <i>(x[0],...,x[n-1])</i>.
		 */
		virtual double eval( double *x ) const = 0;

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
		virtual void derive( double *grad , double *x ) const;
	};

}


#endif /* THIMBLE_REALFUNCTIONAL_H_ */
