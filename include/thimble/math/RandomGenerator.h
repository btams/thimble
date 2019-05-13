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
 * @file RandomGenerator.h
 *
 * @brief
 *            Provides a portable class for generating pseudo-random
 *            number sequences via a seed.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_RANDOMGENERATOR_H_
#define THIMBLE_RANDOMGENERATOR_H_

#include <stdint.h>

#include <thimble/dllcompat.h>
#include <thimble/security/SHA.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Objects of this class represent states for a generating
	 *            sequences of pseudo-random number.
	 */
	class THIMBLE_DLL RandomGenerator {

	public:

		/**
		 * @brief
		 *           Creates a pseudo-random number generator whose initial
		 *           seed is given by 160-bit value.
		 *
		 * @param seed
		 *           An array containing 20 integers of type
		 *           <code>uint8_t</code>.
		 *
		 * @warning
		 *           If <code>seed</code> is not a well-defined array
		 *           containing 20 integers of type <code>uint8_t</code>
		 *           the program runs into undocumented behavior.
		 */
		RandomGenerator( uint8_t seed[20] );

		/**
		 * @brief
		 *           Generates the next pseudo-random number.
		 *
		 * @details
		 *           The function derives changes its state by replacing
		 *           it with its own \link SHA\endlink hash value. Then
		 *           a 64-bit number is derived from the new state such
		 *           that it looks randomly.
		 *
		 * @return
		 *           A pseudo-random unsigned 64-bit integer.
		 */
		uint64_t rand();

	private:

		/**
		 * @brief
		 *           An object for computing SHA-1 hash values.
		 */
        SHA sha;

        /**
         * @brief
         *           The state of this pseudo-random number generator which,
         *           initially, is assigned with the values of the seed.
         */
		uint8_t hash[20];
	};
}

#endif /* THIMBLE_RANDOMGENERATOR_H_ */
