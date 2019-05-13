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
 * @file RandomGenerator.cpp
 *
 * @brief
 *            Implementation of a portable class for generating pseudo-random
 *            number sequences via a seed as provided by the
 *            'RandomGenerator.h' header.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstring>

#include <thimble/security/SHA.h>
#include <thimble/math/RandomGenerator.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Creates a pseudo-random number generator whose initial
	 *           seed is given by a 160-bit value.
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
	RandomGenerator::RandomGenerator( uint8_t seed[20] ) {

		// Copy the seed to the array hash that represents
		// the state of the random number generator
		memcpy(this->hash,seed,20);
	}

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
	uint64_t RandomGenerator::rand() {

		// Replace hash by its own SHA-1 hash value.
        this->sha.hash(this->hash,this->hash,20);

        // The first 8 bytes form the new 64-bit number (big-endian).
		uint64_t v = 0;
		for ( int j = 0 ; j < 8 ; j++ ) {
			v <<= 8;
			v += this->hash[j];
		}

		return v;
	}
}
