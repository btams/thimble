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
 * @file CTools.h
 *
 * @brief
 *            Provides a class that provides different convenience methods,
 *            e.g., for sorting C arrays.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_CTOOLS_H_
#define THIMBLE_CTOOLS_H_

#include <stdint.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Provides different convenience methods, e.g., for
	 *           sorting C arrays
	 *
	 * @details
	 *           Methods in this class are provided only as they are
	 *           needed for a specific purpose by THIMBLE.
	 */
	class CTools {

	private:

		/**
		 * @brief
		 *           Private standard constructor to avoid that
		 *           objects from this class are created.
		 */
		inline CTools() { }

	public:

		/**
		 * @brief
		 *           Sorts an array of <code>uint32_t</code> in
		 *           ascending order.
		 *
		 * @param array
		 *           Array containing <i>n</i> integers of type
		 *           <code>uint32_t</code>.
		 *
		 * @param n
		 *           Number of integers of type <code>uint32_t</code>
		 *           contained in <code>array</code>.
		 */
		static void sort( uint32_t *array , int n );
	};
}


#endif /* THIMBLE_CTOOLS_H_ */
