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
 *            Implementation of a class that provides different convenience
 *            methods, e.g., for sorting C arrays, as provided by the
 *            'CTools.h' header file.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>

#include <thimble/misc/CTools.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

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
	void CTools::sort( uint32_t *array , int n ) {

		for ( int i = 0 ; i < n-1 ; i++ ) {
			for ( int j = 0 ; j < n - 1 - i ; j++ ) {
				if ( array[j] > array[j+1] ) {
					uint32_t tmp = array[j+1];
					array[j+1] = array[j];
					array[j] = tmp;
				}
			}
		}
	}
}




