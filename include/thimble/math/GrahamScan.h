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
 * @file GrahamScan.h
 *
 * @brief
 *            Provides a mechanism for computing the convex hull containing
 *            a set of two-dimensional coordinates.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_GRAHAMSCAN_H_
#define THIMBLE_GRAHAMSCAN_H_

#include <vector>

#include <thimble/dllcompat.h>

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Provides a static method with which the convex
	 *            hull containing a set of two-dimensional coordinates
	 *            can be computed.
	 */
	class THIMBLE_DLL GrahamScan {

	private:

		/**
		 * @brief
		 *            Private standard constructor.
		 *
		 * @details
		 *            The constructor is provided as private to prevent users
		 *            from creating objects of this class.
		 */
		inline GrahamScan() {
		}

	public:

		/**
		 * @brief
		 *            Computes the convex hull containing a specified list
		 *            of two-dimensional coordinates.
		 *
		 * @param hull
		 *            Output vector that will contain successive
		 *            two-dimensional coordinates defining a convex polygon
		 *            containing the coordinates in <code>list</code>.
		 *
		 * @param list
		 *            Contains the two-dimensional coordinates of which the
		 *            convex hull is computed with this method.
		 */
		static void convexHull
			( std::vector< std::pair<double,double> > & hull ,
			  const std::vector< std::pair<double,double> > & list );
	};
}


#endif /* THIMBLE_GRAHAMSCAN_H_ */
