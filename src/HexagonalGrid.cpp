/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013 Benjamin Tams
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
 * @file HexagonalGrid.cpp
 *
 * @brief
 *            Implements functionalities from 'HexagonalGrid.h' which provides
 *            a mechanism for computing the points of a hexagonal grid that
 *            lay within a specified region.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

#include <thimble/math/HexagonalGrid.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Computes the points of a hexagonal grid centered
	 *            in the specified region and that are of specified
	 *            minimal distance.
	 *
	 * @details
	 *            see 'HexagonalGrid.h'
	 */
	HexagonalGrid::HexagonalGrid
	( double x0 , double y0 , double x1 , double y1 , double lambda ) {

		// Will be frequently used to compute successive hexagonal grid
		// points through accumulation.
		static const double DY =
				tan( 60.0/360.0 * (M_PI+M_PI) ) / 2.0;

		// Check if the arguments are reasonable.
		if ( x0 > x1 || y0 > y1 || lambda <= 0.0 ) {
			cerr << "HexagonalGrid: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// Distance of two grid row's 'y'-coordinate.
		double dy = lambda*DY;

		// Keeps track of the maximal 'x'- and 'y'-coordinates of the
		// temporarily computed grid points. After all points have been
		// computed, this information will be used to center the points
		// in the specified region.
		double xmax = x0, ymax = y0;

		// Two successive rows of hexagonal grid points will not have
		// the same 'x'-coordinate but they are somewhat in between of
		// two grid point's 'x'-coordinate of the previous row. By
		// iterating through the grids rows this boolean value is alternately
		// set 'true' and 'false'; in case it is 'true' the 'x'-coordinates
		// are shifted by an offset of 'lambda*0.5';
		bool shiftX = false;

		// Iteration over the grid rows in the region, i.e. iteration
		// over the 'y'-coordinates.
		for ( double y = y0 ; y <= y1 ; y += dy , shiftX = !shiftX) {

			// Check if update maximal 'y'-coordinate of grid points, which
			// will be used later to center the grid points in the region.
			if ( y > ymax ) {
				ymax = y;
			}

			// Initial 'x'-coordinate of the row which is shifted by
			// '0.5*lambda' for every second row as discussed above.
			double x = x0;
			if ( shiftX ) {
				x += lambda * 0.5;
			}

			// Iteration over grid columns.
			for ( ; x <= x1 ; x += lambda ) {

				// Append the point in the list.
				this->points.push_back( pair<long double,long double>(x,y) );

				// Check if update maximal 'x'-coordinate of grid points, which
				// will be used later to center the grid points in the region.
				if ( x > xmax ) {
					xmax = x;
				}
			}
		}

		// Center the grid points in the specified region.
		double mx = (x1-xmax)*0.5, my = (y1-ymax)*0.5;
		for ( int k = 0 ; k < (int)(this->points.size()) ; k++ ) {
			this->points[k].first  += mx;
			this->points[k].second += my;
		}

	}

}

