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
 * @file HexagonalGrid.h
 *
 * @brief
 *            Provides a mechanism for computing the points of a hexagonal
 *            grid that lay within a specified region.
 *
 * @author Benjamin Tams
 *
 * @see thimble::HexagonalGrid
 */

#ifndef THIMBLE_HEXAGONALGRID_H_
#define THIMBLE_HEXAGONALGRID_H_

#include <vector>

#include <thimble/dllcompat.h>

#ifdef THIMBLE_BUILD_DLL
template struct THIMBLE_DLL std::pair<double,double>;
template class THIMBLE_DLL std::allocator<std::pair<double,double> >;
template class THIMBLE_DLL std::vector<std::pair<double,double>,std::allocator<std::pair<double,double> > >;
#endif

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of this class contain the points of a hexagonal
	 *            grid within a specified region.
	 *
	 * @details
	 *            For example, assume we want to compute the points of a
	 *            hexagonal grid with a minimal distance of \f$\lambda=29\f$
	 *            within the region \f$[0,295]\times[0,559]\f$. Through
	 *            <pre>
	 *             HexagonalGrid grid(0.0,0.0,295.0,559.0,29.0);
	 *            </pre>
	 *            an instance is created that holds the points centered
	 *            in the specified region. The points can be
	 *            accessed via
	 *            <pre>
	 *           std::vector< std::pair<double,double> > p = grid.getPoints();
	 *            </pre>
	 *            which is a vector of <code>double</code> pairs whose
	 *            respective first member correspond to the abscissa
	 *            coordinate and the second member to the ordinate coordinate.
	 */
	class THIMBLE_DLL HexagonalGrid {

	private:

		/**
		 * @brief
		 *            Points of the hexagonal grid centered in the region
		 *            specified on the constructor's call.
		 *
		 * @see getPoints() const
		 */
		std::vector< std::pair<double,double> > points;

	public:

		/**
		 * @brief
		 *            Computes the points of a hexagonal grid centered
		 *            in the specified region and that are of specified
		 *            minimal distance.
		 *
		 * @details
		 *            The points that are computed are centered in the region
		 *            \f$[x0,x1]\times[y0,y1]\f$ and of minimal distance
		 *            \f$\lambda\f$.
		 *
		 * @param x0
		 *            Specifies the minimal <i>x</i>-coordinate of the region.
		 *
		 * @param y0
		 *            Specifies the minimal <i>y</i>-coordinate of the region.
		 *
		 * @param x1
		 *            Specifies the maximal <i>x</i>-coordinate of the region.
		 *
		 * @param y1
		 *            Specifies the maximal <i>y</i>-coordinate of the region.
		 *
		 * @param lambda
		 *            Minimal distance between two distinct hexagonal
		 *            grid points.
		 *
		 * @warning
		 *            If <i>x0>y0</i> or <i>y0>y1</i> or if
		 *            <code>lambda<=0</code> the constructor prints an error
		 *            message to <code>stderr</code> and exits with status
		 *            -1.
		 */
		HexagonalGrid
		( double x0 , double y0 , double x1 , double y1 , double lambda );

		/**
		 * @brief
		 *            Access the points of the hexagonal grid points.
		 *
		 * @return
		 *            The points of the hexagonal grid points.
		 */
		inline const std::vector< std::pair<double,double> > & getPoints() const {
			return this->points;
		}
	};
}

#endif /* THIMBLE_HEXAGONALGRID_H_ */
