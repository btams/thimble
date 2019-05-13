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
 * @file DrawTools.cpp
 *
 * @brief
 *            Implements the drawing operations provided by the
 *            \link DrawTools\endlink class.
 *
 * @author Benjamin Tams
 */

#include "config.h" // We need this to provide a rounding function
#include <cmath>
#include <vector>
#include <algorithm>

#include <thimble/image/DrawTools.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Generic sign function for a primitive signed data
	 *           type.
	 *
	 * @details
	 *           The function can be used to determine the sign of
	 *           a value being of the following primitive signed types:
	 *           <code>int</code>, <code>long</code>, <code>float</code>,
	 *           <code>double</code>, <code>long double</code> etc.
	 *
	 * @param v
	 *           Value of which sign is returned by the function.
	 *
	 * @return
	 *           1 if <code>v</code> is greater than zero; -1 if
	 *           <code>v</code> is smaller than zero; otherwiese,
	 *           if <code>v</code> is equals zero, the function returns
	 *           <code>0</code>.
	 */
	template<class T>
	static inline int signum( const T & v ) {

		if ( v < 0.0 ) {
			return -1;
		}

		if ( v > 0.0 ) {
			return 1;
		}

		return 0;
	}

	/**
	 * @brief
	 *           Draws a line on a binary image.
	 *
	 * @details
	 *           The method implements the <i>Bresenham algorithm</i> and
	 *           draws assigns the pixels of the <i>Bresenham line</i>
	 *           between <code>(y0,x0)</code> and <code>(y1,x1)</code> by
	 *           the value specified by <code>c</code> which is either
	 *           <code>true</code> or <code>false</code>.
	 *
	 * @param image
	 *           Array of <code>bool</code> of size at least
	 *           <code>m*n</code> encoding an image of height
	 *           <code>m</code> and width <code>n</code> in row-major
	 *           order.
	 *
	 * @param m
	 *           Number of rows of <code>image</code>.
	 *
	 * @param n
	 *           Number of columns of <code>image</code>.
	 *
	 * @param y0
	 *           Index of the first pixel's row.
	 *
	 * @param x0
	 *           Index of the first pixel's column.
	 *
	 * @param y1
	 *           Index of the second pixel's row.
	 *
	 * @param x1
	 *           Index of the second pixel's column.
	 *
	 * @param c
	 *           The value <code>true</code> or <code>false</code>
	 *           specifying the color of the Bresenham line.
	 *
	 * @warning
	 *           If <code>m</code> or <code>n</code> are smaller
	 *           than 0, calling this method runs into undocumented
	 *           behavior.
	 *
	 * @warning
	 *           If <code>image</code> has not been allocated to
	 *           hold at least <code>m*n</code>
	 */
	void DrawTools::drawLine
		( bool *image , int m , int n ,
		  int y0 , int x0 , int y1 , int x1 , bool c ) {

		// The method implements the well known 'Bresenham algorithm'.

		int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
		int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
		int err = dx+dy;

		for( ; ; ) {

			DrawTools::drawPoint(image,m,n,y0,x0,c);
			if (x0==x1 && y0==y1) {
		    	break;
		    }
		    int e2 = 2*err;
		    if (e2 >= dy) {
		    	err += dy;
		    	x0 += sx;
		    }
		    if (e2 <= dx) {
		    	err += dx;
		    	y0 += sy;
		    }
		}
	}

	/**
	 * @brief
	 *           Draws the exterior of a polygon on a binary image.
	 *
	 * @details
	 *           Each element of the vector <code>poly</code> defines
	 *           a pixel coordinate <code>(y[i],x[i])</code> where
	 *           <code>y[i]=poly[i].second</code> and
	 *           <code>x[i]=poly[i].first</code>.
	 *
	 *           The method successively draws the Bresenham line between
	 *           <code>(y[i-1],x[i-1])</code> and <code>(y[i],x[i])</code>
	 *           in the color <code>c</code> on the specified pixel data.
	 *           Furthermore, to close the boundary of the polynomial
	 *           the Bresenham line between <code>(y[n-1],x[n-1])</code>
	 *           and <code>(y[0],x[0])</code> is drawn where
	 *           <code>n=poly.size()</code>.
	 *
	 * @param image
	 *           Array of <code>bool</code> of size at least
	 *           <code>m*n</code> encoding an image of height
	 *           <code>m</code> and width <code>n</code> in row-major
	 *           order.
	 *
	 * @param m
	 *           Number of rows of <code>image</code>.
	 *
	 * @param n
	 *           Number of columns of <code>image</code>.
	 *
	 * @param poly
	 *           Edge points defining the specified polygon.
	 *
	 * @param c
	 *           The value <code>true</code> or <code>false</code>
	 *           specifying the color of the polygon's exterior to be
	 *           drawn.
	 *
	 * @warning
	 *           If <code>m</code> or <code>n</code> are smaller
	 *           than 0, calling this method runs into undocumented
	 *           behavior.
	 *
	 * @warning
	 *           If <code>image</code> has not been allocated to
	 */
	void DrawTools::drawPolygon
	( bool *image , int m , int n ,
	  const vector< pair<int,int> > & poly , bool c ) {

		for ( int i0 = 0 ; i0 < (int)poly.size() ; i0++ ) {

			int i1 = i0+1;
			if ( i1 == (int)poly.size() ) {
				i1 = 0;
			}

			int x0 , y0 , x1 , y1;
			x0 = poly[i0].first;
			y0 = poly[i0].second;
			x1 = poly[i1].first;
			y1 = poly[i1].second;

			drawLine(image,m,n,y0,x0,y1,x1,c);
		}
	}

	/**
	 * @brief
	 *           Draws the interior of a polygon on a binary image.
	 *
	 * @details
	 *           Each element of the vector <code>poly</code> defines
	 *           a pixel coordinate <code>(y[i],x[i]) where
	 *           <code>y[i]=poly[i].second</code> and
	 *           <code>x[i]=poly[i].first</code>.
	 *
	 * @param image
	 *           Array of <code>bool</code> of size at least
	 *           <code>m*n</code> encoding an image of height
	 *           <code>m</code> and width <code>n</code> in row-major
	 *           order.
	 *
	 * @param m
	 *           Number of rows of <code>image</code>.
	 *
	 * @param n
	 *           Number of columns of <code>image</code>.
	 *
	 * @param poly
	 *           Edge points defining the specified polygon.
	 *
	 * @param c
	 *           The value <code>true</code> or <code>false</code>
	 *           specifying the color of the polygon's interior to be
	 *           drawn.
	 *
	 * @warning
	 *           If <code>m</code> or <code>n</code> are smaller
	 *           than 0, calling this method runs into undocumented
	 *           behavior.
	 *
	 * @warning
	 *           If <code>image</code> has not been allocated to
	 *           hold at least <code>m*n</code>
	 */
	void DrawTools::fillPolygon
	( bool *image , int m , int n ,
	  const vector< pair<int,int> > & poly , bool c ) {

		if ( poly.size() == 0 ) {
			// Nothing to do
			return;
		}

		vector<int> nodes;

		// Determine the polygon's edges minimal and maximal X-value
		// 'minX' and 'maxX', respectively
		int minX , maxX;
		minX = poly[0].first;
		maxX = poly[0].first;
		for ( int i = 1 ; i < (int)poly.size() ; i++ ) {
			if ( poly[i].first < minX ) {
				minX = poly[i].first;
			}
			if ( poly[i].first > maxX ) {
				maxX = poly[i].first;
			}
		}

		// Process each rows separately
		for ( int x = minX ; x <= maxX ; x++ ) {

			// Collect all y-values after which the color changes.
			// These can be systematically determined by computing each
			// polygon's line piece intersection with the current row
			// defined by 'x'.
			nodes.clear();
			for ( int i = 0 ; i < (int)poly.size() ; i++ ) {
				int i0;
				if ( i == 0 ) {
					i0 = (int)poly.size()-1;
				} else {
					i0 = i-1;
				}
				int xr = poly[i].first-poly[i0].first;
				if ( xr != 0 ) {
					if ( signum(xr) == signum(x-poly[i0].first) &&
						 abs(x-poly[i0].first) <= abs(xr) ) {

						int y = (int)THIMBLE_ROUND
								( (double)(x-poly[i0].first) *
								  (double)(poly[i].second-poly[i0].second) /
								  (double)xr + (double)poly[i0].second );

						bool isContained = false;
						for ( int j = 0 ; j < (int)nodes.size() ; j++ ) {
							if ( nodes[j] == y ) {
								isContained = true;
							}
						}
						if ( !isContained ) {
							nodes.push_back(y);
						}
					}
				}
			}

			// The crossover indices are sorted.
		    sort( nodes.begin() , nodes.end() );

		    // Set pixel in the current rows laying between nodes
		    // of the accessed polygon.
		    for ( int i = 0 ; i+1 < (int)nodes.size() ; i += 2 ) {
		    	for ( int y = nodes[i] ; y <= nodes[i+1] ; y++ ) {
		    		drawPoint(image,m,n,y,x,c);
		    	}
		    }
		}
	}

}



