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
 * @file DrawTools.h
 *
 * @brief
 *            Provides a class implementing elementary drawing operations
 *            as they are needed by the library's main functions.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_DRAWTOOLS_H_
#define THIMBLE_DRAWTOOLS_H_

#include <vector>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Class wrapping static methods for some elementary drawing
	 *            operations.
	 */
	class THIMBLE_DLL DrawTools {

	public:

		/**
		 * @brief
		 *           Set the value of a pixel on a binary image.
		 *
		 * @details
		 *           If the pixel is within the image's region, it
		 *           will be assigned by the specified value
		 *           <code>c</code> which is either <code>true</code>
		 *           or <code>false</code>; otherwise, if the accessed
		 *           pixel is outside the image's region, calling this
		 *           method has no effect.
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
		 * @param y
		 *           Index of the pixel's row being set.
		 *
		 * @param x
		 *           Index of the pixel's column being set.
		 *
		 * @param c
		 *           The value <code>true</code> or <code>false</code>
		 *           that the accessed pixel is assigned with.
		 *
		 * @warning
		 *           If <code>m</code> or <code>n</code> are smaller
		 *           than 0, calling this method runs into undocumented
		 *           behavior.
		 *
		 * @warning
		 *           If <code>image</code> has not been allocated to
		 *           hold at least <code>m*n</code>
		 *
		 */
		static inline void drawPoint
			( bool *image , int m , int n , int y , int x , bool c ) {

			if ( y >= 0 && y < m && x >= 0 && x < n ) {
				image[y*n+x] = c;
			}
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
		static void drawLine
			( bool *image , int m , int n ,
			  int y0 , int x0 , int y1 , int x1 , bool c );

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
		static void drawPolygon
			( bool *image , int m , int n ,
			  const std::vector< std::pair<int,int> > & poly , bool c );

		/**
		 * @brief
		 *           Draws the interior of a polygon on a binary image.
		 *
		 * @details
		 *           Each element of the vector <code>poly</code> defines
		 *           a pixel coordinate <code>(y[i],x[i])</code> where
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
		static void fillPolygon
			( bool *image , int m , int n ,
			  const std::vector< std::pair<int,int> > & poly , bool c );
	};
}


#endif /* THIMBLE_DRAWTOOLS_H_ */
