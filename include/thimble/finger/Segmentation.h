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
 * @file Segmentation.h
 *
 * @brief
 *            Provides a class for estimating a fingerprint image's separation
 *            between foreground and background pixels.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_SEGMENTATION_H_
#define THIMBLE_SEGMENTATION_H_

#include <vector>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            A class for estimating a fingerprint image's separation
	 *            between foreground and background pixels.
	 */
	class THIMBLE_DLL Segmentation {

	public:

		/**
		 * @brief
		 *             Creates an instance for estimating a fingerprint
		 *             image's separation between foreground and background
		 *             pixels.
		 *
		 * @param h
		 *             Controls the radius of circles that are used as
		 *             the structuring elements during internal morphologic
		 *             operations.
		 *
		 * @warning
		 *             If <code>h</code> is negative, the use of the object
		 *             runs into undocumented behavior.
		 */
		Segmentation( int h );

		/**
		 * @brief
		 *             Estimates a polygon that of which inner points are
		 *             the specified fingerprint image's estimated foreground
		 *             pixels.
		 *
		 * @details
		 *             This method calls \link presegment()\endlink which
		 *             estimates a connected region as a coarse separation
		 *             between foreground and background pixels. Then, the
		 *             method \link postsegment()\endlink is called which
		 *             computes the polygon that forms the convex hull of
		 *             the result from the pre-segmentation.
		 *
		 * @param poly
		 *             The output polygon of which inner points are the
		 *             estimated foreground pixels.
		 *
		 * @param intensityImage
		 *             The double-valued intensities of the specified image;
		 *             where the intensity of the pixel at <code>(y,x)</code>
		 *             is stored at <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *             The height of the specified intensity image.
		 *
		 * @param n
		 *             The width of the specified intensity image.
		 *
		 * @warning
		 *             In the following cases, calling this method runs into
		 *             undocumented behavior:
		 *             <ul>
		 *              <li><code>m</code> is negative</li>
		 *              <li><code>n</code> is negative</li>
		 *              <li>
		 *               <code>intensityImage</code> does not contain
		 *               <code>m*n</code> valid <code>double</code> values.
		 *              </li>
		 *             </ul>
		 */
		void estimate
			( std::vector< std::pair<int,int> > & poly ,
			  const double *intensityImage , int m , int n ) const;

		/**
		 * @brief
		 *             Estimates the foreground pixels from the specified
		 *             fingerprint image.
		 *
		 * @details
		 *             This convenience method
		 *             calls \link estimate()\endlink
		 *             and draws the resulting polygon's inner points to
		 *             <code>foregroundImage</code>.
		 *
		 * @param foregroundImage
		 *             Output image that can hold <code>m*n</code> values of
		 *             type <code>bool</code>. If the pixel (y,x) was
		 *             estimated as a foreground pixel, the value
		 *             <code>foregroundImage[y*n+x]</code> is
		 *             <code>true</code> and, otherwise, <code>false</code>.
		 *
		 * @param intensityImage
		 *             The double-valued intensities of the specified image;
		 *             where the intensity of the pixel at <code>(y,x)</code>
		 *             is stored at <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *             The height of the specified intensity image.
		 *
		 * @param n
		 *             The width of the specified intensity image.
		 *
		 * @warning
		 *             In the following cases, calling this method runs into
		 *             undocumented behavior:
		 *             <ul>
		 *              <li><code>m</code> is negative</li>
		 *              <li><code>n</code> is negative</li>
		 *              <li>
		 *               <code>intensityImage</code> does not contain
		 *               <code>m*n</code> valid <code>double</code> values.
		 *              </li>
		 *             </ul>
		 */
		void estimate
			( bool *foregroundImage ,
			  const double *intensityImage , int m , int n ) const;

		/**
		 * @brief
		 *             Estimates a coarse separation between a fingerprint
		 *             image's foreground and background pixels.
		 *
		 * @details
		 *             This method performs <em>Otsu thresholding</em> and
		 *             the result is dilated using a circle of radius as
		 *             specified at creation of this object; then, the largest
		 *             eight-connected region is selected to contain a coarse
		 *             estimation of the fingerprint image's foreground
		 *             pixels.
		 *
		 * @param foregroundImage
		 *             Output image that can hold <code>m*n</code> values of
		 *             type <code>bool</code>. If the pixel (y,x) was
		 *             estimated as a foreground pixel, the value
		 *             <code>foregroundImage[y*n+x]</code> is
		 *             <code>true</code> and, otherwise, <code>false</code>.
		 *
		 * @param intensityImage
		 *             The double-valued intensities of the specified image;
		 *             where the intensity of the pixel at <code>(y,x)</code>
		 *             is stored at <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *             The height of the specified intensity image.
		 *
		 * @param n
		 *             The width of the specified intensity image.
		 *
		 * @warning
		 *             In the following cases, calling this method runs into
		 *             undocumented behavior:
		 *             <ul>
		 *              <li><code>m</code> is negative</li>
		 *              <li><code>n</code> is negative</li>
		 *              <li>
		 *               <code>intensityImage</code> does not contain
		 *               <code>m*n</code> valid <code>double</code> values.
		 *              </li>
		 *             </ul>
		 */
		void presegment
			( bool *foregroundImage ,
			  const double *intensityImage , int m , int n ) const;

		/**
		 * @brief
		 *             Computes the convex hull of an initial fingerprint
		 *             image's foreground estimation as a polygon
		 *             to contain the (final) fingerprint image's
		 *             foreground pixels.
		 *
		 * @param poly
		 *             The output (convec) polygon of which inner points are
		 *             the estimated foreground pixels.
		 *
		 * @param foregroundImage
		 *             An array containing <code>m*n</code> valid values of
		 *             type <code>bool</code> where
		 *             <code>foregroundImage[y*n+x]</code> is
		 *             <code>true</code> if <code>(y,x)</code> was estimated
		 *             as a foreground pixel and, otherwise, if
		 *             <code>(y,x)</code> was estimated as a background pixel,
		 *             the value of <code>foregroundImage[y*n+x]</code> is
		 *             <code>false</code>.
		 *
		 * @param m
		 *             The height of the specified intensity image.
		 *
		 * @param n
		 *             The width of the specified intensity image.
		 */
		void postsegment
			( std::vector< std::pair<int,int> > & poly ,
			  const bool *foregroundImage , int m , int n ) const;

	private:

		/**
		 * @brief
		 *            Controls the size of structuring elements used for
		 *            internal morphologic operations.
		 */
		int h;
	};
}


#endif /* THIMBLE_SEGMENTATION_H_ */
