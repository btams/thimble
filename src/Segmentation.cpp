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
 * @file Segmentation.cpp
 *
 * @brief
 *            Implementation of the class provided by the 'Segmentation.h'
 *            header for estimating a fingerprint image's separation
 *            between foreground and background pixels.
 *
 * @author Benjamin Tams
 */

#include "config.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>

#include <thimble/math/GrahamScan.h>
#include <thimble/image/Otsu.h>
#include <thimble/image/Morphology.h>
#include <thimble/image/Labeler.h>
#include <thimble/image/DrawTools.h>
#include <thimble/finger/Segmentation.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

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
	Segmentation::Segmentation( int h ) {

		this->h = h;
	}

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
	void Segmentation::estimate
		( std::vector< std::pair<int,int> > & poly ,
		  const double *intensityImage , int m , int n ) const {

		bool *tmp = (bool*)malloc( m * n * sizeof(bool) );
		if ( tmp == NULL ) {
			cerr << "Segmentation::estimate: out of memory." << endl;
		}

		presegment(tmp,intensityImage,m,n);
		postsegment(poly,tmp,m,n);

		free(tmp);
	}

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
	void Segmentation::estimate
		( bool *foregroundImage ,
		  const double *intensityImage , int m , int n ) const {

		// Temporarily set all pixels as background, i.e., white (false).
		memset(foregroundImage,0,m*n);

		// Estimate the foreground pixel's surrounding polygon.
		vector< pair<int,int> > poly;
		estimate(poly,intensityImage,m,n);

		// Draw the polygon's inner points as on the foreground image
		// as pixels with value black (true)
		memset(foregroundImage,0,m*n*sizeof(bool));
		DrawTools::fillPolygon(foregroundImage,m,n,poly,true);
	}

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
	void Segmentation::presegment
		( bool *foregroundImage , const double *intensityImage ,
		  int m , int n ) const {

		int histogram[256];
		memset(histogram,0,256*sizeof(int));

		bool *tmp = (bool*)malloc( m * n * sizeof(bool) );
		if ( tmp == NULL ) {
			cerr << "Segmentation::precompute: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Compute histogram
		for ( int i = 0 ; i < m ; i++ ) {
			for ( int j = 0 ; j < n ; j++ ) {
				int k = (int)THIMBLE_ROUND(intensityImage[i*n+j] * 255.0);
				histogram[k] += 1;
			}
		}

		// Apply Otsu's method to the histogram
		double t = Otsu::threshold256(histogram) / 255.0;
		for ( int i = 0 ; i < m ; i++ ) {
			for ( int j = 0 ; j < n ; j++ ) {
				tmp[i*n+j] = (intensityImage[i*n+j]<=t);
			}
		}

		// Dilation
		CircleMask mask(this->h);
		mask.dilation(foregroundImage,tmp,m,n);

		// Select largest eight-connected region
		Labeler labeler(LC_EIGHT_CONNECTED);
		labeler.label(foregroundImage,m,n);
		labeler.get(foregroundImage,0);

		// Erosion
		mask = CircleMask(this->h+this->h);
		mask.dilation(tmp,foregroundImage,m,n);
		mask.erosion(foregroundImage,tmp,m,n);

		free(tmp);
	}

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
	void Segmentation::postsegment
		( std::vector< std::pair<int,int> > & poly ,
		  const bool *foregroundImage , int m , int n ) const {

		// Extract all locations from the image which are border points,
		// i.e. those of which four-connected neighborhood is not set.
		vector< pair<double,double> > hull;
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				if ( foregroundImage[y*n+x] ) {
					if ( y == 0 || x == 0 || y == m-1 || x == n-1 ||
							!foregroundImage[(y-1)*n+x] ||
                            !foregroundImage[(y+1)*n+x] ||
                            !foregroundImage[y*n+x-1] ||
                            !foregroundImage[y*n+x+1] ) {
						hull.push_back(pair<double,double>(x,y));
					}
				}
			}
		}

		// Graham scan to determine the convex hull
		GrahamScan::convexHull(hull,hull);

		poly.clear();
		poly.reserve(hull.size());

		for ( int i = 0 ; i < (int)hull.size(); i++ ) {
			int x , y;
			x = (int)THIMBLE_ROUND(hull[i].first);
			y = (int)THIMBLE_ROUND(hull[i].second);
			poly.push_back( pair<int,int>(x,y) );
		}
	}

}
