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
 * @file Morphology.cpp
 *
 * @brief
 *            Implementation of basic morphological operations and
 *            classes of which objects represent structuring
 *            elements as provided by the 'Morphology.h' header.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/image/Morphology.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Virtual destructor.
	 */
	MorphologicMask::~MorphologicMask() {
	}

	/**
	 * @brief
	 *            Performs the basic morphological erosion operation
	 *            with this structuring element.
	 *
	 * @details
	 *            This method applies the morphological erosion
	 *            operation to an image of height <code>m</code>
	 *            and width <code>n</code> of which value at
	 *            pixel <code>(y,x)</code> is stored at
	 *            <code>in[y*n+x]</code> where <code>true</code>
	 *            encodes a black and <code>false</code> a white
	 *            pixel. The result is stored in the image
	 *            <code>out</code>.
	 *
	 * @param out
	 *            Output binary image that represents the result
	 *            of the erosion.
	 *
	 * @param in
	 *            Input binary image to which the erosion operation
	 *            is applied.
	 *
	 * @param m
	 *            Height of the images represented by
	 *            <code>out</code> and <code>in</code>.
	 *
	 * @param n
	 *            Width of the images represented by
	 *            <code>out</code> and <code>in</code>.
	 *
	 * @warning
	 *            If <code>in</code> does not contain <code>m*n</code>
	 *            values of type <code>bool</code> or if <code>out</code>
	 *            cannot hold <code>m*n</code> values of type
	 *            <code>bool</code> or if <code>m</code> or <code>n</code>
	 *            are negative, the program runs into undocumented
	 *            behavior.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void MorphologicMask::erosion
		( bool *out , const bool *in , int m , int n ) const {

		if ( out == in ) {
			bool *tmp = (bool*)malloc( m * n * sizeof(bool) );
			if ( tmp == NULL ) {
				cerr << "MorphologicMask::erosion: out of memory." << endl;
				exit(EXIT_FAILURE);
			}
			memcpy(tmp,in,m*n*sizeof(bool));
			erosion(out,tmp,m,n);
			free(tmp);
			return;
		}

		memcpy(out,in,m*n*sizeof(bool));

		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				if ( !isContained(in,m,n,y,x) ) {
					out[y*n+x] = false;
				}
			}
		}

	}

	/**
	 * @brief
	 *            Performs the basic morphological erosion operation
	 *            with this structuring element.
	 *
	 * @details
	 *            This method applies the morphological dilation
	 *            operation to an image of height <code>m</code>
	 *            and width <code>n</code> of which value at
	 *            pixel <code>(y,x)</code> is stored at
	 *            <code>in[y*n+x]</code> where <code>true</code>
	 *            encodes a black and <code>false</code> a white
	 *            pixel. The result is stored in the image
	 *            <code>out</code>.
	 *
	 * @param out
	 *            Output binary image that represents the result
	 *            of the image represented by <code>in</code> to
	 *            which the dilation operation has been applied.
	 *
	 * @param in
	 *            Input binary image to which the dilation operation
	 *            is applied.
	 *
	 * @param m
	 *            Height of the images represented by
	 *            <code>out</code> and <code>in</code>.
	 *
	 * @param n
	 *            Width of the images represented by
	 *            <code>out</code> and <code>in</code>.
	 *
	 * @warning
	 *            If <code>in</code> does not contain <code>m*n</code>
	 *            values of type <code>bool</code> or if <code>out</code>
	 *            cannot hold <code>m*n</code> values of type
	 *            <code>bool</code> or if <code>m</code> or <code>n</code>
	 *            are negative, the program runs into undocumented
	 *            behavior.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void MorphologicMask::dilation
		( bool *out , const bool *in , int m , int n ) const {

		if ( out == in ) {
			bool *tmp = (bool*)malloc( m * n * sizeof(bool) );
			if ( tmp == NULL ) {
				cerr << "MorphologicMask::dilation: out of memory." << endl;
				exit(EXIT_FAILURE);
			}
			memcpy(tmp,in,m*n*sizeof(bool));
			dilation(out,tmp,m,n);
			free(tmp);
			return;
		}

		memset(out,0,m*n*sizeof(bool));

		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				if ( in[y*n+x] ) {
					draw(out,m,n,y,x,true);
				}
			}
		}
	}


	/**
	 * @brief
	 *            Creates a filled circle for morphological operations
	 *            with specified radius.
	 *
	 * @param r
	 *            The radius of the circle.
	 */
	CircleMask::CircleMask( int r ) {
		this->r = r;
	}

	/**
	 * @brief
	 *            Draws a filled circle to an image of which center
	 *            is placed at the specified pixel.
	 *
	 * @details
	 *            The circle is drawn to an image of height
	 *            <code>m*n</code> where the value of the pixel
	 *            <code>(y,x)</code> is stored at
	 *            <code>image[y*n+x]</code>.
	 *
	 * @param image
	 *            The image to which this circle is drawn.
	 *
	 * @param m
	 *            The height of the image.
	 *
	 * @param n
	 *            The width of the image.
	 *
	 * @param y
	 *            y-coordinate of the drawn circle's center.
	 *
	 * @param x
	 *            x-coordinate of the drawn circle's center.
	 *
	 * @param c
	 *            The color with which this circle is drawn.
	 *
	 * @warning
	 *            If <code>image</code> does not contain <code>m*n</code>
	 *            values of type <code>bool</code> or if <code>m</code> or
	 *            <code>n</code> is negative, the program runs into
	 *            undocumented behavior.
	 */
	void CircleMask::draw
		( bool *image , int m , int n , int y , int x , bool c ) const {

		for ( int i = -r ; i <= r ; i++ ) {
			if ( i+y >= 0 && i+y < m ) {
				for ( int j = -r ; j <= r ; j++ ) {
					if ( j+x >= 0 && j+x < n ) {
						if ( i*i+j*j <= r*r ) {
							image[(i+y)*n+j+x] = c;
						}
					}
				}
			}
		}
	}

	/**
	 * @brief
	 *            Checks whether this circle with specified center is
	 *            fully contained in a black (true) area of an image.
	 *
	 * @param image
	 *            Row-wise pixel data of the image.
	 *
	 * @param m
	 *            Height of the image.
	 *
	 * @param n
	 *            Width of the image.
	 *
	 * @param y
	 *            y-coordinate of the circle center.
	 *
	 * @param x
	 *            x-coordinate of the circle center.
	 *
	 * @warning
	 *            If <code>image</code> does not contain <code>m*n</code>
	 *            values of type <code>bool</code> or if <code>m</code> or
	 *            <code>n</code> is negative, the program runs into
	 *            undocumented behavior.
	 */
	bool CircleMask::isContained
		( const bool *image , int m , int n , int y , int x ) const {

		for ( int i = -r ; i <= r ; i++ ) {
			if ( i+y >= 0 && i+y < m ) {
				for ( int j = -r ; j <= r ; j++ ) {
					if ( j+x >= 0 && j+x < n ) {
						if ( i*i+j*j <= r*r ) {
							if ( !image[(i+y)*n+j+x] ) {
								return false;
							}
						}
					} else {
						return false;
					}
				}
			} else {
				return false;
			}
		}

		return true;
	}
}




