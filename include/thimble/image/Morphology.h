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
 * @file Morphology.h
 *
 * @brief
 *            Provides classes of which objects represent structuring
 *            elements for applying morphological operations.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_MORPHOLOGY_H_
#define THIMBLE_MORPHOLOGY_H_

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Basic class for objects that represent structuring elements
	 *            with which morpholical operations can be performed.
	 */
	class THIMBLE_DLL MorphologicMask {

	public:

		/**
		 * @brief
		 *            Virtual destructor.
		 */
		virtual ~MorphologicMask();

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
		void erosion
			( bool *out , const bool *in , int m , int n ) const;

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
		void dilation
			( bool *out , const bool *in , int m , int n ) const;

		/**
		 * @brief
		 *            Draws this structuring element to the specified
		 *            image at the specified pixel with specified color.
		 *
		 * @details
		 *            This method is used to draw this structuring element
		 *            at every black pixel on a
		 *            basic \link dilation() dilation operation\endlink.
		 *
		 *            This method is abstract and is implemented for a
		 *            specific structuring element represented by
		 *            objects that have been created from non-abstract
		 *            subclass of this base-class. The implementation of
		 *            this method, together with the implementation of
		 *            the \link isContained()\endlink function, defines this
		 *            structuring element.
		 *
		 *            The structuring element is drawn on the image
		 *            <code>image</code> of height <code>m</code> and
		 *            with <code>n</code> where the pixel at coordinate
		 *            <code>(y,x)</code> is encoded at
		 *            <code>image[y*n+x]</code> where <code>true</code>
		 *            encodes a black and <code>false</code> a white
		 *            pixel.
		 *
		 * @param image
		 *            Image to which this structuring element is drawn.
		 *
		 * @param m
		 *            Height of the image to which this structuring element
		 *            is drawn.
		 *
		 * @param n
		 *            Width of the image to which this structuring element
		 *            is drawn.
		 *
		 * @param y
		 *            y-coordinate of the pixel at which this structuring
		 *            element is drawn.
		 *
		 * @param x
		 *            x-coordinate of the pixel at which this structuring
		 *            element is drawn.
		 *
		 * @param c
		 *            Color with which this structuring element is drawn.
		 */
		virtual void draw
			( bool *image , int m , int n , int y , int x , bool c ) const = 0;

		/**
		 * @brief
		 *            Checks whether this structuring element fully covers
		 *            black pixels at the specified pixel coordinate of
		 *            the specified image.
		 *
		 * @details
		 *            On the basic morphological \link erosion()
		 *            erosion operation\endlink this function is used to
		 *            determine whether an image pixel's value should be
		 *            left <code>true</code> (i.e., if this structuring
		 *            element fully covers black (true) pixels) or if
		 *            should be replaced by a <code>false</code> (white)
		 *            pixel.
		 *
		 *            This function is abstract and is implemented for a
		 *            specific structuring element represented by
		 *            objects that have been created from non-abstract
		 *            subclass of this base-class. The implementation of
		 *            this function, together with the implementation of
		 *            the \link draw()\endlink method, defines this
		 *            structuring element.
		 *
		 * @param image
		 *            Image of type <code>bool</code> of height
		 *            <code>m</code> and width <code>n</code> storing
		 *            the value of the pixel at <code>(y,x)</code> in
		 *            <code>image[y*n+x]</code>.
		 *
		 * @param m
		 *            Height of the image.
		 *
		 * @param n
		 *            Width of the image.
		 *
		 * @param y
		 *            y-coordinate of the pixel.
		 *
		 * @param x
		 *            x-coordinate of the pixel.
		 *
		 * @return
		 *            <code>true</code> if this structuring element fully
		 *            covers black pixel and, otherwise, <code>false</code>.
		 */
		virtual bool isContained
			( const bool *image , int m , int n , int y , int x ) const = 0;
	};

	/**
	 * @brief
	 *            Objects of this class represent filled circles as
	 *            structuring elements in morphological operations.
	 */
	class THIMBLE_DLL CircleMask : public MorphologicMask {

	public:

		/**
		 * @brief
		 *            Creates a filled circle for morphological operations
		 *            with specified radius.
		 *
		 * @param r
		 *            The radius of the circle.
		 *
		 * @warning
		 *            If <code>r</code> is negative, the program runs into
		 *            undocumented behavior.
		 */
		CircleMask( int r );

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
		virtual void draw
			( bool *image , int m , int n , int y , int x , bool c ) const;

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
		virtual bool isContained
			( const bool *image , int m , int n , int y , int x ) const;

	private:

		/**
		 * @brief
		 *            The radius of this circle.
		 */
		int r;
	};
}


#endif /* THIMBLE_MORPHOLOGY_H_ */
