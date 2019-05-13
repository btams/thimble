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
 * @file BitImage.h
 *
 * @brief
 *            Provides a class for reading and writing binary images in
 *            PBM format.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_BITIMAGE_H_
#define THIMBLE_BITIMAGE_H_

#include <stdint.h>
#include <cstdio>

#include <string>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Class representing a binary image.
	 */
	class THIMBLE_DLL BitImage {

	public:
		/**
		 * @brief
		 *           Standard constructor.
		 *
		 * @details
		 *           Constructs an empty image of width 0 and height 0.
		 */
		BitImage();

		/**
		 * @brief
		 *           Constructs a blank image of specified dimension.
		 *
		 * @param height
		 *           The specified height of the image.
		 *
		 * @param width
		 *           The specified width of the image.
		 *
		 * @warning
		 *           If no sufficient memory could be allocated, the
		 *           method prints an error message to <code>stderr</code>
		 *           and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If <code>width</code> or <code>height</code> are smaller
		 *           than 0, the method prints an error message to
		 *           <code>stderr</code> and exits width status 'EXIT_FAILURE'.
		 *
		 */
		BitImage( int height , int width );

		/**
		 * @brief
		 *           Copy constructor.
		 *
		 * @param im
		 *           The image of which the copy is created.
		 *
		 * @details
		 *           Constructs an empty being a copy of <code>im</code>.
		 *
		 * @warning
		 *           If not sufficient memory could be allocated, an error
		 *           message is printed to <code>stderr</code> and the program
		 *           exits with status 'EXIT_FAILURE'.
		 */
		BitImage( const BitImage & im );

		/**
		 * @brief
		 *           Destructor.
		 *
		 * @details
		 *           Frees all the memory allocated by the binary image.
		 */
		~BitImage();

		/**
		 * @brief
		 *           Assignment operator.
		 *
		 * @details
		 *           The image will be set to a copy of <code>im</code>.
		 *
		 * @param im
		 *           The image of which this image becomes a copy of.
		 *
		 * @return
		 *           A reference to this image.
		 *
		 * @warning
		 *           If not sufficient memory could be allocated, an error
		 *           message is printed to <code>stderr</code> and the program
		 *           exits with status 'EXIT_FAILURE'.
		 */
		BitImage &operator=( const BitImage & im);

		/**
		 * @brief
		 *           Access the height of the image.
		 *
		 * @return
		 *           The height of the image.
		 */
		int getHeight() const;

		/**
		 * @brief
		 *           Access the width of the image.
		 *
		 * @return
		 *           The width of the image.
		 */
		int getWidth() const;

		/**
		 * @brief
		 *           Access the value of the pixel at the specified
		 *           coordinate.
		 *
		 * @param y
		 *           Row index of the accessed pixel.
		 *
		 * @param x
		 *           Column index of the accessed pixel.
		 *
		 * @return
		 *           <code>true</code> if the pixel at <code>(y,x)</code>
		 *           is black and <code>false</code> otherwise.
		 *
		 * @warning
		 *           If <code>y</code> or <code>x</code> are smaller than 0
		 *           or if <code>y</code> is greater than or equals the
		 *           image's height or if <code>x</code> is greater than or
		 *           equals the image's width the function prints an error
		 *           message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		bool at( int y , int x ) const;

		/**
		 * @brief
		 *           Sets the pixel at the specified coordinate to black.
		 *
		 * @param y
		 *           Row index of the accessed pixel.
		 *
		 * @param x
		 *           Column index of the accessed pixel.
		 *
		 * @warning
		 *           If <code>y</code> or <code>x</code> are smaller than 0
		 *           or if <code>y</code> is greater than or equals the
		 *           image's height or if <code>x</code> is greater than or
		 *           equals the image's width the function prints an error
		 *           message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		void setAt( int y , int x );

		/**
		 * @brief
		 *           Sets the pixel at the specified coordinate to white.
		 *
		 * @param y
		 *           Row index of the accessed pixel.
		 *
		 * @param x
		 *           Column index of the accessed pixel.
		 *
		 * @warning
		 *           If <code>y</code> or <code>x</code> are smaller than 0
		 *           or if <code>y</code> is greater than or equals the
		 *           image's height or if <code>x</code> is greater than or
		 *           equals the image's width the function prints an error
		 *           message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		void clearAt( int y , int x );

		/**
		 * @brief
		 *           Sets the pixel at the specified coordinate to a specified
		 *           value (black or white)
		 *
		 * @param v
		 *           If <code>true</code> the pixel at <code>(y,x)</code> will
		 *           be set to black; otherwise, if <code>false</code>, the
		 *           pixel will be white.
		 *
		 * @param y
		 *           Row index of the accessed pixel.
		 *
		 * @param x
		 *           Column index of the accessed pixel.
		 *
		 * @warning
		 *           If <code>y</code> or <code>x</code> are smaller than 0
		 *           or if <code>y</code> is greater than or equals the
		 *           image's height or if <code>x</code> is greater than or
		 *           equals the image's width the function prints an error
		 *           message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		void setAt( bool v , int y , int x );

		/**
		 * @brief
		 *           Stores the black-white information of the image in the
		 *           specified array.
		 *
		 * @details
		 *           <code>array</code> will hold the black-white (true-false)
		 *           information of the image in row major order. More
		 *           precisely, the value of the pixel at location
		 *           <code>at(y,x)</code> will be stored in
		 *           <code>array[y*getWidth()+x]</code>.
		 *
		 * @param array
		 *           The array to where the black-white information of the
		 *           image is stored.
		 *
		 * @warning
		 *           If <code>array</code> was not allocated to hold at least
		 *           <code>getHeight()*getWidth()</code> boolean value, calling
		 *           this method runs into undocumented behavior.
		 */
		void toArray( bool *array ) const;

		/**
		 * @brief
		 *           (Re-)Initializes the image by a blank image (white) of
		 *           the specified size.
		 *
		 * @param height
		 *           The specified height of the image.
		 *
		 * @param width
		 *           The specified width of the image.
		 *
		 * @warning
		 *           If no sufficient memory could be allocated, the
		 *           method prints an error message to <code>stderr</code>
		 *           and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If <code>width</code> or <code>height</code> are smaller
		 *           than 0, the method prints an error message to
		 *           <code>stderr</code> and exits width status 'EXIT_FAILURE'.
		 */
		void setSize( int height , int width );

		/**
		 * @brief
		 *           Initializes this image by the portable bitmap stored in
		 *           the specified file.
		 *
		 * @param fileName
		 *           The name of the file from where the bitmap is read.
		 *
		 * @return
		 *           <code>true</code> if the file could successfully be
		 *           read and <code>false</code> otherwise; if
		 *           <code>false</code> the image will be left unchanged.
		 *
		 * @warning
		 *           If not enough memory could be allocated, the function
		 *           prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 */
		bool read( const std::string & fileName );

		/**
		 * @brief
		 *           Initializes this image by the portable bitmap stored
		 *           in the specified input stream.
		 *
		 * @param in
		 *           The input stream from where to read the bitmap data.
		 *
		 * @return
		 *           <code>true</code> if the data could successfully be
		 *           read and <code>false</code> otherwise; if
		 *           <code>false</code> the image will be left unchanged.
		 *
		 * @warning
		 *           If not enough memory could be allocated, the function
		 *           prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 */
		bool read( FILE *in );

		/**
		 * @brief
		 *           Writes the image to the specified file in
		 *           PBM format.
		 *
		 * @param fileName
		 *           The name of the specified file to where the image
		 *           is written.
		 *
		 * @param ascii
		 *           If <code>true</code> the image will be written encoded
		 *           as ASCII; otherwise, if <code>false</code> the image
		 *           is encoded as a binary.
		 *
		 * @return
		 *           <code>true</code> if the image could be successfully
		 *           written; <code>false</code> if an error occurred on
		 *           attempting to write the image.
		 */
		bool write( const std::string & fileName , bool ascii = false ) const;

		/**
		 * @brief
		 *           Writes the image to the specified output stream in
		 *           PBM format.
		 *
		 * @param out
		 *           The output stream to where the image is written.
		 *
		 * @param ascii
		 *           If <code>true</code> the image will be written encoded
		 *           as ASCII; otherwise, if <code>false</code> the image
		 *           is encoded as binary data.
		 *
		 * @return
		 *           <code>true</code> if the image could be successfully
		 *           written; <code>false</code> if an error occurred on
		 *           attempting to write the image.
		 */
		bool write( FILE *out , bool ascii = false ) const;

		/**
		 * @brief
		 *           Swaps the content of two images such the
		 *           the one represents the other.
         *
         * @param im1
         *           Binary image being assigned with <code>im2</code>.
         *
         * @param im2
         *           Binary image being assigned with <code>im1</code>.
		 */
		static void swap( BitImage & im1 , BitImage & im2 );

	private:

		/**
		 * @brief
		 *            The pixel data of the bitmap.
		 */
		uint8_t *data;

		/**
		 * @brief
		 *            The height of the bitmap.
		 */
		int height;

		/**
		 * @brief
		 *            The width of the bitmap.
		 */
		int width;

		/**
		 * @brief
		 *            Initializes this image by the bitmap data (ASCII) stored
		 *            in the specified input stream.
		 *
		 * @param in
		 *            The specified input stream.
		 *
		 * @return
		 *           <code>true</code> if the data could successfully be
		 *           read and <code>false</code> otherwise; if
		 *           <code>false</code> the image will be left unchanged.
		 *
		 * @warning
		 *           If not enough memory could be allocated, the function
		 *           prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 */
		bool readASCII( FILE *in );

		/**
		 * @brief
		 *            Initializes this image by the bitmap data (binary) stored
		 *            in the specified input stream.
		 *
		 * @param in
		 *            The specified input stream.
		 *
		 * @return
		 *           <code>true</code> if the data could successfully be
		 *           read and <code>false</code> otherwise; if
		 *           <code>false</code> the image will be left unchanged.
		 *
		 * @warning
		 *           If not enough memory could be allocated, the function
		 *           prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 */
		bool readBinary( FILE *in );

		/**
		 * @brief
		 *           Writes the image to the specified output stream
		 *           encoded as ASCII data.
		 *
		 * @param out
		 *           The specified output stream.
		 *
		 * @return
		 *           <code>true</code> if the image could be successfully
		 *           written; <code>false</code> if an error occurred on
		 *           attempting to write the image.
		 */
		bool writeASCII( FILE *out ) const;

		/**
		 * @brief
		 *           Writes the image to the specified output stream
		 *           encoded as binary data.
		 *
		 * @param out
		 *           The specified output stream.
		 *
		 * @return
		 *           <code>true</code> if the image could be successfully
		 *           written; <code>false</code> if an error occurred on
		 *           attempting to write the image.
		 */
		bool writeBinary( FILE *out ) const;
	};
}


#endif /* THIMBLE_BITIMAGE_H_ */
