/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013, 2014 Benjamin Tams
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
 * @file GrayImage.h
 *
 * @brief
 *            Provides a class for reading and writing gray-scale images in
 *            PGM format.
 *
 * @details
 * @section sec_grayimages Gray-Scale Images
 *
 * The <a href="http://en.wikipedia.org/wiki/Portable_graymap">PGM format</a>
 * (portable graymap) is an easily usable format for exchanging gray-scale
 * images between different programs and THIMBLE provides a class for
 * representing, reading and writing gray-scale images in PGM format.
 *
 * A gray-scale image is represented by a
 * \link thimble::GrayImage GrayImage\endlink object
 * <pre>
 *  GrayImage image;
 * </pre>
 * The width and height of the image can be accessed via
 * <pre>
 *  int height , width;
 *
 *  height = image.getHeight();
 *  width = image.getWidth();
 * </pre>
 * The intensity of the pixel at the <code>i</code>th row
 * (<code>i=0,...,height-1</code>) and <code>j</code>th column
 * (<code>j=0,...,width-1</code>) can be accessed using the
 * \link thimble::GrayImage::at(int,int)const at()\endlink function,
 * for example by calling
 * <pre>
 *  int p = image.at(i,j);
 * </pre>
 * The value of <code>p</code> in the above example usually varies between
 * 0 (black) and 255 (white). However, the maximal number of pixel intensity
 * can be smaller and (sometimes) even larger than 255 (up to 65535). The bound
 * of pixel intensity associated with <code>image</code> can be accessed
 * with the \link thimble::GrayImage::getMaxValue() getMaxValue()\endlink
 * function, for example,
 * <pre>
 *  int maxvalue = image.getMaxValue();
 * </pre>
 * Consequently, to obtain a normalized representation of the pixel's intensity
 * between 0.0 (black) and 1.0 (white), we may run
 * <pre>
 *  double v = (double)image.at(i,j) / (double)image.getMaxValue();
 * </pre>
 *
 * @subsection sec_grayimages_create Construction and Modification
 *
 * There are multiple ways to create a
 * \link thimble::GrayImage GrayImage\endlink. To create an empty
 * \link thimble::GrayImage GrayImage\endlink <code>image</code> of height
 * <code>height</code> and width <code>width</code> of which all pixels are
 * set to 0 (black) we may run
 * <pre>
 *  GrayImage image(height,width);
 * </pre>
 * Alternatively, we may call the
 * \link thimble::GrayImage::setSize(int,int,int) setSize()\endlink method to
 * overwrite the data of <code>image</code> by an empty image of width
 * <code>width</code> and <code>height</code> of which all pixels will be set
 * to 0 (black)
 * <pre>
 *  GrayImage image;
 *
 *  image.setSize(width,height);
 * </pre>
 * The intensity of a pixel at <code>(i,j)</code> can be specified via
 * the \link thimble::GrayImage::setAt(int,int,int) setAt()\endlink method.
 * For example, if the pixel at row <code>i</code> and column <code>j</code>
 * is sought to be set to <code>p</code> we may run
 * <pre>
 *  image.setAt(p,i,j);
 * </pre>
 * Above <code>p</code> usually denotes an integer between 0 (black) and
 * 255 (white). However, in case the value
 * \link thimble::GrayImage::getMaxValue() getMaxValue\endlink
 * (which is 255 per default but can be specified manually using the
 * \link thimble::GrayImage::GrayImage(int,int,int)
 * GrayImage(int,int,int)\endlink or the
 * \link thimble::GrayImage::setSize(int,int,int)
 * setSize(int,int,int)\endlink method) the bound for the pixel
 * <code>p</code> is 0 and
 * \link thimble::GrayImage::getMaxValue() getMaxValue()\endlink.
 *
 *
 * @subsection sec_grayimages_io In/Out Interface
 *
 * The data of a \link thimble::GrayImage GrayImage\endlink can be read
 * from PGM data via the \link thimble::GrayImage::read() read()\endlink
 * function which returns <code>true</code> if an image has been successfully
 * read and <code>false</code> otherwise. For example, an image from a file
 * <code>in.pgm</code> that is in PGM format can be read via
 * <pre>
 *  if ( !image.read("in.pgm") ) {
 *     cout << "Could not read from 'in.pgm'." << endl;
 *  }
 * </pre>
 * Analogously, the image can be written to a file <code>out.pgm</code>
 * through the member function
 * \link thimble::GrayImage::write() write()\endlink:
 * <pre>
 *  if ( !image.write("out.pgm") ) {
 *     cout << "Could not write to 'out.pgm'." << endl;
 *  }
 * </pre>
 *
 * @see http://en.wikipedia.org/wiki/Portable_graymap
 * @see thimble::GrayImage
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_GRAYIMAGE_H_
#define THIMBLE_GRAYIMAGE_H_

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
	 *            Class representing a gray-scale image.
	 */
	class THIMBLE_DLL GrayImage {

	public:

		/**
		 * @brief
		 *            Standard constructor.
		 *
		 * @details
		 *            Constructs an empty image of width and height 0.
		 */
		GrayImage();

		/**
		 * @brief
		 *            Constructs a black image of specified dimension.
		 *
		 * @details
		 *            The intensities of every pixel is 0 (i.e. black).
		 *
		 * @param height
		 *            The specified height of the image.
		 *
		 * @param width
		 *            The specified width of the image.
		 *
		 * @param maxvalue
		 *            Specifies the maximal value a pixel's intensity
		 *            can have.
		 *
		 * @warning
		 *            If <code>height</code> is smaller than 0 or if
		 *            <code>width</code> is smaller than 0 or if
		 *            <code>maxvalue</code> is smaller than 0 or if
		 *            <code>maxvalue</code> is greater than 65535, an error
		 *            message is written to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory can be allocated, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		GrayImage( int height , int width , int maxvalue = 255 );

		/**
		 * @brief
		 *            Copy constructor.
	     *
		 * @details
		 *            Constructs an image that is a copy of <code>im</code>.
		 *
		 * @param im
		 *            The image of which the copy is created.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		GrayImage( const GrayImage & im );

		/**
		 * @brief
		 *            Destructor.
		 *
		 * @details
		 *            Frees all the memory allocated by the gray-scale image.
		 */
		~GrayImage();

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
		GrayImage & operator=( const GrayImage & im );

		/**
		 * @brief
		 *           Access the maximal intensity value a pixel can have.
		 *
		 * @return
		 *           The maximal intensity value a pixel can have.
		 */
		int getMaxValue() const;

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
		 *            Access the intensity of the specified pixel.
		 *
		 * @param y
		 *            Row index of the accessed pixel.
		 *
		 * @param x
		 *            Column index of the accessed pixel.
		 *
		 * @return
		 *            The intensity of the accessed pixel.
		 *
		 * @warning
		 *            If <code>y</code> or <code>x</code> are smaller than 0
		 *            or if <code>y</code> is greater than or equals the
		 *            image's height or if <code>x</code> is greater than or
		 *            equals the image's width the function prints an error
		 *            message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		int at( int y , int x ) const;

		/**
		 * @brief
		 *           Sets the pixel at the specified coordinate to a specified
		 *           intensity value
		 *
		 * @param v
		 *           The new intensity value of the accessed pixel.
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
		 *           If <code>v</code> is smaller than 0 or if <code>v</code>
		 *           is greater than <code>getMaximalValue()</code>, the
		 *           method prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 */
		void setAt( int v , int y , int x );

		/**
		 * @brief
		 *           Converts the intensity values of the image in an array
		 *           whose values vary between 0.0 and 1.0.
		 *
		 * @details
		 *           The intensities are converted to <code>array</code> in
		 *           row major order. More precisely, the intensity at
		 *           <code>at(y,x)</code> is stored at
		 *           <code>array[y*n+x]</code> where
		 *           <code>n=getWidth()</code>.
		 *
		 * @param array
		 *           The array to where the pixels intensities are converted.
         *
		 * @warning
		 *           If <code>array</code> was not allocated to hold at least
		 *           <code>getHeight()</code>*<code>getWidth()</code> double
		 *           values, calling this method runs into undocumented behavior.
		 */
		void toArray( double *array ) const;

		/**
		 * @brief
		 *           Converts the intensity values of the image in an array
		 *           whose values vary between 0 and 255.
		 *
		 * @details
		 *           The intensities are converted to <code>array</code> in
		 *           row major order. More precisely, the intensity at
		 *           <code>at(y,x)</code> is stored at
		 *           <code>array[y*n+x]</code> where
		 *           <code>n=getWidth()</code>.
		 *
		 * @param array
		 *           The array to where the pixels intensities are converted.
         *
         * @warning
		 *           If <code>array</code> was not allocated to hold at least
		 *           <code>getHeight()</code>*<code>getWidth()</code> double
		 *           values, calling this method runs into undocumented behavior.
		 */
		void toArray( uint8_t *array ) const;

		/**
		 * @brief
		 *            Initializes the image by a black image of specified
		 *            dimension.
		 *
		 * @details
		 *            The intensities of every pixel is 0 (i.e. black).
		 *
		 * @param height
		 *            The specified height of the image.
		 *
		 * @param width
		 *            The specified width of the image.
		 *
		 * @param maxvalue
		 *            Specifies the maximal value a pixel's intensity
		 *            can have.
		 *
		 * @warning
		 *            If <code>height</code> is smaller than 0 or if
		 *            <code>width</code> is smaller than 0 or if
		 *            <code>maxvalue</code> is smaller than 0 or if
		 *            <code>maxvalue</code> is greater than 65535, an error
		 *            message is written to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory can be allocated, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setSize( int height , int width , int maxvalue = 255 );

		/**
		 * @brief
		 *           Initializes this image by the portable gray stored in
		 *           the specified file.
		 *
		 * @param fileName
		 *           The name of the file from where the graymap is read.
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
		 *           Initializes this image by the portable graymap stored
		 *           in the specified input stream.
		 *
		 * @param in
		 *           The input stream from where to read the graymap data.
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
		 *           PGM format.
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
		 *           PGM format.
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

	private:

		/**
		 * @brief
		 *           The data of the graymap image.
		 */
		uint8_t *data;

		/**
		 * @brief
		 *           The maximal intensity value a pixel can have.
		 */
		int maxvalue;

		/**
		 * @brief
		 *            The height of the graymap.
		 */
		int width;

		/**
		 * @brief
		 *            The width of the graymap.
		 */
		int height;

		/**
		 * @brief
		 *            Initializes this image by the graymap data (ASCII) stored
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
		 *            Initializes this image by the graymap data (binary)
		 *            stored in the specified input stream.
		 *
		 * @param in
		 *            The input stream.
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
		bool writeBinary( FILE *out) const;
	};
}

#endif /* THIMBLE_GRAYIMAGE_H_ */
