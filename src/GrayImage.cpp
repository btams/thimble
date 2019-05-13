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
 * @file BitImage.cpp
 *
 * @brief
 *            Implements the functionalities from 'GrayImage.h' which
 *            provides a class for reading and writing gray-scale images
 *            in PGM format.
 *
 * @author Benjamin Tams
 */

#include "config.h"
#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <iostream>

#include <thimble/image/GrayImage.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Standard constructor.
	 *
	 * @details
	 *            Constructs an empty image of width and height 0.
	 */
	GrayImage::GrayImage() {
		this->data = NULL;
		this->maxvalue = 0;
		this->width = 0;
		this->height = 0;
	}

	/**
	 * @brief
	 *            Constructs a black image of specified dimension.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	GrayImage::GrayImage( int width , int height , int maxvalue ) {
		this->data = NULL;
		this->maxvalue = 0;
		this->width = 0;
		this->height = 0;
		setSize(width,height,maxvalue);
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	GrayImage::GrayImage( const GrayImage & im ) {
		this->data = NULL;
		this->maxvalue = 0;
		this->width = 0;
		this->height = 0;
		*this = im;
	}

	/**
	 * @brief
	 *            Destructor.
	 *
	 * @details
	 *            Frees all the memory allocated by the gray-scale image.
	 */
	GrayImage::~GrayImage() {
		free(this->data);
		this->data = NULL;
	}

	/**
	 * @brief
	 *           Assignment operator.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	GrayImage & GrayImage::operator=( const GrayImage & im ) {

		if ( this == &im ) {
			return *this;
		}

		setSize(im.getHeight(),im.getWidth(),im.maxvalue);

		int so = (maxvalue<256)?1:2;

		memcpy(this->data,im.data,width*height*so);

		return *this;
	}

	/**
	 * @brief
	 *           Access the maximal intensity value a pixel can have.
	 *
	 * @return
	 *           The maximal intensity value a pixel can have.
	 */
	int GrayImage::getMaxValue() const {
		return this->maxvalue;
	}

	/**
	 * @brief
	 *           Access the height of the image.
	 *
	 * @return
	 *           The height of the image.
	 */
	int GrayImage::getHeight() const {
		return this->height;
	}

	/**
	 * @brief
	 *           Access the width of the image.
	 *
	 * @return
	 *           The width of the image.
	 */
	int GrayImage::getWidth() const {
		return this->width;
	}

	/**
	 * @brief
	 *            Access the intensity of the specified pixel.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	int GrayImage::at( int y , int x ) const {

		int k = y*width+x;

		if ( this->maxvalue<256 )
			return (int)this->data[k];

		return (((int)data[2*k])<<8)+(int)data[2*k+1];
	}

	/**
	 * @brief
	 *           Sets the pixel at the specified coordinate to a specified
	 *           intensity value
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	void GrayImage::setAt( int v , int y , int x ) {

		if ( y < 0 || x < 0 || y >= this->height || x >= this->width ) {
			cerr << "BitImage::at: Accessed pixel is out of the image's "
				 << "region." << endl;
			exit(EXIT_FAILURE);
		}

		if ( v < 0 || v > maxvalue ) {
			cerr << "GrayImage::setAt: Gray-value is out of range." << endl;
			exit(EXIT_FAILURE);
		}

		int k = y*width+x;

		if ( this->maxvalue<256 ) {
			this->data[k] = (uint8_t)v;
		} else {
			data[2*k] = (uint8_t)(v >> 8);
			data[2*k+1] = (uint8_t)(v & 0xFF);
		}
	}

	/**
	 * @brief
	 *           Converts the intensity values of the image in an array
	 *           whose values vary between 0.0 and 1.0.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	void GrayImage::toArray( double *array ) const {

		int maxValue = getMaxValue();

		for ( int y = 0 ; y < height ; y++ ) {
			for ( int x = 0 ; x < width ; x++ ) {
				array[y*width+x] = (double)at(y,x) / (double)maxValue;
			}
		}
	}

	/**
	 * @brief
	 *           Converts the intensity values of the image in an array
	 *           whose values vary between 0 and 255.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	void GrayImage::toArray( uint8_t *array ) const {

		int maxValue = getMaxValue();

		for ( int y = 0 ; y < height ; y++ ) {
			for ( int x = 0 ; x < width ; x++ ) {
				array[y*width+x] = (uint8_t)THIMBLE_ROUND
					(255.0 * (double)at(y,x) / (double)maxValue);
			}
		}
	}

	/**
	 * @brief
	 *            Initializes the image by a black image of specified
	 *            dimension.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	void GrayImage::setSize( int height , int width , int maxvalue ) {

		free(this->data);

		if ( width < 0 || height < 0 ) {
			cerr << "GrayImage::setSize: Bad arguments: 'height' and 'width' "
					"must be greater than or equals 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( maxvalue <= 0 || maxvalue >= 0x10000 ) {
			cerr << "GrayImage::setSize: Bad argument: "
				 << "'maxvalue' must be between 0 and 65535" << endl;
			exit(EXIT_FAILURE);
		}

		int so = (maxvalue<256)?1:2;
		this->data = (uint8_t*)malloc( width * height * so );

		this->width = width;
		this->height = height;
		this->maxvalue = maxvalue;

		memset(this->data,0,width*height*so);
	}

	/**
	 * @brief
	 *           Initializes this image by the portable gray stored in
	 *           the specified file.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::read( const std::string & fileName ) {

		FILE *in;
		if ( (in = THIMBLE_FOPEN(fileName.c_str(),"rb")) == NULL ) {
			return false;
		}

		bool state = this->read(in);

		fclose(in);

		return state;
	}

	/**
	 * @brief
	 *           Initializes this image by the portable graymap stored
	 *           in the specified input stream.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::read( FILE *in ) {

		if ( fgetc(in) != 'P' ) {
			return false;
		}

		int c;

		c = fgetc(in);

		if ( c != '2' && c != '5' ) {
			return false;
		}

		return (c=='2') ?(readASCII(in)): (readBinary(in));
	}

	/**
	 * @brief
	 *           Writes the image to the specified file in
	 *           PGM format.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::write( const std::string & fileName , bool ascii ) const {

		FILE *out;
		if ( (out=THIMBLE_FOPEN(fileName.c_str(),"wb")) == NULL ) {
			return false;
		}

		bool state = this->write(out,ascii);

		fclose(out);

		return state;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream in
	 *           PGM format.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::write( FILE *out , bool ascii ) const {

		if ( ascii ) {
			return writeASCII(out);
		} else {
			return writeBinary(out);
		}
	}

	/**
	 * @brief
	 *            Returns the next positive integer in decimal representation
	 *            from the specified input stream.
	 *
	 * @details
	 *            The function characters that are no digits; furthermore,
	 *            the subsequent characters of the line are ignored that
	 *            follow the character '#'.
	 *
	 * @param in
	 *            The specified input stream.
	 *
	 * @return
	 *            The next positive integer which is stored in <code>in</code>
	 *            in decimal representation; if the next token is not a valid
	 *            positive integer, the function returns -1.
	 */
	static int nextvalue( FILE *in ) {

		int c , n;

		while ( !isdigit((c=fgetc(in))) ) {
			if ( c < 0 ) {
				return -1;
			}
			if ( c == '#' ) {
				while ( (c=fgetc(in)) != '\n' ) {
					if ( c < 0 ) {
						return -1;
					}
				}
			}
		}

		n = (int)(c-'0');
		while ( isdigit(c=fgetc(in)) ) {
			n *= 10;
			n += c-'0';
		}

		return n;
	}

	/**
	 * @brief
	 *            Initializes this image by the graymap data (ASCII) stored
	 *            in the specified input stream.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	bool GrayImage::readASCII( FILE *in ) {

		int n , m , maxvalue;
		n = nextvalue(in);
		m = nextvalue(in);
		maxvalue = nextvalue(in);

		if ( m < 0 || n < 0 || maxvalue < 0 || maxvalue >= 0x10000 ) {
			return false;
		}

		GrayImage tmp(m,n,maxvalue);

		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				int v = nextvalue(in);
				if ( v < 0 ) {
					return false;
				}
				tmp.setAt(v,y,x);
			}
		}

		*this = tmp;

		return true;
	}

	/**
	 * @brief
	 *            Initializes this image by the graymap data (binary)
	 *            stored in the specified input stream.
	 *
	 * @details
	 *            see 'GrayImage.h'
	 */
	bool GrayImage::readBinary( FILE *in ) {

		int m , n , maxvalue;
		n = nextvalue(in);
		m = nextvalue(in);
		maxvalue = nextvalue(in);

		if ( m < 0 || n < 0 || maxvalue < 0 || maxvalue >= 0x10000 ) {
			return false;
		}

		GrayImage tmp(m,n,maxvalue);

		int so = (maxvalue<256)?1:2;

		size_t size = (size_t)(so*m*n);
		if ( fread(tmp.data,1,size,in) != size ) {
			return false;
		}

		*this = tmp;

		return true;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream
	 *           encoded as ASCII data.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::writeASCII( FILE *out ) const {

		if ( fprintf
				(out,"P2\n"
					 "# Created using THIMBLE\n"
					 "%i %i\n",
				 this->width,this->height) < 0 ) {
			return false;
		}

		int maxValueLength = fprintf(out,"%i",this->maxvalue);
		if ( maxValueLength < 0 ) {
			return false;
		}

		if ( fprintf(out,"\n") < 0 ) {
			return false;
		}

		int n = 0;
		for ( int y = 0 ; y < getHeight() ; y++ ) {
			for ( int x = 0 ; x < getWidth() ; x++ ) {
				int r = fprintf(out,"%i",at(y,x));
				if ( r < 0 ) {
					return false;
				}
				n += r;
				if ( n+maxValueLength+1 > 70 ) {
					if ( fputc('\n',out) < 0 ) {
						return false;
					}
					n = 0;
				} else {
					if ( fputc(' ',out) < 0 ) {
						return false;
					}
					++n;
				}
			}
		}

		return true;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream
	 *           encoded as binary data.
	 *
	 * @details
	 *           see 'GrayImage.h'
	 */
	bool GrayImage::writeBinary( FILE *out ) const {

		int so = (this->maxvalue<256)?1:2;

		if ( fprintf
				(out,"P5\n"
					 "# Created using THIMBLE\n"
					 "%i %i\n%i\n",
				 this->width,this->height,
				 this->maxvalue) < 0 ) {
			return false;
		}

		if ( fwrite(this->data,1,so*this->width*this->height,out) < 0 ) {
			return false;
		}

		return true;
	}
}
