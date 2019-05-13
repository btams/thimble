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
 * @file BitImage.cpp
 *
 * @brief
 *            Implements the functionalities from 'BitImage.h' which
 *            provides a class for reading and writing binary images
 *            in PBM format.
 *
 * @author Benjamin Tams
 */

#include "config.h"
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/image/BitImage.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Standard constructor,
	 *
	 * @details
	 *           Constructs an empty image of width 0 and height 0.
	 */
	BitImage::BitImage() {
		this->data   = NULL;
		this->width  = 0;
		this->height = 0;
	}

	/**
	 * @brief
	 *           Constructs a blank image of specified dimension.
	 *
	 * @details
	 *           see 'BitImage.h'
	 *
	 */
	BitImage::BitImage( int height , int width ) {
		this->data   = NULL;
		this->height = 0;
		this->width  = 0;
		setSize(height,width);
	}

	/**
	 * @brief
	 *           Copy constructor.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	BitImage::BitImage( const BitImage & im ) {
		this->data   = NULL;
		this->width  = 0;
		this->height = 0;

		*this = im;
	}

	/**
	 * @brief
	 *           Destructor.
	 *
	 * @details
	 *           Frees all the memory allocated by the binary image.
	 */
	BitImage::~BitImage() {
		free(this->data);
	}

	/**
	 * @brief
	 *           Assignment operator.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	BitImage &BitImage::operator=( const BitImage & im ) {

		if ( this == &im ) {
			return *this;
		}

		int height , width , size;
		height = im.getHeight();
		width = im.getWidth();
		size = height*width;
		if ( size % 8 ) {
			size /= 8;
			++size;
		} else {
			++size;
		}

		setSize(height,width);

		memcpy(this->data,im.data,size);

		return *this;
	}

	/**
	 * @brief
	 *           Access the height of the image.
	 *
	 * @return
	 *           The height of the image.
	 */
	int BitImage::getHeight() const {
		return this->height;
	}

	/**
	 * @brief
	 *           Access the width of the image.
	 *
	 * @return
	 *           The width of the image.
	 */
	int BitImage::getWidth() const {
		return this->width;
	}

	/**
	 * @brief
	 *           Access the value of the pixel at the specified
	 *           coordinate.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::at( int y , int x ) const {

		if ( y < 0 || x < 0 || y >= this->height || x >= this->width ) {
			cerr << "BitImage::at: Accessed pixel is out of the image's "
				 << "region." << endl;
			exit(EXIT_FAILURE);
		}

		int offset = this->width/8;
		if ( this->width % 8 ) {
			++offset;
		}

		int i , j;
		i = x / 8;
		j = x % 8;

		return (this->data[y*offset+i] & (1<<(7-j)) )?true:false;
	}

	/**
	 * @brief
	 *           Sets the pixel at the specified coordinate to black.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	void BitImage::setAt( int y , int x ) {

		if ( y < 0 || x < 0 || y >= this->height || x >= this->width ) {
			cerr << "BitImage::setAt: Set pixel is out of the image's "
				 << "region." << endl;
			exit(EXIT_FAILURE);
		}

		int offset = this->width/8;
		if ( this->width % 8 ) {
			++offset;
		}

		int i , j;
		i = x / 8;
		j = x % 8;

		this->data[y*offset+i] |= (1<<(7-j));
	}

	/**
	 * @brief
	 *           Sets the pixel at the specified coordinate to white.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	void BitImage::clearAt( int y , int x ) {

		if ( y < 0 || x < 0 || y >= this->height || x >= this->width ) {
			cerr << "BitImage::at: Accessed pixel is out of the image's "
				 << "region." << endl;
			exit(EXIT_FAILURE);
		}

		int offset = this->width/8;
		if ( this->width % 8 ) {
			++offset;
		}

		int i , j;
		i = x / 8;
		j = x % 8;

		this->data[y*offset+i] &= ~(1<<(7-j));

		/*int i , j , k;

		k = y * this->width + x;
		i = k / 8;
		j = k % 8;

		this->data[i] &= ~(1<<(7-j));*/
	}

	/**
	 * @brief
	 *           Sets the pixel at the specified coordinate to a specified
	 *           value (black or white)
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	void BitImage::setAt( bool v , int y , int x ) {

		if ( v ) {
			setAt(y,x);
		} else {
			clearAt(y,x);
		}
	}

	/**
	 * @brief
	 *           Stores the black-white information of the image in the
	 *           specified array.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	void BitImage::toArray( bool *array ) const {

		int m , n;
		m = getHeight();
		n = getWidth();

		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				array[y*n+x] = at(y,x);
			}
		}
	}

	/**
	 * @brief
	 *           (Re-)Initializes the image by a blank image (white) of
	 *           the specified size.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	void BitImage::setSize( int height , int width ) {

		if ( height < 0 || width < 0 ) {
			cerr << "BitImage::setSize: Bad arguments; height and width must"
				 << " be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		free(this->data);

		if ( width == 0 || height == 0 ) {
			this->data = NULL;
		} else {

			int size = width;
			if ( size % 8 ) {
				size /= 8;
				++size;
			} else {
				size /= 8;
			}
			size *= height;

			this->data = (uint8_t*)malloc( size * sizeof(uint8_t) );
			if ( this->data == NULL ) {
				cerr << "BitImage::setSize: Out of memory." << endl;
				exit(EXIT_FAILURE);
			}
			memset(this->data,0,size*sizeof(uint8_t));
		}

		this->height = height;
		this->width  = width;
	}

	/**
	 * @brief
	 *            Initializes this image by the portable bitmap stored
	 *            in the specified input stream.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::read( const std::string & fileName ) {

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
	 *           Initializes this image by the portable bitmap stored
	 *           in the specified input stream.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::read( FILE *in ) {

		if ( fgetc(in) != 'P' ) {
			return false;
		}

		int c;

		c = fgetc(in);
		if ( c == '1') {
			return readASCII(in);
		} else if ( c == '4' ) {
			return readBinary(in);
		} else {
			return false;
		}

		return true;
	}

	/**
	 * @brief
	 *           Writes the image to the specified file in
	 *           PBM format.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::write( const std::string & fileName , bool ascii ) const {

		FILE *out = THIMBLE_FOPEN(fileName.c_str(),"wb");
		if ( out == NULL ) {
			return false;
		}

		bool state = write(out,ascii);

		fclose(out);

		return state;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream in
	 *           PBM format.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::write( FILE *out , bool ascii ) const {
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
	 *            The function ignores characters that are no digits;
	 *            furthermore, the subsequent characters of the line are
	 *            ignored that follow the character '#'.
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
	 *            Returns the next decimal digit
	 *            from the specified input stream.
	 *
	 * @details
	 *            The function ignores characters that are no digits;
	 *            furthermore, the subsequent characters of the line are
	 *            ignored that follow the character '#'.
	 *
	 * @param in
	 *            The specified input stream.
	 *
	 * @return
	 *            The next decimal digit stored in <code>in</code>
	 *            in decimal representation; if no digit could be read
	 *            the function returns -1.
	 */
	static int nextdigit( FILE *in ) {

		int c;

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

		return c-'0';
	}

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
	void BitImage::swap( BitImage & im1 , BitImage & im2 ) {

		uint8_t *data;
		int width , height;

		data = im1.data;
		width = im1.width;
		height = im2.height;

		im1.data = im2.data;
		im1.width = im2.width;
		im1.height = im2.height;

		im2.data = data;
		im2.width = width;
		im2.height = height;
	}

	/**
	 * @brief
	 *            Initializes this image by the bitmap data (ASCII) stored
	 *            in the specified input stream.
	 *
	 * @details
	 *            see 'BitImage.h'
	 */
	bool BitImage::readASCII( FILE *in ) {

		// Read the width and the height
		int width , height;
		width = nextvalue(in);
		height = nextvalue(in);
		if ( width < 0 || height < 0 ) {
			return false;
		}

		// How many bytes 'size'?
		int size = width/8;
		if ( width%8 ) {
			++size;
		}
		size *= height;

		// If the width or height are 0 then size is zero and
		// the image will be initialized correspondingly.
		if ( size == 0 ) {
			free(this->data);
			this->data = NULL;
			this->width = width;
			this->height = height;
			return true;
		}

		BitImage tmp(height,width);

		for ( int y = 0 ; y < height ; y++ ) {
			for ( int x = 0 ; x < width ; x++ ) {

				int p = nextdigit(in);

				if ( p == 1 ) {
					tmp.setAt(y,x);
				} else if ( p != 0 ) {
					return false;
				}
			}
		}

		BitImage::swap(*this,tmp);

		return true;
	}

	/**
	 * @brief
	 *            Initializes this image by the bitmap data (binary) stored
	 *            in the specified input stream.
	 *
	 * @details
	 *            see 'BitImage.h'
	 */
	bool BitImage::readBinary( FILE *in ) {

		// Read the width and the height
		int width , height;
		width  = nextvalue(in);
		height = nextvalue(in);
		if ( width < 0 || height < 0 ) {
			return false;
		}

		// How many bytes 'size'?
		int size = width;
		if ( size % 8 ) {
			size /= 8;
			++size;
		} else {
			size /= 8;
		}
		size *= height;

		if ( size == 0 ) {
			free(this->data);
			this->data = NULL;
			this->width = width;
			this->height = height;
			return true;
		}

		uint8_t *data = (uint8_t*)malloc( size * sizeof(uint8_t) );
		if ( data == NULL ) {
			cerr << "BitImage::read: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		if ( (int)fread(data,1,size,in) != size ) {
			free(data);
			return false;
		}

		free(this->data);
		this->data   = data;
		this->width  = width;
		this->height = height;

		return true;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream
	 *           encoded as ASCII data.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::writeASCII( FILE *out ) const {

		if ( fprintf
				(out,"P1\n"
					 "# Created using THIMBLE\n"
					 "%i %i\n",
				 this->width,this->height ) < 0 ) {
			return false;
		}

		int c = 0;
		int n = std::min(this->width,70);

		for ( int y = 0 ; y < this->height ; y++ ) {
			for ( int x = 0 ; x < this->width ; x++ ) {

				if ( at(y,x) ) {
					if ( fputc('1',out) < 0 ) {
						return false;
					}
				} else {
					if ( fputc('0',out) < 0 ) {
						return false;
					}
				}

				++c;
				c %= n;
				if ( c == 0 ) {
					if ( fputc('\n',out) < 0 ) {
						return false;
					}
				}
			}
		}

		fflush(out);

		return true;
	}

	/**
	 * @brief
	 *           Writes the image to the specified output stream
	 *           encoded as binary data.
	 *
	 * @details
	 *           see 'BitImage.h'
	 */
	bool BitImage::writeBinary( FILE *out ) const {

		if ( fprintf
				(out,"P4\n"
					 "# Created using THIMBLE\n"
					 "%i %i\n",
				 this->width,this->height ) < 0 ) {
			return false;
		}

		int size = this->width;
		if ( size % 8 ) {
			size /= 8;
			++size;
		} else {
			size /= 8;
		}
		size *= this->height;

		if ( fwrite(this->data,1,size,out) < 0 ) {
			return false;
		}

		return true;
	}
}




