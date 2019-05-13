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
 * @file Fingerprint.cpp
 *
 * @brief
 *            Implements mechanisms for estimating certain fingerprint
 *            features as provided by the 'Fingerprint.h' header file.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

#include "config.h"
#include <thimble/image/Orientation.h>
#include <thimble/image/GrayImage.h>
#include <thimble/finger/Segmentation.h>
#include <thimble/finger/TentedArchModel.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/finger/Fingerprint.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Initializes the values of this fingerprint with default
	 *            values.
	 *
	 * @details
	 *            The method is useful to create an empty fingerprint on
	 *            first initialzation where no data must be freed.
	 *            To clear an already initialized fingerprint the method
	 *            \link clear()\endlink should be used.
	 */
	void Fingerprint::first_init() {
		this->dpi = 0;
		this->m = 0;
		this->n = 0;
		this->intensityImage = NULL;
		this->orientationImage = NULL;
		this->foregroundImage = NULL;
		this->drpPtr = NULL;
		this->recordPtr = NULL;
	}

	/**
	 * @brief
	 *            Standard constructor.
	 *
	 * @details
	 *            Creates an empty fingerprint.
	 */
	Fingerprint::Fingerprint() {
		first_init();
	}

	/**
	 * @brief
	 *             Creates a fingerprint specified by its raw gray-scale
	 *             image.
	 *
	 * @param grayImage
	 *             The raw fingerprint image.
	 *
	 * @param dpi
	 *             The resolution in <em>dots per inch</em> at which the
	 *             fingerprint image has been scanned.
	 *
	 * @warning
	 *             If <code>dpi</code> is not greater than 0, an error message
	 *             is printed to <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	Fingerprint::Fingerprint
		( const GrayImage & grayImage , int dpi ) {

		first_init();

		initialize(grayImage,dpi);
	}

	/**
	 * @brief
	 *            Creates a fingerprint specified by its raw gray-scale
	 *            intensity image.
	 *
	 * @param intensityImage
	 *            An array of <i>m*n</i> double values ranging between
	 *            0.0 (black) and 1.0 (white) where the intensity
	 *            at <i>(x,y)</i> is stored at
	 *            <code>intensityImage[y*n+x]</code>.
	 *
	 * @param m
	 *            Height of the fingerprint image.
	 *
	 * @param n
	 *            Width of the fingerprint image.
	 *
	 * @param dpi
	 *            The resolution in <em>dots per inch</em> at which the
	 *            fingerprint image has been scanned.
	 *
	 * @warning
	 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
	 *            are not greater than 0, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	Fingerprint::Fingerprint
		( const double *intensityImage , int m , int n , int dpi ) {

		first_init();

		initialize(intensityImage,m,n,dpi);
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @param fingerprint
	 *            The fingerprint of which a copy is created.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	Fingerprint::Fingerprint
		( const Fingerprint & fingerprint ) {

		first_init();

		assign(fingerprint);
	}

	/**
	 * @brief
	 *           Destructor.
	 *
	 * @details
	 *           Calls \link clear()\endlink which, as a side-effect,
	 *           frees any memory held by this fingerprint.
	 */
	Fingerprint::~Fingerprint() {
		clear();
	}

	/**
	 * @brief
	 *           Clears the data for this fingerprint such that it
	 *           becomes empty.
	 *
	 * @details
	 *           After this method has been called, the function
	 *           \link isEmpty()\endlink returns <code>true</code>.
	 *           More specifically, the result of the functions
	 *           <ul>
	 *            <li>getWidth()</li>
	 *            <li>getHeight()</li>
	 *            <li>getResolution()</li>
	 *           </ul>
	 *           will be 0 right after the fingerprint has been cleared.
	 */
	void Fingerprint::clear() {

		free(this->intensityImage);
		free(this->foregroundImage);

		// For each pixel, free the data consumed by its orientation
		// estimation (if any)
		for ( int y = 0 ; y < this->m ; y++ ) {
			for ( int x = 0 ; x < this->n ; x++ ) {
				delete this->orientationImage[y*this->n+x];
				this->orientationImage[y*this->n+x] = NULL;
			}
		}
		delete[] this->orientationImage;
		delete this->drpPtr;
		delete this->recordPtr;

		// Fill with default value such that this fingerprint
		// becomes empty
		first_init();
	}

	/**
	 * @brief
	 *           Check whether this fingerprint is empty.
	 *
	 * @return
	 *           <code>true</code> if the fingerprint is not initialized
	 *           by a gray-scale image; otherwise, the function returns
	 *           <code>false</code>.
	 */
	bool Fingerprint::isEmpty() const {
		return this->m == 0 && this->n == 0;
	}

	/**
	 * @brief
	 *             Initialize this instance with a new fingerprint
	 *             specified by its raw gray-scale image.
	 *
	 * @param grayImage
	 *             The raw fingerprint image.
	 *
	 * @param dpi
	 *             The resolution in <em>dots per inch</em> at which the
	 *             fingerprint image has been scanned.
	 *
	 * @warning
	 *             If <code>dpi</code> is not greater than 0, an error message
	 *             is printed to <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::initialize( const GrayImage & grayImage , int dpi ) {

		if ( dpi <= 0 ) {
			cerr << "Fingerprint::initialize: "
				 << "resolution must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		// Make the fingerprint empty
		clear();

		int m , n;
		m = grayImage.getHeight();
		n = grayImage.getWidth();

		// Create memory to hold the fingerprint's intensity image ...
		double *intensityImage = (double*)malloc( m * n * sizeof(double) );
		if ( intensityImage == NULL ) {
			cerr << "Fingerprint::initialize: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// ... and get the gray-scale values.
		grayImage.toArray(intensityImage);

		// Fill in the fields needed to represent an initialized
		// fingerprint.
		this->dpi = dpi;
		this->m = m;
		this->n = n;
		this->intensityImage = intensityImage;

		// The orientation image should be an array of NULL
		// pointers on initialization. Therefore, allocate the an
		// array of pointers to an 'Orientation' and ...
		this->orientationImage = new(nothrow) Orientation*[m*n];
		if ( this->orientationImage == NULL ) {
			cerr << "Fingerprint:: initialize: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// ... initialize its elements with 'NULL'
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				this->orientationImage[y*n+x] = NULL;
			}
		}
	}

	/**
	 * @brief
	 *             Initialize this instance with a new fingerprint
	 *             specified by its double-valued gray-scale intensity
	 *             image.
	 *
	 * @param intensityImage
	 *            An array of <i>m*n</i> double values ranging between
	 *            0.0 (black) and 1.0 (white) where the intensity
	 *            at <i>(x,y)</i> is stored at
	 *            <code>intensityImage[y*n+x]</code>.
	 *
	 * @param m
	 *            Height of the fingerprint image.
	 *
	 * @param n
	 *            Width of the fingerprint image.
	 *
	 * @param dpi
	 *            The resolution in <em>dots per inch</em> at which the
	 *            fingerprint image has been scanned.
	 *
	 * @warning
	 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
	 *            are not greater than 0, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::initialize
		( const double *intensityImage , int m , int n , int dpi ) {

		if ( m <= 0 || n <= 0 || dpi <= 0 ) {
			cerr << "Fingerprint::initialize: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		clear();

		// Fill in the fields needed to represent an initialized
		// fingerprint...
		this->dpi = dpi;
		this->m   = m;
		this->n   = n;
		// ... . In particular, allocate an array to hold the intensity
		// image for this fingerprint and ...
		this->intensityImage = (double*)malloc( m * n * sizeof(double) );
		if ( this->intensityImage == NULL ) {
			cerr << "Fingerprint::initialize: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}
		// ... copy the input intensities to the array.
		memcpy
			( this->intensityImage , intensityImage ,
			  m * n * sizeof(double) );

		// The orientation image should be an array of NULL
		// pointers on initialization. Therefore, allocate the an
		// array of pointers to an 'Orientation' and ...
		this->orientationImage = new(nothrow) Orientation*[m*n];
		if ( this->orientationImage == NULL ) {
			cerr << "Fingerprint:: initialize: out of memory." << endl;
		}

		// ... initialize its elements with 'NULL'
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				this->orientationImage[y*n+x] = NULL;
			}
		}
	}

	/**
	 * @brief
	 *             Initialize this instance with a new fingerprint
	 *             specified by its raw 8-bit gray-scale intensity image.
	 *
	 * @param rawImage
	 *            An array of <i>m*n</i> values of type <code>uint8_t</code>
	 *            ranging between 0 (black) and 255 (white) where the
	 *            intensity at <i>(x,y)</i> is stored at
	 *            <code>intensityImage[y*n+x]</code>.
	 *
	 * @param m
	 *            Height of the fingerprint image.
	 *
	 * @param n
	 *            Width of the fingerprint image.
	 *
	 * @param dpi
	 *            The resolution in <em>dots per inch</em> at which the
	 *            fingerprint image has been scanned.
	 *
	 * @warning
	 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
	 *            are not greater than 0, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::initialize
	( const uint8_t *rawImage , int m , int n , int dpi ) {

		if ( m <= 0 || n <= 0 || dpi <= 0 ) {
			cerr << "Fingerprint::initialize: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		clear();

		// Fill in the fields needed to represent an initialized
		// fingerprint ...
		this->dpi = dpi;
		this->m   = m;
		this->n   = n;
		// ... . In particular, allocate an array to hold the intensity
		// image for this fingerprint and ...
		this->intensityImage = (double*)malloc( m * n * sizeof(double) );
		if ( this->intensityImage == NULL ) {
			cerr << "Fingerprint::initialize: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// ... convert the integer pixel intensities between [0,255] to
		// doubles between [0.0,1.0]
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				this->intensityImage[y*n+x] = (double)rawImage[y*n+x]/255.0;
			}
		}

		// The orientation image should be an array of NULL
		// pointers on initialization. Therefore, allocate the an
		// array of pointers to an 'Orientation' and ...
		this->orientationImage = new(nothrow) Orientation*[m*n];
		if ( this->orientationImage == NULL ) {
			cerr << "Fingerprint:: initialize: out of memory." << endl;
		}

		// ... initialize its elements with 'NULL'
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				this->orientationImage[y*n+x] = NULL;
			}
		}
	}

	/**
	 * @brief
	 *            Assignment operator (procedural version).
	 *
	 * @param fingerprint
	 *            The fingerprint assigned to this instance.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::assign( const Fingerprint & fingerprint ) {

		// Nothing to do if self-assignment
		if ( (void*)this == (void*)&fingerprint ) {
			return;
		}

		// If the to-be-copied fingerprint is empty, ...
		if ( fingerprint.isEmpty() ) {
			// ..., make the assigned fingerprint empty, too.
			clear();
			return;
		}

		// Initialize this fingerprint by the intensty image of the
		// to-be-assigned fingerprint.
		initialize
			( fingerprint.intensityImage,fingerprint.m,fingerprint.n,
			  fingerprint.dpi);

		// It remains to fulfill any features already stored by the
		// to-be-assigned fingerprint.

		int m , n;
		m = this->m;
		n = this->n;

		// Copy the foreground image if held by 'fingerprint'.
		if ( fingerprint.foregroundImage != NULL ) {
			this->foregroundImage = (bool*)malloc( m * n * sizeof(bool) );
			if ( this->foregroundImage == NULL ) {
				cerr << "Fingerprint::assign: Out of memory." << endl;
				exit(EXIT_FAILURE);
			}
			memcpy
				( this->foregroundImage , fingerprint.foregroundImage ,
				  m * n * sizeof(bool) );
		}

		// Copy non-NULL orientation estimations held by 'fingerprint'
		for ( int y = 0 ; y < m ; y++ ) {
			for ( int x = 0 ; x < n ; x++ ) {
				if ( fingerprint.orientationImage[y*n+x] != NULL ) {
					this->orientationImage[y*n+x] = new(nothrow) Orientation
							(fingerprint.orientationImage[y*n+x][0]);
					if ( this->orientationImage[y*n+x] == NULL ) {
						cerr << "Fingerprint::assign: Out of memory." << endl;
						exit(EXIT_FAILURE);
					}
				}
			}
		}

		// Copy directed reference point estimation held by 'fingerprint'
		this->has_drp = fingerprint.has_drp;
		if ( this->drpPtr != NULL ) {
			this->drpPtr = new(nothrow) DirectedPoint(*(fingerprint.drpPtr));
			if ( this->drpPtr == NULL ) {
				cerr << "Fingerprint::assign: our of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		if ( fingerprint.recordPtr == NULL ) {
			delete this->recordPtr;
			this->recordPtr = NULL;
		} else {
			if ( this->recordPtr == NULL ) {
				this->recordPtr = new(nothrow) thimble::MinutiaeRecord
						(fingerprint.recordPtr[0]);
				if ( this->recordPtr == NULL ) {
					cerr << "kbeinweg::teil1::ap3::Fingerprint: "
						 << "out of memory." << endl;
					exit(EXIT_FAILURE);
				}
			} else {
				this->recordPtr[0] = fingerprint.recordPtr[0];
			}

		}

	}

	/**
	 * @brief
	 *            Assignment operator.
	 *
	 * @param fingerprint
	 *            The fingerprint assigned to this instance.
	 *
	 * @return
	 *            A reference to this instance.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	Fingerprint & Fingerprint::operator=( const Fingerprint & fingerprint ) {
		assign(fingerprint);
		return *this;
	}

	/**
	 * @brief
	 *            Access the height of the fingerprint image.
	 *
	 * @return
	 *            The height of the fingerprint image.
	 */
	int Fingerprint::getHeight() const {
		return this->m;
	}

	/**
	 * @brief
	 *            Access the width of the fingerprint image.
	 *
	 * @return
	 *            The width of the fingerprint image.
	 */
	int Fingerprint::getWidth() const {
		return this->n;
	}

    /**
     * @brief
     *            Access the resolution specified at which this
     *            fingerprint's image has been scanned.
     * @return
     *            The resolution in <i>dots per inch</i>.
     */
    int Fingerprint::getResolution() const {
        return this->dpi;
    }

	/**
	 * @brief
	 *            Access the intensity image of the original
	 *            fingerprint image.
	 *
	 * @details
	 *            Assume that
	 *            <pre>
	 *             m = getHeight();
	 *             n = getWidth();
	 *            </pre>
	 *            is the height and width of the fingerprint,
	 *            respectively, and
	 *            <pre>
	 *             const double *img = getIntensityImage()
	 *            </pre>
	 *            is the array returned by the function. Then
	 *            <code>img</code> encodes the intensity image as
	 *            an <i>m*n</i> matrix in row-major order. Specifically,
	 *            the intensity at pixel <i>(x,y)</i>
	 *            (where <i>x=0...m-1, y=0...n-1</i>) can be accessed via
	 *            <pre>
	 *             double intensitiy = img[y*n+x];
	 *            </pre>
	 *
	 * @return
	 *            Intensity image of the original fingerprint image.
	 */
	const double *Fingerprint::getIntensityImage() const {
		return this->intensityImage;
	}

	/**
	 * @brief
	 *            Estimates the foreground of the fingerprint and
	 *            stores the result in the specified array.
	 *
	 * @details
	 *            Assume that
	 *            <pre>
	 *             m = getHeight();
	 *             n = getWidth();
	 *            </pre>
	 *            is the height and width of the fingerprint,
	 *            respectively, and
	 *            <pre>
	 *             bool *img = ...;
	 *            </pre>
	 *            is allocated to hold 'm*n' values of type
	 *            <code>bool</code> to which the output of the method
	 *            is stored via
	 *            <pre>
	 *             estimateForegroundImage(img)
	 *            </pre>
	 *            encoded as an <i>m*n</i> matrix in row-major order.
	 *            Specifically, the boolean indicating whether the pixel
	 *            <i>(x,y)</i> (where <i>x=0...m-1, y=0...n-1</i>) belongs
	 *            to the foreground can be accessed via
	 *            <pre>
	 *             bool isForegroundPixel = img[y*n+x];
	 *            </pre>
	 *
	 *
	 * @param foregroundImage
	 *            Array that can hold
	 *            <code>getWidth()*getHeight()</code> values of the
	 *            type <code>bool</code>.
	 *
	 * @warning
	 *            If <code>foregroundImage</code> cannot store
	 *            <code>getWidth()*getHeight()</code> boolean
	 *            values, the method runs into undocumented behavior.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	bool Fingerprint::isForeground( int y , int x ) {

		int m , n;
		m = getHeight();
		n = getWidth();

		if ( y < 0 || x < 0 || y >= m || x >= n ) {
			cerr << "Fingerprint::isForeground: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// Ensures that the foreground is estimated
		getForegroundImage();

		return this->foregroundImage[y*n+x];
	}

	/**
	 * @brief
	 *            Access the foreground image of the fingerprint.
	 *
	 * @details
	 *            Assume that
	 *            <pre>
	 *             m = getHeight();
	 *             n = getWidth();
	 *            </pre>
	 *            is the height and width of the fingerprint,
	 *            respectively, and
	 *            <pre>
	 *             const bool *img = getForegroundImage();
	 *            </pre>
	 *            is the array returned by the function. Then
	 *            <code>img</code> encodes the foreground image as
	 *            an <i>m*n</i> matrix in row-major order.
	 *            Specifically, the boolean indicating whether the pixel
	 *            <i>(x,y)</i> (where <i>x=0...m-1, y=0...n-1</i>) belongs
	 *            to the foreground can be accessed via
	 *            <pre>
	 *             bool isForegroundPixel = img[y*n+x];
	 *            </pre>
	 *
	 *            If the function is called the first time for the
	 *            fingerprint, it causes a call of the
	 *     \link thimble::Fingerprint::estimateForegroundImage(bool*)const
	 *      estimateForegroundImage()\endlink which
	 *            runs an algorithm for fingerprint foreground estimation.
	 *            The result is kept by the fingerprint, thereby
	 *            preventing that the foreground is estimated again.
	 *
	 * @return
	 *            Foreground image of the fingerprint image.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	const bool *Fingerprint::getForegroundImage() {

		// If no foreground has been estimated for this fingerprint, ...
		if ( this->foregroundImage == NULL ) {

			// ..., allocate memory to hold the foreground and ...
			int m , n;
			m = this->m;
			n = this->n;
			this->foregroundImage = (bool*)malloc( m * n * sizeof(bool) );
			if ( this->foregroundImage == NULL ) {
				cerr << "Fingerprint::getForegroundImage: out of memory."
					 << endl;
				exit(EXIT_FAILURE);
			}

			// ... do estimate the foreground.
			estimateForegroundImage(this->foregroundImage);
		}

		return this->foregroundImage;
	}

	/**
	 * @brief
	 *            Access the orientation of the fingerprint at the
	 *            specified pixel.
	 *
	 * @details
	 *            If the orientation at <i>(x,y)</i> is accessed the first
	 *            time, it is estimated using the
	 *            \link estimateOrientationAt(int,int)\endlink
	 *            function and the result is kept stored by the instance.
	 *            Otherwise, if already computed, the function returns a
	 *            reference to the stored orientation.
	 *
	 *            The function calls the low-level
	 *            \link estimateOrientationAt(int,int)\endlink function
	 *            if an orientation is needed to be estimated.
	 *
	 * @param y
	 *            Ordinate varying between 0 and
	 *            \link getHeight()\endlink-1.
	 *
	 * @param x
	 *            Abscissa varying between 0 and
	 *            \link getWidth()\endlink-1.
	 *
	 * @return
	 *            %Orientation at <i>(x,y)</i>.
	 *
	 * @warning
	 *            If <i>(x,y)</i> does not lie between the fingerprint
	 *            image's region, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	const Orientation & Fingerprint::getOrientationAt( int y , int x ) {

		int m , n;
		m = getHeight();
		n = getWidth();

		// Ensure that the accessed pixel lies on the fingerprint image.
		if ( y < 0 || x < 0 || y >= m || x >= n ) {
			cerr << "Fingerprint::getOrientationAt: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// If no orientation estimation is associated for the
		// accessed pixel, ...
		if ( this->orientationImage[y*n+x] == NULL ){
			// ..., estimate the orientation at the pixel and ...
			Orientation v = estimateOrientationAt(y,x);
			// ... store the estimation in the orientation image.
			this->orientationImage[y*n+x] = new(nothrow) Orientation(v);
			if ( this->orientationImage[y*n+x] == NULL ) {
				cerr << "Fingerprint::getOrientationAt: but of memory."
					 << endl;
				exit(EXIT_FAILURE);
			}
		}

		return this->orientationImage[y*n+x][0];
	}

	/**
	 * @brief
	 *            Estimates the foreground of the fingerprint and
	 *            stores the result in the specified array.
	 *
	 * @details
	 *            Assume that
	 *            <pre>
	 *             m = getHeight();
	 *             n = getWidth();
	 *            </pre>
	 *            is the height and width of the fingerprint,
	 *            respectively, and
	 *            <pre>
	 *             bool *img = ...;
	 *            </pre>
	 *            is allocated to hold 'm*n' values of type
	 *            <code>bool</code> to which the output of the method
	 *            is stored via
	 *            <pre>
	 *             estimatedForegroundImage(img)
	 *            </pre>
	 *            encoded as an <i>m*n</i> matrix in row-major order.
	 *            Specifically, the boolean indicating whether the pixel
	 *            <i>(x,y)</i> (where <i>x=0...m-1, y=0...n-1</i>) belongs
	 *            to the foreground can be accessed via
	 *            <pre>
	 *             bool isForegroundPixel = img[y*n+x];
	 *            </pre>
	 *
	 *
	 * @param foregroundImage
	 *            Array that can hold
	 *            <code>getWidth()*getHeight()</code> values of the
	 *            type <code>bool</code>.
	 *
	 * @warning
	 *            If <code>foregroundImage</code> cannot store
	 *            <code>getWidth()*getHeight()</code> boolean
	 *            values, the method runs into undocumented behavior.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::estimateForegroundImage( bool *foregroundImage ) const {

		int h = (int)ceil(5.0/569.0*(double)(this->dpi));

		Segmentation(h).estimate
				(foregroundImage,getIntensityImage(),getHeight(),getWidth());
	}

    /**
     * @brief
     *            Controls the number of gradients used to estimate the
     *            local orientations of the fingerprint with the
     *            gradient method.
     *
     * @details
     *            Unless overwritten by an extension class,
     *            this function is used for estimating local ridge
     *            orientations of the fingerprint with the gradient
     *            method, i.e.,
     *            by \link estimateOrientationAt()\endlink.
     *
     * @return
     *            The function returns
     *            '(int)ceil(8.0/569.0*getResolution())'
     */
    int Fingerprint::getGradientOrientationDistance() const {
        return (int)ceil(8.0/569.0*(this->dpi));
    }

	/**
	 * @brief
	 *            Estimates the orientation of the fingerprint at the
	 *            specified pixel.
	 *
	 * @details
	 *            The orientation is estimated using
	 * \link Orientation::gradient(int,int,const double*,int,int,int)\endlink
	 *            function.
	 *
	 * @param y
	 *            Ordinate varying between 0 and
	 *            \link getHeight()\endlink-1.
	 *
	 * @param x
	 *            Abscissa varying between 0 and
	 *            \link getWidth()\endlink-1.
	 *
	 * @return
	 *            %Orientation at <i>(x,y)</i>.
	 *
	 * @warning
	 *            If <i>(x,y)</i> does not lie between the fingerprint
	 *            image's region, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 */
	Orientation Fingerprint::estimateOrientationAt( int y , int x ) const {

		// How many gradients surrounding the accessed pixel are involved
		// depend on the specified resolution
        int h = getGradientOrientationDistance();

		int m , n;
		m = getHeight();
		n = getWidth();

		// If the accessed pixel's orientation estimation would suffer
		// from edge effects, return a zero orientation
        if ( y < h    || x < h ||
             y+h >= m || x+h >= n ) {
			return Orientation();
		}

		// Orientation estimation via gradients.
        return Orientation::gradient(y,x,getIntensityImage(),m,n,h);
	}

	/**
	 * @brief
	 *            Estimates a directed reference point from this
	 *            fingerprint and returns the result
	 *
	 * @details
	 *            The implementation of the directed reference point
	 *            estimation is performed using the
	 *            \link TentedArchModel\endlink class following the
	 *            procedure as described in
	 *            <ul>
	 *             <li><b>B. Tams (2013)</b>.
	 *             Absolute %Fingerprint Pre-Alignment in Minutiae-Based
	 *             Cryptosystems. <i>Proc. BIOSIG 2013, ser LNI,
	 *             <b>vol. 212</b>, pp. 75--86</i>.</li>
	 *            </ul>
	 *            but can be overridden for a subclass, e.g., by
	 *            a more robust directed reference point estimation.
	 *            <br><br>
	 *            If the result of estimation is invalid, at least one
	 *            of the returned directed reference point's
	 *            <pre>
	 *             DirectedPoint drp = ...
	 *            </pre>
	 *            members
	 *            \link DirectedPoint::x drp.x\endlink,
	 *            \link DirectedPoint::y drp.y\endlink or
	 *            \link Direction::getDirectionAngle()
	 *                  drp.direction.getDirectionAngle()\endlink
	 *            should be equals NAN, i.e., self-comparision results
	 *            in a false boolean expression.
	 *            <br><br>
	 *            Overriding this method affects the behavior (but not
	 *            the meaning) of the
	 *            \link hasDirectedReferencePoint()\endlink and
	 *            \link getDirectedReferencePoint()\endlink function.
	 *
	 * @return
	 *            A reference point attached with a direction.
	 *
	 * @see hasDirectedReferencePoint()
	 * @see getDirectedReferencePoint()
	 *
	 * @warning
	 *            If not enough memory could be provided, an error message
	 *            is printed to stderr and the program exits with status
	 *            'EXIT_FAILURE'.
	 */
	DirectedPoint Fingerprint::estimateDirectedReferencePoint() {

		DirectedPoint drp;

		// Perform the tented arch reference point algorithm; ...
		TentedArchModel model(this->dpi);
		if ( model.estimate(*this) ) {
			// ...; if successful, initialize the attributes of the
			// reference point correspondingly;
			drp.x = model.getComplexCore().real();
			drp.y = model.getComplexCore().imag();
			drp.direction.setDirectionAngle(model.getReferencePointDirection());
		} else {
			// otherwise, make the directed reference point to contain
			// invalid data and ...
			drp.x = THIMBLE_NAN;
			drp.y = THIMBLE_NAN;
			drp.direction.setDirectionAngle(THIMBLE_NAN);
		}

		return drp;
	}

	/**
	 * @brief
	 *            Estimates a minutiae template from this fingerprint
	 *            (Note: this function is not yet implemented).
	 *
	 * @details
	 *            Currently, the implementation of this function does
	 *            nothing else as returning an empty minutiae view.
	 *
	 *            In order to allow programs using the THIMBLE library to
	 *            implement this function, it is provided as a virtual
	 *            function which can be overridden, thereby affecting the
	 *            \link getMinutiaeRecord()\endlink and
	 *            \link getMinutiaeView()\endlink functions.
	 *
	 * @return
	 *            A minutiae view containing the minutiae estimated from
	 *            this fingerprint.
	 *
	 * @see getMinutiaeView()
	 * @see getMinutiaeRecord()
	 *
	 * @attention
	 *            This function is not yet implemented and does nothing
	 *            than returning an empty set of minutiae thereby causing
	 *            the \link getMinutiaeView()\endlink function to return
	 *            a set of empty minutiae (unless initialized manually).
	 */
	MinutiaeView Fingerprint::estimateMinutiae() {

		// Return an empty set of minutiae.
		return MinutiaeView();
	}

	/**
	 * @brief
	 *            Specifies the foreground image manually.
	 *
	 * @details
	 *            The specified <code>foregroundImage</code> should be
	 *            encoded as a
	 *  <code>\link getWidth()\endlink*\link getHeight()\endlink</code>
	 *            matrix in row-major order. More specifically, the value
	 *       at <code>foregroundImage[y*\link getWidth()\endlink+x]</code>
	 *            is <code>true</code> if the pixel at <i>(x,y)</i> is
	 *            specified as the foreground and, otherwise, if specified
	 *            as the background the value of
	 *            <code>foregroundImage[y*getWidth()+x]</code> is
	 *            <code>false</code>.
	 *
	 * @param foregroundImage
	 *            Array containing <code>getWidth()*getHeight()</code>
	 *            boolean values.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::setForegroundImage( const bool *foregroundImage ) {

		// If the outcome of 'getForegroundImage()' is passed to the
		// method, there is nothing to do.
		if ( this->foregroundImage == foregroundImage ) {
			return;
		}

		int m , n;
		m = getHeight();
		n = getWidth();

		// If not data has been allocated to hold the foreground
		// image, ...
		if ( this->foregroundImage == NULL ) {
			// ..., allocate memory; otherwise, the allocated space is reused.
			this->foregroundImage = (bool*)malloc( m * n * sizeof(bool) );
			if ( this->foregroundImage == NULL ) {
				cerr << "Fingerprint::setForegroundImage: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		// Copy the memory.
		memcpy(this->foregroundImage,foregroundImage,m*n*sizeof(bool));
	}

	/**
	 * @brief
	 *            Specifies the orientation of the fingerprint at a
	 *            pixel manually.
	 *
	 * @param orientation
	 *            The orientation assigned to the fingerprint at pixel
	 *            <i>(x,y)</i>.
	 *
	 * @param y
	 *            Ordinate varying between 0 and
	 *            \link getHeight()\endlink-1.
	 *
	 * @param x
	 *            Abscissa varying between 0 and
	 *            \link getWidth()\endlink-1.
	 *
	 * @warning
	 *            If <i>(x,y)</i> does not lie between the fingerprint
	 *            image's region, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::setOrientationAt
		( const Orientation & orientation , int y , int x ) {

		int m , n;
		m = getHeight();
		n = getWidth();

		// Ensure that the accessed pixel is valid.
		if ( y < 0 || x < 0 || y >= m || x >= n ) {
			cerr << "Fingerprint::setOrientationAt: "
				 << "bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		// If no orientation estimation is provided at the
		// specified pixel, ...
		if ( this->orientationImage[y*n+x] == NULL ) {
			// call the copy constructor; ...
			this->orientationImage[y*n+x] =
						new(nothrow) Orientation(orientation);
			if ( this->orientationImage[y*n+x] == NULL ) {
				cerr << "Fingerprint::setOrientationAt: out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			// ...; otherwise, overwrite the present orientation estimation.
			this->orientationImage[y*n+x][0] = orientation;
		}
	}

	/**
	 * @brief
	 *            Access the directed reference point of the fingerprint.
	 *
	 * @details
	 *            If no attempt for directed reference point estimation
	 *            has been run for this fingerprint, the function causes
	 *            an attempt for directed reference point estimation
	 *            by calling the
	 *            \link estimateDirectedReferencePoint()\endlink function
	 *            and keeps track of the result as well as of whether the
	 *            estimation is considered to be successful, thereby
	 *            preventing that the estimation is run another time on
	 *            any further calls of the function.
	 *
	 *            Independently of whether the estimation was successful,
	 *            the function returns a reference to the directed
	 *            reference point but the result should be considered to
	 *            be invalid if the estimation failed.
	 *
	 *            Before the directed reference point is accessed, one
	 *            should test whether the estimation was successful or
	 *            not, for example, by
	 *            <pre>
	 *        %Fingerprint fingerprint = ...;
	 *
	 *        if ( fingerprint.hasDirectedReferencePoint() ) {
	 *           cout << fingerprint.getDirectedReferencePoint() << endl;
	 *        }
	 *            </pre>
	 *            which causes to write the directed reference point
	 *            representation to the standard output stream if the
	 *            estimation was successful.
	 *
	 *            Overriding the virtual
	 *            \link estimateDirectedReferencePoint()\endlink function
	 *            by a possibly more robust method for directed reference point
	 *            estimation may affect the behavior of this method.
	 *
	 * @return
	 *            A constant reference to the result of this fingerprint's
	 *            directed reference point estimation.
	 *
	 * @see estimateDirectedReferencePoint()
	 * @see hasDirectedReferencePoint()
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	const DirectedPoint & Fingerprint::getDirectedReferencePoint() {

		// If no directed reference point estimation has been attempted
		// for this fingerprint, ...
		if ( this->drpPtr == NULL ) {

			this->drpPtr = new(nothrow) DirectedPoint
					(estimateDirectedReferencePoint());
			if ( this->drpPtr == NULL ) {
				cerr << "Fingerprint::getDirectedReferencePoint: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}

			this->has_drp =
			(this->drpPtr->x == this->drpPtr->x) &&
			(this->drpPtr->y == this->drpPtr->y) &&
			(this->drpPtr->direction.getDirectionAngle() ==
			 this->drpPtr->direction.getDirectionAngle());
		}

		// Return the directed reference point.
		return *(this->drpPtr);
	}

	/**
	 * @brief
	 *            Checks whether directed reference point estimation was
	 *            successful for this fingerprint.
	 *
	 * @details
	 *            If for the fingerprint no directed reference point has
	 *            been estimated, the function calls the virtual
	 *            \link estimateDirectedReferencePoint()\endlink function,
	 *            keeps track of the result and of whether the directed
	 *            reference point is valid (i.e., whether none of its
	 *            members equals NAN) and returns <code>true</code> if
	 *            the directed reference point is valid and
	 *            <code>false</code> otherwise. Consequently, if an attempt
	 *            for directed reference point estimation has been
	 *            performed for this fingerprint before, the function
	 *            returns <code>true</code> if the directed reference
	 *            point of the earlier estimation was valid and
	 *            <code>false</code> otherwise.
	 *
	 * @return
	 *            <code>true</code> if directed reference point estimation
	 *            was successful and, otherwise, if the estimation failed,
	 *            the function returns <code>false</code>.
	 *
	 * @see estimateDirectedReferencePoint()
	 * @see getDirectedReferencePoint()
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	bool Fingerprint::hasDirectedReferencePoint() {

		// Ensure that an attempt for directed reference point estimation
		// has been made.
		getDirectedReferencePoint();

		// Return whether the estimation was successful.
		return this->has_drp;
	}

	/**
	 * @brief
	 *            Assigns the directed reference point manually.
	 *
	 * @details
	 *            The method assigns the directed point <code>dp</code> as
	 *            the fingerprint's directed reference point, thereby
	 *            replacing any former estimations or assignments of
	 *            the fingerprint's directed reference point.
	 *
	 *            When a directed reference point is specified manually,
	 *            it is assumed that it is valid, i.e., the function
	 *            \link hasDirectedReferencePoint()\endlink returns
	 *            <code>true</code> after manual directed reference point
	 *            assignment. In particular, the direction and pixel
	 *            coordinates of <code>dp</code> should hold well-defined
	 *            values. Otherwise, the use of this fingerprint instance
	 *            can run into undocumented behavior.
	 *
	 * @param dp
	 *            The valid directed reference point assigned to the
	 *            fingerprint.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::setDirectedReferencePoint( const DirectedPoint & dp ) {

		// If no attempt for directed reference point estimation has been
		// attempted for this fingerprint, ...
		if ( this->drpPtr == NULL ) {
			// ..., call the copy constructor to make a copy of the
			// specified reference point; ...
			this->drpPtr = new(nothrow) DirectedPoint(dp);
			if ( this->drpPtr == NULL ) {
				cerr << "Fingerprint::setDirectedReferencePoint: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			// ...; otherwise, reuse the memory for assignment.
			*(this->drpPtr) = dp;
		}

		// Assume that the specified reference point is valid for this
		// fingerprint.
		this->has_drp = true;
	}

	/**
	 * @brief
	 *            Returns a minutiae record containing the minutiae
	 *            template estimated from this fingerprint.
	 *
	 * @details
	 *            If no attempt for minutiae estimation has been performed
	 *            for this fingerprint, this function creates a
	 *            \link MinutiaeRecord\endlink (with parameters specified
	 *            for this fingerprint) that contains an empty
	 *            \link MinutiaeView\endlink; then the view is assigned
	 *            with the result of the
	 *            \link estimateMinutiae\endlink function; finally, the
	 *            reference to the \link MinutiaeRecord\endlink is returned.
	 *            Otherwise, if a \link MinutiaeRecord\endlink has already
	 *            been created for this fingerprint on a former attempt for
	 *            minutiae estimation, this function returns the reference
	 *            to the already existend \link MinutiaeRecord\endlink.
	 *
	 * @return
	 *            A reference to the minutiae record containing a minutiae
	 *            template estimated from this fingerprint.
	 *
	 * @see estimateMinutiae()
	 * @see getMinutiaeView()
	 */
	const MinutiaeRecord & Fingerprint::getMinutiaeRecord() {

		if ( this->recordPtr == NULL ) {

			int m , n , dpi , ppc;
			m = getWidth();
			n = getHeight();
			dpi = getResolution();
			ppc = (int)(dpi/2.54);

			this->recordPtr = new(nothrow) thimble::MinutiaeRecord
					(m,n,ppc);
			if ( this->recordPtr == NULL ) {
				cerr << "Fingerprint::getMinutiaeRecord(): "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}

			this->recordPtr->addView(0);
			this->recordPtr->getView(0) = estimateMinutiae();
		}

		return *(this->recordPtr);
	}

	/**
	 * @brief
	 *            Returns a minutiae template estimated from this
	 *            fingerprint.
	 *
	 * @details
	 *            If no attempt for minutiae estimation has been performed
	 *            for this fingerprint, this function creates a
	 *            \link MinutiaeRecord\endlink (with parameters specified
	 *            for this fingerprint) that contains an empty
	 *            \link MinutiaeView\endlink; then the view is assigned
	 *            with the result of the
	 *            \link estimateMinutiae\endlink function; finally, the
	 *            reference to the \link MinutiaeView\endlink of the
	 *            \link MinutiaeRecord\endlink is returned.
	 *            Otherwise, if a minutiae template estimation attempt
	 *            has already been performed for this fingerprint, this
	 *            function returns the reference to the already existent
	 *            \link MinutiaeView\endlink.
	 *
	 *            Currently, the implementation of the
	 *            \link estimateMinutiae()\endlink function does nothing
	 *            than returning an empty minutiae template. Consequently,
	 *            by merely using the functionalities provided by this
	 *            class, the user will obtain only empty minutiae
	 *            templates. However, the minutiae templates estimations
	 *            can be specified manually. Therefore, assume that a
	 *            minutiae record has been initialized from the file
	 *            'template.fmr' by
	 *            <pre>
	 *             MinutiaeRecord record;
	 *
	 *             record.read("template.fmr");
	 *            </pre>
	 *            which is assumed to have at least one valid
	 *            \link MinutiaeView\endlink which can be accessed by
	 *            <pre>
	 *             record.getView(0)
	 *            </pre>
	 *            Then the \link MinutiaeView\endlink of this fingerprint
	 *            <pre>
	 *             %Fingerprint fingerprint;
	 *            </pre>
	 *            can be initialized manually via
	 *            <pre>
	 *             fingerprint.getMinutiaeView() = record.getView(0);
	 *            </pre>
	 *
	 * @return
	 *            A reference to the minutiae view containing the minutiae
	 *            set estimated from this fingerprint.
	 *
	 * @see estimateMinutiae()
	 * @see getMinutiaeRecord()
	 */
	MinutiaeView & Fingerprint::getMinutiaeView() {

		getMinutiaeRecord();

		return this->recordPtr->getView(0);
	}

	/**
	 * @brief
	 *            Specify the minutiae for this fingerprint manually.
	 *
	 * @param view
	 *            The minutiae view containing the new minutiae specified
	 *            for this fingerprint object.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void Fingerprint::setMinutiaeView( const MinutiaeView & view ) {

		if  (this->recordPtr == NULL ) {

			int m , n , dpi , ppc;
			m = getWidth();
			n = getHeight();
			dpi = getResolution();
			ppc = (int)(dpi/2.54);

			this->recordPtr = new(nothrow) thimble::MinutiaeRecord
					(m,n,ppc);
			if ( this->recordPtr == NULL ) {
				cerr << "Fingerprint::setMinutiaeRecord(): "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		if ( recordPtr->getViewCount() == 0 ) {
			this->recordPtr->addView(view);
		} else {
			this->recordPtr->getView() = view;
		}
	}

	/**
	 * @brief
	 *            Initialize the fingerprint from a PGM gray-scale
	 *            image file.
	 *
	 * @details
	 *            This function serves as a convenience procedure to
	 *            prevent the user from temporarily creating
	 *            \link GrayImage\endlink objects. More
	 *            specifically, calling the function is equivalent to the
	 *            following:
	 *            <pre>
	 *             %Fingerprint fingerprint = ...; // This instance
	 *             GrayImage pgmImage;
	 *             if ( pgmImage.read(pgmImagePath) ) {
	 *                fingerprint.initialize(pgmImage,dpi);
	 *                return true;
	 *             } else {
	 *                return false;
	 *             }
	 *            </pre>
	 *
	 * @param pgmImagePath
	 *            Specifies the path to the PGM gray-scale fingerprint
	 *            image.
	 *
	 * @param dpi
	 *            The resolution in dots per inch at which the fingerprint
	 *            image has been scanned.
	 *
	 * @return
	 *            <code>true</code> if the fingerprint has been
	 *            successfully initialized; otherwise, the function
	 *            returns <code>false</code>.
	 *
	 * @warning
	 *             If <code>dpi</code> is not greater than 0, an error message
	 *             is printed to <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	bool Fingerprint::fromImageFile( const std::string & pgmImagePath , int dpi ) {

		GrayImage pgmImage;

		if ( !pgmImage.read(pgmImagePath) ) {
			return false;
		}

		this->initialize(pgmImage,dpi);

		return true;
	}
}



