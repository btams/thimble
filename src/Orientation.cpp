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
 * @file Orientation.cpp
 *
 * @brief
 *            Implementation of the functionalities provided by the
 *            'Orientation.h' header file computing orientations from
 *            images and representation of orientations, directions and
 *            two-dimensional point coordinates attached with orientations.
 *            and directions.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include "config.h"
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ostream>
#include <iostream>

#include <thimble/image/Orientation.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Creates an orientation of specified angle and
	 *            coherence (also serves as the standard constructor).
	 *
	 * @details
	 *            Calls \link setOrientationAngle()\endlink and then
	 *            \link setCoherence()\endlink to create a new
	 *            \link Orientation\endlink object.
	 *
	 * @param angle
	 *            Specified orientation angle which should be a value
	 *            from the interval \f$[0,\pi)\f$.
	 *
	 * @param coherence
	 *            The coherence of the orientation indicating how
	 *            reliable the orientation estimation is. Therein,
	 *            a value close to 0 indicates a very unreliable
	 *            estimation and a value close 1 a very reliable
	 *            estimation.
	 *
	 * @warning
	 *            If <code>coherence</code> is not a value from
	 *            \f$[0,1]\f$ an error message will be printed to
	 *            <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 */
	Orientation::Orientation( double orientation , double coherence ) {
		setOrientationAngle(orientation);
		setCoherence(coherence);
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @param orientation
	 *            The orientation of which a copy is created.
	 */
	Orientation::Orientation( const Orientation & v ) {
		this->assign(v);
	}

	/**
	 * @brief
	 *            Updates the orientation angle.
	 *
	 * @details
	 *            The orientation angle should be value between
	 *            \f$[0,\pi)\f$. If <code>angle</code> is a value
	 *            outside of this interval, it is normalized via
	 *            the \link normalize()\endlink function.
	 *
	 *            Replacing the orientation angle affects the
	 *            result of the functions
	 *            <ul>
	 *             <li>\link getOrientationAngle()\endlink</li>
	 *             <li>\link get_dx()\endlink</li>
	 *             <li>\link get_dy()\endlink</li>
	 *             <li>\link get_dx2()\endlink</li>
	 *             <li>\link get_dy2()\endlink</li>
	 *            </ul>
	 *
	 * @param angle
	 *            Specified orientation angle which should be a value
	 *            from the interval \f$[0,\pi)\f$.
	 */
	void Orientation::setOrientationAngle( double orientation ) {

		this->angle = normalize(orientation);

		this->dx = cos(this->angle);
		this->dy = sin(this->angle);

		// The cosines and sines of the doubled orientation angles
		// can be computed via squaring the complex number 'dx+i*dy'
		this->dx2 = this->dx * this->dx - this->dy * this->dy;
		this->dy2 = this->dx * this->dy; this->dy2 += this->dy2;
	}

	/**
	 * @brief
	 *            Updates the coherence.
	 *
	 * @param coherence
	 *            The coherence of the orientation indicating how
	 *            reliable the orientation estimation is. Therein,
	 *            a value close to 0 indicates a very unreliable
	 *            estimation and a value close 1 a very reliable
	 *            estimation.
	 *
	 * @warning
	 *            If <code>coherence</code> is not a value from
	 *            \f$[0,1]\f$ an error message will be printed to
	 *            <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 */
	void Orientation::setCoherence( double coherence ) {

		if ( coherence < 0.0 || coherence > 1.0 ) {
			cerr << "Orientation::setCoherence: Coherence must be between "
				 << "0.0 and 1.0." << endl;
			exit(EXIT_FAILURE);
		}

		this->coherence = coherence;
	}

	/**
	 * @brief
	 *            Assigns this by another orientation.
	 *
	 * @param orientation
	 *            The orientation assigned to this object.
	 */
	void Orientation::assign( const Orientation & v ) {
		this->angle       = v.angle;
		this->coherence   = v.coherence;
		this->dx          = v.dx;
		this->dy          = v.dy;
		this->dx2         = v.dx2;
		this->dy2         = v.dy2;
	}

	/**
	 * @brief
	 *            Assignment operator.
	 *
	 * @details
	 *            Using the '='-operator calls the method
	 *            \link assign()\endlink.
	 *
	 * @param orientation
	 *            The orientation assigned to this object.
	 *
	 * @return
	 *            A reference to this orientation (after assignment).
	 *
	 * @see assign()
	 */
	Orientation & Orientation::operator=( const Orientation & v ) {
		assign(v);
		return *this;
	}

	/**
	 * @brief
	 *            Access the angle of this orientation.
	 *
	 * @return
	 *            A value from the interval \f$[0,\pi)\f$.
	 */
	double Orientation::getOrientationAngle() const { return this->angle; }

	/**
	 * @brief
	 *            Access the coherence of this orientation.
	 *
	 * @return
	 *            A value from the interval \f$[0,1]\f$.
	 */
	double Orientation::getCoherence() const { return this->coherence; }

	/**
	 * @brief
	 *            Access the cosine of the orientation angle.
	 *
	 * @details
	 *            The cosine of the orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            cos(\link getOrientationAngle()\endlink)
	 */
	double Orientation::get_dx() const { return this->dx; }

	/**
	 * @brief
	 *            Access the sine of the orientation angle.
	 *
	 * @details
	 *            The sine of the orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            sin(\link getOrientationAngle()\endlink)
	 */
	double Orientation::get_dy() const { return this->dy; }

	/**
	 * @brief
	 *            Access the cosine of the doubled orientation angle.
	 *
	 * @details
	 *            The cosine of the doubled orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            cos(2.0*\link getOrientationAngle()\endlink)
	 */
	double Orientation::get_dx2() const { return this->dx2; }

	/**
	 * @brief
	 *            Access the sine of the doubled orientation angle.
	 *
	 * @details
	 *            The sine of the doubled orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            sin(2.0*\link getOrientationAngle()\endlink)
	 */
	double Orientation::get_dy2() const { return this->dy2; }

	/**
	 * @brief
	 *            Returns the orientation angle in degree rounded
	 *            to an integer.
	 *
	 * @return
	 *            A value from { 0 , ... , 179 }.
	 */
	int Orientation::getDegree() const {

		int degree;

		degree = (int)THIMBLE_ROUND
				(180.0/M_PI * (double)getOrientationAngle() );

		if ( degree >= 180 ) {
			degree -= 180;
		}

		return degree;
	}

	/**
	 * @brief
	 *            Returns a representative of the specified orientation
	 *            angle from the interval \f$[0,\pi)\f$.
	 *
	 * @details
	 *            Essentially, the function returns the value
	 *            \f[
	 *             angle' = k\cdot\pi+angle
	 *            \f]
	 *            where <i>k</i> is an integer such that
	 *            \f$angle'\in[0,\pi)\f$.
	 *
	 * @param angle
	 *            An (unnormalized) orientation angle.
	 *
	 * @return
	 *            The normalized orientation angle.
	 */
	double Orientation::normalize( double theta ) {

		if ( theta >= 0.0 && theta < M_PI ) {
			// Nothing to normalize
			return theta;
		}

		double dx , dy , normalizedTheta;
		dx = cos(theta);
		dy = sin(theta);

		normalizedTheta = atan2(dy,dx);

		while ( normalizedTheta < 0.0 ) {
			normalizedTheta += M_PI;
		}

		while ( normalizedTheta >= M_PI ) {
			normalizedTheta -= M_PI;
		}

		return normalizedTheta;
	}

	/**
	 * @brief
	 *            Computes a local x-gradient estimation at a
	 *            specified pixel coordinate of a two-dimensional
	 *            double-valued intensity image.
	 *
	 * @details
	 *            The local x-gradient is approximated via \f$3\times 3\f$
	 *            Sobel masks.
	 *
	 * @param y0
	 *            The y-value of the coordinate at which the local
	 *            orientation is estimated.
	 *
	 * @param x0
	 *            The x-value of the coordinate at which the local
	 *            orientation is estimated
	 *
	 * @param image
	 *            The double-valued intensities of the specified image;
	 *            where the intensity of the pixel at <code>(y,x)</code>
	 *            is stored at <code>image[y*n+x]</code>.
	 *
	 * @param m
	 *            The height of the specified image.
	 *
	 * @param n
	 *            The width of the specified image.
	 *
	 * @warning
	 *            In the following cases, the behavior of this function
	 *            is undocumented:
	 *            <ul>
	 *             <li><code>m</code> is negative</li>
	 *             <li><code>n</code> is negative</li>
	 *             <li><code>h</code> is negative</li>
	 *             <li>
	 *              <code>image</code> does not contain
	 *              <code>m*n</code> valid <code>double</code> values.
	 *             </li>
	 *             <li><code>y0</code> is not a value from 0,...,m-1
	 *             <li><code>x0</code> is not a value from 0,...,n-1
	 *            </ul>
	 */
	double Orientation::gradX
		( int y , int x , const double *image , int m , int n ) {


		double dx = 0.0;
		int count = 0;

		if ( x > 0 ) {

			if ( y > 0 ) {
				dx += image[(y-1)*n+x-1];
				++count;
			}

			dx += 2.0 * image[y*n+x-1];
			++count;

			if ( y+1 < m ) {
				dx += image[(y+1)*n+x-1];
				++count;
			}

		}

		if ( x+1 < n ) {

			if ( y > 0 ) {
				dx -= image[(y-1)*n+x+1];
				++count;
			}

			dx -= 2.0*image[y*n+x+1];
			++count;

			if ( y+1 < m ) {
				dx -= image[(y+1)*n+x+1];
				++count;
			}

		}

		return dx / (double)count;
	}

	/**
	 * @brief
	 *            Computes a local y-gradient estimation at a
	 *            specified pixel coordinate of a two-dimensional
	 *            double-valued intensity image.
	 *
	 * @details
	 *            The local y-gradient is approximated via \f$3\times 3\f$
	 *            Sobel masks.
	 *
	 * @param y0
	 *            The y-value of the coordinate at which the local
	 *            orientation is estimated.
	 *
	 * @param x0
	 *            The x-value of the coordinate at which the local
	 *            orientation is estimated
	 *
	 * @param image
	 *            The double-valued intensities of the specified image;
	 *            where the intensity of the pixel at <code>(y,x)</code>
	 *            is stored at <code>image[y*n+x]</code>.
	 *
	 * @param m
	 *            The height of the specified image.
	 *
	 * @param n
	 *            The width of the specified image.
	 *
	 * @warning
	 *            In the following cases, the behavior of this function
	 *            is undocumented:
	 *            <ul>
	 *             <li><code>m</code> is negative</li>
	 *             <li><code>n</code> is negative</li>
	 *             <li><code>h</code> is negative</li>
	 *             <li>
	 *              <code>image</code> does not contain
	 *              <code>m*n</code> valid <code>double</code> values.
	 *             </li>
	 *             <li><code>y0</code> is not a value from 0,...,m-1
	 *             <li><code>x0</code> is not a value from 0,...,n-1
	 *            </ul>
	 */
	double Orientation::gradY
		( int y , int x , const double *image , int m , int n ) {

		double dy = 0.0;
		int count = 0;

		if ( y > 0 ) {

			if ( x > 0 ) {
				dy += image[(y-1)*n+x-1];
				++count;
			}

			dy += 2.0*image[(y-1)*n+x];
			++count;

			if ( x+1 < n ) {
				dy += image[(y-1)*n+x+1];
				++count;
			}

		}

		if ( y+1 < m ) {
			if ( x > 0 ) {
				dy -= image[(y+1)*n+x-1];
				++count;
			}

			dy -= 2.0*image[(y+1)*n+x];
			++count;

			if ( x+1 < n ) {
				dy -= image[(y+1)*n+x+1];
				++count;
			}
		}

		return dy / (double)count;
	}

	/**
	 * @brief
	 *            Estimates the local orientation at the specified
	 *            pixel coordinate from a two-dimensional double-valued
	 *            intensity image.
	 *
	 * @details
	 *            Essentially, this function is implemented as described
	 *            in
	 *            <table border="0" align="center">
	 *             <tr><td>
	 *              [MMJ+09] Maltoni, Maio, Jain, and Prabhakar (2009).
	 *              <i>Handbook of %Fingerprint Recognition</i>, 2nd ed.
	 *              Springer Publishing Company, Incorporated.
	 *             </td></tr>
	 *            </table>
	 *            and calls the functions \link gradX()\endlink and
	 *            \link gradY()\endlink
	 *
	 * @param y0
	 *            The y-value of the coordinate at which the local
	 *            orientation is estimated.
	 *
	 * @param x0
	 *            The x-value of the coordinate at which the local
	 *            orientation is estimated
	 *
	 * @param image
	 *            The double-valued intensities of the specified image;
	 *            where the intensity of the pixel at <code>(y,x)</code>
	 *            is stored at <code>image[y*n+x]</code>.
	 *
	 * @param m
	 *            The height of the specified image.
	 *
	 * @param n
	 *            The width of the specified image.
	 *
	 * @param h
	 *            Controls how many neighboring pixel's gradients are
	 *            computed for estimating the local orientation.
	 *
	 * @return
	 *            The estimated orientation.
	 *
	 * @warning
	 *            In the following cases, the behavior of this function
	 *            is undocumented:
	 *            <ul>
	 *             <li><code>m</code> is negative</li>
	 *             <li><code>n</code> is negative</li>
	 *             <li><code>h</code> is negative</li>
	 *             <li>
	 *              <code>image</code> does not contain
	 *              <code>m*n</code> valid <code>double</code> values.
	 *             </li>
	 *             <li><code>y0</code> is not a value from 0,...,m-1
	 *             <li><code>x0</code> is not a value from 0,...,n-1
	 *            </ul>
	 *
	 * @see gradX()
	 * @see gradY()
	 */
	Orientation Orientation::gradient
		( int y0 , int x0 , const double *image , int m , int n , int h ) {

		// Compute  Gxy, Gxx and Gyy as in Equation (3) in [MMJ+09]
		// on page 104.
		double gxx , gyy , gxy;
		{
			gxy = 0.0;
			gxx = 0.0;
			gyy = 0.0;

			for ( int x = x0-h ; x <= x0+h ; x++ ) {

				if ( x >= 0 && x < n ) {
					for ( int y = y0-h ; y <= y0+h ; y++ ) {

						if ( y >= 0 && y < m ) {

							double dx , dy;
							dx = Orientation::gradX(y,x,image,m,n);
							dy = Orientation::gradY(y,x,image,m,n);
							gxy += dx*dy;
							gxx += dx*dx;
							gyy += dy*dy;
						}

					}
				}

			}
		}

		// Done: Gxy, Gxx and Gyy have been computed.

		// Now, compute the orientation angle from [0,pi) and ...
		double theta = 0.5*(M_PI+atan2(gxy+gxy,gxx-gyy));

		// ... and compute the coherence which can be computed easily
		// from the Gxy, Gxx and Gyy computed above.
		double coherence = gxx-gyy;
		coherence *= coherence;
		coherence += 4.0*gxy*gxy;
		coherence = sqrt(coherence);
		coherence /= gxx+gyy;

		// Due to limited machine precision, the coherence can be
		// slightly above or below the interval [0,1].
		if ( coherence < 0.0 ) {
			coherence = 0.0;
		} else if ( coherence > 1.0 ) {
			coherence = 1.0;
		}

		// Return an orientation object.
		return Orientation(theta,coherence);
	}

	/**
	 * @brief
	 *            Creates a direction of specified angle and
	 *            coherence (also serves as the standard constructor).
	 *
	 * @details
	 *            Calls \link setDirectionAngle()\endlink and then
	 *            \link setCoherence()\endlink to create a new
	 *            \link Direction\endlink object.
	 *
	 * @param angle
	 *            Specified direction angle which should be a value
	 *            from the interval \f$[0,2\pi)\f$.
	 *
	 * @param coherence
	 *            The coherence of the direction indicating how
	 *            reliable the estimation is. Therein,
	 *            a value close to 0 indicates a very unreliable
	 *            estimation and a value close 1 a very reliable
	 *            estimation.
	 *
	 * @warning
	 *            If <code>coherence</code> is not a value from
	 *            \f$[0,1]\f$ an error message will be printed to
	 *            <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 */
	Direction::Direction( double angle , double coherence ) {
		setDirectionAngle(angle);
		setCoherence(coherence);
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @param direction
	 *            The direction of which a copy is created.
	 */
	Direction::Direction( const Direction & direction ) {
		this->angle       = direction.angle;
		this->coherence   = direction.coherence;
		this->dx          = direction.dx;
		this->dy          = direction.dy;
		this->dx2         = direction.dx2;
		this->dy2         = direction.dy2;

		this->s           = direction.s;
	}

	/**
	 * @brief
	 *            Assigns this by another direction.
	 *
	 * @param direction
	 *            The direction assigned to this object.
	 */
	void Direction::assign( const Direction & direction ) {
		this->angle       = direction.angle;
		this->coherence   = direction.coherence;
		this->dx          = direction.dx;
		this->dy          = direction.dy;
		this->dx2         = direction.dx2;
		this->dy2         = direction.dy2;
		this->s           = direction.s;
	}

	/**
	 * @brief
	 *            Assigns this direction by an orientation.
	 *
	 * @details
	 *            Calling this method is equivalent to
	 *            <pre>
	 *   assign(%Direction(orientation.getOrientationAngle(),orientation.getCoherence()))
	 *            </pre>
	 *
	 * @param orientation
	 *            The orientation being assigned to this
	 *            direction object.
	 *
	 * @attention
	 *            This method overrides
	 *            the \link Orientation::assign()\endlink method
	 *            to set a positive sign of this direction.
	 */
	void Direction::assign( const Orientation & orientation ) {
		Orientation::assign(orientation);
		this->s = 1;
	}

	/**
	 * @brief
	 *            Assignment operator.
	 *
	 * @details
	 *            Using the '='-operator calls the method
	 *            \link assign()\endlink.
	 *
	 * @param direction
	 *            The direction assigned to this object.
	 *
	 * @return
	 *            A reference to this direction (after assignment).
	 *
	 * @see assign()
	 */
	Direction & Direction::operator=( const Direction & direction ) {
		assign(direction);
		return *this;
	}

	/**
	 * @brief
	 *            Updates the direction angle.
	 *
	 * @details
	 *            The direction angle should be value between
	 *            \f$[0,2\pi)\f$. If <code>angle</code> is a value
	 *            outside of this interval, it is normalized via
	 *            the \link normalize()\endlink function.
	 *
	 *            Replacing the direction angle affects the
	 *            result of the functions
	 *            <ul>
	 *             <li>\link getDirectionAngle()\endlink</li>
	 *             <li>\link get_dx()\endlink</li>
	 *             <li>\link get_dy()\endlink</li>
	 *             <li>\link get_dx2()\endlink</li>
	 *             <li>\link get_dy2()\endlink</li>
	 *            </ul>
	 *
	 * @param angle
	 *            Specified orientation angle which should be a value
	 *            from the interval \f$[0,\pi)\f$.
	 */
	void Direction::setDirectionAngle( double angle ) {

		double tmpAngle = Direction::normalize(angle);

		if ( tmpAngle < M_PI ) {
			this->s = 1;
		} else {
			this->s = -1;
		}

		setOrientationAngle(tmpAngle);
	}

	/**
	 * @brief
	 *            Access the cosine of the direction angle.
	 *
	 * @details
	 *            The cosine of the orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            cos(\link getDirectionAngle()\endlink)
	 *
	 * @attention
	 *            This method overrides the function
	 *            \link Orientation::get_dx()\endlink and returns
	 *            \link Orientation::get_dx()\endlink if the sign
	 *            of the direction is positive (i.e., if the direction
	 *            angle is in \f$[0,\pi)\f$) and
	 *            -\link Orientation::get_dx()\endlink otherwise if
	 *            the sign is negative (i.e., if the direction
	 *            angle is in \f$[\pi,2\pi)\f$).
	 */
	double Direction::get_dx() const {
		return this->s * Orientation::get_dx();
	}

	/**
	 * @brief
	 *            Access the sine of the direction angle.
	 *
	 * @details
	 *            The sine of the orientation angle is
	 *            pre-computed and stored by this orientation object
	 *            to avoid inefficient re-computations.
	 *
	 * @return
	 *            sin(\link getDirectionAngle()\endlink)
	 *
	 * @attention
	 *            This method overrides the function
	 *            \link Orientation::get_dy()\endlink and returns
	 *            \link Orientation::get_dy()\endlink if the sign
	 *            of the direction is positive (i.e., if the direction
	 *            angle is in \f$[0,\pi)\f$) and
	 *            -\link Orientation::get_dy()\endlink otherwise if
	 *            the sign is negative (i.e., if the direction
	 *            angle is in \f$[\pi,2\pi)\f$).
	 */
	double Direction::get_dy() const {
		return this->s * Orientation::get_dy();
	}

	/**
	 * @brief
	 *            Access the angle of this direction.
	 *
	 * @return
	 *            A value from the interval \f$[0,2\pi)\f$.
	 */
	double Direction::getDirectionAngle() const {
		if ( this->s < 0 ) {
			return M_PI+getOrientationAngle();
		} else {
			return getOrientationAngle();
		}
	}

	/**
	 * @brief
	 *            Returns the direction angle in degree rounded
	 *            to an integer.
	 *
	 * @return
	 *            A value from { 0 , ... , 359 }.
	 *
	 * @attention
	 *            This method overrides the
	 *            \link Orientation::getDegree()\endlink function
	 *            since the degree of a direction can be a value
	 *            in { 0 , ... , 359 } and is not restricted to
	 *            { 0 , ... , 179 }.
	 */
	int Direction::getDegree() const {

		int degree;

		degree = (int)THIMBLE_ROUND
				(360.0/(M_PI+M_PI) * (double)getDirectionAngle() );

		if ( degree >= 360 ) {
			degree -= 360;
		}

		return degree;
	}

	/**
	 * @brief
	 *            Returns a representative of the specified direction
	 *            angle from the interval \f$[0,2\pi)\f$.
	 *
	 * @details
	 *            Essentially, the function returns the value
	 *            \f[
	 *             angle' = k\cdot 2\pi+angle
	 *            \f]
	 *            where <i>k</i> is an integer such that
	 *            \f$angle'\in[0,2\pi)\f$.
	 *
	 * @param angle
	 *            An (unnormalized) direction angle.
	 *
	 * @return
	 *            The normalized direction angle.
	 */
	double Direction::normalize( double angle ) {

		if ( angle >= 0.0 && angle < M_PI+M_PI ) {
			return angle;
		}

		double dx , dy , normalizedAngle;
		dx = cos(angle);
		dy = sin(angle);
		normalizedAngle = atan2(dy,dx);

		while ( normalizedAngle < 0.0 ) {
			normalizedAngle += M_PI+M_PI;
		}

		while ( normalizedAngle >= M_PI+M_PI) {
			normalizedAngle -= M_PI+M_PI;
		}

		return normalizedAngle;
	}

	/**
	 * @brief
	 *             Standard constructor.
	 */
	OrientedPoint::OrientedPoint() {
		this->x = 0.0;
		this->y = 0.0;
	}

	/**
	 * @brief
	 *             Copy constructor.
	 *
	 * @param p
	 *             The oriented point of which a copy is created.
	 */
	OrientedPoint::OrientedPoint( const OrientedPoint & p ) {
		this->x = p.x;
		this->y = p.y;
		this->orientation = p.orientation;
	}

	/**
	 * @brief
	 *             Assignment operator.
	 *
	 * @param p
	 *             The oriented point assigned to this object.
	 *
	 * @return
	 *             A reference to this oriented point (after assignment).
	 */
	OrientedPoint & OrientedPoint::operator=( const OrientedPoint & p ) {
		this->x = p.x;
		this->y = p.y;
		this->orientation = p.orientation;
		return *this;
	}

	/**
	 * @brief
	 *             Standard constructor.
	 */
	DirectedPoint::DirectedPoint() {
		this->x = 0.0;
		this->y = 0.0;
	}

	/**
	 * @brief
	 *             Copy constructor.
	 *
	 * @param p
	 *             The directed point of which a copy is created.
	 */
	DirectedPoint::DirectedPoint( const DirectedPoint & p ) {
		this->x = p.x;
		this->y = p.y;
		this->direction = p.direction;
	}

	/**
	 * @brief
	 *             Assignment operator.
	 *
	 * @param p
	 *             The directed point assigned to this object.
	 *
	 * @return
	 *             A reference to this object (after assignment).
	 */
	DirectedPoint & DirectedPoint::operator=( const DirectedPoint & p ) {
		this->x = p.x;
		this->y = p.y;
		this->direction = p.direction;
		return *this;
	}

	/**
	 * @brief
	 *            Prints a text representation of an orientation to the
	 *            specified output stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions such
	 *            as
	 *            <pre>
	 *             Orientation orientation;
	 *             cout << orientation << endl;
	 *            </pre>
	 *            which does nothing than printing the value
	 *            <code>orientation.getOrientationAngle()</code> in a single
	 *            line of <code>stdout</code>.
	 *
	 * @param out
	 *            The output stream to which
	 *            <code>orientation.getOrientationAngle()</code> is printed.
	 *
	 * @param orientation
	 *            The orientation object of which orientation angle is printed
	 *            to <code>out</code>.
	 *
	 * @return
	 *            The reference to <code>out</code>.
	 */
	ostream &operator<<( ostream &out , const Orientation & orientation ) {
		out << orientation.getOrientationAngle();
		return out;
	}

	/**
	 * @brief
	 *            Prints a text representation of a direction to the
	 *            specified output stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions such
	 *            as
	 *            <pre>
	 *             Direction direction;
	 *             cout << direction << endl;
	 *            </pre>
	 *            which does nothing than printing the value
	 *            <code>direction.getDirectionAngle()</code> in a single
	 *            line of <code>stdout</code>.
	 *
	 * @param out
	 *            The output stream to which
	 *            <code>direction.getDirectionAngle()</code> is printed.
	 *
	 * @param direction
	 *            The direction object of which direction angle is printed
	 *            to <code>out</code>.
	 *
	 * @return
	 *            The reference to <code>out</code>.
	 */
	ostream &operator<<( ostream &out , const Direction & direction ) {
		out << direction.getDirectionAngle();
		return out;
	}

	/**
	 * @brief
	 *            Prints a text representation of an oriented point to the
	 *            specified output stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions such
	 *            as
	 *            <pre>
	 *             OrientedPoint p;
	 *             cout << p << endl;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             cout << p.x << " " << p.y << " " << p.orientation << endl;
	 *            </pre>
	 *
	 * @param out
	 *            The output stream to which <code>p</code> is printed.
	 *
	 * @param p
	 *            The oriented point printed to <code>out</code>.
	 *
	 * @return
	 *            The reference to <code>out</code>.
	 */
	ostream &operator<<( ostream &out , const OrientedPoint & p ) {
		out << p.x << " " << p.y << " " << p.orientation;
		return out;
	}

	/**
	 * @brief
	 *            Prints a text representation of a directed point to the
	 *            specified output stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions such
	 *            as
	 *            <pre>
	 *             DirectedPoint p;
	 *             cout << p << endl;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             cout << p.x << " " << p.y << " " << p.direction << endl;
	 *            </pre>
	 *
	 * @param out
	 *            The output stream to which <code>p</code> is printed.
	 *
	 * @param p
	 *            The directed point printed to <code>out</code>.
	 *
	 * @return
	 *            The reference to <code>out</code>.
	 */
	ostream &operator<<( ostream &out , const DirectedPoint & p ) {
		out << p.x << " " << p.y << " " << p.direction;
		return out;
	}

	/**
	 * @brief
	 *            Reads an orientation from the specified input stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions
	 *            such as
	 *            <pre>
	 *             Orientation orientation;
	 *             cin >> orientation;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             double angle;
	 *             cin >> angle;
	 *             orientation.setOrientationAngle(angle);
	 *            </pre>
	 *
	 * @param in
	 *            The input stream from which the orientation angle
	 *            is read.
	 *
	 * @param orientation
	 *            The orientation to which the read orientation angle is
	 *            stored.
	 *
	 * @return
	 *            A reference to <code>in</code>.
	 */
	istream &operator>>( istream & in , Orientation & orientation ) {
		double angle;
		in >> angle;
		orientation.assign(Orientation(angle));
		return in;
	}

	/**
	 * @brief
	 *            Reads a direction from the specified input stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions
	 *            such as
	 *            <pre>
	 *             Direction direction;
	 *             cin >> direction;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             double angle;
	 *             cin >> angle;
	 *             direction.assign(Direction(angle));
	 *            </pre>
	 *
	 * @param in
	 *            The input stream from which the direction angle
	 *            is read.
	 *
	 * @param direction
	 *            The direction to which the read direction angle is
	 *            stored.
	 *
	 * @return
	 *            A reference to <code>in</code>.
	 */
	istream &operator>>( istream & in , Direction & direction ) {
		double angle;
		in >> angle;
		direction.assign(Direction(angle));
		return in;
	}

	/**
	 * @brief
	 *            Reads an oriented point from the specified input stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions
	 *            such as
	 *            <pre>
	 *             OrientedPoint p;
	 *             cin >> p;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             cin >> p.x >> p.y >> p.orientation;
	 *            </pre>
	 *
	 * @param in
	 *            The input stream from which an oriented point is read.
	 *
	 * @param p
	 *            The oriented point read from <code>in</code>.
	 *
	 * @return
	 *            The reference to <code>in</code>.
	 */
	istream &operator>>( istream & in , OrientedPoint & p ) {
		p.x = THIMBLE_NAN;
		p.y = THIMBLE_NAN;
		in >> p.x >> p.y >> p.orientation;
		return in;
	}

	/**
	 * @brief
	 *            Reads a direction point from the specified input stream.
	 *
	 * @details
	 *            The implementation of this function allows expressions
	 *            such as
	 *            <pre>
	 *             DirectedPoint p;
	 *             cin >> p;
	 *            </pre>
	 *            which is equivalent to
	 *            <pre>
	 *             cin >> p.x >> p.y >> p.direction;
	 *            </pre>
	 *
	 * @param in
	 *            The input stream from which a directed point is read.
	 *
	 * @param p
	 *            The directed point read from <code>in</code>.
	 *
	 * @return
	 *            The reference to <code>in</code>.
	 */
	istream &operator>>( istream & in , DirectedPoint & p ) {

		p.x = THIMBLE_NAN;
		p.y = THIMBLE_NAN;
		p.direction.setDirectionAngle(0.0);

		in >> p.x >> p.y >> p.direction;
		return in;
	}
}
