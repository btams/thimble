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
 * @file Orientation.h
 *
 * @brief
 *            Provides functionalities for computing orientations from
 *            images and representation of orientations, directions and
 *            two-dimensional point coordinates attached with orientations.
 *            and directions.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_ORIENTATION_H_
#define THIMBLE_ORIENTATION_H_

#include <ostream>
#include <istream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of this class represent orientations.
	 *
	 * @details
	 *            An orientation can be represented by an angle in the
	 *            interval \f$[0,\pi)\f$.
	 *
	 * @attention
	 *            Note that an angle in \f$[0,2\pi)\f$ describes a
	 *            direction. The orientation of an oriented pattern
	 *            however is described by an angle from \f$[0,\pi)\f$
	 *            since the angles \f$\phi\f$ and \f$\phi+\pi\f$ both
	 *            define the same (undirected) orientation.
	 */
	class THIMBLE_DLL Orientation {

	protected:

		/**
		 * @brief
		 *            The orientation angle which is a value
		 *            from the interval \f$[0,\pi)\f$.
		 *
		 * @see getOrientationAngle()
		 * @see setOrientationAngle()
		 */
		double angle;

		/**
		 * @brief
		 *            The coherence of this orientation.
		 *
		 * @details
		 *            The coherence is a value from \f$[0,1]\f$ and describes
		 *            the <em>parallelity</em> of the neighboring texture of
		 *            an oriented pattern which can be considered as a
		 *            confidence value for an orientation estimation.
		 *
		 * @see getCoherence()
		 * @see setCoherence()
		 */
		double coherence;

		/**
		 * @brief
		 *            The cosine of the orientation angle.
		 *
		 * @see get_dx()
		 */
		double dx;

		/**
		 * @brief
		 *            The sine of the orientation angle.
		 *
		 * @see get_dy()
		 */
		double dy;

		/**
		 * @brief
		 *            The cosine of the doubled orientation angle.
		 *
		 * @details
		 *            The cosine and sine of the doubled orientation angle
		 *            is useful when, for example, the difference between
		 *            angles are sought to be computed.
		 *
		 * @see get_dx2()
		 */
		double dx2;

		/**
		 * @brief
		 *            The sine of the doubled orientation angle.
		 *
		 * @details
		 *            The cosine and sine of the doubled orientation angle
		 *            is useful when, for example, the difference between
		 *            angles are sought to be computed.
		 *
		 * @see get_dy2()
		 */
		double dy2;

	public:

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
		Orientation( double angle = 0.0 , double coherence = 0.0 );

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @param orientation
		 *            The orientation of which a copy is created.
		 */
		Orientation( const Orientation & orientation );

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
		void setOrientationAngle( double angle );

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
		void setCoherence( double coherence );

		/**
		 * @brief
		 *            Assigns this by another orientation.
		 *
		 * @param orientation
		 *            The orientation assigned to this object.
		 */
		void assign( const Orientation & orientation );

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
		Orientation & operator=( const Orientation & orientation );

		/**
		 * @brief
		 *            Access the angle of this orientation.
		 *
		 * @return
		 *            A value from the interval \f$[0,\pi)\f$.
		 */
		double getOrientationAngle() const;

		/**
		 * @brief
		 *            Access the coherence of this orientation.
		 *
		 * @return
		 *            A value from the interval \f$[0,1]\f$.
		 */
		double getCoherence() const;

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
		double get_dx() const;

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
		double get_dy() const;

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
		double get_dx2() const;

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
		double get_dy2() const;

		/**
		 * @brief
		 *            Returns the orientation angle in degree rounded
		 *            to an integer.
		 *
		 * @return
		 *            A value from { 0 , ... , 179 }.
		 */
		int getDegree() const;

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
		static double normalize( double angle );

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
		static double gradX
			( int y0 , int x0 , const double *image , int m , int n );

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
		static double gradY
			( int y0 , int x0 , const double *image , int m , int n );

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
		static Orientation gradient
			( int y0 , int x0 , const double *image , int m , int n , int h );
	};

	/**
	 * @brief
	 *            Instances of this class represent directions.
	 *
	 * @details
	 *            A direction can be defined as an orientation attached with
	 *            a sign that indicates whether, given an orientation
	 *            angle \f$\varphi\in[0,\pi)\f$, the direction
	 *            equals \f$\varphi\f$ or \f$\varphi+\pi\f$.
	 */
	class THIMBLE_DLL Direction : public Orientation {

	private:

		/**
		 * @brief
		 *            Sign of the orientation that makes it a direction.
		 */
		int s;

	public:

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
		Direction( double angle = 0.0 , double coherence = 0.0 );

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @param direction
		 *            The direction of which a copy is created.
		 */
		Direction( const Direction & direction );

		/**
		 * @brief
		 *            Assigns this by another direction.
		 *
		 * @param direction
		 *            The direction assigned to this object.
		 */
		void assign( const Direction & direction );

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
		void assign( const Orientation & orientation );

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
		Direction &operator=( const Direction & direction );

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
		void setDirectionAngle( double angle );

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
		double get_dx() const;

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
		double get_dy() const;

		/**
		 * @brief
		 *            Access the angle of this direction.
		 *
		 * @return
		 *            A value from the interval \f$[0,2\pi)\f$.
		 */
		double getDirectionAngle() const;

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
		int getDegree() const;

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
		static double normalize( double angle );
	};

	/**
	 * @brief
	 *             Instances of this class represent two-dimensional
	 *             pixel coordinates attached with an orientation.
	 */
	class THIMBLE_DLL OrientedPoint {

	public:

		/**
		 * @brief
		 *             %Orientation attached to this point.
		 */
		Orientation orientation;

		/**
		 * @brief
		 *             X-coordinate.
		 */
		double x;

		/**
		 * @brief
		 *             Y-coordinate.
		 */
		double y;

		/**
		 * @brief
		 *             Standard constructor.
		 */
		OrientedPoint();

		/**
		 * @brief
		 *             Copy constructor.
		 *
		 * @param p
		 *             The oriented point of which a copy is created.
		 */
		OrientedPoint( const OrientedPoint & p );

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
		OrientedPoint & operator=( const OrientedPoint & p );
	};

	/**
	 * @brief
	 *             Instances of this class represent two-dimensional
	 *             pixel coordinates attached with a direction.
	 */
	class THIMBLE_DLL DirectedPoint {
	public:
		/**
		 * @brief
		 *             %Direction attached to this point.
		 */
		Direction direction;

		/**
		 * @brief
		 *             X-coordinate.
		 */
		double x;

		/**
		 * @brief
		 *             Y-coordinate.
		 */
		double y;

		/**
		 * @brief
		 *             Standard constructor.
		 */
		DirectedPoint();

		/**
		 * @brief
		 *             Copy constructor.
		 *
		 * @param p
		 *             The directed point of which a copy is created.
		 */
		DirectedPoint( const DirectedPoint & p );

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
		DirectedPoint & operator=( const DirectedPoint & p );
	};

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
	THIMBLE_DLL std::ostream &operator<<
			( std::ostream & out , const Orientation & orientation );

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
	THIMBLE_DLL std::ostream &operator<<
			( std::ostream & out , const Direction & direction );

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
	THIMBLE_DLL std::ostream &operator<<
			( std::ostream & out , const OrientedPoint & p );

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
	THIMBLE_DLL std::ostream &operator<<
			( std::ostream & out , const DirectedPoint & p );

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
	 *             orientation.assign(Orientation(angle));
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
	THIMBLE_DLL std::istream &operator>>
		( std::istream & in , Orientation & orientation );

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
	THIMBLE_DLL std::istream &operator>>
		( std::istream & in , Direction & direction );

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
	THIMBLE_DLL std::istream &operator>>( std::istream & in , OrientedPoint & p );

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
	THIMBLE_DLL std::istream &operator>>( std::istream & in , DirectedPoint & p );
}

#endif /* THIMBLE_ORIENTATION_H_ */

