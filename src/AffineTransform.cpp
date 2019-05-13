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
 * @file AffineTransform.cpp
 *
 * @brief
 *            Implements functions and methods provided by
 *            'AffineTransform.h' which are related with a class for
 *            representing and computing with two-dimensional real
 *            affine transforms.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include <thimble/math/AffineTransform.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            In case this affine transform is a spatial movement,
	 *            the rotation angle is accessed.
	 *
	 * @details
	 *            see 'AffineTransform.h'
	 */
	double AffineTransform::getRotationAngle() const {

		double theta =
				acos(this->a/sqrt(this->a*this->a+this->c*this->c));

		if ( asin(this->b) < 0.0 ) {
			theta = M_PI+M_PI-theta;
		}

		return theta;
	}

	/**
	 * @brief
	 *            Rotates the affine transform by the
	 *            specified angle.
	 *
	 * @details
	 *            More precisely, the affine transform is multiplied by
	 *            the matrix
	 *            \f[
	 *              \left(
	 *               {\cos(\theta)~-\sin(\theta)}
	 *                   \atop
	 *               {\sin(\theta)~~\cos(\theta)}\right)
	 *            \f]
	 *            such that the affine transform
	 *            \f[
	 *             f:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
	 *             \mapsto
	 *             \left({a~~b}\atop{c~~d}\right)\cdot\left({x\atop y}\right)+
	 *             \left({v\atop w}\right)
	 *            \f]
	 *            is replaced by
	 *            \f[
	 *            f\leftarrow\left(
	 *               {\cos(\theta)~-\sin(\theta)}
	 *                   \atop
	 *               {\sin(\theta)~~\cos(\theta)}\right)\cdot
	 *              f\left({x\atop y}\right).
	 *            \f]
	 *
	 * @param theta
	 *            Specifies the rotation angle.
	 */
	void AffineTransform::rotate( double theta ) {

		double a , b , c , d , v , w , sint , cost;

		// Backup of the members
		a = this->a;
		b = this->b;
		c = this->c;
		d = this->d;
		v = this->v;
		w = this->w;

		// Compute the entries of the multiplication matrix.
	    sint = sin(theta);
	    cost = cos(theta);

	    // Replace the affine transform by a multiplication with
	    // the rotation matrix.

	    this->a =  cost * a - sint * c;
	    this->c =  sint * a + cost * c;

	    this->b =  cost * b - sint * d;
	    this->d =  sint * b + cost * d;

	    this->v =  cost * v - sint * w;
	    this->w =  sint * v + cost * w;
	}

	/**
	 * @brief
	 *            Rotates the affine transform by the negative of the
	 *            specified angle.
	 *
	 * @details
	 *            More precisely, the affine transform is multiplied by
	 *            the matrix
	 *            \f[
	 *              \left(
	 *               {\cos(\theta)~~\sin(\theta)}
	 *                   \atop
	 *              {-\sin(\theta)~\cos(\theta)}\right)
	 *            \f]
	 *            such that the affine transform
	 *            \f[
	 *             f:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
	 *             \mapsto
	 *             \left({a~~b}\atop{c~~d}\right)\cdot\left({x\atop y}\right)+
	 *             \left({v\atop w}\right)
	 *            \f]
	 *            is replaced by
	 *            \f[
	 *            f\leftarrow\left(
	 *               {\cos(\theta)~~\sin(\theta)}
	 *                   \atop
	 *              {-\sin(\theta)~\cos(\theta)}\right)\cdot
	 *              f\left({x\atop y}\right).
	 *            \f]
	 *
	 * @param theta
	 *            Specifies the rotation angle.
	 */
	void AffineTransform::irotate( double theta ) {

		double a , b , c , d , v , w , sint , cost;

		// Backup of the members
		a = this->a;
		b = this->b;
		c = this->c;
		d = this->d;
		v = this->v;
		w = this->w;

		// Compute the entries of the multiplication matrix.
	    sint = sin(theta);
	    cost = cos(theta);

	    // Replace the affine transform by a multiplication with
	    // the rotation matrix.

	    this->a =  cost * a + sint * c;
	    this->c = -sint * a + cost * c;

	    this->b =  cost * b + sint * d;
	    this->d = -sint * b + cost * d;

	    this->v =  cost * v + sint * w;
	    this->w = -sint * v + cost * w;
	}
}
