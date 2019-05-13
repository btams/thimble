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
 * @file AffineTransform.h
 *
 * @brief
 *            Provides a class for representing and computing with
 *            two-dimensional real affine transforms.
 *
 * @author Benjamin Tams
 *
 * @see thimble::AffineTransform
 */

#ifndef THIMBLE_AFFINETRANSFORM_H_
#define THIMBLE_AFFINETRANSFORM_H_

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of this class represent two-dimensional affine
	 *            transforms.
	 *
	 * @details
	 *            Instances members of this class are six real values
	 *            <i>a</i>, <i>b</i>, <i>c</i>, <i>d</i>, <i>v</i>, and
	 *            <i>w</i> represented by <code>double</code>s. Such an
	 *            instance represent the affine transform
	 *            \f[
	 *             f:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
	 *             \mapsto
	 *             \left({a~~b}\atop{c~~d}\right)\cdot\left({x\atop y}\right)+
	 *             \left({v\atop w}\right).
	 *            \f]
	 *            The identity transform
	 *            \f[
	 *             id:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
	 *             \mapsto
	 *             \left({1~~0}\atop{0~~1}\right)\cdot\left({x\atop y}\right)+
	 *             \left({0\atop 0}\right)
	 *            \f]
	 *            can be created through the standard constructor, e.g., via
	 *            <pre>
	 *             AffineTransform id();
	 *            </pre>
	 *            or
	 *            <pre>
	 *             AffineTransform id;
	 *            </pre>
	 *            Modifying an affine transform can be realized by the user
	 *            who has public access to the class members. Furthermore, an
	 *            affine transform can be rotated by an angle \f$\theta\f$ by
	 *            means of multiplication with the matrix
	 *            \f[
	 *             \left(
	 *              {\cos(\theta)~~\sin(\theta)}
	 *                   \atop
	 *              {-\sin(\theta)~\cos(\theta)}\right)
	 *            \f]
	 *            can be performed through
	 *            <pre>
	 *             AffineTransform f = ...;
	 *             f.rotate(theta);
	 *            </pre>
	 *            Very special affine transforms are those that are spatial
	 *            movements, i.e. affine transforms of the shape
	 *            \f[
	 *             f:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
	 *             \mapsto
	 *             \left(
	 *              {\cos(\theta)~~\sin(\theta)}
	 *               \atop{-\sin(\theta)~~\cos(\theta)}
	 *              \right)\cdot\left({x\atop y}\right)+
	 *             \left({v\atop w}\right).
	 *            \f]
	 *            Their rotation angle can be requested via
	 *            <pre>
	 *             AffineTransform f = ...;
	 *             double theta = f.getRotationAngle();
	 *            </pre>
	 *            in the case <i>f</i> does not represent a spatial movement
	 *            the geometric interpretation of the result is undocumented.
	 */
	class THIMBLE_DLL AffineTransform {

	public:

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double a;

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double b;

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double c;

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double d;

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double v;

		/**
		 * @brief
		 *           Public member variable used to encode the affine
		 *           transform.
		 *
		 * @details
		 *           For details see the detailed description of
		 *           \link AffineTransform\endlink.
		 *
		 * @see AffineTransform
		 */
		double w;

		/**
		 * @brief
		 *           Overwrites the affine transform by the identity
		 *           transform.
		 *
		 * @details
		 *           The affine transform will be assigned by
		 *           \f[
		 *             id:{\bf R}^2\rightarrow{\bf R}^2,
		 *             \left({x\atop y}\right)
		 *             \mapsto
		 *             \left({1~~0}\atop{0~~1}\right)\cdot
		 *             \left({x\atop y}\right)+
		 *             \left({0\atop 0}\right).
		 *           \f]
		 */
		inline void setIdentity() {
			this->a = 1.0;
			this->b = 0.0;
			this->c = 0.0;
			this->d = 1.0;
			this->v = 0.0;
			this->w = 0.0;
		}

		/**
		 * @brief
		 *            Assigns the affine transform by the affine transform
		 *            that is encoded by <i>f</i>.
		 *
		 * @param f
		 *            The affine transform that this transform is assigned
		 *            with.
		 */
		inline void assign( const AffineTransform & f ) {
			this->a = f.a;
			this->b = f.b;
			this->c = f.c;
			this->d = f.d;
			this->v = f.v;
			this->w = f.w;
		}

		/**
		 * @brief
		 *            Standard constructor creating the identity transform.
		 *
		 * @details
		 *            The constructor calls \link setIdentity()\endlink.
		 *
		 * @see setIdentity()
		 */
		inline AffineTransform() {
			setIdentity();
		}

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @details
		 *            The constructor calls
		 *            \link assign(const AffineTransform&)\endlink to create
		 *            a copy of the instance encoded by <i>f</i>.
		 *
		 * @param f
		 *            The affine transform of which a copy is created by the
		 *            constructor.
		 */
		inline AffineTransform( const AffineTransform & f ) {
			assign(f);
		}

		/**
		 * @brief
		 *            Assignment operator.
		 *
		 * @details
		 *            Statements like
		 *            <pre>
		 *             AffineTransform f , g;
		 *             ...
		 *             f = g;
		 *            </pre>
		 *            cause <i>f</i> to become a copy of <i>g</i>. The
		 *            implementation of the assignment operator wraps around
		 *            \link assign(const AffineTransform&)\endlink.
		 *
		 * @param f
		 *            The affine transform of which this instance becomes a
		 *            copy of.
		 *
		 * @return
		 *            A reference to this affine transform.
		 */
		inline AffineTransform &operator=( const AffineTransform & f ) {
			assign(f);
			return *this;
		}

		/**
		 * @brief
		 *            In case this affine transform is a spatial movement,
		 *            the rotation angle is accessed.
		 *
		 * @details
		 *            A spatial movement is an affine transform of the shape
		 *            \f[
		 *             f:{\bf R}^2\rightarrow{\bf R}^2, \left({x\atop y}\right)
		 *             \mapsto
		 *             \left(
		 *              {\cos(\theta)~~\sin(\theta)}
		 *               \atop{-\sin(\theta)~~\cos(\theta)}
		 *              \right)\cdot\left({x\atop y}\right)+
		 *             \left({v\atop w}\right).
		 *            \f]
		 *            If this affine transform is of the above the shape with
		 *            \f$\theta\in[0,2\pi)\f$ the function returns
		 *            \f$\theta\f$.
		 *
		 * @return
		 *            The rotation angle of the spatial movement.
		 *
		 * @warning
		 *            If the affine transform is not of the shape of a spatial
		 *            movement, the result and behavior of the function is
		 *            undocumented.
		 */
		double getRotationAngle() const;

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
		void rotate( double theta );

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
		void irotate( double theta );

		/**
		 * @brief
		 *            Computes the evaluation of the affine transform at a
		 *            <code>double</code> pair.
		 *
		 * @details
		 *            If the affine transform is denoted by <i>f</i> the
		 *            values <i>fx</i> and <i>fy</i> change such that
		 *            \f[
		 *            \left({fx\atop fy}\right)\leftarrow
		 *            f\left({x\atop y}\right).
		 *            \f]
		 *
		 * @param fx
		 *            Abscissa coordinate of the evaluation.
		 *
		 * @param fy
		 *            Ordinate coordinate of the evaluation.
		 *
		 * @param x
		 *            Abscissa coordinate of the pair at where we evaluate.
		 *
		 * @param y
		 *            Orindate coordinate of the pair at where we evaluate.
		 */
		inline void eval
		( double & fx , double & fy , double x , double y ) const {
			fx = this->a * x + this->b * y + this->v;
			fy = this->c * x + this->d * y + this->w;
		}
	};
}


#endif /* THIMBLE_AFFINETRANSFORM_H_ */
