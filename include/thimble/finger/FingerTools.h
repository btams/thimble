/*
 *  THIMBLE --- A Research Library Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file FingerTools.h
 *
 * @brief
 *            Provides utility functions related with fingerprints.
 *
 * @details
 *            see the \link thimble::FingerTools\endlink class.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_FINGERTOOLS_H_
#define THIMBLE_FINGERTOOLS_H_

#include <string>

#include <thimble/dllcompat.h>
#include <thimble/math/AffineTransform.h>
#include <thimble/image/Orientation.h>
#include <thimble/finger/MinutiaeRecord.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Class collecting static utility functions related to
	 *            fingerprints.
	 *
	 * @details
	 *            In the current version, there are only functions provided with
	 *            the purpose to align minutiae templates. Its main function
	 *            is
	 *            <code>
	 *             FingerTools::
	 *              align(const MinutiaeView&,const MinutiaeView&,int,double)
	 *            </code>.
	 *            If two minutiae template
	 *            <pre>
	 *             MinutiaeView v , w;
	 *            </pre>
	 *            are given, then the spatial movement which aligns
	 *            <code>b</code> to <code>a</code> can be obtained via
	 *            <pre>
	 *             AffineTransform f = FingerTools::align(v,w);
	 *            </pre>
	 *            The result is an affine transform such that
	 *            <pre>
	 *             MinutiaeView u = FingerTools::eval(f,w);
	 *            </pre>
	 *            is such that <code>u</code> is the aligned version
	 *            of <code>w</code> such that its minutiae well approximate
	 *            the minutiae of <code>v</code>.
	 */
	class THIMBLE_DLL FingerTools {

	public:

		/**
		 * @brief
		 *            Computes the spatial movement as an affine transform
		 *            that maps the minutia <code>b</code> to the minutia
		 *            <code>a</code>.
		 *
		 * @details
		 *            More precisely, write \f$a=(x_a,y_a,\theta_a)\f$ and
		 *            \f$b=(x_b,y_b,\theta_b)\f$ are the input minutiae
		 *            at \f$(x_a,y_a)\f$ and \f$(x_b,y_b)\f$, respectively,
		 *            where <i>a</i> and <i>b</i> are of angle \f$\theta_a\f$
		 *            and \f$\theta_b\f$, respectively. Then the result
		 *            is the affine transform
		 *            \f[
		 *             f:\left({x\atop y}\right)\mapsto
		 *               \left(
		 *                {\cos(\theta_b-\theta_a)~~\sin(\theta_b-\theta_a)}
		 *                 \atop
		 *                {-\sin(\theta_b-\theta_a)~~\cos(\theta_b-\theta_a)}
		 *               \right)
		 *               \left(
		 *                {x-x_b}
		 *                 \atop
		 *                {y-y_b}
		 *               \right)+
		 *               \left(
		 *                {x_a}\atop{y_a}
		 *               \right)
		 *            \f]
		 *            such that
		 *            \f[
		 *             f(x_b,y_b)=(x_a,y_a).
		 *            \f]
		 *            Moreover, the spatial movement is a rotation around the
		 *            origin \f$(x_b,y_b)\f$ of the angle
		 *            \f$\theta_b-\theta_a\f$ such that the angle of <i>b</i>
		 *            is mapped to the angle of <i>a</i> this way.
		 *
		 * @param a
		 *            Minutia to where <i>b</i> is aligned.
		 *
		 * @param b
		 *            Minutia that is aligned to <i>a</i>.
		 *
		 * @return
		 *            Spatial movement that aligns <i>b</i> to <i>a</i>.
		 */
		static AffineTransform align
		( const Minutia & a , const Minutia & b );

		/**
		 * @brief
		 *            Returns the spatial movement of the minutia <i>a</i>
		 *            using the spatial movement <i>f</i>.
		 *
		 * @details
		 *            More precisely, write
		 *            \f[
		 *             f:\left({x\atop y}\right)\mapsto
		 *             \left(
		 *              {\cos(\theta)~~\sin(\theta)}
		 *               \atop
		 *              {-\sin(\theta)~~\cos(\theta)}
		 *             \right)\cdot\left(
		 *              {x\atop y}
		 *             \right)+\left(
		 *              {x_s\atop y_s}
		 *             \right).
		 *            \f]
		 *            Then, if <i>a</i> is at \f$(x_a,y_a)\f$ and of angle
		 *            \f$\theta_a\f$ then the result will be the minutia
		 *            at \f$f(x_a,y_a)\f$ that is of angle
		 *            \f$\theta_a-\theta\f$ modulo \f$2\pi\f$.
		 *            <br><br>
		 *            The type of the returned minutia as well as its quality
		 *            are adopted from <i>a</i>.
		 *
		 * @param f
		 *            The spatial movement given as an affine transform.
		 *
		 * @param a
		 *            The minutia that is moved by <i>f</i>.
		 *
		 * @return
		 *            The minutia <i>a</i> moved by <i>f</i>.
		 *
		 * @warning
		 *            If <i>f</i> does not correspond to a spatial movement,
		 *            i.e. a rotation around an origin, the geometric
		 *            interpretation of the result is undocumented.
		 */
		static Minutia eval
		( const AffineTransform & f , const Minutia & a );

		/**
		 * @brief
		 *            Returns the movement of the given minutiae
		 *            template by the specified spatial movement.
		 *
		 * @details
		 *            The function is a wrapper around
		 *            \link
		 *             FingerTools::eval(const AffineTransform&,const Minutia&)
		 *            \endlink. The <i>i</i>th minutia of the result will
		 *            be the <i>i</i>th minutia of the input moved by
		 *            <i>f</i>.
		 *            <br><br>
		 *         	  <b>Note:</b> for the result we adopt the finger
		 *         	  position, impression type, and finger quality of the
		 *         	  input template.
		 *
		 * @param f
		 *            The spatial movement given as an affine transform.
		 *
		 * @param v
		 *            The minutiae template that is moved by <i>f</i>.
		 *
		 * @return
		 *            The minutiae template <i>v</i> moved by <i>f</i>.
		 */
		static MinutiaeView eval
		( const AffineTransform & f , const MinutiaeView & v );

		/**
		 * @brief
		 *            Returns a measure of dissimilarity between the minutiae
		 *            <i>a</i> and <i>b</i>.
		 *
		 * @details
		 *            More precisely, write \f$a=(x_a,y_a,\theta_a)\f$ and
		 *            \f$b=(x_b,y_b,\theta_b)\f$ for the input minutiae where
		 *            \f$(x_a,y_b)\f$ and \f$(x_b,y_b)\f$ are their respective
		 *            positions and \f$\theta_a\f$ and \f$\theta_b\f$ its
		 *            angles which range in \f$[0,2\pi)\f$. The result of
		 *            the function will be
		 *            \f[
		 *             \|(x_a,y_a)-(x_b,y_b)\|_2+
		 *                 angleWeight\cdot
		 *              \min\{|\theta_a-\theta_b|,|2\pi-\theta_a+\theta_b|\}.
		 *            \f]
		 *            Thus, <i>angleWeight</i> controls how significant the
		 *            minutiae angles are taken into account for measuring
		 *            the minutiae's dissimilarity.
		 *            <br><br>
		 *            In case <i>angleWeight=0</i> the function computes the
		 *            <i>Euclidean distance</i> between the coordinates of
		 *            <i>a</i> and <i>b</i>. Furthermore,
		 *            if <i>angleWeight>0</i> then the result will be 0.0 if
		 *            and only if <i>a</i> and <i>b</i> are equal (by means of
		 *            equal coordinates and angles).
		 *
		 * @param a
		 *            First minutia.
		 *
		 * @param b
		 *            Second minutia.
		 *
		 * @param angleWeight
		 *            Controls how significant the minutiae angles are taken
		 *            into account.
		 *
		 * @return
		 *            &quot;Distance&quot; between <i>a</i> and <i>b</i>.
		 */
		static double dist
		( const Minutia & a , const Minutia & b ,
		  double angleWeight = 11.459 );

		/**
		 * @brief
		 *            Computes a dissimilarity measure between two minutiae
		 *            templates.
		 *
		 * @details
		 *            The result is the sum of the distance between the
		 *            <i>n</i> minutiae pairs from <i>v</i> and <i>w</i> with
		 *            minimal &quot;distance&quot; (by means of \link
		 *             FingerTools::dist(const Minutia&,const Minutia&,double)
		 *            \endlink).
		 *            <br><br>
		 *            More precisely, let <i>v[0],...,v[k-1]</i> and
		 *            <i>w[0],...,w[l-1]</i> be the minutiae of <i>v</i> and
		 *            <i>w</i>, respectively, such that
		 *            <i>
		 *             dist(v[i],w[i],angleWeight) <=
		 *             dist(v[j],w[j],angleWeight)
		 *            </i>
		 *            if <i>i < j</i>. Then the result is
		 *            <pre>
		 *    dist(v[0],w[0],angleWeight)+...+dist(v[m-1],w[m-1],angleWeight)
		 *            </pre>
		 *            where <i>m=min(k,l,n)</i> if <i>n>0</i>; otherwise, if
		 *            <i>n<=0</i> then <i>m=min(k,l)</i>.
		 *
		 * @param v
		 *            First minutiae template.
		 *
		 * @param w
		 *            Second minutiae template.
		 *
		 * @param n
		 *            Controls the Number of minutiae correspondences taken
		 *            into account.
		 *
		 * @param angleWeight
		 *            Controls how significant the minutiae angles are taken
		 *            into account.
		 *
		 * @return
		 *            Dissimilarity between <i>v</i> and <i>w</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated the function
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static double dist
		( const MinutiaeView & v , const MinutiaeView & w ,
		  int n = 5 , double angleWeight = 11.459 );

		/**
		 * @brief
		 *            Determines a spatial movement as an affine transform
		 *            that aligns the minutiae in <code>w</code> to the
		 *            minutiae in <code>v</code>.
		 *
		 * @details
		 *            Among multiple candidates for spatial movements <i>f</i>,
		 *            the function determines the one that minimizes the
		 *            dissimilarity between <i>v</i> and <i>f(w)</i>. Thereby
		 *            the dissimilarity between <i>v</i> and <i>f(w)</i> is
		 *            given by \link
		 * FingerTools::dist(const MinutiaeView&,const MinutiaeView&,int,double)
		 *            \endlink.
		 *            <br><br>
		 *            The candidate alignments are achieved from all minutiae
		 *            correspondences between <i>v</i> and <i>w</i> which
		 *            is of number <i>k*l</i> in total, where
		 *            <i>k=v.getMinutiaeCount()</i> and
		 *            <i>l=w.getMinutiaeCount()</i>: Let <i>a</i> be the
		 *            <i>i</i>th minutia of <i>v</i> and <i>b</i> be the
		 *            <i>j</i>th minutia of <i>w</i>. Then the corresponding
		 *            candidate alignment is given by the movement that maps
		 *            <i>b</i> to <i>a</i>, i.e. by
		 *            <pre>
		 *           AffineTransform f = FingerTools::align(a,b,n,angleWeight)
		 *            </pre>
		 *
		 * @param v
		 *            First minutiae template.
		 *
		 * @param w
		 *            Second minutiae template.
		 *
		 * @param n
		 *            Controls the number of minutiae correspondences that
		 *            are taken into account.
		 *
		 * @param angleWeight
		 *            Controls how significant the minutiae angles are taken
		 *            into account.
		 *
		 * @return
		 *            Spatial movement that aligns the minutiae template
		 *            <i>w</i> to <i>v</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated the function
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static AffineTransform align
		( const MinutiaeView & v , const MinutiaeView & w ,
		  int n = 5 , double angleWeight = 11.459 );

		/**
		 * @brief
		 *            Computes the representation of a minutiae view w.r.t.
		 *            a directed reference point.
		 *
		 * @details
		 *            The function may be, in particular, useful for the problem
		 *            of absolute fingerprint pre-alignment in which the
		 *            minutiae are represented w.r.t. an intrinsic coordinate
		 *            system that may be given by directed reference point.
		 *            For a mechanism to estimate an intrinsic directed reference
		 *            point from a fingerprint, we refer to the
		 *            method \link thimble::Fingerprint::getDirectedReferencePoint()\endlink
		 *
		 * @param v
		 *            The unaligned view.
		 *
		 * @param x
		 *            The x-coordinate of the directed reference point.
		 *
		 * @param y
		 *            The y-coordinate of the directed reference point.
		 *
		 * @param direction
		 *            The direction of the directed reference point.
		 *
		 * @return
		 *            A representation of <code>v</code> w.r.t. the coordinate
		 *            system of origin <i>(x,y)</i> and whose abscissa is of
		 *            direction <i>(cos(direction),sin(direction))</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be provided to perform
		 *            the operations caused by the function, an error message
		 *            is printed to <code>stderr</code> and the program exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static MinutiaeView prealign
		( const MinutiaeView & v , double x , double y , double direction );

		/**
		 * @brief
		 *            Computes the representation of a minutiae view w.r.t.
		 *            a directed point defining the origin of an intrinsic
		 *            coordinate system.
		 *
		 * @details
		 *            Calling this function is equivalent to
		 *            <pre>
		 *             prealign(v,p.x,p.y,p.direction.getDirectionAngle())
		 *            </pre>
		 *
		 * @param v
		 *            The unaligned view.
		 *
		 * @param p
		 *            The directed reference point.
		 *
		 * @return
		 *            A representation of <code>v</code> w.r.t. the coordinate
		 *            system of origin <i>(p.x,p.y)</i> and whose abscissa is of
		 *            direction
		 *            <i>(cos(p.direction.getDirectionAngle()),
		 *               sin(p.direction.getDirectionAngle()))</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be provided to perform
		 *            the operations caused by the function, an error message
		 *            is printed to <code>stderr</code> and the program exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static MinutiaeView prealign
		( const MinutiaeView & v , const DirectedPoint & p );

		/**
		 * @brief
		 *            Computes the representation of a minutia w.r.t.
		 *            a directed reference point.
		 *
		 * @details
		 *            The function may be, in particular, useful for the problem
		 *            of absolute fingerprint pre-alignment in which
		 *            minutiae are represented w.r.t. an intrinsic coordinate
		 *            system that may be given by directed reference point.
		 *            For a mechanism to estimate an intrinsic directed reference
		 *            point from a fingerprint, we refer to
		 *            the \link Fingerprint::getDirectedReferencePoint()\endlink
		 *            method.
		 *
		 * @param m
		 *            The unaligned minutia.
		 *
		 * @param x
		 *            The x-coordinate of the directed reference point.
		 *
		 * @param y
		 *            The y-coordinate of the directed reference point.
		 *
		 * @param direction
		 *            The direction of the directed reference point.
		 *
		 * @return
		 *            A representation of <code>m</code> w.r.t. the coordinate
		 *            system of origin <i>(x,y)</i> and whose abscissa is of
		 *            direction <i>(cos(direction),sin(direction))</i>.
		 */
		static Minutia prealign( const Minutia & m , double x , double y , double direction );

		/**
		 * @brief
		 *            Computes the representation of a minutia w.r.t.
		 *            a directed point defining the origin of an intrinsic
		 *            coordinate system.
		 *
		 * @details
		 *            Calling this function is equivalent to
		 *            <pre>
		 *             prealign(m,p.x,p.y,p.direction.getDirectionAngle())
		 *            </pre>
		 *
		 * @param m
		 *            The unaligned minutia.
		 *
		 * @param p
		 *            The directed reference point.
		 *
		 * @return
		 *            A representation of <code>m</code> w.r.t. the coordinate
		 *            system of origin <i>(p.x,p.y)</i> and whose abscissa is of
		 *            direction
		 *            <i>(cos(p.direction.getDirectionAngle()),
		 *               sin(p.direction.getDirectionAngle()))</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be provided to perform
		 *            the operations caused by the function, an error message
		 *            is printed to <code>stderr</code> and the program exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static Minutia prealign( const Minutia & m , const DirectedPoint & p );
	};
}


#endif /* THIMBLE_FINGERTOOLS_H_ */
