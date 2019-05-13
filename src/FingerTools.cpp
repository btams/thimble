/*
 *  THIMBLE --- A Research Library for Development and Analysis of
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
 * @file FingerTools.cpp
 *
 * @brief
 *            Implements the functionalities from 'FingerTools.h' which
 *            provides utility functions related with fingerprints.
 *
 * @details
 *            see 'FingerTools.h'
 *
 * @author Benjamin Tams
 */
#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <algorithm>

#include <thimble/math/AffineTransform.h>
#include <thimble/math/RigidTransform.h>
#include <thimble/image/Orientation.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/finger/FingerTools.h>

using namespace std;

namespace thimble {

	/**
	 * @brief
	 *            Computes the spatial movement as an affine transform
	 *            that maps the minutia <code>b</code> to the minutia
	 *            <code>a</code>.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	AffineTransform FingerTools::align
	( const Minutia & a , const Minutia & b ) {

		// The alignment initialized as the identity map
		AffineTransform f;

		// Let the rotation matrix be
		//    cos(theta) sin(theta)
		//   -sin(theta) cos(theta)
		// where
		//    theta = b.getAngle()-a.getAngle()
		// which corresponds to the rotation around
		// the origin by the angle theta
		f.irotate(b.getAngle()-a.getAngle());

		// Shift the affine transform as a rotation around
		// the origin (b.getX(),b.getY())
		double x , y;
		f.eval(x,y,b.getX(),b.getY());
		f.v = a.getX()-x;
		f.w = a.getY()-y;

		return f;
	}

	/**
	 * @brief
	 *            Returns the spatial movement of the minutia <i>a</i>
	 *            using the spatial movement <i>f</i>.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	Minutia FingerTools::eval
	( const AffineTransform & f , const Minutia & a ) {

		double x , y , theta;

		// Move the position of the input minutia.
		f.eval(x,y,a.getX(),a.getY());

		// Determine rotation angle.
		theta = a.getAngle()+f.getRotationAngle();

		// Return minutia with parameters corresponding to the moved
		// input minutia; note that the type and the quality of the
		// input is adopted.
		return Minutia(x,y,theta,a.getType(),a.getQuality());
	}

	/**
	 * @brief
	 *            Returns the movement of the given minutiae
	 *            template by the specified spatial movement.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	 MinutiaeView FingerTools::eval
	( const AffineTransform & f , const MinutiaeView & v ) {

		 // First, make a copy of the input template to reserve
		 // memory and to adopt the finger position, impression
		 // type, and finger quality.
		 MinutiaeView w = v;

		 // Remove all minutiae
		 w.removeAllMinutiae();

		 // Successively append the moved minutiae of the input
		 // to the output
		 for ( int i = 0 ; i < v.getMinutiaeCount() ; i++ ) {
			 w.addMinutia(eval(f,v.getMinutia(i)));
		 }

		 return w;
	}

		/**
		 * @brief
		 *            Returns a measure of dissimilarity between the minutiae
		 *            <i>a</i> and <i>b</i>.
		 *
		 * @details
		 *            see 'FingerTools.h'
		 */
	double FingerTools::dist
	( const Minutia & a , const Minutia & b , double angleWeight ) {

		// Access the minutiae's coordinates
		double x0 , y0 , x1 , y1;
		x0 = a.getX();
		y0 = a.getY();
		x1 = b.getX();
		y1 = b.getY();

		// Difference vector between minutiae coordinate
		double dx , dy;
		dx = x0-x1;
		dy = y0-y1;

		// Compute "distance" of the angle on the unity circle.
		double dtheta = a.getAngle()-b.getAngle();
		dtheta = std::min(std::abs(dtheta),std::abs(M_PI+M_PI-dtheta));

		// Combine to the result and return.
		return sqrt(dx*dx+dy*dy)+angleWeight*dtheta;
	}


	/**
	 * @brief
	 *            Computes the dissimilarity between 'v' and 'w'.
	 *
	 * @details
	 *            This function is used by the functions
	 *            \link
	 *  FingerTools::dist(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink
	 *            and
	 *            \link
	 *  FingerTools:align(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink.
	 *            Unlike the user function
	 *            \link
	 *  FingerTools::dist(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink
	 *            it contains an additional parameter <code>distances</code>
	 *            which is required to be an array that can hold at least
	 *            <code>m</code> <code>double</code> values. This is to
	 *            avoid frequent reallocations by the
	 *            \link
	 *  FingerTools:align(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink
	 *            which is imperformant. Because the function can be used to
	 *            implement
	 *            \link
	 *  FingerTools::dist(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink
	 *            it is also used there to follow the minimal principle in
	 *            programming.
	 *            <br><br>
	 *            The result is the sum of the distance between the
     *            <i>m</i> minutiae pairs from <i>v</i> and <i>w</i> with
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
	 *    dist(v[0],w[0],angleWeight)+...+dist(v[m-1],w[m-1],angleWeight).
	 *            </pre>
	 *            Note, that, contrary to
	 *            \link
	 *  FingerTools::dist(const MinutiaeView&v,const MinutiaeView&,int,double)
	 *            \endlink
	 *            here the template <i>v</i> and <i>w</i> are required to
	 *            contain at least <i>m</i> minutiae where <i>m</i> is
	 *            required to be strictly larger than 0.
	 *
	 * @param v
	 *            First minutiae template.
	 *
	 * @param w
	 *            Second minutiae template.
	 *
	 * @param m
	 *            Number of minutiae correspondences between <i>v</i> and
	 *            <i>w</i> that are taken into account.
	 *
	 * @param angleWeight
	 *            Controls how significant the minutiae angles are taken
	 *            into account.
	 *
	 * @param distances
	 *            Pre-allocated array that can hold <i>m</i>
	 *            <code>double</code> values for internal use.
	 *
	 * @return
	 *            Dissimilarity between <i>v</i> and <i>w</i>.
	 *
	 * @warning
	 *            If the prerequisites found in the details are not
	 *            satisfied, the result and the behavior of the function is
	 *            undocumented.
	 */
	static double dist
	( const MinutiaeView & v , const MinutiaeView & w ,
	  int m , double angleWeight , double *distances ) {

		// Initialize the array that keeps track if the 'm'
		// most similar minutiae correspondences betwee 'v'
		// and 'w'
		for ( int i = 0 ; i < m ; i++ ) {
			distances[i] = DBL_MAX;
		}

		// Iterate over all minutiae in 'v'
		for ( int i = 0 ; i < v.getMinutiaeCount() ; i++ ) {

			// Find the most similar minutia from 'w' to the 'i'th of 'v'
			// as well as its distance
			double minDist = DBL_MAX;
			for ( int j = 0 ; j < w.getMinutiaeCount() ; j++ ) {
				double d = FingerTools::dist
						(v.getMinutia(i),w.getMinutia(j),angleWeight);
				if ( d < minDist ) {
					minDist = d;
				}
			}

			// Insert the dissimilarity in 'distances' if it is among the
			// 'm' minimal values
			for ( int k = 0 ; k < m ; k++ ) {
				if ( minDist < distances[k] ) {
					for ( int l = m-1 ; l > k ; l-- ) {
						distances[l] = distances[l-1];
					}
					distances[k] = minDist;
					break;
				}
			}
		}

		// Sum up the values in 'distances' ...
		double d = 0.0;
		for ( int i = 0 ; i < m ; i++ ) {
			if ( distances[i] < DBL_MAX ) {
				d += distances[i];
			}
		}

		// ... which is the result.
		return d;
	}


	/**
	 * @brief
	 *            Computes a similarity measure between two minutiae
	 *            templates.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	double FingerTools::dist
	( const MinutiaeView & v , const MinutiaeView & w ,
	  int n , double angleWeight ) {

		// Number of input minutiae
		int k , l;
		k = v.getMinutiaeCount();
		l = w.getMinutiaeCount();

		// Are 'n' summands possible?
		int m = std::min(k,l);
		if ( n > 0 ) {
			m = std::min(m,n);
		}

		// Special case
		if ( m == 0 ) {
			return 0.0;
		}

		// Allocate memory
		double *distances = (double*)malloc(m*sizeof(double));
		if ( distances == NULL ) {
			cerr << "FingerTools:dist(const MinutiaeView&,"
				 << "const MinutiaeView&,int,double): Out of memory."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// Compute the sum
		double d = thimble::dist(v,w,m,angleWeight,distances);

		// Free memory
		free(distances);

		// Return the result
		return d;
	}

	/**
	 * @brief
	 *            Determines a spatial movement as an affine transform
	 *            that aligns the minutiae in <code>w</code> to the
	 *            minutiae in <code>v</code>.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	AffineTransform FingerTools::align
	( const MinutiaeView & v , const MinutiaeView & w ,
	  int n , double angleWeight ) {

		// Number of input minutiae
		int k , l;
		k = v.getMinutiaeCount();
		l = w.getMinutiaeCount();

		// Are 'n' summands possible?
		int m = std::min(k,l);
		if ( n > 0 ) {
			m = std::min(m,n);
		}

		// Keeps track of the movement that minimizes the dissimilarity
		// between 'v' and the moved 'w'.
		AffineTransform f;

		// Check for special case.
		if ( m == 0 ) {
			return f;
		}

		// Allocate memory
		double *distances = (double*)malloc(m*sizeof(double));
		if ( distances == NULL ) {
			cerr << "FingerTools::align(const MinutiaeView&,"
				 << "const MinutiaeView&,int,double): Out of memory."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// Keeps track of the minimal dissimilarity found.
		double minDist = DBL_MAX;

		// We will use 'u' to save the movement of 'w' under a
		// candidate alignment; we will not use the function
		// 'FingerTools::eval(f,w)' for performance reasons
		// (using 'FingerTools::eval(f,w)' would cause frequent
		//  allocation and frees which is imperformant)
		MinutiaeView u = w;

		// Iterate over all minutiae correspondence
		for ( int i = 0 ; i < v.getMinutiaeCount() ; i++ ) {
			for ( int j = 0 ; j < w.getMinutiaeCount() ; j++ ) {

				// Determine the candidate alignment that moves the 'j'th
				// minutia of 'w' to the 'i'th minutia of 'v'.
				AffineTransform g =
						FingerTools::align(v.getMinutia(i),w.getMinutia(j));

				// Move 'w' under the candidate alignment 'g' and collect
				// the moved minutiae in 'u'
				u.removeAllMinutiae();
				for ( int k = 0 ; k < w.getMinutiaeCount() ; k++ ) {
					Minutia minutia = FingerTools::eval(g,w.getMinutia(k));
					u.addMinutia(minutia);
				}

				// Determine the dissimilarity between 'v' and 'u'
				double d = thimble::dist(v,u,m,angleWeight,distances);

				// If the dissimilarity is smaller, update
				if ( d < minDist ) {
					minDist = d;
					f.assign(g);
				}
			}
		}

		// Free memory, prior returning.
		free(distances);

		// 'f' is the movement that minimizes the dissimilarity between
		// 'v' and 'f(w)' among the iterated candidate alignments.
		return f;
	}

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
	MinutiaeView FingerTools::prealign
	( const MinutiaeView & v , double x , double y , double direction ) {

		AffineTransform f;

		f.v = -x;
		f.w = -y;
		f.irotate(direction);

		return eval(f,v);
	}

	/**
	 * @brief
	 *            Computes the representation of a minutiae view w.r.t.
	 *            a directed point defining the origin of an intrinsic
	 *            coordinate system.
	 *
	 * @details
	 *            see 'FingerTools.h'
	 */
	MinutiaeView FingerTools::prealign
	( const MinutiaeView & v , const DirectedPoint & p ) {

		return prealign(v,p.x,p.y,p.direction.getDirectionAngle());
	}

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
	Minutia FingerTools::prealign( const Minutia & m , double x , double y , double direction ) {

		AffineTransform f;

		f.v = -x;
		f.w = -y;
		f.irotate(direction);

		return eval(f,m);
	}

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
	Minutia FingerTools::prealign( const Minutia & m , const DirectedPoint & p ) {

		return prealign(m,p.x,p.y,p.direction.getDirectionAngle());
	}
}




