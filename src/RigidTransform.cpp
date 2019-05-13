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
 * @file RigidTransform.cpp
 *
 * @brief
 *            Implements the functionality provided by 'RigidTransform.h'
 *            which is to compute the spatial movement that minimizes the
 *            sum of least square distances between two point correspondences.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>

#include <thimble/math/AffineTransform.h>
#include <thimble/math/RigidTransform.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Determines a spatial transform that minimizes
	 *           the sum of squared distances between two given
	 *           two-dimensional points correspondences.
	 *
	 * @details
	 *           see 'RigidTransform.h'
	 */
	AffineTransform RigidTransform::align
		( const double *dx , const double *dy ,
		  const double *mx , const double *my , int n ) {

		AffineTransform f;

		if ( n <= 1 ) {
			return f;
		}

		// Allocate memory for holding the input points shifted w.r.t. their
		// respective origin.
		double *mcx , *mcy , *dcx , *dcy;
		mcx = (double*)malloc( n * sizeof(double) );
		mcy = (double*)malloc( n * sizeof(double) );
		dcx = (double*)malloc( n * sizeof(double) );
		dcy = (double*)malloc( n * sizeof(double) );
		if ( mcx == NULL || mcy == NULL || dcx == NULL || dcy == NULL ) {
			cerr << "RigidTransform::align: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Compute '(mx0,my0)' and '(dx0,dy0)' the barycenters of
		// '{(mx[i],my[i])}' and '{(dx[i],dy[i])}', respectively.
		double mx0 , my0 , dx0 , dy0;
		mx0 = 0.0;
		my0 = 0.0;
		dx0 = 0.0;
		dy0 = 0.0;
		for ( int i = 0 ; i < n ; i++ ) {
			mx0 += mx[i];
			my0 += my[i];
			dx0 += dx[i];
			dy0 += dy[i];
		}
		mx0 /= (double)n;
		my0 /= (double)n;
		dx0 /= (double)n;
		dy0 /= (double)n;

		// Shift the input points w.r.t. their respective origins.
		for ( int i = 0 ; i < n ; i++ ) {
			mcx[i] = mx[i] - mx0;
			mcy[i] = my[i] - my0;
			dcx[i] = dx[i] - dx0;
			dcy[i] = dy[i] - dy0;
		}


		// A = ( 'a' 'b')
		//     ( 'c' 'd')
		// is the correlation matrix of '(mcx[i],mcy[i])_i' and
		// '(dcx[i],dcy[i])_i'.
		double a , b , c , d;
		a = 0.0;
		b = 0.0;
		c = 0.0;
		d = 0.0;
		for ( int i = 0 ; i < n ; i++ ) {
			a += mcx[i] * dcx[i];
			b += mcx[i] * dcy[i];
			c += mcy[i] * dcx[i];
			d += mcy[i] * dcy[i];
		}

		double trmax = -DBL_MAX;

		double agamma , asigma;
		{
			double ad2 , cb2;
			ad2 = a+d;
			ad2 *= ad2;
			cb2 = c-b;
			cb2 *= cb2;

			agamma = ad2 / ( cb2+ad2 );
			asigma = 1.0 - agamma;
			agamma = sqrt(agamma);
			asigma = sqrt(asigma);
		}

		double tr;
		if ( (tr=(a+d)*agamma+(b-c)*asigma) > trmax ) {
			f.a =  agamma;
			f.b = -asigma;
			f.c =  asigma;
			f.d =  agamma;
			trmax = tr;
		}

		if ( (tr=(a+d)*(-agamma)+(b-c)*asigma) > trmax ) {
			f.a = -agamma;
			f.b = -asigma;
			f.c =  asigma;
			f.d = -agamma;
			trmax = tr;
		}

		if ( (tr=(a+d)*agamma+(b-c)*(-asigma)) > trmax ) {
			f.a =  agamma;
			f.b =  asigma;
			f.c = -asigma;
			f.d =  agamma;
			trmax = tr;
		}

		if ( (tr=(a+d)*(-agamma)+(b-c)*(-asigma)) > trmax ) {
			f.a = -agamma;
			f.b =  asigma;
			f.c = -asigma;
			f.d = -agamma;
			trmax = tr;
		}

		// Now, 'f' already has the optimal rotation part.
		// We still need to adopt the optimal translation by mapping
		// the barycenters under rotation.
		f.v = dx0  - f.a * mx0 - f.b * my0;
		f.w = dy0  - f.c * mx0 - f.d * my0;

		// Free allocated memory.
		free(mcx);
		free(mcy);
		free(dcx);
		free(dcy);

		return f;
	}
}



