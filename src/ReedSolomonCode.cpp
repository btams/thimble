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
 * @file ReedSolomonCode.cpp
 *
 * @brief
 *            Implements functionalities provided by 'ReedSolomonCode.h' which
 *            are related with Reed-Solomon error-correcting codes.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/ecc/ReedSolomonCode.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Attempts to decode a (possibly perturbed
	 *            Reed-Solomon code <i>in original view</i> using an
	 *            algorithm by <b>Gao (2002)</b>.
	 *
	 * @details
	 *            see 'ReedSolomonCode.h'
	 */
	bool ReedSolomonCode::gaodecode
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y , int n , int k ) {

		if ( n <= 0 || k <= 0 || k > n ) {
			cerr << "ReedSolomonCode::gaodecode: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		 SmallBinaryFieldPolynomial g0(f.getField());
		 g0.buildFromRoots(x,n);

		 SmallBinaryFieldPolynomial g1(f.getField());
		 g1.interpolate(x,y,n);

		 SmallBinaryFieldPolynomial
		 	 g(f.getField()) ,
		 	 u(f.getField()) ,
		 	 v(f.getField());

		 pgcd(g,u,v,g0,g1, (n+k)&0x1 ? ((n+k-1)>>1)+1 : ( (n+k)>>1) );

		 SmallBinaryFieldPolynomial f1(f.getField()) , r(f.getField());
		 divRem(f1,r,g,v);

		 if ( r.isZero() && f1.deg() < k ) {
			 f = f1;
			 return true;
		 }

		 return false;
	}
}


