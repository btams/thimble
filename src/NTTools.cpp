/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file NTTools.cpp
 *
 * @brief
 *            Implements a variety of functions and methods related to
 *            number theory, as provided by the 'NTTools.h' header, that may
 *            be useful for a variety of purposes.
 *
 * @author Benjamin Tams
 */
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <iostream>

#include <thimble/math/numbertheory/NTTools.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Returns the prime decomposition of a natural number
	 *            having at most 64 bits.
	 *
	 * @details
	 *            The function factors the input number using trial
	 *            division and, therefore, may not be the best choice
	 *            for factoring number having large prime divisors.
	 *
	 *            The function returns a vector of pairs <i>(p,e)</i> such
	 *            that the <i>p</i> are all distinct prime divisors of
	 *            <i>n</i> (sorted increasingly) and the <i>e</i> are their
	 *            respective multiplicities. More specifically, the pairs
	 *            fulfill
	 *            \f[
	 *             n = \prod_{(p,e)}p^e.
	 *            \f]
	 *
	 * @param n
	 *            The number to be factored.
	 *
	 * @return
	 *            A vector of pairs whose first component corresponds
	 *            to a prime <i>p</i> dividing <i>n</i> and the second
	 *            <i>e</i> to its multiplicity.
	 *
	 * @warning
	 *            If <i>n</i> is zero, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            'EXIT_FAILURE'.
	 */
    vector< pair<uint64_t,int> > NTTools::trialDivision( uint64_t n ) {

        if ( n == 0 ) {
            cerr << "NTTools::factor: cannot factor 0." << endl;
            exit(EXIT_FAILURE);
        }

        uint64_t m = (uint64_t)sqrt((double)n);
        vector< pair<uint64_t,int> > list;

        // Iterate through all possible divisors ignoring
        // non-primes
        for ( uint64_t p = 2 ; p <= m ; p++ ) {

            int e = 0;

            // Determine the valuation of 'p' at 'n'
            while ( n % p == 0 ) {
                n /= p;
                ++e;
            }

            // If the multiplicity is greater than 0, the integer
            // 'p' is a prime dividing 'n'.
            if ( e > 0 ) {
                m = (uint64_t)sqrt((double)n);
                list.push_back( pair<uint64_t,int>(p,e) );
            }
        }

        // 'n' must be a prime divisor of single multiplicity if not equals 1.
        if ( n != 1 ) {
            list.push_back( pair<uint64_t,int>(n,1) );
        }

        return list;
    }

    /**
     * @brief
     *            Converts a binary vector into a binary polynomial.
     *
     * @details
     *            Write
     *            \f[
     *             v=(v_0,....,v_{n-1}).
     *            \f]
     *            Then the method sets the polynomial <i>f</i> as follows
     *            \f[
     *             f(X)=v_0+v_1\cdot X+...+v_{n-1}\cdot X^{n-1}
     *            \f]
     *
     * @param v
     *            The vector that is converted to a polynomial.
     *
     * @param f
     *            The polynomial into which the vector <i>v</i> is
     *            converted.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void NTTools::conv( BinaryPolynomial & f , const BinaryVector & v ) {

        // *******************************************************************
        // * WARNING: THIS METHOD USES A 'getData_nonconst()' METHOD WHICH ***
        // * YOU SHOULD ONLY USE ON YOURSELF IF YOU REALLY KNOW WHAT YOU ARE *
        // * DOING ***********************************************************
        // *******************************************************************

        int n = v.getLength();

        // Makes the polynomial the monomial 'X^(n-1)'
        f.setZero();
        f.setCoeff(n-1);

        // Binary copy
        memcpy
        (f.getData_nonconst(),v.getData(),
         (n/32+(n%32?1:0))*sizeof(uint32_t));

        // Since upper coefficients may be zero, the degree could
        // be scmaller
        f.degreeCouldBeSmaller();
    }

    /**
     * @brief
     *            Converts a binary polynomial into a binary vector.
     *
     * @details
     *            Write for the coefficients of the polynomial
     *            \f[
     *             f(X)=f_0+f_1\cdot X+...+f_{n-1}\cdot X^{n-1}
     *            \f]
     *            and for the vector
     *            \f[
     *             v=(v_0,...,v_{n-1}).
     *            \f]
     *            The method performs the assignment
     *            \f[
     *             v_j\leftarrow f_j~~where~~j=0,...,n-1.
     *            \f]
     *            Note that the vector <i>v</i> must have a length greater
     *            than the degree of <i>f</i>.
     *
     * @param v
     *            Vector of length greater than the degree of <i>f</i>
     *            into which the polynomial <i>f</i> is converted.
     *
     * @param f
     *            The polynomial converted to a binary vector
     *            stored in <i>v</i>.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void NTTools::conv( BinaryVector & v , const BinaryPolynomial & f ) {

        // *******************************************************************
        // * WARNING: THIS METHOD USES A 'getData_nonconst()' METHOD WHICH ***
        // * YOU SHOULD ONLY USE ON YOURSELF IF YOU REALLY KNOW WHAT YOU ARE *
        // * DOING ***********************************************************
        // *******************************************************************

        int n = f.deg()+1;

        if ( n > v.getLength() ) {
            cerr << "NTTools::conv: assigned vector has too small length."
                 << endl;
            exit(EXIT_FAILURE);
        }

        v.setZero();

        memcpy
        (v.getData_nonconst(),f.getData(),
         (n/32+(n%32?1:0))*sizeof(uint32_t));
    }
}
