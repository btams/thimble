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
 * @file NTTools.h
 *
 * @brief
 *            Provides a class that provides a variety of static functions
 *            and methods related to number theory that may be useful for a
 *            variety of purposes.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_NTTOOLS_H
#define THIMBLE_NTTOOLS_H

#include <stdint.h>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/BinaryPolynomial.h>
#include <thimble/math/linalg/BinaryVector.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Provides static functions that may be useful for different
     *            purposes with relation to number theory.
     */
    class THIMBLE_DLL NTTools {

    private:

        /**
         * @brief
         *            Private standard constructor.
         *
         * @details
         *            The standard constructor is private to avoid users
         *            to create instance of this tool class that provides only
         *            static functions.
         */
        inline NTTools() {}

    public:

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
        static std::vector< std::pair<uint64_t,int> > trialDivision
            ( uint64_t n );

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
        static void conv( BinaryPolynomial & f , const BinaryVector & v );

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
        static void conv( BinaryVector & v , const BinaryPolynomial & f );
    };

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
    inline void conv( BinaryPolynomial & f , const BinaryVector & v ) {

        NTTools::conv(f,v);
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
    inline void conv( BinaryVector & v , const BinaryPolynomial & f ) {

        NTTools::conv(v,f);
    }
}

#endif /* THIMBLE_NTTOOLS_H */
