/*
 *  THIMBLE --- Research Library for Development and Analysis of
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
 * @file MathTools.h
 *
 * @brief
 *            Provides functions that may be commonly used or a variety of
 *            purposes, especially mathematical ones.
 *
 * @details
 *            The functions are provided as static class functions of the
 *            <code>\link thimble::MathTools\endlink</code> class.
 *
 * @author Benjamin Tams
 *
 * @see thimble::MathTools
 */

#ifndef THIMBLE_MATHTOOLS_H_
#define THIMBLE_MATHTOOLS_H_

#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Provides mathematical functions that may be commonly used
	 *            for a variety of purposes.
	 */
	class THIMBLE_DLL MathTools {

	public:

		/**
		 * @brief
		 *            Rounds a long double value to its nearest integer.
		 *
		 * @param v
		 *            The long double value being rounded.
		 *
		 * @return
		 *            The integer (encoded as a long double) that is closest
		 *            to <code>v</code>.
		 */
		static long double round( long double v );

		/**
		 * @brief
		 *            Rounds a double value to its nearest integer.
		 *
		 * @param v
		 *            The double value being rounded.
		 *
		 * @return
		 *            The integer (encoded as a double) that is closest
		 *            to <code>v</code>.
		 */
		static double round( double v );

		/**
		 * @brief
		 *            Rounds a float to its nearest integer.
		 *
		 * @param v
		 *            The float value being rounded.
		 *
		 * @return
		 *            The integer (encoded as a float) that is closest
		 *            to <code>v</code>.
		 */
		static float round( float v );


		/**
		 * @brief
		 *            Computes the number of bits required to represent the
		 *            given integer.
		 *
		 * @details
		 *            More precisely, the function computes the least positive
		 *            integer <i>m</i> such that \f$ 2^m>n\f$.
		 *
		 * @param n
		 *            integer of which the bit-length is sought
		 *
		 * @return
		 *            the bit-length of the integer <code>n</code>
		 */
		static int numBits( register uint64_t n );

        /**
         * @brief
         *            Computes the Hamming weight of an integer.
         *
         * @details
         *            Write
         *            \f[
         *             n = \sum_{j=0}^{b-1}n_j\cdot 2^j
         *            \f]
         *            where \f$b_j\in\{0,1\}\f$. The function determines
         *            the number of <i>j</i> where \f$n_j\neq 0\f$; more
         *            specifically it returns
         *            \f[
         *             \#\{j=0,...,b-1~|~n_j\neq 0\}.
         *            \f]
         *
         * @param n
         *            The integer of which the Hamming weight is computed.
         *
         * @return
         *            The Hamming weight of <i>n</i>.
         */
        static int hammingWeight( register uint64_t n );


        /**
         * @brief
         *            Determines the scalar product of two vectors of length
         *            being a multiple of 64 having only 0 or 1 as entries.
         *
         * @details
         *            Two vectors of length <code>64*n</code> are passed to
         *            the function and the integral scalar product is
         *            returned. The vector that is represented by the array
         *            <code>v</code> (say; the same holds for the second
         *            vector represented by <code>w</code>) equals the
         *            following:
         *
         *            Write for the <i>i</i>th element in the array
         *            <code>v</code>
         *            \f[
         *             \sum_{j=0}^{63}b_{i,j}\cdot 2^j
         *            \f]
         *            where \f$b_{i,j}\in\{0,1\}\f$. Then <code>(v,n)</code>
         *            represent the vector
         *            \f[
         * (b_{0,0},b_{0,1},...,b_{0,63},b_{1,0},b_{1,1},...,b_{1,63},...,...,
         *  b_{n-1,0},b_{n-1,1},...,b_{n-1,63})
         *            \f]
         *
         * @param v
         *            An array containing <code>n</code> valid 64 bit integers.
         *
         * @param w
         *            An array containing <code>n</code> valid 64 bit integers.
         *
         * @param n
         *            The number of 64 bit integers contained in
         *            <code>v</code> and in <code>w</code>
         *
         * @return
         *            The scalar product between the two vectors.
         */
        static int sprod64
        ( const uint64_t *v , const uint64_t *w , int n );

        /**
         * @brief
         *            Determines the scalar product of two vectors of length
         *            being a multiple of 32 having only 0 or 1 as entries.
         *
         * @details
         *            Two vectors of length <code>32*n</code> are passed to
         *            the function and the integral scalar product is
         *            returned. The vector that is represented by the array
         *            <code>v</code> (say; the same holds for the second
         *            vector represented by <code>w</code>) equals the
         *            following:
         *
         *            Write for the <i>i</i>th element in the array
         *            <code>v</code>
         *            \f[
         *             \sum_{j=0}^{31}b_{i,j}\cdot 2^j
         *            \f]
         *            where \f$b_{i,j}\in\{0,1\}\f$. Then <code>(v,n)</code>
         *            represent the vector
         *            \f[
         * (b_{0,0},b_{0,1},...,b_{0,31},b_{1,0},b_{1,1},...,b_{1,31},...,...,
         *  b_{n-1,0},b_{n-1,1},...,b_{n-1,31})
         *            \f]
         *
         * @param v
         *            An array containing <code>n</code> valid 32 bit integers.
         *
         * @param w
         *            An array containing <code>n</code> valid 32 bit integers.
         *
         * @param n
         *            The number of 32 bit integers contained in
         *            <code>v</code> and in <code>w</code>
         *
         * @return
         *            The scalar product between the two vectors.
         */
        static int sprod32
        ( const uint32_t *v , const uint32_t *w , int n );

        /**
         * @brief
         *            Determines the scalar product of two vectors of length
         *            being a multiple of 16 having only 0 or 1 as entries.
         *
         * @details
         *            Two vectors of length <code>16*n</code> are passed to
         *            the function and the integral scalar product is
         *            returned. The vector that is represented by the array
         *            <code>v</code> (say; the same holds for the second
         *            vector represented by <code>w</code>) equals the
         *            following:
         *
         *            Write for the <i>i</i>th element in the array
         *            <code>v</code>
         *            \f[
         *             \sum_{j=0}^{15}b_{i,j}\cdot 2^j
         *            \f]
         *            where \f$b_{i,j}\in\{0,1\}\f$. Then <code>(v,n)</code>
         *            represent the vector
         *            \f[
         * (b_{0,0},b_{0,1},...,b_{0,15},b_{1,0},b_{1,1},...,b_{1,15},...,...,
         *  b_{n-1,0},b_{n-1,1},...,b_{n-1,15})
         *            \f]
         *
         * @param v
         *            An array containing <code>n</code> valid 16 bit integers.
         *
         * @param w
         *            An array containing <code>n</code> valid 16 bit integers.
         *
         * @param n
         *            The number of 16 bit integers contained in
         *            <code>v</code> and in <code>w</code>
         *
         * @return
         *            The scalar product between the two vectors.
         */
        static int sprod16
        ( const uint16_t *v , const uint16_t *w , int n );

        /**
         * @brief
         *            Determines the scalar product of two vectors of length
         *            being a multiple of 8 having only 0 or 1 as entries.
         *
         * @details
         *            Two vectors of length <code>8*n</code> are passed to
         *            the function and the integral scalar product is
         *            returned. The vector that is represented by the array
         *            <code>v</code> (say; the same holds for the second
         *            vector represented by <code>w</code>) equals the
         *            following:
         *
         *            Write for the <i>i</i>th element in the array
         *            <code>v</code>
         *            \f[
         *             \sum_{j=0}^{7}b_{i,j}\cdot 2^j
         *            \f]
         *            where \f$b_{i,j}\in\{0,1\}\f$. Then <code>(v,n)</code>
         *            represent the vector
         *            \f[
         * (b_{0,0},b_{0,1},...,b_{0,7},b_{1,0},b_{1,1},...,b_{1,7},...,...,
         *  b_{n-1,0},b_{n-1,1},...,b_{n-1,7})
         *            \f]
         *
         * @param v
         *            An array containing <code>n</code> valid 8 bit integers.
         *
         * @param w
         *            An array containing <code>n</code> valid 8 bit integers.
         *
         * @param n
         *            The number of 8 bit integers contained in
         *            <code>v</code> and in <code>w</code>
         *
         * @return
         *            The scalar product between the two vectors.
         */
        static int sprod8
        ( const uint8_t *v , const uint8_t *w , int n );

        /**
         * @brief
         *             Determines the Hamming weight of a bit-vector
         *             having a length that is a multiple of 64 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 64 bit integers that
         *             encodes <code>64*n</code> bits.
         *
         * @param n
         *             The number of valid 64 bit integers contained in
         *             <code>v</code>.
         *
         * @return
         *             The number of 1s encoded in the vector
         *             <code>(v,n)</code>.
         */
        static int hw64
        ( register const uint64_t *v , register int n );

        /**
         * @brief
         *             Determines the Hamming weight of a bit-vector
         *             having a length that is a multiple of 32 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 32 bit integers that
         *             encodes <code>32*n</code> bits.
         *
         * @param n
         *             The number of valid 32 bit integers contained in
         *             <code>v</code>.
         *
         * @return
         *             The number of 1s encoded in the vector
         *             <code>(v,n)</code>.
         */
        static int hw32( const uint32_t *v , int n );

        /**
         * @brief
         *             Determines the Hamming weight of a bit-vector
         *             having a length that is a multiple of 16 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 16 bit integers that
         *             encodes <code>16*n</code> bits.
         *
         * @param n
         *             The number of valid 16 bit integers contained in
         *             <code>v</code>.
         *
         * @return
         *             The number of 1s encoded in the vector
         *             <code>(v,n)</code>.
         */
        static int hw16( const uint16_t *v , int n );

        /**
         * @brief
         *             Determines the Hamming weight of a bit-vector
         *             having a length that is a multiple of 8 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 8 bit integers that
         *             encodes <code>8*n</code> bits.
         *
         * @param n
         *             The number of valid 8 bit integers contained in
         *             <code>v</code>.
         *
         * @return
         *             The number of 1s encoded in the vector
         *             <code>(v,n)</code>.
         */
        static int hw8( const uint8_t *v , int n );

        /**
         * @brief
         *             Determines the Hamming distance between two bit-vectors
         *             having a length that is a multiple of 64 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 64 bit integers that
         *             encodes <code>64*n</code> bits.
         *
         * @param w
         *             An array of <code>n</code> valid 64 bit integers that
         *             encodes <code>64*n</code> bits.
         *
         * @param n
         *             A number of valid 64 bit integers contained in
         *             <code>v</code> and <code>w</code>.
         *
         * @return
         *             The number of bit positions in which the vectors
         *             <code>(v,n)</code> and <code>(w,n)</code> differ.
         */
        static int hd64
        ( register const uint64_t *v , register const uint64_t *w ,
          register int n );

        /**
         * @brief
         *             Determines the Hamming distance between two bit-vectors
         *             having a length that is a multiple of 32 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 32 bit integers that
         *             encodes <code>32*n</code> bits.
         *
         * @param w
         *             An array of <code>n</code> valid 32 bit integers that
         *             encodes <code>32*n</code> bits.
         *
         * @param n
         *             A number of valid 32 bit integers contained in
         *             <code>v</code> and <code>w</code>.
         *
         * @return
         *             The number of bit positions in which the vectors
         *             <code>(v,n)</code> and <code>(w,n)</code> differ.
         */
        static int hd32
        ( const uint32_t *v , const uint32_t *w , int n );

        /**
         * @brief
         *             Determines the Hamming distance between two bit-vectors
         *             having a length that is a multiple of 16 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 16 bit integers that
         *             encodes <code>32*n</code> bits.
         *
         * @param w
         *             An array of <code>n</code> valid 16 bit integers that
         *             encodes <code>16*n</code> bits.
         *
         * @param n
         *             A number of valid 16 bit integers contained in
         *             <code>v</code> and <code>w</code>.
         *
         * @return
         *             The number of bit positions in which the vectors
         *             <code>(v,n)</code> and <code>(w,n)</code> differ.
         */
        static int hd16
        ( const uint16_t *v , const uint16_t *w , int n );

        /**
         * @brief
         *             Determines the Hamming distance between two bit-vectors
         *             having a length that is a multiple of 8 bits.
         *
         * @param v
         *             An array of <code>n</code> valid 8 bit integers that
         *             encodes <code>8*n</code> bits.
         *
         * @param w
         *             An array of <code>n</code> valid 8 bit integers that
         *             encodes <code>8*n</code> bits.
         *
         * @param n
         *             A number of valid 8 bit integers contained in
         *             <code>v</code> and <code>w</code>.
         *
         * @return
         *             The number of bit positions in which the vectors
         *             <code>(v,n)</code> and <code>(w,n)</code> differ.
         */
        static int hd8
        ( const uint8_t *v , const uint8_t *w , int n );

		/**
		 * @brief     Computes the carry-less product of two 32-bit integers
		 *
		 * @details   Given two integers in binary representation their ordinary
		 *            product is computed with a carry which is dismissed for
		 *            the carry-less product.
		 *
		 * @param a
		 *            first factor
		 * @param b
		 *            second factor
		 *
		 * @return    the carry-less product of <code>a</code> and
		 *            <code>b</code>
		 */
		static uint64_t clmul( uint32_t a , uint32_t b );

        /**
         * @brief
         *            Performs exclusive or operations on arrays of 64-bit
         *            integers.
         *
         * @param out
         *            On output the <i>i</i>th 64-bit integer equals the
         *            xor of the input <code>in[i]</code> and
         *            <code>in2[i]</code>.
         *
         * @param in1
         *            First input array containing <code>n</code> valid
         *            integers of 64-bit width.
         *
         * @param in2
         *            Second input array containing <code>n</code> valid
         *            integers of 64-bit width.
         *
         * @param n
         *            Positive integer specifying the number of 64-bit
         *            integers in <code>in1</code> and <code>in2</code>
         *            as well as a number of 64-bit integers that the array
         *            <code>out</code> is capable to hold.
         *
         * @details
         *            The arrays <code>out</code>, <code>in1</code> and
         *            <code>in2</code> are allowed to be equals but not
         *            to overlap at a position different from 0.
         *
         * @warning
         *            If the requirements above are not fulfilled, the method
         *            runs into undocumented behaviour.
         */
        static void mxor64
        ( register uint64_t *out ,
          register const uint64_t *in1 , register const uint64_t *in2 ,
          register int n );

        /**
         * @brief
         *            Performs exclusive or operations on arrays of 32-bit
         *            integers.
         *
         * @param out
         *            On output the <i>i</i>th 32-bit integer equals the
         *            xor of the input <code>in[i]</code> and
         *            <code>in2[i]</code>.
         *
         * @param in1
         *            First input array containing <code>n</code> valid
         *            integers of 32-bit width.
         *
         * @param in2
         *            Second input array containing <code>n</code> valid
         *            integers of 32-bit width.
         *
         * @param n
         *            Positive integer specifying the number of 32-bit
         *            integers in <code>in1</code> and <code>in2</code>
         *            as well as a number of 32-bit integers that the array
         *            <code>out</code> is capable to hold.
         *
         * @details
         *            The arrays <code>out</code>, <code>in1</code> and
         *            <code>in2</code> are allowed to be equals but not
         *            to overlap at a position different from 0.
         *
         * @warning
         *            If the requirements above are not fulfilled, the method
         *            runs into undocumented behaviour.
         */
        static void mxor32
        ( uint32_t *out , const uint32_t *in1 , const uint32_t *in2 , int n );

        /**
         * @brief
         *            Performs exclusive or operations on arrays of 16-bit
         *            integers.
         *
         * @param out
         *            On output the <i>i</i>th 16-bit integer equals the
         *            xor of the input <code>in[i]</code> and
         *            <code>in2[i]</code>.
         *
         * @param in1
         *            First input array containing <code>n</code> valid
         *            integers of 16-bit width.
         *
         * @param in2
         *            Second input array containing <code>n</code> valid
         *            integers of 16-bit width.
         *
         * @param n
         *            Positive integer specifying the number of 16-bit
         *            integers in <code>in1</code> and <code>in2</code>
         *            as well as a number of 16-bit integers that the array
         *            <code>out</code> is capable to hold.
         *
         * @details
         *            The arrays <code>out</code>, <code>in1</code> and
         *            <code>in2</code> are allowed to be equals but not
         *            to overlap at a position different from 0.
         *
         * @warning
         *            If the requirements above are not fulfilled, the method
         *            runs into undocumented behaviour.
         */
        static void mxor16
        ( uint16_t *out , const uint16_t *in1 , const uint16_t *in2 , int n );

        /**
         * @brief
         *            Performs exclusive or operations on arrays of 8-bit
         *            integers.
         *
         * @param out
         *            On output the <i>i</i>th 8-bit integer equals the
         *            xor of the input <code>in[i]</code> and
         *            <code>in2[i]</code>.
         *
         * @param in1
         *            First input array containing <code>n</code> valid
         *            integers of 8-bit width.
         *
         * @param in2
         *            Second input array containing <code>n</code> valid
         *            integers of 8-bit width.
         *
         * @param n
         *            Positive integer specifying the number of 8-bit
         *            integers in <code>in1</code> and <code>in2</code>
         *            as well as a number of 8-bit integers that the array
         *            <code>out</code> is capable to hold.
         *
         * @details
         *            The arrays <code>out</code>, <code>in1</code> and
         *            <code>in2</code> are allowed to be equals but not
         *            to overlap at a position different from 0.
         *
         * @warning
         *            If the requirements above are not fulfilled, the method
         *            runs into undocumented behaviour.
         */
        static void mxor8
        ( uint8_t *out , const uint8_t *in1 , const uint8_t *in2 , int n );


        /**
         * @brief
         *            Tests whether each integer in an array is zero.
         *
         * @param array
         *            Contains <code>n</code> valid 64 bit integers.
         *
         * @param n
         *            The number of 64 bit integers in <code>array</code>.
         *
         * @return
         *            <code>true</code> if all integers in <code>array</code>
         *            are zero; otherwise, if there exists a non-zero integer,
         *            the function returns <code>false</code>.
         */
        static bool zeroTest64
        ( register const uint64_t *array , register int n );

        /**
         * @brief
         *            Tests whether each integer in an array is zero.
         *
         * @param array
         *            Contains <code>n</code> valid 32 bit integers.
         *
         * @param n
         *            The number of 32 bit integers in <code>array</code>.
         *
         * @return
         *            <code>true</code> if all integers in <code>array</code>
         *            are zero; otherwise, if there exists a non-zero integer,
         *            the function returns <code>false</code>.
         */
        static bool zeroTest32
        ( const uint32_t *array , int n );

        /**
         * @brief
         *            Tests whether each integer in an array is zero.
         *
         * @param array
         *            Contains <code>n</code> valid 16 bit integers.
         *
         * @param n
         *            The number of 16 bit integers in <code>array</code>.
         *
         * @return
         *            <code>true</code> if all integers in <code>array</code>
         *            are zero; otherwise, if there exists a non-zero integer,
         *            the function returns <code>false</code>.
         */
        static bool zeroTest16
        ( const uint16_t *array , int n );

        /**
         * @brief
         *            Tests whether each integer in an array is zero.
         *
         * @param array
         *            Contains <code>n</code> valid 8 bit integers.
         *
         * @param n
         *            The number of 8 bit integers in <code>array</code>.
         *
         * @return
         *            <code>true</code> if all integers in <code>array</code>
         *            are zero; otherwise, if there exists a non-zero integer,
         *            the function returns <code>false</code>.
         */
        static bool zeroTest8
        ( const uint8_t *array , int n );


		/**
		 * @brief
		 *            Generates a random 8-bit word.
		 *
		 * @details
		 *            If <code>tryRandom</code> is <code>true</code>
		 *            then the function attempts to use a number generator
		 *            that is cryptographically secure such as
		 *            <code>/dev/urandom</code> on UNIX or <code>rand_s</code>
		 *            on Windows. Otherwise, the functions behaves like a
		 *            wrapper around the standard <code>rand()</code> function.
		 *
		 * @param tryRandom
		 *            specifies whether the function is advised to use
		 *            a cryptographic random generator.
		 *
		 * @return
		 *            A random 8-bit word.
		 *
		 * @warning
		 *            Note, that if we are on neither a UNIX or Windows
		 *            platform but use a generic compiler, the function may
		 *            not be suitable for cryptographic purposes.
		 */
		static uint8_t rand8( bool tryRandom = false );

		/**
		 * @brief
		 *            Generates a random 16-bit word.
		 *
		 * @details
		 *            If <code>tryRandom</code> is <code>true</code>
		 *            then the function attempts to use a number generator
		 *            that is cryptographically secure such as
		 *            <code>/dev/urandom</code> on UNIX or <code>rand_s</code>
		 *            on Windows. Otherwise, the functions behaves like a
		 *            wrapper around the standard <code>rand()</code> function.
		 *
		 * @param tryRandom
		 *            specifies whether the function is advised to use
		 *            a cryptographic random generator.
		 *
		 * @return
		 *            A random 16-bit word.
		 *
		 * @warning
		 *            Note, that if we are on neither a UNIX or Windows
		 *            platform but use a generic compiler, the function may
		 *            not be suitable for cryptographic purposes.
		 */
		static inline uint16_t rand16( bool tryRandom = false ) {
			uint16_t r , l;
			r = rand8(tryRandom);
			r<<=8;
			l = rand8(tryRandom);
			return r+l;
		}

		/**
		 * @brief
		 *            Generates a random 32-bit word.
		 *
		 * @details
		 *            If <code>tryRandom</code> is <code>true</code>
		 *            then the function attempts to use a number generator
		 *            that is cryptographically secure such as
		 *            <code>/dev/urandom</code> on UNIX or <code>rand_s</code>
		 *            on Windows. Otherwise, the functions behaves like a
		 *            wrapper around the standard <code>rand()</code> function.
		 *
		 * @param tryRandom
		 *            specifies whether the function is advised to use
		 *            a cryptographic random generator.
		 *
		 * @return
		 *            A random 32-bit word.
		 *
		 * @warning
		 *            Note, that if we are on neither a UNIX or Windows
		 *            platform but use a generic compiler, the function may
		 *            not be suitable for cryptographic purposes.
		 */
		static inline uint32_t rand32( bool tryRandom = false ) {
			uint32_t r , l;
			r = rand16(tryRandom);
			r<<=16;
			l = rand16(tryRandom);
			return r+l;
		}

		/**
		 * @brief
		 *            Generates a random 64-bit word.
		 *
		 * @details
		 *            If <code>tryRandom</code> is <code>true</code>
		 *            then the function attempts to use a number generator
		 *            that is cryptographically secure such as
		 *            <code>/dev/urandom</code> on UNIX or <code>rand_s</code>
		 *            on Windows. Otherwise, the functions behaves like a
		 *            wrapper around the standard <code>rand()</code> function.
		 *
		 * @param tryRandom
		 *            specifies whether the function is advised to use
		 *            a cryptographic random generator.
		 *
		 * @return
		 *            A random 64-bit word.
		 *
		 * @warning
		 *            Note, that if we are on neither a UNIX or Windows
		 *            platform but use a generic compiler, the function may
		 *            not be suitable for cryptographic purposes.
		 */
		static inline uint64_t rand64( bool tryRandom = false ) {
			uint64_t r , l;
			r = rand32(tryRandom);
			r<<=32;
			l = rand32(tryRandom);
			return r+l;
		}
	};

}


#endif /* THIMBLE_MATHTOOLS_H_ */
