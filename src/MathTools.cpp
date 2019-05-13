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
 * @file MathTools.cpp
 *
 * @brief
 *            Implements the functions provided by 'MathTools.h' which may
 *            be commonly used or a variety of purposes, especially
 *            mathematical ones.
 *
 * @details
 *            see 'MathTools.h'
 *
 * @author Benjamin Tams
 */

#include "config.h"
#ifdef THIMBLE_USE_RAND_S
#define _CRT_RAND_S
#endif
#include <stdint.h>
#include <cmath>
#include <iostream>

#include <thimble/math/MathTools.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

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
	long double MathTools::round( long double v ) {

		return THIMBLE_ROUND(v);
	}

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
	double MathTools::round( double v ) {

		return THIMBLE_ROUND(v);
	}

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
	float MathTools::round( float v ) {

		return THIMBLE_ROUND(v);
	}

	/**
	 * @brief
	 *            Computes the number of bits required to represent the
	 *            given integer.
	 *
	 * @details
	 *            see 'MathTools.h'
	 */
	int MathTools::numBits( register uint64_t n ) {

		register int b = 0;

		while ( n != 0 ) {
			n >>= 1;
			++b;
		}

		return b;
	}

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
    int MathTools::hammingWeight( register uint64_t n ) {

        // Brian Kernighan algorithm

        register int w = 0;

        while ( n ) {

            n &= n-1;
            ++w;
        }

        return w;
    }

	/**
	 * @brief     Computes the carry-less product of two 32-bit integers
	 *
	 * @details   see 'MathTools.h'
	 */
	uint64_t MathTools::clmul( uint32_t a , uint32_t b ) {

#ifdef THIMBLE_GCC_X86_PCLMULQDQ
		uint64_t A , B;

		A = a;
		B = b;

		asm volatile("movups %0,%%xmm1;" : : "m"(A) );
		asm volatile("movups %0,%%xmm2;" : : "m"(B) );
		asm volatile("pclmulqdq $0x00,%xmm1,%xmm2;" );
		asm volatile("movups %%xmm2,%0;" : "=m"(A) );

		return A;
#else
		uint64_t A , B , C;

		if ( a < b ) {
			A = a;
			B = b;
		} else {
			A = b;
			B = a;
		}

		C = 0;

		while ( A ) {
			if ( A & 0x1 )
				C ^= B;
			B <<= 1;
			A >>= 1;
		}

		return C;
#endif
	}


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
    int MathTools::sprod64
    ( register const uint64_t *v , register const uint64_t *w , register int n ) {

        register int hw = 0;

        register int j;

        for ( j = 0 ; j < n ; j++ , v++ , w++ ) {
            hw += hammingWeight((*v)&(*w));
        }

        return hw;
    }

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
    int MathTools::sprod32
    ( const uint32_t *v , const uint32_t *w , int n ) {

        int hw = sprod64((uint64_t*)v,(uint64_t*)w,n/2);

        if ( (n & 0x1) ) {
            hw += hammingWeight(v[n-1]&w[n-1]);
        }

        return hw;
    }

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
    int MathTools::sprod16
    ( const uint16_t *v , const uint16_t *w , int n ) {

        int hw = sprod32((uint32_t*)v,(uint32_t*)w,n/2);

        if ( (n & 0x1) ) {
            hw += hammingWeight(v[n-1]&w[n-1]);
        }

        return hw;
    }

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
    int MathTools::sprod8
    ( const uint8_t *v , const uint8_t *w , int n ) {

        int hw = sprod16((uint16_t*)v,(uint16_t*)w,n/2);

        if ( (n & 0x1) ) {
            hw += hammingWeight(v[n-1]&w[n-1]);
        }

        return hw;
    }

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
    int MathTools::hw64
    ( register const uint64_t *v , register int n ) {

        register int hw;
        register int i;

        hw = 0;

        for ( i = 0 ; i < n ; i++ , v++ ) {
            hw += hammingWeight(v[i]);
        }

        return hw;
    }

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
    int MathTools::hw32( const uint32_t *v , int n ) {

        int hw = hw64((uint64_t*)v,n/2);

        if ( n & 0x1 ) {
            hw += hammingWeight(v[n-1]);
        }

        return hw;
    }

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
    int MathTools::hw16( const uint16_t *v , int n ) {

        int hw = hw32((uint32_t*)v,n/2);

        if ( n & 0x1 ) {
            hw += hammingWeight(v[n-1]);
        }

        return hw;
    }

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
    int MathTools::hw8( const uint8_t *v , int n ) {

        int hw = hw16((uint16_t*)v,n/2);

        if ( n & 0x1 ) {
            hw += hammingWeight(v[n-1]);
        }

        return hw;
    }

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
    int MathTools::hd64
    ( register const uint64_t *v , register const uint64_t *w ,
      register int n ) {

        register int hd = 0;
        register int i;

        for ( i = 0 ; i < n ; i++ ) {
            hd += hammingWeight(v[i]^w[i]);
        }

        return hd;
    }

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
    int MathTools::hd32
    ( const uint32_t *v , const uint32_t *w , int n ) {

        int hd = hd64((uint64_t*)v,(uint64_t*)w,n/2);

        if ( n & 0x1 ) {
            hd += hammingWeight(v[n-1]^w[n-1]);
        }

        return hd;
    }

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
    int MathTools::hd16
    ( const uint16_t *v , const uint16_t *w , int n ) {

        int hd = hd32((uint32_t*)v,(uint32_t*)w,n/2);

        if ( n & 0x1 ) {
            hd += hammingWeight(v[n-1]^w[n-1]);
        }

        return hd;
    }

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
    int MathTools::hd8
    ( const uint8_t *v , const uint8_t *w , int n ) {

        int hd = hd16((uint16_t*)v,(uint16_t*)w,n/2);

        if ( n & 0x1 ) {
            hd += hammingWeight(v[n-1]^w[n-1]);
        }

        return hd;
    }

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
    void MathTools::mxor64
    ( register uint64_t *out ,
      register const uint64_t *in1 , register const uint64_t *in2 ,
      register int n ) {

        register int i;

        for ( i = 0 ; i < n ; i++ , out++ , in1++ , in2++ ) {
            *out = (*in1) ^ (*in2);
        }
    }

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
    void MathTools::mxor32
    ( uint32_t *out , const uint32_t *in1 , const uint32_t *in2 , int n ) {

        mxor64((uint64_t*)out,(uint64_t*)in1,(uint64_t*)in2,n/2);

        if ( (n & 0x1) != 0 ) {
            out[n-1] = in1[n-1] ^ in2[n-1];
        }
    }

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
    void MathTools::mxor16
    ( uint16_t *out , const uint16_t *in1 , const uint16_t *in2 , int n ) {

        mxor32((uint32_t*)out,(uint32_t*)in1,(uint32_t*)in2,n/2);

        if ( (n & 0x1) != 0 ) {
            out[n-1] = in1[n-1] ^ in2[n-1];
        }
    }

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
    void MathTools::mxor8
    ( uint8_t *out , const uint8_t *in1 , const uint8_t *in2 , int n ) {

        mxor16((uint16_t*)out,(uint16_t*)in1,(uint16_t*)in2,n/2);

        if ( (n & 0x1) != 0 ) {
            out[n-1] = in1[n-1] ^ in2[n-1];
        }
    }


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
    bool MathTools::zeroTest64
    ( register const uint64_t *v , register int n ) {

        register int i;

        for ( i = 0 ; i < n ; i++ , v++ ) {
            if ( *v ) {
                return false;
            }
        }

        return true;
    }

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
    bool MathTools::zeroTest32
    ( const uint32_t *v , int n ) {

        if ( !zeroTest64((uint64_t*)v,n/2) ) {
            return false;
        }

        if ( n&0x1 ) {
            if ( v[n-1] ) {
                return false;
            }
        }

        return true;
    }

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
    bool MathTools::zeroTest16
    ( const uint16_t *v , int n ) {

        if ( !zeroTest32((uint32_t*)v,n/2) ) {
            return false;
        }

        if ( n&0x1 ) {
            if ( v[n-1] ) {
                return false;
            }
        }

        return true;
    }

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
    bool MathTools::zeroTest8
    ( const uint8_t *v , int n ) {

        if ( !zeroTest16((uint16_t*)v,n/2) ) {
            return false;
        }

        if ( n&0x1 ) {
            if ( v[n-1] ) {
                return false;
            }
        }

        return true;
    }


	/**
	 * @brief
	 *            Generates a random 8-bit word.
	 *
	 * @details
	 *            see 'MathTools.h'
	 */
	uint8_t MathTools::rand8( bool tryRandom ) {

		uint8_t c;

#ifdef THIMBLE_USE_DEV_URANDOM
		if ( tryRandom ) {
			static FILE *udev = fopen("/dev/urandom","rb");
			if ( udev != NULL ) {
				c = fgetc(udev);
			} else {
				c = rand()&0xFF;
			}
		} else {
			c = rand()&0xFF;
		}
#else
#ifdef THIMBLE_USE_RAND_S
		if ( tryRandom ) {
			unsigned int value = 0;
			rand_s(&value);
			c = value&0xFF;
		} else {
			c = rand()&0xFF;
		}
#else
		c = rand()&0xFF;
#endif
#endif

		return c;
	}
}



