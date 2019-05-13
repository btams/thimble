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
 * @file SHA.cpp
 *
 * @brief
 *            This file implements the functionalities provided by 'SHA.h'
 *            which provides a mechanism for computing the SHA-1 hash value of
 *            data.
 *
 * @details
 *            see 'SHA.h'
 *
 * @author Benjamin Tams
 */
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/security/SHA.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

    /**
     * @brief
     *            Standard constructor.
     */
    SHA::SHA() {
        this->buffer = NULL;
        this->L = 0;
        this->size = 0;
    }

    /**
     * @brief
     *            Destructor.
     */
    SHA::~SHA() {
        free(this->buffer);
    }

    /**
     * @brief
     *            Ensures that this %SHA object can buffer a message of
     *            specified bit length.
     *
     * @details
     *            The value \link L\endlink will be updated indicating
     *            the number of bytes that the
     *            array \link buffer\endlink must be able to store in
     *            order to hash a message of the specified bit length.
     *            If needed, the array \link buffer\endlink will be
     *            reallocated in which case the
     *            integer \link size\endlink is increased.
     *
     * @param l
     *            The number of bits of a message that this %SHA object
     *            must be able to hash.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void SHA::prepareMessageBits( uint64_t l ) {

        // The message will be padded by '100...0' where the number of zeros
        // is equals 'k' which is the smallest non-negative solution to
        // the equation 'l+k+1=448 mod 512'
        uint64_t k;
        if ( (l+1)%512 <= 448 ) {
            k = 448-(l+1)%512;
        } else {
            k = 448+(512-(l+1)%512);
        }

        // The padded message's length will be 'l+k+1+64'. Thus, the number
        // of blocks that are processed is 'N=(l+k+1+64)/512'.
        uint64_t N = (l+k+1+64)/512;

        // Check whether the buffer for the padded message has to be reallocated
        if ( (size_t)(N*64) > size ) {
            size = (size_t)(N*64);
            if ( buffer == NULL )
                buffer =
                   (uint32_t*)malloc( size * sizeof(uint8_t) );
            else
                buffer =
                   (uint32_t*)realloc( buffer , size * sizeof(uint8_t) );
            if ( buffer == NULL ) {
                cerr << __FILE__ << "(" << __LINE__ << "): out of memory."
                        << endl;
                exit(EXIT_FAILURE);
            }
        }

        this->L = N * 16;
    }

    /**
     * @brief
     *            Copies message of bytes to \link buffer\endlink, appends
     *            the bit 1 to the buffer and pads the bit sequence by
     *            0 bits until the bit-length of the buffered message
     *            becomes a multiple of 512.
     *
     * @details
     *            The value \link L\endlink will hold the number of
     *            padded bytes stored in \link buffer\endlink.
     *
     * @param message
     *            An array containing <code>n</code> 8-bit integers,
     *            i.e., bytes, of type <code>uint8_t</code>.
     *
     * @param n
     *            Number of bytes stored in <code>message</code>.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If <code>message</code> does not hold at least
     *            <code>n</code> valid integers of type
     *            <code>uint8_t</code>, the behavior of this
     *            method is undocumented.
     */
    void SHA::padMessage( const uint8_t *message , uint64_t n ) {

        // Length of message in bits
        uint64_t l = 8*n;

        // Ensure that this SHA's buffer is large enough to hash
        // a message of 'l' bits.
        prepareMessageBits(l);

        memset(this->buffer,0,(size_t)(L*sizeof(uint32_t)));

        // Convert the message by intepreting the message's bytes as big-endian
        // 32 bit words.
        for ( uint64_t i = 0 , k = 0 ; i < n ; i += 4 , k++ ) {
            this->buffer[k] = ((uint32_t)(message[i]))<<24;
            if ( i+1 < n ) {
                this->buffer[k] += ((uint32_t)(message[i+1]))<<16;
                if ( i+2 < n ) {
                    this->buffer[k] += ((uint32_t)(message[i+2]))<<8;
                    if ( i+3 < n ) {
                        this->buffer[k] += ((uint32_t)(message[i+3]));
                    }
                }
            }
        }

        // Append the bit '1' to the end of the message
        if ( n%4 == 0 ) {
            this->buffer[n/4] = ((uint32_t)1)<<31;
        } else {
            this->buffer[(n-1)/4] += ((uint32_t)1)<<(31-8*(n%4));
        }

        // Append the 64-bit block that is equal to the number 'l' expressed
        // using a binary representation (see Section 5.1.1. in FIPS 180-2).
        this->buffer[this->L-1] |= l&0xFFFFFFFF;
        this->buffer[this->L-2] |= l>>32;
    }

    /**
     * @brief
     *            Copies a message of 32-bit integers
     *            to \link buffer\endlink, appends the bit 1 to the buffer
     *            and pads the bit sequence by 0 bits until the bit-length
     *            of the buffered message becomes a multiple of 512.
     *
     * @details
     *            The value \link L\endlink will hold the number of
     *            padded bytes stored in \link buffer\endlink.
     *
     * @param message
     *            An array containing <code>n</code> 32-bit integers,
     *            i.e., bytes, of type <code>uint8_t</code>.
     *
     * @param n
     *            Number of 32-bit integers stored in
     *            <code>message</code>.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If <code>message</code> does not hold at least
     *            <code>n</code> valid integers of type
     *            <code>uint32_t</code>, the behavior of this
     *            method is undocumented.
     */
    void SHA::padMessage( const uint32_t *message , uint64_t n ) {

        // Length of message in bits
        uint64_t l = 32*n;

        // Ensure that this SHA's buffer is large enough to hash
        // a message of 'l' bits.
        prepareMessageBits(l);

        // Initialize the buffer and ...
        memcpy(this->buffer,message,(size_t)(n*sizeof(uint32_t)));
        // ... pad the message.
        memset(this->buffer+n,0,(size_t)((L-n)*sizeof(uint32_t)));

        // Append the bit '1' to the end of the message
        this->buffer[n] = ((uint32_t)1)<<31;

        // Append the 64-bit block that is equal to the number 'l' expressed
        // using a binary representation (see Section 5.1.1. in FIPS 180-2).
        this->buffer[this->L-1] |= l&0xFFFFFFFF;
        this->buffer[this->L-2] |= l>>32;
    }

    /**
     * @brief
     *            Circular left shift operation.
     *
     * @details
     *            The operation is defined as
     *            \f$rotl^n(x)=(x << n) or (x >> (w - n))\$f where
     *            \f$w\f$ denotes the width of the integer \f$x\f$, i.e.
     *            (32 bits in our case) and << and >> denote the binary left
     *            shift and right shift operation, respectively.
     *
     * @param x
     *            The integer that is shifted.
     *
     * @param n
     *            The number of of bits to be shifted.
     *
     * @return
     *            The circular left shift of the 32-bit integer <code>x</code>
     *            by <code>n</code> bits.
     */
    inline static uint32_t rotl32( register uint32_t x , int n ) {
        return (x << n) | (x >> (32 - n));
    }

    /**
     * @brief
     *            Computes the hash of the padded message stored by this
     *            %SHA object and stores the result in the specified
     *            output array.
     *
     * @param h
     *            An array that is capable to hold 5 integers of type
     *            <code>uint32_t</code> to which the hash value
     *            of the buffered message is stored.
     */
    void SHA::hash( uint32_t h[5] ) {

        uint32_t H0 , H1 , H2 , H3 , H4 , a , b , c , d , e , f , K , T;

        // Initial values (see Section 5.3.1 in FIPS 180-2)
        H0 = 0x67452301;
        H1 = 0xEFCDAB89;
        H2 = 0x98BADCFE;
        H3 = 0x10325476;
        H4 = 0xC3D2E1F0;

        // ******************************************************************
        // *** BEGIN: SHA-1 Hash Computation (Section 6.2.1 in FIPS 180-2) **
        // ******************************************************************

        // Iterate over the message block
        for ( uint64_t j = 0 ; j < L ; j += 16 ) {

            // Step 1. Prepare the message schedule
            memcpy(W,buffer+j,16*sizeof(uint32_t));
            for ( int i = 16 ; i < 80 ; i++ ) {
                W[i] = rotl32( W[i-3] ^ W[i-8] ^ W[i-14] ^ W[i-16] , 1 );
            }

            // Step 2. Initialize the five working variables with the
            // the (i-1)-st hash value
            a = H0;
            b = H1;
            c = H2;
            d = H3;
            e = H4;

            // Step 3.
            for ( int t = 0 ; t < 80 ; t++ ) {

                if ( t < 20 ) {
                    f = (b&c) ^ ((~b)&d);
                    K = 0x5A827999;
                } else if ( t < 40 ) {
                    f = b ^ c ^ d;
                    K = 0x6ED9EBA1;
                } else if ( t < 60 ) {
                    f = (b&c) ^ (b&d) ^ (c&d);
                    K = 0x8F1BBCDC;
                } else {
                    f = b ^  c ^ d;
                    K = 0xCA62C1D6;
                }

                T = rotl32(a,5) + f + e + K + W[t];
                e = d;
                d = c;
                c = rotl32(b,30);
                b = a;
                a = T;
            }

            // Step 4. Compute the i-th intermediate hash value
            H0 += a;
            H1 += b;
            H2 += c;
            H3 += d;
            H4 += e;
        }

        // ******************************************************************
        // **** END: SHA-1 Hash Computation (Section 6.2.1 in FIPS 180-2) ***
        // ******************************************************************

        // Output
        h[0] = H0;
        h[1] = H1;
        h[2] = H2;
        h[3] = H3;
        h[4] = H4;
    }

    /**
     * @brief
     *            Computes the hash of the padded message stored by this
     *            %SHA object and stores the result in the specified
     *            output array.
     *
     * @param h
     *            An array that is capable to hold 20 integers of type
     *            <code>uint8_t</code> to which the hash value
     *            of the buffered message is stored.
     */
    void SHA::hash( uint8_t h[20] ) {

        uint32_t tmp[5];

        hash(tmp);

        for ( int i = 0 ; i < 5 ; i++ ) {
            for ( int j = 3 ; j >= 0 ; j-- ) {
                h[i*4+j] = tmp[i] & 0xFF;
                tmp[i] >>= 8;
            }
        }
    }

    /**
     * @brief
     *            Computes the SHA-1 hash value of a message.
     *
     * @details
     *            The hash value of the first <code>n</code> bytes
     *            contained in <code>message</code> is computed and
     *            stored in <code>h</code>.
     *
     * @param h
     *            Will contain the SHA-1 hash value of the first <code>n
     *            </code> bytes in <code>message</code>.
     *
     * @param message
     *            Contains at least <code>n</code> valid bytes of which
     *            the SHA hash value is computed.
     *
     * @param n
     *            Number of valid bytes contained in <code>message</code>
     *
     * @warning
     *            If not sufficient memory could be provided to compute
     *            the result an error message is printed to
     *            <code>stderr</code> and the program exits with status
     *            -1.
     *
     * @warning
     *            If <code>message</code> contains no <code>n</code> valid
     *            bytes the method runs into undocumented behavior.
     */
     void SHA::hash( uint32_t h[5] , const uint8_t *message , uint64_t n ) {

         // Pad the message in the buffer ...
         padMessage(message,n);

         // ... and hash from the buffer
         hash(h);
     }

     /**
      * @brief
      *            Computes the SHA-1 hash value of a message.
      *
      * @details
      *            The hash value of the first <code>n</code> integers
      *            of type <code>uint32_t</code> contained in
      *            <code>message</code> are computed and stored in
      *            <code>h</code>.
      *
      * @param h
      *            Will contain the SHA-1 hash value of the first <code>n
      *            </code> 32-bit integers in <code>message</code>.
      *
      * @param message
      *            Contains at least <code>n</code> valid integers of
      *            type <code>uint32_t</code> of which the %SHA hash value
      *            is computed.
      *
      * @param n
      *            Number of valid 32-bit integers contained in
      *            <code>message</code>
      *
      * @warning
      *            If not sufficient memory could be provided to compute
      *            the result an error message is printed to
      *            <code>stderr</code> and the program exits with status
      *            -1.
      *
      * @warning
      *            If <code>message</code> contains no <code>n</code> valid
      *            bytes the method runs into undocumented behavior.
      */
     void SHA::hash( uint32_t h[5] , const uint32_t *message, uint64_t n ) {

         // Pad the message in the buffer ...
         padMessage(message,n);

         // ... and hash from the buffer
         hash(h);
     }

     /**
      * @brief
      *            Computes the SHA-1 hash value of a message.
      *
      * @details
      *            The hash value of the first <code>n</code> bytes
      *            contained in <code>message</code> is computed and
      *            stored in <code>h</code>.
      *
      * @param h
      *            Will contain the SHA-1 hash value of the first <code>n
      *            </code> bytes in <code>message</code>.
      *
      * @param message
      *            Contains at least <code>n</code> valid bytes of which
      *            the SHA hash value is computed.
      *
      * @param n
      *            Number of valid bytes contained in <code>message</code>
      *
      * @warning
      *            If not sufficient memory could be provided to compute
      *            the result an error message is printed to
      *            <code>stderr</code> and the program exits with status
      *            -1.
      *
      * @warning
      *            If <code>message</code> contains no <code>n</code> valid
      *            bytes the method runs into undocumented behavior.
      */
     void SHA::hash( uint8_t h[20] , const uint8_t *message , uint64_t n ) {

         padMessage(message,n);

         hash(h);
     }

     /**
      * @brief
      *            Computes the SHA-1 hash value of a message.
      *
      * @details
      *            The hash value of the first <code>n</code> integers
      *            of type <code>uint32_t</code> contained in
      *            <code>message</code> are computed and stored in
      *            <code>h</code>.
      *
      * @param h
      *            Will contain the SHA-1 hash value of the first <code>n
      *            </code> 32-bit integers in <code>message</code>.
      *
      * @param message
      *            Contains at least <code>n</code> valid integers of
      *            type <code>uint32_t</code> of which the %SHA hash value
      *            is computed.
      *
      * @param n
      *            Number of valid 32-bit integers contained in
      *            <code>message</code>
      *
      * @warning
      *            If not sufficient memory could be provided to compute
      *            the result an error message is printed to
      *            <code>stderr</code> and the program exits with status
      *            -1.
      *
      * @warning
      *            If <code>message</code> contains no <code>n</code> valid
      *            bytes the method runs into undocumented behavior.
      */
     void SHA::hash( uint8_t h[20] , const uint32_t *message , uint64_t n ) {

         padMessage(message,n);

         hash(h);
     }

}





