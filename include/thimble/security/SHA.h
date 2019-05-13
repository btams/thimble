/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file SHA.h
 *
 * @brief
 *            Provides a mechanism for computing the SHA-1 hash value of data.
 *
 * @author Benjamin Tams
 */
#ifndef THIMBLE_SHA_H_
#define THIMBLE_SHA_H_

#include <stdint.h>
#include <cstdlib>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace
 */
namespace thimble {

    /**
     * @brief
     *            Instances of this class represent an environment in
     *            which SHA-1 hash value of data can be computed.
     */
    class THIMBLE_DLL SHA {

    private:

        /**
          * @brief
          *            Contains the padded message and may be reallocated to
          *            ensure it is of sufficient size to perform this. 'size'
          *            always denotes the number of bytes that 'buffer' can
          *            store.
          */
        uint32_t *buffer;

        /**
         * @brief
         *             The number of bytes stored in the array
         *             \link buffer\endlink that forms the message of which
         *             hash value is computed.
         */
        uint64_t L;

        /**
         * @brief
         *             The size in bytes of the array \link buffer\endlink.
         */
        size_t size;

        /**
         * @brief
         *            Will be used as the words of the message schedule.
         */
        uint32_t W[80];


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
        void prepareMessageBits( uint64_t l );

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
        void padMessage( const uint8_t *message , uint64_t n );

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
        void padMessage( const uint32_t *message , uint64_t n );

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
        void hash( uint32_t h[5] );

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
        void hash( uint8_t h[20] );

    public:

        /**
         * @brief
         *            Standard constructor.
         */
        SHA();

        /**
         * @brief
         *            Destructor.
         */
        ~SHA();

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
        void hash( uint32_t h[5] , const uint8_t *message , uint64_t n );

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
        void hash( uint32_t h[5] , const uint32_t *message, uint64_t n );

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
        void hash( uint8_t h[20] , const uint8_t *message , uint64_t n );

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
        void hash( uint8_t h[20] , const uint32_t *message , uint64_t n );
    };
}

#endif /* THIMBLE_SHA_H_ */

