/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file AES.h
 *
 * @brief
 *            Provides mechanisms for encrypting and decrypting data
 *            using the <em>advanced encryption standard</em>.
 *
 * @see <b>National Institute of Standards and Technology (2001).</b>
 *      <i>
 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
 *       (AES).
 *      </i>
 *      <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf" target="_blank">
 *      available online</a>.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_AES_H_
#define THIMBLE_AES_H_

#include <stdint.h>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Generates a 128 bit AES key which can be used
	 *            for data encryption and data decryption.
	 *
	 * @details
	 *            <h2>High-Level Example</h2>
	 *            For an example on how the class can be used, we show
	 *            how it can encrypt a random 64 byte array using the password
	 *            given by the C-string <code>"secret key"</code>.
	 *            Therefore, let
	 *            <pre>
	 *             uint8_t msg[64];
	 *             for ( int i = 0 ; i < 64 ; i++ ) {
	 *                msg[i] = rand()&0xFF;
	 *             }
	 *            </pre>
	 *            be the random 64 byte array. An AES key is constructed
	 *            via
	 *            <pre>
	 *             AES128 aes("secret key");
	 *            </pre>
	 *            which can be used to create the encrypted message
	 *            by
	 *            <pre>
	 *             uint8_t encMsg[64];
	 *
	 *             aes.encrypt(encMsg,msg,64);
	 *            </pre>
	 *            Using the same key, the message can be decrypted by
	 *            <pre>
	 *             uint8_t decMsg[64];
	 *
	 *             aes.decrypt(decMsg,encMsg,64);
	 *            </pre>
	 *            which yields the same byte sequence as given by
	 *            the original message <code>msg</code>. More specifically,
	 *            <pre>
	 *             cout << memcmp(msg,decMsg,64) << endl;
	 *            </pre>
	 *            prints 0 to <code>stdout</code>.
	 *
	 *            <h2>Low-Level Example</h2>
	 *            In the standard specification (see below for a reference)
	 *            an encryption and decryption example with the plain text
	 *            <code>
	 *             00112233445566778899aabbccddeeff
	 *            </code>
	 *            using the 128 bit key
	 *            <code>
	 *             000102030405060708090a0b0c0d0e0f
	 *            </code>
	 *            is given. In this section we show how the example can
	 *            be reproduced using the implementation of the class.
	 *            Let
	 *            <pre>
	 * uint8_t msg[] = {0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88,0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xff};
	 *            </pre>
	 *            be the message and
	 *            <pre>
	 * uint8_t key[] = {0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f};
	 *            </pre>
	 *            be the key. The key is expanded on constructing an instance
	 *            of this class via
	 *            <pre>
	 * AES128 aes(key);
	 *            </pre>
	 *            The message is used to initialize the state of the AES via
	 *            <pre>
	 * memcpy(aes.getState(),msg,16);
	 *            </pre>
	 *            and the state is encrypted by running
	 *            <pre>
	 * aes.cipher();
	 *            </pre>
	 *            Finally the snippet
	 *            <pre>
	 * cout << std::hex;
	 * for ( int j = 0 ; j < 16 ; j++ )
	 *    cout << (int)aes.getState()[j];
	 * cout << std::dec << endl;
	 *            </pre>
	 *            prints the message
	 *            <code>
	 *             69c4e0d86a7b0430d8cdb78070b4c55a
	 *            </code>
	 *            to <code>stdout</code> as reported on page 36 of the standard.
	 *            The method
	 *            <pre>
	 * aes.invCipher();
	 *            </pre>
	 *            causes the state to contain the same byte sequence
	 *            as the original message <code>msg</code>.
	 *
	 * @see <b>National Institute of Standards and Technology (2001).</b>
     *      <i>
     *       FIPS PUB 197: Announcing the Advanced Encryption Standard
     *       (AES).
     *      </i>
     *      <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf" target="_blank">
     *      available online</a>.
	 */
	class THIMBLE_DLL AES128 {

	public:

		/**
		 * @brief
		 *           Creates an AES key using an array of 16 bytes which
		 *           define a 128 bit value.
		 *
		 * @param key
		 *           An array of 16 bytes defining the 128 bit key.
		 *
		 */
		AES128( const uint8_t key[16] );

		/**
		 * @brief
		 *           Derives a 128 bit AES key with the help of a
		 *           C-string defining a password.
		 *
		 * @details
		 *           The constructor computes the 20 byte
         *           \link SHA SHA-1 hash value\endlink of the
		 *           first <code>strlen(password)</code> bytes in
		 *           <code>password</code> and uses the first 16 bytes to
		 *           initialize the 128 bit AES key.
		 *
		 * @param password
		 *           A <code>'\0'</code> terminated C-string defining
		 *           a password.
		 */
        AES128( const char *password = "" );

		/**
		 * @brief
		 *           Copy constructor.
		 *
		 * @details
		 *           Creates a copy of the input AES key.
		 *
		 * @param aes
		 *           AES key of which a copy is constructed.
		 */
		AES128( const AES128 & aes );

		/**
		 * @brief
		 *           Destructor
		 */
		~AES128();

		/**
		 * @brief
		 *           Assignment operator
		 *
		 * @details
		 *           Sets this instance to a copy of the AES key
		 *           specified by <code>aes</code>.
		 *
		 * @param aes
		 *           The AES key assigned to this AES key.
		 *
		 * @return
		 *           A reference to this AES key after assignment.
		 */
		AES128 &operator=(const AES128 & aes );

		/**
		 * @brief
		 *           Encrypt data using this AES key.
		 *
		 * @details
		 *           The input data <code>in</code> is divided into 16 bytes
		 *           (i.e., 128 bit) blocks and encrypted successively in
		 *           <i>
		 * <a href="http://en.wikipedia.org/wiki/Block_cipher_mode_of_operation"
		 * target="_blank">
		 *           cipher-block chaining</a> (with zero initialization vector)
		 *           </i> mode by storing the output 16 byte blocks
		 *           successively to the output <code>out</code>.
		 *
		 *           Note that <code>out</code> and <code>in</code> are
		 *           allowed to be equal
		 *           (i.e., <code>out==in</code> can be <code>true</code>) but should
		 *           not overlap in their remaining first <i>n</i> bytes.
		 *
		 * @param in
		 *           Array of <i>n</i> well-defined bytes containing the
		 *           unencrypted data.
		 *
		 * @param out
		 *           Array that can hold <i>n</i> bytes to which the
		 *           encrypted data is written.
		 *
		 * @param n
		 *           A non-negative integer being a multiple of 16 defining
		 *           the size of the to-be-encrypted data.
		 *
		 * @warning
		 *           If <i>n</i> is not a multiple of 16, an error message
		 *           is printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
		 */
		void encrypt( uint8_t *out , const uint8_t *in , int n );

		/**
		 * @brief
		 *           Decrypt data using this AES key.
		 *
		 * @details
		 *           The input data <code>in</code> is divided into 16 bytes
		 *           (i.e., 128 bit) blocks and decrypted successively by
		 *           assuming that the data was encrypted in
		 *           <i>
		 * <a href="http://en.wikipedia.org/wiki/Block_cipher_mode_of_operation"
		 * target="_blank">
		 *           cipher-block chaining</a>
		 *           </i> mode (with zero initialization vector).
		 *           Each decrypted block is successively written
		 *           to <code>out</code>.
		 *
		 *           Note that <code>out</code> and <code>in</code> are
		 *           allowed to be equal
		 *           (i.e., <code>out==in</code> can be <code>true</code>) but
		 *           should not overlap in their remaining first <i>n</i>
		 *           bytes.
		 *
		 * @param in
		 *           Array of <i>n</i> well-defined bytes containing the
		 *           encrypted data.
		 *
		 * @param out
		 *           Array that can hold <i>n</i> bytes to which the
		 *           decrypted data is written.
		 *
		 * @param n
		 *           A non-negative integer being a multiple of 16 defining
		 *           the size of the to-be-encrypted data.
		 *
		 * @warning
		 *           If <i>n</i> is not a multiple of 16, an error message
		 *           is printed to <code>stderr</code> and the program
         *           exits with status 'EXIT_FAILURE'.
		 */
		void decrypt( uint8_t *out , const uint8_t *in , int n );

		/**
		 * @brief
		 *           Access the 16 bytes of this AES's state.
		 *
		 * @details
		 *           The advanced encryption standard works with
		 *           128 bit blocks of data which are (temporarily) copied
		 *           to the state. Using the AES key, the state is
		 *           <em>shuffled</em> on encryption and decryption to
		 *           produce the encrypted and decrypted data, respectively.
		 *
		 * @return
		 *           Array containing 16 bytes.
		 */
		uint8_t *getState();

		/**
		 * @brief
		 *           Access the 16 bytes of this AES's state
		 *           (constant version).
		 *
		 * @see getState()
		 *
		 * @return
		 *           Constant array of 16 bytes.
		 */
		const uint8_t* getState() const;

		/**
		 * @brief
		 *           <em>Shuffle</em> the state to produce the encrypted
		 *           state.
		 *
		 * @details
		 *           This method performs operations on the initial
		 *           \link getState() state\endlink such that it will
		 *           contain the encrypted \link getState() state\endlink
		 *           after the method has finished its execution.
		 *
		 *           This method is the implementation of the <em>Cipher</em>
		 *           method specified by the standard.
		 *
		 * @see Section 5.1, page 13 in
		 *      <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void cipher();

		/**
		 * @brief
		 *           <em>Unshuffle</em> the state to produce the decrypted
		 *           state.
		 *
		 * @details
		 *           This method performs operations on the encrypted
		 *           \link getState() state\endlink such that it will
		 *           contain the decrypted \link getState() state\endlink
		 *           after the method has finished its execution.
		 *
		 *           This method implements the <em>Inverse Cipher</em> method
		 *           specified by the standard.
		 *
		 * @see Section 5.3, page 20 in
		 *      <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void invCipher();

		/**
		 * @brief
		 *           Assign a new encryption/decryption key to this
		 *           instance.
		 *
		 * @details
		 *           This method is used by the constructors
		 *           \link AES128(const uint8_t[16])\endlink and
		 *           \link AES128(const char *password)\endlink
		 *           to create the AES key. If this instance is sought
		 *           to present another AES key, this method can be used
		 *           for initialization.
		 *
		 *           This method implements the <em>Key Expansion</em>
		 *           procedure specified by the standard.
		 *
		 * @param key
		 *           Array of 16 bytes specifing the 128 bit AES key.
		 *
		 * @see Section 5.2, page 19 in
		 *      <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void keyExpansion( const uint8_t key[16] );

		/**
		 * @brief
		 *           Performs the <em>AddRoundKey</em> transformation
		 *           specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void addRoundKey( int round );

		/**
		 * @brief
		 *            Performs the <em>SubBytes()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void subBytes();

		/**
		 * @brief
		 *            Performs the <em>ShiftRows()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void shiftRows();

		/**
		 * @brief
		 *            Performs the <em>MixColumns()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void mixColumns();

		/**
		 * @brief
		 *            Performs the <em>InvSubBytes()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void invSubBytes();

		/**
		 * @brief
		 *            Performs the <em>InvShiftRows()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void invShiftRows();

		/**
		 * @brief
		 *            Performs the <em>InvMixColumns()</em> transformation
		 *            specified by the standard.
		 *
		 * @see <b>National Institute of Standards and Technology (2001).</b>
		 *      <i>
		 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
		 *       (AES).
		 *      </i>
		 * <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf"
		 *      target="_blank">
		 *      available online</a>.
		 */
		void invMixColumns();

	private:
		/**
		 * @brief
		 *            Holds the AES key after expansion.
		 *
		 * @details
		 *            A 128-bit AES key is given by an 16 byte array. This
		 *            array is expanded to one-dimensional array of 44 words
		 *            of length 32 bit using the
		 *            \link keyExpansion(const uint8_t[16])\endlink method
		 *            which is called automatically by the constructors
		 *            \link AES128(const uint8_t[16])\endlink and
		 *            \link AES128(const char *password)\endlink
		 *            on initialization.
		 */
		uint32_t w[4*(10+1)];

		/**
		 * @brief
		 *            The internal state on which the AES algorithms operate.
		 *
		 * @details
		 *            The state is a two-dimensional 4*4 array of bytes.
		 *            The array \link state\endlink is to be interpreted as
		 *            a matrix in column-major order, i.e.,
		 *            \link state\endlink represents the matrix
		 *            \f[
		 *             \left(\begin{array}{cccc}
		 *              state[0]&state[4]&state[8]&state[12]\\
		 *              state[1]&state[5]&state[8]&state[13]\\
		 *              state[2]&state[6]&state[8]&state[14]\\
		 *              state[3]&state[7]&state[8]&state[15]
		 *             \end{array}\right).
		 *            \f]
		 */
		uint8_t state[16];
	};

}


#endif /* THIMBLE_AES_H_ */
