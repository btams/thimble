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
 * @file AES.cpp
 *
 * @brief
 *            Implements mechanisms for encrypting and decrypting data
 *            with the <em>advanced encryption standard</em> that are
 *            provided by the 'ProtectedMinutiaeTemplate.h' header.
 *
 * @author Benjamin Tams
 */
#include <stdint.h>
#include <cstring>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryPolynomial.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/security/SHA.h>
#include <thimble/security/AES.h>


using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Finite field in which internal AES computations are performed.
	 *
	 * @details
	 *           The standard works with bytes each representing an element
	 *           in the finite field defined by the polynomial
	 *           \f$x^8+x^4+x^3+x+1\f$ which is represented by the integer
	 *           \f$2^8+2^4+2^3+2+1=283\f$.
	 *
	 * @see Section 4.2, Equation (4.1), page 10 in
	 *      <b>National Institute of Standards and Technology (2001).</b>
	 *      <i>
	 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
	 *       (AES).
	 *      </i>
	 *      <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf" target="_blank">
	 *      available online</a>.
	 */
	static SmallBinaryField _GF(SmallBinaryPolynomial(283));

	/**
	 * @brief
	 *           Substitution values for the byte <code>xy</code> in
	 *           (hexadecimal format).
	 *
	 * @see Figure 7 in
	 *      <b>National Institute of Standards and Technology (2001).</b>
	 *      <i>
	 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
	 *       (AES).
	 *      </i>
	 *      <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf" target="_blank">
	 *      available online</a>.
	 */
	static const uint8_t _SBOX[256] = {
	0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
	0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
	0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
	0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
	0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
	0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
	0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
	0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
	0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
	0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
	0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
	0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
	0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
	0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
	0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
	0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16};

	/**
	 * @brief
	 *           Inverse substitution values for the byte <code>xy</code> in
	 *           (hexadecimal format).
	 *
	 * @see Figure 14 in
	 *      <b>National Institute of Standards and Technology (2001).</b>
	 *      <i>
	 *       FIPS PUB 197: Announcing the Advanced Encryption Standard
	 *       (AES).
	 *      </i>
	 *      <a href="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf" target="_blank">
	 *      available online</a>.
	 */
	static const uint8_t _ISBOX[256] = {
	0x52, 0x09, 0x6A, 0xD5, 0x30, 0x36, 0xA5, 0x38, 0xBF, 0x40, 0xA3, 0x9E, 0x81, 0xF3, 0xD7, 0xFB,
	0x7C, 0xE3, 0x39, 0x82, 0x9B, 0x2F, 0xFF, 0x87, 0x34, 0x8E, 0x43, 0x44, 0xC4, 0xDE, 0xE9, 0xCB,
	0x54, 0x7B, 0x94, 0x32, 0xA6, 0xC2, 0x23, 0x3D, 0xEE, 0x4C, 0x95, 0x0B, 0x42, 0xFA, 0xC3, 0x4E,
	0x08, 0x2E, 0xA1, 0x66, 0x28, 0xD9, 0x24, 0xB2, 0x76, 0x5B, 0xA2, 0x49, 0x6D, 0x8B, 0xD1, 0x25,
	0x72, 0xF8, 0xF6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xD4, 0xA4, 0x5C, 0xCC, 0x5D, 0x65, 0xB6, 0x92,
	0x6C, 0x70, 0x48, 0x50, 0xFD, 0xED, 0xB9, 0xDA, 0x5E, 0x15, 0x46, 0x57, 0xA7, 0x8D, 0x9D, 0x84,
	0x90, 0xD8, 0xAB, 0x00, 0x8C, 0xBC, 0xD3, 0x0A, 0xF7, 0xE4, 0x58, 0x05, 0xB8, 0xB3, 0x45, 0x06,
	0xD0, 0x2C, 0x1E, 0x8F, 0xCA, 0x3F, 0x0F, 0x02, 0xC1, 0xAF, 0xBD, 0x03, 0x01, 0x13, 0x8A, 0x6B,
	0x3A, 0x91, 0x11, 0x41, 0x4F, 0x67, 0xDC, 0xEA, 0x97, 0xF2, 0xCF, 0xCE, 0xF0, 0xB4, 0xE6, 0x73,
	0x96, 0xAC, 0x74, 0x22, 0xE7, 0xAD, 0x35, 0x85, 0xE2, 0xF9, 0x37, 0xE8, 0x1C, 0x75, 0xDF, 0x6E,
	0x47, 0xF1, 0x1A, 0x71, 0x1D, 0x29, 0xC5, 0x89, 0x6F, 0xB7, 0x62, 0x0E, 0xAA, 0x18, 0xBE, 0x1B,
	0xFC, 0x56, 0x3E, 0x4B, 0xC6, 0xD2, 0x79, 0x20, 0x9A, 0xDB, 0xC0, 0xFE, 0x78, 0xCD, 0x5A, 0xF4,
	0x1F, 0xDD, 0xA8, 0x33, 0x88, 0x07, 0xC7, 0x31, 0xB1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xEC, 0x5F,
	0x60, 0x51, 0x7F, 0xA9, 0x19, 0xB5, 0x4A, 0x0D, 0x2D, 0xE5, 0x7A, 0x9F, 0x93, 0xC9, 0x9C, 0xEF,
	0xA0, 0xE0, 0x3B, 0x4D, 0xAE, 0x2A, 0xF5, 0xB0, 0xC8, 0xEB, 0xBB, 0x3C, 0x83, 0x53, 0x99, 0x61,
	0x17, 0x2B, 0x04, 0x7E, 0xBA, 0x77, 0xD6, 0x26, 0xE1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0C, 0x7D};

	/**
	 * @brief
	 *            The round constant word array,
	 */
	static const uint32_t _RCON[] = {
	0x1000000 ,0x2000000 ,0x4000000 ,0x8000000 ,
	0x10000000,0x20000000,0x40000000,0x80000000,
	0x1b000000,0x36000000
	};

	/**
	 * @brief
	 *            Concatenates four bytes to an unsigned 32 bit word
	 *            (Big-endian).
	 *
	 * @param b0
	 *            unsigned byte
	 *
	 * @param b1
	 *            unsigned byte
	 *
	 * @param b2
	 *            unsigned byte
	 *
	 * @param b3
	 *            unsigned byte
	 *
	 * @return \f$b3+b2\cdot 256+b1\cdot 256^2+b0\cdot 256^3\f$
	 */
	static uint32_t toWord( uint8_t b0 , uint8_t b1 , uint8_t b2 , uint8_t b3 ) {

		uint32_t w;

		w  = (uint32_t)b0; w <<= 8;
		w += (uint32_t)b1; w <<= 8;
		w += (uint32_t)b2; w <<= 8;
		w += (uint32_t)b3;

		return w;
	}

	/**
	 * @brief
	 *           Splits an unsigned 32-bit word into four bytes.
	 *
	 * @param w
	 *           Input unsigned 32-bit word.
	 *
	 * @param b  Output byte array such that
	 *           \f$w=b[3]+b[2]\cdot 256+b[1]\cdot 256^2+b[0]\cdot 256^3\f$
	 */
	static void toBytes( uint8_t b[4] , uint32_t w ) {
		b[3] = (uint8_t)(w & 0xFF); w >>= 8;
		b[2] = (uint8_t)(w & 0xFF); w >>= 8;
		b[1] = (uint8_t)(w & 0xFF); w >>= 8;
		b[0] = (uint8_t)w;
	}

	/**
	 * @brief
	 *           Substitutes a single byte using the
	 *           \link _SBOX S-box\endlink.
	 *
	 * @param b
	 *           Input byte.
	 *
	 * @return
	 *           Substitution of <i>b</i>.
	 */
	inline static uint8_t SubByte( uint8_t b ) {
		return _SBOX[b];
	}

	/**
	 * @brief
	 *           Performs the inverse substitution of a byte using the
	 *           \link _ISBOX inverse S-box\endlink.
	 *
	 * @details
	 *           This function is inverse to \link SubByte(uint8_t)\endlink.
	 *           More specifically,
	 *           <pre>
	 *            b == InvSubByte(SubByte(b))
	 *           </pre>
	 *           holds <code>true</code> for each <code>uint8_t</code>.
	 *
	 * @param b
	 *           Input byte.
	 *
	 * @return
	 *           Inverse substitution of <i>b</i>.
	 */
	inline static uint8_t InvSubByte( uint8_t b ) {
		return _ISBOX[b];
	}

	/**
	 * @brief
	 *           Returns an unsigned 32-bit word whose four bytes
	 *           are the substituted bytes of the input word.
	 *
	 * @param w
	 *           Input word.
	 *
	 * @return
	 *           \f$b_3'+b_2'\cdot 256+b_1'\cdot 256^2+b_0'\cdot 256^3\f$
	 *           where \f$b_j'=SubByte(b_j)\f$ and
	 *           \$fw=b_3+b_2\cdot 256+b_1\cdot 256^2+b_0\cdot 256^3\f$.
	 */
	static uint32_t SubWord( uint32_t w ) {

		uint8_t b[4];

		toBytes(b,w);

		b[0] = SubByte(b[0]);
		b[1] = SubByte(b[1]);
		b[2] = SubByte(b[2]);
		b[3] = SubByte(b[3]);

		return toWord(b[0],b[1],b[2],b[3]);
	}

	/**
	 * @brief
	 *            Function used in the key expansion routine that takes a
	 *            four-byte word and performs a cyclic permutation.
	 *
	 * @param w
	 *            Input four-byte word
	 *            \f$b_3+b_2\cdot 256+b_1\cdot b_1\cdot 256^2+b_0\cdot 256^3\f$.
	 *
	 * @return
	 *            Cylic permutation of <i>w</i>, i.e.,
	 *            \f$b_0+b_3\cdot 256+b_2\cdot 256^2+b_1\cdot 256^3\f$.
	 */
	static uint32_t RotWord( uint32_t w ) {

		uint8_t b[4];

		toBytes(b,w);

		return toWord(b[1],b[2],b[3],b[0]);
	}

	/**
	 * @brief
	 *           Creates an AES key using an array of 16 bytes which
	 *           define a 128 bit value.
	 *
	 * @param key
	 *           An array of 16 bytes defining the 128 bit key.
	 *
	 */
	AES128::AES128( const uint8_t key[16] ) {
		// Expand the key.
		keyExpansion(key);

		// Initialize the state by zeros to avoid segmentation faults.
		memset(this->state,0,16);
	}

    /**
     * @brief
     *           Derives a 128 bit AES key with the help of a
     *           C-string defining a password.
     *
     * @details
     *           The constructor computes the 20 byte
     *           \link SHA1 SHA-1 hash value\endlink of the
     *           first <code>strlen(password)</code> bytes in
     *           <code>password</code> and uses the first 16 bytes to
     *           initialize the 128 bit AES key.
     *
     * @param password
     *           A <code>'\0'</code> terminated C-string defining
     *           a password.
     */
    AES128::AES128( const char *password ) {

        int n;

		n = strlen(password);

		uint8_t key[16];
		uint8_t hash[20];

		// Hash the password.
        SHA sha;
        sha.hash(hash,(uint8_t*)password,n);

		// Get the first 128 bits for ...
		memcpy(key,hash,16);
		// ... key expansion.
		keyExpansion(key);

		// Initialize the state by zeros to avoid segmentation faults.
		memset(this->state,0,16);
	}

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
	AES128::AES128( const AES128 & aes ) {
		memcpy(this->w    ,aes.w    ,44*sizeof(uint32_t));
		memcpy(this->state,aes.state,16);
	}

	/**
	 * @brief
	 *           Destructor
	 */
	AES128::~AES128() {

		// Just overwrite any data.
		memset(this->w    ,0,44);
		memset(this->state,0,16);
	}

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
	AES128 & AES128::operator=(const AES128 & aes ) {
		memcpy(this->w    ,aes.w    ,44);
		memcpy(this->state,aes.state,16);
		return *this;
	}

	/**
	 * @brief
	 *          XOR the elements of two arrays both containing
	 *          16 bytes.
	 *
	 * @details
	 *          XORs the bytes in <code>out</code> with the
	 *          bytes in <code>in</code> and stores the results
	 *          in <code>out</code>.
	 *
	 *          The method is used to realize the chaining block
	 *          cipher mode on encryption and decryption.
	 *
	 * @param in
	 *          Input array containing 16 bytes.
	 *
	 * @param out
	 *          Input array containing 16 bytes; the array also serves
	 *          as the output.
	 */
    static void mxor8( uint8_t out[16] , const uint8_t in[16]  ) {

		// Interpret the bytes as two 64 bit word which should
		// increase execution on 64 bit machines.

		register uint64_t *_out , *_in;

		_out = (uint64_t*)out;
		_in = (uint64_t*)in;

		*_out ^= *_in;
		++_out; ++_in;

		*_out ^= *_in;
    }

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
	void AES128::encrypt( uint8_t *out , const uint8_t *in , int n ) {

		int i;

		// Ensure data consists of 128 bit blocks.
		if ( n % 16 != 0 ) {
			cerr << "AES128::encrypt: "
				 << "Message must consist of 128 bit blocks." << endl;
            exit(EXIT_FAILURE);
		}

		if ( n <= 0 ) { // Nothing to do
			return;
		}

		// Copy the first 128 bits to the state.
		memcpy(this->state,in,16);

		// Encrypt the first 128 bit and ...
		cipher();

		// ... store the encryption to the output; this is equivalent to
		// perform the encryption of the first block with a zero
		// initialization vector when accounting for the
		// cipher-block chaining mode.
		memcpy(out,this->state,16);

		// Iteration over the remaining blocks
		for ( i = 16 ; i < n ; i += 16 ) {

			in  += 16;
			out += 16;

			// As the state contains the encryption of the previous blocks,
			// we use the state for xoring with the current block, thereby
			// realizing the cipher-block chaining mode.
			mxor8(this->state,in);

			// Cipher the current block and ...
			cipher();

			// ... store the encrpytion to the output.
			memcpy(out,this->state,16);
		}
	}

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
	void AES128::decrypt( uint8_t *out , const uint8_t *in , int n ) {

		int i;

		// These arrays are used to temporarily make backups of encrypted
		// blocks such that we can account for the cipher-block chaining mode
		// even if 'out==in' is 'true'
		uint8_t tmp1[16] , tmp2[16];

		// Ensure data consists of 128 bit blocks.
		if ( n % 16 != 0 ) {
			cerr << "AES128::decrypt: "
				 << "Message must consist of 128 bit blocks." << endl;
            exit(EXIT_FAILURE);
		}

		if ( n <= 0 ) { // Nothing to do
			return;
		}

		// Make a copy of the first 128 bit bytes to account for the
		// cipher-block chaining mode.
		memcpy(tmp1,in,16);

		// Decrypt the first block and ...
		memcpy(this->state,in,16);
		invCipher();
		// ... store the decryption to the output.
		memcpy(out,this->state,16);

		// Iteration over the remaining blocks
		for ( i = 16 ; i < n ; i += 16 ) {

			in  += 16;
			out += 16;

			// Make a copy of the current block to account for the
			// cipher-block chaining mode; note that 'tmp1' still
			// contains the encryption of the previous block.
			memcpy(tmp2,in,16);

			// Decrypt the current block and ...
			memcpy(this->state,in,16);
			invCipher();

			// ... xor with the previous encrypted block to
			// reverse the cipher-block chaining mode.
			memcpy(out,tmp1,16);
			mxor8(out,this->state);

			// Copy the encrypted block to 'tmp1' such that on
			// next iteration, 'tmp1' contains the encryption of the
			// previous block to account for the cipher-block chaining
			// mode.
			memcpy(tmp1,tmp2,16);
		}
	}

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
	uint8_t *AES128::getState() {
		return this->state;
	}

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
	const uint8_t *AES128::getState() const {
		return this->state;
	}

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
	void AES128::cipher() {

		// The following is in accordance with the representation
		// in Figure 5 on page 15 in the standard

		addRoundKey(0);

		for ( int round = 1 ; round < 10 ; round++ ) {
			subBytes();
			shiftRows();
			mixColumns();
			addRoundKey(round);
		}

		subBytes();
		shiftRows();
		addRoundKey(10);
	}

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
	void AES128::invCipher() {

		// The following is in accordance with the presentation
		// of Figure 12 on page 21 of the standard.

		addRoundKey(10);

		for ( int round = 9 ; round >= 1 ; round-- ) {
			invShiftRows();
			invSubBytes();
			addRoundKey(round);
			invMixColumns();
		}

		invShiftRows();
		invSubBytes();
		addRoundKey(0);
	}

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
	void AES128::keyExpansion( const uint8_t key[16] ) {

		// Essentially, the following is in accordance with the presentation
		// of Figure 11 on page 20 of the standard.

		uint32_t temp;
		int i;

		// Unroll the loop constantly performing four iterations in the case
		// of a 128 bit key
		this->w[0] = toWord(key[0],key[1],key[2],key[3]);
		this->w[1] = toWord(key[4],key[5],key[6],key[7]);
		this->w[2] = toWord(key[8],key[9],key[10],key[11]);
		this->w[3] = toWord(key[12],key[13],key[14],key[15]);

		for ( i = 4 ; i < 4*(10+1) ; i++ ) {

			temp = this->w[i-1];

			if ( (i & 0x3) == 0 ) {
				temp = SubWord(RotWord(temp)) ^ _RCON[(i>>2)-1];
			}

			this->w[i] = this->w[i-4] ^ temp;
		}
	}

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
	void AES128::addRoundKey( int round ) {

		uint8_t s[4];

		// Iteration over the columns 'c=0,...,Nb-1=3'
		for ( int c = 0 ; c < 4 ; c++ ) {

			memcpy(s,this->state+4*c,4);

			// Produces the four bytes s_{0,c},s_{1,c},s_{2,c},s_{3,c} as in
			// Equation (5.7) on page 18 of the standard.
			toBytes(s,toWord(s[0],s[1],s[2],s[3])^this->w[round*4+c]);

			memcpy(this->state+4*c,s,4);
		}

	}

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
	void AES128::subBytes() {

		// Substitute all bytes of state using '_SBOX'
		for ( int i = 0 ; i < 16 ; i++ ) {
			this->state[i] = SubByte(this->state[i]);
		}
	}

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
	void AES128::shiftRows() {

		// See Figure 8 on page 17 in the standard.

		uint8_t s[4];

		// Iteration over the last three rows of the state
		for ( int r = 1 ; r < 4 ; r++ ) {

			// Make 's' to contain the shifted row which depends on the
			// row index 'r'
			s[0] = this->state[4*((r+0)%4)+r];
			s[1] = this->state[4*((r+1)%4)+r];
			s[2] = this->state[4*((r+2)%4)+r];
			s[3] = this->state[4*((r+3)%4)+r];

			// Copy 's' to the rows
			this->state[4*0+r] = s[0];
			this->state[4*1+r] = s[1];
			this->state[4*2+r] = s[2];
			this->state[4*3+r] = s[3];
		}
	}

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
	void AES128::mixColumns() {

		// Perform the matrix multiplication as specified by Equation (5.6)
		// on page 18 in the standard.

		uint8_t s[4];

		// Iterate over the columns and ...
		for ( int c = 0 ; c < 4 ; c++ ) {

			// ... multiply column-wise with the matrix by ...

			// ... copying the column to a temporary array ...
			memcpy(s,this->state+4*c,4);

			// ... and updating the column by a matrix-vector
			// multiplication.
			this->state[4*c+0] =
					(uint8_t)(_GF.mul(s[0],0x02)^_GF.mul(s[1],0x03)^s[2]^s[3]);
			this->state[4*c+1] =
					(uint8_t)(_GF.mul(s[1],0x02)^_GF.mul(s[2],0x03)^s[0]^s[3]);
			this->state[4*c+2] =
					(uint8_t)(_GF.mul(s[2],0x02)^_GF.mul(s[3],0x03)^s[0]^s[1]);
			this->state[4*c+3] =
					(uint8_t)(_GF.mul(s[3],0x02)^_GF.mul(s[0],0x03)^s[1]^s[2]);
		}
	}

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
	void AES128::invSubBytes() {

		// Substitute all bytes of state using '_ISBOX'
		for ( int i = 0 ; i < 16 ; i++ ) {
			this->state[i] = InvSubByte(this->state[i]);
		}
	}

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
	void AES128::invShiftRows() {

		// See Figure 13 on page 22 of the standard specification

		uint8_t s[4];

		// Iteration over the last three rows of the state.
		for ( int r = 1 ; r < 4 ; r++ ) {

			// Make a copy of the row and ...
			s[0] = this->state[4*0+r];
			s[1] = this->state[4*1+r];
			s[2] = this->state[4*2+r];
			s[3] = this->state[4*3+r];

			// ... update the row of the inverse shift.
			this->state[4*((r+0)%4)+r] = s[0];
			this->state[4*((r+1)%4)+r] = s[1];
			this->state[4*((r+2)%4)+r] = s[2];
			this->state[4*((r+3)%4)+r] = s[3];
		}
	}

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
	void AES128::invMixColumns() {

		// Perform the matrix multiplication as specified by Equation (5.10)
		// on page 23 in the standard.

		uint8_t s[4];

		// Iterate over the columns and ...
		for ( int c = 0 ; c < 4 ; c++ ) {

			// ... multiply column-wise with the matrix by ...

			// ... copying the column to a temporary array ...
			memcpy(s,this->state+4*c,4);

			// ... and updating the column by a matrix-vector
			// multiplication.
			this->state[4*c+0] =
					(uint8_t)(_GF.mul(0x0e,s[0]) ^ _GF.mul(0x0b,s[1]) ^
					_GF.mul(0x0d,s[2]) ^ _GF.mul(0x09,s[3]));
			this->state[4*c+1] =
					(uint8_t)(_GF.mul(0x09,s[0]) ^ _GF.mul(0x0e,s[1]) ^
					_GF.mul(0x0b,s[2]) ^ _GF.mul(0x0d,s[3]));
			this->state[4*c+2] =
					(uint8_t)(_GF.mul(0x0d,s[0]) ^ _GF.mul(0x09,s[1]) ^
					_GF.mul(0x0e,s[2]) ^ _GF.mul(0x0b,s[3]));
			this->state[4*c+3] =
					(uint8_t)(_GF.mul(0x0b,s[0]) ^ _GF.mul(0x0d,s[1]) ^
					_GF.mul(0x09,s[2]) ^ _GF.mul(0x0e,s[3]));
		}
	}
}


