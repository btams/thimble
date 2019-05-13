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
 * @file ProtectedMinutiaeRecord.cpp
 *
 * @brief
 *            Implementation of a mechanism for generating protected data
 *            from a mult-finger record containing absolutely pre-aligned
 *            minutiae as provided by the 'ProtectedMinutiaeRecord.h' header
 *            file.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <stdint.h>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include <thimble/misc/IOTools.h>
#include <thimble/security/SHA.h>
#include <thimble/security/AES.h>
#include <thimble/ecc/ReedSolomonCode.h>
#include <thimble/ecc/GuruswamiSudanDecoder.h>
#include <thimble/math/Permutation.h>
#include <thimble/math/RandomGenerator.h>
#include <thimble/math/HexagonalGrid.h>
#include <thimble/math/numbertheory/BigInteger.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/finger/ProtectedMinutiaeRecord.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Standard constructor.
	 */
	ProtectedMinutiaeRecord::ProtectedMinutiaeRecord() {
		first_init();
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @param vault
	 *            A reference to the protected minutiae record of which
	 *            a copy is created.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message will be printed to <code>stderr</code> and
	 *            the program exits with status 'EXIT_FAILURE'.
	 */
	ProtectedMinutiaeRecord::ProtectedMinutiaeRecord
	( const ProtectedMinutiaeRecord & vault ) {

		first_init();
		assign(vault);
	}

	/**
	 * @brief
	 *            Destructor.
	 */
	ProtectedMinutiaeRecord::~ProtectedMinutiaeRecord() {
		clear();
		delete this->gfPtr;
	}

	/**
	 * @brief
	 *            Assignment operator.
	 *
	 * @param vault
	 *            A reference to the protected minutiae record of which
	 *            a copy is assigned to this protected minutiae record
	 *            object.
	 *
	 * @return
	 *            A reference to this protected minutiae record object.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message will be printed to <code>stderr</code> and
	 *            the program exits with status 'EXIT_FAILURE'.
	 */
	ProtectedMinutiaeRecord & ProtectedMinutiaeRecord::operator=
	( const ProtectedMinutiaeRecord & vault ) {

		assign(vault);
		return *this;
	}

	/**
	 * @brief
	 *            Assignment operator (procedural version).
	 *
	 * @param vault
	 *            A reference to the protected minutiae record of which
	 *            a copy is assigned to this protected minutiae record
	 *            object.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message will be printed to <code>stderr</code> and
	 *            the program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::assign
	( const ProtectedMinutiaeRecord & vault ) {

		if ( this != &vault ) { // Do only execute if not a self assignment.

			// Copy primitive member variables
			this->width = vault.width;
			this->height = vault.height;
			this->dpi = vault.dpi;
			this->gridDist = vault.gridDist;
			this->s = vault.s;
			this->k = vault.k;
			this->tmax = vault.tmax;
			this->m = vault.m;
			this->slowDownFactor = vault.slowDownFactor;
			this->grid = vault.grid;

			// Copy underlying finite field
			if ( this->gfPtr == NULL && vault.gfPtr != NULL ) {
				this->gfPtr = new (nothrow) SmallBinaryField(vault.gfPtr[0]);
				if ( this->gfPtr == NULL ) {
					cerr << "ProtectedMinutiaeRecord: "
						 << "out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
			} else if ( this->gfPtr != NULL && vault.gfPtr != NULL ){
				this->gfPtr[0] = vault.gfPtr[0];
			} else if ( this->gfPtr != NULL && vault.gfPtr == NULL ) {
				delete this->gfPtr;
				this->gfPtr = NULL;
			}

			this->t = vault.t;

			// Copy vault data
			if ( this->vaultPolynomialData == NULL &&
				vault.vaultPolynomialData != NULL ) {
				int n = vaultDataSize();
				this->vaultPolynomialData =
						(uint8_t*)malloc( n * sizeof(uint8_t) );
				if ( this->vaultPolynomialData == NULL ) {
					cerr << "ProtectedMinutiaeTemplate: out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
				memcpy
				(this->vaultPolynomialData,
				 vault.vaultPolynomialData,n);
			} else if (this->vaultPolynomialData != NULL &&
					    vault.vaultPolynomialData ) {
				int n = vaultDataSize();
				this->vaultPolynomialData =
					(uint8_t*)realloc
					(this->vaultPolynomialData,n*sizeof(uint8_t));
				if ( this->vaultPolynomialData == NULL ) {
					cerr << "ProtectedMinutiaeTemplate: out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
				memcpy
				(this->vaultPolynomialData,
				 vault.vaultPolynomialData,n);
			} else if ( this->vaultPolynomialData != NULL &&
						vault.vaultPolynomialData == NULL ) {
				free(this->vaultPolynomialData);
				this->vaultPolynomialData = NULL;
			}

			// Copy encrypted data
			if ( this->encryptedVaultPolynomialData == NULL &&
				 vault.encryptedVaultPolynomialData != NULL ) {
				int n = vaultDataSize();
				this->encryptedVaultPolynomialData =
						(uint8_t*)malloc( n * sizeof(uint8_t) );
				if ( this->encryptedVaultPolynomialData == NULL ) {
					cerr << "ProtectedMinutiaeView: out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
				memcpy
				(this->encryptedVaultPolynomialData,
				 vault.encryptedVaultPolynomialData,n);
			} else if (this->encryptedVaultPolynomialData != NULL &&
					    vault.encryptedVaultPolynomialData ) {
				int n = vaultDataSize();
				this->encryptedVaultPolynomialData =
					(uint8_t*)realloc
					(this->encryptedVaultPolynomialData,n*sizeof(uint8_t));
				if ( this->encryptedVaultPolynomialData == NULL ) {
					cerr << "ProtectedMinutiaeView: out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
				memcpy
				(this->encryptedVaultPolynomialData,
				 vault.encryptedVaultPolynomialData,n);
			} else if ( this->encryptedVaultPolynomialData != NULL &&
						vault.encryptedVaultPolynomialData == NULL ) {
				free(this->encryptedVaultPolynomialData);
				this->encryptedVaultPolynomialData = NULL;
			}

			// Copy hash value of secret polynomial
			memcpy(this->hash,vault.hash,20);

			// Copy permutation by using its assignment operator
			this->permutation = vault.permutation;
		}
	}

	/**
	 * @brief
	 *            Prepares this protected minutiae record to protect
	 *            multiple absolutely pre-aligned minutiae templates
	 *            of which (un-aligned) minutiae originated from an
	 *            image of the specified width and height that has
	 *            been scanned at the specified resolution.
	 *
	 * @details
	 *            It is actually not necessary to call this method
	 *            before passing an absolutely pre-aligned minutiae
	 *            record to the \link enroll()\endlink function (note,
	 *            that this is different
	 *            for \link ProtectedMinutiaeRecord\endlink objects).
	 *            A \link MinutiaeRecord\endlink already encodes the
	 *            width, height and resolution (in pixel per centimeter)
	 *            of the original fingerprint images and
	 *            this \link inititialize()\endlink method is implicitly
	 *            called by calling the \link enroll()\endlink method.
	 *            However, if the user of this class wants to call
	 *            this object's %quantize() functions, he must call
	 *            this \link initialize()\endlink method first.
	 *
	 * @param width
	 *            The specified pixel width of the original fingerprint
	 *            images.
	 *
	 * @param height
	 *            The specified pixel height of the original fingerprint
	 *            image.
	 *
	 * @param dpi
	 *            The resolution in dots per inch of the original
	 *            fingerprint image.
	 *
	 * @warning
	 *            If the specified width, height or resolution are
	 *            negative or zero, an error message will be printed
	 *            to <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message will be printed to <code>stderr</code> and
	 *            the program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::initialize( int width , int height , int dpi ) {

		if ( width <= 0 || width <= 0 || dpi <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::initialize: bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->width == width && this->height == height && this->dpi == dpi ) {
			return;
		}

		this->width = width;
		this->height = height;
		this->dpi = dpi;

		if ( this->gridDist <= 0 ) {
			this->gridDist = (int)MathTools::round(29.0/569 * (double)dpi);
		}

		updateGrid();
		updateField();
	}

	/**
	 * @brief
	 *            Creates protected data from a record containing
	 *            multiple absolutely pre-aligned minutiae templates.
	 *
	 * @details
	 *            The function passes the input minutiae templates
	 *            through a quantization process and protects the
	 *            quantizations via the improved fuzzy vault scheme
	 *            of which data will be stored by this instance.
	 *
	 *            If not sufficient information could be extracted
	 *            from the input minutiae templates to allow successful
	 *            verification, i.e., a <em>failure to enroll</em>,
	 *            the function returns <code>false</code>;
	 *            otherwise, if successful, the function returns
	 *            <code>true</code>.
	 *
	 * @param record
	 *            A reference to a valid \link MinutiaeRecord\endlink
	 *            storing the multiple to-be-protected absolutely
	 *            pre-aligned minutiae templates.
	 *
	 * @return
	 *            <code>true</code> if the creation was successful;
	 *            otherwise, the function returns <code>false</code>.
	 *
	 * @warning
	 *            If <code>record</code> contains \link MinutiaeView
	 *            minutiae views\endlink suggesting an unknown finger
	 *            position, i.e.,
	 *            where \link MinutiaeView::getFingerPosition()
	 *            getFingerPosition()\endlink ==
	 *            \link UNKNOWN_FINGER\endlink, or if
	 *            <code>record</code> contains \link MinutiaeView
	 *            minutiae views\endlink with
	 *            duplicate \link MinutiaeView::getFingerPosition()
	 *            finger positions\endlink, an error message is printed
	 *            to <code>stderr</code> and the program exits with status
	 *            'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If \link isEnrolled()\endlink is
	 *            <code>true</code>, an error message is printed to
	 *            <code>stderr</code> and the program exits with
	 *            status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If the horizontal resolution of the record is different
	 *            from the vertical resolution, i.e.,
	 *            if \link MinutiaeRecord::getHorizontalResolution()
	 *            record.getHorizontalResolution()\endlink !=
	 *            \link MinutiaeRecord::getVerticalResolution()
	 *            record.getVerticalResolution()\endlink, then an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::enroll( const MinutiaeRecord & record ) {

		// If already enrolled, ...
		if ( isEnrolled() ) {
			// ... print an error message and ...
			cerr << "ProtectedMinutiaeRecord::enroll: "
				 << "already enrolled, clear first." << endl;
            // .. exit with status 'EXIT_FAILURE'.
            exit(EXIT_FAILURE);
		}

		int ppc = record.getHorizontalResolution();

		initialize(record.getWidth(),record.getHeight(),(int)(2.54 * ppc));

		int n = getVaultSize();

		// Temporarily generate secret polynomial and
		SmallBinaryFieldPolynomial f(getField());
		f.random(this->k,true);

		// save its SHA-1 hash value
        SHA().hash(this->hash,f.getData(),f.deg()+1);

        // Derive a record-specific permutation using the
        // SHA-1 hash value as seed
        updatePermutation();

		uint32_t *A = (uint32_t*)malloc( 10 * this->tmax * sizeof(uint32_t) );
		if ( A == NULL ) {
			cerr << "ProtectedMinutiaeRecord::enroll: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Extract the feature set from the record.
		this->t = quantize(A,record);
		// Check if enrollment can be successful.
		if ( this->t < this->k ) {
			free(A);
			return false;
		}

		// Shuffle the feature elements using the record-specific feature
		// permutation process
		for ( int j = 0 ; j < this->t ; j++ ) {
			A[j] = _reorder(A[j]);
		}

		// Supplement the feature set with blending features.
		for ( int j = this->t ; j < 10*(this->tmax) ; j++ ) {
			A[j] = n + MathTools::rand32(true) %
					(getField().getCardinality() - n );
		}
		this->t = 10* (this->tmax);

		// Compute vault polynomial.
		SmallBinaryFieldPolynomial V(getField());
		V.buildFromRoots(A,this->t);
		add(V,V,f);

		packVaultPolynomial(V);

		free(A);

		return true;
	}

	/**
	 * @brief
	 *             Performs a verification procedure to "compare"
	 *             the multiple absolutely pre-aligned minutiae templates
	 *             with the minutiae templates protected by this object.
	 *
	 * @details
	 *             The function passes the input minutiae templates
	 *             through a quantization process from which an instance
	 *             of the Reed-Solomon decoding problem is generated which
	 *             is attempted to be solved using
	 *             the \link ProtectedMinutiaeTemplate::decode()
	 *             decode()\endlink function. If successful, the
	 *             verificationis deemed to be successful and, otherwise,
	 *             failed.
	 *
	 * @param record
	 *             A minutiae record containing absolutely pre-aligned
	 *             query minutiae templates.
	 *
	 * @return
	 *             <code>true</code> if the verification procedure
	 *             was successful; otherwise, if the verification
	 *             procedure was not successful, the function returns
	 *             <code>false</code>.
	 *
	 * @warning
	 *            If <code>record</code> contains \link MinutiaeView
	 *            minutiae views\endlink suggesting an unknown finger
	 *            position, i.e.,
	 *            where \link MinutiaeView::getFingerPosition()
	 *            getFingerPosition()\endlink ==
	 *            \link UNKNOWN_FINGER\endlink, or if
	 *            <code>record</code> contains \link MinutiaeView
	 *            minutiae views\endlink with
	 *            duplicate \link MinutiaeView::getFingerPosition()
	 *            finger positions\endlink, an error message is printed
	 *            to <code>stderr</code> and the program exits with status
	 *            'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If this \link ProtectedMinutiaeRecord\endlink
     *             does not represent a successfully enrolled
     *             decrypted records, i.e.,
     *             if \link ProtectedMinutiaeRecord::isDecrypted()
     *             isDecrypted()\endlink returns <code>false</code>,
     *             the function prints an error message to
     *             <code>stderr</code> and exits with status
     *             'EXIT_FAILURE'.
     *
	 * @warning
	 *             If the horizontal resolution of the record is different
	 *             from the vertical resolution, i.e.,
	 *             if \link MinutiaeRecord::getHorizontalResolution()
	 *             record.getHorizontalResolution()\endlink !=
	 *             \link MinutiaeRecord::getVerticalResolution()
	 *             record.getVerticalResolution()\endlink, then an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::verify( const MinutiaeRecord & record ) const {

		SmallBinaryFieldPolynomial f(getField());
		return open(f,record);
	}

	/**
	 * @brief
	 *             Encrypts this instance using a key.
	 *
	 * @details
	 *             After successful enrollment
	 *             (via  \link enroll()\endlink),
	 *             the \link ProtectedMinutiaeRecord\endlink can
	 *             additionally protected via a (user-specific)
	 *             key derived from a PIN or password. Consequently,
	 *             the security level of the
	 *             \link ProtectedMinutiaeRecord\endlink can be
	 *             improved, though this requires the
	 *             \link ProtectedMinutiaeRecord\endlink to be
	 *             decrypted correctly for a successful verification,
	 *             thereby requiring that the correct can be
	 *             reproduced for successful verification.
	 *
	 * @param key
	 *             The \link AES128\endlink key used to encrypt the
	 *             \link ProtectedMinutiaeRecord\endlink.
	 *
	 * @warning
	 *             If no minutiae template has been enrolled with
	 *             this \link ProtectedMinutiaeRecord\endlink
	 *             (i.e., if \link isEnrolled()\endlink
	 *             returns <code>false</code>)
	 *             or if it contains already encrypted data (i.e.,
	 *             if \link containsEncryptedData()\endlink
	 *             returns <code>true</code>)
	 *             an error message is printed to <code>stderr</code>
	 *             and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @see decrypt()
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::encrypt( const AES128 & key ) {

		// Ensure that this instance does protect a minutiae template.
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeView::encrypt: no vault built." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure that this instance does not already hold encrypted data
		if ( isEncrypted() ) {
			cerr << "ProtectedMinutiaeView::encrypt: already encrypted."
				 << endl;
            exit(EXIT_FAILURE);
		}

		AES128 aes(key);

		int n = vaultDataSize();

		// ENCRYPTION
		aes.encrypt
		(this->vaultPolynomialData,this->vaultPolynomialData,n);

		// Exchange references to encrypted data array
		this->encryptedVaultPolynomialData = this->vaultPolynomialData;
		this->vaultPolynomialData = NULL;
	}

	/**
	 * @brief
	 *             Decrypts this instance using a key.
	 *
	 * @details
	 *             This function is to decrypt a
	 *             \link ProtectedMinutiaeRecord\endlink that has been
	 *             encrypted before using \link encrypt()\endlink.
	 *
	 *             If the \link ProtectedMinutiaeRecord\endlink has been
	 *             encrypted with the key \f$\kappa\f$, then decryption
	 *             with \f$\kappa\f$ will correctly reveal the original
	 *             data; otherwise, if decrypting with a different
	 *             key \f$\kappa'\neq\kappa\f$, this instance will hold
	 *             a candidate for the original data that is (with very
	 *             high probability) different from the original data.
	 *             To date no efficient way for distinguishing incorrect
	 *             candidate data from the correct data is known. Hence,
	 *             the difficulty to break an encrypted
	 *             \link ProtectedMinutiaeRecord\endlink equals the product
	 *             of the difficulty in guessing a  minutiae template being
	 *             of sufficient similarity to the protected minutiae template
	 *             and the difficulty in guessing the correct password.
	 *
	 * @param key
	 *             The \link AES128\endlink key used to decrypt this
	 *             \link ProtectedMinutiaeRecord\endlink.
	 *
	 * @warning
	 *             If this instance does not hold encrypted data
	 *             (i.e., if \link containsEncryptedData()\endlink
	 *             returns <code>false</code>), an error message is
	 *             printed to <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @see encrypt()
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::decrypt( const AES128 & key ) {

		// Ensure that this instance does contain encrypted data.
		if ( !containsEncryptedData() ) {
			cerr << "ProtectedMinutiaeRecord::decrypt: "
				 << "contains no encrypted vault." << endl;
            exit(EXIT_FAILURE);
		}

		AES128 aes(key);

		int n = vaultDataSize();

		// Make sure that the vault polynomial data is capable of holding
		// the encrypted data.
		if ( this->vaultPolynomialData == NULL ) {
			this->vaultPolynomialData =
					(uint8_t*)malloc( n * sizeof(uint8_t) );
			if ( this->vaultPolynomialData == NULL ) {
				cerr << "ProtectedMinutiaeRecord::decrypt: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		aes.decrypt
		(this->vaultPolynomialData,this->encryptedVaultPolynomialData,n);
	}

	/**
	 * @brief
	 *             Uses a record of absolutely pre-aligned minutiae
	 *             to open the vault.
	 *
	 * @details
	 *             This a special version of
	 *             the \link verify()\endlink function (in
	 *             fact, \link verify()\endlink
	 *             wraps around this function) in which the
	 *             secret polynomial that is used to hide the
	 *             minutiae record protected by this instance is revealed
	 *             on successful verification.
	 *
	 * @param f
	 *             Will be set to the secret polynomial on successful
	 *             verification; otherwise the content will be left
	 *             unchanged
	 *
	 * @param record
	 *             The query record of absolutely pre-aligned
	 *             minutiae template.
	 *
	 * @return
	 *             <code>true</code> if the verification procedure was
	 *             successful; otherwise, if the verification was not
	 *             successful, the function returns <code>false</code>.
	 *
	 * @warning
	 *             If this \link ProtectedMinutiaeRecord\endlink
     *             does not represent a successfully enrolled
     *             (and decrypted) record, i.e.,
     *             if \link isEnrolled()\endlink
     *             or \link isDecrypted()\endlink return
     *             <code>false</code>, the function
     *             prints an error message to <code>stderr</code> and
     *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>record</code> contains \link MinutiaeView
	 *             minutiae views\endlink suggesting an unknown finger
	 *             position, i.e.,
	 *             where \link MinutiaeView::getFingerPosition()
	 *             getFingerPosition()\endlink ==
	 *             \link UNKNOWN_FINGER\endlink, or if
	 *             <code>record</code> contains \link MinutiaeView
	 *             minutiae views\endlink with
	 *             duplicate \link MinutiaeView::getFingerPosition()
	 *             finger positions\endlink, an error message is printed
	 *             to <code>stderr</code> and the program exits with status
	 *             'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If the horizontal resolution of the record is different
	 *             from the vertical resolution, i.e.,
	 *             if \link MinutiaeRecord::getHorizontalResolution()
	 *             record.getHorizontalResolution()\endlink !=
	 *             \link MinutiaeRecord::getVerticalResolution()
	 *             record.getVerticalResolution()\endlink, then an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::open
	( SmallBinaryFieldPolynomial & f , const MinutiaeRecord & record ) const {

		// Allocate memory to temporarily hold the feature set.
		uint32_t *B = (uint32_t*)malloc( 10 * this->tmax * sizeof(uint32_t) );
		if ( B == NULL ) {
			cerr << "ProtectedMinutiaeRecord::open: out of memory." << endl;
            exit(EXIT_FAILURE);
		}

		// Extract the feature set and ...
		int t = quantize(B,record);

		// ... attempt to open with the feature set.
		bool success = open(f,B,t);

		// Free temporarily allocated memory.
		free(B);

		return success;
	}

	/**
	 * @brief
	 *             Uses a feature set as a query to open the vault.
	 *
	 * @details
	 *             This is a variant of the \link open()\endlink
	 *             function in which a set of quantized features (e.g.,
	 *             derived from a records absolutely pre-aligned minutiae) is
	 *             used to open the vault. Therein, a quantized minutia is
	 *             encoded as an element of the finite field that can be
	 *             accessed via the \link getField()\endlink function.
     *
     *             To obtain the feature set of a record with absolutely
     *             pre-aligned minutiae templates,
     *             the \link quantize()\endlink method can be used.
	 *
	 * @param f
	 *             Will be set to the secret polynomial on successful
	 *             verification; otherwise the content will be left
	 *             unchanged
	 *
	 * @param B
	 *             The feature set.
	 *
	 * @param t
	 *             Number of elements in the feature set.
	 *
	 * @return
	 *             <code>true</code> if the verification procedure was
	 *             successful; otherwise, if the verification was not
	 *             successful, the function returns <code>false</code>.
	 *
	 * @warning
	 *             If this \link ProtectedMinutiaeRecord\endlink
     *             does not represent a successfully enrolled (and
     *             decrypted) protected minutiae template, i.e., if
     *             if \link isEnrolled()\endlink
     *             or \link isDecrypted()\endlink return
     *             <code>false</code>, the function prints an error
     *             message to <code>stderr</code> and exits with status
     *             'EXIT_FAILURE'.
     *
     * @warning
     *             If <code>B</code> does not contain distinct integers
     *             of type <code>uint32_t</code> being elements of the
     *             finite field represented by \link getField()\endlink,
     *             the function runs into undocumented behavior.
     *
     * @warning
     *             If <code>B</code> does not contain at least
     *             <code>t</code> well-defined distinct feature elements,
     *             the function runs into undocumented behavior.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::open
	( SmallBinaryFieldPolynomial & f , const uint32_t *B, int t ) const {

		// Ensure that this instance does protect a feature set and ...
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::open: "
				 << "no minutiae template protected by this view." << endl;
            exit(EXIT_FAILURE);
		}

		// ... contains a decrypted polynomial.
		if ( isEncrypted() ) {
			cerr << "ProtectedMinutiaeRecord::open: "
				 << "vault is encrypted; decrypt first." << endl;
            exit(EXIT_FAILURE);
		}

		// Allocate memory to hold the set of unlocking pairs '{(x[j],y[j])}'
		uint32_t  *x , *y;
		x = (uint32_t*)malloc( t * sizeof(uint32_t) );
		y = (uint32_t*)malloc( t * sizeof(uint32_t) );
		if ( x == NULL || y == NULL ) {
			cerr << "ProtectedMinutiaeRecord::open: "
				 << "Out of memory." << endl;
            exit(EXIT_FAILURE);
		}

		bool success = false;

		// Iterate of candidates of slow-down values until
		// decoding is successful or the whole slow-down range
		// has been tested
		for ( BigInteger slowDownVal = 0 ;
			  BigInteger::compare(slowDownVal,this->slowDownFactor) < 0 ;
			  add(slowDownVal,slowDownVal,1) ) {

			SmallBinaryFieldPolynomial V = unpackVaultPolynomial(slowDownVal);

			// Build unlocking set and ...
			for ( int j = 0 ; j < t ; j++ ) {

				// ... don't forget to apply the permutation process
				x[j] = _reorder(B[j]);

				y[j] = V.eval(x[j]);
			}

			// Attempt to decode the unlocking set
			success = decode(f,x,y,t,this->k,this->hash,this->m);

			if ( success ) {
				break;
			}
		}

		free(x);
		free(y);

		return success;
	}

	/**
	 * @brief
	 *             Implementation of a Guruswami-Sudan algorithm based
	 *             decoder with which the verification process of a
	 *             this multiple minutiae template based fuzzy vault is
	 *             performed.
	 *
	 * @details
	 *             This program tries to find a polynomial <i>f</i>
	 *             with hash value equals <code>hash</code> from the
	 *             finite field pairs <i>(x[i],y[i])</i>. In the first
	 *             step, a classical Reed-Solomon decoder is used
	 *             and if it successfully output the desired polynomial,
	 *             the function returns <code>true</code>. Otherwise,
	 *             a attempt to find the desired polynomial with a
	 *             Guruswami-Sudan algorithm is performed successfully
	 *             iterating multplicities from 1 to <i>m</i> until
	 *             the desired polynomial could be decoded. If the
	 *             desired polynomial could be found, it will be
	 *             stored in <code>f</code> and the function returns
	 *             <code>true</code>; otherwise, the function
	 *             returns <code>false</code>.
	 *
	 * @see GuruswamiSudanDecoder
	 * @see ReedSolomonCode
	 *
	 * @param f
	 *             will contain the recovery of the secret polynomial on
	 *             successful decoding; otherwise, will be left unchanged.
	 *
	 * @param x
	 *             Abscissa coordinates
	 *
	 * @param y
	 *             Ordinate coordinates
	 *
	 * @param t
	 *             Number of valid finite field pairs <i>(x[j],y[j])</i>.
	 *
	 * @param k
	 *             Length of the secret polynomial whose degree is smaller
	 *             than <i>k</i>.
	 *
	 * @param hash
	 *             160-bit hash value of the correct secret polynomial
	 *             <i>f</i>, that has been originally computed via an
	 *             equivalent to
	 *             <code>
	 *              SHA().hash(hash,f.getData(),f.deg()+1)
	 *             </code>
	 *
	 * @param m
	 *             Maximal bound of multiplicity with which a
	 *             Guruswami-Sudan algorithm will is run during a
	 *             decoding attempt.
	 *
	 * @return
	 *             <code>true</code> if the decoding attempt was successful;
	 *             otherwise, the function returns <code>false</code>.
	 *
	 * @warning
	 *             If the <i>x[i]</i> contain elements not encoding a
	 *             element in the finite field represented
	 *             by \link SmallBinaryFieldPolynomial::getField()
	 *             f.getField()\endlink or of <i>x[i]</i> contain
	 *             duplicates, the decoder runs into undocumented
	 *             behavior.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::decode
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y ,
	  int t , int k , const uint8_t hash[20] , int m ) const {

		if ( k > t ) {
			return false;
		}

		SHA sha;
		uint8_t _hash[20];
		SmallBinaryFieldPolynomial _f(f.getField());

		if ( ReedSolomonCode::decode(_f,x,y,t,k) ) {
			sha.hash(_hash,_f.getData(),_f.deg()+1);
			if ( memcmp(hash,_hash,20) == 0 ) {
				f = _f;
				return true;
			}
		}

		for ( int m0 = 1 ; m0 <= m ; m0++ ) {

			GuruswamiSudanDecoder dec;

			dec.decode(x,y,t,k,m0,getField());

			for ( int j = 0 ; j < (int)dec.getDecodedList().size() ; j++ ) {

				_f = dec.getDecodedList().at(j);

				sha.hash(_hash,_f.getData(),_f.deg()+1);
				if ( memcmp(hash,_hash,20) == 0 ) {
					f = _f;
					return true;
				}
			}
		}

		return false;
	}

	/**
	 * @brief
	 *             Computes the quantization of an (absolutely pre-aligned)
	 *             minutia stemming from a finger with specified position.
	 *
	 * @param minutia
	 *             The absolutely pre-aligned minutia of which the
	 *             quantization is computed.
	 *
	 * @param L
	 *             The position of the finger from which the absolutely
	 *             pre-aligned minutia has been measured.
	 *
	 * @return
	 *             The quantization of <code>minutia</code> encoded as a
	 *             finite field element.
	 *
	 * @warning
	 *             If this object has not been initialized (either
	 *             explicitly via \link initialize()\endlink or
	 *             implicitly via \link enroll()\endlink), i.e.,
	 *             \link isInitialized()\endlink returns
	 *             <code>false</code>, an error message is printed
	 *             to <code>stderr</code> and the program exits with
	 *             status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>L</code> equals \link UNKNOWN_FINGER\endlink,
	 *             an error message is printed to <code>stderr</code> and
	 *             the program exits with status 'EXIT_FAILURE'.
	 */
	uint32_t ProtectedMinutiaeRecord::quantize
	( const thimble::Minutia & minutia , FINGER_POSITION_T L ) const {

		if ( !isInitialized() ) {
			cerr << "ProtectedMinutiaeRecord::quantize: "
				 << "not initialized." << endl;
			exit(EXIT_FAILURE);
		}

		if ( L == UNKNOWN_FINGER ) {
			cerr << "ProtectedMinutiaeRecord::quantize: finger position "
				 << "must not be unknown." << endl;
			exit(EXIT_FAILURE);
		}

		// Access relevant data of the minutia
		double x , y , theta;
		x = minutia.getX();
		y = minutia.getY();
		theta = minutia.getAngle();

		// 'i' will be used to quantize the minutia's
		// position while 'j' quantizes the minutia's
		// angle
		int i = -1 , j = -1;

		// Variable helping to save the squared distance to the closest
		// grid points yet found; its corresponding index will be
		// stored in 'i'.
		double minDist = DBL_MAX;

		// *******************************************************************
		// ********* START: Quantization of the minutia's positions **********
		// *******************************************************************

		// Iteration over the grid points
		for ( int l = 0 ; l < (int)(this->grid.size()) ; l++ ) {

			// Offset of the grid point to the position of the minutia
			double dx , dy;
			dx = this->grid[l].first  - x;
			dy = this->grid[l].second - y;

			// The squared distance of the current grid point to the minutia.
			// For reasons of efficiency we do not compute the square root.
			double dist = dx*dx+dy*dy;

			// Check whether current grid point is closer to the minutia
			// than the grid points before; if true, update;
			if ( dist < minDist ) {
				minDist = dist;
				i = l;
			}
		}

		// *******************************************************************
		// ********** END: Quantization of the minutia's positions ***********
		// *******************************************************************


		// *******************************************************************
		// *********** START: Quantization of the minutia's angle ************
		// *******************************************************************
		j = (int)floor(theta / (M_PI+M_PI) * this->s);
		// *******************************************************************
		// *********** END: Quantization of the minutia's angle ************
		// *******************************************************************

		// Combine the position's quantization with the angle's quantization
		// and return the result.
		uint32_t a = (uint32_t)i+(uint32_t)j*(uint32_t)(this->grid.size());

		return (L-1)+10*a;
	}

	/**
	 * @brief
	 *             Computes the quantization set of an
	 *             (absolutely pre-aligned) minutia template.
	 *
	 * @details
	 *             The function extracts at
	 *             most \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             quantizations from <code>view</code> (preferring
	 *             those having the highest quality estimation;
	 *             see \link MinutiaeView::getQuality()\endlink).
	 *             Consequently, the array should be able to contain
	 *             at least \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             integers of type <code>uint32_t</code>.
	 *
	 * @param A
	 *             Array in which the quantization elements will be
	 *             stored.
	 *
	 * @param view
	 *             The absolutely pre-aligned minutiae template of which
	 *             the quantization set is extracted.
	 *
	 * @return
	 *             The number of (distinct) quantization stored
	 *             in the output array <code>A</code>.
	 *
	 * @warning
	 *             If this object has not been initialized (either
	 *             explicitly via \link initialize()\endlink or
	 *             implicitly via \link enroll()\endlink), i.e.,
	 *             \link isInitialized()\endlink returns
	 *             <code>false</code>, an error message is printed
	 *             to <code>stderr</code> and the program exits with
	 *             status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If \link MinutiaeView::getFingerPosition()
	 *             view.getFingerPosition()\endlink
	 *             equals \link UNKNOWN_FINGER\endlink,
	 *             an error message is printed to <code>stderr</code> and
	 *             the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>A</code> is not able to store at
	 *             least \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             or at least \link MinutiaView::getMinutiaeCount()
	 *             view.getMinutiaeCount()\endlink integers of type
	 *             <code>uint32_t</code>, the program runs into
	 *             undocumented behavior.
	 */
	int ProtectedMinutiaeRecord::quantize
	( uint32_t *A , const MinutiaeView & view ) const {

		// Copy the minutiae template and ..
		MinutiaeView w = view;
		// ... sort their minutiae w.r.t. their estimated qualities.
		w.sortWithRespectToMinutiaeQuality();

		int t = 0; // Count of minutiae quantizations

		// Iterate over the minutiae in 'w' until the minutiae quantization
		// count reaches the bound 'tmax'
		for ( int j = 0 ; j < w.getMinutiaeCount() && t < this->tmax ; j++ ) {

			// Quantize the current minutia and ...
			uint32_t q = quantize(w.getMinutia(j),w.getFingerPosition());

			// ... determine whether the quantization is already contained
			// in the output array; ...
			bool alreadyContained = false;

			for ( int _j = 0 ; _j < t ; _j++ ) {
				if ( A[_j] == q ) {
					alreadyContained = true;
					break;
				}
			}

			// ...; if not, append it in the output array.
			if ( !alreadyContained ) {
				A[t++] = q;
			}
		}

		// Return the number of differend minutiae quantizations extraced
		// from 'view'.
		return t;
	}

	/**
	 * @brief
	 *             Computes the quantization set of a minutiae record
	 *             containing (absolutely pre-aligned) minutia template.
	 *
	 * @details
	 *             The function extracts at
	 *             most \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             (distinct) quantizations from each minutiae template
	 *             stored in <code>record</code> (preferring
	 *             those having the highest quality estimation;
	 *             see \link Minutia::getQuality()\endlink).
	 *             Consequently, the array should be able to contain
	 *             at least \link MinutiaeRecord::getViewCount()
	 *             record.getViewCount()\endlink
	 *             * \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             integers of type <code>uint32_t</code>.
	 *
	 * @param A
	 *             Array in which the quantization elements will be
	 *             stored.
	 *
	 * @param record
	 *             Record containing the absolutely pre-aligned minutiae
	 *             templates from which the quantization set is extracted.
	 *
	 * @return
	 *             The number of (distinct) quantization stored
	 *             in the output array <code>A</code>.
	 *
	 * @warning
	 *             If this object has not been initialized (either
	 *             explicitly via \link initialize()\endlink or
	 *             implicitly via \link enroll()\endlink), i.e.,
	 *             \link isInitialized()\endlink returns
	 *             <code>false</code>, an error message is printed
	 *             to <code>stderr</code> and the program exits with
	 *             status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If there is a minutiae template in <code>record</code>
	 *             specifying an unknown fingerposition, i.e.,
	 *             if \link MinutiaeView::getFingerPosition()
	 *             view.getFingerPosition()\endlink
	 *             equals \link UNKNOWN_FINGER\endlink,
	 *             an error message is printed to <code>stderr</code> and
	 *             the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>record</code> contains minutiae templates
	 *             specifying duplicate finger positions, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>A</code> is not able to store at
	 *             least \link MinutiaeRecord::getViewCount()
	 *             record.getViewCount()\endlink
	 *             * \link getMaxGenuineFeaturesPerFinger()\endlink
	 *             or at least \link MinutiaeRecord::getViewCount()
	 *             record.getViewCount()\endlink
	 *             * \link MinutiaeView::getMinutiaeCount()
	 *             view.getMinutiaeCount()\endlink integers of type
	 *             <code>uint32_t</code>, the program runs into
	 *             undocumented behavior.
	 */
	int ProtectedMinutiaeRecord::quantize
	( uint32_t *A , const MinutiaeRecord & record ) const {

		// Ensure that the horizontal and vertical resolution are equal
		if ( record.getHorizontalResolution() !=
			 record.getVerticalResolution() ) {
			cerr << "ProtectedMinutiaeRecord::quantize: horizontal resolution "
				 << "must not be different from vertical resolution." << endl;
			exit(EXIT_FAILURE);
		}

		// Ensure that there are no duplicate finger positions.
		for ( int i = 0 ; i < record.getViewCount() ; i++ ) {
			for ( int j = i+1 ; j < record.getViewCount() ; j++ ) {
				if ( record.getView(i).getFingerPosition() ==
					 record.getView(j).getFingerPosition() ) {
					cerr << "ProtectedMinutiaeRecord::quantize: "
						 << "duplicate finger positions." << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		int t = 0;

		for ( int i = 0 ; i < record.getViewCount() ; i++ ) {
			t += quantize(A+t,record.getView(i));
		}

		return t;
	}

	/**
	 * @brief
	 *             Evaluate the reordering of a quantized feature which
	 *             is associated with an (enrolled) protected minutiae
	 *             template.
	 *
	 * @details
	 *             To prevent the application of extended Euclidean algorithm
	 *             based record multiplicity attack, each protected minutiae
	 *             record is attached with a public random permutation
	 *             process acting on the quantized feature universe.
	 *             The reordering of such a quantized feature (encoded as
	 *             a finite field element) can be evaluated via this
	 *             function.
	 *
	 *             For further details on the attack that is prevented via
	 *             the implemented reordering process we refer to
	 *             <ul>
	 *              <li>
	 *               <b>J. Merkle and B. Tams (2013)</b>.
	 *               Security of the Improved Fuzzy Vault Scheme in the
	 *               Presence of Record Multiplicity.
	 *               <i>CoRR
	 *                <a href="http://arxiv.org/abs/1312.5225"
	 *                 target="_blank">
	 *                 abs/1312.5225
	 *                </a>
	 *               </i>.
	 *              </li>
	 *             </ul>
     *
	 * @param a
	 *             The feature element of which the reordering is
	 *             evaluated.
	 *
	 * @return
	 *             The reordering of <i>a</i>.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink returns
	 *             <code>false</code>, an error message may be printed to
	 *             <code>stderr</code> and the program exits with status
	 *             -1.
	 *
	 * @warning
	 *             If <i>a</i> does not encode a valid element in the
	 *             finite field that can be access via
	 *             \link getField()\endlink, an error message is printed
	 *             to <code>stderr</code> and the program exits with
	 *             status 'EXIT_FAILURE'.
	 */
	uint32_t ProtectedMinutiaeRecord::reorder( uint32_t a ) const {

		// Ensure the presence of a random permutation process
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::reorder: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		// Call reordering without check.
		return _reorder(a);
	}

	/**
	 * @brief
	 *             Changes the distance of the hexagonal grid used for
	 *             minutiae quantization.
	 *
	 * @details
	 *             Changing the distance of the hexagonal grid can cause
	 *             an update of the \link getField()
	 *             underlying finite field\endlink.
	 *
	 * @param gridDist
	 *             New distance of the hexagonal grid.
	 *
	 * @see getGridDist()
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink returns <code>true</code>,
	 *             the function prints an error message to
	 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>gridDist</code> is smaller than or equals 0,
	 *             the function prints an error message to
	 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setGridDist( int gridDist ) {

		if ( gridDist <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setGridDist: "
				 << "argument must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::setGridDist: "
				 << "parameter cannot be changed after enrollment." << endl;
			exit(EXIT_FAILURE);
		}

		this->gridDist = gridDist;
	}

	/**
	 * @brief
	 *             Change the number of possible quantizations for
	 *             fingerprint minutiae angles.
	 *
	 * @details
	 *             Changing the number of possible quantizations of the
	 *             involved minutiae angles can cause an update of
	 *             the \link getField() underlying finite field\endlink.
	 *
	 * @param s
	 *             The new number of possible quantizations of fingerprint
	 *             minutiae angles.
	 *
	 * @see getNumAngleQuanta()
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>true</code>, the
	 *             method prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <i>s</i> is smaller than or equals 0, the function
	 *             prints an error message to <code>stderr</code> and exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setNumAngleQuanta( int s ) {

		if ( s <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setNumAngleQuanta: "
				 << "argument must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::setNumAngleQuanta: "
				 << "parameter cannot be changed after enrollment." << endl;
			exit(EXIT_FAILURE);
		}

		this->s = s;

		if ( isInitialized() ) {
			updateField();
		}
	}

	/**
	 * @brief
	 *             Change the size of the secret polynomial being used
	 *             to obfuscate the minutiae template's quantization in
	 *             order to protect it.
	 *
	 * @param k
	 *             The new size of the secret polynomial being used to
	 *             obfuscate the minutiae template's quantization.
	 *
	 * @see getSecretSize()
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>true</code>, the
	 *             method prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <i>k</i> is smaller than or equals 0, the method
	 *             prints an error message to <code>stderr</code> and
	 *             exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setSecretSize( int k ) {

		if ( k <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setSecretSize: "
				 << "argument must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::setSecretSize: "
				 << "parameter cannot be changed after enrollment." << endl;
			exit(EXIT_FAILURE);
		}

		this->k = k;
	}

	/**
	 * @brief
	 *             Changes the maximal number of minutiae quantizations
	 *             extracted each of the to-be-protected record's
	 *             absolutely pre-aligned minutiae templates during
	 *             enrollment.
	 *
	 * @param tmax
	 *             The new maximal number of minutiae quantizations
	 *             extracted from each record's absolutely pre-aligned
	 *             minutiae template.
	 *
	 * @see getMaxGenuineFeatures()
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>true</code>,
	 *             the method prints an error message to
	 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>tmax</code> is smaller than or equals 0, the
	 *             method prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setMaxGenuineFeaturesPerFinger( int tmax ) {

		if ( tmax <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setMaxGenuineFeaturesPerFinger: "
				 << "argument must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::setMaxGenuineFeaturesPerFinger: "
				 << "parameter cannot be changed after enrollment." << endl;
			exit(EXIT_FAILURE);
		}

		this->tmax = tmax;
	}

	/**
	 * @brief
	 *            Sets the maximal multiplicity used during a decoding
	 *            attempt with a Guruswami-Sudan algorithm.
	 *
	 * @details
	 *            The specified multiplicity parameter will be used
	 *            by the \link verify()\endlink and \link open()\endlink
	 *            functions that pass it as input to
	 *            the \link decode()\endlink function.
	 *
	 * @see getMultiplicity()
	 *
	 * @param m
	 *            The maximal multiplicity used in Guruswami-Sudan
	 *            algorithm based decoding attempt.
	 *
	 *
	 * @warning
	 *            If <code>m</code> is smaller than 1, an error message
	 *            will be printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setMultiplicity( int m ) {

		if ( m <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setMultiplicity: "
				 << "argument must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		this->m = m;
	}

	/**
	 * @brief
	 *            Specifies the slow-down factor to artificially
	 *            increase the time for an impostor verification
	 *            attempt.
	 *
	 * @details
	 *            The time needed for an impostor to perform
	 *            off-line attacks is multiplied by the slow-down
	 *            factor, thereby increasing security while, on the
	 *            other hand, also increasing response time of genuine
	 *            verification attempts---a trade-off thus but with
	 *            the capability of maintaining a relative security
	 *            measure with the advent of more and more powerful
	 *            computers.
	 *
	 * @param slowDownFactor
	 *            The specified slow-down factor for this protected
	 *            minutiae record.
	 *
	 * @see getSlowDownFactor()
	 *
	 * @warning
	 *            If this protected minutiae record object already
	 *            contains protected data, i.e.,
	 *            if \link isEnrolled()\endlink returns <code>true</code>,
	 *            an error message will be printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>slowDownFactor</code> encodes an integer
	 *            smaller than 1, an error message will be printed
	 *            to <code>stderr</code> and the program exits
	 *            with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be allocated, then an
	 *            error message will be printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::setSlowDownFactor
	( const BigInteger & slowDownFactor ) {

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::setSlowDownFactor: "
				 << "parameter cannot be changed after enrollment." << endl;
			exit(EXIT_FAILURE);
		}

		if ( slowDownFactor.sign() <= 0 ) {
			cerr << "ProtectedMinutiaeRecord::setSlowDownFactor: "
				 << "parameter must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		this->slowDownFactor = slowDownFactor;
	}

	/**
	 * @brief
	 *             Access the distance of the hexagonal grid used
	 *             for minutiae quantization.
	 *
	 * @return
	 *             The distance of the hexagonal grid used for minutiae
	 *             quantization.
	 *
	 * @see setGridDist(int)
	 */
	int ProtectedMinutiaeRecord::getGridDist() const {

		return this->gridDist;
	}

	/**
	 * @brief
	 *             Access the number of possible quantizations for
	 *             fingerprint minutiae angles.
	 *
	 * @return
	 *             The number of possible quantizations of fingerprint
	 *             minutiae angles.
	 *
	 * @see setNumAngleQuanta(int)
	 */
	int ProtectedMinutiaeRecord::getNumAngleQuanta() const {

		return this->s;
	}

	/**
	 * @brief
	 *             Access the size of the secret polynomial being used
	 *             to obfuscate the minutiae template's quantization in
	 *             order to protect it.
	 *
	 * @return
	 *             The size of the secret polynomial being used to
	 *             obfuscate the minutiae template's quantization.
	 *
	 * @see setSecretSize(int)
	 */
	int ProtectedMinutiaeRecord::getSecretSize() const {

		return this->k;
	}

	/**
	 * @brief
	 *             Access the maximal number of minutiae quantizations
	 *             extracted from the to-be-protected record's
	 *             minutiae templates on enrollment.
	 *
	 * @return
	 *             The maximal number of minutiae quantizations extracted
	 *             from the to-be-protected record's minutiae templates on
	 *             enrollment.
	 *
	 * @see setMaxGenuineFeaturesPerFinger(int)
	 */
	int ProtectedMinutiaeRecord::getMaxGenuineFeaturesPerFinger() const {

		return this->tmax;
	}

	/**
	 * @brief
	 *             Accesses the maximal multiplicity bound used by this
	 *             protected minutiae record object's during decoding
	 *             with a Guruswami-Sudan algorithm.
	 *
	 * @return
	 *             The maximal multiplicity bound used by this
	 *             protected minutiae record object's during decoding
	 *             with a Guruswami-Sudan algorithm.
	 *
	 * @see setMultiplicity(int)
	 */
	int ProtectedMinutiaeRecord::getMultiplicity() const {

		return this->m;
	}

	/**
	 * @brief
	 *             Access the resolution of the fingerprint images from
	 *             which the minutiae contained in the protected minutiae
	 *             record's minutiae template have been measured.
	 *
	 * @return
	 *             The resolution of the fingerprint images from which the
	 *             minutiae contained in the protected minutiae  record's
	 *             minutiae template have been measured.
	 */
	int ProtectedMinutiaeRecord::getResolution() const {

		return this->dpi;
	}

	/**
	 * @brief
	 *            Returns the slow-down factor used to artifically
	 *            increase the time for a verification attempt on
	 *            impostor verification.
	 *
	 * @details
	 *            Unless manually specified via
	 *            the \link setSlowDownFactor()\endlink, the result
	 *            of this function is 1.
	 *
	 *            The time needed for an impostor to perform
	 *            off-line attacks is multiplied by the slow-down
	 *            factor, thereby increasing security while also
	 *            increasing response time of genuine verification
	 *            attempts---a trade-off thus but with the capability
	 *            of maintaining a relative security measure with
	 *            the advent of more and more powerful computers,
	 *
	 * @return
	 *            The slow-down factor of this protected minutiae
	 *            template.
	 *
	 * @see setSlowDownFactor()
	 */
	const BigInteger & ProtectedMinutiaeRecord::getSlowDownFactor() const {

		return this->slowDownFactor;
	}

	/**
	 * @brief
	 *             Access the vector of hexagonal grid points centered in
	 *             a region in which absolutely pre-aligned minutiae
	 *             can occur.
	 *
	 * @return
	 *             The vector of hexagonal grid points covering the region
	 *             in which absolutely pre-aligned minutiae can occur.
	 */
	const vector< pair<double,double> > & ProtectedMinutiaeRecord::getGrid() const {

		return this->grid;
	}

	/**
	 * @brief
	 *             Returns the number of hexagonal grid points centered
	 *             in a region in which absolutel pre-algined minutiae can
	 *             occur.
	 *
	 * @return
	 *             (int)\link getGrid()\endlink .size()
	 */
	int ProtectedMinutiaeRecord::getGridSize() const {

		return (int)(getGrid().size());
	}

	/**
	 * @brief
	 *             Returns the size of the vault.
	 *
	 * @details
	 *             The size of the vault is the number of possible
	 *             minutiae quantizations.
	 *
	 * @return
	 *             \link getGridSize()\endlink
	 *             * \link getNumAngleQuanta()\endlink
	 */
	int ProtectedMinutiaeRecord::getVaultSize() const {

		return 10 * getGridSize() * getNumAngleQuanta();
	}

	/**
	 * @brief
	 *             Returns the number in bytes that is needed to represent
	 *             the data of the vault polynomial stored by this
	 *             protected minutiae template.
	 *
	 * @details
	 *             This function returns the number of bytes being stored
	 *             in the array returned by \link getVaultData()\endlink
	 *             and is a multiple of 16.
	 *
	 * @return
	 *             The number in bytes of this object's vault data.
	 *
	 * @see getVaultData()
	 */
	int ProtectedMinutiaeRecord::getVaultDataSize() const {

		if ( !isEnrolled() ) {
			return 0;
		}

		return vaultDataSize();
	}

	/**
	 * @brief
	 *             Returns an array in which the data for the vault
	 *             polynomial of this object is stored.
	 *
	 * @details
	 *             Let
	 *             \f[
	 *              V(X)=X^t+\sum_{j=0}^{t-1} V_j\cdot X^j
	 *             \f]
	 *             be the vault polynomial that protects the feature set
	 *             \f$A\f$ (derived from the minutiae template protected
	 *             by this object). Write \f$b_{0},...,b_{m-1}\f$
	 *             be the bit sequence resulting from a concatenation of
	 *             the coefficients \f$V_0,...,V_{t-1}\f$. Let
	 *             \f$b_{m},...,b_{M-1}\f$ be random bits such that
	 *             \f$M=96\cdot\lceil m/96\rceil\f$, i.e.,
	 *             \f$ b_{0},...,b_{m-1}\f$ are padded with random bits
	 *             such that the total number of bits becomes a multiple
	 *             of 96. An \link AES128 AES key\endlink is then derived
	 *             from a random positive integer smaller
	 *             than \link getSlowDownFactor()\endlink and used to
	 *             "encrypt" the bits \f$b_0,...,b_{M-1}\f$. It is important
	 *             to note that, if the slow-down factor is equals 1, then
	 *             there is only one possibility to generate the AES key,
	 *             i.e., derived from the integer 0. In such we obtain the
	 *             "encrypted" bit sequence
	 *             \f[
	 *              (b'_0,...,b'_{M-1})
	 *             \f]
	 *             which is stored byte-wisely in the array returned by
	 *             this function.
	 *
	 *             The reason why the polynomial data needs to be padded
	 *             such that it can be passed through AES functions is the
	 *             following. If \link getSlowDownFactor()\endlink is
	 *             larger than 1, then the \link AES128 AES key\endlink
	 *             through which the padded vault polynomial bits are
	 *             passed serve as a slow-down mechanism: On a false
	 *             verification attempt, the impostor has to iterate
	 *             through all \link AES128 AES keys\endlink derived from
	 *             integers smaller than \link getSlowDownFactor()\endlink
	 *             before he (most likely) gets rejected thereby improving
	 *             security against false-accept attack while, on the other
	 *             hand, also increasing verification time on a genuine
	 *             verification attempt---a trade-off thus. However, as
	 *             computers become more and more efficient, the slow-down
	 *             factor may be increased and thus a relative security
	 *             measure can be maintained.
	 *
	 *             Also note that the vault data can be further
	 *             encrypted using a user-specific password, for which
	 *             purpose the method \link encrypt()\endlink
	 *             and \link decrypt()\endlink can be used for encryption
	 *             and decryption, respectively.
	 *
	 * @return
	 *             The array containing the vault polynomial's (possibly
	 *             encrypted) data stored by this object; the result is
	 *             <code>NULL</code> in case \link isEnrolled()\endlink
	 *             returns <code>false/<code>.
	 *
	 * @see getVaultDataSize()
	 * @see setSlowDownFactor()
	 * @see getSlowDownFactor()
	 * @see encrypt()
	 * @see decrypt()
	 */
	const uint8_t * ProtectedMinutiaeRecord::getVaultData() const {

		if ( this->vaultPolynomialData != NULL ) {
			return this->vaultPolynomialData;
		} else if ( this->encryptedVaultPolynomialData != NULL ) {
			return this->encryptedVaultPolynomialData;
		}

		return NULL;
	}

	/**
	 * @brief
	 *            Derives the vault polynomial from the data stored by
	 *            this object.
	 *
	 * @details
	 *            If \link getSlowDownFactor()\endlink equals 1 (which, by
	 *            default, is the case), then this function can be called
	 *            without specifying an argument thereby returning the
	 *            correct vault polynomial.
	 *
	 *            Otherwise, the correct vault polynomial has been
	 *            encrypted using a random secret key derived from an
	 *            integer, called the <i>slow-down value</i> smaller
	 *            than \link getSlowDownFactor()\endlink and then the
	 *            correct vault polynomial will be returned if the correct
	 *            slow-down value has been passed to this function.
	 *
	 *            This function is used by the \link open()\endlink
	 *            and implicitly by the \link verify()\endlink
	 *            functions in which slow-down values starting from
	 *            zero is successively increased and in each iteration
	 *            a decoding attempt is performed until the vault
	 *            has been successfully be opened or all slow-down
	 *            values smaller than \link getSlowDownFactor()\endlink
	 *            have been tested.
	 *
	 * @param slowDownVal
	 *            A guess for the slow-down value.
	 *
	 * @return
	 *            A candidate for the correct vault polynomial (which is
	 *            correct if the guess for the slow-down value is correct).
	 *
	 * @warning
	 *            If this object does not contain (a candidate) for decrypted
	 *            vault polynomial data, i.e.,
	 *            if \link isDecrypted()\endlink returns <code>false</code>,
	 *            then an error message will be printed to
	 *            <code>stderr</code> and the program exits with status
	 *            'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>slowDownVal</code> is negative or greater than
	 *            or equals \link getSlowDownFactor()\endlink, then an
	 *            error message will be printed to <code>stderr</code> and
	 *            the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be allocated, an error message
	 *            will be printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 */
	SmallBinaryFieldPolynomial ProtectedMinutiaeRecord::unpackVaultPolynomial
	( const BigInteger & slowDownValue ) const {


		if ( !isDecrypted() ) {
			cerr << "ProtectedMinutiaeRecord::unpackVaultPolynomial: "
				 << "no decrypted vault data." << endl;
			exit(EXIT_FAILURE);
		}

		if ( slowDownValue.sign() < 0 ||
			 BigInteger::compare(slowDownValue,getSlowDownFactor()) >= 0 ) {
			cerr << "ProtectedMinutiaeRecord::unpackVaultPolynomial: "
				 << "slow-down value must be positive and smaller than the "
				 << "slow-down factor" << endl;
			exit(EXIT_FAILURE);
		}

		AES128 aes = deriveKey(slowDownValue);

		int t , d , n;
		t = this->t;
		d = this->gfPtr->getDegree();
		n = vaultDataSize();

		uint8_t *data = (uint8_t*)malloc( n * sizeof(uint8_t) );
		uint32_t *coeffs = (uint32_t*)malloc( t * sizeof(uint32_t ) );
		if ( data == NULL || coeffs == NULL ) {
			cerr << "ProtectedMinutiaeRecord::createVaultPolynomialCandidate: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		aes.decrypt(data,this->vaultPolynomialData,n);

		// Unpack decrypted data into coefficient vector
		split_into_bit_vectors(coeffs,data,t,d);

		SmallBinaryFieldPolynomial V(getField());
		V.setCoeff(t,1);
		for ( int j = 0 ; j < t ; j++ ) {
			V.setCoeff(j,coeffs[j]);
		}

		free(data);
		free(coeffs);

		return V;
	}

	/**
	 * @brief
	 *             Access the SHA-1 hash value of the secret key
	 *             generated on enrollment used to obfuscate the
	 *             correct minutiae record's quantization.
	 *
	 * @return
	 *             The SHA-1 hash value of the secret key used to
	 *             obfuscate the correct minutiae template's quantization.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink returns
	 *             <code>false</code>, an error message is printed to
	 *             <code>stderr</code> and the program exits with status
	 *             -1.
	 */
	const uint8_t *ProtectedMinutiaeRecord::getHash() const {

		// Hash value of secret polynomial is extracted on enrollment
		if  ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::getHash(): "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		return this->hash;
	}

	/**
	 * @brief
	 *             Access the finite field in which computations with
	 *             the improved fuzzy vault scheme are performed to
	 *             protect multiple minutiae templates.
	 *
	 * @returns
	 *             A reference to the \link SmallBinaryField
	 *             finite field\endlink in which computations with the
	 *             improved fuzzy vault scheme are performed with this
	 *             object.
	 */
	const SmallBinaryField & ProtectedMinutiaeRecord::getField() const {

		if ( this->gfPtr == NULL ) {
			cerr << "ProtectedMinutiaeRecord::getField: "
				 << "not initialized." << endl;
			exit(EXIT_FAILURE);
		}

		return this->gfPtr[0];
	}

	/**
	 * @brief
	 *            Returns a boolean value that indicates whether this
	 *            object has been initialized such that it can protect an
	 *            absolutely pre-aligned minutiae template.
	 *
	 * @details
	 *            After an object of this class has been created via the
	 *            standard
	 *            constructor \link ProtectedMinutiaeRecord()\endlink
	 *            the width, height, and resolution to which a
	 *            to-be-protected minutiae record corresponds to is
	 *            unknown but need to be known on enrollment and can
	 *            be specified via the \link initialize()\endlink
	 *            method; the \link initialize()\endlink method is
	 *            called implicitly during a call of
	 *            the \link enroll()\endlink; however, if we want to
	 *            make use of the \link quantize()\endlink function
	 *            without generating actual protected minutiae record
	 *            data, we have to call \link initialize()\endlink
	 *            explicitly. This function returns whether a call of
	 *            the \link initialize()\endlink method has been made and
	 *            is used by the \link quantize()\endlink function to
	 *            produce error handling if the fmethod was not
	 *            called.
	 *
	 * @return
	 *            <code>true</code> if this object has been initialized
	 *            such that \link quantize()\endlink can be called safely;
	 *            otherwise, this function returns <code>false</code>.
	 */
	bool ProtectedMinutiaeRecord::isInitialized() const {

		return this->width > 0 && this->height > 0 &&
			   this->gridDist > 0 && this->s > 0;
	}

	/**
	 * @brief
	 *            Access whether this object does hold
	 *            actual protected data of a minutiae record.
	 *
	 * @return
	 *            <code>true</code> if this object  does protect
	 *            a minutiae record; otherwise <code>false</code>.
	 *
	 * @see enroll()
	 */
	bool ProtectedMinutiaeRecord::isEnrolled() const {

		return this->vaultPolynomialData != NULL ||
			   this->encryptedVaultPolynomialData != NULL;
	}

	/**
	 * @brief
	 *            Determine whether this instance does hold a
	 *            protected minutiae template that is additionally
	 *            protected by a (user-specific) key.
	 *
	 * @details
	 *            If this function returns <code>true</code>, the
	 *            \link ProtectedMinutiaeRecord\endlink has been
	 *            encrypted via \link encrypt()\endlink;
	 *            otherwise, the function returns <code>false</code>.
	 *            If this \link ProtectedMinutiaeRecord\endlink
	 *            has been decrypted via \link decrypt()\endlink
	 *            such that it represents a candidate
	 *            for a decrypted \link ProtectedMinutiaeRecord\endlink,
	 *            the function returns <code>false</code>.
	 *
	 * @return
	 *            <code>true</code> if this
	 *            \link ProtectedMinutiaeRecord\endlink does protect
	 *            a minutiae template being additionally protected
	 *            by an \link AES128\endlink key and does not represent
	 *            a candidate for a decrypted
	 *            \link ProtectedMinutiaeRecord\endlink; otherwise, the
	 *            function returns <code>false</code>.
	 *
	 * @see encrypt()
	 */
	bool ProtectedMinutiaeRecord::isEncrypted() const {

		return this->encryptedVaultPolynomialData != NULL &&
				this->vaultPolynomialData == NULL;
	}

	/**
	 * @brief
	 *            Determine whether this instance holds a decrypted
	 *            candidate for a protected minutiae template.
	 *
	 * @details
	 *            Any encrypted \link ProtectedMinutiaeRecord\endlink
	 *            can be decrypted via an \link AES128\endlink key and,
	 *            subsequently, contains a candidate for the decrypted
	 *            \link ProtectedMinutiaeRecord\endlink.
	 *
	 *            Note that, if this
	 *            \link ProtectedMinutiaeRecord\endlink
	 *            has not been decrypted via \link decrypt()\endlink
	 *            but contains enrolled data, the function returns
	 *            <code>true</code>.
	 *
	 * @return
	 *            <code>true</code> if this
	 *            \link ProtectedMinutiaeRecord\endlink contains a
	 *            candidate for the unencrypted
	 *            \link ProtectedMinutiaeRecord\endlink; otherwise,
	 *            the result will be <code>false</code>.
	 *
	 * @see decrypt()
	 */
	bool ProtectedMinutiaeRecord::isDecrypted() const {

		return this->vaultPolynomialData != NULL;
	}

	/**
	 * @brief
	 *            Determine whether this instance does hold
	 *            encrypted data of a protected minutiae template.
	 *
	 * @details
	 *            Note that a \link ProtectedMinutiaeRecord\endlink
	 *            keeps the encrypted data even after being decrypted
	 *            via \link decrypt()\endlink
	 *
	 * @return
	 *            <code>true</code> if the
	 *            \link ProtectedMinutiaeRecord\endlink contains
	 *            encrypted data; otherwise, the function returns
	 *            <code>false</code>.
	 *
	 * @see encrypt()
	 */
	bool ProtectedMinutiaeRecord::containsEncryptedData() const {

		return this->encryptedVaultPolynomialData != NULL;
	}

	/**
	 * @brief
	 *             Clears all data that is related with a minutiae
	 *             record protected by this instance (if any).
	 *
	 * @details
	 *             Any \link ProtectedMinutiaeRecord\endlink is empty
	 *             (cleared) by default until a
	 *             \link MinutiaeRecord minutiae record\endlink is passed
	 *             through the \link enroll()\endlink function returning
	 *             <code>true</code>. If the user wants to reuse this
	 *             \link ProtectedMinutiaeRecord\endlink objectto protect
	 *             another \link MinutiaeView minutiae template\endlink
	 *             via \link enroll()\endlink, he must clear it first via
	 *             this method.
	 */
	void ProtectedMinutiaeRecord::clear() {

		// Free data
		free(this->vaultPolynomialData);
		free(this->encryptedVaultPolynomialData);

		// Set NULL points
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;

		this->t = -1;

		memset(this->hash,0,20);
	}

	/**
	 * @brief
	 *             Determines the size in bytes needed to store this
	 *             protected minutiae record.
	 *
	 * @return
	 *             The size in bytes needed to store the data for this
	 *             protected minutiae template.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>false</code>, this
	 *             function prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeRecord::getSizeInBytes() const {

		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::getSizeInBytes: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		int size = 0;

		size += 10; // size for the header
		size += 2; // size to encode 'width'
		size += 2; // size to encode 'height'
		size += 2; // size to encode 'dpi'
		size += 1; // size to encode 'gridDist'
		size += 1; // size to encode 's'
		size += 1; // size to encode 'k'
		size += 1; // size to encode 'tmax'
		size += 4; // size to encode 'm'
		size += 4 + this->slowDownFactor.getSizeInBytes(); // size to encode the slow-down factor
		size += 4; // size to encode generator polynomial of finite field
		size += 2; // size to encode 't'
		size += 1; // size to encode whether the vault is encrypted or not
		size += vaultDataSize(); // size to encode the vault data
		size += 20; // size to encode hash value of the secret polynomial.

		return size;
	}

	/**
	 * @brief
	 *             Stores data defining this protected minutiae template to
	 *             the specified byte array.
	 *
	 * @param data
	 *             The byte array to which the data for this protected
	 *             minutiae record is written.
	 *
	 * @return
	 *             The function returns \link getSizeInBytes()\endlink.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>false</code>, this
	 *             function prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>data</code> cannot hold at least
	 *             \link getSizeInBytes()\endlink bytes, the function runs
	 *             into undocumented behavior.
	 */
	int ProtectedMinutiaeRecord::toBytes( uint8_t *data ) const {

		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeRecord::toBytes: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		// Temporary variable
		uint32_t tmp;

		// Count that keeps track of the written bytes and that
		// will be equals getSizeInBytes() after all data have been
		// written and is finally the returned value of this
		// function
		int offset = 0;



		// Write header
		memcpy(data,"PMR140822",10);
		offset = 10;

		// Write 'width'
		tmp = (uint32_t)(this->width);
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 2;

		// Write 'height'
		tmp = (uint32_t)(this->height);
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 2;

		// Write 'dpi'
		tmp = (uint32_t)(this->dpi);
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 2;

		// Write 'gridDist'
		data[offset+0] = (uint8_t)(this->gridDist);
		offset += 1;

		// Write 's'
		data[offset+0] = (uint8_t)(this->s);
		offset += 1;

		// Write 'k'
		data[offset+0] = (uint8_t)(this->k);
		offset += 1;

		// Write 'tmax'
		data[offset+0] = (uint8_t)(this->tmax);
		offset += 1;

		// Write 'm'
		tmp = (uint32_t)(this->m);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write data size for 'slowDownFactor'
		tmp = (uint32_t)(this->slowDownFactor.getSizeInBytes());
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;
		// Write data of 'slowDownFactor'
		this->slowDownFactor.toBytes(data+offset);
		offset += this->slowDownFactor.getSizeInBytes();

		// Write finite field generator polynomial
		tmp = (uint32_t)(this->gfPtr->getDefiningPolynomial().rep);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write 't'
		tmp = (uint32_t)(this->t);
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 2;

		// Write encrypted/unencrypted vault data
		int n = vaultDataSize();
		if ( containsEncryptedData() ) {
			data[offset] = 1;
			offset++;
			memcpy(data+offset,this->encryptedVaultPolynomialData,n);
		} else {
			data[offset] = 0;
			offset++;
			memcpy(data+offset,this->vaultPolynomialData,n);
		}
		offset += n;

		// Write hash of secret polynomial
		memcpy(data+offset,this->hash,20);
		offset += 20;

		return offset;
	}

	/**
	 * @brief
	 *             Initialize this protected minutiae record by the
	 *             specified bytes.
	 *
	 * @details
	 *             If this protected minutiae records has been successfully
	 *             initialized by the specified bytes, the function returns
	 *             the number of relevant bytes; otherwise, if an error
	 *             occurred, this protected minutiae template will be left
	 *             unchanged and the result of the function will be -1.
	 *
	 * @param data
	 *             Contains <code>size</code> well-defined bytes of type
	 *             <code>uint8_t</code>.
	 *
	 * @param size
	 *             The number of well-defined bytes contained in
	 *             <code>data</code>.
	 *
	 * @return
	 *             The number of bytes (<=size) that were read from
	 *             <code>data</code> to initialize this protected minutiae
	 *             template; otherwise, if the first bytes of
	 *             <code>data</code> do not specify a well-defined protected
	 *             minutiae template, the function returns -1.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeRecord::fromBytes( const uint8_t *data , int size ) {

		// Keeps track of the read data and is finally the returned result
		// if all data has been successfully read
		int offset = 0;

		if ( size < offset+10 ) {
			return -1;
		}

		// Temporarily stores the read data; if successfully read, then
		// this object will be initialized using the static 'swap'
		// method
		ProtectedMinutiaeRecord tmp;

		{ // Check if header is valid
			if ( data[0] != 'P' || data[1] != 'M' || data[2] != 'R' ||
				 data[9] != '\0') {
				return -1;
			}

			for ( int j = 3 ; j < 9 ; j++ ) {
				if ( !isdigit(data[j]) ) {
					return -1;
				}
			}
		}

		offset += 10;

		// Read and check 'width'
		if ( size < offset+2 ) { return -1; }
		tmp.width = (int)data[offset]; tmp.width <<= 8; tmp.width += (int)data[offset+1];
		if ( tmp.width == 0 ) { return -1; }
		offset += 2;

		// Read 'height'
		if ( size < offset+2 ) { return -1; }
		tmp.height = (int)data[offset]; tmp.height <<= 8; tmp.height += (int)data[offset+1];
		if ( tmp.height == 0 ) { return -1; }
		offset += 2;

		// Read and check 'dpi'
		if ( size < offset+2 ) { return -1; }
		tmp.dpi = (int)data[offset]; tmp.dpi <<= 8; tmp.dpi += data[offset+1];
		if ( tmp.dpi < 300 || tmp.dpi > 1000 ) { return -1; }
		offset += 2;

		// Read and check 'gridDist'
		if ( size < offset+1 ) { return -1; }
		tmp.gridDist = (int)data[offset];
		if ( tmp.gridDist == 0 ) { return -1; }
		offset += 1;

		// Read and check 's'
		if ( size < offset+1 ) { return -1; }
		tmp.s = (int)data[offset];
		if ( tmp.s == 0 ) { return -1; }
		offset += 1;

		// Read and check 'k'
		if ( size < offset+1 ) { return -1; }
		tmp.k = (int)data[offset];
		if ( tmp.k == 0 ) { return -1; }
		offset += 1;

		// Read 'tmax' and check
		if ( size < offset+1 ) { return -1; }
		tmp.tmax = (int)data[offset];
		if ( tmp.tmax == 0 || tmp.k > tmp.tmax ) { return -1; }
		offset += 1;

		// Read 'm'
		if ( size < offset+4 ) { return -1; }
		tmp.m = (int)data[offset+0]; tmp.m <<= 8;
		tmp.m += (int)data[offset+1]; tmp.m <<= 8;
		tmp.m += (int)data[offset+2]; tmp.m <<= 8;
		tmp.m += (int)data[offset+3];
		offset += 4;

		// Read slow-down factor
		{
			if ( size < offset+4 ) { return -1; }
			int sizeOfSlowDownFactor;
			sizeOfSlowDownFactor  = (int)data[offset+0]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+1]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+2]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+3];
			offset += 4;
			if ( size < offset+sizeOfSlowDownFactor ) { return -1; }
			tmp.slowDownFactor.fromBytes(data+offset,sizeOfSlowDownFactor);
			offset += sizeOfSlowDownFactor;
		}

		// Build the grid
		tmp.updateGrid();

		// Read the finite field
		{
			if ( size < offset+4 ) { return -1; }

			uint32_t rep;
			rep  = (uint32_t)data[offset+0]; rep <<= 8;
			rep += (uint32_t)data[offset+1]; rep <<= 8;
			rep += (uint32_t)data[offset+2]; rep <<= 8;
			rep += (uint32_t)data[offset+3];
			offset += 4;


			if ( this->gfPtr != NULL && this->gfPtr->getDefiningPolynomial().rep == rep ) {
				tmp.gfPtr = new(nothrow) SmallBinaryField(this->gfPtr[0]);
			} else {
				tmp.gfPtr = new(nothrow) SmallBinaryField(SmallBinaryPolynomial(rep));
			}
			if ( tmp.gfPtr == NULL ) {
				cerr << "ProtectedMinutiaeRecord::fromBytes: "
					 << "out of memory." << endl;
                exit(EXIT_FAILURE);
			}

			if ( (int)(tmp.gfPtr->getCardinality()) <
				 10 * this->s * (int)tmp.grid.size() ) {
				return -1;
			}
		}

		// Read 't'
		if ( size < offset+2 ) { return -1; }
		tmp.t = (int)data[offset]; tmp.t <<= 8; tmp.t += (int)data[offset+1];
		offset += 2;

		// *******************************************************************
		// ******** BEGIN: Read vault data ***********************************
		// *******************************************************************

		// Test whether there is enough data provided.
		int n = tmp.vaultDataSize();
		if ( size < offset+n+1 ) {
			return -1;
		}

		// Allocate memory for reading encrypted/unencrypted vault data
		uint8_t *vaultData = (uint8_t*)malloc( n * sizeof(uint8_t) );
		if ( vaultData == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::fromBytes: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Read data
		memcpy(vaultData,data+offset+1,n);
		if ( data[offset] != 0 ) {
			tmp.encryptedVaultPolynomialData = vaultData;
		} else {
			tmp.vaultPolynomialData = vaultData;
		}
		offset += n+1;

		// *******************************************************************
		// ********** End: Read vault data ***********************************
		// *******************************************************************

		// Read 'hash'
		if ( size < offset+20 ) { return -1; }
		memcpy(tmp.hash,data+offset,20);
		offset += 20;

		// Update the permutation which is selected pseudo-randomly using
		// the hash value as seed.
		tmp.updatePermutation();

		// The data has been successfully read. Thus, swap this objects'
		// content with the read record and ...
		swap(*this,tmp);

		// ... return the number of read bytes.
		return offset;
	}

	/**
	 * @brief
	 *             Write data defining this protected minutiae record to
	 *             the specified file pointer.
	 *
	 * @param out
	 *             Specifies the file to which this protected minutiae
	 *             template is written.
	 *
	 * @return
	 *             The function returns the number of bytes that have been
	 *             successfully written to <code>out</code>.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>false</code>, this
	 *             function prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeRecord::write( FILE *out ) const {

		uint8_t *data;
		int size , wsize;

		// Initialize byte array ...
		size = getSizeInBytes();
		data = (uint8_t*)malloc( size * sizeof(uint8_t) );
		if ( data == NULL ) {
			cerr << "ProtectedMinutiaeRecord::write: "
				 << "out of memory." << endl;
            exit(EXIT_FAILURE);
		}

		// ... in which this protected minutiae template is packed.
		toBytes(data);

		// Write the byte array to 'out'.
		wsize = fwrite(data,sizeof(uint8_t),size,out);

		free(data);

		return wsize;
	}

	/**
	 * @brief
	 *             Writes data defining this protected minutiae record to
	 *             the specified file.
	 *
	 * @param file
	 *             Path to the file to which this protected minutiae
	 *             record is written.
	 *
	 * @return
	 *             The function returns <code>true</code> if this protected
	 *             minutiae record has been successfully written
	 *             to the specified file; otherwise, if an error occurred,
	 *             the function returns <code>false</code>.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink is <code>false</code>, this
	 *             function prints an error message to <code>stderr</code>
	 *             and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::write( const std::string & file ) const {

		FILE *out;
		bool success;

		// Open a FILE pointer to the specified path ...
		out = IOTools::fopen(file.c_str(),"wb");
		if ( out == NULL ) {
			return false;
		}

		// ... and attempt to write this protected minutiae tempalte.
		write(out);

		// Check for errors ...
		if ( ferror(out) ) {
			success = false;
		} else {
			success = true;
		}

		fclose(out);

		// ... and return a boolean indicating whether this protected
		// minutiae template has been written successfully to the
		// specified file path.
		return success;
	}

	/**
	 * @brief
	 *             Initialize this protected minutiae record from the
	 *             specified file pointer.
	 *
	 * @details
	 *             If this protected minutiae record has been successfully
	 *             initialized, the function returns the number read bytes;
	 *             otherwise, if an error occurred, this protected minutiae
	 *             record will be left unchanged and the result of the
	 *             function will be -1.
	 *
	 * @param in
	 *             The file pointer from which the protected minutiae
	 *             record is read.
	 *
	 * @return
	 *             The number of bytes that were successfully read from
	 *             <code>in</code> to initialize this protected minutiae
	 *             record; otherwise, if an error occurred,
	 *             the function returns -1.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeRecord::read( FILE *in ) {

		int offset = 0;

		int c;
		ProtectedMinutiaeRecord tmp;

		{ // Check if header is existent and valid
			uint8_t header[10];
			for ( int j = 0 ; j < 10 ; j++ ) {
				if ( (c=fgetc(in)) < 0 ) { return -1; }
				header[j] = (uint8_t)c;
			}
			if ( header[0] != 'P' || header[1] != 'M' || header[2] != 'R' ||
				 header[9] != '\0' ) {
				return-1;
			}
			for ( int j = 3 ; j < 9 ; j++ ) {
				if ( !isdigit(header[j]) ) {
					return -1;
				}
			}
		}

		offset += 10;

		// Read and check 'width'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.width = c; tmp.width <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.width += c;
		if ( tmp.width == 0 ) { return -1; }
		offset += 2;

		// Read 'height'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.height = c; tmp.height <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.height += c;
		if ( tmp.height == 0 ) { return -1; }
		offset += 2;

		// Read and check 'dpi'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.dpi = c; tmp.dpi <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.dpi += c;
		if ( tmp.dpi < 300 || tmp.dpi > 1000 ) { return -1; }
		offset += 2;

		// Read and check 'gridDist'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.gridDist = c;
		if ( tmp.gridDist == 0 ) { return -1; }
		offset += 1;

		// Read and check 's'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.s = c;
		if ( tmp.s == 0 ) { return -1; }
		offset += 1;

		// Read and check 'k'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.k = c;
		if ( tmp.k == 0 ) { return -1; }
		offset += 1;

		// Read 'tmax' and check
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.tmax = c;
		if ( tmp.tmax == 0 || tmp.k > tmp.tmax ) { return -1; }
		offset += 1;

		// Read 'm'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.m  = c; tmp.m <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.m += c; tmp.m <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.m += c; tmp.m <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.m += c;
		offset += 4;

		// Read slow-down factor
		{
			int size;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			size  = c; size <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			size += c; size <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			size += c; size <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			size += c;
			offset += 4;

			uint8_t *data = (uint8_t*)malloc( size * sizeof(uint8_t) );
			if ( data == NULL ) {
				cerr << "ProtectedMinutiaeTemplate::read: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
			if ( fread(data,1,size,in) != (size_t)size ) {
				free(data);
				return -1;
			}
			offset += size;

			tmp.slowDownFactor.fromBytes(data,size);
			free(data);

			if ( tmp.slowDownFactor.sign() <= 0 ) {
				return -1;
			}
		}

		// Build the grid
		tmp.updateGrid();

		// Read the finite field
		{
			uint32_t rep;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			rep  = c; rep <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			rep += c; rep <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			rep += c; rep <<= 8;
			if ( (c=fgetc(in)) < 0 ) { return -1; }
			rep += c;
			offset += 4;

			if ( this->gfPtr != NULL && this->gfPtr->getDefiningPolynomial().rep == rep ) {
				tmp.gfPtr = new(nothrow) SmallBinaryField(this->gfPtr[0]);
			} else {
				tmp.gfPtr = new(nothrow) SmallBinaryField(SmallBinaryPolynomial(rep));
			}
			if ( tmp.gfPtr == NULL ) {
				cerr << "ProtectedMinutiaeTemplate::read: out of memory." << endl;
                exit(EXIT_FAILURE);
			}

			if ( (int)(tmp.gfPtr->getCardinality()) < this->s * (int)tmp.grid.size() ) {
				return -1;
			}
		}

		// Read 't'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.t = c; tmp.t <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.t += c;
		if ( tmp.t == 0 ) { return -1; }
		offset += 2;

		// *******************************************************************
		// ********************** BEGIN: Read data ***************************
		// *******************************************************************

		int n = tmp.vaultDataSize();

		uint8_t *vaultData = (uint8_t*)malloc( n * sizeof(uint8_t) );
		if ( vaultData == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::read: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		if ( (c=fgetc(in)) < 0 ) {
			free(vaultData);
			return -1;
		}
		offset += 1;

		if ( fread(vaultData,1,n,in) != (size_t) n ) {
			free(vaultData);
			return -1;
		}
		offset += n;

		if ( c != 0 ) {
			tmp.encryptedVaultPolynomialData = vaultData;
		} else {
			tmp.vaultPolynomialData = vaultData;
		}

		// *******************************************************************
		// ********************** END: Read data ***************************
		// *******************************************************************

		// Read 'hash'
		if ( fread(tmp.hash,1,20,in) != 20 ) {
			return -1;
		}
		offset += 20;

		tmp.updatePermutation();

		swap(*this,tmp);

		return offset;
	}

	/**
	 * @brief
	 *             Initialize this protected minutiae record from the
	 *             data contained in the specified file path.
	 *
	 * @details
	 *             If this protected minutiae record has been successfully
	 *             initialized, the function returns <code>true</code>;
	 *             otherwise, if an error occurred, this protected minutiae
	 *             record will be left unchanged and the result of the
	 *             function will be <code>false</code>.
	 *
	 * @param file
	 *             String representation of the path to the file from
	 *             which the protected minutiae record is read.
	 *
	 * @return
	 *             <code>true</code> if this protected minutiae record
	 *             has been successfully initialized by the data contained
	 *             in the specified file; otherwise, the result is
	 *             <code>false</code>.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeRecord::read( const std::string & file ) {

		bool success;

		// Open a FILE pointer to the specified path ...
		FILE *in = IOTools::fopen(file.c_str(),"rb");
		if ( in == NULL ) {
			return false;
		}

		// ... and attempt to initialize this protected minutiae template
		// by the data from the file.
		if ( read(in) >= 0 ) {
			success = true;
		} else {
			success = false;
		}

		fclose(in);

		// Return a boolean indicating whether this instance has successfully
		// been initialized by the data contained in the specified file.
		return success;
	}

	/**
	 * @brief
	 *             Swaps the content of two
	 *             \link ProtectedMinutiaeRecord\endlink.
	 *
	 * @details
	 *             After invocation, <code>vault1</code> holds the content of
	 *             the <code>vault2</code> at input and vice versa.
	 *
	 * @param vault1
	 *             A reference to a
	 *             \link ProtectedMinutiaeRecord\endlink
	 *
	 * @param vault2
	 *             A reference to a
	 *             \link ProtectedMinutiaeRecord\endlink
	 */
	void ProtectedMinutiaeRecord::swap
	( ProtectedMinutiaeRecord & vault1 , ProtectedMinutiaeRecord & vault2 ) {

		// Copies of the member variables of 'view1' into temporary
		// member variables
		int width = vault1.width;
		int height = vault1.height;
		int dpi = vault1.dpi;
		int gridDist = vault1.gridDist;
		int s = vault1.s;
		int k = vault1.k;
		int tmax = vault1.tmax;
		int m = vault1.m;
		std::vector< std::pair<double,double> > grid = vault1.grid;
		SmallBinaryField *gfPtr = vault1.gfPtr;
		int t = vault1.t;
		uint8_t *vaultPolynomialData = vault1.vaultPolynomialData;
		uint8_t *encryptedVaultPolynomialData = vault1.encryptedVaultPolynomialData;

		// Copy the members of 'view2' into 'view1'
		vault1.width = vault2.width;
		vault1.height = vault2.height;
		vault1.dpi = vault2.dpi;
		vault1.gridDist = vault2.gridDist;
		vault1.s = vault2.s;
		vault1.k = vault2.k;
		vault1.tmax = vault2.tmax;
		vault1.m = vault2.m;
		vault1.grid = vault2.grid;
		vault1.gfPtr = vault2.gfPtr;
		vault1.t =vault2.t;
		vault1.vaultPolynomialData = vault2.vaultPolynomialData;
		vault1.encryptedVaultPolynomialData = vault2.encryptedVaultPolynomialData;

		// Copy temporary member variables to 'view2'
		vault2.width = width;
		vault2.height = height;
		vault2.dpi = dpi;
		vault2.gridDist = gridDist;
		vault2.s = s;
		vault2.k = k;
		vault2.tmax = tmax;
		vault2.m = m;
		vault2.grid = grid;
		vault2.gfPtr = gfPtr;
		vault2.t = t;
		vault2.vaultPolynomialData = vaultPolynomialData;
		vault2.encryptedVaultPolynomialData = encryptedVaultPolynomialData;

		// SPECIAL CASE: Swap the 'slowDownFactor' fields using
		// the 'BigInteger::swap' method.
        vault1.slowDownFactor.swap(vault2.slowDownFactor);

		// SPECIAL CASE: swap the secret polynomials' hash values
		// via 'memcpy'
		uint8_t hash[20];
		memcpy(hash,vault1.hash,20);
		memcpy(vault1.hash,vault2.hash,20);
		memcpy(vault2.hash,hash,20);

		// SPECIAL CASE: Swap the 'permutation' fields using
		// the 'Permutation::swap' method.
        Permutation::swap(vault1.permutation,vault2.permutation);
	}

	/**
	 * @brief
	 *            Initializes this object's members such that it
	 *            represents an empty uninitialized protected minutiae
	 *            record.
	 *
	 * @details
	 *            This method is called by the standard and copy
	 *            constructor.
	 */
	void ProtectedMinutiaeRecord::first_init() {

		this->width = -1;
		this->height = -1;
		this->dpi = -1;

		this->gridDist = -1;
		this->s = 6;
		this->k = 25;
		this->tmax = 44;
		this->m = 3;
		this->slowDownFactor = BigInteger(1);

		this->gfPtr = NULL;
		this->t = -1;
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;

		memset(this->hash,0,20);
	}

	/**
	 * @brief
	 *            Update the vector of hexagonal grid points after
	 *            corresponding parameters have changed or specified
	 *            the very first time.
	 *
	 * @details
	 *            If the member variables \link width\endlink,
	 *            \link height\endlink, or \link gridDist\endlink
	 *            have changed, the
	 *            \link grid vector of hexagonal grid point\endlink
	 *            needs an update which is performed by calling
	 *            this method.
	 */
	void ProtectedMinutiaeRecord::updateGrid() {

		this->grid.clear();

		int maxRadius = (int)ceil(ceil
				(sqrt((double)(width*width+height*height)))+0.5*gridDist);

		// Create the points of a hexagonal grid of minimal distance 'gridDist'
		HexagonalGrid hexaGrid
			(-maxRadius,-maxRadius,maxRadius,maxRadius,this->gridDist);

		// Keep only those hexagonal grid points which are of a distance
		// more than those of which absolutely pre-aligned minutiae can occur.
		for ( int j = 0 ; j < (int)hexaGrid.getPoints().size() ; j++ ) {

			double dx , dy;
			dx = hexaGrid.getPoints().at(j).first;
			dy = hexaGrid.getPoints().at(j).second;

			if ( dx*dx+dy*dy <= maxRadius*maxRadius ) {
				this->grid.push_back(hexaGrid.getPoints().at(j));
			}
		}
	}

	/**
	 * @brief
	 *            Updates the finite field such that it can encode
	 *            each possible minutia's quantization.
	 *
	 * @details
	 *            After the \link grid hexagonal grid points\endlink
	 *            or the
	 *            \link s number of minutiae angle quantizations\endlink
	 *            have been specified, the finite field might need an
	 *            update which can be performed by calling this method.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeRecord::updateField() {

		// minimal size the field can have to be able to encode every minutia's
		// quantization
		int n = getVaultSize();

		int degree = 0;
		{ // Determine the minimal degree of the binary field
			int tmp = n;
			while ( tmp != 0 ) {
				++degree;
				tmp >>= 1;
			}
		}

		// Double the size of the field to allow the incorporation of
		// blending features
		degree += 1;

		if ( this->gfPtr == NULL || this->gfPtr->getDegree() != degree ) {

			// Delete the old field
			delete this->gfPtr;
			this->gfPtr = NULL;

			// Create new binary field of the desired degree
			this->gfPtr = new(nothrow) SmallBinaryField(degree);
			if ( this->gfPtr == NULL ) {
				cerr << "ProtectedMinutiaeRecord::updateField: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	/**
	 * @brief
	 *            Generate a public random user-specific permutation
	 *            process acting on the feature universe.
	 *
	 * @details
	 *            On enrollment, the feature set extracted from a
	 *            minutiae template have to be passed through a public
	 *            random user-specific permutation process to prevent
	 *            certain record multiplicity attacks
	 *            (see \link EEAAttack\endlink). The method generates
	 *            this permutation process by updating the
	 *            \link permutation\endlink field using a
	 *            \link RandomGenerator random generator\endlink whose
	 *            seed is derived from the hash value of the secret
	 *            polynomial (i.e., the \link hash\endlink field).
	 *            In such the permutation is randomly generated as the
	 *            secret polynomial, and thus, \link hash\endlink has
	 *            been randomly generated and is public such that on
	 *            verification the query feature sets can be passed
	 *            through the permutation process, too.
	 *
	 * @see permutation
	 */
	void ProtectedMinutiaeRecord::updatePermutation() {

		int n = getVaultSize();
        this->permutation.setDimension(n);

		RandomGenerator gen(this->hash);

		for ( int x0 = 0 ; x0 < n ; x0++ ) {
			int x1 = (int)(gen.rand() % n);
			this->permutation.exchange(x0,x1);
		}
	}

	/**
	 * @brief
	 *             Evaluate the reordering of a quantized feature which
	 *             is associated with an (enrolled) protected minutiae
	 *             template.
	 *
	 * @details
	 *             This function is not public and intended for internal
	 *             purposes. It is a variant of the
	 *             \link reorder() reorder(uint32_t)\endlink function but
	 *             does not ensure that \link isEnrolled()\endlink is
	 *             <code>true</code> thereby assuming that the
	 *             \link permutation\endlink field is initialized
	 *             properly. Otherwise, the function runs into
	 *             undocumented behavior.
	 *
	 * @param a
	 *             The feature element of which the reordering is
	 *             evaluated.
	 *
	 * @return
	 *             The reordering of <i>a</i>.
	 *
	 * @warning
	 *             If <i>a</i> does not encode a valid element in the
	 *             finite field that can be access via
	 *             \link getField()\endlink, an error message is printed
	 *             to <code>stderr</code> and the program exits with
     *             status 'EXIT_FAILURE'.
	 */
	uint32_t ProtectedMinutiaeRecord::_reorder( uint32_t a ) const {
		return (uint32_t)(this->permutation.eval((int)a));
	}

	/**
	 * @brief
	 *             Passes the specified polynomial to the
	 *             array \link vaultPolynomialData\endlink.
	 *
	 * @details
	 *             Let
	 *             \f[
	 *              V(X)=X^t+\sum_{j=0}^{t-1} V_j\cdot X^j
	 *             \f]
	 *             be the vault polynomial passed to this function.
	 *             Write \f$b_{0},...,b_{m-1}\f$
	 *             be the bit sequence resulting from a concatenation of
	 *             the coefficients \f$V_0,...,V_{t-1}\f$. Let
	 *             \f$b_{m},...,b_{M-1}\f$ be random bits such that
	 *             \f$M=96\cdot\lceil m/96\rceil\f$, i.e.,
	 *             \f$ b_{0},...,b_{m-1}\f$ are padded with random bits
	 *             such that the total number of bits becomes a multiple
	 *             of 96. An \link AES128 AES key\endlink is then derived
	 *             from a random positive integer smaller
	 *             than \link getSlowDownFactor()\endlink and used to
	 *             "encrypt" the bits \f$b_0,...,b_{M-1}\f$. In such we
	 *             obtain the "encrypted" bit sequence
	 *             \f[
	 *              (b'_0,...,b'_{M-1})
	 *             \f]
	 *             which is stored byte-wisely in
	 *             \link vaultPolynomialData\endlink
	 *
	 * @param V
	 *             The polynomial packed
	 *             into \link vaultPolynomialData\endlink.
	 *
	 *
	 * @warning
	 *             The field \link vaultPolynomialData\endlink will be
	 *             overwritten by newly allocated memory; thus, before
	 *             calling this method it must be ensured that
	 *             possibly allocated memory to
	 *             which \link vaultPolynomialData\endlink refers is freed;
	 *             otherwise, the program can cause memory leaks during
	 *             its runtime.
	 *
	 * @warning
	 *             The behavior of this method is undocumented
	 *             if \link t\endlink and \link gfPtr\endlink are not
	 *             initialized properly.
	 *
	 * @warning
	 *             If not enough memory can be provided, an error message
	 *             is printed <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'
	 *
	 * @warning
	 *             The polynomial represented by <code>V</code> should
	 *             be defined over the field \link getField()\endlink;
	 *             otherwise, the method runs into undocumented
	 *             behavior.
	 */
	void ProtectedMinutiaeRecord::packVaultPolynomial
	( const SmallBinaryFieldPolynomial & V ) {

		// Random slow-down value
		BigInteger slowDownVal; slowDownVal.random(this->slowDownFactor,true);

		// Derive AES key from slow-down value.
		AES128 aes = deriveKey(slowDownVal);

		int t , d , n;
		t = this->t;
		d = this->gfPtr->getDegree();
		n = vaultDataSize();

		// *******************************************************************
		// ***** BEGIN: Pack coefficients of vault polynomial into ***********
		// ***** encryption array ********************************************
		// *******************************************************************

		this->vaultPolynomialData = (uint8_t*)malloc
				( n * sizeof(uint8_t) );
		uint32_t *coeffs = (uint32_t*)malloc( t * sizeof(uint32_t) );
		if ( this->vaultPolynomialData == NULL || coeffs == NULL ) {
			cerr << "ProtectedMinutiaeView::slowDown: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		for ( int j = 0 ; j < t ; j++ ) {
			coeffs[j] = V.getCoeff(j);
		}

		concat_bit_vectors
			(this->vaultPolynomialData,coeffs,t,d,n);

		free(coeffs);

		// *******************************************************************
		// ***** END: Pack coefficients of vault polynomial into *************
		// ***** encryption array *******************************************
		// *******************************************************************

		// ENCRYPTION
		aes.encrypt
		(this->vaultPolynomialData,this->vaultPolynomialData,n);
	}

	/**
	 * @brief
     *            Determines how many bytes are needed to hold
	 *            the data of the packed vault polynomial.
	 *
	 * @details
	 *            The function is used to determine how many bytes
	 *            are needed to initialize an array of 128 bit blocks
	 *            that can hold the data of an instance of the improved
	 *            fuzzy vault scheme (encoded as a monic polynomial)
	 *            such that it can be protected via an
	 *            \link AES128 AES key\endlink.
	 *
	 * @return
	 *            Number of bytes needed to hold the data of a monic
	 *            polynomial of specified degree having coefficients
	 *            in a binary finite field of specified degree such
	 *            that the byte array is a multiple of 16 bytes
	 *            (i.e., 128 bits).
	 */
	int ProtectedMinutiaeRecord::vaultDataSize() const {

		int t , d;
		t = this->t;
		d = this->gfPtr->getDegree();

		int b = t * d;

		int n = b / 8 + (b%8?1:0);

		if ( n % 16 != 0 ) {
			n += 16-n%16;
		}

		return n;
	}

	/**
	 * @brief
	 *            Derives an AES key from the specified integer.
	 *
	 * @details
	 *            The bytes of the specified integer are extracted
	 *            via the \link BigInteger::toBytes() x.toBytes()\endlink
	 *            method and then its \link SHA SHA hash value\endlink
	 *            is computed; the first 16 bytes of the value
	 *            are used to define the \link AES128 AES key\endlink
	 *            returned by this function.
	 *
	 * @param x
	 *            The integer of which an AES key is computed with
	 *            this function.
	 *
	 * @return
	 *            The AES key derived from <code>x</code>.
	 *
	 * @warning
	 *             If not enough memory can be provided, an error message
	 *             is printed <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'
	 */
	AES128 ProtectedMinutiaeRecord::deriveKey( const BigInteger & x ) {

		uint8_t *array = (uint8_t*)malloc( x.getSizeInBytes() );
		if ( array == NULL ) {
			cerr << "createAESKey: out of memory."
				  << endl;
			exit(EXIT_FAILURE);
		}

		x.toBytes(array);
		uint8_t hash[20];
		SHA().hash(hash,array,x.getSizeInBytes());

		AES128 aes(hash);

		free(array);

		return aes;
	}

	/**
	 * @brief
	 *           Store the concatenation of a sequence of <i>d</i>-bit
	 *           vectors in a byte array and supplement it with random
	 *           bits.
	 *
	 * @details
	 *           The method is used to convert an instance of the improved
	 *           fuzzy vault scheme (encoded by <i>t</i>
	 *           elements of a binary finite field of degree <i>d</i>)
	 *           to a byte array that can be encrypted via an
	 *           \link AES128 AES key\endlink. The random bits are needed
	 *           to make the unencrypted array indistinguishable from falsely
	 *           decrypted data.
	 *
	 * @param bits
	 *           The concatenation of the bit vectors.
	 *
	 * @param vecs
	 *           Array of <i>t</i> 32-bit word each encoding a <i>d</i>-bit
	 *           vector.
	 *
	 * @param t
	 *           Number of 32-bit words contained in <code>array</code>.
	 *
	 * @param d
	 *           Number of significant bits of each 32-bit vector.
	 *
	 * @param n
	 *           Controls how many bits in <code>bits</code> are
	 *           supplemented by random values. More precisely, the latest
	 *           <i>8*n-t*d</i> bits in <code>bits</code> are selected
	 *           randomly.
	 */
	void ProtectedMinutiaeRecord::concat_bit_vectors
	( uint8_t *bits , const uint32_t *vecs , int t , int d , int n ) {

		int idx = 0; //Index of the active output byte starting from 0
		int count = 0; // Count of the number of bits

		// Set the output bytes initially to '0'.
		memset(bits,0,n);

		// Iteration over the array of 32-bit vector
		for( int i = 0 ; i < t ; i++ ) {

			// Iterate through the first 'd' bits of each vector starting
			// from the most significant down to the least significant.
			for ( int j = d-1 ; j >= 0 ; j-- ) {

				// Multiply the active byte by 2 ...
				bits[idx] <<= 1;

				// ... and if the currently considered bit is set ...
				if ( vecs[i] & (1<<j) ) {
					// ... increment.
					bits[idx] += 1;
				}

				// Check by bit count whether the index of the active
				// output byte to be changed.
				++count;
				if ( count % 8 == 0 ) {
					++idx;
				}
			}
		}

		// Supplement the least significant bits of the active
		// output byte by random bits.
		for ( ; count % 8 != 0 ; count++ ) {
			bits[idx] <<= 1;
			if ( (MathTools::rand8(true) & 0x1) != 0 ) {
				bits[idx] += 1;
			}
		}

		// Supplement the remaining bytes by random bytes
		++idx;
		for ( ; idx < n ; idx++ ) {
			bits[idx] = MathTools::rand8(true);
		}
	}

	/**
	 * @brief
	 *            Deconcatenate a bit sequence into <i>d</i>-bit strings
	 *            stored in 32 bit words.
	 *
	 * @details
	 *            The method is used to re-obtain the candidate of
	 *            the improved fuzzy vault instance data from a byte
	 *            array or after decryption.
	 *
	 * @param vecs
	 *            Output vector that will contain <i>t</i> 32-bit words
	 *            of <i>d</i> significant bits.
	 *
	 * @param bits
	 *            Byte array encoding the <i>t*d</i> bit sequence that
	 *            is deconcatenated into <code>vecs</code>.
	 *
	 * @param t
	 *            Number of <i>d</i>-bit strings stored in the bit sequence
	 *            <code>bits</code>.
	 *
	 * @param d
	 *            Positive integer smaller than or equals 32 controlling
	 *            the length of the bit strings encoded in the sequence
	 *            <code>bits</code>.
	 */
	void ProtectedMinutiaeRecord::split_into_bit_vectors
	( uint32_t *vecs , const uint8_t *bits , int t , int d  ) {

		int idx = 0; //Active byte index
		int l = 7; // Active bit index

		// Iteration over the 'd'-bit bitstrings.
		for ( int i = 0 ; i < t ; i++ ) {

			// Initialize bitstring with '0'.
			vecs[i] = 0;

			// Iterate over the bitstring's position
			for ( int j = d-1 ; j >= 0 ; j-- ) {

				// Multiply active bitstring by 2 ...
				vecs[i] <<= 1;
				// ... and if the active bit index is set, ...
				if ( bits[idx] & (1<<l) ) {
					// ... increment.
					vecs[i] += 1;
				}

				// Check whether active byte index to be changed.
				--l;
				if ( l < 0 ) {
					++idx;
					l = 7;
				}
			}
		}
	}

}

