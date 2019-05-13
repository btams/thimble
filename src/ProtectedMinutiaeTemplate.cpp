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
 * @file ProtectedMinutiaeTemplate.cpp
 *
 * @brief
 *            Implements functions and methods provided by
 *            'ProtectedMinutiaeTemplate.h' which, in particular, provides a
 *            class for representing protected minutiae templates of a single
 *            finger that are absolutely pre-algined w.r.t. a directed reference
 *            point.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <stdint.h>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <vector>
#include <iostream>

#include "config.h"
#include <thimble/math/HexagonalGrid.h>
#include <thimble/math/Permutation.h>
#include <thimble/math/RandomGenerator.h>
#include <thimble/security/SHA.h>
#include <thimble/security/AES.h>
#include <thimble/security/FuzzyVaultTools.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/finger/ProtectedMinutiaeTemplate.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {


	/**
	 * @brief
	 *            Standard constructor.
	 *
	 * @details
	 *            Creates an empty vault in which no member fields
	 *            are set. The created object will not be useful
	 *            for protect an absolutely pre-aligned minutiae
	 *            before initialization, e.g., via
	 *            <pre>
	 *             ProtectedMinutiaeTemplate vault;
	 *             vault = ProtectedMinutiaeTemplate(width,height);
	 *            </pre>
	 *            in which the width and the height of the fingerprint
	 *            image is specified.
	 */
	ProtectedMinutiaeTemplate::ProtectedMinutiaeTemplate() {
		this->width = -1;
		this->height = -1;
		this->fingerPosition = UNKNOWN_FINGER;
		this->dpi = -1;
		this->gridDist = -1;
		this->s = -1;
		this->k = -1;
		this->tmax = -1;
		this->D = -1;
		this->gfPtr = NULL;
		this->t = -1;
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;
		memset(this->hash,0,20);
		this->is_initialized = false;
	}

	/**
	 * @brief
	 *             Constructor.
	 *
	 * @details
	 *             The constructor generates parameters needed to properly
	 *             quantize absolutely pre-aligned minutiae templates that
	 *             have been estimated from a fingerprint image of the
	 *             specified dimension and that has been scanned at the
	 *             specified resolution. If the finger's position is known,
	 *             it can be specified as well.
	 *
	 * @param width
	 *             The width of the fingerprint image.
	 *
	 * @param height
	 *             The height of the fingerprint image.
	 *
	 * @param dpi
	 *             The resolution in (dots per inch) at which
	 *             the fingerprint image has been scanned.
	 *
	 * @param fingerPosition
	 *             The position of the finger (if known).
	 *
	 * @attention
	 *             If <code>width</code> or <code>height</code> are
	 *             not greater than 0 or if <code>dpi</code> does
	 *             not range between 300 and 1000, an error message
	 *             is printed to <code>stderr</code> and the program
     *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	ProtectedMinutiaeTemplate::ProtectedMinutiaeTemplate
	( int width , int height , int dpi , FINGER_POSITION_T fingerPosition ) {

		first_init(width,height,dpi,fingerPosition);
	}

	/**
	 * @brief
	 *             Copy constructor
	 *
	 * @param vault
	 *             The protected minutiae template of which
	 *             a copy is constructed.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	ProtectedMinutiaeTemplate::ProtectedMinutiaeTemplate
	( const ProtectedMinutiaeTemplate & vault ) {

		this->width = -1;
		this->height = -1;
		this->fingerPosition = UNKNOWN_FINGER;
		this->dpi = -1;
		this->gridDist = -1;
		this->s = -1;
		this->k = -1;
		this->tmax = -1;
		this->D = -1;
		this->gfPtr = NULL;
		this->t = -1;
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;
		memset(this->hash,0,20);
		this->is_initialized = false;

		*this = vault;
	}

	/**
	 * @brief
	 *             Destructor
	 *
	 * @details
	 *             Frees all the memory helt by the protected minutiae
	 *             template.
	 */
	ProtectedMinutiaeTemplate::~ProtectedMinutiaeTemplate() {

		clear();
		delete this->gfPtr;
	}

	/**
	 * @brief
	 *             Initializes this object such that it can protected
	 *             absolutely pre-aligned minutiae templates.
	 *
	 * @details
	 *             This method generates parameters needed to properly
	 *             quantize absolutely pre-aligned minutiae templates that
	 *             have been estimated from a fingerprint image of the
	 *             specified dimension and that has been scanned at the
	 *             specified resolution. If the finger's position is known,
	 *             it can be specified as well.
	 *
	 *             Any data hold by this object will be deleted and
	 *             replaced by new data.
	 *
	 * @param width
	 *             The width of the fingerprint image.
	 *
	 * @param height
	 *             The height of the fingerprint image.
	 *
	 * @param dpi
	 *             The resolution in (dots per inch) at which
	 *             the fingerprint image has been scanned.
	 *
	 * @param fingerPosition
	 *             The position of the finger (if known).
	 *
	 * @attention
	 *             If <code>width</code> or <code>height</code> are
	 *             not greater than 0 or if <code>dpi</code> does
	 *             not range between 300 and 1000, an error message
	 *             is printed to <code>stderr</code> and the program
	 *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::initialize
	( int width , int height ,
	  int dpi , FINGER_POSITION_T fingerPosition ) {

		clear();
		delete this->gfPtr;

		first_init(width,height,dpi,fingerPosition);
	}

	/**
	 * @brief
	 *             Assignment operator
	 *
	 * @details
	 *             Sets the \link ProtectedMinutiaeTemplate\endlink to a
	 *             copy of <code>vault</code>.
	 *
	 * @param vault
	 *             A \link ProtectedMinutiaeTemplate\endlink of which this
	 *             instance becomes a copy of.
	 *
	 * @return
	 *             A reference to this
	 *             \link ProtectedMinutiaeTemplate\endlink
	 *             (after the assignment side-effect).
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	ProtectedMinutiaeTemplate & ProtectedMinutiaeTemplate::operator=
			( const ProtectedMinutiaeTemplate & vault ) {

		if ( this != &vault ) { // Do only execute if not a self assignment.

			// Copy primitive member variables
			this->width = vault.width;
			this->height = vault.height;
			this->fingerPosition = vault.fingerPosition;
			this->dpi = vault.dpi;
			this->gridDist = vault.gridDist;
			this->s = vault.s;
			this->k = vault.k;
			this->tmax = vault.tmax;
			this->D = vault.D;
			this->t = vault.t;
			this->slowDownFactor = vault.slowDownFactor;

			// Adopt assignment operator of the 'std::vector' class
			this->grid = vault.grid;


			// Copy underlying finite field
			if ( this->gfPtr == NULL && vault.gfPtr != NULL ) {
				this->gfPtr = new (nothrow) SmallBinaryField(vault.gfPtr[0]);
				if ( this->gfPtr == NULL ) {
					cerr << "ProtectedMinutiaeTemplate: out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
			} else if ( this->gfPtr != NULL && vault.gfPtr != NULL ){
				this->gfPtr[0] = vault.gfPtr[0];
			} else if ( this->gfPtr != NULL && vault.gfPtr == NULL ) {
				delete this->gfPtr;
				this->gfPtr = NULL;
			}

			// Copy vault down data
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
					cerr << "ProtectedMinutiaeTemplate: out of memory." << endl;
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
					cerr << "ProtectedMinutiaeTemplate: out of memory." << endl;
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

			this->is_initialized = vault.is_initialized;
		}

		return *this;
	}

	/**
	 * @brief
	 *             Swaps the content of two
	 *             \link ProtectedMinutiaeTemplate\endlink.
	 *
	 * @details
	 *             After invocation, <code>vault1</code> holds the content of
	 *             the <code>vault2</code> at input and vice versa.
	 *
	 * @param vault1
	 *             A reference to a valid
	 *             \link ProtectedMinutiaeTemplate\endlink
	 *
	 * @param vault2
	 *             A reference to a valid
	 *             \link ProtectedMinutiaeTemplate\endlink
	 */
	void ProtectedMinutiaeTemplate::swap
	( ProtectedMinutiaeTemplate & vault1 , ProtectedMinutiaeTemplate & vault2 ) {

		// Copies of the member variables of 'vault1' into temporary
		// member variables
		int width = vault1.width;
		int height = vault1.height;
		FINGER_POSITION_T fingerPosition = vault1.fingerPosition;
		int dpi = vault1.dpi;
		int gridDist = vault1.gridDist;
		int s = vault1.s;
		int k = vault1.k;
		int tmax = vault1.tmax;
		int D = vault1.D;
		std::vector< std::pair<double,double> > grid = vault1.grid;
		SmallBinaryField *gfPtr = vault1.gfPtr;
		int t = vault1.t;
		uint8_t *vaultPolynomialData = vault1.vaultPolynomialData;
		uint8_t *encryptedVaultPolynomialData =
				vault1.encryptedVaultPolynomialData;
		bool is_initialized = vault1.is_initialized;

		// Copy the members of 'vault2' into 'vault1'
		vault1.width = vault2.width;
		vault1.height = vault2.height;
		vault1.fingerPosition = vault2.fingerPosition;
		vault1.dpi = vault2.dpi;
		vault1.gridDist = vault2.gridDist;
		vault1.s = vault2.s;
		vault1.k = vault2.k;
		vault1.tmax = vault2.tmax;
		vault1.D = vault2.D;
		vault1.grid = vault2.grid;
		vault1.gfPtr = vault2.gfPtr;
		vault1.t = vault2.t;
		vault1.vaultPolynomialData = vault2.vaultPolynomialData;
		vault1.encryptedVaultPolynomialData =
				vault2.encryptedVaultPolynomialData;
		vault1.is_initialized = vault2.is_initialized;


		// Copy temporary member variables to 'vault2'
		vault2.width = width;
		vault2.height = height;
		vault2.fingerPosition = fingerPosition;
		vault2.dpi = dpi;
		vault2.gridDist = gridDist;
		vault2.s = s;
		vault2.k = k;
		vault2.tmax = tmax;
		vault2.D = D;
		vault2.grid = grid;
		vault2.gfPtr = gfPtr;
		vault2.t = t;
		vault2.vaultPolynomialData = vaultPolynomialData;
		vault2.encryptedVaultPolynomialData =
				encryptedVaultPolynomialData;
		vault2.is_initialized = is_initialized;

		// SPECIAL CASE: swap the secret polynomials' hash values
		// via 'memcpy'
		uint8_t hash[20];
		memcpy(hash,vault1.hash,20);
		memcpy(vault1.hash,vault2.hash,20);
		memcpy(vault2.hash,hash,20);
		// SPECIAL CASE: Swap the 'permutation' fields using
		// the 'Permutation::swap' method.
        Permutation::swap(vault1.permutation,vault2.permutation);
		// SPECIAL CASE: Swap the 'slowDownFactor' fields using
		// the 'BigInteger::swap' method.
        vault1.slowDownFactor.swap(vault2.slowDownFactor);
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
	 *            constructor \link ProtectedMinutiaeTemplate()\endlink
	 *            the width, height, resolution to which a to-be-protected
	 *            minutiae view corresponds to are unknown. However, this
	 *            information must be known for the vault's attributed to
	 *            be initialized properly and the member is_initialized
	 *            equals <code>false</code> as long as they have not
	 *            specified and, otherwise, <code>true</code>.
	 *
	 * @return
	 *            <code>true</code> if essential information about the
	 *            fingerprint images from which (to-be-)protected minutiae
	 *            templates have been extracted and, otherwise,
	 *            <code>false</code>.
	 */
	bool ProtectedMinutiaeTemplate::isInitialized() const {

		return this->is_initialized;
	}

	/**
	 * @brief
	 *            Determine whether this instance does hold a
	 *            protected minutiae template.
	 *
	 * @return
	 *            <code>true</code> if this
	 *            \link ProtectedMinutiaeTemplate\endlink does protect
	 *            a minutiae template; otherwise <code>false</code>.
	 *
	 * @see ProtectedMinutiaeTemplate::enroll()
	 */
	bool ProtectedMinutiaeTemplate::isEnrolled() const {
		return this->t >= 0;
	}

	/**
	 * @brief
	 *            Determine whether this instance does hold a
	 *            protected minutiae template that is additionally
	 *            protected by a (user-specific) key.
	 *
	 * @details
	 *            If this function returns <code>true</code>, the
	 *            \link ProtectedMinutiaeTemplate\endlink has been
	 *            encrypted via
	 *            \link ProtectedMinutiaeTemplate::encrypt\endlink;
	 *            otherwise, the function returns <code>false</code>.
	 *            If this \link ProtectedMinutiaeTemplate\endlink
	 *            has been decrypted via
	 *     \link ProtectedMinutiaeTemplate::decrypt(const AES128&)\endlink
	 *            such that it represents a candidate
	 *            for a decrypted \link ProtectedMinutiaeTemplate\endlink,
	 *            the function returns <code>false</code>.
	 *
	 * @return
	 *            <code>true</code> if this
	 *            \link ProtectedMinutiaeTemplate\endlink does protect
	 *            a minutiae template, is additionally protected
	 *            by an \link AES128\endlink key and does not represent
	 *            a candidate for a decrypted
	 *            \link ProtectedMinutiaeTemplate\endlink; otherwise, the
	 *            function returns <code>false</code>.
	 *
	 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
	 */
	bool ProtectedMinutiaeTemplate::isEncrypted() const {
		return this->encryptedVaultPolynomialData != NULL &&
				this->vaultPolynomialData == NULL;
	}

	/**
	 * @brief
	 *            Determine whether this instance holds a decrypted
	 *            candidate for a protected minutiae template.
	 *
	 * @details
	 *            Any encrypted \link ProtectedMinutiaeTemplate\endlink
	 *            can be decrypted via an \link AES128\endlink key and,
	 *            subsequently, contains a candidate for the decrypted
	 *            \link ProtectedMinutiaeTemplate\endlink.
	 *
	 *            Note that, if this
	 *            \link ProtectedMinutiaeTemplate\endlink
	 *            has not been decrypted via
	 *    \link ProtectedMinutiaeTemplate::decrypt(const AES128&)\endlink
	 *            but contains enrolled data, the function returns
	 *            <code>true</code>.
	 *
	 * @return
	 *            <code>true</code> if this
	 *            \link ProtectedMinutiaeTemplate\endlink contains a
	 *            candidate for the unencrypted
	 *            \link ProtectedMinutiaeTemplate\endlink; otherwise,
	 *            the result will be <code>false</code>.
	 *
	 * @see ProtectedMinutiaeTemplate::decrypt(const AES128&)
	 */
	bool ProtectedMinutiaeTemplate::isDecrypted() const {
		return this->vaultPolynomialData != NULL;
	}

	/**
	 * @brief
	 *            Determine whether this instance does hold
	 *            encrypted data of a protected minutiae template.
	 *
	 * @details
	 *            Note that a \link ProtectedMinutiaeTemplate\endlink
	 *            keeps the encrypted data even after being decrypted
	 *            via \link ProtectedMinutiaeTemplate::decrypt\endlink
	 *
	 * @return
	 *            <code>true</code> if the
	 *            \link ProtectedMinutiaeTemplate\endlink contains
	 *            encrypted data; otherwise, the function returns
	 *            <code>false</code>.
	 *
	 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
	 */
	bool ProtectedMinutiaeTemplate::containsEncryptedData() const {
		return this->encryptedVaultPolynomialData != NULL;
	}

	/**
	 * @brief
	 *            Create protected data from a minutiae template.
	 *
	 * @details
	 *            The function passes the input minutiae template
	 *            through a quantization process and protects the
	 *            quantization via the improved fuzzy vault scheme
	 *            of which data will be held by this instance.
	 *
	 *            If not sufficient information could be extracted
	 *            from the input minutiae template to allow successful
	 *            verification, the function returns <code>false</code>;
	 *            otherwise, if successul, the function returns
	 *            <code>true</code>.
	 *
	 *
	 * @param view
	 *            A reference to a valid \link MinutiaeView\endlink
	 *            representing the input (absolutely pre-aligned)
	 *            minutiae template.
	 *
	 * @return
	 *            <code>true</code> if the creation was successful;
	 *            otherwise, the function returns <code>false</code>.
	 *
	 * @warning
	 *            If \link isEnrolled()\endlink returns
	 *            <code>true</code>, an error message is printed to
	 *            <code>stderr</code> and the program exits with
     *            status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeTemplate::enroll( const MinutiaeView & view ) {

		// Check if essential information about the original minutiae
		// template, such as width and height, are known.
		if ( !this->is_initialized ) {
			cerr << "ProtectedMinutiaeTemplate::enroll: "
				 << "object not initialized." << endl;
			exit(EXIT_FAILURE);
		}

		// If already enrolled, ...
		if ( isEnrolled() ) {
			// ... print an error message and ...
			cerr << "ProtectedMinutiaeTemplate::enroll: "
				 << "already enrolled, clear first." << endl;
            // .. exit with status 'EXIT_FAILURE'.
            exit(EXIT_FAILURE);
		}

		int n = getVaultSize();

		// Temporarily generate secret polynomial and
		SmallBinaryFieldPolynomial f(getField());
		f.random(this->k,true);

		// save its SHA-1 hash value
        SHA().hash(this->hash,f.getData(),f.deg()+1);

		// Generate user-specific public permutation process
		updatePermutation();

		// Allocate memory that can hold the minutiae template's quantization
		uint32_t *x = (uint32_t*)malloc( this->tmax * sizeof(uint32_t) );
		if ( x == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::enroll: out of memory."
				 << endl;
            exit(EXIT_FAILURE);
		}

		// Extract feature set from input minutiae template
		int t = quantize(x,view);

		// Check whether feature set is large enough to
		// allow successful genuine verification.
		if ( t < this->k ) {
			return false;
		}

		// Apply the random public user-specific permutation process
		// to the feature set.
		for ( int j = 0 ; j < t ; j++ ) {
			x[j] = _reorder(x[j]);
		}

		// Supplement the feature set with blending features.
		for ( int j = t ; j < this->tmax ; j++ ) {
			x[j] = n + MathTools::rand32(true) %
					(getField().getCardinality() - n );
		}

		// The number of feature elements (including blending features)
		t = this->tmax;

		// Build characteristic polynomial of feature set which is the
		// unique monic polynomial having the feature elements as roots.
		SmallBinaryFieldPolynomial V(getField());
		V.buildFromRoots(x,t);

		// Obfuscate the characteristic polynomial by the secret polynomial
		// and vice versa by building the 'vault polynomial'
		add(V,V,f);

		// Save degree of vault polynomial in the corresponding
		// member variable.
		this->t = V.deg();

		// Free memory temporarily holding the minutiae template's
		// quantization
		free(x);

		packVaultPolynomial(V);

		// Successful enrollment
		return true;
	}

	/**
	 * @brief
	 *             Performs a verification procedure to compare
	 *             the input minutiae template with the minutiae
	 *             template protected by this instance.
	 *
	 * @details
	 *             The function passes the input minutiae template
	 *             through a quantization process from which an instance
	 *             of the Reed-Solomon decoding problem is generated which
	 *             is attempted to be solved using the
	 * \link ProtectedMinutiaeTemplate::decode()\endlink
	 *             function. If successful, the verification is deemed
	 *             to be successful and, otherwise, failed.
	 *
	 * @param view
	 *             The (absolutely pre-aligned) query
	 *             minutiae template.
	 *
	 * @return
	 *             <code>true</code> if the verification procedure
	 *             was successful; otherwise, if the verification
	 *             procedure was not successful, the function returns
	 *             <code>false</code>.
	 *
	 * @warning
	 *             If this \link ProtectedMinutiaeTemplate\endlink
     *             does not represent a successfully enrolled
     *             decrypted protected minutiae template, i.e., if
     *      \link ProtectedMinutiaeTemplate::isDecrypted()\endlink
     *             returns <code>false</code>, the function
     *             prints an error message to <code>stderr</code> and
     *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeTemplate::verify( const MinutiaeView & view ) const {
		SmallBinaryFieldPolynomial f(getField());
		return open(f,view);
	}

	/**
	 * @brief
	 *             Encrypts this instance using a key.
	 *
	 * @details
	 *             After successful enrollment (via
	 *             \link ProtectedMinutiaeTemplate::enroll()\endlink),
	 *             the \link ProtectedMinutiaeTemplate\endlink can
	 *             additionally protected via a (user-specific)
	 *             key derived from a PIN or password. Consequently,
	 *             the security level of the
	 *             \link ProtectedMinutiaeTemplate\endlink can be
	 *             improved, though this requires the
	 *             \link ProtectedMinutiaeTemplate\endlink to be
	 *             decrypted correctly for a successful verification,
	 *             thereby requiring that the correct can be
	 *             reproduced for successful verification.
	 *
	 * @param key
	 *             The \link AES128\endlink key used to encrypt the
	 *             \link ProtectedMinutiaeTemplate\endlink.
	 *
	 * @warning
	 *             If no minutiae template has been enrolled with
	 *             this \link ProtectedMinutiaeTemplate\endlink
	 *             (i.e., if
	 * \link ProtectedMinutiaeTemplate::isEnrolled()\endlink
	 *             returns <code>false</code>)
	 *             or if it contains already encrypted data (i.e., if
	 * \link ProtectedMinutiaeTemplate::containsEncryptedData()\endlink
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
	void ProtectedMinutiaeTemplate::encrypt( const AES128 & key ) {

		// Ensure that this instance does protect a minutiae template.
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::encrypt: no vault built." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure that this instance does not already hold encrypted data
		if ( isEncrypted() ) {
			cerr << "ProtectedMinutiaeTemplate::encrypt: already encrypted."
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
	 *             \link ProtectedMinutiaeTemplate\endlink that has been
	 *             encrypted before using
	 *             \link ProtectedMinutiaeTemplate::encrypt(const AES128&)\endlink.
	 *
	 *             If the \link ProtectedMinutiaeTemplate\endlink has been
	 *             encrypted with the key \f$\kappa\f$, then decryption
	 *             with \f$\kappa\f$ will correctly reveal the original
	 *             data; otherwise, if decrypting with a different
	 *             key \f$\kappa'\neq\kappa\f$, this instance will hold
	 *             a candidate for the original data that is (with very
	 *             high probability) different from the original data.
	 *             To date no efficient way for distinguishing incorrect
	 *             candidate data from the correct data is known. Hence,
	 *             the difficulty to break an encrypted
	 *             \link ProtectedMinutiaeTemplate\endlink equals the product
	 *             of the difficulty in guessing a  minutiae template being
	 *             of sufficient similarity to the protected minutiae template
	 *             and the difficulty in guessing the correct password.
	 *
	 * @param key
	 *             The \link AES128\endlink key used to decrypt this
	 *             \link ProtectedMinutiaeTemplate\endlink.
	 *
	 * @warning
	 *             If this instance does not hold encrypted data
	 *             (i.e., if
	 *    \link ProtectedMinutiaeTemplate::containsEncryptedData()\endlink
	 *             returns <code>false</code>), an error message is
	 *             printed to <code>stderr</code> and the program exits
	 *             with status 'EXIT_FAILURE'.
	 *
	 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::decrypt( const AES128 & key ) {

		// Ensure that this instance does contain encrypted data.
		if ( !containsEncryptedData() ) {
			cerr << "ProtectedMinutiaeTemplate::decrypt: "
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
				cerr << "ProtectedMinutiaeTemplate::decrypt: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		aes.decrypt
		(this->vaultPolynomialData,this->encryptedVaultPolynomialData,n);
	}

	/**
	 * @brief
	 *             Use a minutiae template as a query to open the
	 *             vault.
	 *
	 * @details
	 *             This a special version of
	 *             the \link verify()\endlink function (in
	 *             fact, \link verify()\endlink
	 *             wraps around this function) in which the
	 *             secret polynomial that is used to hide the minutiae template
	 *             protected by this instance is revealed on successful
	 *             verification.
	 *
	 * @param f
	 *             Will be set to the secret polynomial on successful
	 *             verification; otherwise the content will be left
	 *             unchanged
	 *
	 * @param view
	 *             The (absolutely pre-aligned) query minutiae template.
	 *
	 * @return
	 *             <code>true</code> if the verification procedure was
	 *             successful; otherwise, if the verification was not
	 *             successful, the function returns <code>false</code>.
	 *
	 * @warning
	 *             If this \link ProtectedMinutiaeTemplate\endlink
     *             does not represent a successfully enrolled
     *             decrypted protected minutiae template, i.e.,
     *             if \link ProtectedMinutiaeTemplate::isDecrypted()\endlink
     *             returns <code>false</code>, the function
     *             prints an error message to <code>stderr</code> and
     *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeTemplate::open
	( SmallBinaryFieldPolynomial & f , const MinutiaeView & view ) const {

		// Allocate memory to temporarily hold the feature set.
		uint32_t *x = (uint32_t*)malloc( this->tmax * sizeof(uint32_t) );
		if ( x == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::open: out of memory." << endl;
            exit(EXIT_FAILURE);
		}

		// Extract the feature set and ...
		int t = quantize(x,view);

		// ... attempt to open with the feature set.
		bool success = open(f,x,t);

		// Free temporarily allocated memory.
		free(x);

		return success;
	}

	/**
	 * @brief
	 *             Use a quantized minutiae template as a query to open the
	 *             vault.
	 *
	 * @details
	 *             This is a variant of the \link open()\endlink
	 *             function in which a set of quantized minutiae is
	 *             used to open the vault. Therein, a quantized minutia is
	 *             encoded as an element of the finite field that can be
	 *             accessed via \link getField()\endlink.
     *
     *             To obtain the feature set of an absolutely pre-aligned
     *             minutiae template, the \link quantize()\endlink method
     *             can be used.
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
	 *             If this \link ProtectedMinutiaeTemplate\endlink
     *             does not represent a successfully enrolled and
     *             decrypted protected minutiae template, i.e.,
     *             if \link isDecrypted()\endlink
     *             returns <code>false</code>, the function
     *             prints an error message to <code>stderr</code> and
     *             exits with status 'EXIT_FAILURE'.
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
	bool ProtectedMinutiaeTemplate::open
	( SmallBinaryFieldPolynomial & f , const uint32_t *B , int t ) const {

		// Ensure that this instance does protect a feature set and ...
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::open: "
				 << "no minutiae template protected by this view." << endl;
            exit(EXIT_FAILURE);
		}

		// ... contains a decrypted polynomial.
		if ( isEncrypted() ) {
			cerr << "ProtectedMinutiaeTemplate::open: "
				 << "vault is encrypted; decrypt first." << endl;
            exit(EXIT_FAILURE);
		}

		// Allocate memory to hold the set of unlocking pairs '{(x[j],y[j])}'
		uint32_t  *x , *y;
		x = (uint32_t*)malloc( t * sizeof(uint32_t) );
		y = (uint32_t*)malloc( t * sizeof(uint32_t) );
		if ( x == NULL || y == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::open: "
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
			success = decode(f,x,y,t,this->k,this->hash,this->D);

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
	 *             Implementation of the <i>randomized decoder</i> to
	 *             solve an instance of the Reed-Solomon decoding problem.
	 *
	 * @details
	 *             The function randomly selects <i>k</i> distinct
	 *             pairs <i>(x[j],y[j])</i>, computes its interpolation
	 *             polynomial, checks if its 160 bit hash value agrees
	 *             with <i>hash</i>, and, if true, assigns <i>f</i> by
	 *             the interpolation polynomial returning
	 *             <code>true</code>; otherwise, the function repeats
	 *             the procedure up to <i>D</i> times and returns
	 *             <code>false</code> if none of them yielded the
	 *             correct hash value.
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
	 * @param D
	 *             Number of decoding iterations
	 *
	 * @return
	 *             <code>true</code> if the decoding attempt was successful;
	 *             otherwise, the function returns <code>false</code>.
	 */
	bool ProtectedMinutiaeTemplate::decode
	( SmallBinaryFieldPolynomial & f , const uint32_t *x , const uint32_t *y ,
	  int t , int k , const uint8_t hash[20] , int D ) const {

		// Check whether unlocking set is of size that can possibly
		// be decoded.
		if ( t < k ) {
			return false;
		}

		// Essentially, the randomized decoder consists
		// of the first 'D' steps of a randomized brute-force
		// attack.
		return FuzzyVaultTools::bfattack(f,x,y,t,k,hash,D);
	}

	/**
	 * @brief
	 *             Computes the quantization of an (absolutely pre-aligned)
	 *             minutia.
	 *
	 * @param minutia
	 *             The minutia of which the quantization is computed.
	 *
	 * @return
	 *             The quantization of <code>minutia</code> encoded as a
	 *             finite field element.
	 */
	uint32_t ProtectedMinutiaeTemplate::quantize( const Minutia & minutia ) const {

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
		return (uint32_t)i+(uint32_t)j*(uint32_t)(this->grid.size());
	}

	/**
	 * @brief
	 *             Compute the quantization of an (absolutely pre-aligned)
	 *             minutiae template encoded as a set of finite field
	 *             elements.
	 *
	 * @details
	 *             At most the first <code>tmax</code> quantizations of
	 *             the minutiae in <code>view</code> are computed and
	 *             added to <code>array</code> (if not contained already).
	 *             The order of minutiae is given by their estimated
	 *             quality (if any), thereby guaranteeing that the minutiae
	 *             of the highest estimated quality are quantized first and
	 *             so on.
	 *
	 *             The function returns the number <i>t</i> of distinct
	 *             minutiae  quantizations that have been put into
	 *             <code>array</code> which depends on the size of the input
	 *             minutiae template
	 *             (see \link MinutiaeView::getMinutiaeCount()\endlink)
	 *             , <code>tmax</code> and on the count how many minutiae
	 *             had equal quantizations. Thereby, the function
	 *             guarantees that the first <i>t</i> elements of
	 *             <code>array</code> are distinct.
	 *
	 * @param array
	 *             Will contain at most <code>tmax</code>
	 *             finite field elements encoding the first minutiae'
	 *             quantization in <code>view</code>.
	 *
	 * @param view
	 *             The minutiae template of which the quantization set
	 *             is computed.
	 *
	 * @param tmax
	 *             Maximal number of minutiae quantizations to be extracted
	 *             from <code>view</code>.
	 *
	 * @return
	 *             The number of minutiae quantizations extracted from
	 *             <code>view</code>.
	 *
	 * @warning
	 *             If <code>tmax</code> is not positive, an error message
	 *             is printed to <code>stderr</code> and the program exits
     *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>array</code> cannot hold <code>tmax</code>
	 *             elements of type <code>uint32_t</code> calling this
	 *             function runs into undocumented behavior.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeTemplate::quantize
	( uint32_t *array , const MinutiaeView & view  , int tmax ) const {

		// Ensure that the bound on the maximal number of minutiae
		// quantizations is positive.
		if ( tmax < 0  ) {
			cerr << "ProtectedMinutiaeTemplate::quantize: "
				 << "bound on maximal number of minutiae quantizations must "
				 << "be positive." << endl;
            exit(EXIT_FAILURE);
		}

		// Copy the minutiae template and ..
		MinutiaeView w = view;
		// ... sort their minutiae w.r.t. their estimated qualities.
		w.sortWithRespectToMinutiaeQuality();

		int t = 0; // Count of minutiae quantizations

		// Iterate over the minutiae in 'w' until the minutiae quantization
		// count reaches the bound 'tmax'
		for ( int j = 0 ; j < w.getMinutiaeCount() && t < tmax ; j++ ) {

			// Quantize the current minutia and ...
			uint32_t q = quantize(w.getMinutia(j));

			// ... determine whether the quantization is already contained
			// in the output array; ...
			bool alreadyContained = false;

			for ( int _j = 0 ; _j < t ; _j++ ) {
				if ( array[_j] == q ) {
					alreadyContained = true;
					break;
				}
			}

			// ...; if not, append it in the output array.
			if ( !alreadyContained ) {
				array[t++] = q;
			}
		}

		// Return the number of differend minutiae quantizations extraced
		// from 'view'.
		return t;
	}

	/**
	 * @brief
	 *             Compute the quantization of an (absolutely pre-aligned)
	 *             minutiae template encoded as a set of finite field
	 *             elements.
	 *
	 * @details
	 *             At most the first <code>tmax</code> quantizations of
	 *             the minutiae in <code>view</code> are computed and
	 *             added to <code>array</code> (if not contained already).
	 *             The order of minutiae is given by their estimated
	 *             quality (if any), thereby guaranteeing that the minutiae
	 *             of the highest estimated quality are quantized first and
	 *             so on.
	 *
	 *             The function returns the number <i>t</i> of distinct
	 *             minutiae  quantizations that have been put into
	 *             <code>array</code> which depends on the size of the input
	 *             minutiae template
	 *             (see \link MinutiaeView::getMinutiaeCount()\endlink)
	 *             , <code>tmax</code> and on the count how many minutiae
	 *             had equal quantizations. Thereby, the function
	 *             guarantees that the first <i>t</i> elements of
	 *             <code>array</code> are distinct.
	 *
	 * @param array
	 *             Will contain at most \link getMaxGenuineFeatures()\endlink
	 *             finite field elements encoding the first minutiae'
	 *             quantization in <code>view</code>.
	 *
	 * @param view
	 *             The minutiae template of which the quantization set
	 *             is computed.
	 *
	 * @return
	 *             The number of minutiae quantizations extracted from
	 *             <code>view</code>.
	 *
	 * @warning
	 *             If <code>tmax</code> is not positive, an error message
	 *             is printed to <code>stderr</code> and the program exits
     *             with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <code>array</code> cannot hold
     *             \link getMaxGenuineFeatures()\endlink
	 *             elements of type <code>uint32_t</code>, calling this
	 *             function runs into undocumented behavior.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeTemplate::quantize
	( uint32_t *array , const MinutiaeView & view ) const {
		return quantize(array,view,getMaxGenuineFeatures());
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
	 *             template is attached with a public random permutation
	 *             process acting on the quantized feature universe.
	 *             The reordering of such a quantized feature (encoded as
	 *             a finite field element) can be evaluated via this
	 *             function.
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
	 *             'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If <i>a</i> does not encode a valid element in the
	 *             finite field that can be access via
	 *             \link getField()\endlink, an error message is printed
	 *             to <code>stderr</code> and the program exits with
     *             status 'EXIT_FAILURE'.
	 */
	uint32_t ProtectedMinutiaeTemplate::reorder( uint32_t a ) const {

		// Ensure the presence of a random permutation process
		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::reorder: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		// Call reordering without check.
		return _reorder(a);
	}

	/**
	 * @brief
	 *             Access the width of the fingerprint images related
	 *             to this protected minutiae template.
	 *
	 * @return
	 *             The width of the fingerprint images related to this
	 *             protected minutiae template.
	 *
	 * @see setDimension(int,int)
	 */
	int ProtectedMinutiaeTemplate::getWidth() const {
		return this->width;
	}

	/**
	 * @brief
	 *             Access the height of the fingerprint images related
	 *             to this protected minutiae template.
	 *
	 * @return
	 *             The height of the fingerprint images related to this
	 *             protected minutiae template.
	 *
	 * @see setDimension(int,int)
	 */
	int ProtectedMinutiaeTemplate::getHeight() const {
		return this->height;
	}

	/**
	 * @brief
	 *             Access the finger position of this protected
	 *             minutiae template.
	 *
	 * @return
	 *             The finger position of this protected minutiae
	 *             template.
	 *
	 * @see setFingerPosition(FINGER_POSITION_T)
	 */
	FINGER_POSITION_T ProtectedMinutiaeTemplate::getFingerPosition() const {
		return this->fingerPosition;
	}

	/**
	 * @brief
	 *             Access the value indicating the assumed resolution
	 *             in <i>dots per inch</i> at which the fingerprint
	 *             images related with this instance have been scanned.
	 *
	 * @details
	 *             The resolution in dots per inch.
	 *
	 * @see setResolution(int)
	 */
	int ProtectedMinutiaeTemplate::getResolution() const {
		return this->dpi;
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
	int ProtectedMinutiaeTemplate::getGridDist() const {
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
	int ProtectedMinutiaeTemplate::getNumAngleQuanta() const {
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
	int ProtectedMinutiaeTemplate::getSecretSize() const {
		return this->k;
	}

	/**
	 * @brief
	 *             Access the maximal number of minutiae quantizations
	 *             extracted from the to-be-protected minutiae template
	 *             on enrollment.
	 *
	 * @return
	 *             The maximal number of minutiae quantizations extracted
	 *             from the to-be-protected minutiae template on
	 *             enrollment.
	 *
	 * @see setMaxGenuineFeatures(int)
	 */
	int ProtectedMinutiaeTemplate::getMaxGenuineFeatures() const {
		return this->tmax;
	}

	/**
	 * @brief
	 *             Access the number of possible different minutiae
	 *             quantizations.
	 *
	 * @return
	 *             The number of possible different minutiae quantizations.
	 */
	int ProtectedMinutiaeTemplate::getVaultSize() const {
		return (int)(this->grid.size()) * this->s;
	}

	/**
	 * @brief
	 *             Access the number of hexagonal grid points centered in
	 *             a region in which absolutely pre-aligned minutiae can
	 *             occur.
	 *
	 * @return
	 *             The number of hexagonal gird points covering the
	 *             region in which absolutely pre-aligned minutiae
	 *             can occur.
	 */
	int ProtectedMinutiaeTemplate::getGridSize() const {
		return (int)(this->grid.size());
	}

	/**
	 * @brief
	 *             Access the number of decoding iterations performed
	 *             by the randomized decoder on verification attempts
	 *             with this protected minutiae template.
	 *
	 * @return
	 *             The number of decoding iterations performed by the
	 *             randomized decoder on verification with this protected
	 *             minutiae template.
	 *
	 * @see setNumberOfDecodingIterations(int)
	 */
	int ProtectedMinutiaeTemplate::getNumberOfDecodingIterations() const {
		return this->D;
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
	const vector< pair<double,double> > & ProtectedMinutiaeTemplate::getGrid() const {
		return this->grid;
	}

	/**
	 * @brief
	 *             Access the finite field in which computations with
	 *             the improved fuzzy vault scheme are performed to
	 *             protect a minutiae template.
	 *
	 * @returns
	 *             The \link SmallBinaryField finite field\endlink in which
	 *             computations with the improved fuzzy vault scheme are
	 *             performed to protect a minutiae template.
	 */
	const SmallBinaryField & ProtectedMinutiaeTemplate::getField() const {
		return this->gfPtr[0];
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
	int ProtectedMinutiaeTemplate::getVaultDataSize() const {

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
	 *             returns <code>false</code>.
	 *
	 * @see getVaultDataSize()
	 * @see setSlowDownFactor()
	 * @see getSlowDownFactor()
	 * @see encrypt()
	 * @see decrypt()
	 */
	const uint8_t * ProtectedMinutiaeTemplate::getVaultData() const {

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
	SmallBinaryFieldPolynomial ProtectedMinutiaeTemplate::unpackVaultPolynomial
	( const BigInteger & slowDownValue ) const {


		if ( !isDecrypted() ) {
			cerr << "ProtectedMinutiaeTemplate::unpackVaultPolynomial: "
				 << "no decrypted vault data." << endl;
			exit(EXIT_FAILURE);
		}

		if ( slowDownValue.sign() < 0 ||
			 BigInteger::compare(slowDownValue,getSlowDownFactor()) >= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::unpackVaultPolynomial: "
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
			cerr << "ProtectedMinutiaeTemplate::createVaultPolynomialCandidate: "
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
	const BigInteger & ProtectedMinutiaeTemplate::getSlowDownFactor() const {

		return this->slowDownFactor;
	}

	/**
	 * @brief
	 *             Access the SHA-1 hash value of the secret key
	 *             generated on enrollment used to obfuscate the
	 *             correct minutiae template's quantization.
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
	const uint8_t *ProtectedMinutiaeTemplate::getHash() const {

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
	 *             Change the dimension of the fingerprint images
	 *             related with this protected minutiae template
	 *
	 * @details
	 *             Changing the dimension can cause to update the
	 *             \link getGrid() grid\endlink and the
	 *             \link getField() finite field\endlink related with this
	 *             protected minutiae template.
	 *
	 * @param width
	 *             The width of the fingerprint image.
	 *
	 * @param height
	 *             The height of the fingerprint image.
	 *
	 * @warning
	 *             If \link isEnrolled()\endlink returns
	 *             <code>true</code>, an error message is printed to
	 *             <code>stderr</code> and the program exits with status
	 *             -1.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::setDimension( int width , int height ) {

		// Ensure that this instance is empty
		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setDimension: "
				 << "clear first." << endl;
            exit(EXIT_FAILURE);
		}

		// A change of the dimension may require an update of the grid
		// and of the underlying finite field
		updateGrid();
		updateField();
	}

	/**
	 * @brief
	 *             Change the finger position of this protected
	 *             minutiae template.
	 *
	 * @param fingerPosition
	 *             The new finger position of this protected minutiae
	 *             template.
	 *
	 * @see getFingerPosition()
	 */
	void ProtectedMinutiaeTemplate::setFingerPosition( FINGER_POSITION_T fingerPosition ) {
		this->fingerPosition = fingerPosition;
	}

	/**
	 * @brief
	 *             Change the value indicating the assumed resolution
	 *             in <i>dots per inch</i> at which the fingerprint
	 *             images related with this instance have been scanned.
	 *
	 * @details
	 *             The new resolution in dots per inch.
	 *
	 * @see getResolution(int)
	 *
	 * @warning
	 *             If <code>dpi</code> is smaller than 300 or greater
	 *             than 1000, an error message is printed to
	 *             <code>stderr</code> and the program exits with
     *             status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::setResolution( int dpi ) {

		if ( dpi < 300 || dpi > 1000 ) {
			cerr << "ProtectedMinutiaeTemplate::setResolution: resolution must vary between 300 and 1000 dots per inch." << endl;
            exit(EXIT_FAILURE);
		}

		this->dpi = dpi;
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
	void ProtectedMinutiaeTemplate::setGridDist( int gridDist ) {

		// Ensure that this instance is empty
		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setGridDist: "
				 << "clear first." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure reasonable bounds of the hexagonal grid's distance
		if ( gridDist <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::setGridDist: "
				 << "must be greater than zero." << endl;
            exit(EXIT_FAILURE);
		}

		this->gridDist = gridDist;

		// Update the grid which may require to update the underlying
		// finite field as well
		updateGrid();
		updateField();
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
	void ProtectedMinutiaeTemplate::setNumAngleQuanta( int s ) {

		// Ensure that this instance is empty
		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setNumAngleQuanta: "
				 << "clear first." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure reasonable bounds of quanta size
		if ( s <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::setNumAngleQuanta: "
				 << "must be greater than zero." << endl;
            exit(EXIT_FAILURE);
		}

		this->s = s;

		// The finite field may need an update
		updateField();
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
	void ProtectedMinutiaeTemplate::setSecretSize( int k ) {

		// Ensure that this instance is empty
		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setSecretSize: "
				 << "already enrolled; clear the view first." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure reasonable bounds
		if ( k <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::setSecretSize: "
				 << "Must be greater than zero." << endl;
            exit(EXIT_FAILURE);
		}

		this->k = k;
	}

	/**
	 * @brief
	 *             Change the maximal number of minutiae quantizations
	 *             extracted from the to-be-protected minutiae template
	 *             on enrollment.
	 *
	 * @param tmax
	 *             The new maximal number of minutiae quantizations extracted
	 *             from the to-be-protected minutiae template on
	 *             enrollment.
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
	void ProtectedMinutiaeTemplate::setMaxGenuineFeatures( int tmax ) {

		// Ensure that this instance is empty
		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setMaxGenuineFeatures: "
				 << "clear the view first." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure reasonable bounds
		if ( tmax <= 0 ) {
			cerr << "ProtectedMinutiaeTempalte::setMaxGenuineFeatures: "
				 << "must be greater than zero." << endl;
            exit(EXIT_FAILURE);
		}

		this->tmax = tmax;
	}

	/**
	 * @brief
	 *             Change the number of decoding iterations performed
	 *             by the randomized decoder on verification attempts
	 *             with this protected minutiae template.
	 *
	 * @return
	 *             The new number of decoding iterations that are performed
	 *             by the randomized decoder on verification with this
	 *             protected minutiae template.
	 *
	 * @see getNumberOfDecodingIterations()
	 *
	 * @warning
	 *             If <i>D</i> is smaller than or equals 0, the method
	 *             prints an error message to <code>stderr</code> and
     *             exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::setNumberOfDecodingIterations( int D ) {

		// Ensure reasonable bounds
		if ( D <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::setNumberOfDecodingIterations: "
				 << "must be greater than 0." << endl;
            exit(EXIT_FAILURE);
		}

		this->D = D;
	}

	/**
	 * @brief
	 *            Specifies the slow-down factor manually to artifically
	 *            increase the time for an impostor verification attempt.
	 *
	 * @details
	 *            The time needed for an impostor to perform
	 *            off-line attacks is multiplied by the slow-down
	 *            factor, thereby increasing security while also
	 *            increasing response time of genuine verification
	 *            attempts---a trade-off thus but with the capability
	 *            of maintaining a relative security measure with
	 *            the advent of more and more powerful computers.
	 *
	 * @param slowDownFactor
	 *            The specified slow-down factor for this protected
	 *            minutiae template.
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
	void ProtectedMinutiaeTemplate::setSlowDownFactor( const BigInteger & slowDownFactor ) {

		if ( slowDownFactor.sign() <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate::setSlowDownFactor: argument "
				 << "must not be zero or negative." << endl;
			exit(EXIT_FAILURE);
		}

		if ( isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::setSlowDownFactor: "
				 << "parameter must be specified before enrolment." << endl;
			exit(EXIT_FAILURE);
		}

		this->slowDownFactor = slowDownFactor;
	}

	/**
	 * @brief
	 *             Clears all data that is related with a minutiae
	 *             template protected by this instance (if any).
	 *
	 * @details
	 *             Any \link ProtectedMinutiaeTemplate\endlink is empty
	 *             (cleared) by default until a
	 *             \link MinutiaeView minutiae template\endlink is passed
	 *             through the \link enroll()\endlink function returning
	 *             <code>true</code>. If the user wants to reuse this
	 *             \link ProtectedMinutiaeTemplate\endlink to protect another
	 *             \link MinutiaeView minutiae template\endlink via
	 *             \link enroll()\endlink, he must clear it first via
	 *             this method.
	 */
	void ProtectedMinutiaeTemplate::clear() {

		// Free any data holding the (encrypted) vault.
		free(this->vaultPolynomialData);
		free(this->encryptedVaultPolynomialData);

		// Set NULL pointers.
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;

		// zero 160 bit hash value
		memset(this->hash,0,20);

		// There is no count for the non-leading coefficients
		// in the vault polynomial thus set to -1.
		this->t = -1;
	}

	/**
	 * @brief
	 *             Determine the size in bytes needed to store this
	 *             protected minutiae template.
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
	int ProtectedMinutiaeTemplate::getSizeInBytes() const {

		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::getSizeInBytes: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		int size = 0;

		size += 10; // size for the header
		size += 2; // size to encode 'width'
		size += 2; // size to encode 'height'
		size += 1; // size to encode 'fingerPosition'
		size += 2; // size to encode 'dpi'
		size += 1; // size to encode 'gridDist'
		size += 1; // size to encode 's'
		size += 1; // size to encode 'k'
		size += 1; // size to encode 'tmax'
		size += 4; // size to encode 'D'
		size += 4; // size to encode generator polynomial of finite field
		size += 1; // size to encode 't'
		size += 1; // size to encode whether vault is encrypted or not
		size += 4+this->slowDownFactor.getSizeInBytes(); // size to encode the slow-down factor
		size += vaultDataSize(); // size to encoded vault data.
		size += 20; // size to encode hash value of secret polynomial

		return size;
	}

	/**
	 * @brief
	 *             Store data defining this protected minutiae template to
	 *             the specified byte array.
	 *
	 * @param data
	 *             The byte array to which the data for this protected
	 *             minutiae template is written.
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
	int ProtectedMinutiaeTemplate::toBytes( uint8_t *data ) const {

		if ( !isEnrolled() ) {
			cerr << "ProtectedMinutiaeTemplate::toBytes: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		int offset = 0;

		uint32_t tmp;

		memcpy(data,"PMT140818",10);

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

		// Write 'fingerPosition'
		data[offset+0] = (uint8_t)(this->fingerPosition);
		offset += 1;

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

		// Write 'D'
		tmp = (uint32_t)(this->D);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write finite field generator polynomial
		tmp = (uint32_t)(this->gfPtr->getDefiningPolynomial().rep);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write 't'
		data[offset+0] = (uint8_t)(this->t);
		offset += 1;

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

		// Write seed of secret polynomial
		memcpy(data+offset,this->hash,20);
		offset += 20;

		return offset;
	}

	/**
	 * @brief
	 *             Initialize this protected minutiae template by the
	 *             specified bytes.
	 *
	 * @details
	 *             If this protected minutiae template has been successfully
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
	int ProtectedMinutiaeTemplate::fromBytes( const uint8_t *data , int size ) {

		int offset = 0;

		if ( size < offset+10 ) {
			return -1;
		}

		ProtectedMinutiaeTemplate tmp;

		{ // Check if header is valid
			if ( data[0] != 'P' || data[1] != 'M' || data[2] != 'T' ||
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

		// Check and read 'fingerPosition'
		if ( size < offset+1 ) { return -1; }
		if ( (int)data[offset] > 10 ) { return -1; }
		tmp.fingerPosition = (FINGER_POSITION_T)data[offset];
		offset += 1;

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

		// Read 'D'
		if ( size < offset+4 ) { return -1; }
		tmp.D  = (int)data[offset+0]; tmp.D <<= 8;
		tmp.D += (int)data[offset+1]; tmp.D <<= 8;
		tmp.D += (int)data[offset+2]; tmp.D <<= 8;
		tmp.D += (int)data[offset+3];
		offset += 4;

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
				cerr << "ProtectedMinutiaeTemplate::fromBytes: out of memory." << endl;
                exit(EXIT_FAILURE);
			}

			if ( (int)(tmp.gfPtr->getCardinality()) < this->s * (int)tmp.grid.size() ) {
				return -1;
			}
		}

		// Read 't'
		if ( size < offset+1 ) { return -1; }
		tmp.t = (int)data[offset];
		offset += 1;

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

		// Read 'hash'
		if ( size < offset+20 ) { return -1; }
		memcpy(tmp.hash,data+offset,20);
		offset += 20;

		tmp.updatePermutation();

		swap(*this,tmp);

		return offset;
	}

	/**
	 * @brief
	 *             Write data defining this protected minutiae template to
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
	int ProtectedMinutiaeTemplate::write( FILE *out ) const {

		uint8_t *data;
		int size , wsize;

		// Initialize byte array ...
		size = getSizeInBytes();
		data = (uint8_t*)malloc( size * sizeof(uint8_t) );
		if ( data == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::write: "
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
	 *             Write data defining this protected minutiae template to
	 *             the specified file.
	 *
	 * @param file
	 *             Path to the file to which this protected minutiae
	 *             template is written.
	 *
	 * @return
	 *             The function returns <code>true</code> if this protected
	 *             minutiae template has been successfully be written
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
	bool ProtectedMinutiaeTemplate::write( const std::string & file ) const {

		FILE *out;
		bool success;

		// Open a FILE pointer to the specified path ...
		out = THIMBLE_FOPEN(file.c_str(),"wb");
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
	 *             Initialize this protected minutiae template from the
	 *             specified file pointer.
	 *
	 * @details
	 *             If this protected minutiae template has been successfully
	 *             initialized, the function returns the number read bytes;
	 *             otherwise, if an error occurred, this protected minutiae
	 *             template will be left unchanged and the result of the
	 *             function will be -1.
	 *
	 * @param in
	 *             The file pointer from which the protected minutiae
	 *             template is read.
	 *
	 * @return
	 *             The number of bytes that were successfully read from
	 *             <code>in</code> to initialize this protected minutiae
	 *             template; otherwise, if an error occurred,
	 *             the function returns -1.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	int ProtectedMinutiaeTemplate::read( FILE *in ) {

		int offset = 0;

		int c;
		ProtectedMinutiaeTemplate tmp;

		{ // Check if header is valid
			uint8_t header[10];
			for ( int j = 0 ; j < 10 ; j++ ) {
				if ( (c=fgetc(in)) < 0 ) { return -1; }
				header[j] = (uint8_t)c;
			}
			if ( header[0] != 'P' || header[1] != 'M' || header[2] != 'T' ||
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

		// Check and read 'fingerPosition'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		if ( c > 10 ) { return -1; }
		tmp.fingerPosition = (FINGER_POSITION_T)c;
		offset += 1;

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

		// Read 'D'
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.D  = c; tmp.D <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.D += c; tmp.D <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.D += c; tmp.D <<= 8;
		if ( (c=fgetc(in)) < 0 ) { return -1; }
		tmp.D += c;
		offset += 4;

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
		tmp.t = c;
		offset += 1;

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


		// *******************************************************************
		// ********************** BEGIN: Read data ***************************
		// *******************************************************************

		int n = tmp.vaultDataSize();

		uint8_t *vaultData = (uint8_t*)malloc( n * sizeof(uint8_t) );
		if ( vaultData == NULL ) {
			cerr << "ProtectedMinutiaeTemplate::read: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		if ( (c=fgetc(in)) < 0 ) { free(vaultData); return -1; }
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
	 *             Initialize this protected minutiae template from the
	 *             data contained in the specified file path.
	 *
	 * @details
	 *             If this protected minutiae template has been successfully
	 *             initialized, the function returns <code>true</code>;
	 *             otherwise, if an error occurred, this protected minutiae
	 *             template will be left unchanged and the result of the
	 *             function will be <code>false</code>.
	 *
	 * @param file
	 *             String representation of the path to the file from
	 *             which the protected minutiae template is read.
	 *
	 * @return
	 *             <code>true</code> if this protected minutiae template
	 *             has been successfully initialized by the data contained
	 *             in the specified file; otherwise, the result is
	 *             <code>false</code>.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
	 */
	bool ProtectedMinutiaeTemplate::read( const std::string & file ) {

		bool success;

		// Open a FILE pointer to the specified path ...
		FILE *in = THIMBLE_FOPEN(file.c_str(),"rb");
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
	 *            Convenience method to determine how many elements
	 *            two arrays share.
	 *
	 * @details
	 *            see 'ProtectedMinutiaeTemplate.h'
	 */
	int ProtectedMinutiaeTemplate::overlap
	( const uint32_t *array1 , int t1 , const uint32_t *array2 , int t2 ) {

		// Counts the number of common elements
		int omega = 0;

		// Iterate over all pairs
		for ( int i = 0 ; i < t1 ; i++ ) {
			for ( int j = 0 ; j < t2 ; j++ ) {

				// Increment in case the two elements are equal
				if ( array1[i] == array2[j] ) {
					++omega;
				}
			}
		}

		return omega;
	}

	/**
	 * @brief
	 *            Initializes the member variables of this protected
	 *            minutiae template for the very first time such that
	 *            it becomes ready for \link enroll() enrollment\endlink.
	 *
	 * @details
	 *            The method assumes that the member variables
	 *            \link width\endlink, \link height\endlink,
	 *            \link gridDist\endlink, \link s\endlink have been
	 *            properly initialized already; otherwise, the program
	 *            runs into undocumented behavior.
	 *
	 *            The methods should not be called directly as it is
	 *            called by other initialization methods of this class.
	 *
	 *            The method initializes the grid and the finite field
	 *            using the \link updateGrid()\endlink and
	 *            \link updateField()\endlink functions, respectively.
	 */
	void ProtectedMinutiaeTemplate::first_init() {

		this->t = -1;
		this->gfPtr                        = NULL;
		this->vaultPolynomialData          = NULL;
		this->slowDownFactor               = BigInteger(1);
		this->encryptedVaultPolynomialData = NULL;

		memset(this->hash,0,20);

		updateGrid();
		updateField();

		this->is_initialized = true;
	}

	/**
	 * @brief
	 *             Initializes this object such that it can protected
	 *             absolutely pre-aligned minutiae templates.
	 *
	 * @details
	 *             This method generates parameters needed to properly
	 *             quantize absolutely pre-aligned minutiae templates that
	 *             have been estimated from a fingerprint image of the
	 *             specified dimension and that has been scanned at the
	 *             specified resolution. If the finger's position is known,
	 *             it can be specified as well.
	 *
	 *             This method should be called directly only if no memory
	 *             is reserved by the member of this object such that
	 *             allocation on the attributes does not cause information
	 *             leakage.
	 *
	 * @param width
	 *             The width of the fingerprint image.
	 *
	 * @param height
	 *             The height of the fingerprint image.
	 *
	 * @param dpi
	 *             The resolution in (dots per inch) at which
	 *             the fingerprint image has been scanned.
	 *
	 * @param fingerPosition
	 *             The position of the finger (if known).
	 *
	 * @attention
	 *             If <code>width</code> or <code>height</code> are
	 *             not greater than 0 or if <code>dpi</code> does
	 *             not range between 300 and 1000, an error message
	 *             is printed to <code>stderr</code> and the program
	 *             exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *             If not enough memory could be provided, an error
	 *             message is printed to <code>stderr</code> and the
	 *             program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::first_init
	( int width , int height , int dpi , FINGER_POSITION_T fingerPosition ) {

		// Check whether dimension is okay.
		if ( width <= 0 || height <= 0 ) {
			cerr << "ProtectedMinutiaeTemplate: Invalid image dimensions."
				 << endl;
            exit(EXIT_FAILURE);
		}

		// Check whether resolution is okay
		if ( dpi < 300 || dpi > 1000 ) {
			cerr << "ProtectedMinutiaeTemplate: "
				 << "Resolution must vary between 300 and 1000 dpi." << endl;
            exit(EXIT_FAILURE);
		}

		// Set members specified by the user of the class.
		this->width  = width;
		this->height = height;
		this->dpi = dpi;
		this->fingerPosition = fingerPosition;


		this->gridDist = (int)THIMBLE_ROUND(29.0 / 569.0 * (double)dpi);
		this->s        = 6;
		this->k        = 10;
		this->tmax     = 44;

		// Standard value for the number of decoding iterations performed
		// by the randomized decoder on a verification attempt.
		this->D = 1<<16;

		// Finish initialization
		first_init();
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
	 *            needs an update which is performed when calling
	 *            this method.
	 */
	void ProtectedMinutiaeTemplate::updateGrid() {

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
	 *            Update the finite field such that it can encode
	 *            each possible minutia's quantization.
	 *
	 * @details
	 *            After the \link grid hexagonal grid points\endlink
	 *            or the
	 *            \link s number of minutiae angle quantizations\endlink
	 *            have changed, the finite field might need an update
	 *            which can be performed by calling this method.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
	 */
	void ProtectedMinutiaeTemplate::updateField() {

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

			// Free space used by a possibly non-NULL pointer to a finite field
			delete this->gfPtr;
			this->gfPtr = NULL;

			// Create new binary field of the desired degree
			this->gfPtr = new(nothrow) SmallBinaryField(degree);

			// Check whether sufficient memory could be provided and exit otherwise
			if ( this->gfPtr == NULL ) {
				cerr << "ProtectedMinutiaeTemplate::updateField: out of memory."
					 << endl;
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
	void ProtectedMinutiaeTemplate::updatePermutation() {

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
	uint32_t ProtectedMinutiaeTemplate::_reorder( uint32_t a ) const {
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
	void ProtectedMinutiaeTemplate::packVaultPolynomial
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
			cerr << "ProtectedMinutiaeTemplate::packVaultPolynomial: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		for ( int j = 0 ; j < t ; j++ ) {
			coeffs[j] = V.getCoeff(j);
		}

		concat_bit_vectors(this->vaultPolynomialData,coeffs,t,d,n);

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
	int ProtectedMinutiaeTemplate::vaultDataSize() const {

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
	AES128 ProtectedMinutiaeTemplate::deriveKey( const BigInteger & x ) {

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
	void ProtectedMinutiaeTemplate::concat_bit_vectors
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
	void ProtectedMinutiaeTemplate::split_into_bit_vectors
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
