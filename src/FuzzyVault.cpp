/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
 *
 *  Copyright 2014, 2015 Benjamin Tams
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
 * @file FuzzyVault.cpp
 *
 * @brief
 *            Implementation of a class representing general instances
 *            of an implementation of the improved fuzzy vault scheme by
 *            Dodis et al..
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstring>
#include <iostream>

#include <thimble/misc/CTools.h>
#include <thimble/misc/IOTools.h>
#include <thimble/math/MathTools.h>
#include <thimble/math/RandomGenerator.h>
#include <thimble/math/Permutation.h>
#include <thimble/math/numbertheory/BigInteger.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/ecc/ReedSolomonCode.h>
#include <thimble/ecc/GuruswamiSudanDecoder.h>
#include <thimble/security/AES.h>
#include <thimble/security/SHA.h>
#include <thimble/security/FuzzyVaultTools.h>
#include <thimble/security/FuzzyVault.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	FuzzyVault::FuzzyVault() {
		first_init();
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	FuzzyVault::FuzzyVault( int n , int tmax , int k ) {
		first_init();

		if ( n <= 0 || tmax <= 0 || k <= 0 ||
			 n < tmax || tmax < k ) {
			cerr << "FuzzyVault: invalid arguments." << endl;
			exit(EXIT_FAILURE);
		}

		setVaultSize(n);
		setMaxGenuineFeatures(tmax);
		setSecretSize(k);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	FuzzyVault::FuzzyVault( const FuzzyVault & vault ) {
		first_init();
		assign(vault);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	FuzzyVault::~FuzzyVault() {
		delete this->gfPtr;
		clear();
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::assign( const FuzzyVault & vault ) {

		if ( this != &vault ) {

			// Copy primitive member variables
			this->n = vault.n;
			this->tmax = vault.tmax;
			this->k = vault.k;

			// Copy member for which the =-operator has been
			// implemented.
			this->slowDownFactor = vault.slowDownFactor;
			this->permutation = vault.permutation;

			// Copy underlying finite field
			if ( this->gfPtr == NULL && vault.gfPtr != NULL ) {
				this->gfPtr = new (nothrow) SmallBinaryField(vault.gfPtr[0]);
				if ( this->gfPtr == NULL ) {
					cerr << "FuzzyVault: "
						 << "out of memory." << endl;
                    exit(EXIT_FAILURE);
				}
			} else if ( this->gfPtr != NULL && vault.gfPtr != NULL ){
				this->gfPtr[0] = vault.gfPtr[0];
			} else if ( this->gfPtr != NULL && vault.gfPtr == NULL ) {
				delete this->gfPtr;
				this->gfPtr = NULL;
			}

			// Copy vault data
			if ( this->vaultPolynomialData == NULL &&
				vault.vaultPolynomialData != NULL ) {
				int n = vaultDataSize();
				this->vaultPolynomialData =
						(uint8_t*)malloc( n * sizeof(uint8_t) );
				if ( this->vaultPolynomialData == NULL ) {
					cerr << "FuzzyVault: out of memory." << endl;
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
					cerr << "FuzzyVault: out of memory." << endl;
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
					cerr << "FuzzyVault: out of memory." << endl;
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
					cerr << "FuzzyVault: out of memory." << endl;
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
		}

	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	FuzzyVault & FuzzyVault::operator=( const FuzzyVault & vault ) {
		assign(vault);
		return *this;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::swap( FuzzyVault & vault ) {

		// Backup primitive members of this vault object
		int n = this->n;
		int tmax = this->tmax;
		int k = this->k;
		SmallBinaryField *gfPtr = this->gfPtr;
		uint8_t *vaultPolynomialData = this->vaultPolynomialData;
		uint8_t *encryptedVaultPolynomialData = this->encryptedVaultPolynomialData;

		// Copy primitive members of 'vault' to 'this'
		this->n = vault.n;
		this->tmax = vault.tmax;
		this->k = vault.k;
		this->gfPtr = vault.gfPtr;
		this->vaultPolynomialData = vault.vaultPolynomialData;
		this->encryptedVaultPolynomialData = vault.encryptedVaultPolynomialData;

		// Copy backups to 'vault'
		vault.n = n;
		vault.tmax = tmax;
		vault.k = k;
		vault.gfPtr = gfPtr;
		vault.vaultPolynomialData = vaultPolynomialData;
		vault.encryptedVaultPolynomialData = encryptedVaultPolynomialData;

		// SPECIAL CASE: Swap the 'slowDownFactor' fields using
		// the 'BigInteger::swap' method.
        this->slowDownFactor.swap(vault.slowDownFactor);

		// SPECIAL CASE: swap the secret polynomials' hash values
		// via 'memcpy'
		uint8_t hash[20];
		memcpy(hash,this->hash,20);
		memcpy(this->hash,vault.hash,20);
		memcpy(vault.hash,hash,20);

		// SPECIAL CASE: Swap the 'permutation' fields using
		// the 'Permutation::swap' method.
        Permutation::swap(this->permutation,vault.permutation);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::setVaultSize( int n ) {

		if ( isEnrolled() ) {
			cerr << "FuzzyVault::setVaultSize: "
				 << "already protected data; clear first." << endl;
			exit(EXIT_FAILURE);
		}

		if ( n <= 0 ) {
			cerr << "FuzzyVault::setVaultSize: "
				 << "invalid argument." << endl;
			exit(EXIT_FAILURE);
		}

		this->n = n;

		updateField();
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::setMaxGenuineFeatures( int tmax ) {

		if ( isEnrolled() ) {
			cerr << "FuzzyVault::setMaxGenuineFeatures: "
				 << "already protected data; clear first." << endl;
			exit(EXIT_FAILURE);
		}

		if ( tmax <= 0 ) {
			cerr << "FuzzyVault::setMaxGenuineFeatures: "
				 << "invalid argument." << endl;
			exit(EXIT_FAILURE);
		}

		this->tmax = tmax;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::setSecretSize( int k ) {

		if ( isEnrolled() ) {
			cerr << "FuzzyVault::setSecretSize: "
				 << "already protected data; clear first." << endl;
			exit(EXIT_FAILURE);
		}

		if ( k <= 0 ) {
			cerr << "FuzzyVault::setSecretSize: "
				 << "invalid argument." << endl;
			exit(EXIT_FAILURE);
		}

		this->k = k;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::setSlowDownFactor( const BigInteger & slowDownFactor ) {

		if ( isEnrolled() ) {
			cerr << "FuzzyVault::setSlowDownFactor: "
				 << "already protected data; clear first." << endl;
			exit(EXIT_FAILURE);
		}

		if ( slowDownFactor.sign() <= 0 ) {
			cerr << "FuzzyVault::setSlowDownFactor: "
				 << "invalid argument." << endl;
			exit(EXIT_FAILURE);
		}

		this->slowDownFactor = slowDownFactor;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getVaultSize() const {

		return this->n;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getMaxGenuineFeatures() const {

		return this->tmax;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getSecretSize() const {

		return this->k;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	const BigInteger & FuzzyVault::getSlowDownFactor() const {
		return this->slowDownFactor;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::enroll( const uint32_t *features , int t ) {

		if ( !isInitialized() ) {
			cerr << "FuzzyVault::enroll: "
				 << "parameters not specified; initialize first." << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->n < this->tmax || this->tmax < this->k ) {
			cerr << "FuzzyVault::enroll: "
				 << "invalid parameter configuration." << endl;
			exit(EXIT_FAILURE);
		}

		// Causes the program to exit if feature set is invalid.
		checkFeatureSet(features,t);

		if ( t < this->k ) {
			return false;
		}

		if ( t > this->tmax ) {
			t = this->tmax;
		}

		SmallBinaryFieldPolynomial f(getField());
		f.random(this->k,true);

        SHA().hash(this->hash,f.getData(),f.deg()+1);
        updatePermutation();

		SmallBinaryFieldPolynomial V(getField());

		// Allocate memory.
		uint32_t *A = (uint32_t*)malloc( this->tmax * sizeof(uint32_t) );

		// Shuffle the feature elements using the record-specific feature
		// permutation process
		for ( int j = 0 ; j < t ; j++ ) {
			A[j] = _reorder(features[j]);
		}

		// Supplement the feature set with blending features
		for ( int j = t ; j < this->tmax ; j++ ) {
			A[j] = this->n + MathTools::rand32(true) %
					(getField().getCardinality() - this->n );
		}

		V.buildFromRoots(A,this->tmax);
		V += f;

		packVaultPolynomial(V);

		// Free any allocated memory.
		free(A);

		return true;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::enrol( const uint32_t *features , int t ) {

		return enroll(features,t);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::verify( const uint32_t *queryFeatures , int s ) const {

		SmallBinaryFieldPolynomial f(getField());

		bool state = open(f,NULL,NULL,queryFeatures,s);

		return state;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::encrypt( const AES128 & key ) {

		// Ensure that this instance does protect a minutiae template.
		if ( !isEnrolled() ) {
			cerr << "FuzzyVault::encrypt: no vault built." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure that this instance does not already hold encrypted data
		if ( isEncrypted() ) {
			cerr << "FuzzyVault::encrypt: already encrypted."
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

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::decrypt( const AES128 & key ) {

		// Ensure that this instance does contain encrypted data.
		if ( this->encryptedVaultPolynomialData == NULL ) {
			cerr << "FuzzyVault::decrypt: "
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
				cerr << "FuzzyVault::decrypt: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}

		aes.decrypt
		(this->vaultPolynomialData,this->encryptedVaultPolynomialData,n);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::isInitialized() const {
		return (this->n > 0 && this->tmax > 0 && this->k > 0);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::isEnrolled() const {

		return this->vaultPolynomialData != NULL ||
			   this->encryptedVaultPolynomialData != NULL;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::isEnroled() const {

		return isEnrolled();
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::isEncrypted() const {

		return this->encryptedVaultPolynomialData != NULL;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::isDecrypted() const {

		return this->vaultPolynomialData != NULL;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::clear() {
		free(this->vaultPolynomialData);
		free(this->encryptedVaultPolynomialData);
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	const SmallBinaryField & FuzzyVault::getField() const {

		if ( this->gfPtr == NULL ) {
			cerr << "FuzzyVault::getField: "
				 << "not initialized." << endl;
			exit(EXIT_FAILURE);
		}

		return this->gfPtr[0];
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	const uint8_t * FuzzyVault::getHash() const {
		return this->hash;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::open
	( SmallBinaryFieldPolynomial & f ,
	  uint32_t *features , int *featureSizePtr ,
	  const uint32_t *queryFeatures , int s ) const {

		if ( !isEnrolled() ) {
			cerr << "FuzzyVault::verify: "
				 << "not enrolled." << endl;
			exit(EXIT_FAILURE);
		}

		if ( !isDecrypted() ) {
			cerr << "FuzzyVault::verify: "
				 << "does not contain decrypted vault data." << endl;
			exit(EXIT_FAILURE);
		}

		checkFeatureSet(queryFeatures,s);

		if ( s < getSecretSize() ) {
			return false;
		}

		bool success = false;

		SmallBinaryFieldPolynomial V(getField());
		uint32_t *x , *y;
		x = (uint32_t*)malloc( s * sizeof(uint32_t) );
		y = (uint32_t*)malloc( s * sizeof(uint32_t) );
		if ( x == NULL || y == NULL ) {
			cerr << "FuzzyVault::open: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Iteration over possible slow-down values.
		for ( BigInteger slowDownVal = 0 ;
			  BigInteger::compare(slowDownVal,this->slowDownFactor) < 0 ;
			  add(slowDownVal,slowDownVal,1) ) {

			// Unpack the vault using the current slow-down value
			// as the decryption key.
			V = unpackVaultPolynomial(slowDownVal);

			// Build unlocking set '{ (x[j],y[j]) }'
			for ( int j = 0 ; j < s ; j++ ) {
				// Do not forget to apply the application to the
				// query feature set.
				x[j] = _reorder(queryFeatures[j]);
				y[j] = V.eval(x[j]);
			}

			// Decoding attempt.
			success = decode(f,x,y,s,getSecretSize(),getHash());

			// We are done if the decoding attempt was successful
			if ( success ) {
				break;
			}
		}

		// If the decoding attempt was successful and if
		// the 'features' is non-NULL, ...
		if ( success && features != NULL ) {

			// ..., then recover the entire feature set and store
			// the result in 'features'.

			Permutation inversePermutation;
			inv(inversePermutation,this->permutation);

			// The feature size can be at most 'tmax'
			int featureSize = this->tmax;

			// The roots of 'V-f' encode the feature elements
			(V-f).findRoots(features);

			for ( int i = 0 ; i < featureSize ; i++ ) {

				// A root of 'V-f' can only be a valid feature element
				// if its decimal value is smaller than 'n'; ...
				if ( features[i] < (uint32_t)(this->n) ) {
					features[i] = (uint32_t)inversePermutation.eval
							((int)features[i]);
				} else {
					// ...; otherwise, the root may be a blending
					// element; in this case, the value should be
					// removed from '(features,featureSize)'
					for ( int j = i+1 ; j < featureSize ; j++ ) {
						features[j-1] = features[j];
					}
					--featureSize;
					--i;
				}
			}

			CTools::sort(features,featureSize);

			if ( featureSizePtr != NULL ) {
				*featureSizePtr = featureSize;
			}
		}

		free(x);
		free(y);

		return success;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	SmallBinaryFieldPolynomial FuzzyVault::unpackVaultPolynomial
		( const BigInteger & slowDownValue ) const {

		if ( !isDecrypted() ) {
			cerr << "FuzzyVault::unpackVaultPolynomial: "
				 << "no decrypted vault data." << endl;
			exit(EXIT_FAILURE);
		}

		if ( slowDownValue.sign() < 0 ||
			 BigInteger::compare(slowDownValue,getSlowDownFactor()) >= 0 ) {
			cerr << "FuzzyVault::unpackVaultPolynomial: "
				 << "slow-down value must be positive and smaller than the "
				 << "slow-down factor" << endl;
			exit(EXIT_FAILURE);
		}

		AES128 aes = deriveKey(slowDownValue);

		int t , d , n;
		t = this->tmax;
		d = this->gfPtr->getDegree();
		n = vaultDataSize();

		uint8_t *data = (uint8_t*)malloc( n * sizeof(uint8_t) );
		uint32_t *coeffs = (uint32_t*)malloc( t * sizeof(uint32_t ) );
		if ( data == NULL || coeffs == NULL ) {
			cerr << "FuzzyVault::createVaultPolynomialCandidate: "
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

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	uint32_t FuzzyVault::reorder( uint32_t a ) const {

		// Ensure the presence of a random permutation process
		if ( !isEnrolled() ) {
			cerr << "FuzzyVault::reorder: "
				 << "not enrolled." << endl;
            exit(EXIT_FAILURE);
		}

		// Ensure that 'a' encodes a valid feature element from
		// the feature universe.
		if ( a >= (uint32_t)this->n ) {
			cerr << "FuzzyVault::reorder: "
				 << "feature element is too large." << endl;
			exit(EXIT_FAILURE);
		}

		// Call reordering without check.
		return _reorder(a);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getGuruswamiSudanMultiplicity() const {

		return 1;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getNumDecIts() const {

		return 0;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::decode
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y , int u , int k ,
	  const uint8_t hash[20] ) const {

		if ( k > u ) {
			return false;
		}

		int m = getGuruswamiSudanMultiplicity();

		SHA sha;
		uint8_t _hash[20];
		SmallBinaryFieldPolynomial _f(f.getField());

		if ( m >= 0 ) {
			if ( ReedSolomonCode::decode(_f,x,y,u,k) ) {
				sha.hash(_hash,_f.getData(),_f.deg()+1);
				if ( memcmp(hash,_hash,20) == 0 ) {
					f = _f;
					return true;
				}
			}
		}

		for ( int m0 = 1 ; m0 <= m ; m0++ ) {

			GuruswamiSudanDecoder dec;

			dec.decode(x,y,u,k,m0,f.getField());

			for ( int j = 0 ; j < (int)dec.getDecodedList().size() ; j++ ) {

				_f = dec.getDecodedList().at(j);

				sha.hash(_hash,_f.getData(),_f.deg()+1);
				if ( memcmp(hash,_hash,20) == 0 ) {
					f = _f;
					return true;
				}
			}
		}

		int D = getNumDecIts();

		if ( D > 0 ) {
			return FuzzyVaultTools::bfattack(f,x,y,u,k,hash,D);
		}

		return false;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::getSizeInBytes() const {

		if ( !isEnrolled() ) {
			cerr << "FuzzyVault::getSizeInBytes: "
				 << "does not contained protected data." << endl;
			exit(EXIT_FAILURE);
		}

		int size = 0;

		size += 4; // size for the header string 'FVR'.
		size += 4;  // size for storing the result of 'getSizeInBytes()'

		// a byte encoding 8 bit flags; the 1st bit is reserved to encode
		// whether the vault is encrypted or not
		size += 1;

		size += 4;  // size for 'n'
		size += 4;  // size for 't'
		size += 4;  // size for 'k'
		size += 4;  // size encoding the finite fields defining polynomial

		size += 4 + slowDownFactor.getSizeInBytes();

		// size to encode either 'vaultPolynomialData' or 'encryptedVaultData'
		// depending on whether the flag for an encrypted vault is set or not.
		size += vaultDataSize();

		size += 20; // size needed to store the secret polynomial's hash

		return size;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::toBytes( uint8_t *data ) const {

		if ( !isEnrolled() ) {
			cerr << "FuzzyVault::toBytes: "
				 << "not enrolled." << endl;
			exit(EXIT_FAILURE);
		}

		int sizeInBytes = getSizeInBytes() , offset = 0;
		uint32_t tmp;

		// Write header.
		memcpy(data,"FVR",4);
		offset += 4;

		// Write size in bytes
		tmp = (uint32_t)sizeInBytes;
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write flag byte
		data[offset] = (isEncrypted()?1:0);
		offset += 1;

		// Write 'n'
		tmp = (uint32_t)(this->n);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write 'tmax'
		tmp = (uint32_t)(this->tmax);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write 'k'
		tmp = (uint32_t)(this->k);
		data[offset+3] = tmp & 0xFF; tmp >>= 8;
		data[offset+2] = tmp & 0xFF; tmp >>= 8;
		data[offset+1] = tmp & 0xFF; tmp >>= 8;
		data[offset+0] = tmp & 0xFF;
		offset += 4;

		// Write defining polynomial of finite field
		tmp = (uint32_t)(this->gfPtr->getDefiningPolynomial().rep);
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

		// Write encrypted/unencrypted vault data
		int vds = vaultDataSize();
		if ( isEncrypted() ) {
			memcpy(data+offset,this->encryptedVaultPolynomialData,vds);
		} else {
			memcpy(data+offset,this->vaultPolynomialData,vds);
		}
		offset += vds;

		// Write hash of secret polynomial
		memcpy(data+offset,this->hash,20);
		offset += 20;

		if ( offset != sizeInBytes ) {
			cerr << "FuzzyVault::toBytes: "
				 << "written bytes to no match 'getSizeInBytes()'; "
				 << "probably a bug." << endl;
			exit(EXIT_FAILURE);
		}
	}

	bool FuzzyVault::fromBytes( const uint8_t *data , int dataSize ) {

		FuzzyVault tmpVault;

		// If the data array does not at least hold enough space for
		// the header and the record size, ...
		if ( dataSize < 16 ) {
			// ..., the record cannot be valid.
			return false;

		}

		int offset = 0;

		// Check if header is correct.
		if ( memcmp(data,"FVR",4) != 0 ) {
			return false;
		}
		offset += 4;

		// Determine the size of the record.
		int sizeInBytes = (int)data[offset+0]; sizeInBytes <<= 8;
		sizeInBytes += (int)data[offset+1]; sizeInBytes <<= 8;
		sizeInBytes += (int)data[offset+2]; sizeInBytes <<= 8;
		sizeInBytes += (int)data[offset+3];
		offset += 4;

		// If the size of the record is larger than the size
		// of the specified array, ...
		if ( sizeInBytes > dataSize ) {
			// ..., then the record cannot be valid.
			return false;
		}

		// Read flag byte.
		if ( offset+1 > sizeInBytes ) { return false; }
		uint8_t flag = data[offset];
		offset += 1;

		// Read 'n'
		if ( sizeInBytes < offset+4 ) { return false; }
		tmpVault.n = (int)data[offset+0]; tmpVault.n <<= 8;
		tmpVault.n += (int)data[offset+1]; tmpVault.n <<= 8;
		tmpVault.n += (int)data[offset+2]; tmpVault.n <<= 8;
		tmpVault.n += (int)data[offset+3];
		offset += 4;

		// Read 'tmax'
		if ( sizeInBytes < offset+4 ) { return false; }
		tmpVault.tmax = (int)data[offset+0]; tmpVault.tmax <<= 8;
		tmpVault.tmax += (int)data[offset+1]; tmpVault.tmax <<= 8;
		tmpVault.tmax += (int)data[offset+2]; tmpVault.tmax <<= 8;
		tmpVault.tmax += (int)data[offset+3];
		offset += 4;

		// Read 'k'
		if ( sizeInBytes < offset+4 ) { return false; }
		tmpVault.k = (int)data[offset+0]; tmpVault.k <<= 8;
		tmpVault.k += (int)data[offset+1]; tmpVault.k <<= 8;
		tmpVault.k += (int)data[offset+2]; tmpVault.k <<= 8;
		tmpVault.k += (int)data[offset+3];
		offset += 4;

		// Read the finite field
		{
			if ( sizeInBytes < offset+4 ) { return false; }

			uint32_t rep;
			rep  = (uint32_t)data[offset+0]; rep <<= 8;
			rep += (uint32_t)data[offset+1]; rep <<= 8;
			rep += (uint32_t)data[offset+2]; rep <<= 8;
			rep += (uint32_t)data[offset+3];
			offset += 4;

			if ( this->gfPtr != NULL && this->gfPtr->getDefiningPolynomial().rep == rep ) {
				tmpVault.gfPtr = new(nothrow) SmallBinaryField(this->gfPtr[0]);
			} else {
				tmpVault.gfPtr = new(nothrow) SmallBinaryField(SmallBinaryPolynomial(rep));
			}
			if ( tmpVault.gfPtr == NULL ) {
				cerr << "FuzzyVault::fromBytes: "
					 << "out of memory." << endl;
                exit(EXIT_FAILURE);
			}

			if ( tmpVault.gfPtr->getCardinality() < (uint32_t)(tmpVault.n+tmpVault.n) ) {
				return false;
			}
		}

		// Read slow-down factor
		{
			if ( sizeInBytes < offset+4 ) { return false; }
			int sizeOfSlowDownFactor;
			sizeOfSlowDownFactor  = (int)data[offset+0]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+1]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+2]; sizeOfSlowDownFactor <<= 8;
			sizeOfSlowDownFactor += (int)data[offset+3];
			offset += 4;
			if ( sizeInBytes < offset+sizeOfSlowDownFactor ) { return false; }
			tmpVault.slowDownFactor.fromBytes(data+offset,sizeOfSlowDownFactor);
			offset += sizeOfSlowDownFactor;
		}

		// *******************************************************************
		// ******** BEGIN: Read vault data ***********************************
		// *******************************************************************

		// Test whether there is enough data provided.
		int vds = tmpVault.vaultDataSize();
		if ( sizeInBytes < offset+vds ) {
			return false;
		}

		// Allocate memory for reading encrypted/unencrypted vault data
		uint8_t *vaultData = (uint8_t*)malloc( vds * sizeof(uint8_t) );
		if ( vaultData == NULL ) {
			cerr << "FuzzyVault::fromBytes: "
				 << "out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Read data
		memcpy(vaultData,data+offset,vds);
		if ( (flag&0x1) != 0 ) {
			tmpVault.encryptedVaultPolynomialData = vaultData;
		} else {
			tmpVault.vaultPolynomialData = vaultData;
		}
		offset += vds;

		// *******************************************************************
		// ********** End: Read vault data ***********************************
		// *******************************************************************

		// Read 'hash'
		if ( sizeInBytes < offset+20 ) { return false; }
		memcpy(tmpVault.hash,data+offset,20);
		offset += 20;

		// Update the permutation which is selected pseudo-randomly using
		// the hash value as seed.
		tmpVault.updatePermutation();

		// The data has been successfully read. Thus, swap this objects'
		// content with the read record and ...
		this->swap(tmpVault);

		return true;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::write( FILE *out ) const {

		uint8_t *data;
		int size;

		size = getSizeInBytes();
		data = (uint8_t*)malloc( size );
		if ( data == NULL ) {
			cerr << "FuzzyVault::write: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		toBytes(data);

		// Write the byte array to 'out'.
		int wsize = fwrite(data,sizeof(uint8_t),size,out);

		free(data);

		return wsize == size;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::write( const std::string & filePath ) const {

		FILE *out;
		bool success;

		// Open a FILE pointer to the specified path ...
		out = IOTools::fopen(filePath.c_str(),"wb");
		if ( out == NULL ) {
			return false;
		}

		// ... and attempt to write this fuzzy vault object.
		success = write(out);

		fclose(out);

		return success;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	bool FuzzyVault::read( FILE *in ) {

		bool success = false;
		uint8_t *data;
		int size = 8;

		data = (uint8_t*)malloc( size );
		if ( data == NULL ) {
			cerr << "FuzzyVault::read: out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		if ( fread(data,1,(size_t)size,in) == (size_t)size ) {

			if ( memcmp(data,"FVR",4) == 0 ) {

				// Determine the size of the record.
				size = (int)data[4]; size <<= 8;
				size += (int)data[5]; size <<= 8;
				size += (int)data[6]; size <<= 8;
				size += (int)data[7];

				if ( size > 8 ) {

					// Reallocate memory
					data = (uint8_t*)realloc(data,size);
					if ( data == NULL ) {
						cerr << "FuzzyVault::read: out of memory." << endl;
						exit(EXIT_FAILURE);
					}

					// Read remaining data except the first 16 bytes
					// which have already been read.
					if ( fread(data+8,1,(size_t)(size-8),in) == (size_t)(size-8) ) {

						// Import the fuzzy vault object from the data.
						success = fromBytes(data,size);
					}
				}

			}

		}

		free(data);

		return success;
	}

	bool FuzzyVault::read( const std::string & filePath ) {

		bool success;

		// Open a FILE pointer to the specified path ...
		FILE *in = IOTools::fopen(filePath.c_str(),"rb");
		if ( in == NULL ) {
			return false;
		}

		// ... and attempt to initialize this fuzzy vault object
		// by the data from the file.
		if ( read(in) ) {
			success = true;
		} else {
			success = false;
		}

		fclose(in);

		return success;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	SmallBinaryFieldPolynomial FuzzyVault::createInstance
	( SmallBinaryFieldPolynomial & V ,
	  const uint32_t *features , int t , int k ,
	  bool tryRandom ) {

		int n = (int)V.getField().getCardinality();

		if ( t < 0 || k < 0 || t > n || k > n ) {
			cerr << "FuzzyVault::createInstance: "
				 << "invalid arguments." << endl;
			exit(EXIT_FAILURE);
		}

		SmallBinaryFieldPolynomial f(V.getField());
		f.random(k,tryRandom);

		V.buildFromRoots(features,t);
		V += f;

		return f;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	SmallBinaryFieldPolynomial FuzzyVault::createInstance
	( SmallBinaryFieldPolynomial & V , int t , int k ,
	  bool tryRandom ) {

		int n = (int)V.getField().getCardinality();

		if ( t < 0 || k < 0 || t > n || k > n ) {
			cerr << "FuzzyVault::createInstance: "
				 << "invalid arguments." << endl;
			exit(EXIT_FAILURE);
		}

		SmallBinaryFieldPolynomial f(V.getField());
		f.random(k,tryRandom);

		uint32_t *features = NULL;
		if ( t > 0 ) {
			features = (uint32_t*)malloc( t * sizeof(uint32_t) );
		}

		createFeatureSet(features,t,V.getField(),tryRandom);

		V.buildFromRoots(features,t);
		V += f;

		free(features);

		return f;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::createFeatureSet
	( uint32_t *features , int t , const SmallBinaryField & gf ,
	  bool tryRandom ) {

		uint32_t n = gf.getCardinality();

		if ( t > (int)n || t < 0 ) {
			cerr << "FuzzyVault::createFeatureSet: "
				 << "invalid feature set size." << endl;
			exit(EXIT_FAILURE);
		}

		for ( int i = 0 ; i < t ; i++ ) {

			features[i] = MathTools::rand32(tryRandom) % n;

			bool alreadyContained = false;
			for ( int j = 0 ; j < i ; j++ ) {
				if ( features[j] == features[i] ) {
					alreadyContained = true;
					break;
				}
			}

			if ( alreadyContained ) {
				--i;
				continue;
			}
		}

		CTools::sort(features,t);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::createOverlappingFeatureSets
	( uint32_t *A , int s , uint32_t *B , int t , int omega ,
	  uint32_t n , bool tryRandom) {

		if ( s < 0 || t < 0 || omega < 0 ) {
			cerr << "FuzzyVault: "
				 << "feature set sizes and overlap must not be negative."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( omega > s || omega > t ) {
			cerr << "FuzzyVault: "
				 << "overlap must not be larger than feature set size."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( n <= (uint32_t)s || n <= (uint32_t)t ) {
			cerr << "FuzzyVault: "
				 << "feature sizes must be smaller than universe size."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// Determine common elements
		for ( int i = 0 ; i < omega ; i++ ) {
			uint32_t x =
					(uint32_t)((uint64_t)MathTools::rand32(tryRandom) %
							((uint64_t)n));

			bool alreadyContained = false;
			for ( int j = 0 ; j < i ; j++ ) {
				if ( A[j] == x ) {
					alreadyContained = true;
					break;
				}
			}

			if ( alreadyContained ) {
				--i;
				continue;
			}

			A[i] = x;
			B[i] = x;
		}

		// Pad 'A' with disoint elements
		for ( int i = omega ; i < s ; i++ ) {

			uint32_t x =
					(uint32_t)((uint64_t)MathTools::rand32(tryRandom) %
							((uint64_t)n));

			bool alreadyContained = false;
			for ( int j = 0 ; j < i ; j++ ) {
				if ( A[j] == x ) {
					alreadyContained = true;
					break;
				}
			}

			if ( alreadyContained ) {
				--i;
				continue;
			}

			A[i] = x;
		}

		// Pad 'B' with disjoint elements
		for ( int i = omega ; i < t ; i++ ) {
			uint32_t x =
					(uint32_t)((uint64_t)MathTools::rand32(tryRandom) %
							((uint64_t)n));

			bool alreadyContained = false;
			for ( int j = 0 ; j < s ; j++ ) {
				if ( A[j] == x ) {
					alreadyContained = true;
					break;
				}
			}
			if ( alreadyContained ) {
				--i;
				continue;
			}

			for ( int j = 0 ; j < i ; j++ ) {
				if ( B[j] == x ) {
					alreadyContained = true;
					break;
				}
			}

			if ( alreadyContained ) {
				--i;
				continue;
			}

			B[i] = x;
		}

		// sort feature sets in ascending order
		CTools::sort(A,s);
		CTools::sort(B,t);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::numOverlap
	( uint32_t *A , int s , uint32_t *B , int t ) {

		// Counts the number of common elements
		int omega = 0;

		// Iterate over all pairs
		for ( int i = 0 ; i < s ; i++ ) {
			for ( int j = 0 ; j < t ; j++ ) {

				// Increment in case the two elements are equal
				if ( A[i] == B[j] ) {
					++omega;
				}
			}
		}

		return omega;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::first_init() {
		this->n = 0;
		this->tmax = 0;
		this->k = 0;
		this->gfPtr = NULL;
		this->slowDownFactor = 1;
		this->vaultPolynomialData = NULL;
		this->encryptedVaultPolynomialData = NULL;
		memset(this->hash,0,20);
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::updateField() {

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
				cerr << "FuzzyVault:updateField: "
					 << "out of memory." << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::updatePermutation() {

		int n = getVaultSize();
        this->permutation.setDimension(n);

		RandomGenerator gen(this->hash);

		for ( int x0 = 0 ; x0 < n ; x0++ ) {
			int x1 = (int)(gen.rand() % n);
			this->permutation.exchange(x0,x1);
		}
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::packVaultPolynomial( const SmallBinaryFieldPolynomial & V ) {

		// Random slow-down value
		BigInteger slowDownVal; slowDownVal.random(this->slowDownFactor,true);

		// Derive AES key from slow-down value.
		AES128 aes = deriveKey(slowDownVal);

		int t , d , n;
		t = this->tmax;
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
			cerr << "FuzzyVault::packVaultPolynomialData: out of memory." << endl;
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

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	uint32_t FuzzyVault::_reorder( uint32_t a ) const {
		return (uint32_t)(this->permutation.eval((int)a));
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	int FuzzyVault::vaultDataSize() const {

		int t , d;
		t = this->tmax;
		d = this->gfPtr->getDegree();

		int b = t * d;

		int n = b / 8 + (b%8?1:0);

		if ( n % 16 != 0 ) {
			n += 16-n%16;
		}

		return n;
	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::checkFeatureSet( const uint32_t *features , int t ) const {

		if ( t < 0 ) {
			cerr << "FuzzyVault::checkFeatureSet: "
				 << "size must not be negative." << endl;
			exit(EXIT_FAILURE);
		}

		for ( int i = 0 ; i < t ; i++ ) {

			if ( features[i] >= (uint32_t)(this->n) ) {
				cerr << "FuzzyVault::checkFeatureSet: "
					 << "feature element must not be larger than "
					 << "'getVaultSize()'." << endl;
				exit(EXIT_FAILURE);
			}

			for ( int j = i+1 ; j < t ; j++ ) {

				if ( features[i] == features[j] ) {
					cerr << "FuzzyVault::checkFeatureSet: "
					     << "feature set contains multiple elements." << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

	}

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	AES128 FuzzyVault::deriveKey( const BigInteger & x ) {

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

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::concat_bit_vectors
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

	/*
	 * see 'FuzzyVault.h' for the documentation.
	 */
	void FuzzyVault::split_into_bit_vectors
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
