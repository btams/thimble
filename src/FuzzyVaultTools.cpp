/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013 Benjamin Tams
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
 * @file FuzzyVaultTools.cpp
 *
 * @brief
 *            Provides convenience functionalities related with the fuzzy
 *            vault scheme.
 *
 * @details
 *            see 'FuzzyVaultTools.h'
 *
 * @author Benjamin Tams
 */
#include <stdint.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include <thimble/math/BinomialIterator.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/math/numbertheory/BigInteger.h>
#include <thimble/security/SHA.h>
#include <thimble/security/FuzzyVaultTools.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Choose a subset of specified size from a finite field.
	 *
	 * @details
	 *            The method selects <code>n</code> elements from
	 *            <code>gf</code> and stores them successively in
	 *            <code>x</code>. Before the <code>i</code>-th element
	 *            is updated, it is checked whether it is already contained
	 *            in the list. If <code>true</code> another <code>i</code>-th
	 *            element is selected.
	 *            <br><br>
	 *            As a consequence, this method assumes that the size of the
	 *            finite field specified by <code>gf</code> is larger than
	 *            (or equals) <code>n</code>.
	 *            <br><br>
	 *            If <code>tryRandom</code> is <code>true</code> the method
	 *            is advised to use a random generator with more entropy;
	 *            otherwise, it wraps around the standard <code>rand()</code>
	 *            function.
	 *
	 * @param x
	 *            Will contain <code>n</code> elements of the finite field
	 *            specified by <code>gf</code>.
	 *
	 * @param n
	 *            The size of the randomly selected set.
	 *
	 * @param gf
	 *            The finite field from where elements are selected.
	 *
	 * @param tryRandom
	 *            If <code>true</code> the method is advised to use a
	 *            random generator with more entropy; otherwise, it wraps
	 *            around the standard <code>rand()</code> function.
	 *
	 * @warning
	 *            If <code>x</code> can not store at least <code>n</code>
	 *            elements from the finite field <code>gf</code> the method
	 *            may run into unexpected, undocumented behavior.
	 *
	 * @warning
	 *            If not sufficient memory could be provided, the method
	 *            prints an error message to <code>stderr</code> and exits
	 *            with status 'EXIT_FAILURE'.
	 */
	static void chooseFiniteFieldSetAtRandom
	( uint32_t *x , int n , const SmallBinaryField & gf , bool tryRandom ) {

		uint32_t size = gf.getCardinality();
		uint32_t *elements = (uint32_t*)malloc( size*sizeof(uint32_t));
		if ( elements == NULL ) {
			cerr << "FuzzyVaultTools: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		for ( uint32_t i = 0 ; i < size ; i++ ) {
			elements[i] = i;
		}

		for ( int i = 0 ; i < n ; i++ ) {

			uint32_t j = MathTools::rand32(tryRandom) % (uint32_t)(size-i);

			x[i] = elements[j];
			for ( ; j+1 < size ; j++ ) {
				elements[j] = elements[j+1];
			}
		}

		free(elements);

	}

	/**
	 * @brief
	 *           Chooses <code>k</code> pairwise different integers from the
	 *           range <code>0,...,n-1</code> at random.
	 *
	 * @param indices
	 *           Will contain <code>k</code> pairwise distinct integers
	 *           in the range <code>0,...,n-1</code>.
	 *
	 * @param n
	 *           Specifies the range from where to selected integer.
	 *
	 * @param k
	 *           Specifies the number of integers to be selected.
	 *
	 * @param tryRandom
	 *            If <code>true</code> the method is advised to use a
	 *            random generator with more entropy; otherwise, it wraps
	 *            around the standard <code>rand()</code> function.
	 *
	 * @warning
	 *            If <code>k</code> is greater than or equals <code>n</code>
	 *            the method may never terminate.
	 *
	 * @warning
	 *            If <code>indices</code> cannot store at least <code>k</code>
	 *            integers, the method may run into unexpected, undocumented
	 *            behavior.
	 *
	 * @warning
	 *            If not sufficient memory could be provided, the method
	 *            prints an error message to <code>stderr</code> and exits
	 *            with status 'EXIT_FAILURE'.
	 */
	static void chooseIndicesAtRandom
	( register int *indices , register int n , register int k , bool tryRandom ) {

		register int i;
		register int j;
		register int *range;

		range = (int*)malloc( n * sizeof(int) );
		if ( range == NULL ) {
			cerr << "FuzzyVaultTools: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		for ( i = 0 ; i < n ; i++ ) {
			range[i] = i;
		}

		for ( i = 0 ; i < k ; i++ ) {

			j = MathTools::rand32(tryRandom) % (n-i);
			indices[i] = range[j];

			for ( ; j+1 < n ; j++ ) {
				range[j] = range[j+1];
			}
		}

		free(range);
	}

	/**
	 * @brief
	 *           Chooses <code>k</code> pairwise different integers from the
	 *           range <code>0,...,n-1</code> at random.
	 *
	 * @details
	 *           The method essentially performs the same as the
	 *           \link chooseIndicesAtRandom(int*,int,int,bool)\endlink method
	 *           but assumes <code>tryRandom=false</code>. This enables to use
	 *           <code>rand()</code> directly without frequently considering
	 *           the cases for <code>tryRandom</code>. This may result in a
	 *           more efficient version which is beneficial for, e.g., the
	 *           brute-force attack.
	 *
	 * @param indices
	 *           Will contain <code>k</code> pairwise distinct integers
	 *           in the range <code>0,...,n-1</code>.
	 *
	 * @param n
	 *           Specifies the range from where to selected integer.
	 *
	 * @param k
	 *           Specifies the number of integers to be selected.
	 *
	 * @warning
	 *            If <code>k</code> is greater than or equals <code>n</code>
	 *            the method may never terminate.
	 *
	 * @warning
	 *            If <code>indices</code> cannot store at least <code>k</code>
	 *            integers, the method may run into unexpected, undocumented
	 *            behavior.
	 */
	static void chooseIndicesAtRandom
	( register int *indices , register int n , register int k ) {

		register bool alreadyChosen;
		register int index;
		register int i;
		register int j;

		for ( i = 0 ; i < k ; i++ ) {

			do {
				index = rand() % n;
				alreadyChosen = false;
				for ( j = 0 ; j < i ; j++ ) {
					if ( indices[j] == index ) {
						alreadyChosen = true;
						break;
					}
				}
			} while ( alreadyChosen );

			indices[i] = index;
		}
	}

	/**
	 * @brief
	 *            Creates a random instance of a fuzzy vault with
	 *            specified parameters
	 *
	 * @details
	 *            see 'FuzzyVaultTools.h'
	 */
	SmallBinaryFieldPolynomial FuzzyVaultTools::createRandomInstance
	( uint32_t *x , uint32_t *y , int n , int t , int k ,
	  const SmallBinaryField & gf , bool tryRandom ) {

		if ( (uint64_t)n > gf.getCardinality() ) {
			cerr << "FuzzyVault::createRandomInstance: "
				 << "It is not possible to create a random vault of size "
				 << "larger than the finite field's cardinality." << endl;
			exit(EXIT_FAILURE);
		}

		// Allocate memory for choosing the genuine points
		// and check if this was successful
		int *genuineIndices;
		genuineIndices = (int*)malloc(t*sizeof(int));

		// List in 'x' the elements of a randomly
		// chosen subset of the field that is of size
        // 'n', i.e. a list without multiples.
		chooseFiniteFieldSetAtRandom(x,n,gf,tryRandom);

		// Generate the random vault polynomial
		SmallBinaryFieldPolynomial f(gf);
		f.random(k,tryRandom);

		// Generate the ordinate such that, for now, all
		// of the resulting vault points do not lie
		// on the polynomial's graph
		for ( int i = 0 ; i < n ; i++ ) {

			uint32_t z = f.eval(x[i]);
			do {
				y[i] = gf.random(tryRandom);
			} while ( z == y[i] );
		}

		// Now, select those indices of the vault points
		// that will be genuine points
		chooseIndicesAtRandom(genuineIndices,n,t,tryRandom);

		// Update the genuine vault point's ordinate values
		// by letting the the evaluation of the polynomial
		// at the corresponding x-value
		for ( int i = 0 ; i < t ; i++ ) {
			int j = genuineIndices[i];
			y[j] = f.eval(x[j]);
		}

		// Free the memory
		free(genuineIndices);

		// Output polynomial
		return f;
	}

    /**
     * @brief
     *            Creates a random instance of the improved fuzzy vault scheme.
     *
     * @details
     *            The function generates a random polynomial <i>f</i> of degree
     *            smaller <i>k</i>; furthermore, a product
     *            \f[
     *             \prod_{j=0}^{t-1}(X-x_j)
     *            \f]
     *            where the \f$x_j\f$ are selected to be random and distinct.
     *            Then the function outputs
     *            \f[
     *             V(X)=f(X)+\prod_{j=0}^{t-1}(X-x_j).
     *            \f]
     *            The function also returns <i>f</i> to enable evaluations of
     *            attacks in order to decide whether they were successful or
     *            not.
     *
     * @param V
     *            The random vault polynomial.
     *
     * @param t
     *            The number of random genuine features \f$x_j\f$ protected by
     *            the vault.
     *
     * @param k
     *            The length of the vault's secret polynomial.
     *
     * @param tryRandom
     *            If <code>true</code>, the function is advised to use a
     *            cryptographic number generator; otherwise, the method
     *            wraps around the standard <code>rand()</code> function.
     *
     * @warning
     *            If the relation \f$0\leq k<t\leq n\f$, where \f$n\f$ is
     *            the cardinality of the field over which \f$V\f$ is
     *            defined, the function prints an error message to
     *            <code>stderr</code> and exits with status
     *            'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    SmallBinaryFieldPolynomial FuzzyVaultTools::createRandomInstance
    ( SmallBinaryFieldPolynomial & V , int t , int k , bool tryRandom ) {

        if ( k < 0 || t < 0 ) {
            cerr << "FuzzyVaultTools::createRandomInstance: "
                 << "number of genuine features and length of secret "
                 << "polynomial must be positive." << endl;
            exit(EXIT_FAILURE);
        }

        if ( k >= t ) {
            cerr << "FuzzyVaultTools::createRandomInstance: "
                 << "length of secret polynomial must be smaller than "
                 << "number of genuine features." << endl;
            exit(EXIT_FAILURE);
        }

        if ( (uint64_t)t > (uint64_t)V.getField().getCardinality() ) {

            cerr << "FuzzyVaultTools::createRandomInstance: "
                 << "number of genuine features too large." << endl;
            exit(EXIT_FAILURE);
        }

        // Generate the secret polynomial
        SmallBinaryFieldPolynomial f(V.getField());
        f.random(k,tryRandom);
        // Generate a random subset '{x[j]}' of the field of size 't'
        uint32_t *x = (uint32_t*)malloc( t * sizeof(uint32_t) );
        chooseFiniteFieldSetAtRandom(x,t,V.getField(),tryRandom);
        // Computes the product 'V(X)=(X-x[0])*...*(X-x[t-1])'
        V.buildFromRoots(x,t);
        // Output vault polynomial
        add(V,V,f);

        free(x);

        return f;
    }

	/**
	 * @brief
	 *            Attempts to break an instance of the fuzzy vault scheme.
	 *
	 * @details
	 *            see 'FuzzyVaultTools.h'
	 */
	bool FuzzyVaultTools::bfattack
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y ,
	  int n , int k , const uint32_t hash[5] ,
	  uint64_t maxIts ) {

        SHA sha;

		// Check whether we can choose random points from 'n'
		// points using 'rand()'
		if ( n > RAND_MAX ) {
			cerr << "FuzzyVault::bfattack: The number of vault points must be"
				 << " smaller than or equal RAND_MAX which is"
				 << RAND_MAX << "." << endl;
			exit(EXIT_FAILURE);
		}

		// Check whether the vault is of reasonable parameters
		if ( n <= 0 || k <= 0 || k > n ) {
			cerr << "FuzzyVault::bfattack: The number of vault points must be"
				 << " greater than zero. Furthermore, the size of the secret "
				 << "polynomial must be greater than zero and smaller than "
				 << "(or equal) to the vault's size" << endl;
			exit(EXIT_FAILURE);
		}

		// Keeps track whether a polynomial was yet found or not
		bool state = false;

		// Initialize space for the candidate polynomial
		SmallBinaryFieldPolynomial candidatePolynomial(f.getField());
		candidatePolynomial.ensureCapacity(k);

		// Initalize space for the hash of the candidate polynomial
		uint32_t candidateHash[5];

		uint32_t *a , *b;
		int *indices;

		// Allocate memory to select 'k' random vault
		// points
		a = (uint32_t*)malloc( k * sizeof(uint32_t) );
		b = (uint32_t*)malloc( k * sizeof(uint32_t) );
		indices = (int*)malloc( k * sizeof(uint32_t) );
		if ( a == NULL || b == NULL || indices == NULL ) {
			cerr << "FuzzyVault::bfattack: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Iterate at most 'maxIts' times
		for ( uint64_t it = 0 ; it < maxIts ; it++ ) {

			// Select pairwise different indices in the range
			// '0,...,n-1' and ...
			chooseIndicesAtRandom(indices,n,k);

			// ... set the selected vault points, correspondingly.
			for ( int i = 0 ; i < k ; i++ ) {
				int j = indices[i];
				a[i] = x[j];
				b[i] = y[j];
			}

			// Determine the interpolation polynomial of the selected
			// vault points and ...
			candidatePolynomial.interpolate(a,b,k);
			// ... compute its SHA-1 hash value
            sha.hash
            (candidateHash,
             candidatePolynomial.getData(),candidatePolynomial.deg()+1);

			// Check whether the candidate polynomial's hash value
			// agrees with the hash value of the secret polynomial.
			if ( memcmp(candidateHash,hash,20) == 0 ) {
				// If true, assign 'f', update the 'state' and abort
				// the loop.
				f.assign(candidatePolynomial);
				state = true;
				break;
			}
		}

		// Free memory
		free(a);
		free(b);
		free(indices);

		return state;
	}

	/**
		 * @brief
		 *            Attempts to break an instance of the fuzzy vault scheme.
		 *
		 * @details
		 *            see 'FuzzyVaultTools.h'
		 */
		bool FuzzyVaultTools::bfattack
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int n , int k , const uint8_t hash[20] ,
		  uint64_t maxIts ) {

            SHA sha;

			// Check whether we can choose random points from 'n'
			// points using 'rand()'
			if ( n > RAND_MAX ) {
				cerr << "FuzzyVault::bfattack: The number of vault points must be"
					 << " smaller than or equal RAND_MAX which is"
					 << RAND_MAX << "." << endl;
				exit(EXIT_FAILURE);
			}

			// Check whether the vault is of reasonable parameters
			if ( n <= 0 || k <= 0 || k > n ) {
				cerr << "FuzzyVault::bfattack: The number of vault points must be"
					 << " greater than zero. Furthermore, the size of the secret "
					 << "polynomial must be greater than zero and smaller than "
					 << "(or equal) to the vault's size" << endl;
				exit(EXIT_FAILURE);
			}

			// Keeps track whether a polynomial was yet found or not
			bool state = false;

			// Initialize space for the candidate polynomial
			SmallBinaryFieldPolynomial candidatePolynomial(f.getField());
			candidatePolynomial.ensureCapacity(k);

			// Initalize space for the hash of the candidate polynomial
			uint8_t candidateHash[20];

			uint32_t *a , *b;
			int *indices;

			// Allocate memory to select 'k' random vault
			// points
			a = (uint32_t*)malloc( k * sizeof(uint32_t) );
			b = (uint32_t*)malloc( k * sizeof(uint32_t) );
			indices = (int*)malloc( k * sizeof(uint32_t) );
			if ( a == NULL || b == NULL || indices == NULL ) {
				cerr << "FuzzyVault::bfattack: Out of memory." << endl;
				exit(EXIT_FAILURE);
			}

			// Iterate at most 'maxIts' times
			for ( uint64_t it = 0 ; it < maxIts ; it++ ) {

				// Select pairwise different indices in the range
				// '0,...,n-1' and ...
				chooseIndicesAtRandom(indices,n,k);

				// ... set the selected vault points, correspondingly.
				for ( int i = 0 ; i < k ; i++ ) {
					int j = indices[i];
					a[i] = x[j];
					b[i] = y[j];
				}

				// Determine the interpolation polynomial of the selected
				// vault points and ...
				candidatePolynomial.interpolate(a,b,k);
				// ... compute its SHA-1 hash value
                sha.hash
                (candidateHash,
                 candidatePolynomial.getData(),candidatePolynomial.deg()+1);

				// Check whether the candidate polynomial's hash value
				// agrees with the hash value of the secret polynomial.
				if ( memcmp(candidateHash,hash,20) == 0 ) {
					// If true, assign 'f', update the 'state' and abort
					// the loop.
					f.assign(candidatePolynomial);
					state = true;
					break;
				}
			}

			// Free memory
			free(a);
			free(b);
			free(indices);

			return state;
		}

	/**
	 * @brief
	 *            Attempts to decode a polynomial of degree smaller than
	 *            <i>k</i> that interpolates <i>k</i> given points and
	 *            that is of specified hash value.
	 *
	 * @details
	 *            see 'FuzzyVaultTools.h'
	 */
	bool FuzzyVaultTools::bfdecode
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y ,
	  int n , int k , const uint32_t hash[5] ) {

		if ( n < 0 || k < 0 ) {
			cerr << "FuzzyVaultTools::bfdecode: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		if ( k > n ) {
			return false;
		}

        SHA sha;
		SmallBinaryFieldPolynomial candidatePolynomial(f.getField());
		uint32_t candidateHash[5];

		// Special case: The zero polynomial is not interpolated by any
		// points
		if ( k == 0 ) {
            sha.hash(candidateHash,f.getData(),f.deg()+1);
			return memcmp(candidateHash,hash,20)==0;
		}

		// Keeps track whether a polynomial was yet found or not
		bool state = false;
		uint32_t *a , *b;

		// Allocate memory to select 'k' random vault
		// points
		a = (uint32_t*)malloc( k * sizeof(uint32_t) );
		b = (uint32_t*)malloc( k * sizeof(uint32_t) );
		if ( a == NULL || b == NULL ) {
			cerr << "FuzzyVaultTools::bfdecode: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// The iterator that successfully iterates through all
		// possible states of selected 'k' out of 'n' values
		BinomialIterator<uint32_t> it(n,k);

		do {

			it.select(a,x);
			it.select(b,y);

			// Determine the interpolation polynomial of the selected
			// vault points and ...
			candidatePolynomial.interpolate(a,b,k);
			// ... compute its SHA-1 hash value
            sha.hash
            (candidateHash,
             candidatePolynomial.getData(),candidatePolynomial.deg()+1);

			// Check whether the candidate polynomial's hash value
			// agrees with the hash value of the secret polynomial.
			if ( memcmp(candidateHash,hash,20) == 0 ) {
				// If true, assign 'f', update the 'state' and abort
				// the loop.
				f.assign(candidatePolynomial);
				state = true;
				break;
			}

			// switch to the next state of the iterator; if the iterator
			// has already reached its final state, we are done.
		} while( it.next() );

		// Free memory
		free(a);
		free(b);

		return state;
	}

	/**
	 * @brief
	 *           For a random fuzzy vault of specified parameters,
	 *           the function computes the difficulty for choosing
	 *           the <i>k</i> genuine points in the vault.
	 *
	 * @details
	 *           see 'FuzzyVaultTools.h'
	 *
	 * @warning
	 *           In the following case an error message will be printed
	 *           to <code>stderr</code> and the program will exit with
	 *           status 'EXIT_FAILURE'.
	 *           <ul>
	 *            <li>If either of the parameters is negative.</li>
	 *            <li>If \f$k>t\f$ or if \f$t>n\f$.</li>
	 *            <li>
	 *             If no sufficient memory can be allocated for computing the
	 *             result by this function.
	 *            </li>
	 *           </ul>
	 */
	long double FuzzyVaultTools::bfsecurity( int n , int t , int k ) {

		BigInteger num = BigInteger::binomial(n,k);
		BigInteger den = BigInteger::binomial(t,k);

		return num.toFloat() / den.toFloat();
	}


}



