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
 * @file SmallBinaryField.cpp
 *
 * @brief
 *            Implements the functionalities that are provided by
 *            'SmallBinaryField.h' to provides functionalities
 *            enabling arithmetic in binary galois fields that are of small
 *            cardinality.
 *
 * @author Benjamin Tams
 *
 * @see thimble::SmallBinaryField
 */

#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryPolynomial.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *           Copy constructor.
	 *
	 * @param gf
	 *           The binary finite field of which a copy is constructed.
	 */
	SmallBinaryField::SmallBinaryField( const SmallBinaryField & gf ) {

		this->definingPolynomial = SmallBinaryPolynomial(0);
		this->degree = -1;
		this->size = 0;
		this->generator = 0;
		this->expTable = NULL;
		this->logTable = NULL;
		this->invTable = NULL;

		*this = gf;
	}

	/**
	 * @brief
	 *           Assignment operator.
	 *
	 * @param gf
	 *           The binary finite field of which a copy is assigned
	 *           to this instance.
	 */
	SmallBinaryField & SmallBinaryField::operator=( const SmallBinaryField & gf ) {

		if ( this->getDefiningPolynomial().rep != gf.getDefiningPolynomial().rep ) {

			clear();

			this->definingPolynomial = gf.definingPolynomial;
			this->degree = gf.degree;
			this->size = gf.size;
			this->generator = gf.generator;

			// Allocate memory for the tables
			this->expTable = (uint32_t*)malloc
				( (size_t)(this->size * sizeof(uint32_t)) );
			this->logTable = (uint32_t*)malloc
				( (size_t)(this->size * sizeof(uint32_t)) );
			this->invTable = (uint32_t*)malloc
				( (size_t)(this->size * sizeof(uint32_t)) );

			// Check whether sufficient memory could be allocated
			if ( this->expTable == NULL || this->logTable == NULL ||
				 this->invTable == NULL ) {
				cerr << "SmallBinaryField: out of memory" << endl;
				exit(EXIT_FAILURE);
			}

			memcpy(this->expTable,gf.expTable,this->size * sizeof(uint32_t) );
			memcpy(this->logTable,gf.logTable,this->size * sizeof(uint32_t) );
			memcpy(this->invTable,gf.invTable,this->size * sizeof(uint32_t) );
		}

		return *this;
	}

	/**
	 * @brief
	 *           Swap method.
	 *
	 * @details
	 *           After calling this method, the field <code>gf1</code>
	 *           and <code>gf2</code> encode the fields encoded by
	 *           <code>gf2></code> and <code>gf1</code> on input,
	 *           respectively,
	 *
	 * @param gf1
	 *           First binary field.
	 *
	 * @param gf2
	 *           Second binary field.
	 */
	void SmallBinaryField::swap( SmallBinaryField & gf1 , SmallBinaryField & gf2 ) {

		SmallBinaryPolynomial definingPolynomial;
		int degree;
		uint32_t size;
		uint32_t generator;
		uint32_t *expTable , *logTable , *invTable;

		definingPolynomial = gf1.definingPolynomial;
		degree = gf1.degree;
		size = gf1.size;
		generator = gf1.generator;
		expTable = gf1.expTable;
		logTable = gf1.logTable;
		invTable = gf1.invTable;

		gf1.definingPolynomial = gf2.definingPolynomial;
		gf1.degree = gf2.degree;
		gf1.size = gf2.size;
		gf1.generator = gf2.generator;
		gf1.expTable = gf2.expTable;
		gf1.logTable = gf2.logTable;
		gf1.invTable = gf2.invTable;

		gf2.definingPolynomial = definingPolynomial;
		gf2.degree = degree;
		gf2.size = size;
		gf2.generator = generator;
		gf2.expTable = expTable;
		gf2.logTable = logTable;
		gf2.invTable = invTable;
	}

	/**
	 * @brief    Frees all memory that is allocated by this finite field
	 *           to list the \link expTable\endlink,
	 *           \link logTable\endlink, and \link invTable\endlink
	 *
	 * @details
	 *           see 'SmallBinaryField.h'
	 */
	void SmallBinaryField::clear() {
		free(this->expTable);
		free(this->logTable);
		free(this->invTable);
		this->expTable = NULL;
		this->logTable = NULL;
		this->invTable = NULL;
	}

	/**
	 * @brief    Initializes the finite field with the defining polynomial
	 *           <code>f</code>
	 *
	 * @details
	 *           see 'SmallBinaryField.h'
	 */
	void SmallBinaryField::init( const SmallBinaryPolynomial & f ) {

		// Check whether the polynomial degree is right
		if ( f.deg() > 31 || f.rep < 2 ) {
			cerr << "SmallBinaryField: Defining polynomial must be of degree "
				 << "between 1 and 31" << endl;
			exit(EXIT_FAILURE);
		}

		// Set class members that are already known
		this->definingPolynomial = f;
		this->degree   = f.deg();
		this->size     = 1;
		this->size   <<= this->degree;

		// Allocate memory for the tables
		this->expTable = (uint32_t*)malloc
			( (size_t)(this->size * sizeof(uint32_t)) );
		this->logTable = (uint32_t*)malloc
			( (size_t)(this->size * sizeof(uint32_t)) );
		this->invTable = (uint32_t*)malloc
			( (size_t)(this->size * sizeof(uint32_t)) );

		// Check whether sufficient memory could be allocated
		if ( this->expTable == NULL || this->logTable == NULL ||
			 this->invTable == NULL ) {
			cerr << "SmallBinaryField: Out of memory" << endl;
			exit(EXIT_FAILURE);
		}

		// This for-loop iterates of randomly chosen candidates
		// for the multiplicative group's generator. It is checked
		// within the body, whether a correct generator could be
		// found; if not, the loop is repeated until a correct
		// generator could be generated. If the defining polynomial
		// is reducible, it is not possible to find a generator.
		// Therefore, this loop will not terminate in such a case.
		for(;;) {

			// Generate a random nonzero polynomial
			// representing a candidate for the generator
			// of the fields multiplicative group of
			// unity
			SmallBinaryPolynomial x;
			do {
				x = SmallBinaryPolynomial::rand(this->degree);
			} while ( x.rep == 0 );

			 // The 0-th power will be 1
			this->expTable[0] = 1;
			// Actually, the logarithm of 0 is not defined.
			// But in order to ensure, that every entry is set
			// we set it to zero
			this->logTable[0] = 0;

			// 'y' will contain the 'i'-th power of the generator
			// 'x mod f' during the following for-loop
			SmallBinaryPolynomial y(1);

			uint64_t i;
			for ( i = 1 ; i < this->size ; i++ ) {

				// ' x^i <-- (x^{i-1} * x) mod f'
				y = y.mulMod(x,f);

				// The 'i'-th power of the (alleged) generator
				this->expTable[i] = (uint32_t)y.rep;

				// 'i' is the discrete logarithm of the element
				// represented by 'y mod f' to the base of the
				// generator
				this->logTable[y.rep] = (uint32_t)i;

				// If '(x mod f)^i=1' the element 'x mod f' generates
				// a group of size 'i'.
				if ( y.rep == 1 ) {
					break;
				}
			}

			// Check whether 'i' is equals to the size of the
			// finite fields multiplicative group of unity; if true,
			// 'x' represents a generator and the loop is terminated;
			// otherwise, the loop continues
			if ( i == this->size-1 ) {
				this->generator = (uint32_t)x.rep;
				break;
			}
		}

		// Pre-compute the inverses

		// Actually, the inverse of zero is not defined
		// but in order to ensure every entry is initialized
		// we set it as zero. As a consequence, if the user
		// attempts to compute the inverse of 0 the result will
		// be 0.
		this->invTable[0] = 0;
		for ( uint64_t a = 1 ; a < this->size ; a++ ) {
			// If 'a=x^i' then 'a^(-1)=x^(size-1-i)' where 'x' denotes
			// a fields multiplicative group of unity's generator
			this->invTable[a] = this->expTable[this->size-1-this->logTable[a]];
		}

	}

	/**
	 * @brief
	 *            Returns a constant reference to a field
	 *            with two elements.
	 *
	 * @return
	 *            Constant reference to a field with two elements.
	 */
	const SmallBinaryField & SmallBinaryField::binary() {
		static SmallBinaryField gf(2);
		return gf;

	}
}



