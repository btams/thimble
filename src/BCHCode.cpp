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
 * @file BCHCode.cpp
 *
 * @brief
 *            Implementation of a class for generating binary BCH codes with
 *            extended functionality as provided by the 'BCHCode.h' header.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <iostream>

#include <thimble/math/linalg/BinaryVector.h>
#include <thimble/math/linalg/BinaryMatrix.h>
#include <thimble/math/linalg/LinAlgTools.h>
#include <thimble/math/numbertheory/BinaryPolynomial.h>
#include <thimble/math/numbertheory/NTTools.h>
#include <thimble/ecc/BCHCodeBase.h>
#include <thimble/ecc/BCHCode.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Computes the generator matrix and the check matrix of
	 *            this BCH code.
	 *
	 * @details
	 *            Internally, the error correction algorithmic of BCH
	 *            code is enabled through algebraic properties of
	 *            polynomials. BCH codes for binary polynomials
	 *            are represented by the \link BCHCodeBase\endlink class
	 *            from which instances of this class inherit properties.
	 *            As all properties of the \link BCHCodeBase\endlink
	 *            class have been set, the generator matrix
	 *            \link G\endlink and check matrix \link G\endlink remain
	 *            to be computed. The computation of these matrices of
	 *            performed by this method.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
	void BCHCode::compute_generator_and_check_matrices() {

		int n , k;
		n = getBlockLength();
		k = getDimension();

        // *******************************************************************
        // *** BEGIN: build the generator matrix *****************************
        // *******************************************************************
        this->G.setDimension(n,k);
        {
            BinaryPolynomial m , c;
            for ( int j = k-1 ; j >= 0 ; j-- ) {

                m.setCoeff(j);

                c.assign(m);
                BCHCodeBase::encode(c);

                for ( int i = 0 ; i < n ; i++ ) {
                    if ( c.getCoeff(i) ) {
                        this->G.setAt(i,j);
                    }
                }

                m.clearCoeff(j);
            }
        }
        // *******************************************************************
        // *** END: generator matrix has been built **************************
        // *******************************************************************


        // *******************************************************************
        // *** BEGIN: build the check matrix *********************************
        // *******************************************************************

        // Since the generator matrix is systematic, i.e., of the form
        // 'G=transpose(P|I)', the check matrix equals 'H=(I|P)'.
        this->H.setDimension(n-k,n);
        for ( int i = 0 ; i < n-k ; i++ ) {
            this->H.setAt(i,i);
            for ( int j = 0 ; j < k ; j++ ) {
                if (this->G.getAt(i,j) ) {
                    this->H.setAt(i,j+n-k);
                }
            }
        }

        // *******************************************************************
        // *** END: check matrix has been built ******************************
        // *******************************************************************
	}

	/**
	 * @brief
	 *            Creates a binary BCH code of given block length
	 *            that can tolerate at least the specified number of
	 *            errors.
	 *
	 * @details
	 *            More specifically, a BCH code of block length <i>n</i>
	 *            is created being of maximal dimension <i>k</i> where
	 *            its minimal distance <i>d</i> is greater than or equals
	 *            <i>2*errorTolerance+1</i>.
	 *            <br><br>
	 *            This constructor is wrapped around the constructor
	 *            \link BCHCodeBase(int,int)\endlink.
	 *
	 * @see BCHCodeBase(int,int)
	 *
	 * @param n
	 *            The block length of the BCH code.
	 *
	 * @param numErrors
	 *            A lower bound on the number of errors that the
	 *            code should be able to correct.
	 *
	 * @warning
	 *            If <i>n</i> is negative or even or if the relation
	 *            \f$1\leq 2\cdot errorTolerance+1\leq n\f$ is not
	 *            fulfilled, an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not enough memory could be provided, an error
	 *            message is printed to <code>stderr</code> and the
	 *            program exits with status 'EXIT_FAILURE'.
	 */
    BCHCode::BCHCode( int n , int numErrors ) :
        BCHCodeBase(n,numErrors) {

    	this->compute_generator_and_check_matrices();
    }

    /**
     * @brief
     *            Copy constructor.
     *
     * @details
     *            Creates a BCH code for vectors being an equivalent
     *            for a BCH code for polynomials.
     *
     * @param code
     *            The BCH code (for polynomials) of which a copy
     *            is constructed.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    BCHCode::BCHCode( const BCHCodeBase & code ) :
    		BCHCodeBase(code) {

    	this->compute_generator_and_check_matrices();
    }

    /**
     * @brief
     *            Copy constructor.
     *
     * @param code
     *            The BCH code of which a copy is created.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    BCHCode::BCHCode( const BCHCode & code ) : BCHCodeBase(code) {

    	this->G = code.G;
    	this->H = code.H;
    }

    /**
     * @brief
     *            Destructor.
     *
     * @details
     *            Ensures that all data helt by this instance is
     *            freed.
     */
    BCHCode::~BCHCode() {
        // Every member's class has a destructor. It is actually not
        // necessary to implement the destructor except if we attempt
        // to emphasize that the memory is in fact freed.
    }

    /**
     * @brief
     *            Attempts to round the input vector to its
     *            nearest codeword.
     *
     * @details
     *            If the BCH code contains a word <i>w</i> that differs
     *            in no more than \link getErrorTolerance()\endlink bits
     *            from  <i>c</i>, then the content of <i>c</i>
     *            is replaced by <i>w</i> and the function returns
     *            <code>true</code>. Otherwise, the function is expected
     *            to return <code>false</code> and leaves <i>c</i>
     *            unchanged.
     *
     * @param c
     *            On input and output a vector of length
     *            \link getBlockLength()\endlink.
     *
     * @return
     *            <code>true</code> if <i>c</i> was successfully
     *            rounded to a codeword; otherwise, <code>false</code>.
     *
     * @warning
     *            If <i>c</i> is of degree larger than or equals
     *            \link getBlockLength()\endlink an error message is
     *            printed to <code>stderr</code> and the program exits
     *            with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    bool BCHCode::round( BinaryVector & c ) const {

        int n = getBlockLength();

        if ( c.getLength() != n ) {
            cerr << "BCHCode::round: "
                 << "vector does not match block length." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial cX(n);

        conv(cX,c);

        if ( !BCHCodeBase::round(cX) ) {
            return false;
        }

        conv(c,cX);

        return true;
    }

    /**
     * @brief
     *            Attaches check bits to a message vector such
     *            that the resulting code vector can be transferred
     *            through a noisy channel.
     *
     * @details
     *            The encoding procedure implemented is systematic in the
     *            sense that the last
     *            <i>k=</i>\link getDimension()\endlink coefficients of
     *            the code vector form the coefficients of the input
     *            message vector.
     *
     * @param c
     *            On input, a vector of length
     *            \link getDimension()\endlink; on output, a
     *            vector of length
     *            <i>n=</i>\link getBlockLength()\endlink such that its
     *            last <i>k</i> coefficients form the input message
     *            vector.
     *
     * @warning
     *            If the length of the message vector is different from
     *            <i>k=</i>\link getDimension()\endlink,
     *            an error message is printed to <code>stderr</code> and
     *            the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @see BCHCode::decode()
     * @see BCHCode::round()
     * @see BCHCodeBase::encode()
     */
    void BCHCode::encode( BinaryVector & c ) const {

    	mul(c,this->G,c);
    }

    /**
     * @brief
     *            Rounds the input vector to its nearest codeword
     *            and, if successful, removes the check bits thereby
     *            obtaining the code vector's message polynomial.
     *
     * @details
     *            The implementation of this function is wrapped around
     *            the \link BCHCodeBase::decode()\endlink function.
     *
     * @param m
     *            On input, a vector of length
     *            <i>n=</i>\link getBlockLength()\endlink. On successful
     *            decoding, a vector of length
     *            <i>k=</i>\link getDimension()\endlink; otherwise, if
     *            decoding failed, the content of <i>m</i> is left
     *            unchanged.
     *
     * @return
     *            <code>true</code> if decoding was successful; otherwise,
     *            if the decoding attempt failed, the function returns
     *            <code>false</code>.
     *
     * @warning
     *            If <i>m</i> is of length different from
     *            <i>n=</i>\link getBlockLength()\endlink, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @see BCHCode::encode()
     * @see BCHCode::round()
     * @see BCHCodeBase::decode()
     */
    bool BCHCode::decode( BinaryVector & m ) const {

    	if ( m.getLength() != getBlockLength() ) {
    		cerr << "BCHCode::decode: bad argument, invalid block length"
    			 << endl;
    		exit(EXIT_FAILURE);
    	}

    	BinaryPolynomial tmp;
    	conv(tmp,m);

    	if ( !BCHCodeBase::decode(tmp) ) {
    		return false;
    	}

    	conv(m,tmp);

    	return true;
    }

    /**
     * @brief
     *            Creates a random code word.
     *
     * @param c
     *            On output, a random code word.
     *
     * @param tryRandom
     *            If <code>true</code>, the method uses a cryptographic
     *           number generator if available on the system; otherwise,
     *           the method wraps around the standard <code>rand()</code>
     *           function.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BCHCode::random( BinaryVector & c , bool tryRandom ) const {

    	BinaryVector m(getDimension());
    	m.random(tryRandom);
    	mul(c,this->G,m);
    }

    /**
     * @brief
     *           Assignment operator.
     *
     * @details
     *           Sets this instance to a copy of the specified
     *           BCH code for binary polynomials.
     *
     * @param code
     *           BCH code for binary polynomials
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BCHCode & BCHCode::operator=( const BCHCodeBase & code ) {

    	// Invoke the assignment operator of the 'BCHCodeBase' class ...
    	BCHCodeBase::operator=(code);

    	// ... and compute the generator and check matrix for this
    	// BCH code.
    	this->compute_generator_and_check_matrices();

    	return *this;
    }

    /**
     * @brief
     *           Assignment operator.
     *
     * @param code
     *           The BCH code assigned to this instance.
     *
     * @return
     *           A reference to this BCH code (after assignment).
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BCHCode & BCHCode::operator=( const BCHCode & code ) {

    	// Invoke the assignment operator of the 'BCHCodeBase' class ...
    	BCHCodeBase::operator=(code);

    	// ... and assign the matrices.
    	this->G = code.G;
    	this->H = code.H;

    	return *this;
    }
}
