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
 * @file BCHCode.h
 *
 * @brief
 *            Provides a class for generating binary BCH codes that can
 *            be used for bit error-correction and that provides extended
 *            functionalities as compared to the basic BCH code class
 *            \link thimble::BCHCodeBase BCHCodeBase\endlink.
 * @author Benjamin Tams
 */

#ifndef THIMBLE_BCHCODE_H
#define THIMBLE_BCHCODE_H

#include <thimble/dllcompat.h>
#include <thimble/math/linalg/BinaryVector.h>
#include <thimble/math/linalg/BinaryMatrix.h>
#include <thimble/ecc/BCHCodeBase.h>

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of this class represent binary BCH codes enabling
	 *            a bit error-correcting scheme for binary vectors.
	 *
	 * @details
	 *            The use of this class should be preferred as compared
	 *            to the use of the class \link BCHCodeBase\endlink from
	 *            which all public member functions are inherited. By
	 *            theory BCH codes work with polynomials encoding vectors.
	 *            This class wraps functions for
	 *            \link BinaryVector binary vectors\endlink around the
	 *            functions for
	 *            \link BinaryPolynomial binary polynomials\endlink
	 *            provided by \link BCHCodeBase\endlink.
	 *
	 * @see @see <a href="http://en.wikipedia.org/wiki/BCH_code" target="_blank">http://en.wikipedia.org/wiki/BCH_code</a>
	 */
    class THIMBLE_DLL BCHCode : public BCHCodeBase {

    private:

    	/**
    	 * @brief
    	 *            A generator matrix of this BCH code.
    	 *
    	 * @details
    	 *            This member is equals a generator matrix of
    	 *            this <i>(n,k)</i>-BCH code <i>C</i> in such that
    	 *            \f[
    	 *             C = \{~Gm~|~m\in\{0,1\}^k~\}
    	 *            \f]
    	 *            where <i>n</i> is the
    	 *            \link getBlockLength() block length\endlink and <i>k</i>
    	 *            the \link getDimension() dimension\endlink of <i>C</i>.
    	 *            <br><br>
    	 *            Furthermore, the matrix <i>G</i> is such that it
    	 *            represents a systematic generator matrix which means that
    	 *            the last <i>k</i> bits of <i>Gm</i> agree with the
    	 *            <i>k</i> bits of <i>m</i>.
    	 */
    	BinaryMatrix G;

    	/**
    	 * @brief
    	 *            A check matrix of this BCH code.
    	 *
    	 * @details
    	 *            This member represents a check matrix of this
    	 *            <i>(n,k)</i>-BCH code <i>C</i> in such that for every
    	 *            \f$c\in C\f$ the relation \f$H\cdot c=0\f$ holds.
    	 *            Specifically, <i>H</i> is an \f$(n-k)\times n\f$ matrix
    	 *            of full rank such that \f$H\cdot G=0\f$.
    	 */
        BinaryMatrix H;

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
        void compute_generator_and_check_matrices();

    public:

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
        BCHCode( int n , int numErrors );

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
        BCHCode( const BCHCodeBase & code );

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
        BCHCode( const BCHCode & code );

        /**
         * @brief
         *            Destructor.
         *
         * @details
         *            Ensures that all data helt by this instance is
         *            freed.
         */
        ~BCHCode();

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
        bool round( BinaryVector & c ) const;

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
        void encode( BinaryVector & c ) const;

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
        bool decode( BinaryVector & m ) const;

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
        void random( BinaryVector & c , bool tryRandom = false ) const;

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
        BCHCode & operator=( const BCHCodeBase & code );

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
        BCHCode & operator=( const BCHCode & code );


        /**
         * @brief
         *           Access the binary generator matrix of this
         *           BCH code.
         *
         * @return
         *           A constant reference to the generator matrix of
         *           this BCH code.
         */
        inline const BinaryMatrix & getGeneratorMatrix() const {

            return this->G;
        }


        /**
         * @brief
         *           Access the binary check matrix of this BCH code.
         *
         * @return
         *           A constant reference to the check matrix of this
         *           BCH code.
         */
        inline const BinaryMatrix & getCheckMatrix() const {

            return this->H;
        }

    };
}

#endif /* THIMBLE_BCHCODE_H */
