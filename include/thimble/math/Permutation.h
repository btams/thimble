/*
 *  THIMBLE --- Research Libary for Development and Analysis of
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
 * @file Permutation.h
 *
 * @brief
 *            Provides a class for representation of and computation
 *            with permutations.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_PERMUTATION_H_
#define THIMBLE_PERMUTATION_H_

#include <fstream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Instances of this class represent permutations on
     *            a finite number of elements.
     */
	class THIMBLE_DLL Permutation {

    private:

        /**
         * @brief
         *            The number of elements on which this permutation
         *            operates.
         */
        int n;

        /**
         * @brief
         *            An array of <i>n</i> distinct integers that
         *            encode to which index an index is reordered by this
         *            permutation.
         *
         * @details
         *            Let \f$\pi:\{0,....,n-1\}\rightarrow\{0,...,n-1\}\f$ be
         *            the permutation represented by this instance. Let
         *            \f$y=0,...,n-1\f$ be the element to which
         *            \f$x=0,...,n-1\f$ is mapped. Then <code>data[x]=y</code>.
         */
        int *data;

        /**
         * @brief
         *            Prepares this permutation to encode a permutation of
         *            dimension <i>n</i>.
         *
         * @details
         *            A call of this function causes 1) the field
         *            \link data\endlink to be freed and 2) to be allocated
         *            via <code>malloc</code> to hold <i>n</i> integers of
         *            type <code>int</code> or set to <code>NULL</code> if
         *            <i>n==0</i>. The content of \link data\endlink is not
         *            set to encode a valid permutation; this must be ensured
         *            by the programmer for which reason this method is
         *            declared private.
         *
         *            Furthermore, the argument <i>n</i> must be greater than
         *            or equals 0; otherwise, the program runs into
         *            undocumented behaviour.
         *
         * @param n
         *            The dimension/number of elements on which this
         *            permutation is prepared to operate.
         */
        void prepareDim( int n );

	public:

        /**
         * @brief
         *            Creates the identity permutation to operate on an
         *            <i>n</i> elements.
         *
         * @param n
         *            The number of elements on which this permutation
         *            operates.
         *
         * @warning
         *            If <i>n</i> is negative, an error message is printed
         *            to <code>stderr</code> and the program exits with
         *            status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
		Permutation( int n = 0 );

        /**
         * @brief
         *            Copy constructor.
         *
         * @details
         *            Creates a copy of the specified permutation.
         *
         * @param P
         *            The permutation of which a copy is created.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        Permutation( const Permutation & P );

        /**
         * @brief
         *            Destructor.
         *
         * @details
         *            Frees the data that has been allocated to encode this
         *            permutation.
         */
		~Permutation();

        /**
         * @brief
         *            Assignment operator.
         *
         * @details
         *            Assigns this permutation to a copy of <i>P</i>.
         *
         * @param P
         *            The permutation of which this instance is assigned a
         *            copy of.
         *
         * @return
         *            A reference to this instance (after assignment).
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        Permutation &operator=( const Permutation & P );

        /**
         * @brief
         *            Sets this instance to represent the identity permutation
         *            operating on <i>n</i> elements.
         *
         * @param n
         *            The number of elements on which this permutation
         *            operates.
         *
         * @warning
         *            If <i>n</i> is negative, an error message is printed
         *            to <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void setDimension( int n );

        /**
         * @brief
         *            Returns the number of elements on which this permutation
         *            operates.
         */
        inline int getDimension() const {
            return this->n;
        }

        /**
         * @brief
         *           Evaluates this permutation on the specified index.
         *
         * @details
         *           Write \f$\pi:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ as
         *           the permutation represented by this instance. The function
         *           returns \f$\pi(x)\f$.
         *
         * @param x
         *           The index at which this permutation is evaluated.
         *
         * @return
         *           The evaluation of this permutation at <i>x</i>.
         *
         * @warning
         *           If <i>x</i> is negative or greater than or equals
         *           \link getDimension()\endlink, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         */
		int eval( int x ) const;

        /**
         * @brief
         *           Exchanges the evaluations of this permutation for the
         *           specified arguments.
         *
         * @details
         *           Let \f$\pi:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ be
         *           the permutation represented by this instance. The method
         *           replaces this permutation by the permutation
         *           \f$\pi':\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ where
         *           \f$\pi'(x_0)=\pi(x_1)\f$, \f$\pi'(x_1)=\pi(x_0)\f$ and
         *           if \f$x\neq x_0,x_1\f$, then \f$\pi'(x)=\pi(x)\f$.
         *
         * @param x0
         *           First argument.
         *
         * @param x1
         *           Second argument.
         *
         * @warning
         *           If one of the arguments is negative or greater than or
         *           equals \link getDimension()\endlink, an error message
         *           is printed to <code>stderr</code> and the program exits
         *           with status 'EXIT_FAILURE'.
         */
        void exchange( int x0 , int x1 );

        /**
         * @brief
         *           Replaces this permutation by a random permutation of
         *           operating on the same number of elements.
         *
         * @param tryRandom
         *           If <code>true</code>, the method uses a cryptographic
         *           number generator if available on the system; otherwise,
         *           the method wraps around the standard <code>rand()</code>
         *           function.
         */
        void random( bool tryRandom = false );

        /**
         * @brief
         *           Swaps the content of two permutations such that the one
         *           encodes the other.
         *
         * @param P
         *           First permutation.
         *
         * @param Q
         *           Second permutation.
         */
        static void swap( Permutation & P , Permutation & Q );

        /**
         * @brief
         *           Computes the concatenation of two permutations.
         *
         * @details
         *           More precisely, let
         *           \f$P,Q:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ be the
         *           two permutations of which the concatenation is to be
         *           computed. The method computes \f$R=P\circ Q\f$ which is
         *           the permutation that maps \f$x\f$ to \f$P(Q(x))\f$.
         *
         *           Consider the permutations \f$P\f$ and \f$Q\f$ as
         *           the matrices
         *           \f[
         *            A=(a_{i,j})_{i,j=0,...,n-1}
         *           \f]
         *           and
         *           \f[
         *            B=(b_{i,j})_{i,j=0,...,n-1},
         *           \f]
         *           respectively, where \f$a_{i,j}=1\f$ if \f$P(i)=j\f$ and
         *           \f$a_{i,j}=0\f$ otherwise and \f$b_{i,j}=1\f$ if
         *           \f$Q(i)=j\f$ and \f$b_{i,j}=0\f$ otherwise. The the
         *           concatentation of two permutations can be interpreted
         *           as the matrix multiplication
         *           \f[
         *            C=A\cdot B=(c_{i,j})
         *           \f]
         *           where \f$c_{i,j}=1\f$ if \f$R(i)=P(Q(i))=j\f$ and
         *           \f$c_{i,j}=0\f$ otherwise. Therfore, the name of this
         *           method <code>mul</code> suggests that the concatenation
         *           operation is a multiplication.
         *
         * @param R
         *           On output, the concatentation of <i>P</i> and <i>Q</i>.
         *
         * @param P
         *           Outer permutation.
         *
         * @param Q
         *           Inner permutation.
         *
         * @warning
         *           If <i>P</i> and <i>Q</i> operate on a different number
         *           of elements, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( Permutation & R ,
          const Permutation & P , const Permutation & Q );

        /**
         * @brief
         *           Computes the inverse of a permutation.
         *
         * @details
         *           The inverse of a permutation
         *           \f$P:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ is the
         *           unique \f$R:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$
         *           such that \f$R(P(i))=P(R(i))=i\f$.
         *
         * @param R
         *           On output, the inverse of the permutation <i>P</i>.
         *
         * @param P
         *           The permutation of which the inverse is computed.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void inv
        ( Permutation & R , const Permutation & P );

        /**
         * @brief
         *          Overloaded '*'-operator to allow two permutations
         *          to be concatenated.
         *
         * @details
         *          The '*'-operator allows two permutations to be
         *          concatenated through expressions as
         *          <pre>
         *           Permutation P , Q;
         *           ...
         *           Permutation R = P * Q;
         *          </pre>
         *
         * @param Q
         *          Inner permutation.
         *
         * @return
         *          The concatentation of this permutation as the outer
         *          permutation and <i>Q</i> as the inner permutation, i.e.,
         *          \f$P\circ Q=P(Q(\cdot))\f$.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline Permutation operator*( const Permutation & Q ) {

            Permutation R;
            mul(R,*this,Q);
            return R;
        }
    };

    /**
     * @brief
     *           Swaps the content of two permutations such that the one
     *           encodes the other.
     *
     * @param P
     *           First permutation.
     *
     * @param Q
     *           Second permutation.
     */
    inline void swap( Permutation & P , Permutation & Q ) {
        Permutation::swap(P,Q);
    }

    /**
     * @brief
     *           Computes the concatenation of two permutations.
     *
     * @details
     *           More precisely, let
     *           \f$P,Q:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ be the
     *           two permutations of which the concatenation is to be
     *           computed. The method computes \f$R=P\circ Q\f$ which is
     *           the permutation that maps \f$x\f$ to \f$P(Q(x))\f$.
     *
     *           Consider the permutations \f$P\f$ and \f$Q\f$ as
     *           the matrices
     *           \f[
     *            A=(a_{i,j})_{i,j=0,...,n-1}
     *           \f]
     *           and
     *           \f[
     *            B=(b_{i,j})_{i,j=0,...,n-1},
     *           \f]
     *           respectively, where \f$a_{i,j}=1\f$ if \f$P(i)=j\f$ and
     *           \f$a_{i,j}=0\f$ otherwise and \f$b_{i,j}=1\f$ if
     *           \f$Q(i)=j\f$ and \f$b_{i,j}=0\f$ otherwise. The the
     *           concatentation of two permutations can be interpreted
     *           as the matrix multiplication
     *           \f[
     *            C=A\cdot B=(c_{i,j})
     *           \f]
     *           where \f$c_{i,j}=1\f$ if \f$R(i)=P(Q(i))=j\f$ and
     *           \f$c_{i,j}=0\f$ otherwise. Therfore, the name of this
     *           method <code>mul</code> suggests that the concatenation
     *           operation is a multiplication.
     *
     * @param R
     *           On output, the concatentation of <i>P</i> and <i>Q</i>.
     *
     * @param P
     *           Outer permutation.
     *
     * @param Q
     *           Inner permutation.
     *
     * @warning
     *           If <i>P</i> and <i>Q</i> operate on a different number
     *           of elements, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXUIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( Permutation & R ,
      const Permutation & P , const Permutation & Q ) {
        Permutation::mul(R,P,Q);
    }

    /**
     * @brief
     *           Computes the inverse of a permutation.
     *
     * @details
     *           The inverse of a permutation
     *           \f$P:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ is the
     *           unique \f$R:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$
     *           such that \f$R(P(i))=P(R(i))=i\f$.
     *
     * @param R
     *           On output, the inverse of the permutation <i>P</i>.
     *
     * @param P
     *           The permutation of which the inverse is computed.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void inv
    ( Permutation & R , const Permutation & P ) {
        Permutation::inv(R,P);
    }

    /**
     * @brief
     *           Computes the inverse of a permutation.
     *
     * @details
     *           The inverse of a permutation
     *           \f$P:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$ is the
     *           unique \f$R:\{0,...,n-1\}\rightarrow\{0,...,n-1\}\f$
     *           such that \f$R(P(i))=P(R(i))=i\f$.
     *
     * @param P
     *           The permutation of which the inverse is computed.
     *
     * @return
     *           The inverse of the permutation <i>P</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline Permutation inv( const Permutation & P ) {
        Permutation R;
        Permutation::inv(R,P);
        return R;
    }

    /**
     * @brief
     *           Prints a text representation of a permutation to the
     *           specified output stream.
     *
     * @details
     *           The function prints the following
     *           <pre>
     * [<P(0)> , <P(1)> , ... , <P(n-1)>]
     *           </pre>
     *           where <code><P(i)></code> are replacements for the integers
     *           <code>P.eval(i)</code>. For example, the identity of the
     *           permutation operating on 5 elements is written as
     *           <pre>
     * [0 , 1 , 2 , 3 , 4]
     *           </pre>
     *
     * @param out
     *           The output stream to which a text representation of the
     *           permutation <i>P</i> is written.
     *
     * @param P
     *           The permutation of which a text representation is written
     *           to <code>out</code>.
     *
     * @return
     *           A reference to <code>out</code> after the text representation
     *           has been written.
     */
    THIMBLE_DLL std::ostream & operator<<( std::ostream & out , const Permutation & P );
}


#endif /* THIMBLE_PERMUTATION_H_ */
