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
 * @file BinaryVector.h
 *
 * @brief
 *            Provides a class for representing and computing
 *            with binary vectors of fixed but arbitrary degree.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_BINARYVECTOR
#define THIMBLE_BINARYVECTOR

#include <stdint.h>
#include <fstream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Instances of this class represent vectors with
     *            binary/bit entries.
     */
    class THIMBLE_DLL BinaryVector {

    private:

        /**
         * @brief
         *            Encodes the binary entries of the vector.
         *
         * @see getData()
         */
        uint32_t *data;

        /**
         * @brief
         *            The dimension of the vector which should by kept
         *            greater than or equals 0.
         */
        int length;

        /**
         * @brief
         *            The number of 32-bit words that the array
         *            \link data\endlink can hold.
         *
         * @details
         *            The representation of the binary vector is encoded
         *            by 32 bit blocks. To represent a vector of dimension
         *            <i>n</i> at least \f$\lceil n/32\rceil\f$ integers
         *            of width 32 bits are needed. Thus, the field
         *            should be at least \f$\lceil n/32\rceil\f$.
         */
        int numWords;


    public:

        /**
         * @brief
         *            Standard constructor.
         *
         * @details
         *            Constructs a vector of specified length.
         *
         * @param n
         *            The length of this vector.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        BinaryVector( int n = 0 );

        /**
         * @brief
         *           Copy constructor.
         *
         * @param v
         *           The binary vector of which a copy is constructed.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        BinaryVector( const BinaryVector & v );

        /**
         * @brief
         *            Destructor.
         *
         * @details
         *            Frees the data needed to hold the vector's entries.
         */
        ~BinaryVector();

        /**
         * @brief
         *            Assignment operator (procedural version).
         *
         * @param v
         *            The binary vector assigned to this vector.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void assign( const BinaryVector & v );

        /**
         * @brief
         *            Assignment operator (procedural version).
         *
         * @param v
         *            The binary vector assigned to this vector.
         *
         * @return
         *            A reference to this instance (after assignment).
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryVector &operator=( const BinaryVector & v ) {
            assign(v);
            return *this;
        }

        /**
         * @brief
         *            Compares two binary vectors on equality.
         *
         * @param v
         *            The vector with which this vector is compared.
         *
         * @return
         *            <code>true</code> if this vector equals <i>v</i>;
         *            otherwise, <code>false</code>.
         *
         * @attention
         *            Note that if two vectors are of different length,
         *            they are considered to be different causing
         *            the function to return <code>false</code>.
         */
        bool operator==( const BinaryVector & v ) const;

        /**
         * @brief
         *            Compares two binary vectors on inequality.
         *
         * @param v
         *            The vector with which this vector is compared.
         *
         * @return
         *            <code>false</code> if this vector equals <i>v</i>;
         *            otherwise, <code>true</code>.
         *
         * @attention
         *            Note that if two vectors are of different length,
         *            they are considered to be different causing
         *            the function to return <code>true</code>.
         */
        inline bool operator!=( const BinaryVector & v ) const {
        	return !(*this==v);
        }

        /**
         * @brief
         *            Ensures that this instance is able to represent
         *            a vector of specified length without requiring
         *            reallocation.
         *
         * @param newLength
         *            The ensured length for this vector.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void reserve( int newLength );

        /**
         * @brief
         *            Sets the dimension of the vector to the specified
         *            length.
         *
         * @details
         *            If the length of the vector increases, the entries
         *            of the vector up to the old length remain unchanged
         *            and the new entries are padded with zeros. If the
         *            length of the vector is decreased, the first
         *            <code>length</code> entries are kept.
         *
         * @param length
         *            The new length of the vector.
         *
         * @warning
         *            If length is negative, an error message is printed
         *            to <code>stderr</code> and the program exits with
         *            status 'EXIT_STATUS'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void setLength( int length );

        /**
         * @brief
         *            Accesses the dimension of the vector.
         *
         * @return
         *            The dimension of the vector.
         */
        inline int getLength() const {

            return this->length;
        }

        /**
         * @brief
         *            Access the value at the specified position of
         *            the vector.
         *
         * @param j
         *            The position ranging from
         *            0,...,\link getLength()\endlink-1.
         *
         * @return
         *            <code>true</code> if the <i>j</i>th position
         *            equals 1 and, otherwise, if the <i>j</i>th
         *            position equals 0, <code>false</code>.
         *
         * @warning
         *            If <i>j</i> is negative or if <i>j</i> is greater
         *            than or equals \link getLength()\endlink, an
         *            error message is printed to <code>stderr</code>
         *            and the program exits with status 'EXIT_FAILURE'
         */
        bool getAt( int j ) const;

        /**
         * @brief
         *            Ensures that the value at the specified position
         *            equals 1.
         *
         * @param j
         *            The position ranging from
         *            0,...,\link getLength()\endlink-1.
         *
         * @warning
         *            If <i>j</i> is negative or if <i>j</i> is greater
         *            than or equals \link getLength()\endlink, an
         *            error message is printed to <code>stderr</code>
         *            and the program exits with status 'EXIT_FAILURE'
         */
        void setAt( int j );

        /**
         * @brief
         *            Ensures that the value at the specified position
         *            equals 0.
         *
         * @param j
         *            The position ranging from
         *            0,...,\link getLength()\endlink-1.
         *
         * @warning
         *            If <i>j</i> is negative or if <i>j</i> is greater
         *            than or equals \link getLength()\endlink, an
         *            error message is printed to <code>stderr</code>
         *            and the program exits with status 'EXIT_FAILURE'
         */
        void clearAt( int j );

        /**
         * @brief
         *            Set the value of the vector at a specified
         *            position to a specified value.
         *
         * @param j
         *            The position ranging from
         *            0,...,\link getLength()\endlink-1.
         *
         * @param c
         *            <code>true</code> if the value at <i>j</i> is set
         *            to 1 and <code>false</code> otherwise.
         *
         * @warning
         *            If <i>j</i> is negative or if <i>j</i> is greater
         *            than or equals \link getLength()\endlink, an
         *            error message is printed to <code>stderr</code>
         *            and the program exits with status 'EXIT_FAILURE'
         */
        void setAt( int j , bool c );

        /**
         * @brief
         *            Tests whether the binary vector has only entries equals
         *            0.
         *
         * @return
         *            <code>true</code> if all entries of the vector are 0;
         *            otherwise, if there exists a non-zero entry, the
         *            function returns <code>false</code>.
         */
        bool isZero() const;

        /**
         * @brief
         *            Initialize all entries of this vector with 0.
         */
        void setZero();

        /**
         * @brief
         *           Exchanges the entries of the vector at the two
         *           specified positions.
         *
         * @param j0
         *           First position.
         *
         * @param j1
         *           Second position.
         *
         * @warning
         *           If <i>j0</i> or <i>j1</i> is smaller than zero or
         *           greater than or equals the
         *           \link getLength() length\endlink of the vector, an
         *           error message is printed to <code>stderr</code> and
         *           the program exits with status 'EXIT_FAILURE'.
         */
        void exchange( int j0 , int j1 );

        /**
         * @brief
         *            Initialize the entries of this vector with random
         *            values.
         *
         * @details
         *            Each \link getLength()\endlink entry of the vector is
         *            overwritten with random values.
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
         *           Replaces the entries of this vector such that it becomes
         *           random with the specified Hamming weight.
         *
         * @param hammingWeight
         *           The Hamming weight of this vector after the method has
         *           finished.
         *
         * @param tryRandom
         *           If <code>true</code>, the method is advised to use a
         *           cryptographic number generator; otherwise, the method
         *           wraps around the standard <code>rand()</code> function.
         *
         * @warning
         *           If <code>hammingWeight</code> is negative or greater than
         *           or equals \link getLength()\endlink, an error message is
         *           printed to <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        void wrandom( int hammingWeight , bool tryRandom = false );

        /**
         * @brief
         *           Returns the Hamming weight of the vector.
         *
         * @details
         *           Write for the vector
         *           \f[
         *            v = (v_0,...,v_{n-1})
         *           \f]
         *           where \f$v_j\in\{0,1\}\f$. The <i>Hamming weight</i>
         *           of the vector <i>v</i> is defined as the number of
         *           positions <i>j</i> where \f$v_j\neq 0\f$ which, for
         *           binary vectors, is equivalent to the number of
         *           positions <i>j</i> where \f$v_j=1\f$.
         *
         * @return
         *           The Hamming weight of the vector.
         */
        int hammingWeight() const;

        /**
         * @brief
         *           Initializes the first <i>k</i> of the vector as 1s and
         *           the remaining as zeros.
         *
         * @details
         *           Using the method \link bnext()\endlink, we can systematically
         *           iterate through all vectors of
         *           \link hammingWeight() Hamming weight\endlink exactly equals
         *           <i>k</i>. This method sets this vectors to the first vector
         *           of this sequence, i.e.,
         *           \f[
         *            v=(\underbrace{1,...,1}_{k~\mbox{times}}~,~
         *               \underbrace{0,...,0}_{n-k~\mbox{times}}).
         *           \f]
         *
         * @param k
         *           The number of positions of the vector being set to 1.
         *
         * @see bnext()
         *
         * @warning
         *           If <i>k</i> is negative or greater than
         *           \link getLength()\endlink, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILIRE'.
         */
        void binit( int k );

        /**
         * @brief
         *           Iterates to the ''next'' binary vector of the
         *           same Hamming weight.
         *
         * @details
         *           Assume the length of the vector is <i>n</i> and its
         *           Hamming weight is <i>k</i>. The function selects a
         *           ''next'' vector in a sequence of all \f${k\choose n}\f$
         *           combinations of <i>n</i>-length vectors with Hamming
         *           weight <i>k</i> starting from
         *           \f[
         * (\underbrace{1,...,1}_{k~\mbox{times}}~,~
         *  \underbrace{0,...,0}_{n-k~\mbox{times}})
         *           \f]
         *           to
         *           \f[
         * (\underbrace{0,...,0}_{n-k~\mbox{times}}~,~
         *  \underbrace{1,...,1}_{k~\mbox{times}}).
         *           \f]
         *
         * @return
         *           <code>true</code> if the content of the vector has been
         *           successfully changed to the ''next'' vector of the
         *           same Hamming weight; otherwise, if the vector is already
         *           the ''last'' vector, the result is
         *           <code>false</code>
         *
         * @see binit()
         * @see next()
         */
        bool bnext();

        /**
         * @brief
         *           Iterates to the ''next'' vector in a sequence of distinct
         *           vectors of fixed length.
         *
         * @details
         *           The function can be used to systematically iterate through
         *           all distinct \f$2^n\f$ vectors starting from
         *           \f$(0,...,0)\f$ to \f$(1,...,0)\f$. The sequence of the
         *           iteration has the property that its vectors are
         *           monotonically increasing in their
         *           \link hammingWeight() Hamming weight\endlink which is,
         *           in particular, useful for solving problems related to
         *           <i>maximum likelihood decoding</i>.
         *
         * @return
         *           <code>true</code> if the content has been successfully
         *           changed to the ''next'' vector; otherwise, if the vector
         *           is already the ''last'' vector, the function returns
         *           <code>false</code>.
         *
         * @see bnext()
         */
        bool next();

        /**
         * @brief
         *           Swaps the content of two binary vectors such the
         *           the one represents the other.
         *
         * @param v
         *           Binary vector being assigned the fields of <i>w</i>.
         *
         * @param w
         *           Binary vector being assigned the fields of <i>v</i>.
         */
        static void swap( BinaryVector & v , BinaryVector & w );

        /**
         * @brief
         *           Performs an addition of two binary vectors.
         *
         * @param c
         *           On output, the sum of the input
         *           <code>a</code> and <code>b</code>.
         *
         * @param a
         *           First summand vector.
         *
         * @param b
         *           Second summand vector.
         *
         * @warning
         *           If the length of <code>a</code> and <code>b</code>
         *           are different, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void add
        ( BinaryVector & c ,
          const BinaryVector & a , const BinaryVector & b );

        /**
         * @brief
         *           Performs a subtraction of two binary vectors.
         *
         * @details
         *           Note that subtraction is equivalent to addition for
         *           binary fields.
         *
         * @param c
         *           On output, the difference between the input
         *           <code>a</code> and <code>b</code>.
         *
         * @param a
         *           Minuend.
         *
         * @param b
         *           Subtrahend.
         *
         * @warning
         *           If the length of <code>a</code> and <code>b</code>
         *           are different, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline static void sub
        ( BinaryVector & c ,
          const BinaryVector & a , const BinaryVector & b ) {

            // Subtraction is equivalent to addition in the binary field
            add(c,a,b);
        }

        /**
         * @brief
         *            Determines the Hamming distance between two binary
         *            vectors.
         *
         * @param a
         *            First binary vector.
         *
         * @param b
         *            Second binary vector.
         *
         * @return
         *            The Hamming distance between <code>a</code> and
         *            <code>b</code>.
         *
         * @warning
         *            If <code>a</code> and <code>b</code> have a different
         *            length, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         */
        static int hammingDistance
        ( const BinaryVector & a , const BinaryVector & b );

        /**
         * @brief
         *            Overloaded '+'-operator to perform an addition
         *            of two binary vectors.
         *
         * @details
         *            The '+'-operator allows expressions as
         *            <pre>
         *             BinaryVector a , b , c;
         *             ...
         *             c = a+b;
         *            </pre>
         *            and is wrapped around the \link add()\endlink
         *            function.
         *
         * @param b
         *            Summand.
         *
         * @return
         *            The sum of this vector and <code>b</code>.
         *
         * @warning
         *            If the length of <code>b</code> is different
         *            from the length of this vector, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryVector operator+( const BinaryVector & b ) const {
            BinaryVector c;
            add(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Overloaded '-'-operator to perform a
         *            subtraction/addition of two binary vectors.
         *
         * @details
         *            The '-'-operator allows expressions as
         *            <pre>
         *             BinaryVector a , b , c;
         *             ...
         *             c = a-b;
         *            </pre>
         *
         * @param b
         *            Subtrahend.
         *
         * @return
         *            The difference between this vector and <code>b</code>.
         *
         * @warning
         *            If the length of <code>b</code> is different
         *            from the length of this vector, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryVector operator-( const BinaryVector & b ) const {
            BinaryVector c;
            sub(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Overloaded '+='-operator to add a vector to this
         *            vector inplace.
         *
         * @details
         *            The '+='-operator allows expressions of the
         *            form
         *            <pre>
         *             BinaryVector a , b;
         *             ...
         *             a += b;
         *            </pre>
         *
         * @param b
         *            Summand.
         *
         * @return
         *            A reference to this vector to which the vector
         *            <code>b</code> has been added.
         *
         * @warning
         *            If the length of <code>b</code> is different
         *            from the length of this vector, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryVector & operator+=( const BinaryVector & b ) {
            add(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Overloaded '-='-operator to subtract a vector from this
         *            vector inplace.
         *
         * @details
         *            The '-='-operator allows expressions of the
         *            form
         *            <pre>
         *             BinaryVector a , b;
         *             ...
         *             a -= b;
         *            </pre>
         *
         * @param b
         *            Subtrahend.
         *
         * @return
         *            A reference to this vector from which the vector
         *            <code>b</code> has been subtracted.
         *
         * @warning
         *            If the length of <code>b</code> is different
         *            from the length of this vector, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryVector & operator-=( const BinaryVector & b ) {
            sub(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Access the data array encodes the binary entries of the
         *            vector.
         *
         * @details
         *            If the vector represented by this instance is written as
         *            \f[
         *             (b_0,...,b_{n-1})
         *            \f]
         *            where \f$b_j\in\{0,1\}\f$, then the <i>i</i>th 32 bit
         *            integer, where
         *            \f$i<\lceil n/32\rceil\f$, equals
         *            \f[
         *             \sum_{k=0,...,61,~i\cdot 32+k<n}b_{i\cdot 32+k};
         *            \f]
         *            and \link getLength()\endlink equals
         *            \f$n\geq 0\f$.
         *
         * @return
         *            The data array encoding the binary entries of the
         *            vector.
         */
        inline const uint32_t *getData() const {

            return this->data;
        }

        /**
         * @brief
         *           Returns the array encoding the coefficients of this
         *           vector.
         *
         * @details
         *           The result of this function is the same as the result of
         *           \link getData()\endlink except that the content of the
         *           vector can be changed by changing the content of the
         *           result.
         *
         * @warning
         *           Use the result with caution and perform only changes if
         *           you really know what your are doing; otherwise, you might
         *           experience undocumented behaviour.
         *
         * @return
         *           Array encoding the coefficients of this polynomial.
         */
        inline uint32_t *getData_nonconst() {

            return this->data;
        }
    };

    /**
     * @brief
     *           Performs an addition of two binary vectors.
     *
     * @param c
     *           On output, the sum of the input
     *           <code>a</code> and <code>b</code>.
     *
     * @param a
     *           First summand vector.
     *
     * @param b
     *           Second summand vector.
     *
     * @warning
     *           If the length of <code>a</code> and <code>b</code>
     *           are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    inline void add
    ( BinaryVector & c ,
      const BinaryVector & a , const BinaryVector & b ) {

        BinaryVector::add(c,a,b);
    }

    /**
     * @brief
     *           Performs a subtraction of two binary vectors.
     *
     * @details
     *           Note that subtraction is equivalent to addition for
     *           binary fields.
     *
     * @param c
     *           On output, the difference between the input
     *           <code>a</code> and <code>b</code>.
     *
     * @param a
     *           Minuend.
     *
     * @param b
     *           Subtrahend.
     *
     * @warning
     *           If the length of <code>a</code> and <code>b</code>
     *           are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    inline void sub
    ( BinaryVector & c ,
      const BinaryVector & a , const BinaryVector & b ) {

        // Subtraction is equivalent to addition in the binary field
        BinaryVector::add(c,a,b);
    }

    /**
     * @brief
     *            Determines the Hamming distance between two binary
     *            vectors.
     *
     * @param a
     *            First binary vector.
     *
     * @param b
     *            Second binary vector.
     *
     * @return
     *            The Hamming distance between <code>a</code> and
     *            <code>b</code>.
     *
     * @warning
     *            If <code>a</code> and <code>b</code> have a different
     *            length, an error message is printed to
     *            <code>stderr</code> and the program exits with status
     *            'EXIT_FAILURE'.
     */
    inline int hammingDistance
    ( const BinaryVector & a , const BinaryVector & b ) {

        return BinaryVector::hammingDistance(a,b);
    }

    /**
     * @brief
     *            Prints a representation of a vector to the specified
     *            output stream.
     *
     * @param out
     *            The output stream.
     *
     * @param v
     *            The binary vector being printed to <code>out</code>.
     *
     * @return
     *            A reference to <code>out</code>.
     */
    THIMBLE_DLL std::ostream & operator<<( std::ostream & out , const BinaryVector & v );
}

#endif /* THIMBLE_BINARYVECTOR */
