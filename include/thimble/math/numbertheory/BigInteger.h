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
 * @file BigInteger.h
 *
 * @brief
 *            Provides functionalities for arithmetic of integers of arbitrary
 *            length/precision.
 *
 * @author Benjamin Tams
 */
#ifndef THIMBLE_BIGINTEGER_H_
#define THIMBLE_BIGINTEGER_H_

#include <stdint.h>
#include <iostream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Class that represents integers of arbitrary length/precision.
	 *
	 * @details
	 *            <h3>Creating an Instance of a Big Integer</h3>
	 *            An integer representing 0 can be constructed via
	 *            <pre>
	 *             BigInteger a;
	 *            </pre>
	 *            on initialization, we can also specify the value of the
	 *            integer by a primitive integer in the range between
	 *            -4294967295...4294967295. For example, if we want to create
	 *            an integer initialized by -1 we may run
	 *            <pre>
	 *             BigInteger a = -1;
	 *            </pre>
	 *            or alternatively
	 *            <pre>
	 *             BigInteger a(-1);
	 *            </pre>
	 *            furthermore, after an integer has been created we can
	 *            initialize it by an integer in the range
	 *            -4294967295...4294967295 by simply letting
	 *            <pre>
	 *             a = -1;
	 *            </pre>
	 *            a big integer can be printed to <code>stdout</code> in decimal
	 *            representation via
	 *            <pre>
	 *             cout << a << endl;
	 *            </pre>
	 *            <h3>Arithmetic</h3>
	 *            Addition, subtraction, and multiplication of two big integers
	 *            can be computed via the functions
	 *            <code>
	 *      \link add(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>,
	 *            <code>
	 *      \link sub(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>, and
	 *            <code>
	 *      \link mul(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>, respectively. For example, if two big integers
	 *            <pre>
	 *             BigInteger a , b;
	 *            </pre>
	 *            its sum can be computed via
	 *            <pre>
	 *             BigInteger c;
	 *             add(c,a,b);
	 *            </pre>
	 *            furthermore, it is possible to store the sum of <i>a</i> or
	 *            <i>b</i> in one of those instances by dismissing the old
	 *            values, i.e. after
	 *            <pre>
	 *             add(a,a,b);
	 *            </pre>
	 *            the big integer <i>a</i> will contain the sum of the old
	 *            <i>a</i> and <i>b</i>. Similarly, after
	 *            <pre>
	 *             add(a,a,a);
	 *            </pre>
	 *            the <i>a</i> will be the double value of the former
	 *            <i>a</i>.
	 *            <br><br>
	 *            Subtraction and multiplication work similar.
	 *            <br><br>
	 *            The functions that are related with Euclidean division,
	 *            which are
	 *            <code>
	 *  \link div(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>,
	 *            <code>
	 *  \link rem(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>, and
	 *            <code>
     * \link divRem(BigInteger&,BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>,
	 *            require from the denominator to be non-zero. For example,
	 *            to compute the quotient <i>q</i> of two big integers
	 *            <i>a</i> and <i>b</i>, we may run
	 *            <pre>
	 *             BigInteger q;
	 *             div(q,a,b);
	 *            </pre>
	 *            More precisely, the <i>Euclidean quotient</i> is performed
	 *            which is defined to be such that there exists a positive
	 *            integer <i>r</i> with <i>0<=|r|<|b|</i> and
	 *            \f[
	 *             q\cdot b+r=a;
	 *            \f]
	 *            <i>r</i> is called the remainder which can be computed via
	 *            <pre>
	 *             BigInteger r;
	 *             rem(r,a,b);
	 *            </pre>
	 *            quotient and remainder can both be computed at the same
	 *            time via
	 *            <pre>
	 *             BigInteger q , r;
	 *             divRem(q,r,a,b);
	 *            </pre>
	 *            if in any of the above cases <i>b</i> is zero, i.e. if
	 *            <code>b.isZero()</code> is <code>true</code>, an error
	 *            message will be printed to <code>stderr</code> and causes an
	 *            exit with status 'EXIT_FAILURE'.
	 */
	class THIMBLE_DLL BigInteger {

	private:

		/**
		 * @brief
		 *            Contains the data of the big integer.
		 *
		 * @details
		 *            At any time, the array is able to hold
		 *            <code>capacity</code> unsigned integer of 32 bit length.
		 *            But, only <code>|size|</code> of them are relevant.
		 *            More precisely, the integer that is represented will
		 *            be of absolute value
		 *            \f[
		 *             \sum_{i=0}^{|size|-1}data[i]\cdot 2^{i\cdot 32}
		 *            \f]
		 *            where the integer is positive if <code>size>=0</code>
		 *            and negative if <code>size<0</code>.
		 */
		uint32_t *data;

		/**
		 * @brief
		 *            Specifies the number of relevant <code>uint32_t</code>s
		 *            stored in <code>data</code> as well as the sign of the
		 *            integer.
		 *
		 * @details
		 *            See the <code>\link data\endlink</code> member for
		 *            details.
		 */
		int size;

		/**
		 * @brief
		 *            Specifies how many <code>uint32_t</code>s the array
		 *            <code>\link data\endlink</code> can hold.
		 *
		 * @details
		 *            During computation the value of the instances
		 *            representing big integers change frequently and may
		 *            require the arra <code>\link data\endlink</code> to hold
		 *            a varying number of <code>uint32_t</code>s. This value,
		 *            which always must be larger than or equals
		 *            <code>|size|</code>, indicates how many elements the
		 *            array <code>\link data\endlink</code> can store. If more
		 *            than <code>capacity</code> elements are required, this
		 *            will cause an reallocation of
		 *            <code>\link data\endlink</code>.
		 */
		int capacity;

		/**
		 * @brief
		 *            Normalizes the member <code>\link size\endlink</code> to
		 *            be of minimal absolute value such that the same integer
		 *            is represented by this instances.
		 *
		 * @details
		 *            During computation leading elements in
		 *            <code>\link data\endlink</code> may vanish to 0. In the
		 *            case <code>\link size\endlink</code> still may indicate
		 *            them as to be relevant even if they are acutally not.
		 *            This method makes <code>\link size\endlink</code> such
		 *            that the same integer is represented but such that
		 *            <code>data[|size|-1]</code> is non-zero in case this
		 *            integer is non-zero.
		 */
		void normalize();

	public:

		/**
		 * @brief
		 *            Creates an instance of an big integer intialized by
		 *            the specified value.
		 *
		 * @details
		 *            The specified value must be in the range
		 *            -4294967295...4294967295 even though <i>v</i> can exceed
		 *            the range. If no value is specified, the specified value
		 *            is assumed to be zero. In this way, this constructor
		 *            also serves as the standard constructor allowing
		 *            statements as
		 *            <pre>
		 *             BigInteger a;
		 *            </pre>
		 *
		 * @param v
		 *            Initial value in the range -4294967295...4294967295.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated or if
		 *            <i>v</i> is outside the range
		 *            -4294967295...4294967295 the constructor prints an error
		 *            message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		BigInteger( int64_t v = 0 );

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @details
		 *            The created instance will be the same integer as
		 *            represented by <i>a</i>.
		 *
		 * @param a
		 *            Integer that is copied.
		 */
		BigInteger( const BigInteger & a );

		/**
		 * @brief
		 *            Destructor which frees all the memory that is used
		 *            to represent the integer.
		 */
		~BigInteger();

		/**
		 * @brief
		 *            Returns the number of relevant unsigned 32-bit integers
		 *            needed to represent the integer.
		 *
		 * @details
		 *            More precisely, if the represented integer is
		 *            \f[
		 *             \pm\sum_{i=0}^{n-1}c_i\cdot 2^{i\cdot 32}
		 *            \f]
		 *            with \f$c_i\neq 0\f$ then the result is <i>n</i>.
		 *            If the integer is 0, the result will be 1.
		 *
		 * @return
		 *            Number of relevant unsigned 32-bit integers that
		 *            represents the integer.
		 */
		inline int getNumWords() const {
			return this->size<0?-this->size:this->size;
		}

		/**
		 * @brief
		 *            Returns the number of bits needed to represent the absolute
		 *            value of this integer.
		 *
		 * @details
		 *            More precisely, if the absolute value of this integer is
		 *            <i>a</i> then this function computes the minimal <i>l</i>
		 *            such that \f$2^l>a\f$.
		 *
		 * @return
		 *            Number of bits to represent the absolute value of this
		 *            integer.
		 */
		int numBits() const;

		/**
		 * @brief
		 *            Approximates the binary logarithm of this big integer.
		 *
		 * @return
		 *            Approximation of the binary logarithm of this big
		 *            integer.
		 */
		float log2() const;

		/**
		 * @brief
		 *            Ensures that the integer is capable to hold at least
		 *            <code>numWords</code> unsigned 32-bit integers.
		 *
		 * @details
		 *            After invocation, using the big integer will not cause
		 *            an reallocation, which is imperformant, unless it is
		 *            required to store more than <code>numWords</code>
		 *            unsigned 32-bit integers.
		 *            <br><br>
		 *            Another side effect of the method is that data
		 *            that this integer is capable to store but that is not
		 *            relevant to represent the integer is overwritten by
		 *            zero bytes.
		 *
		 * @param numWords
		 *            Number of unsigned 32-bit integers that this instance
		 *            is guaranteed to hold at least.
		 *
		 * @warning
		 *            If <code>numWords</code> is negative an error message
		 *            is printed to <code>stderr</code> and an exit with
		 *            status 'EXIT_FAILURE' is caused.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated on a
		 *            possible reallocation, the method prints an error
		 *            message to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		void ensureCapacity( int numWords );

		/**
		 * @brief
		 *            Assign the big integer by the integer represented
		 *            by <i>a</i>.
		 *
		 * @details
		 *            After assignment, the integer represented by this
		 *            instance will be the same as represented by <i>a</i>.
		 *
		 * @param a
		 *            Big integer whith that this instance is assigned.
		 *
		 * @return
		 *            A reference to this instance.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated the
		 *            an error message is printed to <code>stderr</code>
		 *            and an exit with status 'EXIT_FAILURE' is caused.
		 */
		BigInteger &operator=( const BigInteger & a );

		/**
		 * @brief
		 *            Swaps the (private) member variables of this instance
		 *            with the member variables of <i>a</i>.
		 *
		 * @details
		 *            As a consequence, after calling this method, the integer
		 *            represented by this instance will be represented by
		 *            <i>a</i> and vice versa.
		 */
		void swap( BigInteger & a );

		/**
		 * @brief
		 *            Determine the sign of the big integer.
		 *
		 * @details
		 *            If the integer represented is negative, the result will
		 *            be -1; if the integer is positive, the result will be
		 *            1; otherwise, if the integer is zero, the result will be
		 *            0.
		 *
		 * @return
		 *            Sign of the integer.
		 */
		inline int sign() const {
			return this->size<0?-1:(this->size==1&&this->data[0]==0?0:1);
		}

		/**
		 * @brief
		 *            Converts this integer to an integer of primitive data
		 *            type.
		 *
		 * @details
		 *            The result agrees with this integer only if it ranges
		 *            between -4294967295...4294967295. More precisely, if
		 *            the integer is denoted by <i>a</i> the result will be
		 *            \f[
		 *             sign(a)\cdot(|a|~mod~2^{32})
		 *            \f]
		 *
		 * @return
		 *            This integer as an integer of primitive data type.
		 */
		inline int64_t toInt() const {
			if ( sign() < 0 ) {
				return -(int64_t)this->data[0];
			} else {
				return this->data[0];
			}
		}

		/**
		 * @brief
		 *            Converts this integer to a long double value.
		 *
		 * @details
		 *            The function attempts in a straightforward manner to
		 *            convert the big integer as a long double value. Any
		 *            possible overflows or incorrectnesses due to the
		 *            system's machine epsilon have to be treated afterwards.
		 *
		 * @return
		 *            This integer as a long double value.
		 */
		long double toFloat() const;

		/**
		 * @brief
		 *            Tests if this integer is zero.
		 *
		 * @details
		 *             If the integer is zero, the function returns
		 *             <code>true</code>; otherwise, if this integer is not
		 *             zero, the result will be <code>false</code>.
		 *
		 * @return
		 *             <code>true</code> if the integer is zero and
		 *             <code>false</code> otherwise.
		 */
		inline bool isZero() const {
			return this->size==1&&this->data[0]==0;
		}

		/**
		 * @brief
		 *            Tests if this integer is on.
		 *
		 * @details
		 *            If the integer is one, the function returns
		 *            <code>true</code>; otherwise, if this integer is not
		 *            one, the result will be <code>false</code>.
		 *
		 * @return
		 *            <code>true</code> if the integer is one and
		 *            <code>false</code> otherwise.
		 */
		inline bool isOne() const {
			return this->size==1&&this->data[0]==1;
		}

		/**
		 * @brief
		 *             Assigns this integer by zero.
		 *
		 * @details
		 *             This method is a convenience methdod. Alternatively,
		 *             we can assign a big integer <i>a</i> by zero via
		 *             <pre>
		 *              BigInteger a;
		 *              ...
		 *              a = 0;
		 *             </pre>
		 *             but using this method is more direct.
		 */
		inline void clear() {
			this->size = 1;
			this->data[0] = 0;
		}

		/**
		 * @brief
		 *             Initializes this object to represent a random integer
		 *             being smaller than the absolute value of the specified
		 *             bound.
		 *
		 * @param bound
		 *             Specifies the bound such that this object is being
		 *             assigned with an integer smaller than the absolute
		 *             value of <code>bound</code>.
		 *
		 * @param tryRandom
		 *             Indicates whether the method is advised to use
		 *             a cryptographic random generator.
		 *
		 * @warning
		 *             If not sufficient memory could be allocated, then
		 *             an error message is printed to <code>stderr</code>
		 *             and an exit with status 'EXIT_FAILURE' is caused.
		 */
		void random( const BigInteger & bound , bool tryRandom = false );

		/**
		 * @brief
		 *             Returns the number of bytes needed to represent this
		 *             integer.
		 *
		 * @details
		 *             Let this integer be equals
		 *             \f[
		 *              \pm\sum_{j=0}^bc_j2^j
		 *             \f]
		 *             where \f$c_j\in\{0,1\}\f$. Then this function returns
		 *             \f$\lceil (b+1)/8 \rceil\f$ which is the number of
		 *             bytes needed to represent this integer including a
		 *             leading bit encoding its sign.
		 *
		 * @return
		 *             The number of bytes needed to represent this integer.
		 */
		int getSizeInBytes() const;

		/**
		 * @brief
		 *             Exports this big integer to byte data.
		 *
		 * @details
		 *             This method write bytes to <code>array</code> such
		 *             that the value
		 *             \f[
		 * (array[s-1]\&127)\cdot 256^{s-1}+\sum_{i=0}^{s-2} array[j]\cdot 256^j
		 *             \f]
		 *             equals the absolute value of this integer and
		 *             <i>s=</i>\link getSizeInBytes()\endlink. If this
		 *             is negative, then the leading bit of
		 *             <code>array[s-1]</code> is set to one.
		 *
		 *             To obtain a BigInteger object back from the exported
		 *             data, we may use the \link fromBytes()\endlink
		 *             method.
		 *
		 * @param array
		 *             The byte array into which the data of this integer
		 *             is exported.
		 *
		 * @see fromBytes()
		 * @see getSizeInBytes()
		 *
		 * @warning
		 *             If <code>array</code> is not capable of holding at
		 *             least \link getSizeInBytes()\endlink bytes, then
		 *             calling this method results into undocumented
		 *             behavior.
		 */
		void toBytes( uint8_t *array ) const;

		/**
		 * @brief
		 *             Initializes this big integer by the data contained
		 *             in the specified array.
		 *
		 * @details
		 *             This object will be initialized such that its
		 *             absolute value will be equals
		 *             \f[
		 * (array[s-1]\&127)\cdot 256^{s-1}+\sum_{i=0}^{s-2} array[j]\cdot 256^j
		 *             \f]
		 *             where <i>s=</i><code>sizeInBytes</code>. If the leading
		 *             bit of <code>array[s-1]</code> is set, then the integer
		 *             will be negative.
		 *
		 *             This method can be considered as the inverse of
		 *             the \link toBytes()\endlink method.
		 *
		 * @param array
		 *             The array containing the data from which this big
		 *             integer is initialized.
		 *
		 * @param sizeInBytes
		 *             The number of valid bytes stored in
		 *             <code>array</code>.
		 *
		 * @see toBytes()
		 *
		 * @warning
		 *             If <code>array</code> does not contain
		 *             <code>sizeInBytes</code> valid bytes, then calling
		 *             this method runs into undocumented behavior.
		 *
		 * @warning
		 *             If <code>sizeInBytes</code> is negative, an error
		 *             message will be printed to <code>stderr</code>
		 *             and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message will be printed to <code>stderr</code>
		 *             and the program exits with status 'EXIT_FAILURE'.
		 */
		void fromBytes( const uint8_t *array , int sizeInBytes );

		/**
		 * @brief
		 *             Initializes this big integer by a string in
		 *             decimal representation.
		 *
		 * @details
		 *             If the specified string <code>str</code> contains
		 *             characters that are no digits (except a possible
		 *             leading '-' or '+'), then the integer will be
		 *             initialized by the first digits only and the
		 *             function returns <code>false</code>; otherwise
		 *             if <code>str</code> represents a valid decimal
		 *             integer, the function returns <code>true</code>.
		 *
		 * @param str
		 *             The string with which this big integer is initialized.
		 *
		 * @return
		 *             <code>true</code> if <code>str</code> represents
		 *             a valid decimal integer; otherwise, the function
		 *             returns <code>false</code>.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message will be printed to <code>stderr</code>
		 *             and the program exits with status 'EXIT_FAILURE'.
		 */
		bool fromString( const std::string & str );

		/**
		 * @brief
		 *            Computes the negative of a big integer.
		 *
		 * @details
		 *            After calling this method, the integer represented by
		 *            <i>a</i> is the negative of the integer that was
		 *            represented by <i>b</i> on input.
		 *
		 * @param a
		 *            Will contain the negative of the integer represented
		 *            by <i>b</i>.
		 *
		 * @param b
		 *            The integer of which the negative is computed.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		inline static void negate
		( BigInteger & a , const BigInteger & b ) {
			a = b;
			if ( !a.isZero() ) {
				a.size = -a.size;
			}
		}

		/**
		 * @brief
		 *            Computes the absolute value of a big integer.
		 *
		 * @details
		 *            After calling this method, the integer represented
		 *            by <i>a</i> will be the absolute value of the integer
		 *            that was represented by <i>b</i> on input. In other
		 *            words: If <i>b</i> is negative, <i>a</i> will be assigned
		 *            by the negative of <i>b</i>; otherwise, <i>a</i> will
		 *            be assigned by <i>b</i>.
		 *
		 * @param a
		 *            Will contain the absolute value of <i>b</i>.
		 *
		 * @param b
		 *            The integer of which the absolute value is computed.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		inline static void abs
		( BigInteger & a , const BigInteger & b ) {
			a = b;
			if ( a.size < 0 ) {
				a.size = -a.size;
			}
		}

		/**
		 * @brief
		 *            Compares the absolute values of two integers.
		 *
		 * @details
		 *            If <i>|a|<|b|</i> the result will be -1; if
		 *            <i>|a|>|b|</i>  the result will be 1; otherwise, if
		 *            <i>|a|=|b|</i> the result will be 0. In other words,
		 *            the result will be the sign of <i>|a|-|b|</i>.
		 *
		 * @param a
		 *            First integer.
		 *
		 * @param b
		 *            Second integer.
		 *
		 * @return
		 *           -1 if <i>|a|<|b|</i>, 1 if  <i>|a|>|b|</i>, or 0
		 *           otherwise if <i>|a|=|b|</i>.
		 */
		static int absCompare
		( const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Compares two big integer.
		 *
		 * @details
		 *            If <i>a<b</i> the result will be -1; if
		 *            <i>a>b</i>  the result will be 1; otherwise, if
		 *            <i>a=b</i> the result will be 0. In other words,
		 *            the result will be the sign of <i>a-b</i>.
		 *
		 * @param a
		 *            First integer.
		 *
		 * @param b
		 *            Second integer.
		 *
		 * @return
		 *           -1 if <i>a<b</i>, 1 if  <i>a>b</i>, or 0
		 *           otherwise if <i>a=b</i>.
		 */
		static int compare
		( const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Performs a binary left shift of a big integer by
		 *            <i>n</i> bits.
		 *
		 * @details
		 *            More precisely, if <i>n>=0</i> the method computes
		 *            <i>a</i> such that
		 *            \f[
		 *             a = b\cdot 2^n.
		 *            \f]
		 *            Otherwise, if <i>n<0</i> then
		 *            \f[
		 *         a = sign(b)\cdot \left\lfloor|b|\cdot {2^n}\right\rfloor.
		 *            \f]
		 *
		 * @param a
		 *            Will contain the binary left shift of <i>b</i> by
		 *            <i>n</i> bits.
		 *
		 * @param b
		 *            The integer of which the binary left shift is computed.
		 *
		 * @param n
		 *            Specifies the number of bits the integer <i>b</i> is
		 *            shifted.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static void leftShift
		( BigInteger & a , const BigInteger & b , int n );

		/**
		 * @brief
		 *            Performs a binary right shift of a big integer by
		 *            <i>n</i> bits.
		 *
		 * @details
		 *            More precisely, if <i>n>=0</i> the method computes
		 *            <i>a</i> such that
		 *            \f[
		 *        a = sign(b)\cdot\left\lfloor\frac{|b|}{2^n}\right\rfloor.
		 *            \f]
		 *            Otherwise, if <i>n<0</i> then
		 *            \f[
		 *             a = b\cdot 2^{-n}
		 *            \f]
		 *
		 * @param a
		 *            Will contain the binary right shift of <i>b</i> by
		 *            <i>n</i> bits.
		 *
		 * @param b
		 *            The integer of which the binary right shift is computed.
		 *
		 * @param n
		 *            Specifies the number of bits the integer <i>b</i> is
		 *            shifted.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static void rightShift
		( BigInteger & a , const BigInteger & b , int n );

		/**
		 * @brief
		 *            Computes the sum of two big integers.
		 *
		 * @details
		 *            After calling this method, the big integer <i>c</i>
		 *            will be the sum of the integer <i>a</i> and <i>b</i>
		 *            on input, i.e.
		 *            \f[
		 *             c = a + b
		 *            \f]
		 *
		 * @param c
		 *            Will contain the sum of <i>a</i> and <i>b</i>.
		 *
		 * @param a
		 *            First summand.
		 *
		 * @param b
		 *            Second summand.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::add(BigInteger&,const BigInteger&a,const BigInteger&b)
		 */
		static void add
		( BigInteger & c , const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the difference between two big integers.
		 *
		 * @details
		 *            After calling this method, the big integer <i>c</i>
		 *            will be the difference between the integers <i>a</i>
		 *            and <i>b</i> on input, i.e.
		 *            \f[
		 *             c = a - b.
		 *            \f]
		 *
		 * @param c
		 *            Will contain the difference between <i>a</i> and <i>b</i>.
		 *
		 * @param a
		 *            Minuend.
		 *
		 * @param b
		 *            Subtrahend.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::sub(BigInteger&,const BigInteger&a,const BigInteger&b)
		 */
		static void sub
		( BigInteger & c , const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the product of two big integers.
		 *
		 * @details
		 *            After calling this method, the big integer <i>c</i>
		 *            will be the difference between the integers <i>a</i>
		 *            and <i>b</i> on input, i.e.
		 *            \f[
		 *             c = a\cdot b
		 *            \f]
		 *
		 * @param c
		 *            Will contain the product of <i>a</i> and <i>b</i>.
		 *
		 * @param a
		 *            First factor.
		 *
		 * @param b
		 *            Second factor.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::mul(BigInteger&,const BigInteger&a,const BigInteger&b)
		 */
		static void mul
		( BigInteger & c , const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the quotient of two big integers.
		 *
		 * @details
		 *            The quotient that is computed is
		 *            \f[
		 *             c = \left\lfloor\frac{a}{b}\right\rfloor.
		 *            \f]
		 *
		 * @param c
		 *            Will contain the quotient of <i>a</i> and <i>b</i>.
		 *
		 * @param a
		 *            Numerator.
		 *
		 * @param b
		 *            Denominator.
		 *
		 * @warning
		 *            If <i>b</i> is zero, the method prints an error message
		 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::div(BigInteger&,const BigInteger&a,const BigInteger&b)
		 */
		static void div
		( BigInteger & c , const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the Euclidean division of two integers.
		 *
		 * @details
		 *            Computes
		 *            \f[
		 *             q = \left\lfloor\frac{a}{b}\right\rfloor
		 *            \f]
		 *            and
		 *            \f[
		 *             r = a-q
		 *            \f]
		 *            which guarantees that \f$0\leq|r|<|b|\f$. Moreover,
		 *            if \f$b>0\f$ then \f$r\geq 0\f$ and if \f$b<0\f$
		 *            then \f$r\leq 0\f$.
		 *
		 * @param q
		 *            Will contain the quotient of <i>a</i> divided by
		 *            <i>b</i>.
		 *
		 * @param r
		 *            Will contain the remainder of <i>a</i> divided by
		 *            <i>b</i>.
		 *
		 * @param a
		 *            Numerator.
		 *
		 * @param b
		 *            Denominator.
		 *
		 * @warning
		 *            If <i>b</i> is zero, the method prints an error message
		 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <i>q</i> and <i>r</i> are of the same reference, this
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::divRem(BigInteger&,BigInteger&r,const BigInteger&,const BigInteger&)
		 */
		static void divRem
		( BigInteger & q , BigInteger & r ,
		  const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the remainder of the Euclidean division of two
		 *            big integers.
		 *
		 * @details
		 *            Computes
		 *            \f[
		 *             r=a-\left\lfloor\frac{a}{b}\right\rfloor
		 *            \f]
		 *            which guarantees that \f$0\leq|r|<|b|\f$. Moreover, if
		 *            \f$b>0\f$ then  \f$r\geq 0\f$ and if \f$b<0\f$ then
		 *            \f$r\leq 0\f$.
		 *
		 * @param c
		 *            Will contain the remainder of <i>a</i> divided by
		 *            <i>b</i>.
		 *
		 * @param a
		 *            Numerator.
		 *
		 * @param b
		 *            Denominator.
		 *
		 * @warning
		 *            If <i>b</i> is zero, the method prints an error message
		 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @see thimble::rem(BigInteger&,const BigInteger&a,const BigInteger&b)
		 */
		inline static void rem
		( BigInteger & c , const BigInteger & a , const BigInteger & b ) {

			BigInteger q;
			divRem(q,c,a,b);
		}

		/**
		 * @brief
		 *            Computes the greatest common divisor of two big
		 *            integers.
		 *
		 * @param g
		 *            On output, the greatest common divisor of <code>a</code>
		 *            and <code>b</code>.
		 *
		 * @param a
		 *            First integer.
		 *
		 * @param b
		 *            Second integer.
		 *
		 * @warning
		 *            If <code>a</code> and <code>b</code> represent zero
		 *            integer, an error message will be printed to
		 *            <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static void gcd
		( BigInteger & g , const BigInteger & a , const BigInteger & b );

		/**
		 * @brief
		 *            Computes the factorial of an integer as a big integer.
		 *
		 * @details
		 *            The result will be the integer as a big integer that is
		 *            equals
		 *            \f[
		 *             n!=\prod_{i=0}^ni.
		 *            \f]
		 *            If <i>n=0</i> the result will be 1.
		 *
		 * @param n
		 *            The integer of which the factorial is computed.
		 *
		 * @return
		 *            The factorial of <i>n</i>.
		 *
		 * @warning
		 *            If <i>n</i> is negative the function prints an error
		 *            message to <code>stdint</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static BigInteger factorial( int n );

		/**
		 * @brief
		 *            Computes the binomial coefficient &quot;<i>n</i> over
		 *            <i>k</i>&quot; as a big integer.
		 *
		 * @details
		 *            The result will be equals
		 *            \f[
		 *             \left({n\atop k}\right)=\frac{n!}{k!\cdot(n-k)!}
		 *             =\prod_{i=1}^k\frac{n+1-j}{j}.
		 *            \f]
		 *
		 * @param n
		 *            Integer as <code>int</code>.
		 *
		 * @param k
		 *            Integer as <code>int</code>.
		 *
		 * @return
		 *            The binomial coefficient &quot;<i>n</i> over <i>k</i>
		 *            &quot; as a big integer.
		 *
		 * @warning
		 *            If <i>n</i> or <i>k</i> are negative the function prints
		 *            an error message to <code>stderr</code> and exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <i>k</i> is greater than <i>n</i> the function prints
		 *            and error message to <code>stderr</code> and exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static BigInteger binomial( int n , int k );

		/**
		 * @brief
		 *            Expands a float from a fraction of big integers.
		 *
		 * @param num
		 *            Numerator.
		 *
		 * @param den
		 *            Denonimator.
		 *
		 * @param numBits
		 *            Number of bits with which the fractional part of
		 *            the ratio is approximated.
		 *
		 * @return
		 *            A float approximating <code>num/den</code>
		 *
		 * @warning
		 *            If <code>den</code> represents the zero number, an
		 *            error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		static long double expand
		( const BigInteger & num , const BigInteger & den ,
		  int numBits = 128 );

		/**
		 * @brief
		 *            Returns the negative of this integer.
		 *
		 * @details
		 *            By providing this function, we may run expressions as
		 *            <pre>
		 *             BigInteger a = ...;
		 *
		 *             BigInteger b = -a;
		 *            </pre>
		 *            to compute the negative integer <code>b</code> of
		 *            <code>a</code>.
		 *
		 * @return
		 *            A big integer representing the negative of this big
		 *            integer.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		inline BigInteger operator-() const {
			BigInteger b;
			negate(b,*this);
			return b;
		}

		/**
		 * @brief
		 *            Returns a copy of this big integer.
		 *
		 * @details
		 *            By providing this function, we may run expressions as
		 *            <pre>
		 *             BigInteger a = ...;
		 *
		 *             BigInteger b = +a;
		 *            </pre>
		 *
		 * @return
		 *            A big integer representing a copy of this big integer.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            method prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 */
		inline BigInteger operator+() const {
			BigInteger b(*this);
			return b;
		}

        /**
         * @brief
         *            Enables that the sum of two big integers can be computed
         *            via the '+'-operator.
         *
         * @param b
         *            Summand.
         *
         * @return
         *            The sum of this integer and <i>b</i>.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger operator+( const BigInteger & b ) const {
            BigInteger c;
            add(c,*this,b);
            return c;
        }

        /**
         * @brief
         *           Enables the the difference between two big integers can be
         *           computed via the '-'-operator.
         *
         * @param b
         *           Subtrahend.
         *
         * @return
         *           The difference between this integer and <i>b</i>.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger operator-( const BigInteger & b ) const {
            BigInteger c;
            sub(c,*this,b);
            return c;
        }


        /**
         * @brief
         *            Enables that the product of two big integers can be
         *            computed via the '*'-operator.
         *
         * @param b
         *            Factor.
         *
         * @return
         *            The product of this integer and <i>b</i>.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger operator*( const BigInteger & b ) const {
            BigInteger c;
            mul(c,*this,b);
            return c;
        }


        /**
         * @brief
         *            Enables that the integral quotient formed by two big
         *            integers can be computed via the '/'-operator.
         *
         * @param b
         *            Denominator.
         *
         * @return
         *            The integral quotient of this integer divided by
         *            <i>b</i>.
         *
         * @warning
         *            If <i>b</i> is zero, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger operator/( const BigInteger & b ) const {
            BigInteger c;
            div(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Enables that the remainder of an integer modulo another
         *            integer can be computed via the '%'-operator.
         *
         * @param b
         *            Denominator.
         *
         * @return
         *            The remainder of this integer modulo <i>b</i>.
         *
         * @warning
         *            If <i>b</i> is zero, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger operator%( const BigInteger & b ) const {
            BigInteger c;
            rem(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Enables that the sum of two integers can be computed
         *            via the '+='-operator.
         *
         * @param b
         *            Summand.
         *
         * @return
         *            A reference to this integer to which <i>b</i> has
         *            been added.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger & operator+=( const BigInteger & b ) {
            add(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Enables that the difference betweem two integers can be
         *            computed via the '-='-operator.
         *
         * @param b
         *            Subtrahend.
         *
         * @return
         *            A reference to this integer from which <i>b</i> has
         *            been subtracted.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger & operator-=( const BigInteger & b ) {
            sub(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Enables that the product of two integers can be computed
         *            via the '*='-operator.
         *
         * @param b
         *            Factor.
         *
         * @return
         *            A reference to this integer to which <i>b</i> has
         *            been multiplied.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger & operator*=( const BigInteger & b ) {
            mul(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Enables that the integral quotient formed by two
         *            integers can be computed via the '+='-operator.
         *
         * @param b
         *            Denominator.
         *
         * @return
         *            A reference to this integer which has been replaced
         *            by the quotient of this instance on input divided
         *            by <i>b</i>.
         *
         * @warning
         *            If <i>b</i> is zero, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger & operator/=( const BigInteger & b ) {
            div(*this,*this,b);
            return *this;
        }

        /**
         * @brief
         *            Enables the the remainder of an integer modulo
         *            another integer can be computed via the '%='-operator.
         *
         * @param b
         *            Denominator.
         *
         * @return
         *            A reference to this integer which has been replaced by
         *            this instance on inputed modulo <i>b</i>.
         *
         * @warning
         *            If <i>b</i> is zero, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         *
         * @warning
         *            If not sufficient memory could be allocated, the
         *            method prints an error message to <code>stderr</code>
         *            and exits with status 'EXIT_FAILURE'.
         */
        inline BigInteger & operator%=( const BigInteger & b ) {
            rem(*this,*this,b);
            return *this;
        }
    };

	/**
	 * @brief
	 *            Computes the sum of two big integers.
	 *
	 * @details
	 *            After calling this method, the big integer <i>c</i>
	 *            will be the sum of the integer <i>a</i> and <i>b</i>
	 *            on input, i.e.
	 *            \f[
	 *             c = a + b
	 *            \f]
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::add(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::add()</code>
	 *            statements.
	 *
	 * @param c
	 *            Will contain the sum of <i>a</i> and <i>b</i>.
	 *
	 * @param a
	 *            First summand.
	 *
	 * @param b
	 *            Second summand.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::add(BigInteger&,const BigInteger&a,const BigInteger&b)
	 */
	inline void add
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {
		BigInteger::add(c,a,b);
	}

	/**
	 * @brief
	 *            Computes the difference between two big integers.
	 *
	 * @details
	 *            After calling this method, the big integer <i>c</i>
	 *            will be the difference between the integers <i>a</i>
	 *            and <i>b</i> on input, i.e.
	 *            \f[
	 *             c = a - b.
	 *            \f]
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::sub(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::sub()</code>
	 *            statements.
	 *
	 * @param c
	 *            Will contain the difference between <i>a</i> and <i>b</i>.
	 *
	 * @param a
	 *            Minuend.
	 *
	 * @param b
	 *            Subtrahend.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::sub(BigInteger&,const BigInteger&a,const BigInteger&b)
	 */
	inline void sub
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {
		BigInteger::sub(c,a,b);
	}

	/**
	 * @brief
	 *            Computes the product of two big integers.
	 *
	 * @details
	 *            After calling this method, the big integer <i>c</i>
	 *            will be the difference between the integers <i>a</i>
	 *            and <i>b</i> on input, i.e.
	 *            \f[
	 *             c = a\cdot b
	 *            \f]
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::mul(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::mul()</code>
	 *            statements.
	 *
	 * @param c
	 *            Will contain the product of <i>a</i> and <i>b</i>.
	 *
	 * @param a
	 *            First factor.
	 *
	 * @param b
	 *            Second factor.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::mul(BigInteger&,const BigInteger&a,const BigInteger&b)
	 */
	inline void mul
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {
		BigInteger::mul(c,a,b);
	}

	/**
	 * @brief
	 *            Computes the quotient of two big integers.
	 *
	 * @details
	 *            The quotient that is computed is
	 *            \f[
	 *             c = \left\lfloor\frac{a}{b}\right\rfloor.
	 *            \f]
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::div(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::div()</code>
	 *            statements.
	 *
	 * @param c
	 *            Will contain the quotient of <i>a</i> and <i>b</i>.
	 *
	 * @param a
	 *            Numerator.
	 *
	 * @param b
	 *            Denominator.
	 *
	 * @warning
	 *            If <i>b</i> is zero, the method prints an error message
	 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::div(BigInteger&,const BigInteger&a,const BigInteger&b)
	 */
	inline void div
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {
		BigInteger::div(c,a,b);
	}

	/**
	 * @brief
	 *            Computes the Euclidean division of two integers.
	 *
	 * @details
	 *            Computes
	 *            \f[
	 *             q = \left\lfloor\frac{a}{b}\right\rfloor
	 *            \f]
	 *            and
	 *            \f[
	 *             r = a-q
	 *            \f]
	 *            which guarantees that \f$0\leq|r|<|b|\f$. Moreover,
	 *            if \f$b>0\f$ then \f$r\geq 0\f$ and if \f$b<0\f$
	 *            then \f$r\leq 0\f$.
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::divRem(BigInteger&,BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::div()</code>
	 *            statements.
	 *
	 * @param q
	 *            Will contain the quotient of <i>a</i> divided by
	 *            <i>b</i>.
	 *
	 * @param r
	 *            Will contain the remainder of <i>a</i> divided by
	 *            <i>b</i>.
	 *
	 * @param a
	 *            Numerator.
	 *
	 * @param b
	 *            Denominator.
	 *
	 * @warning
	 *            If <i>b</i> is zero, the method prints an error message
	 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <i>q</i> and <i>r</i> are of the same reference, this
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::divRem(BigInteger&,BigInteger&r,const BigInteger&,const BigInteger&)
	 */
	inline void divRem
	( BigInteger & q , BigInteger & r ,
	  const BigInteger & a , const BigInteger & b ) {
		BigInteger::divRem(q,r,a,b);
	}

	/**
	 * @brief
	 *            Computes the remainder of the Euclidean division of two
	 *            big integers.
	 *
	 * @details
	 *            Computes
	 *            \f[
	 *             r=a-\left\lfloor\frac{a}{b}\right\rfloor
	 *            \f]
	 *            which guarantees that \f$0\leq|r|<|b|\f$. Moreover, if
	 *            \f$b>0\f$ then  \f$r\geq 0\f$ and if \f$b<0\f$ then
	 *            \f$r\leq 0\f$.
	 *            <br><br>
	 *            This method is a convenience method wrapped around
	 *            <code>
	 * \link BigInteger::rem(BigInteger&,const BigInteger&,const BigInteger&)\endlink
	 *            </code>
	 *            to avoid frequent uses of <code>BigInteger::div()</code>
	 *            statements.
	 *
	 * @param c
	 *            Will contain the remainder of <i>a</i> divided by
	 *            <i>b</i>.
	 *
	 * @param a
	 *            Numerator.
	 *
	 * @param b
	 *            Denominator.
	 *
	 * @warning
	 *            If <i>b</i> is zero, the method prints an error message
	 *            to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be allocated, the
	 *            method prints an error message to <code>stderr</code>
	 *            and exits with status 'EXIT_FAILURE'.
	 *
	 * @see BigInteger::rem(BigInteger&,const BigInteger&a,const BigInteger&b)
	 */
	inline void rem
	( BigInteger & c , const BigInteger & a , const BigInteger & b ) {
		BigInteger::rem(c,a,b);
	}

	/**
	 * @brief
	 *            Prints the big integer in decimal text representation to the
	 *            specified output stream.
	 *
	 * @details
	 *            The implementation of this function allows statements such
	 *            as
	 *            <pre>
	 *             BigInteger a = ...;
	 *             cout << a << endl;
	 *            </pre>
	 *            which prints a big integer <i>a</i> to <code>stdout</code>.
	 *
	 * @param out
	 *            The output stream to where the big integer is printed.
	 *
	 * @param a
	 *            The big integer that is printed to the specifed output
	 *            stream.
	 *
	 * @return
	 *            Returns a reference to the specified output stream.
	 */
	THIMBLE_DLL std::ostream & operator<<
			( std::ostream & out , const BigInteger & a );
}


#endif /* THIMBLE_BIGINTEGER_H_ */
