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
 * @file SmallBinaryFieldPolynomial.h
 *
 * @brief
 *            Provides a class to represent and compute with polynomials
 *            having coefficients in a small binary galois field.
 *
 * @details
 * @section sec_polynomials Polynomials
 *
 * In THIMBLE computations with polynomials having coefficients in a
 * binary field of small degree (16, say) can be performed. In order
 * to create objects that represent these polynomials it is necessary
 * to create an object that presents the finite field to which
 * the polynomials' coefficients belong. Therefore, we assume that
 * a \link thimble::SmallBinaryField SmallBinaryField\endlink object
 * <code>gf</code> has been created, e.g., via
 * <pre>
 *  SmallBinaryField gf(16);
 * </pre>
 * which represents a finite field of size \f$2^{16}\f$; for more details
 * we refer to the documentation about @ref sec_finitefields.
 *
 * To create an \link thimble::SmallBinaryFieldPolynomial
 * SmallBinaryFieldPolynomial\endlink object that represents the zero
 * polynomial over the finite field represented by <code>gf</code> we may
 * run
 * <pre>
 *  SmallBinaryFieldPolynomial f(gf);
 * </pre>
 * constructing a \link thimble::SmallBinaryFieldPolynomial
 * SmallBinaryFieldPolynomial\endlink.
 *
 * The polynomial represented by <code>f</code> can be modified by manually
 * specifying the value of its coefficients via the member method
 * \link thimble::SmallBinaryFieldPolynomial::setCoeff() setCoeff()\endlink.
 * More specifically, if we want to replace the <code>i</code>th coefficient
 * (<code>i>=0</code>) of <code>f</code> by the element represented by
 * <pre>
 *  uint32_t c;
 * </pre>
 * we may run
 * <pre>
 *  f.setCoeff(i,c);
 * </pre>
 * The <code>i</code>th coefficient can be accessed via the
 * \link thimble::SmallBinaryFieldPolynomial::getCoeff(int)const
 * getCoeff(int)\endlink function
 * <pre>
 *  uint32_t c = f.getCoeff(i);
 * </pre>
 *
 * @subsection sec_polynomials_addsub Addition/Subtraction
 *
 * There are multiple equivalent ways of forming the sum of two polynomials
 * with coefficients in a
 * <pre>
 *  SmallBinaryField gf = ...;
 * </pre>
 * For two polynomials
 * <pre>
 *  SmallBinaryFieldPolynomial f(gf) , g(gf);
 * </pre>
 * their sum can be computed using the following equivalent
 * expressions
 * <pre>
 *  SmallBinaryFieldPolynomial h(gf);
 *
 *  h = f + g;
 *  add(h,f,g);
 *  SmallBinaryFieldPolynomial::add(h,f,g);
 * </pre>
 * Similarly, the difference between the polynomials can be formed
 * via
 * <pre>
 *  h = f - g;
 *  sub(h,f,g);
 *  SmallBinaryFieldPolynomial::sub(h,f,g);
 * </pre>
 *
 * @attention In a binary field subtraction is equivalent to addition and
 * it is actually not necessary to distinguish between these two operations
 * in the context of a \link thimble::SmallBinaryFieldPolynomial
 * SmallBinaryFieldPolynomial\endlink; however, in some situations it can be
 * useful to distinguish between addition and subtraction for better
 * readability of program code.
 *
 * @subsection sec_polynomials_mul Multiplication
 *
 * To multiply two polynomials with coefficients in a
 * <pre>
 *  SmallBinaryField gf = ...;
 * </pre>
 * denoted by
 * <pre>
 *  SmallBinaryFieldPolynomial f(gf) , g(gf);
 * </pre>
 * we may run one of the following equivalent expressions
 * <pre>
 *  SmallBinaryFieldPolynomial h(gf);
 *
 *  h = f * g;
 *  mul(h,f,g);
 *  SmallBinaryFieldPolynomial::mul(h,f,g);
 * </pre>
 *
 * @subsection sec_polynomials_divrem Division (with Remainder)
 *
 * To compute the WRITE ME
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_SMALLBINARYFIELDPOLYNOMIAL_H_
#define THIMBLE_SMALLBINARYFIELDPOLYNOMIAL_H_

#include <stdint.h>
#include <cstdlib>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Polynomials with coefficients in a small binary field.
	 *
	 * @details
	 *            Before instances of this class are generated it is necessary
	 *            to specify the coefficient field. Therefore, an instance of
	 *            the \link SmallBinaryField\endlink class has to be generated
	 *            first. The constructors of this class need an instance of
	 *            a \link SmallBinaryField\endlink. For example, the following
	 *            code snipet creates a polynomial with coefficients in a
	 *            finite field of cardinality 2^16:
	 *            <pre>
	 *             SmallBinaryField gf(16);
	 *
	 *             SmallBinaryFieldPolynomial f(gf);
	 *            </pre>
	 */
	class THIMBLE_DLL SmallBinaryFieldPolynomial {

		friend class SmallBinaryFieldBivariatePolynomial;

	private:

		/**
		 * @brief
		 *           Pointer to the underlying finite field.
		 *
		 * @details
		 *           The finite field must be specified before polynomials
		 *           with coefficients in the field can be created.
		 */
		const SmallBinaryField *gfPtr;

		/**
		 * @brief
		 *           Coefficients of the polynomial
		 *
		 * @details
		 *           If the polynomial
		 *           \f$f(X)=c_0+c_1\cdot X+...+c_d\cdot X^d\f$ is represented
		 *           by this instance, then \f$coefficients[i]=c_i\f$.
		 *           <br>
		 *           The coefficients are elements in the field represented
		 *           by the content of \link gfPtr\endlink. Elements of the
		 *           finite field, in turn, are represented as 32-bit integers.
		 *           For details, we refer to \link SmallBinaryField\endlink.
		 */
		uint32_t *coefficients;

		/**
		 * @brief
		 *           The degree of the polynomial
		 *
		 * @details
		 *           The degree species the index of the highest relevant
		 *           coefficient listed in \link coefficients\endlink.
		 */
		int degree;

		/**
		 * @brief
		 *           The number of elements the field
		 *           \link coefficients\endlink can hold.
		 *
		 * @details
		 *           It should be always asserted that this number indicates
		 *           the number of 32-bit word the field
		 *           \link coefficients\endlink can hold. Thus the statement
		 *           <code>
		 *            \link degree\endlink < \link capacity\endlink
		 *           </code>
		 *           always must be <code>true</code>.
		 */
		int capacity;

		/**
		 * @brief
		 *           Normalizes the degree of the polynomial
		 *
		 * @details
		 *           When coefficients of the polynomial become 0 (e.g., due
		 *           to a call of \link setCoeff(int,uint32_t)\endlink or)
		 *           the degree of the polynomial is not necessarily
		 *           equals to \link degree\endlink but can be smaller.
		 *           This method determines the highest non-zero coefficient
		 *           of the polynomial and sets \link degree\endlink to its
		 *           index. In this way, the method ensures that
		 *           \link degree\endlink correctly denotes the degree
		 *           of the polynomial.
		 */
		void normalize();

		/**
		 * @brief
		 *            Computes the product of two polynomials.
		 *
		 * @details
		 *            This methods stores the product of <code>f</code>
		 *            and <code>g</code> in <code>h</code>.
		 *
		 * @param h
		 *            Will contain the product of <code>f</code> and
		 *            <code>g</code>.
		 *
		 * @param f
		 *            First factor polynomial.
		 *
		 * @param g
		 *            Second factor polynomial
		 *
		 * @warning
		 *            If it is not possible to ensure the capacity of
		 *            <code>h</code> to store the result an error message
		 *            is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>h</code>, <code>f</code>, and <code>g</code>
		 *            do not all correspond to the same finite field calling
		 *            this method may run into unexpected behavior.
		 *
		 */
		static void mulUncheck
		( SmallBinaryFieldPolynomial & h ,
		  const SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g );

	public:

		/**
		 * @brief
		 *            Standard constructor which constructs the zero polynomial
		 *            defined over a field with two elements.
		 */
		inline SmallBinaryFieldPolynomial() {
			this->gfPtr = &SmallBinaryField::binary();
			this->coefficients = NULL;
			this->degree = -1;
			this->capacity = 0;
		}

		/**
		 * @brief
		 *           Standard constructor of a polynomial with coefficients
		 *           in the specified binary field.
		 *
		 * @details
		 *           There is no constructor with empty argument list because
		 *           polynomials need to know the field in where it can have
		 *           coefficients.
		 *
		 * @param gf
		 *           The finite field in where the polynomial can have
		 *           coefficients.
		 *
		 * @warning
		 *           If polynomials with coefficients in a finite field
		 *           instance <code>gf</code> have been created but the field
		 *           changes or is destroyed (e.g., by calling its destructor
		 *           explicitly or implicitly) this can result in unexpected
		 *           behavior when the polynomials are used afterwards.
		 *           Therefore, it is recommended to not to change or destroy
		 *           the finite field unless the polynomials are not used
		 *           afterwards.
		 */
		inline SmallBinaryFieldPolynomial( const SmallBinaryField & gf ) {
			this->gfPtr = &gf;
			this->coefficients = NULL;
			this->degree = -1;
			this->capacity = 0;
        }

		/**
		 * @brief
		 *           Assigns this polynomial with the specified polynomial.
		 *
		 * @details
		 *           Sets this polynomial to a copy of the polynomial
         *           <code>f</code>. The method assumes that the assigned
         *           and the to-be-assigned polynomial are related to the
         *           same finite field, i.e., their defining polynomial
         *           must be equal.
		 *
		 * @param f
		 *           The polynomial of which this polynomial will become
		 *           a copy of.
		 *
		 * @return
		 *           A reference to this polynomial.
         *
         * @warning
         *           If the finite field of the assigned and to-be-assigned
         *           polynomial are different, i.e., if their defining
         *           polynomials are different, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
		 * @warning
		 *           If not sufficient memory can be allocated for this
		 *           polynomial to store a copy of <code>f</code> an
		 *           error message is printed to <code>stderr</code>
		 *           an the program exits with status 'EXIT_FAILURE'.
		 */
		void assign( const SmallBinaryFieldPolynomial & f );

		/**
		 * @brief
		 *           Assignment operator.
		 *
		 * @details
		 *           Sets this polynomial to a copy of the polynomial
         *           <code>f</code>. The method assumes that the assigned
         *           and the to-be-assigned polynomial are related to the
         *           same finite field, i.e., their defining polynomial
         *           must be equal.
		 *           <br>
		 *           <br>
		 *           The operator wraps around the
		 *           \link SmallBinaryFieldPolynomial.assign\endlink
		 *           method.
		 *
		 * @param f
		 *           The polynomial of which this polynomial will become
		 *           a copy of.
		 *
		 * @return
		 *           A reference to this polynomial.
		 *
         * @warning
         *           If the finite field of the assigned and to-be-assigned
         *           polynomial are different, i.e., if their defining
         *           polynomials are different, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
		 * @warning
		 *           If not sufficient memory can be allocated for this
		 *           polynomial to store a copy of <code>f</code> an
		 *           error message is printed to <code>stderr</code>
		 *           an the program exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryFieldPolynomial &operator=
				( const SmallBinaryFieldPolynomial & f ) {
			assign(f);
			return *this;
		}

        /**
         * @brief
         *           Copy constructor.
         *
         * @details
         *           Creates a polynomial that is a copy of the polynomial
         *           <code>f</code>.
         *
         * @param f
         *           The polynomial of which the copy is created by this
         *           constructor
         */
        SmallBinaryFieldPolynomial
        ( const SmallBinaryFieldPolynomial & f ) {
            this->gfPtr = f.gfPtr;
            this->coefficients = NULL;
            this->degree = -1;
            this->capacity = 0;
            assign(f);
        }

		/**
		 * @brief
		 *            Swaps this polynomial's content with the content of
		 *            <code>f</code>
		 *
		 * @param f
		 *            Polynomial of which this polynomial changes the content
		 */
		inline void swap( SmallBinaryFieldPolynomial & f ) {

			const SmallBinaryField *gfPtr;
			uint32_t *coefficients;
			int degree , capacity;

			gfPtr              = f.gfPtr;
			coefficients       = f.coefficients;
			degree             = f.degree;
			capacity           = f.capacity;

			f.gfPtr            = this->gfPtr;
			f.coefficients     = this->coefficients;
			f.degree           = this->degree;
			f.capacity         = this->capacity;

			this->gfPtr        = gfPtr;
			this->coefficients = coefficients;
			this->degree       = degree;
			this->capacity     = capacity;
		}

		/**
		 * @brief
		 *           Destructor.
		 *
		 * @details
		 *           The destrutor frees all the memory needed to hold the
		 *           coefficients of the polynomial.
		 *
		 */
		inline ~SmallBinaryFieldPolynomial() {
			free(this->coefficients);
		}

		/**
		 * @brief
		 *           Ensures that the polynomial has enough capacity
		 *           to hold at least <code>newCapacity</code> elements.
		 *
		 * @details
		 *           If the polynomial's coefficients up to index
		 *           <code>i</code> change there will be no reallocation
		 *           of the coefficient array unless <code>i</code> exceeds
		 *           the polynomials capacity. If the polynomial's capacity
		 *           is already larger than <code>newCapacity</code> this
		 *           method will have no effect. Otherwise, the polynomial's
		 *           coefficient array will be reallocated such that it can
		 *           hold at least <code>newCapacity</code> elements.
		 *
		 * @param newCapacity
		 *           The new capacity of the polynomial
		 *
		 * @warning
		 *           If the method fails to allocate enough memory to
		 *           hold at least <code>newCapacity</code> elements an
		 *           error message is printed to <code>stderr</code> and
		 *           the program exits with status 'EXIT_FAILURE'.
		 *
		 */
		void ensureCapacity( int newCapacity );

		/**
		 * @brief
		 *           Access the number of coefficients the polynomial can
		 *           hold.
		 *
		 * @details
		 *           If a coefficient <i>c</i> of the polynomial is set via
		 *           <code>setCoeff(i,c)</code> such that <i>i</i> is smaller
		 *           than the polynomials capacity, there is no need to
		 *           reallocate the coefficient array which would be
		 *           imperformant. Otherwise if <i>i</i> is larger than
		 *           the capacity, setting a coefficient would cause a
		 *           reallocation. Thus, if it is known in advance that the
		 *           polynomial gets along with a certain number of, one may
		 *           ensure its capacity once (<code>ensureCapacity()</code>)
		 *           to suppress imperformant further reallocations.
		 *
		 * @return
		 *           The number of coefficients this polynomial can hold.
		 */
		inline int getCapacity() const {
			return this->capacity;
		}

		/**
		 * @brief
		 *           Access the binary finite field the polynomial
		 *           has its coefficients.
		 *
		 * @return
		 *           A constant reference to this polynomial's underlying
		 *           binary finite field.
		 */
		inline const SmallBinaryField & getField() const {
			return this->gfPtr[0];
		}

		/**
		 * @brief
		 *           Access the degree of the polynomial.
		 *
		 * @details
		 *           If the polynomial is
		 *           \f$f(X)=c_0+c_1\cdot X+...+c_d\cdot X^d\f$ with
		 *           \f$c_d\neq 0\f$ then the result will be \f$c_d\f$.
		 *           If the polynomial is the zero polynomial the returned
		 *           degree will be -1.
		 *
		 * @return
		 *           The degree of the polynomial.
		 */
		inline int deg() const {
			return this->degree;
		}

        /**
         * @brief
         *           Access the array in which the coefficients of the
         *           polynomial are stored.
         *
         * @return
         *           The array of \link deg()\endlink+1 significant
         *           coefficients of the polynomial encoded as unsigned
         *           32 bit integers.
         */
        inline const uint32_t* getData() const {

            return this->coefficients;
        }

		/**
		 * @brief
		 *           Check whether the polynomial is zero.
		 *
		 * @details
		 *           If this polynomial is the zero polynomial the result
		 *           will be <code>true</code> and <code>false</code>
		 *           otherwise.
		 *
		 * @return
		 *           <code>true</code> if the polynomial is the zero
		 *           polynomial and <code>false</code> otherwise.
		 */
		inline bool isZero() const {
			return this->degree < 0;
		}

		/**
		 * @brief
		 *           Access the <i>i</i>-th coefficient of the polynomial.
		 *
		 * @details
		 *           If the polynomial is
		 *           \f$f(X)=c_0+c_1\cdot X+...+c_d\cdot X^d\f$ with
		 *           \f$d\f$ the degree of the polynomial then the
		 *           result will be \f$c_i\f$ if \f$0\leq i\leq d\f$. If
		 *           \f$i<0\f$ or if \f$i>d\f$ the result will be 0.
		 */
		inline uint32_t getCoeff( int i ) const {

			if ( i < 0 || i > this->degree ) {
				return 0;
			}

			return this->coefficients[i];
		}

		/**
		 * @brief
		 *           Determine whether this polynomials is equals to the
		 *           given polynomial.
		 *
		 * @details
		 *           If this polynomial relevant coefficients are the same
		 *           as the coefficients of <code>f</code> the function
		 *           return <code>true</code> and <code>false</code>
		 *           otherwise.
		 *
		 * @return
		 *           <code>true</code> if this polynomial's relevant
		 *           coefficients are the same as the coefficients of
		 *           <code>f</code> and <code>false</code> otherwise.
		 *
		 * @warning
		 *           If this polynomial is defined over a finite field which
		 *           is different to the finite field of <code>f</code>
		 *           (i.e. if the corresponding reference is different) then
		 *           an error message is printed to <code>stderr</code>
		 *           and the program exits with status 'EXIT_FAILURE'.
		 */
		bool equals( const SmallBinaryFieldPolynomial & f) const;

		/**
		 * @brief
		 *            Sets/Replaces the <i>i</i>-th coefficient of the
		 *            polynomial.
		 *
		 * @details
		 *            If the polynomial is
		 *            \f$f(X)=c_0+c_1\cdot X+...+c_i\cdot X^i+...\f$
		 *            then the polynomial is changed by replacing the
		 *            coefficient \f$c_i\f$ by \f$c\f$, i.e.
		 *            \f$c_i\leftarrow c\f$.
		 *
		 * @param i
		 *            The index of the coefficient that is set/replaced
		 *
		 * @param c
		 *            The new <i>i</i>-th coefficient of the polynomial
		 *
		 * @warning
		 *            The method ensures that the coefficient array can
		 *            hold at least <i>i+1</i> elements. Thus, if
		 *            <i>i</i> is too large the method may attempt to
		 *            reallocate an array of which size does not fit in
		 *            memory. If this is the case (or if <i>i</i> is smaller
		 *            than 0), an error message is
		 *            printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		void setCoeff( int i , uint32_t c );

		/**
		 * @brief
		 *           Sets this polynomial to the zero polynomial,
		 *           i.e., \f$f(X)\equiv 0\f$.
		 *
		 * @warning
		 *           If the polynomial could not be ensured to contain
		 *           at least one coefficient an error message is printed
		 *           to <code>stderr</code> and the program exits with
		 *           status 'EXIT_FAILURE'.
		 */
		inline void setZero() {
			ensureCapacity(1);
			this->coefficients[0] = 0;
			this->degree = -1;
		}

		/**
		 * @brief
		 *           Sets this polynomial to the constant polynomial
		 *           equals 1, i.e. \f$f(X)\equiv 1\f$.
		 *
		 * @warning
		 *           If the polynomial could not be ensured to contain
		 *           at least one coefficient an error message is printed
		 *           to <code>stderr</code> and the program exits with
		 *           status 'EXIT_FAILURE'.
		 */
		inline void setOne() {
			ensureCapacity(1);
			this->coefficients[0] = 1;
			this->degree = 0;
		}

		/**
		 * @brief
		 *           Sets this polynomial to the monomial \f$f(X)=X\f$
		 *
		 * @warning
		 *           If the polynomial could not be ensured to contain
		 *           at least two coefficient an error message is printed
		 *           to <code>stderr</code> and the program exits with
		 *           status 'EXIT_FAILURE'.
		 */
		inline void setX() {
			ensureCapacity(2);
			this->coefficients[0] = 0;
			this->coefficients[1] = 1;
			this->degree = 1;
		}

		/**
		 * @brief
		 *            Replaces this polynomial by a random polynomial
		 *            of degree smaller than <code>size</code>
		 *
		 * @details
		 *            Let \f$d=size-1\f$. Then the polynomial will be
		 *            replaced by the polynomial
		 *            \f$f(X)=c_0+c_1\cdot X+...+c_dX^d\f$ where
		 *            the \f$c_i\f$ are chosen uniformly at random.
		 *            As a consequence, the degree of the polynomial
		 *            needs not to be equals to  <code>size-1</code> but it
		 *            smaller than <code>size</code>. Consequently, if
		 *            \f$size\leq 0\f$ then the polynomial is set to the zero
		 *            polynomial.
		 *            <br><br>
		 *            If <code>tryRandom</code> is <code>true</code> the
		 *            method is advised to use a cryptographic random
		 *            generator; otherwise, if <code>tryRandom</code> is
		 *            <code>false</code> the method wraps around the standard
		 *            random generator function <code>rand()</code>.
		 *
		 * @param size
		 *            The number of successive coefficients that are
		 *            randomly generated
		 *
		 * @param tryRandom
		 *            Indicates whether the method is advised to use
		 *            a cryptographic random generator.
		 */
		void random( int size , bool tryRandom = false );

		/**
		 * @brief
		 *            Evaluates the polynomial at a field element.
		 *
		 * @details
		 *            If this polynomial is
		 *            \f$f(X)=\sum_{i=0}^{d-1}c_i\cdot X^i\f$ then the
		 *            function computes
		 *            \f$y\leftarrow \sum_{i=0}^{d-1}c_i\cdot x^i\f$
		 *            and returns the result.
		 *            <br><br>
		 *            The implementation of the function is an implementation
		 *            of Horner's method.
		 *
		 * @param x
		 *            The locator of the finite field
		 *
		 * @return
		 *            The evaluation of the polynomial at <code>x</code>
		 *
		 * @warning
		 *            If <code>x</code> is not a valid element of the
		 *            polynomial's underlying finite field calling this
		 *            function may run into unexpected behavior.
		 */
		uint32_t eval( uint32_t x ) const;

		/**
		 * @brief
		 *            Replaces the polynomial by the polynomial that
		 *            consists of linear factors corresponding to
		 *            having the specified elements as root.
		 *
		 * @details
		 *            The polynomial that is built is
		 *            \f$f(X)=\prod_{i=0}^{n-1}(X-x[i])\f$.
		 *
		 * @param x
		 *            Contains the <code>n</code> roots in the polynomial's
		 *            underlying finite field.
		 *
		 * @param n
		 *            Number of roots contained in <code>x</code>.
		 *
		 * @warning
		 *            If <code>x</code> does not contain at least
		 *            <code>n</code> valid elements in the polynomial's
		 *            underlying binary field then the method may run into
		 *            unexpected behavior.
		 *
		 * @warning
		 *            If not sufficient memory can be allocated to compute
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 */
		void buildFromRoots( const uint32_t *x , int n );

		/**
		 * @brief
		 *            Finds all roots of the polynomial in the coefficient
		 *            field.
		 *
		 * @details
		 *            The function iterates through all finite field elements
		 *            and evaluates them at the polynomial. If zero, the
		 *            root is added to <i>x</i>. The iteration stops if all
		 *            elements have been tested or if the number of found
		 *            roots agrees with the degree of the polynomial, in which
		 *            case no more roots can exist. The function returns the
		 *            number of roots found.
		 *            <br><br>
		 *            If <i>x</i> was not allocated to store all roots, the
		 *            function may run into undefined behavior. Note that
		 *            if <i>x</i> can hold at least as many elements as the
		 *            degree of the polynomial, it is guaranteed that <i>x</i>
		 *            can store all roots.
		 *
		 * @param x
		 *            Output array.
		 *
		 * @return
		 *            Number of roots found.
		 *
		 * @warning
		 *            If not sufficient memory could be allocated, the
		 *            function prints an error message to <code>stderr</code>
		 *            and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If the degree of the polynomial is smaller than 0,
		 *            the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <i>x</i> can not store all the roots of the
		 *            polynomial, the function runs into undefined behavior.
		 */
		int findRoots( uint32_t *x ) const;

		/**
		 * @brief
		 *            Computes the polynomial of minimal degree that
		 *            interpolates the tuples <i>(a[i],b[i])</i>
		 *            for <i>i=0,...,n-1</i>.
		 *
		 * @details
		 *            More precisely, the polynomial is set to the polynomial
		 *            \f$f(X)\f$ of degree \f$<n\f$ such that \f$f(a[i])=b[i]\f$
		 *            for \f$i=0,...,n-1\f$.
		 *            <br><br>
		 *            The implementation of this method corresponds to
		 *            <i>Lagrange interpolation</i>.
		 *
		 * @param a
		 *            Locators.
		 *
		 * @param b
		 *            Values of the polynomial at the locators.
		 *
		 * @param n
		 *            Number of valid locators and values listed in
		 *            <code>a</code> and <code>b</code>.
		 *
		 * @warning
		 *            If the first <code>n</code>elements in <code>a</code>
		 *            are not all pairwise distinct an error message
		 *            is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory can be provided for computing
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>a</code> or <code>b</code> contain invalid
		 *            elements of the polynomial's underlying finite field
		 *            then the method may run into unexpected undocumented
		 *            behavior.
		 */
		void interpolate
		( const uint32_t *a , const uint32_t *b , int n );

        /**
         * @brief
         *           Swaps the content of two polynomials such the
         *           the one represents the other.
         *
         * @param f
         *           Polynomial being assigned the fields of <i>g</i>.
         *
         * @param g
         *           Polynomial being assigned the fields of <i>f</i>.
         */
        inline static void swap
        ( SmallBinaryFieldPolynomial & f , SmallBinaryFieldPolynomial & g ) {

            f.swap(g);
        }

		/**
		 * @brief
		 *            Negates the polynomial <i>f</i> and stores the result
		 *            in <i>h</i>.
		 *
		 * @details
		 *            Because we are over a binary field, running the method
		 *            is equivalent to the expression
		 *            <pre>
		 *             h.assign(f);
		 *            </pre>
		 *            or
		 *            <pre>
		 *             h = f;
		 *            </pre>
		 *            however, for some template-like algorithms that also work
		 *            for non-binary fields, it is nonetheless better to have
		 *            a negation method.
		 *
		 * @param h
		 *            Output polynomial.
		 *
		 * @param f
		 *            Input polynomial.
		 *
		 * @warning
		 *            If not sufficient memory can be allocated to compute
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 */
		inline static void negate
		( SmallBinaryFieldPolynomial & h , const SmallBinaryFieldPolynomial & f ) {
			h.assign(f);
		}

		/**
		 * @brief
		 *            Multiplies the polynomial <i>g(X)</i> by the
		 *            monomial \f$X^n\f$ and stores the result
		 *            in <i>f</i>.
		 *
		 * @details
		 *            If <i>n>=0</i> the method lets
		 *            \f[
		 *             f(X)\leftarrow g(X)\cdot X^n.
		 *            \f]
		 *            If <i>n<0</i> the polynomial <i>f</i> will be set
		 *            as the quotient polynomial of <i>g</i> divided by
		 *            \f$X^n\f$ which is computed with the help of
		 *            <code>rightShift()</code>.
		 *
		 * @param f
		 *            Output polynomial.
		 *
		 * @param g
		 *            Input polynomial.
		 *
		 * @param n
		 *            The degree of the monomial that is multiplied with
		 *            <i>g</i>.
		 *
		 * @warning
		 *            If not sufficient memory could be provided, the method
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static void leftShift
		( SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g , int n );

		/**
		 * @brief
		 *            Computes the quotient of the polynomial <i>g(X)</i>
		 *            divided by the monomial \f$X^n\f$ and stores the result
		 *            in <i>f</i>.
		 *
		 * @details
		 *            For <i>n>=0</i> write
		 *            \f[
		 *             g(X)=\sum_{i=0}^d g_i\cdot X^i.
		 *            \f]
		 *            The method assign the polynomial <i>f</i> by
		 *            \f[
		 *             f(X)\leftarrow\sum_{i=n}^d g_i\cdot X^{i-n}
		 *            \f]
		 *            which corresponds to the result of the quotient that one
		 *            receives from polynomial division <i>g</i> by \f$X^n\f$.
		 *            Otherwise, if <i>n<0</i> the polynomial <i>f</i> is
		 *            assigned as
		 *            \f[
		 *             f(X)\leftarrow g(X)\cdot X^{-n}
		 *            \f]
		 *            using the <code>leftShift()</code> method.
		 *
		 * @param f
		 *            Output polynomial.
		 *
		 * @param g
		 *            Input polynomial.
		 *
		 * @param n
		 *            Number of right-shifted coefficients.
		 *
		 * @warning
		 *            If not sufficient memory could be provided, the method
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 */
		static void rightShift
		( SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g , int n );

		/**
		 * @brief
		 *            Computes the sum of two polynomials.
		 *
		 * @details
		 *            This methods stores the sum of <code>f</code>
		 *            and <code>g</code> in <code>h</code>.
		 *
		 * @param h
		 *            Will contain the sum of <code>f</code> and
		 *            <code>g</code>
		 *
		 * @param f
		 *            first summand
		 *
		 * @param g
		 *            second summand
		 *
		 * @warning
		 *            If it is not possible to ensure the capacity of
		 *            <code>h</code> to store the result an error message
		 *            is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		static void add
		( SmallBinaryFieldPolynomial & h ,
		  const SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g );

		/**
		 * @brief
		 *            Computes the difference of two polynomials.
		 *
		 * @details
		 *            This methods stores the difference of <code>f</code>
		 *            and <code>g</code> in <code>h</code>.
		 *            <br>
		 *            In fact, since the base field is a binary field, there
		 *            is no difference between addition and subtraction.
		 *            However, for the matter of consistency we provide
		 *            a subtraction method here.
		 *
		 * @param h
		 *            Will contain the difference of <code>f</code> and
		 *            <code>g</code>
		 *
		 * @param f
		 *            minuend
		 *
		 * @param g
		 *            subtrahend
		 *
		 * @warning
		 *            If it is not possible to ensure the capacity of
		 *            <code>h</code> to store the result an error message
		 *            is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		inline static void sub
		( SmallBinaryFieldPolynomial & h ,
		  const SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g ) {
			add(h,f,g);
		}

		/**
		 * @brief
		 *            Computes the product of two polynomials.
		 *
		 * @details
		 *            This methods stores the product of <code>f</code>
		 *            and <code>g</code> in <code>h</code>.
		 *
		 * @param h
		 *            Will contain the product of <code>f</code> and
		 *            <code>g</code>.
		 *
		 * @param f
		 *            First factor polynomial.
		 *
		 * @param g
		 *            Second factor polynomial
		 *
		 * @warning
		 *            If it is not possible to ensure the capacity of
		 *            <code>h</code> to store the result an error message
		 *            is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>h</code>, <code>f</code>, and <code>g</code>
		 *            do not all correspond to the same finite field (i.e. if
		 *            the pointer to the finite field instance differ) an error
		 *            message is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		static void mul
		( SmallBinaryFieldPolynomial & h ,
		  const SmallBinaryFieldPolynomial & f ,
		  const SmallBinaryFieldPolynomial & g );

		/**
		 * @brief
		 *            Multiplies a polynomial with a constant element of the
		 *            finite field.
		 *
		 * @details
		 *            This method multiplies the coefficients of the polynomial
		 *            <code>f</code> with the scalar <code>s</code> and stores
		 *            the result in <code>h</code>.
		 *
		 * @param h
		 *            Will contain the coefficients of <code>f</code>
		 *            multiplied with <code>s</code>.
		 *
		 * @param f
		 *            Input polynomial.
		 *
		 * @param s
		 *            Scalar which is an element of the polynomial's
		 *            <code>f</code> underlying finite field.
		 *
		 * @warning
		 *            If not sufficient memory can be provided to store the
		 *            coefficients of <code>f</code> in <code>h</code> an error
		 *            message is printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>s</code> is not a valid element of the
		 *            polynomial's <code>f</code> underlying field, the method
		 *            may run into undocumented behavior.
		 */
		static void mul
		( SmallBinaryFieldPolynomial & h ,
		  const SmallBinaryFieldPolynomial & f ,
		  uint32_t s );

		/**
		 * @brief
		 *            Performs Euclidean division of two polynomial to obtain
		 *            its quotient and remainder.
		 *
		 * @details
		 *            The Euclidean quotient \f$q(X)\f$ and the
		 *            Euclidean remainder \f$r(X)\f$ of two polynomials
		 *            \f$a(X)\f$ and \f$b(X)\f$ are defined to fulfill the
		 *            relation \f$q(X)\cdot b(X)=a(X)+r(X)\f$ where
		 *            \f$\deg(r(X))<\deg(b(X))\f$. If \f$b(X)\neq 0\f$ then
		 *            the Euclidean quotient and remainder are known to exist
		 *            and that they are unique.
		 *
		 * @param a
		 *            numerator
		 *
		 * @param b
		 *            denominator
		 *
		 * @param q
		 *            will contain the quotient of the Euclidean division of
		 *            <code>a</code> divided by <code>b</code>
		 *
		 * @param r
		 *            will contain the remainder of the Euclidean division
		 *            of <code>a</code> divided by <code>b</code>
		 *
		 * @warning
		 *            If <code>b</code> is the zero polynomial or if
		 *            <code>q</code> and <code>r</code> are of the same
		 *            reference then an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *             -1.
		 *
		 * @warning
		 *            If <code>q</code>, <code>r</code>, <code>a</code>, or
		 *            <code>b</code> are defined over different finite
		 *            fields (i.e. if their respective references are
		 *            different) an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If no sufficient memory could be provided to compute
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 */
		static void divRem
		( SmallBinaryFieldPolynomial & q ,
		  SmallBinaryFieldPolynomial & r ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b );

		/**
		 * @brief
		 *            Computes the quotient of two polynomials.
		 *
		 * @details
		 *            This method computes the quotient of <code>a</code>
		 *            divided by <code>b</code> and stores the result in
		 *            <code>q</code>. The quotient (i.e. Euclidean quotient)
		 *            is defined to be the polynomial <code>q</code> such that
		 *            <code>a=q*b+r</code> where the degree of <code>r</code>
		 *            is smaller than the degree of <code>b</code>.
		 *
		 * @param q
		 *            Will contain the result of the polynomial division
		 *            where <code>a</code> is divided by <code>b</code>.
		 *
		 * @param a
		 *            Numerator polynomial.
		 *
		 * @param b
		 *            Denominator polynomial.
		 *
		 * @warning
		 *            If <code>f</code> is the zero polynomial an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be provided to compute
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 */
		inline static void div
		( SmallBinaryFieldPolynomial & q ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b ) {
			SmallBinaryFieldPolynomial r(a.getField());
			divRem(q,r,a,b);
		}

		/**
		 * @brief
		 *            Computes the remainder of a polynomial division.
		 *
		 * @details
		 *            This method computes the remainder of <code>a</code>
		 *            divided by <code>b</code> and stores the result in
		 *            <code>r</code>. The remainder (i.e. Euclidean quotient)
		 *            is defined to be the polynomial <code>r</code> such that
		 *            <code>a=q*b+r</code> where the degree of <code>r</code>
		 *            is smaller than the degree of <code>b</code> and
		 *            <code>q</code> is also a polynomial.
		 *
		 * @param r
		 *            Will contain the remainder of the polynomial division
		 *            where <code>a</code> is divided by <code>b</code>.
		 *
		 * @param a
		 *            Numerator polynomial.
		 *
		 * @param b
		 *            Denominator polynomial.
		 *
		 * @warning
		 *            If <code>f</code> is the zero polynomial an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be provided to compute
		 *            the result an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 */
		inline static void rem
		( SmallBinaryFieldPolynomial & r ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b ) {
			SmallBinaryFieldPolynomial q(a.getField());

			divRem(q,r,a,b);
		}

        /**
         * @brief
         *          Computes the composition of two polynomials with
         *          coefficients in a small binary field.
         *
         * @details
         *          Given two polynomials <i>f</i> and <i>g</i> the composition
         *          <i>h=f(g)</i> is computed.
         *
         * @param h
         *          The composition of the two polynomials <i>f</i> and <i>g</i>
         *          on input, i.e., <i>f(g)</i>.
         *
         * @param f
         *          The outer polynomial.
         *
         * @param g
         *          The inner polynomial.
         *
         * @warning
         *          If not enough memory could be provided, an error
         *          message is printed to <code>stderr</code> and the
         *          program exits with status 'EXIT_FAILURE'.
         */
        static void eval
        ( SmallBinaryFieldPolynomial & h ,
          const SmallBinaryFieldPolynomial & f ,
          const SmallBinaryFieldPolynomial & g );

		/**
		 * @brief
		 *            Computes partial greatest common divisor of two
		 *            polynomials.
		 *
		 * @details
		 *            This methods performs the traditional extended Euclidean
		 *            algorithm, i.e. finds \f$s,t\f$ such that
		 *            \f$g=s\cdot a+t\cdot b\f$, until the degree of
		 *            \f$g\f$ becomes smaller than \f$d\f$.
		 *
		 * @param g
		 *            Will contain the partial greatest common divisor
		 *            of <i>a</i> and <i>b</i>.
		 *
		 * @param s see details
		 * @param t see details
		 * @param a see details
		 * @param b see details
		 * @param d see details
		 */
		static void pgcd
		( SmallBinaryFieldPolynomial & g ,
		  SmallBinaryFieldPolynomial & s ,
		  SmallBinaryFieldPolynomial & t ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b ,
		  int d );

		/**
		 * @brief
		 *            Applies the extended Euclidean algorithm to
		 *            two polynomials and stores the sequence of
		 *            the algorithm in vectors.
		 *
		 * @details
		 *            The method initializes
		 *            <pre>
		 *             g[0] = s[0]*a+t[0]*b
		 *             g[1] = s[1]*a+t[1]*b
		 *            </pre>
		 *            where <code>s[0]=1</code>, <code>t[0]=0</code>
		 *            and <code>s[1]=0</code>, and <code>t[1]=1</code>.
		 *            Then for each <code>i>=2</code> the method computes
		 *            <pre>
		 *             q[i] = g[i-2] rem g[i-1]
		 *            </pre>
		 *            and stores
		 *            <pre>
		 *             g[i] = g[i-2] - q[i]*g[i-1]
		 *             s[i] = s[i-2] - q[i]*s[s-1]
		 *             t[i] = t[i-2] - q[i]*t[i-1]
		 *            </pre>
		 *            in the output vectors passed to this method until
		 *            <code>g.back()</code> is a greatest common divisor
		 *            of <code>a</code> and <code>b</code>. In this way,
		 *            the vectors <code>g, s</code> and <code>t</code>
		 *            have the same length <code>n</code> and for each
		 *            <code>i=0,...,n-1</code> the relation
		 *            <pre>
		 *             g[i] = s[i]*a+t[i]*b
		 *            </pre>
		 *            holds.
		 *
		 * @attention
		 *            Note that, if <code>a</code> and <code>b</code> are
		 *            co-prime, the value <code>g.back()</code> does not
		 *            necessarily need to be equals 1 but is a non-zero
		 *            constant, i.e., <code>g.back.deg()==0</code>.
		 *
		 * @param g
		 *            Output vector of polynomials; see details.
		 *
		 * @param s
		 *            Output vector of polynomials; see details.
		 *
		 * @param t
		 *            Output vector of polynomials; see details.
		 *
		 * @param a
		 *            First input polynomial.
		 *
		 * @param b
		 *            Second input polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>g</code>, <code>s</code> or <code>t</code>
		 *            are vectors containing <code>a</code> or <code>b</code>
		 *            or if <code>g</code>, <code>s</code> and <code>t</code>
		 *            are not of pair-wise different reference, this method
		 *            runs into undocumented behavior.
		 */
		static void xgcd
		( std::vector<SmallBinaryFieldPolynomial> & g ,
		  std::vector<SmallBinaryFieldPolynomial> & s ,
		  std::vector<SmallBinaryFieldPolynomial> & t ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b );

		/**
		 * @brief
		 *            Applies the extended Euclidean algorithm to two
		 *            polynomials.
		 *
		 * @details
		 *            This method computes polynomials <code>g</code>,
		 *            <code>s</code> and <code>t</code> such that
		 *            <code>g</code> is a greatest common divisor of
		 *            <code>a</code> and <code>b</code> where
		 *            <pre>
		 *             g = s * a + t *b.
		 *            </pre>
		 *
		 * @attention
		 *            Note that, if <code>a</code> and <code>b</code> are
		 *            co-prime, the value <code>g</code> does not
		 *            necessarily need to be equals 1 but is a non-zero
		 *            constant, i.e., <code>g.deg()==0</code>.
		 *
		 * @param g
		 *            Output polynomial that is a greatest common divisor
		 *            of <code>a</code> and <code>b</code>.
		 *
		 * @param s
		 *            Output polynomial; see details above.
		 *
		 * @param t
		 *            Output polynomial; see details above.
		 *
		 * @param a
		 *            First input polynomial.
		 *
		 * @param b
		 *            Second input polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>g</code>, <code>s</code> and <code>t</code>
		 *            are not of pair-wise different reference, this method
		 *            runs into undocumented behavior.
		 */
		static void xgcd
		( SmallBinaryFieldPolynomial & g ,
		  SmallBinaryFieldPolynomial & s ,
		  SmallBinaryFieldPolynomial & t ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b );

		/**
		 * @brief
		 *            Computes a greatest common divisor of of two
		 *            polynomials.
		 *
		 * @attention
		 *            Note that, if <code>a</code> and <code>b</code> are
		 *            co-prime, the value <code>g</code> does not
		 *            necessarily need to be equals 1 but is a non-zero
		 *            constant, i.e., <code>g.deg()==0</code>.
		 *
		 * @param g
		 *            Output polynomial that represents a greatest common
		 *            divisor of <code>a</code> and <code>b</code>.
		 *
		 * @param a
		 *            First input polynomial.
		 *
		 * @param b
		 *            Second input polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		static void gcd
		( SmallBinaryFieldPolynomial & g ,
		  const SmallBinaryFieldPolynomial & a ,
		  const SmallBinaryFieldPolynomial & b );

		/**
		 * @brief
		 *            Checks whether two polynomials are co-prime.
		 *
		 * @details
		 *            Two polynomials <code>a</code> and <code>b</code>
		 *            are said to be co-prime if the only non-zero
		 *            polynomials <code>c</code> that divide <code>a</code>
		 *            and <code>b</code> are constants.
		 *
		 * @param a
		 *            First polynomial.
		 *
		 * @param b
		 *            Second polynomial.
		 *
		 * @return
		 *            <code>true</code> if <code>a</code> and <code>b</code>
		 *            are co-prime and, otherwise, <code>false</code>.
		 *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		static bool areCoprime
		( const SmallBinaryFieldPolynomial & a , const SmallBinaryFieldPolynomial & b );

		/**
		 * @brief
		 *            Computes the modular inverse of a polynomial.
		 *
		 * @details
		 *            Given two polynomials \f$x\f$ and \f$m\f$, this function
		 *            determines a polynomial \f$y\f$ with
		 *            \f$\deg(y)<\deg(m)\f$ such that
		 *            \f$(y\cdot x)~mod~m=1\f$. If such a polynomial exists,
		 *            this method will succeed in finding such a polynomial
		 *            and returns <code>true</code>; otherwise, if such a
		 *            polynomial does not exist, <code>y</code> will be
		 *            left unchanged and the function returns
		 *            <code>false</code>.
		 *
		 * @param y
		 *            Output polynomial that forms the inverse of <code>x</code>
		 *            modulo <code>m</code>.
		 *
		 * @param x
		 *            Input polynomial of which the modular inverse is
		 *            computed.
		 *
		 * @param m
		 *            Input modulus polynomial.
		 *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		static bool invMod
		( SmallBinaryFieldPolynomial & y , const SmallBinaryFieldPolynomial & x ,
		  const SmallBinaryFieldPolynomial & m );

		/**
		 * @brief
		 *            Tests this polynomial on equality with another
		 *            polynomial.
		 *
		 * @details
		 *            The implementation of this function allows for
		 *            two polynomials
		 *            <pre>
		 *             SmallBinaryFieldPolynomial a , b;
		 *            </pre>
		 *            to compare them via statements of the form
		 *            <pre>
		 *             if ( a == b ) {
		 *                // then do something
		 *             }
		 *            </pre>
		 *
		 * @param g
		 *            The polynomial with which this polynomial is
		 *            compared.
		 *
		 * @return
		 *            <code>true</code> if this polynomial represents the
		 *            same as the polynomial represented by <code>g</code>
		 *            and, otherwise, <code>false</code>.
		 */
        inline bool operator==( const SmallBinaryFieldPolynomial & g ) const {
            return equals(g);
        }

		/**
		 * @brief
		 *            Tests this polynomial on inequality with another
		 *            polynomial.
		 *
		 * @details
		 *            The implementation of this function allows for
		 *            two polynomials
		 *            <pre>
		 *             SmallBinaryFieldPolynomial a , b;
		 *            </pre>
		 *            to compare them via statements of the form
		 *            <pre>
		 *             if ( a != b ) {
		 *                // then do something
		 *             }
		 *            </pre>
		 *
		 * @param g
		 *            The polynomial with which this polynomial is
		 *            compared.
		 *
		 * @return
		 *            <code>true</code> if this polynomial represents another
		 *            polynomial as the polynomial represented by
		 *            <code>g</code> and, otherwise, <code>false</code>.
		 */
        inline bool operator!=( const SmallBinaryFieldPolynomial & g ) const {
            return !equals(g);
        }

        /**
         * @brief
         *            Computes the sum of this polynomial with another
         *            polynomial and returns the result.
         *
         * @details
         *            The implementation of this function allows to compute
         *            the sum of two polynomials
         *            <pre>
         *             SmallBinaryFieldPolynomial a , b;
         *            </pre>
         *            via statements as
         *            <pre>
         *             SmallBinaryFieldPolynomial c = a + b;
         *            </pre>
         *            and is wrapped around the static
         *            \link SmallBinaryFieldPolynomial::add()\endlink
         *            method.
         *
         * @param g
         *            The polynomial added to this polynomial.
         *
         * @return
         *            The sum of this polynomial and <code>g</code>.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial operator+
        ( const SmallBinaryFieldPolynomial & g ) const {

            SmallBinaryFieldPolynomial h(this->gfPtr[0]);
            add(h,*this,g);
            return h;
        }

        /**
         * @brief
         *            Computes the difference of this polynomial with another
         *            polynomial and returns the result.
         *
         * @details
         *            The implementation of this function allows to compute
         *            the sum of two polynomials
         *            <pre>
         *             SmallBinaryFieldPolynomial a , b;
         *            </pre>
         *            via statements as
         *            <pre>
         *             SmallBinaryFieldPolynomial c = a - b;
         *            </pre>
         *            and is wrapped around the static
         *            \link SmallBinaryFieldPolynomial::sub()\endlink
         *            method.
         *
         * @param g
         *            The polynomial subtracted from this polynomial.
         *
         * @return
         *            The difference between this polynomial and
         *            <code>g</code>.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial operator-
        ( const SmallBinaryFieldPolynomial & g ) const {

        	SmallBinaryFieldPolynomial h(this->gfPtr[0]);
            sub(h,*this,g);
            return h;
        }

        /**
         * @brief
         *            Computes the Euclidean quotient of this polynomial
         *            (as numerator) divided by another polynomial
         *            (as denominator) and returns the result.
         *
         * @details
         *            The implementation of this function allows to compute
         *            the quotient between two polynomials
         *            <pre>
         *             SmallBinaryFieldPolynomial a , b;
         *            </pre>
         *            via statements as
         *            <pre>
         *             SmallBinaryFieldPolynomial c = a / b;
         *            </pre>
         *            and is wrapped around the static
         *            \link SmallBinaryFieldPolynomial::div()\endlink
         *            method.
         *
         * @param g
         *            The denominator polynomial.
         *
         * @return
         *            The quotient of this polynomial divided by
         *            <code>g</code>.
         *
         * @warning
         *            If <code>g</code> is the zero polynomial, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial operator/
        ( const SmallBinaryFieldPolynomial & g ) const {

            SmallBinaryFieldPolynomial h(this->gfPtr[0]);
            div(h,*this,g);
            return h;
        }

        /**
         * @brief
         *            Computes the product of this polynomial with another
         *            polynomial and returns the result.
         *
         * @details
         *            The implementation of this function allows to compute
         *            the product of two polynomials
         *            <pre>
         *             SmallBinaryFieldPolynomial a , b;
         *            </pre>
         *            via statements as
         *            <pre>
         *             SmallBinaryFieldPolynomial c = a * b;
         *            </pre>
         *            and is wrapped around the static
         *            \link SmallBinaryFieldPolynomial::mul()\endlink
         *            method.
         *
         * @param g
         *            The polynomial multiplied to this polynomial.
         *
         * @return
         *            The product of this polynomial and <code>g</code>.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial operator*
        ( const SmallBinaryFieldPolynomial & g ) const {

            SmallBinaryFieldPolynomial h(this->gfPtr[0]);
            mul(h,*this,g);
            return h;
        }

        /**
         * @brief
         *            Computes the Euclidean remainder of this polynomial
         *            (as numerator) divided by another polynomial
         *            (as denominator) and returns the result.
         *
         * @details
         *            The implementation of this function allows to compute
         *            the Euclidean remainders between two polynomials
         *            <pre>
         *             SmallBinaryFieldPolynomial a , b;
         *            </pre>
         *            via statements as
         *            <pre>
         *             SmallBinaryFieldPolynomial c = a % b;
         *            </pre>
         *            and is wrapped around the static
         *            \link SmallBinaryFieldPolynomial::rem()\endlink
         *            method.
         *
         * @param g
         *            The denominator polynomial.
         *
         * @return
         *            The remainder of this polynomial modulo
         *            <code>g</code>.
         *
         * @warning
         *            If <code>g</code> is the zero polynomial, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial operator%
        ( const SmallBinaryFieldPolynomial & g ) const {

            SmallBinaryFieldPolynomial h(this->gfPtr[0]);
            rem(h,*this,g);
            return h;
        }

        /**
         * @brief
         *             Adds a polynomial to this polynomial inplace.
         *
         * @details
         *             The implementation of this functions allows to
         *             increase a polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial a;
         *             </pre>
         *             by another polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial b;
         *             </pre>
         *             using statements as
         *             <pre>
         *              a += b;
         *             </pre>
         *             and is wrapped around the
         *             static \link SmallBinaryFieldPolynomial::add()\endlink
         *             method.
         *
         * @param f
         *             The polynomial added to this polynomial.
         *
         * @return
         *             Reference to this polynomial (after increment).
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial & operator+=
        		( const SmallBinaryFieldPolynomial & f ) {

            add(*this,*this,f);
            return *this;
        }

        /**
         * @brief
         *             Subtracts a polynomial from this polynomial inplace.
         *
         * @details
         *             The implementation of this functions allows to
         *             decrease a polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial a;
         *             </pre>
         *             by another polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial b;
         *             </pre>
         *             using statements as
         *             <pre>
         *              a -= b;
         *             </pre>
         *             and is wrapped around the
         *             static \link SmallBinaryFieldPolynomial::sub()\endlink
         *             method.
         *
         * @param f
         *             The polynomial subtracted from this polynomial.
         *
         * @return
         *             Reference to this polynomial (after decrement).
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial & operator-=
        		( const SmallBinaryFieldPolynomial & f ) {

            sub(*this,*this,f);
            return *this;
        }

        /**
         * @brief
         *             Multiplies a polynomial with this polynomial inplace.
         *
         * @details
         *             The implementation of this functions allows to
         *             multiply a polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial a;
         *             </pre>
         *             by another polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial b;
         *             </pre>
         *             using statements as
         *             <pre>
         *              a *= b;
         *             </pre>
         *             and is wrapped around the
         *             static \link SmallBinaryFieldPolynomial::mul()\endlink
         *             method.
         *
         * @param f
         *             The polynomial multiplied to this polynomial.
         *
         * @return
         *             Reference to this polynomial (after multiplication).
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial & operator*=( const SmallBinaryFieldPolynomial & f ) {
            mul(*this,*this,f);
            return *this;
        }

        /**
         * @brief
         *             Performs inplace polynomial division.
         *
         * @details
         *             The implementation of this functions allows to
         *             divide a polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial a;
         *             </pre>
         *             by another polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial b;
         *             </pre>
         *             using statements as
         *             <pre>
         *              a /= b;
         *             </pre>
         *             and is wrapped around the
         *             static \link SmallBinaryFieldPolynomial::div()\endlink
         *             method.
         *
         * @param f
         *             The denominator polynomial.
         *
         * @return
         *             Reference to this polynomial (after division).
         *
         * @warning
         *            If <code>f</code> represents the zero polynomial, an
         *            error message is written to <code>stderr</code> and
         *            the program exits with status 'EXIT_FAILURE'.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial & operator/=
        		( const SmallBinaryFieldPolynomial & f ) {

            div(*this,*this,f);
            return *this;
        }

        /**
         * @brief
         *             Performs inplace polynomial remaindering.
         *
         * @details
         *             The implementation of this functions allows to
         *             replace a polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial a;
         *             </pre>
         *             by its remainder modulo another polynomial
         *             <pre>
         *              SmallBinaryFieldPolynomial b;
         *             </pre>
         *             using statements as
         *             <pre>
         *              a %= b;
         *             </pre>
         *             and is wrapped around the
         *             static \link SmallBinaryFieldPolynomial::rem()\endlink
         *             method.
         *
         * @param f
         *             The denominator polynomial.
         *
         * @return
         *             Reference to this polynomial (after remaindering).
         *
         * @warning
         *            If <code>f</code> represents the zero polynomial, an
         *            error message is written to <code>stderr</code> and
         *            the program exits with status 'EXIT_FAILURE'.
         *
		 * @warning
		 *            If not enough memory could be provided, an
		 *            error message is printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
         */
        inline SmallBinaryFieldPolynomial & operator%=
        		( const SmallBinaryFieldPolynomial & f ) {

            rem(*this,*this,f);
            return *this;
        }
	};

	/**
	 * @brief
	 *            Negates the polynomial <i>f</i> and stores the result
	 *            in <i>h</i>.
	 *
	 * @details
	 *            Because we are over a binary field, running the method
	 *            is equivalent to the expression
	 *            <pre>
	 *             h.assign(f);
	 *            </pre>
	 *            or
	 *            <pre>
	 *             h = f;
	 *            </pre>
	 *            however, for some template-like algorithms that also work
	 *            for non-binary fields, it is nonetheless better to have
	 *            a negation method.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::negate
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::negate(f,g);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            negate(f,g);
	 *            </pre>
	 *
	 * @param h
	 *            Output polynomial.
	 *
	 * @param f
	 *            Input polynomial.
	 *
	 * @warning
	 *            If not sufficient memory can be allocated to compute
	 *            the result an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 */
	inline void negate
	( SmallBinaryFieldPolynomial & h , const SmallBinaryFieldPolynomial & f ) {
		SmallBinaryFieldPolynomial::negate(h,f);
	}

	/**
	 * @brief
	 *            Multiplies the polynomial <i>g(X)</i> by the
	 *            monomial \f$X^n\f$ and stores the result in <i>f</i>.
	 *
	 * @details
	 *            If <i>n>=0</i> the method lets
	 *            \f[
	 *             f(X)\leftarrow g(X)\cdot X^n.
	 *            \f]
	 *            If <i>n<0</i> the polynomial <i>f</i> will be set
	 *            as the quotient polynomial of <i>g</i> divided by
	 *            \f$X^n\f$ which is computed with the help of
	 *            <code>rightShift()</code>.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::leftShift
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,int)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::leftShift(f,g,n);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            leftShift(f,g,n);
	 *            </pre>
	 *
	 * @param f
	 *            Output polynomial.
	 *
	 * @param g
	 *            Input polynomial.
	 *
	 * @param n
	 *            The degree of the monomial that is multiplied with
	 *            <i>g</i>.
	 *
	 * @warning
	 *            If not sufficient memory could be provided, the method
	 *            prints an error message to <code>stderr</code> and exits
	 *            with status 'EXIT_FAILURE'.
	 */
	inline void leftShift
	( SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g , int n ) {
		SmallBinaryFieldPolynomial::leftShift(f,g,n);
	}

	/**
	 * @brief
	 *            Computes the quotient of the polynomial <i>g(X)</i>
	 *            divided by the monomial \f$X^n\f$ and stores the result
	 *            in <i>f</i>.
	 *
	 * @details
	 *            For <i>n>=0</i> write
	 *            \f[
	 *             g(X)=\sum_{i=0}^d g_i\cdot X^i.
	 *            \f]
	 *            The method assign the polynomial <i>f</i> by
	 *            \f[
	 *             f(X)\leftarrow\sum_{i=n}^d g_i\cdot X^{i-n}
	 *            \f]
	 *            which corresponds to the result of the quotient that one
	 *            receives from polynomial division <i>g</i> by \f$X^n\f$.
	 *            Otherwise, if <i>n<0</i> the polynomial <i>f</i> is
	 *            assigned as
	 *            \f[
	 *             f(X)\leftarrow g(X)\cdot X^{-n}
	 *            \f]
	 *            using the <code>leftShift()</code> method.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::rightShift
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,int)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::rightShift(f,g,n);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            rightShift(f,g,n);
	 *            </pre>
	 *
	 * @param f
	 *            Output polynomial.
	 *
	 * @param g
	 *            Input polynomial.
	 *
	 * @param n
	 *            Number of right-shifted coefficients.
	 *
	 * @warning
	 *            If not sufficient memory could be provided, the method
	 *            prints an error message to <code>stderr</code> and exits
	 *            with status 'EXIT_FAILURE'.
	 */
	inline void rightShift
	( SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g , int n ) {
		SmallBinaryFieldPolynomial::rightShift(f,g,n);
	}

	/**
	 * @brief
	 *            Computes the sum of two polynomials.
	 *
	 * @details
	 *            This methods stores the sum of <code>f</code>
	 *            and <code>g</code> in <code>h</code>.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::add
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::add(h,f,g);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            add(h,f,g);
	 *            </pre>
	 *
	 * @param h
	 *            Will contain the sum of <code>f</code> and
	 *            <code>g</code>
	 *
	 * @param f
	 *            first summand
	 *
	 * @param g
	 *            second summand
	 *
	 * @warning
	 *            If it is not possible to ensure the capacity of
	 *            <code>h</code> to store the result an error message
	 *            is printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 */
	inline void add
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		SmallBinaryFieldPolynomial::add(h,f,g);
	}

	/**
	 * @brief
	 *            Computes the difference of two polynomials.
	 *
	 * @details
	 *            This methods stores the difference of <code>f</code>
	 *            and <code>g</code> in <code>h</code>.
	 *            <br>
	 *            In fact, since the base field is a binary field, there
	 *            is no difference between addition and subtraction.
	 *            However, for the matter of consistency we provide
	 *            a subtraction method here.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::sub
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::sub(h,f,g);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            sub(h,f,g);
	 *            </pre>
	 *
	 * @param h
	 *            Will contain the difference of <code>f</code> and
	 *            <code>g</code>
	 *
	 * @param f
	 *            minuend
	 *
	 * @param g
	 *            subtrahend
	 *
	 * @warning
	 *            If it is not possible to ensure the capacity of
	 *            <code>h</code> to store the result an error message
	 *            is printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 */
	inline void sub
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		SmallBinaryFieldPolynomial::sub(h,f,g);
	}

	/**
	 * @brief
	 *            Computes the product of two polynomials.
	 *
	 * @details
	 *            This methods stores the product of <code>f</code>
	 *            and <code>g</code> in <code>h</code>.
	 *
	 * @param h
	 *            Will contain the product of <code>f</code> and
	 *            <code>g</code>.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::mul
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::mul(h,f,g);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            mul(h,f,g);
	 *            </pre>
	 *
	 * @param f
	 *            First factor polynomial.
	 *
	 * @param g
	 *            Second factor polynomial
	 *
	 * @warning
	 *            If it is not possible to ensure the capacity of
	 *            <code>h</code> to store the result an error message
	 *            is printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>h</code>, <code>f</code>, and <code>g</code>
	 *            do not all correspond to the same finite field (i.e. if
	 *            the pointer to the finite field instance differ) an error
	 *            message is printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 */
	inline void mul
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  const SmallBinaryFieldPolynomial & g ) {

		SmallBinaryFieldPolynomial::mul(h,f,g);
	}

	/**
	 * @brief
	 *            Multiplies a polynomial with a constant element of the
	 *            finite field.
	 *
	 * @details
	 *            This method multiplies the coefficients of the polynomial
	 *            <code>f</code> with the scalar <code>s</code> and stores
	 *            the result in <code>h</code>.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::mul
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 uint32_t)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::mul(h,f,s);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            mul(h,f,s);
	 *            </pre>
	 *
	 * @param h
	 *            Will contain the coefficients of <code>f</code>
	 *            multiplied with <code>s</code>.
	 *
	 * @param f
	 *            Input polynomial.
	 *
	 * @param s
	 *            Scalar which is an element of the polynomial's
	 *            <code>f</code> underlying finite field.
	 *
	 * @warning
	 *            If not sufficient memory can be provided to store the
	 *            coefficients of <code>f</code> in <code>h</code> an error
	 *            message is printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>s</code> is not a valid element of the
	 *            polynomial's <code>f</code> underlying field, the method
	 *            may run into undocumented behavior.
	 */
	inline void mul
	( SmallBinaryFieldPolynomial & h ,
	  const SmallBinaryFieldPolynomial & f ,
	  uint32_t s ) {

		SmallBinaryFieldPolynomial::mul(h,f,s);
	}

	/**
	 * @brief
	 *            Performs Euclidean division of two polynomial to obtain
	 *            its quotient and remainder.
	 *
	 * @details
	 *            The Euclidean quotient \f$q(X)\f$ and the
	 *            Euclidean remainder \f$r(X)\f$ of two polynomials
	 *            \f$a(X)\f$ and \f$b(X)\f$ are defined to fulfill the
	 *            relation \f$q(X)\cdot b(X)=a(X)+r(X)\f$ where
	 *            \f$\deg(r(X))<\deg(b(X))\f$. If \f$b(X)\neq 0\f$ then
	 *            the Euclidean quotient and remainder are known to exist
	 *            and that they are unique.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::divRem
	 *            	(SmallBinaryFieldPolynomial&,SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::divRem(q,r,a,b);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            divRem(q,r,a,b);
	 *            </pre>
	 *
	 * @param a
	 *            numerator
	 *
	 * @param b
	 *            denominator
	 *
	 * @param q
	 *            will contain the quotient of the Euclidean division of
	 *            <code>a</code> divided by <code>b</code>
	 *
	 * @param r
	 *            will contain the remainder of the Euclidean division
	 *            of <code>a</code> divided by <code>b</code>
	 *
	 * @warning
	 *            If <code>b</code> is the zero polynomial or if
	 *            <code>q</code> and <code>r</code> are of the same
	 *            reference then an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *             -1.
	 *
	 * @warning
	 *            If <code>q</code>, <code>r</code>, <code>a</code>, or
	 *            <code>b</code> are defined over different finite
	 *            fields (i.e. if their respective references are
	 *            different) an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 * @warning
	 *            If no sufficient memory could be provided to compute
	 *            the result an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 */
	inline void divRem
	( SmallBinaryFieldPolynomial & q ,
	  SmallBinaryFieldPolynomial & r ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {
		SmallBinaryFieldPolynomial::divRem(q,r,a,b);
	}

	/**
	 * @brief
	 *            Computes the quotient of two polynomials.
	 *
	 * @details
	 *            This method computes the quotient of <code>a</code>
	 *            divided by <code>b</code> and stores the result in
	 *            <code>q</code>. The quotient (i.e. Euclidean quotient)
	 *            is defined to be the polynomial <code>q</code> such that
	 *            <code>a=q*b+r</code> where the degree of <code>r</code>
	 *            is smaller than the degree of <code>b</code>.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::div
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::div(q,a,b);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            div(q,a,b);
	 *            </pre>
	 *
	 * @param q
	 *            Will contain the result of the polynomial division
	 *            where <code>a</code> is divided by <code>b</code>.
	 *
	 * @param a
	 *            Numerator polynomial.
	 *
	 * @param b
	 *            Denominator polynomial.
	 *
	 * @warning
	 *            If <code>f</code> is the zero polynomial an error message
	 *            will be printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be provided to compute
	 *            the result an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 */
	inline void div
	( SmallBinaryFieldPolynomial & q ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial::div(q,a,b);
	}

	/**
	 * @brief
	 *            Computes the remainder of a polynomial division.
	 *
	 * @details
	 *            This method computes the remainder of <code>a</code>
	 *            divided by <code>b</code> and stores the result in
	 *            <code>r</code>. The remainder (i.e. Euclidean quotient)
	 *            is defined to be the polynomial <code>r</code> such that
	 *            <code>a=q*b+r</code> where the degree of <code>r</code>
	 *            is smaller than the degree of <code>b</code> and
	 *            <code>q</code> is also a polynomial.
	 *            <br><br>
	 *            This method is a convenience method wrapped around the
	 *            <code>
	 *              SmallBinaryFieldPolynomial::rem
	 *            	(SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&,
	 *            	 const SmallBinaryFieldPolynomial&)
	 *            </code>
	 *            method. Thus, instead of frequently running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial::rem(q,a,b);
	 *            </pre>
	 *            we may simply run
	 *            <pre>
	 *            rem(q,a,b);
	 *            </pre>
	 *
	 * @param r
	 *            Will contain the remainder of the polynomial division
	 *            where <code>a</code> is divided by <code>b</code>.
	 *
	 * @param a
	 *            Numerator polynomial.
	 *
	 * @param b
	 *            Denominator polynomial.
	 *
	 * @warning
	 *            If <code>f</code> is the zero polynomial an error message
	 *            will be printed to <code>stderr</code> and the program
	 *            exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If not sufficient memory could be provided to compute
	 *            the result an error message is printed to
	 *            <code>stderr</code> and the program exits with status
	 *            -1.
	 *
	 */
	inline void rem
	( SmallBinaryFieldPolynomial & r ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial::rem(r,a,b);
	}

    /**
     * @brief
     *          Computes the composition of two polynomials with
     *          coefficients in a small binary field.
     *
     * @details
     *          Given two polynomials <i>f</i> and <i>g</i> the composition
     *          <i>h=f(g)</i> is computed.
     *
     * @param h
     *          The composition of the two polynomials <i>f</i> and <i>g</i>
     *          on input, i.e., <i>f(g)</i>.
     *
     * @param f
     *          The outer polynomial.
     *
     * @param g
     *          The inner polynomial.
     *
     * @warning
     *          If not enough memory could be provided, an error
     *          message is printed to <code>stderr</code> and the
     *          program exits with status 'EXIT_FAILURE'.
     */
    inline void eval
    ( SmallBinaryFieldPolynomial & h ,
      const SmallBinaryFieldPolynomial & f ,
      const SmallBinaryFieldPolynomial & g ) {

        SmallBinaryFieldPolynomial::eval(h,f,g);
    }

	/**
	 * @brief
	 *            Computes partial greatest common divisor of two
	 *            polynomials.
	 *
	 * @details
	 *            This methods performs the traditional extended Euclidean
	 *            algorithm, i.e. finds \f$s,t\f$ such that
	 *            \f$g=s\cdot a+t\cdot b\f$, until the degree of
	 *            \f$g\f$ becomes smaller than \f$d\f$.
	 *
	 * @param g
	 *            Will contain the partial greatest common divisor
	 *            of <i>a</i> and <i>b</i>.
	 *
	 * @param s see details
	 * @param t see details
	 * @param a see details
	 * @param b see details
	 * @param d see details
	 */
	inline void pgcd
	( SmallBinaryFieldPolynomial & g ,
	  SmallBinaryFieldPolynomial & s ,
	  SmallBinaryFieldPolynomial & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ,
	  int d ) {
		SmallBinaryFieldPolynomial::pgcd(g,s,t,a,b,d);
	}

	/**
	 * @brief
	 *            Applies the extended Euclidean algorithm to
	 *            two polynomials and stores the sequence of
	 *            the algorithm in vectors.
	 *
	 * @details
	 *            The method initializes
	 *            <pre>
	 *             g[0] = s[0]*a+t[0]*b
	 *             g[1] = s[1]*a+t[1]*b
	 *            </pre>
	 *            where <code>s[0]=1</code>, <code>t[0]=0</code>
	 *            and <code>s[1]=0</code>, and <code>t[1]=1</code>.
	 *            Then for each <code>i>=2</code> the method computes
	 *            <pre>
	 *             q[i] = g[i-2] rem g[i-1]
	 *            </pre>
	 *            and stores
	 *            <pre>
	 *             g[i] = g[i-2] - q[i]*g[i-1]
	 *             s[i] = s[i-2] - q[i]*s[s-1]
	 *             t[i] = t[i-2] - q[i]*t[i-1]
	 *            </pre>
	 *            in the output vectors passed to this method until
	 *            <code>g.back()</code> is a greatest common divisor
	 *            of <code>a</code> and <code>b</code>. In this way,
	 *            the vectors <code>g, s</code> and <code>t</code>
	 *            have the same length <code>n</code> and for each
	 *            <code>i=0,...,n-1</code> the relation
	 *            <pre>
	 *             g[i] = s[i]*a+t[i]*b
	 *            </pre>
	 *            holds.
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g.back()</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.back.deg()==0</code>.
	 *
	 * @param g
	 *            Output vector of polynomials; see details.
	 *
	 * @param s
	 *            Output vector of polynomials; see details.
	 *
	 * @param t
	 *            Output vector of polynomials; see details.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>g</code>, <code>s</code> or <code>t</code>
	 *            are vectors containing <code>a</code> or <code>b</code>
	 *            or if <code>g</code>, <code>s</code> and <code>t</code>
	 *            are not of pair-wise different reference, this method
	 *            runs into undocumented behavior.
	 */
	inline void xgcd
	( std::vector<SmallBinaryFieldPolynomial> & g ,
	  std::vector<SmallBinaryFieldPolynomial> & s ,
	  std::vector<SmallBinaryFieldPolynomial> & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {
		SmallBinaryFieldPolynomial::xgcd(g,s,t,a,b);
	}

	/**
	 * @brief
	 *            Applies the extended Euclidean algorithm to two
	 *            polynomials.
	 *
	 * @details
	 *            This method computes polynomials <code>g</code>,
	 *            <code>s</code> and <code>t</code> such that
	 *            <code>g</code> is a greatest common divisor of
	 *            <code>a</code> and <code>b</code> where
	 *            <pre>
	 *             g = s * a + t *b.
	 *            </pre>
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.deg()==0</code>.
	 *
	 * @param g
	 *            Output polynomial that is a greatest common divisor
	 *            of <code>a</code> and <code>b</code>.
	 *
	 * @param s
	 *            Output polynomial; see details above.
	 *
	 * @param t
	 *            Output polynomial; see details above.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 *
	 * @warning
	 *            If <code>g</code>, <code>s</code> and <code>t</code>
	 *            are not of pair-wise different reference, this method
	 *            runs into undocumented behavior.
	 */
	inline void xgcd
	( SmallBinaryFieldPolynomial & g ,
	  SmallBinaryFieldPolynomial & s ,
	  SmallBinaryFieldPolynomial & t ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		SmallBinaryFieldPolynomial::xgcd(g,s,t,a,b);
	}

	/**
	 * @brief
	 *            Computes a greatest common divisor of of two
	 *            polynomials.
	 *
	 * @attention
	 *            Note that, if <code>a</code> and <code>b</code> are
	 *            co-prime, the value <code>g</code> does not
	 *            necessarily need to be equals 1 but is a non-zero
	 *            constant, i.e., <code>g.deg()==0</code>.
	 *
	 * @param g
	 *            Output polynomial that represents a greatest common
	 *            divisor of <code>a</code> and <code>b</code>.
	 *
	 * @param a
	 *            First input polynomial.
	 *
	 * @param b
	 *            Second input polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	inline void gcd
	( SmallBinaryFieldPolynomial & g ,
	  const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {
		SmallBinaryFieldPolynomial::gcd(g,a,b);
	}

	/**
	 * @brief
	 *            Checks whether two polynomials are co-prime.
	 *
	 * @details
	 *            Two polynomials <code>a</code> and <code>b</code>
	 *            are said to be co-prime if the only non-zero
	 *            polynomials <code>c</code> that divide <code>a</code>
	 *            and <code>b</code> are constants.
	 *
	 * @param a
	 *            First polynomial.
	 *
	 * @param b
	 *            Second polynomial.
	 *
	 * @return
	 *            <code>true</code> if <code>a</code> and <code>b</code>
	 *            are co-prime and, otherwise, <code>false</code>.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	inline bool areCoprime
	( const SmallBinaryFieldPolynomial & a ,
	  const SmallBinaryFieldPolynomial & b ) {

		return SmallBinaryFieldPolynomial::areCoprime(a,b);
	}

	/**
	 * @brief
	 *            Computes the modular inverse of a polynomial.
	 *
	 * @details
	 *            Given two polynomials \f$x\f$ and \f$m\f$, this function
	 *            determines a polynomial \f$y\f$ with
	 *            \f$\deg(y)<\deg(m)\f$ such that
	 *            \f$(y\cdot x)~mod~m=1\f$. If such a polynomial exists,
	 *            this method will succeed in finding such a polynomial
	 *            and returns <code>true</code>; otherwise, if such a
	 *            polynomial does not exist, <code>y</code> will be
	 *            left unchanged and the function returns
	 *            <code>false</code>.
	 *
	 * @param y
	 *            Output polynomial that forms the inverse of <code>x</code>
	 *            modulo <code>m</code>.
	 *
	 * @param x
	 *            Input polynomial of which the modular inverse is
	 *            computed.
	 *
	 * @param m
	 *            Input modulus polynomial.
	 *
	 * @warning
	 *            If not enough memory could be provided, an
	 *            error message is printed to <code>stderr</code>
	 *            and the program exits with status 'EXIT_FAILURE'.
	 */
	inline bool invMod
	( SmallBinaryFieldPolynomial & y ,
	  const SmallBinaryFieldPolynomial & x ,
	  const SmallBinaryFieldPolynomial & m ) {
		return SmallBinaryFieldPolynomial::invMod(y,x,m);
	}

	/**
	 * @brief
	 *           Prints a text representation of a polynomial to the specified
	 *           output stream.
	 *
	 * @details
	 *           This functions allows to use statements as
	 *           <pre>
	 *            SmallBinaryFieldPolynomial f = ...; //Initialize
	 *
	 *            cout << f << endl;
	 *           </pre>
	 *
	 * @param out
	 *           The output stream to where the text representation of the
	 *           polynomial is printed
	 *
	 * @param f
	 *           The polynomial of which the text representation is printed to
	 *           <code>out</code>.
	 *
	 * @return
	 *           The reference to <code>out</code>
	 */
	THIMBLE_DLL std::ostream & operator<<
			( std::ostream & out , const SmallBinaryFieldPolynomial & f );
}

#endif /* THIMBLE_SMALLBINARYFIELDPOLYNOMIAL_H_ */
