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
 * @file SmallBinaryField.h
 *
 * @brief
 *            Provides functionalities enabling arithmetic in binary galois
 *            fields that are of small cardinality.
 *
 * @details
 * @section sec_finitefields Finite Fields
 *
 * THIMBLE provides the \link thimble::SmallBinaryField
 * SmallBinaryField\endlink class that enables arithmetic in finite fields
 * being of characteristic 2 having small degree. Theoretically, binary
 * fields of degree up to 31 can be processed. However, to accelerate
 * finite field arithmetic a \link thimble::SmallBinaryField
 * SmallBinaryField\endlink object precomputes tables of size of orders
 * the size of the finite field. Consequently, reasonable finite fields
 * represented by a \link thimble::SmallBinaryField SmallBinaryField\endlink
 * object are of degree up to 20 while a realistic upper bound on the degree
 * is 25 (for a computer with sufficient RAM).
 *
 * To construct a finite field of degree 16, we may run
 * <pre>
 *  SmallBinaryField gf(16);
 * </pre>
 * Elements of the finite field are represented by unsigned 32 bit intgers
 * of type <code>uint32_t</code> provided by the <code>stdint.h</code> header.
 * For a finite field element
 * <pre>
 *  uint32_t a;
 * </pre>
 * to represent a well-defined element from <code>gf</code> it must range
 * between \f$0\f$ and \f$2^{16}-1\f$; more generally, if the finite
 * field is of degree \f$d\f$, a valid element is represented by
 * an (unsigned) integer between \f$0\f$ and \f$2^d-1\f$. Thereby, it is
 * guaranteed that the zero element from a binary field is represented by
 * the integer 0 while the one element is represented by the integer 1.
 *
 * @subsection sec_finitefields_addsub Addition/Subtraction
 *
 * In the following we want to perform arithmetic in a binary finite
 * field
 * <pre>
 *  SmallBinaryField gf(d);
 * </pre>
 * of degree <code>d</code> represented by a \link thimble::SmallBinaryField
 * SmallBinaryField\endlink object. Given two elements
 * <pre>
 *  uint32_t a , b;
 * </pre>
 * their sum can be computed via the \link thimble::SmallBinaryField::add()
 * add()\endlink member function of the \link thimble::SmallBinaryField
 * SmallBinaryField\endlink class
 * <pre>
 *  uint32_t c = gf.add(a,b);
 * </pre>
 * which is equivalent to an exclusive or operation, i.e.,
 * <pre>
 *  uint32_t c = a ^ b;
 * </pre>
 * A subtraction function
 * \link thimble::SmallBinaryField::sub() sub()\endlink is provided as well
 * <pre>
 *  uint32_t c = gf.sub(a,b);
 * </pre>
 * @attention In a binary field subtraction is equivalent to addition and
 * it is actually not necessary to distinguish between these two operations
 * in the context of a \link thimble::SmallBinaryField
 * SmallBinaryField\endlink; however, in some situations it can be useful
 * to distinguish between addition and subtraction for better readability
 * of program code.
 *
 * @subsection sec_finitefields_mul Multiplication
 *
 * In a binary field
 * <pre>
 *  SmallBinaryField gf(d);
 * </pre>
 * multiplication of two elements
 * <pre>
 *  uint32_t a , b;
 * </pre>
 * can be performed with the \link thimble::SmallBinaryField::mul()
 * mul()\endlink member function of a \link thimble::SmallBinaryField
 * SmallBinaryField\endlink object
 * <pre>
 *  uint32_t c = gf.mul(a,b);
 * </pre>
 *
 *
 * @subsection sec_finitefields_div Division
 *
 * In a \link thimble::SmallBinaryField SmallBinaryField\endlink
 * <pre>
 *  SmallBinaryField gf(d);
 * </pre>
 * the division between two elements
 * <pre>
 *  uint32_t a , b;
 * </pre>
 * where <code>b</code> must be non-zero can be formed with the
 * member function \link thimble::SmallBinaryField::div()
 * div()\endlink
 * <pre>
 *  uint32_t c = gf.div(a,b);
 * </pre>
 * such the expression <code>gf.mul(c,b)</code> equals <code>a</code>.
 *
 * @subsection sec_finitefields_inv Inversion
 *
 * To invert an element of a \link thimble::SmallBinaryField
 * SmallBinaryField\endlink
 * <pre>
 *  SmallBinaryField gf(d);
 * </pre>
 * the function \link thimble::SmallBinaryField::inv() inv()\endlink
 * is provided as a member of <code>gf</code>. For example, to invert
 * a non-zero
 * <pre>
 *  uint32_t a;
 * </pre>
 * we may run
 * <pre>
 *  uint32_t b = gf.inv(a);
 * </pre>
 * such that <code>gf.mul(a,b)==1</code>. Note, inversion is
 * equivalent to the expression
 * <pre>
 *  uint32_t b = gf.div(1,a);
 * </pre>
 *
 * @subsection sec_finitefields_neg Negation
 *
 * Even though negation of an element in a binary field does
 * not change the element, it can be useful for better readability. Therefore,
 * a negation function \link thimble::SmallBinaryField::neg() neg()\endlink
 * is provided by a \link thimble::SmallBinaryField SmallBinaryField\endlink
 * object.
 *
 * @author Benjamin Tams
 *
 * @see thimble::SmallBinaryField
 */

#ifndef THIMBLE_SMALLBINARYFIELD_H_
#define THIMBLE_SMALLBINARYFIELD_H_

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryPolynomial.h>

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *           Enables arithmetic of a small binary finite field which are
	 *           fields of power of two size.
	 *
	 * @details
	 *           A binary finite field is defined by an irreducible polynomial
	 *           \f$f(X)\f$ (over the prime field with two elements) that is of
	 *           degree \f$d\geq 1\f$. The corresponding field \f$GF(2^{d})\f$
	 *           is defined to be the smallest field containing the
	 *           field with two elements and in where \f$f(X)\f$
	 *           splits into linear factors. Thus, since finite fields are
	 *           known to be galois, if \f$\alpha\f$ is a root
	 *           of \f$f(X)\f$ then the elements of the field are of
	 *           the shape \f$c_0+c_1\alpha+...+c_{d-1}\alpha^{d-1}\f$ with
	 *           \f$c_i\in\{0,1\}\f$.
	 *           <br>
	 *           Elements of the finite fields that are represented by
	 *           instances of this class are given as <code>uint32_t</code>
	 *           between 0 and (excluding) the finite field's cardinality. If
	 *           an element of the finite field is given as
	 *           \f$c_0+c_1\alpha+...+c_{d-1}\alpha^{d-1}\f$ with
	 *           \f$c_i\in\{0,1\}\f$ then the integer
	 *           \f$c_0+c_1\cdot 2+...+c_{d-1}2^{d-1}\f$ represents the
	 *           element. Consequently, addition of two finite field elements
	 *           corresponds to an exclusive or operation of their
	 *           representing integers.
	 *           <br><br>
	 *           To create instances of this class we need to specify an
	 *           instance of the \link SmallBinaryPolynomial\endlink class,
	 *           which represent binary polynomials of small degree, first.
	 *           Furthermore, the represented polynomial must be irreducible
	 *           i.e. it must not have factors different from constants and
	 *           the polynomial itself. For example, the following code
	 *           generates a binary finite field with \f$2^{16}\f$ elements:
	 *           <pre>
	 * // Access an irreducible polynomial of degree 16
	 * SmallBinaryPolynomial f = SmallBinaryPolynomial::irreducible(16);
	 *
	 * SmallBinaryField field(f);
	 *           </pre>
	 *           The variable <code>field</code> contains exponent-,
	 *           logarithm- and inverse tables which accelerate field
	 *           arithmetic. Alternatively, the same finite field can
	 *           be created via
	 *           <pre>
	 * SmallBinaryField field(16);
	 *           </pre>
	 *           which wraps around the
	 *           \link SmallBinaryPolynomial.irreducible(int)\endlink
	 *           function.
	 *           <br>
	 *           Theoretically, this class can handle binary finite fields
	 *           of degree up to 31 but it is likely that there is not enough
	 *           memory available to hold the computations. Reasonable
	 *           field degrees may vary between 1 and 25.
	 *
	 * @warning
	 *           When instances of this class are created with a reducible
	 *           polynomial the constructor may run into an infinite loop.
	 *           For example, running the following code, does not
	 *           terminate:
	 *           <pre>
	 * // Create the polynomial X^2+1 which is the product
	 * // of X+1 and X+1
	 * SmallBinaryPolynomial f;
	 * f.setCoeff(0);
	 * f.setCoeff(2);
	 * // WARNING: The constructor runs into an infinite loop
	 * SmallBinaryField field(f);
	 *           </pre>
	 *           To prevent infinite loops one must assert that the alleged
	 *           irreducible polynomial used to define the finite field is in
	 *           fact irreducible. The polynomials that are returned by the
	 *           \link SmallBinaryPolynomial.irreducible(int)\endlink are
	 *           guaranteed to be irreducible.
	 */
	class THIMBLE_DLL SmallBinaryField {

	private:

		/**
		 * @brief The defining modulus polynomial of this finite field
		 */
		SmallBinaryPolynomial definingPolynomial;

		/**
		 * @brief The degree of this finite field over <code>GF(2)</code>
		 */
		int degree;

		/**
		 * @brief The cardinality of the finite field
		 */
		uint32_t size;

		/**
		 * @brief Represents a generator of the multiplicative group of unity
		 *        of the finite field.
		 */
		uint32_t generator;

		/**
		 * @brief     Exponent table
		 *
		 * @details
		 *            If <code>i</code> is an integer then
		 *            <code>expTable[i]</code> represents the field element
		 *            that is the <code>i</code>-th power of the generator
		 *            that is represented by \link generator\endlink.
		 */
		uint32_t *expTable;

		/**
		 * @brief     Logarithm table
		 *
		 * @details
		 *            If <code>a</code> is an integer that represents a finite
		 *            field element then <code>logTable[a]</code> is the
		 *            integer <code>i</code> such that the <code>i</code>-th
		 *            power of the \link generator\endlink is equals the field
		 *            element. Therefore, this table is called logarithm table
		 *            because they correspond to discrete logarithms.
		 */
		uint32_t *logTable;

		/**
		 * @brief     Table of multiplicative inverses
		 *
		 * @details
		 *            If <code>a</code> is an integer that represents a
		 *            nonzero finite field element then
		 *            <code>invTable[a]</code> represents its multiplicative
		 *            inverse.
		 */
		uint32_t *invTable;

		/**
		 * @brief    Frees all memory that is allocated by this finite field
		 *           to list the \link expTable\endlink,
		 *           \link logTable\endlink, and \link invTable\endlink
		 */
		void clear();

		/**
		 * @brief    Initializes the finite field with the defining polynomial
		 *           <code>f</code>
		 *
		 * @details
		 *           On initialization, the \link definingPolynomial\endlink
		 *           is set and the \link degree\endlink and the fields
		 *           \link size\endlink (which is <code>2^degree</code>)
		 *           are precomputed.
		 *           <br><br>
		 *           Next, memory is allocated for \link expTable\endlink,
		 *           \link logTable\endlink, and \link invTable\endlink such
		 *           that each of them can store \link size\endlink 32-bit
		 *           integers.
		 *           <br><br>
		 *           A candidate for the fields multiplicative group's
		 *           generator is generated at random and the exponent
		 *           table (\link expTable\endlink) and the logarithm table
		 *           (\link logTable\endlink) entries are filled. If the
		 *           generator turns out to only generate a proper subgroup,
		 *           then this step is repeated until a proper generator
		 *           is found.
		 *           <br><br>
		 *           Finally, the inverse table (\link invTable\endlink)
		 *           is filled.
		 *
		 * @param f
		 *           The defining polynomial of the finite field which
		 *           must be irreducible and of degree between 1 and 31.
		 *
		 * @warning
		 *           The function may run into an infinite loop if
		 *           <code>f</code> is a reducible polynomial. It is
		 *           not checked whether this is the case.
		 *           <br><br>
		 *           If not sufficient memory can be allocated for
		 *           the tables or if this polynomial is of degree smaller
		 *           than 1 or greater than 31 then this functions prints
		 *           an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 */
		void init( const SmallBinaryPolynomial & f );

	public:

		/**
		 * @brief
		 *           Initializes a finite field defined by the specified
		 *           polynomial
		 *
		 * @param f
		 *           The defining polynomial of the finite field which
		 *           must be irreducible and of degree between 1 and 31.
		 *
		 * @warning
		 *           The constructor may run into an infinite loop if
		 *           <code>f</code> is a reducible polynomial. It is
		 *           not checked whether this is the case.
		 *           <br><br>
		 *           If not sufficient memory can be allocated for
		 *           the tables or if this polynomial is of degree smaller
		 *           than 1 or greater than 31 then this functions prints
		 *           an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @see SmallBinaryPolynomial.irreducible(int)
		 */
		inline SmallBinaryField( const SmallBinaryPolynomial & f = 2 ) {
			init(f);
		}

		/**
		 * @brief
		 *           Constructs a binary finite field of specified degree
		 *
		 * @details
		 *           The constructor wraps around the
		 *           \link SmallBinaryPolynomial::irreducible(int)\endlink
		 *           function which is guaranteed to return an irreducible
		 *           polynomial if <code>degree</code> ranges between 1 and
		 *           31. Thus, the constructor, unlike
		 *           \link SmallBinaryField(const SmallBinaryPolynomial&)\endlink
		 *           can not run into an infinite loop. Furthermore, it
		 *           is guaranteed that repeats calls of the constructor
		 *           with the same degree constructs the same finite
		 *           field (i.e. the same isomorphy class).
		 *
		 * @param degree
		 *           The degree of the finite field
		 *
		 * @warning
		 *           If <code>degree</code> does not range between 1 and 31
		 *           the constructor prints an error message to
		 *           <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		inline SmallBinaryField( int degree ) {
			init(SmallBinaryPolynomial::irreducible(degree));
		}

		/**
		 * @brief
		 *           Copy constructor.
		 *
		 * @param gf
		 *           The binary finite field of which a copy is constructed.
		 */
		SmallBinaryField( const SmallBinaryField & gf );

		/**
		 * @brief Destructor that completely frees the memory helt by the
		 *        finite field.
		 */
		inline ~SmallBinaryField() {
			clear();
		}

		/**
		 * @brief
		 *           Assignment operator.
		 *
		 * @param gf
		 *           The binary finite field of which a copy is assigned
		 *           to this instance.
		 */
		SmallBinaryField &operator=( const SmallBinaryField & gf );

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
		static void swap( SmallBinaryField & gf1 , SmallBinaryField & gf2 );

		/**
		 * @brief
		 *            Initializes this finite field by the specified
		 *            polynomial.
		 *
		 * @details
		 *            The data for the old finite field will be dismissed.
		 *
		 * @param f
		 *            The defining polynomial of the finite field which
		 *            must be irreducible and of degree between 1 and 31.
		 *
		 * @warning
		 *           The method may run into an infinite loop if
		 *           <code>f</code> is a reducible polynomial. It is
		 *           not checked whether this is the case.
		 *           <br><br>
		 *           If not sufficient memory can be allocated for
		 *           the tables or if this polynomial is of degree smaller
		 *           than 1 or greater than 31 then this functions prints
		 *           an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @see SmallBinaryPolynomial.irreducible(int)
		 */
		inline void initialize( const SmallBinaryPolynomial & f ) {
			clear();
			init(f);
		}

		/**
		 * @brief Returns the defining polynomial of the finite field
		 *
		 * @return the defining polynomial of the finite field
		 */
		inline SmallBinaryPolynomial getDefiningPolynomial() const {
			return this->definingPolynomial;
		}

		/**
		 * @brief Returns the degree of the field over <code>GF(2)</code>
		 *
		 * @return the degree of the finite field over <code>GF(2)</code>
		 */
		inline int getDegree() const {
			return this->degree;
		}

		/**
		 * @brief Returns the number of elements contained in the finite field.
		 *
		 * @return the size of the finite field
		 */
		inline uint32_t getCardinality() const {
			return this->size;
		}

		/**
		 * @brief Returns a generator of the finite field's multiplicative
		 *        group of unity.
		 *
		 * @return a generator of the finite field
		 */
		inline uint32_t getGenerator() const {
			return this->generator;
		}

		/**
		 * @brief Generates a random element of the finite field
		 *
		 * @param tryRandom
		 *        If <code>tryRandom</code> is <code>true</code> then
		 *        the function attempts to use a cryptographic random
		 *        generator; otherwise, this function wraps around
		 *        the standard <code>rand()</code> function.
		 *
		 */
		inline uint32_t random( bool tryRandom = false ) const {

			return MathTools::rand32(tryRandom)%(this->size);
		}

		/**
		 * @brief
		 *            Computes the sum of two finite field elements
		 *
		 * @details
		 *            The sum of two binary finite field elements in bit
		 *            representation is equivalent to a bitwise
		 *            exclusive or operation.
		 *
		 * @param a
		 *            summand
		 *
		 * @param b
		 *            summand
		 *
		 * @return
		 *            The sum of <code>a</code> and <code>b</code>.
		 */
		inline uint32_t add( uint32_t a , uint32_t b ) const {
			return a^b;
		}

		/**
		 * @brief
		 *            Computes the difference of two finite field elements
		 *
		 * @details
		 *            The difference of two binary finite field elements in
		 *            bit representation is equivalent to a bitwise
		 *            exclusive or operation. In fact, there is no difference
		 *            to addition but for consistency we provide a subtraction
		 *            function.
		 *
		 * @param a
		 *            minuend
		 *
		 * @param b
		 *            subtrahend
		 *
		 * @return
		 *            The difference of <code>a</code> and <code>b</code>.
		 */
		inline uint32_t sub( uint32_t a , uint32_t b ) const {
			return a^b;
		}

		/**
		 * @brief
		 *            Computes the additive negative of a finite field element.
		 *
		 * @details
		 *            For binary fields the negative of an element is just that
		 *            element. But for consistency, we provide a negation
		 *            function.
		 *
		 * @param a
		 *            element of the finite field in integer representation
		 */
		inline uint32_t neg( uint32_t a ) const {
			return a;
		}

		/**
		 * @brief
		 *            Computes the product of two element of the finite field
		 *
		 * @details
		 *            The multiplication of two finite field elements canonically
		 *            corresponds to the multiplication of two polynomials
		 *            (carry-less multiplication) and subsequent reduction
		 *            modulo the defining polynomial. This can become quite time-
		 *            consuming if multiplication is frequently required.
		 *            <br>
		 *            Therefore, this class implements multiplication by accessing
		 *            exponent and logarithm tables which only requires two
		 *            three accesses of arrays. These have been precomputed on
		 *            the initialization of the finite field which is only possible
		 *            for small finite fields (as represented by this class) for
		 *            the tables have to be hold in the memory.
		 *
		 * @param a
		 *            first factor
		 *
		 * @param b
		 *            second factor
		 *
		 * @return
		 *            product of <code>a</code> and <code>b</code>
		 *
		 * @warning
		 *            The function may run into unexpected behaviour if
		 *            <code>a</code> or <code>b</code> are integers larger
		 *            than or equals the finite field's size.
		 */
		inline uint32_t mul( uint32_t a , uint32_t b ) const {

			if ( a==0 || b==0 ) {
				return 0;
			}

			uint64_t i = (uint64_t)(this->logTable[a])+
					     (uint64_t)(this->logTable[b]);

			if ( i >= this->size-1 ) {
				i -= this->size-1;
			}

			return this->expTable[i];
		}

		/**
		 * @brief
		 *            Computes the multiplicative inverse of an element of the
		 *            finite field
		 *
		 * @details
		 *            The computation of the inverse of an element in a finite
		 *            field corresponds to performing an extended Euclidean
		 *            algorithm and subsequent reduction modulo the defining
		 *            polynomial.
		 *            <br>
		 *            For instances of this class the inverses have been
		 *            precomputed and are stored in a table that simply is
		 *            accessed. This makes the effort for computing an inverse
		 *            negligible.
		 *
		 * @param a
		 *            Element of the finite field
		 *
		 * @return
		 *            Multiplicative inverse of <code>a</code> in the
		 *            finite field
		 *
		 * @warning
		 *            The function may run into unexpected behaviour if
		 *            <code>a</code> is integer larger than or equals the
		 *            finite field's size. In case <code>a=0</code> the result
		 *            will be <code>0</code> even though the inverse of 0
		 *            is not defined.
		 */
		inline uint32_t inv( uint32_t a ) const {
			return this->invTable[a];
		}

		/**
		 * @brief
		 *            Computes the quotient of <code>a</code> divided by
		 *            <code>b</code>, i.e. \f$ a\cdot b^{-1}\f$.
		 *
		 * @details
		 *            The function wraps around
		 *            \link inv(uint32_t)const\endlink and
		 *            \link mul(uint32_t,uint32_t)const\endlink.
		 *
		 * @param a
		 *            numerator element
		 *
		 * @param b
		 *            denominator element
		 *
		 * @return
		 *            \f$ a\cdot b^{-1}\f$
		 *
		 * @warning
		 *            The function may run into unexpected behaviour if
		 *            <code>a</code> or <code>b</code> are integers larger
		 *            than or equals the finite field's size.
		 *            If <code>b=0</code> the result will be <code>0</code>
		 *            even though the inverse is of 0 is defined.
		 */
		inline uint32_t div( uint32_t a , uint32_t b ) const {
			return mul(a,inv(b));
		}

		/**
		 * @brief
		 *            Returns a constant reference to a field
		 *            with two elements.
         *
		 * @return
		 *            Constant reference to a field with two elements.
		 */
		static const SmallBinaryField & binary();
	};
}


#endif /* THIMBLE_SMALLBINARYFIELD_H_ */
