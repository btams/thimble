/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file BinaryPolynomial.cpp
 *
 * @brief
 *            Implements a mechanism for representing and computing
 *            with polynomials of arbitrary degree as provided
 *            by the 'BinaryPolynomial.h' header.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>


#include <thimble/math/MathTools.h>
#include <thimble/math/numbertheory/BinaryPolynomial.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

    /**
     * @brief
     *           Standard constructor.
     *
     * @details
     *           If in the degree of the polynomials that this class is used
     *           to represent is known in advance, it can be specified
     *           via the parameter <i>d</i>.
     *
     * @param d
     *           Ensures that the instance can represent polynomials up to
     *           degree <i>d</i> without causing reallocation.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     */
    BinaryPolynomial::BinaryPolynomial( int d ) {

        this->numWords = 1;
        if ( d >= 0 && (d+1)/32+(d%32?1:0) >= 1 ) {
            this->numWords = (d+1)/32+((d+1)%32?1:0);
        }

        this->data = (uint32_t*)malloc( this->numWords * sizeof(uint32_t) );
        if ( this->data == NULL ) {
            cerr << "BinaryPolynomial: out of memory." << endl;
            exit(EXIT_FAILURE);
        }

        memset(this->data,0,this->numWords*sizeof(uint32_t));

        this->degree = -1;
    }

    /**
     * @brief
     *           Copy constructor.
     *
     * @param f
     *           The binary polynomial of which a copy is constructed.
     *
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial::BinaryPolynomial( const BinaryPolynomial & f ) {

        this->numWords = (f.degree+1)/32+((f.degree+1)%32?1:0);
        if ( this->numWords <= 0 ) {
            this->numWords = 1;
        }

        this->data = (uint32_t*)malloc( this->numWords * sizeof(uint32_t) );
        if ( this->data == NULL ) {
            cerr << "BinaryPolynomial: out of memory." << endl;
            exit(EXIT_FAILURE);
        }
        memcpy(this->data,f.data,this->numWords*sizeof(uint32_t));

        this->degree = f.degree;
    }

    /**
     * @brief
     *           Destructor.
     *
     * @details
     *           Frees the memory used by this instances to store
     *           the coefficients of the binary polynomial.
     */
    BinaryPolynomial::~BinaryPolynomial() {
        free(this->data);
    }

    /**
     * @brief
     *           Assignment operator (procedural version).
     *
     * @param f
     *           The binary polynomial of which this instances is assigned
     *           a copy.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::assign( const BinaryPolynomial & f ) {

        // Nothing to do on self-assignment.
        if ( this != &f ) {

            // Make the polynomial zero
            setZero();

            if ( !f.isZero() ) { // Nothing to do if 'f' is zero

                int n = (f.degree+1)/32+((f.degree+1)%32?1:0);

                // Ensure that this instance can hold a polynomial
                // of the degree of 'f'
                ensureDegree(f.degree);

                // Copy the coefficients.
                memcpy(this->data,f.data,n*sizeof(uint32_t));

                // Copy the degree.
                this->degree = f.degree;
            }
        }

    }


    /**
     * @brief
     *           Returns the Hamming weight of this polynomial.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^d c_j\cdot X^j
     *           \f]
     *           where \f$c_j\in\{0,1\}\f$. The function returns
     *           the number of \f$j\f$ where \f$c_j\neq 0\f$. More
     *           specifically, it returns the cardinality
     *           \f[
     *            \#\{j=0,...,d~|~c_j\neq 0\}.
     *           \f]
     *
     * @return
     *           The Hamming weight of this polynomial.
     */
    int BinaryPolynomial::hammingWeight() const {


        return MathTools::hw32
                (this->data,(this->degree+1)/32+((this->degree+1)%32?1:0));
    }

    /**
     * @brief
     *           Clears all data held by this polynomial
     *           such that it becomes identically to 0.
     */
    void BinaryPolynomial::setZero() {

        // The first 'm' elements in 'data' can potentially be non-zero.
        int m = (this->degree+1)/32 + ((this->degree+1)%32?1:0);

        // Overwrite the elements by zeroes.
        memset(this->data,0,m*sizeof(uint32_t));

        // A negative degree indicates a zero polynomial.
        this->degree = -1;
    }

    /**
     * @brief
     *           Access the <i>i</i>th coefficients of the binary
     *           polynomial.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^dc_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The function returns
     *           the coefficient \f$c_i\f$.
     *
     * @param i
     *           The index of the coefficient that is accessed.
     *
     * @return
     *           <code>true</code> if \f$c_i=1\f$; otherwise, if
     *           \f$c_i=0\f$, the function returns <code>false</code>.
     *
     * @warning
     *           If <i>i</i> is less than 0, an error message is
     *           printed to <code>stderr</code> and the program
     *           exits with status 'EXIT_FAILURE'.
     */
    bool BinaryPolynomial::getCoeff( int i ) const {

        if ( i < 0 ) {
            cerr << "BinaryPolynomial::getCoeff: "
                 << "index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        if ( i > this->degree ) {
            return false;
        }

        // Get the correct 32 bit word
        uint32_t c = this->data[i>>5];

        // Which bit in the 32 bit word encodes the coefficient?
        int s = i%32;

        // Return the value of the bit in 32 bit word.
        return c&(1<<s)?true:false;
    }

    /**
     * @brief
     *           Ensures that the <i>i</i>th coefficient of the
     *           polynomial is set to 1.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^dc_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The method ensures
     *           that \f$c_i=1\f$.
     *
     * @param i
     *           The index of the coefficient to ensuring to be
     *           equal 1.
     *
     * @warning
     *           If <i>i</i> is less than 0, an error message is
     *           printed to <code>stderr</code> and the program
     *           exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::setCoeff( int i ) {

        if ( i < 0 ) {
            cerr << "BinaryPolynomial::getCoeff: "
                 << "index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        if ( i > this->degree ) {
            ensureDegree(i);
            this->degree = i;
        }

        int j , s;
        j = i/32;
        s = i%32;
        this->data[j] |= (((uint32_t)1)<<s);
    }

    /**
     * @brief
     *           Ensures that the <i>i</i>th coefficient of the
     *           polynomial is set to 0.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^dc_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The method ensures
     *           that \f$c_i=0\f$.
     *
     * @param i
     *           The index of the coefficient to ensuring to be
     *           equal 1.
     *
     * @waupdateDegreerning
     *           If <i>i</i> is less than 0, an error message is
     *           printed to <code>stderr</code> and the program
     *           exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::clearCoeff( int i ) {

        if ( i < 0 ) {
            cerr << "BinaryPolynomial::getCoeff: "
                 << "index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        if ( i > this->degree ) {
            // Nothing to do because the coefficient is alread zero.
            return;
        }

        int j , s;
        j = i/32;
        s = i%32;
        this->data[j] &= ~(((uint32_t)1)<<s);

        if ( i == this->degree ) {
            // It is only possible that the degree changes if
            // the leading coefficient is cleared.
            degreeCouldBeSmaller();
        }
    }


    /**
     * @brief
     *           Ensures that this polynomial can store at least
     *           <i>d+1/i> coefficients without reallocation.
     *
     * @param d
     *           A degree bound for a polynomial that can
     *           be represented by this instance without
     *           reallocation.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::ensureDegree( int d ) {

        if ( d >= 0 ) {

            int numWords = (d+1)/32+((d+1)%32?1:0);

            if ( numWords > this->numWords ) {

                this->data = (uint32_t*)realloc
                       ( this->data , numWords * sizeof(uint32_t) );
                if ( this->data == NULL ) {
                    cerr << "BinaryPolynomial::ensureCapacity: "
                         << "out of memory." << endl;
                    exit(EXIT_FAILURE);
                }

                // Sets the new coefficients to zero.
                memset
                (this->data+this->numWords,0,
                 (numWords-this->numWords)*sizeof(uint32_t));

                this->numWords = numWords;
            }
        }
    }

    /**
     * @brief
     *           Initializes this instance by a polynmial with the
     *           specified number of random bit coefficients.
     *
     * @details
     *           The degree of the resulting polynomial will be
     *           strictly smaller than <code>numCoefficients</code>.
     *
     * @param numCoefficients
     *           The number of random coefficients.
     *
     * @param tryRandom
     *           If <code>true</code>, the method uses a cryptographic
     *           number generator if available on the system; otherwise,
     *           the method wraps around the standard <code>rand()</code>
     *           function.
     *
     * @warning
     *           If <code>numCoefficients</code> is negative, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::random
        ( int numCoefficients  , bool tryRandom ) {

        if ( numCoefficients < 0 ) {
            cerr << "BinaryPolynomial::rand: "
                 << "number of random coefficients must be greater "
                 << "than or equals 0." << endl;
            exit(EXIT_FAILURE);
        }

        ensureDegree(numCoefficients-1);
        setZero();

        int m0 , m1;
        m0 = numCoefficients/32;
        m1 = numCoefficients%32;

        for ( int j = 0 ; j < m0 ; j++ ) {
            this->data[j] = MathTools::rand32(tryRandom);
        }


        if ( m1 ) {
            this->data[m0] = MathTools::rand32(tryRandom);
            this->data[m0] &= (((uint32_t)1)<<m1)-1;
        }

        this->degree = numCoefficients - 1;
        degreeCouldBeSmaller();
    }

    /**
     * @brief
     *           Replaces this instance by polynomial with the
     *           specified number of random coefficients such that
     *           it has the specified Hamming weight.
     *
     * @param numCoefficients
     *           The number of random coefficients.
     *
     * @param hammingWeight
     *           The Hamming weight of the random polynomial.
     *
     * @param tryRandom
     *           If <code>true</code>, the method uses a cryptographic
     *           number generator if available on the system; otherwise,
     *           the method wraps around the standard <code>rand()</code>
     *           function.
     *
     * @warning
     *           If <code>hammingWeight</code> is negative or greater than
     *           <code>numCoefficients</code>, an error message is printed
     *           to <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If <code>numCoefficients</code> is negative, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::wrandom
    ( int numCoefficients , int hammingWeight , bool tryRandom ) {

        if ( numCoefficients < 0 ) {
            cerr << "BinaryPolynomial::rand2: "
                 << "number of random coefficients must be greater "
                 << "than or equals 0." << endl;
            exit(EXIT_FAILURE);
        }

        if ( hammingWeight < 0 ) {
            cerr << "BinaryPolynomial::rand2: "
                 << "hamming weight must be non-negative." << endl;
            exit(EXIT_FAILURE);
        }

        if  ( hammingWeight > numCoefficients ) {
            cerr << "BinaryPolynomial::rand2: "
                 << "hamming weight must be bounded by the "
                 << "number of coefficients." << endl;
            exit(EXIT_FAILURE);
        }

        setZero();

        if ( numCoefficients > 0 ) {

            int *indices = (int*)malloc
                    ( numCoefficients * sizeof(int));
            if ( indices == NULL ) {
                cerr << "BinaryPolynomial::rand2: out of memory." << endl;
                exit(EXIT_FAILURE);
            }

            int hw = numCoefficients;
            for ( int j = 0 ; j < numCoefficients ; j++ ) {
                indices[j] = j;
            }

            while ( hw > hammingWeight ) {
                // which indices to remove?
                int i = (int)(MathTools::rand64(tryRandom) % (uint64_t)hw);
                for ( int j = i ; j+1 < hw ; j++ ) {
                    indices[j] = indices[j+1];
                }
                --hw;
            }

            // The first 'hw' positions in 'indices' of the polynomial
            // are set to '1'.
            for ( int j = hw-1 ; j >= 0 ; j-- ) {
                setCoeff(indices[j]);
            }

            free(indices);
        }
    }

    /**
     * @brief
     *           Replaces this polynomial by its multipliciation
     *           with \f$X^e\f$.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^dc_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The method replaces the
     *           polynomial by
     *           \f[
     *            X^e\cdot \sum_{j=0}^dc_j\cdot X^j.
     *           \f]
     *
     * @param e
     *           The exponent of the shift.
     *
     * @warning
     *           If <i>e</i> is negative, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::leftShift( int e ) {

        if ( e < 0 ) {
            cerr << "BinaryPolynomial::leftShift: negative exponent."
                 << endl;
            exit(EXIT_FAILURE);
        }

        // Nothing to do if exponent or polynomial is zero
        if ( e > 0 && !isZero() ) {

            ensureDegree(this->degree+e+32);

            // Split the exponent into the number of 32 bit words 'e1'
            // needed to be shifted and into the number of shifts 'e0'
            // that remain.
            int e1 , e0;
            e1 = e/32;
            e0 = e%32;

            // Multiple of 32 shifts
            if ( e1 != 0 ) {

                int size = (this->degree+1)/32+((this->degree+1)%32?1:0);

                for ( int j = 0 ; j < size ; j++ ) {
                    this->data[size+e1-j-1] = this->data[size-j-1];
                }

                memset(this->data,0,e1*sizeof(uint32_t));

                this->degree += e1*32;
            }

            // Shift by the remaining bits
            if ( e0 > 0 ) {

                int size = (this->degree+1)/32+((this->degree+1)%32?1:0);

                uint64_t tmp;
                for ( int j = size-1 ; j >= 0 ; j-- ) {
                    tmp = this->data[j];
                    tmp <<= e0;
                    this->data[j] = tmp & (uint32_t)0xFFFFFFFF;
                    tmp >>= 32;
                    this->data[j+1] ^= (uint32_t)tmp;
                }

                this->degree += e0;
            }
        }
    }

    /**
     * @brief
     *           Essentially, polynomial division by the
     *           monomial \f$X^e\f$.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^dc_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The method replaces the
     *           polynomial by
     *           \f[
     *            \sum_{j=e}^dc_j\cdot X^{j-e}.
     *           \f]
     *
     * @param e
     *           The exponent of the shift.
     *
     * @warning
     *           If <i>e</i> is negative, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::rightShift( int e ) {

        if ( e < 0 ) {
            cerr << "BinaryPolynomial::rightShift: negative exponent."
                 << endl;
            exit(EXIT_FAILURE);
        }

        if ( e > deg() ) {
            setZero();
            return;
        }

        // Nothing to do if exponent or polynomial is zero
        if ( e > 0 && !isZero() ) {

            // Split the exponent into the number of 32 bit words 'e1'
            // needed to be shifted and into the number of shifts 'e0'
            // that remain.
            int e1 , e0;
            e1 = e/32;
            e0 = e%32;

            // Shift multiples of 32
            if ( e1 != 0 ) {

                int size = (this->degree+1)/32+((this->degree+1)%32?1:0);

                for ( int j = 0 ; j < size-e1 ; j++ ) {
                    this->data[j] = this->data[j+e1];
                    this->data[j+e1] = 0;
                }

                this->degree -= e1*32;
            }

            // Shift remaining coefficients/bits
            if ( e0 != 0 && !isZero() ) {

                int size = (this->degree+1)/32+((this->degree+1)%32?1:0);

                uint32_t mask = (e0==31?0xFFFFFFFFu:((1<<(e0+1))-1)) , back;

                this->data[0] >>= e0;

                for ( int j = 1 ; j < size ; j++ ) {
                    back = this->data[j] & mask;
                    back <<= 32-e0;
                    this->data[j-1] ^= back;
                    this->data[j] >>= e0;
                }
            }
        }
    }

    /**
     * @brief
     *           Replaces the polynomial by its derivative.
     *
     * @details
     *           Write the polynomial as
     *           \f[
     *            \sum_{j=0}^d c_j\cdot X^j
     *           \f]
     *           with \f$c_j\in\{0,1\}\f$. The method replaces
     *           the polynomial by the formal derivative
     *           \f[
     *            \sum_{j=1,3,5,...,d} c_j\cdot X^{j-1}.
     *           \f]
     */
    void BinaryPolynomial::derive() {

        static uint32_t evenBits = 2863311530u;

        int m = (this->degree+1)/32+((this->degree+1)%32?1:0);

        for ( int j = 0 ; j < m ; j++ ) {
            this->data[j] &= evenBits;
        }

        rightShift(1);

        degreeCouldBeSmaller();
    }


    /**
     * @brief
     *           Computes the <i>d</i>-reverse of this polynomial.
     *
     * @param d
     *           The <i>d</i> reverse of a polynomial <i>f</i> of degree
     *           \f$n\leq d\f$ is defined as the polynomial
     *           \f[
     *            X^d\cdot \sum_{j=0}^n f_j\cdot X^{d-n+j}.
     *           \f]
     *
     * @return
     *           The <i>d</i>-reverse of this polynomial.
     *
     * @warning
     *           If <i>d</i> is negative or smaller than the
     *           polynomial's degree (i.e., \link deg()\endlink), an
     *           error message is printed to <code>stderr</code> and
     *           the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::reverse( int d ) const {

        if ( d < 0 ) {
            cerr << "BinaryPolynomial::reverse: "
                 << "input must be positive." << endl;
            exit(EXIT_FAILURE);
        }

        int n = deg();

        if ( d < n ) {
            cerr << "BinaryPolynomial::reverse: "
                 << "input must be larger than the polynomial's degree."
                 << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial h(d);

        for ( int j = 0 ; j <= n ; j++ ) {
            h.setCoeff(d-j,getCoeff(j));
        }

        return h;
    }

    /**
     * @brief
     *            Determines the <i>n</i>th cyclotomic polynomial.
     *
     * @details
     *            The <i>n</i>th cyclotomic polynomial is defined as
     *            the following:
     *            \f[
     *   \Phi_n(X):=\prod_{1\leq k<n,\\\gcd(k,n)=1}(X-e^{2\pi ik/n}).
     *            \f]
     *            This function computes this polynomial whose
     *            coefficients are reduced modulo 2.
     *
     *            The implementation of the cyclotomic polynomial
     *            computation follows the description of
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @param n
     *            Index of the cyclotomic polynomial.
     *
     * @return
     *            The <i>n</i>th cyclotomic polynomial.
     *
     * @warning
     *            If <i>n</i> is zero or negative, an error message is
     *            printed to <code>stderr</code> and the program exits
     *            with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::cyclo( int n ) {

        if ( n <= 0 ) {
            cerr << "BinaryPolynomial::binaryCyclotomicPolynomial: "
                 << "input must be positive." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial phi(n) , f(n) , xp(n) , fxp(n) , tmp0(n) , tmp1(n);

        // Step 1 (Algorithm 14.48 in [vzGth])
        f.setCoeff(1);
        f.setCoeff(0);

        // Step 2 (Algorithm 14.48 in [vzGth])
        int m = n , n0 = 1;
        for ( int p = 2 ; p <= m ; p++ ) {

            if ( m % p == 0 ) {

                n0 *= p;

                do {
                    m /= p;
                } while ( m % p == 0 );

                xp.setZero();
                xp.setCoeff(p);
                eval(fxp,f,xp);
                divRem(tmp0,tmp1,fxp,f);
                swap(f,tmp0);
            }
        }

        // Step 3 (Algorithm 14.48 in [vzGth])
        xp.setZero();
        xp.setCoeff(n/n0);
        eval(phi,f,xp);

        return phi;
    }

    /**
     * @brief
     *           Tests two polynomials on equality.
     *
     * @param f
     *           First input polynomial.
     *
     * @param g
     *           Second input polynomial.
     *
     * @return
     *           <code>true</code> if <i>f</i> and <i>g</i> represent
     *           the same polynomials and, otherwise, <code>false<code>.
     */
    bool BinaryPolynomial::areEqual
    ( const BinaryPolynomial & f , const BinaryPolynomial & g ) {

        if ( f.degree != g.degree ) {
            return false;
        }

        int m = (f.degree+1)/32+((f.degree+1)%32?1:0);

        for ( int j = 0 ; j < m ; j++ ) {
            if ( f.data[j] != g.data[j] ) {
                return false;
            }
        }

        return true;
    }


    /**
     * @brief
     *            Low-level function for computing the Hamming distance
     *            between two polynomials.
     *
     * @details
     *            The function assumes that <i>A</i> and <i>B</i> contain
     *            <i>m</i> and <i>n</i> integers of width 32 bit,
     *            respectively, where \f$0\leq m\leq n\$f. If not fulfilled,
     *            the function runs into undocumented behaviour.
     *
     * @param A
     *            Array containing <i>m</i> integers of width 32 bits.
     *
     * @param m
     *            The number of 32 bit integers in <i>A</i>.

     * @param B
     *            Array containing <i>n</i> integers of width 32 bits.
     *
     * @param n
     *            The number of 32 bit integers in <i>B</i>.
     *
     * @return
     *            The Hamming distance between <i>(A,m)</i> and <i>(B,n)</i>.
     */
    static int lowLevelHammingDistance
    ( register const uint32_t *A , register int m ,
      register const uint32_t *B , register int n ) {

        register int j;
        register int w = 0;

        for ( j = 0 ; j < m ; j++ , A++ , B++ ) {
            w += MathTools::hammingWeight((*A)^(*B));
        }

        for ( ; j < n ; j++ , B++ ) {
            w += MathTools::hammingWeight(*B);
        }

        return w;
    }

    /**
     * @brief
     *           Determines the number of coefficients in which two
     *           polynomials differ.
     *
     * @details
     *           Write
     *           \f[
     *            f(X)=\sum_{j=0}^n f_j X^j
     *           \f]
     *           and
     *           \f[
     *            g(X)=\sum_{j=0}^n g_j X^j
     *           \f]
     *           where \f$f_j,g_j\in\{0,1\}\f$. The function returns
     *           the number of <i>j</i> in which the \f$f_j\f$ and
     *           \f$g_j\f$ differ; more specifically, it returns
     *           \f[
     *            \#\{j=0,...,n~|~f_j\neq g_j\}.
     *           \f]
     *
     * @param f
     *           First input polynomial.
     *
     * @param g
     *           Second input polynomial.
     *
     * @return
     *           The hamming distance between <i>f</i> and <i>g</i>.
     */
    int BinaryPolynomial::hammingDistance
    ( const BinaryPolynomial & f , const BinaryPolynomial & g ) {

        int m , n;

        m = (f.degree+1)/32+((f.degree+1)%32?1:0);
        n = (g.degree+1)/32+((g.degree+1)%32?1:0);

        int w;
        if ( m < n ) {
            w = lowLevelHammingDistance(f.data,m,g.data,n);
        } else {
            w = lowLevelHammingDistance(g.data,n,f.data,m);
        }

        return w;
    }

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
    void BinaryPolynomial::swap( BinaryPolynomial & f , BinaryPolynomial & g ) {

        uint32_t *data;
        int numWords , degree;

        data = f.data;
        numWords = f.numWords;
        degree = f.degree;

        f.data = g.data;
        f.numWords = g.numWords;
        f.degree = g.degree;

        g.data = data;
        g.numWords = numWords;
        g.degree = degree;
    }

    /**
     * @brief
     *           Computes the sum of two binary polynomials.
     *
     * @details
     *           The method replaces the polynomial <i>f</i>
     *           by the sum of the two polynomials <i>g</i>
     *           and <i>h</i>.
     *
     * @param f
     *           Sum of <i>g</i> and <i>h</i>.
     *
     * @param g
     *           First summand.
     *
     * @param h
     *           Second summand.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::add
    ( BinaryPolynomial & f ,
      const BinaryPolynomial & g ,
      const BinaryPolynomial & h ) {

        if ( &g == &h ) {
            f.setZero();
            return;
        }

        int d = max(g.degree,h.degree);

        f.ensureDegree(d+1);


        MathTools::mxor32(f.data,g.data,h.data,(d+1)/32+((d+1)%32?1:0));

        f.degree = d+1;

        f.degreeCouldBeSmaller();
    }

    /**
     * @brief
     *           Low-level method for adding the product of two
     *           polynomials, where one is degree not larger than 31,
     *           to another polynomial.
     *
     * @param C
     *           The coefficients of the polynomial to which the
     *           product is added; should contain at least
     *           <i>n+1</i> valid unsigned 32 bit integers.
     *
     * @param A
     *           The coefficients of the first factor polynomial;
     *           should contain at least <i>n</i> valid unsinged 32 bit
     *           integers.
     *
     * @param n
     *           The number of 32 bit words contained in <code>A</code>.
     *
     * @param b
     *           The coefficients of the second factor polynomial being
     *           of degree at most 31.
     */
    static void addMul
    ( register uint32_t *C , register const uint32_t *A , int n , uint32_t b ) {

        register int i;
        register uint32_t carry = 0;
        register uint64_t c;

        for ( i = 0 ; i < n ; i++ , A++ , C++ ) {
            c = MathTools::clmul(*A,b);
            *C ^= (c&((uint32_t)0xFFFFFFFF))^carry;
            carry = c >> 32;
        }

        *C ^= carry;
    }

    /**
     * @brief
     *            Low-level method for computing the product
     *            of two binary polynomials.
     *
     * @details
     *            Stores the product of the polynomial
     *            specified by <i>(A,m)</i> and <i>(B,n)</i>
     *            in the array <i>C</i>.
     *
     * @param C
     *            The array in which the coefficients of the product
     *            is stored; should be able to hold at least
     *            <i>m+n</i> unsigned 32 bit integers.
     *
     * @param A
     *            Coefficients of the first factor polynomial.
     *
     * @param m
     *            Number of unsigned 32 bit words stored in <i>A</i>.
     *
     * @param B
     *            Coefficients of the second factor polynomial.
     *
     * @param n
     *            Number of unsigned 32 bit words stored in <i>B</i>.
     *
     */
    static void plainMul
    ( register uint32_t *C ,
      register const uint32_t *A , register int m ,
      register const uint32_t *B , register int n ) {

        register int j;

        memset(C,0,(m+n)*sizeof(uint32_t));

        for ( j = 0 ; j < m ; j++ , C++ , A++ ) {
            addMul(C,B,n,*A);
        }
    }

    /**
     * @brief
     *           Computes the produce of two binary polynomials.
     *
     * @details
     *           The method replaces the polynomial <i>f</i>
     *           by the product of the two polynomials <i>g</i>
     *           and <i>h</i>.
     *
     * @param f
     *           Product of <i>g</i> and <i>h</i>.
     *
     * @param g
     *           First factor.
     *
     * @param h
     *           Second factor.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::mul
    ( BinaryPolynomial & f ,
      const BinaryPolynomial & g ,
      const BinaryPolynomial & h ) {

        if ( &f == &g ) {
            BinaryPolynomial tg(g);
            mul(f,tg,h);
            return;
        }

        if ( &f == &h ) {
            BinaryPolynomial th(h);
            mul(f,g,th);
            return;
        }

        f.setZero();

        if ( !g.isZero() && !h.isZero() ) {

            int m , n;
            m = (g.degree+1)/32+((g.degree+1)%32?1:0);
            n = (h.degree+1)/32+((h.degree+1)%32?1:0);

            f.ensureDegree(32*(m+n+1));

            plainMul(f.data,g.data,m,h.data,n);

            f.degree = g.degree+h.degree;
        }
    }

    /**
     * @brief
     *           Computes a Euclidean division with remainder.
     *
     * @details
     *           Given two polynomials <i>a</i> and a non-zero <i>b</i>
     *           there exits unique polynomials <i>q</i> and <i>r</i>
     *           such that
     *           \f$
     *            a=q\cdot b+r
     *           \f$
     *           with \f$\deg(r)<\deg(b)\f$. These two polynomials
     *           are determined by the method.
     *
     * @param q
     *           The quotient of the Euclidean division.
     *
     * @param r
     *           The remainder of the Euclidean division.
     *
     * @param a
     *           Numerator.
     *
     * @param b
     *           Denominator.
     *
     * @warning
     *           If <i>b</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If <i>q</i> and <i>r</i> are the same references,
     *           an error message is printed to <code>stderr</code>
     *           and the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::divRem
    ( BinaryPolynomial & q , BinaryPolynomial & r ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        // Preliminary tests

        if ( &q == &r ) {
            cerr << "BinaryPolynomial::divRem: "
                 << "Quotient and remainder must be of different references."
                 << endl;
            exit(EXIT_FAILURE);
        }

        if ( b.isZero() ) {
            cerr << "BinaryPolynomial::divRem: "
                 << "division by zero." << endl;
            exit(EXIT_FAILURE);
        }

        // Check if we have to make a temporary copy before division.
        if ( &b == &q || &b == &r ) {
            BinaryPolynomial tb(b);
            divRem(q,r,a,tb);
            return;
        }

        // School-book division method

        BinaryPolynomial tmp;
        BinaryPolynomial c(b);

        int n , m;
        n = a.deg();
        m = c.deg();
        r.assign(a);
        q.setZero();

        for ( int j = n-m ; j >= 0 ; j-- ) {
            if ( r.deg() == m+j ) {
                q.setCoeff(j);
                tmp.assign(c);
                tmp.leftShift(j);
                sub(r,r,tmp);
            }
        }

        r.degreeCouldBeSmaller();
    }

    /**
     * @brief
     *           Computes the inverse of a polynomial modulo
     *           another polynomial.
     *
     * @details
     *           The method attempts to find a polynomial <i>s</i>
     *           of degree smaller than \f$\deg(b)\f$ such that
     *           \f[
     *            1=c(X)\cdot a(X)~rem~b.
     *           \f]
     *           If such a polynomial exists, the function
     *           returns <code>true</code> and replaces the content
     *           of <i>c</i> correspondingly; otherwise, if the inverse
     *           does not exist, the function returns <code>false</code>
     *           and leaves the content of <i>c</i> unchanged.
     *
     * @param c
     *           Output of the inverse of <i>a</i> modulo <i>b</i>.
     *
     * @param a
     *           Input polynomial of which the inverse is computed.
     *
     * @param b
     *           Input modulus polynomial.
     *
     * @return
     *           <code>true</code> if the inverse have been determined
     *           successfully, in particular, if it exists; otherwise
     *           the function returns <code>false</code>.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> are both zero, the program
     *           prints an error message to <code>stderr</code> and
     *           exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    bool BinaryPolynomial::invMod
    ( BinaryPolynomial & c ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        BinaryPolynomial tmp0 , tmp1 , tmp2;

        xgcd(tmp0,tmp1,tmp2,a,b);

        if ( !tmp0.isOne() ) {
            return false;
        }

        divRem(tmp0,c,tmp1,b);

        return true;
    }

    /**
     * @brief
     *          Computes the greatest common divisor of two binary
     *          polynomials.
     *
     * @param g
     *          The greatest common divisors of <i>a</i> and <i>b</i>.
     *
     * @param a
     *          First polynomial.
     *
     * @param b
     *          Second polynomial.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::gcd
    ( BinaryPolynomial & g ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        BinaryPolynomial g1 , g2 , tmp;

        g1 = b; g = a;

        divRem(tmp,g2,g,g1);

        while ( !g1.isZero() ) {
            divRem(tmp,g2,g,g1);
            swap(g1,g2);
            swap(g,g2);
        }
    }

    /**
     * @brief
     *           Performs the extended Euclidean algorithm.
     *
     * @details
     *           The method determines the greatest common divisor
     *           <i>g</i> of the two polynomials <i>a</i> and <i>b</i>
     *           and co-prime polynomials <i>s</i> and <i>t</i> such that
     *           \f[
     *            g=s\cdot a+t\cdot b.
     *           \f]
     *
     * @param g
     *           Output of the greatest common divisor of <i>a</i>
     *           and <i>b</i>.
     *
     * @param s
     *           Output polynomial as described above.
     *
     * @param t
     *           Output polynomial as described above.
     *
     * @param a
     *           First input polynomial.
     *
     * @param b
     *           Second input polynomial.
     *
     * @warning
     *           If <i>g</i>, <i>s</i> or <i>t</i> are of the same
     *           reference, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> are both zero, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::xgcd
    ( BinaryPolynomial & g ,
      BinaryPolynomial & s , BinaryPolynomial & t ,
      const BinaryPolynomial & a , const BinaryPolynomial & b ) {

        if ( &g == &s || &g == &t || &s == &t ) {
            cerr << "BinaryPolynomial::xgcd: some output polynomials "
                 << "are of same reference." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial
            q , r0 , s0 , t0 , r1 , s1 , t1 , tmp;

        r0.assign(a); s0.setOne() ; t0.setZero();
        r1.assign(b); s1.setZero(); t1.setOne();

        while ( !r1.isZero() ) {

            divRem(q,tmp,r0,r1);

            swap(r0,r1);
            swap(s0,s1);
            swap(t0,t1);

            mul(tmp,q,r0);
            sub(r1,r1,tmp);

            mul(tmp,q,s0);
            sub(s1,s1,tmp);

            mul(tmp,q,t0);
            sub(t1,t1,tmp);
        }

        s.assign(s0);
        t.assign(t0);
        g.assign(r0);
    }

    /**
     * @brief
     *            Computes the partial greatest common divisor of two
     *            polynomials.
     *
     * @details
     *            This methods performs the traditional extended euclidean
     *            algorithm, i.e., determines \f$s,t\f$ such that
     *            \f$g=s\cdot a+t\cdot b\f$, until the degree of
     *            \f$g\f$ becomes smaller than \f$d\f$.
     *
     * @param g
     *            Will contain the partial greatest common divisor
     *            of <i>a</i> and <i>b</i>.
     *
     * @param s
     *            see details
     * @param t
     *            see details
     * @param a
     *            see details
     * @param b
     *            see details
     * @param d
     *            see details
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::pgcd
    ( BinaryPolynomial & g ,
      BinaryPolynomial & s ,
      BinaryPolynomial & t ,
      const BinaryPolynomial & a ,
      const BinaryPolynomial & b ,
      int d ) {

        BinaryPolynomial
            q , r0 , s0 , t0 , r1 , s1 , t1 , tmp;

        r0.assign(a); s0.setOne() ; t0.setZero();
        r1.assign(b); s1.setZero(); t1.setOne();

        while ( !r1.isZero() && r0.deg() >= d ) {

            divRem(q,tmp,r0,r1);

            swap(r0,r1);
            swap(s0,s1);
            swap(t0,t1);

            mul(tmp,q,r0);
            sub(r1,r1,tmp);

            mul(tmp,q,s0);
            sub(s1,s1,tmp);

            mul(tmp,q,t0);
            sub(t1,t1,tmp);
        }

        s.assign(s0);
        t.assign(t0);
        g.assign(r0);
    }

    /**
     * @brief
     *            Computes the <i>(n,n-k)</i>-Pade approximant to the
     *            polynomial.
     *
     * @details
     *            The method determines co-prime polynomials <i>r</i> and
     *            <i>t</i> with \f$\deg(r)<k\f$ and \f$\deg(t)\leq n-k\f$
     *            such that
     *            \f[
     *             \frac{r(X)}{t(X)}=g(X)~mod~X^n
     *            \f]
     *            where \f$t(X)\f$ is not divisible by \f$X\f$.
     *
     *            For more details on <em>Pade approximation</em> we refer
     *            to Section 5.9 in
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @param r
     *            Output numerator polynomial.
     *
     * @param t
     *            Output denominator polynomial.
     *
     * @param g
     *            Input polynomial.
     *
     * @param n
     *            Bounds the degree of <i>r</i>.
     *
     * @param k
     *            Bounds the degree of <i>t</i>, i.e., such that
     *            \f$\deg(t)\leq n-k\f$.
     *
     * @warning
     *            If <i>n</i> or <i>k</i> is negative, the degree
     *            of <i>g</i> is greater than or equals <i>n</i>, or
     *            if <i>k>n</i>, an error message is written to
     *            <code>stderr</code> and the program exits with
     *            status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::pade
    ( BinaryPolynomial & r , BinaryPolynomial & t ,
      const BinaryPolynomial & g , int n , int k ) {

        if ( n < 0 || g.deg() >= n || k > n ) {
            cerr << "BinaryPolynomial::pade: invalid arguments." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial m;
        m.setCoeff(n);

        BinaryPolynomial s;
        pgcd(r,t,s,g,m,k);
    }

    /**
     * @brief
     *           Computes the minimal polynomial of a polynomial modulo
     *           another polynomial.
     *
     * @details
     *           The minimal polynomial of <i>g</i> modulo <i>f</i> is
     *           defined as the non-zero polynomial <i>h</i> of minimal
     *           degree such that
     *           \f[
     *            h(g(X))~mod~f(X)=0.
     *           \f]
     *
     * @param h
     *           Output to contain the minimal polynomial of <i>g</i>
     *           modulo <i>f</i>.
     *
     * @param g
     *           Input polynomial.
     *
     * @param f
     *           Input modulus polynomial; must be of degree greater
     *           than 1.
     *
     * @warning
     *           If <i>f</i> is of degree smaller than 1, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::minPolyMod
    ( BinaryPolynomial & h ,
      const BinaryPolynomial & g , const BinaryPolynomial & f ) {

        int n = f.deg();

        if ( n < 1 ) {
            cerr << "BinaryPolynomial::minPolyMod: "
                 << "modulus must be of degree greater than or equals 1."
                 << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial a(n+n) , b(n+n) , tmp0(n+n) , tmp1(n+n);
        a.setOne();
        b.setOne();

        for ( int j = 1 ; j < n+n ; j++ ) {
            mul(tmp0,b,g);
            divRem(tmp1,b,tmp0,f);
            a.setCoeff(j,b.getCoeff(0));
        }

        BinaryPolynomial s(n) , t(n);
        pade(s,t,a,n+n,n);


        int d;
        if ( s.deg()+1 > t.deg() ) {
            d = s.deg()+1;
        } else {
            d = t.deg();
        }

        h = t.reverse(d);
    }

    /**
     * @brief
     *          Computes the composition of two binary polynomials.
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
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::eval
    ( BinaryPolynomial & h ,
      const BinaryPolynomial & f ,
      const BinaryPolynomial & g ) {

        if ( &h == &f ) {
            BinaryPolynomial tf(f);
            eval(h,tf,g);
            return;
        }

        if ( &h == &g ) {
            BinaryPolynomial tg(g);
            eval(h,f,tg);
            return;
        }

        BinaryPolynomial tmp(f.degree*(g.degree-1));

        h.setZero();
        int d = f.deg();

        // Run Horner's method
        for ( int j = d ; j >= 0 ; j--  ) {
            mul(tmp,h,g);
            swap(h,tmp);
            h.setCoeff(0,f.getCoeff(j)^h.getCoeff(0));
        }
    }

    /**
     * @brief
     *           Computes the composition of two polynomials being
     *           reduced by a modulus.
     *
     * @details
     *           The method performs <i>h(X)=f(g(X))</i> mod <i>m</i>.
     *
     * @param h
     *           Modular composition.
     *
     * @param f
     *           The outer polynomial.
     *
     * @param g
     *           The inner polynomial.
     *
     * @param m
     *           Modulus.
     *
     * @warning
     *           If <i>m</i> is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     *
     */
    void BinaryPolynomial::evalMod
    ( BinaryPolynomial & h , const BinaryPolynomial & f ,
      const BinaryPolynomial & g , const BinaryPolynomial & m ) {

        if ( &h == &f ) {
            BinaryPolynomial tf(f);
            evalMod(h,tf,g,m);
            return;
        }

        if ( &h == &g ) {
            BinaryPolynomial tg(g);
            evalMod(h,f,tg,m);
            return;
        }

        if ( &h == &m ) {
            BinaryPolynomial tm(m);
            evalMod(h,f,g,tm);
            return;
        }

        if ( m.isZero() ) {
            cerr << "BinaryPolynomial::modComp: division by zero." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial
                tmp0(f.degree+g.degree) , tmp1(f.degree+g.degree);

        h.setZero();
        int d = f.deg();

        // Run Horner's method ...
        for ( int j = d ; j >= 0 ; j--  ) {
            mul(tmp0,h,g);
            tmp0.setCoeff(0,f.getCoeff(j)^tmp0.getCoeff(0));
            // ... by incorporating a reduction.
            divRem(tmp1,h,tmp0,m);
        }
    }

    /**
     * @brief
     *          Computes the square root of a polynomial.
     *
     * @details
     *          If the input polynomial <i>f</i> is the square
     *          of another polynomial <i>g</i>, i.e., \f$f=g^2\f$, then
     *          the method will correctly determined <i>g</i>.
     *
     * @param g
     *          On output square root of the input polynomial <i>f</i>.
     *
     * @param f
     *          On input the square of the polynomial <i>g</i>.
     *
     * @see isSquare()
     *
     * @warning
     *          If <i>f</i> is not the square of a polynomial, the
     *          result of the method is undocumented.
     *
     * @warning
     *          If not enough memory could be provided, an error
     *          message is printed to <code>stderr</code> and the
     *          program exits with status 'EXIT_FAILURE'.
     */
    void BinaryPolynomial::squareRoot
        ( BinaryPolynomial & g , const BinaryPolynomial & f ) {

        if ( &g == & f ) {
            BinaryPolynomial tf(f);
            squareRoot(g,tf);
            return;
        }

        int m , n;

        n = f.deg();
        m = n/2;

        g.setZero();

        // If 'g_0+g_1*X+g_2*X^2+...+g_m*X^m' is the square root
        // of 'f', then 'f(X)=g_0+g_1*X^2+g_2*X^4+...+g_m*X^{2m}';
        // in the following we assume that only even coefficients
        // are non-zero and are thus able to revert the square.

        for ( int j = m ; j >= 0 ; j-- ) {
            if ( f.getCoeff(j*2) ) {
                g.setCoeff(j);
            }
        }
    }

    /**
     * @brief
     *          Tests whether this polynomial is the square of
     *          another binary polynomial.
     *
     * @return
     *          <code>true</code> if this polynomial is a square
     *          of another polynomial and <code>false</code> otherwise.
     */
    bool BinaryPolynomial::isSquare() const {

        // The polynomial is a square if and only if it has
        // non-zero coefficients at even positions

        static uint32_t oddBits = 4294967295u;

        int m = (this->degree+1)/32+((this->degree+1)%32?1:0);

        for ( int j = 0 ; j < m ; j++ ) {
            if ( (this->data[j] & oddBits) != 0 ) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief
     *          Determines whether the binary polynomial is irreducible.
     *
     * @details
     *          The implementation of the irreducible test follows
     *          the description of Algorithm 14.36 in
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @return
     *          <code>true</code> if the polynomial is irreducible and
     *          otherwise, if not irreducible, <code>false</code>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    bool BinaryPolynomial::isIrreducible() const {

        int n = deg();

        // The zero polynomial is divisible by each non-zero polynomial
        // and, thus, not irreducible.
        if ( n < 0 ) {
            return false;
        }

        // Units and linear factors cannot be further splitted and are, thus,
        // considered to be irreducible.
        if ( n <= 1) {
            return true;
        }

        // *******************************************
        // *** Step 1 (Algorithm 14.36 in [vzGth]) ***
        // *******************************************
        {
            BinaryPolynomial a , a2 , tmp;
            a.setCoeff(2);
            for ( int j = 0 ; j < n-1 ; j++ ) {
                mul(a2,a,a);
                divRem(tmp,a,a2,*this);
            }
            if ( !a.isX() ) {
                return false;
            }
        }

        BinaryPolynomial b , b2 , tmp , g;


        // *******************************************
        // *** Step 2 (Algorithm 14.36 in [vzGth]) ***
        // *******************************************

        int l = n;
        for ( int t = 2 ; t <= l ; t++ ) {
            if ( l % t == 0 ) {

                // for all prime divisors 't' of 'n' ...
                do {
                    l /= t;
                } while( l % t == 0 );

                // *******************************************
                // *** Step 3 (Algorithm 14.36 in [vzGth]) ***
                // *******************************************

                // ... compute 'b=x^(2^(n/t)) rem f'
                int m = n/t;
                b.setZero();
                b.setCoeff(2);
                for ( int j = 0 ; j < m-1 ; j++ ) {
                    mul(b2,b,b);
                    divRem(tmp,b,b2,*this);
                }

                // Compute 'b(X)-X'
                b.data[0] ^= 2;
                if ( b.degree == 1 ) {
                    b.degreeCouldBeSmaller();
                } else if ( b.degree <= 0 ) {
                    b.degree = 1;
                }

                // Compute 'g=gcd(b(X)-X,f(X))'
                gcd(g,b,*this);

                if ( !g.isOne() ) {
                    return false;
                }
            }
        }


        // *******************************************
        // *** Step 4 (Algorithm 14.36 in [vzGth]) ***
        // *******************************************
        return true;
    }

    /**
     * @brief
     *           Performs <i>Equal Degree Splitting</i> to split
     *           a proper factor from the polynomial as it is
     *           known to be a product of irreducible distinct
     *           polynomials of degree exactly <i>d</i>.
     *
     * @details
     *           For details of the implemented algorithm we refer
     *           to Algorithm 14.8 and Excercise 14.16 in
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @param d
     *           The degree of all irreducible divisors of the polynomial.
     *
     * @return
     *           A proper divisor of the polynomial (of which degree
     *           is a multiple of <i>d</i>).
     *
     * @warning
     *           If the polynomial is not the product of distinct
     *           irreducible polynomials of degree <i>d</i>, the
     *           function runs into undocumented behaviour.
     *
     * @warning
     *           If <i>d</i> is zero, negative or equals the degree
     *           of this polynomial or if <i>d</i> does not divide
     *           the degree, an error message is printed
     *           <code>stderr</code> an the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::equalDegreeSplitting( int d ) const {

        int n = deg();

        if ( d <= 0 || d >= n || n%d != 0 ) {
            cerr << "BinaryPolynomial::equalDegreeSplitting: "
                 << "invalid argument." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial a(n) , g(n) , b(n) , b2(n+n) , T(n) , tmp(n);

        // Iteration over Algorithm 14.8 until a proper factor reveals
        do {

            // ******************************************
            // *** Step 1 (Algorithm 14.8 in [vzGth]) ***
            // ******************************************

            // Choose $a\in\mathbb{F}_q[X]\setminu\mathbb{F}_q$
            // with $\deg a<n$ at random.
            do {
                a.random(n);
            } while( a.deg() <= 0 );

            // ******************************************
            // *** Step 2 (Algorithm 14.8 in [vzGth]) ***
            // ******************************************
            gcd(g,a,*this);
            if ( !g.isOne() ) {
                break;
            }

            // Compute the evaluation 'T' of the 'd'th trace polynomial at 'a'
            // modulo the polynomial 'f' to be factored; see Excercise
            // 14.36 in [vzGth]
            T.setZero();
            for ( int i = 0 ; i < d ; i++ ) {
                mul(b,a,a);
                for ( int j = 0 ; j < i ; j++ ) {
                    mul(b2,b,b);
                    divRem(tmp,b,b2,*this);
                }
                add(T,T,b);
            }

            // Try to split a common factor
            gcd(g,T,*this);

        } while ( g.deg() == 0 || g.deg() == n );

        return g;
    }

    /**
     * @brief
     *           Determines an irreducible divisors from the polynomial
     *           for which it is know that it is a product of irreducible
     *           and distinct polynomials of degree exactly <i>d</i>.
     *
     * @param d
     *           The degree of all irreducible divisors of the polynomial.
     *
     * @return
     *           A proper divisor of the polynomial of degree
     *           <i>d</i> (which is assumed to be irreducible).
     *
     * @warning
     *           If the polynomial is not the product of distinct
     *           irreducible polynomials of degree <i>d</i>, the
     *           function runs into undocumented behaviour.
     *
     * @warning
     *           If <i>d</i> is zero, negative or does not divide
     *           the degree of this polynomial, an error message is
     *           printed <code>stderr</code> an the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::splitIrreducibleEDF( int d ) const {

        BinaryPolynomial g(*this);

        // Determine proper factors of decreasing degree until equals 'd'
        while ( g.deg() != d ) {
            g = g.equalDegreeSplitting(d);
        }

        return g;
    }

    /**
     * @brief
     *           Determines an irreducible divisor from the
     *           polynomial which is assumed to be square-free.
     *
     * @return
     *           An irreducible divisor of this polynomial.
     *
     * @warning
     *           If this polynomial is not square-free, the function
     *           runs into undocumented behaviour.
     *
     * @warning
     *           If the degree of this polynomial is smaller than
     *           1, an error message is printed to <code>stderr</code>
     *           an the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::splitIrreducibleDDF() const {

        int n = deg();
        BinaryPolynomial h , g , tmp;

        h.setX();

        // Perform distinct-degree factorization ...

        for ( int d = 1 ; d < n ; d++ ) {

            mul(tmp,h,h);
            rem(h,tmp,*this);

            h.data[0] ^= 2;
            if ( h.degree == 1 ) {
                h.degreeCouldBeSmaller();
            } else if ( h.degree <= 0 ) {
                h.degree = 1;
            }

            gcd(g,h,*this);

            // ... until a 'd'-distinct-degree factorization
            // reveals.

            if ( !g.isOne() ) {

                // Then determine a proper irreducible factor
                // via Equals-degree factorization
                return g.splitIrreducibleEDF(d);
            }

            h.data[0] ^= 2;
            if ( h.degree == 1 ) {
                h.degreeCouldBeSmaller();
            } else if ( h.degree <= 0 ) {
                h.degree = 1;
            }
        }

        // The polynomial must be irreducible if no non-trivial
        // 'd'-distinct-degree factorization has been found.
        return *this;
    }

    /**
     * @brief
     *           Determines the product of all distinct irreducible
     *           divisors of this polynomial.
     *
     * @return
     *           The product of all distinct and irreducible divisors
     *           of this polynomial.
     *
     * @warning
     *           If the polynomial is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::squareFreeDecomposition() const {

        if ( isZero() ) {
            cerr << "BinaryPolynomial::squareFreeDecomposition: "
                 << "polynomial must be non-zero." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial g(*this) , tmp , h;

        // Splits non-square factor from the polynomial
        while ( g.isSquare() ) {
            squareRoot(tmp,g);
            swap(g,tmp);
        }

        // Determine formal derivative
        tmp.assign(g);
        tmp.derive();

        // The polynomial f/gcd(f,f') is the square-free decomposition
        gcd(tmp,tmp,g);
        div(h,g,tmp);

        return h;
    }

    /**
     * @brief
     *           Determines an irreducible factor of this polynomial.
     *
     * @details
     *           If the polynomial is irreducible, the returned polynomial
     *           will be equals to this polynomial. This method is useful
     *           if a single irreducible divisor of a polynomial is needed.
     *           For irreducibility tests it is more efficient to use
     *           the \link isIrreducible()\endlink function.
     *
     * @return
     *           An irreducible factor of this polynomial.
     *
     * @see isIrreducible()
     *
     * @warning
     *           If the polynomial is zero, an error message is printed
     *           to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryPolynomial BinaryPolynomial::splitIrreducible() const {
        return squareFreeDecomposition().splitIrreducibleDDF();
    }

    /**
     * @brief
     *           Updates the degree of the polynomial.
     *
     * @details
     *           After the coefficients stored in \link data\endlink
     *           change, the %degree of the polynomial can
     *           smaller than indicated by the \link degree\endlink
     *           field. The method updates the degree correspondingly.
     */
    void BinaryPolynomial::degreeCouldBeSmaller() {

        int m = (this->degree+1)/32+((this->degree+1)%32?1:0);

        this->degree = -1;

        for ( int j = m-1 ; j >= 0 ; j-- ) {
            int d = MathTools::numBits(this->data[j])-1;
            if ( d >= 0 ) {
                this->degree = j*32+d;
                break;
            }
        }
    }

    /**
     * @brief
     *            Prints a text representation of a polynomial to the
     *            specified output stream.
     *
     * @param out
     *            The output stream.
     *
     * @param f
     *            The polynomial of which a text representation is
     *            printed to <code>out</code>.
     *
     * @return
     *            A reference to <code>out</code> after the text
     *            representation has been written.
     */
    ostream & operator<<
    ( ostream & out , const BinaryPolynomial & f ) {

        if ( f.isZero() ) {
            out << "0";
        } else {

            int d = f.deg();
            bool alreadyPrinted = false;

            for ( int j = 0 ; j <= d ; j++ ) {
                bool c = f.getCoeff(j);
                if ( c ) {

                    if ( alreadyPrinted ) {
                        out << " + ";
                    }

                    if ( j >= 1 ) {
                        out << "Mod(1,2)*x";
                    } else {
                        out << "Mod(1,2)";
                    }

                    if ( j > 1 ) {
                        out << "^" << j;
                    }

                    alreadyPrinted = true;
                }
            }
        }

        return out;
    }
}
