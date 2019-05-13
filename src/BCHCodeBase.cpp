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
 * @file BCHCodeBase.cpp
 *
 * @brief
 *            Implementation of a class for generating binary BCH codes that
 *            can be used for bit error-correction.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <vector>
#include <iostream>

#include <thimble/math/numbertheory/BinaryPolynomial.h>
#include <thimble/ecc/BCHCodeBase.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Creates a binary BCH code of given length that can
     *            tolerate at least the specified number of bit errors.
     *
     * @details
     *            The constructor takes a positive odd number <i>n</i>
     *            defining the block length of the code and a lower
     *            bound on the designed error tolerance. From the input
     *            a BCH code of block length <i>n</i> is generated that
     *            can tolerate at least the specified number of errors
     *            in a transmitted binary codeword being represented
     *            by \link BinaryPolynomial binary polynomials\endlink
     *            of degree smaller than <i>n</i>.
     *
     *            Parts of the BCH code generation have been implemented
     *            following the description of
     *            <ul>
     *             <li><b>[vzGth]</b> v.z. Gathen and Gerhard (2003).
     *              <i>Modern Computer Algebra</i>. Cambridge University
     *              Press, Cambridge (UK), 2nd edition.
     *             </li>
     *            </ul>
     *
     * @param n
     *            The block length of the BCH code.
     *
     * @param errorTolerance
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
    BCHCodeBase::BCHCodeBase( int n , int errorTolerance ) {

        if ( n < 0 ) {
            cerr << "BCHCode: block length must be positive." << endl;
            exit(EXIT_FAILURE);
        }

        if ( n%2 == 0 ) {
            cerr << "BCHCode: "
                 << "block length must not be a multiple of 2." << endl;
            exit(EXIT_FAILURE);
        }

        if ( errorTolerance < 0 ) {
            cerr << "BCHCode: "
                 << "designed minimal error tolerance must be positive."
                 << endl;
            exit(EXIT_FAILURE);
        }

        // Planned minimal distance
        int delta = 2*errorTolerance+1;

        if ( delta > n ) {
            cerr << "BCHCode: designed minimal distance is greater "
                 << "than block length." << endl;
            exit(EXIT_FAILURE);
        }

        // *******************************************
        // *** Step 1 (Algorithm 14.52 in [vzGth]) ***
        // *******************************************

        // Compute the 'n'th cyclotomic polynomial
        BinaryPolynomial phi = BinaryPolynomial::cyclo(n);

        // *******************************************
        // *** Step 2 (Algorithm 14.52 in [vzGth]) ***
        // *******************************************

        // Let 'f' be an irreducible factor of 'phi'
        this->f = phi.splitIrreducible();

        // *****************************************************
        // * BEGIN: Determine the generator polynomial as the  *
        // * minimal polynomial of all '(X mod f)^j' s.t. j is *
        // * greater or equals the planned minimal distance    *
        // * 'delta'. Note that follow a more straightforward  *
        // * approach as compared to Algorithm 14.52 in        *
        // * [vzGth].                                          *
        // *****************************************************

        this->g.setOne();
        BinaryPolynomial beta , beta_power , gj , tmp;
        beta.setX();
        beta_power.setOne();
        this->powers_of_X_mod_f.push_back(beta_power);
        for ( int j = 1 ; ; j++ ) {

            mul(beta_power,beta_power,beta);
            rem(beta_power,beta_power,this->f);
            this->powers_of_X_mod_f.push_back(beta_power);

            minPolyMod(gj,beta_power,this->f);

            gcd(tmp,this->g,gj);

            div(gj,gj,tmp);

            // Stop iteration as soon as 'j' has reached the
            // minimal distance.
            if ( j >= delta && gj.deg() != 0 ) {
                this->d = j;
                break;
            }

            mul(tmp,this->g,gj);
            swap(this->g,tmp);
        }

        // *****************************************************
        // * END: Generator polynomial and actual minimal      *
        // * distance have been determined                     *
        // *****************************************************

        // Not forget to keep track of the block length ...
        this->n = n;
        // ... and the dimension.
        this->k = n-this->g.deg();
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
    BCHCodeBase::BCHCodeBase( const BCHCodeBase & code ) {

        // Adopt the '='-operator of the members.
        this->n                 = code.n;
        this->k                 = code.k;
        this->d                 = code.d;
        this->g                 = code.g;
        this->f                 = code.f;
        this->powers_of_X_mod_f = code.powers_of_X_mod_f;
    }

    /**
     * @brief
     *            Destructor.
     *
     * @details
     *            Ensures that all data helt by this instance is
     *            freed.
     */
    BCHCodeBase::~BCHCodeBase() {
        // Every member's class has a destructor. It is actually not
        // necessary to implement the destructor except if we attempt
        // to emphasize that the memory is in fact freed.
    }

    /**
     * @brief
     *            Assignment operator.
     *
     * @details
     *            The data of this BCH code is overwritten by the
     *            data in <code>code</code>.
     *
     * @param code
     *            The BCH code assigned to this instance.
     *
     * @return
     *            A reference to this BCH code (after assignment).
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    BCHCodeBase & BCHCodeBase::operator=( const BCHCodeBase & code ) {

        // Adopt the '='-operator of the members. It is actually not
        // necessary to explicitly provide the assignment operator
        // except that we want to emphasize that assignment can safely
        // be performed.
        this->n                 = code.n;
        this->k                 = code.k;
        this->d                 = code.d;
        this->g                 = code.g;
        this->f                 = code.f;
        this->powers_of_X_mod_f = code.powers_of_X_mod_f;

        return *this;
    }

    /**
     * @brief
     *            Attempts to round the input bit sequence to its
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
     *            On input and output a polynomial of degree smaller than
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
    bool BCHCodeBase::round( BinaryPolynomial & c ) const {

        if ( c.deg() >= this->n ) {
            cerr << "BCHCode::round: "
                 << "length of input polynomial too large."
                 << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial r(this->n) , tmp(this->n);
        divRem(tmp,r,c,this->g);

        if ( r.isZero() ) {
            // 'c' is already a codeword; nothing to do
            return true;
        }

        int nu = getErrorTolerance();

        // Compute partial syndromes.
        // TODO: syndrome computation is currently the bottleneck in the
        // implementation
        vector<BinaryPolynomial> S(nu+nu);
        for ( int j = 0 ; j < nu+nu ; j++ ) {
            evalMod(S[j],r,this->powers_of_X_mod_f[j+1],this->f);
        }

        // Find the error-locator polynomial
        vector<BinaryPolynomial> Lambda = BCHCodeBase::berlekampMassey(S,this->f);

        // Find the error-locations which are given by the zeros of Lambda
        vector<int> locs = chienSearch(Lambda);
        if ( locs.size() + 1 != Lambda.size() ) {
            return false;
        }

        // Flip the bits at the error location; as the code is binary, this
        // causes
        for ( int i = 0 ; i < (int)locs.size() ; i++ ) {
            int j = locs.at(i);
            c.setCoeff(j,!c.getCoeff(j));
        }

        // Decoding successul
        return true;
    }

    /**
     * @brief
     *            Attaches checkbits to a message polynomial such
     *            that it becomes its encoding code polynomial.
     *
     * @details
     *            The encoding procedure implemented is systematic in the
     *            sense that the upper
     *            <i>k=</i>\link getDimension()\endlink coefficients of
     *            the code polynomial form the coefficients of the input
     *            message polynomial.
     *
     * @param c
     *            On input, the message polynomial of degree at most
     *            <i>k=</i>\link getDimension()\endlink; on output, a
     *            polynomial of degree smaller than
     *            <i>n=</i>\link getBlockLength()\endlink such that its
     *            upper <i>k</i> coefficients form the input message
     *            polynomial.
     *
     * @warning
     *            If the degree of the input message polynomial is larger
     *            than or equals <i>k=</i>\link getDimension()\endlink,
     *            an error message is printed to <code>stderr</code> and
     *            the program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void BCHCodeBase::encode( BinaryPolynomial & c ) const {

        if ( c.deg() >= this->k ) {
            cerr << "BCHCode::encode: "
                  << "length of message polynomial too large." << endl;
            exit(EXIT_FAILURE);
        }

        BinaryPolynomial p(this->n);

        // Make the upper 'k' coefficient form the input
        // message polynomial
        c.leftShift(this->n-this->k);

        // Check bits
        rem(p,c,this->g);

        // Add check bits
        sub(c,c,p);
    }

    /**
     * @brief
     *            Rounds the input polynomial to its nearest codeword
     *            and, if successful, removes the checkbits thereby
     *            obtaining the code polynomial's message polynomial.
     *
     * @details
     *            Apart from error-correction, this method performs the
     *            inverse of the \link encode()\endlink method.
     *
     * @param m
     *            On input, a polynomial of degree at most
     *            <i>n=</i>\link getBlockLength()\endlink. On successful
     *            decoding, a message polynomial of degree at most
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
     *            If <i>m</i> is of degree larger than or equals
     *            <i>n=</i>\link getBlockLength()\endlink, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    bool BCHCodeBase::decode( BinaryPolynomial & r ) const {

        if ( r.deg() >= this->n ) {
            cerr << "BCHCode::decode: "
                 << "length of received polynomial too large." << endl;
            exit(EXIT_FAILURE);
        }

        if ( !round(r) ) {
            return false;
        }

        r.rightShift(this->n-this->k);

        return true;
    }

    /**
     * @brief
     *           Creates a random code polynomial.
     *
     * @param c
     *           On output, a random code polynomial.
     *
     * @param tryRandom
     *           If <code>true</code>, the method uses a cryptographic
     *           number generator if available on the system; otherwise,
     *           the method wraps around the standard <code>rand()</code>
     *           function.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BCHCodeBase::random
    ( BinaryPolynomial & c , bool tryRandom ) const {

        BinaryPolynomial m(this->k);
        m.random(this->k,tryRandom);
        mul(c,m,this->g);
    }

    /**
     * @brief
     *            Implementation of the Berlekamp-Massey algorithm
     *            for finding the minimal polynomial of a linearly
     *            recurrent sequence in a finite field.
     *
     * @details
     *            If for a sequence \f$(S_j)_{j\geq 0}\f$ in a finite
     *            field \f$F\f$ there exist a non-zero polynomial
     * \f$\Lambda(X)=\lambda_0+\lambda_1\cdot X+...+\lambda_n\cdot X^t\f$
     *            such that
     *            \f[
     *             \sum_{0\leq j\leq t}\lambda_j\cdot S_j=0,
     *            \f]
     *            then \f$\Lambda(X)\f$ is called
     *            <em>characteristic polynomial</em> of <i>S</i>. If, in
     *            addition, the degree of \f$\Lambda(X)\f$ is minimal,
     *            then \f$\Lambda(X)\f$ is called
     *            <em>minimal polynomial</em> of <i>S</i>.
     *
     *            Given the first \f$2t\f$ coefficients of a (linearly
     *            recurrent) sequence in the finite field defined by
     *            <i>f</i>, this function computes its minimal polynomial
     *            \f$\Lambda(X)\f$ by returning a vector
     * <pre>
     *  vector<BinaryPolynomial> Lambda;
     * </pre>
     *            of length at most <i>t+1</i> where the coefficients
     *            \f$\lambda_j\f$ are represented by the
     *            <code>Lambda[j]</code> (modulo <i>f</i>).
     *
     *            This function is used in the \link round()\endlink
     *            function to find the error-locator polynomial on
     *            decoding a BCH code.
     *
     * @param S
     *            Linearly recurrent sequence (modulo <i>f</i>) in the
     *            finite field whose defining polynomial is specified
     *            via the parameter <i>f</i>.
     *
     * @param f
     *            Should be an irreducible polynomial of degree
     *            greater than 1.
     *
     * @return
     *            A vector defining the minimal polynomial of the
     *            linearly recurrent sequence given by <i>S</i> modulo
     *            <i>f</i>.
     *
     * @warning
     *            If <i>f</i> is not an irreducible polynomial of degree
     *            at least 1, the behaviour of the function is
     *            undocumented.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    vector<BinaryPolynomial> BCHCodeBase::berlekampMassey
    ( const std::vector<BinaryPolynomial> & S ,
      const BinaryPolynomial & f ) {

        int s = S.size();

        BinaryPolynomial Delta , tmp0 , tmp1 , tmp2;
        std::vector<BinaryPolynomial> N(s+1);

        //Lambda(X) = 1
        vector<BinaryPolynomial> Lambda;
        Lambda.assign(s+1,BinaryPolynomial());
        Lambda[0].setOne();

        // is equals the degree of 'Lambda'
        // if 'Lambda' contains non-zero higher coefficients,
        // they should be ignored
        int L = 0;

        // T(X) = X
        std::vector<BinaryPolynomial> T(s+1,BinaryPolynomial());
        T[1].setOne();

        for ( int l = 1 ; l <= s ; l++ ) {

            // Delta = S_k+Lambda_1S_{k-1}+...+Lambda_LS_{k-L-1}
            tmp1.setZero();
            for ( int i = 0 ; i <= L ; i++ ) {
                mul(tmp0,Lambda[i],S[l-i-1]);
                add(tmp1,tmp1,tmp0);
            }
            divRem(tmp2,Delta,tmp1,f);

            // don't change 'Lambda' if 'Delta' is zero; otherwise proceed
            // as follows:
            if ( !Delta.isZero() ) {

                // N(X) = Lambda(X)-Delta*T(X); new Delta(X) has 0 discrepancy
                for ( int i = 0 ; i <= s ; i++ ) {
                    mul(tmp0,Delta,T[i]);
                    divRem(tmp2,tmp1,tmp0,f);
                    sub(N[i],Lambda[i],tmp1);
                }

                if ( L < l-L ) {

                    L = l-L;

                    invMod(tmp0,Delta,f);

                    for ( int i = 0 ; i <= L ; i++ ) {
                        mul(tmp1,tmp0,Lambda[i]);
                        divRem(tmp2,T[i],tmp1,f);
                    }
                }

                // Essentially, initializes 'Lambda' by 'N'
                for ( int j = 0 ; j <= L ; j++ ) {
                    swap(Lambda[j],N[j]);
                }
            }

            if ( l < s ) {
                // T(X) <- X*T(X)
                for ( int i = s ; i > 0 ; i-- ) {
                    swap(T[i],T[i-1]);
                }
                T[0].setZero();
            }
        }

        Lambda.erase(Lambda.begin()+L+1,Lambda.end());

        return Lambda;
    }

    /**
     * @brief
     *             Determines the non-zero roots of a polynomial
     *             with coefficients in the field \f$GF(2^n)\f$ defined
     *             by the member <i>\link f\endlink</i>.
     *
     * @details
     *             The functions requires a polynomial
     *             \f[
     * \Lambda(X)=\lambda_0+\lambda_1\cdot X+\lambda_2\cdot X^2+...
     *    +\lambda_t\cdot X^t
     *             \f]
     *             such that \f$\lambda_j\in GF(2^n)\f$ (where
     *             \f$GF(2^n)\f$ is given by polynomial arithmetic modulo
     *             <i>\link f\endlink</i>) and determines its non-zero
     *             roots as follows:
     *             Let \f$\beta=(X~mod~f(X))\in GF(2^n)\f$; then the group
     *             of unity \f$GF(2^n)^\times=\langle\beta\rangle\f$ is
     *             generated by \f$\beta\f$ since a primite
     *             <i>\link n\endlink</i>th root of unity. If
     *             \f$\alpha_1,...,\alpha_\ell\f$ are the roots of
     *             \f$\Lambda(X)\f$, there exists non-negative integers
     *             \f$i_1,...,i_\ell\f$ such that
     *             \f$\alpha_j=\beta^{i_j}\f$. The function returns a vector
     *             containing the distinct exponents \f$i_j\f$
     *             (\f$j=1,...,\ell\f$).
     *
     *             The implementation of the function follows the
     *             <i>
     *              <a href="http://en.wikipedia.org/wiki/Chien_search"
     *               target="_blank">
     *               Chien search
     *              </a>
     *             </i>
     *             approach.
     *
     * @param Lambda
     *             An array of length <i>t+1</i>
     *             (where <i>t</i>\f$\leq\f$<i>\link d\endlink</i>)
     *             such that the entry <code>Lambda[j]</code>
     *             (modulo <i>\link f\endlink</i>) defines
     *             the <i>j</i>th coefficient \f$\lambda_j\f$ of the
     *             polynomial \f$\Lambda(X)\f$.
     *
     * @return
     *             A list of the distinct exponents \f$i_j\f$ of the
     *             non-zero roots \f$\alpha_j=\beta^{i_j}\f$ of the
     *             polynomial \f$\Lambda(X)\f$.
     *
     * @see <ul>
     *       <li>
     *        <b>Chien (1964)</b>.
     *        Cyclic Decoding Procedures for the Bose-Chaudhuri-Hocquenghem
     *        Codes.
     *        <i>IEEE Trans. Information Theory</i>, IT-10 (4): 357--363.
     *       </li>
     *      </ul>
     *
     * @warning
     *             If the requirement on the input <code>Lambda</code>
     *             is violated, a call of this function runs into
     *             undocumented behaviour.
     *
     * @warning
     *             If not enough memory could be provided, an error
     *             message is printed to <code>stderr</code> and the
     *             program exits with status 'EXIT_FAILURE'.
     */
    vector<int> BCHCodeBase::chienSearch
    ( const std::vector<BinaryPolynomial> & Lambda ) const {

        // To collect the error locations
        vector<int> locs;

        // The degree (or an upper bound if leading coefficients are zero)
        // of 'Lambda'
        int t = (int)Lambda.size()-1;

        // Temporary variables.
        BinaryPolynomial tmp0 , tmp1 , sum;


        // Initial table
        vector<BinaryPolynomial> Q(t+1,BinaryPolynomial());
        for ( int j = 0 ; j <= t ; j++ ) {
            mul(tmp0,Lambda[j],this->powers_of_X_mod_f[j]);
            divRem(tmp1,Q[j],tmp0,this->f);
        }

        // Chien search
        for ( int i = 0 ; i < this->n ; i++ ) {

            sum.setZero();
            for ( int j = 0 ; j <= t ; j++ ) {
                add(sum,sum,Q[j]);
            }

            // If the sum is zero, ...
            if ( sum.isZero() ) {

                // ..., we found an error lcoation.
                locs.push_back(this->n-i-1);
                if ( locs.size()+1 == Lambda.size() ) {
                    break;
                }
            }

            if ( i+1 < this->n ) {
                // Update the table
                for ( int j = 1 ; j <= t ; j++ ) {
                    mul(tmp0,Q[j],this->powers_of_X_mod_f[j]);
                    divRem(tmp1,Q[j],tmp0,this->f);
                }
            }
        }

        return locs;
    }
}

