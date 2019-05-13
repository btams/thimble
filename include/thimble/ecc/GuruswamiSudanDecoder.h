/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013 Benjamin Tams
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
 * @file GuruswamiSudanDecoder.h
 *
 * @brief
 *            Provides an implementation of the Guruswami-Sudan list decoding
 *            algorithm for decoding Reed-Solomon codes in original view.
 *
 * @author Benjamin Tams
 *
 * @see thimble::ReedSolomonCode
 * @see thimble::GuruswamiSudanDecoder
 */

#ifndef THIMBLE_GURUSWAMISUDANDECODER_H_
#define THIMBLE_GURUSWAMISUDANDECODER_H_

#include <stdint.h>
#include <ctime>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/math/numbertheory/SmallBinaryFieldBivariatePolynomial.h>

#ifdef THIMBLE_BUILD_DLL
template struct THIMBLE_DLL std::pair<double,double>;
template class THIMBLE_DLL std::allocator< thimble::SmallBinaryFieldPolynomial >;
template class THIMBLE_DLL std::vector< thimble::SmallBinaryFieldPolynomial , std::allocator< thimble::SmallBinaryFieldPolynomial > >;
#endif

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Instances of the class provide functions for performing
	 *            Guruswami-Sudan list decoding of Reed-Solomon codes given
	 *            in original view.
	 *
	 * @details
	 * @section sec_gsdecode Guruswami-Sudan Algorithm
	 *            Objects generated from
	 *            the \link GuruswamiSudanDecoder\endlink class
	 *            can be used to decode Reed-Solomon codes beyond half
	 *            the minimal distance. Assume we are given a set of \f$n\f$
	 *            finite field pairs \f$\{(x[i],y[i])\}\subset{\bf F}\f$
	 *            where the <i>x[i]</i> are all pairwise distinct.
	 *            Let \f$0<k\leq t\leq n\f$ be further integers. The
	 *            Reed-Solomon list decoding problem is to find all
	 *            polynomials \f$f\in{\bf F}[X]\f$ of degree \f$<k\f$ such
	 *            that \f$f(x[i])=y[i]\f$ for at least \f$t\f$ values of
	 *            \f$i\f$. This problem is also known as the <i>polynomial
	 *            reconstruction problem</i> of which assumed hardness the
	 *            <i>fuzzy vault scheme</i> draws its security.
	 *            <br><br>
	 *            If \f$t\geq\lceil (n+k)/2 \rceil\f$ there are efficient
	 *            algorithms which solve the polynomial reconstruction problem
	 *            and it is known that if a solution to the problem exists it
	 *            is even unique. An implementation for the such an algorithm
	 *            is (for example) provided by
	 *            \link
	 * ReedSolomonCode::decode(SmallBinaryFieldPolynomial&,const uint32_t*,const uint32_t*,int,int)
	 *            \endlink. But if \f$t\ll\sqrt{n\cdot(k-1)}\f$ the problem
	 *            is believed to be computationally hard. For the intermediate
	 *            cases where \f$\sqrt{n\cdot(k-1)}<t<\lceil(n+k)/2\rceil\f$
	 *            there are polynomial time algorithms which solve the
	 *            problem. A special class of such Reed-Solomon list decoding
	 *            algorithms are the Guruswami-Sudan list decoders:
	 *            <ul>
	 *             <li>
	 *              <b>M. Sudan (1997)</b>. Decoding of Reed-Solomon Codes
	 *              beyond the Error-Correction Bound. <i>Journal of
	 *              Complexity</i>, 13:180-193.
	 *             </li>
	 *             <li>
	 *              <b>V. Guruswami and M. Sudan (1998)</b>. Improved Decoding
	 *              of Reed-Solomon and Algebraic-Geometric Code.
	 *              <i>IEEE Trans. Intell. Transp. Syst.</i>, 45:1757-1767.
	 *             </li>
	 *            </ul>
	 *            In addition to the finite field and the parameters
	 *            <i>n, t,</i> and <i>k</i> the Guruswami-Sudan decoders
	 *            require an additional integral input, the multiplicity
	 *            <i>m>0</i>, which controls the number of errors <i>e=n-t</i>
	 *            the decoder can tolerate. More precisely, the higher the
	 *            multiplicity, the more errors the decoder can tolerate.
	 *            But note that the number of errors the decoder can tolerate
	 *            converges to
	 *            \f[
	 *             \lfloor n-\sqrt{n\cdot(k-1)}\rfloor
	 *            \f]
	 *            as the multiplicity <i>m</i> grows. Moreover, the higher the
	 *            multiplicity, the higher is the running time of the decoder.
	 *            So, the multiplicity should be chosen carefully.
	 *
	 *            <h2>Example</h2>
	 *            In the following, we give an example of how to use the
	 *            class for Reed-Solomon list decoding. Let
	 *            <pre>
	 *             SmallBinaryField gf(16);
	 *            </pre>
	 *            be a finite field of 2^16 element. For integers
	 *            <pre>
	 *             int n , t , k;
	 *             n = 224;
	 *             t = 48;
	 *             k = 9;
	 *            </pre>
	 *            we can create a random instance of the
	 *            <i>(n,t,k)</i>- polynomial reconstruction problem
	 *            (for simulation purposes) via
	 *            <pre>
	 *             uint32_t *x , *y;
	 *             x = (uint32_t*)malloc( n * sizeof(uint32_t) );
	 *             y = (uint32_t*)malloc( n * sizeof(uint32_t) );
	 *
	 *             FuzzyVaultTools::createRandomInstance(x,y,n,t,k,gf);
	 *            </pre>
	 *            such that there is a polynomial of degree (strictly) smaller
	 *            than <i>k</i> that interpolates (exactly) <i>t</i> of the
	 *            <i>n</i> pairs <i>(x[i],y[i])</i>. Attempting to find the
	 *            polynomial utilizing a classical Reed-Solomon decoder
	 *            <pre>
	 *             SmallBinaryFieldPolynomial f(gf);
	 *             bool success = ReedSolomonCode::decode(f,x,y,n,k);
	 *            </pre>
	 *            would fail, i.e. <code>success</code> will be
	 *            <code>false</code> and <i>f</i> will be left unchanged,
	 *            because the bound \f$t\geq\lceil(n+k)/2\rceil\f$ is not
	 *            satisfied. The implementation of the Guruswami-Sudan
	 *            list decoding, however,
	 *            <pre>
	 *             GuruswamiSudanDecoder gsDecoder;
	 *             int m = 3; //multiplicity
	 *
	 *             bool success = gsDecoder.decode(x,y,n,k,m,gf);
	 *            </pre>
	 *            will succeed (after approximately 0.5 seconds on a
	 *            3.2 Ghz desktop computer), i.e. <code>success</code> will be
	 *            <code>true</code>. The list of decoded polynomials can
	 *            be accessed via
	 *            <pre>
	 *             vector<SmallBinaryFieldPolynomial> list;
	 *
	 *             list = gsDecoder.getDecodedList();
	 *            </pre>
	 *            which will be of size 1 with very high probability even
	 *            though it is possible that it contains more than only
	 *            one polynomial. The list contains polynomials found by the
	 *            decoder interpolating at least <i>k+1</i> pairs
	 *            <i>(x[i],y[i])</i>. Furthermore, the decoded polynomials
	 *            are sorted w.r.t. the number of pairs they interpolate, i.e.
	 *            the first element in the list interpolates the most pairs
	 *            and the last element the fewest pairs. As a consequence, the
	 *            polynomial
	 *            <pre>
	 *             SmallBinaryFieldPolynomial f = list.front();
	 *            </pre>
	 *            will be the correct polynomial with very high probability.
	 */
	class THIMBLE_DLL GuruswamiSudanDecoder {

	public:

		/**
		 * @brief
		 *            Attempts to solve a Reed-Solomon list decoding problem.
		 *
		 * @details
		 *            The algorithm tries to find all polynomials of degree
		 *            \f$<k\f$ with coefficient in the field specified by
		 *            <code>gf</code> that interpolate a reasonable proportion
		 *            of pairs <i>(x[i],y[i])</i>. The number of interpolated
		 *            pairs depends on the specified multiplicity parameter
		 *            <i>m</i> but will be at least
		 *            \f$\lceil\sqrt{n\cdot (k-1)}\rceil\f$.
		 *            <br><br>
		 *            After the function terminated, the list of decoded
		 *            polynomials can be accessed via
		 *            <code>getDecodedList()</code>.
		 *            <br><br>
		 *            The implementation of the decoder proceeds as follows.
		 *            <h2>Interpolation Step</h2>
		 *            First, a non-zero bivariate polynomial
		 *            \f$Q(X,Y)\in{\bf F}(X,Y)\f$ is computed via
		 *            <code>interpolate()</code> which passes through the
		 *            points <i>(x[i],y[i])</i> with multiplicity at least
		 *            <i>m</i> and which is of minimal <i>(1,k-1)</i>-weighted
		 *            degree.
		 *            <h2>Root Step</h2>
		 *            It can be shown that if <i>m</i> is sufficiently large,
		 *            all polynomials \f$f^*(X)\f$ of degree \f$<k\f$ that
		 *            interpolate \f$>\sqrt{n\cdot(k-1)}\f$ pairs
		 *            <i>(x[i],y[i])</i> fulfill \f$Q(X,f^*(X))\equiv 0\f$.
		 *            Thus, next all polynomial \f$f^*\in{\bf F}[X]\f$ of
		 *            degree (strictly) smaller than <i>k</i> are computed
		 *            such that \f$Q(X,f^*(X))\equiv 0\f$. Therefore, the
		 *            class member function <code>roots()</code> is used.
		 *            <br>
		 *            Finally, those polynomials that do not interpolate more
		 *            than \f$\sqrt{n\cdot(k-1)}\f$ pairs are dismissed and
		 *            the remaining are sorted decreasingly w.r.t. the number
		 *            of pairs they interpolate.
		 *
		 * @param x
		 *            Reed-Solomon code locators.
		 *
		 * @param y
		 *            (Erroneous) evaluations of the message polynomials on
		 *            the code locators.
		 *
		 * @param n
		 *            Number of code pairs <i>(x[i],y[i])</i>.
		 *
		 * @param k
		 *            Bound of the length of polynomials that the decoder
		 *            outputs.
		 *
		 * @param m
		 *            Multiplicity parameter.
		 *
		 * @param gf
		 *            Specifies the finite field over which the decoder
		 *            operates.
		 *
		 * @return
		 *            <code>true</code> if at least one polynomial could be
		 *            reconstructed; otherwise <code>false</code>.
		 *
		 * @warning
		 *            If <i>n</i>, <i>k</i> or <i>m</i> are smaller than or
		 *            equals zero or if <i>k</i> is greater than <i>n</i> then
		 *            the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory can be provided, the function
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <i>x</i> or <i>y</i> do not contain at least <i>n</i>
		 *            valid elements in the field specified by <code>gf</code>
		 *            the function runs into undefined behavior.
		 */
		bool decode
		( const uint32_t *x , const uint32_t *y ,
		  int n , int k , int m , const SmallBinaryField & gf );

		/**
		 * @brief
		 *            Constructor.
		 */
		GuruswamiSudanDecoder();

		/**
		 * @brief
		 *            Destructor.
		 */
		~GuruswamiSudanDecoder();

		/**
		 * @brief
		 *            Access the list of recently decoded polynomials.
		 *
		 * @return
		 *            Decoded list.
		 *
		 * @warning
		 *            If the decoder was not used for a decoding attempt
		 *            before, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		const std::vector<SmallBinaryFieldPolynomial> &
		getDecodedList() const;

		/**
		 * @brief
		 *            Access the bivariate polynomial recently computed
		 *            in the interpolation step of the Guruswami-Sudan
		 *            decoder.
		 *
		 * @return
		 *            Result of the last interpolation step running
		 *            <code>decode()</code>.
		 *
		 * @warning
		 *            If the decoder was not used for a decoding attempt
		 *            before, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		const SmallBinaryFieldBivariatePolynomial &
		getBivariatePolynomial() const;

		/**
		 * @brief
		 *            Access the time in seconds the latest interpolation
		 *            step during <code>decode()</code> consumed.
		 *
		 * @return
		 *            The time of the last interpolation step.
		 *
		 * @warning
		 *            If the decoder was not used for a decoding attempt
		 *            before, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		double getInterpolationSecs() const;

		/**
		 * @brief
		 *            Access the time in seconds the latest root
		 *            step during <code>decode()</code> consumed.
		 *
		 * @return
		 *            The time of the last root step.
		 *
		 * @warning
		 *            If the decoder was not used for a decoding attempt
		 *            before, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		double getRootSecs() const;

		/**
		 * @brief
		 *            Access the time in seconds decoding attempt running
		 *            <code>decode()</code> consumed.
		 *
		 * @return
		 *            The overall time of the last decoding attempt.
		 *
		 * @warning
		 *            If the decoder was not used for a decoding attempt
		 *            before, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		double getOverallSecs() const;

		/**
		 * @brief
		 *            Implementation of the interpolation step in the
		 *            Guruswami-Sudan algorithm.
		 *
		 * @details
		 *            The function returns a non-zero bivariate polynomial
		 *            \f$Q(X,Y)\f$ with the property that the pairs
		 *            <i>(x[i],y[i])</i> pass through \f$Q\f$ with
		 *            multiplicity at least <i>m</i> and such that \f$Q\f$
		 *            is of minimal <i>(1,k-1)</i>-weighted degree.
		 *            <br><br>
		 *            The implementation of the function follows the
		 *            description of
		 *            <ul>
		 *             <li>
		 *              <b>Trifonov, P. (2010)</b>. Efficient interpolation in
		 *              the Guruswami-Sudan Algorithm. <i>IEEE Transactions on
		 *              Information Theory</i>, 56(9):4341-4349.
		 *             </li>
		 *            </ul>
		 *
		 * @param x
		 *            Reed-Solomon code locators.
		 *
		 * @param y
		 *            (Erroneous) evaluations of the message polynomials on
		 *            the code locators.
		 *
		 * @param n
		 *            Number of code pairs <i>(x[i],y[i])</i>.
		 *
		 * @param k
		 *            Bound of the length of polynomials that the decoder
		 *            outputs.
		 *
		 * @param m
		 *            Multiplicity parameter.
		 *
		 * @param gf
		 *            Specifies the finite field over which the decoder
		 *            operates.
		 *
		 * @return
		 *            A bivariate non-zero polynomial of minimal
		 *            <i>(1,k-1)</i>-weighted degree having the input pairs
		 *            <i>(x[i],y[i])</i> as roots with multiplicity at least
		 *            <i>m</i>.
		 *
		 * @warning
		 *            If <i>n</i>, <i>k</i> or <i>m</i> are smaller than or
		 *            equals zero or if <i>k</i> is greater than <i>n</i> then
		 *            the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory can be provided, the function
		 *            prints an error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <i>x</i> or <i>y</i> do not contain at least <i>n</i>
		 *            valid elements in the field specified by <code>gf</code>
		 *            the function runs into undefined behavior.
		 */
		SmallBinaryFieldBivariatePolynomial interpolate
		( const uint32_t *x , const uint32_t *y ,
		  int n , int k , int m , const SmallBinaryField & gf ) const;

		/**
		 * @brief
		 *            Implementation of the root step in the Guruswami-Sudan
		 *            list decoding algorithm.
		 *
		 * @details
		 *            Given a non-zero bivariate polynomial <i>Q</i>,
		 *            the function returns all polynomials \f$f^*\f$ of degree
		 *            (strictly) smaller than <i>k</i> such that
		 *            \f$Q(X,f^*(X))\equiv 0\f$.
		 *            <br><br>
		 *            The implementation of the function follows the
		 *            description of
		 *            <ul>
		 *             <li>
		 *              <b>McEliece, R. (2003)</b>. The Guruswami-Sudan
		 *              Decoding Algorithm for Reed-Solomon Codes.
		 *             </li>
		 *            </ul>
		 *            Also see
		 *            <ul>
		 *             <li>
		 *              <b>R. Roth and G. Ruckenstein (2000)</b>. Efficient
		 *              Decoding of Reed–Solomon Codes beyond Half the Minimum
		 *              Distance, <i>IEEE Trans. Inform. Theory</i>, vol. 46,
		 *              no. 1, pp. 246-257.
		 *             </li>
		 *            </ul>
		 *
		 * @param Q
		 *            Bivariate polynomial of which the roots are output.
		 *
		 * @param k
		 *            Size of the output roots.
		 *
		 * @return
		 *            Vector of polynomials that are roots of <i>Q</i>.
		 *
		 * @warning
		 *            If <i>Q</i> is zero or if <i>k</i> is smaller than or
		 *            equals 0, the function prints an error message to
		 *            <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not sufficient memory could be provided, the function
		 *            prints an error message to <code>stderr</code> and
		 *            exits with status 'EXIT_FAILURE'.
		 */
		std::vector<SmallBinaryFieldPolynomial> roots
		( const SmallBinaryFieldBivariatePolynomial & Q , int k ) const;

	private:
		/**
		 * @brief
		 *            The CPU time for the last run of the interpolation step
		 *            in the decoder.
		 *
		 * @details
		 *            The time (in seconds) can be accessed via
		 *            <code>getInterpolationSecs()</code>.
		 *
		 * @see getInterpolationSecs() const
		 */
		clock_t interpolationTime;

		/**
		 * @brief
		 *            The CPU time for the last run of the root step
		 *            in the decoder.
		 *
		 * @details
		 *            The time (in seconds) can be accessed via
		 *            <code>getRootSecs()</code>.
		 *
		 * @see getRootSecs() const
		 */
		clock_t rootTime;

		/**
		 * @brief
		 *            The CPU time for the last run the entire steps of
		 *            the decoder (i.e. interpolation and root step plus
		 *            minor intermediate steps).
		 *
		 * @details
		 *            The time (in seconds) can be accessed via
		 *            <code>getOverallSecs()</code>.
		 *
		 * @see getOverallSecs() const
		 */
		clock_t overallTime;

		/**
		 * @brief
		 *           Pointer to the bivariate polynomial as the result
		 *           of the last interpolation step.
		 *
		 * @details
		 *           If the member is <code>NULL</code> this indicates that
		 *           the no attempt to decode a Reed-Solomon list decoding
		 *           problem using this instance has been performed by now.
		 *           <br><br>
		 *           The latest interpolated polynomial can be accessed via
		 *           <code>getBivariatePolynomial()</code>.
		 *
		 * @see getBivariatePolynomial()
		 */
		SmallBinaryFieldBivariatePolynomial *Qptr;

		/**
		 * @brief
		 *           Contains all polynomials that have been decoded from
		 *           a last attempt to solve a Reed-Solomon list decoding
		 *           problem.
		 *
		 * @details
		 *           A constant reference to this member can be accessed
		 *           via <code>getDecodedList()</code>.
		 *
		 * @see getDecodedList()
		 */
		std::vector<SmallBinaryFieldPolynomial> decodedList;
	};
}


#endif /* THIMBLE_GURUSWAMISUDANDECODER_H_ */
