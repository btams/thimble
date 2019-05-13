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
 * @file ReedSolomonCode.h
 *
 * @brief
 *            Provides functionalities related to Reed-Solomon
 *            error-correcting codes.
 *
 * @author Benjamin Tams
 *
 * @see thimble::ReedSolomonCode
 */

#ifndef THIMBLE_REEDSOLOMONCODE_H_
#define THIMBLE_REEDSOLOMONCODE_H_

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Provides functionalities related to the Reed-Solomon
	 *            error-correcting codes.
	 *
	 * @details
	 * @section sec_rsdecode Reed-Solomon Decoder
	 *            The class \link ReedSolomonCode\endlink provides a mechanism
	 *            for decoding Reed-Solomon codes. Assume that we are given a
	 *            finite field
	 *            <pre>
	 *             SmallBinaryField gf = ...;
	 *            </pre>
	 *            Furthermore, assume that we are given a set of pairs
	 *            \f$\{(x[i],y[i])\}\f$ in the finite field
	 *            <pre>
	 *             uint32_t *x , *y;
	 *            </pre>
	 *            of length <i>n</i> where the <i>x[i]</i> are all pairwise
	 *            distinct (note that this assumes that <i>n</i> is smaller
	 *            than or equals the finite field's cardinality). Moreover,
	 *            let <i>t</i> pairs <i>(x[i],y[i])</i> lay on a common
	 *            polynomial of degree (strictly) smaller than <i>k</i>.
	 *            By running
	 *            <pre>
	 *             SmallBinaryFieldPolynomial f(gf);
	 *             bool success = ReedSolomonCode::decode(f,x,y,n,k);
	 *            </pre>
	 *            we may attempt to find the polynomial. Now, if
	 *            \f$t\geq\left\lceil\frac{n+k}{2}\right\rceil\f$ the variable
	 *            <code>success</code> will be <code>true</code> in which case
	 *            <i>f</i> is the polynomial interpolating <i>t</i> of the
	 *            pairs <i>(x[i],y[i])</i>. Otherwise, <code>success</code>
	 *            is <code>false</code> and <i>f</i> will be left unchanged.
	 *            <br><br>
	 *            There may be further functions/methods provided by the
	 *            class in the future but currently, there is only one decoder
	 *            of Reed-Solomon codes provided which is the one proposed in
	 *            <ul>
	 *             <li>
	 *              <b>Gao, S. (2002)</b>. A new algorithm for decoding
	 *              reed-solomon codes. In <i>Communications, Information and
	 *              Networks Security, V. Bhargava, H. V. Poor, V. Takrokh,
	 *              and S. Yoon</i>, pages 55-68. Kluwer.
	 *             </li>
	 *            </ul>
	 *            which is implemented in the
	 *            <code>ReedSolomonCode::gaodecode()</code> function.
	 *            <code>ReedSolomonCode::decode()</code> is just a wrapper
	 *            around <code>ReedSolomonCode::gaodecode()</code>. It is meant
	 *            to serve as a mechanism for choosing the most appropriate
	 *            decoder among those provided by the class in the future.
	 */
	class THIMBLE_DLL ReedSolomonCode {
	public:

		/**
		 * @brief
		 *            Attempts to decode a (possibly) perturbed
		 *            Reed-Solomon code <i>in original view</i> using an
		 *            algorithm by <b>Gao (2002)</b>.
		 *
		 * @details
		 *            If there is a polynomial <i>f</i> with coefficients
		 *            in the specified finite field that is of degree
		 *            (strictly) smaller than <i>k</i> such that
		 *            <i>f(x[i])=y[i]</i> for at least
		 *            \f$t\geq\left\lceil\frac{(n+k)}{2}\right\rceil\f$
		 *            pairs <i>i=0,...,n-1</i> then the polynomial <i>f</i>
		 *            is unique and will be discovered by the function. In
		 *            this case the function returns <code>true</code>.
		 *            Otherwise, if no such polynomial exists, the function
		 *            returns <code>false</code>.
		 *
		 * @param f
		 *            On success, <i>f</i> will contain the decoded message
		 *            polynomial.
		 *
		 * @param x
		 *            Locators of the Reed-Solomon code which should be all
		 *            pairwise distinct.
		 *
		 * @param y
		 *            Received vector of the (erroneous) Reed-Solomon
		 *            code.
		 *
		 * @param n
		 *            Length of the Reed-Solomon code.
		 *
		 * @param k
		 *            Size of the Reed-Solomon code.
		 *
		 * @return
		 *            <code>true</code> if decoding the received vector
		 *            was successful and <code>false</code> otherwise.
		 *
		 * @warning
		 *            In the following cases the function prints an
		 *            error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 *            <ul>
		 *             <li>
		 *              If <i>n<=0</i> or if <i>k<=0</i>.
		 *             </li>
		 *             <li>
		 *              If <i>k>n</i>.
		 *             </li>
		 *             <li>
		 *              If <i>x</i> does not contain elements that are
		 *              pairwise distinct.
		 *             </li>
		 *            </ul>
		 *
		 * @warning
		 *            If the first <i>n</i> elements of <i>x</i> or <i>y</i>
		 *            are not valid <code>uint32_t</code>s the behavior of the
		 *            function is undefined.
		 */
		static bool gaodecode
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y , int n , int k );

		/**
		 * @brief
		 *            Attempts to decode a (possibly perturbed)
		 *            Reed-Solomon code <i>in original view</i>.
		 *
		 * @details
		 *            &quot;Original view&quot; means that a Reed-Solomon code
		 *            is given as a set of points that lay on its message
		 *            polynomial.
		 *            <br><br>
		 *            Currently, the function just passes the arguments to
		 *            <code>gaodecode()</code> and returns its result.
		 *            This, may change in the future.
		 *
		 * @param f
		 *            On success, <i>f</i> will contain the decoded message
		 *            polynomial.
		 *
		 * @param x
		 *            Locators of the Reed-Solomon code which should be all
		 *            pairwise distinct.
		 *
		 * @param y
		 *            Received vector of the (erroneous) Reed-Solomon
		 *            code.
		 *
		 * @param n
		 *            Length of the Reed-Solomon code.
		 *
		 * @param k
		 *            Size of the Reed-Solomon code.
		 *
		 * @return
		 *            <code>true</code> if decoding the received vector
		 *            was successful and <code>false</code> otherwise.
		 *
		 * @warning
		 *            In the following cases the function prints an
		 *            error message to <code>stderr</code> and exits
		 *            with status 'EXIT_FAILURE'.
		 *            <ul>
		 *             <li>
		 *              If <i>n<=0</i> or if <i>k<=0</i>.
		 *             </li>
		 *             <li>
		 *              If <i>k>n</i>.
		 *             </li>
		 *             <li>
		 *              If <i>x</i> does not contain elements that are
		 *              pairwise distinct.
		 *             </li>
		 *            </ul>
		 *
		 * @warning
		 *            If the first <i>n</i> elements of <i>x</i> or <i>y</i>
		 *            are not valid <code>uint32_t</code>s the behavior of the
		 *            function is undefined.
		 */
		inline static bool decode
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y , int n , int k ) {
			return gaodecode(f,x,y,n,k);
		}
	};
}


#endif /* THIMBLE_REEDSOLOMONCODE_H_ */
