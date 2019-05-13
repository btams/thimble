/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
 *
 *  Copyright 2014, 2015 Benjamin Tams
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
 * @file ProtectedMinutiaeRecord.h
 *
 * @brief
 *            Provides a class for generating and representing protected data
 *            from a mult-finger record containing absolutely pre-aligned
 *            minutiae.
 *
 * @details
 * @section sec_mffv Minutiae Template Protection for Multiple Fingerprints
 *
 * Objects generated from the \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink class
 * can be used to protect absolutely pre-aligned minutiae templates from multiple
 * fingerprint stored in a \link thimble::MinutiaeRecord MinutiaeRecord\endlink.
 * The implementation of \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink follows the description of
 * <ul>
 *  <li>
 *   <b>B. Tams (2015)</b>. Unlinkable Minutiae-Based Fuzzy Vault for
 *   Multiple Fingerprints, <i>IET Biometrics</i>, to appear
 *   (<a href="http://www.stochastik.math.uni-goettingen.de/preprints/ffv.pdf">
 *    Preprint</a>).
 *  </li>
 * </ul>
 * and its interface is very similar to the interface of objects generated
 * from the the \link thimble::ProtectedMinutiaeTemplate\endlink class. We
 * therefore refer to the tutorial @ref sec_ffv and recommend to get familiar
 * with how single absolutely pre-aligned minutiae templates can be protected
 * with THIMBLE first. In the following we describe the difference.
 *
 * @subsection sec_mffv_record Enrollment/Verification with Minutiae Records
 *
 * To protect absolutely pre-aligned minutiae templates from multiple
 * fingerprints the \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink class makes use of
 * the \link thimble::MinutiaeRecord MinutiaeRecord\endlink class which
 * can store multiple minutiae templates encoded by objects generated from
 * the \link thimble::MinutiaeView MinutiaeView\endlink. These minutiae
 * templates must be absolutely pre-aligned (e.g., see the tutorial
 * @ref sec_ffv_prealign).
 *
 * Furthermore, the minutiae templates stored
 * in the record must not encode an unknown finger position, i.e.,
 * their member function \link thimble::MinutiaeView::getFingerPosition()
 * getFingerPosition()\endlink must not return \link thimble::UNKNOWN_FINGER
 * UNKNOWN_FINGER\endlink. Moreover, the minutiae records must contain
 * minutiae templates indicating pair-wise different finger positions.
 *
 * Another difference of the \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink class is that it does not need to
 * be initialized before enrollment while
 * the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class requires initialization of
 * the fingerprint image's width and height prior enrollment: An object
 * from the \link thimble::MinutiaeRecord MinutiaeRecord\endlink class
 * encodes the width, height, and even scanning resolution from the
 * original fingerprint images; consequently, they are adopted by the
 * enrollment and verification functions of
 * a \link thimble::ProtectedMinutiaeRecord ProtectedMinutiaeRecord\endlink
 * object.
 *
 * @subsection sec_mffv_enrol Enrollment
 *
 * Given a minutiae record
 * <pre>
 *    MinutiaeRecord reference;
 * </pre>
 * we may generate protected data from it via
 * <pre>
 *    ProtectedMinutiaeRecord vault;
 *    vault.enroll(record)
 * </pre>
 * which return <code>true</code> on successful enrollment and otherwise
 * <code>false</code>
 *
 * @subsection sec_mffv_verify Verification
 *
 * Given a query record
 * <pre>
 *     MinutiaeRecord query;
 * </pre>
 * we can verify its authenticity versus a successfully enrolled
 * <pre>
 *     ProtectedMinutiaeRecord vault;
 * </pre>
 * by running
 * <pre>
 *     bool verified = vault.verify(query);
 *
 *     if ( verified ) {
 *        cout << "ACCEPT" << endl;
 *     } else {
 *        cout << "REJECT" << endl;
 *     }
 * </pre>
 *
 * @subsection sec_mffv_other Other Functionalities
 *
 * Most of the other functions that are provided by objects of
 * the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class (e.g., writing and reading the
 * protected data and @ref sec_ffv_password and @ref sec_ffv_slowdown)
 * are provided for objects from the \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink class in the same way.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_PROTECTEDMINUTIAERECORD_H_
#define THIMBLE_PROTECTEDMINUTIAERECORD_H_

#include <stdint.h>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/security/SHA.h>
#include <thimble/math/Permutation.h>
#include <thimble/math/numbertheory/BigInteger.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/finger/MinutiaeRecord.h>

#ifdef THIMBLE_BUILD_DLL
template struct THIMBLE_DLL std::pair<double,double>;
template class THIMBLE_DLL std::allocator<std::pair<double,double> >;
template class THIMBLE_DLL std::vector<std::pair<double,double>,std::allocator<std::pair<double,double> > >;
#endif

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            A class providing mechanisms for generating
	 *            protected data from multi-finger records containing
	 *            absolutely pre-aligned minutiae templates.
	 */
	class THIMBLE_DLL ProtectedMinutiaeRecord {

	public:

		/**
		 * @brief
		 *            Standard constructor.
		 */
		ProtectedMinutiaeRecord();

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @param vault
		 *            A reference to the protected minutiae record of which
		 *            a copy is created.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 */
		ProtectedMinutiaeRecord( const ProtectedMinutiaeRecord & vault );

		/**
		 * @brief
		 *            Destructor.
		 */
		~ProtectedMinutiaeRecord();

		/**
		 * @brief
		 *            Assignment operator.
		 *
		 * @param vault
		 *            A reference to the protected minutiae record of which
		 *            a copy is assigned to this protected minutiae record
		 *            object.
		 *
		 * @return
		 *            A reference to this protected minutiae record object.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 */
		ProtectedMinutiaeRecord &operator=
		( const ProtectedMinutiaeRecord & vault );

		/**
		 * @brief
		 *            Assignment operator (procedural version).
		 *
		 * @param vault
		 *            A reference to the protected minutiae record of which
		 *            a copy is assigned to this protected minutiae record
		 *            object.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 */
		void assign( const ProtectedMinutiaeRecord & vault );

		/**
		 * @brief
		 *            Prepares this protected minutiae record to protect
		 *            multiple absolutely pre-aligned minutiae templates
		 *            of which (un-aligned) minutiae originated from an
		 *            image of the specified width and height that has
		 *            been scanned at the specified resolution.
		 *
		 * @details
		 *            It is actually not necessary to call this method
		 *            before passing an absolutely pre-aligned minutiae
		 *            record to the \link enroll()\endlink function (note,
		 *            that this is different
		 *            for \link ProtectedMinutiaeRecord\endlink objects).
		 *            A \link MinutiaeRecord\endlink already encodes the
		 *            width, height and resolution (in pixel per centimeter)
		 *            of the original fingerprint images and
		 *            this \link initialize()\endlink method is implicitly
		 *            called by calling the \link enroll()\endlink method.
		 *            However, if the user of this class wants to call
		 *            this object's %quantize() functions, he must call
		 *            this \link initialize()\endlink method first.
		 *
		 * @param width
		 *            The specified pixel width of the original fingerprint
		 *            images.
		 *
		 * @param height
		 *            The specified pixel height of the original fingerprint
		 *            image.
		 *
		 * @param dpi
		 *            The resolution in dots per inch of the original
		 *            fingerprint image.
		 *
		 * @warning
		 *            If the specified width, height or resolution are
		 *            negative or zero, an error message will be printed
		 *            to <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 */
		void initialize( int width , int height , int dpi = 500 );

		/**
		 * @brief
		 *            Creates protected data from a record containing
		 *            multiple absolutely pre-aligned minutiae templates.
		 *
		 * @details
		 *            The function passes the input minutiae templates
		 *            through a quantization process and protects the
		 *            quantizations via the improved fuzzy vault scheme
		 *            of which data will be stored by this instance.
		 *
		 *            If not sufficient information could be extracted
		 *            from the input minutiae templates to allow successful
		 *            verification, i.e., a <em>failure to enroll</em>,
		 *            the function returns <code>false</code>;
		 *            otherwise, if successful, the function returns
		 *            <code>true</code>.
		 *
		 * @param record
		 *            A reference to a valid \link MinutiaeRecord\endlink
		 *            storing the multiple to-be-protected absolutely
		 *            pre-aligned minutiae templates.
		 *
		 * @return
		 *            <code>true</code> if the creation was successful;
		 *            otherwise, the function returns <code>false</code>.
		 *
		 * @warning
		 *            If <code>record</code> contains \link MinutiaeView
		 *            minutiae views\endlink suggesting an unknown finger
		 *            position, i.e.,
		 *            where \link MinutiaeView::getFingerPosition()
		 *            getFingerPosition()\endlink ==
		 *            \link UNKNOWN_FINGER\endlink, or if
		 *            <code>record</code> contains \link MinutiaeView
		 *            minutiae views\endlink with
		 *            duplicate \link MinutiaeView::getFingerPosition()
		 *            finger positions\endlink, an error message is printed
		 *            to <code>stderr</code> and the program exits with status
		 *            'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If \link isEnrolled()\endlink is
		 *            <code>true</code>, an error message is printed to
		 *            <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If the horizontal resolution of the record is different
		 *            from the vertical resolution, i.e.,
		 *            if \link MinutiaeRecord::getHorizontalResolution()
		 *            record.getHorizontalResolution()\endlink !=
		 *            \link MinutiaeRecord::getVerticalResolution()
		 *            record.getVerticalResolution()\endlink, then an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool enroll( const MinutiaeRecord & record );

		/**
		 * @brief
		 *             Performs a verification procedure to "compare"
		 *             the multiple absolutely pre-aligned minutiae templates
		 *             with the minutiae templates protected by this object.
		 *
		 * @details
		 *             The function passes the input minutiae templates
		 *             through a quantization process from which an instance
		 *             of the Reed-Solomon decoding problem is generated which
		 *             is attempted to be solved using
		 *             the \link ProtectedMinutiaeTemplate::decode()
		 *             decode()\endlink function. If successful, the
		 *             verificationis deemed to be successful and, otherwise,
		 *             failed.
		 *
		 * @param record
		 *             A minutiae record containing absolutely pre-aligned
		 *             query minutiae templates.
		 *
		 * @return
		 *             <code>true</code> if the verification procedure
		 *             was successful; otherwise, if the verification
		 *             procedure was not successful, the function returns
		 *             <code>false</code>.
		 *
		 * @warning
		 *            If <code>record</code> contains \link MinutiaeView
		 *            minutiae views\endlink suggesting an unknown finger
		 *            position, i.e.,
		 *            where \link MinutiaeView::getFingerPosition()
		 *            getFingerPosition()\endlink ==
		 *            \link UNKNOWN_FINGER\endlink, or if
		 *            <code>record</code> contains \link MinutiaeView
		 *            minutiae views\endlink with
		 *            duplicate \link MinutiaeView::getFingerPosition()
		 *            finger positions\endlink, an error message is printed
		 *            to <code>stderr</code> and the program exits with status
		 *            'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If this \link ProtectedMinutiaeRecord\endlink
	     *             does not represent a successfully enrolled
	     *             decrypted records, i.e.,
	     *             if \link ProtectedMinutiaeRecord::isDecrypted()
	     *             isDecrypted()\endlink returns <code>false</code>,
	     *             the function prints an error message to
	     *             <code>stderr</code> and exits with status
	     *             'EXIT_FAILURE'.
	     *
		 * @warning
		 *             If the horizontal resolution of the record is different
		 *             from the vertical resolution, i.e.,
		 *             if \link MinutiaeRecord::getHorizontalResolution()
		 *             record.getHorizontalResolution()\endlink !=
		 *             \link MinutiaeRecord::getVerticalResolution()
		 *             record.getVerticalResolution()\endlink, then an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool verify( const MinutiaeRecord & record ) const;

		/**
		 * @brief
		 *             Encrypts this instance using a key.
		 *
		 * @details
		 *             After successful enrollment
		 *             (via  \link enroll()\endlink),
		 *             the \link ProtectedMinutiaeRecord\endlink can
		 *             additionally protected via a (user-specific)
		 *             key derived from a PIN or password. Consequently,
		 *             the security level of the
		 *             \link ProtectedMinutiaeRecord\endlink can be
		 *             improved, though this requires the
		 *             \link ProtectedMinutiaeRecord\endlink to be
		 *             decrypted correctly for a successful verification,
		 *             thereby requiring that the correct can be
		 *             reproduced for successful verification.
		 *
		 * @param key
		 *             The \link AES128\endlink key used to encrypt the
		 *             \link ProtectedMinutiaeRecord\endlink.
		 *
		 * @warning
		 *             If no minutiae template has been enrolled with
		 *             this \link ProtectedMinutiaeRecord\endlink
		 *             (i.e., if \link isEnrolled()\endlink
		 *             returns <code>false</code>)
		 *             or if it contains already encrypted data (i.e.,
		 *             if \link containsEncryptedData()\endlink
		 *             returns <code>true</code>)
		 *             an error message is printed to <code>stderr</code>
		 *             and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void encrypt( const AES128 & key );

		/**
		 * @brief
		 *             Decrypts this instance using a key.
		 *
		 * @details
		 *             This function is to decrypt a
		 *             \link ProtectedMinutiaeRecord\endlink that has been
		 *             encrypted before using \link encrypt()\endlink.
		 *
		 *             If the \link ProtectedMinutiaeRecord\endlink has been
		 *             encrypted with the key \f$\kappa\f$, then decryption
		 *             with \f$\kappa\f$ will correctly reveal the original
		 *             data; otherwise, if decrypting with a different
		 *             key \f$\kappa'\neq\kappa\f$, this instance will hold
		 *             a candidate for the original data that is (with very
		 *             high probability) different from the original data.
		 *             To date no efficient way for distinguishing incorrect
		 *             candidate data from the correct data is known. Hence,
		 *             the difficulty to break an encrypted
		 *             \link ProtectedMinutiaeRecord\endlink equals the product
		 *             of the difficulty in guessing a  minutiae template being
		 *             of sufficient similarity to the protected minutiae template
		 *             and the difficulty in guessing the correct password.
		 *
		 * @param key
		 *             The \link AES128\endlink key used to decrypt this
		 *             \link ProtectedMinutiaeRecord\endlink.
		 *
		 * @warning
		 *             If this instance does not hold encrypted data
		 *             (i.e., if \link containsEncryptedData()\endlink
		 *             returns <code>false</code>), an error message is
		 *             printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @see encrypt()
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void decrypt( const AES128 & key );

		/**
		 * @brief
		 *             Uses a record of absolutely pre-aligned minutiae
		 *             to open the vault.
		 *
		 * @details
		 *             This a special version of
		 *             the \link verify()\endlink function (in
		 *             fact, \link verify()\endlink
		 *             wraps around this function) in which the
		 *             secret polynomial that is used to hide the
		 *             minutiae record protected by this instance is revealed
		 *             on successful verification.
		 *
		 * @param f
		 *             Will be set to the secret polynomial on successful
		 *             verification; otherwise the content will be left
		 *             unchanged
		 *
		 * @param record
		 *             The query record of absolutely pre-aligned
		 *             minutiae template.
		 *
		 * @return
		 *             <code>true</code> if the verification procedure was
		 *             successful; otherwise, if the verification was not
		 *             successful, the function returns <code>false</code>.
		 *
		 * @warning
		 *             If this \link ProtectedMinutiaeRecord\endlink
	     *             does not represent a successfully enrolled
	     *             (and decrypted) record, i.e.,
	     *             if \link isEnrolled()\endlink
	     *             or \link isDecrypted()\endlink return
	     *             <code>false</code>, the function
	     *             prints an error message to <code>stderr</code> and
	     *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>record</code> contains \link MinutiaeView
		 *             minutiae views\endlink suggesting an unknown finger
		 *             position, i.e.,
		 *             where \link MinutiaeView::getFingerPosition()
		 *             getFingerPosition()\endlink ==
		 *             \link UNKNOWN_FINGER\endlink, or if
		 *             <code>record</code> contains \link MinutiaeView
		 *             minutiae views\endlink with
		 *             duplicate \link MinutiaeView::getFingerPosition()
		 *             finger positions\endlink, an error message is printed
		 *             to <code>stderr</code> and the program exits with status
		 *             'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If the horizontal resolution of the record is different
		 *             from the vertical resolution, i.e.,
		 *             if \link MinutiaeRecord::getHorizontalResolution()
		 *             record.getHorizontalResolution()\endlink !=
		 *             \link MinutiaeRecord::getVerticalResolution()
		 *             record.getVerticalResolution()\endlink, then an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool open
		( SmallBinaryFieldPolynomial & f ,
		  const MinutiaeRecord & record ) const;

		/**
		 * @brief
		 *             Uses a feature set as a query to open the vault.
		 *
		 * @details
		 *             This is a variant of the \link open()\endlink
		 *             function in which a set of quantized features (e.g.,
		 *             derived from a records absolutely pre-aligned minutiae) is
		 *             used to open the vault. Therein, a quantized minutia is
		 *             encoded as an element of the finite field that can be
		 *             accessed via the \link getField()\endlink function.
         *
         *             To obtain the feature set of a record with absolutely
         *             pre-aligned minutiae templates,
         *             the \link quantize()\endlink method can be used.
		 *
		 * @param f
		 *             Will be set to the secret polynomial on successful
		 *             verification; otherwise the content will be left
		 *             unchanged
		 *
		 * @param B
		 *             The feature set.
		 *
		 * @param t
		 *             Number of elements in the feature set.
		 *
		 * @return
		 *             <code>true</code> if the verification procedure was
		 *             successful; otherwise, if the verification was not
		 *             successful, the function returns <code>false</code>.
		 *
		 * @warning
		 *             If this \link ProtectedMinutiaeRecord\endlink
	     *             does not represent a successfully enrolled (and
	     *             decrypted) protected minutiae template, i.e., if
	     *             if \link isEnrolled()\endlink
	     *             or \link isDecrypted()\endlink return
	     *             <code>false</code>, the function prints an error
	     *             message to <code>stderr</code> and exits with status
	     *             'EXIT_FAILURE'.
	     *
	     * @warning
	     *             If <code>B</code> does not contain distinct integers
	     *             of type <code>uint32_t</code> being elements of the
	     *             finite field represented by \link getField()\endlink,
	     *             the function runs into undocumented behavior.
	     *
	     * @warning
	     *             If <code>B</code> does not contain at least
	     *             <code>t</code> well-defined distinct feature elements,
	     *             the function runs into undocumented behavior.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool open
		( SmallBinaryFieldPolynomial & f , const uint32_t *B , int t ) const;

		/**
		 * @brief
		 *             Implementation of a Guruswami-Sudan algorithm based
		 *             decoder with which the verification process of a
		 *             this multiple minutiae template based fuzzy vault is
		 *             performed.
		 *
		 * @details
		 *             This program tries to find a polynomial <i>f</i>
		 *             with hash value equals <code>hash</code> from the
		 *             finite field pairs <i>(x[i],y[i])</i>. In the first
		 *             step, a classical Reed-Solomon decoder is used
		 *             and if it successfully output the desired polynomial,
		 *             the function returns <code>true</code>. Otherwise,
		 *             a attempt to find the desired polynomial with a
		 *             Guruswami-Sudan algorithm is performed successfully
		 *             iterating multplicities from 1 to <i>m</i> until
		 *             the desired polynomial could be decoded. If the
		 *             desired polynomial could be found, it will be
		 *             stored in <code>f</code> and the function returns
		 *             <code>true</code>; otherwise, the function
		 *             returns <code>false</code>.
		 *
		 * @see GuruswamiSudanDecoder
		 * @see ReedSolomonCode
		 *
		 * @param f
		 *             will contain the recovery of the secret polynomial on
		 *             successful decoding; otherwise, will be left unchanged.
		 *
		 * @param x
		 *             Abscissa coordinates
		 *
		 * @param y
		 *             Ordinate coordinates
		 *
		 * @param t
		 *             Number of valid finite field pairs <i>(x[j],y[j])</i>.
		 *
		 * @param k
		 *             Length of the secret polynomial whose degree is smaller
		 *             than <i>k</i>.
		 *
		 * @param hash
		 *             160-bit hash value of the correct secret polynomial
		 *             <i>f</i>, that has been originally computed via an
		 *             equivalent to
		 *             <code>
		 *              SHA().hash(hash,f.getData(),f.deg()+1)
		 *             </code>
		 *
		 * @param m
		 *             Maximal bound of multiplicity with which a
		 *             Guruswami-Sudan algorithm will is run during a
		 *             decoding attempt.
		 *
		 * @return
		 *             <code>true</code> if the decoding attempt was successful;
		 *             otherwise, the function returns <code>false</code>.
		 *
		 * @warning
		 *             If the <i>x[i]</i> contain elements not encoding a
		 *             element in the finite field represented
		 *             by \link SmallBinaryFieldPolynomial::getField()
		 *             f.getField()\endlink or of <i>x[i]</i> contain
		 *             duplicates, the decoder runs into undocumented
		 *             behavior.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool decode
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int t , int k , const uint8_t hash[20] , int m ) const;

		/**
		 * @brief
		 *             Computes the quantization of an (absolutely pre-aligned)
		 *             minutia stemming from a finger with specified position.
		 *
		 * @param minutia
		 *             The absolutely pre-aligned minutia of which the
		 *             quantization is computed.
		 *
		 * @param L
		 *             The position of the finger from which the absolutely
		 *             pre-aligned minutia has been measured.
		 *
		 * @return
		 *             The quantization of <code>minutia</code> encoded as a
		 *             finite field element.
		 *
		 * @warning
		 *             If this object has not been initialized (either
		 *             explicitly via \link initialize()\endlink or
		 *             implicitly via \link enroll()\endlink), i.e.,
		 *             \link isInitialized()\endlink returns
		 *             <code>false</code>, an error message is printed
		 *             to <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>L</code> equals \link UNKNOWN_FINGER\endlink,
		 *             an error message is printed to <code>stderr</code> and
		 *             the program exits with status 'EXIT_FAILURE'.
		 */
		uint32_t quantize
		( const Minutia & minutia , FINGER_POSITION_T L ) const;

		/**
		 * @brief
		 *             Computes the quantization set of an
		 *             (absolutely pre-aligned) minutia template.
		 *
		 * @details
		 *             The function extracts at
		 *             most \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             quantizations from <code>view</code> (preferring
		 *             those having the highest quality estimation;
		 *             see \link Minutia::getQuality()\endlink).
		 *             Consequently, the array should be able to contain
		 *             at least \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             integers of type <code>uint32_t</code>.
		 *
		 * @param A
		 *             Array in which the quantization elements will be
		 *             stored.
		 *
		 * @param view
		 *             The absolutely pre-aligned minutiae template of which
		 *             the quantization set is extracted.
		 *
		 * @return
		 *             The number of (distinct) quantization stored
		 *             in the output array <code>A</code>.
		 *
		 * @warning
		 *             If this object has not been initialized (either
		 *             explicitly via \link initialize()\endlink or
		 *             implicitly via \link enroll()\endlink), i.e.,
		 *             \link isInitialized()\endlink returns
		 *             <code>false</code>, an error message is printed
		 *             to <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If \link MinutiaeView::getFingerPosition()
		 *             view.getFingerPosition()\endlink
		 *             equals \link UNKNOWN_FINGER\endlink,
		 *             an error message is printed to <code>stderr</code> and
		 *             the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>A</code> is not able to store at
		 *             least \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             or at least \link MinutiaeView::getMinutiaeCount()
		 *             view.getMinutiaeCount()\endlink integers of type
		 *             <code>uint32_t</code>, the program runs into
		 *             undocumented behavior.
		 */
		int quantize( uint32_t *A , const MinutiaeView & view ) const;

		/**
		 * @brief
		 *             Computes the quantization set of a minutiae record
		 *             containing (absolutely pre-aligned) minutia template.
		 *
		 * @details
		 *             The function extracts at
		 *             most \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             (distinct) quantizations from each minutiae template
		 *             stored in <code>record</code> (preferring
		 *             those having the highest quality estimation;
		 *             see \link Minutia::getQuality()\endlink).
		 *             Consequently, the array should be able to contain
		 *             at least \link MinutiaeRecord::getViewCount()
		 *             record.getViewCount()\endlink
		 *             * \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             integers of type <code>uint32_t</code>.
		 *
		 * @param A
		 *             Array in which the quantization elements will be
		 *             stored.
		 *
		 * @param record
		 *             Record containing the absolutely pre-aligned minutiae
		 *             templates from which the quantization set is extracted.
		 *
		 * @return
		 *             The number of (distinct) quantization stored
		 *             in the output array <code>A</code>.
		 *
		 * @warning
		 *             If this object has not been initialized (either
		 *             explicitly via \link initialize()\endlink or
		 *             implicitly via \link enroll()\endlink), i.e.,
		 *             \link isInitialized()\endlink returns
		 *             <code>false</code>, an error message is printed
		 *             to <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If there is a minutiae template in <code>record</code>
		 *             specifying an unknown fingerposition, i.e.,
		 *             if \link MinutiaeView::getFingerPosition()
		 *             view.getFingerPosition()\endlink
		 *             equals \link UNKNOWN_FINGER\endlink,
		 *             an error message is printed to <code>stderr</code> and
		 *             the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>record</code> contains minutiae templates
		 *             specifying duplicate finger positions, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>A</code> is not able to store at
		 *             least \link MinutiaeRecord::getViewCount()
		 *             record.getViewCount()\endlink
		 *             * \link getMaxGenuineFeaturesPerFinger()\endlink
		 *             or at least \link MinutiaeRecord::getViewCount()
		 *             record.getViewCount()\endlink
		 *             * \link MinutiaeView::getMinutiaeCount()
		 *             view.getMinutiaeCount()\endlink integers of type
		 *             <code>uint32_t</code>, the program runs into
		 *             undocumented behavior.
		 */
		int quantize( uint32_t *A , const MinutiaeRecord & record ) const;

		/**
		 * @brief
		 *             Evaluate the reordering of a quantized feature which
		 *             is associated with an (enrolled) protected minutiae
		 *             template.
		 *
		 * @details
		 *             To prevent the application of extended Euclidean algorithm
		 *             based record multiplicity attack, each protected minutiae
		 *             record is attached with a public random permutation
		 *             process acting on the quantized feature universe.
		 *             The reordering of such a quantized feature (encoded as
		 *             a finite field element) can be evaluated via this
		 *             function.
		 *
		 *             For further details on the attack that is prevented via
		 *             the implemented reordering process we refer to
		 *             <ul>
		 *              <li>
		 *               <b>J. Merkle and B. Tams (2013)</b>.
		 *               Security of the Improved Fuzzy Vault Scheme in the
		 *               Presence of Record Multiplicity.
		 *               <i>CoRR
		 *                <a href="http://arxiv.org/abs/1312.5225"
		 *                 target="_blank">
		 *                 abs/1312.5225
		 *                </a>
		 *               </i>.
		 *              </li>
		 *             </ul>
	     *
		 * @param a
		 *             The feature element of which the reordering is
		 *             evaluated.
		 *
		 * @return
		 *             The reordering of <i>a</i>.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink returns
		 *             <code>false</code>, an error message may be printed to
		 *             <code>stderr</code> and the program exits with status
		 *             -1.
		 *
		 * @warning
		 *             If <i>a</i> does not encode a valid element in the
		 *             finite field that can be access via
		 *             \link getField()\endlink, an error message is printed
		 *             to <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 */
		uint32_t reorder( uint32_t a ) const;

		/**
		 * @brief
		 *             Changes the distance of the hexagonal grid used for
		 *             minutiae quantization.
		 *
		 * @details
		 *             Changing the distance of the hexagonal grid can cause
		 *             an update of the \link getField()
		 *             underlying finite field\endlink.
		 *
		 * @param gridDist
		 *             New distance of the hexagonal grid.
		 *
		 * @see getGridDist()
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink returns <code>true</code>,
		 *             the function prints an error message to
		 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>gridDist</code> is smaller than or equals 0,
		 *             the function prints an error message to
		 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void setGridDist( int gridDist );

		/**
		 * @brief
		 *             Change the number of possible quantizations for
		 *             fingerprint minutiae angles.
		 *
		 * @details
		 *             Changing the number of possible quantizations of the
		 *             involved minutiae angles can cause an update of
		 *             the \link getField() underlying finite field\endlink.
		 *
		 * @param s
		 *             The new number of possible quantizations of fingerprint
		 *             minutiae angles.
		 *
		 * @see getNumAngleQuanta()
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>true</code>, the
		 *             method prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <i>s</i> is smaller than or equals 0, the function
		 *             prints an error message to <code>stderr</code> and exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void setNumAngleQuanta( int s );

		/**
		 * @brief
		 *             Change the size of the secret polynomial being used
		 *             to obfuscate the minutiae template's quantization in
		 *             order to protect it.
		 *
		 * @param k
		 *             The new size of the secret polynomial being used to
		 *             obfuscate the minutiae template's quantization.
		 *
		 * @see getSecretSize()
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>true</code>, the
		 *             method prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <i>k</i> is smaller than or equals 0, the method
		 *             prints an error message to <code>stderr</code> and
		 *             exits with status 'EXIT_FAILURE'.
		 */
		void setSecretSize( int k );

		/**
		 * @brief
		 *             Changes the maximal number of minutiae quantizations
		 *             extracted each of the to-be-protected record's
		 *             absolutely pre-aligned minutiae templates during
		 *             enrollment.
		 *
		 * @param tmax
		 *             The new maximal number of minutiae quantizations
		 *             extracted from each record's absolutely pre-aligned
		 *             minutiae template.
		 *
		 * @see getMaxGenuineFeatures()
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>true</code>,
		 *             the method prints an error message to
		 *             <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>tmax</code> is smaller than or equals 0, the
		 *             method prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 */
		void setMaxGenuineFeaturesPerFinger( int tmax );

		/**
		 * @brief
		 *            Sets the maximal multiplicity used during a decoding
		 *            attempt with a Guruswami-Sudan algorithm.
		 *
		 * @details
		 *            The specified multiplicity parameter will be used
		 *            by the \link verify()\endlink and \link open()\endlink
		 *            functions that pass it as input to
		 *            the \link decode()\endlink function.
		 *
		 * @see getMultiplicity()
		 *
		 * @param m
		 *            The maximal multiplicity used in Guruswami-Sudan
		 *            algorithm based decoding attempt.
		 *
		 *
		 * @warning
		 *            If <code>m</code> is smaller than 1, an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		void setMultiplicity( int m );

		/**
		 * @brief
		 *            Specifies the slow-down factor to artificially
		 *            increase the time for an impostor verification
		 *            attempt.
		 *
		 * @details
		 *            The time needed for an impostor to perform
		 *            off-line attacks is multiplied by the slow-down
		 *            factor, thereby increasing security while, on the
		 *            other hand, also increasing response time of genuine
		 *            verification attempts---a trade-off thus but with
		 *            the capability of maintaining a relative security
		 *            measure with the advent of more and more powerful
		 *            computers.
		 *
		 * @param slowDownFactor
		 *            The specified slow-down factor for this protected
		 *            minutiae record.
		 *
		 * @see getSlowDownFactor()
		 *
		 * @warning
		 *            If this protected minutiae record object already
		 *            contains protected data, i.e.,
		 *            if \link isEnrolled()\endlink returns <code>true</code>,
		 *            an error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>slowDownFactor</code> encodes an integer
		 *            smaller than 1, an error message will be printed
		 *            to <code>stderr</code> and the program exits
		 *            with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be allocated, then an
		 *            error message will be printed to <code>stderr</code>
		 *            and the program exits with status 'EXIT_FAILURE'.
		 */
		void setSlowDownFactor( const BigInteger & slowDownFactor );

		/**
		 * @brief
		 *             Access the distance of the hexagonal grid used
		 *             for minutiae quantization.
		 *
		 * @return
		 *             The distance of the hexagonal grid used for minutiae
		 *             quantization.
		 *
		 * @see setGridDist(int)
		 */
		int getGridDist() const;

		/**
		 * @brief
		 *             Access the number of possible quantizations for
		 *             fingerprint minutiae angles.
		 *
		 * @return
		 *             The number of possible quantizations of fingerprint
		 *             minutiae angles.
		 *
		 * @see setNumAngleQuanta(int)
		 */
		int getNumAngleQuanta() const;

		/**
		 * @brief
		 *             Access the size of the secret polynomial being used
		 *             to obfuscate the minutiae template's quantization in
		 *             order to protect it.
		 *
		 * @return
		 *             The size of the secret polynomial being used to
		 *             obfuscate the minutiae template's quantization.
		 *
		 * @see setSecretSize(int)
		 */
		int getSecretSize() const;

		/**
		 * @brief
		 *             Access the maximal number of minutiae quantizations
		 *             extracted from the to-be-protected record's
		 *             minutiae templates on enrollment.
		 *
		 * @return
		 *             The maximal number of minutiae quantizations extracted
		 *             from the to-be-protected record's minutiae templates on
		 *             enrollment.
		 *
		 * @see setMaxGenuineFeaturesPerFinger(int)
		 */
		int getMaxGenuineFeaturesPerFinger() const;

		/**
		 * @brief
		 *             Accesses the maximal multiplicity bound used by this
		 *             protected minutiae record object's during decoding
		 *             with a Guruswami-Sudan algorithm.
		 *
		 * @return
		 *             The maximal multiplicity bound used by this
		 *             protected minutiae record object's during decoding
		 *             with a Guruswami-Sudan algorithm.
		 *
		 * @see setMultiplicity(int)
		 */
		int getMultiplicity() const;

		/**
		 * @brief
		 *             Access the resolution of the fingerprint images from
		 *             which the minutiae contained in the protected minutiae
		 *             record's minutiae template have been measured.
		 *
		 * @return
		 *             The resolution of the fingerprint images from which the
		 *             minutiae contained in the protected minutiae  record's
		 *             minutiae template have been measured.
		 */
		int getResolution() const;

		/**
		 * @brief
		 *            Returns the slow-down factor used to artifically
		 *            increase the time for a verification attempt on
		 *            impostor verification.
		 *
		 * @details
		 *            Unless manually specified via
		 *            the \link setSlowDownFactor()\endlink, the result
		 *            of this function is 1.
		 *
		 *            The time needed for an impostor to perform
		 *            off-line attacks is multiplied by the slow-down
		 *            factor, thereby increasing security while also
		 *            increasing response time of genuine verification
		 *            attempts---a trade-off thus but with the capability
		 *            of maintaining a relative security measure with
		 *            the advent of more and more powerful computers,
		 *
		 * @return
		 *            The slow-down factor of this protected minutiae
		 *            template.
		 *
		 * @see setSlowDownFactor()
		 */
		const BigInteger & getSlowDownFactor() const;

		/**
		 * @brief
		 *             Access the vector of hexagonal grid points centered in
		 *             a region in which absolutely pre-aligned minutiae
		 *             can occur.
		 *
		 * @return
		 *             The vector of hexagonal grid points covering the region
		 *             in which absolutely pre-aligned minutiae can occur.
		 */
		const std::vector< std::pair<double,double> > & getGrid() const;

		/**
		 * @brief
		 *             Returns the number of hexagonal grid points centered
		 *             in a region in which absolutel pre-algined minutiae can
		 *             occur.
		 *
		 * @return
		 *             (int)\link getGrid()\endlink .size()
		 */
		int getGridSize() const;

		/**
		 * @brief
		 *             Returns the size of the vault.
		 *
		 * @details
		 *             The size of the vault is the number of possible
		 *             minutiae quantizations.
		 *
		 * @return
		 *             \link getGridSize()\endlink
		 *             * \link getNumAngleQuanta()\endlink
		 */
		int getVaultSize() const;

		/**
		 * @brief
		 *             Returns the number in bytes that is needed to represent
		 *             the data of the vault polynomial stored by this
		 *             protected minutiae template.
		 *
		 * @details
		 *             This function returns the number of bytes being stored
		 *             in the array returned by \link getVaultData()\endlink
		 *             and is a multiple of 16.
		 *
		 * @return
		 *             The number in bytes of this object's vault data.
		 *
		 * @see getVaultData()
		 */
		int getVaultDataSize() const;

		/**
		 * @brief
		 *             Returns an array in which the data for the vault
		 *             polynomial of this object is stored.
		 *
		 * @details
		 *             Let
		 *             \f[
		 *              V(X)=X^t+\sum_{j=0}^{t-1} V_j\cdot X^j
		 *             \f]
		 *             be the vault polynomial that protects the feature set
		 *             \f$A\f$ (derived from the minutiae template protected
		 *             by this object). Write \f$b_{0},...,b_{m-1}\f$
		 *             be the bit sequence resulting from a concatenation of
		 *             the coefficients \f$V_0,...,V_{t-1}\f$. Let
		 *             \f$b_{m},...,b_{M-1}\f$ be random bits such that
		 *             \f$M=96\cdot\lceil m/96\rceil\f$, i.e.,
		 *             \f$ b_{0},...,b_{m-1}\f$ are padded with random bits
		 *             such that the total number of bits becomes a multiple
		 *             of 96. An \link AES128 AES key\endlink is then derived
		 *             from a random positive integer smaller
		 *             than \link getSlowDownFactor()\endlink and used to
		 *             "encrypt" the bits \f$b_0,...,b_{M-1}\f$. It is important
		 *             to note that, if the slow-down factor is equals 1, then
		 *             there is only one possibility to generate the AES key,
		 *             i.e., derived from the integer 0. In such we obtain the
		 *             "encrypted" bit sequence
		 *             \f[
		 *              (b'_0,...,b'_{M-1})
		 *             \f]
		 *             which is stored byte-wisely in the array returned by
		 *             this function.
		 *
		 *             The reason why the polynomial data needs to be padded
		 *             such that it can be passed through AES functions is the
		 *             following. If \link getSlowDownFactor()\endlink is
		 *             larger than 1, then the \link AES128 AES key\endlink
		 *             through which the padded vault polynomial bits are
		 *             passed serve as a slow-down mechanism: On a false
		 *             verification attempt, the impostor has to iterate
		 *             through all \link AES128 AES keys\endlink derived from
		 *             integers smaller than \link getSlowDownFactor()\endlink
		 *             before he (most likely) gets rejected thereby improving
		 *             security against false-accept attack while, on the other
		 *             hand, also increasing verification time on a genuine
		 *             verification attempt---a trade-off thus. However, as
		 *             computers become more and more efficient, the slow-down
		 *             factor may be increased and thus a relative security
		 *             measure can be maintained.
		 *
		 *             Also note that the vault data can be further
		 *             encrypted using a user-specific password, for which
		 *             purpose the method \link encrypt()\endlink
		 *             and \link decrypt()\endlink can be used for encryption
		 *             and decryption, respectively.
		 *
		 * @return
		 *             The array containing the vault polynomial's (possibly
		 *             encrypted) data stored by this object; the result is
		 *             <code>NULL</code> in case \link isEnrolled()\endlink
		 *             returns <code>false</code>.
		 *
		 * @see getVaultDataSize()
		 * @see setSlowDownFactor()
		 * @see getSlowDownFactor()
		 * @see encrypt()
		 * @see decrypt()
		 */
		const uint8_t *getVaultData() const;

		/**
		 * @brief
		 *            Derives the vault polynomial from the data stored by
		 *            this object.
		 *
		 * @details
		 *            If \link getSlowDownFactor()\endlink equals 1 (which, by
		 *            default, is the case), then this function can be called
		 *            without specifying an argument thereby returning the
		 *            correct vault polynomial.
		 *
		 *            Otherwise, the correct vault polynomial has been
		 *            encrypted using a random secret key derived from an
		 *            integer, called the <i>slow-down value</i> smaller
		 *            than \link getSlowDownFactor()\endlink and then the
		 *            correct vault polynomial will be returned if the correct
		 *            slow-down value has been passed to this function.
		 *
		 *            This function is used by the \link open()\endlink
		 *            and implicitly by the \link verify()\endlink
		 *            functions in which slow-down values starting from
		 *            zero is successively increased and in each iteration
		 *            a decoding attempt is performed until the vault
		 *            has been successfully be opened or all slow-down
		 *            values smaller than \link getSlowDownFactor()\endlink
		 *            have been tested.
		 *
		 * @param slowDownVal
		 *            A guess for the slow-down value.
		 *
		 * @return
		 *            A candidate for the correct vault polynomial (which is
		 *            correct if the guess for the slow-down value is correct).
		 *
		 * @warning
		 *            If this object does not contain (a candidate) for decrypted
		 *            vault polynomial data, i.e.,
		 *            if \link isDecrypted()\endlink returns <code>false</code>,
		 *            then an error message will be printed to
		 *            <code>stderr</code> and the program exits with status
		 *            'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If <code>slowDownVal</code> is negative or greater than
		 *            or equals \link getSlowDownFactor()\endlink, then an
		 *            error message will be printed to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be allocated, an error message
		 *            will be printed to <code>stderr</code> and the program
		 *            exits with status 'EXIT_FAILURE'.
		 */
		SmallBinaryFieldPolynomial unpackVaultPolynomial
		( const BigInteger & slowDownVal = BigInteger(0) ) const;

		/**
		 * @brief
		 *             Access the SHA-1 hash value of the secret key
		 *             generated on enrollment used to obfuscate the
		 *             correct minutiae record's quantization.
		 *
		 * @return
		 *             The SHA-1 hash value of the secret key used to
		 *             obfuscate the correct minutiae template's quantization.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink returns
		 *             <code>false</code>, an error message is printed to
		 *             <code>stderr</code> and the program exits with status
		 *             -1.
		 */
		const uint8_t *getHash() const;

		/**
		 * @brief
		 *             Access the finite field in which computations with
		 *             the improved fuzzy vault scheme are performed to
		 *             protect multiple minutiae templates.
		 *
		 * @returns
		 *             A reference to the \link SmallBinaryField
		 *             finite field\endlink in which computations with the
		 *             improved fuzzy vault scheme are performed with this
		 *             object.
		 */
		const SmallBinaryField & getField() const;

		/**
		 * @brief
		 *            Returns a boolean value that indicates whether this
		 *            object has been initialized such that it can protect an
		 *            absolutely pre-aligned minutiae template.
		 *
		 * @details
		 *            After an object of this class has been created via the
		 *            standard
		 *            constructor \link ProtectedMinutiaeRecord()\endlink
		 *            the width, height, and resolution to which a
		 *            to-be-protected minutiae record corresponds to is
		 *            unknown but need to be known on enrollment and can
		 *            be specified via the \link initialize()\endlink
		 *            method; the \link initialize()\endlink method is
		 *            called implicitly during a call of
		 *            the \link enroll()\endlink; however, if we want to
		 *            make use of the \link quantize()\endlink function
		 *            without generating actual protected minutiae record
		 *            data, we have to call \link initialize()\endlink
		 *            explicitly. This function returns whether a call of
		 *            the \link initialize()\endlink method has been made and
		 *            is used by the \link quantize()\endlink function to
		 *            produce error handling if the fmethod was not
		 *            called.
		 *
		 * @return
		 *            <code>true</code> if this object has been initialized
		 *            such that \link quantize()\endlink can be called safely;
		 *            otherwise, this function returns <code>false</code>.
		 */
		bool isInitialized() const;

		/**
		 * @brief
		 *            Access whether this object does hold
		 *            actual protected data of a minutiae record.
		 *
		 * @return
		 *            <code>true</code> if this object  does protect
		 *            a minutiae record; otherwise <code>false</code>.
		 *
		 * @see enroll()
		 */
		bool isEnrolled() const;

		/**
		 * @brief
		 *            Determine whether this instance does hold a
		 *            protected minutiae template that is additionally
		 *            protected by a (user-specific) key.
		 *
		 * @details
		 *            If this function returns <code>true</code>, the
		 *            \link ProtectedMinutiaeRecord\endlink has been
		 *            encrypted via \link encrypt()\endlink;
		 *            otherwise, the function returns <code>false</code>.
		 *            If this \link ProtectedMinutiaeRecord\endlink
		 *            has been decrypted via \link decrypt()\endlink
		 *            such that it represents a candidate
		 *            for a decrypted \link ProtectedMinutiaeRecord\endlink,
		 *            the function returns <code>false</code>.
		 *
		 * @return
		 *            <code>true</code> if this
		 *            \link ProtectedMinutiaeRecord\endlink does protect
		 *            a minutiae template being additionally protected
		 *            by an \link AES128\endlink key and does not represent
		 *            a candidate for a decrypted
		 *            \link ProtectedMinutiaeRecord\endlink; otherwise, the
		 *            function returns <code>false</code>.
		 *
		 * @see encrypt()
		 */
		bool isEncrypted() const;

		/**
		 * @brief
		 *            Determine whether this instance holds a decrypted
		 *            candidate for a protected minutiae template.
		 *
		 * @details
		 *            Any encrypted \link ProtectedMinutiaeRecord\endlink
		 *            can be decrypted via an \link AES128\endlink key and,
		 *            subsequently, contains a candidate for the decrypted
		 *            \link ProtectedMinutiaeRecord\endlink.
		 *
		 *            Note that, if this
		 *            \link ProtectedMinutiaeRecord\endlink
		 *            has not been decrypted via \link decrypt()\endlink
		 *            but contains enrolled data, the function returns
		 *            <code>true</code>.
		 *
		 * @return
		 *            <code>true</code> if this
		 *            \link ProtectedMinutiaeRecord\endlink contains a
		 *            candidate for the unencrypted
		 *            \link ProtectedMinutiaeRecord\endlink; otherwise,
		 *            the result will be <code>false</code>.
		 *
		 * @see decrypt()
		 */
		bool isDecrypted() const;

		/**
		 * @brief
		 *            Determine whether this instance does hold
		 *            encrypted data of a protected minutiae template.
		 *
		 * @details
		 *            Note that a \link ProtectedMinutiaeRecord\endlink
		 *            keeps the encrypted data even after being decrypted
		 *            via \link decrypt()\endlink
		 *
		 * @return
		 *            <code>true</code> if the
		 *            \link ProtectedMinutiaeRecord\endlink contains
		 *            encrypted data; otherwise, the function returns
		 *            <code>false</code>.
		 *
		 * @see encrypt()
		 */
		bool containsEncryptedData() const;

		/**
		 * @brief
		 *             Clears all data that is related with a minutiae
		 *             record protected by this instance (if any).
		 *
		 * @details
		 *             Any \link ProtectedMinutiaeRecord\endlink is empty
		 *             (cleared) by default until a
		 *             \link MinutiaeRecord minutiae record\endlink is passed
		 *             through the \link enroll()\endlink function returning
		 *             <code>true</code>. If the user wants to reuse this
		 *             \link ProtectedMinutiaeRecord\endlink objectto protect
		 *             another \link MinutiaeView minutiae template\endlink
		 *             via \link enroll()\endlink, he must clear it first via
		 *             this method.
		 */
		void clear();

		/**
		 * @brief
		 *             Determines the size in bytes needed to store this
		 *             protected minutiae record.
		 *
		 * @return
		 *             The size in bytes needed to store the data for this
		 *             protected minutiae template.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>false</code>, this
		 *             function prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 */
		int getSizeInBytes() const;

		/**
		 * @brief
		 *             Stores data defining this protected minutiae template to
		 *             the specified byte array.
		 *
		 * @param data
		 *             The byte array to which the data for this protected
		 *             minutiae record is written.
		 *
		 * @return
		 *             The function returns \link getSizeInBytes()\endlink.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>false</code>, this
		 *             function prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>data</code> cannot hold at least
		 *             \link getSizeInBytes()\endlink bytes, the function runs
		 *             into undocumented behavior.
		 */
		int toBytes( uint8_t *data ) const;

		/**
		 * @brief
		 *             Initialize this protected minutiae record by the
		 *             specified bytes.
		 *
		 * @details
		 *             If this protected minutiae records has been successfully
		 *             initialized by the specified bytes, the function returns
		 *             the number of relevant bytes; otherwise, if an error
		 *             occurred, this protected minutiae template will be left
		 *             unchanged and the result of the function will be -1.
		 *
		 * @param data
		 *             Contains <code>size</code> well-defined bytes of type
		 *             <code>uint8_t</code>.
		 *
		 * @param size
		 *             The number of well-defined bytes contained in
		 *             <code>data</code>.
		 *
		 * @return
		 *             The number of bytes (<=size) that were read from
		 *             <code>data</code> to initialize this protected minutiae
		 *             template; otherwise, if the first bytes of
		 *             <code>data</code> do not specify a well-defined protected
		 *             minutiae template, the function returns -1.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		int fromBytes( const uint8_t *data , int size );

		/**
		 * @brief
		 *             Write data defining this protected minutiae record to
		 *             the specified file pointer.
		 *
		 * @param out
		 *             Specifies the file to which this protected minutiae
		 *             template is written.
		 *
		 * @return
		 *             The function returns the number of bytes that have been
		 *             successfully written to <code>out</code>.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>false</code>, this
		 *             function prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		int write( FILE *out ) const;

		/**
		 * @brief
		 *             Writes data defining this protected minutiae record to
		 *             the specified file.
		 *
		 * @param file
		 *             Path to the file to which this protected minutiae
		 *             record is written.
		 *
		 * @return
		 *             The function returns <code>true</code> if this protected
		 *             minutiae record has been successfully written
		 *             to the specified file; otherwise, if an error occurred,
		 *             the function returns <code>false</code>.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink is <code>false</code>, this
		 *             function prints an error message to <code>stderr</code>
		 *             and exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool write( const std::string & file ) const;

		/**
		 * @brief
		 *             Initialize this protected minutiae record from the
		 *             specified file pointer.
		 *
		 * @details
		 *             If this protected minutiae record has been successfully
		 *             initialized, the function returns the number read bytes;
		 *             otherwise, if an error occurred, this protected minutiae
		 *             record will be left unchanged and the result of the
		 *             function will be -1.
		 *
		 * @param in
		 *             The file pointer from which the protected minutiae
		 *             record is read.
		 *
		 * @return
		 *             The number of bytes that were successfully read from
		 *             <code>in</code> to initialize this protected minutiae
		 *             record; otherwise, if an error occurred,
		 *             the function returns -1.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		int read( FILE *in );

		/**
		 * @brief
		 *             Initialize this protected minutiae record from the
		 *             data contained in the specified file path.
		 *
		 * @details
		 *             If this protected minutiae record has been successfully
		 *             initialized, the function returns <code>true</code>;
		 *             otherwise, if an error occurred, this protected minutiae
		 *             record will be left unchanged and the result of the
		 *             function will be <code>false</code>.
		 *
		 * @param file
		 *             String representation of the path to the file from
		 *             which the protected minutiae record is read.
		 *
		 * @return
		 *             <code>true</code> if this protected minutiae record
		 *             has been successfully initialized by the data contained
		 *             in the specified file; otherwise, the result is
		 *             <code>false</code>.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool read( const std::string & file );

		/**
		 * @brief
		 *             Swaps the content of two
		 *             \link ProtectedMinutiaeRecord\endlink.
		 *
		 * @details
		 *             After invocation, <code>vault1</code> holds the content of
		 *             the <code>vault2</code> at input and vice versa.
		 *
		 * @param vault1
		 *             A reference to a
		 *             \link ProtectedMinutiaeRecord\endlink
		 *
		 * @param vault2
		 *             A reference to a
		 *             \link ProtectedMinutiaeRecord\endlink
		 */
		static void swap
		( ProtectedMinutiaeRecord & vault1 , ProtectedMinutiaeRecord & vault2 );

	private:

		/**
		 * @brief
		 *            The width of the fingerprint images that this protect
		 *            minutiae record protects.
		 *
		 * @see initialize()
		 */
		int width;

		/**
		 * @brief
		 *            The height of the fingerprint images that this protect
		 *            minutiae record protects.
		 *
		 * @see initialize()
		 */
		int height;

		/**
		 * @brief
		 *            The resolution at which the fingerprint images protected
		 *            by this object have been scanned.
		 *
		 * @see initialize()
		 */
		int dpi;

		/**
		 * @brief
		 *            The distance of the hexagonal grid used for minutiae
		 *            quantization.
		 *
		 * @details
		 *            The value of this field by default depends on the
		 *            resolution specified via the \link initialize()\endlink
		 *            method and chosen 25 for 500 dots per inch.
		 *
		 * @see getGridDist()
		 * @see setGridDist(int)
		 */
		int gridDist;

		/**
		 * @brief
		 *            The number of possible quantizations for fingerprint
		 *            minutiae angles.
		 *
		 * @details
		 *            The value of this field is set to 6 by default.
		 *
		 * @see setNumAngleQuanta(int)
		 * @see getNumAngleQuanta()
		 */
		int s;

		/**
		 * @brief
		 *            The size of the secret polynomial being used to obfuscate
		 *            the minutiae records's quantization.
		 *
		 * @details
		 *            The value of this field is set to 25 by default.
		 *
		 * @see setSecretSize(int)
		 * @see getSecretSize()
		 */
		int k;

		/**
		 * @brief
		 *            The maximal number of minutiae quantizations extracted
		 *            from each record's to-be-protected minutiae template
		 *            on enrollment.
		 *
		 * @details
		 *            The value of this field is set to 44 by default.
		 *
		 * @see setMaxGenuineFeatures(int)
		 * @see getMaxGenuineFeatures()
		 */
		int tmax;

		/**
		 * @brief
		 *            The maximal multiplicity used on opening the
		 *            vault via the Guruswami-Sudan algorithm based
		 *            decoder with which this vault is opened at a
		 *            verification attempt.
		 *
		 * @see setMultiplicity(int)
		 * @see getMultiplicity()
		 */
		int m;

		/**
		 * @brief
		 *            The slow-down factor associated with this protected
		 *            minutiae record.
		 *
		 * @see getSlowDownFactor()
		 * @see setSlowDownFactor()
		 */
		BigInteger slowDownFactor;

		/**
		 * @brief
		 *            Vector of hexagonal grid points used for
		 *            minutiae quantization.
		 *
		 * @details
		 *            The state of this field is controlled by the object and
		 *            can be controlled implicitly by specifying the
		 *            grid distance (e.g., via \link setGridDist()\endlink)
		 *            and the fingerprint image dimensions
		 *            (e.g., via \link initialize()\endlink).
		 *
		 *            A constant reference to this field can be accessed
		 *            via \link getGrid()\endlink.
		 *
		 * @see getGrid()
		 */
		std::vector< std::pair<double,double> > grid;

		/**
		 * @brief
		 *            Pointer to the finite field used for template
		 *            protection via the improved fuzzy vault scheme.
		 *
		 * @details
		 *            On constructing an instance this field should first be
		 *            set to <code>NULL</code> or a valid pointer to a
		 *            \link SmallBinaryField finite field\endlink and, after the
		 *            member variables
		 *            \link width\endlink, \link height\endlink,
		 *            \link gridDist\endlink, and \link s\endlink have been
		 *            specified a proper finite field is updated via
		 *            \link updateField()\endlink.
		 *
		 * @see updateField()
		 * @see getField()
		 */
		SmallBinaryField *gfPtr;

		/**
		 * @brief
		 *            The size of the feature set protected by this instance.
		 *
		 * @details
		 *            On \link enroll() enrollment\endlink a
		 *            \link MinutiaeRecord minutiae record\endlink is used to
		 *            extract a feature set via \link quantize()\endlink
		 *            which is protected using the improved fuzzy vault scheme.
		 *            This field is actually deemed to be equals the size of the
		 *            protected feature set. However, in this implementation,
		 *            the feature set is supplemented by a certain number of
		 *            blending feature elements such that it reaches a size
		 *            of 10 *\link tmax\endlink. Consequently, it is
		 *            impossible to derive information about the (say) quality
		 *            of the protected minutiae record just from observing
		 *            <i>t</i>.
		 */
		int t;

		/**
		 * @brief
		 *            Array containing the (decrypted) vault data.
		 *
		 * @see getVaultData()
		 * @see unpackVaultPolynomial()
		 */
		uint8_t *vaultPolynomialData;

		/**
		 * @brief
		 *            Encrypted vault data.
		 *
		 * @details
		 *            After \link enroll() enrollment\endlink the protected
		 *            minutiae template can additionally be encrypted via
		 *            a (say) password or PIN to improve security. If the
		 *            vault is additionally protected by a
		 *            \link AES128 key\endlink, this member will be
		 *            non-<code>NULL</code> and <code>NULL</code> otherwise.
		 *
		 * @see encrypt()
		 * @see decrypt()
		 * @see getVaultData()
		 */
		uint8_t *encryptedVaultPolynomialData;

		/**
		 * @brief
		 *            SHA-1 hash value of the secret polynomial protecting
		 *            the enrolled minutiae record's quantization.
		 *
		 * @details
		 *            On \link enroll() enrollment\endlink,
		 *            a \link MinutiaeRecord minutiae record\endlink is bound
		 *            to a \link SmallBinaryFieldPolynomial
		 *            secret polynomial\endlink of which any data is dismissed
		 *            except of its SHA-1 hash value which is stored in this
		 *            field.
		 *
		 * @see getHash()
		 */
		uint8_t hash[20];

		/**
		 * @brief
		 *            A random public permutation process used to shuffle
		 *            protected feature sets in order to avoid
		 *            record multiplicity attacks.
		 *
		 * @details
		 *            Without preventions, an adversay may determine related
		 *            from non-related protected records via the improved
		 *            fuzzy vault scheme and, even worse, he may break them
		 *            entirely with quite a high probability using an attack
		 *            based on the Euclidean algorithm (see
		 *            \link EEAAttack\endlink). To prevent successful
		 *            application of the attack, each protected minutiae
		 *            record is attached with a random public
		 *            record-specific permutation process acting on the
		 *            feature sets which is represented by this field.
		 *
		 *            This field is updated on each enrollment via
		 *            \link updatePermutation()\endlink. The permutation
		 *            can be accessed via the \link reorder()\endlink
		 *            function.
		 *
		 * @see reorder()
		 * @see _reorder()
		 * @see updatePermutation()
		 */
		Permutation permutation;

		/**
		 * @brief
		 *            Initializes this object's members such that it
		 *            represents an empty uninitialized protected minutiae
		 *            record.
		 *
		 * @details
		 *            This method is called by the standard and copy
		 *            constructor.
		 */
		void first_init();

		/**
		 * @brief
		 *            Update the vector of hexagonal grid points after
		 *            corresponding parameters have changed or specified
		 *            the very first time.
		 *
		 * @details
		 *            If the member variables \link width\endlink,
		 *            \link height\endlink, or \link gridDist\endlink
		 *            have changed, the
		 *            \link grid vector of hexagonal grid point\endlink
		 *            needs an update which is performed by calling
		 *            this method.
		 */
		void updateGrid();

		/**
		 * @brief
		 *            Updates the finite field such that it can encode
		 *            each possible minutia's quantization.
		 *
		 * @details
		 *            After the \link grid hexagonal grid points\endlink
		 *            or the
		 *            \link s number of minutiae angle quantizations\endlink
		 *            have been specified, the finite field might need an
		 *            update which can be performed by calling this method.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void updateField();

		/**
		 * @brief
		 *            Generate a public random user-specific permutation
		 *            process acting on the feature universe.
		 *
		 * @details
		 *            On enrollment, the feature set extracted from a
		 *            minutiae template have to be passed through a public
		 *            random user-specific permutation process to prevent
		 *            certain record multiplicity attacks
		 *            (see \link EEAAttack\endlink). The method generates
		 *            this permutation process by updating the
		 *            \link permutation\endlink field using a
		 *            \link RandomGenerator random generator\endlink whose
		 *            seed is derived from the hash value of the secret
		 *            polynomial (i.e., the \link hash\endlink field).
		 *            In such the permutation is randomly generated as the
		 *            secret polynomial, and thus, \link hash\endlink has
		 *            been randomly generated and is public such that on
		 *            verification the query feature sets can be passed
		 *            through the permutation process, too.
		 *
		 * @see permutation
		 */
		void updatePermutation();

		/**
		 * @brief
		 *             Evaluate the reordering of a quantized feature which
		 *             is associated with an (enrolled) protected minutiae
		 *             template.
		 *
		 * @details
		 *             This function is not public and intended for internal
		 *             purposes. It is a variant of the
		 *             \link reorder() reorder(uint32_t)\endlink function but
		 *             does not ensure that \link isEnrolled()\endlink is
		 *             <code>true</code> thereby assuming that the
		 *             \link permutation\endlink field is initialized
		 *             properly. Otherwise, the function runs into
		 *             undocumented behavior.
		 *
		 * @param a
		 *             The feature element of which the reordering is
		 *             evaluated.
		 *
		 * @return
		 *             The reordering of <i>a</i>.
		 *
		 * @warning
		 *             If <i>a</i> does not encode a valid element in the
		 *             finite field that can be access via
		 *             \link getField()\endlink, an error message is printed
		 *             to <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 */
		uint32_t _reorder( uint32_t a ) const;

		/**
		 * @brief
		 *             Passes the specified polynomial to the
		 *             array \link vaultPolynomialData\endlink.
		 *
		 * @details
		 *             Let
		 *             \f[
		 *              V(X)=X^t+\sum_{j=0}^{t-1} V_j\cdot X^j
		 *             \f]
		 *             be the vault polynomial passed to this function.
		 *             Write \f$b_{0},...,b_{m-1}\f$
		 *             be the bit sequence resulting from a concatenation of
		 *             the coefficients \f$V_0,...,V_{t-1}\f$. Let
		 *             \f$b_{m},...,b_{M-1}\f$ be random bits such that
		 *             \f$M=96\cdot\lceil m/96\rceil\f$, i.e.,
		 *             \f$ b_{0},...,b_{m-1}\f$ are padded with random bits
		 *             such that the total number of bits becomes a multiple
		 *             of 96. An \link AES128 AES key\endlink is then derived
		 *             from a random positive integer smaller
		 *             than \link getSlowDownFactor()\endlink and used to
		 *             "encrypt" the bits \f$b_0,...,b_{M-1}\f$. In such we
		 *             obtain the "encrypted" bit sequence
		 *             \f[
		 *              (b'_0,...,b'_{M-1})
		 *             \f]
		 *             which is stored byte-wisely in
		 *             \link vaultPolynomialData\endlink
		 *
		 * @param V
		 *             The polynomial packed
		 *             into \link vaultPolynomialData\endlink.
		 *
		 *
		 * @warning
		 *             The field \link vaultPolynomialData\endlink will be
		 *             overwritten by newly allocated memory; thus, before
		 *             calling this method it must be ensured that
		 *             possibly allocated memory to
		 *             which \link vaultPolynomialData\endlink refers is freed;
		 *             otherwise, the program can cause memory leaks during
		 *             its runtime.
		 *
		 * @warning
		 *             The behavior of this method is undocumented
		 *             if \link t\endlink and \link gfPtr\endlink are not
		 *             initialized properly.
		 *
		 * @warning
		 *             If not enough memory can be provided, an error message
		 *             is printed <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'
		 *
		 * @warning
		 *             The polynomial represented by <code>V</code> should
		 *             be defined over the field \link getField()\endlink;
		 *             otherwise, the method runs into undocumented
		 *             behavior.
		 */
		void packVaultPolynomial( const SmallBinaryFieldPolynomial & V );

		/**
		 * @brief Determines how many bytes are needed to hold the data of
		 * the packed vault polynomial.
		 *
		 * @details The function is used to determine how many bytes are
		 * needed to initialize an array of 128 bit blocks that can hold the
		 * data of this fuzzy vault object's vault polynomial being of degree
		 * smaller than \link t\endlink such that it can be protected with
		 * an \link AES128 AES key\endlink.
		 *
		 * @return Number of bytes needed to hold the data of a monic
		 * polynomial of degree smaller than \link t\endlink having
		 * coefficients in the field \link gfPtr\endlink such that the byte
		 * array is a multiple of 16 bytes (i.e., 128 bits).
		 */
		int vaultDataSize() const;

		/**
		 * @brief
		 *            Derives an AES key from the specified integer.
		 *
		 * @details
		 *            The bytes of the specified integer are extracted
		 *            via the \link BigInteger::toBytes() x.toBytes()\endlink
		 *            method and then its \link SHA SHA hash value\endlink
		 *            is computed; the first 16 bytes of the value
		 *            are used to define the \link AES128 AES key\endlink
		 *            returned by this function.
		 *
		 * @param x
		 *            The integer of which an AES key is computed with
		 *            this function.
		 *
		 * @return
		 *            The AES key derived from <code>x</code>.
		 *
		 * @warning
		 *             If not enough memory can be provided, an error message
		 *             is printed <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'
		 */
		static AES128 deriveKey( const BigInteger & x );

		/**
		 * @brief
		 *           Store the concatenation of a sequence of <i>d</i>-bit
		 *           vectors in a byte array and supplement it with random
		 *           bits.
		 *
		 * @details
		 *           The method is used to convert an instance of the improved
		 *           fuzzy vault scheme (encoded by <i>t</i>
		 *           elements of a binary finite field of degree <i>d</i>)
		 *           to a byte array that can be encrypted via an
		 *           \link AES128 AES key\endlink. The random bits are needed
		 *           to make the unencrypted array indistinguishable from falsely
		 *           decrypted data.
		 *
		 * @param bits
		 *           The concatenation of the bit vectors.
		 *
		 * @param vecs
		 *           Array of <i>t</i> 32-bit word each encoding a <i>d</i>-bit
		 *           vector.
		 *
		 * @param t
		 *           Number of 32-bit words contained in <code>array</code>.
		 *
		 * @param d
		 *           Number of significant bits of each 32-bit vector.
		 *
		 * @param n
		 *           Controls how many bits in <code>bits</code> are
		 *           supplemented by random values. More precisely, the latest
		 *           <i>8*n-t*d</i> bits in <code>bits</code> are selected
		 *           randomly.
		 */
		static void concat_bit_vectors
		( uint8_t *bits , const uint32_t *vecs , int t , int d , int n );

		/**
		 * @brief
		 *            Deconcatenate a bit sequence into <i>d</i>-bit strings
		 *            stored in 32 bit words.
		 *
		 * @details
		 *            The method is used to re-obtain the candidate of
		 *            the improved fuzzy vault instance data from a byte
		 *            array or after decryption.
		 *
		 * @param vecs
		 *            Output vector that will contain <i>t</i> 32-bit words
		 *            of <i>d</i> significant bits.
		 *
		 * @param bits
		 *            Byte array encoding the <i>t*d</i> bit sequence that
		 *            is deconcatenated into <code>vecs</code>.
		 *
		 * @param t
		 *            Number of <i>d</i>-bit strings stored in the bit sequence
		 *            <code>bits</code>.
		 *
		 * @param d
		 *            Positive integer smaller than or equals 32 controlling
		 *            the length of the bit strings encoded in the sequence
		 *            <code>bits</code>.
		 */
		static void split_into_bit_vectors
		( uint32_t *vecs , const uint8_t *bits , int t , int d  );
	};
}


#endif /* THIMBLE_PROTECTEDMINUTIAERECORD_H_ */
