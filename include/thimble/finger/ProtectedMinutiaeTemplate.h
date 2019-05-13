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
 * @file ProtectedMinutiaeTemplate.h
 *
 * @brief
 *            Provides a class for representing protected minutiae templates
 *            of a single finger that are absolutely pre-algined w.r.t. a
 *            directed reference point.
 *
 * @details
 * @section sec_ffv Minutiae Template Protection
 *
 * THIMBLE provides the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class for generating protected data from
 * a minutiae template that has been pre-aligned w.r.t. a
 * fingerprint's coordinate system estimation. From the protected data it
 * should be hard to reveal information about the protected fingerprint
 * while, given a second fingerprint that matches with the protected
 * fingerprint, it is still possible to verify that they stem from the
 * same finger, similar as a cryptographic password hash.
 * The so-called <em>fuzzy vault scheme</em> is a possible
 * solution that is currently considered for providing this kind of
 * <em>biometric template protection</em> for the fingerprint modality.
 * The \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink implements a fuzzy vault scheme
 * for "hashing" fingerprint minutiae templates following the ideas
 * presented in
 * <ul>
 *  <li><b>B. Tams (2013)</b>.
 *  Absolute %Fingerprint Pre-Alignment in Minutiae-Based
 *  Cryptosystems. <i>Proc. BIOSIG 2013, ser LNI,
 *  <b>vol. 212</b>, pp. 75--86</i>.</li>
 *  <li><b>B. Tams, P. Mihailescu, and A. Munk (2015)</b>.
 *  Security Considerations in Minutiae-Based Fuzzy Vault,
 *  <i>IEEE Trans. Inf. Forensics Security</i>, vol. 10, no. 5, May 2015
 *  (<a href="http://www.stochastik.math.uni-goettingen.de/preprints/ffv.pdf">
 *    Preprint
 *   </a>).</li>
 * </ul>
 *
 * @subsection sec_ffv_multi Minutiae Template Protection for Multiple Fingerprints
 * THIMBLE also provides a mechanism for protecting minutiae templates
 * estimated from multiple fingerprints using objects generated from
 * the \link thimble::ProtectedMinutiaeRecord ProtectedMinutiaeRecord\endlink
 * class. This implementation follows the description of
 * <ul>
 *  <li>
 *   <b>B. Tams (2015)</b>. Unlinkable Minutiae-Based Fuzzy Vault for
 *   Multiple Fingerprints, <i>IET Biometrics</i>, to appear
 *   (<a href="http://www.stochastik.math.uni-goettingen.de/preprints/ffv.pdf">
 *    Preprint</a>).
 *  </li>
 * </ul>
 * Its interface is very similar to the interface of
 * the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class and we therefore recommend
 * to get familiar with using the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class first. For a tutorial on
 * the difference in using the \link thimble::ProtectedMinutiaeRecord
 * ProtectedMinutiaeRecord\endlink class see the section on
 * @ref sec_mffv.
 *
 * In the following, we describe how to
 * generate objects of the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class that hold protected minutiae
 * data. In the context of the class, fingerprint minutiae are assumed
 * to be pre-aligned and we therefore devote it a brief section first.
 *
 * @subsection sec_ffv_prealign Absolutely Pre-Aligned Minutiae
 *
 * Assume that we are given
 * a minutiae template of a fingerprint
 * <pre>
 *      MinutiaeView view;
 * </pre>
 * encoded as a \link thimble::MinutiaeView MinutiaeView\endlink object;
 * for more details on how minutiae templates are processed in THIMBLE,
 * we refer to the tutorial @ref sec_minutiae. The coordinate system may
 * be encoded by the coordinate of its origin
 * <pre>
 *      double x , y;
 * </pre>
 * and the angle of the abscissa (say) axis
 * <pre>
 *      double abscissaAngle;
 * </pre>
 * Then, the pre-aligned minutiae can be computed via
 * <pre>
 *      MinutiaeView prealignedView = FingerTools::prealign(view,x,y,abscissaAngle);
 * </pre>
 * For a mechanism to estimate a fingerprint's coordinate system
 * (or equivalently a fingerprint's reference point with a direction), we
 * refer to the tutorials @ref sec_drp and @ref sec_tarp.
 *
 * @subsection sec_ffv_enrol Enrolment
 *
 * Assume that we are given an absolutely pre-aligned minutiae template
 * <pre>
 *      MinutiaeView reference;
 * </pre>
 * that has been measured from a fingerprint image of width
 * <pre>
 *      int width;
 * </pre>
 * and height
 * <pre>
 *      int height
 * </pre>
 * scanned at a resolution of
 * <pre>
 *      int dpi;
 * </pre>
 * dots per inch. We may generate protected data from the reference
 * template via
 * <pre>
 *      ProtectedMinutiaeTemplate vault(width,height,dpi);
 *
 *      vault.enroll(reference);
 * </pre>
 * using the \link thimble::ProtectedMinutiaeTemplate::enroll()
 * enroll()\endlink function of objects generated by
 * the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class. Note
 * that \link thimble::ProtectedMinutiaeTemplate::enroll()
 * enroll()\endlink returns <code>true</code> if protected data
 * could be successfully generated and, otherwise, it returns
 * <code>false</code>.
 *
 * After successful enrollment, the protected data can be stored in a file
 * path
 * <pre>
 *      char *path;
 * </pre>
 * using
 * the \link thimble::ProtectedMinutiaeTemplate::write(const std::string&)const write()\endlink
 * function, i.e.,
 * <pre>
 *      vault.write(path);
 * </pre>
 * which returns <code>true</code> if the data has been successfully written
 * to the file path and, otherwise, <code>false</code>.
 *
 * To read a protected minutiae template from a file path, we may call
 * the \link thimble::ProtectedMinutiaeTemplate::read(const std::string&) read()\endlink
 * function which, again, returns <code>true</code> if the data has been
 * read from the path and otherwise <code>false</code>. For example, to
 * read a protected minutiae template we may run
 * <pre>
 *      ProtectedMinutiaeTemplate vault;
 *      vault.read(path);
 * </pre>
 *
 * @subsection sec_ffv_verify Verification
 *
 * Given a protected minutiae template
 * <pre>
 *      ProtectedMinutiaeTemplate vault;
 * </pre>
 * at which an absolutely pre-aligned reference template has been successfully
 * enroled and an absolutey pre-aligned query template
 * <pre>
 *      MinutiaeView query;
 * </pre>
 * we may use it test whether it is of sufficient agreement with the
 * protected reference template using
 * the \link thimble::ProtectedMinutiaeTemplate::verify() verify()\endlink
 * function, i.e., via
 * <pre>
 *      bool verified = vault.verify(query);
 * </pre>
 * The boolean <code>verified</code> will be equals <code>true</code>
 * if the minutiae contained in <code>query</code> are of sufficient
 * agreement with the minutiae protected by the <code>vault</code>.
 * Otherwise, if <code>verified</code> is <code>false</code> there exist
 * only a few (at best) minutiae in <code>query</code> agreeing
 * with the protected minutiae deemed to be a negative verification
 * result.
 *
 * @subsection sec_ffv_password Improving Security via Password
 *
 * A successfully enrolled protected minutiae template can additionally
 * be protected using a password/PIN. Encrypting minutiae template
 * can be performed using
 * the \link thimble::ProtectedMinutiaeTemplate::encrypt() encrypt()\endlink
 * method. For example, using the PIN "1234" we may run
 * <pre>
 *    vault.encrypt("1234");
 * </pre>
 * On verification, the same PIN must be provided and the vault may be
 * decrypted via
 * <pre>
 *    vault.decrypt("1234");
 * </pre>
 * and then we can perform a verification attempt
 * <pre>
 *    bool verified = vault.verify(query);
 * </pre>
 * If the false password is used to decrypt the vault, e.g.,
 * <pre>
 *    vault.decrypt("4321");
 * </pre>
 * then a verification attempt attempt
 * <pre>
 *    bool verified = vault.verify(query);
 * </pre>
 * will most likely return <code>false</code>.
 *
 * Internally, the password encryption/decryption mechanism is realized
 * via an implementation of the \link thimble::AES128 avanced encryption
 * algorithm\endlink.
 *
 * It is important to note that the encryption and decryption interface
 * provided by objects of the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class is such that it is not possible
 * to distinguish falsely from correctly decrypted templates: Bit blocks
 * containing the to-be-encrypted data of which only the first encode the
 * data are padded with random bits; on verification these bits are
 * dismissed. In such, the additional security provided by the strength of
 * the passwords/PINs multiplies with the potentially low security provided
 * by a protected minutiae template.
 *
 * @subsection sec_ffv_slowdown Maintaining Security via Slow-Down Functions
 *
 * There is also a slow-down mechanism to artificially increase verification
 * time on impostor verification attempts thereby increasing the difficulty
 * of an intruder to run off-line attack.
 *
 * To artificially increase the time for a verification we may specify a
 * slow-down factor with
 * the \link thimble::ProtectedMinutiaeTemplate::setSlowDownFactor()
 * setSlowDownFactor()\endlink method before enrollment. For example, if
 * we want to artificially slow down the verification process by the
 * factor 8, we may run
 * <pre>
 *    vault.setSlowDownFactor(8);
 *
 *    vault.enroll(reference);
 * </pre>
 * Then the time needed for running
 * <pre>
 *    bool verified = vault.verify(query);
 * </pre>
 * is increased by the factor 8 on a reject.
 *
 * Internally, the slow-down mechanism implemented in
 * the \link thimble::ProtectedMinutiaeTemplate
 * ProtectedMinutiaeTemplate\endlink class exploits a mechanism for encrypting
 * the protected templates similar as in the section described above
 * @ref sec_ffv_password. Therein an \link thimble::AES128 AES key\endlink
 * is derived from a randomly chosen positive integer smaller than the
 * specified slow-down factor called slow-down value; the key is then used
 * to decrypt the vault. On verification, an impostor has no chance to iterate
 * through all possible unknown slow-down values of which at most one will
 * could possibly yield to an accept given the query template is of reasonable
 * agreement with the protected reference template.
 *
 * It is important to note that the time needed for genuine verification
 * is increased by the slow-down mechanism as well. Therefore, the
 * slow-down factor must be chosen carefully. However, the slow-down mechanism
 * has the potential to maintain absolute vault security over time with which
 * computer that an attacker can use to break the vaults become more and more
 * powerful.
 *
 * @section sec_thimble_fjfx Fully Automatic Minutiae-Based Fuzzy Vault with THIMBLE and FingerJetFX OSE
 *
 * To implement a fully automatic we may first implement a subclass
 * <code>FJFXFingerprint</code> of the \link thimble::Fingerprint\endlink
 * class in which minutiae are estimated with Digital Persona's FingerJetFX
 * OSE. Therefore, see the tutorial @ref sec_fjfx_fingerprint.
 *
 * @subsection sec_thimble_fjfx_enrollment Enrollment
 *
 * On enrollment, we may read the to-be-protected fingerprint by running
 * <pre>
 *    FJFXFingerprint reference;
 *    bool success;
 *
 *    // Read fingerprint image from 'fingerprint.pgm'
 *    success = reference.fromImageFile("fingerprint.pgm");
 *    if ( !success ) {
 *        cerr << "ERROR: Could not read 'fingerprint.pgm'" << endl;
 *        exit(EXIT_FAILURE);
 *    }
 *
 *    // Ensure that a directed reference point can be estimated from
 *    // the reference fingerprint.
 *    if ( !reference.hasDirectedReferencePoint() ) {
 *        cerr << "FAILURE TO ALIGN" << endl;
 *        exit(EXIT_FAILURE);
 *    }
 *
 *    // Absolute Fingerprint Pre-Alignment
 *    MinutiaeView prealignedView = FingerTools::prealign
 *       (reference.getMinutiaeView(),reference.getDirectedReferencePoint());
 *
 *    // Initialize a vault object
 *    ProtectedMinutiaeTemplate vault(reference.getWidth(),reference.getHeight());
 *
 *    // Enrollment
 *    success = vault.enroll(prealignedView);
 *
 *    // If enrollment was successful, ...
 *    if ( success ) {
 *
 *        // ..., attempt to write the vault.
 *        success = vault.write("vault.dat");
 *        if ( !success ) {
 *            cerr << "ERROR: could not write 'vault.dat'" << endl;
 *            exit(EXIT_FAILURE);
 *        }
 *
 *    } else {
 *
 *        cerr << "FAILURE TO ENROL" << endl;
 *        exit(EXIT_FAILURE);
 *
 *    }
 * </pre>
 *
 * @subsection sec_thimble_fjfx_verification Verification
 *
 * We may implement the verification process by reading a query fingerprint
 * from 'query.pgm', say, and and verify its authenticity versus the protected
 * minutiae template read from 'vault.dat' as follows:
 * <pre>
 *    FJFXFingerprint query;
 *    bool success;
 *
 *    // Read query fingerprint
 *    success = query.fromImageFile("query.pgm");
 *    if ( !success ) {
 *         cerr << "ERROR: Could not read 'query.pgm'" << endl;
 *         exit(EXIT_FAILURE);
 *    }
 *
 *    // Ensure that a directed reference point can be estimated from
 *    // the query fingerprint.
 *    if ( !query.hasDirectedReferencePoint() ) {
 *        cerr << "FAILURE TO ALIGN" << endl;
 *        exit(EXIT_FAILURE);
 *    }
 *
 *    // Absolute Fingerprint Pre-Alignment
 *    MinutiaeView prealignedView = FingerTools::prealign
 *       (query.getMinutiaeView(),query.getDirectedReferencePoint());
 *
 *    // Read protected minutiae template from 'vault.dat'
 *    ProtectedMinutiaeTemplate vault;
 *    success = vault.read("vault.dat");
 *
 *    // Verification attempt
 *    bool verified = vault.verify(prealignedView);
 *
 *    if ( verified ) {
 *        cout << "ACCEPT" << endl;
 *    } else {
 *        cout << "REJECT" << endl;
 *    }
 *
 * </pre>
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_PROTECTEDMINUTIAETEMPLATE_H_
#define THIMBLE_PROTECTEDMINUTIAETEMPLATE_H_

#include <stdint.h>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/security/AES.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/math/Permutation.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/math/numbertheory/BigInteger.h>

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
	 *            protected from unprotected minutiae templates with
	 *            the <em>improved fuzzy vault scheme</em>.
	 *
	 *
	 * @see <b>B. Tams, P. Mihailescu, A. Munk (2014).</b>
	 *      <i>Security-Improved Minutiae-Based Fuzzy Vault</i>.
	 *      submitted.
	 */
	class THIMBLE_DLL ProtectedMinutiaeTemplate {

	public:

		/**
		 * @brief
		 *            Standard constructor.
		 *
		 * @details
		 *            Creates an empty vault in which no member fields
		 *            are set. The created object will not be useful
		 *            for protect an absolutely pre-aligned minutiae
		 *            before initialization, e.g., via
		 *            <pre>
		 *             ProtectedMinutiaeTemplate vault;
		 *             vault = ProtectedMinutiaeTemplate(width,height);
		 *            </pre>
		 *            in which the width and the height of the fingerprint
		 *            image is specified.
		 */
		ProtectedMinutiaeTemplate();

		/**
		 * @brief
		 *             Constructor.
		 *
		 * @details
		 *             The constructor generates parameters needed to properly
		 *             quantize absolutely pre-aligned minutiae templates that
		 *             have been estimated from a fingerprint image of the
		 *             specified dimension and that has been scanned at the
		 *             specified resolution. If the finger's position is known,
		 *             it can be specified as well.
		 *
		 * @param width
		 *             The width of the fingerprint image.
		 *
		 * @param height
		 *             The height of the fingerprint image.
		 *
		 * @param dpi
		 *             The resolution in (dots per inch) at which
		 *             the fingerprint image has been scanned.
		 *
		 * @param fingerPosition
		 *             The position of the finger (if known).
		 *
		 * @attention
		 *             If <code>width</code> or <code>height</code> are
		 *             not greater than 0 or if <code>dpi</code> does
		 *             not range between 300 and 1000, an error message
		 *             is printed to <code>stderr</code> and the program
		 *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		ProtectedMinutiaeTemplate
		( int width , int height ,
		  int dpi = 500 , FINGER_POSITION_T fingerPosition = UNKNOWN_FINGER );

		/**
		 * @brief
		 *             Copy constructor
		 *
		 * @param vault
		 *             The protected minutiae template of which
		 *             a copy is constructed.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		ProtectedMinutiaeTemplate
		( const ProtectedMinutiaeTemplate & vault );

		/**
		 * @brief
		 *             Destructor
		 *
		 * @details
		 *             Frees all the memory helt by the protected minutiae
		 *             template.
		 */
		~ProtectedMinutiaeTemplate();

		/**
		 * @brief
		 *             Initializes this object such that it can protected
		 *             absolutely pre-aligned minutiae templates.
		 *
		 * @details
		 *             This method generates parameters needed to properly
		 *             quantize absolutely pre-aligned minutiae templates that
		 *             have been estimated from a fingerprint image of the
		 *             specified dimension and that has been scanned at the
		 *             specified resolution. If the finger's position is known,
		 *             it can be specified as well.
		 *
		 *             Any data hold by this object will be deleted and
		 *             replaced by new data.
		 *
		 * @param width
		 *             The width of the fingerprint image.
		 *
		 * @param height
		 *             The height of the fingerprint image.
		 *
		 * @param dpi
		 *             The resolution in (dots per inch) at which
		 *             the fingerprint image has been scanned.
		 *
		 * @param fingerPosition
		 *             The position of the finger (if known).
		 *
		 * @attention
		 *             If <code>width</code> or <code>height</code> are
		 *             not greater than 0 or if <code>dpi</code> does
		 *             not range between 300 and 1000, an error message
		 *             is printed to <code>stderr</code> and the program
		 *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void initialize
		( int width , int height ,
		  int dpi = 500 , FINGER_POSITION_T fingerPosition = UNKNOWN_FINGER );

		/**
		 * @brief
		 *             Assignment operator
		 *
		 * @details
		 *             Sets the \link ProtectedMinutiaeTemplate\endlink to a
		 *             copy of <code>vault</code>.
		 *
		 * @param vault
		 *             A \link ProtectedMinutiaeTemplate\endlink of which this
		 *             instance becomes a copy of.
		 *
		 * @return
		 *             A reference to this
		 *             \link ProtectedMinutiaeTemplate\endlink
		 *             (after the assignment side-effect).
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		ProtectedMinutiaeTemplate &operator=
				( const ProtectedMinutiaeTemplate & vault );

		/**
		 * @brief
		 *             Swaps the content of two
		 *             \link ProtectedMinutiaeTemplate\endlink.
		 *
		 * @details
		 *             After invocation, <code>vault1</code> holds the content of
		 *             the <code>vault2</code> at input and vice versa.
		 *
		 * @param vault1
		 *             A reference to a valid
		 *             \link ProtectedMinutiaeTemplate\endlink
		 *
		 * @param vault2
		 *             A reference to a valid
		 *             \link ProtectedMinutiaeTemplate\endlink
		 */
		static void swap
		( ProtectedMinutiaeTemplate & vault1 ,
		  ProtectedMinutiaeTemplate & vault2 );

		/**
		 * @brief
		 *            Returns a boolean value that indicates whether this
		 *            object has been initialized such that it can protect an
		 *            absolutely pre-aligned minutiae template.
		 *
		 * @details
		 *            After an object of this class has been created via the
		 *            standard
		 *            constructor \link ProtectedMinutiaeTemplate()\endlink
		 *            the width, height, resolution to which a to-be-protected
		 *            minutiae view corresponds to are unknown. However, this
		 *            information must be known for the vault's attributed to
		 *            be initialized properly and the member is_initialized
		 *            equals <code>false</code> as long as they have not
		 *            specified and, otherwise, <code>true</code>.
		 *
		 * @return
		 *            <code>true</code> if essential information about the
		 *            fingerprint images from which (to-be-)protected minutiae
		 *            templates have been extracted and, otherwise,
		 *            <code>false</code>.
		 */
		bool isInitialized() const;

		/**
		 * @brief
		 *            Determine whether this instance does hold a
		 *            protected minutiae template.
		 *
		 * @return
		 *            <code>true</code> if this
		 *            \link ProtectedMinutiaeTemplate\endlink does protect
		 *            a minutiae template; otherwise <code>false</code>.
		 *
		 * @see ProtectedMinutiaeTemplate::enroll()
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
		 *            \link ProtectedMinutiaeTemplate\endlink has been
		 *            encrypted via
		 *            \link ProtectedMinutiaeTemplate::encrypt\endlink;
		 *            otherwise, the function returns <code>false</code>.
		 *            If this \link ProtectedMinutiaeTemplate\endlink
		 *            has been decrypted via
		 *     \link ProtectedMinutiaeTemplate::decrypt(const AES128&)\endlink
		 *            such that it represents a candidate
		 *            for a decrypted \link ProtectedMinutiaeTemplate\endlink,
		 *            the function returns <code>false</code>.
		 *
		 * @return
		 *            <code>true</code> if this
		 *            \link ProtectedMinutiaeTemplate\endlink does protect
		 *            a minutiae template, is additionally protected
		 *            by an \link AES128\endlink key and does not represent
		 *            a candidate for a decrypted
		 *            \link ProtectedMinutiaeTemplate\endlink; otherwise, the
		 *            function returns <code>false</code>.
		 *
		 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
		 */
		bool isEncrypted() const;

		/**
		 * @brief
		 *            Determine whether this instance holds a decrypted
		 *            candidate for a protected minutiae template.
		 *
		 * @details
		 *            Any encrypted \link ProtectedMinutiaeTemplate\endlink
		 *            can be decrypted via an \link AES128\endlink key and,
		 *            subsequently, contains a candidate for the decrypted
		 *            \link ProtectedMinutiaeTemplate\endlink.
		 *
		 *            Note that, if this
		 *            \link ProtectedMinutiaeTemplate\endlink
		 *            has not been decrypted via
		 *    \link ProtectedMinutiaeTemplate::decrypt(const AES128&)\endlink
		 *            but contains enrolled data, the function returns
		 *            <code>true</code>.
		 *
		 * @return
		 *            <code>true</code> if this
		 *            \link ProtectedMinutiaeTemplate\endlink contains a
		 *            candidate for the unencrypted
		 *            \link ProtectedMinutiaeTemplate\endlink; otherwise,
		 *            the result will be <code>false</code>.
		 *
		 * @see ProtectedMinutiaeTemplate::decrypt(const AES128&)
		 */
		bool isDecrypted() const;

		/**
		 * @brief
		 *            Determine whether this instance does hold
		 *            encrypted data of a protected minutiae template.
		 *
		 * @details
		 *            Note that a \link ProtectedMinutiaeTemplate\endlink
		 *            keeps the encrypted data even after being decrypted
		 *            via \link ProtectedMinutiaeTemplate::decrypt\endlink
		 *
		 * @return
		 *            <code>true</code> if the
		 *            \link ProtectedMinutiaeTemplate\endlink contains
		 *            encrypted data; otherwise, the function returns
		 *            <code>false</code>.
		 *
		 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
		 */
		bool containsEncryptedData() const;

		/**
		 * @brief
		 *            Create protected data from a minutiae template.
		 *
		 * @details
		 *            The function passes the input minutiae template
		 *            through a quantization process and protects the
		 *            quantization via the improved fuzzy vault scheme
		 *            of which data will be held by this instance.
		 *
		 *            If not sufficient information could be extracted
		 *            from the input minutiae template to allow successful
		 *            verification, the function returns <code>false</code>;
		 *            otherwise, if successul, the function returns
		 *            <code>true</code>.
		 *
		 *
		 * @param view
		 *            A reference to a valid \link MinutiaeView\endlink
		 *            representing the input (absolutely pre-aligned)
		 *            minutiae template.
		 *
		 * @return
		 *            <code>true</code> if the creation was successful;
		 *            otherwise, the function returns <code>false</code>.
		 *
		 * @warning
		 *            If \link isEnrolled()\endlink returns
		 *            <code>true</code>, an error message is printed to
		 *            <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool enroll( const MinutiaeView & view );

		/**
		 * @brief
		 *             Performs a verification procedure to compare
		 *             the input minutiae template with the minutiae
		 *             template protected by this instance.
		 *
		 * @details
		 *             The function passes the input minutiae template
		 *             through a quantization process from which an instance
		 *             of the Reed-Solomon decoding problem is generated which
		 *             is attempted to be solved using the
		 * \link ProtectedMinutiaeTemplate::decode()\endlink
		 *             function. If successful, the verification is deemed
		 *             to be successful and, otherwise, failed.
		 *
		 * @param view
		 *             The (absolutely pre-aligned) query
		 *             minutiae template.
		 *
		 * @return
		 *             <code>true</code> if the verification procedure
		 *             was successful; otherwise, if the verification
		 *             procedure was not successful, the function returns
		 *             <code>false</code>.
		 *
		 * @warning
		 *             If this \link ProtectedMinutiaeTemplate\endlink
	     *             does not represent a successfully enrolled
	     *             decrypted protected minutiae template, i.e., if
	     *      \link ProtectedMinutiaeTemplate::isDecrypted()\endlink
	     *             returns <code>false</code>, the function
	     *             prints an error message to <code>stderr</code> and
	     *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool verify( const MinutiaeView & view ) const;

		/**
		 * @brief
		 *             Encrypts this instance using a key.
		 *
		 * @details
		 *             After successful enrollment (via
		 *             \link ProtectedMinutiaeTemplate::enroll()\endlink),
		 *             the \link ProtectedMinutiaeTemplate\endlink can
		 *             additionally protected via a (user-specific)
		 *             key derived from a PIN or password. Consequently,
		 *             the security level of the
		 *             \link ProtectedMinutiaeTemplate\endlink can be
		 *             improved, though this requires the
		 *             \link ProtectedMinutiaeTemplate\endlink to be
		 *             decrypted correctly for a successful verification,
		 *             thereby requiring that the correct can be
		 *             reproduced for successful verification.
		 *
		 * @param key
		 *             The \link AES128\endlink key used to encrypt the
		 *             \link ProtectedMinutiaeTemplate\endlink.
		 *
		 * @warning
		 *             If no minutiae template has been enrolled with
		 *             this \link ProtectedMinutiaeTemplate\endlink
		 *             (i.e., if
		 * \link ProtectedMinutiaeTemplate::isEnrolled()\endlink
		 *             returns <code>false</code>)
		 *             or if it contains already encrypted data (i.e., if
		 * \link ProtectedMinutiaeTemplate::containsEncryptedData()\endlink
		 *             returns <code>true</code>)
		 *             an error message is printed to <code>stderr</code>
		 *             and the program exits with status 'EXIT_FAILURE'.
		 *
		 * @see decrypt()
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
		 *             \link ProtectedMinutiaeTemplate\endlink that has been
		 *             encrypted before using
		 *             \link ProtectedMinutiaeTemplate::encrypt(const AES128&)\endlink.
		 *
		 *             If the \link ProtectedMinutiaeTemplate\endlink has been
		 *             encrypted with the key \f$\kappa\f$, then decryption
		 *             with \f$\kappa\f$ will correctly reveal the original
		 *             data; otherwise, if decrypting with a different
		 *             key \f$\kappa'\neq\kappa\f$, this instance will hold
		 *             a candidate for the original data that is (with very
		 *             high probability) different from the original data.
		 *             To date no efficient way for distinguishing incorrect
		 *             candidate data from the correct data is known. Hence,
		 *             the difficulty to break an encrypted
		 *             \link ProtectedMinutiaeTemplate\endlink equals the product
		 *             of the difficulty in guessing a  minutiae template being
		 *             of sufficient similarity to the protected minutiae template
		 *             and the difficulty in guessing the correct password.
		 *
		 * @param key
		 *             The \link AES128\endlink key used to decrypt this
		 *             \link ProtectedMinutiaeTemplate\endlink.
		 *
		 * @warning
		 *             If this instance does not hold encrypted data
		 *             (i.e., if
		 *    \link ProtectedMinutiaeTemplate::containsEncryptedData()\endlink
		 *             returns <code>false</code>), an error message is
		 *             printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @see ProtectedMinutiaeTemplate::encrypt(const AES128&)
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void decrypt( const AES128 & key );

		/**
		 * @brief
		 *             Use a minutiae template as a query to open the
		 *             vault.
		 *
		 * @details
		 *             This a special version of
		 *             the \link verify()\endlink function (in
		 *             fact, \link verify()\endlink
		 *             wraps around this function) in which the
		 *             secret polynomial that is used to hide the minutiae template
		 *             protected by this instance is revealed on successful
		 *             verification.
		 *
		 * @param f
		 *             Will be set to the secret polynomial on successful
		 *             verification; otherwise the content will be left
		 *             unchanged
		 *
		 * @param view
		 *             The (absolutely pre-aligned) query minutiae template.
		 *
		 * @return
		 *             <code>true</code> if the verification procedure was
		 *             successful; otherwise, if the verification was not
		 *             successful, the function returns <code>false</code>.
		 *
		 * @warning
		 *             If this \link ProtectedMinutiaeTemplate\endlink
	     *             does not represent a successfully enrolled
	     *             decrypted protected minutiae template, i.e.,
	     *             if \link ProtectedMinutiaeTemplate::isDecrypted()\endlink
	     *             returns <code>false</code>, the function
	     *             prints an error message to <code>stderr</code> and
	     *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		bool open( SmallBinaryFieldPolynomial & f , const MinutiaeView & view ) const;

		/**
		 * @brief
		 *             Use a quantized minutiae template as a query to open the
		 *             vault.
		 *
		 * @details
		 *             This is a variant of the \link open()\endlink
		 *             function in which a set of quantized minutiae is
		 *             used to open the vault. Therein, a quantized minutia is
		 *             encoded as an element of the finite field that can be
		 *             accessed via \link getField()\endlink.
         *
         *             To obtain the feature set of an absolutely pre-aligned
         *             minutiae template, the \link quantize()\endlink method
         *             can be used.
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
		 *             If this \link ProtectedMinutiaeTemplate\endlink
	     *             does not represent a successfully enrolled and
	     *             decrypted protected minutiae template, i.e.,
	     *             if \link isDecrypted()\endlink
	     *             returns <code>false</code>, the function
	     *             prints an error message to <code>stderr</code> and
	     *             exits with status 'EXIT_FAILURE'.
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
		 *             Implementation of the <i>randomized decoder</i> to
		 *             solve an instance of the Reed-Solomon decoding problem.
		 *
		 * @details
		 *             The function randomly selects <i>k</i> distinct
		 *             pairs <i>(x[j],y[j])</i>, computes its interpolation
		 *             polynomial, checks if its 160 bit hash value agrees
		 *             with <i>hash</i>, and, if true, assigns <i>f</i> by
		 *             the interpolation polynomial returning
		 *             <code>true</code>; otherwise, the function repeats
		 *             the procedure up to <i>D</i> times and returns
		 *             <code>false</code> if none of them yielded the
		 *             correct hash value.
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
		 * @param D
		 *             Number of decoding iterations
		 *
		 * @return
		 *             <code>true</code> if the decoding attempt was successful;
		 *             otherwise, the function returns <code>false</code>.
		 */
		bool decode
		( SmallBinaryFieldPolynomial & f , const uint32_t *x , const uint32_t *y ,
		  int t , int k , const uint8_t hash[20] , int D ) const;

		/**
		 * @brief
		 *             Computes the quantization of an (absolutely pre-aligned)
		 *             minutia.
		 *
		 * @param minutia
		 *             The minutia of which the quantization is computed.
		 *
		 * @return
		 *             The quantization of <code>minutia</code> encoded as a
		 *             finite field element.
		 */
		uint32_t quantize( const Minutia & minutia ) const;

		/**
		 * @brief
		 *             Compute the quantization of an (absolutely pre-aligned)
		 *             minutiae template encoded as a set of finite field
		 *             elements.
		 *
		 * @details
		 *             At most the first <code>tmax</code> quantizations of
		 *             the minutiae in <code>view</code> are computed and
		 *             added to <code>array</code> (if not contained already).
		 *             The order of minutiae is given by their estimated
		 *             quality (if any), thereby guaranteeing that the minutiae
		 *             of the highest estimated quality are quantized first and
		 *             so on.
		 *
		 *             The function returns the number <i>t</i> of distinct
		 *             minutiae  quantizations that have been put into
		 *             <code>array</code> which depends on the size of the input
		 *             minutiae template
		 *             (see \link MinutiaeView::getMinutiaeCount()\endlink)
		 *             , <code>tmax</code> and on the count how many minutiae
		 *             had equal quantizations. Thereby, the function
		 *             guarantees that the first <i>t</i> elements of
		 *             <code>array</code> are distinct.
		 *
		 * @param array
		 *             Will contain at most <code>tmax</code>
		 *             finite field elements encoding the first minutiae'
		 *             quantization in <code>view</code>.
		 *
		 * @param view
		 *             The minutiae template of which the quantization set
		 *             is computed.
		 *
		 * @param tmax
		 *             Maximal number of minutiae quantizations to be extracted
		 *             from <code>view</code>.
		 *
		 * @return
		 *             The number of minutiae quantizations extracted from
		 *             <code>view</code>.
		 *
		 * @warning
		 *             If <code>tmax</code> is not positive, an error message
		 *             is printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>array</code> cannot hold <code>tmax</code>
		 *             elements of type <code>uint32_t</code> calling this
		 *             function runs into undocumented behavior.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		int quantize( uint32_t *array , const MinutiaeView & view , int tmax ) const;

		/**
		 * @brief
		 *             Compute the quantization of an (absolutely pre-aligned)
		 *             minutiae template encoded as a set of finite field
		 *             elements.
		 *
		 * @details
		 *             At most the first <code>tmax</code> quantizations of
		 *             the minutiae in <code>view</code> are computed and
		 *             added to <code>array</code> (if not contained already).
		 *             The order of minutiae is given by their estimated
		 *             quality (if any), thereby guaranteeing that the minutiae
		 *             of the highest estimated quality are quantized first and
		 *             so on.
		 *
		 *             The function returns the number <i>t</i> of distinct
		 *             minutiae  quantizations that have been put into
		 *             <code>array</code> which depends on the size of the input
		 *             minutiae template
		 *             (see \link MinutiaeView::getMinutiaeCount()\endlink)
		 *             , <code>tmax</code> and on the count how many minutiae
		 *             had equal quantizations. Thereby, the function
		 *             guarantees that the first <i>t</i> elements of
		 *             <code>array</code> are distinct.
		 *
		 * @param array
		 *             Will contain at most \link getMaxGenuineFeatures()\endlink
		 *             finite field elements encoding the first minutiae'
		 *             quantization in <code>view</code>.
		 *
		 * @param view
		 *             The minutiae template of which the quantization set
		 *             is computed.
		 *
		 * @return
		 *             The number of minutiae quantizations extracted from
		 *             <code>view</code>.
		 *
		 * @warning
		 *             If <code>tmax</code> is not positive, an error message
		 *             is printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If <code>array</code> cannot hold
         *             \link getMaxGenuineFeatures()\endlink
		 *             elements of type <code>uint32_t</code>, calling this
		 *             function runs into undocumented behavior.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		int quantize( uint32_t *array , const MinutiaeView & view ) const;

		/**
		 * @brief
		 *             Evaluate the reordering of a quantized feature which
		 *             is associated with an (enrolled) protected minutiae
		 *             template.
		 *
		 * @details
		 *             To prevent the application of extended Euclidean algorithm
		 *             based record multiplicity attack, each protected minutiae
		 *             template is attached with a public random permutation
		 *             process acting on the quantized feature universe.
		 *             The reordering of such a quantized feature (encoded as
		 *             a finite field element) can be evaluated via this
		 *             function.
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
		 *             'EXIT_FAILURE'.
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
		 *             Access the width of the fingerprint images related
		 *             to this protected minutiae template.
		 *
		 * @return
		 *             The width of the fingerprint images related to this
		 *             protected minutiae template.
		 *
		 * @see setDimension(int,int)
		 */
		int getWidth() const;

		/**
		 * @brief
		 *             Access the height of the fingerprint images related
		 *             to this protected minutiae template.
		 *
		 * @return
		 *             The height of the fingerprint images related to this
		 *             protected minutiae template.
		 *
		 * @see setDimension(int,int)
		 */
		int getHeight() const;

		/**
		 * @brief
		 *             Access the finger position of this protected
		 *             minutiae template.
		 *
		 * @return
		 *             The finger position of this protected minutiae
		 *             template.
		 *
		 * @see setFingerPosition(FINGER_POSITION_T)
		 */
		FINGER_POSITION_T getFingerPosition() const;

		/**
		 * @brief
		 *             Access the value indicating the assumed resolution
		 *             in <i>dots per inch</i> at which the fingerprint
		 *             images related with this instance have been scanned.
		 *
		 * @details
		 *             The resolution in dots per inch.
		 *
		 * @see setResolution(int)
		 */
		int getResolution() const;

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
		 *             extracted from the to-be-protected minutiae template
		 *             on enrollment.
		 *
		 * @return
		 *             The maximal number of minutiae quantizations extracted
		 *             from the to-be-protected minutiae template on
		 *             enrollment.
		 *
		 * @see setMaxGenuineFeatures(int)
		 */
		int getMaxGenuineFeatures() const;

		/**
		 * @brief
		 *             Access the number of possible different minutiae
		 *             quantizations.
		 *
		 * @return
		 *             The number of possible different minutiae quantizations.
		 */
		int getVaultSize() const;

		/**
		 * @brief
		 *             Access the number of hexagonal grid points centered in
		 *             a region in which absolutely pre-aligned minutiae can
		 *             occur.
		 *
		 * @return
		 *             The number of hexagonal gird points covering the
		 *             region in which absolutely pre-aligned minutiae
		 *             can occur.
		 */
		int getGridSize() const;

		/**
		 * @brief
		 *             Access the number of decoding iterations performed
		 *             by the randomized decoder on verification attempts
		 *             with this protected minutiae template.
		 *
		 * @return
		 *             The number of decoding iterations performed by the
		 *             randomized decoder on verification with this protected
		 *             minutiae template.
		 *
		 * @see setNumberOfDecodingIterations(int)
		 */
		int getNumberOfDecodingIterations() const;

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
		 *             Access the finite field in which computations with
		 *             the improved fuzzy vault scheme are performed to
		 *             protect a minutiae template.
		 *
		 * @returns
		 *             The \link SmallBinaryField finite field\endlink in which
		 *             computations with the improved fuzzy vault scheme are
		 *             performed to protect a minutiae template.
		 */
		const SmallBinaryField & getField() const;

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
		 *             Access the SHA-1 hash value of the secret key
		 *             generated on enrollment used to obfuscate the
		 *             correct minutiae template's quantization.
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
		 *             Change the dimension of the fingerprint images
		 *             related with this protected minutiae template
		 *
		 * @details
		 *             Changing the dimension can cause to update the
		 *             \link getGrid() grid\endlink and the
		 *             \link getField() finite field\endlink related with this
		 *             protected minutiae template.
		 *
		 * @param width
		 *             The width of the fingerprint image.
		 *
		 * @param height
		 *             The height of the fingerprint image.
		 *
		 * @warning
		 *             If \link isEnrolled()\endlink returns
		 *             <code>true</code>, an error message is printed to
		 *             <code>stderr</code> and the program exits with status
		 *             -1.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void setDimension( int width , int height );

		/**
		 * @brief
		 *             Change the finger position of this protected
		 *             minutiae template.
		 *
		 * @param fingerPosition
		 *             The new finger position of this protected minutiae
		 *             template.
		 *
		 * @see getFingerPosition()
		 */
		void setFingerPosition( FINGER_POSITION_T fingerPosition );

		/**
		 * @brief
		 *             Change the value indicating the assumed resolution
		 *             in <i>dots per inch</i> at which the fingerprint
		 *             images related with this instance have been scanned.
		 *
		 * @details
		 *             The new resolution in dots per inch.
		 *
		 * @see getResolution(int)
		 *
		 * @warning
		 *             If <code>dpi</code> is smaller than 300 or greater
		 *             than 1000, an error message is printed to
		 *             <code>stderr</code> and the program exits with
		 *             status 'EXIT_FAILURE'.
		 */
		void setResolution( int dpi );

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
		 *             Change the maximal number of minutiae quantizations
		 *             extracted from the to-be-protected minutiae template
		 *             on enrollment.
		 *
		 * @param tmax
		 *             The new maximal number of minutiae quantizations extracted
		 *             from the to-be-protected minutiae template on
		 *             enrollment.
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
		void setMaxGenuineFeatures( int tmax );

		/**
		 * @brief
		 *             Change the number of decoding iterations performed
		 *             by the randomized decoder on verification attempts
		 *             with this protected minutiae template.
		 *
		 * @return
		 *             The new number of decoding iterations that are performed
		 *             by the randomized decoder on verification with this
		 *             protected minutiae template.
		 *
		 * @see getNumberOfDecodingIterations()
		 *
		 * @warning
		 *             If <i>D</i> is smaller than or equals 0, the method
		 *             prints an error message to <code>stderr</code> and
		 *             exits with status 'EXIT_FAILURE'.
		 */
		void setNumberOfDecodingIterations( int D );

		/**
		 * @brief
		 *            Specifies the slow-down factor manually to artifically
		 *            increase the time for an impostor verification attempt.
		 *
		 * @details
		 *            The time needed for an impostor to perform
		 *            off-line attacks is multiplied by the slow-down
		 *            factor, thereby increasing security while also
		 *            increasing response time of genuine verification
		 *            attempts---a trade-off thus but with the capability
		 *            of maintaining a relative security measure with
		 *            the advent of more and more powerful computers.
		 *
		 * @param slowDownFactor
		 *            The specified slow-down factor for this protected
		 *            minutiae template.
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
		 *             Clears all data that is related with a minutiae
		 *             template protected by this instance (if any).
		 *
		 * @details
		 *             Any \link ProtectedMinutiaeTemplate\endlink is empty
		 *             (cleared) by default until a
		 *             \link MinutiaeView minutiae template\endlink is passed
		 *             through the \link enroll()\endlink function returning
		 *             <code>true</code>. If the user wants to reuse this
		 *             \link ProtectedMinutiaeTemplate\endlink to protect another
		 *             \link MinutiaeView minutiae template\endlink via
		 *             \link enroll()\endlink, he must clear it first via
		 *             this method.
		 */
		void clear();

		/**
		 * @brief
		 *             Determine the size in bytes needed to store this
		 *             protected minutiae template.
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
		 *             Store data defining this protected minutiae template to
		 *             the specified byte array.
		 *
		 * @param data
		 *             The byte array to which the data for this protected
		 *             minutiae template is written.
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
		 *             Initialize this protected minutiae template by the
		 *             specified bytes.
		 *
		 * @details
		 *             If this protected minutiae template has been successfully
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
		 *             Write data defining this protected minutiae template to
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
		 *             Write data defining this protected minutiae template to
		 *             the specified file.
		 *
		 * @param file
		 *             Path to the file to which this protected minutiae
		 *             template is written.
		 *
		 * @return
		 *             The function returns <code>true</code> if this protected
		 *             minutiae template has been successfully be written
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
		 *             Initialize this protected minutiae template from the
		 *             specified file pointer.
		 *
		 * @details
		 *             If this protected minutiae template has been successfully
		 *             initialized, the function returns the number read bytes;
		 *             otherwise, if an error occurred, this protected minutiae
		 *             template will be left unchanged and the result of the
		 *             function will be -1.
		 *
		 * @param in
		 *             The file pointer from which the protected minutiae
		 *             template is read.
		 *
		 * @return
		 *             The number of bytes that were successfully read from
		 *             <code>in</code> to initialize this protected minutiae
		 *             template; otherwise, if an error occurred,
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
		 *             Initialize this protected minutiae template from the
		 *             data contained in the specified file path.
		 *
		 * @details
		 *             If this protected minutiae template has been successfully
		 *             initialized, the function returns <code>true</code>;
		 *             otherwise, if an error occurred, this protected minutiae
		 *             template will be left unchanged and the result of the
		 *             function will be <code>false</code>.
		 *
		 * @param file
		 *             String representation of the path to the file from
		 *             which the protected minutiae template is read.
		 *
		 * @return
		 *             <code>true</code> if this protected minutiae template
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
		 *            Convenience method to determine how many elements
		 *            two arrays share.
		 *
		 * @details
		 *            The function computes the sum
		 *            \f[
		 *             \sum_{i=0}^{t1-1}\sum_{j=0}^{t2-1}
		 *                \delta(array1[i],array2[j])
		 *            \f]
		 *            where \f$\delta(a,b)=1\f$ if \f$a=b\f$ and
		 *            \f$\delta(a,b)\f$ if \f$a\neq b\f$.
		 *            <br><br>
		 *            Consequently, if <code>array1</code> and
		 *            <code>array2</code> both contain pairwise distinct
		 *            elements if we let
		 *            \f$A=\{array1[0],...,array1[t1-1]\}\f$
		 *            and \f$B=\{array2[0],...,array2[t2-1]\}\f$ the result
		 *            will be \f$|A\cap B|\f$.
		 *
		 * @param array1
		 *            First array.
		 *
		 * @param t1
		 *            Size of first array.
		 *
		 * @param array2
		 *            Second array.
		 *
		 * @param t2
		 *            Size of second array.
		 *
		 * @return
		 *            The number of elements which are contained
		 *            simultaneously contained in <code>array1</code>
		 *            and <code>array2</code>.
		 *
		 * @warning
		 *            If <code>array1</code> and <code>array2</code> do not
		 *            contain at least <code>t1</code> and <code>t2</code>
		 *            valid <code>uint32_t</code>, respectively, the
		 *            function runs into undocumented behavior.
		 */
		static int overlap
		( const uint32_t *array1 , int t1 , const uint32_t *array2 , int t2 );

	private:

		/**
		 * @brief
		 *            The width of the fingerprint images related with this
		 *            protected minutiae template.
		 *
		 * @details
		 *            The value of this field must be specified on
		 *            constructing a valid instance via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The function
		 *            \link getWidth()\endlink serves as the getter for this
		 *            field and can be modified via the
		 *            \link setDimension(int,int)\endlink function.
		 *
		 * @see ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)
		 * @see getWidth()
		 * @see setDimension(int,int)
		 */
		int width;

		/**
		 * @brief
		 *            The height of the fingerprint images related with this
		 *            protected minutiae template.
		 *
		 * @details
		 *            The value of this field is must be specified on
		 *            constructing a valid instance via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The function
		 *            \link getHeight()\endlink serves as the getter for this
		 *            field and can be modified via the
		 *            \link setDimension(int,int)\endlink function.
		 *
		 * @see ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)
		 * @see getHeight()
		 * @see setDimension(int,int)
		 */
		int height;

		/**
		 * @brief
		 *            The position of the finger from which the minutiae
		 *            template protected by this instance has been estimated.
		 *
		 * @details
		 *            The value of this field can (but does not need to) be
		 *            specified explicitly on constructing an instance
		 *            (UNKNOWN_FINGER per default) via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The function
		 *            \link getFingerPosition()\endlink and method
		 *            \link setFingerPosition(FINGER_POSITION_T)\endlink serve
		 *            as the getter and setter of this field, respectively.
		 *
		 * @see ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)
		 * @see getFingerPosition()
		 * @see setFingerPosition(FINGER_POSITION_T)
		 */
		FINGER_POSITION_T fingerPosition;

		/**
		 * @brief
		 *            The resolution in dots per inch at which the fingerprint
		 *            images related with this instance have been scanned.
		 *
		 * @details
		 *            The value of this field can (but does not need to) be specified
		 *            explicitly on constructing an instance (500 per default) via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The function
		 *            \link getResolution()\endlink and method
		 *            \link setResolution(int)\endlink serve as the
		 *            getter and setter of this field, respectively.
		 *
		 * @see ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)
		 * @see getResolution()
		 * @see setResolution(int)
		 */
		int dpi;

		/**
		 * @brief
		 *            The distance of the hexagonal grid used for minutiae
		 *            quantization.
		 *
		 * @details
		 *            The value of this field is specified implicitly on
		 *            constructing an instance (via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink)
		 *            in dependence of the specified
		 *            resolution and is 25 pixels for a default resolution of
		 *            500 dots per inch. The function
		 *            \link getGridDist()\endlink and method
		 *            \link setGridDist(int)\endlink serve as the getter and
		 *            setter of this field, respectively.
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
		 *            The value of this field is automatically set to 6 per
		 *            default on constructing an instance via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            Its setter and getter are the method
		 *            \link setNumAngleQuanta(int)\endlink and the function
		 *            \link getNumAngleQuanta()\endlink, respectively.
		 *
		 * @see setNumAngleQuanta(int)
		 * @see getNumAngleQuanta()
		 */
		int s;

		/**
		 * @brief
		 *            The size of the secret polynomial being used to obfuscate
		 *            the minutiae template's quantization in order to protect
		 *            it.
		 *
		 * @details
		 *            The value of this field is automatically set to 12
		 *            per default when an instance is constructed via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            Its setter and getter are the method
		 *            \link setSecretSize(int)\endlink and function
		 *            \link getSecretSize()\endlink, respectively.
		 *
		 * @see setSecretSize(int)
		 * @see getSecretSize()
		 */
		int k;

		/**
		 * @brief
		 *            The maximal number of minutiae quantizations extracted
		 *            from the to-be-protected minutiae template on enrollment.
		 *
		 * @details
		 *            The value of this field is automatically set to 44
		 *            per default when an instance is constructed via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The setter and getter of this field are the method
		 *            \link setMaxGenuineFeatures(int)\endlink and the function
		 *            \link getMaxGenuineFeatures()\endlink.
		 *
		 * @see setMaxGenuineFeatures(int)
		 * @see getMaxGenuineFeatures()
		 */
		int tmax;

		/**
		 * @brief
		 *            The number of decoding iterations performed by the
		 *            randomized decoder on a verification attempt with this
		 *            protected minutiae template.
		 *
		 * @details
		 *            The value of this field is automatically set to
		 *            \f$2^{16}\f$ per default when an instance is constructed
		 *            via
		 * \link ProtectedMinutiaeTemplate(int,int,int,FINGER_POSITION_T)\endlink.
		 *            The setter and getter of this field are
		 *            \link setNumberOfDecodingIterations(int)\endlink and
		 *            \link getNumberOfDecodingIterations()\endlink,
		 *            respectively.
		 *
		 * @see setNumberOfDecodingIterations(int)
		 * @see getNumberOfDecodingIterations()
		 */
		int D;

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
		 *            (e.g., via \link setDimension(int,int)\endlink).
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
		 */
		SmallBinaryField *gfPtr;

		/**
		 * @brief
		 *            The size of the feature set protected by this instance.
		 *
		 * @details
		 *            On \link enroll() enrollment\endlink a
		 *            \link MinutiaeView minutiae template\endlink is used to
		 *            extract a feature set via \link quantize()\endlink
		 *            which is protected using the improved fuzzy vault scheme.
		 *            This field is actually deemed to be of the size of the
		 *            protected feature set. However, in this implementation,
		 *            the feature set is supplemented by a certain number of
		 *            blending feature elements such that it reaches a size
		 *            of \link tmax\endlink. Consequently, it is harder to
		 *            derive information about the (say) quality of the
		 *            protected minutiae template just from the intensity
		 *            of <i>t</i>.
		 */
		int t;

		/**
		 * @brief
		 *            The slow-down factor of this object.
		 *
		 * @see getSlowDownFactor()
		 * @see setSlowDownFactor()
		 */
		BigInteger slowDownFactor;

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
		 *            the enrolled minutiae template's quantization.
		 *
		 * @details
		 *            On \link enroll() enrollment\endlink,
		 *            a \link MinutiaeView minutiae template\endlink is bound
		 *            to a \link SmallBinaryFieldPolynomial
		 *            secret polynomial\endlink of which any data is dismissed
		 *            except of its SHA-1 hash value which is stored in this
		 *            field.
		 *
		 *            The getter for this field is the
		 *            \link getHash()\endlink function.
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
		 *            from non-related protected templates via the improved
		 *            fuzzy vault scheme and, even worse, he may break them
		 *            entirely with quite a high probability using an attack
		 *            based on the Euclidean algorithm (see
		 *            \link EEAAttack\endlink). To prevent successful
		 *            application of the attack, each protected minutiae
		 *            template is attached with a random public
		 *            record-specific permutation process acting on the
		 *            feature sets which is represented by this field.
		 *
		 *            This field is updated on each enrollment via
		 *            \link updatePermutation()\endlink. The permutation
		 *            can be accessed via the \link reorder()\endlink
		 *            function.
		 *
		 * @see reorder()
		 */
		Permutation permutation;

		/**
		 * @brief
		 *            Boolean value that indicates whether this object has
		 *            been initialized such that it can protect an absolutely
		 *            pre-aligned minutiae template.
		 *
		 * @details
		 *            After an object of this class has been created via the
		 *            standard
		 *            constructor \link ProtectedMinutiaeTemplate()\endlink
		 *            the width, height, resolution to which a to-be-protected
		 *            minutiae view corresponds to are unknown. However, this
		 *            information must be known for the vault's attributed to
		 *            be initialized properly and the member is_initialized
		 *            equals <code>false</code> as long as they have not
		 *            specified and, otherwise, <code>true</code>.
		 *
		 * @see isInitialized()
		 */
		bool is_initialized;

		/**
		 * @brief
		 *            Initializes the member variables of this protected
		 *            minutiae template for the very first time such that
		 *            it becomes ready for \link enroll() enrollment\endlink.
		 *
		 * @details
		 *            The method assumes that the member variables
		 *            \link width\endlink, \link height\endlink,
		 *            \link gridDist\endlink, \link s\endlink have been
		 *            properly initialized already; otherwise, the program
		 *            runs into undocumented behavior.
		 *
		 *            The methods should not be called directly as it is
		 *            called by other initialization methods of this class.
		 *
		 *            The method initializes the grid and the finite field
		 *            using the \link updateGrid()\endlink and
		 *            \link updateField()\endlink functions, respectively.
		 */
		void first_init();

		/**
		 * @brief
		 *             Initializes this object such that it can protected
		 *             absolutely pre-aligned minutiae templates.
		 *
		 * @details
		 *             This method generates parameters needed to properly
		 *             quantize absolutely pre-aligned minutiae templates that
		 *             have been estimated from a fingerprint image of the
		 *             specified dimension and that has been scanned at the
		 *             specified resolution. If the finger's position is known,
		 *             it can be specified as well.
		 *
		 *             This method should be called directly only if no memory
		 *             is reserved by the member of this object such that
		 *             allocation on the attributes does not cause information
		 *             leakage.
		 *
		 * @param width
		 *             The width of the fingerprint image.
		 *
		 * @param height
		 *             The height of the fingerprint image.
		 *
		 * @param dpi
		 *             The resolution in (dots per inch) at which
		 *             the fingerprint image has been scanned.
		 *
		 * @param fingerPosition
		 *             The position of the finger (if known).
		 *
		 * @attention
		 *             If <code>width</code> or <code>height</code> are
		 *             not greater than 0 or if <code>dpi</code> does
		 *             not range between 300 and 1000, an error message
		 *             is printed to <code>stderr</code> and the program
		 *             exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void first_init
		( int width , int height , int dpi , FINGER_POSITION_T fingerPosition );

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
		 *            needs an update which is performed when calling
		 *            this method.
		 */
		void updateGrid();

		/**
		 * @brief
		 *            Update the finite field such that it can encode
		 *            each possible minutia's quantization.
		 *
		 * @details
		 *            After the \link grid hexagonal grid points\endlink
		 *            or the
		 *            \link s number of minutiae angle quantizations\endlink
		 *            have changed, the finite field might need an update
		 *            which can be performed by calling this method.
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
		 * @brief
	     *            Determines how many bytes are needed to hold
		 *            the data of the packed vault polynomial.
		 *
		 * @details
		 *            The function is used to determine how many bytes
		 *            are needed to initialize an array of 128 bit blocks
		 *            that can hold the data of an instance of the improved
		 *            fuzzy vault scheme (encoded as a monic polynomial)
		 *            such that it can be protected via an
		 *            \link AES128 AES key\endlink.
		 *
		 * @return
		 *            Number of bytes needed to hold the data of a monic
		 *            polynomial of specified degree having coefficients
		 *            in a binary finite field of specified degree such
		 *            that the byte array is a multiple of 16 bytes
		 *            (i.e., 128 bits).
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


#endif /* THIMBLE_PROTECTEDMINUTIAEVIEW_H_ */
