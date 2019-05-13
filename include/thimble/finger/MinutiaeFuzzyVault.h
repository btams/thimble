/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file MinutiaeFuzzyVault.h
 *
 * @brief
 *            Provides an implementation of a minutiae fuzzy vault.
 *
 * @details
 *            see \link thimble::MinutiaeFuzzyVault\endlink
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_MINUTIAEFUZZYVAULT_H_
#define THIMBLE_MINUTIAEFUZZYVAULT_H_

#include <stdint.h>
#include <vector>

#include <thimble/dllcompat.h>
#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/finger/MinutiaeRecord.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           An implementation of the minutiae fuzzy vault.
	 *
	 * @details
	 *           To create an instance of the class one must specify the
	 *           fingerprint image's width and height. For example the code
	 *           <pre>
	 *            MinutiaeFuzzyVault vault(296,560);
	 *           </pre>
	 *           creates a vault to protect the minutiae whose horizontal
	 *           and vertical position vary between 0,...,295 and
	 *           0,...,559, respectively. To protect the minutiae of a finger
	 *           given as
	 *           <pre>
	 *            MinutiaeView minutiae;
	 *           </pre>
	 *           we may run
	 *           <pre>
	 *            bool success = vault.enroll(minutiae);
	 *           </pre>
	 *           If the enrollment was successful, the variable
	 *           <code>success</code> will be <code>true</code>. Otherwise
	 *           it will be <code>false</code> in which case <code>vault</code>
	 *           does not protect the template specified by
	 *           <code>minutiae</code>. For the following, we assume that the
	 *           enrollment was successful.
	 *           <br>
	 *           Now, as aligned query minutiae
	 *           <pre>
	 *            MinutiaeView alignedQuery;
	 *           </pre>
	 *           are given we can match against the protected template
	 *           <code>vault</code> via
	 *           <pre>
	 *            SmallBinaryFieldPolynomial f(vault.getField());
	 *            bool success = vault.open(f,alignedQuery);
	 *           </pre>
	 *           If successful, <code>success</code> will be <code>true</code>.
	 *           In this case, <code>f</code> contains a secret that was
	 *           randomly generated on enrollment and binded to the enrolled
	 *           minutiae. Otherwise, if <code>success</code> is
	 *           <code>false</code> the query minutiae could not be used
	 *           to successfully authenticate versus the protected template
	 *           <code>vault</code>. In this case, <code>f</code> will be left
	 *           unchanged.
	 *           <br><br>
	 *           The implementation mainly follows the description of
	 *           <ul>
	 *            <li>
	 *             <b>K. Nandakumar, A. K. Jain, and S.Pankanti (2007)</b>,
	 *             &quot;Fingerprint-based fuzzy vault: Implementation
	 *             and performance.&quot; <i>IEEE Trans. Inf. Forensics
	 *             Security</i>, vol. 2, no. 4, pp. 744-757.
	 *            </li>
	 *           </ul>
	 *           with some minor modifications. For further details we refer
	 *           to
	 *           <ul>
	 *            <li>
	 *            <b>Tams, B. (2013)</b>. Attacks and Countermeasures in
	 *            Fingerprint Based Biometric Cryptosystems,
	 *            <a href="http://arxiv.org/abs/1304.7386">arXiv:1304.7386</a>.
	 *            </li>
	 *           </ul>
	 *           We did not implement an automatic alignment step but assume
	 *           that the query minutiae templates are already aligned to the
	 *           vault which has to be implement manually by the user of the
	 *           class. Currently, proposals which provide automatic alignment
	 *           from a security perspective are not solved in a satisfactory
	 *           manner.
	 *           <br><br>
	 */
	class THIMBLE_DLL MinutiaeFuzzyVault {

	private:

		/**
		 * @brief
		 *           Indicates whether the library is advised to use a
		 *           cryptographically secure random generator.
		 *
		 * @details
		 *           This member can be specified when the constructor of the
		 *           class is called.
		 *
		 * @see MinutiaeFuzzyVault(int,int,bool)
		 */
		bool tryRandom;

		/**
		 * @brief
		 *           Weight specifying how significant minutiae angles are
		 *           taken into account when the similarity of two minutiae
		 *           is measured.
		 *
		 * @details
		 *           The default value of this member is 11.459.
		 *           <br><br>
		 *           This member can be modified using the
		 *           <code>setAngleWeight(double)</code> method. It can be
		 *           accessed via the <code>getAngleWeight()</code> function.
		 *
		 * @see getAngleWeight()
		 * @see setAngleWeight(double)
		 */
		double angleWeight;

		/**
		 * @brief
		 *           Horizontal dimension in where the abscissa coordinate
		 *           of protected minutiae vary.
		 *
		 * @details
		 *           This member can be specified when the constructor of the
		 *           class is called. It can be accessed via the
		 *           <code>getWidth()</code> function.
		 *
		 * @see MinutiaeFuzzyVault(int,int,bool)
		 * @see getWidth()
		 */
		int width;

		/**
		 * @brief
		 *           Vertical dimension where the ordinate coordinate of
		 *           protected minutiae vary.
		 *
		 * @details
		 *           This member can be specified when the constructor of the
		 *           class is called. It can be accessed via the
		 *           <code>getHeight()</code> function.
		 *
		 * @see MinutiaeFuzzyVault(int,int,bool)
		 * @see getHeight()
		 */
		int height;

		/**
		 * @brief
		 *           Size of the vault, i.e. number of vault minutiae/points
		 *           (chaff and genuine).
		 *
		 * @details
		 *           The default value of this member is 224.
		 *           <br><br>
		 *           This member can be accessed using the
		 *           <code>getVaultSize()</code> function and be modified
		 *           using the <code>setVaultSize(int)</code> method.
		 *
		 * @see getVaultSize()
		 * @see setVaultSize(int)
		 */
		int n;

		/**
		 * @brief
		 *           Lower bound on the number of well-separated genuine
		 *           minutiae required for successful vault enrollment,
		 *           i.e. lower bound on the size of the genuine vault points.
		 *
		 * @details
		 *           The default value of this member is 18.
		 *           <br><br>
		 *           This member can be modified using the
		 *           <code>setMinimalGenuineSetSize(int)</code> method and be
		 *           accessed via <code>getMinimaGenuineSetSize()</code>.
		 *
		 * @see setMinimalGenuineSetSize(int)
		 * @see getMinimalGenuineSetSize()
		 */
		int tmin;

		/**
		 * @brief
		 *           Bound on the maximal number of well-separated genuine
		 *           minutiae that are protected, i.e. upper bound on the
		 *           size of the genuine vault points.
		 *
		 * @details
		 *           The default value of this member is 24.
		 *           <br><br>
		 *           This member can be modified using the
		 *           <code>setMaximalGenuineSetSize(int)</code> method and be
		 *           accessed via <code>getMaximalGenuineSetSize()</code>.
		 *
		 * @see setMaximalGenuineSetSize(int)
		 * @see getMaximalGenuineSetSize()
		 */
		int tmax;

		/**
		 * @brief
		 *            Size of the secret polynomial which is binded to the
		 *            protected templates.
		 *
		 * @details
		 *            On enrollment, the secret polynomial is generated at
		 *            random and it can not be specified the user. The default
		 *            value of its number of coefficients, i.e. <i>k</i>, that
		 *            are generated is 9.
		 *            <br><br>
		 *            This member can be accessed via the
		 *            <code>getSecretSize()</code> function and be modified by
		 *            <code>setSecretSize(int)</code>.
		 *
		 * @see getSecretSize()
		 * @see setSecretSize(int)
		 */
		int k;

		/**
		 * @brief
		 *            Controls the minimal distance two vault minutiae
		 *            are allowed to have.
		 *
		 * @details
		 *            The default value of this member is 25.0 and can be
		 *            modified by
		 *            <code>setMinimalInterVaultDistance(double)</code>.
		 *            Furthermore, it can be accessed via
		 *            <code>getMinimalInterVaultDistance()</code>.
		 *
		 * @see getMinimalInterVaultDistance()
		 * @see setMinimalInterVaultDistance(double)
		 */
		double minInterVaultDistance;

		/**
		 * @brief
		 *            Controls the maximal distance of matching query minutiae
		 *            to vault minutiae.
		 *
		 * @details
		 *            The default value of this member is 30.0 and can be
		 *            modified by
		 *            <code>setMaximalMatchDistance(double)</code>.
		 *            Furthermore, it can be accessed via
		 *            <code>getMaximalMatchDistance()</code>.
		 *
		 * @see getMaximalMatchDistance()
		 * @see setMaximalMatchDistance(double)
		 */
		double maxMatchDistance;

		/**
		 * @brief
		 *           Will contain all vault minutiae sorted w.r.t.
		 *           lexicographical order.
		 *
		 * @details
		 *           Within the vector the genuine minutiae are dispersed
		 *           among the chaff minutiae. Because the minutiae are sorted
		 *           w.r.t. lexicographical order, it should be computationally
		 *           hard to distinguish genuine from chaff.
		 *           <br><br>
		 *           A constant reference to this member can be accessed by
		 *           <code>getVaultMinutiae()</code>.
		 *
		 * @see getVaultMinutiae()
		 */
		std::vector<Minutia> vaultMinutiae;

		/**
		 * @brief
		 *           List of abscissa coordinates of the vault points
		 *           corresponding to vault minutiae.
		 *
		 * @details
		 *           There is a one-to-one correspondence between vault
		 *           minutiae and vault points. More precisely, the vault
		 *           minutia <code>vautlMinutiae[i]</code> corresponds to the
		 *           vault point <code>(vaultX[i],vaultY[i])</code>.
		 *           <br><br>
		 *           The abscissa value <code>vaultX[i]</code> is equals
		 *           the finite field element that is encoded by the integer
		 *           <code>i</code>. Therefore, the presence of the array is
		 *           somewhat redundant be it is nonetheless here for the
		 *           reason of consistency.
		 *           <br><br>
		 *           The array is created automatically on enrollment and
		 *           can be accessed via <code>getVaultAbscissas()</code>.
		 *
		 * @see getVaultAbscissas()
		 */
		uint32_t *vaultX;

		/**
		 * @brief
		 *           List of ordinate coordinates of the vault points
		 *           corresponding to vault minutiae.
		 *
		 * @details
		 *           There is a one-to-one correspondence between vault
		 *           minutiae and vault points. More precisely, the vault
		 *           minutia <code>vautlMinutiae[i]</code> corresponds to the
		 *           vault point <code>(vaultX[i],vaultY[i])</code>.
		 *           <br><br>
		 *           The ordinate binds the genuine minutiae to a secret
		 *           polynomial <i>f</i> over the finite field given by
		 *           <code>gfPtr</code>. More precisely, for a
		 *           genuine vault minutia <code>vaultMinutia[i]</code> we
		 *           have that <code>vaultY[i]=f(vaultX[i])</code> while for
		 *           a chaff minutia <code>vaultMinutia[j]</code> the
		 *           corresponding ordinate is random such that
		 *           <code>vaultY[j]!=f(vaultX[j])</code>.
		 *           <br><br>
		 *           The array is created automatically on enrollment and
		 *           can be accessed via <code>getVaultOrdinates()</code>.
		 *
		 * @see getVaultOrdinates()
		 */
		uint32_t *vaultY;

		/**
		 * @brief
		 *           Pointer to the instance of the finite field over which
		 *           the vault is constructed.
		 *
		 * @details
		 *           The instance of pointed by this member is a binary finite
		 *           field which is isomorphic to the field with 2^16 elements.
		 *           More specifically, the field is equals the field that is
		 *           constructed via <code>SmallBinaryField gf(16)</code>.
		 *           A constant reference can be accessed via
		 *           <code>getField()</code>.
		 *
		 * @see getField()
		 */
		SmallBinaryField *gfPtr;

		/**
		 * @brief
		 *           Contains the SHA-1 hash value of the secret polynomial
		 *           which is generated  and binded to the protected template
		 *           on enrollment.
		 *
		 * @see getHashValue()
		 */
		uint32_t hash[5];

	public:

		/**
		 * @brief
		 *           Constructor which creates an instance of the minutiae
		 *           fuzzy vault to protect minutiae template whose
		 *           coordinates vary within the specified spatial region.
		 *
		 * @param width
		 *           Specifies the range in where a minutia's horizontal
		 *           position should vary.
		 *
		 * @param height
		 *           Specifies the range in where a minutia's vertical
		 *           position should vary.
		 *
		 * @param tryRandom
		 *           If <code>true</code> the instance is advised to use a
		 *           random generator that is cryptographically secure
		 *           (e.g., see <code>MathTools::rand8()</code>).
		 *
		 * @warning
		 *           If <code>width</code> or <code>height</code> is smaller
		 *           than or equals 0 the constructor prints an error message
		 *           to <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 */
		MinutiaeFuzzyVault
		( int width , int height , bool tryRandom = true );

		/**
		 * @brief
		 *           Destructor which frees all memory helt by the minutiae
		 *           fuzzy vault.
		 */
		~MinutiaeFuzzyVault();

		/**
		 * @brief
		 *           Protects the given minutiae template using the instance.
		 *
		 * @details
		 *           1. The minutiae of the view are sorted w.r.t. their
		 *           respective quality.
		 *           <br><br>
		 *           2. Then only those minutiae are selected that keep at least
		 *           \f$\delta_1\f$=<code>getMinimalInterVaultDistance()</code>
		 *           distance to each other minutiae in the view. Thereby the
		 *           distance of two minutiae <code>Minutia minutia1</code>
		 *           and <code>Minutia minutia2</code> is meaured by
		 *           <code>FingerTools::dist(minutia1,minutia2,w)</code> where
		 *           <code>w=getAngleWeight()</code>.
		 *           <br><br>
		 *           3. If not at least
		 *           \f$t_{\min}\f$=<code>getMinimalGenuineSetSize()</code> can
		 *           be selected this way, the function terminates by returning
		 *           <code>false</code> in which case the instance will be left
		 *           unchanged (i.e. if the vault already protects a template
		 *           it remains protected); otherwise, it continues as follows.
		 *           <br><br>
		 *           4. Let
		 *           \f$t\leq t_{\max}\f$=<code>getMaximalGenuineSetSize()</code>
		 *           be the number of selected minutiae. Then <i>n-t</i>
		 *           random chaff minutiae that lay within the specified region
		 *           are generated such that they keep mutual distance
		 *           of \f$\delta_1\f$ to each other as well as to selected
		 *           minutiae.
		 *           <br><br>
		 *           5. The selected genuine minutiae and chaff minutiae are
		 *           put in a list which is sorted w.r.t. lexicographical.
		 *           <br><br>
		 *           6. A polynomial \f$f\in{\bf F}\f$ where \f${\bf F}\f$=
		 *           <code>getField()</code> of degree smaller than <i>k</i>=
		 *           <code>getSecretSize()</code> is generated at random whose
		 *           SHA-1 hash value is saved by the instance.
		 *           <br><br>
		 *           7. As list of vault points
		 *           \f$(x_0,y_0),...,(x_{n-1},y_{n-1})\f$ is generated such that
		 *           \f$x_0,...,x_{n-1}\f$ are elements of \f${\bf F}\f$ that
		 *           are encoded by the integers \f$0,...,n-1\f$, respectively;
		 *           furthermore, if <i>i</i> is the index of a genuine vault
		 *           minutia, the corresponding ordinate will be assigned as
		 *           \f$y_i=f(x_i)\f$; otherwise, if <i>j</i> is the index of a
		 *           chaff minutia \f$y_j\f$ is randomly chosen such that
		 *           \f$y_j\neq f(x_j)\f$. In this way, there is a one-to-one
		 *           correspondence between genuine vault points and genuine
		 *           vault points as well as between chaff points and
		 *           chaff minutiae.
		 *           <br><br>
		 *           8. Finally, the function returns <code>true</code> to
		 *           indicate that the minutiae template has been protected
		 *           by the instance successfully.
		 *
		 * @param view
		 *           Minutiae template that is protected.
		 *
		 * @return
		 *           <code>true</code> if <code>view</code> could be
		 *           successfully protected; otherwise, if no suffciently many
		 *           well-separated minutiae could be selected from
		 *           <code>view</code> the functions returns <code>false</code>.
		 *
		 * @warning
		 *           If the parameters of the instance are not carefully set,
		 *           in particular if
		 *           <code>getMinimalInterVaultDistance()</code> is too large,
		 *           the generation of sufficiently many chaff points may
		 *           become impossible which causes an infinite loop.
		 *
		 * @warning
		 *           If the vault parameters have been initialized such that
		 *           they make no sense, the function prints an error message
		 *           to <code>stderr</code> and exits with status 'EXIT_FAILURE'. This
		 *           will be the case if the relation
		 *           \f$k\leq t_{\min}\leq t_{\max}\leq n\f$ is not satisfied
		 *           where
		 *           <i>k</i>=<code>getSecretSize()</code>,
		 *           \f$t_{\min}\f$=<code>getMinimalGenuineSetSize()</code>,
		 *           \f$t_{\max}\f$=<code>getMaximalGenuineSetSize()</code>,
		 *           and <i>n</i>=<code>getVaultSize()</code>.
		 *
		 * @warning
		 *           If not enough memory could be allocated the functions
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 */
		bool enroll( const MinutiaeView & view );

		/**
		 * @brief
		 *           Attempts to open the vault using an aligned query minutiae
		 *           template.
		 *
		 * @details
		 *           1. In the same way as on enrollment, the function first
		 *           selects at most the first
		 *           \f$t_{\max}\f$=<code>getMaximalGenuineSetSize</code>
		 *           minutiae from the aligned query template that keep a
		 *           mutual distance of \f$\delta_1\f$=
		 *           <code>getMinimalInterVaultDistance()</code> (see
		 *           <code>enroll()</code>).
		 *           <br><br>
		 *           2. For each selected query minutiae its nearest vault
		 *           minutia is determined.
		 *           <br><br>
		 *           3. If the distance is at most \f$\delta_2\f$=
		 *           <code>getMaximalMatchDistance()</code> the corresponding
		 *           vault point is included in the unlocking set.
		 *           <br><br>
		 *           4. Then the function attempts to decoded the secret
		 *           polynomial from the unlocking set. Therefore, the
		 *           unlocking points as well as the SHA-1 hash value of the
		 *           secret polynomial are passed to the <code>decode()</code>
		 *           function.
		 *           <br><br>
		 *           5. If the polynomial could be successfully decoded the
		 *           reference <code>f</code> will be assigned by the decoded
		 *           polynomial and the function returns <code>true</code>;
		 *           otherwise, the function returns <code>false</code>.
		 *
		 * @param f
		 *           If opening the vault is successful, the polynomial will
		 *           be equals to the polynomial that is protected by the
		 *           vault.
		 *
		 * @param view
		 *           The aligned query template used to open the vault.
		 *
		 * @return
		 *           <code>true</code> if opening the vault was successful
		 *           and <code>false</code> otherwise.
		 *
		 * @warning
		 *           If the instance does not protect a template, the function
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If not sufficient memory can be allocated the function
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 */
		bool open
			( SmallBinaryFieldPolynomial & f ,
			  const MinutiaeView & view ) const;

		/**
		 * @brief
		 *           Attempts to decode a polynomial given unlocking points.
		 *
		 * @details
		 *           The implementation of the method is a wrapper around the
		 *           <code>FuzzyVaultTools::bfdecode()</code> function.
		 *           <br><br>
		 *           More precisely, let
		 *           \f${\bf U}=\{(x[0],y[0]),...,(x[t-1],y[t-1])\}\f$ be
		 *           the unlocking set. Then the decoder iterates through all
		 *           possible choices of <i>k</i>=<code>getSecretSize()</code>
		 *           subsets of the unlocking sets. Then its corresponding
		 *           candidate polynomial \f$f^*(X)\f$ of degree smaller
		 *           <i>k</i> is computed which interpolates the points in the
		 *           subset. The SHA-1 hash value of the candidate polynomial
		 *           is then compared with the hash value of the
		 *           correct polynomial; if both are equal, the function
		 *           assigns <i>f</i> with the candidate polynomial and
		 *           returns <code>true</code>; otherwise, if no candidate
		 *           polynomial with the hash of the correct polynomial can be
		 *           found, <i>f</i> will be left unchanged and the function
		 *           returns <code>false</code>.
		 *
		 * @param f
		 *           If decoding was successful, the polynomial will be equals
		 *           the decoded polynomial.
		 *
		 * @param x
		 *           Successive abscissas of the unlocking set.
		 *
		 * @param y
		 *           Successive ordinates of the unlocking set.
		 *
		 * @param t
		 *           Size of the unlocking set.
		 *
		 * @return
		 *           <code>true</code> if decoding was successful and
		 *           <code>false</code> otherwise.
		 *
		 * @warning
		 *           If the instance does not protect a template the function
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If the first <i>t</i> entries of <i>x</i>
		 *           are not pairwise distinct an error message will be
		 *           printed to <code>stderr</code> and an exit
		 *           with status 'EXIT_FAILURE' is caused.
		 *
		 * @warning
		 *           If no sufficient memory could be provided the function
		 *           prints an error message to <code>stderr</code> and
		 *           exits with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If <i>x</i> or if <i>y</i> do not contain at least
		 *           <i>t</i> valid elements in the finite field over which
		 *           <i>f</i> is defined, the function runs into unexpected
		 *           behavior.
		 */
		bool decode
		( SmallBinaryFieldPolynomial & f ,
		  const uint32_t *x , const uint32_t *y ,
		  int t ) const;

		/**
		 * @brief
		 *           Sorts the given minutiae w.r.t. their respective quality
		 *           and dismisses all minutiae that do not keep sufficient
		 *           mutual distance.
		 *
		 * @details
		 *           This function is used in the <code>enroll()</code> and
		 *           <code>open()</code> methods to select minutiae. It works
		 *           as follows.
		 *           <br><br><br>
		 *           1. The minutiae are first sorted w.r.t. their quality. If
		 *           no quality was assigned to the minutiae, i.e. if their
		 *           quality is 0, the relative ordering of the minutiae
		 *           is kept for the sorting algorithm used is stable.
		 *           <br><br>
		 *           2. Afterwards, all minutiae are dismissed that do not
		 *           keep at least
		 *           \f$\delta_1\f$=<code>getMinimalInterVaultDistance()</code>
		 *           distance to the other minutiae. The distance of two
		 *           minutiae will be determined by the <code>dist()</code>
		 *           funtion.
		 *           <br><br>
		 *           3. Finally, the type and the quality attribute of the
		 *           minutiae are set to <code>UNKNOWN</code> and 0,
		 *           respectively, to avoid that genuine minutiae can be
		 *           distinguished from chaff minutiae just because they have
		 *           known minutia type or specified quality.
		 *           <br><br>
		 *           4. The resulting minutiae are returned by the function.
		 *           <br><br>
		 *           <b>Note:</b> for the result we adopt the finger position,
		 *           impression type, and finger quality of the input
		 *           template.
		 *
		 * @param view
		 *           The template from where the minutiae are selected.
		 *
		 * @return
		 *           All minutiae (sorted w.r.t. their respective quality)
		 *           from <code>view</code> that keep sufficient distance to
		 *           all other minutiae in <code>view</code>.
		 */
		MinutiaeView select
		( const MinutiaeView & view ) const;

		/**
		 * @brief
		 *           Measures the distance of two minutiae.
		 *
		 * @details
		 *           The function is used for the <code>select()</code>
		 *           function to select well-separated minutiae. It is a
		 *           wrapper around <code>FingerTools::dist(a,b,w)</code>
		 *           where <code>w=getAngleWeight()</code>.
		 *
		 * @param a
		 *           First minutia.
		 *
		 * @param b
		 *           Second minutia.
		 *
		 * @return
		 *           Distance of <code>a</code> to <code>b</code> by
		 *           accounting for their position and angle.
		 */
		double dist
		( const Minutia & a , const Minutia & b ) const;

		/**
		 * @brief
		 *           Generates a candidate for a chaff minutia.
		 *
		 * @details
		 *           The result will be a random minutia with abscissa
		 *           coordinate in the range 0..<code>getWidth()</code>,
		 *           ordinate coordinate in the range
		 *           0..<code>getHeight()</code>, and random angle. Thereby
		 *           the angle is generated from 256 quantas to avoid that
		 *           chaff minutia can be distinguished from genuine minutia
		 *           just because genuine minutia from ISO templates have
		 *           angle in the quanta but chaff have not. Furthermore,
		 *           the generated minutia is of unknown type and not
		 *           specified quality such as minutiae that are returned
		 *           by the minutiae selection function
		 *           <code>select()</code>.
		 *
		 * @return
		 *           Random chaff minutia.
		 */
		Minutia generateChaffMinutia() const;

		/**
		 * @brief
		 *           Specifies how significant the minutiae angles are taken
		 *           into account when the similarity of two minutiae is
		 *           measured.
		 *
		 * @details
		 *           The specified value will be used by
		 *           <code>dist(const Minutia&,const Minutia&)</code> which is
		 *           a wrapper around
		 *           <code>FingerTools::dist(const Minutia&a,const Minutia&,
		 *           angleWeight)</code>.
		 *           <br><br>
		 *           If not set, the value is assumed as 11.549 per default.
		 *           Note that in Nandakumar et al (2007) the weight was
		 *           chosen as 0.2 but there it is w.r.t. angular measure
		 *           while ours is w.r.t. radian which is equivalent.
		 *
		 * @param angleWeight
		 *           Specifies how significant the minutiae angles are taken
		 *           into account when the similarity of two minutiae is
		 *           measured.
		 *
		 * @warning
		 *           If the instances protects a minutiae template which has
		 *           not been cleared or if <code>angleWeight</code> is
		 *           smaller than or equals 0.0, the method prints an
		 *           error message to <code>stderr</code> and exits with
		 *           status 'EXIT_FAILURE'.
		 *
		 * @see getAngleWeight()
		 */
		void setAngleWeight( double angleWeight );

		/**
		 * @brief
		 *           Specifies the size of the vault.
		 *
		 * @details
		 *           If <i>t</i> minutiae have been selected on enrollment,
		 *           <i>n-t</i> chaff points are generated such that the vault
		 *           reaches a size of <i>n</i>.
		 *           <br><br>
		 *           If not set, the value is assumed as 224 per default.
		 *
		 * @param n
		 *           The specified size the vault has after successful
		 *           enrollment.
		 *
		 * @warning
		 *           If the vault already protects a minutiae template or
		 *           if <i>i</i> is smaller than or equals 0, the method
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'. Note, we can ensure that the vault does
		 *           not protect a template by a call of <code>clear()</code>.
		 *
		 * @see getVaultSize()
		 */
		void setVaultSize( int n );

		/**
		 * @brief
		 *           Specifies how many well-separated genuine minutiae are
		 *           required in a minutiae template to successfully be
		 *           enrolled.
		 *
		 * @details
		 *           If not at least <code>tmin</code> well-separeted minutiae
		 *           can be selected from a minutiae template passed to
		 *           <code>enroll()</code> it will return <code>false</code>.
		 *           <br><br>
		 *           If not set, the value is assumed as 18 per default.
		 *
		 * @param tmin
		 *           Minimal number of well-separated minutiae that are
		 *           required to successfull enroll.
		 *
		 * @warning
		 *           If the vault already protects a minutiae template or
		 *           if <i>tmin</i> is smaller than or equals 0, the method
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'. Note, we can ensure that the vault does
		 *           not protect a template by a call of <code>clear()</code>.
		 *
		 * @see getMinimalGenuineSetSize()
		 */
		void setMinimalGenuineSetSize( int tmin );

		/**
		 * @brief
		 *           Specifies how many well-separated genuine minutiae are
		 *           at most selected from a minutiae template that is
		 *           passed to <code>enroll()</code>
		 *
		 * @details
		 *           On enrollment, only the first <code>tmax</code> that
		 *           are of best quality are selected as genuine minutiae.
		 *           <br><br>
		 *           If not set, the value is assumed as 24 per default.
		 *
		 * @param tmax
		 *           Maximal number of well-separated minutiae that are
		 *           selected on enrollment.
		 *
		 * @warning
		 *           If the vault already protects a minutiae template or
		 *           if <i>tmax</i> is smaller than or equals 0, the method
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'. Note, we can ensure that the vault does
		 *           not protect a template by a call of <code>clear()</code>.
		 *
		 * @see getMaximalGenuineSetSize()
		 */
		void setMaximalGenuineSetSize( int tmax );

		/**
		 * @brief
		 *           Specifies the length of the secret polynomial that is
		 *           randomly generated on enrollment for being binded to the
		 *           genuine minutiae.
		 *
		 * @details
		 *           The polynomial that is generated has <i>k</i>
		 *           coefficients in the field given by
		 *           <code>getField()</code>. Thus, the secret polynomial will
		 *           be of degree strictly smaller than <i>k</i>.
		 *           <br><br>
		 *           If not set, the value is assumed as 9 per default.
		 *
		 * @param k
		 *           The size of the secret polynomial that is generated on
		 *           enrollment.
		 *
		 * @warning
		 *           If the vault already protects a minutiae template or
		 *           if <i>k</i> is smaller than or equals 0, the method
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'. Note, we can ensure that the vault does
		 *           not protect a template by a call of <code>clear()</code>.
		 *
		 * @see getSecretSize()
		 */
		void setSecretSize( int k );

		/**
		 * @brief
		 *           Specifies the minimal distance between vault minutiae.
		 *
		 * @details
		 *           For each two vault minutiae <i>a</i> and <i>b</i> on
		 *           enrollment genuine and vault minutiae are selected and
		 *           generated such that
		 *           <pre>
		 *           dist(a,b) >= interVaultDistance
		 *           </pre>
		 *           <br><br>
		 *           If not set, the value is assumed as 25.0 per default.
		 *
		 * @param interVaultDistance
		 *           Minimal distance between distinct vault minutiae.
		 *
		 * @warning
		 *           If the vault already protects a minutiae template or
		 *           if <code>interVaultDistance</code> is smaller than or
		 *           equals 0, the method prints an error message to
		 *           <code>stderr</code> and exits with status 'EXIT_FAILURE'. Note, we
		 *           can ensure that the vault does not protect a template
		 *           by a call of <code>clear()</code>.
		 *
		 * @see getMinimalInterVaultDistance()
		 */
		void setMinimalInterVaultDistance
			( double interVaultDistance );

		/**
		 * @brief
		 *           Specifies how much two minutiae are allowed to differ
		 *           such that they are considered to match.
		 *
		 * @details
		 *           On authentication (see <code>enroll()</code>) the query
		 *           minutiae are used to extract vault points. Therein, for
		 *           each query minutia its nearest vault minutia is extracted.
		 *           If their distance (by means of <code>dist()</code>) is
		 *           smaller than or equals <code>matchDistance</code> the
		 *           corresponding vault point is included in the unlocking
		 *           set.
		 *           <br><br>
		 *           If not set, the value is assumed as 30.0 per default.
		 *
		 * @param matchDistance
		 *           The maximal distance between two matching minutiae.
		 *
		 * @warning
		 *           If <code>interVaultDistance</code> is smaller than or
		 *           equals 0, the method prints an error message to
		 *           <code>stderr</code> and exits with status 'EXIT_FAILURE'. Note, we
		 *           can ensure that the vault does not protect a template
		 *           by a call of <code>clear()</code>.
		 *
		 * @see getMaximalMatchDistance()
		 */
		void setMaximalMatchDistance
			( double matchDistance );

		/**
		 * @brief
		 *           Ensures that the vault dismisses all data related to
		 *           a protected minutiae template, if any.
		 *
		 * @details
		 *           The method clears all vault minutiae and frees the data
		 *           holding the corresponding vault points. Furthermore, the
		 *           hash value of the secret polynomial is overwritten by
		 *           zero bytes.
		 *           <br><br>
		 *           The methods which are used to specify vault parameters
		 *           require that the vault does not protect a template. If
		 *           it does, with this method we can assert that the vault
		 *           is at the state in where no template is enrolled.
		 */
		void clear();

		/**
		 * @brief
		 *           Access the SHA-1 hash value of the secret polynomial
		 *           that is binded to the genuine vault minutiae.
		 *
		 * @details
		 *           On enrollment, a secret polynomial <i>f</i> is generated
		 *           which will be finally dismissed. But to allow safe
		 *           recovery on genuine authentication, its SHA-1 hash value
		 *           (i.e. the output of <code>f.hashValue(hash)</code>) is
		 *           saved.
		 *
		 * @return
		 *           The secret polynomial's hash value.
		 *
		 * @warning
		 *           If no minutiae template was successfully enrolled, the
		 *           function prints an error message to <code>stderr</code>
		 *           and exits with status 'EXIT_FAILURE',
		 */
		const uint32_t* getHashValue() const;

		/**
		 * @brief
		 *           Access the width specifying the range of a chaff
		 *           minutia's abscissa coordinate.
		 *
		 * @details
		 *           Chaff minutiae are generated such that their
		 *           corresponding abscissa value is between
		 *           0..<code>width-1</code>.
		 *           <br><br>
		 *           The height must be specified on construction of an
		 *           instance of the class.
		 *
		 * @return
		 *           The width specifying the range of a chaff minutia's
		 *           abscissa coordinate.
		 *
		 * @see MinutiaeFuzzyVault(int,int,bool)
		 * @see getHeight()
		 */
		inline int getWidth() const { return this->width; }

		/**
		 * @brief
		 *           Access the height specifying the range of a chaff
		 *           minutia's ordinate coordinate.
		 *
		 * @details
		 *           Chaff minutiae are generated such that their
		 *           corresponding ordinate value is between
		 *           0..<code>height-1</code>.
		 *           <br><br>
		 *           The height must be specified on construction of an
		 *           instance of the class.
		 *
		 * @return
		 *           The height specifying the range of a chaff minutia's
		 *           ordinate coordinate.
		 *
		 * @see MinutiaeFuzzyVault(int,int,bool)
		 * @see getWidth()
		 */
		inline int getHeight() const { return this->height; }

		/**
		 * @brief
		 *           Access the weight that controls how significant the
		 *           minutiae angles are taken into account when the
		 *           similarity of two minutiae is measured.
		 *
		 * @details
		 *           The returned value is used by
		 *           <code>dist(const Minutia&,const Minutia&)</code> which is
		 *           a wrapper around
		 *           <code>FingerTools::dist(const Minutia&a,const Minutia&,
		 *           angleWeight)</code>.
		 *           <br><br>
		 *           Unless set using <code>setAngleWeight(double)</code>,
		 *           the function returns 11.549 per default.
		 *
		 * @return
		 *           Weight controlling the significance of minutiaes angles
		 *           on similarity measurement.
		 *
		 * @see setAngleWeight(double)
		 */
		inline double getAngleWeight() const { return this->angleWeight; }

		/**
		 * @brief
		 *           Returns the size the vault will have after successful
		 *           enrollment.
		 *
		 * @details
		 *           Unless set using <code>setVaultSize(int)</code>, the
		 *           function returns 224 per default.
		 *
		 * @return
		 *           The size of the vault.
		 *
		 * @see setVaultSize(int)
		 */
		inline int getVaultSize() const { return this->n; }

		/**
		 * @brief
		 *           Returns the minimal number of well-separated minutiae
		 *           that are required from a minutiae template to
		 *           successfully enroll.
		 *
		 * @details
		 *           Unless set using
		 *           <code>setMinimalGenuineSetSize(int)</code>, the function
		 *           returns 18 per default.
		 *
		 * @return
		 *           Lower bound on number of genuine vault minutiae.
		 *
		 * @see setMinimalGenuineSetSize(int)
		 */
		inline int getMinimalGenuineSetSize() const { return this->tmin; }

		/**
		 * @brief
		 *           Returns the maximal number of well-separated minutiae
		 *           that are selected from a minutiae template as genuine
		 *           vault minutiae.
		 *
		 * @details
		 *           Unless set using
		 *           <code>setMaximalGenuineSetSize(int)</code>, the function
		 *           returns 24 per default.
		 *
		 * @return
		 *           Upper bound on number of genuine vault minutiae.
		 *
		 * @see setMaximalGenuineSetSize(int)
		 */
		inline int getMaximalGenuineSetSize() const { return this->tmax; }

		/**
		 * @brief
		 *           Access the number of coefficients of the
		 *           secret polynomial which is binded to the genuine vault
		 *           points on enrollment.
		 *
		 * @details
		 *           Unless set using <code>setSecretSize(int)</code>, the
		 *           function returns 9 per default.
		 *
		 * @return
		 *           The size of the secret polynomial that is generated on
		 *           enrollment.
		 *
		 * @see setSecretSize(int)
		 */
		inline int getSecretSize() const { return this->k; }

		/**
		 * @brief
		 *           Access the minimal guaranteed distance between vault
		 *           minutiae.
		 *
		 * @details
		 *           Unless set using
		 *           <code>setMinimalInterVaultDistance(int)</code>, the
		 *           function returns 25.0 per default.
		 *
		 * @return
		 *           Minimal distance between distinct vault minutiae.
		 *
		 * @see setMinimalInterVaultDistance(int)
		 */
		inline double getMinimalInterVaultDistance() const {
			return this->minInterVaultDistance;
		}

		/**
		 * @brief
		 *           Returns the maximal between two matching such that they
		 *           are considered to match.
		 *
		 * @details
		 *           The distance between two minutiae is computed using
		 *           <code>dist()</code>.
		 *           <br><br>
		 *           Unless set using <code>setMaximalMatchDistance()</code>,
		 *           the function returns 30.0 per default.
		 *
		 * @return
		 *           The maximal distance between two matching minutiae.
		 *
		 * @see setMaximalMatchDistance(int)
		 */
		inline double getMaximalMatchDistance() const {
			return this->maxMatchDistance;
		}

		/**
		 * @brief
		 *           Access a constant reference to the finite field over that
		 *           the vault is constructed.
		 *
		 * @details
		 *           The finite field can not be changed in any way and will
		 *           be the finite field that is constructed via
		 *           <code>SmallBinaryField gf(16)</code> which is of size
		 *           2^16=65563. The size of the finite field  allow the vault
		 *           to be of size up to 65536 which is sufficient for the
		 *           protection of fingerprint minutiae templates.
		 *
		 * @return
		 *           The vault's underlying finite field.
		 */
		inline const SmallBinaryField & getField() const {
			return this->gfPtr[0];
		}

		/**
		 * @brief
		 *           Access the vault minutiae.
		 *
		 * @return
		 *           A constant reference to the vector of vault minutiae.
		 *
		 * @see getVaultAbscissas()
		 * @see getVaultOrdinates()
		 */
		inline const std::vector<Minutia> & getVaultMinutiae() const {
			return this->vaultMinutiae;
		}

		/**
		 * @brief
		 *           Access the vault points abscissas.
		 *
		 * @details
		 *           After successful enrollment, a list of vault
		 *           minutiae \f$T\f$ as well as a list of vault
		 *           points \f${\bf V}\subset{\bf F}\times{\bf F}\f$ is
		 *           generated. There is a one-to-one correspondence between
		 *           vault minutiae and vault points. For our
		 *           implementation let
		 *           <pre>
		 *            MinutiaeFuzzyVault vault = ... //Enrolled vault
		 *            vector<Minutia> vaultMinutiae = vault.getVaultMinutiae();
		 *            uint32_t *x, *y;
		 *            x = vault.getVaultAbscissas();
		 *            y = vault.getVaultOrdinates();
		 *           </pre>
		 *           then the vault point that corresponds to the
		 *           <i>i</i>th vault minutiae, i.e.
		 *           <code>vaultMinutiae[i]</code> is the vault point
		 *           <code>(x[i],y[i])</code>. Furthermore, the <i>i</i>th
		 *           vault point's abscissa <code>x[i]</code> is the finite
		 *           field element that is encoded by the integer <i>i</i>.
		 *
		 * @return
		 *           The vault points abscissas.
		 *
		 * @see getVaultMinutiae()
		 * @see getVaultOrdinates()
		 */
		inline const uint32_t *getVaultAbscissas() const {
			return this->vaultX;
		}

		/**
		 * @brief
		 *           Access the vault points ordinates.
		 *
		 * @return
		 *           The vault points ordinates.
		 *
		 * @see getVaultMinutiae()
		 * @see getVaultAbscissas()
		 */
		inline const uint32_t *getVaultOrdinates() const {
			return this->vaultY;
		}
	};
}

#endif /* THIMBLE_MINUTIAEFUZZYVAULT_H_ */
