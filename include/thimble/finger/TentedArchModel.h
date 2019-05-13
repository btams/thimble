/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
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
 * @file TentedArchModel.h
 *
 * @brief
 *            Provides a class that represents the ridge flow of
 *            tented arch type fingerprints and provides a mechanism
 *            for fitting them to fingerprint orientation field
 *            estimations.
 *
 * @section sec_tarp Absolute Fingerprint Pre-Alignment
 *
 * Given a fingerprint encoded as an object of the \link thimble::Fingerprint
 * Fingerprint\endlink class and a minutiae template as an object generated
 * from \link thimble::MinutiaeView MinutiaeView\endlink
 * <pre>
 *     Fingerprint fingerprint; // Minutiae template
 *     MinutiaeView view; // Minutiae template
 * </pre>
 * the easiest way to absolutely pre-aligned the minutiae template
 * is to represent their minutiae with respect to coordinate system
 * derived from the fingerprint's directed reference point estimation.
 * The directed reference point can be estimated as follows
 * <pre>
 *    DirectedPoint drp;
 *
 *    if ( !fingerprint.hasDirectedReferencePoint() ) {
 *       cerr << "Directed reference point estimation failed." << endl;
 *       exit(EXIT_FAILURE);
 *    }
 *
 *    drp = fingerprint.getDirectedReferencePoint();
 * </pre>
 * and the absolutely pre-aligned minutiae can be computed via
 * <pre>
 *    MinutiaeView prealignedView = FingerTools::prealign(view,drp);
 * </pre>
 * The critical step is the estimation of the fingerprint's directed reference
 * point which is implemented by the
 * virtual \link thimble::Fingerprint::estimateDirectedReferencePoint()
 * fingerprint.estimateDirectedReferencePoint()\endlink function around
 * which the functions \link thimble::Fingerprint::hasDirectedReferencePoint()
 * fingerprint.hasDirectedReferencePoint()\endlink
 * and \link thimble::Fingerprint::getDirectedReferencePoint()
 * fingerprint.getDirectedReferencePoint()\endlink are wrapped.
 *
 * Internally,
 * the \link thimble::Fingerprint::estimateDirectedReferencePoint()
 * fingerprint.estimateDirectedReferencePoint()\endlink function uses
 * an object of the \link thimble::TentedArchModel TentedArchModel\endlink
 * class which implements the methods presented in
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
 * @author Benjamin Tams
 */

#ifndef THIMBLE_TENTEDARCHMODEL_H_
#define THIMBLE_TENTEDARCHMODEL_H_

#include <vector>
#include <complex>

#include <thimble/dllcompat.h>
#include <thimble/math/AffineTransform.h>
#include <thimble/image/Orientation.h>
#include <thimble/finger/Fingerprint.h>

#ifdef THIMBLE_BUILD_DLL
template<> class THIMBLE_DLL std::complex<double>;
#endif

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Represents the ridge flow of a tented arch type fingerprint.
	 *
	 * @details
	 *            A tented arch can be constructed via
	 *            <pre>
	 *             %TentedArchModel model;
	 *            </pre>
	 *            with origin (0,0) and direction of longitudinal axis along
	 *            the coordinate system's ordinate. Each tented arch has a
	 *            core and a delta placed on the longitudinal axis whose
	 *            absolute coordinates (as complex numbers <i>x+i*y</i>) can
	 *            be accessed via
	 *            <pre>
	 *             model.\link getComplexCore()\endlink;
	 *            </pre>
	 *            and
	 *            <pre>
	 *             model.\link getComplexDelta()\endlink;
	 *            </pre>
	 *            respectively.
	 *            The representation of the tented arch presented by
	 *            <code>model</code> can be moved by a translation
	 *            \f$z\mapsto z+s\f$ via
	 *            <pre>
	 *             model.\link translate(const std::complex<double>&) translate(s)\endlink;
	 *            </pre>
	 *            To obtain the representation of the model rotated by an angle
	 *            \f$\theta\in[0,\pi)\f$ around the arch's core, we may run
	 *            <pre>
	 *             model.\link rotate(double) rotate(theta)\endlink;
	 *            </pre>
	 *
	 *            <h2>
	 *             Directed Reference Point Estimation using the
	 *             %TentedArchModel class
	 *            </h2>
	 *            Let the orientation field of a fingerprint be given by a
	 *            vector of sample points constituted with an orientation, i.e.,
	 *            <pre>
	 *             vector<\link OrientedPoint\endlink> samplePoints = ...;
	 *            </pre>
	 *            The translation of an initial tented arch
	 *            <pre>
	 *             %TentedArchModel model;
	 *            </pre>
	 *            can be adjusted via
	 *            <pre>
	 *             model.\link fitTranslation() fitTranslation(samplePoints)\endlink;
	 *            </pre>
	 *            such that the ridge flow in its core's neighborhood well
	 *            matches the ridge flow of the orientation field.
	 *            Correspondingly, the tented arch's rotation can be adjusted
	 *            around its \link getComplexCore() core\endlink via
	 *            <pre>
	 *             model.\link fitRotation() fitRotation(samplePoints)\endlink;
	 *            </pre>
	 *            Both the fitting functions may be useful for minimizing the
	 *            cost function given
	 *            <pre>
	 *             double cost = model.\link distance() distance(samplePoints)\endlink
	 *            </pre>
	 *            over the rotation and translation parameters. The adjusted
	 *            tented arche's core
	 *            <pre>
	 *             double x, y;
	 *             x = model.getComplexCore().real();
	 *             y = model.getComplexCore().imag();
	 *            </pre>
	 *            may serve as a reference point and the angle
	 *            <pre>
	 *             double angle = model.\link getReferencePointDirection()\endlink;
	 *            </pre>
	 *            may define its direction defining the origin of the
	 *            abscissa of a directed reference point estimation.
	 *
	 *            <h2>
	 *             Absolute %Fingerprint Pre-Alignment
	 *            </h2>
	 *            Given a directed reference point <code>(x,y,angle)</code> as
	 *            above of a fingerprint's orientation field, the features
	 *            (minutiae, say) can be shifted to the polar coordinate system
	 *            defined by the reference point. For example, if
	 *            <pre>
	 *             MinutiaeView minutiae = ...;
	 *            </pre>
	 *            contain the minutiae of the fingerprint, their presentation
	 *            w.r.t. to the new coordinate system can be computed
	 *            by
	 *            <pre>
	 *             MinutiaeView prealigned = \link FingerTools::prealign() FingerTools::prealign(minutiae,x,y,angle)\endlink;
	 *            </pre>
	 *
	 * @see <ul>
	 *       <li>
	 *        <b>Huckemann, Hotz, Munk (2008).</b>
	 *        <i>
	 *         Global Models for the %Orientation Field of Fingerprints:
	 *         An Approach Based on Quadratic Differentials.
	 *        </i>
	 *        IEEE Trans. Pattern Anal. Machine Intell., vol. 30, no. 9,
	 *        pp. 1507--1519
	 *       </li>
	 *       <li>
	 *        <b>Tams (2013).</b>
	 *        <i>
	 *         Absolute %Fingerprint Pre-Alignment in Minutiae-Based
	 *         Cryptosystems.
	 *        </i>
	 *        in Proc. of BIOSIG, pp. 75--86.
	 *       </li>
	 *       <li>
	 *        <b>Tams, Mihailescu, and Munk (2014).</b>
	 *        <i>
	 *         Security-Improved Minutiae-Based Fuzzy Vault.
	 *        </i>
	 *        in preparation.
	 *       </li>
	 *      </ul>
	 */
	class THIMBLE_DLL TentedArchModel {

	private:
		/**
		 * @brief
		 *            Controls the distance of the rectangular grid at which
		 *            the orientations of a fingerprint's orientation field
		 *            are measured.
		 *
		 * @details
		 *            If the dimension of the fingerprint image is
		 *            <code>width*height</code>, then the points at which
		 *            the orientations are measures are placed on a rectangular
		 *            grid of distance <code>2*halfGridDist+1</code> in which
		 *            the points have a distance to the boundary at of at
		 *            least \link halfGridDist\endlink to allow
		 *            that the orientations can be estimated by the gradient
		 *            method without suffering from boundary effects.
		 *
		 * @see getHalfGridDistance()
		 * @see setHalfGridDistance(int)
		 */
        int halfGridDist;

		/**
		 * @brief
		 *            Essentially, the value of this field is constantly 1.
		 *
		 * @details
		 *            The ridge flow of a tented arch is, essentially, the
		 *            ridge flow of an arch influenced by a core and a delta
		 *            both lying on the longitudinal axis of the arch. An
		 *            arch can ibe modeled as the complex function
		 *            \f[
		 *             \psi(z)=\lambda^2\cdot(z^2-R^2)^{2n},
		 *            \f]
		 *            in general, where <i>n=1</i> is practice.
		 */
		int n;

		/**
		 * @brief
		 *            The stretching factor of the (tented) arch.
		 *
		 * @details
		 *            The field specifies the real parameter \f$\lambda\f$
		 *            for the function
		 *            \f[
		 *             \psi(z)=\lambda^2\cdot(z^2-R^2)^{2n}
		 *            \f]
		 *            defining the tented arch model's underlying
		 *            arch.
		 *
		 * @see getLambda()
		 * @see setLambda(double)
		 */
		double lambda;

		/**
		 * @brief
		 *            Specifies the distance of the poles from arche's
		 *            origina along the coordinate system's abscissa.
		 *
		 * @details
		 *            The field specifies the real parameter \f$R\f$
		 *            for the function
		 *            \f[
		 *             \psi(z)=\lambda^2\cdot(z^2-R^2)^{2n}
		 *            \f]
		 *            defining the tented arch model's underlying
		 *            arch.
		 *
		 * @see getPoleDistance()
		 * @see setPoleDistance(double)
		 */
		double poleDistance;

		/**
		 * @brief
		 *            Defines the distance of the tented arch's core
		 *            from the origin along the longitudinal axis.
		 *
		 * @details
		 *            A tented arch is modeled via the complex function
		 *            \f[
		 *             \tau(z)=\psi(z)\cdot\frac{z^2+d_{core}}{z^2+d_{delta}}
		 *            \f]
		 *            where \f$d_{core}\f$ and \f$d_{delta}\f$ specifies
		 *            the distance of the core and delta, respectively. The
		 *            field <code>coreDistance</code> specifies the real
		 *            parameter \f$d_{core}\f$.
		 *
		 * @see getCoreDistance()
		 * @see setCoreDistance(double)
		 */
		double coreDistance;

		/**
		 * @brief
		 *            Defines the distance of the tented arch's core
		 *            from the origin along the longitudinal axis.
		 *
		 * @details
		 *            A tented arch is modeled via the complex function
		 *            \f[
		 *             \tau(z)=\psi(z)\cdot\frac{z^2+d_{core}}{z^2+d_{delta}}
		 *            \f]
		 *            where \f$d_{core}\f$ and \f$d_{delta}\f$ specifies
		 *            the distance of the core and delta, respectively. The
		 *            field <code>coreDistance</code> specifies the real
		 *            parameter \f$d_{delta}\f$.
		 *
		 * @see getDeltaDistance()
		 * @see setDeltaDistance(double)
		 */
		double deltaDistance;

		/**
		 * @brief
		 *            Specifies the distance from a tented arch's core
		 *            at which a fingerprint orientation field estimations
		 *            are taken into account with the highest weight when
		 *            the tented arch model is fitted to an orientation field
		 *            estimation.
		 *
		 * @details
		 *            Whenever the rotation/translation of the tented
		 *            arch model \f$\tau_{\alpha,\beta}(z)\f$ (with complex
		 *            rotation part \f$\alpha\f$ where \f$|\alpha|=1\f$,
		 *            complex translation part \f$\beta\f$ and complex core
		 *            \f$\gamma_{\alpha,\beta}\f$) is adjusted
		 *            to an orientation field estimation \f$\{(z_j,v_j)\}\f$
		 *            of a fingerprint (where
		 *            \f$v_j=\cos(2\varphi_j)+i\cdot\sin(2\varphi)\f$,
		 *            \f$\varphi_j\in[0,\pi)\f$), the cost function
		 *            \f[
		 *             \kappa(\alpha,\beta) =
		 * \sum_j\exp\left(\frac{|z-\gamma_{\alpha,\beta}|-\rho)^2}{2\cdot\sigma^2}\right)\cdot
		 * \left|\frac{\tau_{\alpha,\beta}(z_j)}{|\tau_{\alpha,\beta}(z_j)|}-v_j\right|^2
		 *            \f]
		 *            is attempted to be minimized over the rotation part
		 *            \f$\alpha\f$ and the translation part \f$\beta\f$.
		 *            <br><br>
		 *            The field <code>rho</code> specifies the real parameter
		 *            \f$\rho\f$ in the cost function
		 *            \f$\kappa(\cdot,\cdot)\f$.
		 *
		 * @see The cost function \f$\kappa(\alpha,\beta)\f$ can be evaluated
		 *      using the function
		 *      \link distance()\endlink.
		 * @see getRho()
		 * @see setRho(double)
		 */
		double rho;

		/**
		 * @brief
		 *            Specifies the standard deviation of a gaussian tapered
		 *            in the to-be-minimized cost function on fitting the
		 *            tented arch to an orientation field estimation.
		 *
		 * @details
		 *            Whenever the rotation/translation of the tented
		 *            arch model \f$\tau_{\alpha,\beta}(z)\f$ (with complex
		 *            rotation part \f$\alpha\f$ where \f$|\alpha|=1\f$,
		 *            complex translation part \f$\beta\f$ and complex core
		 *            \f$\gamma_{\alpha,\beta}\f$) is adjusted
		 *            to an orientation field estimation \f$\{(z_j,v_j)\}\f$
		 *            of a fingerprint (where
		 *            \f$v_j=\cos(2\varphi_j)+i\cdot\sin(2\varphi)\f$,
		 *            \f$\varphi_j\in[0,\pi)\f$), the cost function
		 *            \f[
		 *             \kappa(\alpha,\beta) =
		 * \sum_j\exp\left(\frac{|z-\gamma_{\alpha,\beta}|-\rho)^2}{2\cdot\sigma^2}\right)\cdot
		 * \left|\frac{\tau_{\alpha,\beta}(z_j)}{|\tau_{\alpha,\beta}(z_j)|}-v_j\right|^2
		 *            \f]
		 *            is attempted to be minimized over the rotation part
		 *            \f$\alpha\f$ and the translation part \f$\beta\f$.
		 *            <br><br>
		 *            The field <code>sigma</code> specifies the real parameter
		 *            \f$\sigma\f$ in the cost function
		 *            \f$\kappa(\cdot,\cdot)\f$.
		 *
		 * @see The cost function \f$\kappa(\alpha,\beta)\f$ can be evaluated
		 *      using the function
		 *      \link distance()\endlink.
		 * @see getSigma()
		 * @see setSigma(double)
		 */
		double sigma;

		/**
		 * @brief
		 *            Determines the maximal number of different initial
		 *            tented arch models for which  minimization attempts
		 *            are conducted on finding a good fit of the tented
		 *            arch model to a fingerprint orientation field
		 *            estimation.
		 *
		 * @details
		 *            The value of the field is chosen as 20 per default.
		 *
		 * @see getMaxInitialCandidates()
		 * @see setMaxInitialCandidates(int)
		 */
		int maxInitialCandidates;

		/**
		 * @brief
		 *            Specifies the precision with which
		 *            <i>Jacobian matrices</i> for numerical differentiation
		 *            are approximated.
		 *
		 * @details
		 *            The value of the field is chosen as 1e-6 and it is not
		 *            intended that its value can be changed manually.
		 *
		 * @see getDerivPrecision()
		 */
		double derivPrecision;

		/**
		 * @brief
		 *            Specifies the precision with which the cost function
		 *            is minimized when fitting the tented arch to a
		 *            fingerprint orientation field estimation.
		 *
		 * @details
		 *            The value of the field is chosen as 1e-6 and it is not
		 *            intended that its value can be changed manually.
		 *
		 * @see getMinimizerPrecision()
		 */
		double minimizerPrecision;

		/**
		 * @brief
		 *            Specifies the maximal bound of iterations run by a
		 *            steepest descent method to minimize the cost function
		 *            when the tented arch is fitted to a fingerprint
		 *            orientation field estimation.
		 *
		 * @details
		 *            The value of the field is chosen as 10000 and it is not
		 *            intended that its value can be changed manually.
		 */
		int maxIterations;

		/**
		 * @brief
		 *            Specifies the number of rerunning the rotation and then
		 *            the translation fit when an initial tented arch is
		 *            fitted to a fingerprint orientation field estimation.
		 *
		 * @details
		 *            The value of the field is chosen as 10 and it is not
		 *            intended that its value can be changed manually.
		 */
		int numRefinements;

		/**
		 * @brief
		 *             The rotation part of the tented arch model.
		 *
		 * @details
		 *             A (fitted) tented arch is modeled by the complex
		 *             function
		 *             \f[
		 *   \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *             \f]
		 *             where \f$|\alpha|=1\f$ is the complex rotation
		 *             part and specified by the field <code>alpha</code>.
		 *
		 * @see getAlpha()
		 */
		std::complex<double> alpha;

		/**
		 * @brief
		 *             The translation part of the tented arch model.
		 *
		 * @details
		 *             A (fitted) tented arch is modeled by the complex
		 *             function
		 *             \f[
		 *   \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *             \f]
		 *             where \f$|\alpha|=1\f$ is the complex rotation
		 *             part. The field <code>beta</code> specifies the complex
		 *             translation part \f$\beta\f$ of this tented
		 *             arch model.
		 *
		 * @see getBeta()
		 */
		std::complex<double> beta;

		/**
		 * @brief
		 *             Implements a real functional to be minimized on
		 *             fitting the translation part.
		 *
		 * @details
		 *             The class is declared here but defined internally.
		 *
		 * @see fitTranslation()
		 * @see RealFunctional
		 */
		class TranslationParameterFunction;

		/**
		 * @brief
		 *             Implements a real functional to be minimized on
		 *             fitting the rotation part.
		 *
		 * @details
		 *             The class is declared here but defined internally.
		 *
		 * @see fitRotation()
		 * @see RealFunctional
		 */
		class RotationParameterFunction;

	public:

		/**
		 * @brief
		 *            Creates a tented arch model w.r.t. to the standard
		 *            coordinate system.
		 *
		 * @details
		 *            Some of the private parameters (unless constant)
		 *            controlling the ridge flow of this tented arch are set
		 *            in dependence of the specified parameter
		 *            <code>dpi</code> by scaling them correspondingly.
		 *
		 * @param dpi
		 *            Specifies the resolution in dots per inch of the
		 *            fingerprints to which this tented arch relates.
		 *
		 * @warning
		 *            If <code>dpi</code> is not greater than zero, an
		 *            error message is written to <code>stderr</code> and
		 *            the program exits with status 'EXIT_FAILURE'.
		 */
		TentedArchModel( int dpi = 500 );

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @param model
		 *            The tented arch model of which is copy is created.
		 */
		TentedArchModel( const TentedArchModel & model );

		/**
		 * @brief
		 *            Assignment operator.
		 *
		 * @param model
		 *            The tented arch model of which this tented arch
		 *            is assigned a copy.
		 *
		 * @return
		 *            A reference to this tented arch model.
		 */
		TentedArchModel & operator=( const TentedArchModel & model );

		/**
		 * @brief
		 *            Assignment operator (procedural version).
		 *
		 * @param model
		 *            The tented arch model of which this tented arch
		 *            is assigned a copy.
		 */
		void assign( const TentedArchModel & model );

		/**
		 * @brief
		 *            Access the reference point estimation after the
		 *            tented arch has been fitted to an orientation field.
		 *
		 * @details
		 *            This method is a convenience method returning the result
		 *            of \link getComplexCore()\endlink whose real part
		 *            is considered as the reference point's abscissa and its
		 *            imaginary part as its ordinate coordinate.
		 *
		 * @return
		 *            \link getComplexCore()\endlink
		 */
		std::complex<double> getComplexReferencePoint() const;

		/**
		 * @brief
		 *            Access the direction of the directed reference point
		 *            estimation given by this tented arch model after it has
		 *            been successfully fitted to an orientation field.
		 *
		 * @return
		 *            An angle in radian (i.e., value in the interval
		 *            \f$[0,2\pi)\f$) specifying the direction.
		 */
		double getReferencePointDirection() const;

		/**
		 * @brief
		 *            Access the affine transform that aligns fingerprint
		 *            features (minutiae, say) to the coordinate system
		 *            specified by the directed reference point of this
		 *            tented arch model.
		 *
		 * @details
		 *            If the core of this tented arch is placed at
		 *            \f$(x_0,y_0)\f$ and \f$\theta\in[0,2\pi)\f$ is the angle
		 *            returned by \link getReferencePointDirection()\endlink,
		 *            then the returned affine transform represents the
		 *            isometry
		 *            \f[
		 *            \left(\begin{array}{c}x\\y\end{array}\right)\mapsto
		 *            \left(\begin{array}{cc}
		 *             \cos(\theta)&\sin(\theta)\\
		 *            -\sin(\theta)&\cos(\theta)
		 *                  \end{array}\right)\cdot
		 *            \left(\begin{array}{c}
		 *             x-x_0\\
		 *             y-y_0
		 *            \end{array}\right)
		 *            \f]
		 *            such that it maps the coordinate system encoded
		 *            by the directed reference point to the standard
		 *            coordinate system.
		 *
		 * @return
		 *            An affine transform defining the prealigning isometry.
		 */
		AffineTransform getPrealigningTransform() const;

		/**
		 * @brief
		 *            Access the origin of this tented arch's coordinate
		 *            system.
		 *
		 * @details
		 *            This function
		 *            returns the complex value
		 *            \f[
		 *             \omega=-(\beta\cdot\alpha^{-1})
		 *            \f]
		 *            which represents the origin of this tented arch's
		 *            coordinate system. For more details on the the
		 *            tented arch's parameters, we refer to the
		 *            documentation of the \link eval()\endlink function.
		 *
		 * @return
		 *            The origin of this tented arch's coordinate system
		 *            encoded by a complex number.
		 *
		 * @attention
		 *            The origin being return by this function should not
		 *            be confused with the origin of the intrinsic coordinate
		 *            system defined by the directed reference point of which
		 *            origin serves the (by convention) core which can be
		 *            accessed via \link getComplexCore()\endlink and
		 *            \link getComplexReferencePoint()\endlink.
 		 */
		std::complex<double> getComplexOrigin() const;

		/**
		 * @brief
		 *            Computes the absolute position of this tented arch
		 *            model's core as a complex number.
		 *
		 * @details
		 *            For details on the parameters that control the tented
		 *            arch, we refer to the documentation of the
		 *            \link eval()\endlink function.
		 *
		 * @return
		 *            \f$
		 *             \gamma_{\alpha,\beta}=-\alpha^{-1}\cdot
		 *             (\beta+i\cdot d_{core})
		 *            \f$.
		 *
		 * @attention
		 *            Note that in Tams, Mihailescu, Munk (2014), the position
		 *            of the complex core equals
		 *            \f$\alpha^{-1}\cdot(i\cdot d_{core}-\beta)\f$ in a
		 *            mathematical coordinate system. In an image coordinate
		 *            system, however, where the origin is the left upper
		 *            corner, some signs require to be switched.
		 *
		 * @see eval()
		 */
		std::complex<double> getComplexCore() const;

		/**
		 * @brief
		 *             Computes the absolute position of this tented arch
		 *             model's delta as a complex number.
		 *
		 * @details
		 *             For details the parameters that control the tented
		 *             arch, we refere to the documentation of the
		 *             \link eval()\endlink function.
		 *
		 * @return
		 *            \f$
		 *             \gamma_{\alpha,\beta}=-\alpha^{-1}\cdot
		 *             (\beta+i\cdot d_{delta})
		 *            \f$.
		 *
		 * @see eval()
		 */
		std::complex<double> getComplexDelta() const;

		/**
		 * @brief
		 *            Access the rotation part of the tented arch model.
		 *
		 * @details
		 *            The tented arch model is described by a complex function
		 *            of the form
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            with \f$|\alpha|=1\f$. The function returns the complex
		 *            \f$\alpha\f$.
		 *
		 * @return
		 *            \f$\alpha\f$ with \f$|\alpha|=1\f$.
		 */
		std::complex<double> getAlpha() const;

		/**
		 * @brief
		 *            Access the translation part of the tented arch model.
		 *
		 * @details
		 *            The tented arch model is described by a complex function
		 *            of the form
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            with \f$|\alpha|=1\f$. The function returns the complex
		 *            \f$\beta\f$.
		 *
		 * @return
		 *            The translation part \f$\beta\f$ a complex number.
		 */
		std::complex<double> getBeta() const;

		/**
		 * @brief
		 *            Changes the translation part of the tented arch model
		 *            such that the specified complex becomes the model
		 *            coordinate system's origin.
		 *
		 * @details
		 *            The parameter \f$\beta\f$ for the complex function
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            is replaced by
		 *            \f[
		 *             -\alpha^{-1}\cdot\omega
		 *            \f]
		 *            where \f$\omega\f$ denotes the specified complex.
		 *
		 * @param origin
		 *            The new origin \f$\omega\f$ of this model's coordinate
		 *            system.
		 */
		void setOrigin( const std::complex<double> & origin );

		/**
		 * @brief
		 *            Moves the tented arch by the specified translation.
		 *
		 * @details
		 *            The method moves the tented arch by the
		 *            translation \f$z\mapsto z+s\f$. More specifically,
		 *            the parameter \f$\beta\f$ in the function
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            is replaced by
		 *            \f[
		 *             \beta'=\beta+s.
		 *            \f]
		 *
		 * @param s
		 *            The translation by which the tented arch is moved.
		 */
		void translate( const std::complex<double> & s );

		/**
		 * @brief
		 *            Rotates the tented arch around the core by the
		 *            specified angle.
		 *
		 * @details
		 *            Updates rotation part and translation part such that
		 *            the tented arch's core remains fixed. Specifically,
		 *            the parameters \f$\alpha\f$ and \f$\beta\f$ for the
		 *            function
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            are replaced by
		 *            \f[
		 *             \alpha'=r_\theta\cdot\alpha
		 *            \f]
		 *            and
		 *            \f[
		 * \beta'=r_\theta\cdot(i\cdot d_{core}+\beta)-i\cdot d_{core},
		 *            \f]
		 *            respectively, where
		 *            \f$d_{core}=\f$\link getCoreDistance()\endlink and
		 *            \f$r_\theta=\cos(\theta)+i\cdot\sin(\theta)\f$ denotes
		 *            the complex describing the rotation by the angle
		 *            \f$\theta\f$.
		 *
		 * @param theta
		 *            Specifies the rotation angle \f$\theta\in[0,2\pi)\f$.
		 */
		void rotate( double theta );

		/**
		 * @brief
		 *            Replaces this model by a tented arch that well
		 *            approximates the orientation field of a fingerprint.
		 *
		 * @details
		 *            Essentially, this procedure implements the entire
		 *            tented arch fitting approach as considered in
		 *            <ul>
		 *             <li>
		 *              <b>Tams (2013).</b>
		 *              <i>
		 *               Absolute %Fingerprint Pre-Alignment in Minutiae-Based
		 *               Cryptosystems.
		 *              </i>
		 *              in Proc. of BIOSIG, pp. 75--86.
		 *             </li>
		 *             <li>
		 *              <b>Tams, Mihailescu, Munk (2014).</b>
		 *              <i>
		 *               Security-Improved Minutiae-Based Fuzzy Vault.
		 *              </i>
		 *              in preparation.
		 *             </li>
		 *            </ul>
		 *            More specifically, if <i>h =</i>
		 * \link Fingerprint::getGradientOrientationDistance()
		 * fingerprint.getGradientOrientationDistance()\endlink,
		 *            <i>g = 2 * </i>\link getHalfGridDistance()\endlink + 1,
		 *            <i>M =</i>
		 * \link Fingerprint::getHeight() fingerprint.getHeight()\endlink,
		 *            and <i>N =</i>
		 * \link Fingerprint::getWidth() fingerprint.getWidth()\endlink, then
		 *            the orientation field is estimated at each pixel
		 *            <i>(x,y)</i> where <i>x=h+j*g, y=h+j'*g</i>, and
		 *            <i>j,j'=0,1,2,3...</i> such that <i>h<=x<N-h</i> and
		 *            <i>h<=y<M<h</i>. Furthermore, the orientations are
		 *            measured using the
		 * \link Fingerprint::getOrientationAt(int,int)
		 * fingerprint.getOrientationAt(y,x)\endlink function, thereby
		 *            building a vector containing
		 *            \link OrientedPoint OrientedPoints\endlink defining
		 *            an orientation field estimation.
		 *            <br><br>
		 *            Furthermore, a grid on the fingerprint's foreground is
		 *            build that is arranged concurrently to the grid at which
		 *            orientations have been estimated which are used as as
		 *            candidate list of initial tented arch's cores. More
		 *            specifically, a vector of pairs containing the
		 *            coordinates <i>(x,y)</i> where <i>x=h+(g-1)/2+j*g</i>,
		 *            <i>y=h+(g-1)/2+j'*g</i>, and <i>j,j'=0,1,2,3,...</i>
		 *            such that <i>x<N</i> and <i>y<M</i>.
		 *            <br><br>
		 *            The vectors and the fingerprint's foreground image are
		 *            passed to the lower level
		 * \link TentedArchModel::estimate(const std::vector<OrientedPoint>&,const bool*,int,int,const std::vector< std::pair<double,double> >&)
		 * estimate() \endlink
		 *            procedure whose result is returned, thereby searching
		 *            for an initial model, fitting its rotation and
		 *            translation iteratively, and, if failed, repeating
		 *            the fit up to 20 times.
		 *
		 * @param fingerprint
		 *            The fingerprint to which orientation field the tented
		 *            arch model is fitted.
		 *
		 * @return
		 *            <code>true</code> if the tented arch has been
		 *            successfully fitted; otherwise, if the fit yielded no
		 *            tented arch having a core on the fingerprint's
		 *            foreground, the function returns <code>false</code>.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool estimate( Fingerprint & fingerprint );

		/**
		 * @brief
		 *            Replaces this model by a tented arch that well
		 *            approximates the orientation field of a fingerprint
		 *            (low-level function)
		 *
		 * @details
		 *            This function
		 *            <ol>
		 *             <li>
		 *              selects a tented arch (0,1) longitudinal axis
		 *              having an element from <code>grid</code> as core
		 *              for which the \link distance()\endlink to
		 *              <code>samples</code> is minimal;
		 *             </li>
		 *             <li>
		 *              updates the rotation via \link fitRotation()\endlink.
		 *             </li>
		 *             <li>
		 *              updates the translation via
		 *              \link fitTranslation()\endlink;
		 *             </li>
		 *             <li>
		 *              repeats steps 2 and 3 until convergence.
		 *             </li>
		 *            </ol>
		 *            If successfully converged, this tented arch is replaced
		 *            by the adjusted tented arch and the function returns
		 *            <code>true</code>. Otherwise, if looping over steps 2
		 *            and 3 do not cause the tented arch selected in Step 1 to
		 *            converge, the second minimal core is selected in Step 1
		 *            and it is looked up whether looping over steps 2 and 3
		 *            converges. If true, this instance is replaced by the
		 *            content of the secondly tested tented arch; otherwise,
		 *            if false, another candidate (at
		 *            most \link getMaxInitialCandidates()\endlink;
		 *            20 per default), is tried to be adjusted to
		 *            the orientation field. If none of these candidates yield
		 *            a tented arch of which core lays on the fingerprint
		 *            foreground, the content of this instance is left
		 *            unchanged and the function returns <code>false</code>.
		 *
		 * @param samples
		 *            Contains the orientations and positions which build
		 *            the orientation field.
		 *
		 * @param foregroundImage
		 *            A two-dimensional image of height <i>m</i> and width
		 *            <i>n</i> where the pixel at <i>(x,y)</i> is encoded
		 *            at <code>foregroundImage[y*n+x]</code> where the value
		 *            <code>true</code> indicates that the pixel <i>(x,y)</i>
		 *            belongs to the foreground and, otherwise, to the
		 *            background.
		 *
		 * @param m
		 *            The height of the fingerprint image.
		 *
		 * @param n
		 *            The width of the fingerprint image.
		 *
		 * @param grid
		 *            The grid from which the cores for the candidates of
		 *            the initial tented arch model are chosen.
		 *
		 * @return
		 *            <code>true</code> if the tented arch has been
		 *            successfully fitted; otherwise, if the attempt resulted
		 *            in no tented arch having a core on the fingerprint's
		 *            foreground, the function returns <code>false</code>.
		 *
		 * @warning
		 *            If <code>m</code> or <code>n</code> is smaller than or
		 *            equals 0, an error message is printed to
		 *            <code>stderr</code> and the program exits with
		 *            status 'EXIT_FAILURE'.
		 */
		bool estimate
		( const std::vector<OrientedPoint> & samples ,
		  const bool *foregroundImage , int m , int n ,
		  const std::vector< std::pair<double,double> > & grid );

		/**
		 * @brief
		 *            Fits this tented arch model's translation to well
		 *            fit the specified orientation field around the model's
		 *            core.
		 *
		 * @details
		 *            Changes over the
		 *            \link getBeta() translation part\endlink of this model
		 *            may change the value of
		 *            \link distance() distance(samples)\endlink.
		 *            This can be modeled as a
		 *            \link RealFunctional two-dimensional functional\endlink
		 *            of which a local minimum corresponds to a fit of this
		 *            model by adding a suitable translation part. Using a
		 *            \link SteepestDescent steepest-descent method\endlink
		 *            this function tries to find such a minimum and, if
		 *            successfully converged, replaces this model by a
		 *            corresponding fit and returns <code>true</code>;
		 *            otherwise, if minimizing the functional, did not
		 *            converge, the properties of this model is left unchanged
		 *            and this function returns <code>false</code>.
		 *
		 *
		 * @param samples
		 *            A vector of positions attached with orientations
		 *            building the orientation field.
		 *
		 * @return
		 *            <code>true</code> if this model converged and,
		 *            otherwise, <code>false</code>.
		 */
		bool fitTranslation
		( const std::vector<OrientedPoint> & samples );

		/**
		 * @brief
		 *           Fits this tented arch model's rotation to well
		 *           fit the specified orientation field around the model's
		 *           core.
		 *
		 * @details
		 *           Small rotations of this model around its
		 *           \link getComplexCore() core\endlink causes small changes
		 *           of \link distance() distance(samples)\endlink.
		 *           This can be modeled as a
		 *           \link RealFunctional one-dimensional functional\endlink
		 *           over the rotation angle of which a local minimum
		 *           corresponds to a fit of this model by a suitable rotation
		 *           around its \link getComplexCore() core\endlink.
		 *           Using a
		 *           \link SteepestDescent steepest-descent method\endlink
		 *           this function tries to find such a minimum and, if
		 *           successfully converged, replaces this model by a
		 *           corresponding fit and returns <code>true</code>;
		 *           otherwise, if minimizing the functional, did not
		 *           converge, the properties of this model is left unchanged
		 *           and this function returns <code>false</code>.
		 *
		 * @param samples
		 *            A vector of positions attached with orientations
		 *            building the orientation field.
		 *
		 * @return
		 *            <code>true</code> if this model converged and,
		 *            otherwise, <code>false</code>.
		 */
		bool fitRotation
		( const std::vector<OrientedPoint> & samples );

		/**
		 * @brief
		 *            Computes the (weighted) distance of this tented arch
		 *            model to the specified orientation field estimation.
		 *
		 * @details
		 *            Let
		 *            \f[
		 *    \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            be the tented arch model represented by this instance
		 *            with \link getBeta() translation part\endlink
		 *            \f$\beta\f$ and \link getAlpha() rotation part\endlink
		 *            \f$\alpha\f$ (see \link eval()\endlink for more
		 *            details). Furthermore, let \f$\{(z_j,v_j)\}\f$ denote
		 *            the orientation field estimation passed through this
		 *            function via the <code>samples</code> argument: Here
		 *            \f$v_j=\cos(2\varphi_j)+i\cdot\sin(2\varphi_j)\f$,
		 *            \f$\phi_j\in[0,\pi)\f$, denotes the orientation attached
		 *            to the complex position \f$z_j\f$. This function returns
		 *            \f[
		 *             \kappa(\alpha,\beta):=\sum_j w(z_j)\cdot\left|
		 *             \frac{\tau_{\alpha,\beta}(z_j)}
		 *             {|\tau_{\alpha,\beta}(z_j)|}-v_j\right|^2
		 *            \f]
		 *            where the
		 *            \f$w(z_j)\f$=\link weight() weight(\f$z_j\f$)\endlink.
		 *            <br><br>
		 *            This function is used by the functions
		 *            \link fitTranslation()\endlink and
		 *            \link fitRotation()\endlink to minimize its value
		 *            over the translation part and rotation angle,
		 *            respectively.
		 *
		 * @param samples
		 *            Defines the orientation field as described in the
		 *            details where \f$z_j\f$=samples[j].x+i*samples[j].y
		 *            and \f$\varphi_j\f$=
		 *            \link Orientation::getOrientationAngle()
		 *             samples[j].orientation.getOrientationAngle()
		 *            \endlink.
		 *
		 * @return
		 *            \f$\kappa(\alpha,\beta)\f$
		 */
		double distance( const std::vector<OrientedPoint> & samples ) const;

		/**
		 * @brief
		 *            Computes the (weighted) distance of this tented arch
		 *            model to the specified orientation field estimation
		 *            (quick version).
		 *
		 * @details
		 *            This function returns the same as the function
		 *            \link distance()\endlink except that those positions
		 *            from the orientation field estimation are ignored
		 *            for which the weight with which they go into the sum
		 *            \f$\kappa(\alpha,\beta)\f$ is very small. More
		 *            precisely, only those positions being closer than
		 *            \f$2\sigma+\rho\f$
		 *            (where \f$\sigma=\f$\link getSigma()\endlink and
		 *             \f$\rho\f$=\link getRho()\endlink) to this model's
		 *            \link getComplexCore() core\endlink are taken into
		 *            accounting assuming that the others are close to zero.
		 *            This yields to a faster evaluation of the distance.
		 *            Yet, being less accurate.
		 *
		 * @param samples
		 *            see \link distance()\endlink
		 *
		 * @return
		 *            see \link distance()\endlink
		 */
        	double qdistance( const std::vector<OrientedPoint> & samples ) const;

		/**
		 * @brief
		 *            Returns the weight with which an orientation at the
		 *            specified coordinate is taken into account.
		 *
		 * @details
		 *            On fitting this model to an orientation field estimation,
		 *            e.g., by calling the methods \link estimate()\endlink,
		 *            \link fitTranslation()\endlink or
		 *            \link fitRotation()\endlink, the function
		 *            \link distance()\endlink is minimized. Those orientations
		 *            estimated being close to a neighborhood around this
		 *            tented arch's core (\link getComplexCore()\endlink) are
		 *            taken into account with the highest weight.
		 *            Specifically, the weight with which an orientation
		 *            measured at <i>z</i> is taken into account equals
		 *            \f[
		 *          \exp\left(\frac{(|z-\gamma|-\rho)^2}{2\cdot\sigma}\right)
		 *            \f]
		 *            being returned by this function where
		 *            \f$\gamma\f$=\link getComplexCore()\endlink,
		 *            \f$\sigma\f$=\link getSigma()\endlink and
		 *            \f$\rho\f$=\link getRho()\endlink.
		 * @param z
		 *            The coordinate (as a complex number) of which
		 *            orientation's weight is returned by this function.
		 *
		 * @return
		 *            The weight with which this orientation at the specified
		 *            coordinate is taken into account.
		 */
		double weight( const std::complex<double> & z ) const;

		/**
		 * @brief
		 *            Returns the weight with which an orientation at the
		 *            specified coordinate is taken into account.
		 *
		 * @details
		 *            Essentially, this function performs as
		 *            the \link weight()\endlink function except that it
		 *            requires the complex
		 *            core \link getComplexCore()\endlink to be
		 *            specified manually. Invoking the
		 *            function \link getComplexCore()\endlink requires one
		 *            complex addition and multiplication of which repetition
		 *            can be saved.
		 *
		 * @param z
		 *            The coordinate (as a complex number) of which
		 *            orientation's weight is returned by this function.
		 *
		 * @param gamma
		 *            The coordinate of this tented arch's core (as a complex
		 *            number) specified manually.
		 *
		 * @return
		 *            The weight with which this orientation at the specified
		 *            coordinate is taken into account.
		 */
		double weight
		( const std::complex<double> & z ,
		  const std::complex<double> & gamma ) const;

		/**
		 * @brief
		 *            Evaluation of the complex function defining this tented
		 *            arch model.
		 *
		 * @details
		 *            The complex function defining this tented arch
		 *            \f[
		 * \tau_{\alpha,\beta}(z)=\alpha^{-2}\cdot\tau(\alpha\cdot z+\beta)
		 *            \f]
		 *            where
		 *           \f[
		 * \tau(z)=\psi(z)\cdot\frac{z^2+d_{core}^2}{z^2+d_{delta}^2}
		 *           \f]
		 *           is the tented arch w.r.t. the standard coordinate system
		 *           and
		 *           \f[
		 * \psi(z)=\lambda^2\cdot(z^2-R^2)^2
		 *           \f]
		 *          defines the underlying arch.
		 *          <br><br>
		 *          The involved parameters are
		 *          <ul>
 		 *           <li>
 		 *            \f$\alpha=\f$\link getAlpha()\endlink
 		 *           </li>
		 *           <li>
		 *            \f$\beta=\f$\link getBeta()\endlink
		 *           </li>
		 *           <li>
		 *            \f$d_{core}=\f$\link getCoreDistance()\endlink
		 *           </li>
		 *           <li>
		 *            \f$d_{delta}=\f$\link getDeltaDistance()\endlink
		 *           </li>
		 *           <li>
		 *            \f$\lambda=\f$\link getLambda()\endlink
		 *           </li>
		 *           <li>
		 *            \f$R=\f$\link getPoleDistance()\endlink
		 *           </li>
		 *          </ul>
		 *
		 * @param z
		 *          The complex at which \f$\tau_{\alpha,\beta}(\cdot)\f$
		 *          is evaluated.
		 *
		 * @return
		 *          The evaluation of
		 *          \f$\tau_{\alpha,\beta}(\cdot)\f$ at <i>z</i>.
		 */
		std::complex<double> eval( const std::complex<double> & z ) const;

		/**
		 * @brief
		 *          Accesses the precision with which gradients
		 *          are approximated.
		 *
		 * @details
		 *          The value of this field may be 1e-6 and is used for
		 *          computing gradients while minimizing functionals
		 *          during rotation and translation part fitting.
		 *
		 * @return
		 *          Precision with which gradients are approximated.
		 *
		 * @see fitTranslation()
		 * @see fitRotation()
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; the corresponding private member is
		 *          accessed directly when the fitting methods are run and is
		 *          not intended to be changed manually.
		 */
		double getDerivPrecision() const;

		/**
		 * @brief
		 *          Accesses the precision with which cost function
		 *          are minimized during rotation and translation part
		 *          fitting.
		 *
		 * @details
		 *          The value of this field will be 1e-6.
		 *
		 * @return
		 *          Precision with which cost functions are minimized.
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; the corresponding private member is
		 *          accessed directly when the fitting methods are run and is
		 *          not intended to be changed manually.
		 */
		double getMinimizerPrecision() const;

		/**
		 * @brief
		 *          Access the distance of the rectangular grid at which the
		 *          orientations of a fingerprint's orientation field are
		 *          measured when this model is fitted to a fingerprint image.
		 *
		 * @details
		 *          If the dimension of the fingerprint image is
		 *          <code>width*height</code>, then the points at which
		 *          the orientations are measures are placed on a rectangular
		 *          grid of distance <code>2*halfGridDist+1</code> in which
		 *          the points have a distance to the boundary at of at
		 *          least <code>halfGridDist</code> to allow
		 *          that the orientations can be estimated by the gradient
		 *          method without suffering from boundary effects.
		 *
		 * @return
		 *          The distance of the rectangular grid at which the
		 *          orientations of a fingerprint's orientation field are
		 *          measured when this model is fitted to a fingerprint image.
		 *
		 * @see setHalfGridDistance(int)
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; the corresponding private member is
		 *          accessed directly when the fitting methods are run.
		 */
		int getHalfGridDistance() const;

		/**
		 * @brief
		 *            Returns the maximal number of different initial models
		 *            for which fitting attempts are performed on adjusting
		 *            this model to a fingerprint's orientation field
		 *            estimation.
		 *
		 * @details
		 *            The return value will be 20 per default and can be changed
		 *            via the \link setMaxInitialCandidates(int)\endlink method.
		 *
		 * @return
		 *            see brief description.
		 *
		 * @attention
		 *            Do not override this function in order to change the
		 *            behavior of a fit; the corresponding private member is
		 *            accessed directly when the fitting methods are run.
		 */
        	int getMaxInitialCandidates() const;

		/**
		 * @brief
		 *          Sets the distance of the rectangular grid at which the
		 *          orientations of a fingerprint's orientation field are
		 *          measured when this model is fitted to a fingerprint image
		 *          to the specified value.
		 *
		 * @details
		 *          For more details of the meaning of this value, we refer
		 *          to the documentation of the function
		 *          \link getHalfGridDistance()\endlink.
		 *
		 * @param halfGridDist
		 *          The distance of the rectangular grid at which the
		 *          orientations of a fingerprint's orientation field are
		 *          measured when this model is fitted to a fingerprint image
		 *          to the specified value.
		 *
		 * @see getHalfGridDistance()
		 *
		 * @warning
		 *          If <code>halfGridDist</code> is negative (the value 0
		 *          is explicitly allowed), an error message is printed to
		 *          <code>stderr</code> and the program exits with status
		 *          'EXIT_FAILURE'.
		 */
		void setHalfGridDistance( int halfGridDist );

		/**
		 * @brief
		 *          Sets the maximal number of different initial models
		 *          for which fitting attempts are performed on adjusting
		 *          this model to a fingerprint's orientation field
		 *          estimation.
		 *
		 * @param maxInitialCandidates
		 *          The specified maximal number of different initial
		 *          models.
		 *
		 * @warning
		 *          If <code>maxInitialCandidates</code> is smaller than
		 *          1, an error message is printed to <code>stderr</code>
		 *          and the program exits with status 'EXIT_FAILURE'.
		 */
		void setMaxInitialCandidates( int maxInitialCandidates );

		/**
		 * @brief
		 *          Access the stretching factor of this tented arch
		 *          model.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function returns the value of \f$\lambda\f$.
		 *
		 * @return
		 *          The stretching factor of this tented arch model.
		 *
		 * @see setLambda(double)
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; use \link setLambda(double)\endlink
		 *          instead. The corresponding private member is
		 *          accessed directly when methods for fitting this model
		 *          are called.
		 */
		double getLambda() const;

		/**
		 * @brief
		 *          Access the distance of the poles of this tented arch
		 *          model.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function returns the value of \f$R\f$.
		 *
		 * @return
		 *          The pole distance of this tented arch model.
		 *
		 * @see setPoleDistance(double)
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; use
		 *          \link setPoleDistance(double)\endlink instead. The
		 *          corresponding private member is accessed directly when
		 *          methods for fitting this model are called.
		 */
		double getPoleDistance() const;

		/**
		 * @brief
		 *          Access the position of this tented arch's delta
		 *          on its longitudinal axis.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function returns the value of \f$d_{delta}\f$.
		 *
		 * @return
		 *          The position of this tented arch's delta
		 *          on its longitudinal axis.
		 *
		 * @see setDeltaDistance(double)
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; use
		 *          \link setDeltaDistance(double)\endlink instead.
		 *          The corresponding private member is accessed directly
		 *          when methods for fitting this model are called.
		 */
		double getDeltaDistance() const;

		/**
		 * @brief
		 *          Access the position of this tented arch's core
		 *          on its longitudinal axis.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function returns the value of \f$d_{core}\f$.
		 *
		 * @return
		 *          The position of this tented arch's core
		 *          on its longitudinal axis.
		 *
		 * @see setCoreDistance(double)
		 */
		double getCoreDistance() const;

		/**
		 * @brief
		 *          Access the distance from this tented arch's core
		 *          at which orientations are accounted with the highest
		 *          weight.
		 *
		 * @details
		 *          An orientation measured at a complex <i>z</i> is accounted
		 *          with the weight
		 *          \f[
		 *          \exp\left(\frac{(|z-\gamma|-\rho)^2}{2\cdot\sigma}\right)
		 *          \f]
		 *          when the model is fitted to an orientation field
		 *          estimation (see \link weight()\endlink). This function
		 *          returns the value for \f$\rho\f$ in the above formula.
		 *
		 * @return
		 *          The distance from this tented arch's core at which
		 *          orientations are accounted with the highest weight.
		 *
		 * @see setRho(double)
		 * @see weight()
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; use \link setRho(double)\endlink
		 *          instead. The corresponding private member is
		 *          accessed directly when methods for fitting this model
		 *          are called.
		 */
		double getRho() const;

		/**
		 * @brief
		 *          Access the Gaussian standard deviation used to compute the
		 *          weight with which orientations are accounted on fitting the
		 *          model to an orientation field estimation.
		 *
		 * @details
		 *          An orientation measured at a complex <i>z</i> is accounted
		 *          with the weight
		 *          \f[
		 *          \exp\left(\frac{(|z-\gamma|-\rho)^2}{2\cdot\sigma}\right)
		 *          \f]
		 *          when the model is fitted to an orientation field
		 *          estimation (see \link weight()\endlink). This function
		 *          returns the value for \f$\sigma\f$ in the above formula.
		 *
		 * @return
		 *          The distance from this tented arch's core at which
		 *          orientations are accounted with the highest weight.
		 *
		 * @see setSigma(double)
		 * @see weight()
		 *
		 * @attention
		 *          Do not override this function in order to change the
		 *          behavior of a fit; use \link setSigma(double)\endlink
		 *          instead. The corresponding private member is
		 *          accessed directly when methods for fitting this model
		 *          are called.
		 */
		double getSigma() const;

		/**
		 * @brief
		 *          Replaces the stretching factor of this tented arch
		 *          model by the specified value.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function replaces the value of \f$\lambda\f$
		 *          by the value specified by <code>lambda</code>.
		 *
		 * @param lambda
		 *          The new value for \f$\lambda\f$.
		 *
		 * @see getLambda()
		 */
		void setLambda( double lambda);

		/**
		 * @brief
		 *          Replaces the distance of the poles of this tented arch
		 *          model by the specified value.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function replaces the value for \f$R\f$
		 *          by the value specified by <code>poleDistance</code>.
		 *
		 * @param poleDistance
		 *          The new value for \f$R\f$.
		 *
		 * @see getPoleDistance()
		 */
		void setPoleDistance( double poleDistance );

		/**
		 * @brief
		 *          Replaces the position of this model's delta on its
		 *          longitudinal axis by the specified value.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function replaces the value of \f$d_{delta}\f$
		 *          by the value specified by <code>deltaDistance</code>.
		 *
		 * @param deltaDistance
		 *          The value specifying \f$d_{delta}\f$.
		 *
		 * @see getDeltaDistance()
		 */
		void setDeltaDistance( double deltaDistance );

		/**
		 * @brief
		 *          Replaces the position of this model's core on its
		 *          longitudinal axis by the specified value.
		 *
		 * @details
		 *          This tented arch is a (possibly moved) complex function
		 *          of the form
		 *          \f[
		 *           \tau(z)=\lambda^2\cdot(z^2-R^2)^2\cdot
		 *           \frac{z^2+d_{core}^2}{z^2+d_{delta}^2}.
		 *          \f]
		 *          This function replaces the value of \f$d_{core}\f$
		 *          by the value specified by <code>coreDistance</code>.
		 *
		 * @param coreDistance
		 *          The value specifying \f$d_{core}\f$.
		 *
		 * @see getCoreDistance()
		 */
		void setCoreDistance( double coreDistance );

		/**
		 * @brief
		 *          Replaces the distance from this tented arch's core
		 *          at which orientations are accounted with the highest
		 *          weight.
		 *
		 * @details
		 *          An orientation measured at a complex <i>z</i> is accounted
		 *          with the weight
		 *          \f[
		 *          \exp\left(\frac{(|z-\gamma|-\rho)^2}{2\cdot\sigma}\right)
		 *          \f]
		 *          when the model is fitted to an orientation field
		 *          estimation (see \link weight()\endlink). This method
		 *          serves as a setter for the value of \f$\rho\f$.
		 *
		 * @param rho
		 *          The new value for \f$\rho\f$.
		 *
		 * @see getRho()
		 */
		void setRho( double rho );

		/**
		 * @brief
		 *          Replaces the Gaussian standard deviation used to compute the
		 *          weight with which orientations are accounted on fitting the
		 *          model to an orientation field estimation.
		 *
		 * @details
		 *          An orientation measured at a complex <i>z</i> is accounted
		 *          with the weight
		 *          \f[
		 *          \exp\left(\frac{(|z-\gamma|-\rho)^2}{2\cdot\sigma}\right)
		 *          \f]
		 *          when the model is fitted to an orientation field
		 *          estimation (see \link weight()\endlink). This method
		 *          serves as a setter for the value of \f$\sigma\f$.
		 *
		 * @param sigma
		 *          The new value for \f$\sigma\f$.
		 *
		 * @see getSigma()
		 *
		 * @attention
		 *          Note that if the specified value is very close or even
		 *          equals to 0.0, subsequent fitting attempts may fail due
		 *          to limited machine accuracy.
		 */
		void setSigma( double sigma );
	};
}


#endif /* THIMBLE_TENTEDARCHMODEL_H_ */
