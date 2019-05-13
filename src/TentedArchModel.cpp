/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
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
 * @file TentedArchModel.cpp
 *
 * @brief
 *            Implements the mechanisms provided by the 'TentedArchModel.h'
 *            header file which defines a class for representing the ridge
 *            flow of a tented arch type fingerprint and for them to
 *            fingerprint orientation field estimations.
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include "config.h"
#include <cmath>
#include <cfloat>
#include <complex>
#include <vector>
#include <algorithm>
#include <iostream>

#include <thimble/math/numerical/RealFunctional.h>
#include <thimble/math/numerical/NumericalDifferentiation.h>
#include <thimble/math/numerical/SteepestDescent.h>
#include <thimble/image/Orientation.h>
#include <thimble/finger/Fingerprint.h>
#include <thimble/finger/TentedArchModel.h>

#include <thimble/image/GrayImage.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Creates a tented arch model w.r.t. the standard
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
	TentedArchModel::TentedArchModel( int dpi ) {

		// Ensure that the specified resolution is greater than 0.
		if ( dpi <= 0 ) {
			cerr << "TentedArchModel:: "
				 << "resolution must be greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		// Set parameters in dependence of the specified resolution. The
		// specific values have been specified during a parameter tuning
		// with fingerprints scanned at 569 dpi whose final result
		// can be found in 'Tams, Mihailescu, Munk (2014)'; note that in
		// 'Tams (2013)' an intermediate result is reported as the training
		// was very time consuming and was not finished at submitting
		// the manuscript 'Tams (2013)'
		this->halfGridDist
		         = (int)ceil((3.0/569.0)*(double)dpi);
		this->n             = 1;
		this->lambda        = 1.81;
		this->poleDistance  = (175.0/569.0) * (double)dpi;
		this->coreDistance  = (160.0/569.0) * (double)dpi;
		this->deltaDistance = (22.0/569.0)  * (double)dpi;
		this->rho           = (45.0/569.0)  * (double)dpi;
		this->sigma         = (12.0/569.0)  * (double)dpi;

		// Some constant fields to avoid magic number and to let them
		// become variable in the future easily.
		this->maxInitialCandidates = 20;
		this->derivPrecision       = 1e-6;
		this->minimizerPrecision   = 1e-6;
		this->maxIterations        = 10000;
		this->numRefinements       = 10;

		// Initial rotation and translation part.
		this->alpha = complex<double>(1.0,0.0);
		this->beta = complex<double>(0.0,0.0);
	}

	/**
	 * @brief
	 *            Copy constructor.
	 *
	 * @param model
	 *            The tented arch model of which is copy is created.
	 */
	TentedArchModel::TentedArchModel( const TentedArchModel & model ) {
		// Use the procedural assignment method
		assign(model);
	}

	/**
	 * @brief
	 *            Assignment operator.
	 *
	 * @param model
	 *            The tented arch model of which this tented arch
	 *            is assigned a copy.
	 */
	TentedArchModel & TentedArchModel::operator=
			( const TentedArchModel & model ) {
		assign(model);
		return *this;
	}

	/**
	 * @brief
	 *            Assignment operator (procedural version).
	 *
	 * @param model
	 *            The tented arch model of which this tented arch
	 *            is assigned a copy.
	 */
	void TentedArchModel::assign( const TentedArchModel & model ) {

		// Copy all primitive members ...
		this->halfGridDist             = model.halfGridDist;
		this->n                        = model.n;
		this->lambda                   = model.lambda;
		this->poleDistance             = model.poleDistance;
		this->coreDistance             = model.coreDistance;
		this->deltaDistance            = model.deltaDistance;
		this->rho                      = model.rho;
		this->sigma                    = model.sigma;
		this->maxInitialCandidates     = model.maxInitialCandidates;
		this->derivPrecision           = model.derivPrecision;
		this->minimizerPrecision       = model.minimizerPrecision;
		this->maxIterations            = model.maxIterations;
		this->numRefinements           = model.numRefinements;

		// .. and use the =-operator of the non-primitive
		// 'std::complex' template-class.
		this->alpha = model.alpha;
		this->beta  = model.beta;
	}

	/**
	 * @brief
	 *            Access the reference point estimation after the
	 *            tented arch has been fitted to an orientation field.
	 *
	 * @details
	 *            This method is a convenience method return the result
	 *            of the \link getComplexCore()\endlink whose real part
	 *            is considered as the reference point's abscissa and its
	 *            imaginary part as its ordinate coordinate.
	 *
	 * @return
	 *            \link getComplexCore()\endlink
	 *
	 * @see getComplexCore()
	 */
	complex<double> TentedArchModel::getComplexReferencePoint() const {
		return getComplexCore();
	}

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
	double TentedArchModel::getReferencePointDirection() const {

		// The argument of the rotation part's conjugate (and not of
		// unconjugated rotation part) is used as the angle to be in
		// accordance of the convention of angles in image processing
		// in which they are interpreted clock-wise
		// (and not anti-clockwise as typically interpreted mathematically).
		double angle = arg(conj(this->alpha));

		// Ensure the angle is in the interval '[0,2\pi)'
		while ( angle < 0.0 ) {
			angle += M_PI+M_PI;
		}
		while ( angle >= M_PI+M_PI ) {
			angle -= M_PI+M_PI;
		}

		return angle;
	}

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
	AffineTransform TentedArchModel::getPrealigningTransform() const {

		complex<double> refPoint = getComplexReferencePoint();
		double angle = getReferencePointDirection();

		AffineTransform f;
		f.v = -real(refPoint);
		f.w = -imag(refPoint);
		f.irotate(angle);

		return f;
	}

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
	complex<double> TentedArchModel::getComplexOrigin() const {
        // Note that, since |alpha|=1, the conjugation operations is
        // equivalent to inversion.
		return -(this->beta * conj(this->alpha));
	}

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
	complex<double> TentedArchModel::getComplexCore() const {
        complex<double> tmp(
          -real(this->beta),
          -(imag(this->beta)+this->coreDistance));
        tmp *= conj(this->alpha);
        return tmp;
	}

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
	complex<double> TentedArchModel::getComplexDelta() const {
        complex<double> tmp(
          -real(this->beta),
          -(imag(this->beta)+this->deltaDistance));
        tmp *= conj(this->alpha);
        return tmp;
	}

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
    std::complex<double> TentedArchModel::getAlpha() const {
        return this->alpha;
    }

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
    std::complex<double> TentedArchModel::getBeta() const {
        return this->beta;
    }

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
     void TentedArchModel::setOrigin( const std::complex<double> & origin ) {
		this->beta = -(this->alpha * origin);
     }

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
     void TentedArchModel::translate( const std::complex<double> & s ) {
		this->beta += s;
     }

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
     void TentedArchModel::rotate( double theta ) {

		double c , s;

		c = cos(theta);
		s = sin(theta);

		complex<double> r(c,s);

        this->alpha *= r;

        complex<double> tmp(0.0,coreDistance);
        this->beta = r*(tmp+this->beta)-tmp;
    }

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
    bool TentedArchModel::estimate
	( Fingerprint & fingerprint ) {

        int h = fingerprint.getGradientOrientationDistance();

        int M , N;
        M = fingerprint.getHeight();
        N = fingerprint.getWidth();

        // Estimate the orientation field
		vector<OrientedPoint> samples;
        int g = halfGridDist+halfGridDist+1;
        for ( int y = h ; y+h < M ; y += g ) {
            for ( int x = h ; x+h < N ; x += g ) {
				OrientedPoint p;
				p.x = (double)x;
				p.y = (double)y;
				p.orientation = fingerprint.getOrientationAt(y,x);
				samples.push_back(p);
			}
		}

        // Build the grid for the cores of candidates of inital
        // tented arches.
		vector< pair<double,double> > grid;
        for ( int y = h+g/2 ; y < M ; y += g ) {
            for ( int x = h+g/2 ; x < N ; x += g ) {
				if ( fingerprint.isForeground(y,x) ) {
					grid.push_back( pair<double,double>(x,y) );
				}
			}
		}

        // The final fit is performed by the lower level procedure
        return estimate(samples,fingerprint.getForegroundImage(),M,N,grid);
     }

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
	bool TentedArchModel::estimate
		( const vector<OrientedPoint> & samples ,
		  const bool *foregroundImage , int m , int n ,
		  const vector< pair<double,double> > & grid ) {

		if ( m <= 0 || n <= 0 ) {
			cerr << "TentedArchModel::estimate: invalid image dimension."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( grid.size() == 0 ) {
			return false;
		}

		TentedArchModel tmp(*this);

		vector< pair<double,double> > sortedGrid;
		sortedGrid.reserve(grid.size());

		{
			vector< double > sortedDistances;
			sortedDistances.reserve(grid.size());
			for ( int i = 0 ; i < (int)grid.size() ; i++ ) {

				tmp.alpha = complex<double>(1.0,0.0);
				tmp.beta = complex<double>(-grid[i].first,-grid[i].second-this->coreDistance);

                double dist = tmp.qdistance(samples);

				if ( sortedGrid.size() == 0 ||
					 sortedDistances.back() <= dist) {

					sortedGrid.push_back(grid[i]);
					sortedDistances.push_back(dist);

				} else {
					for ( int j = 0 ; j < (int)sortedGrid.size() ; j++ ) {

						if ( sortedDistances[j] > dist ) {

							sortedGrid.insert
								(sortedGrid.begin()+j,grid[i]);

							sortedDistances.insert
								(sortedDistances.begin()+j,dist);

							break;
						}
					}
				}

			}

		}

		// Iterate over the sorted grid points.
		for ( int j = 0 ;
			  j < this->maxInitialCandidates && sortedGrid.size() > 0 ;
			  j++ ) {

			tmp.alpha = complex<double>(1.0,0.0);
			tmp.beta = complex<double>
			(-sortedGrid.front().first,
			 -sortedGrid.front().second-this->coreDistance);

			sortedGrid.erase(sortedGrid.begin());

			{ // Check whether everything converged.
				bool converged = true;

				for ( int i = 0 ; i < this->numRefinements ; i++ ) {
					if ( !tmp.fitRotation(samples) ) {
						converged = false;
						break;
					}
					if ( !tmp.fitTranslation(samples) ) {
						converged = false;
						break;
					}
				}

				if ( !converged ) {
					continue;
				}
			}

			// Return only if the core is within the fingerprint's foreground.
			// If not, continue with the next grid point.
			complex<double> core = tmp.getComplexCore();
			int x , y;
			x = (int)THIMBLE_ROUND(real(core));
			y = (int)THIMBLE_ROUND(imag(core));
			if ( x >= 0 && y >= 0 && x < n && y < m &&
				 foregroundImage[y*n+x] ) {

				assign(tmp);
				return true;
			}

		}

		// Not successful determining a valid coordinate system.
		return false;
	}

	/**
	 * @brief
	 *             Implements a real functional to be minimized on
	 *             fitting the translation part.
	 *
	 * @see fitTranslation()
	 * @see RealFunctional
	 */
	class TentedArchModel::TranslationParameterFunction :
			public RealFunctional {

	public:

		/**
		 * @brief
		 *            This field will be a copy of the tented arch model which
		 *            is needed to evaluate how well the model moved by a
		 *            candidate shift approximates the orientation field
		 *            estimation.
		 */
		TentedArchModel copy;

		/**
		 * @brief
		 *            A pointer to the orientation field estimation to
		 *            which the tented arch model is fitted.
		 */
		const vector<OrientedPoint> *samplesPtr;

		/**
		 * @brief
		 *            Constructor of the parameter function.
		 *
		 * @param model
		 *            The tented arch model of which the translation is
		 *            updated to better fit the orientation field estimation
		 *            given by <code>samples</code>.
		 *
		 * @param samples
		 *            The orientation field to which <code>model</code>'s
		 *            translation part is fitted.
		 */
		TranslationParameterFunction
		( const TentedArchModel & model ,
		  const vector<OrientedPoint> & samples );

		/**
		 * @brief
		 *            Returns 2 which is the dimension of the translation
		 *            part.
		 *
		 * @return
		 *            The integer 2.
		 */
		virtual int getPreimageDimension() const;

		/**
		 * @brief
		 *            Evaluates the model specified on construction of this
		 *            parameter function moved by the specified translation
		 *            vector.
		 *
		 * @param x
		 *            An array containing two valid double values
		 *            <code>x[0]</code> and <code>x[1]</code>.
		 *
		 * @return
		 *            The \link distance() distance\endlink of the tented arch
		 *            model specified by the field <code>model</code> moved by
		 *            the shift <i>(x[0],x[1])</i> to the orientation field
		 *            specified by the field <code>samplesPtr</code>.
		 */
		virtual double eval( double *x ) const;

		/**
		 * @brief
		 *            Approximates the gradient of this functional
		 *            at the specified values.
		 *
		 * @details
		 *            The implementation of this method is needed to find
		 *            a fit of the translation via a numerical steepest
		 *            descent method and is wrapped around the
		 *            \link NumericalDifferentiation::jacobi()\endlink
		 *            method approximating the gradient using
		 *            the \link eval()\endlink function of this class.
		 *            For more details, see the documentation
		 *            of the \link SteepestDescent\endlink class.
		 *
		 * @param grad
		 *            A two-dimensional array to contain the gradient of this
		 *            functional.
		 *
		 * @param x
		 *            A two-dimensional array containing the values at which
		 *            the gradient is approximated.
		 */
		virtual void derive( double *grad , double *x ) const;
	};

	/**
	 * @brief
	 *            Constructor of the parameter function.
	 *
	 * @param model
	 *            The tented arch model of which the translation is
	 *            updated to better fit the orientation field estimation
	 *            given by <code>samples</code>.
	 *
	 * @param samples
	 *            The orientation field to which <code>model</code>'s
	 *            translation part is fitted.
	 */
	TentedArchModel::TranslationParameterFunction::
	TranslationParameterFunction
	( const TentedArchModel & model ,
	  const vector<OrientedPoint> & samples ) {

		 this->copy       = model;
		 this->samplesPtr = &samples;
	}

	/**
	 * @brief
	 *            Returns 2 which is the dimension of the translation
	 *            part.
	 *
	 * @return
	 *            The integer 2.
	 */
	int TentedArchModel::TranslationParameterFunction::
	getPreimageDimension() const {
		return 2;
	}

	/**
	 * @brief
	 *            Evaluates the model specified on construction of this
	 *            parameter function moved by the specified translation
	 *            vector.
	 *
	 * @param x
	 *            An array containing two valid double values
	 *            <code>x[0]</code> and <code>x[1]</code>.
	 *
	 * @return
	 *            The \link distance() distance\endlink of the tented arch
	 *            model specified by the field <code>model</code> moved by
	 *            the shift <i>(x[0],x[1])</i> to the orientation field
	 *            specified by the field <code>samplesPtr</code>.
	 */
	double TentedArchModel::TranslationParameterFunction::
	eval( double *x ) const {

		TentedArchModel tmp(this->copy);

		tmp.translate(complex<double>(x[0],x[1]));

		return tmp.distance(*samplesPtr);
	}

	/**
	 * @brief
	 *            Approximates the gradient of this functional
	 *            at the specified values.
	 *
	 * @details
	 *            The implementation of this method is needed to find
	 *            a fit of the translation via a numerical steepest
	 *            descent method and is wrapped around the
	 *            \link NumericalDifferentiation::jacobi()\endlink
	 *            method approximating the gradient using
	 *            the \link eval()\endlink function of this class.
	 *            For more details, see the documentation
	 *            of the \link SteepestDescent\endlink class.
	 *
	 * @param grad
	 *            A two-dimensional array to contain the gradient of this
	 *            functional.
	 *
	 * @param x
	 *            A two-dimensional array containing the values at which
	 *            the gradient is approximated.
	 */
	void TentedArchModel::TranslationParameterFunction::derive
		( double *grad , double *x ) const {

		NumericalDifferentiation::jacobi
			(grad,*this,x,this->copy.getDerivPrecision());
	}

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
	bool TentedArchModel::fitTranslation
	( const vector<OrientedPoint> & samples ) {

		TranslationParameterFunction f(*this,samples);

		double x[2];
		x[0] = 0.0;
		x[1] = 0.0;

		SteepestDescent minimizer
			(this->maxIterations,this->minimizerPrecision);

		SD_STATE_T state = minimizer.minimize(x,f,x);
		if ( state == SD_DIVERGED ) {
			return false;
		}

		translate(complex<double>(x[0],x[1]));

		return true;
	}

	/**
	 * @brief
	 *             Implements a real functional to be minimized on
	 *             fitting the rotation part.
	 *
	 * @see fitRotation()
	 * @see RealFunctional
	 */
	class TentedArchModel::RotationParameterFunction : public RealFunctional {

	public:

		/**
		 * @brief
		 *            This field will be a copy of the tented arch model which
		 *            is needed to evaluate how well the rotated model
		 *            approximates the orientation field estimation.
		 */
		TentedArchModel copy;

		/**
		 * @brief
		 *            A pointer to the orientation field estimation to
		 *            which the tented arch model is fitted.
		 */
		const vector<OrientedPoint> *samplesPtr;

		/**
		 * @brief
		 *            Constructor of the parameter function.
		 *
		 * @param model
		 *            The tented arch model of which an update of the rotation
		 *            part can be computed by minimizing this parameter function.
		 *
		 * @param samples
		 *            The orientation field to which <code>model</code>'s
		 *            rotation part is fitted.
		 */
		RotationParameterFunction
			( const TentedArchModel & model ,
			  const vector<OrientedPoint> & samples );

		/**
		 * @brief
		 *            Returns 1 which is the dimension of the rotation angle.
		 *
		 * @return
		 *            The integer 1.
		 */
		virtual int getPreimageDimension() const;

		/**
		 * @brief
		 *            Evaluates the model specified on construction of this
		 *            parameter function rotated by the specified angle.
		 *
		 * @param x
		 *            An array containing on valid double value
		 *            <code>x[0]</code>.
		 *
		 * @return
		 *            The \link distance() distance\endlink of the tented arch
		 *            model specified by the field <code>model</code> rotated
		 *            by the angle <code>x[0]</code> around the
		 *            \link getComplexCore() core\endlink.
		 */
		virtual double eval( double *x ) const;

		/**
		 * @brief
		 *            Approximates the derivation of this function.
		 *
		 * @details
		 *            The implementation of this method is needed to find
		 *            a fit of the translation via a numerical steepest
		 *            descent method and is wrapped around the
		 *            \link NumericalDifferentiation::jacobi()\endlink
		 *            method approximating the derivation using
		 *            the \link eval()\endlink function of this class.
		 *            For more details, see the documentation
		 *            of the \link SteepestDescent\endlink class.
		 *
		 * @param grad
		 *            A pointer to a double value to contain the derivation.
		 *
		 * @param x
		 *            A pointer to a double value specifying the value at
		 *            which the derivation of this function is approximated.
		 */
		virtual void derive( double *grad , double *x ) const;
	};

	/**
	 * @brief
	 *            Constructor of the parameter function.
	 *
	 * @param model
	 *            The tented arch model of which an update of the rotation
	 *            part can be computed by minimizing this parameter function.
	 *
	 * @param samples
	 *            The orientation field to which <code>model</code>'s
	 *            rotation part is fitted.
	 */
	TentedArchModel::RotationParameterFunction::RotationParameterFunction
	( const TentedArchModel & model ,
	  const vector<OrientedPoint> & samples ) {

		this->copy       = model;
		this->samplesPtr = &samples;
	}

	/**
	 * @brief
	 *            Returns 1 which is the dimension of the rotation angle.
	 *
	 * @return
	 *            The integer 1.
	 */
	int TentedArchModel::RotationParameterFunction::
	getPreimageDimension() const {

		return 1;
	}

	/**
	 * @brief
	 *            Evaluates the model specified on construction of this
	 *            parameter function rotated by the specified angle.
	 *
	 * @param x
	 *            An array containing on valid double value
	 *            <code>x[0]</code>.
	 *
	 * @return
	 *            The \link distance() distance\endlink of the tented arch
	 *            model specified by the field <code>model</code> rotated
	 *            by the angle <code>x[0]</code> around the
	 *            \link getComplexCore() core\endlink.
	 */
	double TentedArchModel::RotationParameterFunction::
	eval( double *x ) const {

		TentedArchModel tmp(this->copy);

		tmp.rotate(x[0]);

		return tmp.distance(*(this->samplesPtr));
	}

	/**
	 * @brief
	 *            Approximates the derivation of this function.
	 *
	 * @details
	 *            The implementation of this method is needed to find
	 *            a fit of the translation via a numerical steepest
	 *            descent method and is wrapped around the
	 *            \link NumericalDifferentiation::jacobi()\endlink
	 *            method approximating the derivation using
	 *            the \link eval()\endlink function of this class.
	 *            For more details, see the documentation
	 *            of the \link SteepestDescent\endlink class.
	 *
	 * @param grad
	 *            A pointer to a double value to contain the derivation.
	 *
	 * @param x
	 *            A pointer to a double value specifying the value at
	 *            which the derivation of this function is approximated.
	 */
	void TentedArchModel::RotationParameterFunction::derive
		( double *grad , double *x ) const {

		NumericalDifferentiation::jacobi
			(grad,*this,x,this->copy.getDerivPrecision());
	}

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
	bool TentedArchModel::fitRotation
	( const std::vector<OrientedPoint> & samples ) {

		RotationParameterFunction f(*this,samples);

		double x[1];
		x[0] = 0.0;

		SteepestDescent minimizer
			(this->maxIterations,this->minimizerPrecision);

		SD_STATE_T state = minimizer.minimize(x,f,x);
        if ( state == SD_DIVERGED ) {
			return false;
		}

		rotate(x[0]);

		return true;
	}

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
	double TentedArchModel::distance
		( const vector<OrientedPoint> & samples ) const {

        complex<double> gamma = getComplexCore();

		int l = (int)samples.size();
        double kappa = 0.0;

		for ( int i = 0 ; i < l ; i++ ) {

			complex<double> z
				(samples[i].x,samples[i].y);
            complex<double> v
            	(samples[i].orientation.get_dx2(),
            	 samples[i].orientation.get_dy2());

            complex<double> q = eval(z);
            q = q / abs(q);

            double w = weight(z,gamma);
            kappa += w * norm(q-v);
		}

        return kappa;
	}

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
    double TentedArchModel::qdistance
        ( const vector<OrientedPoint> & samples ) const {

        complex<double> gamma = getComplexCore();
        double R = 2.0*(this->sigma)+this->rho;
        double R2 = R*R;

        int l = (int)samples.size();
        double kappa = 0.0;

        for ( int i = 0 ; i < l ; i++ ) {

            if ( fabs(real(gamma)-samples[i].x) <= R &&
            	 fabs(imag(gamma)-samples[i].y) <= R ) {

                complex<double> z(samples[i].x,samples[i].y);

                if ( norm(z-gamma) <= R2 ) {

                    complex<double> v
                    	(samples[i].orientation.get_dx2(),
                    	 samples[i].orientation.get_dy2());

                    complex<double> q = eval(z);
                    q = q / abs(q);

                    double w = weight(z,gamma);
                    kappa += w * norm(q-v);
                }
            }
        }

        return kappa;
    }

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
    double TentedArchModel::weight( const std::complex<double> & z ) const {

        double x = abs(z-getComplexCore());

        double w = exp
            (-0.5*(x-this->rho)*(x-this->rho)/(this->sigma*this->sigma));

        return w;
    }

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
    double TentedArchModel::weight
    ( const std::complex<double> & z ,
      const std::complex<double> & gamma ) const {

        double x = abs(z-gamma);

        double w = exp
            (-0.5*(x-this->rho)*(x-this->rho)/(this->sigma*this->sigma));

        return w;
    }

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
	complex<double> TentedArchModel::eval( const std::complex<double> & z ) const {

		complex<double> y , q;

		y = this->alpha*z+this->beta;
		q = y*y;
		q = complex<double>(real(q)-this->poleDistance*this->poleDistance,imag(q));
		q = pow(q,this->n+this->n);
		q = complex<double>(real(q),imag(q)*this->lambda*this->lambda);

		complex<double> y2 = y*y;

		q *= y2+this->coreDistance*this->coreDistance;
		q /= y2+this->deltaDistance*this->deltaDistance;

		if ( y.imag() >= 0.0 ) {
			q = complex<double>(abs(q),0.0);
		}

		q *= conj(this->alpha * this->alpha);

		return q;
	}

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
	double TentedArchModel::getDerivPrecision() const {
		return this->derivPrecision;
	}

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
	double TentedArchModel::getMinimizerPrecision() const {
		return this->minimizerPrecision;
	}

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
    int TentedArchModel::getHalfGridDistance() const {
        return this->halfGridDist;
    }

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
    int TentedArchModel::getMaxInitialCandidates() const {
    	return this->maxInitialCandidates;
    }

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
    void TentedArchModel::setHalfGridDistance( int halfGridDist ) {

    	if ( halfGridDist < 0 ) {
    		cerr << "TentedArchModel::setHalfGridDistance: "
    			 << "argument must be non-negative." << endl;
    		exit(EXIT_FAILURE);
    	}

    	this->halfGridDist = halfGridDist;
    }

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
    void TentedArchModel::setMaxInitialCandidates
    ( int maxInitialCandidates ) {

    	if ( maxInitialCandidates <= 0 ) {
    		cerr << "TentedArchModel::setMaxInitialCandidates: "
    			 << "argument must be greater than zero." << endl;
    		exit(EXIT_FAILURE);
    	}

    	this->maxInitialCandidates = maxInitialCandidates;
    }

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
    double TentedArchModel::getLambda() const {
        return this->lambda;
    }

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
    double TentedArchModel::getPoleDistance() const {
        return this->poleDistance;
    }

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
    double TentedArchModel::getDeltaDistance() const {
        return this->deltaDistance;
    }

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
    double TentedArchModel::getCoreDistance() const {
        return this->coreDistance;
    }

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
    double TentedArchModel::getRho() const {
        return this->rho;
    }

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
    double TentedArchModel::getSigma() const {
        return this->sigma;
    }

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
    void TentedArchModel::setLambda( double lambda) {
		this->lambda = lambda;
    }

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
     void TentedArchModel::setPoleDistance( double poleDistance ) {
		this->poleDistance = poleDistance;
     }

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
    void TentedArchModel::setDeltaDistance( double deltaDistance ) {
		this->deltaDistance = deltaDistance;
    }

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
    void TentedArchModel::setCoreDistance( double coreDistance ) {
		this->coreDistance = coreDistance;
    }

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
    void TentedArchModel::setRho( double rho ) {
		this->rho = rho;
    }

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
    void TentedArchModel::setSigma( double sigma ) {
		this->sigma = sigma;
    }
}
