/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file Fingerprint.h
 *
 * @brief
 *            Provides a class that represents a fingerprint and
 *            mechanisms for estimating certain fingerprint features.
 *
 * @details
 * @section sec_fingerprints Fingerprints
 *
 * A <em>fingerprint</em> is defined to be the traces that the ridges of a
 * human's finger leave on a surface. Typically, these traces are measured
 * by a device (e.g., a scanner) and transformed into digital images which
 * can be processed by computers. An interface to digital gray-scale
 * images is provided by the \link thimble::GrayImage GrayImage\endlink
 * class (see documentation on @ref sec_grayimages).
 *
 * The \link thimble::Fingerprint Fingerprint\endlink class which can be
 * created from a \link thimble::GrayImage GrayImage\endlink represents
 * fingerprints in the context of the THIMBLE library. In its current
 * version it provides member functions and member methods for
 * <ul>
 *  <li>@ref sec_of</li>
 *  <li>@ref sec_seg</li>
 *  <li>@ref sec_drp</li>
 * </ul>
 * There is even a very basic mechanism for processing with
 * minutiae templates (although, in the present version, not for
 * minutiae estimation):
 * <ul>
 *  <li>@ref sec_minutiae_estimation</li>
 * </ul>
 *
 * @subsection sec_of Orientation Field Estimation
 *
 * Local orientations, i.e., the angles of the ridge
 * flow, can be estimated from a
 * <pre>
 *  Fingerprint fingerprint;
 *
 *  fingerprint.fromImageFile("fingerprint.pgm");
 * </pre>
 * that was read from the file <code>fingerprint.pgm</code>, say, using
 * the convenience function
 * \link thimble::Fingerprint::fromImageFile fromImageFile()\endlink.
 * The interface for accessing a fingerprint's local orientation is
 * provided by the \link thimble::Fingerprint::getOrientationAt(int,int)
 * getOrientationAt(int,int)\endlink function. More specifically,
 * to estimate the orientation at the <code>i</code>th row and
 * <code>j</code>th column of a fingerprint
 * (<code>i=0,...,fingerprint.getHeight()-1</code>,
 *  <code>j=0,...,fingerprint.getWidth()-1</code>) we may call
 * <pre>
 *  Orientation orientation = fingerprint.getOrientationAt(i,j);
 * </pre>
 * where an \link thimble::Orientation Orientation\endlink object is
 * returned representing the local orientation of a pattern. The
 * \link thimble::Orientation::getOrientationAngle()
 * getOrientationAngle()\endlink function
 * can be used to access the angle of an \link thimble::Orientation
 * Orientation\endlink object
 * <pre>
 *  double angle = orientation.getOrientationAngle();
 * </pre>
 * in the range \f$[0,\pi)\f$. If <code>angle</code> is not a valid
 * <code>double</code> value between \f$[0,\pi)\f$ the orientation
 * estimation at the corresponding pixel should be considered as being
 * invalid (e.g., due to a local orientation of being non-existent
 * or border effects).
 *
 * Sometimes, it may be useful to obtain
 * a degree representation of the orientation between 0 and 179.
 * Therefore, we may make use of the convenience function
 * \link thimble::Orientation::getDegree() getDegree()\endlink
 * from an \link thimble::Orientation Orientation\endlink object
 * <pre>
 *  int degree = orientation.getDegree();
 * </pre>
 * Again, if <code>degree</code> is a value different from the
 * interval [0,179] the estimation at the pixel is considered of having
 * failed.
 *
 * Consider the following example.
 *
 * <table border="0">
 *  <tr>
 *   <td>\image html a_fingerprint.png %Fingerprint</td>
 *   <td>\image html a_of.png %Orientation Field</td>
 *  </tr>
 * </table>
 *
 * The above two images visualize the result of an orientation field
 * estimation using the
 * \link thimble::Fingerprint::getOrientationAt(int,int)
 * getOrientationAt(int,int)\endlink function of a
 * \link thimble::Fingerprint Fingerprint\endlink object.
 * The left image shows a fingerprint (being contained in
 * <a href="http://bias.csr.unibo.it/fvc2000/download.asp"
 *  target="_blank">
 *  FVC 2000 DB2_B</a>) from which a \link thimble::Fingerprint
 * Fingerprint\endlink object can be created; the right image visualizes
 * each of the fingerprint's pixel orientation estimated with the
 * \link thimble::Fingerprint::getOrientationAt(int,int)
 * getOrientationAt(int,int)\endlink function for each of the
 * fingerprint's pixel. The values of an orientation vary between
 * 0 and 179 degrees of which values are visualized in the 8-bit
 * orientation field image shown above.
 *
 * @subsection sec_seg Foreground Estimation
 *
 * On an image containing a fingerprint there exist region in which
 * no traces of a finger's ridges are shown. An interface fingerprint
 * segmentation, that is the process of dividing
 * a fingerprint image into regions that contain the fingerprint
 * (foreground; occasionally called <em>region of interest</em>)
 * and regions that do not show the fingerprint (background), is provided
 * by the \link thimble::Fingerprint::isForeground(int,int)
 * isForeground(int,int)\endlink function.
 *
 * For example, let <code>fingerprint</code> be a
 * \link thimble::Fingerprint Fingerprint\endlink object that has been
 * created from the PGM file <code>fingerprint.pgm</code> (say), i.e.,
 * <pre>
 *  Fingerprint fingerprint;
 *
 *  fingerprint.fromImageFile("fingerprint.pgm");
 * </pre>
 * To estimate whether the pixel at the <code>i</code>th row
 * (<code>i=0,...,fingerprint.getHeight()-1</code>) and
 * <code>i</code>th column
 * (<code>j=0,...,fingerprint.getWidth()-1</code>) of the fingerprint
 * belongs to the foreground, we may use the
 * \link thimble::Fingerprint::isForeground(int,int)
 * isForeground(int,int)\endlink function:
 * <pre>
 *  bool b = fingerprint.isForeground(int,int);
 * </pre>
 * If <code>(i,j)</code> has been estimated as a foreground pixel, the
 * value of <code>b</code> will be <code>true</code>; otherwise, if
 * <code>(i,j)</code> was estimated as a background pixel, the
 * value of <code>b</code> will be <code>false</code>.
 *
 * Consider the following example.
 *
 * <table border="0" target="_blank">
 *  <tr>
 *   <td>\image html b_fingerprint.png %Fingerprint</td>
 *   <td>\image html b_seg.png %Foreground Estimation</td>
 *  </tr>
 * </table>
 *
 * The left image shows a fingerprint (being contained in the
 * <a href="http://bias.csr.unibo.it/fvc2002/download.asp"
 *  target="_blank">FVC 2002 DB1_B</a>) from which a
 * \link thimble::Fingerprint Fingerprint\endlink can be created. The right
 * image visualizes the result of each pixels foreground/background estimation
 * using the \link thimble::Fingerprint::isForeground(int,int)
 * isForeground(int,int)\endlink function: If a pixel was estimated as
 * foreground (<code>true</code>) it is indicated by a black color and,
 * otherwise, by a white color.
 *
 * @subsection sec_drp Directed Reference Point Estimation
 *
 * \link thimble::Fingerprint Fingerprint\endlink objects do provide a
 * mechansim for estimating a fingerprint's reference that are constituted
 * with a direction, i.e., <em>directed reference point</em>. Most methods
 * that have been proposed do only account for estimating a fingerprint's
 * reference point (without direction). However, the presence of a robust
 * direction may be very useful, for example, if the minutiae estimated
 * from a fingerprint are required to be absolutely pre-aligned w.r.t. to
 * the Cartesian coordinate system that can be derived from a directed
 * reference point.
 *
 * Assume that a \link thimble::Fingerprint Fingerprint\endlink object has
 * been created from a PGM file <code>fingerprint.pgm</code> say, e.g.,
 * with the help of the convenience method
 * \link thimble::Fingerprint::fromImageFile() fromImageFile()\endlink:
 * <pre>
 *  Fingerprint fingerprint;
 *
 *  fingerprint.fromImageFile("fingerprint.pgm");
 * </pre>
 * Then, an estimation of the fingerprint's directed reference point can
 * be accessed via
 * <pre>
 *  DirectedPoint drp = fingerprint.getDirectedReferencePoint();
 * </pre>
 * using the \link thimble::Fingerprint::getDirectedReferencePoint()
 * getDirectedReferencePoint()\endlink
 * method returning a reference to a \link thimble::DirectedPoint
 * DirectedPoint\endlink object. The abscissa coordinate of the
 * directed reference point is given by
 * \link thimble::DirectedPoint::x drp.x\endlink and its
 * ordinate coordinate by
 * \link thimble::DirectedPoint::y drp.y\endlink.
 * <pre>
 *  double x , y;
 *
 *  x = drp.x;
 *  y = drp.y;
 * </pre>
 * The angle of the direction \link thimble::DirectedPoint::direction
 * drp.direction\endlink, which is given as a \link thimble::Direction
 * Direction\endlink object,
 * can be accessed using the \link thimble::Direction::getDirectionAngle()
 * getDirectionAngle()\endlink function:
 * <pre>
 *  double angle = drp.direction.getDirectionAngle();
 * </pre>
 * It can happen that directed reference point estimation fails due to
 * the fingerprint's center being not visible or due to (pre-)processing
 * errors. In this case, the attributes of the \link thimble::DirectedPoint
 * DirectedPoint\endlink object returned by the
 * \link thimble::Fingerprint::getDirectedReferencePoint()
 * getDirectedReferencePoint()\endlink function are equals the NAN
 * expression. Alternatively and more conveniently, the function
 * \link thimble::Fingerprint::hasDirectedReferencePoint()
 * hasDirectedReferencePoint()\endlink of a \link thimble::Fingerprint
 * Fingerprint\endlink returns <code>true</code> if directed reference point
 * estimation was successful and, otherwise, <code>false</code>. Consequently,
 * if we want to output a fingerprint's directed reference point only if
 * estimation was successful, we may run
 * <pre>
 *  if ( fingerprint.hasDirectedReferencePoint() ) {
 *
 *     // Print a text representation of the directed reference point
 *     cout << fingerprint.getDirectedReferencePoint() << endl;
 *  }
 * </pre>
 * Below examples for results of the
 * \link thimble::Fingerprint::getDirectedReferencePoint()
 * getDirectedReferencePoint()\endlink for a finger's fingerprints
 * (contained in
 * <a href="http://bias.csr.unibo.it/fvc2002/download.asp"
 *  target="_blank">FVC 2002 DB2_B</a>) are visualized by red
 * circles surrounding the reference points' coordinates and red lines for
 * their direction angle; the captions list the coordinate and direction of
 * each directed reference point estimation; furthermore, the Cartesian
 * coordinate system that can be derived from a directed reference point is
 * visualized in the examples by a thick yellow line (ordinate) and thin
 * yellow line (abscissa).
 * <table>
 *  <tr>
 *   <td>\image html c_drp1.png (x y angle)=(178.43 220.84 6.09)</td>
 *   <td>\image html c_drp2.png (x y angle)=(185.06 249.68 6.14)</td>
 *   <td>\image html c_drp3.png (x y angle)=(181.91 226.60 6.14)</td>
 *   <td>\image html c_drp4.png (x y angle)=(189.63 238.74 6.05)</td>
 *  </tr>
 *  <tr>
 *   <td>\image html c_drp5.png (x y angle)=(174.38 78.52 6.16)</td>
 *   <td>\image html c_drp6.png (x y angle)=(180.00 380.74 6.05)</td>
 *   <td>\image html c_drp7.png (x y angle)=(153.88 450.82 6.05)</td>
 *   <td>\image html c_drp8.png (x y angle)=(162.69 388.70 6.12)</td>
 *  </tr>
 * </table>
 *
 * @subsection sec_minutiae_estimation Estimating Minutiae Templates
 *
 * Object generated from the \link thimble::Fingerprint Fingerprint\endlink
 * class provide the function \link thimble::Fingerprint::getMinutiaeView()
 * getMinutiaeView()\endlink which may suggest that THIMBLE comes along
 * with an automatic method for estimating the minutiae from a fingerprint
 * image. However, the \link thimble::Fingerprint::getMinutiaeView()
 * getMinutiaeView()\endlink returns nothing but an empty minutiae template:
 * THIMBLE currently contains no automatic method for minutiae estimation.
 * However, if there is an external method that estimates minutiae and
 * stores them in ISO 19794-2:2005 format they can be read with the help
 * of the \link thimble::MinutiaeRecord MinutiaeRecord\endlink class and
 * then it is possible to read specify the minutiae template manually:
 * <pre>
 *    bool success;
 *    MinutiaeRecord record;
 *
 *    success = record.read("minutiae.fmr"); // Read minutiae from the file 'minutiae.fmr' in ISO 19794-2:2005 format
 *    if ( !success ) {
 *       cerr << "Could not read 'minutiae.fmr'." << endl;
 *       exit(EXIT_FAILURE);
 *    }
 *
 *    Fingerprint fingerprint;
 *
 *    success = fingerprint.fromImageFile("fingerprint.pgm");
 *    if ( !success ) {
 *       cerr << "Could not read 'fingerprint.pgm'." << endl;
 *       exit(EXIT_FAILURE);
 *    }
 *
 *    // Manual specification the minutiae template computed by an external algorithm
 *    fingerprint.getMinutiaeView() = record.getView();
 * </pre>
 *
 * @subsection sec_fjfx_fingerprint Minutiae Estimation with Digital Persona's FingerJETFX OSE library
 *
 * Another alternative is to extend the \link thimble::Fingerprint
 * Fingerprint\endlink class in which the minutiae estimation is performed
 * by an external C++ library.
 * The \link thimble::Fingerprint::getMinutiaeView()
 * getMinutiaeView()\endlink is wrapped around
 * the \link thimble::Fingerprint::estimateMinutiae()\endlink function which
 * is called only by the first calls of
 * the \link thimble::Fingerprint::getMinutiaeView()
 * getMinutiaeView()\endlink function thereby preventing unnecessary
 * repeated minutiae estimations for the same fingerprint object. We can use
 * Digital Persona's FingerJetFX OSE C++ library
 * (<a href="http://digitalpersona.com/fingerjetfx" target="_blank">
 *  http://digitalpersona.com/fingerjetfx
 *  </a>) to implement an automatic minutiae estimation for fingerprint
 *  objects.
 *
 * We assume that the FJFX library is accessible within the program through
 * the 'fjfx.h' header library and our preambel should contain the following
 * includes
 * <pre>
 *    \#include <fjfx.h>
 *    \#include <thimble/all.h>
 * </pre>
 * Furthermore, we assume as usual that we are using the following
 * namespaces
 * <pre>
 *    using namespace std;
 *    using namespace thimble;
 * </pre>
 * The header for the new subclass may be of the form:
 * <pre>
 *     class FJFXFingerprint : public Fingerprint {
 *
 *     public:
 *
 *        // Overrides the empty minutiae estimation method from the 'Fingerprint' class
 *        virtual MinutiaeView estimateMinutiae();
 *
 *     };
 *  </pre>
 *  Then the implementation of the function may look as follows:
 *  <pre>
 *     // Implementation of the 'virtual Minutiae estimateMinutiae()' function
 *     // estimating a fingerprint image's minutiae using Digital Persona's
 *     // FingerJetFX OSE library
 *     MinutiaeView FJFXFingerprint::estimateMinutiae() {
 *
 *     // The fingerprint's image dimension
 *     int m , n;
 *     m = getHeight();
 *     n = getWidth();
 *
 *     // Allocate memory for raw pixel data of the fingerprint image
 *     uint8_t *raw_image = (uint8_t*)malloc( m * n );
 *     if ( raw_image == NULL ) {
 *        cerr << "FJFXFingerprint::estimateMinutiae: out of memory." << endl;
 *        exit(EXIT_FAILURE);
 *     }
 *
 *     // Convert the intensities from [0.0,1.0] to 8-bit values in [0,255]
 *     for ( int y = 0 ; y < m ; y++ ) {
 *        for ( int x = 0 ; x < n ; x++ ) {
 *           raw_image[y*n+x] = (uint8_t)(getIntensityImage()[y*n+x] * 255.0);
 *        }
 *     }
 *
 *     // fingerprint minutiae data
 *     uint8_t fmd[FJFX_FMD_BUFFER_SIZE];
 *     unsigned int size_of_fmd = FJFX_FMD_BUFFER_SIZE;
 *
 *     // ***********************************************************************
 *     // ****** CALL OF DIGITAL PERSONA'S FINGERJETFX OSE LIBRARY FUNCTION *****
 *     // ***********************************************************************
 *     // ****** Remark: the function allocates some memory of size of orders ***
 *     // ****** the pixel data that is, however, not freed afterwards **********
 *     // ***********************************************************************
 *     int code = fjfx_create_fmd_from_raw
 *     (raw_image,getResolution(),m,n,FJFX_FMD_ISO_19794_2_2005,fmd,&size_of_fmd);
 *
 *     MinutiaeRecord record;
 *
 *     if ( code == FJFX_SUCCESS ) {
 *
 *        // If the data output by Digital Persona's FingerJetFX OSE function
 *        // cannot be read by THIMBLE's MinutiaeRecord class, we would expect
 *        // it to be a bug in THIMBLE; however, an error should not occur and
 *        // we have not experienced any for the present version of THIMBLE
 *        // as of now.
 *        if ( !record.fromBytes(fmd) ) {
 *           cerr << "ERROR: probably a bug in THIMBLE." << endl;
 *           exit(EXIT_FAILURE);
 *        }
 *
 *     } else {
 *
 *        // If, for which reason ever, FJFX has not output a valid minutiae
 *        // estimate, we use an empty view instead which might result in
 *        // a treatment of a 'failure to enrolment' later.
 *        record.addView(MinutiaeView());
 *     }
 *
 *     // Memories, you are free
 *     free(raw_image);
 *
 *     // Return the view
 *     return record.getView();
 *   }
 *  </pre>
 *  Now, objects that are generated from the <code>FJFXFingerprint</code>
 *  class provide an automatic method for minutiae estimation:
 *  <pre>
 *     FJFXFingerprint fingerprint;
 *
 *     // Read the fingerprint image
 *     bool success = fingerprint.fromImageFile("fingerprint.pgm");
 *     if ( !success ) {
 *        cerr << "Could not read 'fingerprint.pgm'." << endl;
 *        exit(EXIT_FAILURE);
 *     }
 *
 *     // Access the non-empty minutiae template
 *     MinutiaeView view = fingerprint.getMinutiaeView();
 *
 *     // Print a text representation of the minutiae template to
 *     // 'cout'
 *     cout << view << endl;
 *  </pre>
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_FINGERPRINT_H_
#define THIMBLE_FINGERPRINT_H_

#include <thimble/dllcompat.h>
#include <thimble/image/Orientation.h>
#include <thimble/image/GrayImage.h>
#include <thimble/finger/MinutiaeRecord.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *           Instances of this class represent fingerprints.
	 */
	class THIMBLE_DLL Fingerprint {

	private:

		/**
		 * @brief
		 *            Holds the resolution in dots per inch at which this
		 *            fingerprint has been scanned.
		 */
		int dpi;

		/**
		 * @brief
		 *            The height in pixels of the fingerprint image.
		 */
		int m;

		/**
		 * @brief
		 *            The width in pixels of the fingerprint image.
		 */
		int n;

		/**
		 * @brief
		 *            Array representing the intensities of the original
		 *            fingerprint image.
		 *
		 * @details
		 *            The intensities of the fingerprint image are stored
		 *            in the array encoding an
		 *            <i>\link m\endlink*\link n\endlink</i>-matrix in row-major
		 *            order, i.e., the intensity at pixel <i>(x,y)</i> where
		 *            <i>x=0,...,\link n\endlink-1</i> and
		 *            <i>y=0,...,\link m\endlink-1</i> is stored at
		 *            <code>intensityImage[y*n+x]</code>.
		 *
		 *            The intensity of a pixel varies between 0.0 (black) and
		 *            1.0 (white).
		 */
		double *intensityImage;

		/**
		 * @brief
		 *            Array encoding for each pixel whether it belongs to
		 *            the fingerprint's foreground or background.
		 *
		 * @details
		 *            The array is <code>NULL</code> if the no foreground
		 *            of the fingerprint has been estimated. Otherwise,
		 *            the boolean foreground information for the pixel
		 *            <i>(x,y)</i> is stored at
		 *            <code>foregroundImage[y*n+x]</code> which is
		 *            <code>true</code> if <i>(x,y)</i> belongs to the
		 *            foreground and <code>false</code> otherwise.
		 *
		 *            The foreground is computed automatically during the
		 *            first call of the \link isForeground(int,int)\endlink
		 *            function.
		 *
		 * @see isForeground(int,int)
		 */
		bool *foregroundImage;

		/**
		 * @brief
		 *            Array encoding the local ridge flow of a pixel on
		 *            the fingerprint image.
		 *
		 * @details
		 *            Initially, the array contains <i>m*n</i> pointers
		 *            to an instance of the \link Orientation\endlink class
		 *            Initially, these pointers are <code>NULL</code> and
		 *            non-<code>NULL</code> only if an orientation has been
		 *            already computed for the corresponding pixel. If
		 *            non-<code>NULL</code>, the orientation at pixel
		 *            <i>(x,y)</i> can be accessed via
		 *            <code>orientationImage[y*n+x][0]</code>.
		 *
		 *            When calling
		 *            \link getOrientationAt(int,int) getOrientation(x,y)\endlink,
		 *            the function checks whether
		 *            <code>orientationImage[y*n+x]</code> is NULL and, if true,
		 *            initialize the entry with a pointer to the pixel's
		 *            orientation estimation to which a reference is returned
		 *            by the function; otherwise, if
		 *            <code>orientationImage[y*n+x]</code> is
		 *            non-<code>NULL</code>, it is assumed that the orientation
		 *            at <i>(x,y)</i> has already been computed and a reference
		 *            to the orientation is returned, thereby avoiding that the
		 *            orientation at <i>(x,y)</i> is estimated another time.
		 */
		Orientation **orientationImage;

		/**
		 * @brief
		 *            Indicates whether this fingerprint contains a valid
		 *            estimation of a directed reference point.
		 *
		 * @details
		 *            The value of this field is <code>false</code> unless
		 *            an attempt for estimating the fingerprint's directed
		 *            reference point has been made via
		 *            \link getDirectedReferencePoint()\endlink. If an
		 *            attempt for estimating the directed reference point
		 *            was successful, the value of the field is set
		 *            <code>true</code> and the field \link drpPtr\endlink
		 *            references to the directed reference point; otherwise,
		 *            if the estimation failed, the value of this field
		 *            is left <code>false</code> and \link drpPtr\endlink
		 *            refers to an instance of the DirectedPoint class but
		 *            is nonetheless assumed to be invalid.
		 */
		bool has_drp;

		/**
		 * @brief
		 *            Reference to the directed reference point estimation of
		 *            the fingerprint.
		 *
		 * @details
		 *            The value of this field is <code>NULL</code> unless an
		 *            attempt for estimating the fingerprint's directed
		 *            reference point via
		 *            \link getDirectedReferencePoint()\endlink has been
		 *            performed. After directed reference point estimation the
		 *            value of the field will be non-<code>NULL</code>
		 *            independent from whether the estimation is considered to
		 *            be successful if not. If successful, the value of
		 *            \link has_drp\endlink will be set to <code>true</code>
		 *            in addition and left <code>false</code> otherwise.
		 */
		DirectedPoint *drpPtr;

		/**
		 * @brief
		 *            Pointer to a minutiae record containing this
		 *            fingerprint's minutiae template.
		 *
		 * @details
		 *            If no minutiae extraction has been performed for this
		 *            fingerprint, the value of this field remains
		 *            <code>NULL</code>. Otherwise, this field contains
		 *            a pointer to a valid \link MinutiaeRecord\endlink
		 *            containing the minutiae template associated with
		 *            this fingerprint.
		 *
		 * @see getMinutiaeRecord()
		 * @see getMinutiaeView()
		 */
		MinutiaeRecord *recordPtr;

		/**
		 * @brief
		 *            Initializes the values of this fingerprint with default
		 *            values.
		 *
		 * @details
		 *            The method is useful to create an empty fingerprint on
		 *            first initialzation where no data must be freed.
		 *            To clear an already initialized fingerprint the method
		 *            \link clear()\endlink should be used.
		 */
		void first_init();

	public:

		/**
		 * @brief
		 *            Standard constructor.
		 *
		 * @details
		 *            Creates an empty fingerprint.
		 */
		Fingerprint();

		/**
		 * @brief
		 *             Creates a fingerprint specified by its raw gray-scale
		 *             image.
		 *
		 * @param grayImage
		 *             The raw fingerprint image.
		 *
		 * @param dpi
		 *             The resolution in <em>dots per inch</em> at which the
		 *             fingerprint image has been scanned.
		 *
		 * @warning
		 *             If <code>dpi</code> is not greater than 0, an error message
		 *             is printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		Fingerprint( const GrayImage & grayImage , int dpi = 500 );

		/**
		 * @brief
		 *            Creates a fingerprint specified by its raw gray-scale
		 *            intensity image.
		 *
		 * @param intensityImage
		 *            An array of <i>m*n</i> double values ranging between
		 *            0.0 (black) and 1.0 (white) where the intensity
		 *            at <i>(x,y)</i> is stored at
		 *            <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *            Height of the fingerprint image.
		 *
		 * @param n
		 *            Width of the fingerprint image.
		 *
		 * @param dpi
		 *            The resolution in <em>dots per inch</em> at which the
		 *            fingerprint image has been scanned.
		 *
		 * @warning
		 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
		 *            are not greater than 0, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		Fingerprint
		( const double *intensityImage , int m , int n , int dpi = 500 );

		/**
		 * @brief
		 *            Copy constructor.
		 *
		 * @param fingerprint
		 *            The fingerprint of which a copy is created.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		Fingerprint( const Fingerprint & fingerprint );

		/**
		 * @brief
		 *           Destructor.
		 *
		 * @details
		 *           Calls \link clear()\endlink which, as a side-effect,
		 *           frees any memory held by this fingerprint.
		 */
		virtual ~Fingerprint();

		/**
		 * @brief
		 *           Clears the data for this fingerprint such that it
		 *           becomes empty.
		 *
		 * @details
		 *           After this method has been called, the function
		 *           \link isEmpty()\endlink returns <code>true</code>.
		 *           More specifically, the result of the functions
		 *           <ul>
		 *            <li>getWidth()</li>
		 *            <li>getHeight()</li>
		 *            <li>getResolution()</li>
		 *           </ul>
		 *           will be 0 right after the fingerprint has been cleared.
		 */
		virtual void clear();

		/**
		 * @brief
		 *           Check whether this fingerprint is empty.
		 *
		 * @return
		 *           <code>true</code> if the fingerprint is not initialized
		 *           by a gray-scale image; otherwise, the function returns
		 *           <code>false</code>.
		 */
		bool isEmpty() const;

		/**
		 * @brief
		 *             Initialize this instance with a new fingerprint
		 *             specified by its raw gray-scale image.
		 *
		 * @param grayImage
		 *             The raw fingerprint image.
		 *
		 * @param dpi
		 *             The resolution in <em>dots per inch</em> at which the
		 *             fingerprint image has been scanned.
		 *
		 * @warning
		 *             If <code>dpi</code> is not greater than 0, an error message
		 *             is printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *             If not enough memory could be provided, an error
		 *             message is printed to <code>stderr</code> and the
		 *             program exits with status 'EXIT_FAILURE'.
		 */
		void initialize( const GrayImage & grayImage , int dpi = 500 );

		/**
		 * @brief
		 *             Initialize this instance with a new fingerprint
		 *             specified by its double-valued gray-scale intensity
		 *             image.
		 *
		 * @param intensityImage
		 *            An array of <i>m*n</i> double values ranging between
		 *            0.0 (black) and 1.0 (white) where the intensity
		 *            at <i>(x,y)</i> is stored at
		 *            <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *            Height of the fingerprint image.
		 *
		 * @param n
		 *            Width of the fingerprint image.
		 *
		 * @param dpi
		 *            The resolution in <em>dots per inch</em> at which the
		 *            fingerprint image has been scanned.
		 *
		 * @warning
		 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
		 *            are not greater than 0, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void initialize
		( const double *intensityImage , int m , int n , int dpi = 500 );

		/**
		 * @brief
		 *             Initialize this instance with a new fingerprint
		 *             specified by its raw 8-bit gray-scale intensity image.
		 *
		 * @param rawImage
		 *            An array of <i>m*n</i> values of type <code>uint8_t</code>
		 *            ranging between 0 (black) and 255 (white) where the
		 *            intensity at <i>(x,y)</i> is stored at
		 *            <code>intensityImage[y*n+x]</code>.
		 *
		 * @param m
		 *            Height of the fingerprint image.
		 *
		 * @param n
		 *            Width of the fingerprint image.
		 *
		 * @param dpi
		 *            The resolution in <em>dots per inch</em> at which the
		 *            fingerprint image has been scanned.
		 *
		 * @warning
		 *            If <code>dpi</code>, <code>m</code> or <code>n</code>
		 *            are not greater than 0, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void initialize
		( const uint8_t *rawImage , int m , int n , int dpi = 500 );

		/**
		 * @brief
		 *            Assignment operator (procedural version).
		 *
		 * @param fingerprint
		 *            The fingerprint assigned to this instance.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		virtual void assign( const Fingerprint & fingerprint );

		/**
		 * @brief
		 *            Assignment operator.
		 *
		 * @param fingerprint
		 *            The fingerprint assigned to this instance.
		 *
		 * @return
		 *            A reference to this instance.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		Fingerprint &operator=( const Fingerprint & fingerprint );

		/**
		 * @brief
		 *            Access the height of the fingerprint image.
		 *
		 * @return
		 *            The height of the fingerprint image.
		 */
		int getHeight() const;

		/**
		 * @brief
		 *            Access the width of the fingerprint image.
		 *
		 * @return
		 *            The width of the fingerprint image.
		 */
		int getWidth() const;

        /**
         * @brief
         *            Access the resolution specified at which this
         *            fingerprint's image has been scanned.
         * @return
         *            The resolution in <i>dots per inch</i>.
         */
        int getResolution() const;

		/**
		 * @brief
		 *            Access the intensity image of the original
		 *            fingerprint image.
		 *
		 * @details
		 *            Assume that
		 *            <pre>
		 *             m = getHeight();
		 *             n = getWidth();
		 *            </pre>
		 *            is the height and width of the fingerprint,
		 *            respectively, and
		 *            <pre>
		 *             const double *img = getIntensityImage()
		 *            </pre>
		 *            is the array returned by the function. Then
		 *            <code>img</code> encodes the intensity image as
		 *            an <i>m*n</i> matrix in row-major order. Specifically,
		 *            the intensity at pixel <i>(x,y)</i>
		 *            (where <i>x=0...m-1, y=0...n-1</i>) can be accessed via
		 *            <pre>
		 *             double intensitiy = img[y*n+x];
		 *            </pre>
		 *
		 * @return
		 *            Intensity image of the original fingerprint image.
		 */
		const double *getIntensityImage() const;

		/**
		 * @brief
		 *            Check whether a pixel is assigned to the fingerprint's
		 *            foreground.
		 *
		 * @param y
		 *            ordinate varying between 0 and
		 *            \link getHeight()\endlink-1.
		 *
		 * @param x
		 *            abscissa varying between 0 and
		 *            \link getWidth()\endlink-1.
		 *
		 *            If the function is called the first time for the
		 *            fingerprint, it causes a call of the
		 *     \link thimble::Fingerprint::estimateForegroundImage(bool*)const
		 *      estimateForegroundImage()\endlink which
		 *            runs an algorithm for fingerprint foreground estimation.
		 *            The result is kept by the fingerprint, thereby
		 *            preventing that the foreground is estimated again.
		 *
		 * @warning
		 *            If <i>(x,y)</i> does not lie between the fingerprint
		 *            image's region, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool isForeground( int y , int x );

		/**
		 * @brief
		 *            Access the foreground image of the fingerprint.
		 *
		 * @details
		 *            Assume that
		 *            <pre>
		 *             m = getHeight();
		 *             n = getWidth();
		 *            </pre>
		 *            is the height and width of the fingerprint,
		 *            respectively, and
		 *            <pre>
		 *             const bool *img = getForegroundImage();
		 *            </pre>
		 *            is the array returned by the function. Then
		 *            <code>img</code> encodes the foreground image as
		 *            an <i>m*n</i> matrix in row-major order.
		 *            Specifically, the boolean indicating whether the pixel
 	 	 *            <i>(x,y)</i> (where <i>x=0...m-1, y=0...n-1</i>) belongs
 	 	 *            to the foreground can be accessed via
		 *            <pre>
		 *             bool isForegroundPixel = img[y*n+x];
		 *            </pre>
		 *
		 *            If the function is called the first time for the
		 *            fingerprint, it causes a call of the
		 *     \link thimble::Fingerprint::estimateForegroundImage(bool*)const
		 *      estimateForegroundImage()\endlink which
		 *            runs an algorithm for fingerprint foreground estimation.
		 *            The result is kept by the fingerprint, thereby
		 *            preventing that the foreground is estimated again.
		 *
		 * @return
		 *            Foreground image of the fingerprint image.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		const bool *getForegroundImage();

		/**
		 * @brief
		 *            Access the orientation of the fingerprint at the
		 *            specified pixel.
		 *
		 * @details
		 *            If the orientation at <i>(x,y)</i> is accessed the first
		 *            time, it is estimated using the
		 *            \link estimateOrientationAt(int,int)const
		 *            estimateOrientationAt(int,int)\endlink
		 *            function and the result is kept stored by the instance.
		 *            Otherwise, if already computed, the function returns a
		 *            reference to the stored orientation.
		 *
		 *            The function calls the low-level
		 *            \link estimateOrientationAt(int,int)const
		 *            estimateOrientationAt(int,int)\endlink function
		 *            if an orientation is needed to be estimated.
		 *
		 * @param y
		 *            Ordinate varying between 0 and
		 *            \link getHeight()\endlink-1.
		 *
		 * @param x
		 *            Abscissa varying between 0 and
		 *            \link getWidth()\endlink-1.
		 *
		 * @return
		 *            %Orientation at <i>(x,y)</i>.
		 *
		 * @warning
		 *            If <i>(x,y)</i> does not lie between the fingerprint
		 *            image's region, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		const Orientation & getOrientationAt( int y , int x );

		/**
		 * @brief
		 *            Estimates the foreground of the fingerprint and
		 *            stores the result in the specified array.
		 *
		 * @details
		 *            Assume that
		 *            <pre>
		 *             m = getHeight();
		 *             n = getWidth();
		 *            </pre>
		 *            is the height and width of the fingerprint,
		 *            respectively, and
		 *            <pre>
		 *             bool *img = ...;
		 *            </pre>
		 *            is allocated to hold 'm*n' values of type
		 *            <code>bool</code> to which the output of the method
		 *            is stored via
		 *            <pre>
		 *             estimateForegroundImage(img)
		 *            </pre>
 	 	 *            encoded as an <i>m*n</i> matrix in row-major order.
 	 	 *            Specifically, the boolean indicating whether the pixel
 	 	 *            <i>(x,y)</i> (where <i>x=0...m-1, y=0...n-1</i>) belongs
 	 	 *            to the foreground can be accessed via
		 *            <pre>
		 *             bool isForegroundPixel = img[y*n+x];
		 *            </pre>
		 *
		 *
		 * @param foregroundImage
		 *            Array that can hold
		 *            <code>getWidth()*getHeight()</code> values of the
		 *            type <code>bool</code>.
		 *
		 * @warning
		 *            If <code>foregroundImage</code> cannot store
		 *            <code>getWidth()*getHeight()</code> boolean
		 *            values, the method runs into undocumented behavior.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		virtual void estimateForegroundImage( bool *foregroundImage ) const;

        /**
         * @brief
         *            Controls the number of gradients used to estimate the
         *            local orientations of the fingerprint with the
         *            gradient method.
         *
         * @details
         *            Unless overwritten by an extension class,
         *            this function is used for estimating local ridge
         *            orientations of the fingerprint with the gradient
         *            method, i.e.,
         *            by \link estimateOrientationAt()\endlink.
         *
         * @return
         *            The function returns
         *            '(int)ceil(8.0/569.0*getResolution())'
         */
        virtual int getGradientOrientationDistance() const;

		/**
		 * @brief
		 *            Estimates the orientation of the fingerprint at the
		 *            specified pixel.
		 *
		 * @details
		 *            The orientation is estimated using
		 * \link Orientation::gradient(int,int,const double*,int,int,int)\endlink
		 *            function.
		 *
		 * @param y
		 *            Ordinate varying between 0 and
		 *            \link getHeight()\endlink-1.
		 *
		 * @param x
		 *            Abscissa varying between 0 and
		 *            \link getWidth()\endlink-1.
		 *
		 * @return
		 *            %Orientation at <i>(x,y)</i>.
		 *
		 * @warning
		 *            If <i>(x,y)</i> does not lie between the fingerprint
		 *            image's region, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 */
		virtual Orientation estimateOrientationAt( int y , int x ) const;

		/**
		 * @brief
		 *            Estimates a directed reference point from this
		 *            fingerprint and returns the result
		 *
		 * @details
		 *            The implementation of the directed reference point
		 *            estimation is performed using the
		 *            \link TentedArchModel\endlink class following the
		 *            procedure as described in
		 *            <ul>
		 *             <li><b>B. Tams (2013)</b>.
		 *             Absolute %Fingerprint Pre-Alignment in Minutiae-Based
		 *             Cryptosystems. <i>Proc. BIOSIG 2013, ser LNI,
		 *             <b>vol. 212</b>, pp. 75--86</i>.</li>
		 *            </ul>
		 *            but can be overridden for a subclass, e.g., by
		 *            a more robust directed reference point estimation.
		 *            <br><br>
		 *            If the result of estimation is invalid, at least one
		 *            of the returned directed reference point's
		 *            <pre>
		 *             DirectedPoint drp = ...
		 *            </pre>
		 *            members
		 *            \link DirectedPoint::x drp.x\endlink,
		 *            \link DirectedPoint::y drp.y\endlink or
		 *            \link Direction::getDirectionAngle()
		 *                  drp.direction.getDirectionAngle()\endlink
		 *            should be equals NAN, i.e., self-comparision results
		 *            in a false boolean expression.
		 *            <br><br>
		 *            Overriding this method affects the behavior (but not
		 *            the meaning) of the
		 *            \link hasDirectedReferencePoint()\endlink and
		 *            \link getDirectedReferencePoint()\endlink function.
		 *
		 * @return
		 *            A reference point attached with a direction.
		 *
		 * @see hasDirectedReferencePoint()
		 * @see getDirectedReferencePoint()
		 *
		 * @warning
		 *            If not enough memory could be provided, an error message
		 *            is printed to stderr and the program exits with status
		 *            'EXIT_FAILURE'.
		 */
		virtual DirectedPoint estimateDirectedReferencePoint();

		/**
		 * @brief
		 *            Estimates a minutiae template from this fingerprint
		 *            (Note: this function is not yet implemented).
		 *
		 * @details
		 *            Currently, the implementation of this function does
		 *            nothing else as returning an empty minutiae view.
		 *
		 *            In order to allow programs using the THIMBLE library to
		 *            implement this function, it is provided as a virtual
		 *            function which can be overridden, thereby affecting the
		 *            \link getMinutiaeRecord()\endlink and
		 *            \link getMinutiaeView()\endlink functions.
		 *
		 * @return
		 *            A minutiae view containing the minutiae estimated from
		 *            this fingerprint.
		 *
		 * @see getMinutiaeView()
		 * @see getMinutiaeRecord()
		 *
		 * @attention
		 *            This function is not yet implemented and does nothing
		 *            than returning an empty set of minutiae thereby causing
		 *            the \link getMinutiaeView()\endlink function to return
		 *            a set of empty minutiae (unless initialized manually).
		 */
		virtual MinutiaeView estimateMinutiae();

		/**
		 * @brief
		 *            Specifies the foreground image manually.
		 *
		 * @details
		 *            The specified <code>foregroundImage</code> should be
		 *            encoded as a
		 *  <code>\link getWidth()\endlink*\link getHeight()\endlink</code>
		 *            matrix in row-major order. More specifically, the value
		 *       at <code>foregroundImage[y*\link getWidth()\endlink+x]</code>
		 *            is <code>true</code> if the pixel at <i>(x,y)</i> is
		 *            specified as the foreground and, otherwise, if specified
		 *            as the background the value of
		 *            <code>foregroundImage[y*getWidth()+x]</code> is
		 *            <code>false</code>.
		 *
		 * @param foregroundImage
		 *            Array containing <code>getWidth()*getHeight()</code>
		 *            boolean values.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setForegroundImage( const bool *foregroundImage );

		/**
		 * @brief
		 *            Specifies the orientation of the fingerprint at a
		 *            pixel manually.
		 *
		 * @param orientation
		 *            The orientation assigned to the fingerprint at pixel
		 *            <i>(x,y)</i>.
		 *
		 * @param y
		 *            Ordinate varying between 0 and
		 *            \link getHeight()\endlink-1.
		 *
		 * @param x
		 *            Abscissa varying between 0 and
		 *            \link getWidth()\endlink-1.
		 *
		 * @warning
		 *            If <i>(x,y)</i> does not lie between the fingerprint
		 *            image's region, an error message is printed to
		 *            <code>stderr</code> and the program exits with status
		 *            -1.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setOrientationAt
			( const Orientation & orientation , int y , int x );

		/**
		 * @brief
		 *            Access the directed reference point of the fingerprint.
		 *
		 * @details
		 *            If no attempt for directed reference point estimation
		 *            has been run for this fingerprint, the function causes
		 *            an attempt for directed reference point estimation
		 *            by calling the
		 *            \link estimateDirectedReferencePoint()\endlink function
		 *            and keeps track of the result as well as of whether the
		 *            estimation is considered to be successful, thereby
		 *            preventing that the estimation is run another time on
		 *            any further calls of the function.
		 *
		 *            Independently of whether the estimation was successful,
		 *            the function returns a reference to the directed
		 *            reference point but the result should be considered to
		 *            be invalid if the estimation failed.
		 *
		 *            Before the directed reference point is accessed, one
		 *            should test whether the estimation was successful or
		 *            not, for example, by
		 *            <pre>
		 *        %Fingerprint fingerprint = ...;
		 *
		 *        if ( fingerprint.hasDirectedReferencePoint() ) {
		 *           cout << fingerprint.getDirectedReferencePoint() << endl;
		 *        }
		 *            </pre>
		 *            which causes to write the directed reference point
		 *            representation to the standard output stream if the
		 *            estimation was successful.
		 *
		 *            Overriding the virtual
		 *            \link estimateDirectedReferencePoint()\endlink function
		 *            by a possibly more robust method for directed reference point
		 *            estimation may affect the behavior of this method.
		 *
		 * @return
		 *            A constant reference to the result of this fingerprint's
		 *            directed reference point estimation.
		 *
		 * @see estimateDirectedReferencePoint()
		 * @see hasDirectedReferencePoint()
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		const DirectedPoint & getDirectedReferencePoint();

		/**
		 * @brief
		 *            Checks whether directed reference point estimation was
		 *            successful for this fingerprint.
		 *
		 * @details
		 *            If for the fingerprint no directed reference point has
		 *            been estimated, the function calls the virtual
		 *            \link estimateDirectedReferencePoint()\endlink function,
		 *            keeps track of the result and of whether the directed
		 *            reference point is valid (i.e., whether none of its
		 *            members equals NAN) and returns <code>true</code> if
		 *            the directed reference point is valid and
		 *            <code>false</code> otherwise. Consequently, if an attempt
		 *            for directed reference point estimation has been
		 *            performed for this fingerprint before, the function
		 *            returns <code>true</code> if the directed reference
		 *            point of the earlier estimation was valid and
		 *            <code>false</code> otherwise.
		 *
		 * @return
		 *            <code>true</code> if directed reference point estimation
		 *            was successful and, otherwise, if the estimation failed,
		 *            the function returns <code>false</code>.
		 *
		 * @see estimateDirectedReferencePoint()
		 * @see getDirectedReferencePoint()
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool hasDirectedReferencePoint();

		/**
		 * @brief
		 *            Assigns the directed reference point manually.
		 *
		 * @details
		 *            The method assigns the directed point <code>dp</code> as
		 *            the fingerprint's directed reference point, thereby
		 *            replacing any former estimations or assignments of
		 *            the fingerprint's directed reference point.
		 *
		 *            When a directed reference point is specified manually,
		 *            it is assumed that it is valid, i.e., the function
		 *            \link hasDirectedReferencePoint()\endlink returns
		 *            <code>true</code> after manual directed reference point
		 *            assignment. In particular, the direction and pixel
		 *            coordinates of <code>dp</code> should hold well-defined
		 *            values. Otherwise, the use of this fingerprint instance
		 *            can run into undocumented behavior.
		 *
		 * @param dp
		 *            The valid directed reference point assigned to the
		 *            fingerprint.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setDirectedReferencePoint( const DirectedPoint & dp );

		/**
		 * @brief
		 *            Returns a minutiae record containing the minutiae
		 *            template estimated from this fingerprint.
		 *
		 * @details
		 *            If no attempt for minutiae estimation has been performed
		 *            for this fingerprint, this function creates a
		 *            \link MinutiaeRecord\endlink (with parameters specified
		 *            for this fingerprint) that contains an empty
		 *            \link MinutiaeView\endlink; then the view is assigned
		 *            with the result of the
		 *            \link estimateMinutiae\endlink function; finally, the
		 *            reference to the \link MinutiaeRecord\endlink is returned.
		 *            Otherwise, if a \link MinutiaeRecord\endlink has already
		 *            been created for this fingerprint on a former attempt for
		 *            minutiae estimation, this function returns the reference
		 *            to the already existend \link MinutiaeRecord\endlink.
		 *
		 * @return
		 *            A reference to the minutiae record containing a minutiae
		 *            template estimated from this fingerprint.
		 *
		 * @see estimateMinutiae()
		 * @see getMinutiaeView()
		 */
		const MinutiaeRecord & getMinutiaeRecord();

		/**
		 * @brief
		 *            Returns a minutiae template estimated from this
		 *            fingerprint.
		 *
		 * @details
		 *            If no attempt for minutiae estimation has been performed
		 *            for this fingerprint, this function creates a
		 *            \link MinutiaeRecord\endlink (with parameters specified
		 *            for this fingerprint) that contains an empty
		 *            \link MinutiaeView\endlink; then the view is assigned
		 *            with the result of the
		 *            \link estimateMinutiae\endlink function; finally, the
		 *            reference to the \link MinutiaeView\endlink of the
		 *            \link MinutiaeRecord\endlink is returned.
		 *            Otherwise, if a minutiae template estimation attempt
		 *            has already been performed for this fingerprint, this
		 *            function returns the reference to the already existent
		 *            \link MinutiaeView\endlink.
		 *
		 *            Currently, the implementation of the
		 *            \link estimateMinutiae()\endlink function does nothing
		 *            than returning an empty minutiae template. Consequently,
		 *            by merely using the functionalities provided by this
		 *            class, the user will obtain only empty minutiae
		 *            templates. However, the minutiae templates estimations
		 *            can be specified manually. Therefore, assume that a
		 *            minutiae record has been initialized from the file
		 *            'template.fmr' by
		 *            <pre>
		 *             MinutiaeRecord record;
		 *
		 *             record.read("template.fmr");
		 *            </pre>
		 *            which is assumed to have at least one valid
		 *            \link MinutiaeView\endlink which can be accessed by
		 *            <pre>
		 *             record.getView(0)
		 *            </pre>
		 *            Then the \link MinutiaeView\endlink of this fingerprint
		 *            <pre>
		 *             %Fingerprint fingerprint;
		 *            </pre>
		 *            can be initialized manually via
		 *            <pre>
		 *             fingerprint.getMinutiaeView() = record.getView(0);
		 *            </pre>
		 *
		 * @return
		 *            A reference to the minutiae view containing the minutiae
		 *            set estimated from this fingerprint.
		 *
		 * @see estimateMinutiae()
		 * @see getMinutiaeRecord()
		 */
		MinutiaeView & getMinutiaeView();

		/**
		 * @brief
		 *            Specify the minutiae for this fingerprint manually.
		 *
		 * @param view
		 *            The minutiae view containing the new minutiae specified
		 *            for this fingerprint object.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		void setMinutiaeView( const MinutiaeView & view );

		/**
		 * @brief
		 *            Initialize the fingerprint from a PGM gray-scale
		 *            image file.
		 *
		 * @details
		 *            This function serves as a convenience procedure to
		 *            prevent the user from temporarily creating
		 *            \link GrayImage\endlink objects. More
		 *            specifically, calling the function is equivalent to the
		 *            following:
		 *            <pre>
		 *             %Fingerprint fingerprint = ...; // This instance
		 *             GrayImage pgmImage;
		 *             if ( pgmImage.read(pgmImagePath) ) {
		 *                fingerprint.initialize(pgmImage,dpi);
		 *                return true;
		 *             } else {
		 *                return false;
		 *             }
		 *            </pre>
		 *
		 * @param pgmImagePath
		 *            Specifies the path to the PGM gray-scale fingerprint
		 *            image.
		 *
		 * @param dpi
		 *            The resolution in dots per inch at which the fingerprint
		 *            image has been scanned.
		 *
		 * @return
		 *            <code>true</code> if the fingerprint has been
		 *            successfully initialized; otherwise, the function
		 *            returns <code>false</code>.
		 *
		 * @warning
		 *             If <code>dpi</code> is not greater than 0, an error message
		 *             is printed to <code>stderr</code> and the program exits
		 *             with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *            If not enough memory could be provided, an error
		 *            message is printed to <code>stderr</code> and the
		 *            program exits with status 'EXIT_FAILURE'.
		 */
		bool fromImageFile
		( const std::string & pgmImagePath , int dpi = 500 );
	};
}

#endif /* THIMBLE_FINGERPRINT_H_ */
