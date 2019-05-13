/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file MinutiaeRecord.h
 *
 * @brief
 *            Provides classes representing minutiae, minutiae templates, and
 *            minutiae records in ISO 19794-2:2005 format.
 *
 * @details
 *
 * @section sec_minutiae Minutiae Templates
 *
 * Assume that a fingerprint's minutiae template is stored in the file
 * <code>in.fmr</code> in ISO 19794-2:2005 format. Then the class
 * \link thimble::MinutiaeRecord MinutiaeRecord\endlink can be used as an
 * interface to read the content of the file for further use in a program.
 * This can be implemented as follows
 * <pre>
 *  MinutiaeRecord record;
 *
 *  if ( !record.read("in.fmr") ) {
 *     cout << "could not read from 'in.fmr'." << endl;
 *  }
 * </pre>
 * Similarly, the record can written in a file <code>out.fmr</code>
 * via
 * <pre>
 *  if ( !record.write("out.fmr") ) {
 *     cout << "could not write to 'out.fmr'." << endl;
 *  }
 * </pre>
 * A \link thimble::MinutiaeRecord MinutiaeRecord\endlink usually contains
 * (at least) one \link thimble::MinutiaeView MinutiaeView\endlink that can
 * be accessed using the
 * \link thimble::MinutiaeRecord::getView(int)
 * getView()\endlink function via
 * <pre>
 *  MinutiaeView view = record.getView();
 * </pre>
 * Each view essentially represents a minutiae template that has been
 * estimated from a fingerprint. The number of minutiae contained in the
 * template can be accessed using the
 * \link thimble::MinutiaeView::getMinutiaeCount()
 * getMinutiaeCount()\endlink function:
 * <pre>
 *  int n = view.getMinutiaeCount();
 * </pre>
 * The \link thimble::MinutiaeView::getMinutia() getMinutia(int)\endlink
 * function may be used to access the <i>i</i>th minutia where
 * <i>i=0,...,n-1</i>:
 * <pre>
 *  Minutia minutia = view.getMinutia(i);
 * </pre>
 * To access the abscissa coordinate at which the
 * \link thimble::Minutia Minutia\endlink is placed, one may run
 * <pre>
 *  int x = minutia.getX();
 * </pre>
 * Analogous, the ordinate coordinate can be accessed by
 * <pre>
 *  int y = minutia.getY();
 * </pre>
 * Finally, to obtain the minutia angle, as a value between \f$[0,2\pi)\f$,
 * we may run
 * <pre>
 *  double angle = minutia.getAngle();
 * </pre>
 *
 * @see thimble::MinutiaeRecord
 * @see thimble::MinutiaeView
 * @see thimble::Minutia
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_MINUTIAERECORD_H_
#define THIMBLE_MINUTIAERECORD_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief Enumerates all relevant types for a minutia.
	 */
	typedef enum {

		/**
		 * @brief Code for an unknown minutia type
		 */
		UNKNOWN_MINUTIA_TYPE     = -1 ,

		/**
		 * @brief Code for a minutia ending
		 */
		ENDING_MINUTIA_TYPE      =  1 ,

		/**
		 * @brief Code for a minutia bifurcation
		 */
		BIFURCATION_MINUTIA_TYPE =  3

	} MINUTIA_TYPE_T;

	/**
	 * @brief Enumerates relevant finger position.
	 *
	 * @details
	 *        The codes for finger positions are defined by the
	 *        ANSI/NIST-ITL 1-2000
	 *		  <i>Data Format for the Interchange of Fingerprint
	 *		  Information</i>.
	 *        There are plain finger codes ranging from 11 to 14 but these
	 *        codes are not relevant here.
	 */
	typedef enum {

		/**
		 * @brief Finger position code for an unknown finger which is 0
		 */
		UNKNOWN_FINGER = 0 ,

		/**
		 * @brief Finger position code for a right thumb which is 1
		 */
		RIGHT_THUMB = 1 ,

		/**
		 * @brief Finger position code for a right index finger which is 2
		 */
		RIGHT_INDEX_FINGER = 2 ,

		/**
		 * @brief Finger position code for a right middle finger which is 3
		 */
		RIGHT_MIDDLE_FINGER = 3 ,

		/**
		 * @brief Finger position code for a ring middle finger which is 4
		 */
		RIGHT_RING_FINGER = 4 ,

		/**
		 * @brief Finger position code for a right little finger which is 5
		 */
		RIGHT_LITTLE_FINGER = 5 ,

		/**
		 * @brief Finger position code for a left thumb which is 6
		 */
		LEFT_THUMB = 6 ,

		/**
		 * @brief Finger position code for a left index finger which is 7
		 */
		LEFT_INDEX_FINGER = 7 ,

		/**
		 * @brief Finger position code for a left middle finger which is 8
		 */
		LEFT_MIDDLE_FINGER = 8 ,

		/**
		 * @brief Finger position code for a left ring finger which is 9
		 */
		LEFT_RING_FINGER = 9 ,

		/**
		 * @brief Finger position code for a left little finger which is 10
		 */
		LEFT_LITTLE_FINGER = 10

	} FINGER_POSITION_T;

	/**
	 * @brief Enumerate codes for a finger view's impression type.
	 *
	 * @details
	 *        The impression type corresponds to the type of capture device
	 *        used to to scan the image from where the minutiae were extraced.
	 *        <br>
	 *        There are other codes between 4 and 7 which correspond to latent
	 *        impressions but which are not relevant here.
	 */
	typedef enum {

		/**
		 * Code for a live-scan plain impression which is 0
		 */
		LIVESCAN_PLAIN = 0 ,

		/**
		 * Code for a live-scan rolled impression which is 1
		 */
		LIVESCAN_ROLLED = 1 ,

		/**
		 * Code for a non-live-scan plain impression which is 2
		 */
		NONLIVESCAN_PLAIN = 2 ,

		/**
		 * Code for a non-live-scan rolled impression which is 3
		 */
		NONLIVESCAN_ROLLED = 3 ,

		/**
		 * Code for a swipe impression which is 8
		 */
		SWIPE = 8

	} FINGER_IMPRESSION_TYPE_T;

	/**
	 * @brief
	 *      A class that represents a single minutia.
	 *
	 * @details
	 *      An instance contains a minutia's position, its angle, a measure
	 *      of its quality, and its determined type.
	 *
	 * @author Benjamin Tams
	 */
	class THIMBLE_DLL Minutia {

		friend class MinutiaeView;
		friend class MinutiaeRecord;

	private:

		/**
		 * @brief Horizontal position of the minutia
		 */
		double x;

		/**
		 * @brief Vertical position of the minutia
		 */
		double y;

		/**
		 * @brief
		 *      The angle of the minutia.
		 * @details
		 *      This is a value between 0.0 and (excluding)
		 *      <code>2.0*PI</code>.
		 */
		double angle;

		/**
		 * @brief
		 *      Encodes the type of the minutia.
		 * @details
		 * 		The type of a minutia can be either an ending, a bifurcation,
		 * 		or unknown which are encoded by the constant values
		 *      \link ENDING_MINUTIA_TYPE\endlink,
		 *      \link BIFURCATION_MINUTIA_TYPE\endlink, or
		 *      \link UNKNOWN_MINUTIA_TYPE\endlink, respectively.
		 */
		MINUTIA_TYPE_T typ;

		/**
		 * @brief The quality the minutia is obeyed with.
		 *
		 * @details This is a value between 1 and 100 where 1 corresponds
		 *          to the lowest and 100 to the highest quality. In case there
		 *          is was no quality determined for this minutiae, the value
		 *          will be 0.
		 */
		int quality;

		/**
		 * @brief Initializes this minutia as specified.
		 *
		 * @warning
		 *      If <code>quality</code> is not between 0 and 100
		 *      then the program outputs an error message on
		 *      <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @param x
		 *          horizontal position
		 * @param y
		 *          vertical position
		 * @param angle
		 *          angle of the minutia
		 * @param typ
		 *          type of the minutia
		 * @param quality
		 *          quality of the minutia which should lay between 1
		 *          and 100 and be 0 only if no quality was
		 *          determined for this minutia.
		 */
		void init
		( double x , double y , double angle , MINUTIA_TYPE_T typ ,
		  int quality );

	public:

		/**
		 * @brief
		 *             Standard constructor.
		 */
		inline Minutia() {
			init(0.0,0.0,0.0,UNKNOWN_MINUTIA_TYPE,0);
		}

		/**
		 * @brief
		 *      Creates a minutia at the specified position and that is
		 *      of the specified angle.
		 *
		 * @details
		 *      If <code>angle</code> is outside the range 0.0 and
		 *      <code>2*PI</code> then this minutia will be set with the
		 *      corresponding representative (modulo <code>2*PI</code>) of
		 *      the specified angle.
		 *
		 * @param x
		 *          horizontal position of the minutia
		 * @param y
		 *          vertical position of the minutia
		 * @param angle
		 *          the minutia's angle between 0.0 and
		 *          (excluding) <code>2*PI</code>
		 */
		inline Minutia
		( double x , double y , double angle ) {
			init(x,y,angle,UNKNOWN_MINUTIA_TYPE,0);
		}

		/**
		 * @brief
		 *      Creates a minutia of specified position, angle and type.
		 *
		 * @details
		 *      If <code>angle</code> is outside the range 0.0 and
		 *      <code>2*PI</code> then this minutia will be set with the
		 *      corresponding representative (modulo <code>2*PI</code>) of
		 *      the specified angle.
		 *
		 * @param x
		 *          horizontal position of the minutia
		 * @param y
		 *          vertical position of the minutia
		 * @param angle
		 *          the minutia's angle between 0.0 and
		 *          (excluding) <code>2*PI</code>
		 * @param typ
		 *          the minutia's type which is either
		 *          \link UNKNOWN_MINUTIA_TYPE\endlink,
		 *          \link ENDING_MINUTIA_TYPE\endlink, or
		 *          \link BIFURCATION_MINUTIA_TYPE\endlink
		 */
		inline Minutia
		( double x , double y , double angle , MINUTIA_TYPE_T typ ) {
			init(x,y,angle,typ,0);
		}

		/**
		 * @brief
		 *      Creates a minutia of specified position, angle, type, and
		 *      quality.
		 *
		 * @details
		 *      If <code>angle</code> is outside the range 0.0 and
		 *      <code>2*PI</code> then this minutia will be set with the
		 *      corresponding representative (modulo <code>2*PI</code>) of
		 *      the specified angle.
		 *
		 * @warning
		 *      If <code>quality</code> is not between 0 and 100
		 *      then the program outputs an error message on
		 *      <code>stderr</code> and exits with status 'EXIT_FAILURE'.
		 *
		 * @param x
		 *          horizontal position of the minutia
		 * @param y
		 *          vertical position of the minutia
		 * @param angle
		 *          the minutia's angle between 0.0 and
		 *          (excluding) <code>2*PI</code>
		 * @param typ
		 *          the minutia's type which is either
		 *          \link UNKNOWN_MINUTIA_TYPE\endlink,
		 *          \link ENDING_MINUTIA_TYPE\endlink, or
		 *          \link BIFURCATION_MINUTIA_TYPE\endlink
		 * @param quality
		 *          the minutia's quality which should lay between 1 and 100
		 *          and be 0 only if there was no quality measured for this
		 *          minutia.
		 */
		inline Minutia
		( double x , double y , double angle , MINUTIA_TYPE_T typ , int quality ) {
			init(x,y,angle,typ,quality);
		}

		/**
		 * @brief
		 *           Assignment operator of the \link Minutia\endlink class.
		 *
		 * @details
		 *           Sets this minutia to a copy of the given
		 *           <code>minutia</code>.
		 *
		 * @param minutia
		 *           The minutia of which this minutia will become a copy of.
		 */
		inline Minutia &operator=( const Minutia & minutia ) {
			this->x       = minutia.x;
			this->y       = minutia.y;
			this->angle   = minutia.angle;
			this->typ     = minutia.typ;
			this->quality = minutia.quality;
			return *this;
		}

		/**
		 * @brief
		 *          Copy constructor of the \link Minutia\endlink class
		 *
		 * @details
		 *          Creates an instance of the \link Minutia\endlink class that
		 *          is a copy of the given <code>minutia</code>.
		 *
		 * @param minutia
		 *          The minutia of which the copy is created by the
		 *          constructor
		 */
		inline Minutia( const Minutia & minutia ) {
			*this = minutia;
		}


		/**
		 * @brief Accesses the minutia's horizontal position
		 *
		 * @return horizontal position of the minutia
		 */
		inline double getX() const { return this->x; }

		/**
		 * @brief Accesses the minutia's vertical position
		 *
		 * @return vertical position of the minutia
		 */
		inline double getY() const { return this->y; }

		/**
		 * @brief Accesses the minutia's angle.
		 *
		 * @details
		 *       The angle is a value between 0.0 and (excluding)
		 *       <code>2.0*PI</code>
		 *
		 * @return angle of the minutia
		 */
		inline double getAngle() const { return this->angle; }

		/**
		 * @brief Accesses the type of the minutia
		 *
		 * @details
		 *        Possible results are
		 *        \link UNKNOWN_MINUTIA_TYPE\endlink (-1),
		 *        \link ENDING_MINUTIA_TYPE\endlink (1), and
		 *        \link BIFURCATION_MINUTIA_TYPE\endlink (3)
		 *
		 * @return the type of the minutia
		 */
		inline MINUTIA_TYPE_T getType() const { return this->typ; }

		/**
		 * @brief Accesses the quality of the minutia.
		 *
		 * @details
		 *       The quality is a value between 1 and 100 where
		 *       1 specifies the lowest and 100 the highest
		 *       quality. If no quality is associated with the minutia
		 *       the result will be 0.
		 *
		 * @return the measured quality of the minutia
		 */
		inline int getQuality() const { return this->quality; }
	};

	/**
	 * @brief
	 *           Checks whether two minutiae are equal.
	 *
	 * @param a
	 *           First minutia.
	 *
	 * @param b
	 *           Second minutia.
	 *
	 * @return
	 *           <code>true</code> if <code>a</code> has equal members
	 *           as <code>b</code> and otherwise <code>false</code>.
	 */
	inline bool operator==( const Minutia & a , const Minutia & b ) {
		if ( a.getX() != b.getX() ) {
			return false;
		}
		if ( a.getY() != b.getY() ) {
			return false;
		}
		if ( a.getAngle() != b.getAngle() ) {
			return false;
		}
		if ( a.getType() != b.getType() ) {
			return false;
		}
		if ( a.getQuality() != b.getQuality() ) {
			return false;
		}
		return true;
	}

	/**
	 * @brief
	 *           Checks whether two minutiae are not equal.
	 *
	 * @param a
	 *           First minutia.
	 *
	 * @param b
	 *           Second minutia.
	 *
	 * @return
	 *           <code>false</code> if <code>a</code> has equal members
	 *           as <code>b</code> and otherwise <code>true</code>.
	 */
	inline bool operator!=( const Minutia & a , const Minutia & b ) {
		return !(a==b);
	}

	/**
	 * @brief Prints a text representation of the given minutia to an output
	 *        stream.
	 *
	 * @param out
	 *         The output stream to where the text representation of
	 *         <code>minutia</code> is printed
	 * @param minutia
	 *         The minutia which text representation is printed to the
	 *         output stream
	 *
	 * @return <code>out</code>
	 */
	THIMBLE_DLL std::ostream & operator<<
			( std::ostream & out , const Minutia & minutia );


}

#ifdef THIMBLE_BUILD_DLL
template class THIMBLE_DLL std::allocator<thimble::Minutia>;
template class THIMBLE_DLL std::vector<thimble::Minutia,std::allocator<thimble::Minutia> >;
#endif

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief Represents a single view in a general
	 *        \link MinutiaeRecord\endlink.
	 *
	 * @details
	 *        Instances of this class list all measured minutiae of a finger
	 *        view. Furthermore, information about the finger's position, i.e.
	 *        either right/left thumb, right/left index finger, and so on are
	 *        encoded. Also, information about the finger's impression types
	 *        are encoded, e.g., if it is either a live-scan, non-live-scan,
	 *        or swipe. It also contains a value that measures the
	 *        fingerprint's overall quality. Finally, it contains a field of
	 *        the corresponding view number which stands for an index in
	 *        different record views.
	 *
	 * @author Benjamin Tams
	 */
	class THIMBLE_DLL MinutiaeView {

		friend class MinutiaeRecord;

	private:

		/**
		 * @brief The finger position of the finger this view corresponds to.
		 */
		FINGER_POSITION_T fingerPosition;

		/**
		 * @brief The index number of the view.
		 *
		 * @details
		 *        This is to distinguish multiple fingers in a record.
		 *        If only one view is recorded the view should have index 0.
		 */
		int viewNumber;

		/**
		 * The impression type of the view.
		 */
		FINGER_IMPRESSION_TYPE_T impressionType;

		/**
		 * @brief The quality of the view.
		 *
		 * @details
		 *       The value should be between 0 and 100 where 0 corresponds to
		 *       the lowest and 100 to the highest quality.
		 */
		int fingerQuality;

		/**
		 * @brief Lists all minutiae in the record view
		 */
		std::vector<Minutia> minutiae;

		/**
		 * @brief
		 *            Determines how much size in bytes this view would
		 *            require when stored.
		 *
		 * @return
		 *            Size in bytes this view needs on storage
		 */
		inline int getSizeInBytes() const {

			// Each minutia requires 6 bytes in the standard. Thus,
			// the total size in bytes is 6 times the number of minutiae,
			// plus the size of the header.

			return 6 * this->minutiae.size() + 4;
		}

		/**
		 * @brief
		 *            Writes the minutiae view to the specified output
		 *            stream.
		 *
		 * @details
		 *            The function attempts to write the minutiae view to the
		 *            specified output stream. If successfully written, the
		 *            function returns <code>true</code>; otherwise,
		 *            <code>false</code> is returned. The latter occurs if
		 *            an error occurred or if any of the minutiae coordinates
		 *            is negative or can not be recorded in 14 bits.
		 *
		 * @param out
		 *            The binary output stream to where the view is
		 *            written.
		 *
		 * @return
		 *            <code>true</code> if the view could be successfully
		 *            written; otherwise, e.g., if an error occurred, the
		 *            result is <code>false</code>.
		 *
		 * @warning
		 *            If <code>out</code> was not successfully opened as
		 *            a binary output stream before, the behavior of the
		 *            function is undocumented.
		 */
		bool write( FILE *out ) const;

		/**
		 * @brief Reads a minutiae view's data from the specified
		 *        <code>FILE</code> and stores the result in this view.
		 *
		 * @details
		 *       The data of the view will be overwritten during the
		 *       call.
		 *
		 * @param in
		 *       The <code>FILE</code> from where the view will be read.
		 *
		 * @return
		 *       If the view could be successfully read then the
		 *       function returns <code>true</code>; otherwise, if
		 *       an error occurred the result will be <code>false</code>
		 *
		 * @warning
		 *       If <code>in</code> was not successfully opened as a binary
		 *       input stream before, the behavior of the function is
		 *       undocumented.
		 */
		bool read( FILE *in );

		/**
		 * @brief Reads a minutiae view from the specified data
		 *        and overwrites the data of this view by the read data.
		 *
		 * @param data
		 *        The data of the minutiae view being read
		 *
		 * @return
		 *        If an error occurred, the result will be -1; otherwise,
		 *        the function returns the number of bytes read from
		 *        <code>dat</code>.
		 *
		 * @warning
		 *        If <code>dat</code> does not contain data of a valid
		 *        minutiae view, the function may run into undocumented
		 *        behavior.
		 */
		int read( const void *data );

	public:

		/**
		 * @brief Constructor for an empty view with default values.
		 *
		 * @details
		 *        It is constructed an empty view of unknown finger type
		 *        (\link UNKNOWN_FINGER\endlink) with view number 0
		 *        of a plain live scan (\link LIVESCAN_PLAIN\endlink).
		 *        Furthermore, the finger quality is set to the maximal valid
		 *        value 100 per default.
		 */
		inline MinutiaeView() {
			this->fingerPosition = UNKNOWN_FINGER;
			this->viewNumber     = 0;
			this->impressionType = LIVESCAN_PLAIN;
			this->fingerQuality  = 100;
		}

		/**
		 * @brief Constructor for an empty view with default values but
		 *        specified view number.
		 *
		 * @details
		 *        It is constructed an empty view of an unknown finger type
		 *        (\link UNKNOWN_FINGER\endlink)  of a plain live scan
		 *        (\link LIVESCAN_PLAIN\endlink). Furthermore, the finger
		 *        quality is set to the maximal valid value 100 per default.
		 *
		 * @warning
		 *        The view number should be >=0 and be representable as a
		 *        four bit integer
		 */
		inline MinutiaeView( int viewNumber ) {

			this->fingerPosition = UNKNOWN_FINGER;
			this->viewNumber     = viewNumber;
			this->impressionType = LIVESCAN_PLAIN;
			this->fingerQuality  = 100;
		}

		/**
		 * @brief Sets this view to a hard copy of the given <code>view</code>.
		 *
		 * @param view
		 *         The view of which this view will become a copy of.
		 */
		inline MinutiaeView &operator=( const MinutiaeView & view ) {
			this->fingerPosition = view.fingerPosition;
			this->viewNumber     = view.viewNumber;
			this->impressionType = view.impressionType;
			this->fingerQuality  = view.fingerQuality;
			this->minutiae       = view.minutiae;
			return *this;
		}

		/**
		 * @brief
		 *          Copy constructor of the \link MinutiaeView\endlink class
		 *
		 * @details
		 *          Creates an instance of the \link MinutiaeView\endlink class
		 *          that is a hard copy of the given <code>view</code>.
		 *
		 * @param view
		 *          The view of which the copy is created by the
		 *          constructor
		 */
		inline MinutiaeView( const MinutiaeView & view ) {
			*this = view;
		}

		/**
		 * @brief Accesses the number of minutiae that are contained in the
		 *        view.
		 *
		 * @return
		 *       the number of minutiae contained in the view
		 */
		inline int getMinutiaeCount() const { return this->minutiae.size(); }

		/**
		 * @brief Accesses the <i>i</i>-th minutia of the view.
		 *
		 * @param i index of the <i>i</i>-th minutia
		 *
		 * @return constant reference to the <i>i</i>-th minutia in the view
		 *
		 * @warning
		 *         if <code>i<0</code> or
		 *         <code>i>=\link getMinutiaeCount()\endlink</code>
		 *         an error message is printed to <code>stderr</code> and
		 *         the program exits with status 'EXIT_FAILURE'.
		 */
		inline const Minutia & getMinutia( int i ) const {

			if ( i < 0 || i >= (int)this->minutiae.size() ) {
				std::cerr << "Accessed minutia is out of bounds" << std::endl;
				exit(EXIT_FAILURE);
			}

			return this->minutiae.at(i);
		}

		/**
		 * @brief
		 *          Removes all minutiae that are contained in the view.
		 *
		 * @details
		 *          After calling the method the result of
		 *          \link getMinutiaeCount()\endlink will be 0.
		 */
		inline void removeAllMinutiae() {
			this->minutiae.clear();
		}

		/**
		 * @brief
		 *          Adds a minutia to the view.
		 *
		 * @param minutia
		 *          The minutia that is added to the view.
		 */
		inline void addMinutia( const Minutia & minutia ) {
			this->minutiae.push_back(minutia);
		}

		/**
		 * @brief
		 *          Removes the <i>i</i>-th minutia from the view.
		 *
		 * @param i
		 *          The index of the minutia that is removed.
		 *
		 * @warning
		 *         if <code>i</code> is negative or if
		 *         <code>i>=\link getMinutiaeCount()\endlink</code>
		 *         an error message is printed to <code>stderr</code> and
		 *         the program exits with status 'EXIT_FAILURE'.
		 */
		inline void removeMinutia( int i ) {

			if ( i < 0 || i >= (int)this->minutiae.size() ) {
				std::cerr << "Accessed minutia is out of bounds" << std::endl;
				exit(EXIT_FAILURE);
			}

			this->minutiae.erase(this->minutiae.begin()+i);
		}

		/**
		 * @brief
		 *           Inserts a minutia at the specified position.
		 *
		 * @param i
		 *           The index in the view of the inserted minutia.
		 *
		 * @param minutia
		 *           The minutia that is inserted.
		 *
		 * @warning
		 *         if <code>i</code> is negative or if
		 *         <code>i>\link getMinutiaeCount()\endlink</code>
		 *         an error message is printed to <code>stderr</code> and
		 *         the program exits with status 'EXIT_FAILURE'.
		 */
		inline void insertMinutia( int i , const Minutia & minutia ) {

			if ( i < 0 || i > (int)this->minutiae.size() ) {
				std::cerr << "Index is out of bounds" << std::endl;
				exit(EXIT_FAILURE);
			}

			this->minutiae.insert(this->minutiae.begin()+i,minutia);
		}

		/**
		 * @brief
		 *            Reserves the specified capacity for the view to hold
		 *            minutiae.
		 *
		 * @param capacity
		 *            The capacity in terms of number of minutiae that is
		 *            reserved.
		 */
		inline void ensureCapacity( int capacity ) {
			if ( capacity > 0 ) {
				this->minutiae.reserve(capacity);
			}
		}

		/**
		 * @brief
		 *            Sorts the minutiae in the view w.r.t. their
		 *            associated quality.
		 *
		 * @details
		 *            The minutiae with the lowest index
		 *            (i.e. <code>getMinutia(0)</code>) will correspond
		 *            to the minutiae of highest quality while the minutiae
		 *            with highest index
		 *            (i.e. <code>getMinutia(getMinutiaCount()-1)</code>) will
		 *            correspond to the minutiae with worst quality in the
		 *            view.
		 */
		void sortWithRespectToMinutiaeQuality();

		/**
		 * @brief Accesses the finger position of the view.
		 *
		 * @return the view's finger position
		 */
		inline FINGER_POSITION_T getFingerPosition() const {
			return this->fingerPosition;
		}

		/**
		 * @brief Accesses the impression type of the view.
		 *
		 * @details The impression type of a view corresponds to the
		 *          type of capture device used to scan the image
		 *          from where the minutiae are extracted.
		 *
		 * @return the view's impression type
		 */
		inline FINGER_IMPRESSION_TYPE_T getImpressionType() const {
			return this->impressionType;
		}

		/**
		 * @brief Access the quality of the view.
		 *
		 * @details
		 *        Valid values for a finger quality are integers between
		 *        0 and 100 where 0 corresponds the lowest and 100 to the
		 *        highest quality.
		 *
		 * @return the overall quality of the view
		 *
		 * @see setFingerPosition(FINGER_POSITION_T)
		 */
		inline int getFingerQuality() const { return this->fingerQuality; }

		/**
		 * @brief Sets the finger position of the view.
		 *
		 * @param fingerPosition
		 *             the finger's position
		 *
		 * @see getFingerPosition()
		 */
		inline void setFingerPosition( FINGER_POSITION_T fingerPosition ) {
			this->fingerPosition = fingerPosition;
		}

		/**
		 * @brief    Sets the impression type of the finger view.
		 *
		 * @details
		 *           The impression type of a view corresponds to the
		 *           type of capture device used to scan the image
		 *           from where the minutiae are extracted.
		 *
		 * @param impressionType
		 *           the new impression type
		 */
		inline void setImpressionType
		( FINGER_IMPRESSION_TYPE_T impressionType ) {
			this->impressionType = impressionType;
		}

		/**
		 * @brief    Specify the quality of the view.
		 *
		 * @details
		 *           Valid qualities are values between 0 and 100.
		 *
		 * @param fingerQuality
		 *           the specified finger quality
		 *
		 * @warning
		 *           if <code>fingerQuality<0</code> or
		 *           <code>fingerQuality>100</code> an error message is
		 *           printed to <code>stderr</code> and the program exits
		 *           with status 'EXIT_FAILURE'.
		 */
		inline void setQuality( int fingerQuality ) {

			if ( fingerQuality < 0 || fingerQuality > 100 ) {
				std::cerr <<
					"The finger's quality must range between 0 and 100"
						 << std::endl;
				exit(EXIT_FAILURE);
			}

			this->fingerQuality = fingerQuality;
		}
	};

	/**
	 * @brief Prints a text representation of the given minutiae view to an
	 *        output stream.
	 *
	 * @param out
	 *         The output stream to where the text representation of
	 *         <code>view</code> is printed
	 * @param view
	 *         The minutiae view whose text representation is printed to the
	 *         output stream
	 *
	 * @return <code>out</code>
	 */
	THIMBLE_DLL std::ostream & operator<<
			( std::ostream & out , const MinutiaeView & view );

}

#ifdef THIMBLE_BUILD_DLL
template class THIMBLE_DLL std::allocator<thimble::MinutiaeView>;
template class THIMBLE_DLL std::vector<thimble::MinutiaeView,std::allocator<thimble::MinutiaeView> >;
#endif

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Provides functionalities for reading  and writing finger
	 *            minutiae data in ISO 19794-2:2005 format
	 */
	class THIMBLE_DLL MinutiaeRecord {

	private:

		/**
		 * @brief
		 *            Contains four bits to indicate that the capture
		 *            equipment used to capture the original image was
		 *            compliant with a standard certification method for
		 *            such an equipment. For further information we refer
		 *            to the specification of the standard.
		 */
		bool captureEquipmentCertifications[4];

		/**
		 * @brief
		 *            These twelve bits are used to identify the type or
		 *            model of the capture device used to acquire the original
		 *            biometric sample. Unreported capture devices are
		 *            indicated when all bits are set to <code>false</code>.
		 *            For further informations we refer to the specification
		 *            of the standard.
		 */
		bool captureDeviceTypeID[12];

		/**
		 * @brief
		 *            The horizontal size (width) of the fingerprint images
		 *            this record corresponds to.
		 */
		int sizeX;

		/**
		 * @brief
		 *            The vertical size (height) of the fingerprint images
		 *            this record corresponds to.
		 */
		int sizeY;

		/**
		 * @brief
		 *            The horizontal resolution of the fingerprint images
		 *            this record corresponds to. The resolution corresponds
		 *            to the number of pixels per centimeter (<i>ppc</i>).
		 */
		int resX;

		/**
		 * @brief
		 *            The vertical resolution of the fingerprint images this
		 *            record corresponds to. The resolution corresponds to the
		 *            number of pixels per centimeter (<i>ppc</i>).
		 */
		int resY;

		/**
		 * @brief
		 *            This field is reserved for future use, and to align
		 *            the end of the record header on a long-word (four byte)
		 *            boundary. For the 2005 version of the standard, this
		 *            field shall be set to zero.
		 *            Also, see the specification of the standard.
		 */
		int reservedByte;

		/**
		 * @brief
		 *           Lists all record views that are associated with this
		 *           record.
		 */
		std::vector<MinutiaeView> views;

		/**
		 * @brief
		 *        Initializes this record with default values
		 *
		 * @details
		 *        In particular every entry of the fields
		 *        \link captureEquipmentCertifications\endlink and
		 *        \link captureDeviceTypeID\endlink
		 *        are set to <code>false</code>.
		 *        <br>
		 *        The fingerprint image's dimension are set to
		 *        <code>\link sizeX\endlink=0</code> and
		 *        <code>\link sizeY\endlink=0</code>.
		 *        <br>
		 *        The resolutions are set to
		 *        <code>\link resX\endlink=197</code> and
		 *        <code>\link resY\endlink=197</code> pixels per centimeter.
		 *        <br>
		 *        The field <code>\link reservedByte\endlink</code>
		 *        is set to 0.
		 *        <br>
		 *        Finally, all views in the record are removed, i.e.
		 *        <code>\link views\endlink.clear()</code>
		 *
		 */
		void init();


	public:

		/**
		 * @brief
		 *          Standard constructor.
		 *
		 * @details
		 *          Creates an empty minutiae record for fingerprints
		 *          of zero width and height corresponding to a resolution
		 *          of 197 pixels per centimeter (i.e., 500 dots per inch).
		 *
		 *          After an empty minutiae record has been created its
		 *          parameters its content should be initialized from an
		 *          external file or data in ISO 19794-2:2005 format using one
		 *          of its read members.
		 *
		 * @see read(FILE*)
		 * @see read(const void*)
		 * @see read(const std::string&)
		 */
		MinutiaeRecord();

		/**
		 * @brief Constructs an empty minutiae record with specified width,
		 *        height, horizontal-, and vertical resolution.
		 *
		 * @param width
		 *            determines in which range the horizontal position of record
		 *            minutiae are allowed
		 * @param height
		 *            determines in which range the vertical position of record
		 *            minutiae are allowed
		 * @param ppcX
		 *            indicates the horizontal resolution in pixels per centimeter
		 * @param ppcY
		 *            indicates the vertical resolution in pixels per centimeter
		 * @warning
		 *            The constructor causes the print of an error message
		 *            to <code>stderr</code> and an exit with status 'EXIT_FAILURE'
		 *            if <code>width</code>, <code>height</code>,
		 *            <code>ppcX</code>, or <code>ppcY</code> is smaller
		 *            than or equal 0 or if they can not be encoded by two
		 *            bytes.
		 */
		MinutiaeRecord(int width, int height, int ppcX, int ppcY);

		/**
		 * @brief Constructs an empty minutiae record with specified width,
		 *        height, and resolution.
		 *
		 * @param width
		 *            determines in which range the horizontal position of record
		 *            minutiae are allowed
		 * @param height
		 *            determines in which range the vertical position of record
		 *            minutiae are allowed
		 * @param ppc
		 *            indicates the resolution in pixels per centimeter
		 *            (horizontal and vertical)
		 *
		 * @warning
		 *            The constructor causes the print of an error message
		 *            to <code>stderr</code> and an exit with status 'EXIT_FAILURE'
		 *            if <code>width</code>, <code>height</code>, or
		 *            <code>ppc</code> is smaller than or equal 0 or if they
		 *            can not be encoded by two bytes.
		 */
		MinutiaeRecord(int width, int height, int ppc);

		/**
		 * @brief Constructs an empty minutiae record with specified width and
		 *        height.
		 *
		 * @details
		 *        The horizontal and vertical resolution are both set as 197
		 *        pixels per centimeter.
		 *
		 * @param width
		 *            determines in which range the horizontal position of
		 *            record minutiae are allowed
		 * @param height
		 *            determines in which range the vertical position of record
		 *            minutiae are allowed
		 *
		 * @warning
		 *            The constructor causes the print of an error message
		 *            to <code>stderr</code> and an exit with status 'EXIT_FAILURE'
		 *            if <code>width</code> or <code>height</code> is smaller
		 *            than or equal 0 or if they can not be encoded by two
		 *            bytes.
		 */
		MinutiaeRecord(int width, int height);


		/**
		 * @brief Sets this record to a hard copy of the given
		 *        <code>record</code>.
		 *
		 * @param record
		 *         The record of which this record will become a copy of.
		 */
		MinutiaeRecord &operator=( const MinutiaeRecord & record );

		/**
		 * @brief
		 *          Copy constructor.
		 *
		 * @details
		 *          Creates an instance of the \link MinutiaeRecord\endlink class
		 *          that is a hard copy of the given <code>record</code>.
		 *
		 * @param record
		 *          The record of which the copy is created by the
		 *          constructor
		 */
		inline MinutiaeRecord( const MinutiaeRecord & record ) {
			*this = record;
		}

		/**
		 * @brief Accesses the horizontal size in where record minutiae can
		 *        range.
		 *
		 * @return horizontal size
		 */
		inline int getWidth() const { return this->sizeX; }

		/**
		 * @brief Accesses the vertical size in where record minutiae can
		 *        range.
		 *
		 * @return vertical size
		 */
		inline int getHeight() const { return this->sizeY; }

		/**
		 * @brief Accesses the horizontal resolution associated with this
		 *        record.
		 *
		 * @details
		 *        The result is the horizontal resolution in pixels per
		 *        centimeter of the capture device of the fingerprint image
		 *        from where the record's minutiae have been extracted.
		 *
		 * @return horizontal resolution
		 */
		inline int getHorizontalResolution() const { return this->resX; }

		/**
		 * @brief Accesses the vertical resolution associated with this
		 *        record.
		 *
		 * @details
		 *        The result is the vertical resolution in pixels per
		 *        centimeter of the capture device of the fingerprint image
		 *        from where the record's minutiae have been extracted.
		 *
		 * @return vertical resolution
		 */
		inline int getVerticalResolution() const { return this->resY; }

		/**
		 * @brief Accesses the number of views in the record.
		 *
		 * @return the number of views in the record
		 */
		inline int getViewCount() const { return this->views.size(); }

		/**
		 * @brief
		 *             Access the view in the record of given view number.
		 *
		 * @param viewNumber
		 *             view number of the accessed view
		 * @return
		 *             the view with indicated view number or <code>null</code>
		 *             if no view with <code>viewNumber</code> is contained
		 *             in the record
		 *
		 * @warning
		 *             if no view with the specified view number exists
		 *             in the record an error is printed to
		 *             <code>stderr</code> and the program exits with status
		 *             -1.
		 */
		MinutiaeView & getView( int viewNumber = 0 );

		/**
		 * @brief
		 *             Access the view in the record of given view number
		 *             (constant version).
		 *
		 * @param viewNumber
		 *             view number of the accessed view
		 * @return
		 *             the view with indicated view number or <code>null</code>
		 *             if no view with <code>viewNumber</code> is contained
		 *             in the record
		 *
		 * @warning
		 *             if no view with the specified view number exists
		 *             in the record an error is printed to
		 *             <code>stderr</code> and the program exits with status
		 *             -1.
		 */
		const MinutiaeView & getView( int viewNumber ) const;

		/**
		 * @brief
		 *             Adds the specified view to the minutiae record.
		 *
		 * @details
		 *             The added view's view number may change such that
		 *             it indicates its index position in the record's
		 *             view list.
		 *
		 * @param view
		 *             The view that is added to the record.
		 */
		inline void addView( const MinutiaeView & view = MinutiaeView() ) {
			this->views.push_back(view);
			this->views.back().viewNumber = (int)(this->views.size()-1);
		}

		/**
		 * @brief
		 *            Determine the number of bytes required to store
		 *            the record.
		 *
		 * @return
		 *            The number of bytes required to store the record.
		 */
		int getSizeInBytes() const;

		/**
		 * @brief
		 *            Writes the minutiae record to the specified binary
		 *            output stream.
		 *
		 * @details
		 *            The funtion attempts to write the minutiae record to the
		 *            specified output stream. If successfully written, the
		 *            function returns <code>true</code> and
		 *            <code>false</code> otherwise. The latter may occur if an
		 *            error occurred or if any of the minutiae coordinates is
		 *            negative or can not be reocrded in 14 bits.
		 *
		 * @param out
		 *            The binary output stream to where the record is written.
		 *
		 * @return
		 *            <code>true</code> if the view could be successfully
		 *            written; otherwise, the result is <code>false</code>.
		 *
		 * @warning
		 *            If <code>out</code> was not successfully opened as a
		 *            binary output stream before, the behavior of the
		 *            function is undefined.
		 */
		bool write( FILE *out ) const;

		/**
		 * @brief
		 *            Writes the minutiae record to the file specified by
		 *            the given file name.
		 *
		 * @details
		 *            The function is a wrapper around
		 *            <code>write(FILE*) const</code> and forwards its
		 *            return result.
		 *
		 * @param fileName
		 *            The file name of the file to where the record is
		 *            written.
		 *
		 * @return
		 *            <code>true</code> if the view could be successfully
		 *            written; otherwise, the result is <code>false</code>.
		 *
		 * @see write(FILE*) const
		 */
		bool write( const std::string & fileName ) const;

		/**
		 * @brief    Reads a minutiae record from the specified
		 *           <code>FILE</code> and returns the result.
		 *
		 * @details
		 *          Attempts to read a minutiae record from <code>in</code>
		 *          which must be a valid non-null <code>FILE</code>
		 *          which was opened before with the reading flag. If
		 *          a valid minutiae record has been read from <code>in</code>,
		 *          the function returns <code>true</code> and this instance
		 *          represents the data read from <code>in</code>; otherwise
		 *          the function returns <code>false</code> and this instance
		 *          will be left unchanged.
		 *
		 * @param in
		 *          The <code>FILE</code> from where to read the record
		 *
		 * @return
		 *          <code>true</code> if the minutiae record has been
		 *          successfully read from <code>in</code>; otherwise, the
		 *          function returns <code>false</code>.
		 */
		bool read( FILE *in );

		/**
		 * @brief    Reads a minutiae record from the specified
		 *           data array.
		 *
		 * @param data
		 *           Contains the data of the specified minutiae record
		 *
		 * @return
		 *          <code>true</code> if the minutiae record has been
		 *          successfully read from <code>data</code>; otherwise, the
		 *          function returns <code>false</code>.
		 */
		bool fromBytes( const void *data );

		/**
		 * @brief    Reads a minutiae record from a file specified
		 *           by its path.
		 *
		 * @details
		 *           This function is a wrapper around
		 *           \link read(FILE*)\endlink. It opens a
		 *           <code>FILE</code> from the path given by
		 *           <code>fileName</code> and attempts to read the
		 *           corresponding minutiae record. Afterwards the opened
		 *           <code>FILE</code> is closed. If the opened
		 *           <code>FILE</code> was <code>NULL</code> or if the
		 *           call of \link read(FILE*)\endlink caused an error
		 *           the function returns <code>false</code> and leaves
		 *           this instance unchanged; otherwise, this instance
		 *           data is updated by the read data and the function
		 *           returns <code>true</code>.
		 *
		 * @param fileName
		 *          The file name from where to read the record
		 *
		 * @return
		 *          <code>true</code> if the minutiae record has been
		 *          successfully read from <code>in</code>; otherwise, the
		 *          function returns <code>false</code>.
		 */
		bool read( const std::string & fileName );
	};
}

#endif /* THIMBLE_MINUTIAERECORD_H_ */
