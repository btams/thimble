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
 * @file Labeler.h
 *
 * @brief
 *            Provides a class for dividing two-dimensional binary images into
 *            their connected components.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_LABELER_H_
#define THIMBLE_LABELER_H_

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Enumerates codes encoding the connectivity type
	 *            in the context of a \link thimble::Labeler Labeler\endlink.
	 */
	typedef enum {

		/**
		 * @brief
		 *            Code indicating the four-connectivity.
		 */
		LC_FOUR_CONNECTED ,


		/**
		 * @brief
		 *            Code indicating eight-connectivity.
		 */
		LC_EIGHT_CONNECTED

	} LABEL_CONNECTIVITY_T;

	/**
	 * @brief
	 *            An instance of this class can be used to divide a
	 *            two-dimensional binary image into their connected
	 *            components.
	 */
	class THIMBLE_DLL Labeler {

	public:

		/**
		 * @brief
		 *            Creates a labeler object for dividing a two-dimensional
		 *            binary image into their connected components.
		 *
		 * @param connectivity
		 *            The connectivity type of the binary images.
		 */
		Labeler( LABEL_CONNECTIVITY_T connectivity = LC_EIGHT_CONNECTED );

		/**
		 * @brief
		 *            Destructor.
		 */
		~Labeler();

		/**
		 * @brief
		 *            Accesses the connectivity of this labeler.
		 *
		 * @return
		 *            The connectivity code with which this object has been
		 *            created.
		 */
		LABEL_CONNECTIVITY_T getConnectivity() const;

		/**
		 * @brief
		 *            Divides the specified binary image into their
		 *            connected black pixel components.
		 *
		 * @details
		 *            The result of the separation will be hold by this
		 *            object; any existing resulting hold by this object
		 *            will be overwritten.
		 *
		 * @param binaryImage
		 *            Contains <code>m*n</code> values of type
		 *            <code>bool</code> where if <code>(y,x)</code> is
		 *            a black pixel, then <code>binaryImage[y*n+x]</code>
		 *            equals <code>true</code> and, otherwise, if
		 *            <code>(y,x)</code> is white, then
		 *            <code>binaryImage[y*n+x]</code> equals
		 *            <code>false</code>.
		 *
		 * @param m
		 *            Height of the binary image.
		 *
		 * @param n
		 *            Width of the binary image.
		 *
		 * @warning
		 *            In the following cases, the behavior of this method
		 *            is undocumented:
		 *            <ul>
		 *             <li>
		 *              <code>binaryImage</code> does not contain
		 *              <code>m*n</code> valid <code>bool</code> values from
		 *              <code>{true,false}</code>.
		 *             </li>
		 *             <li><code>m</code> is negative</li>
		 *             <li><code>n</code> is negative</li>
		 *            </ul>
		 */
		void label( const bool *binaryImage , int m , int n );

		/**
		 * @brief
		 *            Returns the height of the latest binary image for
		 *            which this labeler has been applied.
		 *
		 * @return
		 *            The height of the latest binary image for which
		 *            this labeler has been applied.
		 */
		int getHeight() const;

		/**
		 * @brief
		 *            Returns the width of the latest binary image for
		 *            which this labeler has been applied.
		 *
		 * @return
		 *            The width of the latest binary image for which
		 *            this labeler has been applied.
		 */
		int getWidth() const;

		/**
		 * @brief
		 *            Accesses the number of connected components hold by this
		 *            labeler.
		 *
		 * @details
		 *            After a call of \link label()\endlink applied to a
		 *            non-empty image (i.e., an image that contains at least
		 *            one black pixel), the result of this function will be
		 *            greater than zero.
		 *
		 * @return
		 *            The number of connected components hold by this labeler.
		 */
		int getLabelCount() const;

		/**
		 * @brief
		 *            Draws the pixels of the specified component's
		 *            index to the specified binary image.
		 *
		 * @details
		 *            The array <code>binaryImage</code> of size at
		 *            least \link getHeight()\endlink*\link getWidth()\endlink
		 *            will have entries <code>false</code> except for those
		 *            pixels <code>(y,x)</code> that belong to the specified
		 *            component, i.e., the value
		 *            <code>binaryImage[y*\link getWidth()\endlink+x]</code>
		 *            equals <code>true</code> iff <code>(y,x)</code> belongs
		 *            to the <code>labelIndex</code>'th component.
		 *
		 *            The array <code>binaryImage</code> will contain a
		 *            non-empty image if <code>labelIndex=0,
		 *            ..., \link getLabelCount()\endlink-1</code>. If
		 *            <code>labelIndex=-1</code>, then
		 *            <code>binaryImage</code> forms the image of which pixel
		 *            are set <code>true</code> that belong to none
		 *            component. Otherwise, if <code>labelIndex <= -2</code>
		 *            or if <code>labelIndex >= getLabelCount()</code> the
		 *            image <code>binaryImage</code> will be empty.
		 *
		 * @param binaryImage
		 *            An array that can
		 *            hold \link getHeight()\endlink*\link getWidth()\endlink
		 *            values of type <code>bool</code> to which the pixels
		 *            of the specified components are drawn.
		 *
		 * @param labelIndex
		 *            The index of component of which pixels are drawn to
		 *            <code>binaryImage</code>.
		 */
		void get( bool *binaryImage , int labelIndex ) const;

		/**
		 * @brief
		 *            Counts the number of pixels that belong to the specified
		 *            component.
		 *
		 * @details
		 *            <b>NOTE:</b> A labeler object stores the connected
		 *            component's sorted w.r.t. their pixel size. More
		 *            precisely, the largest component is specified by the
		 *            label index 0 where the smallest non-empty component's
		 *            index equals \link getLabelCount()\endlink-1.
		 *
		 * @param labelIndex
		 *            The index of the component of which pixel size is
		 *            counted.
		 *
		 * @return
		 *            The number of pixels belonging to the
		 *            <code>labelIndex</code>'th component.
		 */
		int labelSize( int labelIndex ) const;

	private:

		/**
		 * @brief
		 *            The connectivity of this labeler.
		 *
		 * @details
		 *            The connectivity type, which can be
		 *            either \link LC_FOUR_CONNECTED four-connected\endlink
		 *            or \link LC_EIGHT_CONNECTED eight-connected\endlink, is
		 *            specified on construction of a \link Labeler\endlink
		 *            object via
		 *            the \link Labeler(LABEL_CONNECTIVITY_T)\endlink
		 *            constructor.
		 *
		 * @see getConnectivity()
		 * @see Labeler(LABEL_CONNECTIVITY_T)
		 */
		LABEL_CONNECTIVITY_T connectivity;

		/**
		 * @brief
		 *            The height of latest image to which this labeler
		 *            has been applied.
		 *
		 * @details
		 *            Will be 0 as long as this labeler object has not been
		 *            applied to an image.
		 *
		 * @see getHeight()
		 */
		int m;

		/**
		 * @brief
		 *            The width of latest image to which this labeler
		 *            has been applied.
		 *
		 * @details
		 *            Will be 0 as long as this labeler object has not been
		 *            applied to an image.
		 *
		 * @see getWidth()
		 */
		int n;

		/**
		 * @brief
		 *            Stores the components' indices at the pixel to which
		 *            the pixel belongs.
		 *
		 * @details
		 *            This array can store \link m\endlink * \link n\endlink
		 *            values of type <code>int</code>. If <code>(y,x)</code>
		 *            is a pixel that belongs to the <i>j</i>th label index,
		 *            then <code>labelImage[y*\link n\endlink+x]==j</code>.
		 *            Otherwise, if <code>(y,x)</code> does not belong to any
		 *            non-empty component, then
		 *            <code>labelImage[y*\link n\endlink+x]==-1</code>.
		 */
		int *labelImage;

		/**
		 * @brief
		 *            The number of different components stored by this
		 *            labeler object.
		 *
		 * @details
		 *            If <code>c</code> is the maximum over all
		 *            <code>labelImage[y*\link n\endlink+x]</code> where
		 *            <code>y=0,...,\link m\endlink-1</code> and
		 *            <code>y=0,...,\link n\endlink-1</code>, then
		 *            <code>labelCount=c+1</code>.
		 *
		 * @see getLabelCount()
		 */
		int labelCount;

		/**
		 * @brief
		 *            Initially iterates through all pixels of the specified
		 *            image the first time.
		 *
		 * @details
		 *            The function assumes that the specified binary image
		 *            is of dimension \link m\endlink * \link n\endlink and
		 *            that the array \link labelImage\endlink contains
		 *            \link m\endlink * \link n\endlink integers equals
		 *            -1.
		 *
		 *            After the function has terminated, the
		 *            members \link labelImage\endlink
		 *            and \link labelCount\endlink will be initialized such
		 *            that \link labelImage\endlink contains a separation
		 *            between \link labelCount\endlink components covering
		 *            all black pixels in <code>binaryImage</code> but which
		 *            do not need to be disjoint.
		 *
		 *            A full separation between the components can be
		 *            determined during subsequent repeated calls of
		 *            the \link passForward()\endlink
		 *            and \link passBackward()\endlink functions until
		 *            the state of the object remains unchanged.
		 *
		 * @return
		 *            <code>true</code> if the specified image is non-empty
		 *            and, otherwise, if the specified image does only
		 *            have <code>false</code> values, the result of this
		 *            function will be <code>false</code>.
		 */
		bool firstPass( const bool *binaryImage );

		/**
		 * @brief
		 *            Forward iteration through all pixels of the current
		 *            component separation in which each considered pixel's
		 *            label index is updated to the smallest index of its
		 *            neighboring pixels.
		 *
		 * @details
		 *            The function iterates through all pixels row-wise,
		 *            starting from (0,0) and ending at
		 *            (\link m\endlink-1,\link n\endlink-1).
		 *
		 * @return
		 *            <code>true</code> if the component separation was
		 *            changed during the iteration and, otherwise,
		 *            <code>false</code>.
		 */
		bool passForward();

		/**
		 * @brief
		 *            Backward iteration through all pixels of the current
		 *            component separation in which each considered pixel's
		 *            label index is updated to the smallest index of its
		 *            neighboring pixels.
		 *
		 * @details(\link m\endlink-1,\link n\endlink-1)
		 *            The function iterates through all pixels row-wise,
		 *            starting from (\link m\endlink-1,\link n\endlink-1) and
		 *            ending at (0,0).
		 *
		 * @return
		 *            <code>true</code> if the component separation was
		 *            changed during the iteration and, otherwise,
		 *            <code>false</code>.
		 */
		bool passBackward();

		/**
		 * @brief
		 *            Inserts a new pixel to the current connected
		 *            component separation.
		 *
		 * @details
		 *            This function is called by
		 *            the \link firstPass()\endlink function and
		 *            called for each black pixel.
		 *
		 * @param count
		 *            The number of different component indices contained
		 *            in the current connected component separation.
		 *
		 * @param binaryImage
		 *            Array that contains \link m\endlink * \link n\endlink
		 *            values of type <code>bool</code>.
		 *
		 * @param i
		 *            y-coordinate of the new pixel.
		 *
		 * @param j
		 *            x-coordinate of the new pixel.
		 *
		 * @return
		 *            The new number of different component indices
		 *            which can be either <code>count</code> or
		 *            <code>count+1</code>.
		 */
		int newLabel
			( int count , const bool *binaryImage , int i , int j ) const;

		/**
		 * @brief
		 *            Sets the pixel's label index to the minimum of its
		 *            neighboring labeled pixel's indices.
		 *
		 * @details
		 *            This function is used by the \link passForward()\endlink
		 *            and \link passBackward()\endlink functions and called
		 *            only if the specified pixel belongs to a component,
		 *            i.e., if the
		 *            entry \link labelImage\endlink[i*\link n\endlink+j] is
		 *            greater than or equals 0.
		 *
		 * @param i
		 *            y-coordinate of the re-labeled pixel.
		 *
		 * @param j
		 *            x-coordinate of the re-labeled pixel.
		 *
		 * @return
		 *            <code>true</code> if the specified pixel's label has
		 *            been updated and, otherwise, if the specified pixel's
		 *            label remains unchanged, the function returns
		 *            <code>false</code>.
		 */
		bool relabel( int i , int j );

		/**
		 * @brief
		 *            This function returns the smallest index of that
		 *            component being stored by this labeler object that is
		 *            greater than the specified component's label index.
		 *
		 * @details
		 *            After the function \link label()\endlink has finished
		 *            with separating a binary image's connected components,
		 *            the component's label indices do not need to cover the
		 *            elements of an integer
		 *            interval { 0 , 1 , 2 , ... , n-1 }. This function is
		 *            used in the \link label()\endlink function to iterate
		 *            through the index
		 *            sequence \f$\{\ell_0,...,\ell_{n-1}\}\f$ covering
		 *            the component's label indices where \f$\ell_i\geq i\f$.
		 *            In each iteration the  replacements of the
		 *            indices \f$\ell_i\leftarrow i\f$ is performed such that
		 *            the the indices fully cover the
		 *            interval { 0 , 1 , 2 , ... , n-1 }.
		 *
		 * @param labelIndex
		 *            The index of a labeled component.
		 *
		 * @return
		 *            The smallest index of that component hold by this
		 *            labeler that is strictly greater than the specified
		 *            label index if existent; otherwise, the function
		 *            returns 'INT_MAX'.
		 */
		int nextLabel( int labelIndex ) const;

		/**
		 * @brief
		 *            Sorts the disjoint components hold by this labeler
		 *            w.r.t. their pixel size.
		 *
		 * @details
		 *            This method assumes that a full separation between
		 *            a binary image's connected components has been performed
		 *            and that the component's indices cover the
		 *            interval { 0 , 1 , 2 , ... , n-1 }. This methods
		 *            updates the elements in \link labelImage\endlink
		 *            such that the pixels with value <i>i</i> belong
		 *            to the <i>i</i>th largest component
		 *            (<i>i=0,...,n-1</i>).
		 */
		void sortLabels();
	};
}


#endif /* THIMBLE_LABELER_H_ */
