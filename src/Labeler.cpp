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
 * @file Labeler.cpp
 *
 * @brief
 *            Implementation of the class provided in the 'Labeler.h' header
 *            for dividing two-dimensional binary images into their connected
 *            components.
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <climits>
#include <vector>
#include <algorithm>
#include <iostream>

#include <thimble/image/Labeler.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Creates a labeler object for dividing a two-dimensional
	 *            binary image into their connected components.
	 *
	 * @param connectivity
	 *            The connectivity type of the binary images.
	 */
	Labeler::Labeler( LABEL_CONNECTIVITY_T connectivity ) {

		this->connectivity = connectivity;

		this->m            = 0;
		this->n            = 0;
		this->labelImage   = NULL;
		this->labelCount   = 0;
	}

	/**
	 * @brief
	 *            Destructor.
	 */
	Labeler::~Labeler() {
		free(this->labelImage);
	}

	/**
	 * @brief
	 *            Accesses the connectivity of this labeler.
	 *
	 * @return
	 *            The connectivity code with which this object has been
	 *            created.
	 */
	LABEL_CONNECTIVITY_T Labeler::getConnectivity() const {

		return this->connectivity;
	}

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
	void Labeler::label( const bool *binaryImage , int m , int n ) {

		free(this->labelImage);

		// initialize an empty label image whose entries are all -1
		this->m = m;
		this->n = n;
		this->labelImage = (int*)malloc( m * n * sizeof(int) );
		if ( this->labelImage == NULL ) {
			cerr << "Labeler::label: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				this->labelImage[i*n+j] = -1;
			}
		}

		//Pass through the image the first time
		bool hasComponents = firstPass(binaryImage);

		this->labelCount = 0;

		if ( hasComponents ) {

			// If the first pass determined non-empty components we
			// alternately pass through the image forward and backward
			// as long as there are changes in the labels.
			bool forward, backward;
			do {
				forward = passForward();
				backward = passBackward();
			} while (forward || backward);
			// Since the last passes through the image did not causes changes
			// of the labels all connected components must have been examined

			// Next, we count the number of different labels and make their
			// make sure that the label indices range from 0 to the 'number
			// of different labels'-1

			// Variable that keeps track of the successive label index
			// of connected components founds
			int k = 0;
			for ( int l = -1 ; ; k++ ) {

				// Determine the next label index that is greater than 'l'.
				l = nextLabel(l);
				if (l == INT_MAX) {
					// If the label index is 'Integer.MAX_VALUE' then there
					// are no more labels.
					break;
				}

				if (k != l) {
					// If the label index is greater than 'k' then replace
					// all entries equals the label index by 'k'.
					for (int i = 0; i < m; i++) {
						for (int j = 0; j < n; j++) {
							if ( this->labelImage[i*n+j] == l) {
								 this->labelImage[i*n+j] = k;
							}
						}
					}
				}
			}

			// The number of connected components is equals to the latest
			// updated label index plus '1'.
			this->labelCount = k + 1;
		}

		// Sort the labels by their corresponding component size in
		// descending order
		sortLabels();
	}

	/**
	 * @brief
	 *            Returns the height of the latest binary image for
	 *            which this labeler has been applied.
	 *
	 * @return
	 *            The height of the latest binary image for which
	 *            this labeler has been applied.
	 */
	int Labeler::getHeight() const {

		return this->m;
	}

	/**
	 * @brief
	 *            Returns the width of the latest binary image for
	 *            which this labeler has been applied.
	 *
	 * @return
	 *            The width of the latest binary image for which
	 *            this labeler has been applied.
	 */
	int Labeler::getWidth() const {

		return this->n;
	}

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
	int Labeler::getLabelCount() const {

		return this->labelCount;
	}

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
	void Labeler::get( bool *binaryImage , int labelIndex ) const {

		for ( int i = 0 ; i < m ; i++ ) {
			for ( int j = 0 ; j < n ; j++ ) {
				binaryImage[i*n+j] = (labelImage[i*n+j]==labelIndex);
			}
		}
	}

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
	int Labeler::labelSize( int labelIndex ) const {

		// Counts the number of positions '(i,j)' in where
		// 'labelImage[i][j]==labelIndex'
		int count = 0;

		// Iterate over all pixels
		for (int i = 0; i < m ; i++) {
			for (int j = 0; j < n ; j++) {

				if ( labelImage[i*n+j] == labelIndex) {
					// If the current pixel's value in 'labelImage'
					// is equals to 'labelIndex' increment the count
					++count;
				}
			}
		}

		// Return the count
		return count;
	}

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
	bool Labeler::firstPass( const bool *binaryImage ) {

		int m , n;
		m = this->m;
		n = this->n;

		// Variable that holds the number of connected pixel areas for
		// that a common connected area could be determined. Note,
		// 'count' is not the final number of connected components.
		// More precisely, 'count' will be greater or equals to the
		// final number of connected components.
		int count = 0;

		// Variable which indicates whether there have been found at least
		// one non-empty connected component
		bool changed = false;

		// Iterate over each pixel
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {

				if ( binaryImage[i*n+j] ) {

					// If a pixel's value is 'true' there must be a connected
					// component to where the pixel '(i,j)' belongs to.
					// Parts of this component are firstly labeled now.
					count = newLabel(count,binaryImage, i, j);
					// Now, count holds the number of connected pixel areas
					// for that a connected area could be determined by now.

					// The labels have changed which is the case if and only
					// if a pixel with value 'true' exists.
					changed = true;
				}
			}
		}

		// Return whether there are non-empty components or not.
		return changed;
	}

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
	bool Labeler::passForward() {

		// Stores whether at least one pixel's label changed
		bool changed = false;

		// Iterate over all rows from top to bottom
		for (int i = 0; i < m; i++) {

			// Iterate over all columns from left to right
			for (int j = 0; j < n; j++) {

				if (labelImage[i*n+j] >= 0) {
					// If the pixel '(i,j)' belongs to a component, i.e. if
					// 'labelImage[i][j]>=0', check whether to relabel.
					changed |= relabel(i, j);
				}
			}
		}

		// Is 'true' if at least one pixel was relabeled
		return changed;
	}

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
	bool Labeler::passBackward() {

		// Stores whether at least one pixel's label changed
		bool changed = false;

		// Iterate over all rows from bottom to top
		for (int i = m - 1; i >= 0; i--) {

			// Iterate over all columns from right to left
			for (int j = n - 1; j >= 0; j--) {
				if (labelImage[i*n+j] >= 0) {

					// If the pixel '(i,j)' belongs to a component, i.e. if
					// 'labelImage[i][j]>=0', check whether to relabel.
					changed |= relabel(i, j);
				}
			}
		}

		// Is 'true' if at least one pixel was relabeled
		return changed;
	}

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
	int Labeler::newLabel
		( int count , const bool *binaryImage , int i , int j ) const {

		int maxRowIndex , maxColIndex;
		// maximal index that accesses a valid row
		maxRowIndex = this->m-1;
		// maximal index that accesses a valid column
		maxColIndex = this->n-1;

		// Determine those neighbors relative to the pixel '(i,j)' that
		// also belong to a component. Note that if the pixel '(i,j)'
		// encompasses the image's border then those neighbors that would
		// exceed the image's bound are set to 'false'.
		// Note, the eight-connected neighbors which are
		// 'northeast', 'southeast', 'southwest', and 'northwest'
		// can only be 'true' if 'this.connectivity' indicates that the
		// labeled image is eight-connected
		bool north, east, south, west,
		northeast, southeast, southwest, northwest;
		north = (i > 0 ? binaryImage[(i-1)*n+j] : false);
		east = (j < maxColIndex ? binaryImage[i*n+j+1] : false);
		south = (i < maxRowIndex ? binaryImage[(i+1)*n+j] : false);
		west = (j > 0 ? binaryImage[i*n+j-1] : false);
		northeast = (this->connectivity == LC_EIGHT_CONNECTED && i > 0
				&& j < maxColIndex ? binaryImage[(i-1)*n+j+1] : false);
		southeast = (this->connectivity == LC_EIGHT_CONNECTED
				&& i < maxRowIndex &&
				j < maxColIndex ? binaryImage[(i+1)*n+j + 1] : false);
		southwest = (this->connectivity == LC_EIGHT_CONNECTED
				&& i < maxRowIndex && j > 0 ? binaryImage[(i+1)*n+j-1] : false);
		northwest = (this->connectivity == LC_EIGHT_CONNECTED && i > 0
				&& j > 0 ? binaryImage[(i-1)*n+j-1] : false);

		int labelIndex , neighborLabelIndex;

		// label index of pixel '(i,j)'
		labelIndex = count;

		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its north neighbor. Labeling by the north neighbor is only
		// conceived if the neighbor corresponds to a smaller label index than
		// the pixel '(i,j)' does
		if (north &&
			(neighborLabelIndex = labelImage[(i-1)*n+j]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its east neighbor.
		if (east &&
			(neighborLabelIndex = labelImage[i*n+j+1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its south neighbor.
		if (south &&
			(neighborLabelIndex = labelImage[(i+1)*n+j]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its west neighbor.
		if (west &&
			(neighborLabelIndex = labelImage[i*n+j-1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its northeast neighbor.
		if (northeast &&
			(neighborLabelIndex = labelImage[(i-1)*n+j+1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its southeast neighbor.
		if (southeast &&
			(neighborLabelIndex = labelImage[(i+1)*n+j+1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its southwest neighbor.
		if (southwest &&
			(neighborLabelIndex = labelImage[(i+1)*n+j-1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}
		// Checks whether the pixel '(i,j)' should be labeled by the label
		// of its northwest neighbor.
		if (northwest &&
			(neighborLabelIndex = labelImage[(i-1)*n+j-1]) >= 0 &&
			neighborLabelIndex < labelIndex) {

			labelIndex = neighborLabelIndex;
		}

		// Label the pixel '(i,j)' by the minimal label index of neighbor
		// corresponding to an examined connected area. If no neighbor
		// corresponds to an examined connected area then 'labelIndex' is
		// equals to 'count'
		labelImage[i*n+j] = labelIndex;

		if ( labelIndex == count ) {
			// If no neighbors correspond to an examined connected area
			// then we included a further label index which requires to
			// increase the count
			return count + 1;
		}

		// return the number of examined connected areas
		return count;
	}

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
	bool Labeler::relabel( int i , int j ) {

		int maxRowIndex , maxColIndex;
		// maximal index that accesses a valid row
		maxRowIndex = m-1;
		// maximal index that accesses a valid column
		maxColIndex = n-1;

		// Determine those neighbors relative to the pixel '(i,j)' that
		// also belong to a component. Note that if the pixel '(i,j)'
		// encompasses the image's border then those neighbors that would
		// exceed the image's bound are set to 'false'.
		// Note, the eight-connected neighbors which are
		// 'northeast', 'southeast', 'southwest', and 'northwest'
		// can only be 'true' if
		// 'this.connectivity' indicates that the labeled image was
		// eight-connected
		bool north, east, south, west, northeast,
	    southeast, southwest, northwest;
		north = (i > 0 ? labelImage[(i-1)*n+j] >= 0 : false);
		east = (j < maxColIndex ? labelImage[i*n+j + 1] >= 0 : false);
		south = (i < maxRowIndex ? labelImage[(i+1)*n+j] >= 0 : false);
		west = (j > 0 ? labelImage[i*n+j-1] >= 0 : false);
		northeast =
				(this->connectivity == LC_EIGHT_CONNECTED && i > 0 &&
				 j < maxColIndex ? labelImage[(i-1)*n+j+1] >= 0 : false);
		southeast =
				(this->connectivity == LC_EIGHT_CONNECTED &&
				 i < maxRowIndex &&
				 j < maxColIndex ? labelImage[(i+1)*n+j+1] >= 0
				: false);
		southwest =
				(this->connectivity == LC_EIGHT_CONNECTED &&
				 i < maxRowIndex &&
				 j > 0 ? labelImage[(i+1)*n+j-1] >= 0 : false);
		northwest =
				(this->connectivity == LC_EIGHT_CONNECTED &&
				i > 0 &&
				j > 0 ? labelImage[(i-1)*n+j-1] >= 0 : false);

		int labelIndex , neighborLabelIndex;

		// label index of pixel '(i,j)'
		labelIndex = labelImage[i*n+j];

		// keeps track whether labels change
		bool changed = false;

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its north neighbor. Relabeling is only conceived if the neighbor
		// corresponds to a smaller label index than the pixel '(i,j)' does
		if (north &&
			(neighborLabelIndex = labelImage[(i-1)*n+j]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its east neighbor.
		if (east &&
			(neighborLabelIndex = labelImage[i*n+j+1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its south neighbor.
		if (south &&
			(neighborLabelIndex = labelImage[(i+1)*n+j]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its west neighbor.
		if (west &&
			(neighborLabelIndex = labelImage[i*n+j-1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its northeast neighbor.
		if (northeast &&
		    (neighborLabelIndex = labelImage[(i-1)*n+j+1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its southeast neighbor.
		if (southeast &&
			(neighborLabelIndex = labelImage[(i+1)*n+j+1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its southwest neighbor.
		if (southwest &&
			(neighborLabelIndex = labelImage[(i+1)*n+j-1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// Checks whether the pixel '(i,j)' should be relabeled by the label
		// of its northwest neighbor.
		if (northwest &&
			(neighborLabelIndex = labelImage[(i-1)*n+j-1]) < labelIndex) {
			labelIndex = neighborLabelIndex;
			changed = true;
		}

		// If '(i,j)' was found to be relabeled by the label of at least one
		// of its neighbors, i.e. if 'changed' is 'true', then the label
		// of pixel '(i,j)' is set to the examined label index
		if (changed) {
			labelImage[i*n+j] = labelIndex;
		}

		// Return whether the pixel '(i,j)' was relabeled
		return changed;
	}

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
	int Labeler::nextLabel( int labelIndex ) const {

		// Stores the label index found which is the smallest among those
		// greater than 'labelIndex'.
		int nextLabelIndex = INT_MAX;

		// Iterate over each pixel's coordinate
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {

				// label index of the current pixel
				int currentLabelIndex = labelImage[i*n+j];

				if ( currentLabelIndex < nextLabelIndex &&
					 currentLabelIndex > labelIndex) {

					// If the current label index is strictly greater than
					// 'labelIndex' and smaller than the smallest candidate
					// update the next label index
					nextLabelIndex = currentLabelIndex;
				}
			}
		}

		// return the next label index
		return nextLabelIndex;
	}

	/**
	 * @brief
	 *            Comparator class for comparing two pairs of integers
	 *            w.r.t. to their second entry.
	 *
	 * @details
	 *            Objects of this class are used to compare two components
	 *            stored by a \link Labeler\endlink object where the
	 *            first entry of a pair holds the label index and the
	 *            second the pixel size of that component.
	 */
	class LabelComparator {

	public:

		/**
		 * @brief
		 *             Compares two components.
		 *
		 * @details
		 *             A component is described as a pair of integers where
		 *             the first entry equals its label index stored by
		 *             a Labeler object and the second equals the component's
		 *             pixel size.
		 *
		 * @param label1
		 *             First component.
		 *
		 * @param label2
		 *             Second component.
		 *
		 * @return
		 *             <code>true</code> if the pixel size of
		 *             <code>label1</code> is greater than the pixel size of
		 *             <code>label2</code>; otherwise, this function returns
		 *             <code>false</code>.
		 */
		inline bool operator()
		(const pair<int,int> & label1 , const pair<int,int> & label2 ) {

			return label1.second > label2.second;
		}
	};

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
	void Labeler::sortLabels() {

		// Create a vector of pairs holding a pair of integers,
		// where the first entries hold a component's label index
		// and the second the component's pixel size.
		vector< pair<int,int> > labels(this->labelCount);
		for ( int k = 0 ; k < this->labelCount ; k++ ) {
			int size = labelSize(k);
			labels.push_back( pair<int,int>(k,size) );
		}

		// Sort the components w.r.t. their pixel size.
		sort(labels.begin(),labels.end(),LabelComparator());

		// Replace label indices their index + the number of different
		// components. The shift by the number of different
		// components is performed temporarily to guarantee a separation
		// between disjoint components at any time.
		for ( int l = 0 ; l < labelCount ; l++ ) {

			int k = labels[l].first;

			for ( int i = 0 ; i < m ; i++ ) {
				for ( int j = 0 ; j < n ; j++ ) {
					if ( labelImage[i*n+j] == k ) {
						labelImage[i*n+j] = l+labelCount;
					}
				}
			}
		}

		// Shift back by the number of different components such that
		// the indices cover the interval { 0 , 1 , 2 , ... , n-1 }.
		for ( int i = 0 ; i < m ; i++ ) {
			for ( int j = 0 ; j < n ; j++ ) {
				if ( labelImage[i*n+j] >= labelCount ) {
					labelImage[i*n+j] -= labelCount;
				}
			}
		}

	}

}




