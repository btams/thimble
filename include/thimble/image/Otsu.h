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
 * @file Otsu.h
 *
 * @brief
 *            Provides an implementation of %Otsu's thresholding method
 *            for image segmentation
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_OTSU_H_
#define THIMBLE_OTSU_H_

#include <stdint.h>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Class that provides a static function that implements
	 *            %Otsu's method from a 256-sized histogram.
	 */
	class THIMBLE_DLL Otsu {

	public:

		/**
		 * @brief
		 *            Performs %Otsu thresholding to separate two classes
		 *            from an 8-bit source with specified histogram.
		 *
		 * @details
		 *            For example, the 8-bit source can be a 8-bit
		 *            gray-scale intensity image in which the number <i>c</i>
		 *            of pixels with intensity <i>i</i> are counted; this
		 *            number is then contained in the histogram such that
		 *            <code>histogram[i]=c</code>.
		 *
		 * @param histogram
		 *            Array containing 256 positive integers (>=0) that
		 *            describes a histogram of an 8-bit source.
		 *
		 * @return
		 *            A threshold that minmizes the intra-class varaiance
		 *            between two classes from an 8-bit source with
		 *            the specified histogram.
		 */
		static int threshold256( const int *histogram );
	};
}


#endif /* THIMBLE_OTSU_H_ */
