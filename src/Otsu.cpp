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
 * @file Otsu.cpp
 *
 * @brief
 *            Implementation of Otsu's thresholding method
 *            for image segmentation.
 *
 * @author Benjamin Tams
 */
#include <iostream>
#include <cfloat>

#include <thimble/image/Otsu.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

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
	int Otsu::threshold256( const int *histogram ) {

		 // Number of intensities
		static const int NUM_INTENSITIES = 256;

		// Number of intensities in the source (e.g., image), i.e.,
		// sum_{j=0}^{255} histogram[j]
		double sum;

		// sum_{j=0}^{255} t * histogram[j]
		double sum_t;

		// sum_{j=0}^t histogram[j] for each t = 0 ,... , 255
		double sum_histogram;

		// sum_{j=0}^t j * histogram[j] for each t = 0 , ... , 255
		double sum_t_histogram;

		// Keeps track of the maximal objective quotient between the between
		// variance and the inner variance
		double Qmax = -DBL_MAX;

		// Keeps track of the optimal threshold that maximizes the
		// objective quotient
		int thres = -1;

		// Determine the sums as above
		sum = 0.0;
		sum_t = 0.0;
		for ( int j = 0 ; j < NUM_INTENSITIES ; j++ ) {
			sum   += (double)histogram[j];
			sum_t += (double)j * (double)histogram[j];
		}

		sum_histogram = 0.0;
		sum_t_histogram = 0.0;

		// Iterate through all pixel intensities as a candidate for the
		// optimal threshold, ...
		for ( int t = 0 ; t < NUM_INTENSITIES ; t++ ) {

			sum_histogram += (double)histogram[t];
			sum_t_histogram += (double)t * (double)histogram[t];

			// inner variance
			double varInner =
					(double)sum_histogram * (double)(sum - sum_histogram);

			// ... compute the objective quotient and ...
			double Q;
			if ( varInner == 0 ) {
				Q = 0.0;
			} else {
				// Between variance
				double varBetween =
					((double) sum_histogram / (double)sum) * sum_t -
						sum_t_histogram;
				Q = (varBetween * varBetween) / varInner;
			}

			// ... and if it is larger then those before, ...
			if ( Q > Qmax ) {
				// ..., update.
				Qmax = Q;
				thres = t;
			}

		}

		// Return the optimal threshold.
		return thres;
	}
}



