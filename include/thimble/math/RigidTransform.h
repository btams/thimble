/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013 Benjamin Tams
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
 * @file RigidTransform.h
 *
 * @brief
 *            Provides a functionality for computing the spatial movement that
 *            minimizes the sum of least square distances between two point
 *            correspondences.
 *
 * @author Benjamin Tams
 *
 * @see thimble::RigidTransform
 */

#ifndef THIMBLE_RIGIDTRANSFORM_H_
#define THIMBLE_RIGIDTRANSFORM_H_

#include <thimble/dllcompat.h>
#include <thimble/math/AffineTransform.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Provides a static class function that can compute a spatial
	 *            movement minimizing the sum of squared distances between two
	 *            two-dimensional point correspondences.
	 */
	class THIMBLE_DLL RigidTransform {

	public:

		/**
		 * @brief
		 *           Determines a spatial movmement that minimizes
		 *           the sum of squared distances between two given
		 *           two-dimensional points correspondences.
		 *
		 * @details
		 *           The spatial movement (encoded as an affine transform)
		 *           \f$f\f$ is returned such that
		 *           \f[
		 *            \sum_{i=0}^{n-1}\|(f(mx[i],my[i])-(dx[i],dy[i]))\|_2
		 *           \f]
		 *           becomes minimal.
		 *
		 * @param dx
		 *           The abscissa coordinates of the reference points.
		 *
		 * @param dy
		 *           The ordinate coordinates of the reference points.
		 *
		 * @param mx
		 *           The abscissa coordinates of the points that are moved
		 *           to the reference points.
		 *
		 * @param my
		 *           The ordinate coordinates of the points that are moved
		 *           to the reference points.
		 *
		 * @param n
		 *           The number of input point correspondences.
		 *
		 * @return
		 *           The spatial movement of the points \f$(mx[i],my[i])\f$
		 *           that minimizes the sum of least square distance between
		 *           \f$(dx[i],dy[i])\f$.
		 *
		 * @attention
		 *           If <i>n<=1</i> the result will be the identity transform.
		 *
		 * @warning
		 *           If not sufficient memory could be allocated, the function
		 *           prints an error message to <code>stderr</code> and exits
		 *           with status 'EXIT_FAILURE'.
		 *
		 * @warning
		 *           If <code>dx</code>, <code>dy</code>, <code>mx</code>, or
		 *           <code>my</code> do not contain at least <i>n</i> valid
		 *           <code>double</code> values, the behavior of the function
		 *           is undocumented.
		 */
		static AffineTransform align
		( const double *dx , const double *dy ,
		  const double *mx , const double *my ,
		  int n );
	};
}


#endif /* THIMBLE_RIGIDTRANSFORM_H_ */
