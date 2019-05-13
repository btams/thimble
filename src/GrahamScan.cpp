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
 * @file GrahamScan.cpp
 *
 * @brief
 *            Implementation of a mechanism for computing the convex hull
 *            containing a set of two-dimensional coordinates as provided
 *            by the 'GrahamScan.h' header file.
 *
 * @author Benjamin Tams
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#include <thimble/math/GrahamScan.h>

using namespace std;

/**
 * @brief
 *            The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Swaps the content of two two-dimensional coordinates.
	 *
	 * @param p1
	 *            The first coordinate with which the content of
	 *            <code>p2</code> is exchanged.
	 *
	 * @param p2
	 *            The second coordinate with which the content of
	 *            <code>p1</code> is exchanged.
	 */
	static inline void swap
	( pair<double,double> & p1 , pair<double,double> & p2 ) {

		double first , second;

		first  = p1.first;
		second = p1.second;

		p1.first  = p2.first;
		p1.second = p2.second;

		p2.first  = first;
		p2.second = second;
	}

	/**
	 * @brief
	 *            Computes the cross-product of two two-dimensional
	 *            coordinates.
	 *
	 * @param p1
	 *            First coordinate.
	 *
	 * @param p2
	 *            Second coordinate.
	 *
	 * @return
	 *            The cross-product of <code>p1</code> and <code>p2</code>.
	 */
	static inline double cross
	( const pair<double,double> & p1 ,
	  const pair<double,double> & p2 ) {

		return p1.first * p2.second - p1.second * p2.first;
	}

	/**
	 * @brief
	 *            Computes the Manhattan norm of a two-dimensional
	 *            coordinate.
	 *
	 * @param p
	 *            The two-dimensional coordiante of which the Manhattan
	 *            norm is computed.
	 *
	 * @return
	 *            The Manhattan norm of <code>p</code>.
	 */
	static inline double mnorm( const pair<double,double> & p ) {

		return abs(p.first) + abs(p.second);
	}

	/**
	 * @brief
	 *            Computes the Manhattan distance between two two-dimensional
	 *            coordinates.
	 *
	 * @param p1
	 *            The first two-dimensional coordinate.
	 *
	 * @param p2
	 *            The second two-dimensional coordinate.
	 *
	 * @return
	 *            The Manhattan distance between <code>p1</code> and
	 *            <code>p2</code>.
	 */
	static double mdist
	( const pair<double,double> & p1, const pair<double,double> & p2) {

		return abs(p1.first - p2.first) + abs(p1.second - p2.second);
	}

	/**
	 * @brief
	 *            Tests whether <code>p1</code> is between the line
	 *            segment <code>p2</code> and <code>p3</code>
	 *            (in the Manhattan norm).
	 *
	 * @param p1
	 *            First coordinate.
	 *
	 * @param p2
	 *            Second coordinate.
	 *
	 * @param p3
	 *            Third coordinate.
	 *
	 * @return
	 *            <code>true</code> if <code>p1</code> is between
	 *            <code>p2</code> and <code>p3</code> and
	 *            <code>false</code> otherwise.
	 */
	static inline bool isBetween
		( const pair<double,double> & p1 ,
		  const pair<double,double> & p2 ,
		  const pair<double,double> & p3 ) {

		return mdist(p2, p3) >= mdist(p1, p2) + mdist(p1, p3);
	}

	/**
	 * @brief
	 *            Determines an element from <code>list</code> that
	 *            is an element of its convex hull.
	 *
	 * @details
	 *            More precisely, the function determines the index of
	 *            the element with the smallest second component with the
	 *            smallest first element.
	 *
	 * @param list
	 *            The list of the coordinate of which a pivot is determined.
	 *
	 * @return
	 *            The index of the pivot element in <code>list</code>.
	 */
	static int pivot( const vector< pair<double,double> > & list );

	/**
	 * @brief
	 *            Adds a constant shift to the elements in a list of
	 *            two-dimensional coordinates.
	 *
	 * @param list
	 *            The elements to which the two-dimensional vector
	 *            <code>q</code> is added.
	 *
	 * @param q
	 *            The constant shift added inplace to the elements in
	 *            <code>list</code>.
	 */
	static void move
		( vector< pair<double,double> > & list ,
		  const pair<double,double> & q );


	/**
	 * @brief
	 *            A class of which object serves as a comparator of
	 *            two two-dimensional coordinates to determine the
	 *            coordinate that the smaller angle.
	 */
	class GrahamScanComparator {
	public:

		/**
		 * @brief
		 *            Tests whether <code>p1</code> forms a smaller
		 *            angle with the abscissa axis than <code>p2</code>
		 *
		 * @param p1
		 *            First coordinate.
		 *
		 * @param p2
		 *            Second coordinate.
		 *
		 * @return
		 *            <code>true</code> if <code>p1</code> forms a
		 *            smaller angle with the abscissa axis than
		 *            <code>p2</code>; otherwise, the function returns
		 *            <code>false</code>.
		 */
		inline bool operator()
		(const pair<double,double> & p1,const pair<double,double> & p2 ) {
			double f = cross(p1,p2);
			return f > 0.0 || (f==0.0&&mnorm(p1)>mnorm(p2));
		}
	};

	/**
	 * @brief
	 *            Tests whether the triangle formed by three two-dimensional
	 *            coordinates is convex.
	 *
	 * @param p1
	 *            First edge of the triangle.
	 *
	 * @param p2
	 *            Second edge of the triangle.
	 *
	 * @param p3
	 *            Third edge of the triangle.
	 *
	 * @return
	 *            <code>true</code> if the triangle is convex; otherwise, if
	 *            the triangle is not convex, the function returns
	 *            <code>false</code>.
	 */
	static bool isConvex
	( const pair<double,double> & p1, const pair<double,double> & p2,
	  const pair<double,double> & p3 );

	/**
	 * @brief
	 *            Computes the convex hull containing a specified list
	 *            of two-dimensional coordinates.
	 *
	 * @param hull
	 *            Output vector that will contain successive
	 *            two-dimensional coordinates defining a convex polygon
	 *            containing the coordinates in <code>list</code>.
	 *
	 * @param list
	 *            Contains the two-dimensional coordinates of which the
	 *            convex hull is computed with this method.
	 */
	void GrahamScan::convexHull
		( vector< pair<double,double> > & hull ,
		  const vector< pair<double,double> > & list ) {

		int n = list.size();

		// First, make a hard copy of the points in 'list'.
		hull = list;

		// From now on, if possible, we will compute the convex hull in-place

		// A triangle is always convex; thus, if fewer than three points are
		// given we just return the hard copy of 'list'.
		if (n <= 3) {
			return;
		}

		// Determine from where to start the Graham scan
		int pivotIndex = pivot(hull);

		// Ensure that the initial point is the first point in the list.
		// Therefore, we exchange the first point in 'hull' with the initial
		// point
		swap(hull[0],hull[pivotIndex]);

		// Make a hard copy of the initial point
		pair<double,double> q = hull[0];

		// Negate 'q'. Note, if 'q' were no hard copy of the
		// initial point, the negation would also affect the initial point in
		// 'hull'.
		q.first  = -q.first;
		q.second = -q.second;

		// Translate all elements of 'hull' by 'q'. Afterwards, 'hull' will
		// held the elements of the former 'hull' with respect to the initial
		// point which was returned by 'pivot(List<Point2D.Double>)' above.
		// In particular, the first element will have 'x'- and 'y'-coordinate
		// equals to '0'.
		move(hull, q);

		// Sort the array w.r.t. their angles and distance as following
		// GrahamScanComparator'
		sort(hull.begin(),hull.end(),GrahamScanComparator());

		// Add the initial point back to the points in 'hull'. Therefore,
		// 'q' is negated which will be equals the initial point.
		// Subsequently, the elements of 'hull' are moved by using the
		// initial point 'q' as the translation vector. Note, the first
		// element of 'hull' will be shifted from the origin back to the
		// initial point
		q.first  = -q.first;
		q.second = -q.second;
		move(hull, q);

		// *******************************************************************
		// ************** BEGIN: the actual Graham scan **********************
		// *******************************************************************

		int i = 3, k = 3;

		// In each iteration we make sure that the elements with index
		// '0' to 'i' will define a convex polygon.
		while ( k < n ) {

			// 'k' is the index of the element that we currently consider.
			// Therefore we exchange it be the element 'i'.
			swap(hull[i],hull[k]);

			// Makes sure that the edge currently considered is convex
			while (!isConvex(hull[i-1], hull[i-2],hull[i])) {

				swap(hull[i-1], hull[i]);
				if ( i > 2 ) {
					--i;
				}
			}

			// Now, the 'k'-the element is processed and we continue with
			// the next, if any.

			++k;
			++i;
		}

		// *******************************************************************
		// ************** END: the actual Graham scan is finished ************
		// *******************************************************************

		// The elements in 'hull' with index '0' to 'i' will be the successive
		// points defining the polygon whose area is the convex hull of the
		// points in 'list'.

		// Remove all points from 'hull' of index larger than 'i'
		hull.erase(hull.begin()+i,hull.end());
	}

	/**
	 * @brief
	 *            Determines an element from <code>list</code> that
	 *            is an element of its convex hull.
	 *
	 * @details
	 *            More precisely, the function determines the index of
	 *            the element with the smallest second component with the
	 *            smallest first element.
	 *
	 * @param list
	 *            The list of the coordinate of which a pivot is determined.
	 *
	 * @return
	 *            The index of the pivot element in <code>list</code>.
	 */
	static int pivot( const vector< pair<double,double> > & list ) {

		int minIdx = -1;
		int size = list.size();

		if ( size <= 0) {
			// if 'list' contains no element the returned index is '-1'
			return minIdx;
		}

		double minX , minY;

		// initialize search with the first element in the list
		minX = list[0].first;
		minY = list[0].second;
		minIdx = 0;

		// iterate over all remaining elements in the list
		for ( int i = 1; i < size; i++) {

			if ( list[i].second < minY
					|| (list[i].second == minY && list[i].first < minX)) {

				// If current element has smaller 'y'-coordinate than the
				// elements before or if it has minimal 'y'-coordinate but
				// with smaller 'x'-coordinate then the minimal index is
				// set to the current index. Furthermore, the 'y'- and
				// 'x'-coordinates of the current element is saved.
				minY = list[i].second;
				minX = list[i].first;
				minIdx = i;
			}
		}

		// After all elements in the list have been considered
		// the element with smallest 'y'-coordinate (or given multiple
		// elements with smallest 'y' exist then among them the element
		// with smallest 'x'-coordinate) has been examined and its index
		// is held by 'minIdx'

		// return the index of the pivot element
		return minIdx;
	}

	/**
	 * @brief
	 *            Adds a constant shift to the elements in a list of
	 *            two-dimensional coordinates.
	 *
	 * @param list
	 *            The elements to which the two-dimensional vector
	 *            <code>q</code> is added.
	 *
	 * @param q
	 *            The constant shift added inplace to the elements in
	 *            <code>list</code>.
	 */
	static void move
		( vector< pair<double,double> > & list ,
		  const pair<double,double> & q ) {

		//Successively, iterate over all elements in the list
		for (int i = 0; i < (int)list.size() ; i++) {

			// translate each entry
			list[i].first  += q.first;
			list[i].second += q.second;
		}
	}

	/**
	 * @brief
	 *            Tests whether the triangle formed by three two-dimensional
	 *            coordinates is convex.
	 *
	 * @param p1
	 *            First edge of the triangle.
	 *
	 * @param p2
	 *            Second edge of the triangle.
	 *
	 * @param p3
	 *            Third edge of the triangle.
	 *
	 * @return
	 *            <code>true</code> if the triangle is convex; otherwise, if
	 *            the triangle is not convex, the function returns
	 *            <code>false</code>.
	 */
	static bool isConvex
	( const pair<double,double> & p1, const pair<double,double> & p2,
	  const pair<double,double> & p3 ) {


		// compute the vectors 'p2-p1' and 'p3-p1'
		pair<double,double> q1, q2;
		q1.first  = p2.first  - p1.first;
		q1.second = p2.second - p1.second;
		q2.first  = p3.first  - p1.first;
		q2.second = p3.second - p1.second;

		// compute the cross product
		// (according to 'cross(Point2D.Double,Point2D.Double)'
		double f = cross(q1,q2);

		// the edge '(p1,p2,p3)' is convex if the cross product
		// (according to 'cross(Point2D.Double,Point2D.Double)' of
		// 'p2-p1' and 'p3-p1' is '<0'; if '>0' the edge is not convex;
		// if the cross product is identical to '0' then the edge is
		// convex if 'p1' lies between the line segment with respective
		// edge points 'p2' and 'p3'
		return f < 0 || (f == 0 && !isBetween(p1, p2, p3));
	}

}




