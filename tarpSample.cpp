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
 * @file tarpSample.cpp
 *
 * @brief
 *            Implements a small sample program that uses the THIMBLE
 *            library in which a directed reference point is estimated
 *            from a 500DPI fingerprint image that is then output to
 *            'stdout'
 *
 * @author Benjamin Tams
 */

#include <cstdlib>
#include <iostream>

#include <thimble/all.h>

// Exit status if no directed reference could be estimated
#define EXIT_FTA (EXIT_FAILURE+1)

using namespace std;
using namespace thimble;

/**
 * @brief
 *           Reads a fingerprint image specified by the argument
 *           list and attempts to estimate directed reference point
 *           from the image that is printed to 'stdout' on success.
 *
 * @details
 *           The program exits with status 'EXIT_FAILURE' if an
 *           error occurred and with  status 'EXIT_FTA' if no
 *           directed reference could be estimated. Otherwise, the
 *           main function returns 'true'
 *
 * @return
 *           'true' on success; otherwise, the program exits with
 *           a corresponding error state.
 */
int main ( int argc , char *args[] ) {

	if ( argc != 2 ) {
		cerr << "USAGE: " << args[0] << " <pgmFile>" << endl;
		cerr << "where <pgmFile> is the PGM file containing " << endl
			 << "      a 500DPI fingerprint image." << endl;
		cerr << "The program prints a directed reference point (format "
			 << "'x y radianAngle') estimated from <pgmFile> to 'stdout'."
			 << endl;
		exit(EXIT_FAILURE);
	}

	Fingerprint fingerprint;
	bool success;

	// Read fingerprint image
	success = fingerprint.fromImageFile(args[1]);
	if ( !success ) {
		cerr << "ERROR: could not read PGM file from '"
			 << args[1] << "'" << endl;
		exit(EXIT_FAILURE);
	}

	// Ensure that the fingerprint contains a valid directed reference point
	// estimation
	if ( !fingerprint.hasDirectedReferencePoint() ) {
		cerr << "FAILURE TO ALIGN" << endl;
		exit(EXIT_FTA);
	}

	// Output the directed reference point estimation
	cout.precision(21);
	cout << fingerprint.getDirectedReferencePoint() << endl;

	return 0;
}
