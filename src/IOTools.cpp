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
 * @file IOTools.cpp
 *
 * @brief
 *            Implements functionalites of the class \link thimble::IOTools
 *            IOTools\endlink provided some static convenience methods
 *            related to reading and writing data as provided by the
 *            'IOTools.h' header file.
 *
 * @author Benjamin Tams
 */

#include <cstdio>
#include <cstdlib>

#include "config.h"
#include <thimble/misc/IOTools.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

	/**
	 * @brief
	 *            Portable version of the standard %fopen() function
	 *            for opening a file for reading/writing.
	 *
	 * @details
	 *            In UNIX systems files are typically opened via
	 *            the standard %fopen() command. In Windows, however,
	 *            the use of %fopen() is deprecated and one should use
	 *            the %fsopen() command instead. This function wraps
	 *            around these two functions depending the operating
	 *            system on which THIMBLE is used.
	 *
	 * @param path
	 *            C-String specifying the path to the file to be
	 *            opened.
	 *
	 * @param mode
	 *            C-String specifying the mode with which the file
	 *            is attempted to be opened.
	 *
	 * @return
	 *            The pointer to a FILE successfully opened or NULL
	 *            if some failure occured.
	 */
	FILE *IOTools::fopen( const char *path , const char *mode ) {

		return THIMBLE_FOPEN(path,mode);
	}
}

