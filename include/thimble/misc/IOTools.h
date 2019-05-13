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
 * @file IOTools.h
 *
 * @brief
 *            Provides a class with some static convenience methods related
 *            to reading and writing data.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_IOTOOLS_H_
#define THIMBLE_IOTOOLS_H_

#include <cstdio>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {


	/**
	 * @brief
	 *           Provides some static convenience methods related to reading
	 *           and writing data.
	 *
	 * @details
	 *           This class will be constituted with new functionalities
	 *           only if needed during the development of THIMBLE.
	 */
	class THIMBLE_DLL IOTools {

	private:

		/**
		 * @brief
		 *            Private standard constructor to prevent that instances
		 *            of this class are created.
		 */
		inline IOTools();

	public:

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
		static FILE *fopen( const char *path , const char *mode );
	};
}


#endif /* THIMBLE_IOTOOLS_H_ */
