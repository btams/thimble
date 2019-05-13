/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2014, 2015 Benjamin Tams
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
 * @file dllcompat.h
 *
 * @brief
 *            This file provides a makro 'THIMBLE_DLL' with which classes
 *            and functions can automatically be exported for Win32 DLLs.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_DLLCOMPAT_H_
#define THIMBLE_DLLCOMPAT_H_

#ifdef _MSC_VER

// Comment out if we do not want to build a DLL library
// #define THIMBLE_BUILD_DLL

#ifdef THIMBLE_BUILD_DLL
#define THIMBLE_DLL __declspec(dllexport)
//#define THIMBLE_DLL __declspec(dllimport)
#else
#define THIMBLE_DLL
#endif

#else

#define THIMBLE_DLL

#endif /* _MSC_VER */


#endif /* THIMBLE_DLLCOMPAT_H_ */
