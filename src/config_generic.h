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
 * @file config_generic.h
 *
 * @brief
 *            Provides defs for the library source file to
 *            successfully on platforms that are compatible to
 *            a current C++ standard.
 *
 * @author Benjamin Tams
 */
#ifndef THIMBLE_CONFIG_GENERIC_H_
#define THIMBLE_CONFIG_GENERIC_H_

/*
 * We use a makro controlling the interface to open a file to
 * avoid warnings when Microsoft Visual C++ Express 2010 is
 * used to compile the library while keeping the library compatible
 * with compilers compatible with current C++ standards.
 */
#define THIMBLE_FOPEN(str,mode) \
	::fopen((str),(mode))

/*
 * In Linux (GCC-4.6.3) there is a round function but this is not defined in
 * the C++ standard. In fact, for example in Windows (Visual C++ 2010
 * Express), there is no such round function and we have to wrap around
 * 'std::ceil' and 'std::floor'. Therefore, we define a makro for rounding
 * here.
 */
#define THIMBLE_ROUND(v) \
	((std::ceil((v))-(v))<((v)-std::floor((v)))?std::ceil((v)):std::floor((v)))

#define THIMBLE_FSCANF fscanf

#define THIMBLE_NAN (0.0/0.0)

#endif /* THIMBLE_CONFIG_GENERIC_H_ */
