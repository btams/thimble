/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013, 2015 Benjamin Tams
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
 * @file config_windows_msc.h
 *
 * @brief
 *            Provides defs for the library source files to
 *            successfully compile as a WIN32 library with Microsoft's
 *            Visual C++ compilers.
 *
 * @details
 *            We successfully tested compilation using Microsoft
 *            Visual C++ Express 2010 and 2012.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_CONFIG_WINDOWS_MSC_H_
#define THIMBLE_CONFIG_WINDOWS_MSC_H_

#include <share.h>

/*
 * For Microsoft Visual C++ the cryptographic random number
 * generator is 'rand_s' but which is not available for other
 * platforms. If this is defined, the library is advised to
 * use 'rand_s' to produce the result of
 * 'thimble::MathTools::rand8(true)' and thus, implicitly, for
 *  'thimble::MathTools::rand16(true)',
 *  'thimble::MathTools::rand32(true)', and
 *  'thimble::MathTools::rand64(true)'.
 */
#define THIMBLE_USE_RAND_S

/*
 * In Linux (GCC-4.6.3) there is a round function but this is not defined in
 * the C++ standard. In fact, for example in Windows (Visual C++ 2010
 * Express), there is no such round function and we have to wrap around
 * 'std::ceil' and 'std::floor'. Therefore, we define a makro for rounding
 * here.
 */
#define THIMBLE_ROUND(v) \
	((std::ceil((v))-(v))<((v)-std::floor((v)))?std::ceil((v)):std::floor((v)))

/*
 * We use a makro controlling the interface to open a file to
 * avoid warnings when Microsoft Visual C++ Express 2010 is
 * used to compile the library while keeping the library compatible
 * with compilers compatible with current C++ standards.
 */
#define THIMBLE_FOPEN(str,mode) \
	_fsopen((str),(mode),_SH_DENYWR)

#define THIMBLE_FSCANF fscanf_s

#define THIMBLE_NAN (std::numeric_limits<double>::quiet_NaN())

#endif /* THIMBLE_CONFIG_WINDOWS_MSC_H_ */
