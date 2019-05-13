/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013, 2014 Benjamin Tams
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
 * @file config.h
 *
 * @brief
 *            Provides the defs for the library source file to successfully
 *            compile on this platform.
 *
 * @details
 *            This header examines the type of this platform and to
 *            include a corresponding header file on compilation. If the type
 *            can not be examined, a generic header file is included that
 *            should work for compilers that are compatible with the current
 *            C++ standard.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_IMPL_CONFIG_H_
#define THIMBLE_IMPL_CONFIG_H_

#if defined(__unix__) || defined(__unix) || \
	(defined(__APPLE__) && defined(__MACH__))
#include "config_unix.h"
#else

#ifdef _MSC_VER
#include "config_windows_msc.h"
#else
#include "config_generic.h"
#endif

#endif

#endif /* THIMBLE_IMPL_CONFIG_H_ */
