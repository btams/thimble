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
 * @file config_unix.h
 *
 * @brief
 *            Provides defs for the library source file to
 *            successfully compile on UNIX platforms.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_CONFIG_UNIX_H_
#define THIMBLE_CONFIG_UNIX_H_

/*
 * If this is defined, this will cause the function
 * 'thimble::MathTools::rand8(true)' to generate
 *  its result from '/dev/urandom'. Because the function
 *  is implicitly invoked when
 *  'thimble::MathTools::rand16(true)',
 *  'thimble::MathTools::rand32(true)', or
 *  'thimble::MathTools::rand64(true)' are called, this
 *  affects the entire random generation related with the
 *  'thimble::MathTools' class.
 */
#define THIMBLE_USE_DEV_URANDOM

/*
 * In Linux (GCC-4.6.3) there is a round function but this is not defined in
 * the C++ standard. In fact, for example in Windows (Visual C++ 2010
 * Express), there is no such round function and we have to wrap around
 * 'std::ceil' and 'std::floor'. Therefore, we define a makro for rounding
 * here.
 */
#define THIMBLE_ROUND(v) ::round((v))

/*
 * If 'THIMBLE_GCC_X86_PCLMULQDQ' is defined, the compiler
 * implements the 'thimble::MathTools::clmul(uint32_t,uint32_t)
 * using the 'PCLMULQDQ' from the "AES instruction set" which is
 * faster than implementing the function using generic methods.
 * A consequence from this is, that multiplication of
 * 'thimble::SmallBinaryPolynomial' becomes faster.
 *
 * However, we should only define this if we are absolutely
 * certain that the compiled library will be used on a machine
 * with "AES instruction set". Therefore, we leave it undefined.
 */
//#define THIMBLE_GCC_X86_PCLMULQDQ

/*
 * We use a makro controlling the interface to open a file to
 * avoid warnings when Microsoft Visual C++ Express 2010 is
 * used to compile the library while keeping the library compatible
 * with compilers compatible with current C++ standards.
 */
#define THIMBLE_FOPEN(str,mode) \
	::fopen((str),(mode))

#define THIMBLE_FSCANF fscanf

#define THIMBLE_NAN NAN

#endif /* THIMBLE_CONFIG_UNIX_H_ */
