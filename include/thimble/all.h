/*
 *  THIMBLE --- A Library for Research, Development, and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
 *
 *  Copyright 2013, 2014, 2015 Benjamin Tams
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
 * @file include/thimble/all.h
 *
 * @brief
 *            Provides all the functionalities that are provided by
 *            the THIMBLE library.
 *
 * @author Benjamin Tams
 */

/**
 * @mainpage A Research Library for Analysis and Proof-of-Concept Development of Fingerprint-Based Biometric Cryptosystems
 *
 * @section sec_toc Contents
 * <ul>
 *  <li>
 *   @ref sec_intro
 *  </li>
 *  <li>
 *   @ref sec_interface
 *  </li>
 *  <li>
 *   @ref sec_overview
 *  </li>
 *  <li>
 *   @ref sec_compile_unix
 *  </li>
 *  <li>
 *   @ref sec_compile_windows
 *  </li>
 *  <li>
 *   @ref sec_knownproblems
 *  </li>
 *  <li>
 *   @ref sec_changelog
 *  </li>
 * </ul>
 *
 * @section sec_intro Introduction
 *
 * THIMBLE is a portable C++ library providing data structures for the
 * protection of fingerprint minutiae templates. Many sub-routines
 * (such as decoders of Reed-Solomon codes) are provided as well
 * making THIMBLE a useful tool for research, analysis, and
 * development of biometric template protection --- not just restricted to
 * fingerprint minutiae templates.
 *
 * THIMBLE was implemented with the intention of being highly portable:
 * It compiles on Linux and Windows without fitting the code.
 * Furthermore, it comes without any external dependencies --- except
 * the C++ standard library.
 *
 * I publish THIMBLE in the hope that it will be useful to the biometric
 * community; in this way, the implementations that I have already
 * implemented can reused by other researchers.
 *
 * For me, as a single person, it is pretty much work to implement a
 * C++ library that is accessible by others and to provide a proper
 * documentation in English which is not my native language. There are,
 * most likely, typos, errors, and inconsistencies contained in the
 * documentation and I appreciate that you understand if you encounter some
 * problems with THIMBLE due to errors in the documentation. So if you
 * discover misleading flaws in the documentation (or even programming bugs),
 * I would be most grateful if you report them to me and in return I will try
 * to help you out with your programming problem resulting from a wrong
 * documentation in THIMBLE.
 *
 * @warning THIMBLE has been implemented with all due caution and care. Yet,
 * neither me nor any other person can be made responsible for any damages
 * caused by using THIMBLE. By using THIMBLE you agree to this policy.
 *
 * @section sec_interface Programming Interface
 *
 * Each of THIMBLE's provided function, method, and data structure can be
 * provided by the <code>thimble/all.h</code> header. Furthermore, each of
 * THIMBLE's function, method, and data structure is wrapped into the
 * namespace <code>thimble</code>. Consequently, an
 * &quot;allround no-worry&quot; way to access THIMBLE's interface in a
 * program is to write
 * <pre>
 *  \#include <thimble/all.h>
 *
 *  using namespace std;
 *  using namespace thimble;
 * </pre>
 * before implementing functions. <b>NOTE:</b> During the documentation
 * of the THIMBLE library we assume the program code has included all of
 * THIMBLE's functionalities and that it uses the namespaces
 * <code>std</code> and <code>thimble</code> as above.
 *
 * @section sec_overview Overview of THIMBLE's Modules
 *
 * THIMBLE's main functionalities include
 * <ul>
 *  <li>
 *   @ref sec_minutiae
 *  </li>
 *  <li>
 *   @ref sec_grayimages
 *  </li>
 *  <li>
 *   @ref sec_fingerprints
 *  </li>
 *  <li>
 *   @ref sec_tarp
 *  </li>
 *  <li>
 *   @ref sec_ffv
 *  </li>
 * </ul>
 * Some lower level functionalities are provided as well and
 * may be useful when one wants to implement biometric template protection
 *.In particular, in order to implement a fuzzy vault, important ingredients
 *.are finite fields and polynomials with coefficients over these finite
 *.fields. For both, brief tutorials on how THIMBLE provides
 * them are given:
 * <ul>
 *  <li>
 *   @ref sec_finitefields
 *  </li>
 *  <li>
 *   @ref sec_polynomials
 *  </li>
 * </ul>
 * Furthermore, THIMBLE can be used to implement programs requiring
 * the following ingredients:
 * <ul>
 *  <li>
 *   @ref sec_rsdecode
 *  </li>
 *  <li>
 *   @ref sec_gsdecode for List-Decoding Reed-Solomon Codes
 *  <li>
 *   \link thimble::BinaryVector Binary Vectors\endlink
 *   and \link thimble::BinaryMatrix Binary Matrices\endlink
 *  </li>
 *  <li>
 *   \link thimble::BCHCode Binary BCH Codes\endlink (in particular useful
 *   to perform simulations with a binary fuzzy commitment scheme)
 *  </li>
 *  <li>
 *   \link thimble::AES128 Advanced Encryption Algorithm\endlink
 *  </li>
 *  <li>
 *   \link thimble::SHA Secure Hash Algorithm\endlink
 *  </li>
 *  <li>
 *   \link thimble::BigInteger Big Integers\endlink
 *  </li>
 *  <li>
 *   \link thimble::SteepestDescent Numerical Minimization\endlink
 *   of \link thimble::RealFunctional
 *   Real Functionals\endlink / \link thimble::RealFunction
 *   Real Functions\endlink
 *  </li>
 * </ul>
 * Furthermore, THIMBLE can be used to reproduce the
 * attack presented in <i><b>Merkle and Tams (2013).</b>
 * Security of the Improved Fuzzy Vault Scheme in the Presence of Record
 * Multiplicity, <a href="http://arxiv.org/abs/1312.5225" target="_blank">
 * arXiv:1312.5225</a></i>:
 * <ul>
 *  <li>
 *   @ref sec_arm_improved_fv
 *  </li>
 * </ul>
 * For more of THIMBLE's functionalities, we refer the reader of this
 * DOXYGEN documentation to the class/file list.
 *
 * @section sec_compile_unix Compilation on UNIX
 *
 * NOTE: I tested the compilation and installation as described here with
 * Ubuntu Linux 14.04, but it may also work for other UNIX like systems like
 * MinGW/MSYS etc..
 *
 * To obtain the source code of THIMBLE, download the latest version of
 * THIMBLE from
 * <a href="http://www.stochastik.math.uni-goettingen.de/biometrics/thimble">
 * here</a> to the local folder <code>$FOLDER</code> say. Then extract it in
 * the directory by doing
 * <pre>
 *  cd $FOLDER
 *  tar xvf thimble-yyyy.mm.dd.tar.gz
 * </pre>
 * where 'yyyy', 'mm', and 'dd' are delimiters for the year, month, and day,
 * respectively, the version was released.
 *
 * @subsection sec_make_compilation Compilation
 *
 * To compile THIMBLE run
 * <pre>
 *  cd $FOLDER/thimble-yyyy.mm.dd
 *  make
 * </pre>
 * If successful, this should produce the file
 * <pre>
 *  libthimble.a
 * </pre>
 * as well as the sample binary
 * <pre>
 *  tarpSample
 * </pre>
 * To test whether compilation was successful we can run
 * <pre>
 *  ./tarpSample fingerprint.pgm
 * </pre>
 * (there should exist the file 'fingerprint.pgm' along with the binary).
 * This may output the following directed reference point estimation of
 * the fingerprint image:
 * <pre>
 *  136.851968181911871625 228.042596947062747859 0.173504083499781386335
 * </pre>
 *
 *
 * @subsection sec_usingthimble Using THIMBLE in a Program
 *
 * If we want to write an executable program using THIMBLE whose source code
 * <code>foo.cpp</code> is of the form
 * <pre>
 * \#include <thimble/all.h>
 *
 * using namespace std;
 * using namespace thimble;
 *
 * int main() {
 * 	// Implementation
 * 	return 0;
 * }
 * </pre>
 * we can compile it via <code>g++</code>. Therefore, we need to specify
 * THIMBLE's include-directory and the directory where THIMBLE's static
 * library besides; furthermore, we need to link the executable with the
 * library <code>libthimble.a</code> as well as with the global standard
 * math library. These tasks can be realized via the
 * command
 * <pre>
 *  g++ -I$FOLDER/thimble/include foo.cpp -o foo -L$FOLDER/thimble-yyyy.mm.dd/ -lthimble -lm
 * </pre>
 * which should produce the executable file
 * <pre>
 *  foo
 * </pre>
 * in the current work direction.
 *
 * @subsection sec_install Installation
 *
 * If we want to store the THIMBLE library in a global place, we may
 * copy the include-directory to <code>/usr/local/include/</code> and the
 * static library to <code>/usr/local/lib/</code>. This can be realized
 * via the command
 * <pre>
 *  make install
 * </pre>
 * as root in the directory $FOLDER/thimble. WARNING: This will replace the
 * folder <code>/usr/local/include/thimble</code> (if it already exists) by
 * deleting the entire old content; furthermore, the file
 * <code>/usr/local/lib/libthimble.a</code> (if any) will be replaced.
 * Deinstallation can be performed manually via the
 * <pre>
 *  make uninstall
 * </pre>
 * command (again as root).
 *
 * After putting the THIMBLE library in a global place, we can compile
 * the files like <code>foo.cpp</code> from above more conveniently, via
 * <pre>
 *  g++ foo.cpp -o foo -lthimble -lm
 * </pre>
 * without specifying THIMBLE's include-directory and the place where
 * the static library besides. Furthermore, it is available to all users
 * of the system.
 *
 * @subsection sec_doxy Build the Documentation
 *
 * To build this documentation we may run
 * <pre>
 *  make doxy
 * </pre>
 * This will create the documentation in
 * <code>$FOLDER/thimble-yyyy.mm.dd/doc/html</code>.
 * Opening <code>$FOLDER/thimble-yyyy.mm.dd/doc/html/index.html</code> in a
 * web-browser will get us to the main-page of the documentation.
 *
 * Building the documentation requires an installation of
 * <code>doxygen</code> and <code>doxygen-latex</code> which can easily
 * be installed via <code>apt</code> for example on Ubuntu:
 * <pre>
 *  sudo apt-get install doxygen doxygen-latex
 * </pre>
 *
 * @subsection sec_cleaning Cleaning Up
 *
 * To clean all the object files, the static library, and the binaries
 * we may run
 * <pre>
 *  make clean
 * </pre>
 * in the <code>$FOLDER/thimble-yyyy.mm.dd/</code>. If we also want to delete
 * the documentation possibly generated by <code>doxygen</code>, we can run
 * <pre>
 *  make deldoxy
 * </pre>
 * All the files that have been generate by <code>make</code> in
 * <code>$FOLDER/thimble-yyyy.mm.dd/</code> can be removed via
 * <pre>
 *  make clobber
 * </pre>
 * Note, that this will not remove the installation of THIMBLE in the global
 * directories. If we also want to remove the installation, we
 * should run
 * <pre>
 * make uninstall
 * </pre>
 * as root.
 *
 * @section sec_compile_windows Compilation on Windows
 *
 * I also tested compilation of the library as well compilation and running
 * the sample executable in Windows. I will describe a procedure on how this
 * can be realized for Microsoft Visual C++ Express 2010. This procedure may be
 * adopted for compilation with Microsoft Visual Studio 2012 Express with some
 * minor but intuitive modifications.
 *
 * First download the current version of <code>thimble-yyyy.mm.dd.zip</code> from
 * <a href="http://www.stochastik.math.uni-goettingen.de/biometrics/thimble">
 * here</a> and extract the content to the folder
 * <pre>C:\\...\\thimble-yyyy.mm.dd.</pre>
 * Now, open Visual C++ Express and
 * <pre>
 * Click "File->New->Project"
 *
 * 		Choose "Win32 Project"
 * 		and let the project "Name" be 'thimble'
 * 		Click "OK" button
 *
 * 		Click "Next" button
 *
 * 		Select "Empty project"
 * 		Choose "Application type" 'Static library'
 * 		Deselect "Precompiled headers"
 * 		Click "Finish" button
 * </pre>
 * Next, copy all files that are contained in the extracted
 * <pre>
 *  C:\\...\\thimble-yyyy.mm.dd\\src\\
 * </pre>
 * to the
 * <pre>
 * Source Folder
 * </pre>
 * in the Visual C++ Express window. Afterwards, we must specify the folder
 * containing the library's include-files:
 * <pre>
 *  Click "Project->thimble Properties..."
 *	Select "Configuration Properties->C/C++->General"
 *	Edit "Additional Include Directories"
 *		Add the folder 'C:\\...\\thimble-yyyy.mm.dd\\include' to the list
 *		Click "OK" button
 *	Click "Apply" button and then "OK" button
 * </pre>
 * Finally, to compile the library we may build the project which can be
 * achieved by a right-click on the folder <code>thimble</code> in the Solution
 * explorer and then clicking <code>Build</code>. This will create the file
 * <pre>
 * C:\\...\\Visual Studio 2010\\Projects\\thimble\\Debug\\thimble.lib
 * </pre>
 *
 * @subsection sec_windows_executable Compiling Windows Executables
 *
 * I will describe the procedure of compiling executables that use THIMBLE
 * at the example of the file
 * <pre>C:\\...\\thimble-yyyy.mm.dd\\tarpSample.cpp.</pre>
 * The steps are similar to the steps for compiling the library:
 * <pre>
 *  Click "File->New->Project"
 *
 *    Choose "Win32 Console Application" and let the project "Name" be 'bfattack'
 *    Click on "OK" button
 *
 *      Click "Next" button"
 *
 *      Choose "Application type" 'Console application'
 *      Select "Empty project"
 *      Deselect "Precompiled header"
 *      Click "Finish" button
 * </pre>
 * Copy the file
 * <pre>
 *  C:\\...\\thimble-yyyy.mm.dd\\tarpSample.cpp
 * </pre>
 * to the
 * <pre>
 *  Source Files
 * </pre>
 * folder in the Visual C++ Express window. Again, we need to specify a
 * folder to THIMBLE's include-files:
 * <pre>
 *  Click "Project->bfattack Properties..."
 *      Select "Configuration Properties->C/C++->General"
 *      Edit "Additional Include Directories"
 *      Add the folder 'C:\\...\\thimble-yyyy.mm.dd\\include' to the list
 * </pre>
 * In addition we have to specify the path to the folder containing the static
 * library <code>thimble.lib</code>:
 * <pre>
 *  	Select "Configuration Properties->Linker->General"
 *  	Edit "Additional Library Directories"
 *  	Add the "C:\\...\\Visual Studio 2010\\Projects\\thimble\\Debug" folder to the list
 * </pre>
 * Also, we have to link the library <code>thimble.lib</code>:
 * <pre>
 *  	Select "Configuration Properties->Linker->Input"
 *  	Edit "Additional Dependencies"
 *  	Add the line 'thimble.lib'
 *  	Click "OK" button
 * </pre>
 * Finally, we apply the changes via
 * <pre>
 *  Click "Apply" button and the "OK" button
 * </pre>
 * Now, we can build the executable by a right-click on the folder
 * <code>thimble</code> in the Solution explorer and then clicking
 * <code>Build</code>. This will create the file
 * <pre>
 *  C:\\...\\Visual Studio 2010\\Projects\\bfattack\\Debug\\tarpSample.exe
 * </pre>
 *
 * @subsection sec_running_windows_executables Running Windows Executables
 *
 * Open the command prompt, e.g., by clicking on the Start menu and entering
 * <pre>cmd</pre> Then cd to the folder
 * <pre>C:\\...\\Visual Studio 2010\\Projects\\bfattack\\Debug\\</pre>
 * by typing
 * <pre>
 *  cd C:\\...\\Visual Studio 2010\\Projects\\bfattack\\Debug
 * </pre>
 * and pressing ENTER afterwards. Now, run
 * <pre>
 *  tarpSample.exe fingerprint.pgm
 * </pre>
 * (make sure that a fingerprint image file in PGM format
 * <code>fingerprint.pgm</code> exists along with the executable
 * <code>tarpSample.exe</code>); this should output either output a directed
 * reference point estimation of the form
 * <pre>
 *  136.851968181911871625 228.042596947062747859 0.173504083499781386335
 * </pre>
 * or
 * <pre>
 *  FAILURE TO ALIGN
 * </pre>
 *
 * @attention
 *   We may also want to build Release versions of the library and the
 *   executables. Therefore, we have to repeat the just-described compilation
 *   steps but by selecting the Release configuration in Visual C++ Express's
 *   toolbar. Note that if we run the Debug versions of the executables that
 *   they are significantly slower than their Release counterparts.
 *
 * @section sec_knownproblems Known Problems
 *
 * <ul>
 *  <li>
 *   On some <code>Intel Atom</code> architectures numerical results can
 *   be unrobust due to low precision of the FPU. In particular,
 *   directed reference point estimation
 *   via \link thimble::Fingerprint::hasDirectedReferencePoint()\endlink
 *   may return <code>false</code> in way much too many cases as compared
 *   to architectures with more powerful FPU.
 *  </li>
 * </ul>
 *
 * @section sec_changelog Change-Log
 *
 * @subsection sec_changelog_current THIMBLE-2015.06.01
 *
 * <ul>
 *  <li>
 *   Added and tested support for MinGW compilers to compile THIMBLE in a
 *   UNIX-way on Windows platforms.
 *  </li>
 *  <li>
 *   New class \link thimble::FuzzyVault\endlink of which objects represent
 *   general instances of an improved fuzzy vault scheme.
 *  </li>
 *  <li>
 *   New class \link thimble::CTools\endlink that provides convenience methods
 *   which are added as they are needed for a specific purpose in THIMBLE.
 *  </li>
 *  <li>
 *   Added the \link thimble::BigInteger::gcd()\endlink method for
 *   computing the greatest common divisor between two \link thimble::BigInteger
 *   BigIntegers\endlink.
 *  </li>
 *  <li>
 *   Added the \link thimble::BigInteger::log2()\endlink function for
 *   approximating the binary logarithm of a \link thimble::BigInteger
 *   BigInteger\endlink.
 *  </li>
 *  <li>
 *   Added the static \link thimble::BigInteger::expand()\endlink function
 *   for approximating the fraction of two \link thimble::BigInteger
 *   BigIntegers\endlink as a <code>long double</code>.
 *  </li>
 *  <li>
 *   Fixed a performance issue in
 *   the \link thimble::BigInteger::binomial()\endlink function by accounting
 *   for the equality \f${n\choose k}={n\choose n-k}\f$.
 *  </li>
 *  <li>
 *   Fixed a bug in the \link thimble::Fingerprint::setMinutiaeView()\endlink
 *   method which simply appended the
 *   specified \link thimble::MinutiaeView\endlink to
 *   the \link thimble::MinutiaeRecord\endlink object instead of overwriting
 *   its first \link thimble::MinutiaeView\endlink
 *  </li>
 *   <li>
 *   Fixed a bug that concerned the support for Win32 DLLs in which the operator
 *   << for writing text representations of
 *   a \link thimble::Minutia\endlink
 *   and \link thimble::MinutiaeView\endlink were not exported.
 *  </li>
 * </ul>
 *
 * @subsection sec_changelog_20140907 THIMBLE-2014.09.07
 *
 * <ul>
 *  <li>
 *   Added support for the library to be compiled as a Win32 DLL.
 *   To enable Win32 DLL support uncomment the line
 *   <code>// \#define THIMBLE_BUILD_DLL</code> in the file 'dllcompat.h'
 *  </li>
 *  <li>
 *   Added the static \link thimble::SmallBinaryField::binary()\endlink
 *   function return a constant reference to a finite field of size 2.
 *  </li>
 *  <li>
 *   The 'thimble::leftShift' and 'thimble::rightShift' methods operating
 *   on objects of the \link thimble::SmallBinaryFieldPolynomial\endlink
 *   class has been falsely provided as static functions. This has been
 *   fixed.
 *  </li>
 *  <li>
 *   Added \link thimble::Fingerprint::setMinutiaeView()
 *   setMinutiaeView()\endlink for manual specification of
 *   minutiae templates for a \link thimble::Fingerprint Fingerprint\endlink
 *   object.
 *  </li>
 *  <li>
 *   Changed the use of <code>std::complex<double></code> objects
 *   within \link thimble::TentedArchModel TentedArchModel\endlink which fixes
 *   compilation errors with some older <code>gcc</code> versions.
 *  </li>
 * </ul>
 *
 * @subsection sec_changelog_20140824 THIMBLE-2014.08.24
 *
 * Significant changes. Do not expect backwards compability.
 *
 * @subsection sec_changelog_20130430 THIMBLE-2013.04.30
 * First release.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_ALL_H_
#define THIMBLE_ALL_H_

#include <thimble/ecc/all.h>
#include <thimble/finger/all.h>
#include <thimble/image/all.h>
#include <thimble/math/all.h>
#include <thimble/misc/all.h>
#include <thimble/security/all.h>

#endif /* THIMBLE_ALL_H_ */
