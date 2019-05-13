/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint Based Biometric Cryptosystems.
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
 * @file LinAlgTools.h
 *
 * @brief
 *            Provides a class and utility functions/methods that might be
 *            useful for purposes related to linear algebra.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_LINALGTOOLS_H
#define THIMBLE_LINALGTOOLS_H

#include <thimble/dllcompat.h>
#include <thimble/math/linalg/BinaryVector.h>
#include <thimble/math/linalg/BinaryMatrix.h>
#include <thimble/math/Permutation.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Class implementing static functions that might
     *            be useful for purposes with relation to linear
     *            algebra.
     */
    class THIMBLE_DLL LinAlgTools {

    private:

       /**
        * @brief
        *            Private standard constructor.
        *
        * @details
        *            The standard constructor is private to avoid users
        *            to create instance of this tool class that provides only
        *            static functions.
        */
        inline LinAlgTools() { }

    public:

        /**
         * @brief
         *           Performs a matrix vector multiplication over the
         *           binary field.
         *
         * @details
         *           The method performs
         *           \f[
         *            y \leftarrow A\cdot x
         *           \f]
         *           where the vectors are interpreted as column vectors.
         *
         * @param y
         *           Output vector in which the result of the multiplication
         *           is stored.
         *
         * @param A
         *           Binary matrix.
         *
         * @param x
         *           Binary vector.
         *
         * @warning
         *           If the
         *           \link BinaryMatrix::numCols() number of columns\endlink
         *           of <i>A</i> and the
         *           \link BinaryVector::getLength() length\endlink of
         *           <i>x</i> are different, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( BinaryVector & y ,
          const BinaryMatrix & A , const BinaryVector & x );

        /**
         * @brief
         *           Multiplies a permutation matrix with a binary matrix.
         *
         * @param B
         *           On output, the input matrix <i>A</i> but with reordered
         *           rows as specified by <i>P</i>.
         *
         * @param P
         *           The specified permutation matrix being multiplied with
         *           <i>A</i>.
         *
         * @param A
         *           Binary matrix to which rows the permutation <i>P</i> is
         *           applied.
         *
         * @warning
         *           If the dimension of <i>P</i> is different from the number
         *           of rows of <i>A</i>, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A );

        /**
         * @brief
         *           Multiplies the inverse of a permutation matrix with a
         *           binary matrix.
         *
         * @param B
         *           On output, the input matrix <i>A</i> but with reordered
         *           rows as specified by the inverse of <i>P</i>.
         *
         * @param P
         *           The specified permutation matrix of which inverse is
         *           multiplied with <i>A</i>.
         *
         * @param A
         *           Binary matrix to which rows the inverse permutation
         *           <i>P</i> is applied.
         *
         * @warning
         *           If the dimension of <i>P</i> is different from the number
         *           of rows of <i>A</i>, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void imul
        ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A );

        /**
         * @brief
         *           Multiplies a permutation matrix with a binary vector.
         *
         * @param y
         *           On output, the input vector <i>y</i> but with reordered
         *           entries as specified by <i>P</i>.
         *
         * @param P
         *           The specified permutation matrix being multiplied with
         *           <i>x</i>.
         *
         * @param x
         *           Binary vector of which entries are reordered by the
         *           permutation <i>P</i>.
         *
         * @warning
         *           If the dimension of <i>P</i> is different from the length
         *           of <i>x</i>, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( BinaryVector & y , const Permutation & P ,const BinaryVector & x );

        /**
         * @brief
         *           Computes a basis of the image of a matrix.
         *
         * @param A
         *           The matrix of which image's basis is computed.
         *
         * @param B
         *           On output, the column vectors of <i>B</i> form a
         *           basis of the image of <i>A</i>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void image( BinaryMatrix & B , const BinaryMatrix & A );

        /**
         * @brief
         *           Computes a basis of the kernel of a binary matrix.
         *
         * @param B
         *           On output, the column vectors of <i>B</i> form a
         *           basis of the kernel of <i>A</i>.
         *
         * @param A
         *           The input matrix of which kernel's basis is computed.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void kernel( BinaryMatrix & B , const BinaryMatrix & A );

        /**
         * @brief
         *           Solves a linear equation.
         *
         * @details
         *           The method attempts to find a vector <i>x</i> such that
         *           \f$A\cdot x=y\f$ holds. If solvable, the function sets
         *           <i>x</i> to a valid solution and returns <code>true</code>;
         *           otherwise, if the equation cannot be solved, the vector
         *           <i>x</i> will be left unchanged and the function returns
         *           <code>false</code>.
         *
         * @param x
         *           Will be set to a solution of the linear system, if any.
         *
         * @param A
         *           Encodes the relation of the linear equation.
         *
         * @param y
         *           Encodes the data of the linear equation.
         *
         * @return
         *           <code>true</code> if a valid solution exists and
         *           <code>false</code> otherwise.
         *
         * @warning
         *           If <i>y</i> has length different from the number of columns
         *           of <i>A</i>, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static bool solve
        ( BinaryVector & x , const BinaryMatrix & A ,
          const BinaryVector & y );

        /**
         * @brief
         *           Checks whether a vector is an element of the kernel
         *           of a matrix.
         *
         * @param v
         *           The vector for which it is tested whether it is an
         *           element of the kernel of <i>A</i>.
         *
         * @param A
         *           The matrix for which it is testet whether <i>v</i> is
         *           an element of its kernel.
         *
         * @return
         *           <code>true</code> if <i>v</i> is an element of the
         *           kernel of <i>A</i>; otherwise, <code>false</code>.
         *
         * @warning
         *           If the
         *           \link BinaryMatrix::numCols() number of columns\endlink
         *           of <i>A</i> and the
         *           \link BinaryVector::getLength() length\endlink of
         *           <i>v</i> are different, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         */
        static bool isInKernel
        ( const BinaryVector & v , const BinaryMatrix & A );


    };


    /**
     * @brief
     *           Performs a matrix vector multiplication over the
     *           binary field.
     *
     * @details
     *           The method performs
     *           \f[
     *            y \leftarrow A\cdot x
     *           \f]
     *           where the vectors are interpreted as column vectors.
     *
     * @param y
     *           Output vector in which the result of the multiplication
     *           is stored.
     *
     * @param A
     *           Binary matrix.
     *
     * @param x
     *           Binary vector.
     *
     * @warning
     *           If the
     *           \link BinaryMatrix::numCols() number of columns\endlink
     *           of <i>A</i> and the
     *           \link BinaryVector::getLength() length\endlink of
     *           <i>x</i> are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( BinaryVector & y ,
      const BinaryMatrix & A , const BinaryVector & x ) {

        LinAlgTools::mul(y,A,x);
    }

    /**
     * @brief
     *           Multiplies a permutation matrix with a binary matrix.
     *
     * @param B
     *           On output, the input matrix <i>A</i> but with reordered
     *           rows as specified by <i>P</i>.
     *
     * @param P
     *           The specified permutation matrix being multiplied with
     *           <i>A</i>.
     *
     * @param A
     *           Binary matrix to which rows the permutation <i>P</i> is
     *           applied.
     *
     * @warning
     *           If the dimension of <i>P</i> is different from the number
     *           of rows of <i>A</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A ) {

        LinAlgTools::mul(B,P,A);
    }

    /**
     * @brief
     *           Multiplies the inverse of a permutation matrix with a
     *           binary matrix.
     *
     * @param B
     *           On output, the input matrix <i>A</i> but with reordered
     *           rows as specified by the inverse of <i>P</i>.
     *
     * @param P
     *           The specified permutation matrix of which inverse is
     *           multiplied with <i>A</i>.
     *
     * @param A
     *           Binary matrix to which rows the inverse permutation
     *           <i>P</i> is applied.
     *
     * @warning
     *           If the dimension of <i>P</i> is different from the number
     *           of rows of <i>A</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void imul
    ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A ) {

        LinAlgTools::imul(B,P,A);
    }

    /**
     * @brief
     *           Multiplies a permutation matrix with a binary vector.
     *
     * @param y
     *           On output, the input vector <i>y</i> but with reordered
     *           entries as specified by <i>P</i>.
     *
     * @param P
     *           The specified permutation matrix being multiplied with
     *           <i>x</i>.
     *
     * @param x
     *           Binary vector of which entries are reordered by the
     *           permutation <i>P</i>.
     *
     * @warning
     *           If the dimension of <i>P</i> is different from the length
     *           of <i>x</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( BinaryVector & y , const Permutation & P ,const BinaryVector & x ) {

        LinAlgTools::mul(y,P,x);
    }

    /**
     * @brief
     *           Computes a basis of the image of a matrix.
     *
     * @details
     *           The implementation performs a Gaussian elimination to
     *           compute a matrix of full rank whose column vectors
     *           form a basis for the image of <i>A</i> and stores
     *           the result in <i>B</i>.
     *
     * @param A
     *           The matrix of which image's basis is computed.
     *
     * @param B
     *           On output, the column vectors of <i>B</i> form a
     *           basis of the image of <i>A</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void image( BinaryMatrix & B , const BinaryMatrix & A ) {

        LinAlgTools::image(B,A);
    }

    /**
     * @brief
     *           Computes a basis of the image of a matrix.
     *
     * @param A
     *           The matrix of which image's basis is computed.
     *
     * @return
     *           A matrix of which column vectors form a
     *           basis of the image of <i>A</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix image( const BinaryMatrix & A ) {

        BinaryMatrix B;
        LinAlgTools::image(B,A);
        return B;
    }

    /**
     * @brief
     *           Computes a basis of the kernel of a binary matrix.
     *
     * @param B
     *           On output, the column vectors of <i>B</i> form a
     *           basis of the kernel of <i>A</i>.
     *
     * @param A
     *           The input matrix of which kernel's basis is computed.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void kernel( BinaryMatrix & B , const BinaryMatrix & A ) {

        LinAlgTools::kernel(B,A);
    }

    /**
     * @brief
     *           Computes a basis of the kernel of a binary matrix.
     *
     * @param A
     *           The input matrix of which kernel's basis is computed.
     *
     * @return
     *           A binary matrix of which column vectors form a
     *           basis of the kernel of <i>A</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix kernel( const BinaryMatrix & A ) {

        BinaryMatrix B;
        LinAlgTools::kernel(B,A);
        return B;
    }

    /**
     * @brief
     *           Solves a linear equation.
     *
     * @details
     *           The method attempts to find a vector <i>x</i> such that
     *           \f$A\cdot x=y\f$ holds. If solvable, the function sets
     *           <i>x</i> to a valid solution and returns <code>true</code>;
     *           otherwise, if the equation cannot be solved, the vector
     *           <i>x</i> will be left unchanged and the function returns
     *           <code>false</code>.
     *
     * @param x
     *           Will be set to a solution of the linear system, if any.
     *
     * @param A
     *           Encodes the relation of the linear equation.
     *
     * @param y
     *           Encodes the data of the linear equation.
     *
     * @return
     *           <code>true</code> if a valid solution exists and
     *           <code>false</code> otherwise.
     *
     * @warning
     *           If <i>y</i> has length different from the number of columns
     *           of <i>A</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline bool solve
    ( BinaryVector & x , const BinaryMatrix & A ,
      const BinaryVector & y ) {

    	return LinAlgTools::solve(x,A,y);
    }

    /**
     * @brief
     *           Checks whether a vector is an element of the kernel
     *           of a matrix.
     *
     * @param v
     *           The vector for which it is tested whether it is an
     *           element of the kernel of <i>A</i>.
     *
     * @param A
     *           The matrix for which it is testet whether <i>v</i> is
     *           an element of its kernel.
     *
     * @return
     *           <code>true</code> if <i>v</i> is an element of the
     *           kernel of <i>A</i>; otherwise, <code>false</code>.
     *
     * @warning
     *           If the
     *           \link BinaryMatrix::numCols() number of columns\endlink
     *           of <i>A</i> and the
     *           \link BinaryVector::getLength() length\endlink of
     *           <i>v</i> are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     */
    inline bool isInKernel
    ( const BinaryVector & v , const BinaryMatrix & A ) {

    	return LinAlgTools::isInKernel(v,A);
    }

    /**
     * @brief
     *           Performs a matrix vector multiplication over the
     *           binary field.
     *
     * @details
     *           The method performs
     *           \f[
     *            y \leftarrow A\cdot x
     *           \f]
     *           and returns <i>y</i> (here the vectors are intepreted as
     *           column vectors)
     *
     * @param A
     *           Binary matrix.
     *
     * @param x
     *           Binary vector.
     *
     * @return
     *           \f$A\cdot x\f$.
     *
     * @warning
     *           If the
     *           \link BinaryMatrix::numCols() number of columns\endlink
     *           of <i>A</i> and the
     *           \link BinaryVector::getLength() length\endlink of
     *           <i>x</i> are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryVector operator*
    ( const BinaryMatrix & A , const BinaryVector & x ) {

        BinaryVector y;
        LinAlgTools::mul(y,A,x);
        return y;
    }

    /**
     * @brief
     *           Applies a permutation to a binary vector.
     *
     * @param P
     *           The permutation applied to <i>x</i>.
     *
     * @param x
     *           The binary vector to which <i>P</i> is applied.
     *
     * @return
     *           A vector whose entries form the reordered entries
     *           of <i>x</i> as specified by <i>P</i>.
     *
     * @warning
     *           If the dimension of <i>P</i> is different from the
     *           length of <i>x</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryVector operator*
    ( const Permutation & P , const BinaryVector & x ) {

        BinaryVector y;
        LinAlgTools::mul(y,P,x);
        return y;
    }

    /**
     * @brief
     *           Mutliplies a permutation matrix with a binary matrix.
     *
     * @param P
     *           The specified permutation matrix being multiplied with
     *           <i>A</i>.
     *
     * @param A
     *           Binary matrix to which rows the permutation <i>P</i> is
     *           applied.
     *
     * @return
     *           The product <i>P*A</i>.
     *
     * @warning
     *           If the dimension of <i>P</i> is different from the number
     *           of rows of <i>A</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix operator*
    ( const Permutation & P , const BinaryMatrix & A ) {

        BinaryMatrix B;
        LinAlgTools::mul(B,P,A);
        return B;
    }
}

#endif /* THIMBLE_LINALGTOOLS_H */
