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
 * @file BinaryMatrix.h
 *
 * @brief
 *            Provides a class for representing and computing
 *            with binary matrices of arbitrary dimension.
 *
 * @author Benjamin Tams
 */

#ifndef THIMBLE_BINARYMATRIX_H
#define THIMBLE_BINARYMATRIX_H

#include <stdint.h>
#include <fstream>

#include <thimble/dllcompat.h>

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Instances of this class represent matrices with
     *            binary/bit entries.
     */
    class THIMBLE_DLL BinaryMatrix {

    private:

        /**
         * @brief
         *            The data containing the entries of this matrix.
         *
         * @see getData()
         */
        uint32_t *data;

        /**
         * @brief
         *            The number of rows of this matrix.
         */
        int m;

        /**
         * @brief
         *            The number of columns of this matrix.
         */
        int n;

        /**
         * @brief
         *            The number of 32-bit integers needed to encode an entire
         *            row of the matrix.
         */
        int N;

        /**
         * @brief
         *            Changes the attributes of this matrix such that
         *            it becomes of the specified dimension.
         *
         * @details
         *            If the change of the dimension causes an reallocation,
         *            the content of the data will not be changed furthermore.
         *            The programmer should ensure that the data will be
         *            initialized with valid data. Therefore, the method is
         *            declared private. However, the public
         *            \link setDimension(int,int)\endlink method can be used
         *            in order to ensure that the matrix' entries are
         *            initialized as zeros.
         *
         * @param m
         *            Positive integer specifying the new number of rows.
         *
         * @param n
         *            Positive integer specifying the new number of columns.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void prepareDim( int m , int n );

    public:

        /**
         * @brief
         *            Creates a binary zero matrix of the specified dimension.
         *
         * @details
         *            This constructor also serves as the standard constructor
         *            which, by omitting to specify the parameters <i>m</i>
         *            and <i>n</i>, creates an empty matrix of dimension
         *            \f$0\times 0\f$.
         *
         * @param m
         *            The number of rows of the matrix.
         *
         * @param n
         *            The number of columns of the matrix.
         *
         * @warning
         *            If <i>m</i> or <i>n</i> is negative, an error message
         *            is printed to <code>stderr</code> and the program exits
         *            with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        BinaryMatrix( int m = 0 , int n = 0 );


        /**
         * @brief
         *            Copy constructor.
         *
         * @param mat
         *            The matrix of which a copy is created.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        BinaryMatrix( const BinaryMatrix & mat );

        /**
         * @brief
         *            Destructor.
         *
         * @details
         *            Frees all data needed to represent the matrix.
         */
        ~BinaryMatrix();

        /**
         * @brief
         *            Assignment operator (procedural version).
         *
         * @details
         *            Replaces the data held by this matrix by the
         *            data of <code>mat</code>.
         *
         * @param mat
         *            The matrix assigned to this matrix.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void assign( const BinaryMatrix & mat );

        /**
         * @brief
         *           Assignment operator.
         *
         * @details
         *            Replaces the data held by this matrix by the
         *            data of <code>mat</code>.
         *
         * @param mat
         *            The matrix assigned to this matrix.
         *
         * @return
         *            A reference to this matrix (after assignment).
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryMatrix & operator=( const BinaryMatrix & mat ) {

            assign(mat);
            return *this;
        }

        /**
         * @brief
         *            Returns the number of rows of this matrix.
         *
         * @return
         *            The number of rows of this matrix.
         */
        inline int numRows() const {

            return this->m;
        }

        /**
         * @brief
         *            Returns the number of columns of this matrix.
         *
         * @return
         *            The number of columns of this matrix.
         */
        inline int numCols() const {

            return this->n;
        }

        /**
         * @brief
         *            Changes the attributes of this matrix such that
         *            it becomes of the specified dimension.
         *
         * @details
         *            All entries of the matrix are overwritten such that
         *            it becomes the zero <i>m*n</i>-matrix.
         *
         * @param m
         *            New number of rows.
         *
         * @param n
         *            New number of columns.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void setDimension( int m , int n );

        /**
         * @brief
         *            Access the value at the specified position of
         *            the matrix.
         *
         * @param i
         *            Index of accessed row.
         *
         * @param j
         *            Index of accessed column.
         *
         * @return
         *            <code>true</code> if the accessed value equals
         *            1; otherwise <code>false</code>.
         *
         * @warning
         *            If <i>i</i> or <i>j</i> is negative or if <i>i</i>
         *            of <i>j</i> is greater than \link numRows()\endlink
         *            or \link numCols()\endlink, respectively, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        bool getAt( int i , int j ) const;

        /**
         * @brief
         *            Ensures that the value at the specified position is
         *            equals 1.
         *
         * @param i
         *            Index of accessed row.
         *
         * @param j
         *            Index of accessed column.
         *
         * @warning
         *            If <i>i</i> or <i>j</i> is negative or if <i>i</i>
         *            of <i>j</i> is greater than \link numRows()\endlink
         *            or \link numCols()\endlink, respectively, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void setAt( int i , int j );

        /**
         * @brief
         *            Ensures that the value at the specified position is
         *            equals 0.
         *
         * @param i
         *            Index of accessed row.
         *
         * @param j
         *            Index of accessed column.
         *
         * @warning
         *            If <i>i</i> or <i>j</i> is negative or if <i>i</i>
         *            of <i>j</i> is greater than \link numRows()\endlink
         *            or \link numCols()\endlink, respectively, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        void clearAt( int i , int j );

        /**
         * @brief
         *            Set the value of the matrix at a specified
         *            position to a specified value.
         *
         * @param i
         *            Row position.
         *
         * @param j
         *            Column position.
         *
         * @param c
         *            <code>true</code> if the entry at <i>(i,j)</i> is
         *            to be set to 1; otherwise, if <i>(i,j)</i> is set
         *            to 0, the value of <code>c</code> is <code>false</code>.
         *
         * @warning
         *            If <i>i</i> or <i>j</i> is negative or if <i>i</i>
         *            of <i>j</i> is greater than \link numRows()\endlink
         *            or \link numCols()\endlink, respectively, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline void setAt( int i , int j , bool c ) {

            if ( c ) {
                setAt(i,j);
            } else {
                clearAt(i,j);
            }
        }

        /**
         * @brief
         *            Sets all entries of the matrix to 0.
         */
        void setZero();

        /**
         * @brief
         *            Tests whether this matrix has only zero entries.
         *
         * @return
         *            <code>true</code> if each entry of the matrix is
         *            zero; otherwise, if there exists a non-zero entry,
         *            the function returns <code>false</code>.
         */
        bool isZero() const;

        /**
         * @brief
         *            Sets this matrix to the identity matrix.
         *
         * @warning
         *            If the number of rows is different from the number
         *            of columns, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         */
        void setIdentity();

        /**
         * @brief
         *            Exchange two rows of the binary matrix.
         *
         * @param i0
         *            The index of the first row replaced by the second
         *            row.
         *
         * @param i1
         *            The index of the second row replaced by the first
         *            row.
         *
         * @warning
         *            If <i>i0</i> or <i>i1</i> is negative or if <i>i0</i>
         *            or <i>i1</i> is greater than or equals the number
         *            of rows of this matrix, an error message is printed to
         *            <code>stderr</code> and the program exits with status
         *            'EXIT_FAILURE'.
         */
        void exchangeRows( int i0 , int i1 );

        /**
         * @brief
         *            Exchange two columns of the binary matrix.
         *
         * @param j0
         *            The index of the first column replaced by the second
         *            column.
         *
         * @param j1
         *            The index of the second column replaced by the first
         *            column.
         *
         * @warning
         *            If <i>j0</i> or <i>j1</i> is negative or if <i>j0</i>
         *            or <i>j1</i> is greater than or equals the number
         *            of columns of this matrix, an error message is printed
         *            to <code>stderr</code> and the program exits with
         *            status 'EXIT_FAILURE'.
         */
        void exchangeCols( int j0 , int j1 );

        /**
         * @brief
         *            Adds a row of the matrix to another row.
         *
         * @param i0
         *            The index of the row added to the <i>i1</i>th row.
         *
         * @param i1
         *            The index of the row to which the <i>i0</i>th row
         *            is added.
         *
         * @param j
         *            Indicates that the entries of the <i>i0</i>th row
         *            at column index smaller <i>j</i> are already zero
         *            and do therfore not need to be added to the <i>i1</i>th
         *            row. If there are entries in the <i>i0</i>th row and
         *            column index smaller <i>j</i>, the result of the
         *            method is undocumented.
         *
         * @warning
         *            If <i>i0</i>, <i>i1</i> or <i>j</i> is negative or
         *            if <i>i0</i> or <i>i1</i> is greater than or equals
         *            the number of rows of this matrix or if <i>j</i>
         *            is greater or equals the number of columns, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILIURE'.
         */
        void addRow( int i0 , int i1 , int j = 0 );

        /**
         * @brief
         *            Makes the matrix an upper triangular matrix using
         *            elementary row and column operations.
         *
         * @return
         *            The rank of the matrix.
         */
        int gauss();

        /**
         * @brief
         *            Returns the rank of the matrix.
         *
         * @return
         *            The rank of the matrix.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        int rank() const;

        /**
         * @brief
         *            Overwrites the coefficient of this matrix with
         *            random values.
         *
         * @param tryRandom
         *           If <code>true</code>, the method uses a cryptographic
         *           number generator if available on the system; otherwise,
         *           the method wraps around the standard <code>rand()</code>
         *           function.
         */
        void random( bool tryRandom = false );

        /**
         * @brief
         *            Access the data containing the entries of this matrix.
         *
         * @details
         *            Let <i>m</i> be the number of rows and <i>n</i> the
         *            number of columns of this matrix. The array contains
         *            \f$m\cdot N\f$ integers of width 32 bit where
         *            \f$N=\lceil n/32\rceil\f$. The entries of the <i>i</i>th
         *            row are encoded by the \f$\lceil m/32\rceil\f$ integers
         *            <pre>
         *             data[i*N+0],data[i*N+1],...,data[i*N+ceil(m/32)-1].
         *            </pre>
         *            Consequently, if we write the 32-bit integer
         *            <code>data[i*N+k]</code> as
         *            \f[
         *             \sum_{\ell=0}^{31}b(i,k,\ell)\cdot 2^\ell
         *            \f]
         *            where \f$b_\ell\in\{0,1\}\f$, the boolean entry of the
         *            matrix at row <i>i</i> and column <i>j</i> is given by
         *            \f$b(i,\lfloor j/32\rfloor,j\%32)\f$.
         *
         *            The field can be <code>NULL</code> if the number of rows
         *            or the number of columns is 0.
         *
         * @return
         *            The data encoding the coefficients of the matrix.
         */
        inline const uint32_t *getData() const {

            return this->data;
        }

        /**
         * @brief
         *            Access the data containing the entries of this matrix.
         *
         * @details
         *            The result of this function is the same as the result of
         *            \link getData()\endlink except that the content of the
         *            matrix can be changed by changing the content of the
         *            result.
         *
         * @warning
         *            Use the result with caution and perform only changes if
         *            you really know what your are doing; otherwise, you might
         *            experience undocumented behaviour.
         *
         * @return
         *            The data encoding the coefficients of the matrix.
         */
        inline uint32_t *getData_nonconst() {
            return this->data;
        }

        /**
         * @brief
         *           Swaps the content of two matrices such the
         *           the one represents the other.
         *
         * @param a
         *           Matrix assigned to <i>b</i>.
         *
         * @param b
         *           Matrix assigned to <i>a</i>.
         */
        static void swap( BinaryMatrix & a , BinaryMatrix & b );

        /**
         * @brief
         *           Addition of two binary matrices.
         *
         * @param c
         *           Output matrix being the sum of the input
         *           <i>a</i> and <i>b</i>.
         *
         * @param a
         *           First summand.
         *
         * @param b
         *           Second summand.
         *
         * @warning
         *           If <i>a</i> and <i>b</i> have a different number of rows
         *           or a different number of columns, an error message is
         *           printed to <code>stderr</code> and the program exits
         *           with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void add
        ( BinaryMatrix & c ,
          const BinaryMatrix & a , const BinaryMatrix & b );

        /**
         * @brief
         *           Subtraction of two binary matrices.
         *
         * @details
         *           For binary matrices, subtraction is equivalent to
         *           addition and, thus, it is is actually not necessary
         *           to provide a subtraction method explicitly. However, in
         *           some situations it can be useful to distinguish addition
         *           from subtraction for better readability of code.
         *
         * @param c
         *           Output matrix being the difference between the input
         *           <i>a</i> and <i>b</i>.
         *
         * @param a
         *           Minuend matrix.
         *
         * @param b
         *           Subtrahend matrix.
         *
         * @warning
         *           If <i>a</i> and <i>b</i> have a different number of rows
         *           or a different number of columns, an error message is
         *           printed to <code>stderr</code> and the program exits
         *           with status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline static void sub
        ( BinaryMatrix & c ,
          const BinaryMatrix & a , const BinaryMatrix & b ) {

            // Addition is equivalent to subtraction for binary matrices
            BinaryMatrix::add(c,a,b);
        }

        /**
         * @brief
         *           Determines the transpose of a binary matrix.
         *
         * @param b
         *           On output, the transpose of the input matrix <i>a</i>.
         *
         * @param a
         *           The input matrix.
         *
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void transpose
        ( BinaryMatrix & b , const BinaryMatrix & a );

        /**
         * @brief
         *           Computes the product of a matrix times the transpose
         *           of another matrix.
         *
         * @details
         *           The method performs the matrix multiplication
         *           \f[
         *            c = a\cdot b^\top.
         *           \f]
         *
         * @param c
         *           On output, the product of the input <i>a</i> times
         *           the transpose of the input <i>b</i>.
         *
         * @param a
         *           First factor matrix.
         *
         * @param b
         *           Transpose of the second factor matrix.
         *
         * @warning
         *           If the number of columns of <i>a</i> is different from the
         *           number of columns of <i>b</i>, an error message is printed
         *           to <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul_t
        ( BinaryMatrix & c ,
          const BinaryMatrix & a , const BinaryMatrix & b );

        /**
         * @brief
         *           Performs a matrix multipliciation.
         *
         * @details
         *           The method performs the matrix multipliciation
         *           \f[
         *            c = a \cdot b.
         *           \f]
         *
         * @param c
         *           One output, the product of the input <i>a</i> times
         *           the input <i>b</i>.
         *
         * @param a
         *           First factor matrix.
         *
         * @param b
         *           Second factor matrix.
         *
         * @warning
         *           If the number of columns of <i>a</i> is different from the
         *           number of rows of <i>b</i>, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void mul
        ( BinaryMatrix & c ,
          const BinaryMatrix & a , const BinaryMatrix & b );

        /**
         * @brief
         *           Negates a binary matrix bit-wisely.
         *
         * @param B
         *           On output, the bitwise negation of <i>A</i>.
         *
         * @param A
         *           The matrix of which the bits are negated and stored
         *           in <i>B</i>.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void lnot
        ( BinaryMatrix & B , const BinaryMatrix & A );

        /**
         * @brief
         *           Cuts off a submatrix of a binary matrix.
         *
         * @param A0
         *           On output, the submatrix of <i>A</i> of dimension
         *           \f$m_0\times n_0\f$ starting from \f$(i_0,j_0)\f$.
         *
         * @param A
         *           The matrix of which a submatrix is extracted.
         *
         * @param i0
         *           The initial row index.
         *
         * @param j0
         *           The initial column index.
         *
         * @param m0
         *           The number of rows of the submatrix.
         *
         * @param n0
         *           The number of columns of the submatrix.
         *
         * @warning
         *           If \f$(i_0,j_0)\f$ or \f$(i_0+m_0),j_0+n_0)\f$ is not a
         *           valid position in <i>A</i> or if \f$m_0\f$ or \f$n_0\f$
         *           is negative, an error message is printed to
         *           <code>stderr</code> and the program exits with status
         *           'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void submatrix
        ( BinaryMatrix & A0 , const BinaryMatrix & A ,
          int i0 , int j0 , int m0 , int n0 );

        /**
         * @brief
         *           Concatenates the columns of two matrices.
         *
         * @details
         *           The function's input are two matrices \f$A\f$ and
         *           \f$B\f$ and sets
         *           \f[
         *            C=(A|B).
         *           \f]
         *
         * @param C
         *           On output, the concatenated matrix \f$(A|B)\f$.
         *
         * @param A
         *           The first columns of the concatenated matrix.
         *
         * @param B
         *           The last columns of the concatenated matrix.
         *
         * @warning
         *           If <i>A</i> and <i>B</i> have a different number
         *           of rows, an error message is printed to
         *           <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void concatCols
        ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B );

        /**
         * @brief
         *           Concatenates the rows of two matrices.
         *
         * @details
         *           The function's input are two matrices \f$A\f$ and
         *           \f$B\f$ and sets
         *           \f[
         *            C=\left(\begin{array}{c}A\\B\end{array}\right)
         *           \f]
         *
         * @param C
         *           On output, the concatenation of the rows of <i>A</i>
         *           an <i>B</i>.
         *
         * @param A
         *           The first rows of the concatenated matrix.
         *
         * @param B
         *           The last rows of the concatenated matrix.
         *
         *
         * @warning
         *           If <i>A</i> and <i>B</i> have a different number
         *           of columns, an error message is printed to
         *           <code>stderr</code> and the program exits with
         *           status 'EXIT_FAILURE'.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        static void concatRows
        ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B );

        /**
         * @brief
         *           Allows that the negative of a binary matrix
         *           can be computed via the '-'-operator.
         *
         * @details
         *           Computing the negative of a binary matrix is not
         *           really useful. However, there are cases in which
         *           the existence of the negation operator is useful,
         *           e.g., when one wants to implement templates that
         *           work for binary and non-binary matrices.
         *
         * @return
         *           The negative of this matrix.
         *
         * @warning
         *           If not enough memory could be provided, an error
         *           message is printed to <code>stderr</code> and the
         *           program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryMatrix operator-() const {
            return *this; // Just returns a copy of this matrix
        }

        /**
         * @brief
         *            Overloaded '+'-operator to perform an addition
         *            of two binary matrices.
         *
         * @details
         *            The '+'-operator allows expressions as
         *            <pre>
         *             BinaryMatrix a , b , c;
         *             ...
         *             c = a+b;
         *            </pre>
         *            and is wrapped around the \link add()\endlink
         *            function.
         *
         * @param b
         *            Summand.
         *
         * @return
         *            The sum of this matrix and <code>b</code>.
         *
         * @warning
         *            If the dimension of <code>b</code> is different
         *            from the dimension of this matrix, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryMatrix operator+( const BinaryMatrix & b ) const {

            BinaryMatrix c;
            add(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Overloaded '-'-operator to perform a
         *            subtraction/addition of two binary matrices.
         *
         * @details
         *            The '-'-operator allows expressions as
         *            <pre>
         *             BinaryMatrix a , b , c;
         *             ...
         *             c = a-b;
         *            </pre>
         *
         * @param b
         *            Subtrahend.
         *
         * @return
         *            The difference between this matrix and <code>b</code>.
         *
         * @warning
         *            If the dimension of <code>b</code> is different
         *            from the dimension of this matrix, an error message
         *            is printed to <code>stderr</code> and the program
         *            exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryMatrix operator-( const BinaryMatrix & b ) const {

            BinaryMatrix c;
            sub(c,*this,b);
            return c;
        }

        /**
         * @brief
         *            Overloaded '*'-operator to perform a binary matrix
         *            multiplication.
         *
         * @details
         *            The '*'-operator allows expressions as
         *            <pre>
         *             BinaryMatrix a , b , c;
         *             ...
         *             c = a*b;
         *            </pre>
         *            and is wrapped around the \link mul()\endlink.
         *
         * @param b
         *            Factor.
         *
         * @return
         *            The product of this matrix times <code>b</code>.
         *
         * @warning
         *            If the number of columns of this matrix is different
         *            from the number of rows of <code>b</code>, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         *
         * @warning
         *            If not enough memory could be provided, an error
         *            message is printed to <code>stderr</code> and the
         *            program exits with status 'EXIT_FAILURE'.
         */
        inline BinaryMatrix operator*( const BinaryMatrix & b ) const {

            BinaryMatrix c;
            mul(c,*this,b);
            return c;
        }

        /**
         * @brief
         *           Implementation of the += operator.
         *
         * @details
         *           The function performs
         *           <pre>
         *            add(*this,*this,mat).
         *           </pre>
         *           For more information see the
         *           \link BinaryMatrix::add()\endlink function.
         *
         * @param mat
         *           The matrix added to this matrix.
         *
         * @return
         *           A reference to this instance.
         */
        inline BinaryMatrix & operator+=( const BinaryMatrix & mat ) {

            add(*this,*this,mat);
            return *this;
        }

        /**
         * @brief
         *           Implementation of the -= operator.
         *
         * @details
         *           The function performs
         *           <pre>
         *            sub(*this,*this,mat)
         *           </pre>
         *           or, equivalently,
         *           <pre>
         *            add(*this,*this,mat).
         *           </pre>
         *           For more information see the
         *           \link BinaryMatrix::sub()\endlink function.
         *
         * @param mat
         *           The matrix added to/subtracted from this matrix.
         *
         * @return
         *           A reference to this instance.
         */
        inline BinaryMatrix & operator-=( const BinaryMatrix & mat ) {

            sub(*this,*this,mat);
            return *this;
        }

        /**
         * @brief
         *           Implementation of the *= operator.
         *
         * @details
         *           The function performs
         *           <pre>
         *            mul(*this,*this,mat).
         *           </pre>
         *           For more information see the
         *           \link BinaryMatrix::mul()\endlink function.
         *
         * @param mat
         *           The matrix with which this matrix multiplied.
         *
         * @return
         *           A reference to this instance.
         */
        inline BinaryMatrix & operator*=( const BinaryMatrix & mat ) {

            mul(*this,*this,mat);
            return *this;
        }
    };

    /**
     * @brief
     *           Swaps the content of two matrices such the
     *           the one represents the other.
     *
     * @param a
     *           Matrix assigned to <i>b</i>.
     *
     * @param b
     *           Matrix assigned to <i>a</i>.
     */
    inline void swap( BinaryMatrix & a , BinaryMatrix & b ) {

        BinaryMatrix::swap(a,b);
    }

    /**
     * @brief
     *           Addition of two binary matrices.
     *
     * @param c
     *           Output matrix being the sum of the input
     *           <i>a</i> and <i>b</i>.
     *
     * @param a
     *           First summand.
     *
     * @param b
     *           Second summand.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> have a different number of rows
     *           or a different number of columns, an error message is
     *           printed to <code>stderr</code> and the program exits
     *           with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void add
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        BinaryMatrix::add(c,a,b);
    }

    /**
     * @brief
     *           Subtraction of two binary matrices.
     *
     * @details
     *           For binary matrices, subtraction is equivalent to
     *           addition and, thus, it is is actually not necessary
     *           to provide a subtraction method explicitly. However, in
     *           some situations it can be useful to distinguish addition
     *           from subtraction for better readability of code.
     *
     * @param c
     *           Output matrix being the difference between the input
     *           <i>a</i> and <i>b</i>.
     *
     * @param a
     *           Minuend matrix.
     *
     * @param b
     *           Subtrahend matrix.
     *
     * @warning
     *           If <i>a</i> and <i>b</i> have a different number of rows
     *           or a different number of columns, an error message is
     *           printed to <code>stderr</code> and the program exits
     *           with status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void sub
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        BinaryMatrix::sub(c,a,b);
    }

    /**
     * @brief
     *           Determines the transpose of a binary matrix.
     *
     * @param b
     *           On output, the transpose of the input matrix <i>a</i>.
     *
     * @param a
     *           The input matrix.
     *
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void transpose
    ( BinaryMatrix & b , const BinaryMatrix & a ) {

        BinaryMatrix::transpose(b,a);
    }

    /**
     * @brief
     *           Determines the transpose of a binary matrix.
     *
     * @param A
     *           The input matrix of which transposed matrix is returned.
     *
     * @return
     *           The transpose of the input matrix <i>A</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix transpose( const BinaryMatrix & A ) {

        BinaryMatrix AT;
        BinaryMatrix::transpose(AT,A);
        return AT;
    }

    /**
     * @brief
     *           Computes the product of a matrix times the transpose
     *           of another matrix.
     *
     * @details
     *           The method performs the matrix multiplication
     *           \f[
     *            c = a\cdot b^\top.
     *           \f]
     *
     * @param c
     *           On output, the product of the input <i>a</i> times
     *           the transpose of the input <i>b</i>.
     *
     * @param a
     *           First factor matrix.
     *
     * @param b
     *           Transpose of the second factor matrix.
     *
     * @warning
     *           If the number of columns of <i>a</i> is different from the
     *           number of columns of <i>b</i>, an error message is printed
     *           to <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul_t
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        BinaryMatrix::mul_t(c,a,b);
    }

    /**
     * @brief
     *           Performs a matrix multipliciation.
     *
     * @details
     *           The method performs the matrix multipliciation
     *           \f[
     *            c = a \cdot b.
     *           \f]
     *
     * @param c
     *           One output, the product of the input <i>a</i> times
     *           the input <i>b</i>.
     *
     * @param a
     *           First factor matrix.
     *
     * @param b
     *           Second factor matrix.
     *
     * @warning
     *           If the number of columns of <i>a</i> is different from the
     *           number of rows of <i>b</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void mul
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        BinaryMatrix::mul(c,a,b);
    }

    /**
     * @brief
     *           Negates a binary matrix bit-wisely.
     *
     * @param B
     *           On output, the bitwise negation of <i>A</i>.
     *
     * @param A
     *           The matrix of which the bits are negated and stored
     *           in <i>B</i>.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void lnot
    ( BinaryMatrix & B , const BinaryMatrix & A ) {

    	BinaryMatrix::lnot(B,A);
    }

    /**
     * @brief
     *           Cuts off a submatrix of a binary matrix.
     *
     * @param A0
     *           On output, the submatrix of <i>A</i> of dimension
     *           \f$m_0\times n_0\f$ starting from \f$(i_0,j_0)\f$.
     *
     * @param A
     *           The matrix of which a submatrix is extracted.
     *
     * @param i0
     *           The initial row index.
     *
     * @param j0
     *           The initial column index.
     *
     * @param m0
     *           The number of rows of the submatrix.
     *
     * @param n0
     *           The number of columns of the submatrix.
     *
     * @warning
     *           If \f$(i_0,j_0)\f$ or \f$(i_0+m_0),j_0+n_0)\f$ is not a
     *           valid position in <i>A</i> or if \f$m_0\f$ or \f$n_0\f$
     *           is negative, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void submatrix
    ( BinaryMatrix & A0 , const BinaryMatrix & A ,
      int i0 , int j0 , int m0 , int n0 ) {

        BinaryMatrix::submatrix(A0,A,i0,j0,m0,n0);
    }

    /**
     * @brief
     *           Cuts off a submatrix of a binary matrix.
     *
     * @param A
     *           The matrix of which a submatrix is extracted.
     *
     * @param i0
     *           The initial row index.
     *
     * @param j0
     *           The initial column index.
     *
     * @param m0
     *           The number of rows of the submatrix.
     *
     * @param n0
     *           The number of columns of the submatrix.
     *
     * @return
     *           The submatrix of <i>A</i> of dimension
     *           \f$m_0\times n_0\f$ starting from \f$(i_0,j_0)\f$.
     *
     * @warning
     *           If \f$i_0,j_0,m_0\f$, or \f$n_0\f$ is negative or if
     *           the accessed submatrix does not fit in the matrix
     *           <i>A</i>, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix submatrix
    ( const BinaryMatrix & A ,
      int i0 , int j0 , int m0 , int n0 ) {

        BinaryMatrix A0;
        BinaryMatrix::submatrix(A0,A,i0,j0,m0,n0);
        return A0;
    }


    /**
     * @brief
     *           Concatenates the columns of two matrices.
     *
     * @details
     *           The function's input are two matrices \f$A\f$ and
     *           \f$B\f$ and sets
     *           \f[
     *            C=(A|B).
     *           \f]
     *
     * @param C
     *           On output, the concatenated matrix \f$(A|B)\f$.
     *
     * @param A
     *           The first columns of the concatenated matrix.
     *
     * @param B
     *           The last columns of the concatenated matrix.
     *
     * @warning
     *           If <i>A</i> and <i>B</i> have a different number
     *           of rows, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void concatCols
    ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B ) {

        BinaryMatrix::concatCols(C,A,B);
    }

    /**
     * @brief
     *           Concatenates the rows of two matrices.
     *
     * @details
     *           The function's input are two matrices \f$A\f$ and
     *           \f$B\f$ and sets
     *           \f[
     *            C=\left(\begin{array}{c}A\\B\end{array}\right)
     *           \f]
     *
     * @param C
     *           On output, the concatenation of the rows of <i>A</i>
     *           an <i>B</i>.
     *
     * @param A
     *           The first rows of the concatenated matrix.
     *
     * @param B
     *           The last rows of the concatenated matrix.
     *
     *
     * @warning
     *           If <i>A</i> and <i>B</i> have a different number
     *           of columns, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline void concatRows
    ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B ) {

        BinaryMatrix::concatRows(C,A,B);
    }

    /**
     * @brief
     *           Concatenates the columns of two matrices.
     *
     * @details
     *           The function's input are two matrices \f$A\f$ and
     *           \f$B\f$ and the function returns the matrix
     *           \f[
     *            (A|B).
     *           \f]
     *
     * @param A
     *           The first columns of the concatenated matrix.
     *
     * @param B
     *           The last columns of the concatenated matrix.
     *
     * @return
     *           \f$(A|B)\f$.
     *
     * @warning
     *           If <i>A</i> and <i>B</i> have a different number
     *           of rows, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix concatCols
    (  const BinaryMatrix & A , const BinaryMatrix & B ) {

        BinaryMatrix C;
        BinaryMatrix::concatCols(C,A,B);
        return C;
    }


    /**
     * @brief
     *           Concatenates the rows of two matrices.
     *
     * @details
     *           The function's input are two matrices \f$A\f$ and
     *           \f$B\f$ and the function returns the matrix
     *           \f[
     *            \left(\begin{array}{c}A\\B\end{array}\right)
     *           \f]
     *
     * @param A
     *           The first rows of the concatenated matrix.
     *
     * @param B
     *           The last rows of the concatenated matrix.
     *
     * @return
     *           See details above.
     *
     * @warning
     *           If <i>A</i> and <i>B</i> have a different number
     *           of columns, an error message is printed to
     *           <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    inline BinaryMatrix concatRows
    ( const BinaryMatrix & A , const BinaryMatrix & B ) {

        BinaryMatrix C;
        BinaryMatrix::concatRows(C,A,B);
        return C;
    }

    /**
     * @brief
     *            Prints a representation of a binary matrix to the specified
     *            output stream.
     *
     * @param out
     *            The output stream.
     *
     * @param mat
     *            The binary matrix being printed to <code>out</code>.
     *
     * @return
     *            A reference to <code>out</code>.
     */
    THIMBLE_DLL std::ostream & operator<<
    ( std::ostream & out , const BinaryMatrix & mat );
}

#endif /* THIMBLE_BINARYMATRIX_H */
