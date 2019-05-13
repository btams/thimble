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
 * @file BinaryMatrix.cpp
 *
 * @brief
 *            Implementation of a class for representation and computation
 *            of binary matrices of arbitrary dimension as provided by the
 *            'BinaryMatrix.h' header.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include <thimble/math/MathTools.h>
#include <thimble/math/linalg/BinaryMatrix.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

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
    BinaryMatrix::BinaryMatrix( int m , int n ) {

        if ( m < 0 || n < 0 ) {
            cerr << "BinaryMatrix: negative dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        this->m = m;
        this->n = n;
        this->N = n/32+(n%32?1:0);

        if ( this->m > 0 && this->N > 0) {
            this->data = (uint32_t*)malloc
                    ( this->m * this->N * sizeof(uint32_t) );
            if ( this->data == NULL ) {
                cerr << "BinaryMatrix: out of memory." << endl;
                exit(EXIT_FAILURE);
            }
            memset(this->data,0,this->m * this->N * sizeof(uint32_t));
        } else {
           this->data = NULL;
        }
    }

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
    BinaryMatrix::BinaryMatrix( const BinaryMatrix & mat ) {

        this->m = mat.m;
        this->n = mat.n;
        this->N = mat.N;

        if ( mat.data != NULL ) {
            this->data = (uint32_t*)malloc
                    ( this->m * this->N * sizeof(uint32_t ));
            if ( this->data == NULL ) {
                cerr << "BinaryMatrix: out of memory." << endl;
                exit(EXIT_FAILURE);
            }
            memcpy(this->data,mat.data,this->m * this->N * sizeof(uint32_t) );
        } else {
            this->data = NULL;
        }
    }

    /**
     * @brief
     *            Destructor.
     *
     * @details
     *            Frees all data needed to represent the matrix.
     */
    BinaryMatrix::~BinaryMatrix() {
        free(this->data);
    }

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
    void BinaryMatrix::assign( const BinaryMatrix & mat ) {

        // Nothing to do on self-assignment
        if ( this != &mat ) {

            this->m = mat.m;
            this->n = mat.n;
            this->N = mat.N;

            if ( mat.data == NULL ) {

                free(this->data);
                this->data = NULL;

            } else {

                if ( this->data == NULL ) {
                    this->data = (uint32_t*)malloc
                            ( this->m * this->N * sizeof(uint32_t) );
                } else {
                    this->data = (uint32_t*)realloc( this->data ,
                            this->m * this->N * sizeof(uint32_t) );
                }

                if ( this->data == NULL ) {
                    cerr << "BinaryMatrix::assign: "
                         << "out of memory." << endl;
                    exit(EXIT_FAILURE);
                }

                memcpy
                (this->data,mat.data,this->m * this->N * sizeof(uint32_t));
            }

        }
    }

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
    void BinaryMatrix::prepareDim( int m , int n ) {

        if ( m == 0 || n == 0 ) {

            this->m = m;
            this->n = n;
            this->N = n/32+(n%32?1:0);

            free(this->data);
            this->data = NULL;

        } else {

            if ( this->m != m || this->n != n  ) {

                int N = n/32+(n%32?1:0);

                if ( this->data == NULL ) {
                    this->data = (uint32_t*)malloc
                             ( m * N * sizeof(uint32_t) );
                } else {
                    this->data = (uint32_t*)realloc( this->data ,
                             m * N * sizeof(uint32_t) );
                }

                if ( this->data == NULL ) {
                    cerr << "BinaryMatrix::setDim: "
                         << "out of memory." << endl;
                    exit(EXIT_FAILURE);
                }

                this->m = m;
                this->n = n;
                this->N = N;
            }
        }
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
    void BinaryMatrix::setDimension( int m , int n ) {

        if ( m < 0 || n < 0 ) {
            cerr << "BinaryMatrix: negative dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        prepareDim(m,n);

        // Sets all entries to zero.
        memset(this->data,0,this->m * this->N * sizeof(uint32_t) );
    }

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
    bool BinaryMatrix::getAt( int i , int j ) const {

        if ( i < 0 || j < 0 || i >= this->m || j >= this->n ) {
            cerr << "BinaryMatrix::getAt: invalid position." << endl;
            exit(EXIT_FAILURE);
        }


        int j0 , j1;

        j1 = j/32;
        j0 = j%32;

        return ((this->data[i*this->N+j1]) & (((uint32_t)1)<<j0))?true:false;
    }

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
    void BinaryMatrix::setAt( int i , int j ) {

        if ( i < 0 || j < 0 || i >= this->m || j >= this->n ) {
            cerr << "BinaryMatrix::setAt: invalid position." << endl;
            exit(EXIT_FAILURE);
        }

        int j0 , j1;

        j1 = j/32;
        j0 = j%32;

        this->data[i*this->N + j1] |= ((uint32_t)1)<<j0;
    }

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
    void BinaryMatrix::clearAt( int i , int j ) {

        if ( i < 0 || j < 0 || i >= this->m || j >= this->n ) {
            cerr << "BinaryMatrix::clearAt: invalid position." << endl;
            exit(EXIT_FAILURE);
        }

        int j0 , j1;

        j1 = j/32;
        j0 = j%32;

        this->data[i*this->N + j1] &= ~(((uint32_t)1)<<j0);
    }

    /**
     * @brief
     *            Sets all entries of the matrix to 0.
     */
    void BinaryMatrix::setZero() {

        memset(this->data,0,this->m*this->N*sizeof(uint32_t));
    }

    /**
     * @brief
     *            Tests whether this matrix has only zero entries.
     *
     * @return
     *            <code>true</code> if each entry of the matrix is
     *            zero; otherwise, if there exists a non-zero entry,
     *            the function returns <code>false</code>.
     */
    bool BinaryMatrix::isZero() const {

        return MathTools::zeroTest32(this->data,this->m*this->N);
    }

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
    void BinaryMatrix::setIdentity() {

        if ( this->m != this->n ) {
            cerr << "BinaryMatrix::setIdentity: not a square matrix." << endl;
            exit(EXIT_FAILURE);
        }

        setZero();

        for ( int i = 0 ; i < this->n ; i++ ) {
            setAt(i,i);
        }
    }

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
    void BinaryMatrix::exchangeRows( int i0 , int i1 ) {

        if ( i0 < 0 || i1 < 0 || i0 >= this->m || i1 >= this->m ) {
            cerr << "BinaryMatrix::exchangeRows: "
                 << "invalid rows indices." << endl;
            exit(EXIT_FAILURE);
        }

        if ( i0 != i1 ) {

            int N;
            uint32_t *row1 , *row2;
            N = this->N;
            row1 = this->data+i0*N;
            row2 = this->data+i1*N;


            int _N;
            uint64_t *_row1 , *_row2;
            _N = N/2;
            _row1 = (uint64_t*)row1;
            _row2 = (uint64_t*)row2;

            { // Exchange 64 bit blocks
                uint64_t tmp;
                for ( int i = 0 ; i < _N ; i++ , _row1++ , _row2++ ) {
                    tmp = *_row1;
                    *_row1 = *_row2;
                    *_row2 = tmp;
                }
            }

            // Exchange leading 32 bit block (if any)
            if ( N&0x1 ) {
                uint32_t tmp = row1[N-1];
                row1[N-1] = row2[N-1];
                row2[N-1] = tmp;
            }
        }
    }

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
    void BinaryMatrix::exchangeCols( int j0 , int j1 ) {

        if ( j0 < 0 || j1 < 0 || j0 >= this->n || j1 >= this->n ) {
            cerr << "BinaryMatrix::exchangeRows: "
                 << "invalid rows indices." << endl;
            exit(EXIT_FAILURE);
        }

        if ( j0 != j1 ) {

            // I do not know how to exchange the columns inplace
            // different from exchanging each entry separately.
            // Exchanging rows is likely to much more efficient.

            bool c;
            for ( int i = 0 ; i < this->m ; i++ ) {
                c = getAt(i,j0);
                setAt(i,j0,getAt(i,j1));
                setAt(i,j1,c);
            }

        }
    }

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
    void BinaryMatrix::addRow( int i0 , int i1 , int j ) {

        if ( i0 < 0 || i1 < 0 || j < 0 ||
             i0 >= this->m || i1 >= this->m || j >= this->n ) {

            cerr << "BinaryMatrix::addRow: invalid arguments." << endl;
            exit(EXIT_FAILURE);
        }

        int N0 = j/32;

        uint32_t *row0 , *row1;
        row0 = this->data+ i0*this->N + N0;
        row1 = this->data+ i1*this->N + N0;

        MathTools::mxor32(row1,row0,row1,this->N-N0);
    }

    /**
     * @brief
     *            Makes the matrix an upper triangular matrix using
     *            elementary row and column operations.
     *
     * @return
     *            The rank of the matrix.
     */
    int BinaryMatrix::gauss() {

        int n;
        n = min(numRows(),numCols());

        int k;
        for ( k = 0 ; k < n ; k++ ) {

            bool p = false;

            // Attempt to find a pivot element.
            int i0 = k , j0 = k;
            for ( int j = k ; j < numCols() ; j++ ) {
                for ( int i = k ; i < numRows() ; i++ ) {
                    if ( getAt(i,j) ) { // Found a pivot
                        p = true;
                        i0 = i;
                        j0 = j;
                        break;
                    }
                }
                if ( p ) {
                    break;
                }
            }

            if ( !p ) {
                // If no pivot has been been found, we are done.
                break;
            }

            if ( k != j0 ) {
                exchangeCols(k,j0);
            }

            if ( k != i0 ) {
                exchangeRows(k,i0);
            }

            // Now, the element at (k,k) is non-zero and
            // we can eliminate the elements in further rows of the
            // 'k'th column.
            for ( int i = k+1 ; i < numRows() ; i++ ) {
                if ( getAt(i,k) ) {
                    addRow(k,i,k);
                }
            }
        }

        return k;
    }

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
    int BinaryMatrix::rank() const {

    	BinaryMatrix tmp(*this);

    	return tmp.gauss();
    }

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
    void BinaryMatrix::random( bool tryRandom ) {

        int m , N;

        m = this->m;
        N = this->N;

        for ( int i = 0 ; i < m ; i++ ) {
            for ( int j = 0 ; j < N ; j++ ) {
                this->data[i*N+j] = MathTools::rand32(tryRandom);
            }
        }

        if ( this->n % 32 > 0 ) {

            uint32_t lb;
            lb = 1;
            lb <<= (this->n%32);
            lb--;

            for ( int i = 0 ; i < m ; i++ ) {
                this->data[i*N+N-1] &= lb;
            }
        }
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
    void BinaryMatrix::swap( BinaryMatrix & a , BinaryMatrix & b ) {

        uint32_t *data;
        int m , n , N;

        data = a.data;
        m    = a.m;
        n    = a.n;
        N    = a.N;

        a.data = b.data;
        a.m    = b.m;
        a.n    = b.n;
        a.N    = b.N;

        b.data = data;
        b.m    = m;
        b.n    = n;
        b.N    = N;
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
    void BinaryMatrix::add
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        if ( a.m != b.m || a.n != b.n ) {
            cerr << "BinaryMatrix::add: incompatible dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        c.prepareDim(a.m,a.n);

        MathTools::mxor32(c.data,a.data,b.data, c.m*c.N );
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
    void BinaryMatrix::transpose
    ( BinaryMatrix & b , const BinaryMatrix & a ) {

        if ( &b == &a ) {
            BinaryMatrix tb;
            transpose(tb,a);
            swap(b,tb);
            return;
        }

        int m , n;
        m = a.numRows();
        n = a.numCols();

        b.setDimension(n,m);

        for ( int i = 0 ; i < m ; i++ ) {
            for ( int j = 0 ; j < n ; j++ ) {
                b.setAt(j,i,a.getAt(i,j));
            }
        }

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
    void BinaryMatrix::mul_t
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        if ( &c == &a || &c == &b ) {
            BinaryMatrix tc;
            mul_t(tc,a,b);
            swap(c,tc);
            return;
        }

        if ( a.numCols() != b.numCols() ) {
            cerr << "BinaryMatrix::mul_t: incompatible dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        int m , n , N;
        m = a.numRows();
        n = b.numRows();
        N = a.N;

        c.setDimension(m,n);

        uint32_t *A , *B;
        for ( int i = 0 ; i < m ; i++ ) {
            A = a.data + i * N;
            for ( int j = 0 ; j < n ; j++ ) {
                B = b.data + j * N;
                if ( (MathTools::sprod32(A,B,N) & 0x1) ) {
                    c.setAt(i,j);
                }
            }
        }
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
    void BinaryMatrix::mul
    ( BinaryMatrix & c ,
      const BinaryMatrix & a , const BinaryMatrix & b ) {

        BinaryMatrix tb;
        transpose(tb,b);
        mul_t(c,a,tb);
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
    void BinaryMatrix::lnot
    ( BinaryMatrix & B , const BinaryMatrix & A ) {

    	if ( &B == &A ) {
    		BinaryMatrix tB;
    		lnot(tB,A);
    		swap(B,tB);
    		return;
    	}

    	int m , n;
    	m = A.numRows();
    	n = A.numCols();

    	B.setDimension(m,n);
    	for ( int i = 0 ; i < m ; i++ ) {
    		for ( int j = 0 ; j < n ; j++ ) {
    			if ( !A.getAt(i,j) ) {
    				B.setAt(i,j);
    			}
    		}
    	}
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
    void BinaryMatrix::submatrix
    ( BinaryMatrix & A0 , const BinaryMatrix & A ,
      int i0 , int j0 , int m0 , int n0 ) {

        if ( &A0 == &A ) {
            BinaryMatrix tA(A);
            submatrix(A0,tA,i0,j0,m0,n0);
            return;
        }

        if ( i0 < 0 || j0 < 0 || m0 < 0 || n0 < 0 ||
             i0+m0 > A.m || j0+n0 > A.n ) {
            cerr << "BinaryMatrix::submatrix: invalid arguments." << endl;
            exit(EXIT_FAILURE);
        }

        A0.setDimension(m0,n0);
        for ( int i = 0 ; i < m0 ; i++ ) {
            for ( int j = 0 ; j < n0 ; j++ ) {
                if ( A.getAt(i0+i,j0+j) ) {
                    A0.setAt(i,j);
                }
            }
        }
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
    void BinaryMatrix::concatCols
    ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B ) {

        if ( &C == &A || &C == &B ) {
            BinaryMatrix tC;
            concatCols(tC,A,B);
            swap(C,tC);
            return;
        }

        if ( A.m != B.m ) {
            cerr << "BinaryMatrix::concatCols: "
                 << "matrices must have the same number of rows." << endl;
            exit(EXIT_FAILURE);
        }

        int m , n0 , n1;
        m = A.m;
        n0 = A.n;
        n1 = B.n;

        C.setDimension(m,n0+n1);

        // The first matrix can be copied efficiently but ...
        for ( int i = 0 ; i < m ; i++ ) {
            memcpy(C.data+i*C.N,A.data+i*A.N,A.N*sizeof(uint32_t));
        }

        // ... the second matrix's entries are copied separately.
        for ( int i = 0 ; i < m ; i++ ) {
            for ( int j = 0 ; j < n1 ; j++ ) {
                if ( B.getAt(i,j) ) {
                    C.setAt(i,n0+j);
                }
            }
        }
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
    void BinaryMatrix::concatRows
    ( BinaryMatrix & C , const BinaryMatrix & A , const BinaryMatrix & B ) {

        if ( &C == &A || &C == &B ) {
            BinaryMatrix tC;
            concatRows(tC,A,B);
            swap(C,tC);
            return;
        }

        if ( A.n != B.n ) {
            cerr << "BinaryMatrix::concatRows: "
                 << "matrices must have the same number of columns." << endl;
            exit(EXIT_FAILURE);
        }

        int n , m0 , m1 , N;
        n = A.n;
        m0 = A.m;
        m1 = B.m;
        N = A.N;

        C.prepareDim(m0+m1,n);

        for ( int i = 0 ; i < m0 ; i++ ) {
            memcpy(C.data+i*N,A.data+i*N,N*sizeof(uint32_t));
        }

        for ( int i = 0 ; i < m1 ; i++ ) {
            memcpy(C.data+(m0+i)*N,B.data+i*N,N*sizeof(uint32_t));
        }
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
    ostream & operator<<
    ( ostream & out , const BinaryMatrix & mat ) {

        out << "[";

        for ( int i = 0 ; i < mat.numRows() ; i++ ) {

            if ( i > 0 ) {
                out << "," << endl << " ";
            }

            out << "[";

            for ( int j = 0 ; j < mat.numCols() ; j++ ) {
                if ( mat.getAt(i,j) ) {
                    out << "1";
                } else {
                    out << "0";
                }
                if (j+1<mat.numCols() ) {
                    out << ",";
                }
            }

            out << "]";

        }

        out << "]";

        return out;
    }

}
