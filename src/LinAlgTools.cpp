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
 * @file LinAlgTools.cpp
 *
 * @brief
 *            Implementation of a class and utility functions/methods that may
 *            be useful for purposes related to linear algebra.
 *
 * @author Benjamin Tams
 */

#include <cstring>
#include <iostream>

#include <thimble/math/MathTools.h>
#include <thimble/math/linalg/BinaryVector.h>
#include <thimble/math/linalg/BinaryMatrix.h>
#include <thimble/math/linalg/LinAlgTools.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

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
     *           where the vectors are intepreted as column vectors.
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
    void LinAlgTools::mul
    (BinaryVector & y , const BinaryMatrix & A , const BinaryVector & x ) {

        if ( &y == &x ) {
            BinaryVector tx(x);
            mul(y,A,tx);
            return;
        }

        if ( x.getLength() != A.numCols() ) {
            cerr << "LinAlgTools::mul: "
                 << "invalid dimension in matrix-vector multiplication."
                 << endl;
            exit(EXIT_FAILURE);
        }

        int m , N;
        m = A.numRows();
        N = A.numCols()/32+(A.numCols()%32?1:0);

        y.setLength(m);
        y.setZero();

        const uint32_t *_A , *_x;
        _A = A.getData();
        _x = x.getData();

        // Runs 'm' scalar products
        for ( int i = 0 ; i < m ; i++ , _A += N ) {
            if ( (MathTools::sprod32(_A,_x,N)&0x1) ) {
                y.setAt(i);
            }
        }
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
    void LinAlgTools::mul
    ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A ) {

        if ( &B == &A ) {
            BinaryMatrix tA(A);
            mul(B,P,tA);
            return;
        }

        if ( P.getDimension() != A.numRows() ) {
            cerr << "LinAlgTools::mul: incompatible dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        int m , n , N;
        m = A.numRows();
        n = A.numCols();
        N = n/32+(n%32?1:0);

        B.setDimension(m,n);

        const uint32_t *in = A.getData();
        uint32_t *out = B.getData_nonconst();

        // Copy the rows of 'A' to 'B' but with reordered index as
        // specified by 'P'
        for ( int i = 0 ; i < m ; i++ ) {

            int j = P.eval(i);

            memcpy(out+j*N,in+i*N, N * sizeof(uint32_t));
        }
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
    void LinAlgTools::imul
    ( BinaryMatrix & B , const Permutation & P , const BinaryMatrix & A ) {

        if ( &B == &A ) {
            BinaryMatrix tA(A);
            mul(B,P,tA);
            return;
        }

        if ( P.getDimension() != A.numRows() ) {
            cerr << "LinAlgTools::imul: incompatible dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        int m , n , N;
        m = A.numRows();
        n = A.numCols();
        N = n/32+(n%32?1:0);

        B.setDimension(m,n);

        const uint32_t *in = A.getData();
        uint32_t *out = B.getData_nonconst();

        // Copy the rows of 'A' to 'B' but with reordered index as
        // specified by the inverse of 'P'
        for ( int i = 0 ; i < m ; i++ ) {

            int j = P.eval(i);

            memcpy(out+i*N,in+j*N, N * sizeof(uint32_t));
        }
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
    void LinAlgTools::mul
    ( BinaryVector & y , const Permutation & P ,const BinaryVector & x ) {

        if ( &y == &x ) {
            BinaryVector tx(x);
            mul(y,P,tx);
            return;
        }

        int n = P.getDimension();

        if ( n != x.getLength() ) {
            cerr << "LinAlgTools::mul: incompatible dimensions." << endl;
            exit(EXIT_FAILURE);
        }

        y.setLength(n);
        y.setZero();

        for ( int i = 0 ; i < n ; i++ ) {
            if ( x.getAt(i) ) {
                int j = P.eval(i);
                y.setAt(j);
            }
        }
    }

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
    void LinAlgTools::image( BinaryMatrix & B , const BinaryMatrix & A ) {

        // Implementation follows Algorithm 2.3.2 in
        // H. Cohen (1996). A Course in Computational Algebraic Number Theory.
        // Springer-Verlag

        BinaryMatrix M(A);

        int m , n;
        m = A.numRows();
        n = A.numCols();

        // Step 1. [Initialize]

        int r = 0;

        int *c = (int*)malloc(m*sizeof(int));
        if ( c == NULL ) {
            cerr << "LinAlgTools::image: out of memory." << endl;
            exit(EXIT_FAILURE);
        }
        for ( int i = 0 ; i < m ; i++ ) {
            c[i] = -1;
        }

        for (int k = 0 ; k < n ; k++ ) {

            // Step 2. [Scan column]
            int j = -1;
            for ( int i = 0 ; i < m ; i++ ) {
                if ( c[i] < 0 && M.getAt(i,k) ) {
                    j = i;
                    break;
                }
            }

            if ( j < 0 ) {
                ++r;
            } else {

                // Step 3. [Eliminate]
                for ( int i = 0 ; i < m ; i++ ) {
                    if ( i != j ) {
                        if ( M.getAt(i,k) ) {
                            M.addRow(j,i);
                        }
                    }
                }

                c[j] = k;
            }
        }

        // Step 5. [Output image]
        B.setDimension(m,n-r);
        for ( int j = 0 , l = 0 ; j < m ; j++ ) {
            if ( c[j] >= 0 ) {
                for ( int i = 0 ; i < m ; i++ ) {
                    if ( A.getAt(i,c[j]) ) {
                        B.setAt(i,l);
                    }
                }
                ++l;
            }
        }

        free(c);
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
    void LinAlgTools::kernel( BinaryMatrix & B , const BinaryMatrix & A ) {

        // Make a copy of 'A' in which Gauss elimination is performed inplace
        BinaryMatrix C(A);

        // Permutation applied on the columns of 'A' to bring it in upper
        // triangular form.
        Permutation P(C.numCols());

        int n;
        n = min(C.numRows(),C.numCols());

        int k;
        for ( k = 0 ; k < n ; k++ ) {

            bool p = false;

            // Attempt to find a pivot element.
            int i0 = k , j0 = k;
            for ( int j = k ; j < C.numCols() ; j++ ) {
                for ( int i = k ; i < C.numRows() ; i++ ) {
                    if ( C.getAt(i,j) ) { // Found a pivot
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
                C.exchangeCols(k,j0);
                P.exchange(k,j0);
            }

            if ( k != i0 ) {
                C.exchangeRows(k,i0);
            }

            // Now, the element at (k,k) is non-zero and
            // we can eliminate the elements in further rows of the
            // 'k'th column.
            for ( int i = k+1 ; i < C.numRows() ; i++ ) {
                if ( C.getAt(i,k) ) {
                    C.addRow(k,i,k);
                }
            }
        }

        BinaryMatrix PB;
        PB.setDimension(C.numCols(),C.numCols()-k);

        for ( int l = 0 ; l < PB.numCols() ; l++ ) {

            PB.setAt(l+k,l);

            for ( int i = k-1 ; i >= 0 ; i-- ) {

                bool c = 0;
                for ( int j = i+1 ; j <= l+k ; j++ ) {
                    if ( PB.getAt(j,l) && C.getAt(i,j) ) {
                        c ^= 1;
                    }
                }

                if ( c ) {
                    PB.setAt(i,l);
                }
            }
        }

        mul(B,P,PB);

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
    bool LinAlgTools::solve
    ( BinaryVector & x , const BinaryMatrix & A ,
      const BinaryVector & y ) {

    	if ( y.getLength() != A.numRows() ) {
    		cerr << "LinAlgTools::solve: invalid dimensions." << endl;
    		exit(EXIT_FAILURE);
    	}

        // Make a copy of 'A' in which Gauss elimination is performed inplace
        BinaryMatrix C(A);

        // Make a copy of 'y' apply exchanges to its rows on performing
        // Gaussian elimination.
        BinaryVector z(y);

        // Permutation applied on the columns of 'A' to bring it in upper
        // triangular form.
        Permutation P(C.numCols());

        int n;
        n = min(C.numRows(),C.numCols());

        int k;
        for ( k = 0 ; k < n ; k++ ) {

            bool p = false;

            // Attempt to find a pivot element.
            int i0 = k , j0 = k;
            for ( int j = k ; j < C.numCols() ; j++ ) {
                for ( int i = k ; i < C.numRows() ; i++ ) {
                    if ( C.getAt(i,j) ) { // Found a pivot
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
                C.exchangeCols(k,j0);
                P.exchange(k,j0);
            }

            if ( k != i0 ) {
            	z.exchange(k,i0);
                C.exchangeRows(k,i0);
            }

            // Now, the element at (k,k) is non-zero and
            // we can eliminate the elements in further rows of the
            // 'k'th column.
            bool zk = z.getAt(k);
            for ( int i = k+1 ; i < C.numRows() ; i++ ) {
                if ( C.getAt(i,k) ) {
                	z.setAt(i,zk^z.getAt(i));
                    C.addRow(k,i,k);
                }
            }
        }

        // Now 'C' is an upper triangular matrix and 'k' the rank of 'C'
        // and thus of 'A'.

        // Check if the system can be solved.
        for ( int j = k ; j < z.getLength() ; j++ ) {
        	// If any of the entries of 'z' at position higher
        	// than the rank of 'C', is set, 'z' is not an element
        	// of the image of 'C' and thus the system cannot be solved.
        	if ( z.getAt(j) ) {
        		return false;
        	}
        }

        int N = C.numCols()/32+(C.numCols()%32?1:0);

        x.setLength(C.numCols());

        for ( int j = k - 1 ; j >= 0 ; j-- ) {

        	bool s = (MathTools::sprod32
        			  (C.getData()+j*N,x.getData(),N)&0x1)?true:false;

        	bool zj = z.getAt(j);

        	if ( (s && !zj) || (!s && zj) ) {
        		x.setAt(j);
        	}
        }

        mul(x,P,x);

        return true;
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
    bool LinAlgTools::isInKernel
    ( const BinaryVector & v , const BinaryMatrix & A ) {

        if ( v.getLength() != A.numCols() ) {
            cerr << "LinAlgTools::isInKernel: "
                 << "invalid dimension."
                 << endl;
            exit(EXIT_FAILURE);
        }

        int m , N;
        m = A.numRows();
        N = A.numCols()/32+(A.numCols()%32?1:0);

        const uint32_t *_A , *_v;
        _A = A.getData();
        _v = v.getData();

        // Runs 'm' scalar products
        for ( int i = 0 ; i < m ; i++ , _A += N ) {
            if ( (MathTools::sprod32(_A,_v,N)&0x1) ) {
                return false;
            }
        }

        return true;
    }
}
