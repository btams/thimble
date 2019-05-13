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
 * @file BinaryVector.cpp
 *
 * @brief
 *            Implementation of a class for the represention of and
 *            computation with binary vectors of fixed but arbitrary
 *            degree. The class is provided by the 'BinaryVector.h'
 *            header.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <thimble/math/MathTools.h>
#include <thimble/math/linalg/BinaryVector.h>

using namespace std;

/**
 * @brief The library's namespace.
 */
namespace thimble {

    /**
     * @brief
     *            Standard constructor.
     *
     * @details
     *            Constructs a vector of specified length.
     *
     * @param n
     *            The length of this vector.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    BinaryVector::BinaryVector( int n ) {

        this->data = NULL;
        this->length = 0;
        this->numWords = 0;

        setLength(n);
    }


    /**
     * @brief
     *           Copy constructor.
     *
     * @param v
     *           The binary vector of which a copy is constructed.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    BinaryVector::BinaryVector( const BinaryVector & v ) {
        this->data = NULL;
        this->length = 0;
        this->numWords = 0;
        assign(v);
    }



    /**
     * @brief
     *            Destructor.
     *
     * @details
     *            Frees the data needed to hold the vector's entries.
     */
    BinaryVector::~BinaryVector() {
        free(this->data);
    }



    /**
     * @brief
     *            Assignment operator (procedural version).
     *
     * @param v
     *            The binary vector assigned to this vector.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::assign( const BinaryVector & v ) {

        if ( this != &v ) {

            int n = v.length/32+(v.length%32?1:0);

            reserve(v.length);
            memcpy(this->data,v.data,n*sizeof(uint32_t));
            this->length = v.length;
        }
    }

    /**
     * @brief
     *            Compares two binary vectors on equality.
     *
     * @param v
     *            The vector with which this vector is compared.
     *
     * @return
     *            <code>true</code> if this vector equals <i>v</i>;
     *            otherwise, <code>false</code>.
     *
     * @attention
     *            Note that if two vectors are of different length,
     *            they are considered to be different causing
     *            the function to return <code>false</code>.
     */
    bool BinaryVector::operator==( const BinaryVector & v ) const {

    	if ( this->length != v.length ) {
    		return false;
    	}

    	if ( this != &v ) {

			int N = this->length/32 + (this->length%32?1:0);

			for ( int j = 0 ; j < N ; j++ ) {
				if ( this->data[j] != v.data[j] ) {
					return false;
				}
			}
    	}

    	return true;
    }

    /**
     * @brief
     *            Ensures that this instance is able to represent
     *            a vector of specified length without requiring
     *            reallocation.
     *
     * @param length
     *            The reserved length for this vector.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::reserve( int length ) {

        if ( length > 0 ) {

            int numWords = length/32+(length%32?1:0);

            if ( numWords > this->numWords ) {

                if ( this->numWords > 0 ) {
                    this->data = (uint32_t*)realloc
                            (this->data,numWords*sizeof(uint32_t));
                } else {
                    this->data = (uint32_t*)malloc
                            (numWords*sizeof(uint32_t));
                }

                if ( this->data == NULL ) {
                    cerr << "BinaryVector::ensureCapacity: "
                         << "out of memory." << endl;
                    exit(EXIT_FAILURE);
                }

                this->numWords = numWords;
            }
        }
    }

    /**
     * @brief
     *            Sets the dimension of the vector to the specified
     *            length.
     *
     * @details
     *            If the length of the vector increases, the entries
     *            of the vector up to the old length remain unchanged
     *            and the new entries are padded with zeros. If the
     *            length of the vector is decreased, the first
     *            <code>length</code> entries are kept.
     *
     * @param newLength
     *            The new length of the vector.
     *
     * @warning
     *            If length is negative, an error message is printed
     *            to <code>stderr</code> and the program exits with
     *            status 'EXIT_STATUS'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::setLength( int newLength ) {

        if ( newLength < 0 ) {
            cerr << "BinaryVector::setLength: must be positive." << endl;
            exit(EXIT_FAILURE);
        }

        // Reserve
        reserve(newLength);

        int oldLength , oldNumWords , newNumWords;

        oldLength = this->length;
        oldNumWords = oldLength/32+(oldLength%32?1:0);
        newNumWords = newLength/32+(newLength%32?1:0);

        // Initialize the new words (if any) as zeros.
        if ( newNumWords > oldNumWords ) {
            memset(this->data+oldNumWords,0,(newNumWords-oldNumWords)*sizeof(uint32_t));
        }

        // Pad the last bits of the latest word with zeros
        if ( newLength % 32 != 0 ) {
            uint32_t lb;
            lb = 1;
            lb <<= newLength%32;
            lb--;
            this->data[newNumWords-1] &= lb;
        }

        this->length = newLength;
    }


    /**
     * @brief
     *            Access the value at the specified position of
     *            the vector.
     *
     * @param j
     *            The position ranging from
     *            0,...,\link getLength()\endlink-1.
     *
     * @return
     *            <code>true</code> if the <i>j</i>th position
     *            equals 1 and, otherwise, if the <i>j</i>th
     *            position equals 0, <code>false</code>.
     *
     * @warning
     *            If <i>j</i> is negative or if <i>j</i> is greater
     *            than or equals \link getLength()\endlink, an
     *            error message is printed to <code>stderr</code>
     *            and the program exits with status 'EXIT_FAILURE'
     */
    bool BinaryVector::getAt( int j ) const {

        if ( j < 0 || j >= this->length ) {
            cerr << "BinaryVector::getAt: index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        int j0 , j1;

        j0 = j/32;
        j1 = j%32;

        return (this->data[j0]&((uint32_t)1<<j1))?true:false;
    }


    /**
     * @brief
     *            Ensures that the value at the specified position
     *            equals 1.
     *
     * @param j
     *            The position ranging from
     *            0,...,\link getLength()\endlink-1.
     *
     * @warning
     *            If <i>j</i> is negative or if <i>j</i> is greater
     *            than or equals \link getLength()\endlink, an
     *            error message is printed to <code>stderr</code>
     *            and the program exits with status 'EXIT_FAILURE'
     */
    void BinaryVector::setAt( int j ) {

        if ( j < 0 || j >= this->length ) {
            cerr << "BinaryVector::setAt: index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        int j0 , j1;

        j0 = j/32;
        j1 = j%32;

        this->data[j0] |= (((uint32_t)1)<<j1);
    }

    /**
     * @brief
     *            Tests whether the binary vector has only entries equals
     *            0.
     *
     * @return
     *            <code>true</code> if all entries of the vector are 0;
     *            otherwise, if there exists a non-zero entry, the
     *            function returns <code>false</code>.
     */
    bool BinaryVector::isZero() const {
    	int s = this->length/32+(this->length%32?1:0);
    	return MathTools::zeroTest32(this->data,s);
    }

    /**
     * @brief
     *            Initialize all entries of this vector with 0.
     */
    void BinaryVector::setZero() {
        int n = this->length / 32 + (this->length%32?1:0);
        memset(this->data,0,n*sizeof(uint32_t));
    }

    /**
     * @brief
     *           Exchanges the entries of the vector at the two
     *           specified positions.
     *
     * @param j0
     *           First position.
     *
     * @param j1
     *           Second position.
     *
     * @warning
     *           If <i>j0</i> or <i>j1</i> is smaller than zero or
     *           greater than or equals the
     *           \link getLength() length\endlink of the vector, an
     *           error message is printed to <code>stderr</code> and
     *           the program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::exchange( int j0 , int j1 ) {

    	int n = getLength();

    	if ( j0 < 0 || j1 < 0 || j0 >= n || j1 >= n ) {
    		cerr << "BinaryVector::exchange: invalid arguments." << endl;
    		exit(EXIT_FAILURE);
    	}

    	if ( j0 != j1 ) {
    		bool c0 , c1;
    		c0 = getAt(j0);
    		c1 = getAt(j1);
    		setAt(j0,c1);
    		setAt(j1,c0);
    	}
    }

    /**
     * @brief
     *            Ensures that the value at the specified position
     *            equals 0.
     *
     * @param j
     *            The position ranging from
     *            0,...,\link getLength()\endlink-1.
     *
     * @warning
     *            If <i>j</i> is negative or if <i>j</i> is greater
     *            than or equals \link getLength()\endlink, an
     *            error message is printed to <code>stderr</code>
     *            and the program exits with status 'EXIT_FAILURE'
     */
    void BinaryVector::clearAt( int j ) {

        if ( j < 0 || j >= this->length ) {
            cerr << "BinaryVector::clearAt: index out of bounds." << endl;
            exit(EXIT_FAILURE);
        }

        int j0 , j1;

        j0 = j/32;
        j1 = j%32;

        this->data[j0] &= ~(((uint32_t)1)<<j1);
    }


    /**
     * @brief
     *            Set the value of the vector at a specified
     *            position to a specified value.
     *
     * @param j
     *            The position ranging from
     *            0,...,\link getLength()\endlink-1.
     *
     * @param c
     *            <code>true</code> if the value at <i>j</i> is set
     *            to 1 and <code>false</code> otherwise.
     *
     * @warning
     *            If <i>j</i> is negative or if <i>j</i> is greater
     *            than or equals \link getLength()\endlink, an
     *            error message is printed to <code>stderr</code>
     *            and the program exits with status 'EXIT_FAILURE'
     */
    void BinaryVector::setAt( int j , bool c ) {

        if ( c ) {
            setAt(j);
        } else {
            clearAt(j);
        }
    }


    /**
     * @brief
     *            Initialize the entries of this vector with random
     *            values.
     *
     * @details
     *            Each \link getLength()\endlink entry of the vector is
     *            overwritten with random values.
     *
     * @param tryRandom
     *           If <code>true</code>, the method uses a cryptographic
     *           number generator if available on the system; otherwise,
     *           the method wraps around the standard <code>rand()</code>
     *           function.
     */
    void BinaryVector::random( bool tryRandom ) {

        int n = this->length/32+(this->length%32?1:0);

        for ( int j = 0 ; j < n ; j++ ) {
            this->data[j] = MathTools::rand32(tryRandom);
        }

        // Set the most significant bits of the last word to
        // zero.
        if ( this->length % 32 != 0 ) {
            uint32_t lb;
            lb = 1;
            lb <<= this->length%32;
            lb--;
            this->data[n-1] &= lb;
        }
    }

    /**
     * @brief
     *           Replaces the entries of this vector such that it becomes
     *           random with the specified Hamming weight.
     *
     * @param hammingWeight
     *           The Hamming weight of this vector after the method has
     *           finished.
     *
     * @param tryRandom
     *           If <code>true</code>, the method is advised to use a
     *           cryptographic number generator; otherwise, the method
     *           wraps around the standard <code>rand()</code> function.
     *
     * @warning
     *           If <code>hammingWeight</code> is negative or greater than
     *           or equals \link getLength()\endlink, an error message is
     *           printed to <code>stderr</code> and the program exits with
     *           status 'EXIT_FAILURE'.
     *
     * @warning
     *           If not enough memory could be provided, an error
     *           message is printed to <code>stderr</code> and the
     *           program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::wrandom( int hammingWeight , bool tryRandom ) {

        if ( hammingWeight < 0 || hammingWeight > this->length ) {
            cerr << "BinaryVector::wrandom: invalid Hamming weight." << endl;
            exit(EXIT_FAILURE);
        }

        setZero();

        int n = this->length;

        int *indices = (int*)malloc( n * sizeof(int) );
        if ( indices == NULL ) {
            cerr << "BinaryVector::wrandom: out of memory." << endl;
            exit(EXIT_FAILURE);
        }
        for ( int j = 0 ; j < n ; j++ ) {
            indices[j] = j;
        }

        int hw = n;

        while ( hw > hammingWeight ) {
            // which indices to remove?
            int i = (int)(MathTools::rand64(tryRandom) % (uint64_t)hw);
            for ( int j = i ; j+1 < hw ; j++ ) {
                indices[j] = indices[j+1];
            }
            --hw;
        }

        // Now 'hw==hammingWeight'.

        // The first 'hw' positions in 'indices' of the vector
        // are set to '1'.
        for ( int j = hw-1 ; j >= 0 ; j-- ) {
            setAt(indices[j]);
        }

        free(indices);
    }

    /**
     * @brief
     *           Returns the Hamming weight of the vector.
     *
     * @details
     *           Write for the vector
     *           \f[
     *            v = (v_0,...,v_{n-1})
     *           \f]
     *           where \f$v_j\in\{0,1\}\f$. The <i>Hamming weight</i>
     *           of the vector <i>v</i> is defined as the number of
     *           positions <i>j</i> where \f$v_j\neq 0\f$ which, for
     *           binary vectors, is equivalent to the number of
     *           positions <i>j</i> where \f$v_j=1\f$.
     *
     * @return
     *           The Hamming weight of the vector.
     */
    int BinaryVector::hammingWeight() const {

    	int hw;
    	int s;
    	uint32_t *data;
    	int j;

    	hw = 0;
    	s = (this->length/32) + (this->length%32?1:0);
    	data = this->data;

    	for ( j = 0 ; j < s ; j++ , data++ ) {
    		hw += MathTools::hammingWeight(*data);
    	}

    	return hw;
    }

    /**
     * @brief
     *           Initializes the first <i>k</i> of the vector as 1s and
     *           the remaining as zeros.
     *
     * @details
     *           Using the method \link bnext()\endlink, we can systematically
     *           iterate through all vectors of
     *           \link hammingWeight() Hamming weight\endlink exactly equals
     *           <i>k</i>. This method sets this vectors to the first vector
     *           of this sequence, i.e.,
     *           \f[
     *            v=(\underbrace{1,...,1}_{k~\mbox{times}}~,~
     *               \underbrace{0,...,0}_{n-k~\mbox{times}}).
     *           \f]
     *
     * @param k
     *           The number of positions of the vector being set to 1.
     *
     * @see bnext()
     *
     * @warning
     *           If <i>k</i> is negative or greater than
     *           \link getLength()\endlink, an error message is printed
     *           to <code>stderr</code> and the program exits with status
     *           'EXIT_FAILIRE'.
     */
    void BinaryVector::binit( int k ) {

    	int n = getLength();

    	if ( k < 0 || k > n ) {
    		cerr << "BinaryVector::binit: invalid argument." << endl;
    		exit(EXIT_FAILURE);
    	}

    	setZero();

    	for ( int j = 0 ; j < k ; j++ ) {
    		setAt(j);
    	}
    }

    /**
     * @brief
     *           Iterates to the ''next'' binary vector of the
     *           same Hamming weight.
     *
     * @details
     *           Assume the length of the vector is <i>n</i> and its
     *           Hamming weight is <i>k</i>. The function selects a
     *           ''next'' vector in a sequence of all \f${k\choose n}\f$
     *           combinations of <i>n</i>-length vectors with Hamming
     *           weight <i>k</i> starting from
     *           \f[
     * (\underbrace{1,...,1}_{k~\mbox{times}}~,~
     *  \underbrace{0,...,0}_{n-k~\mbox{times}})
     *           \f]
     *           to
     *           \f[
     * (\underbrace{0,...,0}_{n-k~\mbox{times}}~,~
     *  \underbrace{1,...,1}_{k~\mbox{times}}).
     *           \f]
     *
     * @return
     *           <code>true</code> if the content of the vector has been
     *           successfully changed to the ''next'' vector of the
     *           same Hamming weight; otherwise, if the vector is already
     *           the ''last'' vector, the result is
     *           <code>false</code>
     *
     * @see binit()
     * @see next()
     */
    bool BinaryVector::bnext() {

    	int n , i;

    	n = getLength();

    	i = 0;

    	for ( int j = 0 ; j < n ; j++ ) {

    		if ( getAt(j) ) {

    			if ( j+1 >= n ) {
    				return false;
    			}

    			clearAt(j);

    			if ( getAt(j+1) ) {
    				setAt(i);
    				++i;
    			} else {
    				setAt(j+1);
    				return true;
    			}

    		}

    	}

    	return false;
    }

    /**
     * @brief
     *           Iterates to the ''next'' vector in a sequence of distinct
     *           vectors of fixed length.
     *
     * @details
     *           The function can be used to systematically iterate through
     *           all distinct \f$2^n\f$ vectors starting from
     *           \f$(0,...,0)\f$ to \f$(1,...,0)\f$. The sequence of the
     *           iteration has the property that its vectors are
     *           monotonically increasing in their
     *           \link hammingWeight() Hamming weight\endlink which is,
     *           in particular, useful for solving problems related to
     *           <i>maximum likelihood decoding</i>.
     *
     * @return
     *           <code>true</code> if the content has been successfully
     *           changed to the ''next'' vector; otherwise, if the vector
     *           is already the ''last'' vector, the function returns
     *           <code>false</code>.
     *
     * @see bnext()
     */
    bool BinaryVector::next() {

    	if ( bnext() ) {
    		return true;
    	}

    	int k = hammingWeight()+1;

    	if ( k > getLength() ) {
    		return false;
    	}

    	binit(k);

    	return true;
    }

    /**
     * @brief
     *           Swaps the content of two binary vectors such the
     *           the one represents the other.
     *
     * @param v
     *           Binary vector being assigned the fields of <i>w</i>.
     *
     * @param w
     *           Binary vector being assigned the fields of <i>v</i>.
     */
    void BinaryVector::swap( BinaryVector & v , BinaryVector & w ) {

        uint32_t *data;
        int length , numWords;

        data = v.data;
        length = v.length;
        numWords = v.numWords;

        v.data = w.data;
        v.length = w.length;
        v.numWords = w.numWords;

        w.data = data;
        w.length = length;
        w.numWords = numWords;
    }

    /**
     * @brief
     *           Performs an addition of two binary vectors.
     *
     * @param c
     *           On output, the sum of the input
     *           <code>a</code> and <code>b</code>.
     *
     * @param a
     *           First summand vector.
     *
     * @param b
     *           Second summand vector.
     *
     * @warning
     *           If the length of <code>a</code> and <code>b</code>
     *           are different, an error message is printed to
     *           <code>stderr</code> and the program exits with status
     *           'EXIT_FAILURE'.
     *
     * @warning
     *            If not enough memory could be provided, an error
     *            message is printed to <code>stderr</code> and the
     *            program exits with status 'EXIT_FAILURE'.
     */
    void BinaryVector::add
    ( BinaryVector & c , const BinaryVector & a , const BinaryVector & b ) {

        int n = a.length;

        if ( n != b.length ) {
            cerr << "BinaryVector::add: incompatible vector lengths." << endl;
            exit(EXIT_FAILURE);
        }

        c.reserve(n);

        MathTools::mxor32(c.data,a.data,b.data,n/32+(n%32?1:0));

        c.length = n;
    }

    /**
     * @brief
     *            Determines the Hamming distance between two binary
     *            vectors.
     *
     * @param a
     *            First binary vector.
     *
     * @param b
     *            Second binary vector.
     *
     * @return
     *            The Hamming distance between <code>a</code> and
     *            <code>b</code>.
     *
     * @warning
     *            If <code>a</code> and <code>b</code> have a different
     *            length, an error message is printed to
     *            <code>stderr</code> and the program exits with status
     *            'EXIT_FAILURE'.
     */
    int BinaryVector::hammingDistance
    ( const BinaryVector & a , const BinaryVector & b ) {

        if ( a.length != b.length ) {
            cerr << "BinaryVector::hammingDistance: "
                 << "vectors have different length." << endl;
            exit(EXIT_FAILURE);
        }

        return MathTools::hd32(a.data,b.data,a.length/32+(a.length%32?1:0));
    }

    /**
     * @brief
     *            Prints a representation of a vector to the specified
     *            output stream.
     *
     * @param out
     *            The output stream.
     *
     * @param v
     *            The binary vector being printed to <code>out</code>.
     *
     * @return
     *            A reference to <code>out</code>.
     */
    ostream & operator<<( ostream & out , const BinaryVector & v ) {

        out << "[";

        for ( int j = 0 ; j < v.getLength() ; j++ ) {
            if ( v.getAt(j) ) {
                out << "1";
            } else {
                out << "0";
            }
        }

        out << "]";

        return out;
    }

}
