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
 * @file GuruswamiSudanDecoder.cpp
 *
 * @brief
 *            Implements the functionalities provided by
 *            'GuruswamiSudanDecoder.h' which is to enable a Guruswami-Sudan
 *            list decoding algorithm for decoding Reed-Solomon codes in
 *            original view.
 *
 * @author Benjamin Tams
 */

#include <stdint.h>
#include <climits>
#include <vector>
#include <algorithm>
#include <iostream>

#include <thimble/math/MathTools.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/math/numbertheory/SmallBinaryFieldBivariatePolynomial.h>
#include <thimble/ecc/GuruswamiSudanDecoder.h>

using namespace std;

namespace thimble {

	class DecodedPolynomialComparator {

	private:
		const uint32_t *x;
		const uint32_t *y;
		int n;

	public:
		inline DecodedPolynomialComparator
		( const uint32_t *x , const uint32_t *y , int n ) {

			this->x = x;
			this->y = y;
			this->n = n;
		}

		bool operator()
		( const SmallBinaryFieldPolynomial & f1 ,
		  const SmallBinaryFieldPolynomial & f2 ) const;
	};

	bool DecodedPolynomialComparator::operator()
	( const SmallBinaryFieldPolynomial & f1 ,
	  const SmallBinaryFieldPolynomial & f2 ) const {

		int hd1 , hd2;

		hd1 = 0;
		hd2 = 0;

		for ( int i = 0 ; i < n ; i++ ) {
			if ( f1.eval(this->x[i]) != this->y[i] ) {
				hd1++;
			}
			if ( f2.eval(this->x[i]) != this->y[i] ) {
				hd2++;
			}
		}

		return hd1 < hd2;
	}

	/**
	 * @brief
	 *            Attempts to solve a Reed-Solomon list decoding problem.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	bool GuruswamiSudanDecoder::decode
	( const uint32_t *x , const uint32_t *y ,
	  int n , int k , int m , const SmallBinaryField & gf ) {


		this->overallTime = clock();

		if ( n <= 0 || k <= 0 || m <= 0 || k > n ) {
			cerr << "GuruswamiSudanDecoder::decode: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		this->decodedList.clear();
		delete this->Qptr;
		this->Qptr = NULL;

		this->interpolationTime = clock();
		this->Qptr =
				new SmallBinaryFieldBivariatePolynomial
					(interpolate(x,y,n,k,m,gf));
		this->interpolationTime = clock()-this->interpolationTime;

		this->rootTime = clock();
		vector<SmallBinaryFieldPolynomial> ps = roots(*(this->Qptr),k);
		this->rootTime = clock() - this->rootTime;

		for ( int i = 0 ; i < (int)ps.size() ; i++ ) {

			int hammingDistance = 0;
			for ( int j = 0 ; j < n ; j++ ) {
				uint32_t c = ps[i].eval(x[j]);
				if ( c != y[j] ) {
					++hammingDistance;
				}
			}

			if ( hammingDistance <= n-k-1 ) {
				this->decodedList.push_back(ps[i]);
			}

		}

		sort
		(this->decodedList.begin(),this->decodedList.end(),
		 DecodedPolynomialComparator(x,y,n));

		this->overallTime = clock()-this->overallTime;

		return this->decodedList.size() > 0;
	}

	/**
	 * @brief
	 *            Constructor.
	 */
	GuruswamiSudanDecoder::GuruswamiSudanDecoder() {
		this->interpolationTime = 0;
		this->rootTime = 0;
		this->overallTime = 0;
		this->Qptr = NULL;
	}

	/**
	 * @brief
	 *            Destructor.
	 */
	GuruswamiSudanDecoder::~GuruswamiSudanDecoder() {
		delete this->Qptr;
	}

	/**
	 * @brief
	 *            Access the list of recently decoded polynomials.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	const vector<SmallBinaryFieldPolynomial> &
	GuruswamiSudanDecoder::getDecodedList() const {

		if ( this->Qptr == NULL ) {
			cerr << "GuruswamiSudanDecoder: Missing recent decoding attempt."
				 << endl;
			exit(EXIT_FAILURE);
		}

		return this->decodedList;
	}

	/**
	 * @brief
	 *            Access the bivariate polynomial recently computed
	 *            in the interpolation step of the Guruswami-Sudan
	 *            decoder.
	 *
	 * @return
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	const SmallBinaryFieldBivariatePolynomial &
	GuruswamiSudanDecoder::getBivariatePolynomial() const {

		if ( this->Qptr == NULL ) {
			cerr << "GuruswamiSudanDecoder: Missing recent decoding attempt."
				 << endl;
			exit(EXIT_FAILURE);
		}

		return *(this->Qptr);
	}

	/**
	 * @brief
	 *            Access the time in seconds the latest interpolation
	 *            step during <code>decode()</code> consumed.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	double GuruswamiSudanDecoder::getInterpolationSecs() const {

		if ( this->Qptr == NULL ) {
			cerr << "GuruswamiSudanDecoder: Missing recent decoding attempt."
				 << endl;
			exit(EXIT_FAILURE);
		}

		return
			(double)(this->interpolationTime) / (double)(CLOCKS_PER_SEC);
	}

	/**
	 * @brief
	 *            Access the time in seconds the latest root
	 *            step during <code>decode()</code> consumed.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	double GuruswamiSudanDecoder::getRootSecs() const {

		if ( this->Qptr == NULL ) {
			cerr << "GuruswamiSudanDecoder: Missing recent decoding attempt."
				 << endl;
			exit(EXIT_FAILURE);
		}

		return (double)(this->rootTime) / (double)(CLOCKS_PER_SEC);
	}

	/**
	 * @brief
	 *            Access the time in seconds decoding attempt running
	 *            <code>decode()</code> consumed.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	double GuruswamiSudanDecoder::getOverallSecs() const {

		if ( this->Qptr == NULL ) {
			cerr << "GuruswamiSudanDecoder: Missing recent decoding attempt."
				 << endl;
			exit(EXIT_FAILURE);
		}

		return (double)(this->overallTime) / (double)(CLOCKS_PER_SEC);
	}

	inline static int WeightedDegreeCompare
	( const pair<int,int> & w1 , const pair<int,int> & w2 , int a , int b ) {

		long wd1 , wd2;

		wd1 = a*w1.first+b*w1.second;
		wd2 = a*w2.first+b*w2.second;

		if ( wd1 < wd2 ) {
			return -1;
		} else if ( wd1 > wd2 ) {
			return 1;
		} else if ( w1.second < w2.second ) {
			return -1;
		} else if ( w1.second > w2.second ) {
			return 1;
		} else if ( w1.first < w2.first ) {
			return -1;
		} else if ( w1.first > w2.first ) {
			return 1;
		} else {
			return 0;
		}
	}

	static void TrifonovRand
	( SmallBinaryFieldBivariatePolynomial & f ,
	  const vector<SmallBinaryFieldBivariatePolynomial> & B ,
	  const SmallBinaryField & gf ) {

		SmallBinaryFieldBivariatePolynomial g(gf);

		f.setZero();

		int n = B.size();
		for ( int i = 0 ; i < n ; i++ ) {
			mul(g,B[i],gf.random());
			add(f,f,g);
		}
	}

	static int TrifonovDelta
	( const vector<SmallBinaryFieldBivariatePolynomial> & A , int k ) {

		pair<int,int> d;
		int dlt = 0;

		int n = A.size();
		for ( int i = 0 ; i < n ; i++ ) {
			d = A[i].deg(1,k);
			dlt += d.first;
		}

		return dlt;
	}

	static void TrifonovReduce
	( std::vector<SmallBinaryFieldBivariatePolynomial> & S ,
	  const SmallBinaryFieldBivariatePolynomial & P , int k ,
	  const SmallBinaryField & gf ) {

		SmallBinaryFieldBivariatePolynomial tmp(gf);

		int i = S.size();

		S.push_back(P);

		for (;;) {

			std::pair<int,int> di = S[i].deg(1,k);

			int dxi , dy;
			dxi = di.first;
			dy  = di.second;

			std::pair<int,int> dj;
			int j = -1;

			for ( int l = 0 ; l < i ; l++ ) {
				dj = S[l].deg(1,k);
				if ( di.second == dj.second ) {
					j = l;
					break;
				}
			}

			if ( j < 0 ) {
				break;
			}

			int dxj;
			dxj = dj.first;

			uint32_t ci , cj;
			ci = S[i].getCoeff(dxi,dy);
			cj = S[j].getCoeff(dxj,dy);

			if ( dxi <= dxj ) {

				uint32_t c = gf.div(cj,ci);
				mul(tmp,S[i],c);
				leftShiftX(tmp,tmp,dxj-dxi);
				sub(tmp,S[j],tmp);
				S[j].swap(S[i]);
				S[i].swap(tmp);

			} else {

				uint32_t c = gf.div(ci,cj);
				mul(tmp,S[j],c);
				leftShiftX(tmp,tmp,dxi-dxj);
				sub(S[i],S[i],tmp);

			}
		}

		if ( S.back().isZero() ) {
			S.pop_back();
		}
	}

	static void TrifonovMerge
	( vector<SmallBinaryFieldBivariatePolynomial> & Q ,
	  const vector<SmallBinaryFieldBivariatePolynomial> & P ,
	  const vector<SmallBinaryFieldBivariatePolynomial> & S ,
	  int delta0 , int k , const SmallBinaryField & gf ) {


		if ( &Q == &P || &Q == &S ) {
			std::vector<SmallBinaryFieldBivariatePolynomial> tQ;
			TrifonovMerge(tQ,P,S,delta0,k,gf);
			Q = tQ;
			return;
		}

		SmallBinaryFieldBivariatePolynomial tmp1(gf) , tmp2(gf) , tmp3(gf);

		int u , v;
		u = (int)(P.size()-1);
		v = (int)(S.size()-1);

		Q.assign(u+v+1,tmp1);

		for ( int i = 0 ; i <= u+v ; i++ ) {

			std::pair<int,int> dmin;
			bool chosen = false;

			for ( long j = 0 ; j <= v ; j++ ) {

				if ( i-j >= 0 && i-j <= u ) {

					mul(tmp1,P[i-j],S[j]);

					std::pair<int,int> d = tmp1.deg(1,k);
					if ( !chosen || WeightedDegreeCompare(d,dmin,1,k) < 0 ) {
						dmin = d;
						Q[i] = tmp1;
						chosen = true;
					}

				}

			}
		}

		int delta = TrifonovDelta(Q,k);

		while( delta > delta0 ){

			TrifonovRand(tmp1,P,gf);
			TrifonovRand(tmp2,S,gf);
			mul(tmp3,tmp1,tmp2);

			TrifonovReduce(Q,tmp3,k,gf);

			delta = TrifonovDelta(Q,k);
		}

	}

	static SmallBinaryFieldBivariatePolynomial TrifonovInterpolate
	( const uint32_t *x , const uint32_t *y ,
	  int n , int k , int r , const SmallBinaryField & gf ) {

		SmallBinaryFieldBivariatePolynomial Q(gf);
		SmallBinaryFieldPolynomial T(gf);
		std::vector<SmallBinaryFieldBivariatePolynomial> G;

		{
			SmallBinaryFieldPolynomial phi(gf);
			phi.buildFromRoots(x,n);
			G.push_back(SmallBinaryFieldBivariatePolynomial(phi));

			for ( int i = 0 ; i < n ; i++ ) {

				SmallBinaryFieldPolynomial monom(gf);
				monom.setX();
				monom.setCoeff(0,gf.neg(x[i]));

				SmallBinaryFieldPolynomial tmp(gf);
				div(tmp,phi,monom);

				uint32_t den = 1;
				for ( int j = 0 ; j < n ; j++ ) {
					if ( j != i ) {
						den = gf.mul(den,gf.sub(x[i],x[j]));
					}
				}
				den = gf.inv(den);
				den = gf.mul(den,y[i]);
				mul(tmp,tmp,den);
				add(T,T,tmp);
			}
		}

		SmallBinaryFieldPolynomial ONE(gf); ONE.setCoeff(0,1);
		SmallBinaryFieldBivariatePolynomial tmp1(gf), tmp2(gf);

		int j = 0;

		for (;;) {

			// Let 'tmp1(X,Y)<-Y-T(X)' and note that for binary fields
			// there is no need to negate 'T(X)'.
			tmp1.setY();
			tmp1.setCoeffY(0,T);
			leftShiftY(tmp1,tmp1,j);

			TrifonovReduce(G,tmp1,k,gf);
			j = G.size()-1;

			tmp1.setZero();
			tmp1.setCoeffY(j,ONE);
			leadTerm(tmp2,G[j],1,k);
			if ( tmp2.equals(tmp1) ) {
				break;
			}

		}


		std::vector<SmallBinaryFieldBivariatePolynomial> B , C;
		B = G;

		int m = thimble::MathTools::numBits(r)-1;
		int R = 1;

		for ( int j = m - 1 ; j >= 0 ; j-- ) {
			R <<= 1;
			TrifonovMerge(C,B,B,n*(R*(R+1))/2,k,gf);
			B.swap(C);
			if ( r & (0x1<<j) ) {
				R++;
				TrifonovMerge(C,B,G,n*(R*(R+1))/2,k,gf);
				B.swap(C);
			}
		}

		Q = B.front();
		pair<int,int> dmin = Q.deg(1,k);

		for ( int i = 1 ; i < (int)B.size() ; i++ ) {
			std::pair<int,int> d = B[i].deg(1,k);
			if ( WeightedDegreeCompare(d,dmin,1,k) < 0 ) {
				Q = B[i];
				dmin = d;
			}
		}

		return Q;
	}

	static vector<SmallBinaryFieldPolynomial> RothRuckensteinRoots
	( const SmallBinaryFieldBivariatePolynomial & Q , int k , int i ,
	  const SmallBinaryField & gf ) {

		vector<SmallBinaryFieldPolynomial> roots;

	    int r = 0;
	    SmallBinaryFieldBivariatePolynomial Myx(gf);
	    swapXY(Myx,Q);

	    while( Myx.degX(0) < 0 ) {
	    	rightShiftY(Myx,Myx,1);
	    	++r;
	    }

	    int numGamma = Myx.degX(0);
	    uint32_t *gamma = NULL;
	    if ( numGamma > 0 ) {
	    	gamma = (uint32_t*)malloc(numGamma*sizeof(uint32_t));
	    	if ( gamma == NULL ) {
	    		cerr << "RothRuckensteinRoots: Out of memory." << endl;
	    		exit(EXIT_FAILURE);
	    	}
	    	if ( numGamma == 1 ) {
	    		gamma[0] = gf.div
	    			(Myx.getCoeffY(0).getCoeff(0),Myx.getCoeffY(0).getCoeff(1));
	    	} else {
	    		numGamma = Myx.getCoeffY(0).findRoots(gamma);
	    	}
	    }

	    SmallBinaryFieldBivariatePolynomial Mxy(gf);
    	SmallBinaryFieldBivariatePolynomial nQ(gf);
    	SmallBinaryFieldBivariatePolynomial y(gf);
	    SmallBinaryFieldPolynomial phi(gf);

	    rightShiftX(Mxy,Q,r);

	    for (int j = 0 ; j < numGamma ; j++ ) {

	    	phi.setZero();
	    	phi.setCoeff(i,gamma[j]);

	    	if ( i == k ) {
	    		roots.push_back(phi);
	    		break;
	    	}

	    	y.setY();
	    	y.setCoeff(0,0,gamma[j]);
	    	evalY(nQ,Mxy,y);

	    	y.setZero();
	    	y.setCoeff(1,1,1);
	    	evalY(nQ,nQ,y);

	    	vector<SmallBinaryFieldPolynomial> next_roots;
	    	//Recursive call
	    	next_roots = RothRuckensteinRoots(nQ,k,i+1,gf);

	    	for ( int l = 0 ; l < (int)next_roots.size() ; l++ ) {
	    		roots.push_back(SmallBinaryFieldPolynomial(gf));
	    		add(roots.back(),next_roots[l],phi);
	    	}

	    }

	    free(gamma);

	    return roots;
	}

	SmallBinaryFieldBivariatePolynomial GuruswamiSudanDecoder::interpolate
	( const uint32_t *x , const uint32_t *y ,
	  int n , int k , int m , const SmallBinaryField & gf ) const {

		if ( n <= 0 || k <= 0 || m <= 0 || k > n ) {
			cerr << "GuruswamiSudanDecoder::interpolate: Bad arguments." << endl;
			exit(EXIT_FAILURE);
		}

		return TrifonovInterpolate(x,y,n,k-1,m,gf);
	}

	/**
	 * @brief
	 *            Implementation of the root step in the Guruswami-Sudan
	 *            list decoding algorithm.
	 *
	 * @details
	 *            see 'GuruswamiSudanDecoder.h'
	 */
	vector<SmallBinaryFieldPolynomial> GuruswamiSudanDecoder::roots
	( const SmallBinaryFieldBivariatePolynomial & Q , int k ) const {

		if ( Q.isZero() ) {
			cerr << "GuruswamiSudanDecoder::roots: "
				 << "Bivariate polynomial is zero." << endl;
			exit(EXIT_FAILURE);
		}

		if ( k <= 0 ) {
			cerr << "GuruswamiSudanDecoder::roots: "
				 << "Root size must be greater than zero." << endl;
			exit(EXIT_FAILURE);
		}

		return RothRuckensteinRoots(Q,k-1,0,Q.getField());
	}


}




