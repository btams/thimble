/*
 *  THIMBLE --- Research Library for Development and Analysis of
 *  Fingerprint-Based Biometric Cryptosystems.
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
 * @file MinutiaeFuzzyVault.cpp
 *
 * @brief
 *           Implements the functionalities that are provided by
 *           'MinutiaeFuzzyVault.h' which provides an implementation
 *           of a minutiae fuzzy vault.
 *
 * @details
 *           see 'MinutiaeFuzzyVault.h'
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <stdint.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <iostream>

#include <thimble/math/numbertheory/SmallBinaryField.h>
#include <thimble/math/numbertheory/SmallBinaryFieldPolynomial.h>
#include <thimble/security/SHA.h>
#include <thimble/security/FuzzyVaultTools.h>
#include <thimble/finger/MinutiaeRecord.h>
#include <thimble/finger/FingerTools.h>

#include <thimble/finger/MinutiaeFuzzyVault.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *           Constructor which creates an instance of the minutiae
	 *           fuzzy vault to protect minutiae template whose
	 *           coordinates vary within the specified spatial region.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	MinutiaeFuzzyVault::MinutiaeFuzzyVault
	( int width , int height , bool tryRandom) {

		// Check whether the specified region has positive ranges. ...
		if ( width <= 0 || height <= 0 ) {

			// ... . If not, verbose and exit
			cerr << "MinutiaeFuzzyVault: Fingerprint image must be of "
				 <<" width and height both greater than 0." << endl;
			exit(EXIT_FAILURE);
		}

		// Adopt the specified parameters
		this->width  = width;
		this->height = height;
		this->tryRandom   = tryRandom;

		// *******************************************************************
		// *** BEGIN: Set default values *************************************
		// *******************************************************************

		this->angleWeight = 11.459;

		this->n    = 224;
		this->tmin = 18;
		this->tmax = 24;
		this->k    = 9;

		this->minInterVaultDistance = 25.0;
		this->maxMatchDistance      = 30.0;

		// *******************************************************************
		// *** END: Set default values ***************************************
		// *******************************************************************

		// Initialize the vault points as NULL
		this->vaultX = NULL;
		this->vaultY = NULL;

		// The finite field will be the finite field of 2^16 elements
		this->gfPtr = new SmallBinaryField(16);

		// Initialize the zero hash value
		memset(this->hash,0,20);
	}

	/**
	 * @brief
	 *           Destructor which frees all memory helt by the minutiae
	 *           fuzzy vault.
	 */
	MinutiaeFuzzyVault::~MinutiaeFuzzyVault() {
		delete this->gfPtr;
		clear();
	}

	/**
	 * @brief
	 *           Comparator class which is used to sort the vault minutiae
	 *           w.r.t. lexicographical order using the <code>sort()</code>
	 *           method from the STL.
	 */
	class LexicographicalMinutiaComparator {
	public:

		/**
		 * @brief
		 *            Returns <code>true</code> if and only if <code>a</code>
		 *            is lexicographically smaller than <code>b</code>.
		 *
		 * @return
		 *            <code>true</code> if <code>a</code> is lexicographically
		 *            smaller than <code>b</code> and <code>false</code>
		 *            otherwise.
		 */
		inline bool operator()( const Minutia & a , const Minutia & b ) {

			// The abscissa is defined to be the most relevant coordinate
			// for the lexicographical order.
			if ( a.getX() < b.getX() ) {
				return true;
			} else if ( a.getX() > b.getX() ) {
				return false;
			}

			// The ordinate is defined to be the second most relevant
			// coordinate for the lexicographical order.
			if ( a.getY() < b.getY() ) {
				return true;
			} else if ( a.getY() > b.getY() ) {
				return false;
			}

			// If the position of both minutiae is equals, the angle
			// determines the order.
			if ( a.getAngle() < b.getAngle() ) {
				return true;
			} else if ( a.getAngle() > b.getAngle() ) {
				return false;
			}

			// The other attributes of the minutiae are not relevant
			// for the vault construction.
			return false;
		}
	};

	/**
	 * @brief
	 *           Protects the given minutiae template using the instance.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	bool MinutiaeFuzzyVault::enroll( const MinutiaeView & view ) {

		// Check whether the specified vault parameters make sense
		if ( this->k > this->tmin || this->k > this->tmax ||
			 this->k > this->n || this->tmin > this->tmax ||
			 this->tmin > this->n || this->tmax > this->n ) {

			cerr << "MinutiaeFuzzyVault::enroll: Bad vault parameters."
				 << endl;
			exit(EXIT_FAILURE);
		}

		// Select minutiae: The minutiae are sorted w.r.t. their
		// respective quality and then all of those minutiae
		// are removed that do not keep a distance of at least
		// 'getMinimalInterVaultDistance()' to other minutiae; The minutiae
		// that remain are returned.
		MinutiaeView selectedView = select(view);

		// If not at least 'tmin' could be selected, the minutiae template
		// will not be protected; in this case, we return 'false';
		if ( selectedView.getMinutiaeCount() < this->tmin ) {
			return false;
		}

		// We will protect at most the 'tmax' minutiae that are of
		// the best quality. Note, that 'selectedView' is sorted w.r.t.
		// the minutiae quality, such that the first minutiae is of best
		// and the last minutiae is of worst quality.
		int t = min(selectedView.getMinutiaeCount(),this->tmax);

		// Dismiss the data related with a template protected by this instance
		// if any
		clear();

		// Reserve memory
		this->vaultMinutiae.reserve(this->n);
		this->vaultX = (uint32_t*)malloc( this->n * sizeof(uint32_t) );
		this->vaultY = (uint32_t*)malloc( this->n * sizeof(uint32_t) );
		if ( this->vaultX == NULL || this->vaultY == NULL ) {
			cerr << "MinutiaeFuzzyVault::enroll: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// Temporarily, let the first 't' vault minutiae be the first 't'
		// selected minutiae that are of best quality.
		for ( int i = 0 ; i < t ; i++ ) {
			this->vaultMinutiae.push_back(selectedView.getMinutia(i));
		}

		// Loop 'n-t' times to generate 'n-t' chaff minutiae
		for ( int j = t ; j < this->n ; j++ ) {

			// Loop until a randomly chosen chaff minutia
			// fulfills all requirements. More precisely, a chaff minutia
			// must keep mutual distance to all chaff minutiae already
			// generated as well as to genuine minutiae
			for(;;) {

				// First, generate a candidate for the chaff minutia
				Minutia chaff = generateChaffMinutia();

				// Determine the minimal distance over all vault minutiae
				// included in the vector 'vaultMinutiae'
				double dmin = DBL_MAX;
				for ( int i = 0 ; i < j ; i++ ) {
					double d = dist(this->vaultMinutiae[i],chaff);
					if ( d < dmin ) {
						dmin = d;
					}
				}

				// If the minimal distance is at least 'minInterVaultDistance'
				// the chaff minutia is valid and the loop stops; otherwise
				// the loop continues until a valid chaff minutia could be
				// generated
				if ( dmin > this->minInterVaultDistance ) {
					this->vaultMinutiae.push_back(chaff);
					break;
				}
			}
		}

		// Now, as all vault minutiae are contained in 'vaultMinutiae' we need
		// to disperse the first 't' minutiae in the rest of the vector to such
		// that genuine minutiae can not be immediately be determined. An
		// elegant way to guarantee this is to sort the minutiae w.r.t. their
		// lexicographical order to reach an arrangement of the vault points
		// that can be computed easily at any time, no matter how the vault
		// minutiae have been shuffled.
		sort
			(this->vaultMinutiae.begin(),this->vaultMinutiae.end(),
			 LexicographicalMinutiaComparator());

		// Now, generate a random polynomial of degree '<k'...
		SmallBinaryFieldPolynomial f(getField());
		f.random(this->k,this->tryRandom);

		// ... and save its hash SHA-1 hash value
        SHA().hash(this->hash,f.getData(),f.deg()+1);

		// Build the vault points
		for ( int i = 0 ; i < this->n ; i++ ) {

			// Determine whether the 'i'th vault minutia is genuine
			bool isGenuine = false;
			for ( int j = 0 ; j < t ; j++ ) {
				if ( selectedView.getMinutia(j) == this->vaultMinutiae[i] ) {
					isGenuine = true;
					break;
				}
			}

			// The abscissa of the vault point corresponding to the 'i'th
			// vault minutia is the field element that is encoded by the
			// integer 'i'.
			this->vaultX[i] = (uint32_t)i;

			uint32_t z = f.eval(this->vaultX[i]);

			if ( isGenuine ) {
				// If the vault point is genuine then set the ordinate to
				// the evaluation of the polynomial at the abscissa... ;
				this->vaultY[i] = z;
			} else {
				// ...; otherwise, choose a random ordinate value that is
				// different from the evaluation of the polynomial at the
				// abscissa.
				uint32_t y;
				do {
					y = this->gfPtr->random(this->tryRandom);
				} while ( y == z );
				this->vaultY[i] = y;
			}
		}

		//The template is now successfully protected
		return true;
	}

	/**
	 * @brief
	 *           Attempts to open the vault using an aligned query minutiae
	 *           template.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	bool MinutiaeFuzzyVault::open
		( SmallBinaryFieldPolynomial & f , const MinutiaeView & view ) const {

		// Check whether the instance protects a template. ...
		if ( this->vaultMinutiae.size() == 0 ) {
			// ... . If not, verbose and exit.
			cerr << "MinutiaeFuzzyVault::open: The vault does not protect a "
			     << "template. Enroll first." << endl;
			exit(EXIT_FAILURE);
		}

		// Keeps track of whether decoding was successful
		bool state = false;

		// Select minutiae in the same way as on enrollment: The minutiae are
		// sorted w.r.t. their respective quality and then all of those
		// minutiae are removed that do not keep a distance of at least
		// 'getMinimalInterVaultDistance()' to other minutiae; The minutiae
		// that remain are returned.
		MinutiaeView selectedView = select(view);

		// We only consider at most the first 'tmax' selected minutiae that
		// are of the best quality
		int t = min(selectedView.getMinutiaeCount(),this->tmax);

		// Reserve memory to save the unlocking points; the size
		// of the unlocking points can be at most 't<=tmax'
		uint32_t *x , *y;
		x = (uint32_t*)malloc( t * sizeof(uint32_t));
		y = (uint32_t*)malloc( t * sizeof(uint32_t));
		if ( x == NULL || y == NULL ) {
			cerr << "MinutiaeFuzzyVault::open: Out of memory." << endl;
			exit(EXIT_FAILURE);
		}

		// 's' will be the size of the unlocking set
		int s = 0;
		// Extract unlocking points by iterating over all query minutiae.
		for ( int i = 0 ; i < t ; i++ ) {

			// Determine the index and the distance of the vault minutiae
			// that is of minimal distance to the query minutia
			int minIndex = -1;
			double minDistance = DBL_MAX;
			for ( int j = 0 ; j < (int)(this->vaultMinutiae.size()) ; j++ ) {
				double d = dist(selectedView.getMinutia(i),vaultMinutiae[j]);
				if ( d < minDistance ) {
					minIndex = j;
					minDistance = d;
				}
			}

			// If vault minutia matches the query minutia ...
			if ( minDistance <= this->maxMatchDistance ) {

				// ... check whether the corresponding unlocking point
				// is already contained in the unlocking set. ...
				bool alreadyContained = false;
				for ( int j = 0 ; j < s ; j++ ) {
					if ( x[j] == this->vaultX[minIndex] ) {
						alreadyContained = true;
						break;
					}
				}

				// ... . If not, append the point at the end of the
				// unlocking set.
				if ( !alreadyContained ) {
					x[s] = this->vaultX[minIndex];
					y[s] = this->vaultY[minIndex];
					++s;
				}

			}
		}

		// Pass the unlocking set to the decoder and keep track of whether
		// decoding was successful to ...
		state = decode(f,x,y,s);

		// ... free the data that was used to hold the unlocking points ...
		free(x);
		free(y);

		// ... and return afterwards.
		return state;
	}

	/**
	 * @brief
	 *           Attempts to decode a polynomial given unlocking points.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	bool MinutiaeFuzzyVault::decode
	( SmallBinaryFieldPolynomial & f ,
	  const uint32_t *x , const uint32_t *y ,
	  int t ) const {

		// Check whether the instance protects a minutiae template. ...
		if ( this->vaultMinutiae.size() == 0 ) {
			// ... . If not, verbose and exit.
			cerr << "MinutiaeFuzzyVault::decode: Instance does not protect "
				 << "a template. Enroll first." << endl;
			exit(EXIT_FAILURE);
		}

		// The decoder is a wrapper around the brute-force decoder of the
		// FuzzyVaultTools-class.
		return FuzzyVaultTools::bfdecode(f,x,y,t,this->k,this->hash);
	}

	/**
	 * @brief
	 *           Sorts the given minutiae w.r.t. their respective quality
	 *           and dismisses all minutiae that do not keep sufficient
	 *           mutual distance.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	MinutiaeView MinutiaeFuzzyVault::select
	( const MinutiaeView & view ) const {

		// Let 'tmp' be the input minutiae sorted w.r.t. their respective
		// quality. Note that the used sorting algorithm is stable such that
		// the relative ordering of minutiae of equal quality will not
		// change.
		MinutiaeView tmp = view;
		tmp.sortWithRespectToMinutiaeQuality();

		// In 'selectedView' we collect the minutiae from 'tmp' that keep
		// sufficient distance from other minutiae in 'tmp'. We adopt the
		// finger position, impression type, and finger quality.
		MinutiaeView selectedView;
		selectedView.ensureCapacity(tmp.getMinutiaeCount());
		selectedView.setFingerPosition(tmp.getFingerPosition());
		selectedView.setImpressionType(tmp.getImpressionType());
		selectedView.setQuality(tmp.getFingerQuality());

		// Iterate over the sorted minutiae
		for ( int i = 0 ; i < tmp.getMinutiaeCount() ; i++ ) {

			// Determine the minimal distance of the current minutiae
			// to the other minutiae
			double dmin = DBL_MAX;
			for ( int j = 0 ; j < tmp.getMinutiaeCount() ; j++ ) {
				if ( i != j ) {
					double d = dist(tmp.getMinutia(i),tmp.getMinutia(j));
					if ( d < dmin ) {
						dmin = d;
					}
				}
			}

			// If the current minutia is sufficiently separated from the other
			// minutiae ...
			if ( dmin >= this->minInterVaultDistance ) {
				// ... adopt its relevant data and ...
				Minutia minutia
					(tmp.getMinutia(i).getX(),
					 tmp.getMinutia(i).getY(),
					 tmp.getMinutia(i).getAngle());
				// ... select it
				selectedView.addMinutia(minutia);
			}
		}

		return selectedView;
	}

	/**
	 * @brief
	 *           Measures the distance of two minutiae.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	double MinutiaeFuzzyVault::dist
	( const Minutia & a , const Minutia & b ) const {
		return FingerTools::dist(a,b,this->angleWeight);
	}

	/**
	 * @brief
	 *           Generates a candidate for a chaff minutia.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	Minutia MinutiaeFuzzyVault::generateChaffMinutia() const {

		Minutia chaff(

				// Random abscissa in the specified range
				MathTools::rand32(this->tryRandom)%this->width,

				// Random ordinate in the specified range
				MathTools::rand32(this->tryRandom)%this->height,

				// Random angle from 256 quanta
				MathTools::rand8 (this->tryRandom) / 256.0 * (M_PI+M_PI)

		);

		return chaff;
	}

	/**
	 * @brief
	 *           Specifies how significant the minutiae angles are taken
	 *           into account when the similarity of two minutiae is
	 *           measured.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setAngleWeight( double angleWeight ) {

		if ( angleWeight <= 0.0 ) {
			cerr << "MinutiaeFuzzyVault::setAngleWeight: "
				 << "Must be greater than or equals 0.0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0 ) {
			cerr << "MinutiaeFuzzyVault::setAngleWeight: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->angleWeight = angleWeight;
	}

	/**
	 * @brief
	 *           Specifies the size of the vault.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setVaultSize( int n ) {

		if ( n <= 0 ) {
			cerr << "MinutiaeFuzzyVault::setVaultSize: "
				 << "Must be greater than 0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0 ) {
			cerr << "MinutiaeFuzzyVault::setVaultSize: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->n = n;
	}

	/**
	 * @brief
	 *           Specifies how many well-separated genuine minutiae are
	 *           required in a minutiae template to successfully be
	 *           enrolled.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setMinimalGenuineSetSize( int tmin ) {

		if ( tmin <= 0 ) {
			cerr << "MinutiaeFuzzyVault::setMinimalGenuineSetSize: "
				 << "Must be greater than 0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0) {
			cerr << "MinutiaeFuzzyVault::setMinimalGenuineSetSize: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->tmin = tmin;
	}

	/**
	 * @brief
	 *           Specifies how many well-separated genuine minutiae are
	 *           at most selected from a minutiae template that is
	 *           passed to <code>enroll()</code>
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setMaximalGenuineSetSize( int tmax ) {

		if ( tmax <= 0 ) {
			cerr << "MinutiaeFuzzyVault::setMaximalGenuineSetSize: "
				 << "Must be greater than 0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0 ) {
			cerr << "MinutiaeFuzzyVault::setMaximalGenuineSetSize: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->tmax = tmax;
	}

	/**
	 * @brief
	 *           Specifies the length of the secret polynomial that is
	 *           randomly generated on enrollment for being binded to the
	 *           genuine minutiae.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setSecretSize( int k ) {

		if ( k <= 0 ) {
			cerr << "MinutiaeFuzzyVault::setSecretSetSize: "
				 << "Must be greater than 0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0) {
			cerr << "MinutiaeFuzzyVault::setSecretSize: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->k = k;
	}

	/**
	 * @brief
	 *           Specifies the minimal distance between vault minutiae.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setMinimalInterVaultDistance
	( double minInterVaultDistance ) {

		if ( minInterVaultDistance <= 0.0 ) {
			cerr << "MinutiaeFuzzyVault::setInterVaultDistance: "
				 << "Must be greater than or equals 0.0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		if ( this->vaultMinutiae.size() != 0 ) {
			cerr << "MinutiaeFuzzyVault::setInterVaultDistance: "
				 << "Already enrolled. Clear first."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->minInterVaultDistance = minInterVaultDistance;
	}

	/**
	 * @brief
	 *           Specifies how much two minutiae are allowed to differ
	 *           such that they are considered to match.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::setMaximalMatchDistance
	( double matchDistance ) {

		if ( matchDistance <= 0.0 ) {
			cerr << "MinutiaeFuzzyVault::setMatchDistance: "
				 << "Must be greater than or equals 0.0."
				 << endl;
			exit(EXIT_FAILURE);
		}

		this->maxMatchDistance = matchDistance;
	}

	/**
	 * @brief
	 *           Ensures that the vault dismisses all data related to
	 *           a protected minutiae template, if any.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	void MinutiaeFuzzyVault::clear() {
		this->vaultMinutiae.clear();
		free(this->vaultX);
		free(this->vaultY);
		this->vaultX = NULL;
		this->vaultY = NULL;
		memset(this->hash,0,20);
	}

	/**
	 * @brief
	 *           Access the SHA-1 hash value of the secret polynomial
	 *           that is binded to the genuine vault minutiae.
	 *
	 * @details
	 *           see 'MinutiaeFuzzyVault.h'
	 */
	const uint32_t* MinutiaeFuzzyVault::getHashValue() const {

		if ( this->vaultMinutiae.size() == 0 ) {
			std::cerr
			 << "MinutiaeFuzyyVault::getHashValue: "
			 << "No finger was yet enrolled to the vault."
			 << std::endl;
			exit(EXIT_FAILURE);
		}

		return this->hash;
	}
}



