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
 * @file MinutiaeRecord.cpp
 *
 * @brief
 *            Implements the functionalities from 'MinutiaeRecord.h' which
 *            provides classes representing minutiae, minutiae templates, and
 *            minutiae records in the ISO 19794-2:2005 standard.
 *
 * @details
 *            see 'MinutiaeRecord.h'
 *
 * @author Benjamin Tams
 */

#define _USE_MATH_DEFINES
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iostream>

#include "config.h"

#include <thimble/finger/MinutiaeRecord.h>

using namespace std;

/**
 * @brief The library's namespace
 */
namespace thimble {

	/**
	 * @brief
	 *            Initializes this minutia as specified.
	 *
	 * @details
	 *            see 'MinutiaRecord.h'
	 *
	 */
	void Minutia::init
	( double x , double y , double angle , MINUTIA_TYPE_T typ ,
	  int quality ) {

		this->x = x;
		this->y = y;

		if ( angle < 0.0 || angle >= M_PI+M_PI ) {

			double dx , dy;
			dx = cos(angle);
			dy = sin(angle);

			this->angle = atan2(dy,dx);

			while ( this->angle < 0.0 ) {
				this->angle += M_PI+M_PI;
			}

			while ( this->angle >= M_PI+M_PI ) {
				this->angle -= M_PI+M_PI;
			}

		} else {
			this->angle = angle;
		}

		this->typ = typ;

		if ( quality < 0 || quality > 100 ) {
			cerr << "A minutia's quality must range between 1 and 100 and "
				 << "should be 0 only to indicate that no quality is "
				 << "associated with the minutia." << endl;
			exit(EXIT_FAILURE);
		}

		this->quality = quality;
	}

	/**
	 * @brief Prints a text representation of the given minutia to an output
	 *        stream.
	 *
	 * @details
	 *            see 'MinutiaeRecord.h'
	 */
	ostream & operator<<( ostream & out , const Minutia & minutia ) {
		out << minutia.getX() << " "
			<< minutia.getY() << " "
			<< minutia.getAngle() / (M_PI+M_PI) * 360.0 << " "
			<< minutia.getType() << " "
			<< minutia.getQuality();
		return out;
	}

	/**
	 * @brief Class used to sort minutiae w.r.t. their associated quality.
	 */
	class MinutiaQualityComparator {
	public:

		/**
		 * @brief
		 *           Returns <code>true</code> if the minutiae <code>a</code>
		 *           is associated with a higher quality than <code>b</code>;
		 *           otherwise <code>false</code> is returned.
		 *
		 * @return
		 *           <code>true</code> if <code>a</code> is associated with
		 *           a higher quality than <code>b</code> and
		 *           <code>false</code> otherwise.
		 */
		inline bool operator()
		( const Minutia & a , const Minutia & b ) const {
			return a.getQuality() > b.getQuality();
		}
	};

	/**
	 * @brief
	 *            Sorts the minutiae in the view w.r.t. their
	 *            associated quality.
	 *
	 * @details
	 *            see 'MinutiaeRecord.h'
	 */
	void MinutiaeView::sortWithRespectToMinutiaeQuality() {

		// We use the sort implementation from the STL and the
		// MinutiaeQualityComparator class to implement the method

		stable_sort(this->minutiae.begin(),this->minutiae.end(),
					MinutiaQualityComparator());
	}

	/**
	 * @brief
	 *            Writes the minutiae view to the specified output
	 *            stream.
	 *
	 * @details
	 *            see 'MinutiaeRecord.h'
	 */
	bool MinutiaeView::write( FILE *out ) const {

		unsigned char header[4];
		unsigned char minutiaData[6];

		header[0] = (unsigned char) this->fingerPosition;
		header[1] = (unsigned char) (this->viewNumber + 0x10 * this->impressionType);
		header[2] = (unsigned char) this->fingerQuality;
		header[3] = (unsigned char) this->minutiae.size();

		if ( fwrite(header,1,4,out) != 4 ) {
			return false;
		}

		for (int i = 0; i < (int)(this->minutiae.size()) ; i++) {

			Minutia minutia = this->minutiae.at(i);

			switch (minutia.getType()) {
			case ENDING_MINUTIA_TYPE:
				minutiaData[0] = 64;
				break;
			case BIFURCATION_MINUTIA_TYPE:
				minutiaData[0] = 128;
				break;
			default:
				minutiaData[0] = 0;
				break;
			}

			int x, y;
			double theta;
			x = (int)THIMBLE_ROUND(minutia.x);
			y = (int)THIMBLE_ROUND(minutia.y);
			theta = minutia.angle;

			if ( x < 0 || y < 0 || x >= (1<<14) || y >= (1<<14) ) {
				return false;
			}

			minutiaData[0] += (unsigned char) ((x >> 8) & 63);
			minutiaData[1] = (unsigned char) (x & 0xFF);
			minutiaData[2] = (unsigned char) ((y >> 8) & 63);
			minutiaData[3] = (unsigned char) (y & 0xFF);

			minutiaData[4] = (unsigned char) THIMBLE_ROUND(theta / (M_PI + M_PI)
					* 256.0);
			minutiaData[5] = (unsigned char) minutia.quality;

			if ( fwrite(minutiaData,1,6,out) != 6 ) {
				return false;
			}
		}

		return true;
	}

	/**
	 * @brief Reads a minutiae view's data from the specified
	 *        <code>FILE</code> and stores the result in this view.
	 *
	 * @details
	 *       see 'MinutiaeRecord.h'
	 */
	bool MinutiaeView::read( FILE *in ) {

		unsigned char header[4];

		if ( fread(header,sizeof(unsigned char),4,in) != 4 )
			return false;

		int impressionType = (header[1]>>4)&0xF;

		if ( impressionType != 0 && impressionType != 1 &&
			 impressionType != 2 && impressionType != 3 &&
			 impressionType != 4 && impressionType != 8 ) {
			return false;
		}

		if ( header[0] > 10 || header[2] > 100 || impressionType ) {
			return false;
		}

		this->fingerPosition = (FINGER_POSITION_T)header[0];
		this->viewNumber = header[1]&0xF;
		this->impressionType = (FINGER_IMPRESSION_TYPE_T)impressionType;
		this->fingerQuality = header[2];

		int numMinutiae = header[3];

		this->minutiae.assign(numMinutiae,Minutia());

		unsigned char minutiaData[6];
		for ( int i = 0 ; i < numMinutiae ; i++ ) {

			if ( fread(minutiaData,sizeof(unsigned char),6,in) != 6 ) {
				return false;
			}

			switch( (minutiaData[0]>>6) & 0x3 ) {
			case 1:
				this->minutiae[i].typ = ENDING_MINUTIA_TYPE;
				break;
			case 2:
				this->minutiae[i].typ = BIFURCATION_MINUTIA_TYPE;
				break;
			default:
				this->minutiae[i].typ = UNKNOWN_MINUTIA_TYPE;
				break;
			}

			uint16_t x , y;
			x = 0x100*(minutiaData[0]&((1<<6)-1)) + minutiaData[1];
			y = 0x100*(minutiaData[2]&((1<<6)-1)) + minutiaData[3];

			double angle = (M_PI+M_PI) * ( (double)minutiaData[4] / 256.0 );

			uint8_t quality = minutiaData[5];

			this->minutiae[i].x = x;
			this->minutiae[i].y = y;
			this->minutiae[i].angle = angle;
			this->minutiae[i].quality = quality;
		}

		return true;
	}

	/**
	 * @brief Reads a minutiae view from the specified data
	 * and overwrites the data of this view by the read data.
	 *
	 * @details
	 *       see 'MinutiaeRecord.h'
	 */
	int MinutiaeView::read( const void *data ) {

		const unsigned char *dat = (unsigned char*)data;

		unsigned char header[4];
		memcpy(header,dat,4*sizeof(unsigned char));

		int impressionType = (header[1]>>4)&0xF;

		if ( impressionType != 0 && impressionType != 1 &&
			 impressionType != 2 && impressionType != 3 &&
			 impressionType != 4 && impressionType != 8 ) {
			return false;
		}

		if ( header[0] > 10 || header[2] > 100 || impressionType ) {
			return false;
		}

		this->fingerPosition = (FINGER_POSITION_T)header[0];
		this->viewNumber = header[1]&0xF;
		this->impressionType = (FINGER_IMPRESSION_TYPE_T)impressionType;
		this->fingerQuality = header[2];

		int numMinutiae = header[3];

		this->minutiae.assign(numMinutiae,Minutia());

		dat += 4;

		unsigned char minutiaData[6];
		for ( int i = 0 ; i < numMinutiae ; i++ , dat += 6 ) {

			memcpy(minutiaData,dat,6*sizeof(unsigned char));

			switch( (minutiaData[0]>>6) & 0x3 ) {
			case 1:
				this->minutiae[i].typ = ENDING_MINUTIA_TYPE;
				break;
			case 2:
				this->minutiae[i].typ = BIFURCATION_MINUTIA_TYPE;
				break;
			default:
				this->minutiae[i].typ = UNKNOWN_MINUTIA_TYPE;
				break;
			}

			uint16_t x , y;
			x = 0x100*(minutiaData[0]&((1<<6)-1)) + minutiaData[1];
			y = 0x100*(minutiaData[2]&((1<<6)-1)) + minutiaData[3];

			double angle = (M_PI+M_PI) * ( (double)minutiaData[4] / 256.0 );

			uint8_t quality = minutiaData[5];

			this->minutiae[i].x = x;
			this->minutiae[i].y = y;
			this->minutiae[i].angle = angle;
			this->minutiae[i].quality = quality;
		}

		return true;
	}

	/**
	 * @brief Prints a text representation of the given minutiae view to an
	 *        output stream.
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	ostream & operator<<( ostream & out , const MinutiaeView & view ) {

		int n = view.getMinutiaeCount();

		out << n;

		for ( int k = 0 ; k < n ; k++ ) {
			out << endl << view.getMinutia(k);
		}

		return out;
	}

	/**
	 * @brief Initializes this record with default values
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 *
	 */
	void MinutiaeRecord::init() {

		memset(this->captureEquipmentCertifications,0,4*sizeof(bool));
		memset(this->captureDeviceTypeID,0,12*sizeof(bool));

		this->sizeX = 0;
		this->sizeY = 0;
		this->resX = 197;
		this->resY = 197;
		this->reservedByte = 0;

		this->views.clear();
	}

	/**
	 * @brief
	 *          Standard constructor.
	 *
	 * @details
	 *          Creates an empty minutiae record for fingerprints
	 *          of zero width and height corresponding to a resolution
	 *          of 197 pixels per centimeter (i.e., 500 dots per inch).
	 *
	 *          After an empty minutiae record has been created its
	 *          parameters its content should be initialized from an
	 *          external file or data in ISO 19794-2:2005 format using one
	 *          of its read members.
	 *
	 * @see read(FILE*)
	 * @see read(const void*)
	 * @see read(const std::string&)
	 */
	MinutiaeRecord::MinutiaeRecord() {

		init();
	}

	/**
	 * @brief Constructs an empty minutiae record with specified width,
	 *        height, horizontal-, and vertical resolution.
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	MinutiaeRecord::MinutiaeRecord
	(int width, int height, int ppcX, int ppcY) {

		init();

		if (width <= 0 || width > 0xFFFF ) {
			cerr << "MinutiaeRecord: "
				"The width must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (height <= 0) {
			cerr << "MinutiaeRecord: "
				"The height must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (ppcX <= 0) {
			cerr << "MinutiaeRecord: "
				"The horizontal resolution must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (ppcY <= 0) {
			cerr << "MinutiaeRecord: "
				"The vertical resolution must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		this->sizeX = width;
		this->sizeY = height;
		this->resX = ppcX;
		this->resY = ppcY;
	}

	/**
	 * @brief Constructs an empty minutiae record with specified width,
	 *        height, and resolution.
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	MinutiaeRecord::MinutiaeRecord(int width, int height, int ppc) {

		init();

		if (width <= 0 || width > 0xFFFF ) {
			cerr << "MinutiaeRecord: "
				"The width must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (height <= 0) {
			cerr << "MinutiaeRecord: "
				"The height must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (ppc <= 0) {
			cerr << "MinutiaeRecord: "
				"The resolution must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		this->sizeX = width;
		this->sizeY = height;
		this->resX = ppc;
		this->resY = ppc;
	}

	/**
	 * @brief Constructs an empty minutiae record with specified width and
	 *        height.
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	MinutiaeRecord::MinutiaeRecord(int width, int height) {

		init();

		if (width <= 0 || width > 0xFFFF ) {
			cerr << "MinutiaeRecord: "
				"The width must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		if (height <= 0) {
			cerr << "MinutiaeRecord: "
				"The height must be greater than zero and "
				"be representable by two bytes" << endl;
			exit(EXIT_FAILURE);
		}

		this->sizeX = width;
		this->sizeY = height;
		this->resX = 197;
		this->resY = 197;
	}

	/**
	 * @brief
	 *             Access the view in the record of given view number.
	 *
	 * @details
	 *             see 'MinutiaeRecord.h'
	 */
	MinutiaeView & MinutiaeRecord::getView( int viewNumber ) {

		int k;

		for ( k = 0; k < (int)(this->views.size()) ; k++ ) {
			if ( this->views[k].viewNumber == viewNumber ) {
				break;
			}
		}

		if ( k >= (int)(this->views.size()) ) {
			cerr << "MinutiaeRecord: "
				"No view with specified view number was found in the record"
				 << endl;
			exit(EXIT_FAILURE);
		}

		return this->views[k];
	}

	/**
	 * @brief
	 *        Access the view in the record of given view number
	 *        (constant version).
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	const MinutiaeView & MinutiaeRecord::getView( int viewNumber ) const {

		int k;

		for ( k = 0; k < (int)(this->views.size()) ; k++ ) {
			if ( this->views[k].viewNumber == viewNumber ) {
				break;
			}
		}

		if ( k >= (int)(this->views.size()) ) {
			cerr << "MinutiaeRecord: "
				"No view with specified view number was found in the record"
				 << endl;
			exit(EXIT_FAILURE);
		}

		return this->views[k];
	}

	/**
	 * @brief Sets this record to a hard copy of the given
	 *        <code>record</code>.
	 *
	 * @details
	 *        see 'MinutiaeRecord.h'
	 */
	MinutiaeRecord & MinutiaeRecord::operator=( const MinutiaeRecord & record ) {

		memcpy
		 (this->captureEquipmentCertifications,
		  record.captureEquipmentCertifications,
		  4*sizeof(bool));

		memcpy
		 ( this->captureDeviceTypeID,record.captureDeviceTypeID,12*sizeof(bool));

		this->sizeX = record.sizeX;
		this->sizeY = record.sizeY;
		this->resX  = record.resX;
		this->resY  = record.resY;
		this->reservedByte = record.reservedByte;

		this->views = record.views;

		return *this;
	}

	/**
	 * @brief
	 *            Determine the number of bytes required to store
	 *            the record.
	 *
	 * @details
	 *            see 'MinutiaeRecord.h'
	 */
	int MinutiaeRecord::getSizeInBytes() const {
		int size = 24;
		for (int i = 0; i < (int)(this->views.size()) ; i++) {
			size += this->views.at(i).getSizeInBytes();
		}
		return size;
	}

	/**
	 * @brief
	 *            Writes the minutiae record to the specified binary
	 *            output stream.
	 *
	 * @details
	 *            see 'MinutiaeRecord.h'
	 */
	bool MinutiaeRecord::write( FILE *out ) const {

		unsigned char header[24];

		// First four bytes correspond to format identifier
		header[0] = 'F';
		header[1] = 'M';
		header[2] = 'R';
		header[3] = '\0';

		// Next four bytes correspond to version number
		header[4] = ' ';
		header[5] = '2';
		header[6] = '0';
		header[7] = '\0';

		// Next four bytes correspond to length of template in bytes
		int size = getSizeInBytes();
		header[8] = (unsigned char) ((size >> 24) & 0xFF);
		header[9] = (unsigned char) ((size >> 16) & 0xFF);
		header[10] = (unsigned char) ((size >> 8) & 0xFF);
		header[11] = (unsigned char) ((size >> 0) & 0xFF);

		header[12] = 0;
		header[13] = 0;

		if (this->captureEquipmentCertifications[0]) {
			header[12] += (unsigned char) 0x1;
		}
		if (this->captureEquipmentCertifications[1]) {
			header[12] += (unsigned char) 0x2;
		}
		if (this->captureEquipmentCertifications[2]) {
			header[12] += (unsigned char) 0x4;
		}
		if (this->captureEquipmentCertifications[3]) {
			header[12] += (unsigned char) 0x8;
		}

		if (this->captureDeviceTypeID[0]) {
			header[12] += (unsigned char) 0x10;
		}
		if (this->captureDeviceTypeID[1]) {
			header[12] += (unsigned char) 0x20;
		}
		if (this->captureDeviceTypeID[2]) {
			header[12] += (unsigned char) 0x40;
		}
		if (this->captureDeviceTypeID[3]) {
			header[12] += (unsigned char) 0x80;
		}
		if (this->captureDeviceTypeID[4]) {
			header[13] += (unsigned char) 0x1;
		}
		if (this->captureDeviceTypeID[5]) {
			header[13] += (unsigned char) 0x2;
		}
		if (this->captureDeviceTypeID[6]) {
			header[13] += (unsigned char) 0x4;
		}
		if (this->captureDeviceTypeID[7]) {
			header[13] += (unsigned char) 0x8;
		}
		if (this->captureDeviceTypeID[8]) {
			header[13] += (unsigned char) 0x10;
		}
		if (this->captureDeviceTypeID[9]) {
			header[13] += (unsigned char) 0x20;
		}
		if (this->captureDeviceTypeID[10]) {
			header[13] += (unsigned char) 0x40;
		}
		if (this->captureDeviceTypeID[11]) {
			header[13] += (unsigned char) 0x80;
		}

		header[14] = (unsigned char) (this->sizeX >> 8);
		header[15] = (unsigned char) (this->sizeX);
		header[16] = (unsigned char) (this->sizeY >> 8);
		header[17] = (unsigned char) (this->sizeY);
		header[18] = (unsigned char) (this->resX >> 8);
		header[19] = (unsigned char) (this->resX);
		header[20] = (unsigned char) (this->resY >> 8);
		header[21] = (unsigned char) (this->resY);
		header[22] = (unsigned char) (this->views.size());
		header[23] = (unsigned char) (this->reservedByte);

		if ( fwrite(header,1,24,out) != 24 ) {
			return false;
		}

		for (int i = 0; i < (int)(this->views.size()) ; i++) {
			if ( !this->views.at(i).write(out) ) {
				return false;
			}
		}

		return true;
	}

	bool MinutiaeRecord::write( const std::string & fileName ) const {

		FILE *out;

		if ( (out=THIMBLE_FOPEN(fileName.c_str(),"wb")) == NULL ) {
			return false;
		}

		MinutiaeRecord record;

		bool state = write(out);

		fclose(out);

		return state;
	}

	/**
	 * @brief    Reads a minutiae record from the specified
	 *           <code>FILE</code> and returns the result.
	 *
	 * @details
	 *          see 'MinutiaeRecord.h'
	 */
	bool MinutiaeRecord::read( FILE *in ) {

		MinutiaeRecord record;

		//Read the header
		unsigned char header[24];
		if ( fread(header,sizeof(unsigned char),24,in) != 24 ) {
			return false;
		}

		// First four bytes correspond to format identifier
		if ( header[0] != 'F' || header[1] != 'M' || header[2] != 'R' ||
			 header[3] != '\0' ) {
			return false;
		}

		// Next four bytes correspond to version number
		if ( header[4] != ' ' || header[5] != '2' || header[6] != '0' ||
			 header[7] != '\0' ) {
			return false;
		}

		//int size = header[11]+0x100*(header[10]+0x100*(header[9]+0x100*header[8]));

		record.captureEquipmentCertifications[0] = (header[12]&0x1?true:false);
		record.captureEquipmentCertifications[1] = (header[12]&0x2?true:false);
		record.captureEquipmentCertifications[2] = (header[12]&0x4?true:false);
		record.captureEquipmentCertifications[3] = (header[12]&0x8?true:false);

		record.captureDeviceTypeID[0]  = (header[12]&0x10?true:false);
		record.captureDeviceTypeID[1]  = (header[12]&0x20?true:false);
		record.captureDeviceTypeID[2]  = (header[12]&0x40?true:false);
		record.captureDeviceTypeID[3]  = (header[12]&0x80?true:false);
		record.captureDeviceTypeID[4]  = (header[13]&0x1?true:false);
		record.captureDeviceTypeID[5]  = (header[13]&0x2?true:false);
		record.captureDeviceTypeID[6]  = (header[13]&0x4?true:false);
		record.captureDeviceTypeID[7]  = (header[13]&0x8?true:false);
		record.captureDeviceTypeID[8]  = (header[13]&0x10?true:false);
		record.captureDeviceTypeID[9]  = (header[13]&0x20?true:false);
		record.captureDeviceTypeID[10] = (header[13]&0x40?true:false);
		record.captureDeviceTypeID[11] = (header[13]&0x80?true:false);

		record.sizeX = (int)header[15]+0x100*(int)header[14];
		record.sizeY = (int)header[17]+0x100*(int)header[16];
		record.resX = (int)header[19]+0x100*(int)header[18];
		record.resY = (int)header[21]+0x100*(int)header[20];
		record.reservedByte = header[23];

		int numViews = header[22];

		record.views.assign(numViews,MinutiaeView());
		for ( int i = 0 ; i < numViews ; i++ ) {
			if ( !record.views[i].read(in) ) {
				return false;
			}
		}

		*this = record;

		return true;
	}

	/**
	 * @brief    Reads a minutiae record from the specified
	 *           data array.
	 *
	 * @param data
	 *           Contains the data of the specified minutiae record
	 *
	 * @return
	 *          <code>true</code> if the minutiae record has been
	 *          successfully read from <code>data</code>; otherwise, the
	 *          function returns <code>false</code>.
	 */
	bool MinutiaeRecord::fromBytes( const void *data ) {

		const unsigned char *dat = (unsigned char*)data;

		MinutiaeRecord record;

		//Read the header
		unsigned char header[24];
		header[0] = (unsigned char)dat[0];
		header[1] = (unsigned char)dat[1];
		header[2] = (unsigned char)dat[2];
		header[3] = (unsigned char)dat[3];

		memcpy(header,dat,24*sizeof(unsigned char));

		// First four bytes correspond to format identifier
		if ( header[0] != 'F' || header[1] != 'M' || header[2] != 'R' ||
			 header[3] != '\0' ) {
			return false;
		}

		// Next four bytes correspond to version number
		if ( header[4] != ' ' || header[5] != '2' || header[6] != '0' ||
			 header[7] != '\0' ) {
			return false;
		}

		//int size = header[11]+0x100*(header[10]+0x100*(header[9]+0x100*header[8]));

		record.captureEquipmentCertifications[0] = (header[12]&0x1?true:false);
		record.captureEquipmentCertifications[1] = (header[12]&0x2?true:false);
		record.captureEquipmentCertifications[2] = (header[12]&0x4?true:false);
		record.captureEquipmentCertifications[3] = (header[12]&0x8?true:false);

		record.captureDeviceTypeID[0]  = (header[12]&0x10?true:false);
		record.captureDeviceTypeID[1]  = (header[12]&0x20?true:false);
		record.captureDeviceTypeID[2]  = (header[12]&0x40?true:false);
		record.captureDeviceTypeID[3]  = (header[12]&0x80?true:false);
		record.captureDeviceTypeID[4]  = (header[13]&0x1?true:false);
		record.captureDeviceTypeID[5]  = (header[13]&0x2?true:false);
		record.captureDeviceTypeID[6]  = (header[13]&0x4?true:false);
		record.captureDeviceTypeID[7]  = (header[13]&0x8?true:false);
		record.captureDeviceTypeID[8]  = (header[13]&0x10?true:false);
		record.captureDeviceTypeID[9]  = (header[13]&0x20?true:false);
		record.captureDeviceTypeID[10] = (header[13]&0x40?true:false);
		record.captureDeviceTypeID[11] = (header[13]&0x80?true:false);

		record.sizeX = (int)header[15]+0x100*(int)header[14];
		record.sizeY = (int)header[17]+0x100*(int)header[16];
		record.resX = (int)header[19]+0x100*(int)header[18];
		record.resY = (int)header[21]+0x100*(int)header[20];
		record.reservedByte = header[23];

		int numViews = header[22];

		dat += 24;

		record.views.assign(numViews,MinutiaeView());
		for ( int i = 0 ; i < numViews ; i++ ) {

			int state = record.views[i].read((void*)dat);
			if ( state < 0 ) {
				return false;
			}

			dat += state;

		}

		*this = record;

		return true;
	}

	/**
	 * @brief    Reads a minutiae record from a file specified
	 *           by its path and returns the result.
	 *
	 * @details
     *           see 'MinutiaeRecord.h'
	 */
	bool MinutiaeRecord::read( const string & fileName ) {

		FILE *in;
		
		if ( (in=THIMBLE_FOPEN(fileName.c_str(),"rb")) == NULL ) {
			return false;
		}

		MinutiaeRecord record;
		if ( !record.read(in) ) {
			return false;
		}

		*this = record;

		fclose(in);

		return true;
	}
}



