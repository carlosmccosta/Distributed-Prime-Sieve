#include "BidimensionalMatrix.h"

BidimensionalMatrix::BidimensionalMatrix() :
	numberColumns(2), numberLines(2), matrixData(NULL), memoryHandle(NULL) {}

BidimensionalMatrix::BidimensionalMatrix( unsigned int numberOfColumns, unsigned int numberOfLines ) :
	numberColumns(numberOfColumns), numberLines(numberOfLines), matrixData(NULL), memoryHandle(NULL) {}


BidimensionalMatrix::~BidimensionalMatrix()
{
	releaseMemoryOfMatrixData();
}

bool BidimensionalMatrix::initializeMatrix( double defaultValue )
{
	if (allocateMemoryForMatrixData()) {
		for (unsigned int line = 0; line < numberLines; ++line) {
			for (unsigned int column = 0; column < numberColumns; ++column) {
				putValue(column, line, defaultValue);
			}
		}
	} else {
		return false;
	}

	return true;
}

bool BidimensionalMatrix::initializeMatrixFromFile( string initializationFile ) {
	return true;
}

bool BidimensionalMatrix::initializeMatrixWithSequenceOnLines() {
	if (allocateMemoryForMatrixData()) {
		for (unsigned int line = 0; line < numberLines; ++line) {
			double lineValue = line + 1;
			for (unsigned int column = 0; column < numberColumns; ++column) {
				putValue(column, line, lineValue);
			}
		}
	} else {
		return false;
	}

	return true;
}

void BidimensionalMatrix::putValue( unsigned int matrixColum, unsigned int matrixLine, double value ) {
	matrixData[numberColumns*matrixLine + matrixColum] = value;
}

double BidimensionalMatrix::getValue( unsigned int matrixColumn, unsigned int matrixLine ) {
	return matrixData[numberColumns*matrixLine + matrixColumn];
}

bool BidimensionalMatrix::allocateMemoryForMatrixData() {
	memoryHandle = GlobalAlloc(GMEM_FIXED, numberColumns * numberLines * sizeof(double));
	if (memoryHandle == NULL) {
		return false;
	}

	matrixData = (double*)GlobalLock(memoryHandle);
	if (matrixData == NULL) {
		return false;
	}	

	return true;
}

bool BidimensionalMatrix::releaseMemoryOfMatrixData() {
	if (memoryHandle != NULL) {
		GlobalUnlock(memoryHandle);
		GlobalFree(memoryHandle);
	} else {
		return false;
	}

	return true;
}

bool BidimensionalMatrix::validateResultOfDefaultMatrixInitialization(unsigned int numberLinesRightMatrix) {
	double sumOfPowers = (numberLinesRightMatrix * (numberLinesRightMatrix + 1)) / 2;
	unsigned int columnPos, linePos;


	for (linePos = 0; linePos < numberLines; ++linePos) {
		for (columnPos = 0; columnPos < numberColumns; ++columnPos) {
			if (getValue(columnPos, linePos) != sumOfPowers) {
				return false;
			}
		}
	}

	return true;
}

bool BidimensionalMatrix::exportMatrixToFile( string filename ) {
	ofstream outputStream;
	outputStream.open(filename);

	if (outputStream.is_open()) {
		unsigned int columnPos, linePos;

		for (linePos = 0; linePos < numberLines; ++linePos) {
			for (columnPos = 0; columnPos < numberColumns; ++columnPos) {
				outputStream << getValue(columnPos, linePos);
				if (columnPos == numberColumns-1) {
					outputStream << endl;
				} else {
					outputStream << " ";
				}
			}
		}

		outputStream.close();
		return true;
	} else {
		return false;
	}
}
