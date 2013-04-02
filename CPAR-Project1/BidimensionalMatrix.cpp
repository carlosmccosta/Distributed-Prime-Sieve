#include "BidimensionalMatrix.h"

BidimensionalMatrix::BidimensionalMatrix() :
	numberColumns(2), numberLines(2), matrixData(NULL), memoryHandle(NULL) {}

BidimensionalMatrix::BidimensionalMatrix( unsigned int numberOfColumns, unsigned int numberOfLines ) :
	numberColumns(numberOfColumns), numberLines(numberOfLines), matrixData(NULL), memoryHandle(NULL) {}


BidimensionalMatrix::~BidimensionalMatrix() {
	releaseMemoryOfMatrixData();
}

bool BidimensionalMatrix::initializeMatrix( double defaultValue ) {
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
	ifstream inputStream;
	inputStream.open(initializationFile);

	if (inputStream.is_open()) {
		string fileLine;
		getline(inputStream, fileLine);
		istringstream iSStrAttrs(fileLine);
		iSStrAttrs >> numberColumns >> numberLines;

		if (allocateMemoryForMatrixData()) {
			unsigned int columnPos = 0, linePos = 0;
			while (getline(inputStream, fileLine)) {
				istringstream iSStr(fileLine);
				double fileNumber;
				while (iSStr >> fileNumber) {
					putValue(columnPos, linePos, fileNumber);
					++columnPos;
				}
				++linePos;
				columnPos = 0;
			}

		}
		inputStream.close();
		return true;
	} else {
		cout << "    -> Loading of input file failed!\n\n";
		return false;
	}
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



bool BidimensionalMatrix::allocateMemoryForMatrixData() {
#ifdef _WIN32
	memoryHandle = GlobalAlloc(GMEM_FIXED, numberColumns * numberLines * sizeof(double));
	if (memoryHandle == NULL) {
		return false;
	}

	matrixData = (double*)GlobalLock(memoryHandle);
	if (matrixData == NULL) {
		GlobalFree(memoryHandle);
		return false;
	}	

	return true;
#else
	matrixData = (double*)malloc(numberColumns * numberLines * sizeof(double));
	if (matrixData == NULL) {
		return false;
	} else {
		return true;
	}

#endif
}

bool BidimensionalMatrix::releaseMemoryOfMatrixData() {
#ifdef _WIN32
	if (memoryHandle != NULL) {
		GlobalUnlock(memoryHandle);
		GlobalFree(memoryHandle);
	} else {
		return false;
	}

	return true;
#else
	free matrixData;
	return true;
#endif
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


bool BidimensionalMatrix::validateResultFromFile(string expectedResultMatrixFilename) {
	ifstream inputStream;
	inputStream.open(expectedResultMatrixFilename);

	if (inputStream.is_open()) {
		unsigned int columnPos, linePos;

		double newValFromFile;
		for (linePos = 0; linePos < numberLines; ++linePos) {
			for (columnPos = 0; columnPos < numberColumns; ++columnPos) {
				inputStream >> newValFromFile;

				if (newValFromFile != getValue(columnPos, linePos)) {
					return false;
				}
			}
		}

		inputStream.close();
		return true;
	} else {
		cout << "    -> Loading of result file failed!\n\n";
		return false;
	}
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
		cout << "    -> Export of file " << filename << " failed!\n\n";
	}
}

