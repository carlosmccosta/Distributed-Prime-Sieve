#pragma once

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

using std::cout;
using std::string;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::endl;

class BidimensionalMatrix {
public:
	BidimensionalMatrix();
	BidimensionalMatrix(unsigned int numberOfColumns, unsigned int numberOfLines);
	virtual ~BidimensionalMatrix();

	bool initializeMatrix(double defaultValue = 1.0);
	bool initializeMatrixFromFile(string initializationFile);
	bool initializeMatrixWithSequenceOnLines();

	bool allocateMemoryForMatrixData();
	bool releaseMemoryOfMatrixData();

	bool validateResultOfDefaultMatrixInitialization(unsigned int numberLinesRightMatrix); //sum of powers
	bool validateResultFromFile(string expectedResultMatrixFilename);
	bool exportMatrixToFile(string filename);


	inline double getValue(unsigned int matrixColumn, unsigned int matrixLine) {
		return matrixData[numberColumns * matrixLine + matrixColumn];
	}

	inline void putValue(unsigned int matrixColumn, unsigned int matrixLine, double value) {
		matrixData[numberColumns * matrixLine + matrixColumn] = value;
	}

	inline void addValue(unsigned int matrixColumn, unsigned int matrixLine, double value) {
		matrixData[numberColumns * matrixLine + matrixColumn] += value;
	}

	inline void multiplyValue(unsigned int matrixColumn, unsigned int matrixLine, double value) {
		matrixData[numberColumns * matrixLine + matrixColumn] *= value;
	}

	inline unsigned int getNumberColumns() const { return numberColumns; }
	inline unsigned int getNumberLines() const { return numberLines; }

private:
	unsigned int numberColumns, numberLines;
	double* matrixData;

#ifdef _WIN32
	HGLOBAL memoryHandle;
#endif

};
