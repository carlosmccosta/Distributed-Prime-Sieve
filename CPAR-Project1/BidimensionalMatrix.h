#pragma once

#include <windows.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

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

	double getValue(unsigned int matrixColumn, unsigned int matrixLine);
	void putValue(unsigned int matrixColumn, unsigned int matrixLine, double value);
	void addValue(unsigned int matrixColumn, unsigned int matrixLine, double value);
	void multiplyValue(unsigned int matrixColumn, unsigned int matrixLine, double value);

	bool validateResultOfDefaultMatrixInitialization(unsigned int numberLinesRightMatrix); //sum of powers
	bool validateResultFromFile(string expectedResultMatrixFilename);
	bool exportMatrixToFile(string filename);

	unsigned int getNumberColumns() const { return numberColumns; }
	unsigned int getNumberLines() const { return numberLines; }

private:
	unsigned int numberLines, numberColumns;
	double* matrixData;
	HGLOBAL memoryHandle;
};
