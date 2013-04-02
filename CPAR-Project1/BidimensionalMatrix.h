#pragma once

#include <windows.h>
#include <string>
#include <fstream>

using std::string;
using std::ofstream;
using std::endl;

class BidimensionalMatrix
{
public:
	BidimensionalMatrix();
	BidimensionalMatrix(unsigned int numberOfColumns, unsigned int numberOfLines);
	virtual ~BidimensionalMatrix();

	bool initializeMatrix(double defaultValue = 1.0);
	bool initializeMatrixFromFile(string initializationFile);
	bool initializeMatrixWithSequenceOnLines();

	bool allocateMemoryForMatrixData();
	bool releaseMemoryOfMatrixData();

	void putValue(unsigned int matrixColum, unsigned int matrixLine, double value);
	double getValue(unsigned int matrixColum, unsigned int matrixLine);

	bool validateResultOfDefaultMatrixInitialization(unsigned int numberLinesRightMatrix); //sum of powers
	bool exportMatrixToFile(string filename);

	unsigned int getNumberColumns() const { return numberColumns; }
	unsigned int getNumberLines() const { return numberLines; }
	

private:
	unsigned int numberLines, numberColumns;
	double* matrixData;
	HGLOBAL memoryHandle;
};

