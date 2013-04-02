#include "MatrixMultiplicationLine.h"

MatrixMultiplicationLine::MatrixMultiplicationLine() {}
MatrixMultiplicationLine::~MatrixMultiplicationLine() {}

shared_ptr<BidimensionalMatrix> MatrixMultiplicationLine::performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix ) {
	unsigned int numberColumnsLeftMatrix = leftMatrix.getNumberColumns();
	unsigned int numberLinesLeftMatrix = leftMatrix.getNumberLines();
	unsigned int numberColumnsRightMatrix = rightMatrix.getNumberColumns();

	shared_ptr<BidimensionalMatrix> resultMatrix = make_shared<BidimensionalMatrix>(numberColumnsRightMatrix, numberLinesLeftMatrix);
	resultMatrix->allocateMemoryForMatrixData();

	unsigned int columnOfLeftMatrix, lineOfLeftMatrix, columnOfRightMatrix;
	double parcialResultMatrixCellValue;

	performanceTimer.reset();
	performanceTimer.start();

	for (lineOfLeftMatrix = 0; lineOfLeftMatrix < numberLinesLeftMatrix; ++lineOfLeftMatrix) {
		for (columnOfLeftMatrix = 0; columnOfLeftMatrix < numberColumnsLeftMatrix; ++columnOfLeftMatrix) {
			double firstMatrixLineValue = leftMatrix.getValue(columnOfLeftMatrix, lineOfLeftMatrix);
			for (columnOfRightMatrix = 0; columnOfRightMatrix < numberColumnsRightMatrix; ++columnOfRightMatrix) {
				parcialResultMatrixCellValue = firstMatrixLineValue * rightMatrix.getValue(columnOfRightMatrix, columnOfLeftMatrix);
				resultMatrix->addValue(columnOfRightMatrix, lineOfLeftMatrix, parcialResultMatrixCellValue);
			}
		}
	}

	performanceTimer.stop();

	return resultMatrix;
}
