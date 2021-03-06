#include "MatrixMultiplicationBasic.h"

MatrixMultiplicationBasic::MatrixMultiplicationBasic() {}
MatrixMultiplicationBasic::~MatrixMultiplicationBasic() {}

shared_ptr<BidimensionalMatrix> MatrixMultiplicationBasic::performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix, unsigned int blockColumnSize, unsigned int blockLineSize ) {
	unsigned int numberColumnsLeftMatrix = leftMatrix.getNumberColumns();
	unsigned int numberLinesLeftMatrix = leftMatrix.getNumberLines();
	unsigned int numberColumnsRightMatrix = rightMatrix.getNumberColumns();

	shared_ptr<BidimensionalMatrix> resultMatrix = make_shared<BidimensionalMatrix>(numberColumnsRightMatrix, numberLinesLeftMatrix);
	BidimensionalMatrix* resultMatrixPtr = resultMatrix.get();
	resultMatrixPtr->allocateMemoryForMatrixData();

	unsigned int lineOfLeftMatrix, columnOfRightMatrix, positionOnLineAndColumn;
	double resultMatrixCellValue;

	performanceTimer.reset();
	performanceTimer.start();

	for (lineOfLeftMatrix = 0; lineOfLeftMatrix < numberLinesLeftMatrix; ++lineOfLeftMatrix) {
		for (columnOfRightMatrix = 0; columnOfRightMatrix < numberColumnsRightMatrix; ++columnOfRightMatrix) {
			resultMatrixCellValue = 0;
			for (positionOnLineAndColumn = 0; positionOnLineAndColumn < numberColumnsLeftMatrix; ++positionOnLineAndColumn) {
				resultMatrixCellValue += (leftMatrix.getValue(positionOnLineAndColumn, lineOfLeftMatrix) * rightMatrix.getValue(columnOfRightMatrix, positionOnLineAndColumn));
			}
			resultMatrixPtr->putValue(columnOfRightMatrix, lineOfLeftMatrix, resultMatrixCellValue);
		}
	}

	performanceTimer.stop();

	return resultMatrix;
}
