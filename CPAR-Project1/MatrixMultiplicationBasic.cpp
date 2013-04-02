#include "MatrixMultiplicationBasic.h"
#include "BidimensionalMatrix.h"


MatrixMultiplicationBasic::MatrixMultiplicationBasic() {}
MatrixMultiplicationBasic::~MatrixMultiplicationBasic() {}


shared_ptr<BidimensionalMatrix> MatrixMultiplicationBasic::performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix ) {
	unsigned int numberColumnsLeftMatrix = leftMatrix.getNumberColumns();
	unsigned int numberLinesLeftMatrix = leftMatrix.getNumberLines();
	unsigned int numberColumnsRightMatrix = rightMatrix.getNumberColumns();

	shared_ptr<BidimensionalMatrix> resultMatrix = make_shared<BidimensionalMatrix>(numberColumnsRightMatrix, numberLinesLeftMatrix);
	resultMatrix->allocateMemoryForMatrixData();

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
			resultMatrix->putValue(columnOfRightMatrix, lineOfLeftMatrix, resultMatrixCellValue);
		}
	}

	performanceTimer.stop();

	return resultMatrix;
}


