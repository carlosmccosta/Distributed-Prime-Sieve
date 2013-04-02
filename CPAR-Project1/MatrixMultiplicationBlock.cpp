#include "MatrixMultiplicationBlock.h"

MatrixMultiplicationBlock::MatrixMultiplicationBlock() {}
MatrixMultiplicationBlock::~MatrixMultiplicationBlock() {}

shared_ptr<BidimensionalMatrix> MatrixMultiplicationBlock::performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix, unsigned int blockColumnSize, unsigned int blockLineSize ) {
	unsigned int numberColumnsLeftMatrix = leftMatrix.getNumberColumns();
	unsigned int numberLinesLeftMatrix = leftMatrix.getNumberLines();
	unsigned int numberColumnsRightMatrix = rightMatrix.getNumberColumns();

	shared_ptr<BidimensionalMatrix> resultMatrix = make_shared<BidimensionalMatrix>(numberColumnsRightMatrix, numberLinesLeftMatrix);
	BidimensionalMatrix* resultMatrixPtr = resultMatrix.get();
	resultMatrixPtr->allocateMemoryForMatrixData();

	unsigned int columnOfLeftMatrix, lineOfLeftMatrix, columnOfRightMatrix;
	unsigned int beginBlockColumnOfLeftMatrix, beginBlockLineOfLeftMatrix, beginBlockColumnOfRightMatrix;
	unsigned int endBlockColumnOfLeftMatrix = blockColumnSize, endBlockLineOfLeftMatrix = blockLineSize, endBlockColumnOfRightMatrix = blockColumnSize;
	double parcialResultMatrixCellValue;

	performanceTimer.reset();
	performanceTimer.start();

	for (beginBlockLineOfLeftMatrix = 0;
		beginBlockLineOfLeftMatrix < numberLinesLeftMatrix;
		beginBlockLineOfLeftMatrix += blockLineSize) {
		for (beginBlockColumnOfLeftMatrix = 0;
			beginBlockColumnOfLeftMatrix < numberColumnsLeftMatrix;
			beginBlockColumnOfLeftMatrix += blockColumnSize) {
			for (beginBlockColumnOfRightMatrix = 0;
				beginBlockColumnOfRightMatrix < numberColumnsRightMatrix;
				beginBlockColumnOfRightMatrix += blockColumnSize) {

				for (lineOfLeftMatrix = beginBlockLineOfLeftMatrix;
					lineOfLeftMatrix < endBlockLineOfLeftMatrix;
					++lineOfLeftMatrix) {
					for (columnOfLeftMatrix = beginBlockColumnOfLeftMatrix;
						columnOfLeftMatrix < endBlockColumnOfLeftMatrix;
						++columnOfLeftMatrix) {
						double firstMatrixLineValue = leftMatrix.getValue(columnOfLeftMatrix, lineOfLeftMatrix);
						for (columnOfRightMatrix = beginBlockColumnOfRightMatrix;
							columnOfRightMatrix < endBlockColumnOfRightMatrix;
							++columnOfRightMatrix) {
							parcialResultMatrixCellValue = firstMatrixLineValue * rightMatrix.getValue(columnOfRightMatrix, columnOfLeftMatrix);
							resultMatrixPtr->addValue(columnOfRightMatrix, lineOfLeftMatrix, parcialResultMatrixCellValue);
						}
					}
				}
				endBlockLineOfLeftMatrix = min(numberLinesLeftMatrix, endBlockLineOfLeftMatrix + blockLineSize);
				endBlockColumnOfLeftMatrix = min(numberColumnsLeftMatrix, endBlockColumnOfLeftMatrix + blockColumnSize);
				endBlockColumnOfRightMatrix = min(numberColumnsRightMatrix, endBlockColumnOfRightMatrix + blockColumnSize);

			}
		}
	}


	performanceTimer.stop();

	return resultMatrix;
}
