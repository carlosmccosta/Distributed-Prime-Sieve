#include "MatrixMultiplicationBlock.h"

MatrixMultiplicationBlock::MatrixMultiplicationBlock() {}
MatrixMultiplicationBlock::~MatrixMultiplicationBlock() {}

shared_ptr<BidimensionalMatrix> MatrixMultiplicationBlock::performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix, unsigned int blockColumnSize, unsigned int blockLineSize ) {
	unsigned int numberColumnsLeftMatrix = leftMatrix.getNumberColumns();
	unsigned int numberLinesLeftMatrix = leftMatrix.getNumberLines();
	unsigned int numberColumnsRightMatrix = rightMatrix.getNumberColumns();
	unsigned int numberLinesRightMatrix = rightMatrix.getNumberLines();

	if (blockColumnSize > numberColumnsLeftMatrix || blockColumnSize > numberColumnsRightMatrix) {
		blockColumnSize = min(numberColumnsLeftMatrix, numberColumnsRightMatrix);
		cout << "\n\n    -> Maximum block column size exceeded! Changed to limit: " << blockColumnSize << "\n";
	}

	if (blockLineSize > numberLinesLeftMatrix || blockLineSize > numberLinesRightMatrix) {
		blockLineSize = min(numberLinesLeftMatrix, numberLinesRightMatrix);
		cout << "    -> Maximum block line size exceeded! Changed to limit: " << blockColumnSize << "\n\n";
	}

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
		beginBlockLineOfLeftMatrix = endBlockLineOfLeftMatrix, endBlockLineOfLeftMatrix = min(numberLinesLeftMatrix, endBlockLineOfLeftMatrix + blockLineSize)) {
			for (beginBlockColumnOfLeftMatrix = 0;
				beginBlockColumnOfLeftMatrix < numberColumnsLeftMatrix;
				beginBlockColumnOfLeftMatrix = endBlockColumnOfLeftMatrix, endBlockColumnOfLeftMatrix = min(numberColumnsLeftMatrix, endBlockColumnOfLeftMatrix + blockColumnSize)) {
					for (beginBlockColumnOfRightMatrix = 0;
						beginBlockColumnOfRightMatrix < numberColumnsRightMatrix;
						beginBlockColumnOfRightMatrix = endBlockColumnOfRightMatrix, endBlockColumnOfRightMatrix = min(numberColumnsRightMatrix, endBlockColumnOfRightMatrix + blockColumnSize)) {

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
					}	
			}
	}


	performanceTimer.stop();

	return resultMatrix;
}
