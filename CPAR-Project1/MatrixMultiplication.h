#pragma once

#include "BidimensionalMatrix.h"
#include "PerformanceTimer.h"
#include <memory>

using std::shared_ptr;

class MatrixMultiplication {
public:
	MatrixMultiplication();
	virtual ~MatrixMultiplication();
	virtual shared_ptr<BidimensionalMatrix> performMultiplication(BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix, unsigned int blockColumnSize = 1, unsigned int blockLineSize = 1) = 0;
	
	PerformanceTimer getPerformanceTimer() const { return performanceTimer; }

protected:
	PerformanceTimer performanceTimer;
};
