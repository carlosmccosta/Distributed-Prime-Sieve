#pragma once
#include "BidimensionalMatrix.h"
#include "PerformanceTimer.h"

#include <memory>

using std::shared_ptr;

class MatrixMultiplication
{
public:
	MatrixMultiplication();
	virtual ~MatrixMultiplication();
	virtual shared_ptr<BidimensionalMatrix> performMultiplication(BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix) = 0;
	
	PerformanceTimer getPerformanceTimer() const { return performanceTimer; }

protected:
	PerformanceTimer performanceTimer;
};

