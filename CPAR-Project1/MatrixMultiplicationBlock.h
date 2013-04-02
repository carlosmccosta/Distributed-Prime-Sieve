#pragma once

#include "MatrixMultiplication.h"
#include "BidimensionalMatrix.h"

using std::make_shared;

class MatrixMultiplicationBlock : public MatrixMultiplication {
public:
	MatrixMultiplicationBlock();
	virtual ~MatrixMultiplicationBlock();

	virtual shared_ptr<BidimensionalMatrix> performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix, unsigned int blockColumnSize = 1, unsigned int blockLineSize = 1 );
};
