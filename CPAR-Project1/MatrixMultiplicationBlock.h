#pragma once

#include "MatrixMultiplication.h"

class MatrixMultiplicationBlock : public MatrixMultiplication {
public:
	MatrixMultiplicationBlock();
	virtual ~MatrixMultiplicationBlock();

	virtual shared_ptr<BidimensionalMatrix> performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix );
};
