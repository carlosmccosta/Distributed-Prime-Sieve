#pragma once

#include "MatrixMultiplication.h"
#include "BidimensionalMatrix.h"

using std::make_shared;

class MatrixMultiplicationLine : public MatrixMultiplication {
public:
	MatrixMultiplicationLine();
	virtual ~MatrixMultiplicationLine();

	virtual shared_ptr<BidimensionalMatrix> performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix );
};
