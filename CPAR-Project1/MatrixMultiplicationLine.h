#pragma once
#include "MatrixMultiplication.h"
class MatrixMultiplicationLine :
	public MatrixMultiplication
{
public:
	MatrixMultiplicationLine();
	virtual ~MatrixMultiplicationLine();

	virtual shared_ptr<BidimensionalMatrix> performMultiplication( BidimensionalMatrix& leftMatrix, BidimensionalMatrix& rightMatrix );

};

