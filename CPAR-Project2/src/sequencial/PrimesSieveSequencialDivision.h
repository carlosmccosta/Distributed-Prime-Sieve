#pragma once

#include "../PrimesSieve.h"

class PrimesSieveSequencialDivision: public PrimesSieve {
	public:
		PrimesSieveSequencialDivision() {}
		virtual ~PrimesSieveSequencialDivision() {}

		vector<int>& computePrimes(size_t maxRange);
};

