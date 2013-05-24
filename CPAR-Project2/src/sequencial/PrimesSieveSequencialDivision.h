#pragma once

#include "../PrimesSieve.h"

#include <cmath>

using std::sqrt;

class PrimesSieveSequencialDivision: public PrimesSieve {
	public:
		PrimesSieveSequencialDivision() {
		}
		virtual ~PrimesSieveSequencialDivision() {
		}
		
		void computePrimes(size_t maxRange);
};

