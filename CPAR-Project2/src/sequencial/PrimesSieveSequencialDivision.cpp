/*
 * PrimesSieveSequencialDivision.cpp
 *
 *  Created on: May 24, 2013
 *      Author: carloscosta
 */

#include "PrimesSieveSequencialDivision.h"


vector<int>& PrimesSieveSequencialDivision::computePrimes(size_t maxRange) {
	primesBitset = vector<bool>(getNumberBitsToStore(maxRange));




	return primesValues;
}
