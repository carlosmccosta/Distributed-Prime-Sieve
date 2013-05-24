#include "PrimesSieveSequencialDivision.h"

vector<size_t>& PrimesSieveSequencialDivision::computePrimes(size_t maxRange) {
	initPrimesCompositeBitset(maxRange);
	size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
	
	performanceTimer.reset();
	performanceTimer.start();
	
	for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
		// for each number not marked as composite (prime number)
		if (!getPrimesBitsetValue(primeNumber)) {
			//use it to calculate his composites
			for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= maxRange; compositeNumber += 2) {
				if ((compositeNumber % primeNumber) == 0) {
					setPrimesBitsetValue(compositeNumber, true);
				}
			}
		}
	}
	
	performanceTimer.stop();
	
	return extractPrimesFromBitset();
}
