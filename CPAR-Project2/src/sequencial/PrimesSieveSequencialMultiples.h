#pragma once

#include "../PrimesSieve.h"

#include <cmath>

using std::sqrt;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiples: public PrimesSieve<FlagsContainer> {
	public:
		PrimesSieveSequencialMultiples() {
		}
		
		virtual ~PrimesSieveSequencialMultiples() {
		}
		
		void computePrimes(size_t maxRange) {
			this->template initPrimesCompositeBitset(maxRange);
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			
			this->template performanceTimer.reset();
			this->template performanceTimer.start();

			for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (!(this->template getPrimesBitsetValue(primeNumber))) {
					//use it to calculate his composites
					for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= maxRange; compositeNumber += primeNumber) {
						this->template setPrimesBitsetValue(compositeNumber, true);
					}
				}
			}
			
			this->template performanceTimer.stop();
			this->template primesValues.clear();
		}
	};

