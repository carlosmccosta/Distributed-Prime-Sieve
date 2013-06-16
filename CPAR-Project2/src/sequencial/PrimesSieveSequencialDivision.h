#pragma once

#include "../PrimesSieve.h"

#include <cmath>

using std::sqrt;

template<typename FlagsContainer>
class PrimesSieveSequencialDivision: public PrimesSieve<FlagsContainer> {
	public:
		PrimesSieveSequencialDivision() {
		}
		
		virtual ~PrimesSieveSequencialDivision() {
		}
		
		void computePrimes(size_t maxRange) {
			this->template primesValues.clear();
			this->template initPrimesBitset(maxRange);
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			
			this->template performanceTimer.reset();
			this->template performanceTimer.start();

			for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= maxRange; compositeNumber += 2) {
						if ((compositeNumber % primeNumber) == 0) {
							this->template setPrimesBitsetValue(compositeNumber, false);
						}
					}
				}
			}
			
			this->template performanceTimer.stop();
		}
	};

