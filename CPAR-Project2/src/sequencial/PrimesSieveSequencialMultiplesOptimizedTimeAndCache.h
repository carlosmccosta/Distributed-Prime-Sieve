#pragma once

#include "../PrimesSieve.h"

#include <cmath>

using std::sqrt;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedTimeAndCache: public PrimesSieve<FlagsContainer> {
	public:
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache() {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedTimeAndCache() {
		}
		
		void computePrimes(size_t maxRange) {
			this->template initPrimesBitset(maxRange);
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			
			this->template performanceTimer.reset();
			this->template performanceTimer.start();



			for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= maxRange; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}
				}
			}
			
			this->template performanceTimer.stop();
			this->template primesValues.clear();
		}
	};

