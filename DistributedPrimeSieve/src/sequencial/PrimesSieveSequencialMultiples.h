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
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template getPrimesValues().clear();
			this->template initPrimesBitset(maxRange);
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (!this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= maxRange; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, true);
					}
				}
			}

			this->template getPerformanceTimer().stop();
		}
	};

