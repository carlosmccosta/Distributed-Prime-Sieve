#pragma once

#include "PrimesSieveSequencialMultiplesOptimized.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>

using std::min;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedSpaceAndCache: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	public:
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache(size_t blockSizeInBytes = 64 * 1024) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}

		virtual ~PrimesSieveSequencialMultiplesOptimizedSpaceAndCache() {
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			for (size_t primesIndex = 0; primesIndex < blockIndexBegin; ++primesIndex) {
				if (this->template getPrimesBitset()[primesIndex]) {
					size_t primeNumber = this->template getNumberAssociatedWithBitsetPosition(primesIndex);
					size_t primeDoubled = primeNumber << 1;
					size_t compositeNumber = PrimesUtils::closestPrimeMultiple(primeNumber, blockBeginNumber);
					while (compositeNumber < blockBeginNumber || (compositeNumber % 2 == 0)) {
						compositeNumber += primeNumber;
					}

					for (; compositeNumber < blockEndNumber; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}
				}
			}
		}

		void calculatePrimesInBlock(size_t primeNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t maxPrimeNumberSearch = min(maxRangeSquareRoot + 1, blockEndNumber);

			for (; primeNumber < maxPrimeNumberSearch; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}
				}
			}
		}

		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template initPrimesBitset(maxRange);
		}
};

