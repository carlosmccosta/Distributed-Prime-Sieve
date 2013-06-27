#pragma once

#include "PrimesSieveSequencialMultiplesOptimized.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>
#include <utility>

using std::min;
using std::pair;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	protected:
		vector<pair<size_t, size_t> > _sievingPrimes;

	public:
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache(size_t blockSizeInBytes = 64 * 1024) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache() {
		}
		
		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			size_t sievingPrimesSize = _sievingPrimes.size();
			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				size_t primeDoubled = _sievingPrimes[sievingPrimesIndex].first;
				size_t primeMultiple = _sievingPrimes[sievingPrimesIndex].second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeDoubled) {
					this->template setPrimesBitsetValue(primeMultiple, false);
				}

				_sievingPrimes[sievingPrimesIndex].second = primeMultiple;
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

					_sievingPrimes.push_back(pair<size_t, size_t>(primeDoubled, compositeNumber));
				}
			}
		}

		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template initPrimesBitset(maxRange);

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.reserve(numberSievingPrimes);
		}
};

