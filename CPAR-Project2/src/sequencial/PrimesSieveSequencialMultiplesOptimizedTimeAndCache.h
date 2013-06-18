#pragma once

#include "PrimesSieveSequencialMultiplesOptimized.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>

using std::min;
using std::max;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedTimeAndCache: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	public:
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache(size_t blockSizeInBytes = sizeof(size_t) * 1024) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes / sizeof(size_t)) {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedTimeAndCache() {
		}
		
		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			this->template resetPrimesBitsetBlock();

			size_t primesFoundSize = this->template getPrimesValues().size();
			// primesFoundIndex = 1 because in position 0 is prime number 2, and since there are no even numbers, it is not necessary to sieve number 2
			for (size_t primesFoundIndex = 1; primesFoundIndex < primesFoundSize; ++primesFoundIndex) {
				size_t primeNumber = this->template getPrimesValues()[primesFoundIndex];
				size_t primeDoubled = primeNumber << 1;
				size_t compositeNumber = PrimesUtils::closestPrimeMultiple(primeNumber, blockBeginNumber);
				compositeNumber = max(compositeNumber, primeNumber * primeNumber);
				while (compositeNumber < blockBeginNumber || (compositeNumber % 2 == 0)) {
					compositeNumber += primeNumber;
				}

				for (; compositeNumber < blockEndNumber; compositeNumber += primeDoubled) {
					this->template setPrimesBitsetValueBlock(compositeNumber, blockBeginNumber, false);
				}
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t maxPrimeNumberSearch = min(maxRangeSquareRoot+1, blockEndNumber);
			size_t primeNumber = blockBeginNumber;

			for (; primeNumber < maxPrimeNumberSearch; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValueBlock(primeNumber, blockBeginNumber)) {
					this->template getPrimesValues().push_back(primeNumber);

					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValueBlock(compositeNumber, blockBeginNumber, false);
					}
				}
			}

			size_t numberToStartCulling = maxPrimeNumberSearch;
			if (maxPrimeNumberSearch % 2 == 0) {
				numberToStartCulling = maxPrimeNumberSearch + 1;
			}

			this->template collectRemainingPrimes(numberToStartCulling, blockBeginNumber, blockEndNumber);
		}

		void collectRemainingPrimes(size_t numberToStartCulling, size_t blockBeginNumber, size_t blockEndNumber) {
			size_t cullingIndex;
			if (numberToStartCulling >= blockBeginNumber) {
				cullingIndex = (numberToStartCulling - blockBeginNumber) >> 1;
			} else {
				cullingIndex = this->template getBitsetPositionToNumber(blockBeginNumber) % this->template getPrimesBitset().size();
			}
			size_t endIndex = this->template getPrimesBitset().size();
			for (; cullingIndex < endIndex; ++cullingIndex) {
				if (this->template getPrimesBitset()[cullingIndex]) {
					size_t numberToCull = blockBeginNumber + (cullingIndex << 1);
					if (numberToCull >= blockEndNumber)
						break;

					this->template getPrimesValues().push_back(numberToCull);
				}
			}
		}

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
			this->template getPrimesValues().push_back(2);
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			return this->template getPrimesValues();
		}

		virtual void initPrimesBitSetBlock(size_t maxRange) {
			size_t numberElements = min(this->template getNumberBitsToStore(maxRange), this->template getBlockSizeInElements());
			this->template initPrimesBitsetBlock(numberElements);
		}
	};

