#pragma once

#include "PrimesSieveSequencialMultiplesOptimized.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>

using std::min;
using std::max;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedTimeAndCache: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	protected:
		vector<size_t> _primesMultiples;

	public:
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache(size_t blockSizeInBytes = 204800) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedTimeAndCache() {
		}
		
		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			this->template resetPrimesBitsetBlock();

			size_t primesMultiplesSize = _primesMultiples.size();
			vector<size_t>& primesValues = this->template getPrimesValues();
			// primesFoundIndex = 1 because in position 0 is prime number 2, and since there are no even numbers, it is not necessary to sieve number 2
			for (size_t primesMultiplesIndex = 1; primesMultiplesIndex < primesMultiplesSize; ++primesMultiplesIndex) {
				size_t primeMultiple = _primesMultiples[primesMultiplesIndex];
				size_t primeDoubled = primesValues[primesMultiplesIndex] << 1;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeDoubled) {
					this->template setPrimesBitsetValueBlock(primeMultiple, blockBeginNumber, false);
				}

				_primesMultiples[primesMultiplesIndex] = primeMultiple;
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t maxPrimeNumberSearch = min(maxRangeSquareRoot + 1, blockEndNumber);
			size_t primeNumber = blockBeginNumber;

			vector<size_t>& primesValues = this->template getPrimesValues();

			for (; primeNumber < maxPrimeNumberSearch; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValueBlock(primeNumber, blockBeginNumber)) {
					primesValues.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValueBlock(compositeNumber, blockBeginNumber, false);
					}
					_primesMultiples.push_back(compositeNumber);
				}
			}

			size_t numberToStartCulling = maxPrimeNumberSearch;
			if (maxPrimeNumberSearch % 2 == 0) {
				numberToStartCulling = maxPrimeNumberSearch + 1;
			}

			this->template collectRemainingPrimes(numberToStartCulling, blockBeginNumber, blockEndNumber, maxRangeSquareRoot);
		}

		void collectRemainingPrimes(size_t numberToStartCulling, size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t cullingIndex;
			if (numberToStartCulling >= blockBeginNumber) {
				cullingIndex = (numberToStartCulling - blockBeginNumber) >> 1;
			} else {
				cullingIndex = this->template getBitsetPositionToNumber(blockBeginNumber) % this->template getPrimesBitset().size();
			}

			size_t endIndex = this->template getPrimesBitset().size();
			vector<size_t>& primesValues = this->template getPrimesValues();
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			for (; cullingIndex < endIndex; ++cullingIndex) {
				if (primesBitset[cullingIndex]) {
					size_t numberToCull = blockBeginNumber + (cullingIndex << 1);
					if (numberToCull >= blockEndNumber)
						break;

					primesValues.push_back(numberToCull);
					if (numberToCull < maxRangeSquareRoot) {
						_primesMultiples.push_back(numberToCull * numberToCull);
					}
				}
			}
		}

		virtual void clearPrimesValues() {
			size_t numberPrimesInRange = this->template getNumberOfPrimesInRange(this->template getMaxRange());
			this->template resetPrimesValues(numberPrimesInRange);
			this->template getPrimesValues().push_back(2);

			_primesMultiples = vector<size_t>(numberPrimesInRange);
			_primesMultiples.push_back(2);
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			return this->template getPrimesValues();
		}

		virtual void initPrimesBitSetBlock(size_t maxRange) {
			size_t numberElements = min(this->template getNumberBitsToStore(maxRange), this->template getBlockSizeInElements());
			this->template initPrimesBitsetBlock(numberElements);
		}
};

