#pragma once

#include "PrimesSieveSequencialMultiplesOptimized.h"
#include "../lib/PrimesUtils.h"
#include "../WheelFactorization.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <limits>

using std::min;
using std::max;
using std::pair;
using std::endl;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	protected:
		vector<pair<size_t, size_t> > _sievingPrimes;
		WheelType wheelSieve;

	public:
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel(size_t blockSizeInBytes = 64 * 1024) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel() {
		}
		
		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			size_t sievingPrimesSize = _sievingPrimes.size();

			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				size_t primeMultiple = _sievingPrimes[sievingPrimesIndex].first;
				size_t primeMultipleIncrement = _sievingPrimes[sievingPrimesIndex].second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					this->template setPrimesBitsetValue(primeMultiple, false);
				}

				_sievingPrimes[sievingPrimesIndex].first = primeMultiple;
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t maxPrimeNumberSearch = min(maxRangeSquareRoot + 1, blockEndNumber);

			size_t primeNumber = blockBeginNumber;
			if (wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = wheelSieve.getNextPossiblePrime(primeNumber);
			}

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}

					_sievingPrimes.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template initPrimesBitset(maxRange);

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual size_t getBlockBeginNumber() {
			return wheelSieve.getFirstPrimeToSieve();
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			vector<size_t>& primesValues = this->template getPrimesValues();

			primesValues.clear();
			primesValues.push_back(2);
			primesValues.push_back(3);
			primesValues.push_back(5);
			primesValues.push_back(7);

			size_t maxRange = this->template getMaxRange();
			for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (this->template getPrimesBitsetValue(possiblePrime)) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
		}

		virtual void savePrimes(ostream& outputStream) {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() <= 2) {
				outputStream << 2 << endl;
				outputStream << 3 << endl;
				outputStream << 5 << endl;
				outputStream << 7 << endl;

				size_t maxRange = this->template getMaxRange();
				for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (this->template getPrimesBitsetValue(possiblePrime)) {
						outputStream << possiblePrime << endl;
					}
				}
			} else {
				size_t iSize = primesValues.size();
				for (size_t i = 0; i < iSize; ++i) {
					outputStream << primesValues[i] << endl;
				}
			}
		}

		virtual size_t getNumberPrimesFound() {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() >= 2)
				return primesValues.size();

			size_t primesFound = 4;
			size_t maxRange = this->template getMaxRange();
			for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (this->template getPrimesBitsetValue(possiblePrime)) {
					++primesFound;
				}
			}

			return primesFound;
		}
};

