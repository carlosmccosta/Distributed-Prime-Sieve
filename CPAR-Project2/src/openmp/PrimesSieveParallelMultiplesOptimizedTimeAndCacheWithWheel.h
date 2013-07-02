#pragma once

#include "PrimesSieveParallelMultiplesOptimized.h"
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
class PrimesSieveParallelMultiplesOptimizedTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimized<FlagsContainer> {
	protected:
		vector<pair<size_t, size_t> > _sievingPrimes;
		WheelType wheelSieve;

	public:
		PrimesSieveParallelMultiplesOptimizedTimeAndCacheWithWheel(size_t blockSizeInBytes = 64 * 1024) :
				PrimesSieveParallelMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}
		
		virtual ~PrimesSieveParallelMultiplesOptimizedTimeAndCacheWithWheel() {
		}
		
		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			size_t sievingPrimesSize = _sievingPrimes.size();
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				pair<size_t, size_t> primeCompositeInfo = _sievingPrimes[sievingPrimesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					primesBitset[primeMultiple] = false;
				}

				_sievingPrimes[sievingPrimesIndex].first = primeMultiple;
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot) {
			size_t maxPrimeNumberSearch = min(maxRangeSquareRoot + 1, blockEndNumber);
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			size_t primeNumber = blockBeginNumber;
			if (wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = wheelSieve.getNextPossiblePrime(primeNumber);
			}

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (primesBitset[primeNumber]) {
					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						primesBitset[compositeNumber] = false;
					}

					_sievingPrimes.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template setMaxRange(maxRange);
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(maxRange + 1);

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual size_t getBlockBeginNumber() {
			return wheelSieve.getFirstPrimeToSieve();
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			vector<size_t>& primesValues = this->template getPrimesValues();
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			primesValues.clear();
			primesValues.push_back(2);
			primesValues.push_back(3);
			primesValues.push_back(5);
			primesValues.push_back(7);

			size_t maxRange = this->template getMaxRange();
			for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (primesBitset[possiblePrime]) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
		}

		virtual void savePrimes(ostream& outputStream) {
			vector<size_t>& primesValues = this->template getPrimesValues();
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			if (primesValues.size() <= 2) {
				outputStream << 2 << endl;
				outputStream << 3 << endl;
				outputStream << 5 << endl;
				outputStream << 7 << endl;

				size_t maxRange = this->template getMaxRange();
				for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (primesBitset[possiblePrime]) {
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
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			if (primesValues.size() >= 2)
				return primesValues.size();

			size_t primesFound = 4;
			size_t maxRange = this->template getMaxRange();
			for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (primesBitset[possiblePrime]) {
					++primesFound;
				}
			}

			return primesFound;
		}
};

