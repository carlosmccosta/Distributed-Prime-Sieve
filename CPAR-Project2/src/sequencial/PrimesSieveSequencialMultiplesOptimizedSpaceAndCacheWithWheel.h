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
class PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel: public PrimesSieveSequencialMultiplesOptimized<FlagsContainer> {
	protected:
		vector<pair<size_t, size_t> > _sievingPrimes;
		WheelType wheelSieve;

	public:
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel(size_t blockSizeInBytes = 64 * 1024) :
				PrimesSieveSequencialMultiplesOptimized<FlagsContainer>(blockSizeInBytes * 8) {
		}

		virtual ~PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel() {
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) {
			size_t sievingPrimesSize = _sievingPrimes.size();
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				size_t primeMultiple = _sievingPrimes[sievingPrimesIndex].first;
				size_t primeMultipleIncrement = _sievingPrimes[sievingPrimesIndex].second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					wheelSieve.setBitsetPositionToNumber(primesBitset, primeMultiple, false);
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
				size_t positionOnBitset = wheelSieve.getBitsetPositionToNumber(primeNumber);
				if (primesBitset[positionOnBitset]) {
					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						wheelSieve.setBitsetPositionToNumber(primesBitset, compositeNumber, false);
					}

					_sievingPrimes.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template setMaxRange(maxRange);
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(wheelSieve.getNumberBitsToStore(maxRange));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.clear();
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

			size_t possiblePrime = 7;
			if (wheelSieve.getNumberPrimesSievedByTheWheel() == 4) {
				primesValues.push_back(7);
				possiblePrime = 11;
			}

			size_t iSize = primesBitset.size();
			for (size_t i = 1; i < iSize; ++i) {   // position 0 has number 1 of spoke 1
				if (primesBitset[i]) {
					primesValues.push_back(possiblePrime);
				}

				possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
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

				size_t possiblePrime = 7;
				if (wheelSieve.getNumberPrimesSievedByTheWheel() == 4) {
					outputStream << 7 << endl;
					possiblePrime = 11;
				}

				size_t iSize = primesBitset.size();
				for (size_t i = 1; i < iSize; ++i) {   // position 2 has prime number 11
					if (primesBitset[i]) {
						outputStream << possiblePrime << endl;
					}

					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
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

			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t primesFound = wheelSieve.getNumberPrimesSievedByTheWheel();
			size_t iSize = primesBitset.size();
			for (size_t i = 1; i < iSize; ++i) { // in position 0 is spoke 1 of wheel and number 1 is not prime
				if (primesBitset[i]) {
					++primesFound;
				}
			}

			return primesFound;
		}
};

