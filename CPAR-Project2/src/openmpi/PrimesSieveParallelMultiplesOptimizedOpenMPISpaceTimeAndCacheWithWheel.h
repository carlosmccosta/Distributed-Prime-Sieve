#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"
#include "../lib/PrimesUtils.h"
#include "../WheelFactorization.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <limits>
#include <vector>

using std::vector;
using std::min;
using std::max;
using std::pair;
using std::endl;
using std::cout;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType> {
	protected:
		vector<size_t> _sievingPrimes;
		WheelType _wheelSieve;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, bool sendResultsToRoot = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, sendResultsToRoot, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel() {
		}

		virtual void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			sievingMultiples.clear();
			size_t sievingPrimesSize = _sievingPrimes.size();
			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				size_t primeNumber = _sievingPrimes[sievingPrimesIndex];
				size_t primeMultiple = PrimesUtils::closestPrimeMultiple(primeNumber, blockBeginNumber);
				size_t primeMultipleIncrement = primeNumber << 1;

				if (primeMultiple < blockBeginNumber || primeMultiple == primeNumber) {
					primeMultiple += primeNumber;
				}

				if (primeMultiple % 2 == 0) {
					primeMultiple += primeNumber;
				}

				sievingMultiples.push_back(pair<size_t, size_t>(primeMultiple, primeMultipleIncrement));
			}
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t sievingMultiplesSize = sievingMultiples.size();
			for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
				pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					this->template setPrimesBitsetValueMPI(primeMultiple, false);
				}

				sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxPrimeNumberSearch = blockEndNumber;
			if (maxPrimeNumberSearch >= maxRangeSquareRoot) {
				maxPrimeNumberSearch = maxRangeSquareRoot + 1;
			}

			size_t primeNumber = blockBeginNumber;
			if (_wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber);
			}

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValueMPI(primeNumber)) {
					_sievingPrimes.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						this->template setPrimesBitsetValueMPI(compositeNumber, false);
					}
					sievingMultiples.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual void initPrimesBitSetSizeForRoot(size_t maxRange, size_t blockEndNumber) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(blockEndNumber));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual void initPrimesBitSetSizeForRootWithAllValues(size_t maxRange) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRange));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual void initPrimesBitSetSizeForSievingPrimes(size_t maxRangeSquareRoot) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStoreSievingPrimes(maxRangeSquareRoot));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange(maxRangeSquareRoot);
			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual void initPrimesBitSetSizeForSieving(size_t blockSize) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize((blockSize >> 1) + 1);
		}

		virtual inline size_t getNumberBitsToStore(size_t maxRange) {
			return ((maxRange - this->template getStartSieveNumber()) >> 1) + 1;
		}

		virtual inline size_t getNumberBitsToStoreSievingPrimes(size_t maxRange) {
			return ((maxRange - this->template getBlockBeginNumber()) >> 1) + 1;
		}

		virtual inline size_t getBlockBeginNumber() {
			return _wheelSieve.getFirstPrimeToSieve();
		}

		inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return (position << 1) + this->template getStartSieveNumber();
		}

		inline bool getPrimesBitsetValueMPI(size_t number) {
			return this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)];
		}

		inline void setPrimesBitsetValueMPI(size_t number, bool newValue) {
			this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)] = newValue;
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			vector<size_t>& primesValues = this->template getPrimesValues();

			primesValues.clear();
			primesValues.push_back(2);
			primesValues.push_back(3);
			primesValues.push_back(5);
			primesValues.push_back(7);

			size_t possiblePrime = this->template getStartSieveNumber();
			size_t maxRange = this->template getMaxRange();
			for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (this->template getPrimesBitsetValueMPI(possiblePrime)) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
		}

		virtual void savePrimes(ostream& outputStream) {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() <= 2) {

				if (this->template getProcessId() == 0) {
					outputStream << 2 << endl;
					outputStream << 3 << endl;
					outputStream << 5 << endl;
					outputStream << 7 << endl;
				}

				size_t possiblePrime = this->template getStartSieveNumber();
				if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t maxRange = this->template getMaxRange();
				for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (this->template getPrimesBitsetValueMPI(possiblePrime)) {
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
			size_t primesFound = this->template getPrimesCount();
			// avoid recomputation
			if (primesFound != 0) {
				return primesFound;
			}

			vector<size_t>& primesValues = this->template getPrimesValues();
			if (primesValues.size() >= 2) {
				return primesValues.size();
			}

			size_t possiblePrime = this->template getStartSieveNumber();

			if (possiblePrime == this->template getBlockBeginNumber()) {
				primesFound = _wheelSieve.getNumberPrimesSievedByTheWheel();
			} else {
				primesFound = 0;
			}

			if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
				possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
			}
			size_t maxRange = this->template getMaxRange();
			for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (this->template getPrimesBitsetValueMPI(possiblePrime)) {
					++primesFound;
				}
			}

			this->template setPrimesCount(primesFound);
			return primesFound;
		}

		inline const vector<size_t>& getSievingPrimes() const {
			return _sievingPrimes;
		}

		inline void setSievingPrimes(const vector<size_t>& sievingPrimes) {
			_sievingPrimes = sievingPrimes;
		}

		inline WheelType& getWheelSieve() {
			return _wheelSieve;
		}

		inline void setWheelSieve(WheelType wheelSieve) {
			this->_wheelSieve = wheelSieve;
		}
};

