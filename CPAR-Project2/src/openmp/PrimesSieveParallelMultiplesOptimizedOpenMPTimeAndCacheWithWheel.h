#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMP.h"
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

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMP<FlagsContainer> {
	protected:
		vector<size_t> _sievingPrimes;
		WheelType wheelSieve;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel(size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0) :
				PrimesSieveParallelMultiplesOptimizedOpenMP<FlagsContainer>(blockSizeInBytes * 8, numberOfThreads) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel() {
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
//			cout << "init sievingMultiples in block [" << blockBeginNumber << ", " << blockEndNumber << "]" << endl;
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			size_t sievingMultiplesSize = sievingMultiples.size();
			for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
				pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					primesBitset[primeMultiple] = false;
				}

				sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
			}
		}

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxPrimeNumberSearch = blockEndNumber;
			if (maxPrimeNumberSearch >= maxRangeSquareRoot) {
				maxPrimeNumberSearch = maxRangeSquareRoot + 1;
			}
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			size_t primeNumber = blockBeginNumber;
			if (wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = wheelSieve.getNextPossiblePrime(primeNumber);
			}

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (primesBitset[primeNumber]) {
					_sievingPrimes.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						primesBitset[compositeNumber] = false;
					}
					sievingMultiples.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual inline size_t getNumberBitsToStore(size_t maxRange) {
			return (maxRange) + 1;
		}


		virtual void initPrimesBitSetSize(size_t maxRange) {
			this->template setMaxRange(maxRange);
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRange));

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
			int maxNumberThreads = omp_get_max_threads();
			size_t numberPrimesToCheckInBlock = maxRange / maxNumberThreads;

			#pragma omp parallel for \
			default(shared) \
			firstprivate(maxRange, maxNumberThreads, numberPrimesToCheckInBlock) \
			schedule(guided) \
			reduction(+: primesFound) \
			num_threads(maxNumberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < maxNumberThreads; ++threadBlockNumber) {
				size_t possiblePrime;
				if (threadBlockNumber == 0) {
					possiblePrime = 11;
				} else {
					possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock;
					if (!wheelSieve.isNumberPossiblePrime(possiblePrime)) {
						possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
					}
				}

				size_t nextPossiblePrimeNumberEndBlock;
				if (threadBlockNumber != maxNumberThreads - 1) {
					nextPossiblePrimeNumberEndBlock = (threadBlockNumber + 1) * numberPrimesToCheckInBlock;
				} else {
					nextPossiblePrimeNumberEndBlock = maxRange + 1;
				}

				primesFound = 0;
				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (primesBitset[possiblePrime]) {
						++primesFound;
					}
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}
			}

			return primesFound;
		}
};

