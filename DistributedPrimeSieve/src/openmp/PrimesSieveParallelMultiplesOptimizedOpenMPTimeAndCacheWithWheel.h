#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMP.h"
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
class PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMP<FlagsContainer, WheelType> {
	protected:
		WheelType _wheelSieve;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel(size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0) :
				PrimesSieveParallelMultiplesOptimizedOpenMP<FlagsContainer, WheelType>(blockSizeInBytes * 8, numberOfThreads) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel() {
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			FlagsContainer& primesBitset = this->template getPrimesBitset();

			size_t sievingMultiplesSize = sievingMultiples.size();

			// worse performance with omp parallel because of false sharing
			/*	#pragma omp parallel for \
				default(shared) \
				schedule(guided) \
				firstprivate(sievingMultiplesSize, blockEndNumber)*/
			for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
				pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					primesBitset[primeMultiple] = true;
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
			if (_wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber);
			}

			vector<size_t>& _sievingPrimes = this->template getSievingPrimes();

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (!primesBitset[primeNumber]) {
					_sievingPrimes.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						primesBitset[compositeNumber] = true;
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
			vector<size_t>& _sievingPrimes = this->template getSievingPrimes();

			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual size_t getBlockBeginNumber() {
			return _wheelSieve.getFirstPrimeToSieve();
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
			for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (!primesBitset[possiblePrime]) {
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
				for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (!primesBitset[possiblePrime]) {
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
			size_t minNumberPrimesPerThread = 100;
			int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));
			size_t numberPrimesToCheckInBlock = maxRange / numberThreads;

#pragma omp parallel for \
			default(shared) \
			firstprivate(maxRange, numberThreads, numberPrimesToCheckInBlock) \
			schedule(guided) \
			reduction(+: primesFound) \
			num_threads(numberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
				size_t possiblePrime;

				possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + 11;
				if (!_wheelSieve.isNumberPossiblePrime(possiblePrime)) {
					possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t nextPossiblePrimeNumberEndBlock = min((threadBlockNumber + 1) * numberPrimesToCheckInBlock + 11, maxRange + 1);

				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (!primesBitset[possiblePrime]) {
						++primesFound;
					}
					possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
				}
			}

			return primesFound;
		}

		inline WheelType& getWheelSieve() {
			return _wheelSieve;
		}
};
