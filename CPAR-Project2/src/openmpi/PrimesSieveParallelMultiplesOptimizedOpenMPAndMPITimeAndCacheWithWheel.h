#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel.h"

#include <omp.h>

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPAndMPITimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<FlagsContainer, WheelType> {
	protected:
		size_t _numberOfThreads;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPITimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot), _numberOfThreads(
						numberOfThreads) {
		}
		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPITimeAndCacheWithWheel() {
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			int _processID = this->template getProcessId();
			cout << "    --> Removing composites in process with rank " << _processID << " in [" << processBeginBlockNumber << ", " << (processEndBlockNumber - 1) << "]" << endl;

			const size_t blockSizeInElements = this->template getBlockSizeInElements();
			const size_t numberBlocks = ceil((double) (processEndBlockNumber - processBeginBlockNumber) / (double) blockSizeInElements);
			size_t numberThreadsToUse = omp_get_max_threads();
			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					numberThreadsToUse = numberBlocks;
				} else {
					numberThreadsToUse = _numberOfThreads;
				}
			}

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

#pragma omp parallel for \
			default(shared) \
			firstprivate(processBeginBlockNumber, processEndBlockNumber, blockSizeInElements, sievingMultiples, priviousBlockNumber) \
			schedule(guided, 128) \
			num_threads(numberThreadsToUse)
			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockBeginNumber = blockNumber * blockSizeInElements + processBeginBlockNumber;
				size_t blockEndNumber = blockBeginNumber + blockSizeInElements;

				if (blockEndNumber > processEndBlockNumber) {
					blockEndNumber = processEndBlockNumber;
				}

				if (blockNumber == 0 && this->template getProcessId() == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
		}

		virtual size_t getNumberPrimesFound() {
			size_t primesFound = this->template getPrimesCount();
			// avoid recomputation
			if (primesFound != 0) {
				return primesFound;
			}

			vector<size_t>& primesValues = this->template getPrimesValues();
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			if (primesValues.size() >= 2) {
				return primesValues.size();
			}

			size_t startSieveNumber = this->template getStartSieveNumber();
			if (startSieveNumber == this->template getBlockBeginNumber()) {
				primesFound = this->template getWheelSieve().getNumberPrimesSievedByTheWheel();
			} else {
				primesFound = 0;
			}

			size_t maxRange = this->template getMaxRange();
			int maxNumberThreads = omp_get_max_threads();
			size_t minNumberPrimesPerThread = 100;
			int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));
			size_t numberPrimesToCheckInBlock = (maxRange - startSieveNumber) / numberThreads;

			WheelType& wheelSieve = this->template getWheelSieve();

#pragma omp parallel for \
			default(shared) \
			firstprivate(maxRange, numberThreads, numberPrimesToCheckInBlock, startSieveNumber) \
			schedule(static) \
			reduction(+: primesFound) \
			num_threads(numberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
				size_t possiblePrime;

				possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + startSieveNumber;
				if (!(wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t nextPossiblePrimeNumberEndBlock;
				if (threadBlockNumber == numberThreads - 1) {
					nextPossiblePrimeNumberEndBlock = maxRange + 1;
				} else {
					nextPossiblePrimeNumberEndBlock = (threadBlockNumber + 1) * numberPrimesToCheckInBlock + startSieveNumber;
				}

				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (primesBitset[possiblePrime]) {
						++primesFound;
					}
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}
			}

			this->template setPrimesCount(primesFound);
			return primesFound;
		}

		size_t getNumberOfThreads() const {
			return _numberOfThreads;
		}

		void setNumberOfThreads(size_t numberOfThreads) {
			_numberOfThreads = numberOfThreads;
		}
};

