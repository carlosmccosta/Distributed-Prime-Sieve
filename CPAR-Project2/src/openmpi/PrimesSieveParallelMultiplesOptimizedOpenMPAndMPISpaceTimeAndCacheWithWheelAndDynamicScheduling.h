#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling: public PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true, bool sendPrimesCountToRoot =
				true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling() {
		}

		virtual void computePrimes(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();

			performanceTimer.reset();
			performanceTimer.start();

			this->template setPrimesCount(0);
			this->template clearPrimesValues();

			int _processID = this->template getProcessId();
			int _numberProcesses = this->template getNumberProcesses();

			size_t processStartBlockNumber = max(this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange), this->template getBlockBeginNumber());
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);

			if (processStartBlockNumber % 2 == 0) {
				++processStartBlockNumber;
			}

			if (_processID == _numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}

			this->template setStartSieveNumber(this->template getBlockBeginNumber());
			this->template setMaxRange(processEndBlockNumber - 1); // processEndBlockNumber is max + 1 to use < operator instead of <=

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			vector<pair<size_t, size_t> > sievingMultiples;

			if (_processID == 0) {
				if (this->template isSendResultsToRoot()) {
					this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
				} else {
					this->template initPrimesBitSetSizeForRoot(maxRange, processEndBlockNumber);
				}
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}

			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);
			if (_processID != 0) {
				this->template initPrimesBitSetSizeForSieving(processEndBlockNumber - processStartBlockNumber);
			}

			++processEndBlockNumber; // force check of block limit

			this->template setStartSieveNumber(processStartBlockNumber);
			this->template removeComposites(processStartBlockNumber, processEndBlockNumber, sievingMultiples);

			if (_processID != 0) {
				performanceTimer.stop();
				cout << "    > Finish sieving in process with rank " << _processID << " in " << performanceTimer.getElapsedTimeFormated() << endl;
			}

			this->template syncProcesses(maxRange);
		}
};

