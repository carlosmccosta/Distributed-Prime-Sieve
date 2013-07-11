#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel.h"

#include <vector>
#include <utility>
#include <string>

using std::string;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling: public PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType> {
	protected:
		size_t _processStartBlockNumber;
		size_t _processEndBlockNumber;
		string _outputResultsFilename;
		string _resultsConfirmationFile;
		size_t _dynamicSchedulingBlockSizeInElements;

		size_t _numerBlocksSieved;
		vector<vector<pair<size_t, size_t> > > _blockDistribution;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true, bool sendPrimesCountToRoot =
				true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, sendPrimesCountToRoot), _processStartBlockNumber(0), _processEndBlockNumber(
						0), _dynamicSchedulingBlockSizeInElements(1048576), _numerBlocksSieved(0) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling() {
		}

		virtual void computePrimes(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();

			performanceTimer.reset();
			performanceTimer.start();

			this->template setPrimesCount(0);
			this->template clearPrimesValues();

			int processID = this->template getProcessId();

			if (processID == 0) {
				this->template initBlockDistribution(maxRange);
			}

			// init primes used in sieving
			this->template setStartSieveNumber(this->template getBlockBeginNumber());
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			vector<pair<size_t, size_t> > sievingMultiples;
			if (processID == 0) {
				this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);

			// sieve all assigned blocks
			while (this->template getNewBlockFromRoot()) {
				if (processID != 0) {
					this->template initPrimesBitSetSizeForSieving(_processEndBlockNumber - _processStartBlockNumber);
				}
				this->template removeComposites(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);

				if (processID != 0) {
					this->template sendSievedBlockToRoot();
				}
			}

			performanceTimer.stop();
			cout << "    > Finish sieving " << _numerBlocksSieved << " blocks in process with rank " << processID << " in " << performanceTimer.getElapsedTimeFormated() << endl;

			if (processID == 0) {
				cout << "\n    >>>>> Finish sieving in all processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;
			}
		}

		virtual void initBlockDistribution(size_t maxRange) {
			int numberProcesses = this->template getNumberProcesses();
			_blockDistribution = vector<vector<pair<size_t, size_t> > >(numberProcesses, vector<pair<size_t, size_t> >());

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			const size_t numberBlocks = ceil((double) (maxRange - maxRangeSquareRoot) / (double) _dynamicSchedulingBlockSizeInElements);

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockBegin = blockNumber * _dynamicSchedulingBlockSizeInElements + maxRangeSquareRoot;
				size_t blockEnd = blockBegin + _dynamicSchedulingBlockSizeInElements;

				if (blockEnd > maxRange) {
					blockEnd = maxRange;
				}

//				_blockDistribution.push_back(pair<size_t, size_t>(blockBegin, blockEnd));

			}
		}

		virtual bool getNewBlockFromRoot() {

			this->template setMaxRange(_processEndBlockNumber - 1); // processEndBlockNumber is max + 1 to use < operator instead of <=
			++_processEndBlockNumber; // force check of block limit
			this->template setStartSieveNumber(_processStartBlockNumber);

			return false;
		}

		virtual void blocksAvailableServer() {

		}

		virtual void sendSievedBlockToRoot() {

		}

		virtual bool outputResults() {
			if (_outputResultsFilename == "stdout") {
				cout << "\n\n=============================================  Computed primes  =============================================\n\n";
				this->template printPrimesToConsole();
				cout << "\n" << endl;
			} else if (_outputResultsFilename != "") {
				int numberProcesses = this->template getNumberProcesses();
				bool sendResultsToRoot = this->template isSendResultsToRoot();

				if (sendResultsToRoot || numberProcesses == 1) {
					cout << "\n    > Exporting results to file " << _outputResultsFilename << "...";
				}
				this->template savePrimesToFile(_outputResultsFilename);
				if (sendResultsToRoot || numberProcesses == 1) {
					cout << "\n    --> Export to file " << _outputResultsFilename << " finished!" << endl;
				}
			} else {
				return false;
			}

			return true;
		}

};

