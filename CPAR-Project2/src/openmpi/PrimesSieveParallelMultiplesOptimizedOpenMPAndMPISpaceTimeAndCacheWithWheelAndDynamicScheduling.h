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
		size_t _dynamicSchedulingBlockSizeInElements;

		unsigned long long _numberSegmentsSieved;
		vector<vector<pair<size_t, size_t> > > _blockDistribution;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode,
						sendPrimesCountToRoot), _processStartBlockNumber(0), _processEndBlockNumber(0), _dynamicSchedulingBlockSizeInElements(1048576), _numberSegmentsSieved(0) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling() {
		}

		virtual void computePrimes(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			bool _sendResultsToRoot = this->template isSendResultsToRoot();

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
				if (_sendResultsToRoot) {
					this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
				}
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);

			// sieve all assigned blocks
			while (this->template getNewBlockFromRoot()) {
				if (processID != 0 || (processID == 0 && !_sendResultsToRoot)) {
					this->template initPrimesBitSetSizeForSieving(_processEndBlockNumber - _processStartBlockNumber);
				}
				this->template removeComposites(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);

				if (processID != 0) {
					this->template sendSievedBlockToRoot();
				}
			}

			performanceTimer.stop();
			cout << "    > Finish sieving " << _numberSegmentsSieved << " blocks in process with rank " << processID << " in " << performanceTimer.getElapsedTimeFormated() << endl;

			if (processID == 0) {
				cout << "\n    >>>>> Finish sieving in all processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;
			}
		}

		virtual void initBlockDistribution(size_t maxRange) {
			size_t numberProcesses = (size_t) this->template getNumberProcesses();
			_blockDistribution = vector<vector<pair<size_t, size_t> > >(numberProcesses, vector<pair<size_t, size_t> >());

			size_t dynamicSchedulingBlockSizeInElements = _dynamicSchedulingBlockSizeInElements;
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			size_t numberBlocks = ceil((double) (maxRange - maxRangeSquareRoot) / (double) _dynamicSchedulingBlockSizeInElements);
			size_t numberBlocksPerProcess = numberBlocks / numberProcesses;

			size_t numberBlocksAdded = 0;

			size_t blockNumber = 0;
#pragma omp parallel for \
				default(shared) \
				firstprivate(numberProcesses, numberBlocksPerProcess, dynamicSchedulingBlockSizeInElements, maxRangeSquareRoot) \
				reduction(+: numberBlocksAdded)
			for (size_t processNumber = 0; processNumber < numberProcesses; ++processNumber) {
				for (size_t processBlockNumber = 0; processBlockNumber < numberBlocksPerProcess; ++processBlockNumber) {
					size_t blockBeginNumber = (processNumber * numberBlocksPerProcess + blockNumber) * dynamicSchedulingBlockSizeInElements + maxRangeSquareRoot;
					size_t blockEndNumber = blockBeginNumber + dynamicSchedulingBlockSizeInElements;

					if (blockBeginNumber % 2 == 0) {
						++blockBeginNumber;
					}

					if (blockEndNumber > maxRange) {
						blockEndNumber = maxRange;
					}
					_blockDistribution[processNumber].push_back(pair<size_t, size_t>(blockBeginNumber, blockEndNumber));
					++numberBlocksAdded;
				}
			}

			// add any missing blocks to last process because of division rounding
			size_t lastProcessNumberIndex = numberProcesses - 1;
			for (; numberBlocksAdded < numberBlocks; ++numberBlocksAdded) {
				size_t blockBegin = numberBlocksAdded * dynamicSchedulingBlockSizeInElements + maxRangeSquareRoot;
				size_t blockEnd = blockBegin + dynamicSchedulingBlockSizeInElements;

				if (blockEnd > maxRange) {
					blockEnd = maxRange;
				}
				_blockDistribution[lastProcessNumberIndex].push_back(pair<size_t, size_t>(blockBegin, blockEnd));
			}
		}

		virtual bool getNewBlockFromRoot() {
			++_numberSegmentsSieved;
			unsigned long long segmentRange[2];

			MPI_Status status;

			MPI_Send(&_numberSegmentsSieved, 1, MPI_UNSIGNED_LONG_LONG, 0, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD);
			MPI_Recv(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, 0, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
				_processStartBlockNumber = segmentRange[0];
				_processEndBlockNumber = segmentRange[1];

				if (_processStartBlockNumber == 0 && _processEndBlockNumber == 0) {
					return false;
				}

				this->template setStartSieveNumber(_processStartBlockNumber);
				this->template setMaxRange(_processEndBlockNumber - 1); // processEndBlockNumber is max + 1 to use < operator instead of <=
				return true;
			} else {
				return false;
			}
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
					cout << "\n    > Exporting results to file " << _outputResultsFilename << "..." << endl;
				}

				PerformanceTimer performanceTimer;
				performanceTimer.reset();
				performanceTimer.start();
				if (this->template savePrimesToFile(_outputResultsFilename)) {
					performanceTimer.stop();
					cout << "    --> Export to file " << _outputResultsFilename << " finished in " << performanceTimer.getElapsedTimeFormated() << endl;
				} else {
					cerr << "    --> Export to file failed!" << endl;
				}
			} else {
				cerr << "\n    --> Export failed!" << endl;
				return false;
			}

			return true;
		}

};

