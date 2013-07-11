#pragma once

#include "../PrimesSieve.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <mpi.h>

using std::sqrt;
using std::min;
using std::max;
using std::pair;

template<typename FlagsContainer>
class PrimesSieveParallelMultiplesOptimizedOpenMPI: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;
		int _processID;
		int _numberProcesses;
		bool _sendResultsToRoot;
		bool _sendPrimesCountToRoot;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPI(size_t maxRange, size_t blockSizeInElements = 16 * 1024, bool sendResultsToRoot = true, bool sendPrimesCountToRoot = true) :
				_blockSizeInElements(blockSizeInElements), _sendResultsToRoot(sendResultsToRoot), _sendPrimesCountToRoot(sendPrimesCountToRoot) {
			int flag;
			MPI_Initialized(&flag);
			if (!flag) {
				MPI_Init(NULL, NULL);
			}

			MPI_Comm_rank(MPI_COMM_WORLD, &_processID);
			MPI_Comm_size(MPI_COMM_WORLD, &_numberProcesses);

			size_t processStartBlockNumber = max(this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange), this->template getBlockBeginNumber());
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);

			if (_processID == _numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}

			this->template setStartSieveNumber(processStartBlockNumber);
			this->template setMaxRange(processEndBlockNumber - 1);
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPI() {
//			int flag;
//			MPI_Finalized(&flag);
//			if (!flag) {
//				MPI_Finalize();
//			}
		}

		virtual void computePrimes(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();

			performanceTimer.reset();
			performanceTimer.start();

			this->template setPrimesCount(0);
			this->template clearPrimesValues();

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
				if (_sendResultsToRoot) {
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

			++processEndBlockNumber;

			this->template setStartSieveNumber(processStartBlockNumber);
			this->template removeComposites(processStartBlockNumber, processEndBlockNumber, sievingMultiples);

			if (_processID != 0) {
				performanceTimer.stop();
				cout << "    > Finish sieving in process with rank " << _processID << " in " << performanceTimer.getElapsedTimeFormated() << endl;
			}

			if (_processID == 0) {
				if (_numberProcesses > 1) {
					cout << "    > Finish sieving in root process in " << performanceTimer.getElapsedTimeFormated() << endl;
					this->template waitForAllProcessesInGroupToFinishSieving();
					performanceTimer.stop();
					cout << "    > Finish sieving in all processes in " << performanceTimer.getElapsedTimeFormated() << endl;

					this->template setStartSieveNumber(this->template getBlockBeginNumber());
					if (_sendPrimesCountToRoot) {
						this->template collectPrimesCountFromProcessGroup();
					}

					if (_sendResultsToRoot) {
						this->template setMaxRange(maxRange);
						this->template collectResultsFromProcessGroup(maxRange);
					}
				} else {
					performanceTimer.stop();
				}
			} else {
				this->template sendFinishSievingMessageToRoot();

				if (_sendPrimesCountToRoot) {
					this->template sendPrimesCountToRoot();
				}

				if (_sendResultsToRoot) {
					this->template sendResultsToRootProcess(maxRange);
				}
			}
		}

		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumberMPI(maxRangeSquareRoot);
			size_t blockBeginNumber = this->template getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t blockIndexBegin = this->template getBitsetPositionToNumberMPI(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + _blockSizeInElements, maxIndexRangeSquareRoot + 1);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;
				if (blockIndexEnd > maxIndexRangeSquareRoot) {
					blockIndexEnd = maxIndexRangeSquareRoot;
				}

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + processBeginBlockNumberIndex;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				if (blockIndexEnd > processEndBlockNumberIndex) {
					blockIndexEnd = processEndBlockNumberIndex;
				}

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockNumber == 0 && _processID == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void sendFinishSievingMessageToRoot() {
			double elapsedTimeMicroSec = this->template getPerformanceTimer().getElapsedTimeInMicroSec();
			MPI_Send(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}

		virtual void waitForAllProcessesInGroupToFinishSieving() {
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				double elapsedTimeMicroSec;
				MPI_Status status;
				MPI_Recv(&elapsedTimeMicroSec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR != MPI_SUCCESS) {
					cout << "    --> MPI_Recv detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
		}

		virtual void sendPrimesCountToRoot() {
			unsigned long long primesFound = (unsigned long long)this->template getNumberPrimesFound();
			MPI_Send(&primesFound, 1, MPI_LONG_LONG, 0, 2, MPI_COMM_WORLD);
		}

		virtual void collectPrimesCountFromProcessGroup() {
			PerformanceTimer countingPrimesTimer;
			countingPrimesTimer.reset();
			countingPrimesTimer.start();

			// compute primes in this process
			size_t numberPrimesFound = this->template getNumberPrimesFound();
			cout << "    --> Process " << 0 << " counted " << numberPrimesFound << " primes" << endl;

			// update primes count with the parcial count from the remaining processes
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				unsigned long long primesCount;
				MPI_Status status;
				MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					this->template incrementPrimesCount((size_t) primesCount);
					cout << "    --> Process " << status.MPI_SOURCE << " counted " << primesCount << " primes" << endl;
				} else {
					cout << "    --> MPI_Recv detected the following error code: " << status.MPI_ERROR << endl;
				}
			}

			numberPrimesFound = this->template getNumberPrimesFound(); // get updated value
			cout << "    --> Counted " << numberPrimesFound << " primes on " << _numberProcesses << " processes " << " in " << countingPrimesTimer.getElapsedTimeFormated() << endl;
			countingPrimesTimer.stop();
		}

		virtual void sendResultsToRootProcess(size_t maxRange) {
			cout << "    > Sending results from process with rank " << _processID << " to root process..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t processStartBlockNumber = this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange);
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);

			if (processStartBlockNumber % 2 == 0) {
				++processStartBlockNumber;
			}

			if (_processID == _numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}
			size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;

			MPI_Send(&primesBitset[0], blockSize, MPI_UNSIGNED_CHAR, 0, 3, MPI_COMM_WORLD);
			cout << "    --> Finished sending results from process with rank " << _processID << " to root process" << endl;
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
			cout << "\n    > Collecting results from other processes..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				cout << "    > Probing for results..." << endl;
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, _numberProcesses, maxRange);
					size_t processEndBlockNumber = this->template getProcessEndBlockNumber(processID, _numberProcesses, maxRange);

					if (processStartBlockNumber % 2 == 0) {
						++processStartBlockNumber;
					}

					if (processID == _numberProcesses - 1) {
						processEndBlockNumber = maxRange + 1;
					}
					size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					cout << "    --> Collecting results from process with rank " << processID << endl;
					MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
					cout << "    --> Finished collecting results from process with rank " << processID << endl;
				} else {
					cout << "    --> MPI_Probe detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
			cout << "    --> Finished collecting all results\n" << endl;
		}

		virtual void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) = 0;

		virtual void initPrimesBitSetSizeForRoot(size_t maxRange, size_t blockEndNumber) = 0;
		virtual void initPrimesBitSetSizeForRootWithAllValues(size_t maxRange) = 0;
		virtual void initPrimesBitSetSizeForSievingPrimes(size_t maxRangeSquareRoot) = 0;
		virtual void initPrimesBitSetSizeForSieving(size_t maxRange) = 0;

		inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return (position << 1) + this->template getStartSieveNumber();
		}

		inline size_t getProcessStartBlockNumber(size_t processID, size_t numberProcesses, size_t maxRange) {
			return processID * maxRange / numberProcesses;
		}

		inline size_t getProcessEndBlockNumber(size_t processID, size_t numberProcesses, size_t maxRange) {
			return ((processID + 1) * maxRange / numberProcesses);
		}

		inline size_t getProcessBlockSize(size_t processID, size_t numberProcesses, size_t maxRange) {
			return (this->template getProcessEndBlockNumber(processID, numberProcesses, maxRange) - this->template getProcessStartBlockNumber(processID, numberProcesses, maxRange));
		}

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
		}

		virtual inline size_t getBlockBeginNumber() {
			return this->template getNumberAssociatedWithBitsetPositionMPI(0);
		}

		inline size_t getBlockSizeInBytes() {
			return _blockSizeInElements / 8;
		}

		inline void setBlockSizeInBytes(size_t blockSize) {
			_blockSizeInElements = blockSize * 8;
		}

		inline size_t getBlockSizeInElements() const {
			return _blockSizeInElements;
		}

		inline void setBlockSizeInElements(size_t blockSizeInElements) {
			_blockSizeInElements = blockSizeInElements;
		}

		inline int getNumberProcesses() const {
			return _numberProcesses;
		}

		inline int getProcessId() const {
			return _processID;
		}
	};

