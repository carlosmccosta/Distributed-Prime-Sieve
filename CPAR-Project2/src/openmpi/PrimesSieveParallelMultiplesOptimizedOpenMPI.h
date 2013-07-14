#pragma once

#include "../PrimesSieve.h"
#include "../WheelFactorization.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <mpi.h>

using std::sqrt;
using std::min;
using std::max;
using std::pair;

enum MessageTags {
	MSG_NODE_SIEVING_FINISHED = 0, MSG_NODE_PRIMES_FOUND_COUNT, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MSG_REQUEST_NEW_SEGMET, MSG_ASSIGN_NEW_SEGMET
};

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPI: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;
		int _processID;
		int _numberProcesses;
		bool _sendResultsToRoot;
		bool _countNumberOfPrimesOnNode;
		bool _sendPrimesCountToRoot;

		vector<size_t> _sievingPrimes;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPI(size_t maxRange, size_t blockSizeInElements = 16 * 1024, bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true) :
				_blockSizeInElements(blockSizeInElements), _sendResultsToRoot(sendResultsToRoot), _countNumberOfPrimesOnNode(countNumberOfPrimesOnNode), _sendPrimesCountToRoot(sendPrimesCountToRoot) {
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
		}

		virtual void computePrimes(size_t maxRange) {
			PerformanceTimer& totalPerformanceTimer = this->template getPerformanceTimer();

			totalPerformanceTimer.reset();
			totalPerformanceTimer.start();

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
					this->template initPrimesBitSetSizeForRoot(maxRange, processEndBlockNumber - 1);
				}
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}

			// compute sieving primes
			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.reset();
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);
			computeSievingPrimesTimer.stop();
			cout << "    --> Computed sieving primes in process with rank " << _processID << " in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;

			// remove composites
			if (_processID == 0) {
				this->template removeComposites(maxRangeSquareRoot, processEndBlockNumber, sievingMultiples);
			} else {
				this->template setStartSieveNumber(processStartBlockNumber);
				this->template initPrimesBitSetSizeForSieving(processEndBlockNumber - processStartBlockNumber);
				this->template removeComposites(processStartBlockNumber, processEndBlockNumber, sievingMultiples);
			}

			if (_processID != 0) {
				totalPerformanceTimer.stop();
			}

			cout << "    > Finish sieving in process with rank " << _processID << " in " << totalPerformanceTimer.getElapsedTimeFormated() << endl;

			this->template syncProcesses(maxRange);
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

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber >= maxRangeSquareRoot) {
					blockEndNumber = maxRangeSquareRoot + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			cout << "    --> Removing composites in process with rank " << _processID << " in [" << processBeginBlockNumber << ", " << (processEndBlockNumber - 1)<< "]" << endl;

			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + processBeginBlockNumberIndex;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber > processEndBlockNumber) {
					blockEndNumber = processEndBlockNumber;
				}

				if (blockNumber == 0 && _processID == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void syncProcesses(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();

			if (_processID == 0) {
				if (_numberProcesses > 1) {
					this->template waitForAllProcessesInGroupToFinishSieving();
					performanceTimer.stop();
					cout << "\n    >>>>> Finish sieving in all processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;

					if (_sendPrimesCountToRoot) {
						this->template collectPrimesCountFromProcessGroup();
					} else if (_countNumberOfPrimesOnNode) {
						this->template countPrimesInNode();
					}

					if (_sendResultsToRoot) {
						this->template setMaxRange(maxRange);
						this->template collectResultsFromProcessGroup(maxRange);

						if (_countNumberOfPrimesOnNode) {
							//force recount with results from other processes
							this->template setPrimesCount(0);
							this->template countPrimesInNode();
						}
					}
				} else {
					performanceTimer.stop();
					if (_countNumberOfPrimesOnNode) {
						this->template countPrimesInNode();
					}
				}
			} else {
				this->template sendFinishSievingMessageToRoot();

				if (_sendPrimesCountToRoot) {
					this->template sendPrimesCountToRoot();
				} else if (_countNumberOfPrimesOnNode) {
					this->template countPrimesInNode();
				}

				if (_sendResultsToRoot) {
					this->template sendResultsToRootProcess(maxRange);
				}
			}
		}

		virtual void sendFinishSievingMessageToRoot() {
			double elapsedTimeMicroSec = this->template getPerformanceTimer().getElapsedTimeInMicroSec();
			MPI_Send(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 0, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD);
		}

		virtual void waitForAllProcessesInGroupToFinishSieving() {
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				double elapsedTimeMicroSec;
				MPI_Status status;
				MPI_Recv(&elapsedTimeMicroSec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR != MPI_SUCCESS) {
					cout << "    --> MPI_Recv detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
		}

		virtual size_t countPrimesInNode() {
			PerformanceTimer countingPrimesTimer;
			countingPrimesTimer.reset();
			countingPrimesTimer.start();

			size_t startPossiblePrime = this->template getStartSieveNumber();
			size_t maxRange = this->template getMaxRange();
			cout << "    --> Process with rank " << _processID << " counting primes in [" << startPossiblePrime << ", " << maxRange << "]" << endl;

			// compute primes in this process
			PerformanceTimer performanceTimer;
			performanceTimer.reset();
			performanceTimer.start();
			size_t numberPrimesFound = this->template getNumberPrimesFound();
			performanceTimer.stop();
			cout << "    --> Process " << _processID << " counted " << numberPrimesFound << " primes in [" << startPossiblePrime << ", " << maxRange << "] in " << performanceTimer.getElapsedTimeFormated() << endl;
			return numberPrimesFound;
		}

		virtual void sendPrimesCountToRoot() {
			unsigned long long numberPrimesFound = (unsigned long long)this->template countPrimesInNode();
			MPI_Send(&numberPrimesFound, 1, MPI_LONG_LONG, 0, MSG_NODE_PRIMES_FOUND_COUNT, MPI_COMM_WORLD);
		}

		virtual void collectPrimesCountFromProcessGroup() {
			PerformanceTimer countingPrimesTimer;
			countingPrimesTimer.reset();
			countingPrimesTimer.start();

			size_t numberPrimesFound = this->template countPrimesInNode();

			// update primes count with the partial count from the remaining processes
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				unsigned long long primesCount;
				MPI_Status status;
				MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_PRIMES_FOUND_COUNT, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					this->template incrementPrimesCount((size_t) primesCount);
				} else {
					cout << "    --> MPI_Recv detected the following error code: " << status.MPI_ERROR << endl;
				}
			}

			numberPrimesFound = this->template getNumberPrimesFound(); // get updated value
			cout << "\n    >>>>> Counted " << numberPrimesFound << " primes on " << _numberProcesses << " processes in " << countingPrimesTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;
			countingPrimesTimer.stop();
		}

		virtual void sendResultsToRootProcess(size_t maxRange) {
			cout << "    > Sending results from process with rank " << _processID << " to root process..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t blockSize = this->template getProcessBitsetSize(_processID, _numberProcesses, maxRange);

			MPI_Send(&primesBitset[0], blockSize, MPI_UNSIGNED_CHAR, 0, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD);
			cout << "    --> Finished sending results from process with rank " << _processID << " to root process" << endl;
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
			cout << "\n    > Collecting results from other processes..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				cout << "    > Probing for results..." << endl;
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, _numberProcesses, maxRange);
					size_t blockSize = this->template getProcessBitsetSize(_processID, _numberProcesses, maxRange);
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					cout << "    --> Collecting results from process with rank " << processID << endl;
					MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD, &status);
					cout << "    --> Finished collecting results from process with rank " << processID << endl;
				} else {
					cout << "    --> MPI_Probe detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
			cout << "    --> Finished collecting all results\n" << endl;
		}

		virtual bool savePrimesToFile(string filename) {
			if (!_sendResultsToRoot && _numberProcesses > 1) {
				stringstream sstream;
				size_t startPossiblePrime = this->template getStartSieveNumber();
				size_t maxRange = this->template getMaxRange();
				sstream << filename.substr(0, filename.rfind("."));
				sstream << "_primes_in[" << startPossiblePrime << ", " << maxRange << "]";
				sstream << filename.substr(filename.rfind("."));

				filename = sstream.str();
			}

			ofstream outputStream(filename.c_str());

			if (outputStream.is_open()) {
				if (!_sendResultsToRoot && _numberProcesses > 1) {
					cout << "\n    > Exporting partial results to file " << filename << "...";
				}
				this-> template savePrimes(outputStream);
				if (!_sendResultsToRoot && _numberProcesses > 1) {
					cout << "\n    --> Export partial results to file " << filename << " finished!" << endl;
				}
				return true;
			}

			return false;
		}

//		virtual bool checkPrimesFromFile(string filename) {
//			WheelType& wheelSieve = this->template getWheelSieve();
//			FlagsContainer& primesBitset = this->template getPrimesBitset();
//			size_t startSieveNumber = this->template getStartSieveNumber();
//			return wheelSieve.checkPrimesFromFileMPI(primesBitset, startSieveNumber, filename);
//		}

		void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
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

		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) = 0;

		virtual void initPrimesBitSetSizeForRoot(size_t maxRange, size_t maxNumberToStore) = 0;
		virtual void initPrimesBitSetSizeForRootWithAllValues(size_t maxRange) = 0;
		virtual void initPrimesBitSetSizeForSievingPrimes(size_t maxRangeSquareRoot) = 0;
		virtual void initPrimesBitSetSizeForSieving(size_t blockSize) = 0;
		virtual WheelType& getWheelSieve() = 0;

		virtual size_t getBitsetPositionToNumberMPI(size_t number) = 0;
		virtual size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) = 0;
		virtual size_t getProcessBitsetSize(size_t processID, size_t numberProcesses, size_t maxRange) = 0;

		virtual inline size_t getBlockBeginNumber() {
			return WheelType().getFirstPrimeToSieve();
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

		inline bool isSendPrimesCountToRoot() const {
			return _sendPrimesCountToRoot;
		}

		inline void setSendPrimesCountToRoot(bool sendPrimesCountToRoot) {
			_sendPrimesCountToRoot = sendPrimesCountToRoot;
		}

		inline bool isSendResultsToRoot() const {
			return _sendResultsToRoot;
		}

		inline void setSendResultsToRoot(bool sendResultsToRoot) {
			_sendResultsToRoot = sendResultsToRoot;
		}

		inline bool isCountNumberOfPrimesOnNode() const {
			return _countNumberOfPrimesOnNode;
		}

		inline void setCountNumberOfPrimesOnNode(bool countNumberOfPrimesOnNode) {
			_countNumberOfPrimesOnNode = countNumberOfPrimesOnNode;
		}

		inline void setNumberProcesses(int numberProcesses) {
			_numberProcesses = numberProcesses;
		}

		inline vector<size_t>& getSievingPrimes() {
			return _sievingPrimes;
		}

		inline void setSievingPrimes(const vector<size_t>& sievingPrimes) {
			_sievingPrimes = sievingPrimes;
		}
	};

