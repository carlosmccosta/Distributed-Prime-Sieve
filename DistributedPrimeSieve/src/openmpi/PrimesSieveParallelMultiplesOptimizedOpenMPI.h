#pragma once

#include "../PrimesSieve.h"
#include "../WheelFactorization.h"
#include "../lib/PrimesUtils.h"

#include "../lib/Settings.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <mpi.h>

using std::endl;
using std::sqrt;
using std::min;
using std::max;
using std::pair;

enum MessageTags {
	MSG_NODE_SIEVING_FINISHED = 0,
	MSG_NODE_PRIMES_COUNT_FOUND,
	MSG_NODE_COMPUTATION_RESULTS_SEGMENT,
	MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE,
	MSG_REQUEST_NEW_SEGMET,
	MSG_ASSIGN_NEW_SEGMET
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
		WheelType _wheelSieve;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPI(size_t maxRange, size_t blockSizeInElements = 16 * 1024, bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true,
				bool sendPrimesCountToRoot = true) :
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


		virtual void collectResultsFromProcessGroup(size_t maxRange) = 0;
		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) = 0;
		virtual size_t getNumberBitsToStoreBlock(size_t blockSize) = 0;


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
					PerformanceTimer memoryAllocationTimer;
					memoryAllocationTimer.start();
					this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
					memoryAllocationTimer.stop();
					cout << "    --> Allocated memory for " << maxRange << " numbers in " << memoryAllocationTimer.getElapsedTimeFormated() << endl;
				} else {
					this->template initPrimesBitSetSizeForRoot(maxRange, processEndBlockNumber - 1);
				}
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}

			// compute sieving primes
			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);
			computeSievingPrimesTimer.stop();
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in process with rank " << _processID << " in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;

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

			cout << "    > Finished sieving in process with rank " << _processID << " in " << totalPerformanceTimer.getElapsedTimeFormated() << endl;

			this->template syncProcesses(maxRange);
		}


		void syncProcesses(size_t maxRange) {
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();

			if (_processID == 0) {
				if (_numberProcesses > 1) {
					this->template waitForAllProcessesInGroupToFinishSieving();
					performanceTimer.stop();
					cout << "\n    >>>>> Finished sieving in all " << _numberProcesses << " processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;

					if (_sendPrimesCountToRoot) {
						this->template collectPrimesCountFromProcessGroup();
					} else if (_countNumberOfPrimesOnNode) {
						this->template countPrimesInNode();
					}

					if (_sendResultsToRoot) {
						this->template setMaxRange(maxRange);
						this->template collectResultsFromProcessGroup(maxRange);

						if (_countNumberOfPrimesOnNode && !_sendPrimesCountToRoot) {
							//force recount with results from other processes
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
					this->template sendPrimesCountToCollector(0);
				} else if (_countNumberOfPrimesOnNode) {
					this->template countPrimesInNode();
				}

				if (_sendResultsToRoot) {
					this->template sendResultsToRootProcess(maxRange);
				}
			}
		}

		MPI_Status receiveSievingDataMPI(FlagsContainer& primesBitset, size_t positionToStoreResults, size_t blockSize, int source, int tag) {
			MPI_Status status;
			cerr << "\n\n!!!!!Missing implementation for this type of bitset container !!!!!" << endl << endl;
			return status;
		}

		void sendSievingDataMPI(FlagsContainer& primesBitset, size_t startPositionOfResults, size_t blockSize, int destination, int tag) {
			cerr << "\n\n!!!!!Missing implementation for this type of bitset container !!!!!" << endl << endl;
		}

		inline void sendFinishSievingMessageToRoot() {
			double elapsedTimeMicroSec = this->template getPerformanceTimer().getElapsedTimeInMicroSec();
			MPI_Send(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 0, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD);
		}

		void waitForAllProcessesInGroupToFinishSieving() {
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				double elapsedTimeMicroSec;
				MPI_Status status;
				MPI_Recv(&elapsedTimeMicroSec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR != MPI_SUCCESS) {
					cout << "    --> MPI_Recv detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
		}

		size_t countPrimesInNode() {
			this->template setPrimesCount(0);
#					ifdef DEBUG_OUTPUT
					size_t startPossiblePrime = this->template getStartSieveNumber();
					size_t maxRange = this->template getMaxRange();
					cout << "    --> Process with rank " << _processID << " counting primes in [" << startPossiblePrime << ", " << maxRange << "]" << endl;

					// compute primes in this process
					PerformanceTimer performanceTimer;
					performanceTimer.start();
#					endif

					size_t numberPrimesFound = this->template getNumberPrimesFound();

#					ifdef DEBUG_OUTPUT
					performanceTimer.stop();
					cout << "    --> Process " << _processID << " counted " << numberPrimesFound << " primes in [" << startPossiblePrime << ", " << maxRange << "] in " << performanceTimer.getElapsedTimeFormated() << endl;
#					endif
			return numberPrimesFound;
		}

		void sendPrimesCountToCollector(int sendPrimesCountToRoot) {
			unsigned long long numberPrimesFound = (unsigned long long)this->template countPrimesInNode();
#			ifdef DEBUG_OUTPUT
			cout << "    --> Sending " << numberPrimesFound << " primes count to node " << sendPrimesCountToRoot << " from " << _processID << endl;
#			endif
			MPI_Send(&numberPrimesFound, 1, MPI_LONG_LONG, sendPrimesCountToRoot, MSG_NODE_PRIMES_COUNT_FOUND, MPI_COMM_WORLD);
		}

		void collectPrimesCountFromProcessGroup() {
			PerformanceTimer countingPrimesTimer;
			countingPrimesTimer.start();

			size_t numberPrimesFound = this->template countPrimesInNode();

			// update primes count with the partial count from the remaining processes
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
				unsigned long long primesCount;
				MPI_Status status;
				MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_PRIMES_COUNT_FOUND, MPI_COMM_WORLD, &status);

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

		void sendResultsToRootProcess(size_t maxRange) {
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t blockSize = this->template getProcessBitsetSize(_processID, _numberProcesses, maxRange);

			this->template sendSievingDataMPI(primesBitset, 0, blockSize, 0, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
		}

		bool savePrimesToFile(string filename) {
			if (!_sendResultsToRoot && _numberProcesses > 1) {
				stringstream sstream;
				size_t startPossiblePrime = this->template getStartSieveNumber();

				if (startPossiblePrime == this->template getWheelSieve().getFirstPrimeToSieve()) {
					startPossiblePrime = 2;
				}

				size_t maxRange = this->template getMaxRange();
				sstream << filename.substr(0, filename.rfind("."));
				sstream << "_primes_in[" << startPossiblePrime << ", " << maxRange << "]";
				sstream << filename.substr(filename.rfind("."));

				filename = sstream.str();
			}

			ofstream outputStream(filename.c_str());

			if (outputStream.is_open()) {
#				ifdef DEBUG_OUTPUT
				int processID = this->template getProcessId();
				if (!_sendResultsToRoot && _numberProcesses > 1) {
					cout << "\n    --> Exporting partial results to file " << filename << " in process " << processID << "..." << endl;
				}
#				endif

				this-> template savePrimes(outputStream);

#				ifdef DEBUG_OUTPUT
				if (!_sendResultsToRoot && _numberProcesses > 1) {
					cout << "    --> Export partial results to file " << filename << " in process " << processID << " finished!" << endl;
				}
#				endif
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
			size_t sievingPrimesSize = _sievingPrimes.size();
			sievingMultiples.clear();
			sievingMultiples.reserve(sievingPrimesSize);
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

		void initPrimesBitSetSizeForRoot(size_t maxRange, size_t maxNumberToStore) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxNumberToStore));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));

			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		void initPrimesBitSetSizeForRootWithAllValues(size_t maxRange) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRange));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));

			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		void initPrimesBitSetSizeForSievingPrimes(size_t maxRangeSquareRoot) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRangeSquareRoot));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange(maxRangeSquareRoot);

			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		void initPrimesBitSetSizeForSieving(size_t blockSize) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStoreBlock(blockSize));
		}

		size_t getProcessBitsetSize(size_t processID, size_t numberProcesses, size_t maxRange) {
			size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, numberProcesses, maxRange);
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(processID, numberProcesses, maxRange);

			if (processStartBlockNumber % 2 == 0) {
				++processStartBlockNumber;
			}

			if (processID == numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}
			return this->template getNumberBitsToStoreBlock(processEndBlockNumber - processStartBlockNumber);
		}

		virtual inline size_t getNumberBitsToStore(size_t maxRange) {
			return (this->template getNumberBitsToStoreBlock((maxRange + 1) - this->template getStartSieveNumber()));
		}

		virtual inline int getProcessIDWithFirstPrimesBlock() {
			return 0;
		}

		virtual inline size_t getBlockBeginNumber() {
			return _wheelSieve.getFirstPrimeToSieve();
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

		inline void clearPrimesValues() {
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

		inline WheelType& getWheelSieve() {
			return _wheelSieve;
		}

		inline void setWheelSieve(WheelType wheelSieve) {
			this->_wheelSieve = wheelSieve;
		}
	};

