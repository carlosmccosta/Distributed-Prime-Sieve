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
#			ifdef DEBUG_OUTPUT
			cout << "    --> Removing composites in process with rank " << _processID << " in [" << processBeginBlockNumber << ", " << (processEndBlockNumber - 1) << "]" << endl;
#			endif

			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockEndNumber = -1;

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
				} else if (sievingMultiples.empty() || blockBeginNumber != priviousBlockEndNumber) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockEndNumber = blockEndNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void syncProcesses(size_t maxRange) {
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

		virtual MPI_Status receiveSievingDataMPI(FlagsContainer& primesBitset, size_t positionToStoreResults, size_t blockSize, int source, int tag) {
			MPI_Status status;
			cerr << "\n\n!!!!!Missing implementation for this type of bitset container !!!!!" << endl << endl;
			return status;
		}

		virtual void sendSievingDataMPI(FlagsContainer& primesBitset, size_t startPositionOfResults, size_t blockSize, int destination, int tag) {
			cerr << "\n\n!!!!!Missing implementation for this type of bitset container !!!!!" << endl << endl;
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

		virtual void sendPrimesCountToCollector(int sendPrimesCountToRoot) {
			unsigned long long numberPrimesFound = (unsigned long long)this->template countPrimesInNode();
#			ifdef DEBUG_OUTPUT
			cout << "    --> Sending " << numberPrimesFound << " primes count to node " << sendPrimesCountToRoot << " from " << _processID << endl;
#			endif
			MPI_Send(&numberPrimesFound, 1, MPI_LONG_LONG, sendPrimesCountToRoot, MSG_NODE_PRIMES_COUNT_FOUND, MPI_COMM_WORLD);
		}

		virtual void collectPrimesCountFromProcessGroup() {
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

		virtual void sendResultsToRootProcess(size_t maxRange) {
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t blockSize = this->template getProcessBitsetSize(_processID, _numberProcesses, maxRange);

			this->template sendSievingDataMPI(primesBitset, 0, blockSize, 0, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
			cout << "\n    > Collecting results from other processes..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
#				ifdef DEBUG_OUTPUT
				cout << "    > Probing for results..." << endl;
#				endif
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, _numberProcesses, maxRange);
					if (processStartBlockNumber % 2 == 0) {
						++processStartBlockNumber;
					}

					size_t blockSize = this->template getProcessBitsetSize(processID, _numberProcesses, maxRange);
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					this->template receiveSievingDataMPI(primesBitset, positionToStoreResults, blockSize, status.MPI_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
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

		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t sievingMultiplesSize = sievingMultiples.size();
			for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
				pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					this->template setPrimesBitsetValueMPI(primeMultiple, true);
				}

				sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
			}
		}

		virtual void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
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
				if (!this->template getPrimesBitsetValueMPI(primeNumber)) {
					_sievingPrimes.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						this->template setPrimesBitsetValueMPI(compositeNumber, true);
					}
					sievingMultiples.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
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
				if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
		}

		virtual void savePrimes(ostream& outputStream) {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() <= 2) {
				size_t possiblePrime = this->template getStartSieveNumber();

				if (possiblePrime == this->template getWheelSieve().getFirstPrimeToSieve()) {
					outputStream << 2 << endl;
					outputStream << 3 << endl;
					outputStream << 5 << endl;
					outputStream << 7 << endl;
				}

				if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t maxRange = this->template getMaxRange();
				for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
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
				if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
					++primesFound;
				}
			}

			this->template setPrimesCount(primesFound);
			return primesFound;
		}

		virtual void initPrimesBitSetSizeForRoot(size_t maxRange, size_t maxNumberToStore) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxNumberToStore));

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
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRangeSquareRoot));

			size_t numberSievingPrimes = this->template getNumberOfPrimesInRange(maxRangeSquareRoot);

			_sievingPrimes.clear();
			_sievingPrimes.reserve(numberSievingPrimes);
		}

		virtual void initPrimesBitSetSizeForSieving(size_t blockSize) {
			this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStoreBlock(blockSize));
		}

		virtual inline size_t getProcessBitsetSize(size_t processID, size_t numberProcesses, size_t maxRange) {
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

		// bitset specific memory management
		virtual size_t getNumberBitsToStoreBlock(size_t blockSize) = 0;
		virtual size_t getBitsetPositionToNumberMPI(size_t number) = 0;
		virtual size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) = 0;
		// end bitset specific memory management

		virtual int getProcessIDWithFirstPrimesBlock() {
			return 0;
		}

		virtual inline size_t getBlockBeginNumber() {
			return _wheelSieve.getFirstPrimeToSieve();
		}

		inline bool getPrimesBitsetValueMPI(size_t number) {
			return this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)];
		}

		inline void setPrimesBitsetValueMPI(size_t number, bool newValue) {
			this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)] = newValue;
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

		inline WheelType& getWheelSieve() {
			return _wheelSieve;
		}

		inline void setWheelSieve(WheelType wheelSieve) {
			this->_wheelSieve = wheelSieve;
		}
	};

