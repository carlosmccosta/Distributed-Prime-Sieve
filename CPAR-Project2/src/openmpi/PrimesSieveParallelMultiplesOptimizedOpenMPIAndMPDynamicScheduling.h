#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP.h"

#include <vector>
#include <utility>
#include <string>

using std::string;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling: public PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType> {
	protected:
		size_t _processStartBlockNumber;
		size_t _processEndBlockNumber;
		string _outputResultsFilename;
		size_t _dynamicSchedulingSegmentSizeInElements;

		unsigned long long _numberSegmentsSieved;
		size_t _numberSegmentsCreated;
		vector<vector<pair<size_t, size_t> > > _segmentsDistribution;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, string outputResultsFilename = "") :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot), _processStartBlockNumber(
						0), _processEndBlockNumber(0), _outputResultsFilename(outputResultsFilename), _dynamicSchedulingSegmentSizeInElements(1048576), _numberSegmentsSieved(0), _numberSegmentsCreated(0) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling() {
		}

		virtual void computePrimes(size_t maxRange) {
			size_t numberProcesses = (size_t) this->template getNumberProcesses();

			if (numberProcesses < 2) {
				cerr << "\n   !!!!! This algorithm must have a MPI_Comm_size of at least 2 processes !!!!!\n" << endl;
				exit(EXIT_FAILURE);
			}

			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			bool _sendResultsToRoot = this->template isSendResultsToRoot();
			bool _countNumberOfPrimesOnNode = this->template isCountNumberOfPrimesOnNode();
			bool _sendPrimesCountToRoot = this->template isSendPrimesCountToRoot();

			performanceTimer.reset();
			performanceTimer.start();

			this->template setPrimesCount(0);
			this->template clearPrimesValues();

			int _processID = this->template getProcessId();
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			// init primes used in sieving
			this->template setStartSieveNumber(this->template getBlockBeginNumber());
			this->template setMaxRange(maxRangeSquareRoot);

			if (_processID == 0) {
				if (_sendResultsToRoot) {
					this->template setMaxRange(maxRange);
					this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
				} else {
					this->template initPrimesBitSetSizeForRoot(maxRangeSquareRoot, maxRangeSquareRoot);
				}
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			}

			if (_processID == 0) {
				this->template setupSegmentDistribution(maxRange);
				this->template setupRootSegmentManagementTasks(maxRangeSquareRoot);
			} else {
				vector<pair<size_t, size_t> > sievingMultiples;
				this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);

				size_t previousProcessEndBlockNumber = 0;
				// sieve all assigned blocks
				while (this->template getNewSegmentFromRoot()) {
					this->template initPrimesBitSetSizeForSieving(_processEndBlockNumber - _processStartBlockNumber);

					// to avoid computing sieving multiples when the segment are contiguous
					if (previousProcessEndBlockNumber != _processStartBlockNumber) {
						this->template computeSievingMultiples(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);
						previousProcessEndBlockNumber = _processEndBlockNumber;
					}

					this->template removeComposites(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);
					++_numberSegmentsSieved;

					if (_sendPrimesCountToRoot) {
						this->template sendPrimesCountToRoot();
					} else if (_countNumberOfPrimesOnNode) {
						this->template countPrimesInNode();
					}

					if (_sendResultsToRoot) {
						this->template sendSievedBlockToRoot();
					} else {
						this->template outputResults();
					}

					cout << "    -> Finish sieving segment [" << _processStartBlockNumber << ", " << _processEndBlockNumber << "]" << endl;
				}

				performanceTimer.stop();
				cout << "    > Finish sieving " << _numberSegmentsSieved << " segments of " << _dynamicSchedulingSegmentSizeInElements << " numbers each, in process with rank " << _processID << " in "
						<< performanceTimer.getElapsedTimeFormated() << endl;
			}
		}

		void calculateSievingMultiples(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			int _processID = this->template getProcessId();

			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.reset();
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);
			computeSievingPrimesTimer.stop();
			size_t startSievingNumber = this->template getStartSieveNumber();
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in [" << startSievingNumber << ", " << maxRangeSquareRoot << "], on process with rank " << _processID << " in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;
		}

		virtual void setupSegmentDistribution(size_t maxRange) {
			size_t numberProcesses = (size_t) this->template getNumberProcesses() - 1; // process 0 acts only as a management node
			_segmentsDistribution = vector<vector<pair<size_t, size_t> > >(numberProcesses, vector<pair<size_t, size_t> >());

			size_t dynamicSchedulingBlockSizeInElements = _dynamicSchedulingSegmentSizeInElements;
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			size_t numberBlocks = ceil((double) (maxRange - maxRangeSquareRoot) / (double) dynamicSchedulingBlockSizeInElements);
			size_t numberBlocksPerProcess = numberBlocks / numberProcesses;

			size_t numberBlocksAdded = 0;

			size_t blockNumber = 0;

			#pragma omp parallel for \
				default(shared) \
				firstprivate(numberProcesses, numberBlocksPerProcess, dynamicSchedulingBlockSizeInElements, maxRangeSquareRoot) \
				reduction(+: numberBlocksAdded)
			for (size_t processPosition = 0; processPosition < numberProcesses; ++processPosition) {
				for (size_t processBlockNumber = 0; processBlockNumber < numberBlocksPerProcess; ++processBlockNumber) {
					size_t blockBeginNumber = (processPosition * numberBlocksPerProcess + blockNumber) * dynamicSchedulingBlockSizeInElements + maxRangeSquareRoot;
					size_t blockEndNumber = blockBeginNumber + dynamicSchedulingBlockSizeInElements;

					if (blockBeginNumber % 2 == 0) {
						++blockBeginNumber;
					}

					if (blockEndNumber >= maxRange) {
						blockEndNumber = maxRange + 1;
					}
					_segmentsDistribution[processPosition].push_back(pair<size_t, size_t>(blockBeginNumber, blockEndNumber));
					++numberBlocksAdded;
				}
			}

			// add any missing blocks to last process because of division rounding
			size_t lastProcessNumberIndex = numberProcesses - 1;
			for (; numberBlocksAdded < numberBlocks; ++numberBlocksAdded) {
				size_t blockBegin = numberBlocksAdded * dynamicSchedulingBlockSizeInElements + maxRangeSquareRoot;
				size_t blockEnd = blockBegin + dynamicSchedulingBlockSizeInElements;

				if (blockEnd >= maxRange) {
					blockEnd = maxRange + 1;
				}
				_segmentsDistribution[lastProcessNumberIndex].push_back(pair<size_t, size_t>(blockBegin, blockEnd));
				++numberBlocksAdded;
			}

			_numberSegmentsCreated = numberBlocksAdded;
		}

		virtual void setupRootSegmentManagementTasks(size_t maxRangeSquareRoot) {
			#pragma omp parallel default(shared)
			{
				#pragma omp single nowait
				{
					#pragma omp task default(shared)
					this->template startServerOfSegmentsAvailableTask();

					if (this->template isSendPrimesCountToRoot()) {
						#pragma omp task default(shared)
						{
							vector<pair<size_t, size_t> > sievingMultiples;
							this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);

							this->template startCollectorOfNodesPrimesCountTask();
						}
					}

					if (this->template isSendResultsToRoot()) {
						#pragma omp task default(shared)
						this->template startCollectorOfSievingDataResultsTask();
					}
				}
			}
		}

		virtual bool getNewSegmentFromRoot() {
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

		virtual void startServerOfSegmentsAvailableTask() {
			cout << "    --> Root server task for available segments started..." << endl;
			size_t numberSegmentsRemaining = _numberSegmentsCreated;
			while (numberSegmentsRemaining > 0) {
				MPI_Status status;
				unsigned long long numberSegmentsSievedByNode;
				MPI_Recv(&numberSegmentsSievedByNode, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD, &status);

				unsigned long long segmentRange[2];
				if (status.MPI_ERROR == MPI_SUCCESS) {
					pair<size_t, size_t> newSegmentRange;
					size_t processPositionToTakeSegment = status.MPI_SOURCE - 1; // process 0 isn't present in segment distribution because only acts as manager

					if (_segmentsDistribution[processPositionToTakeSegment].empty()) {
						processPositionToTakeSegment = this->template getProcessPositionWithMoreSegmentsRemaining();
					}

					newSegmentRange = _segmentsDistribution[processPositionToTakeSegment].back();
					_segmentsDistribution[processPositionToTakeSegment].pop_back();

					segmentRange[0] = newSegmentRange.first;
					segmentRange[1] = newSegmentRange.second;

					MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD);
					--numberSegmentsRemaining;
				} else {
					cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving request for a new segment from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				}
			}

			this->template endSievingNodes();
		}

		virtual void endSievingNodes() {
			cout << "    --> All segments processed. Terminating tasks..." << endl;
			size_t numberProcessesToTerminate = (size_t) this->template getNumberProcesses() - 1;
			for (size_t numberNodesTerminated = 0; numberNodesTerminated < numberProcessesToTerminate; ++numberNodesTerminated) {
				MPI_Status status;
				unsigned long long numberSegmentsSievedByNode;
				MPI_Recv(&numberSegmentsSievedByNode, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					unsigned long long segmentRange[2];
					segmentRange[0] = 0;
					segmentRange[1] = 0;

					MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD);
				} else {
					cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving request for a new segment from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				}
			}

			size_t numberProcesses = (size_t) this->template getNumberProcesses();
			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			cout << "\n    >>>>> Finish sieving in all " << numberProcesses << " processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;
		}

		size_t getProcessPositionWithMoreSegmentsRemaining() {
			size_t processPositionToTakeSegment = 0;
			size_t maxNumberBlocksOnProcess = _segmentsDistribution[0].size();

			size_t numberProcessPositions = _segmentsDistribution.size();
			for (size_t processPosition = 1; processPosition < numberProcessPositions; ++processPosition) {
				size_t currentProcessBlockSize = _segmentsDistribution[processPosition].size();
				if (currentProcessBlockSize > maxNumberBlocksOnProcess) {
					maxNumberBlocksOnProcess = currentProcessBlockSize;
					processPositionToTakeSegment = processPosition;
				}
			}

			return processPositionToTakeSegment;
		}

		virtual void sendSievedBlockToRoot() {
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t blockSize = _processEndBlockNumber - _processStartBlockNumber;

			// send results range
			unsigned long long segmentRange[2];
			segmentRange[0] = (unsigned long long) _processStartBlockNumber;
			segmentRange[1] = (unsigned long long) _processEndBlockNumber;
			MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, 0, MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE, MPI_COMM_WORLD);

			// send block sieving data
			this->template sendSievingDataMPI(primesBitset, 0, blockSize, 0, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
		}

		virtual void receiveSievedBlockFromNode(int nodeSource) {
			// received block range
			unsigned long long segmentRange[2];
			MPI_Status status;
			MPI_Recv(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, nodeSource, MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
				size_t blockSize = this->template getNumberBitsToStoreBlock(segmentRange[1] - segmentRange[0]);
				size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(segmentRange[0]);

				// receive sieving data
				FlagsContainer& primesBitset = this->template getPrimesBitset();
				this->template receiveSievingDataMPI(primesBitset, positionToStoreResults, blockSize, nodeSource, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
			} else {
				cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving sieving segment results data from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
			}
		}

		virtual void startCollectorOfSievingDataResultsTask() {
			cout << "    --> Root task for collecting nodes segments results started..." << endl;
			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				#pragma omp task default(shared)
				{
					MPI_Status status;
					MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);
					if (status.MPI_ERROR == MPI_SUCCESS) {
						this->template receiveSievedBlockFromNode(status.MPI_SOURCE);
					} else {
						cout << "    --> MPI_Probe detected the following error code when probing for segment sieving data results from node " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
					}
				}
			}
		}

		virtual void startCollectorOfNodesPrimesCountTask() {
			cout << "    --> Root task for collecting nodes segments primes count started..." << endl;

			// root node only has sieving primes
			size_t numberSievingPrimes = this->template getSievingPrimes().size();
			size_t numberPrimesSievedByTheWheel = this->template getWheelSieve().getNumberPrimesSievedByTheWheel();
			size_t numberPrimesFound = numberSievingPrimes + numberPrimesSievedByTheWheel;
			this->template setPrimesCount(numberPrimesFound);

			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				unsigned long long primesCount;
				MPI_Status status;
				MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_PRIMES_FOUND_COUNT, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					this->template incrementPrimesCount((size_t) primesCount);
				} else {
					cout << "    --> MPI_Recv detected the following error code when receiving segment primes count from " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
				}
			}

			int _numberProcesses = this->template getNumberProcesses();
			numberPrimesFound = this->template getPrimesCount(); // get updated value
			cout << "\n    >>>>> Counted " << numberPrimesFound << " primes on " << _numberProcesses << " processes <<<<<\n" << endl;
		}

		virtual bool outputResults() {
			if (_outputResultsFilename == "stdout") {
				cout << "\n\n=============================================  Computed primes  =============================================\n\n";
				this->template printPrimesToConsole();
				cout << "\n" << endl;
				return true;
			} else if (_outputResultsFilename != "") {
				bool sendResultsToRoot = this->template isSendResultsToRoot();

				if (sendResultsToRoot) {
					cout << "\n    > Exporting results to file " << _outputResultsFilename << "..." << endl;
				}

				PerformanceTimer performanceTimer;
				performanceTimer.reset();
				performanceTimer.start();
				if (this->template savePrimesToFile(_outputResultsFilename)) {
					performanceTimer.stop();

					if (sendResultsToRoot) {
						cout << "    --> Export to file " << _outputResultsFilename << " finished in " << performanceTimer.getElapsedTimeFormated() << endl;
					}
					return true;
				} else {
					cerr << "    !!!!! Export to file " << _outputResultsFilename << "failed !!!!!" << endl;
				}
			}
			return false;
		}
};

