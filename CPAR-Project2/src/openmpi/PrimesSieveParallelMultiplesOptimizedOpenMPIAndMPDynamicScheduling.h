#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP.h"
#include "../lib/Settings.h"

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
		size_t _dynamicSchedulingNumberSegments;

		unsigned long long _numberSegmentsSieved;
		size_t _numberSegmentsCreated;
		vector<vector<pair<size_t, size_t> > > _segmentsDistribution;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, size_t dynamicSchedulingSegmentSizeInElements = 1048576, size_t dynamicSchedulingNumberSegments = 0,
				string outputResultsFilename = "") :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode,
						sendPrimesCountToRoot), _processStartBlockNumber(0), _processEndBlockNumber(0), _outputResultsFilename(outputResultsFilename), _dynamicSchedulingSegmentSizeInElements(
						dynamicSchedulingSegmentSizeInElements), _dynamicSchedulingNumberSegments(dynamicSchedulingNumberSegments), _numberSegmentsSieved(0), _numberSegmentsCreated(0) {
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

			if (_dynamicSchedulingNumberSegments != 0) {
				_dynamicSchedulingSegmentSizeInElements = maxRange / _dynamicSchedulingNumberSegments;
			} else {
				_dynamicSchedulingNumberSegments = ceil((double) (maxRange - maxRangeSquareRoot) / (double) _dynamicSchedulingSegmentSizeInElements);
			}

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

				this->template setupSegmentDistribution(maxRange);
				this->template setupRootSegmentManagementTasks(maxRangeSquareRoot);
//				this->template setupRootSegmentManagementSections(maxRangeSquareRoot);
			} else {
				this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);

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

#ifdef DEBUG_OUTPUT
					cout << "    -> Finish sieving segment [" << _processStartBlockNumber << ", " << (_processEndBlockNumber - 1) << "]" << endl;
#endif
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
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in [" << startSievingNumber << ", " << maxRangeSquareRoot << "], on process with rank " << _processID << " in "
					<< computeSievingPrimesTimer.getElapsedTimeFormated() << endl;
		}

		virtual void setupSegmentDistribution(size_t maxRange) {
			size_t numberProcesses = (size_t) this->template getNumberProcesses() - 1; // process 0 acts only as a management node
			_segmentsDistribution = vector<vector<pair<size_t, size_t> > >(numberProcesses, vector<pair<size_t, size_t> >());

			size_t maxRangeSquareRootOffset = (size_t) sqrt(maxRange) + 1;
			size_t numberSegmentsPerProcess = _dynamicSchedulingNumberSegments / numberProcesses;

			size_t numbersegmentsAdded = 0;
			size_t dynamicSchedulingBlockSizeInElements = _dynamicSchedulingSegmentSizeInElements;

			/*#pragma omp parallel for \
				default(shared) \
				firstprivate(numberProcesses, numberBlocksPerProcess, dynamicSchedulingBlockSizeInElements, maxRangeSquareRoot) \
				reduction(+: numberBlocksAdded)*/
			for (size_t processPosition = 0; processPosition < numberProcesses; ++processPosition) {
				for (size_t processSegmentNumber = 0; processSegmentNumber < numberSegmentsPerProcess; ++processSegmentNumber) {
					size_t segmentBeginNumber = (processPosition * numberSegmentsPerProcess + processSegmentNumber) * dynamicSchedulingBlockSizeInElements + maxRangeSquareRootOffset;
					size_t segmentEndNumber = segmentBeginNumber + dynamicSchedulingBlockSizeInElements;

					if (segmentBeginNumber % 2 == 0) {
						++segmentBeginNumber;
					}

					if (segmentEndNumber >= maxRange) {
						segmentEndNumber = maxRange + 1;
					}
					_segmentsDistribution[processPosition].push_back(pair<size_t, size_t>(segmentBeginNumber, segmentEndNumber));
#ifdef DEBUG_OUTPUT
					cout << "    --> Added segment [" << segmentBeginNumber << ", " << (segmentEndNumber - 1) << "] to list of probable segments to process " << (processPosition + 1) << endl;
#endif
					++numbersegmentsAdded;
				}
			}

			// add any missing blocks to last process because of division rounding
			size_t lastProcessNumberIndex = numberProcesses - 1;
			for (; numbersegmentsAdded < _dynamicSchedulingNumberSegments; ++numbersegmentsAdded) {
				size_t segmentBeginNumber = numbersegmentsAdded * dynamicSchedulingBlockSizeInElements + maxRangeSquareRootOffset;
				size_t segmentEndNumber = segmentBeginNumber + dynamicSchedulingBlockSizeInElements;

				if (segmentEndNumber >= maxRange) {
					segmentEndNumber = maxRange + 1;
				}
				_segmentsDistribution[lastProcessNumberIndex].push_back(pair<size_t, size_t>(segmentBeginNumber, segmentEndNumber));
#ifdef DEBUG_OUTPUT
				cout << "    --> Added segment [" << segmentBeginNumber << ", " << (segmentEndNumber - 1) << "] to list of probable segments to process " << (lastProcessNumberIndex + 1) << endl;
#endif
			}

			_numberSegmentsCreated = numbersegmentsAdded;
			cout << "    --> Root process initialized " << _numberSegmentsCreated << " segments to sieve" << endl;
		}

		virtual void setupRootSegmentManagementTasks(size_t maxRangeSquareRoot) {
#pragma omp parallel default(shared)
			{
#pragma omp single nowait
				{
#pragma omp task default(shared)
					this->template startServerOfSegmentsAvailable();

					if (this->template isSendPrimesCountToRoot()) {
#pragma omp task default(shared)
						{
							vector<pair<size_t, size_t> > sievingMultiples;
							this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);

							this->template startCollectorOfNodesPrimesCount();
						}
					}

					if (this->template isSendResultsToRoot()) {
#pragma omp task default(shared)
						this->template startCollectorOfSievingDataResultsTask();
					}
				}

#pragma omp taskwait
			}
		}

		virtual void setupRootSegmentManagementSections(size_t maxRangeSquareRoot) {
#pragma omp parallel sections default(shared)
			{
#pragma omp section
				{
					this->template startServerOfSegmentsAvailable();
				}

#pragma omp section
				{
					if (this->template isSendPrimesCountToRoot()) {
						vector<pair<size_t, size_t> > sievingMultiples;
						this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);
						this->template startCollectorOfNodesPrimesCount();
					}
				}

#pragma omp section
				{
					if (this->template isSendResultsToRoot()) {
						this->template startCollectorOfSievingDataResults();
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

#ifdef DEBUG_OUTPUT
				int processID = this->template getProcessId();
				cout << "    --> Process " << processID << " received segment [" << segmentRange[0] << ", " << (segmentRange[1] - 1) << "]" << endl;
#endif
				return true;
			} else {
				return false;
			}
		}

		virtual void startServerOfSegmentsAvailable() {
			cout << "    --> Root server for available segments started..." << endl;
			size_t numberSegmentsRemaining = _numberSegmentsCreated;
			while (numberSegmentsRemaining > 0) {
				MPI_Status status;
				unsigned long long numberSegmentsSievedByNode;
				MPI_Recv(&numberSegmentsSievedByNode, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					pair<size_t, size_t> newSegmentRange;
					size_t processPositionToTakeSegment = status.MPI_SOURCE - 1; // process 0 isn't present in segment distribution because only acts as manager

					if (_segmentsDistribution[processPositionToTakeSegment].empty() || processPositionToTakeSegment >= _segmentsDistribution.size()) {
						processPositionToTakeSegment = this->template getProcessPositionWithMoreSegmentsRemaining();
					}

					newSegmentRange = _segmentsDistribution[processPositionToTakeSegment].back();
					_segmentsDistribution[processPositionToTakeSegment].pop_back();

					unsigned long long segmentRange[2];
					segmentRange[0] = newSegmentRange.first;
					segmentRange[1] = newSegmentRange.second;

#ifdef DEBUG_OUTPUT
					cout << "    --> Assigning segment [" << segmentRange[0] << ", " << (segmentRange[1] - 1) << "] to process " << status.MPI_SOURCE << endl;
#endif
					MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD);
					--numberSegmentsRemaining;
				} else {
					cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving request for a new segment from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				}

#ifdef DEBUG_OUTPUT
				cout << "    --> " << numberSegmentsRemaining << " remaining segments" << endl;
#endif
			}

			this->template endSievingNodes();
		}

		virtual void endSievingNodes() {
			cout << "    --> All segments processed. Terminating auxiliary nodes..." << endl;
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
			size_t blockSize = this->template getNumberBitsToStoreBlock(_processEndBlockNumber - _processStartBlockNumber);

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

		virtual void startCollectorOfSievingDataResults() {
			cout << "    --> Root collector for nodes segments results started..." << endl;
			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					this->template receiveSievedBlockFromNode(status.MPI_SOURCE);
				} else {
					cout << "    --> MPI_Probe detected the following error code when probing for segment sieving data results from node " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
				}

			}
		}

		virtual void startCollectorOfSievingDataResultsTask() {
			cout << "    --> Root collector for nodes segments results started..." << endl;

#pragma omp for \
				schedule(dynamic)
			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					this->template receiveSievedBlockFromNode(status.MPI_SOURCE);
				} else {
					cout << "    --> MPI_Probe detected the following error code when probing for segment sieving data results from node " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
				}
			}
		}

		virtual void startCollectorOfNodesPrimesCount() {
			cout << "    --> Root collector for nodes segments primes count started..." << endl;

			// root node only has sieving primes
			size_t numberSievingPrimes = this->template getSievingPrimes().size();
			size_t numberPrimesSievedByTheWheel = this->template getWheelSieve().getNumberPrimesSievedByTheWheel();
			size_t numberPrimesFound = numberSievingPrimes + numberPrimesSievedByTheWheel;
			this->template setPrimesCount(numberPrimesFound);

			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				this->template receivePrimesCountFromNode(MPI_ANY_SOURCE);
			}

			int _numberProcesses = this->template getNumberProcesses();
			numberPrimesFound = this->template getPrimesCount(); // get updated value
			cout << "\n    >>>>> Counted " << numberPrimesFound << " primes on " << _numberProcesses << " processes <<<<<\n" << endl;
		}

		virtual size_t receivePrimesCountFromNode(int nodeID) {
			unsigned long long primesCount;
			MPI_Status status;
			MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_PRIMES_COUNT_FOUND, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
#pragma omp critical
				this->template incrementPrimesCount((size_t) primesCount);
			} else {
				cout << "    --> MPI_Recv detected the following error code when receiving segment primes count from " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
			}

			return (size_t) primesCount;
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

