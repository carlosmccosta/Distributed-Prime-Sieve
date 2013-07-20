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
		bool _outputOnlyLastSegment;
		size_t _dynamicSchedulingSegmentSizeInElements;
		size_t _dynamicSchedulingNumberSegments;

		int _mpiThreadSupport;
		int _processIDWithFirstPrimesBlock;
		size_t _numberSievingProcesses;

		unsigned long long _numberSegmentsSieved;
		size_t _numberSegmentsCreated;
		vector<vector<pair<size_t, size_t> > > _segmentsDistribution;

		size_t _globalMaxRange;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, size_t dynamicSchedulingSegmentSizeInElements = 1048576, size_t dynamicSchedulingNumberSegments = 0,
				string outputResultsFilename = "", bool outputOnlyLastSegment = false) :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType>(maxRange, blockSizeInElements, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode,
						sendPrimesCountToRoot), _processStartBlockNumber(0), _processEndBlockNumber(0), _outputResultsFilename(outputResultsFilename), _outputOnlyLastSegment(outputOnlyLastSegment), _dynamicSchedulingSegmentSizeInElements(
						dynamicSchedulingSegmentSizeInElements), _dynamicSchedulingNumberSegments(dynamicSchedulingNumberSegments), _mpiThreadSupport(MPI_THREAD_SINGLE), _processIDWithFirstPrimesBlock(
						0), _numberSievingProcesses(1), _numberSegmentsSieved(0), _numberSegmentsCreated(0), _globalMaxRange(11) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling() {
		}

		virtual void computePrimes(size_t maxRange) {
			if (!(this->template checkifThereIsSufficientNumberOfProcess())) {
				cerr << "!!!! Insufficient number of processes in work group !!!!!" << endl;
				cerr << "    --> Check with ompi_info if current MPI implementation supports MPI_THREAD_MULTIPLE" << endl;
				cerr << "    --> In ompi_info should appear the following line: Thread support: posix (mpi: yes, progress: no)" << endl;
				cerr << " If your implementation doesn't provide MPI_THREAD_MULTIPLE, it will be necessary:" << endl;
				cerr << "    --> 1 process to dynamic schedule the blocks for each sieving process" << endl;
				cerr << "    --> 1 Additional process if flag --sendPrimesCountToRoot is set" << endl;
				cerr << "    --> 1 Additional process if flag --sendResultsToRoot is set" << endl;
				cerr << "    --> The remaining process are used to sieve the primes segments\n" << endl;
				cerr << " Although the current implementation is prepared for systems with support for MPI_THREAD_MULTIPLE, since it isn't tested yet, it was decided to be left deactivated" << endl;
				cerr
						<< " Implementation tested in a Ubuntu 13.04 64 bits, with Open MPI: 1.4.5, Open MPI SVN revision: r25905, Open MPI release date: Feb 10, 2012, Thread support: posix (mpi: no, progress: no)\n\n"
						<< endl;

				exit(EXIT_FAILURE);
			}

			_globalMaxRange = maxRange;

			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			performanceTimer.reset();
			performanceTimer.start();

			this->template setPrimesCount(this->template getWheelSieve().getNumberPrimesSievedByTheWheel());
			this->template clearPrimesValues();

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			this->template setStartSieveNumber(this->template getBlockBeginNumber());
			this->template setMaxRange(maxRangeSquareRoot);

			if (_dynamicSchedulingNumberSegments != 0) {
				_dynamicSchedulingSegmentSizeInElements = ceil((double) maxRange / (double) _dynamicSchedulingNumberSegments);
			} else {
				_dynamicSchedulingNumberSegments = ceil((double) (maxRange - maxRangeSquareRoot) / (double) _dynamicSchedulingSegmentSizeInElements);
			}

			size_t numberProcesses = (size_t) this->template getNumberProcesses();
			int _processID = this->template getProcessId();
			bool _sendResultsToRoot = this->template isSendResultsToRoot();
			bool _countNumberOfPrimesOnNode = this->template isCountNumberOfPrimesOnNode();
			bool _sendPrimesCountToRoot = this->template isSendPrimesCountToRoot();

			_numberSievingProcesses = numberProcesses - 1;
			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				_processIDWithFirstPrimesBlock = 1;
				if (_sendResultsToRoot) {
					--_numberSievingProcesses;
				}
				if (_sendPrimesCountToRoot) {
					--_numberSievingProcesses;
				}

				if (!_sendResultsToRoot && _sendPrimesCountToRoot) {
					_processIDWithFirstPrimesBlock = 2;
				}
			}

			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				if (_processID == 0) {
					this->template setupSegmentDistribution(maxRange);
					this->template startServerOfSegmentsAvailable(_numberSegmentsCreated, _segmentsDistribution, _numberSievingProcesses);
				} else if (_sendResultsToRoot && _processID == 1) {
					this->template setMaxRange(maxRange);
					this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
					vector<pair<size_t, size_t> > sievingMultiples;
					this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);
					this->template startCollectorOfSievingDataResults(_dynamicSchedulingNumberSegments);
					this->template outputResults();
					double elapsedTimeMicroSec = this->template getPerformanceTimer().getElapsedTimeInMicroSec();
					MPI_Send(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 0, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD);
				} else if ((_sendResultsToRoot && _sendPrimesCountToRoot && _processID == 2) || (!_sendResultsToRoot && _sendPrimesCountToRoot && _processID == 1)) {
					this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
					vector<pair<size_t, size_t> > sievingMultiples;
					this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);
					this->template setPrimesCount(this->template startCollectorOfNodesPrimesCount(_dynamicSchedulingNumberSegments));
					double elapsedTimeMicroSec = this->template getPerformanceTimer().getElapsedTimeInMicroSec();
					MPI_Send(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 0, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD);
				} else {
					this->template startSievingSegmentsInSievingProcesses(maxRangeSquareRoot, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
				}
			} else {
				if (_processID == 0) {
					if (_sendResultsToRoot) {
						this->template setMaxRange(maxRange);
						this->template initPrimesBitSetSizeForRootWithAllValues(maxRange);
					} else {
						this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
					}

					this->template setupSegmentDistribution(maxRange);
//					this->template setupRootSegmentManagementTasks(maxRangeSquareRoot);
					this->template setupRootSegmentManagementSections(maxRangeSquareRoot, _numberSievingProcesses);
				} else {
					this->template startSievingSegmentsInSievingProcesses(maxRangeSquareRoot, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
				}
			}
		}

		bool checkifThereIsSufficientNumberOfProcess() {
			size_t numberProcesses = (size_t) this->template getNumberProcesses();

			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				size_t minimumNumberOfProcesses = 2;

				if (this->template isSendResultsToRoot()) {
					++minimumNumberOfProcesses;
				}

				if (this->template isSendPrimesCountToRoot()) {
					++minimumNumberOfProcesses;
				}

				return numberProcesses >= minimumNumberOfProcesses;
			} else {
				return numberProcesses >= 2;
			}
		}

		void startSievingSegmentsInSievingProcesses(size_t maxRangeSquareRoot, bool sendResultsToRoot, bool countNumberOfPrimesOnNode, bool sendPrimesCountToRoot) {
			this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);

			vector<pair<size_t, size_t> > sievingMultiples;
			this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);

			size_t previousProcessEndBlockNumber = 0;

			int processIDToSendResults = 0;
			int processIDToSendPrimesCount = 0;

			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				processIDToSendResults = 1;

				if (sendResultsToRoot) {
					processIDToSendPrimesCount = 2;
				} else {
					processIDToSendPrimesCount = 1;
				}

			}

			// save sieving primes to partial file
			if (!sendResultsToRoot) {
				this->template outputResults();
			}

			// sieve all assigned blocks
			while (this->template getNewSegmentFromRoot()) {
				this->template initPrimesBitSetSizeForSieving(_processEndBlockNumber - _processStartBlockNumber);

				// to avoid computing sieving multiples when the segments are contiguous
				if (previousProcessEndBlockNumber != _processStartBlockNumber) {
					this->template computeSievingMultiples(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);
					previousProcessEndBlockNumber = _processEndBlockNumber;
				}

				this->template removeComposites(_processStartBlockNumber, _processEndBlockNumber, sievingMultiples);
				++_numberSegmentsSieved;

				if (sendPrimesCountToRoot) {
					this->template sendPrimesCountToCollector(processIDToSendPrimesCount);
				} else if (countNumberOfPrimesOnNode) {
					this->template countPrimesInNode();
				}

				if (sendResultsToRoot) {
					this->template sendSievedBlockDataToCollector(processIDToSendResults);
				} else {
					this->template outputResults();
				}

#				ifdef DEBUG_OUTPUT
				cout << "    --> Finish sieving segment [" << _processStartBlockNumber << ", " << (_processEndBlockNumber - 1) << "]" << endl;
#				endif
			}

			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			performanceTimer.stop();
			int _processID = this->template getProcessId();
			cout << "\n\n    > Finish sieving " << _numberSegmentsSieved << " segments of " << _dynamicSchedulingSegmentSizeInElements << " numbers each, in process with rank " << _processID << " in "
					<< performanceTimer.getElapsedTimeFormated() << endl;
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
					<< computeSievingPrimesTimer.getElapsedTimeFormated() << endl << endl;
		}

		virtual void setupSegmentDistribution(size_t maxRange) {
			size_t numberProcesses = (size_t) this->template getNumberProcesses() - 1; // process 0 acts only as a management node
			_segmentsDistribution = vector<vector<pair<size_t, size_t> > >(numberProcesses, vector<pair<size_t, size_t> >());

			size_t maxRangeSquareRootOffset = (size_t) sqrt(maxRange) + 1;
			size_t numberSegmentsPerProcess = _dynamicSchedulingNumberSegments / numberProcesses;

			size_t numbersegmentsAdded = 0;
			size_t dynamicSchedulingBlockSizeInElements = _dynamicSchedulingSegmentSizeInElements;

#			ifdef DEBUG_OUTPUT
			cout << "    > Initializing segment distribution in root node (" << _dynamicSchedulingNumberSegments << " segments in " << numberProcesses << " processes)" << endl;
#			endif

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

#					ifdef DEBUG_OUTPUT
					cout << "    --> Added segment [" << segmentBeginNumber << ", " << (segmentEndNumber - 1) << "] to list of probable segments to process " << (processPosition + 1) << endl;
#					endif

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

#				ifdef DEBUG_OUTPUT
				cout << "    --> Added segment [" << segmentBeginNumber << ", " << (segmentEndNumber - 1) << "] to list of probable segments to process " << (lastProcessNumberIndex + 1) << endl;
#				endif
			}

			_numberSegmentsCreated = numbersegmentsAdded;
			cout << "    --> Root process initialized " << _numberSegmentsCreated << " segments to sieve\n" << endl;
		}

		virtual void setupRootSegmentManagementTasks(size_t maxRangeSquareRoot, size_t numberSievingProcesses) {
			size_t numberSegmentsCreated = _numberSegmentsCreated;
			vector<vector<pair<size_t, size_t> > >& segmentsDistribution = _segmentsDistribution;

#			pragma omp parallel \
				default(shared) \
				firstprivate(numberSegmentsCreated) \
				shared(segmentsDistribution)
			{
#				pragma omp single nowait
				{
#					pragma omp task \
						default(shared) \
						firstprivate(numberSegmentsCreated) \
						shared(segmentsDistribution)
					this->template startServerOfSegmentsAvailable(numberSegmentsCreated, segmentsDistribution, numberSievingProcesses);

					if (this->template isSendPrimesCountToRoot()) {
#						pragma omp task default(none) firstprivate(numberSegmentsCreated)
						{
//							vector<pair<size_t, size_t> > sievingMultiples;
//							this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);
							this->template startCollectorOfNodesPrimesCount(numberSegmentsCreated);
						}
					}

					if (this->template isSendResultsToRoot()) {
#						pragma omp task default(none)
						this->template startCollectorOfSievingDataResultsTask();
					}
				}

#				pragma omp taskwait
			}
		}

		virtual void setupRootSegmentManagementSections(size_t maxRangeSquareRoot, size_t numberSievingProcesses) {
			size_t numberSegmentsCreated = _numberSegmentsCreated;
			vector<vector<pair<size_t, size_t> > >& segmentsDistribution = _segmentsDistribution;
			size_t numberPrimesFound;

#			pragma omp parallel sections \
				default(shared) \
				firstprivate(numberSegmentsCreated) \
				shared(segmentsDistribution, numberPrimesFound)
			{
#				pragma omp section
				{
					this->template startServerOfSegmentsAvailable(numberSegmentsCreated, segmentsDistribution, numberSievingProcesses);
				}

#				pragma omp section
				{
					if (this->template isSendPrimesCountToRoot()) {
						vector<pair<size_t, size_t> > sievingMultiples;
						this->template calculateSievingMultiples(maxRangeSquareRoot, sievingMultiples);
						numberPrimesFound = this->template startCollectorOfNodesPrimesCount(numberSegmentsCreated);
					}
				}

#				pragma omp section
				{
					if (this->template isSendResultsToRoot()) {
						this->template startCollectorOfSievingDataResults(numberSegmentsCreated);
					}
				}
			}

			this->template setPrimesCount(numberPrimesFound);
		}

		virtual bool getNewSegmentFromRoot() {
			unsigned long long segmentRange[2];

#			ifdef DEBUG_OUTPUT
			cout << "\n    > Requesting new segment from root node..." << endl;
#			endif

			MPI_Status status;
			MPI_Send(&_numberSegmentsSieved, 1, MPI_UNSIGNED_LONG_LONG, 0, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD);
			MPI_Recv(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, 0, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD, &status);
//			MPI_Sendrecv(&_numberSegmentsSieved, 1, MPI_UNSIGNED_LONG_LONG, 0, MSG_REQUEST_NEW_SEGMET, &segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, 0, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
				_processStartBlockNumber = segmentRange[0];
				_processEndBlockNumber = segmentRange[1];

				if (_processStartBlockNumber == 0 && _processEndBlockNumber == 0) {
					return false;
				}

				this->template setStartSieveNumber(_processStartBlockNumber);
				this->template setMaxRange(_processEndBlockNumber - 1); // processEndBlockNumber is max + 1 to use < operator instead of <=

#				ifdef DEBUG_OUTPUT
				int processID = this->template getProcessId();
				cout << "    --> Process " << processID << " received segment [" << segmentRange[0] << ", " << (segmentRange[1] - 1) << "]" << endl;
#				endif
				return true;
			} else {
				return false;
			}
		}

		virtual void startServerOfSegmentsAvailable(size_t numberSegmentsCreated, vector<vector<pair<size_t, size_t> > >& segmentsDistribution, size_t numberSievingProcesses) {
			cout << "    > Root server for available segments started..." << endl;
			size_t numberSegmentsRemaining = numberSegmentsCreated;
			while (numberSegmentsRemaining > 0) {
				MPI_Status status;
				unsigned long long numberSegmentsSievedByNode;

#				ifdef DEBUG_OUTPUT
				cout << "    --> Ready to send new segment..." << endl;
#				endif

				MPI_Recv(&numberSegmentsSievedByNode, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					pair<size_t, size_t> newSegmentRange;
					size_t processPositionToTakeSegment = status.MPI_SOURCE - 1; // process 0 isn't present in segment distribution because only acts as manager

					if (segmentsDistribution[processPositionToTakeSegment].empty() || processPositionToTakeSegment >= segmentsDistribution.size()) {
						processPositionToTakeSegment = this->template getProcessPositionWithMoreSegmentsRemaining();
					}

					newSegmentRange = segmentsDistribution[processPositionToTakeSegment].back();
					segmentsDistribution[processPositionToTakeSegment].pop_back();

					unsigned long long segmentRange[2];
					segmentRange[0] = newSegmentRange.first;
					segmentRange[1] = newSegmentRange.second;

#					ifdef DEBUG_OUTPUT
					cout << "    --> Assigning segment [" << segmentRange[0] << ", " << (segmentRange[1] - 1) << "] to process " << status.MPI_SOURCE << endl;
					cout << "    --> " << numberSegmentsRemaining << " remaining segments" << endl;
#					endif

					MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD);
					--numberSegmentsRemaining;
				} else {
					cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving request for a new segment from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				}
			}

			this->template endSievingNodes(numberSievingProcesses);
		}

		virtual void endSievingNodes(size_t numberSievingProcesses) {
			cout << "\n    > All segments processed. Terminating auxiliary nodes..." << endl;
			for (size_t numberNodesTerminated = 0; numberNodesTerminated < numberSievingProcesses; ++numberNodesTerminated) {
				MPI_Status status;
				unsigned long long numberSegmentsSievedByNode;

#				ifdef DEBUG_OUTPUT
				cout << "    --> Waiting for sieve node to terminate..." << endl;
#				endif

				MPI_Recv(&numberSegmentsSievedByNode, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_REQUEST_NEW_SEGMET, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					unsigned long long segmentRange[2];
					segmentRange[0] = 0;
					segmentRange[1] = 0;

					MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, MSG_ASSIGN_NEW_SEGMET, MPI_COMM_WORLD);

#					ifdef DEBUG_OUTPUT
					cout << "    --> Node " << status.MPI_SOURCE << " finished" << endl;
#					endif
				} else {
					cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving request for a new segment from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				}
			}

			PerformanceTimer& performanceTimer = this->template getPerformanceTimer();
			cout << "\n    >>>>> Finish sieving in " << _numberSievingProcesses << " sieving processes in " << performanceTimer.getElapsedTimeFormated() << " <<<<<\n" << endl;

			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				int processIDForPrimesCountCollector = 1;
				if (this->template isSendResultsToRoot()) {
					cout << "    --> Waiting for results collector to terminate..." << endl;
					double elapsedTimeMicroSec;
					MPI_Status status;
					MPI_Recv(&elapsedTimeMicroSec, 1, MPI_DOUBLE, 1, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD, &status);
					++processIDForPrimesCountCollector;
					cout << "    --> Results collector finished\n" << endl;
				}

				if (this->template isSendPrimesCountToRoot()) {
					cout << "    --> Waiting for primes count collector to terminate..." << endl;
					double elapsedTimeMicroSec;
					MPI_Status status;
					MPI_Recv(&elapsedTimeMicroSec, 1, MPI_DOUBLE, processIDForPrimesCountCollector, MSG_NODE_SIEVING_FINISHED, MPI_COMM_WORLD, &status);
					cout << "    --> Primes count collector finished" << endl;
				}
			}
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

		virtual void sendSievedBlockDataToCollector(int processIDToSendResults) {
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t blockSize = this->template getNumberBitsToStoreBlock(_processEndBlockNumber - _processStartBlockNumber);

			// send results range
			unsigned long long segmentRange[2];
			segmentRange[0] = (unsigned long long) _processStartBlockNumber;
			segmentRange[1] = (unsigned long long) _processEndBlockNumber;
			MPI_Send(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, processIDToSendResults, MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE, MPI_COMM_WORLD);

			// send block sieving data
			this->template sendSievingDataMPI(primesBitset, 0, blockSize, processIDToSendResults, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
		}

		virtual size_t receiveSievedBlockFromNode() {
			// received block range
			unsigned long long segmentRange[2];
			MPI_Status status;

#			ifdef DEBUG_OUTPUT
			cout << "    --> Ready to receive segment results from sieving node..." << endl;
#			endif

			MPI_Recv(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
				size_t blockSize = this->template getNumberBitsToStoreBlock(segmentRange[1] - segmentRange[0]);
				size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(segmentRange[0]);

				// receive sieving data
				FlagsContainer& primesBitset = this->template getPrimesBitset();
				this->template receiveSievingDataMPI(primesBitset, positionToStoreResults, blockSize, status.MPI_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
				return blockSize;
			} else {
				cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving sieving segment results data from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				return 0;
			}
		}

		virtual void startCollectorOfSievingDataResults(size_t numberSegmentsCreated) {
			size_t numberBytesReceived = 0;
			cout << "    > Collector for nodes segments results started..." << endl;
			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				numberBytesReceived += this->template receiveSievedBlockFromNode();
			}

			cout << "    --> Results collector received " << numberBytesReceived << " bytes from " << numberSegmentsCreated << " segments" << endl;
		}

		virtual void startCollectorOfSievingDataResultsTask() {
			cout << "    > Collector for nodes segments results started..." << endl;

#			pragma omp for \
				schedule(dynamic)
			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < _numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				this->template receiveSievedBlockFromNode();
			}
		}

		virtual size_t startCollectorOfNodesPrimesCount(size_t numberSegmentsCreated) {
			cout << "    > Collector for nodes segments primes count started..." << endl;

			// root node only has sieving primes
			size_t numberSievingPrimes = this->template getSievingPrimes().size();
			size_t numberPrimesSievedByTheWheel = this->template getWheelSieve().getNumberPrimesSievedByTheWheel();
			size_t numberPrimesFound = numberSievingPrimes + numberPrimesSievedByTheWheel;

			for (size_t numberSegmentsPrimesCountReceived = 0; numberSegmentsPrimesCountReceived < numberSegmentsCreated; ++numberSegmentsPrimesCountReceived) {
				numberPrimesFound += this->template receivePrimesCountFromNode(MPI_ANY_SOURCE);
			}

			cout << "\n    >>>>> Counted " << numberPrimesFound << " primes on " << _numberSievingProcesses << " sieving processes <<<<<\n" << endl;
			return numberPrimesFound;
		}

		virtual size_t receivePrimesCountFromNode(int nodeID) {
			unsigned long long primesCount;
			MPI_Status status;

#			ifdef DEBUG_OUTPUT
			cout << "    --> Ready to receive primes count from node..." << endl;
#			endif

			MPI_Recv(&primesCount, 1, MPI_UNSIGNED_LONG_LONG, nodeID, MSG_NODE_PRIMES_COUNT_FOUND, MPI_COMM_WORLD, &status);

#			ifdef DEBUG_OUTPUT
			cout << "    --> Received " << primesCount << " primes count from node " << status.MPI_SOURCE << endl;
#			endif

			if (status.MPI_ERROR == MPI_SUCCESS) {
				return (size_t) primesCount;
			} else {
				cerr << "    --> MPI_Recv detected the following error code when receiving segment primes count from " << status.MPI_SOURCE << ": " << status.MPI_ERROR << endl;
				return 0;
			}
		}

		virtual bool outputResults() {
			if (_outputOnlyLastSegment) {
				size_t maxRangeSegment = this->template getMaxRange();

				if (maxRangeSegment != _globalMaxRange) {
					return false;
				}
			}

			int _processID = this->template getProcessId();
			bool _sendResultsToRoot = this->template isSendResultsToRoot();
			bool _sendPrimesCountToRoot = this->template isSendPrimesCountToRoot();

			if (_mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				if (_processID == 0 || (_processID == 1 && !_sendResultsToRoot && _sendPrimesCountToRoot) || (_processID == 2 && _sendResultsToRoot && _sendPrimesCountToRoot)) {
					return false;
				}
			}

			if (_outputResultsFilename == "stdout") {
				cout << "\n\n=============================================  Computed primes  =============================================\n\n";
				this->template printPrimesToConsole();
				cout << "\n" << endl;
				return true;
			} else if (_outputResultsFilename != "") {
				bool sendResultsToRoot = this->template isSendResultsToRoot();

				int processID = this->template getProcessId();

				if (sendResultsToRoot) {
					cout << "\n    > Exporting results to file " << _outputResultsFilename << " in process " << processID << "..." << endl;
				}

				PerformanceTimer performanceTimer;
				performanceTimer.reset();
				performanceTimer.start();
				if (this->template savePrimesToFile(_outputResultsFilename)) {
					performanceTimer.stop();

					if (sendResultsToRoot) {
						cout << "    --> Export to file " << _outputResultsFilename << " in process " << processID << " finished in " << performanceTimer.getElapsedTimeFormated() << endl;
					}
					return true;
				} else {
					cerr << "    !!!!! Export to file " << _outputResultsFilename << "failed !!!!!" << endl;
				}
			}
			return false;
		}

		int getMpiThreadSupport() const {
			return _mpiThreadSupport;
		}

		void setMpiThreadSupport(int mpiThreadSupport) {
			_mpiThreadSupport = mpiThreadSupport;
		}

		virtual int getProcessIDWithFirstPrimesBlock() {
			return _processIDWithFirstPrimesBlock;
		}
};

