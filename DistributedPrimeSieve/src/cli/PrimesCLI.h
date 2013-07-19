#pragma once

#include "../lib/ConsoleInput.h"
#include "../WheelFactorization.h"
#include "../sequencial/PrimesSieveSequencialDivision.h"
#include "../sequencial/PrimesSieveSequencialMultiples.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceAndCache.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedTimeAndCache.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel.h"
#include "../openmp/PrimesSieveParallelMultiplesOptimizedOpenMPSpaceTimeAndCacheWithWheel.h"
#include "../openmp/PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPI.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPSpaceTimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel.h"
#include "../openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel.h"

#include <cmath>
#include <string>
#include <omp.h>
#include <mpi.h>

using std::string;

class PrimesCLI {
	protected:
		PrimesSieve<PrimesFlagsContainer>* _primesSieve;
		PrimesSieve<PrimesFlagsContainerMPI>* _primesSieveMPI;

		int _algorithmToUse;
		size_t _primesMaxRange;
		size_t _cacheBlockSize;
		size_t _dynamicSchedulingSegmentSizeInElements;
		size_t _dynamicSchedulingNumberSegments;
		size_t _numberOfThreadsToUseInSieving;
		string _outputResultsFilename;
		string _resultsConfirmationFile;
		bool _countNumberOfPrimesOnNode;
		bool _sendPrimesCountToRoot;
		bool _sendResultsToRoot;

		int _mpiThreadSupport;

		string _programName;
		PerformanceTimer performanceTimer;

	public:
		PrimesCLI() :
				_primesSieve(NULL), _primesSieveMPI(NULL), _algorithmToUse(13), _primesMaxRange(7920), _cacheBlockSize(16384), _dynamicSchedulingSegmentSizeInElements(1048576), _dynamicSchedulingNumberSegments(
						0), _numberOfThreadsToUseInSieving(0), _outputResultsFilename(""), _resultsConfirmationFile(""), _countNumberOfPrimesOnNode(false), _sendPrimesCountToRoot(false), _sendResultsToRoot(
						false), _mpiThreadSupport(MPI_THREAD_SINGLE), _programName("PrimeSieve") {
		}

		virtual ~PrimesCLI() {
			delete _primesSieve;
			_primesSieve = NULL;
			delete _primesSieveMPI;
			_primesSieveMPI = NULL;
		}

		void startInteractiveCLI();
		bool computePrimes();
		size_t countNumberOfPrimes();
		bool checkPrimesFromFile();
		bool outputResults();
		void showTotalComputationTimte();

		bool parseCLIParameters(int argc, char** argv);
		void showCurrentConfiguration();
		void showUsage(string message = "");
		void showVersion();
		void showProgramHeader();

		void startTimer() {
			performanceTimer.reset();
			performanceTimer.start();
		}

		void stopTimer() {
			performanceTimer.stop();
		}

		string getElapsedTimeFormated() {
			return performanceTimer.getElapsedTimeFormated();
		}

		int getAlgorithmToUse() const {
			return _algorithmToUse;
		}

		void setAlgorithmToUse(int algorithmToUse) {
			_algorithmToUse = algorithmToUse;
		}

		size_t getBlockSize() const {
			return _cacheBlockSize;
		}

		void setBlockSize(size_t blockSize) {
			_cacheBlockSize = blockSize;
		}

		bool isCountNumberOfPrimes() const {
			return _countNumberOfPrimesOnNode;
		}

		void setCountNumberOfPrimes(bool countNumberOfPrimes) {
			_countNumberOfPrimesOnNode = countNumberOfPrimes;
		}

		size_t getNumberOfThreadsToUseInSieving() const {
			return _numberOfThreadsToUseInSieving;
		}

		void setNumberOfThreadsToUseInSieving(size_t numberOfThreadsToUseInSieving) {
			_numberOfThreadsToUseInSieving = numberOfThreadsToUseInSieving;
		}

		const string& getOutputResultsFilename() const {
			return _outputResultsFilename;
		}

		void setOutputResultsFilename(const string& outputResultsFilename) {
			_outputResultsFilename = outputResultsFilename;
		}

		size_t getPrimesMaxRange() const {
			return _primesMaxRange;
		}

		void setPrimesMaxRange(size_t primesMaxRange) {
			_primesMaxRange = primesMaxRange;
		}

		const string& getResultsConfirmationFile() const {
			return _resultsConfirmationFile;
		}

		void setResultsConfirmationFile(const string& resultsConfirmationFile) {
			_resultsConfirmationFile = resultsConfirmationFile;
		}

		const string& getProgramName() const {
			return _programName;
		}

		void setProgramName(const string& programName) {
			_programName = programName;
		}

		int getMpiThreadSupport() const {
			return _mpiThreadSupport;
		}

		void setMpiThreadSupport(int mpiThreadSupport) {
			_mpiThreadSupport = mpiThreadSupport;
		}
};
