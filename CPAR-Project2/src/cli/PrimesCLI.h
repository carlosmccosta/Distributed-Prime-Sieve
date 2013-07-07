#pragma once

#include "../lib/ConsoleInput.h"
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

#include <cmath>
#include <string>
#include <omp.h>

using std::string;

class PrimesCLI {
	protected:
		PrimesSieve<vector<bool> >* _primesSieve;

		int _algorithmToUse;
		size_t _primesMaxRange;
		size_t _blockSize;
		size_t _numberOfThreadsToUseInSieving;
		string _outputResultsFilename;
		string _resultsConfirmationFile;
		bool _countNumberOfPrimes;

	public:
		PrimesCLI() :
			_primesSieve(NULL),
			_algorithmToUse(13),
			_primesMaxRange(pow(2, 32)),
			_blockSize(16384),
			_numberOfThreadsToUseInSieving(0),
			_outputResultsFilename(""),
			_resultsConfirmationFile(""),
			_countNumberOfPrimes(false) {}

		virtual ~PrimesCLI() {
			delete _primesSieve;
		}

		void startInteractiveCLI();
		bool computePrimes();
		size_t countNumberOfPrimes();
		bool checkPrimesFromFile();
		bool outputResults();

		bool parseCLIParameters(int argc, char** argv);
		void showUsage(string programName, string message = "");
		void showVersion();
		void showProgramHeader();


		int getAlgorithmToUse() const {
			return _algorithmToUse;
		}

		void setAlgorithmToUse(int algorithmToUse) {
			_algorithmToUse = algorithmToUse;
		}

		size_t getBlockSize() const {
			return _blockSize;
		}

		void setBlockSize(size_t blockSize) {
			_blockSize = blockSize;
		}

		bool isCountNumberOfPrimes() const {
			return _countNumberOfPrimes;
		}

		void setCountNumberOfPrimes(bool countNumberOfPrimes) {
			_countNumberOfPrimes = countNumberOfPrimes;
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
};
