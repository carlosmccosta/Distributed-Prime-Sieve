#include "PrimesCLI.h"

void outOfMemmoryHandler() {
	cerr << "\n\n##### Unable to allocate memory!!! #####\n\n" << endl;

	int flagFinalize;
	MPI_Finalized(&flagFinalize);
	if (!flagFinalize) {
		MPI_Finalize();
	}
	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv) {
	std::set_new_handler(outOfMemmoryHandler);

	PrimesCLI primesCLI;

	int flagInit;
	MPI_Initialized(&flagInit);
	if (!flagInit) {
		MPI_Init(&argc, &argv);
	}

	if (argc == 1) {
		primesCLI.startInteractiveCLI();
	} else {
		primesCLI.showProgramHeader();
		if (primesCLI.parseCLIParameters(argc, argv)) {
			primesCLI.showCurrentConfiguration();
			if (primesCLI.computePrimes()) {
				primesCLI.countNumberOfPrimes();
				primesCLI.checkPrimesFromFile();
				primesCLI.outputResults();
			} else {
				cout << "    --> Failed to allocate memory for primes computation!\n" << endl;
			}
		}
	}

	int flagFinalize;
	MPI_Finalized(&flagFinalize);
	if (!flagFinalize) {
		MPI_Finalize();
	}

	return 0;
}

void PrimesCLI::startInteractiveCLI() {
	do {
		showProgramHeader();

		cout << "  1 - Single processor implementation (using division to check for primes)\n";
		cout << "  2 - Single processor implementation (using multiples to check for primes)\n";
		cout << "  3 - Single processor implementation (using block search with bitset with all numbers)\n";
		cout << "  4 - Single processor implementation (using block search with bitset only with numbers in the block)\n";
		cout << "  5 - Single processor implementation (using block search with bitset with all numbers optimized for time)\n";
		cout << "  6 - Single processor implementation (using block search with bitset with all numbers optimized for space and with modulo 30 wheel factorization)\n";
		cout << "  7 - Single processor implementation (using block search with bitset with all numbers optimized for space and with modulo 210 wheel factorization)\n";
		cout << "  8 - Single processor implementation (using block search with bitset with all numbers optimized for time and space and with modulo 30 wheel factorization)\n";
		cout << "  9 - Single processor implementation (using block search with bitset with all numbers optimized for time and space and with modulo 210 wheel factorization)\n";
		cout << " 10 - Single processor implementation (using block search with bitset with all numbers optimized for time and with modulo 30 wheel factorization)\n";
		cout << " 11 - Single processor implementation (using block search with bitset with all numbers optimized for time and with modulo 210 wheel factorization)\n";
		cout << " 12 - OpenMP implementation optimized for space and time and with modulo 210 wheel\n";
		cout << " 13 - OpenMP implementation optimized for time and with modulo 210 wheel\n";
		cout << " 14 - OpenMPI implementation optimized for space and time and with modulo 210 wheel\n\n";
		cout << "  0 - Exit\n\n\n" << endl;

		_algorithmToUse = ConsoleInput::getInstance()->getIntCin("  >>> Implementation to use [0, 14]: ", "    -> Insert one of the listed algorithms!\n", 0, 15);

		if (_algorithmToUse == 0) {
			break;
		}

		bool inputRangeInBits = ConsoleInput::getInstance()->getYesNoCin("\n   ## Max primes search range in bits (2^n)? (No for direct max search range specification) (Y/N): ");

		if (inputRangeInBits) {
			_primesMaxRange = ConsoleInput::getInstance()->getIntCin("    # n: ", "Range must be [4, 64]", 4, 65);
			_primesMaxRange = (size_t) pow(2, _primesMaxRange);
		} else {
			_primesMaxRange = ConsoleInput::getInstance()->getIntCin("    # Max search range: ", "Range must be >= 11", 11);
		}

		if (_algorithmToUse > 2) {
			_blockSize = ConsoleInput::getInstance()->getIntCin("    # Block size in bytes: ", "Block size must be > 4", 5);
		}

		if (_algorithmToUse == 12 || _algorithmToUse == 13) {
			_numberOfThreadsToUseInSieving = ConsoleInput::getInstance()->getIntCin("    # Number of threads to use in sieving (0 to let openMP decide): ", "Number of threads must be >= 0");
		}

		cout << "   ## Output result to file (filename, stdout or empty to avoid output): ";
		_outputResultsFilename = ConsoleInput::getInstance()->getLineCin();

		cout << "   ## Confirm results from file (empty to skip confirmation): ";
		_resultsConfirmationFile = ConsoleInput::getInstance()->getLineCin();

		computePrimes();

		_countNumberOfPrimes = ConsoleInput::getInstance()->getYesNoCin("\n   ## Count primes found? (Y/N): ");
		countNumberOfPrimes();

		checkPrimesFromFile();
		outputResults();

		cout << endl << endl;

		ConsoleInput::getInstance()->getUserInput();
		ConsoleInput::getInstance()->clearConsoleScreen();
	} while (_algorithmToUse != 0);
}

bool PrimesCLI::computePrimes() {
	delete _primesSieve;
	_primesSieve = NULL;
	delete _primesSieveMPI;
	_primesSieveMPI = NULL;
	bool validAlgorithmToUse = true;

	switch (_algorithmToUse) {
		case 1: {
			_primesSieve = new PrimesSieveSequencialDivision<vector<bool> >();
			break;
		}

		case 2: {
			_primesSieve = new PrimesSieveSequencialMultiples<vector<bool> >();
			break;
		}

		case 3: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> >(_blockSize);
			break;
		}

		case 4: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<bool> >(_blockSize);
			break;
		}

		case 5: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> >(_blockSize);
			break;
		}

		case 6: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_blockSize);
			break;
		}

		case 7: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_blockSize);
			break;
		}

		case 8: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_blockSize);
			break;
		}

		case 9: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_blockSize);
			break;
		}

		case 10: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_blockSize);
			break;
		}

		case 11: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_blockSize);
			break;
		}

		case 12: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPSpaceTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_blockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 13: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_blockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 14: {
			_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte>(_primesMaxRange, _blockSize);
			break;
		}

		default: {
			validAlgorithmToUse = false;
			break;
		}
	}

	if (validAlgorithmToUse) {
		if ((_algorithmToUse < 14 && _primesSieve == NULL) || (_algorithmToUse >= 14 && _primesSieveMPI == NULL)) {
			return false;
		}

		size_t processStartBlockNumber = (_algorithmToUse > 13 ? _primesSieveMPI->getStartSieveNumber() : _primesSieve->getStartSieveNumber());
		size_t processEndBlockNumber = (_algorithmToUse > 13 ? _primesSieveMPI->getMaxRange() : _primesMaxRange);
		if (_algorithmToUse > 13) {
			int processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char> >*) _primesSieveMPI)->getProcessId();
			cout << "\n    > Process with rank " << processRank;
		}
		cout << "\n    > Computing primes from " << processStartBlockNumber << " to " << processEndBlockNumber << "..." << endl;
		(_algorithmToUse > 13 ? _primesSieveMPI->computePrimes(_primesMaxRange) : _primesSieve->computePrimes(_primesMaxRange));

		cout << "    --> Finished ";
		if (_algorithmToUse > 13) {
			int processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char> >*) _primesSieveMPI)->getProcessId();
			cout << "process with rank " << processRank << " ";
		}
		cout << "in " << (_algorithmToUse > 13 ? _primesSieveMPI->getPerformanceTimer().getElapsedTimeFormated() : _primesSieve->getPerformanceTimer().getElapsedTimeFormated());
		if (_algorithmToUse == 12 || _algorithmToUse == 13) {
			if (_numberOfThreadsToUseInSieving != 0) {
				cout << " using " << ((PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads() << " threads";
			} else {
				cout << " using at most " << omp_get_num_procs() << " processors and " << omp_get_max_threads() << " threads";
			}
		}
		cout << "\n" << endl;
		return true;
	} else {
		return false;
	}

}

size_t PrimesCLI::countNumberOfPrimes() {
	if (_countNumberOfPrimes && (_algorithmToUse > 13 ? (((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char> >*) _primesSieveMPI)->getProcessId() == 0) : true)) {
		cout << "    > Counting number of primes found..." << endl;
		PerformanceTimer countingPrimesTimer;
		countingPrimesTimer.reset();
		countingPrimesTimer.start();
		size_t numberPrimesFound = (_algorithmToUse > 13 ? _primesSieveMPI->getNumberPrimesFound() : _primesSieve->getNumberPrimesFound());
		countingPrimesTimer.stop();
		cout << "    --> Computed " << numberPrimesFound << " primes in " << countingPrimesTimer.getElapsedTimeFormated() << endl;

		return numberPrimesFound;
	}

	return 0;
}

bool PrimesCLI::checkPrimesFromFile() {
	if (_resultsConfirmationFile != "" && (_algorithmToUse > 13 ? (((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char> >*) _primesSieveMPI)->getProcessId() == 0) : true)) {
		cout << "\n    > Validating computed primes with result file supplied...\n";
		bool validationResult = (_algorithmToUse > 13 ? _primesSieveMPI->checkPrimesFromFile(_resultsConfirmationFile) : _primesSieve->checkPrimesFromFile(_resultsConfirmationFile));

		if (validationResult) {
			cout << "    --> Computed primes are correct!\n\n";
			return true;
		} else {
			cout << "    --> Computed primes are different from the ones in supplied file!\n\n";
			return false;
		}
	}

	return false;
}

bool PrimesCLI::outputResults() {
	if (_algorithmToUse > 13 ? (((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char> >*) _primesSieveMPI)->getProcessId() == 0) : true) {
		if (_outputResultsFilename == "stdout") {
			cout << "\n\n=============================================  Computed primes  =============================================\n\n";
			(_algorithmToUse > 13 ? _primesSieveMPI->printPrimesToConsole() : _primesSieve->printPrimesToConsole());
			cout << "\n" << endl;
		} else if (_outputResultsFilename != "") {
			cout << "\n    > Exporting results to file " << _outputResultsFilename << "...";
			(_algorithmToUse > 13 ? _primesSieveMPI->savePrimesToFile(_outputResultsFilename) : _primesSieve->savePrimesToFile(_outputResultsFilename));
			cout << "\n    --> Export finished!" << endl;
		} else {
			return false;
		}
	} else {
		return false;
	}

	return true;
}

bool PrimesCLI::parseCLIParameters(int argc, char** argv) {
	if (argc <= 1) {
		return false;
	}

	for (int argNumber = 1; argNumber < argc; ++argNumber) {
		string argSelector(argv[argNumber]);

		if (argSelector == "--help") {
			showUsage(argv[0]);
			return false;
		} else if (argSelector == "--version") {
			showVersion();
			return false;
		}

		if (++argNumber >= argc) {
			showUsage(argv[0], "  >>> Missing argument for --<argSelector>");
			return false;
		}

		string argValue(argv[argNumber]);

		if (argSelector == "--algorithm") {
			stringstream sstream(argValue);
			int algorithm;
			if (!(sstream >> algorithm) || (algorithm < 1 || algorithm > 14)) {
				showUsage(argv[0], "  >>> Invalid algorithm selector! Must be a number [1, 13]");
				return false;
			} else {
				_algorithmToUse = algorithm;
			}
		} else if (argSelector == "--maxRange") {
			stringstream sstream(argValue);
			size_t range;
			if (!(sstream >> range) || (range < 11)) {
				showUsage(argv[0], "  >>> Invalid primes max range! Max range must be >= 11");
				return false;
			} else {
				_primesMaxRange = range;
			}
		} else if (argSelector == "--blockSize") {
			stringstream sstream(argValue);
			size_t blockSize;
			if (!(sstream >> blockSize) || (blockSize < 4)) {
				showUsage(argv[0], "  >>> Invalid block size! Block size must be >= 4");
				return false;
			} else {
				_blockSize = blockSize;
			}
		} else if (argSelector == "--numberThreads") {
			stringstream sstream(argValue);
			size_t numThreads;
			if (!(sstream >> numThreads) || (numThreads < 0)) {
				showUsage(argv[0], "  >>> Invalid number of threads! Number of threads must be >= 0 (0 to use default)");
				return false;
			} else {
				_numberOfThreadsToUseInSieving = numThreads;
			}
		} else if (argSelector == "--outputResult") {
			_outputResultsFilename = argValue;
		} else if (argSelector == "--checkResult") {
			_resultsConfirmationFile = argValue;
		} else if (argSelector == "--countPrimes") {
			if (argValue == "Y" || argValue == "y") {
				_countNumberOfPrimes = true;
			} else if (argValue == "N" || argValue == "n") {
				_countNumberOfPrimes = false;
			} else {
				showUsage(argv[0], "  >>> Invalid primes count flag! Flag must be Y or N");
				return false;
			}
		} else {
			return false;
		}
	}

	return true;
}

void PrimesCLI::showCurrentConfiguration() {
	cout << " >>> Current configuration:\n";
	cout << "\t --algorithm      -> " << _algorithmToUse << "\n";
	cout << "\t --maxRange       -> " << _primesMaxRange << "\n";
	cout << "\t --blockSize      -> " << _blockSize << "\n";
	cout << "\t --numberThreads  -> " << _numberOfThreadsToUseInSieving << "\n";
	cout << "\t --outputResult   -> " << _outputResultsFilename << "\n";
	cout << "\t --checkResult    -> " << _resultsConfirmationFile << "\n";
	cout << "\t --countPrimes    -> " << (_countNumberOfPrimes ? "Y" : "N") << "\n" << endl;
}

void PrimesCLI::showUsage(string programName, string message) {
	if (message != "") {
		cout << message << "\n" << endl;
	}

	cout << " >>> Usage:" << endl;
	cout << programName << " [--algorithm <number>] [--maxRange <number>] [--blockSize <number>] [--numberThreads <number>] [--outputResult <filename>] [--checkResult <filename>] [--countPrimes <Y/N>] [--help] [--version]" << endl;
	cout << "\t --algorithm      -> number in [1, 13]" << endl;
	cout << "\t --maxRange       -> number >= 11 (default 2^32)" << endl;
	cout << "\t --blockSize      -> block size in bytes >= 4 (default 16384)" << endl;
	cout << "\t --numberThreads  -> number threads to use in sieving >= 0 (default 0 -> let algorithm choose the best number of threads)" << endl;
	cout << "\t --outputResult   -> filename of file to output results (default doesn't output results)" << endl;
	cout << "\t --checkResult    -> filename of file with primes to check the algorithm result (default doesn't check algorithm result)" << endl;
	cout << "\t --countPrimes    -> Y or N to count the primes computed (default N)" << endl;
	cout << "\t --help           -> show program usage" << endl;
	cout << "\t --version        -> show program version" << endl;
	cout << "\n\tWith no arguments starts interactive command line interface\n\n" << endl;

}

void PrimesCLI::showVersion() {
	cout << "Version 1.0 developed for Parallel computing (4th year, 2nd semester, MIEIC, FEUP)" << endl;
	cout << "Author: Carlos Miguel Correia da Costa (carlos.costa@fe.up.pt / carloscosta.cmcc@gmail.com)" << endl;
	cout << "Released 08/07/2013\n\n" << endl;
}

void PrimesCLI::showProgramHeader() {
	cout << "#######################################################################################\n";
	cout << "  >>>               Primes generator - Sieve of Eratosthenes                      <<<  \n";
	cout << "#######################################################################################\n\n";
}
