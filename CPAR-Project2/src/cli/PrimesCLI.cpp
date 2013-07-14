#include "PrimesCLI.h"

void outOfMemmoryHandler() {
	int flag;
	MPI_Initialized(&flag);
	if (flag) {
		int processID;
		MPI_Comm_rank(MPI_COMM_WORLD, &processID);
		cerr << "\n\n##### Unable to allocate memory in process with rank " << processID << " !!! #####\n\n" << endl;

		int flagFinalize;
		MPI_Finalized(&flagFinalize);
		if (!flagFinalize) {
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
			MPI_Finalize();
		}
	} else {
		cerr << "\n\n##### Unable to allocate memory!!! #####\n\n" << endl;
	}

	exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
	std::set_new_handler(outOfMemmoryHandler);

	PrimesCLI primesCLI;
	primesCLI.startTimer();

	int flagInit;
	MPI_Initialized(&flagInit);
	if (!flagInit) {
		MPI_Init(&argc, &argv);
	}

	primesCLI.setProgramName(argv[0]);

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

	primesCLI.showTotalComputationTimte();

	return 0;
}

void PrimesCLI::startInteractiveCLI() {
	do {
		showProgramHeader();

		cout << "  1 - Single processor implementation using modulo division to cross of composites\n";
		cout << "  2 - Single processor implementation using primes multiples to cross of composites\n";
		cout << "  3 - Single processor implementation using block search with bitset with all even numbers\n";
		cout << "  4 - Single processor implementation using block search with bitset containing only the even numbers in the block\n";
		cout << "  5 - Single processor implementation using block search with bitset with all even numbers optimized for time\n";
		cout << "  6 - Single processor implementation using block search with bitset with only possible primes numbers optimized for space and with modulo 30 wheel factorization\n";
		cout << "  7 - Single processor implementation using block search with bitset with only possible primes numbers optimized for space and with modulo 210 wheel factorization\n";
		cout << "  8 - Single processor implementation using block search with bitset with all even numbers optimized for time and space and with modulo 30 wheel factorization\n";
		cout << "  9 - Single processor implementation using block search with bitset with all even numbers optimized for time and space and with modulo 210 wheel factorization\n";
		cout << " 10 - Single processor implementation using block search with bitset with all numbers optimized for time and with modulo 30 wheel factorization\n";
		cout << " 11 - Fastest single processor implementation using block search with bitset with all numbers optimized for time and with modulo 210 wheel factorization\n";
		cout << " 12 - OpenMP implementation using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 13 - Fastest OpenMP implementation using block search with bitset with all numbers optimized for time and with modulo 210 wheel\n";
		cout << " 14 - OpenMPI implementation using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 15 - Hybrid implementation with OpenMPI and OpenMP using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 16 - Fastest hybrid implementation with OpenMPI and OpenMP using block search with bitset with all even numbers optimized for space and time with modulo 210 wheel and with dynamic scheduling\n\n";
		cout << " 17 - Command line help\n";
		cout << " 18 - About\n";
		cout << "  0 - Exit\n\n\n" << endl;

		_algorithmToUse = ConsoleInput::getInstance()->getIntCin("  >>> Option [0, 18]: ", "    -> Insert one of the listed algorithms!\n", 0, 19);

		if (_algorithmToUse == 0) {
			break;
		} else if (_algorithmToUse == 17) {
			ConsoleInput::getInstance()->clearConsoleScreen();
			showProgramHeader();
			showUsage();

		} else if (_algorithmToUse == 18) {
			ConsoleInput::getInstance()->clearConsoleScreen();
			showProgramHeader();
			showVersion();
		} else {
			bool inputRangeInBits = ConsoleInput::getInstance()->getYesNoCin("\n   ## Max primes search range in bits (2^n)? (No for direct max search range specification) (Y/N): ");

			if (inputRangeInBits) {
				_primesMaxRange = ConsoleInput::getInstance()->getIntCin("    # n: ", "Range must be [4, 64]", 4, 65);
				_primesMaxRange = (size_t) pow(2, _primesMaxRange);
			} else {
				_primesMaxRange = ConsoleInput::getInstance()->getIntCin("    # Max search range: ", "Range must be >= 11", 11);
			}

			if (_algorithmToUse > 2) {
				_cacheBlockSize = ConsoleInput::getInstance()->getIntCin("    # Cache block size in bytes: ", "Block size must be >= 512", 512);
			}

			if (_algorithmToUse == 16) {
				_cacheBlockSize = ConsoleInput::getInstance()->getIntCin("    # Block size in elements to split primes domain for mpi dynamic scheduling: ", "Block size must be >= 100", 100);
			}

			if (_algorithmToUse == 12 || _algorithmToUse == 13 || _algorithmToUse == 15 || _algorithmToUse == 16) {
				_numberOfThreadsToUseInSieving = ConsoleInput::getInstance()->getIntCin("    # Number of threads to use in sieving (0 to let openMP decide): ", "Number of threads must be >= 0");
			}

//			if (_algorithmToUse > 13) {
//				_sendPrimesCountToRoot = ConsoleInput::getInstance()->getYesNoCin("\n   # Perform primes computation on all processes and send it back to root node? (Y/N): ");
//				_sendResultsToRoot = ConsoleInput::getInstance()->getYesNoCin("\n   # Send primes results to root node? (Y/N): ");
//			}

			cout << "   ## Output result to file (filename, stdout or empty to avoid output): ";
			_outputResultsFilename = ConsoleInput::getInstance()->getLineCin();

			cout << "   ## Confirm results from file (empty to skip confirmation): ";
			_resultsConfirmationFile = ConsoleInput::getInstance()->getLineCin();

			if (!_sendPrimesCountToRoot) {
				_countNumberOfPrimesOnNode = ConsoleInput::getInstance()->getYesNoCin("\n   ## Count primes found? (Y/N): ");
			}

			computePrimes();
			countNumberOfPrimes();

			checkPrimesFromFile();
			outputResults();
		}

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
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> >(_cacheBlockSize);
			break;
		}

		case 4: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<bool> >(_cacheBlockSize);
			break;
		}

		case 5: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> >(_cacheBlockSize);
			break;
		}

		case 6: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 7: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 8: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 9: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 10: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 11: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 12: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPSpaceTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_cacheBlockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 13: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>(_cacheBlockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 14: {
			_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize, _sendResultsToRoot,
					_countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			break;
		}

		case 15: {
			_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize, _numberOfThreadsToUseInSieving,
					_sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			break;
		}

		case 16: {
			_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling<vector<unsigned char>, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize,
					_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
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

		cout << "\n    > Computing primes";
		if (_algorithmToUse > 13) {
			int processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
			cout << " in process with rank " << processRank;
		}
		cout << " from " << processStartBlockNumber << " to " << processEndBlockNumber << "..." << endl;
		(_algorithmToUse > 13 ? _primesSieveMPI->computePrimes(_primesMaxRange) : _primesSieve->computePrimes(_primesMaxRange));

		if (_algorithmToUse < 14) {
			cout << "    --> Finished sieving in " << _primesSieve->getPerformanceTimer().getElapsedTimeFormated();
			if (_algorithmToUse == 12 || _algorithmToUse == 13) {
				cout << " using at most " << omp_get_num_procs() << " processors and ";
				if (_numberOfThreadsToUseInSieving != 0) {
					cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads();
				} else {
					cout << omp_get_max_threads();
				}
				cout << " threads";
			}
			cout << endl;
		}
		return true;
	} else {
		return false;
	}

}

size_t PrimesCLI::countNumberOfPrimes() {
	if (_countNumberOfPrimesOnNode) {
		if (_algorithmToUse < 14) {
			size_t startPossiblePrime = _primesSieve->getStartSieveNumber();
			size_t maxRange = _primesSieve->getMaxRange();

			cout << "    --> Counting primes in [" << startPossiblePrime << ", " << maxRange << "]" << endl;

			PerformanceTimer countingPrimesTimer;
			countingPrimesTimer.reset();
			countingPrimesTimer.start();
			size_t numberPrimesFound = _primesSieve->getNumberPrimesFound();
			countingPrimesTimer.stop();
			cout << "    --> Counted " << numberPrimesFound << " primes in [" << startPossiblePrime << ", " << maxRange << "]" << " in " << countingPrimesTimer.getElapsedTimeFormated() << endl;

			return numberPrimesFound;
		} else {
			return _primesSieveMPI->getPrimesCount();
		}
	}

	return 0;
}

bool PrimesCLI::checkPrimesFromFile() {
	if (_resultsConfirmationFile != "" && (_algorithmToUse > 13 ? (((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId() == 0) : true)) {
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
	if (_algorithmToUse == 16) {
		return ((PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheelAndDynamicScheduling<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->outputResults();
	}

	if (_algorithmToUse > 13 ? (((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId() == 0 || !_sendResultsToRoot) : true) {
		if (_outputResultsFilename == "stdout") {
			cout << "\n\n=============================================  Computed primes  =============================================\n\n";
			(_algorithmToUse > 13 ? _primesSieveMPI->printPrimesToConsole() : _primesSieve->printPrimesToConsole());
			cout << "\n" << endl;
		} else if (_outputResultsFilename != "") {
			int numberProcesses = (_algorithmToUse > 13 ? ((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getNumberProcesses() : 1);

			if (_sendResultsToRoot || numberProcesses == 1) {
				cout << "\n    > Exporting results to file " << _outputResultsFilename << "..." << endl;
			}

			PerformanceTimer performanceTimer;
			performanceTimer.reset();
			performanceTimer.start();
			if ((_algorithmToUse > 13 ? _primesSieveMPI->savePrimesToFile(_outputResultsFilename) : _primesSieve->savePrimesToFile(_outputResultsFilename))) {
				performanceTimer.stop();
				cout << "    --> Export to file " << _outputResultsFilename << " finished in " << performanceTimer.getElapsedTimeFormated() << endl;
			} else {
				cerr << "    --> Export to file failed!" << endl;
			}
		} else {
			return false;
		}
	} else {
		return false;
	}

	return true;
}

void PrimesCLI::showTotalComputationTimte() {
	cout << "\n\n  > Finished program in ";
	if (_algorithmToUse > 13) {
		int processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
		cout << "process with rank " << processRank << " in ";
	}
	stopTimer();
	cout << getElapsedTimeFormated();
	if (_algorithmToUse == 12 || _algorithmToUse == 13 || _algorithmToUse == 15 || _algorithmToUse == 16) {
		cout << " using at most " << omp_get_num_procs() << " processors and using ";
		if (_numberOfThreadsToUseInSieving != 0) {
			if (_algorithmToUse == 15 || _algorithmToUse == 16) {
				cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte>*) _primesSieveMPI)->getNumberOfThreads();
			} else {
				cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<vector<bool>, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads();
			}
		} else {
			cout << omp_get_max_threads();
		}
		cout << " threads";
	}
	cout << endl;
}

bool PrimesCLI::parseCLIParameters(int argc, char** argv) {
	if (argc <= 1) {
		return false;
	}

	for (int argNumber = 1; argNumber < argc; ++argNumber) {
		string argSelector(argv[argNumber]);

		if (argSelector == "--help") {
			showUsage();
			return false;
		} else if (argSelector == "--version") {
			showVersion();
			return false;
		}

		if (++argNumber >= argc) {
			showUsage("  >>> Missing argument for --<argSelector>");
			return false;
		}

		string argValue(argv[argNumber]);

		if (argSelector == "--algorithm") {
			stringstream sstream(argValue);
			int algorithm;
			if (!(sstream >> algorithm) || (algorithm < 1 || algorithm > 16)) {
				showUsage("  >>> Invalid algorithm selector! Must be a number [1, 16]");
				return false;
			} else {
				_algorithmToUse = algorithm;
			}
		} else if (argSelector == "--maxRange") {
			stringstream sstream(argValue);
			size_t range;
			if (!(sstream >> range) || (range < 11)) {
				showUsage("  >>> Invalid primes max range! Max range must be >= 11");
				return false;
			} else {
				_primesMaxRange = range;
			}
		} else if (argSelector == "--cacheBlockSize") {
			stringstream sstream(argValue);
			size_t blockSize;
			if (!(sstream >> blockSize) || (blockSize < 512)) {
				showUsage("  >>> Invalid cache block size! Block size must be >= 512");
				return false;
			} else {
				_cacheBlockSize = blockSize;
			}
		} else if (argSelector == "--dynamicSchedulingBlockSize") {
			stringstream sstream(argValue);
			size_t scheduleBlockSize;
			if (!(sstream >> scheduleBlockSize) || (scheduleBlockSize < 100)) {
				showUsage("  >>> Invalid dynamic scheduling block size! Block size must be >= 100");
				return false;
			} else {
				_dynamicSchedulingBlockSizeInElements = scheduleBlockSize;
			}
		} else if (argSelector == "--numberThreads") {
			stringstream sstream(argValue);
			size_t numThreads;
			if (!(sstream >> numThreads) || (numThreads < 0)) {
				showUsage("  >>> Invalid number of threads! Number of threads must be >= 0 (0 to use default)");
				return false;
			} else {
				_numberOfThreadsToUseInSieving = numThreads;
			}
		} else if (argSelector == "--outputResult") {
			_outputResultsFilename = argValue;
		} else if (argSelector == "--checkResult") {
			_resultsConfirmationFile = argValue;
		} else if (argSelector == "--countPrimesInNode") {
			if (argValue == "Y" || argValue == "y") {
				_countNumberOfPrimesOnNode = true;
			} else if (argValue == "N" || argValue == "n") {
				_countNumberOfPrimesOnNode = false;
			} else {
				showUsage("  >>> Invalid --countPrimesInNode flag! Flag must be Y or N");
				return false;
			}
		} else if (argSelector == "--sendResultsToRoot") {
			if (argValue == "Y" || argValue == "y") {
				_sendResultsToRoot = true;
			} else if (argValue == "N" || argValue == "n") {
				_sendResultsToRoot = false;
			} else {
				showUsage("  >>> Invalid --sendResultsToRoot flag! Flag must be Y or N");
				return false;
			}
		} else if (argSelector == "--sendPrimesCountToRoot") {
			if (argValue == "Y" || argValue == "y") {
				_sendPrimesCountToRoot = true;
			} else if (argValue == "N" || argValue == "n") {
				_sendPrimesCountToRoot = false;
			} else {
				showUsage("  >>> Invalid --sendPrimesCountToRoot flag! Flag must be Y or N");
				return false;
			}
		} else {
			return false;
		}
	}

	return true;
}

void PrimesCLI::showCurrentConfiguration() {
	cout << "  >>> Current configuration:\n";
	cout << "\t --algorithm                   -> " << _algorithmToUse << "\n";
	cout << "\t --maxRange                    -> " << _primesMaxRange << "\n";
	cout << "\t --cacheBlockSize              -> " << _cacheBlockSize << "\n";
	cout << "\t --dynamicSchedulingBlockSize  -> " << _dynamicSchedulingBlockSizeInElements << "\n";
	cout << "\t --numberThreads               -> " << _numberOfThreadsToUseInSieving << "\n";
	cout << "\t --outputResult                -> " << _outputResultsFilename << "\n";
	cout << "\t --checkResult                 -> " << _resultsConfirmationFile << "\n";
	cout << "\t --countPrimesInNode           -> " << (_countNumberOfPrimesOnNode ? "Y" : "N") << "\n";
	cout << "\t --sendPrimesCountToRoot       -> " << (_sendPrimesCountToRoot ? "Y" : "N") << "\n";
	cout << "\t --sendResultsToRoot           -> " << (_sendResultsToRoot ? "Y" : "N") << "\n";
	cout << endl;
}

void PrimesCLI::showUsage(string message) {
	if (message != "") {
		cout << message << "\n" << endl;
	}

	cout << " >>> Usage:" << endl;
	cout << _programName;
	cout << "\n  [--algorithm <number>]";
	cout << "\n  [--maxRange <number>]";
	cout << "\n  [--cacheBlockSize <number>]";
	cout << "\n  [--dynamicSchedulingBlockSize <number>]";
	cout << "\n  [--numberThreads <number>]";
	cout << "\n  [--outputResult <filename>]";
	cout << "\n  [--checkResult <filename>]";
	cout << "\n  [--countPrimes <Y/N>]";
	cout << "\n  [--help]";
	cout << "\n  [--version]" << endl;
	cout << "\t --algorithm                   -> number in [1, 16]" << endl;
	cout << "\t --maxRange                    -> number >= 11 (default 2^32)" << endl;
	cout << "\t --cacheBlockSize              -> block size in bytes >= 512 to optimize cache hit rate (default 16384)" << endl;
	cout << "\t --dynamicSchedulingBlockSize  -> block size in elements >= 100 to split the primes domain in blocks to perform dynamic allocation in mpi (default 1048576)" << endl;
	cout << "\t --numberThreads               -> number threads to use in sieving >= 0 (default 0 -> let algorithm choose the best number of threads)" << endl;
	cout << "\t --outputResult                -> filename of file to output results (default doesn't output results)" << endl;
	cout << "\t --checkResult                 -> filename of file with primes to check the algorithm result in root node (default doesn't check algorithm result)" << endl;
	cout << "\t --countPrimesInNode           -> Y/N to count the primes computed each node (default N)" << endl;
	cout << "\t --sendPrimesCountToRoot       -> Y/N to to count the number of primes found in each mpi process and send the result to the root node (default Y and overrides the --countPrimesInNode flag if set to N)"
			<< endl;
	cout << "\t --sendResultsToRoot           -> Y/N to send the computation results to the root node (default Y)" << endl;
	cout << "\t --help                        -> show program usage" << endl;
	cout << "\t --version                     -> show program version" << endl;
	cout << "\n\tWith no arguments starts interactive command line interface\n\n" << endl;

}

void PrimesCLI::showVersion() {
	cout << "Version 1.0 developed for Parallel Computing (4th year, 2nd semester, MIEIC, FEUP)" << endl;
	cout << "Author: Carlos Miguel Correia da Costa (carlos.costa@fe.up.pt / carloscosta.cmcc@gmail.com)\n\n" << endl;
}

void PrimesCLI::showProgramHeader() {
	cout << "#######################################################################################\n";
	cout << "  >>>                           Distributed primes sieve                          <<<  \n";
	cout << "#######################################################################################\n\n";
}
