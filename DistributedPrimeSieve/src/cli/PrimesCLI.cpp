#include "PrimesCLI.h"

void outOfMemmoryHandler() {
	int flag;
	MPI_Initialized(&flag);
	if (flag) {
		int processID;
		MPI_Comm_rank(MPI_COMM_WORLD, &processID);
		cerr << "\n\n!!!!! Unable to allocate memory in process with rank " << processID << " !!!!!\n\n" << endl;

		int flagFinalize;
		MPI_Finalized(&flagFinalize);
		if (!flagFinalize) {
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
			MPI_Finalize();
		}
	} else {
		cerr << "\n\n!!!!! Unable to allocate memory !!!!!\n\n" << endl;
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
		int threadLevel;
//		MPI_Init(&argc, &argv);
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadLevel);

		primesCLI.setMpiThreadSupport(threadLevel);

//		if (threadLevel < MPI_THREAD_MULTIPLE) {
//			cerr << "\n\n!!!!! MPI implementation used does not provide MPI_THREAD_MULTIPLE support (provided: " << threadLevel << ") !!!!!\n\n" << endl;
//			exit(EXIT_FAILURE);
//		}
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
				primesCLI.showTotalComputationTimte();
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

		cout << "  1 - Single processor implementation using modulo division to cross of composites\n";
		cout << "  2 - Single processor implementation using primes multiples to cross of composites\n";
		cout << "  3 - Single processor implementation using block search with bitset with all even numbers\n";
		cout << "  4 - Single processor implementation using block search with bitset containing only the even numbers in the block\n";
		cout << "  5 - Single processor implementation using block search with bitset with all even numbers optimized for time\n";
		cout << "  6 - Single processor implementation using block search with bitset with only possible primes numbers optimized for space and with modulo 30 wheel factorization\n";
		cout << "  7 - Single processor implementation using block search with bitset with only possible primes numbers optimized for space and with modulo 210 wheel factorization\n";
		cout << "  8 - Single processor implementation using block search with bitset with all even numbers optimized for time and space and with modulo 30 wheel factorization\n";
		cout << "  9 - Fastest single processor implementation using block search with bitset with all even numbers optimized for time and space and with modulo 210 wheel factorization\n";
		cout << " 10 - Single processor implementation using block search with bitset with all numbers optimized for time and with modulo 30 wheel factorization\n";
		cout << " 11 - Single processor implementation using block search with bitset with all numbers optimized for time and with modulo 210 wheel factorization\n";
		cout << " 12 - Fastest OpenMP implementation using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 13 - OpenMP implementation using block search with bitset with all numbers optimized for time and with modulo 210 wheel\n";
		cout << " 14 - OpenMPI implementation using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 15 - OpenMPI implementation using block search with bitset with all numbers optimized for time and with modulo 210 wheel\n";
		cout << " 16 - Hybrid implementation with OpenMPI and OpenMP using block search with bitset with all even numbers optimized for space and time and with modulo 210 wheel\n";
		cout << " 17 - Hybrid implementation with OpenMPI and OpenMP using block search with bitset with all numbers optimized for time and with modulo 210 wheel\n";
		cout
				<< " 18 - Fastest hybrid implementation with OpenMPI and OpenMP using block search with bitset with all even numbers optimized for space and time with modulo 210 wheel and with dynamic scheduling\n";
		cout << " 19 - Hybrid implementation with OpenMPI and OpenMP using block search with bitset with all numbers optimized for time with modulo 210 wheel and with dynamic scheduling\n\n";
		cout << " 20 - Command line help\n";
		cout << " 21 - About\n";
		cout << "  0 - Exit\n\n\n" << endl;

		_algorithmToUse = ConsoleInput::getInstance()->getIntCin("  >>> Option [0, 21]: ", "    -> Insert one of the listed algorithms!\n", 0, 22);

		if (_algorithmToUse == 0) {
			break;
		} else if (_algorithmToUse == 20) {
			ConsoleInput::getInstance()->clearConsoleScreen();
			showProgramHeader();
			showUsage();

		} else if (_algorithmToUse == 21) {
			ConsoleInput::getInstance()->clearConsoleScreen();
			showProgramHeader();
			showVersion();
		} else {
			bool inputRangeInBits = ConsoleInput::getInstance()->getYesNoCin("\n   ## Max primes search range in bits (2^n)? (No for direct max search range specification) (Y/N): ");

			if (inputRangeInBits) {
				_primesMaxRange = ConsoleInput::getInstance()->getIntCin("    # n: ", "Range must be [4, 64]", 4, 65);
				_primesMaxRange = (size_t) pow(2, _primesMaxRange);
			} else {
				_primesMaxRange = ConsoleInput::getInstance()->getNumberCin("    # Max search range: ", "Range must be >= 11", (size_t) 11);
			}

			if (_algorithmToUse >= 3) {
				_cacheBlockSize = ConsoleInput::getInstance()->getIntCin("    # Cache block size in bytes: ", "Block size must be >= 512", 512);
			}

//			if (_algorithmToUse == 16) {
//				_cacheBlockSize = ConsoleInput::getInstance()->getIntCin("    # Block size in elements to split primes domain for mpi dynamic scheduling: ", "Block size must be >= 100", 100);
//			}

			if (_algorithmToUse >= 12 || _algorithmToUse <= 17) {
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
			_primesSieve = new PrimesSieveSequencialDivision<PrimesFlagsContainer>();
			break;
		}

		case 2: {
			_primesSieve = new PrimesSieveSequencialMultiples<PrimesFlagsContainer>();
			break;
		}

		case 3: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<PrimesFlagsContainer>(_cacheBlockSize);
			break;
		}

		case 4: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCache<PrimesFlagsContainer>(_cacheBlockSize);
			break;
		}

		case 5: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<PrimesFlagsContainer>(_cacheBlockSize);
			break;
		}

		case 6: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<PrimesFlagsContainer, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 7: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 8: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 9: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 10: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo30Wheel>(_cacheBlockSize);
			break;
		}

		case 11: {
			_primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_cacheBlockSize);
			break;
		}

		case 12: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_cacheBlockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 13: {
			_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_cacheBlockSize, _numberOfThreadsToUseInSieving);
			break;
		}

		case 14: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize,
						_sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize, _sendResultsToRoot,
						_countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			}
			break;
		}

		case 15: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize,
						_sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize, _sendResultsToRoot,
						_countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			}
			break;
		}

		case 16: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			}
			break;
		}

		case 17: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot);
			}
			break;
		}

		case 18: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange,
						_cacheBlockSize, _numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot, _dynamicSchedulingSegmentSizeInElements,
						_dynamicSchedulingNumberSegments, _outputResultsFilename);
//				((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->setMpiThreadSupport(_mpiThreadSupport);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot, _dynamicSchedulingSegmentSizeInElements,
						_dynamicSchedulingNumberSegments, _outputResultsFilename, _outputOnlyLastSegment);
//				((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieve)->setMpiThreadSupport(_mpiThreadSupport);
			}
			break;
		}

		case 19: {
			if (_sendResultsToRoot) {
				_primesSieveMPI = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>(_primesMaxRange,
						_cacheBlockSize, _numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot, _dynamicSchedulingSegmentSizeInElements,
						_dynamicSchedulingNumberSegments, _outputResultsFilename);

//				((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->setMpiThreadSupport(_mpiThreadSupport);
			} else {
				_primesSieve = new PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>(_primesMaxRange, _cacheBlockSize,
						_numberOfThreadsToUseInSieving, _sendResultsToRoot, _countNumberOfPrimesOnNode, _sendPrimesCountToRoot, _dynamicSchedulingSegmentSizeInElements,
						_dynamicSchedulingNumberSegments, _outputResultsFilename, _outputOnlyLastSegment);
//				((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieve)->setMpiThreadSupport(_mpiThreadSupport);
			}
			break;
		}

		default: {
			validAlgorithmToUse = false;
			break;
		}
	}

	if (validAlgorithmToUse) {
		if ((_algorithmToUse <= 13 && _primesSieve == NULL) || (_algorithmToUse >= 14 && _sendResultsToRoot && _primesSieveMPI == NULL)) {
			return false;
		}

		size_t processStartBlockNumber = ((_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->getStartSieveNumber() : _primesSieve->getStartSieveNumber());
		size_t processEndBlockNumber = ((_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->getMaxRange() : _primesMaxRange);

		if (_algorithmToUse <= 17) {
			cout << "\n    > Computing primes";
		}

		if (_algorithmToUse >= 14 && _algorithmToUse <= 17) {
			int processRank;
			if (_sendResultsToRoot) {
				processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
			} else {
				processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getProcessId();
			}
			cout << " in process with rank " << processRank;
		}

		if (_algorithmToUse <= 17) {
			cout << " from " << processStartBlockNumber << " to " << processEndBlockNumber << "..." << endl;
		}

		((_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->computePrimes(_primesMaxRange) : _primesSieve->computePrimes(_primesMaxRange));

		if (_algorithmToUse <= 13) {
			cout << "    --> Finished sieving in " << _primesSieve->getPerformanceTimer().getElapsedTimeFormated();
			if (_algorithmToUse == 12 || _algorithmToUse == 13) {
				cout << " using at most " << omp_get_num_procs() << " processors and ";
				if (_numberOfThreadsToUseInSieving != 0) {
					cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads();
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
		if (_algorithmToUse <= 13) {
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
			if (_sendResultsToRoot) {
				return _primesSieveMPI->getPrimesCount();
			} else {
				return _primesSieve->getPrimesCount();
			}
		}
	}

	return 0;
}

bool PrimesCLI::checkPrimesFromFile() {
	if (_resultsConfirmationFile != "") {
		int processID = 0;
		if (_algorithmToUse >= 14) {
			if (_sendResultsToRoot) {
				processID = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
			} else {
				processID = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieve)->getProcessId();
			}

			if ((_algorithmToUse >= 18) && _mpiThreadSupport < MPI_THREAD_MULTIPLE) {
				if (processID != 1) {
					return false;
				}
			} else {
				if (processID != 0) {
					return false;
				}
			}
		}

		cout << "\n    > Validating computed primes with file " << _resultsConfirmationFile;
		if (_algorithmToUse >= 14) {
			cout << " in process " << processID;
		}
		cout << "... " << endl;

		bool validationResult = (
				(_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->checkPrimesFromFile(_resultsConfirmationFile) : _primesSieve->checkPrimesFromFile(_resultsConfirmationFile));

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
	if (_algorithmToUse >= 18) {
		return false;
//
//		if (_sendResultsToRoot) {
//			return ((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->outputResults();
//		} else {
//			return ((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->outputResults();
//		}
	}

	if (_algorithmToUse >= 14) {
		if (_sendResultsToRoot) {
			if (((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId() != 0 && _sendResultsToRoot) {
				return false;
			}
		} else {
			if (((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getProcessId() != 0 && _sendResultsToRoot) {
				return false;
			}
		}
	}

	if (_outputResultsFilename == "stdout") {
		cout << "\n\n=============================================  Computed primes  =============================================\n\n";
		((_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->printPrimesToConsole() : _primesSieve->printPrimesToConsole());
		cout << "\n" << endl;
	} else if (_outputResultsFilename != "") {
		int numberProcesses = 1;
		int processID = 0;
		if (_algorithmToUse >= 14) {
			if (_sendResultsToRoot) {
				numberProcesses = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getNumberProcesses();
				processID = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
			} else {
				numberProcesses = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getNumberProcesses();
				processID = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieve)->getProcessId();
			}
		} else {
			numberProcesses = 1;
		}

		if (_sendResultsToRoot || numberProcesses == 1) {
			cout << "\n    > Exporting results to file " << _outputResultsFilename;
		}

		if (_algorithmToUse >= 14) {
			cout << " in process " << processID;
		}

		cout << "..." << endl;

		PerformanceTimer performanceTimer;
		performanceTimer.reset();
		performanceTimer.start();
		if ((_algorithmToUse >= 14 && _sendResultsToRoot) ? _primesSieveMPI->savePrimesToFile(_outputResultsFilename) : _primesSieve->savePrimesToFile(_outputResultsFilename)) {
			performanceTimer.stop();
			cout << "    --> Export to file " << _outputResultsFilename << " finished in " << performanceTimer.getElapsedTimeFormated() << endl;
		} else {
			cerr << "    !!!!! Export to file " << _outputResultsFilename << " failed !!!!!" << endl;
		}
	} else {
		return false;
	}

	return true;
}

void PrimesCLI::showTotalComputationTimte() {
	cout << "\n\n  > Finished program in ";
	if (_algorithmToUse >= 14) {
		int processRank;
		if (_sendResultsToRoot) {
			processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getProcessId();
		} else {
			processRank = ((PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getProcessId();
		}
		cout << "process with rank " << processRank << " in ";
	}
	stopTimer();
	cout << getElapsedTimeFormated();
	if (_algorithmToUse == 12 || _algorithmToUse == 13 || _algorithmToUse >= 16) {
		cout << " using at most " << omp_get_num_procs() << " processors and using ";
		if (_numberOfThreadsToUseInSieving != 0) {
			if (_algorithmToUse >= 16) {
				if (_sendResultsToRoot) {
					cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>*) _primesSieveMPI)->getNumberOfThreads();
				} else {
					cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads();
				}
			} else {
				cout << ((PrimesSieveParallelMultiplesOptimizedOpenMPTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel>*) _primesSieve)->getNumberOfThreads();
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
			if (!(sstream >> algorithm) || (algorithm < 1 || algorithm > 19)) {
				showUsage("  >>> Invalid algorithm selector! Must be a number [1, 19]");
				return false;
			} else {
				_algorithmToUse = algorithm;
			}
		} else if (argSelector == "--maxRangeInBits") {
			stringstream sstream(argValue);
			size_t rangeInBits;
			if (!(sstream >> rangeInBits) || (rangeInBits < 4 || rangeInBits > 64)) {
				showUsage("  >>> Invalid primes max range in bits! Max range must be >= 4 and <= 64");
				return false;
			} else {
				_primesMaxRange = (size_t) pow(2, rangeInBits);
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
			if (!(sstream >> blockSize) || (blockSize < 128)) {
				showUsage("  >>> Invalid cache block size! Block size must be >= 128");
				return false;
			} else {
				_cacheBlockSize = blockSize;
			}
		} else if (argSelector == "--dynamicSchedulingSegmentSize") {
			stringstream sstream(argValue);
			size_t dynamicScheduleBlockSize;
			if (!(sstream >> dynamicScheduleBlockSize) || (dynamicScheduleBlockSize < 256)) {
				showUsage("  >>> Invalid scheduling segment size for mpi! Block size must be >= 256");
				return false;
			} else {
				_dynamicSchedulingSegmentSizeInElements = dynamicScheduleBlockSize;
			}
		} else if (argSelector == "--dynamicSchedulingNumberSegments") {
			stringstream sstream(argValue);
			size_t dynamicScheduleNumberSegments;
			if (!(sstream >> dynamicScheduleNumberSegments) || (dynamicScheduleNumberSegments < 1)) {
				showUsage("  >>> Invalid number of segments for dynamic scheduling mpi! Block size must be >= 1");
				return false;
			} else {
				_dynamicSchedulingNumberSegments = dynamicScheduleNumberSegments;
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
		} else if (argSelector == "--outputOnlyLastSegment") {
			if (argValue == "Y" || argValue == "y") {
				_outputOnlyLastSegment = true;
			} else if (argValue == "N" || argValue == "n") {
				_outputOnlyLastSegment = false;
			} else {
				showUsage("  >>> Invalid --outputOnlyLastSegment flag! Flag must be Y or N");
				return false;
			}
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
			showUsage("  >>> Invalid arguments found!");
			return false;
		}
	}

	return true;
}

void PrimesCLI::showCurrentConfiguration() {
	cout << "  >>> Current configuration:\n";
	cout << "\t --algorithm                       -> " << _algorithmToUse << "\n";
	cout << "\t --maxRange                        -> " << _primesMaxRange << "\n";
	cout << "\t --cacheBlockSize                  -> " << _cacheBlockSize << "\n";
	cout << "\t --dynamicSchedulingSegmentSize    -> " << _dynamicSchedulingSegmentSizeInElements << "\n";
	cout << "\t --dynamicSchedulingNumberSegments -> " << _dynamicSchedulingNumberSegments << "\n";
	cout << "\t --numberThreads                   -> " << _numberOfThreadsToUseInSieving << "\n";
	cout << "\t --outputResult                    -> " << _outputResultsFilename << "\n";
	cout << "\t --outputOnlyLastSegment           -> " << _outputOnlyLastSegment << "\n";
	cout << "\t --checkResult                     -> " << _resultsConfirmationFile << "\n";
	cout << "\t --countPrimesInNode               -> " << (_countNumberOfPrimesOnNode ? "Y" : "N") << "\n";
	cout << "\t --sendPrimesCountToRoot           -> " << (_sendPrimesCountToRoot ? "Y" : "N") << "\n";
	cout << "\t --sendResultsToRoot               -> " << (_sendResultsToRoot ? "Y" : "N") << "\n";
	cout << endl;
}

void PrimesCLI::showUsage(string message) {
	if (message != "") {
		cout << message << "\n" << endl;
	}

	cout << " >>> Usage:" << endl;
	cout << "  " << _programName;
	cout << "\n    [--algorithm <number>]";
	cout << "\n    [--maxRangeInBits <number>]";
	cout << "\n    [--maxRange <number>]";
	cout << "\n    [--cacheBlockSize <number>]";
	cout << "\n    [--dynamicSchedulingSegmentSize <number>]";
	cout << "\n    [--dynamicSchedulingNumberSegments <number>]";
	cout << "\n    [--numberThreads <number>]";
	cout << "\n    [--outputResult <filename>]";
	cout << "\n    [--outputOnlyLastSegment <Y/N>]";
	cout << "\n    [--checkResult <filename>]";
	cout << "\n    [--countPrimesInNode <Y/N>]";
	cout << "\n    [--sendPrimesCountToRoot <Y/N>]";
	cout << "\n    [--sendResultsToRoot <Y/N>]";
	cout << "\n    [--help]";
	cout << "\n    [--version]" << endl;
	cout << "\t --algorithm                        -> number in [1, 19]" << endl;
	cout << "\t --maxRangeInBits                   -> number >= 4  and <= 64 (used to set maxRange using the number of bits instead of direct range -> 2^n)" << endl;
	cout << "\t --maxRange                         -> number >= 11 and <= 2^64 (default 2^32)" << endl;
	cout << "\t --cacheBlockSize                   -> block size in bytes >= 128 to optimize cache hit rate (default 16384; used in --algorithm >= 3)" << endl;
	cout
			<< "\t --dynamicSchedulingSegmentSize     -> block size in elements >= 256 to split the primes domain in blocks to perform dynamic allocation in mpi (default 1048576; used in --algorithm 16 and 19)"
			<< endl;
	cout
			<< "\t --dynamicSchedulingNumberSegments  -> number segments to use in mpi dynamic scheduling >= 1 (overrides dynamicSchedulingSegmentSize if set, default not used; applies to --algorithm 16 and 19)"
			<< endl;
	cout << "\t --numberThreads                    -> number threads to use in sieving >= 0 (default 0 -> let algorithm choose the best number of threads; used in --algorithm 12, 13, 15, 16, 18, 19)"
			<< endl;
	cout << "\t --outputResult                     -> filename of file to output results (default doesn't output results; used in all algorithms)" << endl;
	cout
			<< "\t --outputOnlyLastSegment            -> Y/N to output only the primes found on the last segment (default N; used only in mpi with dynamic scheduling - algorithms 18 and 19 - and only applies when --sendResultsToRoot is false)"
			<< endl;
	cout << "\t --checkResult                      -> filename of file with primes to check the algorithm result in root node (default doesn't check algorithm result; used in all algorithms)" << endl;
	cout << "\t --countPrimesInNode                -> Y/N to count the primes computed each node (default N; used in all algorithms)" << endl;
	cout
			<< "\t --sendPrimesCountToRoot            -> Y/N to to count the number of primes found in each mpi process and send the result to the root node (default Y and overrides the --countPrimesInNode flag if set to N; used in --algorithm >= 14)"
			<< endl;
	cout << "\t --sendResultsToRoot                -> Y/N to send the computation results to the root node (default N; used in --algorithm >= 14)" << endl;
	cout << "\t --help                             -> show program usage" << endl;
	cout << "\t --version                          -> show program version" << endl;
	cout << "\n\tWith no arguments starts interactive command line interface\n\n" << endl;

}

void PrimesCLI::showVersion() {
	cout << "Version 1.0 developed in 2013 for Parallel Computing (4th year, 2nd semester, MIEIC, FEUP)" << endl;
	cout << "Author: Carlos Miguel Correia da Costa (carlos.costa@fe.up.pt / carloscosta.cmcc@gmail.com)\n\n" << endl;
}

void PrimesCLI::showProgramHeader() {
	cout << "#######################################################################################\n";
	cout << "  >>>                           Distributed primes sieve                          <<<  \n";
	cout << "#######################################################################################\n\n";
}
