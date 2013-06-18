#include "PrimesCLI.h"

PrimesCLI::PrimesCLI(void) {
}
PrimesCLI::~PrimesCLI(void) {
}

int main() {
	ConsoleInput::getInstance()->clearConsoleScreen();
	
	int option = 1;
	bool validOption = true;
	
	do {
		ConsoleInput::getInstance()->clearConsoleScreen();
		cout << "###############   Parallel computing - Project 2 - Carlos Costa   #####################\n";
		cout << "  >>>               Primes generator - Sieve of Eratosthenes                      <<<  \n";
		cout << "#######################################################################################\n\n\n";
		
		cout << " 1 - Single processor implementation (using division to check for primes)\n";
		cout << " 2 - Single processor implementation (using multiples to check for primes)\n";
		cout << " 3 - Single processor implementation (using block search with bitset with all numbers)\n";
		cout << " 4 - Single processor implementation (using block search with bitset only with numbers in the block)\n";
		cout << " 5 - OpenMP implementation\n";
		cout << " 6 - OpenMPI implementation\n\n";
		cout << " 0 - Exit\n\n\n" << endl;
		
		option = ConsoleInput::getInstance()->getIntCin("  >>> Option [0,6]: ", "    -> Insert one of the listed options!\n", 0, 7);
		
		if (option == 0) {
			break;
		}
		
		bool inputRangeInBits = ConsoleInput::getInstance()->getYesNoCin("\n   ## Max primes search range in bits (2^n)? (No for direct max search range specification) (Y/N): ");
		size_t inputMaxRange;
		
		if (inputRangeInBits) {
			inputMaxRange = ConsoleInput::getInstance()->getIntCin("    # n: ", "Range must be [0, 32]", 0, 33);
			inputMaxRange = (size_t) pow(2, inputMaxRange);
		} else {
			inputMaxRange = ConsoleInput::getInstance()->getIntCin("    # Max search range: ", "Range must be > 2", 3);
		}
		
		size_t blockSize = 0;
		if (option == 3 || option == 4) {
			blockSize = ConsoleInput::getInstance()->getIntCin("    # Block size in bytes: ", "Block size must be > 4", 5);
		}

		cout << "   ## Output result to file (filename, stdout or empty to avoid output): ";
		string resultFilename = ConsoleInput::getInstance()->getLineCin();
		
		cout << "   ## Confirm results from file (empty to skip confirmation): ";
		string resultConfirmationFilename = ConsoleInput::getInstance()->getLineCin();
		
		validOption = true;
		PrimesSieve<vector<bool> >* primesSieve;

		switch (option) {
			case 1: {
				primesSieve = new PrimesSieveSequencialDivision<vector<bool> >();
				break;
			}
				
			case 2: {
				primesSieve = new PrimesSieveSequencialMultiples<vector<bool> >();
				break;
			}
				
			case 3: {
				primesSieve = new PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> >(blockSize);
				break;
			}
				
			case 4: {
				primesSieve = new PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<bool> >(blockSize);
				break;
			}
				
			case 5: {

				break;
			}
				
			case 6: {
				
				break;
			}
				
			default: {
				validOption = false;
				break;
			}
		}
		
		if (validOption) {
			cout << "\n    > Computing primes from 2 to " << inputMaxRange << "..." << endl;
			primesSieve->computePrimes(inputMaxRange);
			primesSieve->extractPrimesFromBitset();
			cout << "    --> Computed " << primesSieve->getNumberPrimesFound() << " primes in " << primesSieve->getPerformanceTimer().getElapsedTimeFormated() << "\n" << endl;
			
			if (!primesSieve->getNumberPrimesFound() != 1) {
				bool validationResult;
				
				if (resultConfirmationFilename != "") {
					cout << "    > Validating computed primes with result file supplied...\n";
					validationResult = primesSieve->checkPrimesFromFile(resultConfirmationFilename);
					
					if (validationResult) {
						cout << "    --> Computed primes are correct!\n\n";
					} else {
						cout << "    --> Computed primes are different from the ones in supplied file!\n\n";
					}
				}
				
				if (resultFilename == "stdout") {
					cout << "\n=============================================  Computed primes  =============================================\n\n";
					primesSieve->printPrimesToConsole();
					cout << endl;
				} else if (resultFilename != "") {
					cout << "    > Exporting results to file " << resultFilename << "...";
					primesSieve->savePrimesToFile(resultFilename);
					cout << "\n    --> Export finished!\n";
				}
			}
			
			delete primesSieve;
			
			cout << endl << endl;
			ConsoleInput::getInstance()->getUserInput();
		}
		
	} while (option != 0);
	
	return 0;
}
