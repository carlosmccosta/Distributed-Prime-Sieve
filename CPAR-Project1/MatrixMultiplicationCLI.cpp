#include "MatrixMultiplicationCLI.h"

MatrixMultiplicationCLI::MatrixMultiplicationCLI(void) {}
MatrixMultiplicationCLI::~MatrixMultiplicationCLI(void) {}


int main() {
	ConsoleInput::getInstance()->clearConsoleScreen();

	int option = 1;
	bool roundLoaded = false;

	do {
		ConsoleInput::getInstance()->clearConsoleScreen();
		cout << "#########################   Parallel computing - Project 1   ##########################\n";
		cout << "  >>>                           Matrix multiplication                             <<<  \n";
		cout << "#######################################################################################\n\n\n";

		cout << " 1 - Basic matrix multiplication\n";
		cout << " 2 - Line matrix multiplication\n";
		cout << " 3 - Block matrix multiplication\n\n";
		cout << " 0 - Exit\n\n\n" << endl;


		option = ConsoleInput::getInstance()->getIntCin("  >>> Option: ", "Insert one of the listed options!\n", 0, 4);

		if (option == 0) {
			break;
		}

		bool validOption = true;
		bool loadFromFiles = false;

		BidimensionalMatrix *leftMatrix, *rightMatrix;
		MatrixMultiplication* matrixMultAlgorithm;

		string leftMatrixFilename;
		string rightMatrixFilename;
		string resultFile = "";
		unsigned int blockColumnSize = 1, blockLineSize = 1;



		loadFromFiles = ConsoleInput::getInstance()->getYesNoCin("\n   ## Load matrices from files? (Y/N): ");
		if (loadFromFiles) {
			cout << "    # Left matrix file: ";
			leftMatrixFilename = ConsoleInput::getInstance()->getLineCin();

			cout << "    # Right matrix file: ";
			rightMatrixFilename = ConsoleInput::getInstance()->getLineCin();

			bool loadResultFile = ConsoleInput::getInstance()->getYesNoCin("\n   ## Load result file to compare? (Y/N): ");
			if (loadResultFile) {
				cout << "    # Result matrix file: ";
				resultFile = ConsoleInput::getInstance()->getLineCin();
			}

			leftMatrix = new BidimensionalMatrix();
			leftMatrix->initializeMatrixFromFile(leftMatrixFilename);
			rightMatrix = new BidimensionalMatrix();
			rightMatrix->initializeMatrixFromFile(rightMatrixFilename);
		} else {
			cout << "    > Using default initialization...\n\n";
			unsigned int numberOfColumnsOfLeftMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of columns of left matrix: ", "Insert a number > 1!\n", 2, INT_MAX);
			unsigned int numberOfLinesOfLeftMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of lines of left matrix: ", "Insert a number > 1!\n", 2, INT_MAX);
			unsigned int numberOfColumnsOfRightMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of columns of right matrix: ", "Insert a number > 1!\n", 2, INT_MAX);


			leftMatrix = new BidimensionalMatrix(numberOfColumnsOfLeftMatrix, numberOfLinesOfLeftMatrix);
			leftMatrix->initializeMatrix();
			rightMatrix = new BidimensionalMatrix(numberOfColumnsOfRightMatrix, numberOfColumnsOfLeftMatrix);
			rightMatrix->initializeMatrixWithSequenceOnLines();
		}


		switch (option) {
		case 1: {
			matrixMultAlgorithm = new MatrixMultiplicationBasic();
			break;
				}

		case 2: {
			matrixMultAlgorithm = new MatrixMultiplicationLine();
			break;
				}

		case 3: {
			blockColumnSize = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Block column size: ", "Insert a number > 1!\n", 2, INT_MAX);
			blockLineSize = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Block line size: ", "Insert a number > 1!\n", 2, INT_MAX);
			matrixMultAlgorithm = new MatrixMultiplicationBlock();
			break;
				}

		default: {
			validOption = false;
			break;
				 }
		}

		if (validOption) {
			cout << endl << "\n    > Multiplying matrices...";
			shared_ptr<BidimensionalMatrix> resultMatrix = matrixMultAlgorithm->performMultiplication(*leftMatrix, *rightMatrix, blockColumnSize, blockLineSize);
			cout << "\n    -> Multiplication finished in " << matrixMultAlgorithm->getPerformanceTimer().getElapsedTimeInSec() <<  " seconds\n" << endl;

			bool validationResult;

			if (loadFromFiles) {				
				cout << "    > Validating result matrix with result file supplied...\n";
				validationResult = resultMatrix->validateResultFromFile(resultFile);
			} else {
				cout << "    > Validating result matrix from default initialization (result should be sum of powers)...\n";
				validationResult = resultMatrix->validateResultOfDefaultMatrixInitialization(rightMatrix->getNumberLines());
			}


			if (validationResult) {
				cout << "    -> Result matrix validated!\n\n";
			} else {
				cout << "    -> Result matrix incorrect!\n\n";
			}

			int exportType = ConsoleInput::getInstance()->getIntCin("   ## Export matrices? (1-Only result, 2-All, 0-None):  ", "Insert the number of the option!\n", 0, 3);


			if (exportType == 2) {
				leftMatrix->exportMatrixToFile("leftMatrix.txt");
				cout << "    > Exporting left matrix to leftMatrix.txt...\n";
				rightMatrix->exportMatrixToFile("rightMatrix.txt");
				cout << "    > Exporting right matrix to rightMatrix.txt...\n";
			}

			if (exportType == 1 || exportType == 2) {
				cout << "    > Exporting result matrix to resultMatrix.txt...\n";
				resultMatrix->exportMatrixToFile("resultMatrix.txt");
			}

			cout << "    -> Export finished!\n";

			delete leftMatrix;
			delete rightMatrix;
			resultMatrix.reset();
			delete matrixMultAlgorithm;
			cout << endl << endl;
			ConsoleInput::getInstance()->getUserInput();
		}

	} while (option != 0);

	return 0;
}
