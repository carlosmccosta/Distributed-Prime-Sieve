#include "MatrixMultiplicationCLI.h"




MatrixMultiplicationCLI::MatrixMultiplicationCLI(void)
{
}


MatrixMultiplicationCLI::~MatrixMultiplicationCLI(void)
{
}



int main( int argc, char *argv[] )
{
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
		bool validOption = true;
		bool loadFromFiles = false;

		BidimensionalMatrix leftMatrix, rightMatrix;
		MatrixMultiplication* matrixMultAlgorithm;

		switch (option) {
		case 1: {
			loadFromFiles = ConsoleInput::getInstance()->getYesNoCin("   ## Load matrices from files? (Y/N): ");
			if (loadFromFiles) {

			} else {
				cout << "    > Using default initialization...\n\n";
				unsigned int numberOfColumnsOfLeftMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of columns of left matrix: ", "Insert a number > 1!\n", 0, INT_MAX);
				unsigned int numberOfLinesOfLeftMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of lines of left matrix: ", "Insert a number > 1!\n", 0, INT_MAX);
				unsigned int numberOfColumnsOfRightMatrix = (unsigned int)ConsoleInput::getInstance()->getIntCin("   ## Number of columns of right matrix: ", "Insert a number > 1!\n", 0, INT_MAX);


				leftMatrix = BidimensionalMatrix(numberOfColumnsOfLeftMatrix, numberOfLinesOfLeftMatrix);
				leftMatrix.initializeMatrix();
				rightMatrix = BidimensionalMatrix(numberOfColumnsOfRightMatrix, numberOfColumnsOfLeftMatrix);
				rightMatrix.initializeMatrixWithSequenceOnLines();

				matrixMultAlgorithm = new MatrixMultiplicationBasic();
			}

			break;
				}

		case 2: {

			break;
				}

		case 3: {

			break;
				}

		default: {
			validOption = false;
			break;
				 }
		}

		if (validOption) {
			cout << endl << "\n    > Multiplying matrices...";
			shared_ptr<BidimensionalMatrix> resultMatrix = matrixMultAlgorithm->performMultiplication(leftMatrix, rightMatrix);
			cout << "\n    -> Multiplication finished in " << matrixMultAlgorithm->getPerformanceTimer().getElapsedTimeInSec() <<  " seconds\n" << endl;

			if (loadFromFiles) {				
				bool exportResult = ConsoleInput::getInstance()->getYesNoCin("   # Export result matrix? (Y/N): ");
			} else {
				cout << "    > Validating result matrix from default initialization (result should be sum of powers)...\n";
				bool validationResult = resultMatrix->validateResultOfDefaultMatrixInitialization(rightMatrix.getNumberLines());

				if (validationResult) {
					cout << "    -> Result matrix validated!\n\n";

					int exportType = ConsoleInput::getInstance()->getIntCin("   ## Export matrices? (1-Only result, 2-All, 0-None):  ", "Insert the number of the option!\n", 0, INT_MAX);


					if (exportType == 2) {
						leftMatrix.exportMatrixToFile("leftMatrix.txt");
						cout << "    > Exporting left matrix to leftMatrix.txt...\n";
						rightMatrix.exportMatrixToFile("rightMatrix.txt");
						cout << "    > Exporting right matrix to rightMatrix.txt...\n";
					}

					if (exportType == 1 || exportType == 2) {
						cout << "    > Exporting result matrix to resultMatrix.txt...\n";
						resultMatrix->exportMatrixToFile("resultMatrix.txt");
					}

					cout << "    -> Export finished!\n";
				} else {
					cout << "    -> Result matrix incorrect!\n";
				}

			}

			delete matrixMultAlgorithm;
			cout << endl << endl;
			ConsoleInput::getInstance()->getUserInput();
		}

	} while (option != 0);

	return 0;
}
