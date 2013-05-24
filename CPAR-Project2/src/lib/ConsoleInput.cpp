#include "ConsoleInput.h"


ConsoleInput::ConsoleInput() {}
ConsoleInput::~ConsoleInput() {}


void ConsoleInput::getUserInput() {
	cout << "Press ENTER to continue..." << endl;

	cin.clear();

	cin.sync();
	string temp;
	getline(cin, temp);

	cin.clear();
	cin.sync();
}


string ConsoleInput::getLineCin() {
	cin.clear();
	cin.sync();

	string input;
	getline(cin, input);

	cin.clear();
	cin.sync();

	return input;
}


void ConsoleInput::clearConsoleScreen() {
	for (size_t i = 0; i < 80; ++i) {
		cout << "\n";
	}

	cout << endl;
}


int ConsoleInput::getIntCin(const char* message, const char* errorMessage, int min, int size) {
	int number;
	do {
		cout << message;


		if (!(cin >> number)) {
			cin.clear();
			cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

			cout << errorMessage << endl;
		} else {
			if (number >= min && number < size)
				break;
			else
				cout << errorMessage << endl;
		}

	} while (true);

	cin.clear();
	cin.sync();

	return number;
}



bool ConsoleInput::getYesNoCin(const char* message, const char* errorMessage) {
	bool stop = false;
	bool incorrectOption;
	string option;

	do {
		cout << message;

		option = getLineCin();
		if ((option == "Y") || (option == "y")) {
			stop = true;
			incorrectOption = false;
		}
		else if ((option == "N") || (option == "n")) {
			stop = false;
			incorrectOption = false;
		}
		else {
			cout << errorMessage;
			incorrectOption = true;
		}
	} while (incorrectOption);

	cin.clear();
	cin.sync();

	return stop;
}

ConsoleInput* ConsoleInput::getInstance() {
	if (instance == NULL) {
		instance = new ConsoleInput();
	}

	return instance;
}

ConsoleInput* ConsoleInput::instance = NULL;
