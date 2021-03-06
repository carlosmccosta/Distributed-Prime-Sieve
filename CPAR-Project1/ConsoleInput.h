#pragma once

#include <iostream>
#include <string>
#include <limits>


using std::string;
using std::getline;
using std::cin;
using std::cout;
using std::endl;


class ConsoleInput {
public:
	static ConsoleInput* getInstance();

	void getUserInput();
	string getLineCin();
	void clearConsoleScreen();
	int getIntCin(const char* message, const char* errorMessage, int min = 0, int size = INT_MAX);
	bool getYesNoCin(const char* message, const char* errorMessage = "    -> Incorrect answer! Insert Y or N!\n");

private:
	ConsoleInput(void);
	virtual ~ConsoleInput(void);
	static ConsoleInput* instance;
};
