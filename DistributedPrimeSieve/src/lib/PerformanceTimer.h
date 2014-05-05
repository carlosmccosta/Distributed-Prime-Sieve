#pragma once

#include "TimeUtils.h"

#ifdef _WIN32   // Windows system specific
#include <windows.h>
#else           // Unix based system specific
#include <sys/time.h>
#include <ctime>
#endif

#include <stdlib.h>
#include <string>
#include <sstream>

using std::string;
using std::stringstream;

class PerformanceTimer {
	public:
		PerformanceTimer();
		virtual ~PerformanceTimer();

		void start();                             // start timer
		void stop();                              // stop the timer
		void reset();
		unsigned long long getClockCounts();
		double getCPUTimeInSec();
		double getElapsedTimeInSec();             // get elapsed time in second
		double getElapsedTimeInMilliSec();        // get elapsed time in milli-second
		double getElapsedTimeInMicroSec();        // get elapsed time in micro-second
		string getElapsedTimeFormated();

		bool updateState();


	private:
		double elapsedTimeMicroSec;               // starting time in micro-second
		bool stopped;                             // stop flag
#ifdef _WIN32
		LARGE_INTEGER frequencyWin;               // ticks per second
		LARGE_INTEGER startCountWin;
		LARGE_INTEGER endCountWin;
#else
		timeval startCount;
		timeval endCount;
		clock_t clockTicksStart;
		clock_t clockTicksEnd;
#endif
		
		void calculateElapsedTimeMicroSec();
};
