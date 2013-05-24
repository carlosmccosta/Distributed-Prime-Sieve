#pragma once

#ifdef _WIN32   // Windows system specific
#include <windows.h>
#else           // Unix based system specific
#include <sys/time.h>
#endif

#include <stdlib.h>

class PerformanceTimer {
	public:
		PerformanceTimer();
		virtual ~PerformanceTimer();

		void start();                             // start timer
		void stop();                              // stop the timer
		void reset();
		double getElapsedTimeInSec();               // get elapsed time in second
		double getElapsedTimeInMilliSec();          // get elapsed time in milli-second
		double getElapsedTimeInMicroSec();          // get elapsed time in micro-second
		
	private:
		double startTimeInMicroSec;                 // starting time in micro-second
		double endTimeInMicroSec;                   // ending time in micro-second
		int stopped;                             // stop flag 
#ifdef _WIN32
		LARGE_INTEGER frequency;                    // ticks per second
		LARGE_INTEGER startCount;
		LARGE_INTEGER endCount;
#else
		timeval startCount;
		timeval endCount;
#endif
		
		double getElapsedTime(double timeConversionOffset);
};
