#include "PerformanceTimer.h"

PerformanceTimer::PerformanceTimer() {
	reset();
}

PerformanceTimer::~PerformanceTimer() {
}

void PerformanceTimer::start() {
	stopped = 0;
#ifdef _WIN32
	QueryPerformanceCounter(&startCount);
#else
	gettimeofday(&startCount, NULL);
#endif
}

void PerformanceTimer::stop() {
	stopped = 1;

#ifdef _WIN32
	QueryPerformanceCounter(&endCount);
#else
	gettimeofday(&endCount, NULL);
#endif
}

void PerformanceTimer::reset() {
#ifdef _WIN32
	QueryPerformanceFrequency(&frequency);
	startCount.QuadPart = 0;
	endCount.QuadPart = 0;
#else
	startCount.tv_sec = startCount.tv_usec = 0;
	endCount.tv_sec = endCount.tv_usec = 0;
#endif
	stopped = 0;
	startTimeInMicroSec = 0;
	endTimeInMicroSec = 0;
}

double PerformanceTimer::getElapsedTimeInSec() {
	return this->getElapsedTime(1.0);
}

double PerformanceTimer::getElapsedTimeInMilliSec() {

	return this->getElapsedTime(1000.0);
}

double PerformanceTimer::getElapsedTimeInMicroSec() {
	return this->getElapsedTime(1000000.0);
}

double PerformanceTimer::getElapsedTime(double timeConversionOffset) {
#ifdef _WIN32
	if(!stopped) {
		QueryPerformanceCounter(&endCount);
	}

	startTimeInMicroSec = startCount.QuadPart * (timeConversionOffset / frequency.QuadPart);
	endTimeInMicroSec = endCount.QuadPart * (timeConversionOffset / frequency.QuadPart);
#else
	if (!stopped)
		gettimeofday(&endCount, NULL);

	startTimeInMicroSec = (startCount.tv_sec * timeConversionOffset) + startCount.tv_usec;
	endTimeInMicroSec = (endCount.tv_sec * timeConversionOffset) + endCount.tv_usec;
#endif

	return endTimeInMicroSec - startTimeInMicroSec;
}
