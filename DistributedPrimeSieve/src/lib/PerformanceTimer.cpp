#include "PerformanceTimer.h"

PerformanceTimer::PerformanceTimer() {
	reset();
}

PerformanceTimer::~PerformanceTimer() {
}

void PerformanceTimer::start() {
	stopped = 0;
#ifdef _WIN32
	QueryPerformanceCounter(&startCountWin);
#else
	gettimeofday(&startCount, NULL);
	clockTicksStart = clock();
#endif
}

void PerformanceTimer::stop() {
	stopped = true;
	
#ifdef _WIN32
	QueryPerformanceCounter(&endCountWin);
#else
	gettimeofday(&endCount, NULL);
	clockTicksEnd = clock();
#endif
	
	calculateElapsedTimeMicroSec();
}

void PerformanceTimer::reset() {
#ifdef _WIN32
	QueryPerformanceFrequency(&frequencyWin);
	startCountWin.QuadPart = 0;
	endCountWin.QuadPart = 0;
#else
	startCount.tv_sec = startCount.tv_usec = 0;
	endCount.tv_sec = endCount.tv_usec = 0;
	clockTicksStart = 0;
	clockTicksEnd = 0;
#endif
	
	stopped = false;
	elapsedTimeMicroSec = 0;
}


unsigned long long PerformanceTimer::getClockCounts() {
	updateState();

#ifdef _WIN32
	return (unsigned long long) (endCountWin - startCountWin);
#else
	return (unsigned long long) (clockTicksEnd - clockTicksStart);
#endif
}

double PerformanceTimer::getCPUTimeInSec() {
	updateState();

#ifdef _WIN32
	return getElapsedTimeInSec(); // todo for windows
#else
	double cpuTicks = clockTicksEnd - clockTicksStart;
	return cpuTicks / (double)CLOCKS_PER_SEC;
#endif
}


double PerformanceTimer::getElapsedTimeInSec() {
	return getElapsedTimeInMicroSec() / 1000000.0;
}

double PerformanceTimer::getElapsedTimeInMilliSec() {
	
	return getElapsedTimeInMicroSec() / 1000.0;
}

double PerformanceTimer::getElapsedTimeInMicroSec() {
	updateState();
	
	return elapsedTimeMicroSec;
}

void PerformanceTimer::calculateElapsedTimeMicroSec() {
#ifdef _WIN32
	elapsedTimeMicroSec = (endCountWin.QuadPart - startCountWin.QuadPart) * 1000000.0 / frequencyWin.QuadPart;
#else
	elapsedTimeMicroSec = (double)((endCount.tv_sec - startCount.tv_sec) * 1000000 + (endCount.tv_usec - startCount.tv_usec));
#endif
	
}

string PerformanceTimer::getElapsedTimeFormated() {
	stringstream timeFormated;
	timeFormated << "[ Wall time: " << TimeUtils::formatSecondsToDate(getElapsedTimeInSec());

#ifdef __linux__

	timeFormated << " | CPU time: " << TimeUtils::formatSecondsToDate(getCPUTimeInSec());
#endif

	timeFormated << " ]";

	return timeFormated.str();
}


bool PerformanceTimer::updateState() {
	if (!stopped) {
#ifdef _WIN32
		QueryPerformanceCounter(&endCountWin);
#else
		gettimeofday(&endCount, NULL);
		clockTicksEnd = clock();
#endif
		calculateElapsedTimeMicroSec();
		return true;
	} else {
		return false;
	}
}
