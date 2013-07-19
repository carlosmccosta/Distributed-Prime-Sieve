#pragma once

#include "lib/PerformanceTimer.h"

#include <stddef.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <fstream>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::deque;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;


//typedef vector<bool> PrimesBitset;
//typedef vector<unsigned char> PrimesBitset;

template<typename FlagsContainer>
class PrimesSieve {
	protected:
		size_t _startSieveNumber;
		size_t _maxRange;
		FlagsContainer _primesBitset;
		vector<size_t> _primesValues;
		PerformanceTimer _performanceTimer;
		size_t _primesCount;

	public:
		PrimesSieve() :
				_startSieveNumber(11), _maxRange(7920), _primesCount(0) {
		}

		virtual ~PrimesSieve() {
		}

		/**
		 * Compute the number of bits to store in the bitset according to max range
		 * @param maxRange Maximum number to include in the primes search [3, maxRange]
		 * @return number of bits that the bitset must contain
		 */
		virtual inline size_t getNumberBitsToStore(size_t maxRange) {
			size_t blockSize = (maxRange + 1) - 3;
			if (blockSize % 2 == 0) {
				return (blockSize >> 1);
			} else {
				return ((blockSize >> 1) + 1);
			}
		}

		/**
		 * Prime counting function that gives an upper bound of the number of primes that are less or equal than range
		 * @param maxNumber
		 * @return
		 */
		static inline size_t getNumberOfPrimesInRange(size_t range) {
			return ((30 * log(113) / 113) * range / log(range));   // 30*log(113)/113 = 1.2550587
		}

		inline size_t getBitsetPositionToNumber(size_t number) {
			return (number - 3) >> 1;
		}

		inline size_t getNumberAssociatedWithBitsetPosition(size_t position) {
			return (position << 1) + 3;
		}

		inline bool getPrimesBitsetValue(size_t number) {
			return _primesBitset[getBitsetPositionToNumber(number)];
		}

		inline void setPrimesBitsetValue(size_t number, bool newValue) {
			_primesBitset[getBitsetPositionToNumber(number)] = newValue;
		}

		inline bool getPrimesBitsetValueBlock(size_t number, size_t blockBeginNumber) {
			return _primesBitset[(number - blockBeginNumber) >> 1];
		}

		inline void setPrimesBitsetValueBlock(size_t number, size_t blockBeginNumber, bool newValue) {
			_primesBitset[(number - blockBeginNumber) >> 1] = newValue;
		}

		virtual void computePrimes(size_t maxRange) = 0;

		virtual vector<size_t>& extractPrimesFromBitset() {
			_primesValues.clear();
			_primesValues.push_back(2);
			size_t iSize = _primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (!_primesBitset[i]) {
					_primesValues.push_back(getNumberAssociatedWithBitsetPosition(i));
				}
			}

			return _primesValues;
		}

		bool checkComputedPrimes(const vector<size_t>& expectedPrimes) {
			if (_primesValues.size() != expectedPrimes.size() || expectedPrimes.empty())
				return false;

			if (!_primesBitset.empty()) {
				extractPrimesFromBitset();
			}

			size_t iSize = _primesValues.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (_primesValues[i] != expectedPrimes[i])
					return false;
			}

			return true;
		}

		virtual bool checkPrimesFromFile(string filename) {
			if (!_primesBitset.empty() && _primesValues.empty()) {
				extractPrimesFromBitset();
			}

			ifstream inputStream(filename.c_str());

			if (inputStream.is_open()) {
				size_t numberRead;
				size_t iSize = _primesValues.size();
				size_t i = 0;
				while (inputStream >> numberRead) {
					if (i >= iSize || (numberRead != _primesValues[i])) {
						return false;
					}
					++i;
				}

				if (i != iSize) {
					return false;
				} else {
					return true;
				}
			} else {
				cerr << "    !!!!! File " << filename << " is not available !!!!!" << endl;
			}

			return false;
		}

		virtual void savePrimes(ostream& outputStream) {
			if (_primesValues.size() <= 2) {
				outputStream << 2 << endl;
				size_t iSize = _primesBitset.size();
				for (size_t i = 0; i < iSize; ++i) {
					if (!_primesBitset[i]) {
						outputStream << getNumberAssociatedWithBitsetPosition(i) << endl;
					}
				}
			} else {
				size_t iSize = _primesValues.size();
				for (size_t i = 0; i < iSize; ++i) {
					outputStream << _primesValues[i] << endl;
				}
			}
		}

		virtual bool savePrimesToFile(string filename) {
			ofstream outputStream(filename.c_str());

			if (outputStream.is_open()) {
				savePrimes(outputStream);
				return true;
			}

			return false;
		}

		virtual void printPrimesToConsole() {
			savePrimes(cout);
		}

		virtual void initPrimesBitset(size_t maxRange) {
			this->_maxRange = maxRange;
			_primesBitset = FlagsContainer(this->template getNumberBitsToStore(maxRange), false);
//			size_t iSize = _primesBitset.size();
//			for (size_t i = 0; i < iSize; ++i) {
//				_primesBitset[i] = true;
//			}
		}

		virtual void initPrimesBitSetSize(size_t newBitsetSize) {
			_primesBitset = FlagsContainer(newBitsetSize, false);
//			size_t iSize = _primesBitset.size();
//			for (size_t i = 0; i < iSize; ++i) {
//				_primesBitset[i] = true;
//			}
		}

		void resetPrimesBitsetBlock() {
			size_t iSize = _primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				_primesBitset[i] = false;
			}
		}

		inline size_t getMaxRange() const {
			return _maxRange;
		}

		inline void setMaxRange(size_t maxRange) {
			_maxRange = maxRange;
		}

		inline PerformanceTimer& getPerformanceTimer() {
			return _performanceTimer;
		}

		inline FlagsContainer& getPrimesBitset() {
			return _primesBitset;
		}

		inline vector<size_t>& getPrimesValues() {
			return _primesValues;
		}

		inline void resetPrimesValues(size_t newSize = 0) {
			_primesValues.clear();
			_primesValues.reserve(newSize);
		}

		virtual size_t getNumberPrimesFound() {
			if (_primesValues.size() >= 2)
				return _primesValues.size();

			size_t primesFound = 1;   // prime number 2 isn't in _primesBitset
			size_t iSize = _primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (!_primesBitset[i]) {
					++primesFound;
				}
			}

			return primesFound;
		}

		inline size_t getStartSieveNumber() const {
			return _startSieveNumber;
		}

		inline void setStartSieveNumber(size_t startSieveNumber) {
			_startSieveNumber = startSieveNumber;
		}

		inline size_t getPrimesCount() const {
			return _primesCount;
		}

		inline void setPrimesCount(size_t primesCount) {
			_primesCount = primesCount;
		}

		inline void incrementPrimesCount(size_t count) {
			_primesCount += count;
		}
};

