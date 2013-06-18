#pragma once

#include "lib/PerformanceTimer.h"

#include <stddef.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <fstream>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::deque;
using std::string;
using std::ofstream;
using std::ifstream;

//typedef vector<bool> PrimesBitset;
//typedef vector<unsigned char> PrimesBitset;

template<typename FlagsContainer>
class PrimesSieve {
	protected:
		size_t maxRange;
		FlagsContainer primesBitset;
		vector<size_t> primesValues;
		PerformanceTimer performanceTimer;

	public:
		PrimesSieve() :
				maxRange(0) {
		}
		
		virtual ~PrimesSieve() {
		}
		
		/**
		 * Compute the number of bits to store in the bitset according to max range
		 * @param maxRange Maximum number to include in the primes search [3, maxRange]
		 * @return number of bits that the bitset must contain
		 */
		static inline size_t getNumberBitsToStore(size_t maxRange) {
			return ((maxRange - 3) >> 1) + 1;
		}
		
		static inline size_t getBitsetPositionToNumber(size_t number) {
			return (number - 3) >> 1;
		}
		
		static inline size_t getNumberAssociatedWithBitsetPosition(size_t position) {
			return (position << 1) + 3;
		}
		
		inline bool getPrimesBitsetValue(size_t number) {
			return primesBitset[getBitsetPositionToNumber(number)];
		}

		inline void setPrimesBitsetValue(size_t number, bool newValue) {
			primesBitset[getBitsetPositionToNumber(number)] = newValue;
		}

		inline bool getPrimesBitsetValueBlock(size_t number, size_t blockBeginNumber) {
			return primesBitset[(number - blockBeginNumber) >> 1];
		}

		inline void setPrimesBitsetValueBlock(size_t number, size_t blockBeginNumber, bool newValue) {
			primesBitset[(number - blockBeginNumber) >> 1] = newValue;
		}

		virtual void computePrimes(size_t maxRange) = 0;

		virtual vector<size_t>& extractPrimesFromBitset() {
			primesValues.clear();
			primesValues.push_back(2);
			size_t iSize = primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (primesBitset[i]) {
					primesValues.push_back(getNumberAssociatedWithBitsetPosition(i));
				}
			}
			
			return primesValues;
		}
		
		bool checkComputedPrimes(const vector<size_t>& expectedPrimes) {
			if (primesValues.size() != expectedPrimes.size() || expectedPrimes.empty())
				return false;
			
			if (!primesBitset.empty()) {
				extractPrimesFromBitset();
			}
			
			size_t iSize = primesValues.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (primesValues[i] != expectedPrimes[i])
					return false;
			}
			
			return true;
		}
		
		bool checkPrimesFromFile(string filename) {
			if (!primesBitset.empty() && primesValues.empty()) {
				extractPrimesFromBitset();
			}
			
			ifstream inputStream(filename.c_str());
			
			if (inputStream.is_open()) {
				size_t numberRead;
				size_t iSize = primesValues.size();
				size_t i = 0;
				while (inputStream >> numberRead) {
					if (i >= iSize || (numberRead != primesValues[i])) {
						return false;
					}
					++i;
				}
				
				return true;
			} else {
				cerr << "    -> File " << filename << " is not available!" << endl;
			}
			
			return false;
		}
		
		bool savePrimesToFile(string filename) {
			ofstream outputStream(filename.c_str());
			
			if (outputStream.is_open()) {
				if (primesValues.empty()) {
					outputStream << 2 << endl;
					size_t iSize = primesBitset.size();
					for (size_t i = 0; i < iSize; ++i) {
						if (primesBitset[i]) {
							outputStream << getNumberAssociatedWithBitsetPosition(i) << endl;
						}
					}
				} else {
					size_t iSize = primesValues.size();
					for (size_t i = 0; i < iSize; ++i) {
						outputStream << primesValues[i] << endl;
					}
				}
				return true;
			}
			
			return false;
		}
		
		void printPrimesToConsole() {
			if (primesValues.empty()) {
				cout << 2 << endl;
				size_t iSize = primesBitset.size();
				for (size_t i = 0; i < iSize; ++i) {
					if (primesBitset[i]) {
						cout << getNumberAssociatedWithBitsetPosition(i) << endl;
					}
				}
			} else {
				size_t iSize = primesValues.size();
				for (size_t i = 0; i < iSize; ++i) {
					cout << primesValues[i] << endl;
				}
			}
		}
		
		void initPrimesBitset(size_t maxRange) {
			this->maxRange = maxRange;
			primesBitset = FlagsContainer(getNumberBitsToStore(maxRange));
			size_t iSize = primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				primesBitset[i] = true;
			}
		}
		
		void initPrimesBitsetBlock(size_t blockSize) {
			this->maxRange = blockSize;
			primesBitset = FlagsContainer(blockSize);
			size_t iSize = primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				primesBitset[i] = true;
			}
		}

		void resetPrimesBitsetBlock() {
			size_t iSize = primesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				primesBitset[i] = true;
			}
		}

		inline size_t getMaxRange() const {
			return maxRange;
		}
		
		inline PerformanceTimer& getPerformanceTimer() {
			return performanceTimer;
		}
		
		inline FlagsContainer& getPrimesBitset() {
			return primesBitset;
		}
		
		inline vector<size_t>& getPrimesValues() {
			return primesValues;
		}
		
		inline size_t getNumberPrimesFound() const {
			return primesValues.size();
		}
};
