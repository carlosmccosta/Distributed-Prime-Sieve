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
		FlagsContainer primesCompositesBitset;
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
		
		static inline size_t getBitSetPositionToNumber(size_t number) {
			return (number - 3) >> 1;
		}
		
		static inline size_t getNumberAssociatedWithBitSetPosition(size_t number) {
			return (number << 1) + 3;
		}
		
		virtual void computePrimes(size_t maxRange) = 0;

		vector<size_t>& extractPrimesFromBitset() {
			primesValues.clear();
			primesValues.push_back(2);
			size_t iSize = primesCompositesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				if (!primesCompositesBitset[i]) {
					primesValues.push_back(getNumberAssociatedWithBitSetPosition(i));
				}
			}
			
			return primesValues;
		}
		
		bool checkComputedPrimes(const vector<size_t>& expectedPrimes) {
			if (primesValues.size() != expectedPrimes.size() || expectedPrimes.empty())
				return false;
			
			if (!primesCompositesBitset.empty()) {
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
			if (!primesCompositesBitset.empty() && primesValues.empty()) {
				extractPrimesFromBitset();
			}
			
			ifstream inputStream(filename.c_str());
			
			if (inputStream.is_open()) {
				size_t numberRead;
				size_t iSize = primesValues.size();
				size_t i = 0;
				while (inputStream >> numberRead) {
					if (i >= iSize || numberRead != primesValues[i]) {
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
				outputStream << 2 << endl;
				if (primesValues.empty()) {
					size_t iSize = primesCompositesBitset.size();
					for (size_t i = 0; i < iSize; ++i) {
						if (!primesCompositesBitset[i]) {
							outputStream << getNumberAssociatedWithBitSetPosition(i) << endl;
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
			cout << 2 << endl;
			if (primesValues.empty()) {
				size_t iSize = primesCompositesBitset.size();
				for (size_t i = 0; i < iSize; ++i) {
					if (!primesCompositesBitset[i]) {
						cout << getNumberAssociatedWithBitSetPosition(i) << endl;
					}
				}
			} else {
				size_t iSize = primesValues.size();
				for (size_t i = 0; i < iSize; ++i) {
					cout << primesValues[i] << endl;
				}
			}
		}
		
		inline bool getPrimesBitsetValue(size_t number) {
			return primesCompositesBitset[getBitSetPositionToNumber(number)];
		}
		
		inline void setPrimesBitsetValue(size_t number, bool newValue) {
			primesCompositesBitset[getBitSetPositionToNumber(number)] = newValue;
		}
		
		void initPrimesCompositeBitset(size_t maxRange) {
			this->maxRange = maxRange;
			primesCompositesBitset = FlagsContainer(getNumberBitsToStore(maxRange));
			size_t iSize = primesCompositesBitset.size();
			for (size_t i = 0; i < iSize; ++i) {
				primesCompositesBitset[i] = false;
			}
		}
		
		size_t getMaxRange() const {
			return maxRange;
		}
		
		PerformanceTimer& getPerformanceTimer() {
			return performanceTimer;
		}
		
		const FlagsContainer& getPrimesCompositesBitset() const {
			return primesCompositesBitset;
		}
		
		const vector<size_t>& getPrimesValues() const {
			return primesValues;
		}
		
		size_t getNumberPrimesFound() const {
			return primesValues.size();
		}
		
};

