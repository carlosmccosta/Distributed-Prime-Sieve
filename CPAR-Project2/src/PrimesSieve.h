#pragma once

#include "lib/PerformanceTimer.h"

#include <stddef.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;

class PrimesSieve {
	protected:
		size_t maxRange;
		vector<bool> primesCompositesBitset;
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
		vector<size_t>& extractPrimesFromBitset();

		bool checkComputedPrimes(const vector<size_t>& expectedPrimes);
		bool checkPrimesFromFile(string filename);
		bool savePrimesToFile(string filename);
		void printPrimesToConsole();

		inline bool getPrimesBitsetValue(size_t number) {
			return primesCompositesBitset[getBitSetPositionToNumber(number)];
		}
		
		inline void setPrimesBitsetValue(size_t number, bool newValue) {
			primesCompositesBitset[getBitSetPositionToNumber(number)] = newValue;
		}
		
		void initPrimesCompositeBitset(size_t maxRange);

		size_t getMaxRange() const {
			return maxRange;
		}
		
		PerformanceTimer& getPerformanceTimer() {
			return performanceTimer;
		}
		
		const vector<bool>& getPrimesCompositesBitset() const {
			return primesCompositesBitset;
		}
		
		const vector<size_t>& getPrimesValues() const {
			return primesValues;
		}
		
		size_t getNumberPrimesFound() const {
			return primesValues.size();
		}
};

