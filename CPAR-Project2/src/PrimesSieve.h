#pragma once

#include <stddef.h>
#include <vector>
#include <string>

using std::vector;
using std::string;

class PrimesSieve {
	protected:
		size_t maxRange;
		vector<bool> primesBitset;
		vector<int> primesValues;

	public:
		PrimesSieve() :
				maxRange(0) {
		}
		virtual ~PrimesSieve() {
		}

		virtual vector<int>& computePrimes(size_t maxRange) = 0;
		bool checkPrimesFromFile(string filename);
		bool savePrimesToFile(string filename);

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

		inline bool getPrimesBitsetValue(size_t number) {
			return primesBitset[getBitSetPositionToNumber(number)];
		}

		inline void setPrimesBitsetValue(size_t number, bool newValue) {
			primesBitset[getBitSetPositionToNumber(number)] = newValue;
		}
};

