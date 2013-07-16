#pragma once

#include <stddef.h>
#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ifstream;

#include <boost/dynamic_bitset.hpp>

//typedef boost::dynamic_bitset<> PrimesFlagsContainer;
typedef vector<bool> PrimesFlagsContainer;
//typedef vector<unsigned char> PrimesFlagsContainer;
typedef vector<unsigned char> PrimesFlagsContainerMPI;


struct WheelElement {
		unsigned char wheelPositionIndex;
		unsigned char nextPossiblePrimeIncrement;
};

extern const WheelElement wheel30Elements[30];
extern const WheelElement wheel210Elements[210];

template<typename FlagsContainer, size_t WheelSize, size_t NumberOfPossiblePrimesPerWheel, size_t FirstPrimeToSieve, size_t NumberPrimesSievedByTheWheel, const WheelElement* WheelElements>
class WheelFactorization {
	protected:
		size_t _wheelSize;
		size_t _numberOfPossiblePrimesPerWheel;
		size_t _firstPrimeToSieve;
		size_t _numberPrimesSievedByTheWheel;
		const WheelElement* _wheelElements;

	public:
		WheelFactorization() :
				_wheelSize(WheelSize), _numberOfPossiblePrimesPerWheel(NumberOfPossiblePrimesPerWheel), _firstPrimeToSieve(FirstPrimeToSieve), _numberPrimesSievedByTheWheel(NumberPrimesSievedByTheWheel), _wheelElements(WheelElements) {
		}

		virtual ~WheelFactorization() {
		}

		inline size_t getNextPossiblePrime(size_t number) {
			return number + _wheelElements[number % _wheelSize].nextPossiblePrimeIncrement;
		}

		inline size_t getBitsetPositionToNumber(size_t number) {
			return ((number / _wheelSize) * _numberOfPossiblePrimesPerWheel) + _wheelElements[number % _wheelSize].wheelPositionIndex;
		}

		void setBitsetPositionToNumber(FlagsContainer& flagsContainer, size_t number, bool newValue) {
			unsigned char positionInsideWheel = _wheelElements[number % _wheelSize].wheelPositionIndex;
			if (positionInsideWheel != 0xFF) {
				flagsContainer[((number / _wheelSize) * _numberOfPossiblePrimesPerWheel) + positionInsideWheel] = newValue;
			}
		}

		size_t getBitsetPositionToNumberWithCheck(size_t number) {
			unsigned char positionInsideWheel = _wheelElements[number % _wheelSize].wheelPositionIndex;
			if (positionInsideWheel == 0xFF) {
				return std::numeric_limits<std::size_t>::max();
			}

			return ((number / _wheelSize) * _numberOfPossiblePrimesPerWheel) + positionInsideWheel;
		}

		bool isNumberPossiblePrime(size_t number) {
			if (number > 0) {
				return (_wheelElements[(number - 1) % _wheelSize].nextPossiblePrimeIncrement != (_wheelElements[number % _wheelSize].nextPossiblePrimeIncrement + 1));
			}
			return false;
		}

		inline size_t getNumberBitsToStore(size_t maxRange) {
			size_t numberBitsToStore = getBitsetPositionToNumberWithCheck(maxRange);
			if (numberBitsToStore == std::numeric_limits<std::size_t>::max()) {
				size_t nextPossiblePrime = getNextPossiblePrime(maxRange);
				numberBitsToStore = getBitsetPositionToNumber(nextPossiblePrime);
			} else {
				++numberBitsToStore;
			}

			return numberBitsToStore;
		}

		virtual bool checkPrimesFromFileMPI(FlagsContainer& primesBitset, size_t startSieveNumber, string filename) {
			ifstream inputStream(filename.c_str());

			if (inputStream.is_open()) {
				size_t numberRead;

				size_t primesBitsetSize = primesBitset.size();

				size_t possiblePrime = startSieveNumber;
				size_t possiblePrimePosition = getBitsetPositionToNumberMPI(possiblePrime, startSieveNumber);

				if (!(isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = getNextPossiblePrime(possiblePrime);
				}

				size_t numberPrimesFromFile = 0;
				while (inputStream >> numberRead) {
					++numberPrimesFromFile;

					if (numberRead < getFirstPrimeToSieve()) {
						continue;
					}

					if (possiblePrimePosition >= primesBitsetSize || primesBitset[possiblePrimePosition] == false) {
						return false;
					}
					possiblePrime = getNextPossiblePrime(possiblePrime);
					possiblePrimePosition = getBitsetPositionToNumberMPI(possiblePrime, startSieveNumber);
				}

//				if (possiblePrimePosition < primesBitsetSize && primesBitset[possiblePrimePosition]) {
//					return false;
//				} else {
//					return true;
//				}
			} else {
				cerr << "    !!!!! File " << filename << " is not available !!!!!" << endl;
			}

			return false;
		}

		inline size_t getBitsetPositionToNumberMPI(size_t number, size_t startSieveNumber) {
			return (number - startSieveNumber) >> 1;
		}

		inline size_t getWheelSize() const {
			return _wheelSize;
		}

		size_t getFirstPrimeToSieve() const {
			return _firstPrimeToSieve;
		}

		size_t getNumberOfPossiblePrimesPerWheel() const {
			return _numberOfPossiblePrimesPerWheel;
		}

		WheelElement* getWheelElements() const {
			return _wheelElements;
		}

		size_t getNumberPrimesSievedByTheWheel() const {
			return _numberPrimesSievedByTheWheel;
		}
};


/// 3rd wheel, skips multiples of 2, 3 and 5
typedef WheelFactorization<PrimesFlagsContainer, (size_t) 30, (size_t) 8, (size_t) 7, (size_t) 3, wheel30Elements> Modulo30Wheel;
//typedef WheelFactorization<FlagsContainerMPI, (size_t) 30, (size_t) 8, (size_t) 7, (size_t) 3, wheel30Elements> Modulo30Wheel;

/// 4th wheel, skips multiples of 2, 3, 5 and 7
typedef WheelFactorization<PrimesFlagsContainer, (size_t) 210, (size_t) 48, (size_t) 11, (size_t) 4, wheel210Elements> Modulo210Wheel;
//typedef WheelFactorization<FlagsContainerMPI, (size_t) 210, (size_t) 48, (size_t) 11, (size_t) 4, wheel210Elements> Modulo210Wheel;
typedef WheelFactorization<PrimesFlagsContainerMPI, (size_t) 210, (size_t) 48, (size_t) 11, (size_t) 4, wheel210Elements> Modulo210WheelByte;
