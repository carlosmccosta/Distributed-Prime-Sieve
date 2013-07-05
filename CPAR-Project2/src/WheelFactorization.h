#pragma once

#include <stddef.h>
#include <limits>
#include <vector>

using std::vector;

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
				_wheelSize(WheelSize),
				_numberOfPossiblePrimesPerWheel(NumberOfPossiblePrimesPerWheel),
				_firstPrimeToSieve(FirstPrimeToSieve),
				_numberPrimesSievedByTheWheel(NumberPrimesSievedByTheWheel),
				_wheelElements(WheelElements) {
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
				return _wheelElements[(number - 1) % _wheelSize].wheelPositionIndex != _wheelElements[number % _wheelSize].wheelPositionIndex;
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
typedef WheelFactorization<vector<bool>, (size_t)30, (size_t)8, (size_t)7, (size_t)3, wheel30Elements> Modulo30Wheel;

/// 4th wheel, skips multiples of 2, 3, 5 and 7
typedef WheelFactorization<vector<bool>, (size_t)210, (size_t)48, (size_t)11, (size_t)4, wheel210Elements> Modulo210Wheel;
