#pragma once

#include <stddef.h>
#include <limits>

struct WheelElement {
		unsigned char wheelPositionIndex;
		unsigned char nextPossiblePrimeIncrement;

};

extern const WheelElement wheel30Elements[30];
extern const WheelElement wheel210Elements[210];

template<size_t WheelSize, size_t NumberOfPossiblePrimesPerWheel, size_t FirstPrimeToSieve, const WheelElement* WheelElements>
class WheelFactorization {
	protected:
		size_t _wheelSize;
		size_t _numberOfPossiblePrimesPerWheel;
		size_t _firstPrimeToSieve;
		const WheelElement* _wheelElements;

	public:
		WheelFactorization() :
				_wheelSize(WheelSize),
				_numberOfPossiblePrimesPerWheel(NumberOfPossiblePrimesPerWheel),
				_firstPrimeToSieve(FirstPrimeToSieve),
				_wheelElements(WheelElements) {
		}

		virtual ~WheelFactorization() {
		}

		inline size_t getNextPossiblePrime(size_t number) {
			return number + _wheelElements[number % _wheelSize].nextPossiblePrimeIncrement;
		}

		inline size_t getBitsetPositionToNumber(size_t number) {
			size_t wheelElementPosition = number % _wheelSize;
			size_t positionInsideWheel = _wheelElements[wheelElementPosition].wheelPositionIndex;

			if (positionInsideWheel == 0xFF) {
				return std::numeric_limits<std::size_t>::max();
			}

			size_t wheelIndex = ((number / _wheelSize) * _numberOfPossiblePrimesPerWheel);
			return wheelIndex + positionInsideWheel;
		}

		inline size_t getNumberBitsToStore(size_t maxRange) {
			size_t numberBitsToStore = getBitsetPositionToNumber(maxRange);
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
};

/// 3rd wheel, skips multiples of 2, 3 and 5
typedef WheelFactorization<(size_t)30, (size_t)8, (size_t)7, wheel30Elements> Modulo30Wheel;

/// 4th wheel, skips multiples of 2, 3, 5 and 7
typedef WheelFactorization<(size_t)210, (size_t)48, (size_t)11, wheel210Elements> Modulo210Wheel;
