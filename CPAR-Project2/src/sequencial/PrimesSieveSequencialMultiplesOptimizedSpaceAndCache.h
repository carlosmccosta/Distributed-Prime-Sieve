#pragma once

#include "../PrimesSieve.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>

using std::sqrt;
using std::min;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimizedSpaceAndCache: public PrimesSieve<FlagsContainer> {
	private:
		size_t _blockSize;

	public:
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache(size_t blockSize = 256 * 1024) :
				_blockSize(blockSize * 8) {
		}
		
		virtual ~PrimesSieveSequencialMultiplesOptimizedSpaceAndCache() {
		}
		
		void computePrimes(size_t maxRange) {
			this->template initPrimesBitset(maxRange);
			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			
			this->template performanceTimer.reset();
			this->template performanceTimer.start();

			size_t numberBlocks = ceil((double) maxRange / (double) _blockSize);
			size_t blockBegin = 0;
			size_t blockEnd = min(_blockSize, maxRange);
			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				if (blockNumber == 0) {
					calculateFirstBlock(blockEnd, maxRangeSquareRoot);
				} else {
					blockBegin = blockEnd;
					blockEnd = blockBegin + _blockSize;
					if (blockEnd > maxRange)
						blockEnd = maxRange;

					this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBegin, blockEnd);

					size_t primeNumber = this->template getNumberAssociatedWithBitsetPosition(blockBegin);
					size_t maxPrimeNumberSearch = min(maxRangeSquareRoot, blockEnd);

					for (; primeNumber <= maxPrimeNumberSearch; primeNumber += 2) {
						// for each number not marked as composite (prime number)
						if (this->template getPrimesBitsetValue(primeNumber)) {
							//use it to calculate his composites
							size_t primeDoubled = primeNumber << 1;
							for (size_t compositeNumber = primeNumber; compositeNumber <= blockEnd; compositeNumber += primeDoubled) {
								this->template setPrimesBitsetValue(compositeNumber, false);
							}
						}
					}
				}
			}

			this->template performanceTimer.stop();
			this->template primesValues.clear();
		}

		void calculateFirstBlock(size_t blockSize, size_t maxRangeSquareRoot) {
			for (size_t primeNumber = 3; primeNumber <= maxRangeSquareRoot; primeNumber += 2) {
				// for each number not marked as composite (prime number)
				if (this->template getPrimesBitsetValue(primeNumber)) {
					//use it to calculate his composites
					size_t primeDoubled = primeNumber << 1;
					for (size_t compositeNumber = primeNumber * primeNumber; compositeNumber <= blockSize; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}
				}
			}
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBegin, size_t blockEnd) {
			for (size_t primesIndex = 0; primesIndex < blockBegin; ++primesIndex) {
				if (this->template getPrimesBitset()[primesIndex]) {
					size_t primeNumber = this->template getNumberAssociatedWithBitsetPosition(primesIndex);
					size_t primeDoubled = primeNumber << 1;
					size_t compositeNumber = PrimesUtils::closestPrimeMultiple(primeNumber, this->template getNumberAssociatedWithBitsetPosition(primesIndex));
					if (primeNumber == compositeNumber)
						compositeNumber += primeDoubled;

					for (; compositeNumber <= blockEnd; compositeNumber += primeDoubled) {
						this->template setPrimesBitsetValue(compositeNumber, false);
					}
				}
			}
		}

		inline size_t getBlockSize() {
			return _blockSize;
		}

		inline void setBlockSize(size_t blockSize) {
			_blockSize = blockSize;
		}
	};

