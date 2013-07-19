#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true,
				bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel() {
		}

		inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			if (blockSize % 2 == 0) {
				return (blockSize >> 1);
			} else {
				return ((blockSize >> 1) + 1);
			}
		}

		inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return (position << 1) + this->template getStartSieveNumber();
		}
};

