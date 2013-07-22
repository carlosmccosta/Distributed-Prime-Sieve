#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true,
				bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel() {
		}

		virtual inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			return blockSize;
		}

		virtual inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber());
		}

		virtual inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return position + this->template getStartSieveNumber();
		}
};

