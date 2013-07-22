#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPTimeAndCacheWithWheel() {
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

