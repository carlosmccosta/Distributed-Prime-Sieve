#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, size_t dynamicSchedulingSegmentSizeInElements = 1048576, size_t dynamicSchedulingNumberSegments = 0, string outputResultsFilename = "") :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, numberOfThreads, sendResultsToRoot, countNumberOfPrimesOnNode,
						sendPrimesCountToRoot, dynamicSchedulingSegmentSizeInElements, dynamicSchedulingNumberSegments, outputResultsFilename) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingTimeAndCacheWithWheel() {
		}


		inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			return blockSize;
		}

		inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber());
		}

		inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return position + this->template getStartSieveNumber();
		}
};
