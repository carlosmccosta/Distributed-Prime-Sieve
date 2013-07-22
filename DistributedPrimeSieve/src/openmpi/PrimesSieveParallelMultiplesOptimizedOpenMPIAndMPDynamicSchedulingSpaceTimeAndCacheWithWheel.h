#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0,
				bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, size_t dynamicSchedulingSegmentSizeInElements = 1048576,
				size_t dynamicSchedulingNumberSegments = 0, string outputResultsFilename = "", bool outputOnlyLastSegment = false) :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, numberOfThreads, sendResultsToRoot,
						countNumberOfPrimesOnNode, sendPrimesCountToRoot, dynamicSchedulingSegmentSizeInElements, dynamicSchedulingNumberSegments, outputResultsFilename, outputOnlyLastSegment) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel() {
		}

		virtual inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			if (blockSize % 2 == 0) {
				return (blockSize >> 1);
			} else {
				return ((blockSize >> 1) + 1);
			}
		}

		virtual inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		virtual inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return (position << 1) + this->template getStartSieveNumber();
		}
};

