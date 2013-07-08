#pragma once

#include "../PrimesSieve.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <mpi.h>

using std::sqrt;
using std::min;
using std::pair;

template<typename FlagsContainer>
class PrimesSieveParallelMultiplesOptimizedOpenMPI: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;
		int _processID;
		int _numberProcesses;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPI(size_t blockSizeInElements = 128 * 1024) :
				_blockSizeInElements(blockSizeInElements) {
			_processID = 0;
			_numberProcesses = 1;
//			MPI_Comm_rank(MPI_COMM_WORLD, &_processID);
//			MPI_Comm_size(MPI_COMM_WORLD, &_numberProcesses);
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPI() {
		}

		void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template clearPrimesValues();
			if (_processID == 0) {
				this->template initPrimesBitSetSize(maxRange);
			} else {
				this->template initPrimesBitSetSize(this->template getProcessBlockSize(_processID, _numberProcesses, maxRange));
			}

			// adjustment of maxRange in order to calculate the primes <= maxRange instead of < maxRange (since algorithm computes [begin, maxRange[
			++maxRange;

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			vector<pair<size_t, size_t> > sievingMultiples;
			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);

			size_t processStartBlockNumber = this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange);
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);

			if (_processID == _numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}
			this->template removeComposites(processStartBlockNumber, processEndBlockNumber, sievingMultiples);

			if (_processID == 0) {
				this->template collectResultsFromProcessGroup(maxRange);
			} else {
				this->template sendResultsToRootProcess(maxRange);
			}

			this->template getPerformanceTimer().stop();
		}

		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumber(maxRangeSquareRoot);
			size_t blockBeginNumber = getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t blockIndexBegin = this->template getBitsetPositionToNumber(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + _blockSizeInElements, maxIndexRangeSquareRoot);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;
				if (blockIndexEnd > maxIndexRangeSquareRoot) {
					blockIndexEnd = maxIndexRangeSquareRoot;
				}

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t processEndBlockNumberIndex = this->template getNumberBitsToStore(processEndBlockNumber) - 1;
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumber(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + processBeginBlockNumberIndex;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				if (blockIndexEnd > processEndBlockNumberIndex) {
					blockIndexEnd = processEndBlockNumberIndex;
				}

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				if (blockNumber == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		void collectResultsFromProcessGroup(size_t maxRange) {
//			FlagsContainer& primesBitset = this->template getPrimesBitset();
//			MPI_Status status;
//			for (int processID = 1; processID < _numberProcesses; ++processID) {
//				size_t processStartBlockNumber = this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange);
//				size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);
//
//				if (processID == _numberProcesses - 1) {
//					processEndBlockNumber = maxRange + 1;
//				}
//				size_t blockSize = processEndBlockNumber - processStartBlockNumber;
//				size_t positionToStoreResults = this->template getBitsetPositionToNumber(processStartBlockNumber);
//				MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, processID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//			}
		}

		void sendResultsToRootProcess(size_t maxRange) {
//			FlagsContainer& primesBitset = this->template getPrimesBitset();
//			size_t processStartBlockNumber = this->template getProcessStartBlockNumber(_processID, _numberProcesses, maxRange);
//			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(_processID, _numberProcesses, maxRange);
//
//			if (_processID == _numberProcesses - 1) {
//				processEndBlockNumber = maxRange + 1;
//			}
//			size_t blockSize = processEndBlockNumber - processStartBlockNumber;
//
//			MPI_Send(&primesBitset[0], blockSize, MPI_UNSIGNED_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
		}

		virtual void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void initPrimesBitSetSize(size_t maxRange) = 0;


		virtual inline size_t getBitsetPositionToNumber(size_t number) {
			return number;
		}

		virtual inline size_t getNumberAssociatedWithBitsetPosition(size_t position) {
			return position;
		}

		inline size_t getProcessStartBlockNumber(size_t processID, size_t numberProcesses, size_t maxRange) {
			return processID * maxRange / numberProcesses;
		}

		inline size_t getProcessEndBlockNumber(size_t processID, size_t numberProcesses, size_t maxRange) {
			return ((processID + 1) * maxRange / numberProcesses) - 1;
		}

		inline size_t getProcessBlockSize(size_t processID, size_t numberProcesses, size_t maxRange) {
			return (this->template getProcessEndBlockNumber(processID, numberProcesses, maxRange) - this->template getProcessStartBlockNumber(processID, numberProcesses, maxRange));
		}

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
		}

		virtual size_t getBlockBeginNumber() {
			return this->template getNumberAssociatedWithBitsetPosition(0);
		}

		inline size_t getBlockSizeInBytes() {
			return _blockSizeInElements / 8;
		}

		inline void setBlockSizeInBytes(size_t blockSize) {
			_blockSizeInElements = blockSize * 8;
		}

		inline size_t getBlockSizeInElements() const {
			return _blockSizeInElements;
		}

		inline void setBlockSizeInElements(size_t blockSizeInElements) {
			_blockSizeInElements = blockSizeInElements;
		}

		inline int getNumberProcesses() const {
			return _numberProcesses;
		}

		inline int getProcessId() const {
			return _processID;
		}
	};

