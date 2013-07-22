#pragma once

#include "../PrimesSieve.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

#include <omp.h>

using std::vector;
using std::min;
using std::max;
using std::pair;
using std::endl;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;
		size_t _numberOfThreads;
		size_t _segmentSizeInBlocks;
		bool _outputOnlyLastSegment;
		string _outputResultsFilename;

		size_t _globalMaxRange;
		size_t _segmentStartNumber;
		size_t _segmentEndNumber;

		vector<size_t> _sievingPrimes;
		WheelType _wheelSieve;

	public:
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel(size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0, size_t segmentSizeInBlocks = 128 * 1024,
				bool outputOnlyLastSegment = false, string outputResultsFilename = "") :
				_blockSizeInElements(blockSizeInBytes * 8 * 2), _numberOfThreads(numberOfThreads), _segmentSizeInBlocks(segmentSizeInBlocks), _outputOnlyLastSegment(outputOnlyLastSegment), _outputResultsFilename(
						outputResultsFilename), _globalMaxRange(11), _segmentStartNumber(11), _segmentEndNumber(12) {
		}

		virtual ~PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel();

		virtual void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			_globalMaxRange = maxRange;
			this->template clearPrimesValues();

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			this->template setStartSieveNumber(this->template getBlockBeginNumber());
			this->template setMaxRange(maxRangeSquareRoot);

			// compute sieving primes
			vector<pair<size_t, size_t> > sievingMultiples;
			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, maxRange, sievingMultiples);
			computeSievingPrimesTimer.stop();
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;
			this->template outputResults();

			// sieve remaining numbers
			this->template removeComposites(maxRangeSquareRoot, maxRange, sievingMultiples);

			this->template getPerformanceTimer().stop();
		}

		void computeSievingPrimes(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t blockBeginNumber = getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			size_t blockEndNumber = min(blockBeginNumber + _blockSizeInElements, maxRangeSquareRoot + 1);
			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockBeginNumber = blockEndNumber;
				blockEndNumber += _blockSizeInElements;

				if (blockEndNumber >= maxRange) {
					blockEndNumber = maxRange + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocksParallel(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		void removeComposites(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t numberBlocks = ceil((double) ((maxRange + 1) - maxRangeSquareRoot) / (double) blockSizeInElements);
			const size_t numberSegments = numberBlocks / _segmentSizeInBlocks;

			size_t numberThreadsToUse = omp_get_max_threads();
			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					numberThreadsToUse = numberBlocks;
				} else {
					numberThreadsToUse = _numberOfThreads;
				}
			}

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

			for (size_t segmentNumber = 0; segmentNumber < numberSegments; ++segmentNumber) {
				size_t startBlockNumber = segmentNumber * _segmentSizeInBlocks;
				size_t endBlockNumber = startBlockNumber + _segmentSizeInBlocks;

				size_t beginNumber = startBlockNumber * blockSizeInElements + maxRangeSquareRoot;
				size_t endNumber = beginNumber + blockSizeInElements;

				if (endNumber > maxRange) {
					endNumber = maxRange + 1;
				}

				this->template setStartSieveNumber(beginNumber);
				this->template setMaxRange(endNumber);

				this->template initPrimesBitSetSizeForSieving(endNumber - beginNumber);

#				pragma omp parallel for \
					if (numberBlocks > 32) \
					default(shared) \
					num_threads(numberThreadsToUse) \
					shared(sievingMultiplesFirstBlock) \
					firstprivate(maxRange, maxRangeSquareRoot, blockSizeInElements, priviousBlockNumber, sievingMultiples) \
					schedule(guided, 64)
					for (size_t blockNumber = startBlockNumber; blockNumber < endBlockNumber; ++blockNumber) {
						size_t blockBeginNumber = blockNumber * blockSizeInElements + maxRangeSquareRoot;
						size_t blockEndNumber = blockBeginNumber + blockSizeInElements;

						if (blockEndNumber > maxRange) {
							blockEndNumber = maxRange + 1;
						}

						if (blockNumber == 0) {
							sievingMultiples = sievingMultiplesFirstBlock;
						} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
							computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
						}
						priviousBlockNumber = blockNumber;
						this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
					}

					this->template outputResults();
				}

			}

			void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
				sievingMultiples.clear();
				size_t sievingPrimesSize = _sievingPrimes.size();
				for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
					size_t primeNumber = _sievingPrimes[sievingPrimesIndex];
					size_t primeMultiple = PrimesUtils::closestPrimeMultiple(primeNumber, blockBeginNumber);
					size_t primeMultipleIncrement = primeNumber << 1;

					if (primeMultiple < blockBeginNumber || primeMultiple == primeNumber) {
						primeMultiple += primeNumber;
					}

					if (primeMultiple % 2 == 0) {
						primeMultiple += primeNumber;
					}

					sievingMultiples.push_back(pair<size_t, size_t>(primeMultiple, primeMultipleIncrement));
				}
				//					cout << "init sievingMultiples in block [" << blockBeginNumber << ", " << blockEndNumber << "]" << endl;
			}

			void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
				size_t sievingMultiplesSize = sievingMultiples.size();
				for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
					pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
					size_t primeMultiple = primeCompositeInfo.first;
					size_t primeMultipleIncrement = primeCompositeInfo.second;

					for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
						this->PrimesSieve<FlagsContainer>::template setPrimesBitsetValueBlock(primeMultiple, _segmentStartNumber, true);
					}

					sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
				}
			}

			void removeMultiplesOfPrimesFromPreviousBlocksParallel(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
				size_t numberThreadsToUse = omp_get_max_threads();
				size_t numberOfThreads = this->template getNumberOfThreads();

				if (numberOfThreads != 0) {
					numberThreadsToUse = numberOfThreads;
				}

				size_t sievingMultiplesSize = sievingMultiples.size();

#				pragma omp parallel for \
					if (sievingMultiplesSize > 8) \
					default(shared) \
					schedule(guided, 8) \
					num_threads(numberThreadsToUse)
				for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
					pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
					size_t primeMultiple = primeCompositeInfo.first;
					size_t primeMultipleIncrement = primeCompositeInfo.second;

					for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
						this->PrimesSieve<FlagsContainer>::template setPrimesBitsetValueBlock(primeMultiple, _segmentStartNumber, true);
					}

					sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
				}
			}

				void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
					size_t maxPrimeNumberSearch = blockEndNumber;
					if (maxPrimeNumberSearch >= maxRangeSquareRoot) {
						maxPrimeNumberSearch = maxRangeSquareRoot + 1;
					}

					size_t primeNumber = blockBeginNumber;
					if (_wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
						primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber);
					}

					vector<size_t>& _sievingPrimes = this->template getSievingPrimes();

					for (; primeNumber < maxPrimeNumberSearch; primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber)) {
						// for each number not marked as composite (prime number)
						if (!this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(primeNumber, _segmentStartNumber)) {
							_sievingPrimes.push_back(primeNumber);

							//use it to calculate his composites
							size_t primeMultipleIncrement = primeNumber << 1;
							size_t compositeNumber = primeNumber * primeNumber;
							for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
								this->PrimesSieve<FlagsContainer>::template setPrimesBitsetValueBlock(compositeNumber, _segmentStartNumber, true);
							}
							sievingMultiples.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
						}
					}
				}

				void initPrimesBitSetSize(size_t maxRange) {
					this->template initPrimesBitset(maxRange);

					vector<size_t>& _sievingPrimes = this->template getSievingPrimes();
					_sievingPrimes.clear();

					size_t numberSievingPrimes = this->template getNumberOfPrimesInRange((size_t) sqrt(maxRange));
					_sievingPrimes.reserve(numberSievingPrimes);
				}

				void initPrimesBitSetSizeForSievingPrimes(size_t maxRangeSquareRoot) {
					this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStore(maxRangeSquareRoot));

					size_t numberSievingPrimes = this->template getNumberOfPrimesInRange(maxRangeSquareRoot);

					_sievingPrimes.clear();
					_sievingPrimes.reserve(numberSievingPrimes);
				}

				void initPrimesBitSetSizeForSieving(size_t blockSize) {
					this->PrimesSieve<FlagsContainer>::template initPrimesBitSetSize(this->template getNumberBitsToStoreBlock(blockSize));
				}

				inline size_t getBlockBeginNumber() {
					return _wheelSieve.getFirstPrimeToSieve();
				}

				inline size_t getNumberBitsToStore(size_t maxRange) {
					return this->PrimesSieve<FlagsContainer>::template getNumberBitsToStore(maxRange);
				}

				inline size_t getBitsetPositionToNumberOpenMP(size_t number, size_t segmentStartNumber) {
					return (number - segmentStartNumber) >> 1;
				}

				inline size_t getNumberAssociatedWithBitsetPositionOpenMP(size_t position, size_t segmentStartNumber) {
					return (position << 1) + segmentStartNumber;
				}

				vector<size_t>& extractPrimesFromBitset() {
					vector<size_t>& primesValues = this->template getPrimesValues();

					primesValues.clear();
					primesValues.push_back(2);
					primesValues.push_back(3);
					primesValues.push_back(5);
					primesValues.push_back(7);

					size_t maxRange = this->template getMaxRange();
					for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
						if (!this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(possiblePrime, _segmentStartNumber)) {
							primesValues.push_back(possiblePrime);
						}
					}

					return primesValues;
				}

				void savePrimes(ostream& outputStream) {
					vector<size_t>& primesValues = this->template getPrimesValues();

					if (primesValues.size() <= 2) {
						outputStream << 2 << endl;
						outputStream << 3 << endl;
						outputStream << 5 << endl;
						outputStream << 7 << endl;

						size_t maxRange = this->template getMaxRange();
						for (size_t possiblePrime = 11; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
							if (!this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(possiblePrime, _segmentStartNumber)) {
								outputStream << possiblePrime << endl;
							}
						}
					} else {
						size_t iSize = primesValues.size();
						for (size_t i = 0; i < iSize; ++i) {
							outputStream << primesValues[i] << endl;
						}
					}
				}

				size_t getNumberPrimesFound() {
					vector<size_t>& primesValues = this->template getPrimesValues();

					if (primesValues.size() >= 2)
					return primesValues.size();

					size_t primesFound = 4;
					size_t maxRange = this->template getMaxRange();
					int maxNumberThreads = omp_get_max_threads();
					size_t minNumberPrimesPerThread = 100;
					int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));
					size_t numberPrimesToCheckInBlock = maxRange / numberThreads;

#					pragma omp parallel for \
						default(shared) \
						firstprivate(maxRange, numberThreads, numberPrimesToCheckInBlock) \
						reduction(+: primesFound) \
						num_threads(numberThreads)
					for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
						size_t possiblePrime;

						possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + 11;
						if (!_wheelSieve.isNumberPossiblePrime(possiblePrime)) {
							possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
						}

						size_t nextPossiblePrimeNumberEndBlock = min((threadBlockNumber + 1) * numberPrimesToCheckInBlock + 11, maxRange + 1);

						while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
							if (!(this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(possiblePrime, _segmentStartNumber))) {
								++primesFound;
							}
							possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
						}
					}

					return primesFound;
				}

				bool outputResults() {
					if (_outputOnlyLastSegment) {
						size_t maxRangeSegment = this->template getMaxRange();
						if (maxRangeSegment != _globalMaxRange) {
							return false;
						}
					}

					if (_outputResultsFilename == "stdout") {
						cout << "\n\n=============================================  Computed primes  =============================================\n\n";
						this->template printPrimesToConsole();
						cout << "\n" << endl;
						return true;
					} else if (_outputResultsFilename != "") {
						ofstream outputStream(_outputResultsFilename.c_str(), ofstream::out | ofstream::app);
						if (outputStream.is_open()) {
							if (this-> template savePrimes(outputStream)) {
								return true;
							} else {
								cerr << "    !!!!! Export to file " << _outputResultsFilename << "failed !!!!!" << endl;
							}
						}
					}

					return false;
				}

				inline WheelType& getWheelSieve() {
					return _wheelSieve;
				}

				inline void clearPrimesValues() {
					this->template getPrimesValues().clear();
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

				inline size_t getNumberOfThreads() const {
					return _numberOfThreads;
				}

				inline void setNumberOfThreads(size_t numberOfThreads) {
					_numberOfThreads = numberOfThreads;
				}

				inline vector<size_t>& getSievingPrimes() {
					return _sievingPrimes;
				}

				inline void setSievingPrimes(const vector<size_t>& sievingPrimes) {
					_sievingPrimes = sievingPrimes;
				}
			};
