#pragma once

#include "../PrimesSieve.h"
#include "../lib/PrimesUtils.h"
#include "../lib/Settings.h"

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
		bool _countNumberOfPrimes;
		bool _outputOnlyLastSegment;
		string _outputResultsFilename;

		size_t _globalMaxRange;
		size_t _segmentStartNumber;
		size_t _segmentEndNumber;

		size_t _parcialPrimesCount;
		vector<size_t> _sievingPrimes;
		WheelType _wheelSieve;

	public:
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel(size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0, size_t segmentSizeInBlocks = 128 * 1024,
				bool countNumberOfPrimes = false, bool outputOnlyLastSegment = false, string outputResultsFilename = "") :
				_blockSizeInElements(blockSizeInBytes * 8 * 2), _numberOfThreads(numberOfThreads), _segmentSizeInBlocks(segmentSizeInBlocks), _countNumberOfPrimes(countNumberOfPrimes), _outputOnlyLastSegment(
						outputOnlyLastSegment), _outputResultsFilename(outputResultsFilename), _globalMaxRange(11), _segmentStartNumber(11), _segmentEndNumber(12), _parcialPrimesCount(0) {
		}

		virtual ~PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel() {
		}

		virtual void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			_globalMaxRange = maxRange;
			this->template clearPrimesValues();

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			size_t blockBeginNumber = this->template getBlockBeginNumber();

			this->template setStartSieveNumber(blockBeginNumber);
			this->template setMaxRange(maxRangeSquareRoot);

			int numberThreadsToUse = omp_get_max_threads();
			size_t numberOfThreads = this->template getNumberOfThreads();
			if (numberOfThreads != 0) {
				numberThreadsToUse = (int) numberOfThreads;
			}
			omp_set_num_threads(numberThreadsToUse);

			// compute sieving primes
			vector<pair<size_t, size_t> > sievingMultiples;
			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, maxRange, sievingMultiples);
			computeSievingPrimesTimer.stop();
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;

			this->template cleanOutputFile();

			if (_outputResultsFilename != "") {
				cout << "\n    > Exporting results to file " << _outputResultsFilename << endl;;
			}
			this->template outputResults();

			// sieve remaining numbers
			this->template removeComposites(maxRangeSquareRoot, maxRange, sievingMultiples);

			if (_countNumberOfPrimes) {
				this->template setPrimesCount(_parcialPrimesCount);
				cout << "    --> Counted " << _parcialPrimesCount << " primes in [" << blockBeginNumber << ", " << maxRange << "]" << endl;
			}

			this->template getPerformanceTimer().stop();
		}

		void computeSievingPrimes(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t blockBeginNumber = getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			_segmentStartNumber = blockBeginNumber;
			_segmentEndNumber = maxRangeSquareRoot + 1;

			this->template initPrimesBitSetSizeForSievingPrimes(maxRangeSquareRoot);
			size_t blockEndNumber = min(blockBeginNumber + _blockSizeInElements, maxRangeSquareRoot + 1);
			size_t numberBlocks = (size_t) ceil((double) ((maxRangeSquareRoot + 1) - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockBeginNumber = blockEndNumber;
				blockEndNumber += _blockSizeInElements;

				if (blockEndNumber >= maxRangeSquareRoot) {
					blockEndNumber = maxRangeSquareRoot + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocksParallel(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}

			if (_countNumberOfPrimes) {
				_parcialPrimesCount += this->template getNumberPrimesFound();
			}
		}

		void removeComposites(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			++maxRangeSquareRoot; // maxRangeSquareRoot was already sieved

			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t numberBlocks = (size_t) ceil((double) ((maxRange + 1) - maxRangeSquareRoot) / (double) blockSizeInElements);
			const size_t numberSegments = (size_t) ceil((double) numberBlocks / (double) _segmentSizeInBlocks);

			int numberThreadsToUse = omp_get_max_threads();
			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					numberThreadsToUse = (int) numberBlocks;
				} else {
					numberThreadsToUse = (int) _numberOfThreads;
				}
			}

			omp_set_num_threads(numberThreadsToUse);

			size_t startBlock = 0;
			size_t endBlock = min(_segmentSizeInBlocks, numberBlocks);
			size_t numberElementsPerSegment = _blockSizeInElements * _segmentSizeInBlocks;

			for (size_t segmentNumber = 0; segmentNumber < numberSegments; ++segmentNumber) {
				size_t priviousBlockNumber = -1;
				vector<pair<size_t, size_t> > sievingMultiples;

				size_t numberBlocksLeft = endBlock - startBlock;

				_segmentStartNumber = startBlock * blockSizeInElements + maxRangeSquareRoot;
				_segmentEndNumber = _segmentStartNumber + numberElementsPerSegment;

				if (_segmentStartNumber % 2 == 0) {
					++_segmentStartNumber;
				}

				if (_segmentEndNumber > maxRange) {
					_segmentEndNumber = maxRange + 1;
				}

				this->template setStartSieveNumber(_segmentStartNumber);
				this->template setMaxRange(_segmentEndNumber - 1);

				this->template initPrimesBitSetSizeForSieving(_segmentEndNumber - _segmentStartNumber);

#				pragma omp parallel for \
					if (numberBlocksLeft > 8) \
					default(shared) \
					shared(sievingMultiplesFirstBlock) \
					firstprivate(maxRange, maxRangeSquareRoot, blockSizeInElements, priviousBlockNumber, sievingMultiples) \
					schedule(guided, 64)
					for (size_t blockNumber = startBlock; blockNumber < endBlock; ++blockNumber) {
						size_t blockBeginNumber = blockNumber * blockSizeInElements + maxRangeSquareRoot;
						size_t blockEndNumber = blockBeginNumber + blockSizeInElements;

						if (blockEndNumber > maxRange) {
							blockEndNumber = maxRange + 1;
						}

						if (blockNumber == 0) {
							sievingMultiples = sievingMultiplesFirstBlock;
						} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
							this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
						}
						priviousBlockNumber = blockNumber;
						this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
					}

					if (_countNumberOfPrimes) {
						_parcialPrimesCount += this->template getNumberPrimesFound();
					}

					this->template outputResults();

					startBlock = endBlock;
					endBlock += _segmentSizeInBlocks;
					if (endBlock > numberBlocks) {
						endBlock = numberBlocks;
					}
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
				size_t sievingMultiplesSize = sievingMultiples.size();

#				pragma omp parallel for \
					if (sievingMultiplesSize > 4) \
					default(shared) \
					schedule(guided, 4)
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

				virtual inline size_t getNumberBitsToStore(size_t maxRange) {
					return this->template getNumberBitsToStoreBlock((maxRange + 1) - _wheelSieve.getFirstPrimeToSieve());
				}

				inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
					if (blockSize % 2 == 0) {
						return (blockSize >> 1);
					} else {
						return ((blockSize >> 1) + 1);
					}
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
					if (_segmentStartNumber == _wheelSieve.getFirstPrimeToSieve()) {
						primesValues.push_back(2);
						primesValues.push_back(3);
						primesValues.push_back(5);

						if (_wheelSieve.getNumberPrimesSievedByTheWheel() == 4) {
							primesValues.push_back(7);
						}
					}

					size_t possiblePrime = _segmentStartNumber;
					if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
						possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
					}

					size_t maxRange = _segmentEndNumber - 1;
					for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
						if (!this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(possiblePrime, _segmentStartNumber)) {
							primesValues.push_back(possiblePrime);
						}
					}

					return primesValues;
				}

				void savePrimes(ostream& outputStream) {
					vector<size_t>& primesValues = this->template getPrimesValues();

					if (primesValues.size() <= 2) {
						if (_segmentStartNumber == _wheelSieve.getFirstPrimeToSieve()) {
							outputStream << 2 << endl;
							outputStream << 3 << endl;
							outputStream << 5 << endl;

							if (_wheelSieve.getNumberPrimesSievedByTheWheel() == 4) {
								outputStream << 7 << endl;
							}
						}

						size_t maxRange = _segmentEndNumber - 1;
						size_t possiblePrime = _segmentStartNumber;
						if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
							possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
						}

						for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
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

				bool cleanOutputFile() {
					if (_outputResultsFilename != "") {
						ofstream outputStream(_outputResultsFilename.c_str());
						if (outputStream.is_open()) {
							outputStream.close();
							return true;
						}
					}

					return false;
				}

				size_t getNumberPrimesFound() {
					vector<size_t>& primesValues = this->template getPrimesValues();
					if (primesValues.size() >= 2) {
						return primesValues.size();
					}

					size_t primesFound = this->template getPrimesCount();
					if (this->template getPrimesCount() != 0) {
						return primesFound;
					}

					if (_segmentStartNumber == _wheelSieve.getFirstPrimeToSieve()) {
						primesFound = _wheelSieve.getNumberPrimesSievedByTheWheel();
					} else {
						primesFound = 0;
					}

					size_t maxRange = _segmentEndNumber - 1;
//					int maxNumberThreads = omp_get_max_threads();
//					size_t minNumberPrimesPerThread = 100;
//					int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));

					int numberThreadsToUse = omp_get_max_threads();
					if (_numberOfThreads != 0) {
						numberThreadsToUse = (int) _numberOfThreads;
					}

					size_t segmentStartNumber = _segmentStartNumber;
					WheelType wheelSieve = _wheelSieve;

#					ifdef DEBUG_OUTPUT
					cout << "    --> Counting in [" << _segmentStartNumber << ", " << maxRange << "] using " << numberThreadsToUse << " threads" << endl;
#					endif

					size_t numberBlocks = (size_t) ceil((double) ((maxRange + 1) - _segmentStartNumber) / (double) _blockSizeInElements);
					size_t blockSizeInElements = _blockSizeInElements;

#					pragma omp parallel for \
						default(shared) \
						firstprivate(maxRange, blockSizeInElements, numberBlocks, segmentStartNumber, wheelSieve) \
						schedule(guided, 64) \
						reduction(+: primesFound) \
						num_threads(numberThreadsToUse)
					for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
						size_t possiblePrime = blockNumber * blockSizeInElements + segmentStartNumber;
						size_t nextPossiblePrimeNumberEndBlock = min(possiblePrime + blockSizeInElements, maxRange + 1);

						if (!wheelSieve.isNumberPossiblePrime(possiblePrime)) {
							possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
						}

						while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
							if (!(this->PrimesSieve<FlagsContainer>::template getPrimesBitsetValueBlock(possiblePrime, segmentStartNumber))) {
								++primesFound;
							}
							possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
						}
					}

#					ifdef DEBUG_OUTPUT
					cout << "    --> Found " << primesFound << " primes in [" << _segmentStartNumber << ", " << maxRange << "]" << endl;
#					endif

					return primesFound;
				}

				bool outputResults() {
					if (_outputOnlyLastSegment && (_segmentEndNumber - 1)!= _globalMaxRange) {
						return false;
					}

					if (_outputResultsFilename == "stdout") {
						cout << "\n\n=============================================  Computed primes  =============================================\n\n";
						this->template printPrimesToConsole();
						cout << "\n" << endl;
						return true;
					} else if (_outputResultsFilename != "") {
						ofstream outputStream(_outputResultsFilename.c_str(), ofstream::out | ofstream::app);
						if (outputStream.is_open()) {
							this-> template savePrimes(outputStream);
							outputStream.close();
							return true;
						}

						else {
							cerr << "    !!!!! Export to file " << _outputResultsFilename << "failed !!!!!" << endl;
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

				inline bool isOutputOnlyLastSegment() const {
					return _outputOnlyLastSegment;
				}

				inline void setOutputOnlyLastSegment(bool outputOnlyLastSegment) {
					_outputOnlyLastSegment = outputOnlyLastSegment;
				}

				inline const string& getOutputResultsFilename() const {
					return _outputResultsFilename;
				}

				inline void setOutputResultsFilename(const string& outputResultsFilename) {
					_outputResultsFilename = outputResultsFilename;
				}

				inline size_t getSegmentSizeInBlocks() const {
					return _segmentSizeInBlocks;
				}

				inline void setSegmentSizeInBlocks(size_t segmentSizeInBlocks) {
					_segmentSizeInBlocks = segmentSizeInBlocks;
				}
			};
