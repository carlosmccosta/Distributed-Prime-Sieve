#include <boost/test/unit_test.hpp>
#include "openmp/PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel.h"
#include "WheelFactorization.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithModulo210Wheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples;
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(1);
		primesSieveParallelMultiples.computePrimes(546);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/100.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples(1);   // 8 elements block
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks1byte.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(4);
		primesSieveParallelMultiples.computePrimes(7926);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/1000.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks1byte.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples(4);   // 32 elements block
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks4bytes.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(4);
		primesSieveParallelMultiples.computePrimes(7926);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/1000.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks4bytes.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples(16);   // 128 elements block
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks16bytes.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(4);
		primesSieveParallelMultiples.computePrimes(7926);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/1000.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks16bytes.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks128byte) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples(128);   // 1024 elements block
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks128bytes.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(10);
		primesSieveParallelMultiples.computePrimes(1299828);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/100000.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks128bytes.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks512byte) {
		PrimesSieveParallelMultiplesSegmentedOptimizedOpenMPSpaceTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveParallelMultiples(512);   // 4096 elements block
		primesSieveParallelMultiples.setOutputResultsFilename("./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks512bytes.txt");
		primesSieveParallelMultiples.setSegmentSizeInBlocks(20);
		primesSieveParallelMultiples.computePrimes(1299828);

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFiles("./tests/testfiles/100000.txt", "./tests/testresults/openMPSegmentedOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks512bytes.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
