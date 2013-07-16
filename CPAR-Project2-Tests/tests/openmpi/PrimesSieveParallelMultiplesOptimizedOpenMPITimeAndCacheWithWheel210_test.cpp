#include <boost/test/unit_test.hpp>
#include "openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithModulo210Wheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(546, 4);
		primesSieveParallelMultiples.computePrimes(546);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes100.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(7926, 1);   // 8 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks1byte.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(7926, 4);   // 32 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks4bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(7926, 16);   // 128 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks16bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks128byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(1299828, 128);   // 1024 elements block
		primesSieveParallelMultiples.computePrimes(1299828);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes100000_Blocks128bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks512byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte> primesSieveParallelMultiples(1299828, 512);   // 4096 elements block
		primesSieveParallelMultiples.computePrimes(1299828);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/openMPIOptimizedTimeAndCacheWithModulo210WheelPrimes100000_Blocks512bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
