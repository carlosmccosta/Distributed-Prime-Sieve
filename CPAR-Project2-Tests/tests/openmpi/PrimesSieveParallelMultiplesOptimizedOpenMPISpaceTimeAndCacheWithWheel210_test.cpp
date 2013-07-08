#include <boost/test/unit_test.hpp>
#include "openmpi/PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithModulo210Wheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples;
		primesSieveParallelMultiples.computePrimes(546);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples(1);   // 8 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks1byte.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples(4);   // 32 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks4bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples(16);   // 128 elements block
		primesSieveParallelMultiples.computePrimes(7926);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes1000_Blocks16bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks128byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples(128);   // 1024 elements block
		primesSieveParallelMultiples.computePrimes(1299828);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks128bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks512byte) {
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<vector<unsigned char>, Modulo210WheelByte> primesSieveParallelMultiples(512);   // 4096 elements block
		primesSieveParallelMultiples.computePrimes(1299828);
		primesSieveParallelMultiples.savePrimesToFile("tests/testresults/OpenMPIOptimizedSpaceTimeAndCacheWithModulo210WheelPrimes100000_Blocks512bytes.txt");

		BOOST_CHECK(primesSieveParallelMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
