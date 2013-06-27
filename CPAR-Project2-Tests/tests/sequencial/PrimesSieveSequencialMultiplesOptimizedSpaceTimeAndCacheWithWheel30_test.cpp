#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(541);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks10) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples(10 / 8);   // 10 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes1000_Blocks10.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks32) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples(32 / 8);   // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes1000_Blocks32.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks100) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples(100 / 8);   // 100 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes1000_Blocks100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1000) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples(1000 / 8);   // 1000 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes1000_Blocks1000.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks2000) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCacheWithWheel<vector<bool>, Modulo30Wheel> primesSieveSequencialMultiples(2000 / 8);   // 2000 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCacheWithWheelPrimes1000_Blocks2000.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
