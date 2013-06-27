#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithModulo210Wheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(546);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks10) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples(10 / 8);   // 10 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes1000_Blocks10.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks32) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples(32 / 8);   // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes1000_Blocks32.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks100) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples(100 / 8);   // 100 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes1000_Blocks100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1000) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples(1000 / 8);   // 1000 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes1000_Blocks1000.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks2000) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCacheWithWheel<vector<bool>, Modulo210Wheel> primesSieveSequencialMultiples(2000 / 8);   // 2000 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCacheWithModulo210WheelPrimes1000_Blocks2000.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
