#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithModulo210Wheel_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(546);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples(1); // 8 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks1byte.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples(4); // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks4bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples(16); // 128 elements block
		primesSieveSequencialMultiples.computePrimes(7926);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes1000_Blocks16bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks128byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples(128); // 1024 elements block
		primesSieveSequencialMultiples.computePrimes(1299828);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes100000_Blocks128bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute100000PrimesBlocks512byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCacheWithWheel<PrimesFlagsContainer, Modulo210Wheel> primesSieveSequencialMultiples(512); // 4096 elements block
		primesSieveSequencialMultiples.computePrimes(1299828);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCacheWithModulo210WheelPrimes100000_Blocks512bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
