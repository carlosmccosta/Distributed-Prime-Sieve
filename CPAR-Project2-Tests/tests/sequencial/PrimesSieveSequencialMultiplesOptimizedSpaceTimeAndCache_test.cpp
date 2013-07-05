#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache_test)
BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(541);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples(1); // 8 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes1000_Blocks1byte.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples(4); // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes1000_Blocks4bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples(16); // 128 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes1000_Blocks16bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks128byte) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples(128); // 1024 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes1000_Blocks128bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks512byte) {
		PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache<vector<bool> > primesSieveSequencialMultiples(512); // 4096 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceTimeAndCachePrimes1000_Blocks512bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
