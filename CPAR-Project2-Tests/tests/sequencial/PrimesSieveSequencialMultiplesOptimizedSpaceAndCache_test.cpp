#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceAndCache.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedSpaceAndCache_test)
BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> > primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(541);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCachePrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks10) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> > primesSieveSequencialMultiples(64*10/8); // 10 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCachePrimes1000_Blocks10.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks32) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> > primesSieveSequencialMultiples(64*32/8); // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCachePrimes1000_Blocks10.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks100) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> > primesSieveSequencialMultiples(64*100/8); // 100 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCachePrimes1000_Blocks100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1000) {
		PrimesSieveSequencialMultiplesOptimizedSpaceAndCache<vector<bool> > primesSieveSequencialMultiples(64*1000/8); // 1000 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedSpaceAndCachePrimes1000_Blocks1000.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
