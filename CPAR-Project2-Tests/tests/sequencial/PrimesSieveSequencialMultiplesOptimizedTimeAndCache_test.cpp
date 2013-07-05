#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiplesOptimizedTimeAndCache.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiplesOptimizedTimeAndCache_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(541);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes100.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks1byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples(1); // 8 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes1000_Blocks1byte.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks4byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples(4); // 32 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes1000_Blocks4bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks16byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples(16); // 128 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes1000_Blocks16bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks128byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples(128); // 1024 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes1000_Blocks128bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}

	BOOST_AUTO_TEST_CASE(compute1000PrimesBlocks512byte) {
		PrimesSieveSequencialMultiplesOptimizedTimeAndCache<vector<char> > primesSieveSequencialMultiples(512); // 4096 elements block
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/optimizedTimeAndCachePrimes1000_Blocks512bytes.txt");

		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
	}
	BOOST_AUTO_TEST_SUITE_END()
