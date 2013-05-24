#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialMultiples.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialMultiples_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialMultiples<vector<bool> > primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(541);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/v2primes100.txt");
		
		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/100.txt"));
		
	}
	
	BOOST_AUTO_TEST_CASE(compute1000Primes) {
		PrimesSieveSequencialMultiples<vector<bool> > primesSieveSequencialMultiples;
		primesSieveSequencialMultiples.computePrimes(7919);
		primesSieveSequencialMultiples.savePrimesToFile("tests/testresults/v2primes1000.txt");
		
		BOOST_CHECK(primesSieveSequencialMultiples.checkPrimesFromFile("./tests/testfiles/1000.txt"));
		
	}
	BOOST_AUTO_TEST_SUITE_END()
