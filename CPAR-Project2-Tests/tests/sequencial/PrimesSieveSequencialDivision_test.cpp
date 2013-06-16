#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialDivision.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialDivision_test)
	BOOST_AUTO_TEST_CASE(compute100Primes) {
		PrimesSieveSequencialDivision<vector<bool> > primesSieveSequencialDivision;
		primesSieveSequencialDivision.computePrimes(541);
		primesSieveSequencialDivision.savePrimesToFile("tests/testresults/devisionPrimes100.txt");
		
		BOOST_CHECK(primesSieveSequencialDivision.checkPrimesFromFile("./tests/testfiles/100.txt"));
		
	}
	
	BOOST_AUTO_TEST_CASE(compute1000Primes) {
		PrimesSieveSequencialDivision<vector<bool> > primesSieveSequencialDivision;
		primesSieveSequencialDivision.computePrimes(7919);
		primesSieveSequencialDivision.savePrimesToFile("tests/testresults/devisionPrimes1000.txt");
		
		BOOST_CHECK(primesSieveSequencialDivision.checkPrimesFromFile("./tests/testfiles/1000.txt"));
		
	}
	BOOST_AUTO_TEST_SUITE_END()
