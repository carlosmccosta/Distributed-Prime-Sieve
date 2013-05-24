#include <boost/test/unit_test.hpp>
#include "sequencial/PrimesSieveSequencialDivision.h"

BOOST_AUTO_TEST_SUITE(PrimesSieveSequencialDivision_test)

	BOOST_AUTO_TEST_CASE(computePrimes) {
		PrimesSieveSequencialDivision primesSieveSequencialDivision;
		primesSieveSequencialDivision.computePrimes(7920);
		primesSieveSequencialDivision.savePrimesToFile("100primes.txt");

		BOOST_CHECK(primesSieveSequencialDivision.checkPrimesFromFile("./tests/testfiles/100.txt"));

	}
	BOOST_AUTO_TEST_SUITE_END()
