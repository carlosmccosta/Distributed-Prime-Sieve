#include <boost/test/unit_test.hpp>
#include "lib/PrimesUtils.h"

BOOST_AUTO_TEST_SUITE(PrimesUtils_test)
	BOOST_AUTO_TEST_CASE(closestPrimeMultiple) {
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(3, 8), 6);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(3, 9), 9);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(3, 10), 9);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(3, 11), 9);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(3, 12), 12);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(5, 22), 20);


		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(7, 266), 266);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(7, 267), 266);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(7, 268), 266);
		BOOST_CHECK_EQUAL(PrimesUtils::closestPrimeMultiple(7, 269), 266);
	}
	BOOST_AUTO_TEST_SUITE_END()
