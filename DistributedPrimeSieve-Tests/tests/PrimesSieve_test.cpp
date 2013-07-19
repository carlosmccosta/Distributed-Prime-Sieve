#include <boost/test/unit_test.hpp>
#include "PrimesSieve.h"
#include "sequencial/PrimesSieveSequencialMultiples.h"

BOOST_AUTO_TEST_SUITE(PrimesSieve_test)
	BOOST_AUTO_TEST_CASE(getNumberBitsToStore) {
		PrimesSieveSequencialMultiples<vector<bool> > primeSieve;
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(3), 1);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(4), 1);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(5), 2);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(6), 2);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(7), 3);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(8), 3);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(9), 4);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(10), 4);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(11), 5);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(12), 5);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(13), 6);
		BOOST_CHECK_EQUAL(primeSieve.getNumberBitsToStore(14), 6);
	}

	BOOST_AUTO_TEST_CASE(getBitSetPositionToNumber) {
		PrimesSieveSequencialMultiples<vector<bool> > primeSieve;
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(3), 0);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(4), 0);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(5), 1);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(6), 1);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(7), 2);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(8), 2);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(9), 3);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(10), 3);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(11), 4);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(12), 4);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(13), 5);
		BOOST_CHECK_EQUAL(primeSieve.getBitsetPositionToNumber(14), 5);
	}

	BOOST_AUTO_TEST_CASE(getNumberAssociatedWithBitSetPosition) {
		PrimesSieveSequencialMultiples<vector<bool> > primeSieve;
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(0), 3);
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(1), 5);
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(2), 7);
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(3), 9);
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(4), 11);
		BOOST_CHECK_EQUAL(primeSieve.getNumberAssociatedWithBitsetPosition(5), 13);
	}
	BOOST_AUTO_TEST_SUITE_END()
