#include <boost/test/unit_test.hpp>
#include "PrimesSieve.h"

BOOST_AUTO_TEST_SUITE(PrimesSieve_test)

	BOOST_AUTO_TEST_CASE(getNumberBitsToStore) {
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(3), 1);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(4), 1);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(5), 2);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(6), 2);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(7), 3);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(8), 3);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(9), 4);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(10), 4);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(11), 5);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(12), 5);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(13), 6);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberBitsToStore(14), 6);
	}

	BOOST_AUTO_TEST_CASE(getBitSetPositionToNumber) {
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(3), 0);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(4), 0);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(5), 1);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(6), 1);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(7), 2);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(8), 2);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(9), 3);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(10), 3);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(11), 4);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(12), 4);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(13), 5);
		BOOST_CHECK_EQUAL(PrimesSieve::getBitSetPositionToNumber(14), 5);
	}

	BOOST_AUTO_TEST_CASE(getNumberAssociatedWithBitSetPosition) {
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(0), 3);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(1), 5);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(2), 7);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(3), 9);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(4), 11);
		BOOST_CHECK_EQUAL(PrimesSieve::getNumberAssociatedWithBitSetPosition(5), 13);
	}
	BOOST_AUTO_TEST_SUITE_END()
