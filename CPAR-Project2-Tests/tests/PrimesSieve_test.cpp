#include <boost/test/unit_test.hpp>
#include "PrimesSieve.h"

BOOST_AUTO_TEST_SUITE(PrimesSieve_test)
	BOOST_AUTO_TEST_CASE(getNumberBitsToStore) {
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(3), 1);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(4), 1);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(5), 2);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(6), 2);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(7), 3);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(8), 3);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(9), 4);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(10), 4);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(11), 5);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(12), 5);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(13), 6);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberBitsToStore(14), 6);
	}
	
	BOOST_AUTO_TEST_CASE(getBitSetPositionToNumber) {
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(3), 0);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(4), 0);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(5), 1);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(6), 1);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(7), 2);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(8), 2);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(9), 3);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(10), 3);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(11), 4);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(12), 4);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(13), 5);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getBitsetPositionToNumber(14), 5);
	}
	
	BOOST_AUTO_TEST_CASE(getNumberAssociatedWithBitSetPosition) {
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(0), 3);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(1), 5);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(2), 7);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(3), 9);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(4), 11);
		BOOST_CHECK_EQUAL(PrimesSieve<vector<bool> >::getNumberAssociatedWithBitsetPosition(5), 13);
	}
	BOOST_AUTO_TEST_SUITE_END()
