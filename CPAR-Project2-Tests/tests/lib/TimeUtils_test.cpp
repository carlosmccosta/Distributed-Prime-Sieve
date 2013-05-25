#include <boost/test/unit_test.hpp>
#include "lib/TimeUtils.h"

BOOST_AUTO_TEST_SUITE(TimeUtils_test)
	BOOST_AUTO_TEST_CASE(formatSecondsToDate) {
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(0), "0s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(0.002), "0.002s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(1.002), "1.002s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(11.002), "11.002s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(59), "59s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(61), "01m01s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(70), "01m10s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(3599), "59m59s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(3661), "01h01m01s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(76881), "21h21m21s");
		BOOST_CHECK_EQUAL(TimeUtils::formatSecondsToDate(158401), "1d20h00m01s");
	}
	BOOST_AUTO_TEST_SUITE_END()
