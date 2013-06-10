#pragma once

#include <stddef.h>

class PrimesUtils {
	public:
		inline static size_t closestPrimeMultiple(size_t prime, size_t number) {
			return number - (number % prime);
		}
};

