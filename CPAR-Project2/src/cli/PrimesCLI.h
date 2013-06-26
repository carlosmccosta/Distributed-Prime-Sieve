#pragma once

#include "../lib/ConsoleInput.h"
#include "../sequencial/PrimesSieveSequencialDivision.h"
#include "../sequencial/PrimesSieveSequencialMultiples.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceAndCache.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedSpaceTimeAndCache.h"
#include "../sequencial/PrimesSieveSequencialMultiplesOptimizedTimeAndCache.h"

#include <cmath>

class PrimesCLI {
	public:
		PrimesCLI(void);
		virtual ~PrimesCLI(void);
};
