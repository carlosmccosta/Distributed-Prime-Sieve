#include "PrimesSieve.h"

bool PrimesSieve::checkComputedPrimes(const vector<size_t>& expectedPrimes) {
	if (primesValues.size() != expectedPrimes.size())
		return false;
	
	size_t iSize = primesValues.size();
	for (size_t i = 0; i < iSize; ++i) {
		if (primesValues[i] != expectedPrimes[i])
			return false;
	}
	
	return true;
}

bool PrimesSieve::checkPrimesFromFile(string filename) {
	ifstream inputStream(filename.c_str());
	
	if (inputStream.is_open()) {
		size_t numberRead;
		size_t iSize = primesValues.size();
		size_t i = 0;
		while (inputStream >> numberRead) {
			if (i >= iSize || numberRead != primesValues[i]) {
				return false;
			}
			++i;
		}
		
		return true;
	} else {
		cerr << "    -> File " << filename << " is not available!" << endl;
	}
	
	return false;
}

bool PrimesSieve::savePrimesToFile(string filename) {
	ofstream outputStream(filename.c_str());
	
	if (outputStream.is_open()) {
		size_t iSize = primesValues.size();
		for (size_t i = 0; i < iSize; ++i) {
			outputStream << primesValues[i] << endl;
		}
		
		return true;
	}
	
	return false;
}

vector<size_t>& PrimesSieve::extractPrimesFromBitset() {
	primesValues.clear();
	primesValues.push_back(2);
	size_t iSize = primesCompositesBitset.size();
	for (size_t i = 0; i < iSize; ++i) {
		if (!primesCompositesBitset[i]) {
			primesValues.push_back(getNumberAssociatedWithBitSetPosition(i));
		}
	}
	
	return primesValues;
}

void PrimesSieve::printPrimesToConsole() {
	size_t iSize = primesValues.size();
	for (size_t i = 0; i < iSize; ++i) {
		cout << primesValues[i] << endl;
	}
}

void PrimesSieve::initPrimesCompositeBitset(size_t maxRange) {
	this->maxRange = maxRange;
	primesCompositesBitset = vector<bool>(getNumberBitsToStore(maxRange));
	size_t iSize = primesCompositesBitset.size();
	for (size_t i = 0; i < iSize; ++i) {
		primesCompositesBitset[i] = false;
	}
}
