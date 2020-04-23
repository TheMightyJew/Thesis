#pragma once
#ifndef MbbdsBloomFilter_H
#define MbbdsBloomFilter_H

template<typename T, class H>
class MbbdsBloomFilter {
public:

	long long unsigned int M;
	int K;
	unsigned hashOffset;
	std::vector<bool> filter;
	std::uint64_t onesCount;
	std::uint64_t count;
	
	MbbdsBloomFilter(unsigned hashOffset = 0) : onesCount(0), count(0) {
		this->hashOffset = hashOffset;
	}
	MbbdsBloomFilter(long long unsigned int M, int K, unsigned hashOffset = 0) : onesCount(0), count(0) {
		this->hashOffset = hashOffset;
		this->M = M;
		this->K = K;
		filter.resize(M, false);
	}

	void insert(const T& object) {
		for (unsigned i = hashOffset; i < hashOffset + K; ++i) {
			const std::uint64_t b = H::get(object, i*2) % M;
			onesCount += 1 - filter[b];
			filter[b] = true;
		}
		++count;
	}

	bool contains(const T& object) const {
		for (unsigned i = hashOffset; i < hashOffset + K; ++i) {
			const std::uint64_t b = H::get(object, i*2) % M;
			if (!filter[b]) {
				return false;
			}
		}
		return true;
	}
	
	void clear() {
		++hashOffset;
		std::fill(filter.begin(), filter.end(), false);
		onesCount = 0;
		count = 0;
	}

	double getSaturation() const {
		return double(onesCount) / double(M);
	}

	bool isSetEmpty() const {
		return count == 0;
	}

	std::uint64_t getCount() const {
		return count;
	}
	
	std::uint64_t getOnesCount() const {
		return onesCount;
	}
	
	long long unsigned int getM() const {
		return M;
	}
	
	int getK() const {
		return K;
	}
};

#endif

