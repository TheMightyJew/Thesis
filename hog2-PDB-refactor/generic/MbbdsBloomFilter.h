#pragma once
#ifndef MbbdsBloomFilter_H
#define MbbdsBloomFilter_H

template<typename T, class H>
class MbbdsBloomFilter {
public:

	long long unsigned int M;
	unsigned int K;
	unsigned int hashOffset;
	std::vector<bool> filter;
	std::uint64_t onesCount;
	std::uint64_t count;
	
	MbbdsBloomFilter(unsigned int hashOffset = 0) : hashOffset(hashOffset), K(1), onesCount(0), count(0) {
	}
	MbbdsBloomFilter(long long unsigned int M, unsigned int K, unsigned int hashOffset = 0) : M(M), K(K), hashOffset(hashOffset), onesCount(0), count(0) {
		this->filter.resize(M, false);
	}

	std::vector<std::uint64_t> getCellsIDs(const T& object) const{
		std::vector<std::uint64_t> cellsIDs;
		for (unsigned i = this->hashOffset; i < this->hashOffset + this->K; i++) {
			cellsIDs.push_back(H::get(object, i*2) % this->M);
		}
		return cellsIDs;
	}

	void insert(const T& object) {
		std::vector<std::uint64_t> cellsIDs = getCellsIDs(object);
		for (std::uint64_t cellID : cellsIDs) {
			if(!this->filter[cellID]){
				filter[cellID] = true;
				this->onesCount++;
			}
		}
		++this->count;
	}

	bool contains(const T& object) const {
		std::vector<std::uint64_t> cellsIDs = getCellsIDs(object);
		for (std::uint64_t cellID : cellsIDs) {
			if (!this->filter[cellID]) {
				return false;
			}
		}
		return true;
	}
	
	void clear() {
		std::fill(this->filter.begin(), this->filter.end(), false);
		this->onesCount = 0;
		this->count = 0;
	}

	unsigned int getNextOffset(){
		return this->hashOffset + this->K;
	}

	double getSaturation() const {
		return double(this->onesCount) / double(this->M);
	}

	bool isSetEmpty() const {
		return this->count == 0;
	}

	std::uint64_t getCount() const {
		return this->count;
	}
	
	std::uint64_t getOnesCount() const {
		return this->onesCount;
	}
	
	long long unsigned int getM() const {
		return this->M;
	}
	
	unsigned int getK() const {
		return this->K;
	}
};

#endif

