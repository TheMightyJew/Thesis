#ifndef PANCAKE_HASHER_H
#define PANCAKE_HASHER_H

#include "MurmurHash3.h"

template<int N>
struct PancakeHasher {
	static std::uint32_t get(const PancakePuzzleState<N>& pancakePuzzleState, int n) {
		std::uint32_t out[2];
		MurmurHash3_x86_32(pancakePuzzleState.puzzle, N*4, n * 12582917, out);

		return out[n % 2];
	}

	size_t operator()(const PancakePuzzleState<N>& pancakePuzzleState) const {
		return get(pancakePuzzleState, 0);
	}
};
#endif