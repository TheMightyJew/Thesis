#ifndef STP_HASHER_H
#define STP_HASHER_H

#include "MurmurHash3.h"
#include "MNPuzzle.h"

struct STPHasher {
	static std::uint32_t get(const MNPuzzleState<4, 4>& STPState, int n) {
		std::uint32_t out[2];
		MurmurHash3_x86_32(STPState.puzzle, 16*4, n * 12582917, out);

		return out[n % 2];
	}

	size_t operator()(const MNPuzzleState<4, 4>& STPState) const {
		return get(STPState, 0);
	}
};
#endif