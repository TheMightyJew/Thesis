#ifndef GRID_HASHER_H
#define GRID_HASHER_H

#include "MurmurHash3.h"

struct GridHasher {
	static std::uint32_t get(const xyLoc& state, int n) {
		std::uint32_t out[2];
    int loc[2];
    loc[0] = state.x;
    loc[1] = state.y;
		MurmurHash3_x86_32(loc, 2*4, n * 12582917, out);

		return out[n % 2];
	}

	size_t operator()(const xyLoc& state) const {
		return get(state, 0);
	}
};
#endif