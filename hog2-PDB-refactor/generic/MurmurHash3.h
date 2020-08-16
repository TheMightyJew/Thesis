//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

#ifndef MurmurHash3_h
#define MurmurHash3_h

#include <cinttypes>

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32  ( const void * key, int len, std::uint32_t seed, void * out );

void MurmurHash3_x86_128 ( const void * key, int len, std::uint32_t seed, void * out );

void MurmurHash3_x64_128 ( const void * key, int len, std::uint32_t seed, void * out );

//-----------------------------------------------------------------------------

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline std::uint32_t rotl32 ( std::uint32_t x, std::int8_t r )
{
  return (x << r) | (x >> (32 - r));
}

inline std::uint64_t rotl64 ( std::uint64_t x, std::int8_t r )
{
  return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)
#endif

