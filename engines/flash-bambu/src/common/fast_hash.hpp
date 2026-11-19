#pragma once
#include <cstdint>
#include <cstring>
#include <vector>

namespace bambu_next {

inline uint64_t splitmix64(uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline uint64_t hash_combine(uint64_t seed, uint64_t val) noexcept {
    return seed ^ (splitmix64(val) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

inline uint64_t hash_int_pair(int32_t a, int32_t b) noexcept {
    uint64_t k = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
                 static_cast<uint64_t>(static_cast<uint32_t>(b));
    return splitmix64(k);
}

inline uint64_t hash_chain(int32_t seqid, int8_t strand,
                           const int32_t *starts, const int32_t *ends, int32_t n) noexcept {
    uint64_t h = hash_int_pair(seqid, static_cast<int32_t>(strand));
    h = hash_combine(h, static_cast<uint64_t>(n));
    for (int32_t i = 0; i < n; ++i) {
        h = hash_combine(h, hash_int_pair(starts[i], ends[i]));
    }
    return h;
}

} // namespace bambu_next
