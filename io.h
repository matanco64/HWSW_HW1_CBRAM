#pragma once
#include <cstdint>
#include <string>

// Write sigma field as a PPM image (viridis-like colormap, sigma in Q16.16).
// sigma_max_q: the Q16.16 value that maps to the top of the colormap.
void write_ppm(const std::string& path, const int32_t* sigma,
               int N, int32_t sigma_max_q);

// Dump the V field as raw binary (N*N int32_t values, row-major).
void dump_binary(const std::string& path, const int32_t* V, int N);
