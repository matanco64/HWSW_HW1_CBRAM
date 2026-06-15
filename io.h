#pragma once
#include <cstdint>
#include <string>

// Dump the V field as raw binary (N*N int32_t values, row-major).
void dump_binary(const std::string& path, const int32_t* V, int N);
