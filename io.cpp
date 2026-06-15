#include "io.h"
#include <cstdio>

void dump_binary(const std::string& path, const int32_t* V, int N)
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return;
    fwrite(V, sizeof(int32_t), N * N, f);
    fclose(f);
}
