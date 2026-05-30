#include "io.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

// Viridis colormap (25 entries, linearly interpolated)
static const uint8_t VIRIDIS[][3] = {
    {68,  1,  84}, {72,  26, 108}, {71,  47, 125}, {65,  68, 135},
    {57,  86, 140}, {49, 104, 142}, {43, 120, 142}, {38, 136, 141},
    {33, 152, 139}, {34, 167, 133}, {47, 182, 122}, {70, 196, 107},
    {101,207,  89}, {136,216,  70}, {172,220,  52}, {208,221,  38},
    {240,219,  33}, {253,210,  36}, {253,195,  37}, {252,180,  37},
    {252,160,  36}, {251,136,  29}, {244,109,  22}, {231,  82,  18},
    {253, 231,  37}
};
static constexpr int VIRIDIS_N = sizeof(VIRIDIS) / sizeof(VIRIDIS[0]);

static void viridis(float t, uint8_t& r, uint8_t& g, uint8_t& b) {
    t = std::max(0.0f, std::min(1.0f, t));
    float idx = t * (VIRIDIS_N - 1);
    int lo = (int)idx;
    int hi = std::min(lo + 1, VIRIDIS_N - 1);
    float frac = idx - lo;
    r = (uint8_t)(VIRIDIS[lo][0] + frac * (VIRIDIS[hi][0] - VIRIDIS[lo][0]));
    g = (uint8_t)(VIRIDIS[lo][1] + frac * (VIRIDIS[hi][1] - VIRIDIS[lo][1]));
    b = (uint8_t)(VIRIDIS[lo][2] + frac * (VIRIDIS[hi][2] - VIRIDIS[lo][2]));
}

void write_ppm(const std::string& path, const int32_t* sigma,
               int N, int32_t sigma_max_q)
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return;

    fprintf(f, "P6\n%d %d\n255\n", N, N);

    float scale = 1.0f / (float)sigma_max_q;
    uint8_t buf[3];
    for (int i = 0; i < N * N; ++i) {
        float t = (float)sigma[i] * scale;
        viridis(t, buf[0], buf[1], buf[2]);
        fwrite(buf, 1, 3, f);
    }
    fclose(f);
}

void dump_binary(const std::string& path, const int32_t* V, int N)
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return;
    fwrite(V, sizeof(int32_t), N * N, f);
    fclose(f);
}
