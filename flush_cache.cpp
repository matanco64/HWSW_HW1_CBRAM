#include <vector>
#include <cstdlib>
#include <iostream>

int main() {
    constexpr size_t size = 512 * 1024 * 1024; // 512MB
    std::vector<char> buffer(size);

    for (size_t i = 0; i < size; i += 64) {
        buffer[i] = i; // Touch every 64 bytes (typical cache line size)
    }

    std::cout << "Cache flushed by memory thrashing.\n";
    return 0;
}