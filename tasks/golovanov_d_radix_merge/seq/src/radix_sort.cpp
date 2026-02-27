#include "../include/radix_sort.hpp"

#include <cstdint>
#include <cstring>

void RadixSort::sort(std::vector<double>& arr) {
    const size_t n = arr.size();
    if (n == 0) return;

    std::vector<uint64_t> data(n);

    for (size_t i = 0; i < n; ++i) {
        uint64_t bits;
        std::memcpy(&bits, &arr[i], sizeof(double));

        if (bits >> 63) {              // отрицательное
            bits = ~bits;
        } else {                        // положительное
            bits ^= (1ULL << 63);
        }

        data[i] = bits;
    }

    std::vector<uint64_t> buffer(n);

    for (int byte = 0; byte < 8; ++byte) {
        size_t count[256] = {0};

        for (size_t i = 0; i < n; ++i) {
            uint8_t b = (data[i] >> (byte * 8)) & 0xFF;
            count[b]++;
        }

        size_t sum = 0;
        for (int i = 0; i < 256; ++i) {
            size_t tmp = count[i];
            count[i] = sum;
            sum += tmp;
        }

        for (size_t i = 0; i < n; ++i) {
            uint8_t b = (data[i] >> (byte * 8)) & 0xFF;
            buffer[count[b]++] = data[i];
        }

        data.swap(buffer);
    }

    for (size_t i = 0; i < n; ++i) {
        uint64_t bits = data[i];

        if (bits >> 63) {
            bits ^= (1ULL << 63);
        } else {
            bits = ~bits;
        }

        std::memcpy(&arr[i], &bits, sizeof(double));
    }
}