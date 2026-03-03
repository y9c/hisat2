#ifndef CONVERSION_LOOKUP_H
#define CONVERSION_LOOKUP_H
#include <cstdint>
#include <cstring>
#include <cctype>
struct G3NTable {
    uint8_t table[256];
    void init(char f1, char f2, char f3, char f4) {
        std::memset(table, 0, 256);
        if (f1) { table[(uint8_t)f1] = 1; table[(uint8_t)std::tolower((unsigned char)f1)] = 1; }
        if (f2) { table[(uint8_t)f2] = 2; table[(uint8_t)std::tolower((unsigned char)f2)] = 2; }
        if (f3) { table[(uint8_t)f3] = 3; table[(uint8_t)std::tolower((unsigned char)f3)] = 3; }
        if (f4) { table[(uint8_t)f4] = 4; table[(uint8_t)std::tolower((unsigned char)f4)] = 4; }
    }
};
extern G3NTable g3NTable;
#endif
