#ifndef CONVERSION_LOOKUP_H
#define CONVERSION_LOOKUP_H
#include <cstdint>
#include <cstring>
#include <cctype>
struct G3NTable {
    uint8_t table[256];
    void init(char f, char t, char fc, char tc) {
        std::memset(table, 0, 256);
        table[(uint8_t)f] = 1;
        table[(uint8_t)fc] = 2;
        table[(uint8_t)std::tolower((unsigned char)f)] = 1;
        table[(uint8_t)std::tolower((unsigned char)fc)] = 2;
    }
};
extern G3NTable g3NTable;
#endif
