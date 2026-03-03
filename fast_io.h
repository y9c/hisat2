#ifndef FAST_IO_H
#define FAST_IO_H
#include <cstdint>
inline char* fast_append_int(char* buf, int n) {
    if (n == 0) { *buf = '0'; return buf + 1; }
    char* p = buf;
    if (n < 0) { *p++ = '-'; n = -n; }
    char* start = p;
    while (n) {
        *p++ = (n % 10) + '0';
        n /= 10;
    }
    char* end = p;
    p--;
    while (start < p) {
        char tmp = *start;
        *start++ = *p;
        *p-- = tmp;
    }
    return end;
}
#endif
