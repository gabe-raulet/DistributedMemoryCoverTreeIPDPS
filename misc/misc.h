#ifndef MISC_H_
#define MISC_H_

#include "mytypeinfo.h"

template <index_type Index>
Index read_integer(char *str)
{
    double x;
    char *p;

    x = strtod(str, &p);

    if      (toupper(*p) == 'K') x *= (1LL << 10);
    else if (toupper(*p) == 'M') x *= (1LL << 20);
    else if (toupper(*p) == 'G') x *= (1LL << 30);

    return static_cast<Index>(x + 0.499);
}

#endif
