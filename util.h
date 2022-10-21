#ifndef UTIL_H
#define UTIL_H

#include <x86intrin.h>
#include <algorithm>
#include <string>
#include <random>
#include <unordered_map>
#include <climits>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <fstream>
#include "hash.h"

#define TUPLES_LEN 13

struct TIMESTAMP{
    uint8_t array[8];
};

struct TUPLES{
    uint8_t array[TUPLES_LEN];
};

bool operator == (const TUPLES& a, const TUPLES& b){
    return memcmp(a.array, b.array, sizeof(TUPLES)) == 0;
}

namespace std{
    template<>
    struct hash<TUPLES>{
        size_t operator()(const TUPLES& item) const noexcept
        {
            return Hash::BOBHash32((uint8_t*)&item, sizeof(TUPLES), 0);
        }
    };

    string Item (const TUPLES& x){
        string res = "";
        for(int i=0;i<TUPLES_LEN;i++){
            res += to_string(x.array[i]);
        }
        return res;
    }
}

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define hash hash32

typedef TUPLES data_type; // note: according to your dataset
typedef uint32_t count_type; // note: according to your dataset
typedef std::unordered_map< data_type, count_type > HashMap;

inline uint64_t hash64(data_type item, uint32_t seed = 0) {
	return Hash::BOBHash64((uint8_t *)&item, sizeof(data_type), seed);
}

inline uint32_t hash32(data_type item, uint32_t seed = 0) {
	return Hash::BOBHash32((uint8_t *)&item, sizeof(data_type), seed);
}

static std::random_device rd;

count_type Get_Median(count_type result[], uint32_t length) {
	std::sort(result, result + length);
	return (length & 1) ? result[length >> 1]
	                    : (result[length >> 1] + result[(length >> 1) - 1]) / 2;
}

#endif
