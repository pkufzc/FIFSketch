#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"

#define HASH_NUM 3

class CUSketch: public Abstract{
public:
    CUSketch(double _Memory) : Abstract((char *)"CUSketch")
    {
        LENGTH = _Memory / HASH_NUM / sizeof(count_type);
        counters = new count_type*[HASH_NUM];
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            counters[i] = new count_type[LENGTH];
            memset(counters[i], 0, sizeof(count_type) * LENGTH);
        }
    }

    ~CUSketch(){
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            delete [] counters[i];
        }
        delete [] counters;
    }

    void Insert(const data_type key){
        count_type minValue = 0xffffffff;
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            uint32_t pos = hash(key, i) % LENGTH;
            minValue = std::min(minValue, counters[i][pos]);
        }

        count_type tmp = minValue + 1;
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            uint32_t pos = hash(key, i) % LENGTH;
            counters[i][pos] = std::max(counters[i][pos], tmp);
        }
    }

    count_type Query(const data_type key){
        count_type res = 0xffffffff;
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            uint32_t pos = hash(key, i) % LENGTH;
            res = std::min(res, counters[i][pos]);
        }

        return res;
        
    }
    
private:
    count_type **counters;
    uint32_t LENGTH;
};