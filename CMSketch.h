#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"

#define HASH_NUM 3

class CMSketch: public Abstract{
public:
    CMSketch(double _Memory) : Abstract((char *)"CMSketch")
    {
        LENGTH = _Memory / HASH_NUM / sizeof(count_type);
        counters = new count_type*[HASH_NUM];
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            counters[i] = new count_type[LENGTH];
            memset(counters[i], 0, sizeof(count_type) * LENGTH);
        }
    }

    ~CMSketch(){
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            delete [] counters[i];
        }
        delete [] counters;
    }

    void Insert(const data_type key){
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            uint32_t pos = hash(key, i) % LENGTH;
            counters[i][pos]++;
        }
    }

    void Insert_with_f(const data_type key, count_type f=1){
        for (uint32_t i = 0; i < HASH_NUM; i++)
        {
            uint32_t pos = hash(key, i) % LENGTH;
            counters[i][pos] += f;
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