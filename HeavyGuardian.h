#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"

#define HEAVY_CELL_PER_BUCKET 8
#define LIGHT_CELL_PER_BUCKET 32
#define decay 1.08

class HeavyGuardian: public Abstract{
public:
    HeavyGuardian(double _Memory): Abstract((char *)"HeavyGuardian")
    {
        LENGTH = _Memory / sizeof(Bucket);
        buckets = new Bucket[LENGTH];
        memset(buckets, 0, sizeof(Bucket) * LENGTH);
    }

    ~HeavyGuardian(){
        delete [] buckets;
    }

    void Insert(const data_type key){
        uint32_t pos = hash(key) % LENGTH, minPos = 0;
        count_type minCount = 0xffffffff;

        for (uint32_t i = 0; i < HEAVY_CELL_PER_BUCKET; i++)
        {
            if (buckets[pos].heavy_part[i].ID == key)
            {
                buckets[pos].heavy_part[i].count++;
                return;
            }
            if (buckets[pos].heavy_part[i].count == 0)
            {
                buckets[pos].heavy_part[i].ID = key;
                buckets[pos].heavy_part[i].count = 1;
                return;
            }
            if (buckets[pos].heavy_part[i].count < minCount)
            {
                minPos = i;
                minCount = buckets[pos].heavy_part[i].count;
            }
        }

        if(!(rand()%int(pow(decay,minCount))))
        {
            buckets[pos].heavy_part[minPos].count--;
            if (buckets[pos].heavy_part[minPos].count == 0)
            {
                buckets[pos].heavy_part[minPos].ID = key;
                buckets[pos].heavy_part[minPos].count = 1;
            }else{
                uint8_t value = std::min(buckets[pos].light_part[hash(key) % LIGHT_CELL_PER_BUCKET] + 1, 255);
                buckets[pos].light_part[hash(key) % LIGHT_CELL_PER_BUCKET] = value;
            }
            
        }
        
    }

    count_type Query(const data_type key){
        uint32_t pos = hash(key) % LENGTH;
        for (uint32_t i = 0; i < HEAVY_CELL_PER_BUCKET; i++){
            if (buckets[pos].heavy_part[i].ID == key)
            {
                return buckets[pos].heavy_part[i].count;
            }
        }

        return buckets[pos].light_part[hash(key) % LIGHT_CELL_PER_BUCKET];
    }

private:
    struct HeavyCell
    {
        data_type ID;
        count_type count;
    };
    
    struct Bucket
    {
        HeavyCell heavy_part[HEAVY_CELL_PER_BUCKET];
        uint8_t light_part[LIGHT_CELL_PER_BUCKET];    
    };

    Bucket *buckets;

    uint32_t LENGTH;
    
};