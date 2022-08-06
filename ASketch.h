#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"
#include "CMSketch.h"

#define FILTER_SIZE 32

class ASketch : public Abstract{
public:
    ASketch(double _Memory) : Abstract((char *)"ASketch")
    {
        buckets = new Bucket[FILTER_SIZE];
        memset(buckets, 0, sizeof(Bucket) * FILTER_SIZE);
        sketch = new CMSketch(_Memory - sizeof(Bucket) * FILTER_SIZE);
    }

    ~ASketch(){
        delete [] buckets;
        delete sketch;
    }

    void Insert(const data_type key){
        for (uint32_t i = 0; i < FILTER_SIZE; i++)
        {
            if(buckets[i].ID == key){
                buckets[i].new_count++;
                return;
            }
            if(buckets[i].new_count == 0){
                buckets[i].ID = key;
                buckets[i].new_count = 1;
                buckets[i].old_count = 0;
                return;
            }
        }

        sketch->Insert(key);
        count_type estimated_value = sketch->Query(key);
        count_type min_value = 0xffffffff, min_index;
        for (uint32_t i = 0; i < FILTER_SIZE; i++)
        {
            if(buckets[i].new_count < min_value){
                min_index = i;
                min_value = buckets[i].new_count;
            }
        }

        if(estimated_value > min_value){
            count_type tmp = buckets[min_index].new_count - buckets[min_index].old_count;
            if(tmp > 0){
                sketch->Insert_with_f(buckets[min_index].ID, tmp);
            }
            buckets[min_index].ID = key;
            buckets[min_index].new_count = estimated_value;
            buckets[min_index].old_count = estimated_value;
        }
        
    }

    count_type Query(const data_type key){
        for (uint32_t i = 0; i < FILTER_SIZE; i++)
        {
            if(buckets[i].ID == key){
                return buckets[i].new_count;
            }
        }

        return sketch->Query(key);
    }

private:
    struct Bucket
    {
        data_type ID;
        count_type new_count;
        count_type old_count;
    };

    Bucket *buckets;
    CMSketch *sketch;
    
};