#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"

template<uint32_t COUNTER_PER_BUCKET>
class Elastic : public Abstract{
public:
    Elastic(double _MEMORY)
        : Abstract((char *)"Elastic")
    {
        HEAVY_LENGTH = _MEMORY * HEAVY_RATIO / sizeof(Bucket);
        LIGHT_LENGTH = _MEMORY * LIGHT_RATIO / sizeof(uint8_t);

        buckets = new Bucket[HEAVY_LENGTH];
        counters = new uint8_t[LIGHT_LENGTH];
        memset(counters, 0, sizeof(uint8_t) * LIGHT_LENGTH);
        memset(buckets, 0, sizeof(Bucket) * HEAVY_LENGTH);
    }

    ~Elastic(){
        delete [] buckets;
        delete [] counters;
    }

    void Insert(const data_type item) {
        uint32_t pos = hash(item) % HEAVY_LENGTH, minPos = 0;
        count_type minVal = 0xffffffff;

        for (uint32_t i = 0; i < COUNTER_PER_BUCKET; i++){
            if(buckets[pos].ID[i] == item){
                buckets[pos].count[i] += 1;
                return;
            }

            if(buckets[pos].count[i] == 0){
                buckets[pos].ID[i] = item;
                buckets[pos].count[i] = 1;
                return;
            }

            if(buckets[pos].count[i] < minVal){
                minPos = i;
                minVal = buckets[pos].count[i];
            }
        }

        if((buckets[pos].vote + 1) >= minVal * LAMBDA){
            buckets[pos].vote = 0;
            buckets[pos].flags[minPos] = 1;

            Light_Insert(buckets[pos].ID[minPos], buckets[pos].count[minPos]);

            buckets[pos].ID[minPos] = item;
            buckets[pos].count[minPos] = 1;
        }
        else {
            buckets[pos].vote += 1;
            Light_Insert(item);
        }
	}

    count_type Query(const data_type item){
        uint8_t flag = 1;
        count_type result = buckets[hash(item) % HEAVY_LENGTH].Query(item, flag);
        if(flag)
            return result + (count_type)counters[hash(item, 101) % LIGHT_LENGTH];
        else
            return result;
    }

private:
    double HEAVY_RATIO = 0.25;
    double LIGHT_RATIO = 0.75;

    struct Bucket{
        count_type vote;
        uint8_t flags[COUNTER_PER_BUCKET];
        data_type ID[COUNTER_PER_BUCKET];
        count_type count[COUNTER_PER_BUCKET];

        count_type Query(const data_type item, uint8_t& flag) {
            for(uint32_t i = 0; i < COUNTER_PER_BUCKET; i++) {
                if(ID[i] == item) {
                    flag = flags[i];
                    return count[i];
                }
            }
            return 0;
        }
    };

    const uint32_t LAMBDA = 8;
    const count_type MAXNUM = 255;

    uint32_t LIGHT_LENGTH;
    uint32_t HEAVY_LENGTH;

    uint8_t* counters;
    Bucket* buckets;

    void Light_Insert(const data_type item, count_type val = 1) {
        uint32_t position = hash(item, 101) % LIGHT_LENGTH;
        count_type old_val = (count_type)counters[position];
        count_type new_val = (count_type)counters[position] + val;
        new_val = (uint8_t)(new_val > MAXNUM ? MAXNUM : new_val);
        counters[position] = new_val;
    }

};