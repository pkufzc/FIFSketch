#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"
#include "FCMSketch.h"

template<uint32_t COUNTER_PER_BUCKET>
class FCMSketchTopK : public Abstract{
public:
    FCMSketchTopK(double _MEMORY)
        : Abstract((char *)"FCMSketchTopK")
    {
        HEAVY_LENGTH = _MEMORY * HEAVY_RATIO / sizeof(Bucket);

        buckets = new Bucket[HEAVY_LENGTH];
        fcm = new FCMSketch(_MEMORY * LIGHT_RATIO);
        memset(buckets, 0, sizeof(Bucket) * HEAVY_LENGTH);
        
    }

    ~FCMSketchTopK(){
        delete [] buckets;
        delete fcm;
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
            return result + Light_Query(item);
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

    uint32_t HEAVY_LENGTH;

    FCMSketch* fcm;
    Bucket* buckets;

    void Light_Insert(const data_type key, count_type val = 1) { 
        for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
        {
            uint32_t pos = hash(key, i) % (64 * fcm->LENGTH);
            if (fcm->C1[i][pos] <= 254) {
                if (fcm->C1[i][pos] + val <= 254) {
                    fcm->C1[i][pos] = fcm->C1[i][pos] + val;
                    return;
                }
                else {
                    val = val - (254 - fcm->C1[i][pos]);
                    fcm->C1[i][pos] = 255;
                }
            }
            if (fcm->C1[i][pos] == 255) {
                pos = pos >> 3;
                if (fcm->C2[i][pos] <= 65534) {
                    if (fcm->C2[i][pos] + val <= 65534) {
                        fcm->C2[i][pos] += val;
                        return;
                    }
                    else {
                        val = val - (65534 - fcm->C2[i][pos]);
                        fcm->C2[i][pos] = 65535;
                    }
                }
                if (fcm->C2[i][pos] == 65535) {
                    pos = pos >> 3;
                    if (fcm->C3[i][pos] <= 0xffffffff - 1) {
                        fcm->C3[i][pos]=std::min(fcm->C3[i][pos] + val, 0xffffffff);
                    }
                }
            }
        }
    }

    count_type Light_Query(const data_type key) {
        count_type res = 0xffffffff;
        for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
        {
            uint32_t pos = hash(key, i) % (64 * fcm->LENGTH);
            count_type count_query=0;
            if (fcm->C1[i][pos] <= 254) {
                count_query = fcm->C1[i][pos];
            }
            else {
                pos = pos >> 3;
                if (fcm->C2[i][pos] <= 65534) {
                    count_query = 254 + fcm->C2[i][pos];
                }
                else {
                    pos = pos >> 3;
                    if (fcm->C3[i][pos] <= 0xffffffff - 1) {
                        count_query = 254 + 65534 + fcm->C3[i][pos];
                    }
                    else
                        count_query = 254+65534+ 0xffffffff - 1;
                }
            }
            res = std::min(res, count_query);
        }

        return res;

    }

};