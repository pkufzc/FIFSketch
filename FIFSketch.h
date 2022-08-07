#pragma once

#include <iostream>
#include <fstream>
#include "hash.h"
#include "abstract.h"
#include "TowerCU.h"

using namespace std;

// TopK
template<uint32_t COUNTER_PER_BUCKET>
class FIFSketch : public Abstract{
public:
    FIFSketch(double _MEMORY)
        : Abstract((char *)"FIFSketch")
    {
        HEAVY_LENGTH = _MEMORY * HEAVY_RATIO / sizeof(Bucket);;
        
        buckets = new Bucket[HEAVY_LENGTH];
        memset(buckets, 0, sizeof(Bucket) * HEAVY_LENGTH);
        
        towerCU = new TowerCU(_MEMORY * LIGHT_RATIO);

    }

    ~FIFSketch(){
        delete [] buckets;
        delete towerCU;
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

        if(!(rand()%(minVal+1))){
            count_type light_query = towerCU->Query(item);
            towerCU->Insert(buckets[pos].ID[minPos], buckets[pos].count[minPos]);
            buckets[pos].ID[minPos] = item;
            buckets[pos].count[minPos] = light_query + 1;
        }
        else {
            towerCU->Insert(item);
        }
	}

    count_type Query(const data_type item){
        count_type result = buckets[hash(item) % HEAVY_LENGTH].Query(item);
        if(result!=0){
            return result;
        }else{
            return towerCU->Query(item);
        }
    }
    
    void get_distribution(vector<double> &dist){}
    double get_entropy(){return 0.0;}
    int get_cardinality(){return 0;}


private:
    double HEAVY_RATIO = 0.8;
    double LIGHT_RATIO = 0.2;
    
    struct Bucket{
        data_type ID[COUNTER_PER_BUCKET];
        count_type count[COUNTER_PER_BUCKET];

        count_type Query(const data_type item) {
            for(uint32_t i = 0; i < COUNTER_PER_BUCKET; i++) {
                if(ID[i] == item) {
                    return count[i];
                }
            }
            return 0;
        }
    };

    uint32_t HEAVY_LENGTH;

    Bucket* buckets;
    TowerCU* towerCU;

};