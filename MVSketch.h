#pragma once

#include <iostream>
#include <fstream>
#include "hash.h"
#include "abstract.h"

using namespace std;

class MVSketch : public Abstract{
public:
    MVSketch(double _MEMORY): Abstract((char *)"MVSKetch")
    {
        w = _MEMORY / r / sizeof(Bucket);
        sketch = new Bucket*[r];
        for (uint32_t i = 0; i < r; i++)
        {
            sketch[i] = new Bucket[w];
            memset(sketch[i], 0, sizeof(Bucket) * w);
        }
        
    }

    ~MVSketch(){
        for (uint32_t i = 0; i < r; i++)
        {
            delete [] sketch[i];
        }
        delete [] sketch;
    }

    void Insert(const data_type key){
        for (uint32_t i = 0; i < r; i++)
        {
            uint32_t pos = hash(key, i) % w;
            sketch[i][pos].V += 1;
            if (sketch[i][pos].V == 1 && sketch[i][pos].C == 0){
                sketch[i][pos].K = key;
                sketch[i][pos].C = 1;
            }else if (sketch[i][pos].K == key){
                sketch[i][pos].C += 1;
            }else{
                sketch[i][pos].C -= 1;
                if (sketch[i][pos].C < 0)
                {
                    sketch[i][pos].K = key;
                    sketch[i][pos].C = -sketch[i][pos].C;
                }
            }
        }
    }

    count_type Query(const data_type key){
        count_type ret = 0xffffffff;
        for (uint32_t i = 0; i < r; i++)
        {
            uint32_t pos = hash(key, i) % w;
            if (sketch[i][pos].K == key)
            {
                ret = min(ret, (count_type)((sketch[i][pos].V + sketch[i][pos].C) / 2));
            }else{
                ret = min(ret, (count_type)((sketch[i][pos].V - sketch[i][pos].C) / 2));
            }
        }

        return ret;
    }

private:
    struct Bucket{
        data_type K;
        int32_t V;
        int32_t C;
    };

    uint32_t r = 4;
    uint32_t w;
    Bucket** sketch;
};