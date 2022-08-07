#pragma once

#include <iostream>
#include "abstract.h"
#include "CMSketch.h"
#include "Heap.h"

class CMHeap : public Abstract
{
public:
    CMHeap(double _Memory) : Abstract((char *)"CMHeap")
    {
        sketch = new CMSketch(_Memory * SKETCH_RATIO);
        heap = new Heap(_Memory * HEAP_RATIO);
    }

    ~CMHeap(){
        delete sketch;
        delete heap;
    }

    void Insert(const data_type key){
        sketch->Insert(key);
        heap->Insert(key, sketch->Query(key));
    }

    count_type Query(const data_type key){
        count_type ret = heap->Query(key);
        if(!ret){
            return ret;
        }
        return sketch->Query(key);
    }


private:
    const double SKETCH_RATIO = 0.25;
    const double HEAP_RATIO = 0.75;

    CMSketch *sketch;
    Heap *heap;
};