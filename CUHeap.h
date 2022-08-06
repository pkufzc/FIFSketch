#pragma once

#include <iostream>
#include "abstract.h"
#include "CUSketch.h"
#include "Heap.h"

class CUHeap : public Abstract
{
public:
    CUHeap(double _Memory) : Abstract((char *)"CUHeap")
    {
        sketch = new CUSketch(_Memory * SKETCH_RATIO);
        heap = new Heap(_Memory * HEAP_RATIO);
    }

    ~CUHeap(){
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

    CUSketch *sketch;
    Heap *heap;
};