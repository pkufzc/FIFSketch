#pragma once

#include <iostream>

#define GetFourBit(x, pos) ((x >> (4*((pos+1)%2))) & 0xF)

using namespace std;

class TowerCU{
public:
    TowerCU(double _MEMORY){
        LENGTH_2bit = _MEMORY / 2 / sizeof(uint8_t) * 4;
        LENGTH_4bit = _MEMORY / 2 / sizeof(uint8_t) * 2;
        LENGTH = _MEMORY / 2 / sizeof(uint8_t);
        counters_2bit = new uint8_t[LENGTH];
        counters_4bit = new uint8_t[LENGTH];
        memset(counters_2bit, 0, LENGTH * sizeof(uint8_t));
        memset(counters_4bit, 0, LENGTH * sizeof(uint8_t));
    }

    ~TowerCU(){
        delete [] counters_2bit;
        delete [] counters_4bit;
    }

    void Insert(const data_type item, count_type val = 1){
        count_type minVal = 15, value, set_value;
        uint32_t counter_pos, index_in_counter;
        
        for(int i=0;i<2;i++){
            if(i==0){
                counter_pos = hash(item, i) % LENGTH_2bit / 4;
                index_in_counter = hash(item, i) % LENGTH_2bit % 4;
                value = GetTwoBit(counters_2bit[counter_pos],index_in_counter);
                if(value < 3 && value < minVal){
                    minVal = value;
                }
                
            }else{
                counter_pos = hash(item, i) % LENGTH_4bit / 2;
                index_in_counter = hash(item, i) % LENGTH_4bit % 2;
                value = GetFourBit(counters_4bit[counter_pos],index_in_counter);
                if(value < 15 && value < minVal){
                    minVal = value;
                }
                
            }
        }

        if(minVal==15) return;

        for(int i=0;i<2;i++){
            if(i==0){
                counter_pos = hash(item, 0) % LENGTH_2bit / 4;
                index_in_counter = hash(item, 0) % LENGTH_2bit % 4;
                value = GetTwoBit(counters_2bit[counter_pos],index_in_counter);
                set_value = (count_type)max((uint32_t)value, (uint32_t)(minVal + val));
                switch(index_in_counter){
                    case 0:
                        if(set_value >= 3){
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0x3f) + (3 << 6);
                        }else{
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0x3f) + (set_value << 6);
                        }
                        break;
                    case 1:
                        if(set_value >= 3){
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xcf) + (3 << 4);
                        }else{
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xcf) + (set_value << 4);
                        }
                        break;
                    case 2:
                        if(set_value >= 3){
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xf3) + (3 << 2);
                        }else{
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xf3) + (set_value << 2);
                        }
                        break;
                    case 3:
                        if(set_value >= 3){
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xfc) + 3;
                        }else{
                            counters_2bit[counter_pos] = (counters_2bit[counter_pos] & 0xfc) + set_value;
                        }
                        break;
                }
            }else{
                counter_pos = hash(item, 1) % LENGTH_4bit / 2;
                index_in_counter = hash(item, 1) % LENGTH_4bit % 2;
                value = GetFourBit(counters_4bit[counter_pos],index_in_counter);
                set_value = (count_type)max((uint32_t)value, (uint32_t)(minVal + val));
                switch(index_in_counter){
                    case 0:
                        if(set_value >=15){
                            counters_4bit[counter_pos] = (counters_4bit[counter_pos] & 0xf) + (15 << 4);
                        }else{
                            counters_4bit[counter_pos] = (counters_4bit[counter_pos] & 0xf) + (set_value << 4);
                        }
                        break;
                    case 1:
                        if(set_value >=15){
                            counters_4bit[counter_pos] = (counters_4bit[counter_pos] & 0xf0) + 15;
                        }else{
                            counters_4bit[counter_pos] = (counters_4bit[counter_pos] & 0xf0) + set_value;
                        }
                        break;
                }
               
            }
        }
        
    }

    count_type Query(const data_type item){
        count_type minv = 15;
       
        for(int i=0;i<2;i++){
            if(i==0){
                uint32_t counter_pos = hash(item, i) % LENGTH_2bit / 4;
                uint32_t index_in_counter = hash(item, i) % LENGTH_2bit % 4;
                count_type value = GetTwoBit(counters_2bit[counter_pos],index_in_counter);
                if(value==3){
                    continue;
                }else{
                    minv = std::min(minv,value);
                }
               
            }else{
                uint32_t counter_pos = hash(item, i) % LENGTH_4bit / 2;
                uint32_t index_in_counter = hash(item, i) % LENGTH_4bit % 2;
                count_type value = GetFourBit(counters_4bit[counter_pos],index_in_counter);
                if(value==15){
                    continue;
                }else{
                    minv = std::min(minv, value);
                }
               
            }
        }

        if(minv == 15) minv = 255;
       
        return minv;
    }
    

private:
    uint32_t LENGTH_2bit;
    uint32_t LENGTH_4bit;
    uint32_t LENGTH;

    uint8_t* counters_2bit;
    uint8_t* counters_4bit;

    count_type GetTwoBit(uint8_t x,int pos){
        if(pos==0){
            return (x >> 6) & 0x3;
        }else if(pos==1){
            return (x >> 4) & 0x3;
        }else if(pos==2){
            return (x >> 2) & 0x3;
        }else{
            return x & 0x3;
        }
    }
};