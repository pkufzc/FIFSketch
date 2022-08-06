#pragma once

#include <iostream>
#include "abstract.h"
#include "hash.h"
//#include <bitset>
//#include <immintrin.h>
//#include <stdexcept>
//#include <cmath>
//#include <math.h>
//#include <iostream>
//#include <chrono>
//#include <random>
//#include <bitset>
//#include <x86intrin.h>
//using namespace std;

#define NUM_BITS 4
#define PHI 0.77351

class LogFilter {
private:
	const uint32_t threshold = 5;
	const uint32_t d = 3;
	uint32_t w;
	uint32_t f = 0;
	uint32_t *L;
	
public:
	uint32_t timestamp = 0;
	LogFilter(double _Memory) {
		w = _Memory * 8 / NUM_BITS;
		L = new uint32_t[w];//4-bit array
		memset(L, 0, sizeof(uint32_t) * w);
	}

	~LogFilter() {
		delete [] L;
	}

	bool probe(const data_type key) {
		uint32_t v = 0xffffffff;
		for (uint32_t i = 0; i < d; i++) {
			uint32_t value = L[hash(key, 500 + i) % w];
			v = value < v ? value : v;
		}
		return v < threshold;
	}
	
	void insert(const data_type key)
	{
		f++;
		for (uint32_t i = 0; i < d; i++) {
			//uint64_t tmp = ((uint64_t)key) << 32 + timestamp;
			//uint32_t random_num = hash(key, 1000 + i);//Hash::BOBHash64((uint8_t *)&tmp, sizeof(uint64_t), 1000 + i);
			char seeds[sizeof(key) + sizeof(timestamp)];
			memcpy(seeds, &key, sizeof(key));
			memcpy(seeds+sizeof(key), &timestamp, sizeof(timestamp));
			//data_type tmp = key;
			//uint32_t tmpt = timestamp & ((uint32_t*)(tmp.array))[0];
			//memcpy(tmp.array, &tmpt, sizeof(uint32_t));
			uint32_t random_num = (Hash::BOBHash32((const uint8_t *)&seeds, sizeof(seeds), 1000 + i));
			uint32_t temp = 0;
			//bitset<32> bittemp(random_num);
			for (uint32_t i = 31; i >= 0; i--) {
				if (!(random_num & (1<<i))) {
					temp++;
				}
				else
					break;
			}
			uint32_t index = hash(key, 500 + i) % w;
			L[index] = std::max(std::min(temp, threshold), L[index]);
		}
	}


	double query1(const data_type key)
	{
		double frequency = 0.0;
		double sum_value = 0.0;
		for (uint32_t i = 0; i < d; i++) {
			sum_value += L[hash(key, 500 + i) % w];
		}
		sum_value /= d;
		frequency = ((double)(w * d) / (w - d)) * (pow(2, sum_value) / PHI / d - (double)f / w);
		if (frequency > 0) {
			return frequency;
		}
		else {
			return (pow(2, sum_value) / PHI / d);
		}
	}

	double query2()
	{
		double frequency = 0.0;
		frequency = ((double)(w * d) / (w - d)) * (pow(2, threshold) / PHI / d - (double)f / w);
		if (frequency > 0) {
                        return frequency;
                }
                else {
                        return (pow(2, threshold) / PHI / d);
                }
	}
};
