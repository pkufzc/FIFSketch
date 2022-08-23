#pragma once

#include <iostream>
#include "hash.h"
#include "abstract.h"
#define FCMSK_DEPTH 2
//#define FCMSK_WL1 524288
//#define FCMSK_WL2 65536
//#define FCMSK_WL3 8192
typedef uint8_t FCMSK_C1;
typedef uint16_t FCMSK_C2;
typedef uint32_t FCMSK_C3;

class FCMSketch : public Abstract {
public:
	FCMSK_C1** C1;
	FCMSK_C2** C2;
	FCMSK_C3** C3;
	uint32_t LENGTH;

	FCMSketch(double _Memory) : Abstract((char*)"FCMSketch")
	{
		LENGTH = _Memory / FCMSK_DEPTH / 84 / sizeof(FCMSK_C1);
		C1 = new FCMSK_C1 * [FCMSK_DEPTH];
		C2 = new FCMSK_C2 * [FCMSK_DEPTH];
		C3 = new FCMSK_C3 * [FCMSK_DEPTH];
		for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
		{
			C1[i] = new FCMSK_C1[64 * LENGTH];
			memset(C1[i], 0, sizeof(FCMSK_C1) * 64 * LENGTH);
			C2[i] = new FCMSK_C2[8 * LENGTH];
			memset(C2[i], 0, sizeof(FCMSK_C2) * 8 * LENGTH);
			C3[i] = new FCMSK_C3[LENGTH];
			memset(C3[i], 0, sizeof(FCMSK_C3) * LENGTH);
		}
	}
	~FCMSketch() {
		for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
		{
			delete[] C1[i];
			delete[] C2[i];
			delete[] C3[i];
		}
		delete[] C1;
		delete[] C2;
		delete[] C3;
	}
	void Insert(const data_type key) {
		for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
		{
			uint32_t pos = hash(key, i) % (64 * LENGTH);
			if (C1[i][pos] <= 254) {
				C1[i][pos]++;
			}
			if (C1[i][pos] == 255) {
				pos = pos >> 3;
				if (C2[i][pos] <= 65534) {
					C2[i][pos]++;
				}
				if (C2[i][pos] == 65535) {
					pos = pos >> 3;
					if (C3[i][pos] <= 0xffffffff - 1) {
						C3[i][pos]++;
					}
				}
			}
		}
	}
	uint32_t Query(const data_type key) {
		uint32_t res = 0xffffffff;
		for (uint32_t i = 0; i < FCMSK_DEPTH; i++)
		{
			uint32_t pos = hash(key, i) % (64 * LENGTH);
			uint32_t count_query=0;
			if (C1[i][pos] <= 254) {
				count_query = C1[i][pos];
			}
			else {
				pos = pos >> 3;
				if (C2[i][pos] <= 65534) {
					count_query = 254 + C2[i][pos];
				}
				else {
					pos = pos >> 3;
					if (C3[i][pos] <= 0xffffffff - 1) {
						count_query = 254 + 65534 + C3[i][pos];
					}
					else
						count_query=254+65534+ 0xffffffff - 1;
				}
			}
			res = std::min(res, count_query);
		}

		return res;

	}
};